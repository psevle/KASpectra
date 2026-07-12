#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

#include "kaspectra/bh/spectrum.hpp"
#include "kaspectra/constants.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/math/integrate.hpp"
#include "kaspectra/math/interp1d.hpp"

// Additive caching/interpolation layer over bh::dN_dEe. Does not touch
// spectrum.hpp -- dN_dEe and q_pair_spectrum remain available exactly as
// before for callers who want the slow-but-exact answer. Two supported use
// patterns:
//   1) DNdEeTable::over_gamma_p  -- fixed E_e, varying gamma_p
//      (feeds q_pair_spectrum_cached's outer integral)
//   2) DNdEeTable::over_E_e      -- fixed gamma_p, varying E_e
//      (feeds CLI-style sweeps, e.g. examples/bh_spectrum.cpp)
//   3) DNdEeTable2D::build       -- varying BOTH axes
//      (one build amortized across an entire SED sweep of
//      q_pair_spectrum_cached calls; see the class doc below)
//
// Honest cost/benefit: building an N-point table costs N dN_dEe calls
// up front. A single q_pair_spectrum_cached call is roughly a wash against
// plain q_pair_spectrum (its outer quadrature already calls dN_dEe on the
// same order of times as a table build) -- the real win is REPEATED calls
// against one table (e.g. scanning J_p's parameters, or E_p_max cutoffs,
// since dN_dEe never depends on J_p). For an E_e sweep (case 2), the win
// only kicks in once the number of desired sample points exceeds the
// table's build cost (e.g. an example plotting hundreds of points instead
// of a couple dozen).
//
// Usage:
//   io::BlackbodyPhotonField cmb(2.725);
//   auto table = bh::DNdEeTable::over_E_e(1e9 /*gamma_p*/, cmb, 1e-6);
//   double y = table(3.5e5); // cheap lookup, no nested integral

namespace kaspectra::bh {

    // E_e-axis kinematic window at FIXED gamma_p: E_p_min (spectrum.hpp)
    // inverts eps_lo_bound(gamma_p,E_e,m_e)==epsilon_max for gamma_p at fixed
    // E_e (single threshold, eps_lo_bound monotonically decreasing in
    // gamma_p). This is the analogous inversion for E_e at fixed gamma_p --
    // but eps_lo_bound is U-SHAPED in E_e, minimized exactly at
    // E_e = gamma_p*m_e (the spectrum's own peak location) with floor value
    // m_e/gamma_p (set d(eps_lo_bound)/dE_e=0 and solve: g(E_e) =
    // (E_e+a)^2/E_e = E_e+2a+a^2/E_e with a=gamma_p*m_e has g'(E_e)=0 at
    // E_e=a, giving g(a)=4a, eps_lo_min=4a/(4*gamma_p^2)=m_e/gamma_p).
    // So for fixed gamma_p there are generically TWO roots bracketing a
    // finite window [E_e_lo, E_e_hi], from the quadratic
    // E_e^2 + B*E_e + C = 0 with
    //   B = 2*gamma_p*m_e - 4*gamma_p^2*epsilon_max,  C = gamma_p^2*m_e^2.
    struct EeWindow {
        double lo = 0.0;
        double hi = 0.0;
        bool reachable = false;
    };

    inline EeWindow E_e_window(double gamma_p, double epsilon_max, double m_e) {
        const double floor = m_e / gamma_p;
        if (epsilon_max <= floor) return EeWindow{0.0, 0.0, false};

        const double B = 2.0 * gamma_p * m_e - 4.0 * gamma_p * gamma_p * epsilon_max;
        const double C = gamma_p * gamma_p * m_e * m_e;
        const double disc = B * B - 4.0 * C;
        if (disc < 0.0) return EeWindow{0.0, 0.0, false};   // shouldn't happen given the floor check above

        // Stable quadratic roots (Numerical Recipes form) -- B can have
        // either sign on this axis, unlike E_p_min's single-sign shortcut.
        const double sq = std::sqrt(disc);
        const double q = (B >= 0.0) ? -0.5 * (B + sq) : -0.5 * (B - sq);
        const double x1 = q;
        const double x2 = C / q;

        return EeWindow{std::min(x1, x2), std::max(x1, x2), true};
    }

    namespace detail {

        // Benchmarked directly (48-ppd dense dN_dEe reference over a 3-decade
        // band crossing the spectral peak, gamma_p=1e6, warm blackbody;
        // coarser grids simulated by subsampling): max/median interpolation
        // error vs points-per-decade was 3.5%/0.8% at 4, 2.9%/0.2% at 8,
        // 2.5%/0.09% at 12, 1.6%/0.03% at 24. Max error is dominated by
        // curvature near the peak and band edges and barely improves with
        // density, while build cost scales linearly with it -- 12/decade
        // matches the few-percent intrinsic accuracy of dN_dEe itself (see
        // spectrum.hpp's tolerance notes) at half the cost of the previous
        // default of 24.
        constexpr int kDefaultPointsPerDecade = 12;

        // Analytic zero boundaries are backed off by this relative margin
        // before the first/last grid node is placed, matching the existing
        // 1e-9 backoff used for E_minus_hi in spectrum.hpp: drops only a
        // measure-zero sliver right at the kinematic edge, not real accuracy.
        constexpr double kBoundaryMargin = 1e-9;

        // Sentinel for grid nodes where dN_dEe underflows to exact 0.0 despite
        // being kinematically reachable (deep photon-field Wien-tail
        // suppression) -- keeps Interpolator1D's LogLog y>0 precondition
        // satisfied without perturbing anything physically meaningful, since
        // this only occurs where the true value is already negligible.
        constexpr double kYFloor = 1e-300;

    }   // namespace detail

    // Grid density heuristic tied to domain span rather than a fixed count,
    // so callers don't have to guess how many points a given (lo,hi) range
    // needs.
    inline int suggested_n_points(double x_lo, double x_hi,
                                    int points_per_decade = detail::kDefaultPointsPerDecade) {
        if (!(x_hi > x_lo) || !(x_lo > 0.0)) return 4;
        const double decades = std::log10(x_hi / x_lo);
        return std::max(4, static_cast<int>(std::ceil(decades * points_per_decade)) + 1);
    }

    enum class CacheAxis { GammaP, ElectronEnergy };

    // Cached/interpolated stand-in for dN_dEe, valid only for the fixed
    // (f_ph, epsilon_max, and the OTHER axis's fixed value) it was built
    // with. Immutable after construction (no mutable state), so a single
    // instance is safe to share read-only across threads.
    class DNdEeTable {
        public:
            // Default-constructed table is "empty": reachable()==false,
            // operator() always returns 0.0. Only factories below produce a
            // populated table.
            DNdEeTable() = default;

            // Build over gamma_p at FIXED E_e (feeds q_pair_spectrum_cached's
            // outer integral). gamma_p_max is the only bound the caller must
            // supply (typically E_p_max/m_p) -- the lower bound is derived
            // from E_p_min, which is exact.
            static DNdEeTable over_gamma_p(double E_e, const io::PhotonField& f_ph, double epsilon_max,
                                            double gamma_p_max, int n_points = 0,
                                            double abs_tol = 1e-25, double rel_tol = 1e-4,
                                            int max_depth = 20, int panels = 24) {
                using namespace kaspectra::constants;

                DNdEeTable t;
                t.axis_ = CacheAxis::GammaP;
                t.fixed_ = E_e;
                t.two_sided_ = false;

                const double E_p_lo = E_p_min(E_e, epsilon_max);
                if (!std::isfinite(E_p_lo)) return t;   // never reachable at any gamma_p

                const double gamma_p_lo_exact = E_p_lo / m_p;
                const double lo = gamma_p_lo_exact * (1.0 + detail::kBoundaryMargin);
                const double hi = gamma_p_max;
                if (!(lo < hi)) return t;   // requested gamma_p_max doesn't reach threshold

                t.build(lo, hi,
                        [&](double gamma_p) { return dN_dEe(E_e, gamma_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels); },
                        n_points);
                return t;
            }

            // Build over E_e at FIXED gamma_p (feeds CLI-style sweeps).
            // E_e_min/E_e_max default to the FULL analytic reachable window
            // (E_e_window) -- pass narrower values only to clip a sweep to a
            // sub-range you already know you care about.
            static DNdEeTable over_E_e(double gamma_p, const io::PhotonField& f_ph, double epsilon_max,
                                        double E_e_min = 0.0,
                                        double E_e_max = std::numeric_limits<double>::infinity(),
                                        int n_points = 0,
                                        double abs_tol = 1e-25, double rel_tol = 1e-4,
                                        int max_depth = 20, int panels = 24) {
                using namespace kaspectra::constants;

                DNdEeTable t;
                t.axis_ = CacheAxis::ElectronEnergy;
                t.fixed_ = gamma_p;
                t.two_sided_ = true;

                const EeWindow w = E_e_window(gamma_p, epsilon_max, m_e);
                if (!w.reachable) return t;

                const double lo = std::max(E_e_min, w.lo * (1.0 + detail::kBoundaryMargin));
                const double hi = std::min(E_e_max, w.hi * (1.0 - detail::kBoundaryMargin));
                if (!(lo < hi)) return t;   // requested [E_e_min,E_e_max] misses the reachable window

                t.build(lo, hi,
                        [&](double E_e) { return dN_dEe(E_e, gamma_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels); },
                        n_points);
                return t;
            }

            double operator()(double x) const {
                if (!reachable_) return 0.0;
                if (x <= lo_) return 0.0;
                if (two_sided_ && x >= hi_) return 0.0;
                return (*interp_)(x);   // single-sided axis: beyond hi_ falls through to
                                        // Interpolator1D's own clamp-to-boundary behavior
            }

            bool reachable() const { return reachable_; }
            CacheAxis axis() const { return axis_; }
            double fixed_value() const { return fixed_; }
            double domain_lo() const { return lo_; }
            double domain_hi() const { return hi_; }
            std::size_t size() const { return n_; }

        private:
            template <typename F>
            void build(double lo, double hi, F&& f, int n_points) {
                lo_ = lo;
                hi_ = hi;
                reachable_ = true;

                int n = n_points > 0 ? n_points : suggested_n_points(lo, hi);
                n = std::max(n, 4);
                n_ = static_cast<std::size_t>(n);

                std::vector<double> xs(n_), ys(n_);
                const double log_lo = std::log(lo), log_hi = std::log(hi);
                for (std::size_t i = 0; i < n_; ++i) {
                    const double frac = static_cast<double>(i) / static_cast<double>(n_ - 1);
                    xs[i] = std::exp(log_lo + frac * (log_hi - log_lo));
                }

                // Grid points are independent dN_dEe evaluations (seconds each at
                // UHECR scale), each writing only its own slot -- embarrassingly
                // parallel. f is required to be thread-safe: every PhotonField/
                // ProtonSpectrum implementation in io/ is immutable after
                // construction and dN_dEe itself has no shared state. std::async
                // (not raw std::thread) so an exception in any chunk propagates
                // out of get() instead of terminating. Results are bit-identical
                // to the serial loop: same xs, same per-point arithmetic, no
                // reduction order to vary.
                const std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
                const std::size_t n_chunks = std::min(hw, n_);
                std::vector<std::future<void>> chunks;
                chunks.reserve(n_chunks);
                for (std::size_t c = 0; c < n_chunks; ++c) {
                    const std::size_t i_lo = c * n_ / n_chunks;
                    const std::size_t i_hi = (c + 1) * n_ / n_chunks;
                    chunks.push_back(std::async(std::launch::async, [&, i_lo, i_hi] {
                        for (std::size_t i = i_lo; i < i_hi; ++i) {
                            ys[i] = std::max(f(xs[i]), detail::kYFloor);
                        }
                    }));
                }
                for (auto& ch : chunks) ch.get();

                interp_ = std::make_unique<math::Interpolator1D>(std::move(xs), std::move(ys), math::InterpMode::LogLog);
            }

            CacheAxis axis_ = CacheAxis::GammaP;
            double fixed_ = 0.0;
            bool reachable_ = false;
            bool two_sided_ = false;
            double lo_ = 0.0, hi_ = 0.0;
            std::size_t n_ = 0;
            std::unique_ptr<math::Interpolator1D> interp_;
    };

    // Convenience: gamma_p_max derived from the same E_p_max callers already
    // pass to q_pair_spectrum, so building a table needs no new bound to guess.
    inline DNdEeTable make_table_for_q_pair_spectrum(double E_e, const io::PhotonField& f_ph,
                                                        double epsilon_max, double E_p_max,
                                                        int n_points = 0) {
        return DNdEeTable::over_gamma_p(E_e, f_ph, epsilon_max, E_p_max / kaspectra::constants::m_p, n_points);
    }

    // Sibling of q_pair_spectrum (spectrum.hpp) that consults a pre-built
    // DNdEeTable instead of calling dN_dEe directly at every outer quadrature
    // node. table MUST have been built via DNdEeTable::over_gamma_p (or
    // make_table_for_q_pair_spectrum) for this SAME E_e/epsilon_max/f_ph --
    // the checks below catch an axis or E_e mismatch, and an E_p_max that
    // would query beyond the table's built range (both would otherwise
    // silently return a stale clamped/wrong value rather than erroring, the
    // same silent-wrong-answer failure class this project has hit before in
    // math::integrate). They cannot detect a table built against a different
    // f_ph or epsilon_max (not stored on the table, since PhotonField has no
    // equality comparison) -- caller is responsible for keeping those
    // consistent, matching this codebase's "cache is scoped to one
    // field+epsilon_max combination" convention.
    //
    // Real win is REPEATED calls against the same table (varying J_p's own
    // parameters, since dN_dEe never depends on J_p, or re-running at
    // different rel_tol/max_depth for the outer integral) -- a single call
    // is only as fast as the table build was cheap relative to the outer
    // quadrature's own node count (see this file's top-of-file cost/benefit
    // note).
    inline double q_pair_spectrum_cached(double E_e, const io::ProtonSpectrum& J_p, const DNdEeTable& table,
                                            double E_p_max, double epsilon_max,
                                            double abs_tol = 1e-25, double rel_tol = 1e-3,
                                            int max_depth = 10, int panels = 6) {
        using namespace kaspectra::constants;

        if (table.axis() != CacheAxis::GammaP) {
            throw std::invalid_argument("q_pair_spectrum_cached: table must be built via DNdEeTable::over_gamma_p");
        }
        if (std::fabs(table.fixed_value() - E_e) > 1e-9 * std::fabs(E_e)) {
            throw std::invalid_argument("q_pair_spectrum_cached: table was built for a different E_e");
        }

        const double E_p_lo = E_p_min(E_e, epsilon_max);
        if (E_p_lo >= E_p_max) return 0.0;

        // Genuinely reachable per the analytic threshold, but the table
        // itself failed to build (e.g. built earlier with a smaller
        // gamma_p_max that didn't clear E_p_min) -- without this check the
        // integrand below would silently return 0.0 everywhere instead of
        // the true nonzero contribution.
        if (!table.reachable()) {
            throw std::invalid_argument("q_pair_spectrum_cached: table failed to build (unreachable) but "
                                         "[E_p_min(E_e,epsilon_max), E_p_max] is a nonempty window -- rebuild "
                                         "the table with a gamma_p_max covering E_p_max/m_p");
        }

        if (E_p_max / m_p > table.domain_hi() * (1.0 + 1e-6)) {
            throw std::invalid_argument("q_pair_spectrum_cached: E_p_max exceeds the table's built gamma_p range "
                                         "-- rebuild via over_gamma_p/make_table_for_q_pair_spectrum with a "
                                         "gamma_p_max covering E_p_max/m_p");
        }

        auto integrand = [&](double E_p) {
            return J_p(E_p) * table(E_p / m_p);
        };

        return math::integrate_log(integrand, E_p_lo, E_p_max, abs_tol, rel_tol, max_depth, panels);
    }

    // 2D cache over (E_e, gamma_p): a "table of tables" -- one over_gamma_p
    // line per log-spaced E_e node, log-interpolated between adjacent lines at
    // query time. Chosen over a flat rectangular grid because the reachable
    // gamma_p window depends on E_e: each line snaps to ITS OWN exact
    // E_p_min(E_e) threshold, so no grid cell ever blends real values with
    // out-of-window zeros along the gamma_p axis. Along the E_e axis, a node
    // that is unreachable outright (E_p_min infinite, or above gamma_p_max)
    // yields an empty line; interpolation against an empty neighbour falls
    // back from log-log to linear so the value decays continuously to 0
    // instead of blowing up on log(0).
    //
    // Amortization: ONE build serves an entire SED sweep -- every
    // q_pair_spectrum_cached(E_e, ...) call for E_e inside the built band, at
    // any J_p and any E_p_max up to the built bound. Build cost is
    // n_lines * (points per line) dN_dEe calls (each line's evaluation loop is
    // parallel, see DNdEeTable::build).
    class DNdEeTable2D {
        public:
            DNdEeTable2D() = default;

            // Explicitly move-only: lines_ holds move-only DNdEeTable entries,
            // but std::vector declares a copy constructor unconditionally, so
            // is_copy_constructible<DNdEeTable2D> would (wrongly) report true
            // and pybind11/generic code would hard-error inside the vector on
            // instantiation instead of falling back to the move.
            DNdEeTable2D(const DNdEeTable2D&) = delete;
            DNdEeTable2D& operator=(const DNdEeTable2D&) = delete;
            DNdEeTable2D(DNdEeTable2D&&) = default;
            DNdEeTable2D& operator=(DNdEeTable2D&&) = default;

            // E_e band [E_e_min, E_e_max] and proton budget E_p_max are the
            // same numbers an SED sweep already has. epsilon_max is STORED and
            // checked by q_pair_spectrum_cached (unlike the 1D table, whose
            // guards cannot catch it); f_ph remains uncheckable (no equality
            // comparison on PhotonField) -- caller keeps that consistent.
            static DNdEeTable2D build(const io::PhotonField& f_ph, double epsilon_max,
                                        double E_p_max, double E_e_min, double E_e_max,
                                        int n_E_e_lines = 0, int n_points_per_line = 0,
                                        double abs_tol = 1e-25, double rel_tol = 1e-4,
                                        int max_depth = 20, int panels = 24) {
                DNdEeTable2D t;
                t.epsilon_max_ = epsilon_max;
                t.gamma_p_max_ = E_p_max / kaspectra::constants::m_p;
                if (!(E_e_min > 0.0) || !(E_e_min < E_e_max)) return t;

                const int n = n_E_e_lines > 0 ? n_E_e_lines : suggested_n_points(E_e_min, E_e_max);
                t.E_e_nodes_.resize(static_cast<std::size_t>(std::max(n, 2)));
                t.lines_.resize(t.E_e_nodes_.size());

                const double log_lo = std::log(E_e_min), log_hi = std::log(E_e_max);
                for (std::size_t i = 0; i < t.E_e_nodes_.size(); ++i) {
                    const double frac = static_cast<double>(i) / static_cast<double>(t.E_e_nodes_.size() - 1);
                    const double E_e = std::exp(log_lo + frac * (log_hi - log_lo));
                    t.E_e_nodes_[i] = E_e;
                    t.lines_[i] = DNdEeTable::over_gamma_p(E_e, f_ph, epsilon_max, t.gamma_p_max_,
                                                            n_points_per_line, abs_tol, rel_tol,
                                                            max_depth, panels);
                }
                t.built_ = true;
                return t;
            }

            // 0.0 outside the built E_e band (and wherever both bracketing
            // lines are 0/unreachable) -- same "0 beyond domain" convention as
            // the 1D table's two-sided axis.
            double operator()(double E_e, double gamma_p) const {
                if (!built_) return 0.0;
                if (E_e < E_e_nodes_.front() || E_e > E_e_nodes_.back()) return 0.0;

                const auto it = std::upper_bound(E_e_nodes_.begin(), E_e_nodes_.end(), E_e);
                const std::size_t hi = std::min<std::size_t>(
                    static_cast<std::size_t>(it - E_e_nodes_.begin()), E_e_nodes_.size() - 1);
                const std::size_t lo = hi == 0 ? 0 : hi - 1;
                if (lo == hi) return lines_[lo](gamma_p);

                const double y0 = lines_[lo](gamma_p);
                const double y1 = lines_[hi](gamma_p);
                const double frac = (std::log(E_e) - std::log(E_e_nodes_[lo])) /
                                    (std::log(E_e_nodes_[hi]) - std::log(E_e_nodes_[lo]));
                if (y0 > 0.0 && y1 > 0.0) {
                    return std::exp((1.0 - frac) * std::log(y0) + frac * std::log(y1));
                }
                return (1.0 - frac) * y0 + frac * y1;   // near a zero boundary: linear, decays to 0
            }

            bool built() const { return built_; }
            double epsilon_max() const { return epsilon_max_; }
            double gamma_p_max() const { return gamma_p_max_; }
            double E_e_min() const { return built_ ? E_e_nodes_.front() : 0.0; }
            double E_e_max() const { return built_ ? E_e_nodes_.back() : 0.0; }
            std::size_t n_lines() const { return lines_.size(); }

        private:
            bool built_ = false;
            double epsilon_max_ = 0.0;
            double gamma_p_max_ = 0.0;
            std::vector<double> E_e_nodes_;
            std::vector<DNdEeTable> lines_;
    };

    // SED-sweep sibling of q_pair_spectrum_cached: same outer integral, but
    // consulting one shared 2D table for every E_e in the built band instead
    // of needing a fresh 1D table per E_e.
    inline double q_pair_spectrum_cached(double E_e, const io::ProtonSpectrum& J_p, const DNdEeTable2D& table,
                                            double E_p_max, double epsilon_max,
                                            double abs_tol = 1e-25, double rel_tol = 1e-3,
                                            int max_depth = 10, int panels = 6) {
        using namespace kaspectra::constants;

        if (!table.built()) {
            throw std::invalid_argument("q_pair_spectrum_cached: 2D table was never built");
        }
        if (std::fabs(table.epsilon_max() - epsilon_max) > 1e-9 * epsilon_max) {
            throw std::invalid_argument("q_pair_spectrum_cached: 2D table was built for a different epsilon_max");
        }
        if (E_e < table.E_e_min() * (1.0 - 1e-9) || E_e > table.E_e_max() * (1.0 + 1e-9)) {
            throw std::invalid_argument("q_pair_spectrum_cached: E_e outside the 2D table's built band");
        }
        if (E_p_max / m_p > table.gamma_p_max() * (1.0 + 1e-6)) {
            throw std::invalid_argument("q_pair_spectrum_cached: E_p_max exceeds the 2D table's built gamma_p range");
        }

        const double E_p_lo = E_p_min(E_e, epsilon_max);
        if (E_p_lo >= E_p_max) return 0.0;

        auto integrand = [&](double E_p) {
            return J_p(E_p) * table(E_e, E_p / m_p);
        };

        return math::integrate_log(integrand, E_p_lo, E_p_max, abs_tol, rel_tol, max_depth, panels);
    }

}   // namespace kaspectra::bh
