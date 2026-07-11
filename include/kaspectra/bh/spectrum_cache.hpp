#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
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

        // dN_dEe's peak spans ~3 orders of magnitude within ~1.5 decades of
        // its varying axis; 24 points/decade puts ~36 points across that
        // region, keeping log-log-linear interpolation error to a few percent
        // there (see spectrum_cache_test.cpp's table-vs-direct check).
        constexpr int kDefaultPointsPerDecade = 24;

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
                    const double x = std::exp(log_lo + frac * (log_hi - log_lo));
                    const double y = f(x);
                    xs[i] = x;
                    ys[i] = std::max(y, detail::kYFloor);
                }
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

}   // namespace kaspectra::bh
