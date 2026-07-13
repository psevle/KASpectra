#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

namespace kaspectra::math {

    namespace detail {

        // Count of panels that bailed out at max_depth WITHOUT meeting their
        // tolerance -- the signature of every silent-wrong-answer /
        // near-hang quadrature failure found in this project (a NaN-poisoned
        // integrand, a non-Lipschitz cancellation region, a singular edge).
        // thread_local so DNdEeTable's parallel build doesn't race; each
        // thread observes only its own integrations.
        inline long& max_depth_hit_counter() {
            thread_local long count = 0;
            return count;
        }

        inline double simpson(double fa, double fm, double fb, double a, double b) {
            return (b - a) / 6.0 * (fa + 4.0 * fm + fb);
        }

        // Halving rel_tol at each recursion level splits the overall error
        // budget across the growing number of subintervals (standard adaptive
        // Simpson technique) -- but with nothing floored against machine
        // epsilon, rel_tol shrinks geometrically without bound: by depth ~25
        // (2^-25 ~ 3e-8) it is already below DBL_EPSILON (~2.2e-16) once
        // multiplied into a large-magnitude `combined`. Flooring tol at a small
        // multiple of ULP(combined) stops the algorithm from chasing a
        // convergence target tighter than double precision can ever express
        // (see integrate_test.cpp's regression case for this mechanism in
        // isolation, with a synthetic integrand).
        //
        // This does NOT, by itself, bound recursion for every pathological
        // integrand: a real hang was found in bh::interaction_rate at an
        // extreme (but legal) parameter scale, where bh::psi_over_kappa2's
        // [2,4) branch sums O(10-300)-magnitude terms that must cancel almost
        // exactly near kappa==2, and that cancellation is non-Lipschitz at the
        // ULP level (confirmed directly: the returned value fluctuates
        // ~1e-14, sign-flipping, non-monotonic, as kappa creeps toward 2 by
        // successive ULPs). With this floor and the m<=a||m>=b degeneracy
        // bailout below both in place, that specific case still recursed past
        // 20,000,000 evaluations without resolving (confirmed by an
        // instrumented reproduction) -- the interval doesn't reach literal ULP
        // degeneracy until nearly max_depth=50 for that coordinate scale, so
        // the chaotic region gets ~30 full levels of wasted binary recursion
        // before either safeguard here can apply. The fix that actually
        // bounded that case was lowering bh/source.hpp's own max_depth default
        // (see its doc comment) -- both fixes here are still worthwhile
        // general hardening (they cost nothing for well-behaved integrands and
        // do help genuinely precision-limited or genuinely-degenerate cases),
        // just not a complete defense on their own against an integrand that's
        // chaotic well before the interval becomes degenerate.
        constexpr double kEpsilonFloorMultiplier = 4.0;

        // Noise-stagnation bail: a panel whose Richardson error estimate does
        // not SHRINK under refinement is integrating floating-point noise,
        // not structure -- splitting it further fills the entire 2^depth
        // recursion tree learning nothing (measured directly on bh::dN_dEe at
        // its spectral peak: the value is converged to ~4e-5 by depth 6, but
        // cost grows x8 per additional depth level, reaching CPU-hours at the
        // former default depths, because W()'s ~1e9-ULP cancellation noise
        // can never meet a tolerance that halves per level). Simpson on real
        // structure shrinks the estimate x16 per level and a kink confined to
        // one child shrinks it x2, both comfortably beating the factor-0.5
        // gate (the comparison is strict, so exactly-x2 kinks keep refining);
        // noise plateaus fluctuate around x1 and get culled -- a looser 0.75
        // gate measurably let noise pass often enough to keep the growth
        // exponential (~x2.7/level). Applied from depth 8 so early, coarse
        // estimates are never trusted, and counted in max_depth_hit_counter
        // (the tolerance was NOT certified) so it stays observable.
        constexpr int kStagnationMinDepth = 8;
        constexpr double kStagnationFactor = 0.5;

        inline double adaptive_simpson(const std::function<double(double)>& f,
                                        double a, double b, double fa, double fm, double fb,
                                        double whole, double abs_tol, double rel_tol, int depth, int max_depth,
                                        double prev_delta = std::numeric_limits<double>::infinity()) {
            double m = 0.5 * (a + b);

            // Bisection can no longer subdivide this interval in floating point
            // (m coincides with an endpoint) -- any further "recursion" would
            // just re-sample the same or ULP-adjacent points, wasting two more
            // evaluations to learn nothing. `whole` (the enclosing panel's own
            // single-Simpson estimate) is the best available answer here.
            if (m <= a || m >= b) {
                return whole;
            }

            double lm = 0.5 * (a + m);
            double rm = 0.5 * (m + b);
            double flm = f(lm);
            double frm = f(rm);

            double left = simpson(fa, flm, fm, a, m);
            double right = simpson(fm, frm, fb, m, b);
            double combined = left + right;

            const double eps_floor = kEpsilonFloorMultiplier * std::numeric_limits<double>::epsilon() * std::fabs(combined);
            double tol = std::max({abs_tol, rel_tol * std::fabs(combined), eps_floor});

            const double delta = std::fabs(combined - whole);
            const bool converged = delta <= 15.0 * tol;
            const bool stagnant = depth >= kStagnationMinDepth && delta > kStagnationFactor * prev_delta;
            if (depth >= max_depth || converged || stagnant) {
                if (!converged) ++max_depth_hit_counter();
                return combined + (combined - whole) / 15.0;
            }

            return adaptive_simpson(f, a, m, fa, flm, fm, left, abs_tol, rel_tol / 2.0, depth + 1, max_depth, delta) +
                    adaptive_simpson(f, m, b, fm, frm, fb, right, abs_tol, rel_tol / 2.0, depth + 1, max_depth, delta);
        }

    } // namespace detail

    // Diagnostics for the silent-failure class: a panel that exhausts
    // max_depth without converging returns its best estimate with NO error
    // signal, and every bad-numerics episode in this project (the 70-order
    // under-integration, the psi cancellation near-hang, the p_minus=0 NaN
    // grind) manifested exactly this way. Pattern:
    //
    //     math::reset_max_depth_hits();
    //     double v = expensive_integral(...);
    //     if (math::max_depth_hits() > 0) { /* v may be under-converged */ }
    //
    // Per-thread counters (see detail::max_depth_hit_counter), so check on
    // the thread that ran the integration. Zero for every well-behaved
    // integrand; a nonzero count is not necessarily wrong (the estimate may
    // still be adequate) but says the requested tolerance was NOT certified.
    inline long max_depth_hits() { return detail::max_depth_hit_counter(); }
    inline void reset_max_depth_hits() { detail::max_depth_hit_counter() = 0; }

    // A single-panel adaptive Simpson only samples 3 points (a, mid, b) before deciding
    // whether to recurse. If the integrand's entire nonzero support is a narrow sliver
    // that those initial samples miss (e.g. a threshold-peaked integrand evaluated over
    // a domain many decades wide), it can see ~0 everywhere it looks and falsely report
    // convergence, silently returning a result many orders of magnitude too small rather
    // than erroring. Pre-splitting into a fixed number of panels bounds how wide a region
    // any single panel's 3 initial samples must cover, so a narrow feature anywhere in
    // [a,b] is guaranteed to fall inside some panel's own initial sampling.
    constexpr int kDefaultPanels = 32;

    inline double integrate(const std::function<double(double)>& f, double a, double b,
                            double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50,
                            int panels = kDefaultPanels) {
        if (a == b) return 0.0;
        double sign = 1.0;
        if (a > b) { std::swap(a, b); sign = -1.0; }

        const double panel_abs_tol = abs_tol / panels;
        const double step = (b - a) / panels;
        double total = 0.0;
        double pa = a;
        for (int i = 0; i < panels; ++i) {
            double pb = (i == panels - 1) ? b : pa + step;
            double fa = f(pa), fb = f(pb), fm = f(0.5 * (pa + pb));
            double whole = detail::simpson(fa, fm, fb, pa, pb);
            total += detail::adaptive_simpson(f, pa, pb, fa, fm, fb, whole, panel_abs_tol, rel_tol, 0, max_depth);
            pa = pb;
        }
        return sign * total;
    }

    // u = log(E) substitution: integral f(E) dE == integral f(exp(u))*exp(u) du
    // E_p spans ~1 GeV to ~PeV
    inline double integrate_log(const std::function<double(double)>& f, double a, double b,
                                double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50,
                                int panels = kDefaultPanels) {
        auto g = [&f](double u) { double e = std::exp(u); return f(e) * e; };
        return integrate(g, std::log(a), std::log(b), abs_tol, rel_tol, max_depth, panels);
    }

// Usage
//   double result = integrate_log([](double e){ return sigma_inel(e) * J_p(e); },
//                                  1.0 /*GeV*/, 1e6 /*GeV*/);

}   // namespace kaspectra::math