#pragma once

#include <cmath>
#include <functional>

namespace kaspectra::math {

    namespace detail {

        inline double simpson(double fa, double fm, double fb, double a, double b) {
            return (b - a) / 6.0 * (fa + 4.0 * fm + fb);
        }

        inline double adaptive_simpson(const std::function<double(double)>& f,
                                        double a, double b, double fa, double fm, double fb,
                                        double whole, double abs_tol, double rel_tol, int depth, int max_depth) {
            double m = 0.5 * (a + b);
            double lm = 0.5 * (a + m);
            double rm = 0.5 * (m + b);
            double flm = f(lm);
            double frm = f(rm);

            double left = simpson(fa, flm, fm, a, m);
            double right = simpson(fm, frm, fb, m, b);
            double combined = left + right;

            double tol = std::max(abs_tol, rel_tol * std::fabs(combined));

            if (depth >= max_depth || std::fabs(combined - whole) <= 15.0 * tol) {
                return combined + (combined - whole) / 15.0;
            }

            return adaptive_simpson(f, a, m, fa, flm, fm, left, abs_tol, rel_tol / 2.0, depth + 1, max_depth) +
                    adaptive_simpson(f, m, b, fm, frm, fb, right, abs_tol, rel_tol / 2.0, depth + 1, max_depth);
        }

    } // namespace detail

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