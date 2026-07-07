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

    inline double integrate(const std::function<double(double)>& f, double a, double b,
                            double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        if (a == b) return 0.0;
        double sign = 1.0;
        if (a > b) { std::swap(a, b); sign = -1.0; }

        double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
        double whole = detail::simpson(fa, fm, fb, a, b);
        return sign * detail::adaptive_simpson(f, a, b, fa, fm, fb, whole, abs_tol, rel_tol, 0, max_depth);
    }

    // u = log(E) substitution: integral f(E) dE == integral f(exp(u))*exp(u) du
    // E_p spans ~1 GeV to ~PeV
    inline double integrate_log(const std::function<double(double)>& f, double a, double b,
                                double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        auto g = [&f](double u) { double e = std::exp(u); return f(e) * e; };
        return integrate(g, std::log(a), std::log(b), abs_tol, rel_tol, max_depth);
    }

// Usage
//   double result = integrate_log([](double e){ return sigma_inel(e) * J_p(e); },
//                                  1.0 /*GeV*/, 1e6 /*GeV*/);

}   // namespace kaspectra::math