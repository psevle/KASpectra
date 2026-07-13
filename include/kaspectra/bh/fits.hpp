#pragma once

#include <cmath>

namespace kaspectra::bh {
    namespace detail {
        constexpr double zeta3 = 1.2020569031595942854;
        constexpr double pi_bh = 3.14159265358979323846;
    }

    inline double psi_over_kappa2(double kappa) {
        if (kappa < 2.0) return 0.0;
        double psi;
        if (kappa < 4.0) {
            // Chodorowski eq.2.4, restabilized in x = kappa-2. The printed form
            // sums terms of order ~100-300 that cancel EXACTLY through cubic
            // order at threshold (verified analytically: the rational constants
            // sum to 0 over a common denominator of 1440, the log constants
            // give (-19151+1007)/48*ln2 + 189*ln4 == 0, and a 60-digit Taylor
            // expansion shows orders 0-3 vanish identically, psi ~ x^4/16 --
            // matching phi's own explicit (pi/12)x^4 threshold form). Naive
            // evaluation therefore loses ALL significant digits by x ~ 1e-4.
            // Two-regime fix:
            //   x < 1/4:  Taylor series of the printed form (coefficients from
            //             a 60-digit mpmath expansion; a4 = 1/16 and
            //             a5 = -7/160 are exact; convergence radius 2, so
            //             truncation at x^16 is ~6e-12 relative at the switch).
            //   x >= 1/4: the printed form with its (identically zero) constant
            //             block dropped and every term written as a
            //             cancellation-free difference from its kappa=2 value
            //             (log1p, factored rationals); rounding-level agreement
            //             with the printed form at mid-branch, ~6e-12 agreement
            //             with the series at the switch point.
            const double x = kappa - 2.0;
            if (x < 0.25) {
                constexpr double a[13] = {
                    0.0625,                            // 1/16, exact
                    -0.04375,                          // -7/160, exact
                    0.0262369791666666666667,
                    -0.0148670014880952380952,
                    0.0081645965576171875,
                    -0.00439973054108796296296,
                    0.00234190622965494791667,
                    -0.0012355515451142282197,
                    0.000647288892004224989149,
                    -0.000337103391304994240785,
                    0.000174664883386521112351,
                    -0.0000900987121793958875868,
                    0.0000462980523783092697461,
                };
                double s = 0.0;
                for (int i = 12; i >= 0; --i) s = s * x + a[i];
                psi = (2.0 * detail::pi_bh / 3.0) * s * x * x * x * x;
            } else {
                const double q = 4.0 + x;             // == 2 + kappa
                const double q2 = q * q, q3 = q2 * q;
                psi = (2.0 * detail::pi_bh / 3.0) * (
                    - (347.0 / 80.0) * x / (2.0 + x)
                    - (15163.0 / 480.0) * x
                    + (2593.0 / 1920.0) * x * (x + 4.0)
                    - (3904.0 / 9.0) * x * (x * x + 12.0 * x + 48.0) / (64.0 * q3)
                    + (9688.0 / 15.0) * x * (x + 8.0) / (16.0 * q2)
                    - (10676.0 / 15.0) * x / (4.0 * q)
                    + (1007.0 / 48.0) * std::log1p(x / 2.0)
                    + 189.0 * std::log1p(x / 4.0)
                );
            }
        }
        else {
            const double k = kappa;
            const double L = std::log(2.0 * k);
            const double k2 = k * k, k4 = k2 * k2;
            const double L2 = L * L, L3 = L2 * L, L4 = L2 * L2;
            const double pi2 = detail::pi_bh * detail::pi_bh;
            psi = -67.0 / (1728.0 * k4)
                + 29.0 * L / (144.0 * k4)
                + 7.0 / (4.0 * k2)
                  + 3.0 * L / (2.0 * k2)
                - 2.71245
                + 2.0 * (pi2 / 3.0 - 7.0 + 4.0 * detail::zeta3) * L
                + 2.0 * (6.0 - pi2 / 3.0) * L2
                - (4.0 / 3.0) * L3
                + (2.0 / 3.0) * L4
                - (130.0 / 27.0) * k2
                + (14.0 / 9.0) * k2 * L;
        }
        return psi / (kappa * kappa);
    }

    inline double phi_over_kappa2(double kappa) {
        if (kappa < 2.0) return 0.0;
        double phi;
        if (kappa < 25.0) {
            const double x = kappa - 2.0;
            const double c1 = 0.8048, c2 = 0.1459, c3 = 1.137e-3, c4 = -3.879e-6;
            const double x2 = x * x, x4 = x2 * x2;
            const double denom = 1.0 + c1 * x + c2 * x2 + c3 * x2 * x + c4 * x4;
            phi = (detail::pi_bh / 12.0) * x4 / denom;
        }
        else {
            const double k = kappa;
            const double L = std::log(k);
            const double ln2 = std::log(2.0);
            const double pi2 = detail::pi_bh * detail::pi_bh;
            const double d0 = -170.0 + 84.0 * ln2 - 16.0 * ln2 * ln2 + (pi2 / 3.0) * (10.0 - 4.0 * ln2) + 8.0 * detail::zeta3;
            const double d1 = 88.0 - 40.0 * ln2 + 8.0 * ln2 * ln2 - (4.0 / 3.0) * pi2;
            const double d2 = -20.0 + 8.0 * ln2;
            const double d3 = 8.0 / 3.0;
            const double phi_asym = k * (d0 + d1 * L + d2 * L * L + d3 * L * L * L);
            const double f1 = 2.910, f2 = 78.35, f3 = 1837.0;
            const double denom = 1.0 - f1 / k - f2 / (k * k) - f3 / (k * k * k);
            phi = phi_asym / denom;
        }
        return phi / (kappa * kappa);
    }

}   // namespace kaspectra::bh