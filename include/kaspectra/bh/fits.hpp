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
            const double k = kappa;
            psi = (2.0 * detail::pi_bh / 3.0) * (
                347.0 / (40.0 * k)
                - 130903.0 / 1440.0
                - (19151.0 / 48.0) * std::log(2.0)
                - (15163.0 / 480.0) * k
                + (2593.0 / 1920.0) * k * k
                + 3904.0 / (9.0 * (2.0 + k) * (2.0 + k) * (2.0 + k))
                - 9688.0 / (15.0 * (2.0 + k) * (2.0 + k))
                + 10676.0 / (15.0 * (2.0 + k))
                + (1007.0 / 48.0) * std::log(k)
                + 189.0 * std::log(2.0 + k)
            );
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