#pragma once

#include <cmath>
#include <utility>

#include "kaspectra/constants.hpp"

namespace kaspectra::pgamma {

    inline double eta(double epsilon, double E_p) {
        using namespace kaspectra::constants;
        return 4.0 * epsilon * E_p / (m_p * m_p);
    }

    inline double eta_threshold(double R) {
        using namespace kaspectra::constants;
        const double Rr_sum = R + r;
        return Rr_sum * Rr_sum - 1.0;
    }

    inline std::pair<double, double> x_pm(double eta, double R) {
        using namespace kaspectra::constants;

        const double Rr_sum = R + r;
        const double Rr_diff = R - r;
        const double discriminant = (eta + 1.0 - Rr_sum * Rr_sum) * (eta + 1.0 - Rr_diff * Rr_diff);
        const double sqrt_disc = std::sqrt(discriminant); // NaN if disc < 0

        const double numer_common = eta + r * r + 1.0 - R * R;
        const double denom = 2.0 * (1.0 + eta);

        return { (numer_common - sqrt_disc) / denom, (numer_common + sqrt_disc) / denom };
    }

} // namespace kaspectra::pgamma