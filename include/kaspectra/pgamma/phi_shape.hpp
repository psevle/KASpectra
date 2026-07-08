#pragma once

#include <cmath>

namespace kaspectra::pgamma {

    inline double phi_shape(double x, double B, double s, double delta, double psi, double x_minus, double x_plus) {
        if (x <= 0.0 || x_minus <= 0.0 || x_minus >= x_plus) {
            return 0.0;
        }
        if (x >= x_plus) {
            return 0.0;
        }
        if (x < x_minus) {
            return B * std::pow(std::log(2.0), psi);
        }
        const double y = (x - x_minus) / (x_plus - x_minus);
        const double exp_term = std::exp(-s * std::pow(std::log(x / x_minus), delta));
        const double log_term = std::pow(std::log(2.0 / (1.0 + y * y)), psi);
        return B * exp_term * log_term;
    }

} // namespace kaspectra::pgamma