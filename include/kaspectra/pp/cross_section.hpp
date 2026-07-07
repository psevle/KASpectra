#pragma once

#include "kaspectra/constants.hpp"
#include <cmath>

namespace kaspectra::pp {

    // Inelastic pp cross section, arXiv:astro-ph/0606058 eq. 79
    // Input: E_p = proton TOTAL energy [GeV] (not just kinetic)
    // Output: corss section in cm^2
    // L = ln(E_p / 1000 GeV), i.e. E_p expressed in TeV
    // Returns 0 for E_p <= E_threshold_pp (below pion-production threshold,
    // no inelastic/secondary-producing channel is used in this parametrisation)
    // Usage:
    //   double sigma = kaspectra::pp::sigma_inel(10.0); // E_p = 10 GeV, sigma in cm^2
inline double sigma_inel(double E_p) {
    using namespace kaspectra::constants;
    
    if (E_p <= E_threshold_pp) {
        return 0.0;
    }

    const double L = std::log(E_p / 1000.0); // E_p/1 TeV in GeV
    const double ratio = E_threshold_pp / E_p;
    const double bracket = 1.0 - ratio * ratio * ratio * ratio; // (E_th/E_p)^4

    const double sigma_mb = (34.3 + 1.88 * L + 0.25 * L * L) * bracket * bracket;

    constexpr double mb_to_cm2 = 1.0e-27;
    return sigma_mb * mb_to_cm2;
}

} // namespace kaspectra::pp