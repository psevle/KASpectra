#pragma once

#include <cmath>

#include "kaspectra/constants.hpp"

// Parametrised secondary particle spectra from pp interactions.
// Kelner, Aharonian & Bugayov 2006 (arXiv:astro-ph/0606058) eqs. 58, 62, 66
// x = E_secondary / E_p (E_p = incident proton total energy, GeV)
// Valid range: 0.1 TeV <= E_p <= 1e5 TeV, x >~ 1e-3
// Uses assumptions that nu == nubar and e+ == e-

namespace kaspectra::pp {

    namespace detail {

        inline double log_energy_scale(double E_p) {
            return std::log(E_p / 1000.0);
        }

        // Shared functional form for both eq. 58 & 66
        // Same shape, different (B, beta, k) coefficient triples and different
        // arguments (x directly for eq. 58, y = x/0.427 for eq. 66)
        inline double kelner_shape(double u, double B, double beta, double k) {
            const double ub = std::pow(u, beta);
            const double one_minus_ub = 1.0 - ub;
            const double denom = 1.0 + k * ub * one_minus_ub;

            const double term = std::pow(one_minus_ub / denom, 4);
            const double bracket = 1.0 / std::log(u)
                - 4.0 * beta * ub / one_minus_ub
                - 4.0 * k * beta * ub * (1.0 - 2.0 * ub) / denom;
            
            return B * (std::log(u) / u) * term * bracket;
        }

    } // namespace detail

    inline double F_gamma(double x, double E_p) {
        if (x <= 0.0 || x >= 1.0) return 0.0;

        const double L = detail::log_energy_scale(E_p);
        const double B_g = 1.30 + 0.14 * L + 0.011 * L * L;
        const double beta_g = 1.0 / (1.79 + 0.11 * L + 0.008 * L * L);
        const double k_g = 1.0 / (0.801 + 0.049 * L + 0.014 * L * L);

        return detail::kelner_shape(x, B_g, beta_g, k_g);
    }

    inline double F_e(double x, double E_p) {
        if (x <= 0.0 || x>= 1.0) return 0.0;

        const double L = detail::log_energy_scale(E_p);
        const double B_e = 1.0 / (69.5 + 2.65 * L + 0.3 * L * L);
        const double sum = 0.201 + 0.062 * L + 0.00042 * L * L;
        const double beta_e = std::pow(sum, -0.25);
        const double k_e = (0.279 + 0.141 * L + 0.0172 * L * L)
            / (0.3 + (2.3 + L) * (2.3 + L));
        
        const double neg_ln_x = -std::log(x);
        return B_e * std::pow(1.0 + k_e * std::log(x) * std::log(x), 3)
            / (x * (1.0 + 0.3 / std::pow(x, beta_e)))
            * std::pow(neg_ln_x, 5);
    }

    inline double F_numu1(double x, double E_p) {
        if (x <= 0.0 || x >= 0.427) return 0.0;

        const double L = detail::log_energy_scale(E_p);
        const double y = x / 0.427;
        const double B_p = 1.75 + 0.204 * L + 0.010 * L * L;
        const double beta_p = 1.0 / (1.67 + 0.111 * L + 0.0038 * L * L);
        const double k_p = 1.07 - 0.086 * L + 0.002 * L * L;

        return detail::kelner_shape(y, B_p, beta_p, k_p);
    }

    // F_numu2 (second numu, from mu decay)
    inline double F_numu(double x, double E_p) {
        return F_numu1(x, E_p) + F_e(x, E_p);
    }

    // F_nue ~ F_e to <5%, so we stick with this approcimation for now
    inline double F_nue(double x, double E_p) {
        return F_e(x, E_p);
    }

} // namespace kaspectra::pp