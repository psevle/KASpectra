#pragma once

#include <cmath>
#include <utility>

#include "kaspectra/constants.hpp"
#include "kaspectra/pgamma/kinematics.hpp"
#include "kaspectra/pgamma/phi_shape.hpp"
#include "kaspectra/pgamma/phi_tables.hpp"

namespace kaspectra::pgamma {

    namespace detail {
        
        // psi(eta), eq 34 - shared by e+, numubar, nue and numu
        inline double psi_lepton_shared(double eta) {
            return 2.5 + 1.4 * std::log(eta / constants::eta_0);
        }

        // psi(eta), eq 41 - shared by e-/nuebar
        inline double psi_e_minus_family(double eta) {
            const double rho = eta / constants::eta_0;
            if (rho < 4.0) return 0.0;
            return 6.0 * (1.0 - std::exp(1.5 * (4.0 - rho)));
        }

        // x'_=(rho), eq 36
        inline double numu_x_plus(double rho, double x_plus_R1) {
            if (rho < 2.14) {
                return 0.427 * x_plus_R1;
            }
            else if (rho < 10.0) {
                return (0.427 + 0.0729 * (rho - 2.14)) * x_plus_R1;
            }
            else {
                return x_plus_R1;
            }
        }

    } // namespace detail

    inline double Phi_gamma(double x, double eta) {
        if (eta < eta_threshold(1.0)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_I(eta_over_eta0);
        const double psi = 2.5 + 0.4 * std::log(eta_over_eta0);
        const auto [x_minus, x_plus] = x_pm(eta, 1.0);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_minus, x_plus);
    }

    inline double Phi_e_plus(double x, double eta) {
        if (eta < eta_threshold(1.0)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_II_e_plus(eta_over_eta0);
        const double psi = detail::psi_lepton_shared(eta);
        const auto [x_minus_r1, x_plus_r1] = x_pm(eta, 1.0);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_minus_r1 / 4.0, x_plus_r1);
    }

    inline double Phi_numubar(double x, double eta) {
        if (eta < eta_threshold(1.0)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_II_numubar(eta_over_eta0);
        const double psi = detail::psi_lepton_shared(eta);
        const auto [x_minus_r1, x_plus_r1] = x_pm(eta, 1.0);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_minus_r1 / 4.0, x_plus_r1);
    }

    inline double Phi_nue(double x, double eta) {
        if (eta < eta_threshold(1.0)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_II_nue(eta_over_eta0);
        const double psi = detail::psi_lepton_shared(eta);
        const auto [x_minus_r1, x_plus_r1] = x_pm(eta, 1.0);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_minus_r1 / 4.0, x_plus_r1);
    }

    inline double Phi_numu(double x, double eta) {
        if (eta < eta_threshold(1.0)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const double rho = eta_over_eta0;
        const auto params = lookup_table_II_numu(eta_over_eta0);
        const double psi = detail::psi_lepton_shared(eta);
        const auto [x_minus_r1, x_plus_r1] = x_pm(eta, 1.0);

        const double x_minus = 0.427 * x_minus_r1;
        const double x_plus = detail::numu_x_plus(rho, x_plus_r1);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_minus, x_plus);
    }

    inline double Phi_e_minus(double x, double eta) {
        if (eta < eta_threshold(1.0 + constants::r)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_III_e_minus(eta_over_eta0);
        const double psi = detail::psi_e_minus_family(eta);
        const auto [x_min, x_max] = x_pm(eta, 1.0 + constants::r);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_min / 2.0, x_max);
    }

    inline double Phi_nuebar(double x, double eta) {
        if (eta < eta_threshold(1.0 + constants::r)) return 0.0;

        const double eta_over_eta0 = eta / constants::eta_0;
        const auto params = lookup_table_III_nuebar(eta_over_eta0);
        const double psi = detail::psi_e_minus_family(eta);
        const auto [x_min, x_max] = x_pm(eta, 1.0 + constants::r);

        return phi_shape(x, params.B, params.s, params.delta, psi, x_min / 2.0, x_max);
    }

} // namespace kaspectra::pgamma