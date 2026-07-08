#pragma once

#include "kaspectra/constants.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/pgamma/kinematics.hpp"
#include "kaspectra/pgamma/secondary_spectra.hpp"
#include "kaspectra/math/integrate.hpp"
#include "kaspectra/species.hpp"
#include <cmath>
#include <stdexcept>

namespace kaspectra::pgamma {

    using kaspectra::Species;

    inline double Phi_species(Species s, double x, double eta) {
        switch (s) {
            case Species::Gamma:    return Phi_gamma(x, eta);
            case Species::Positron: return Phi_e_plus(x, eta);
            case Species::Electron: return Phi_e_minus(x, eta);
            case Species::NuMu:     return Phi_numu(x, eta);
            case Species::NuMuBar:  return Phi_numubar(x, eta);
            case Species::NuE:      return Phi_nue(x, eta);
            case Species::NuEBar:   return Phi_nuebar(x, eta);
        }
        throw std::invalid_argument("Phi_species: unknown Species");
    }

    namespace detail {

        // R for the threshold each species' Phi_* actually gates on
        // internally: R=1 (single-pion) for gamma/e+/numu/numubar/nue;
        // R=1+r (two-pion, eq.38) for e-/nuebar only. Duplicated here
        // purely as an efficiency device to narrow the eps integration
        // window - Phi_species is safe at any eta regardless.
        inline double R_for_species(Species s) {
            switch (s) {
                case Species::Electron:
                case Species::NuEBar:
                    return 1.0 + kaspectra::constants::r;
                default:
                    return 1.0;
            }
        }
    } // namespace detail

    inline double q_species(Species s, double E, const io::ProtonSpectrum& J_p,
                            const io::PhotonField& f_ph, double E_p_max, double epsilon_max,
                            double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {

        using namespace kaspectra::constants;
        if (E >= E_p_max) return 0.0;

        const double R = detail::R_for_species(s);
        const double eta_th = eta_threshold(R);

        auto outer_integrand = [&](double E_p) {
            const double x = E / E_p;

            const double eps_lo = eta_th * m_p * m_p / (4.0 * E_p);
            if (eps_lo >= epsilon_max) return 0.0;

            auto inner_integrand = [&](double epsilon) {
                return f_ph(epsilon) * Phi_species(s, x, eta(epsilon, E_p));
            };

            const double inner = math::integrate_log(inner_integrand, eps_lo, epsilon_max, abs_tol, rel_tol, max_depth);
            
            return J_p(E_p) * inner / E_p;
        };

        return math::integrate_log(outer_integrand, E, E_p_max, abs_tol, rel_tol, max_depth);
    }

} // namespace kaspectra::pgamma