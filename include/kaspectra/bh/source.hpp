#pragma once

#include <cmath>

#include "kaspectra/constants.hpp"
#include "kaspectra/bh/fits.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/io/photon_field.hpp"
#include "kaspectra/math/integrate.hpp"

namespace kaspectra::bh {

    inline double kappa(double eps, double E_p) {
        using namespace kaspectra::constants;
        return 2.0 * E_p * eps / (m_p * m_e);
    }

    inline double eps_min(double E_p) {
        using namespace kaspectra::constants;
        return m_p * m_e / E_p;
    }

    inline double interaction_rate(double E_p, const io::PhotonField& f_ph, double epsilon_max,
                                    double abs_tol = 1e-10, double rel_tol = 1e-8,
                                    int max_depth = 50, int panels = math::kDefaultPanels) {
        const double eps_lo = eps_min(E_p);
        if (eps_lo >= epsilon_max) return 0.0;

        auto integrand = [&](double eps) {
            return f_ph(eps) * psi_over_kappa2(kappa(eps, E_p));
        };

        return 2.0 * constants::alpha_r0sq_c * math::integrate_log(integrand, eps_lo, epsilon_max, abs_tol, rel_tol, max_depth, panels);
    }

    // -dE_p/dt [GeV/s]. eq.(3.11) (Chodorowski/Zdziarski/Sikora 1992) gives -dgamma/dt
    // proportional to a single power of gamma (via the kappa=2*eps*gamma/m_e substitution's
    // Jacobian); E_p = gamma*m_p introduces no further power when converting to -dE_p/dt.
    // Only ONE factor of (E_p/m_p), not squared.
    inline double energy_loss_rate(double E_p, const io::PhotonField& f_ph, double epsilon_max,
                                    double abs_tol = 1e-10, double rel_tol = 1e-8,
                                    int max_depth = 50, int panels = math::kDefaultPanels) {
        using namespace kaspectra::constants;

        const double eps_lo = eps_min(E_p);
        if (eps_lo >= epsilon_max) return 0.0;

        auto integrand = [&](double eps) {
            return f_ph(eps) * phi_over_kappa2(kappa(eps, E_p));
        };

        return 2.0 * alpha_r0sq_c * m_e * (E_p / m_p) *
            math::integrate_log(integrand, eps_lo, epsilon_max, abs_tol, rel_tol, max_depth, panels);
    }

    inline double q_pair_rate(const io::ProtonSpectrum& J_p, const io::PhotonField& f_ph,
                                double E_p_max, double epsilon_max,
                                double abs_tol = 1e-10, double rel_tol = 1e-8,
                                int max_depth = 50, int panels = math::kDefaultPanels) {
        using namespace kaspectra::constants;

        const double E_p_lo = m_p * m_e / epsilon_max;
        if (E_p_lo >= E_p_max) return 0.0;

        auto integrand = [&](double E_p) {
            return J_p(E_p) * interaction_rate(E_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels);
        };

        return math::integrate_log(integrand, E_p_lo, E_p_max, abs_tol, rel_tol, max_depth, panels);
    }

    inline double q_pair_energy_loss(const io::ProtonSpectrum& J_p, const io::PhotonField& f_ph,
                                        double E_p_max, double epsilon_max,
                                        double abs_tol = 1e-10, double rel_tol = 1e-8,
                                        int max_depth = 50, int panels = math::kDefaultPanels) {
        using namespace kaspectra::constants;

        const double E_p_lo = m_p * m_e / epsilon_max;
        if (E_p_lo >= E_p_max) return 0.0;

        auto integrand = [&](double E_p) {
            return J_p(E_p) * energy_loss_rate(E_p, f_ph, epsilon_max, abs_tol, rel_tol, max_depth, panels);
        };

        return math::integrate_log(integrand, E_p_lo, E_p_max, abs_tol, rel_tol, max_depth, panels);
    }

}   // namespace kaspectra::bh
