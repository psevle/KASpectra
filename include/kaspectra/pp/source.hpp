#pragma once

#include "kaspectra/constants.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/pp/cross_section.hpp"
#include "kaspectra/pp/secondary_spectra.hpp"
#include "kaspectra/pp/delta_approx.hpp"
#include "kaspectra/math/integrate.hpp"
#include <cmath>
#include <stdexcept>

// Usage:
//   kaspectra::io::PowerLawSpectrum J_p(1.0, 2.0);
//   for (double E: {1.0, 10.0, 100.0, 1000.0})
//      double q = kaspectra::pp::q_species(kaspectra::pp::Species::Gamma, E, J_p, 1.0, 1e8, 1e8);

namespace kaspectra::pp {

    // NuMuBar kept only for API symmetry; as per astro-ph/0606058,
    // (nu==nubar, e+==e-) it's numerically identical to NuMu.
    enum class Species { Gamma, ElectronPositron, NuMu, NuMuBar, NuE };

    inline double F_species(Species s, double x, double E_p) {
        switch (s) {
            case Species::Gamma:            return F_gamma(x, E_p);
            case Species::ElectronPositron: return F_e(x, E_p);
            case Species::NuMu:
            case Species::NuMuBar:          return F_numu(x, E_p);
            case Species::NuE:              return F_nue(x, E_p);
        }
        throw std::invalid_argument("q_species: unknown Species");
    }

    // eq. 71/72 E_p_max stands in for the implicit infinity.
    // Lower bound is max(E, E_delta_approx_max), NOT E_threshold_pp: F_species
    // (eq. 58/62/66) is only fitted for E_p >= E_delta_approx_max (100 GeV).
    // Below that, F_e's beta_e = (0.201+0.062L+0.00042L^2)^-0.25 has a negative
    // base for E_p roughly in [1.22, 37] GeV (L in [-144,-3.32]), producing NaN
    // that used to send adaptive Simpson spinning to max_depth on every
    // affected branch - a de facto hang, not just an accuracy loss. Clamping
    // here means low-E secondary spectra (already known to be inaccurate for
    // leptons per q_species's docs) are a truncated integral instead of NaN.
    inline double q_accurate(Species s, double E, const io::ProtonSpectrum& J_p,
                            double n_H, double E_p_max,
                            double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        using namespace kaspectra::constants;
        const double E_lo = std::max(E, E_delta_approx_max);
        if (E_lo >= E_p_max) return 0.0;
        
        auto integrand = [&](double E_p) {
            return sigma_inel(E_p) * J_p(E_p) * F_species(s, E / E_p, E_p) / E_p;
        };

        return c_light * n_H * math::integrate_log(integrand, E_lo, E_p_max, abs_tol, rel_tol, max_depth);
    }

    // n_tilde = accurate(E_delta_approx_max) / unit_delta_approx(E_delta_approx_max)
    // Exposed standalone so callers sweeping many low-E gamma points compute it
    // once instead of paying for two extra integrals on every call
    inline double n_tilde_for(const io::ProtonSpectrum& J_p, double n_H,
                                double E_p_max, double E_pi_max,
                                double K_pi = K_pi_default,
                                double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        using namespace kaspectra::constants;
        const double q_acc = q_accurate(Species::Gamma, E_delta_approx_max, J_p, n_H, E_p_max, abs_tol, rel_tol, max_depth);
        const double q_unit = dN_gamma_dE_gamma_delta_approx(E_delta_approx_max, 1.0, J_p, n_H, E_pi_max, K_pi, abs_tol, rel_tol, max_depth);
        if (q_unit == 0.0) return 0.0;
        return q_acc / q_unit;
    }

    // Convenience stitcher: recompute n_tilde_for() on every low-E call
    inline double q_gamma_full(double E_gamma, const io::ProtonSpectrum& J_p, double n_H,
                                double E_p_max, double E_pi_max,
                                double K_pi = K_pi_default,
                                double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        using namespace kaspectra::constants;
        if (E_gamma >= E_delta_approx_max) {
            return q_accurate(Species::Gamma, E_gamma, J_p, n_H, E_p_max, abs_tol, rel_tol, max_depth);
        }
        const double n_tilde = n_tilde_for(J_p, n_H, E_p_max, E_pi_max, K_pi, abs_tol, rel_tol, max_depth);
        return dN_gamma_dE_gamma_delta_approx(E_gamma, n_tilde, J_p, n_H, E_pi_max, K_pi, abs_tol, rel_tol, max_depth);
    }

    // Only Gamma gets the low-E stitch: astro-ph/0606058 has no delta-function extension
    // for leptons, so other species below 100 GeV are q_accurate() extrapolated
    // outside its fitted validity range - not expected to be accurate there
    inline double q_species(Species s, double E, const io::ProtonSpectrum& J_p, double n_H,
                            double E_p_max, double E_pi_max,
                            double K_pi = K_pi_default,
                            double abs_tol = 1e-10, double rel_tol = 1e-8, int max_depth = 50) {
        if (s == Species::Gamma) {
            return q_gamma_full(E, J_p, n_H, E_p_max, E_pi_max, K_pi, abs_tol, rel_tol, max_depth);
        }
        return q_accurate(s, E, J_p, n_H, E_p_max, abs_tol, rel_tol, max_depth);
    }

}   // namespace kaspectra::pp