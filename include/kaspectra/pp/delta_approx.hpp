#pragma once

#include "kaspectra/constants.hpp"
#include "kaspectra/math/integrate.hpp"
#include "kaspectra/io/proton_spectrum.hpp"
#include "kaspectra/pp/cross_section.hpp"
#include <cmath>

namespace kaspectra::pp {

    // GAMMA-RAY CHANNEL ONLY: arXiv:astro-ph/0606058 works out pi0 -> 2gamma kernel exactly (eq. 78)
    // but gives no explicit lepton-decay kernel for this approximation, so 
    // e+/e-/numu/nue deliberately not implemented here (todo for future)
    //
    // Stitched agains secondary_sepctra.hpp bya. future pp/source.hpp; this
    // header intentionally does not include secondary_spectra.hpp or attempt
    // that mathing itself

    constexpr double K_pi_default = 0.17;

    // eq. 77. n_tilde is NOT a universal constant - obtained by matching this
    // spectrum to the accurate parametrisation at E_delta_approx_max. q_pi (and
    // dN_gamma_dE_gamma_delta_approx below) are exactly linear in n_tilde, so
    // callers compute the unit spectrum (n_tilde=1) once and rescale by a ratio
    // rather than solving for n_tilde iteratively
    inline double q_pi(double E_pi, double n_tilde, const kaspectra::io::ProtonSpectrum& J_p,
                        double n_H, double K_pi = K_pi_default) {
        using namespace kaspectra::constants;
        const double E_p = m_p + E_pi / K_pi;
        return n_tilde * (c_light * n_H / K_pi) * sigma_inel(E_p) * J_p(E_p);
    }

    // eq. 78. E_pi_max stands in for the infinite upper bound; pick it
    // several e-foldings past J_p's own cutoff/support
    // Reuse m_pi_charged
    inline double dN_gamma_dE_gamma_delta_approx(double E_gamma, double n_tilde,
                                                const kaspectra::io::ProtonSpectrum& J_p,
                                                double n_H, double E_pi_max,
                                                double K_pi = K_pi_default,
                                                double abs_tol = 1e-10, double rel_tol = 1e-8,
                                                int max_depth = 50) {
        using namespace kaspectra::constants;
        if (E_gamma <= 0.0) return 0.0;

        const double m_pi2 = m_pi_charged * m_pi_charged;
        const double E_min = E_gamma + m_pi2 / (4.0 * E_gamma);
        if (E_pi_max <= E_min) return 0.0;

        auto integrand = [&](double E_pi) {
            return q_pi(E_pi, n_tilde, J_p, n_H, K_pi) / std::sqrt(E_pi * E_pi - m_pi2);
        };

        return 2.0 * kaspectra::math::integrate(integrand, E_min, E_pi_max, abs_tol, rel_tol, max_depth);
    }

    // Usage:
    //   kaspectra::io::PowerLawSPectrum J_p(1.0, 2.0);
    //   double unit_spec = dN_gamma_dE_gamma_delta_approx(1.0, 1.0, J_p, 1.0, 1e5);
    //   // rescale n_tilde = accurate_spectrum_at_stitch / unit_spec_at_stitch

} // namespace kaspectra::pp