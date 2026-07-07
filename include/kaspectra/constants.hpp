#pragma once

// All energy/mass values are expressed in GeV unless otherwise specified
// Natural units used, except where otherwise specified

namespace kaspectra::constants {
    constexpr double m_p = 0.9382720813; // Proton mass
    constexpr double m_pi_charged = 0.13957039; // Charged pion mass
    constexpr double m_mu = 0.1056583745; // Muon mass

    // Total (rest + kinetic) energy threshold for pi0/pi+ production in pp collisions
    // E_th = m_p + T_th, T_th = 2*m_pi + m_pi^2/(2*m_p) ~ 0.28 GeV kinetic
    // -> E_th ~ 1.22 GeV total. arXiv:astro-ph/0606058 quotes 1.22 GeV;
    // Taken here as total incident-proton energy, not kinetic-only
    constexpr double E_threshold_pp = 1.22;

    constexpr double c_light = 2.99792458e10; // Speed of light in cm/s

    // Boundary in which arXiv:astro-ph/0606058 parametrised high-energy spectra apply
    // instead of low-energy delta-functional approximation
    constexpr double E_delta_approx_max = 100.0;

}   // namespace kaspectra::constants