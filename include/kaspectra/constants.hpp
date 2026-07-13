#pragma once

// All energy/mass values are expressed in GeV unless otherwise specified
// Natural units used, except where otherwise specified

namespace kaspectra::constants {
    constexpr double m_p = 0.9382720813; // Proton mass
    constexpr double m_pi_charged = 0.13957039; // Charged pion mass
    constexpr double m_mu = 0.1056583745; // Muon mass
    constexpr double m_e = 5.10998950e-4; // Electron mass

    // Total (rest + kinetic) energy threshold for pi0/pi+ production in pp collisions
    // E_th = m_p + T_th, T_th = 2*m_pi + m_pi^2/(2*m_p) ~ 0.28 GeV kinetic
    // -> E_th ~ 1.22 GeV total. arXiv:astro-ph/0606058 quotes 1.22 GeV;
    // Taken here as total incident-proton energy, not kinetic-only
    constexpr double E_threshold_pp = 1.22;

    constexpr double c_light = 2.99792458e10; // Speed of light in cm/s

    // Boundary in which arXiv:astro-ph/0606058 parametrised high-energy spectra apply
    // instead of low-energy delta-functional approximation
    constexpr double E_delta_approx_max = 100.0;

    constexpr double pi = 3.14159265358979323846; // pi

    constexpr double k_boltzmann = 8.617333262e-14; // Boltzmann constant, GeV/K

    // (hbar*c)^3 in (GeV*cm)^3. CODATA hbar*c = 1.973269804e-14 GeV*cm;
    // stored as a single precomputed literal rather than cubed at compile
    // time, to avoid compounding rounding across three multiplications of
    // a ~1e-14 value
    constexpr double hbar_c3 = 7.683505e-42;

    // alpha * r0^2 * c in cm^3/s (alpha = fine-structure constant, r0 = classical
    // electron radius). Precomputed literal, same rationale as hbar_c3. Used by
    // bh/ (Chodorowski, Zdziarski & Sikora 1992, ApJ 400, 181).
    constexpr double alpha_r0sq_c = 1.737199184e-17;

    constexpr double r = m_pi_charged / m_p;    // ~0.146
    constexpr double eta_0 = 2.0 * r + r * r;

    // Unit-conversion conveniences for reporting observables in
    // astrophysically conventional units (the library itself works in
    // GeV/cm/s throughout).
    constexpr double cm_per_Mpc = 3.0856775814913673e24;   // IAU 2015 parsec definition x 1e6
    constexpr double seconds_per_year = 3.15576e7;         // Julian year

}   // namespace kaspectra::constants