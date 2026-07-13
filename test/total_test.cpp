#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include <kaspectra/total.hpp>

using namespace kaspectra;
using kaspectra::io::BlackbodyPhotonField;
using kaspectra::io::PowerLawSpectrum;

// Config chosen so every channel is genuinely exercised at bounded cost:
// E = 1e3 GeV secondaries, E_p budget 1e6 GeV, warm blackbody target
// (kT = 1e-7 GeV). Against the real CMB the pair channel at this E_p budget
// is only "reachable" via photons ~4000 kT into the Wien tail -- occupancy
// e^-4000, a true zero in double precision -- so a warmer field is needed
// for the bh term to be a real nonzero number here (same reasoning as the
// warm-field configs in test/bh/spectrum_cache_test.cpp). bh's
// E_p_min(1e3, 2e-6) ~ 1.1e4 GeV sits well inside the budget, and the
// Planckian fast path keeps the cost sub-second.
namespace {
    const PowerLawSpectrum J_p(1.0, 2.0);
    const BlackbodyPhotonField cmb(1e-7 / kaspectra::constants::k_boltzmann);
    constexpr double E = 1e3, E_p_max = 1e6, E_pi_max = 1e6, eps_max = 2e-6;
}

TEST_CASE("q_total_species equals the explicit sum of its channels", "[total]") {
    // gamma: pp + pgamma only (no bh contribution for photons).
    double total_gamma = q_total_species(Species::Gamma, E, J_p, 1.0, cmb, E_p_max, E_pi_max, eps_max);
    double parts_gamma = pp::q_species(Species::Gamma, E, J_p, 1.0, E_p_max, E_pi_max)
                       + pgamma::q_species(Species::Gamma, E, J_p, cmb, E_p_max, eps_max);
    REQUIRE(total_gamma == Catch::Approx(parts_gamma).epsilon(1e-12));

    // electron: pp + pgamma + bh.
    double total_e = q_total_species(Species::Electron, E, J_p, 1.0, cmb, E_p_max, E_pi_max, eps_max);
    double bh_part = bh::q_pair_spectrum(E, J_p, cmb, E_p_max, eps_max);
    double parts_e = pp::q_species(Species::Electron, E, J_p, 1.0, E_p_max, E_pi_max)
                   + pgamma::q_species(Species::Electron, E, J_p, cmb, E_p_max, eps_max)
                   + bh_part;
    CAPTURE(bh_part);
    REQUIRE(bh_part > 0.0);   // config was chosen to make the pair channel real
    REQUIRE(total_e == Catch::Approx(parts_e).epsilon(1e-12));
}

TEST_CASE("q_total_species with n_H == 0 drops the pp channel", "[total]") {
    double total = q_total_species(Species::Gamma, E, J_p, 0.0, cmb, E_p_max, E_pi_max, eps_max);
    double pg_only = pgamma::q_species(Species::Gamma, E, J_p, cmb, E_p_max, eps_max);
    REQUIRE(total == Catch::Approx(pg_only).epsilon(1e-12));
}

TEST_CASE("q_total_species is finite and non-negative for all species", "[total]") {
    for (Species s : {Species::Gamma, Species::Positron, Species::Electron,
                      Species::NuMu, Species::NuMuBar, Species::NuE, Species::NuEBar}) {
        double q = q_total_species(s, E, J_p, 1.0, cmb, E_p_max, E_pi_max, eps_max);
        CAPTURE(static_cast<int>(s), q);
        REQUIRE(std::isfinite(q));
        REQUIRE(q >= 0.0);
    }
}
