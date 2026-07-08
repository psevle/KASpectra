#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/constants.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/pp/source.hpp>

using namespace kaspectra::pp;
using kaspectra::io::PowerLawSpectrum;

TEST_CASE("F_species dispatches to the matching free function", "[pp][source]") {
    double x = 0.1, E = 1000.0;
    REQUIRE(F_species(Species::Gamma, x, E) == F_gamma(x, E));
    REQUIRE(F_species(Species::Positron, x, E) == F_e(x, E));
    REQUIRE(F_species(Species::Electron, x, E) == F_e(x, E));
    REQUIRE(F_species(Species::NuMu, x, E) == F_numu(x, E));
    REQUIRE(F_species(Species::NuMuBar, x, E) == F_species(Species::NuMu, x, E));
    REQUIRE(F_species(Species::NuE, x, E) == F_nue(x, E));
}

TEST_CASE("q_accurate returns zero for a degenerate integration window", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    // E_p_max at or below the effective lower bound (E) collapses the window
    REQUIRE(q_accurate(Species::Gamma, 1000.0, J_p, 1.0, /*E_p_max=*/1000.0) == 0.0);
}

TEST_CASE("q_accurate is non-negative for a sane proton spectrum", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    for (double E : {10.0, 100.0, 1000.0}) {
        REQUIRE(q_accurate(Species::Gamma, E, J_p, 1.0, 1.0e6) >= 0.0);
    }
}

TEST_CASE("n_tilde_for returns zero if the unit delta-approx integral vanishes", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    // E_pi_max tiny enough that the eq.78 integration window is empty
    double n_tilde = n_tilde_for(J_p, 1.0, /*E_p_max=*/1.0e6, /*E_pi_max=*/0.001);
    REQUIRE(n_tilde == 0.0);
}

TEST_CASE("q_gamma_full is continuous across the E_delta_approx_max stitch", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double n_H = 1.0, E_p_max = 1.0e6, E_pi_max = 1.0e6;

    double at_stitch = q_gamma_full(kaspectra::constants::E_delta_approx_max, J_p, n_H, E_p_max, E_pi_max);
    double just_below = q_gamma_full(kaspectra::constants::E_delta_approx_max - 1e-3, J_p, n_H, E_p_max, E_pi_max);
    double just_above = q_gamma_full(kaspectra::constants::E_delta_approx_max + 1e-3, J_p, n_H, E_p_max, E_pi_max);

    REQUIRE(just_below == Catch::Approx(at_stitch).epsilon(1e-4));
    REQUIRE(just_above == Catch::Approx(at_stitch).epsilon(1e-4));
}

TEST_CASE("q_species does not stitch leptons at low energy", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double E = 10.0, n_H = 1.0, E_p_max = 1.0e6, E_pi_max = 1.0e6;
    REQUIRE(q_species(Species::NuMu, E, J_p, n_H, E_p_max, E_pi_max) ==
            q_accurate(Species::NuMu, E, J_p, n_H, E_p_max));
}

TEST_CASE("q_species returns a finite non-negative value for leptons below 100 GeV", "[pp][source]") {
    // Smoke test only: KAB06 gives no low-E extension for leptons, so this is
    // known to be an extrapolation outside the fit's validity range, not an
    // accuracy check.
    PowerLawSpectrum J_p(1.0, 2.0);
    double q = q_species(Species::NuMu, 10.0, J_p, 1.0, 1.0e6, 1.0e6);
    REQUIRE(std::isfinite(q));
    REQUIRE(q >= 0.0);
}

TEST_CASE("q_species routes Gamma through the stitched calculation", "[pp][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double E = 10.0, n_H = 1.0, E_p_max = 1.0e6, E_pi_max = 1.0e6;
    REQUIRE(q_species(Species::Gamma, E, J_p, n_H, E_p_max, E_pi_max) ==
            q_gamma_full(E, J_p, n_H, E_p_max, E_pi_max));
}
