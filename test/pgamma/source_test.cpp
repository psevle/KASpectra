#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/pgamma/source.hpp>
#include <kaspectra/pgamma/kinematics.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/constants.hpp>

using namespace kaspectra::pgamma;
using kaspectra::Species;
using kaspectra::io::PowerLawSpectrum;
using kaspectra::io::BlackbodyPhotonField;

TEST_CASE("Phi_species dispatches to the matching per-species function", "[pgamma][source]") {
    double x = 0.1, eta = 5.0 * kaspectra::constants::eta_0;
    REQUIRE(Phi_species(Species::Gamma, x, eta) == Phi_gamma(x, eta));
    REQUIRE(Phi_species(Species::Positron, x, eta) == Phi_e_plus(x, eta));
    REQUIRE(Phi_species(Species::Electron, x, eta) == Phi_e_minus(x, eta));
    REQUIRE(Phi_species(Species::NuMu, x, eta) == Phi_numu(x, eta));
    REQUIRE(Phi_species(Species::NuMuBar, x, eta) == Phi_numubar(x, eta));
    REQUIRE(Phi_species(Species::NuE, x, eta) == Phi_nue(x, eta));
    REQUIRE(Phi_species(Species::NuEBar, x, eta) == Phi_nuebar(x, eta));
}

TEST_CASE("q_species returns zero for a degenerate integration window (E >= E_p_max)", "[pgamma][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    REQUIRE(q_species(Species::Gamma, 1000.0, J_p, cmb, /*E_p_max=*/1000.0, 1e-9) == 0.0);
    REQUIRE(q_species(Species::Gamma, 1000.0, J_p, cmb, /*E_p_max=*/500.0, 1e-9) == 0.0);
}

TEST_CASE("q_species returns zero when no epsilon window clears the species threshold", "[pgamma][source]") {
    // epsilon_max far too small for eta(eps,Ep) to ever reach eta_threshold,
    // for any Ep up to E_p_max - every outer-integral sample should hit the
    // eps_lo >= epsilon_max early return.
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    double q = q_species(Species::Gamma, 10.0, J_p, cmb, 1.0e6, /*epsilon_max=*/1.0e-20);
    REQUIRE(q == 0.0);
}

TEST_CASE("q_species is finite and non-negative for all 7 species (smoke test, all run to completion)", "[pgamma][source]") {
    // This is the double-integral analog of pp's q_accurate hang bug: the
    // real risk here is a NaN or runaway adaptive-Simpson recursion, not just
    // a wrong value. Every species below must return a finite, non-negative
    // number without hanging.
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    for (Species s : { Species::Gamma, Species::Positron, Species::Electron,
                        Species::NuMu, Species::NuMuBar, Species::NuE, Species::NuEBar }) {
        double q = q_species(s, 1.0e5, J_p, cmb, 1.0e11, 1.0e-9);
        REQUIRE(std::isfinite(q));
        REQUIRE(q >= 0.0);
    }
}

TEST_CASE("q_species stays finite right at each species' own kinematic threshold window", "[pgamma][source]") {
    // Degenerate-edge case: E_p_max chosen so the outer integral's upper
    // bound sits just above E, and epsilon_max chosen so the inner window is
    // extremely tight - exactly the kind of boundary that produced a
    // NaN/hang in pp's q_accurate when eta sat at/near a threshold.
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    for (Species s : { Species::Gamma, Species::Electron, Species::NuMu }) {
        double q = q_species(s, 100.0, J_p, cmb, 101.0, 1.0e-9);
        REQUIRE(std::isfinite(q));
        REQUIRE(q >= 0.0);
    }
}

TEST_CASE("q_species is monotonically non-decreasing in E_p_max (more proton energy budget, at least as much rate)", "[pgamma][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    double q_narrow = q_species(Species::Gamma, 1.0e4, J_p, cmb, 1.0e6, 1.0e-9);
    double q_wide = q_species(Species::Gamma, 1.0e4, J_p, cmb, 1.0e9, 1.0e-9);
    REQUIRE(q_wide >= q_narrow);
}
