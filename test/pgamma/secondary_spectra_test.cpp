#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/pgamma/secondary_spectra.hpp>
#include <kaspectra/pgamma/kinematics.hpp>
#include <kaspectra/constants.hpp>

using namespace kaspectra::pgamma;
using kaspectra::constants::eta_0;
using kaspectra::constants::r;

TEST_CASE("Phi_gamma golden values, independently recomputed", "[pgamma][secondary_spectra]") {
    REQUIRE(Phi_gamma(0.1, 1.5 * eta_0) == Catch::Approx(4.390337440627435e-18).epsilon(1e-6));
    REQUIRE(Phi_gamma(0.1, 5.0 * eta_0) == Catch::Approx(2.4331665839001097e-17).epsilon(1e-6));
    REQUIRE(Phi_gamma(0.1, 15.0 * eta_0) == Catch::Approx(2.3198000258590826e-17).epsilon(1e-6));
}

TEST_CASE("Phi_numu golden values across rho's 3-regime x'_+ (rho<2.14, 2.14<rho<10, rho>10)", "[pgamma][secondary_spectra]") {
    REQUIRE(Phi_numu(0.1, 1.5 * eta_0) == Catch::Approx(3.375593176484207e-19).epsilon(1e-6));
    REQUIRE(Phi_numu(0.1, 5.0 * eta_0) == Catch::Approx(1.1806177353516284e-17).epsilon(1e-6));
    REQUIRE(Phi_numu(0.1, 15.0 * eta_0) == Catch::Approx(1.1895813777682183e-17).epsilon(1e-6));
}

TEST_CASE("Phi_e_minus golden values above the two-pion threshold", "[pgamma][secondary_spectra]") {
    REQUIRE(Phi_e_minus(0.1, 5.0 * eta_0) == Catch::Approx(2.068874328229692e-18).epsilon(1e-6));
    REQUIRE(Phi_e_minus(0.1, 15.0 * eta_0) == Catch::Approx(3.1598614417948953e-18).epsilon(1e-6));
}

TEST_CASE("Phi_gamma/e+/numu/numubar/nue are zero below the single-pion threshold", "[pgamma][secondary_spectra]") {
    double below = eta_threshold(1.0) * 0.9;
    REQUIRE(Phi_gamma(0.1, below) == 0.0);
    REQUIRE(Phi_e_plus(0.1, below) == 0.0);
    REQUIRE(Phi_numubar(0.1, below) == 0.0);
    REQUIRE(Phi_nue(0.1, below) == 0.0);
    REQUIRE(Phi_numu(0.1, below) == 0.0);
}

TEST_CASE("Phi_e_minus/nuebar are zero below the two-pion threshold, even above the single-pion one", "[pgamma][secondary_spectra]") {
    double between = eta_threshold(1.0) * 1.5; // above single-pion, below two-pion (2.14 eta_0)
    REQUIRE(between > eta_threshold(1.0));
    REQUIRE(between < eta_threshold(1.0 + r));
    REQUIRE(Phi_e_minus(0.1, between) == 0.0);
    REQUIRE(Phi_nuebar(0.1, between) == 0.0);
}

TEST_CASE("Phi_e_plus and Phi_electron_minus are genuinely distinct (not aliased, unlike pp)", "[pgamma][secondary_spectra]") {
    double eta = 5.0 * eta_0;
    REQUIRE(Phi_e_plus(0.1, eta) != Catch::Approx(Phi_e_minus(0.1, eta)).epsilon(1e-3));
}

TEST_CASE("all 7 Phi_* functions are non-negative and finite across a range of eta/x", "[pgamma][secondary_spectra]") {
    for (double rho : {1.5, 2.0, 3.0, 5.0, 10.0, 15.0, 50.0}) {
        double eta = rho * eta_0;
        for (double x : {0.01, 0.05, 0.1, 0.3, 0.5, 0.7}) {
            for (double v : { Phi_gamma(x, eta), Phi_e_plus(x, eta), Phi_e_minus(x, eta),
                               Phi_numu(x, eta), Phi_numubar(x, eta), Phi_nue(x, eta), Phi_nuebar(x, eta) }) {
                REQUIRE(std::isfinite(v));
                REQUIRE(v >= 0.0);
            }
        }
    }
}
