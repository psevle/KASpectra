#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/pp/delta_approx.hpp>

using namespace kaspectra::pp;
using kaspectra::io::PowerLawSpectrum;

TEST_CASE("q_pi is exactly linear in n_tilde", "[pp][delta_approx]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double base = q_pi(1.0, 1.0, J_p, /*n_H=*/1.0);
    REQUIRE(q_pi(1.0, 2.5, J_p, 1.0) == Catch::Approx(2.5 * base).epsilon(1e-12));
    REQUIRE(q_pi(1.0, 0.0, J_p, 1.0) == 0.0);
}

TEST_CASE("q_pi is exactly linear in n_H", "[pp][delta_approx]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double base = q_pi(1.0, 1.0, J_p, /*n_H=*/1.0);
    REQUIRE(q_pi(1.0, 1.0, J_p, 4.0) == Catch::Approx(4.0 * base).epsilon(1e-12));
}

TEST_CASE("q_pi is non-negative for physical inputs", "[pp][delta_approx]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    for (double E_pi : {0.1, 1.0, 10.0, 100.0}) {
        REQUIRE(q_pi(E_pi, 1.0, J_p, 1.0) >= 0.0);
    }
}

TEST_CASE("dN_gamma_dE_gamma_delta_approx is linear in n_tilde", "[pp][delta_approx]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    double base = dN_gamma_dE_gamma_delta_approx(1.0, 1.0, J_p, 1.0, /*E_pi_max=*/1.0e5);
    double scaled = dN_gamma_dE_gamma_delta_approx(1.0, 3.0, J_p, 1.0, 1.0e5);
    REQUIRE(scaled == Catch::Approx(3.0 * base).epsilon(1e-6));
}

TEST_CASE("dN_gamma_dE_gamma_delta_approx returns zero for degenerate ranges", "[pp][delta_approx]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    REQUIRE(dN_gamma_dE_gamma_delta_approx(-1.0, 1.0, J_p, 1.0, 1.0e5) == 0.0); // E_gamma <= 0
    REQUIRE(dN_gamma_dE_gamma_delta_approx(1.0, 1.0, J_p, 1.0, /*E_pi_max=*/0.01) == 0.0); // E_pi_max <= E_min
}
