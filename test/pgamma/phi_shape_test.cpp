#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/pgamma/phi_shape.hpp>

using namespace kaspectra::pgamma;

namespace {
    constexpr double B = 1.0e-16, s = 0.2, delta = 2.0, psi = 2.5;
    constexpr double x_minus = 0.01, x_plus = 0.5;
}

TEST_CASE("phi_shape golden values, independently recomputed", "[pgamma][phi_shape]") {
    REQUIRE(phi_shape(0.005, B, s, delta, psi, x_minus, x_plus) ==
            Catch::Approx(4.000033721822121e-17).epsilon(1e-9));
    REQUIRE(phi_shape(0.1, B, s, delta, psi, x_minus, x_plus) ==
            Catch::Approx(1.2254338322221142e-17).epsilon(1e-9));
    REQUIRE(phi_shape(0.49, B, s, delta, psi, x_minus, x_plus) ==
            Catch::Approx(2.87637420163832e-22).epsilon(1e-6));
}

TEST_CASE("phi_shape is continuous at x_minus", "[pgamma][phi_shape]") {
    double below = phi_shape(x_minus - 1e-9, B, s, delta, psi, x_minus, x_plus);
    double at = phi_shape(x_minus, B, s, delta, psi, x_minus, x_plus);
    double expected_below_branch = B * std::pow(std::log(2.0), psi);

    REQUIRE(below == Catch::Approx(expected_below_branch).epsilon(1e-9));
    REQUIRE(at == Catch::Approx(expected_below_branch).epsilon(1e-6));
}

TEST_CASE("phi_shape is zero at and above x_plus", "[pgamma][phi_shape]") {
    REQUIRE(phi_shape(x_plus, B, s, delta, psi, x_minus, x_plus) == 0.0);
    REQUIRE(phi_shape(x_plus + 0.1, B, s, delta, psi, x_minus, x_plus) == 0.0);
}

TEST_CASE("phi_shape returns zero for degenerate/invalid inputs", "[pgamma][phi_shape]") {
    REQUIRE(phi_shape(0.0, B, s, delta, psi, x_minus, x_plus) == 0.0);
    REQUIRE(phi_shape(-1.0, B, s, delta, psi, x_minus, x_plus) == 0.0);
    REQUIRE(phi_shape(0.1, B, s, delta, psi, /*x_minus=*/0.5, /*x_plus=*/0.1) == 0.0); // x_minus >= x_plus
    REQUIRE(phi_shape(0.1, B, s, delta, psi, /*x_minus=*/0.0, x_plus) == 0.0); // x_minus <= 0
}

TEST_CASE("phi_shape is non-negative across its domain", "[pgamma][phi_shape]") {
    for (double x : {0.001, 0.01, 0.05, 0.1, 0.3, 0.49}) {
        REQUIRE(phi_shape(x, B, s, delta, psi, x_minus, x_plus) >= 0.0);
    }
}
