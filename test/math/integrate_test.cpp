#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/math/integrate.hpp>

using kaspectra::math::integrate;
using kaspectra::math::integrate_log;

TEST_CASE("integrate matches closed-form for a polynomial", "[math][integrate]") {
    double result = integrate([](double x) { return x * x; }, 0.0, 3.0);
    REQUIRE(result == Catch::Approx(9.0).epsilon(1e-9)); // integral of x^2 from 0 to 3 is 9
}

TEST_CASE("integrate of a degenerate interval is zero", "[math][integrate]") {
    REQUIRE(integrate([](double x) { return 1.0 / x; }, 5.0, 5.0) == 0.0);
}

TEST_CASE("integrate is antisymmetric under bound swap", "[math][integrate]") {
    auto f = [](double x) { return std::sin(x); };
    double fwd = integrate(f, 0.2, 1.7);
    double rev = integrate(f, 1.7, 0.2);
    REQUIRE(fwd == Catch::Approx(-rev).epsilon(1e-9));
}

TEST_CASE("integrate matches closed-form for a transcendental integrand", "[math][integrate]") {
    double result = integrate([](double x) { return std::exp(x); }, 0.0, 1.0);
    REQUIRE(result == Catch::Approx(std::exp(1.0) - 1.0).epsilon(1e-9));
}

TEST_CASE("integrate_log agrees with integrate via change of variables", "[math][integrate]") {
    auto f = [](double E) { return std::pow(E, -2.0); };
    double direct = integrate(f, 1.0, 1000.0);
    double via_log = integrate_log(f, 1.0, 1000.0);
    REQUIRE(via_log == Catch::Approx(direct).epsilon(1e-6));
}

TEST_CASE("integrate_log handles a wide dynamic range without excessive recursion", "[math][integrate]") {
    auto f = [](double E) { return std::pow(E, -2.5); };
    double result = integrate_log(f, 1.0, 1.0e6);
    REQUIRE(std::isfinite(result));
    REQUIRE(result > 0.0);
}
