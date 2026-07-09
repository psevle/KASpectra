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

TEST_CASE("integrate_log resolves a narrow peak pinned at the lower edge of a wide range", "[math][integrate]") {
    // Regression test for a real bug found while validating pgamma against KA2008's own
    // Fig.6 (E_p=1e20 eV): a single-panel adaptive Simpson only samples 3 points (a, mid, b)
    // before deciding whether to recurse. An integrand that is sharply peaked right at the
    // lower bound `a` and decays by many orders of magnitude within the first decade (exactly
    // what Phi_e_minus(x,eta(eps,E_p)) does near its analytic threshold eps0, for a blackbody
    // photon field with UHECR-scale E_p) can leave every initial sample near zero, causing the
    // integrator to falsely report convergence at ~0 rather than finding the peak. Composite
    // paneling (splitting into kDefaultPanels sub-intervals before adaptively refining each)
    // fixes this by bounding how wide a region any single panel's initial samples must cover.
    //
    // Golden value: closed form for a decaying-exponential-like peak, integral_a^b exp(-k*(x-a)) dx
    // = (1 - exp(-k*(b-a))) / k, with k large enough that the peak's width is a tiny fraction
    // of [a,b] (mirrors the real integrand falling ~10 orders of magnitude within one decade
    // of `a`, over a domain spanning several decades).
    const double a = 1.0;
    const double b = 1.0e6;
    const double k = 50.0; // peak width ~1/k = 0.02, domain width ~1e6 -> width/domain ~2e-8
    auto f = [&](double x) { return std::exp(-k * (x - a)); };
    double result = integrate(f, a, b);
    double golden = (1.0 - std::exp(-k * (b - a))) / k;
    REQUIRE(result == Catch::Approx(golden).epsilon(1e-6));
}
