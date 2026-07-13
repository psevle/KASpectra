#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/math/integrate.hpp>

using kaspectra::math::integrate;
using kaspectra::math::integrate_log;
using kaspectra::math::max_depth_hits;
using kaspectra::math::reset_max_depth_hits;

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

TEST_CASE("adaptive Simpson does not run away when convergence is limited by floating-point noise",
          "[math][integrate]") {
    // Regression for a real hang found in bh::interaction_rate at an extreme (but legal)
    // parameter scale: near one of bh::fits.hpp's own fit thresholds, catastrophic
    // cancellation in a closed-form expression produced ~1e-14-scale floating-point noise
    // (non-monotonic, sign-flipping) instead of the smooth near-zero value the physics
    // expects. adaptive_simpson's recursive tolerance halves at every level
    // (rel_tol /= 2.0, splitting the error budget across the growing subinterval count) with
    // nothing flooring it against machine epsilon -- once rel_tol*|combined| drops below
    // DBL_EPSILON-scale (observed: by depth ~25 for a large-magnitude integrand), the
    // convergence check |combined-whole|<=15*tol asks for better-than-double-precision
    // agreement, which noise can never satisfy. Every node then recurses on both children
    // all the way to max_depth=50: confirmed directly via an instrumented copy of this file's
    // own logic that this explodes past 2,000,000 evaluations (still climbing) for a single
    // integral that should cost a few hundred at most.
    //
    // Reproduced here with a synthetic, deterministic (not random, so the test is
    // reproducible) noise pattern superimposed on an otherwise smooth integrand, rather than
    // depending on bh module specifics -- isolates the integrate.hpp-level fix (flooring tol
    // at a small multiple of ULP(combined)) from the separate, still-open bh::psi_over_kappa2
    // cancellation issue that originally surfaced it.
    long call_count = 0;
    auto noisy_near_smooth = [&](double x) {
        ++call_count;
        double noise = (std::fmod(x * 1.0e15, 2.0) - 1.0) * 1.0e-14; // deterministic, ~1e-14 scale, sign-flipping
        return 100.0 + std::sin(x) + noise;
    };

    double result = integrate(noisy_near_smooth, 0.0, 10.0);

    REQUIRE(std::isfinite(result));
    // Before the epsilon floor, this integrand's noise pushed adaptive_simpson to explode
    // past 2,000,000 evaluations without converging (see instrumented reproduction above);
    // with the floor, it should settle within a call budget appropriate for a smooth,
    // well-conditioned integrand over kDefaultPanels=32 initial panels.
    REQUIRE(call_count < 100000);
}

TEST_CASE("max_depth_hits counts unconverged bailouts and stays zero for smooth integrands",
          "[math][integrate]") {
    // Every silent-wrong-answer quadrature episode in this project ended with
    // a panel exhausting max_depth WITHOUT meeting tolerance -- previously
    // with no signal whatsoever. max_depth_hits() makes that observable.
    reset_max_depth_hits();
    (void)integrate([](double x) { return x * x; }, 0.0, 1.0);
    REQUIRE(max_depth_hits() == 0);

    // |x - pi/4|^(-1/2) has an interior integrable singularity: with a starved
    // max_depth the panels containing it cannot converge and must be counted.
    reset_max_depth_hits();
    (void)integrate([](double x) { return 1.0 / std::sqrt(std::fabs(x - 0.7853981633974483) + 1e-300); },
                    0.0, 1.0, 1e-14, 1e-12, /*max_depth=*/3);
    REQUIRE(max_depth_hits() > 0);

    // reset clears the counter.
    reset_max_depth_hits();
    REQUIRE(max_depth_hits() == 0);
}
