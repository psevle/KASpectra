#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <limits>

#include <kaspectra/bh/spectrum.hpp>
#include <kaspectra/bh/source.hpp>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/constants.hpp>
#include <kaspectra/math/integrate.hpp>

using namespace kaspectra;
using kaspectra::bh::W;
using kaspectra::io::BlackbodyPhotonField;
using kaspectra::io::PowerLawSpectrum;

// Golden values independently recomputed (Python, same closed-form expression
// transcribed directly from Blumenthal 1970, eq.10/11), not copied from this
// header's own output.

TEST_CASE("W is finite and non-negative across near-threshold, mid, and high-k regimes", "[bh][spectrum]") {
    struct Case { double k, E_minus, cos_theta_minus, expected; };
    // (k, E_minus, cos_theta_minus, expected), E_minus at the symmetric midpoint k/2
    const Case cases[] = {
        {2.05, 1.025, -0.9, 5.144954e-02},
        {2.05, 1.025,  0.0, 1.564808e-01},
        {2.05, 1.025,  0.9, 1.396019e-01},
        {5.0,  2.5,   -0.9, 1.117882e-01},
        {5.0,  2.5,    0.0, 9.377224e-01},
        {5.0,  2.5,    0.9, 3.204265e+01},
        {20.0, 10.0,  -0.9, 1.187381e-02},
        {100.0, 50.0,  0.9, 1.232226e+00},
    };
    for (const auto& c : cases) {
        double w = W(c.k, c.E_minus, c.cos_theta_minus);
        CAPTURE(c.k, c.E_minus, c.cos_theta_minus);
        REQUIRE(std::isfinite(w));
        REQUIRE(w >= 0.0);
        REQUIRE(w == Catch::Approx(c.expected).epsilon(1e-5));
    }
}

TEST_CASE("W is finite and non-negative across the full E_minus in [1,k-1] range", "[bh][spectrum]") {
    // Sweep E_minus across 5 fractions of the valid range (including near-boundary
    // corners) at k=2.05 (near threshold), k=5 (mid), k=20 (high), each at
    // cos_theta_minus in {-0.9, 0, 0.9}. This is the domain the eventual lab-frame
    // boost integral will scan over, so it must be well-behaved everywhere in it.
    for (double k : {2.05, 5.0, 20.0}) {
        for (double frac : {0.05, 0.25, 0.5, 0.75, 0.95}) {
            double E_minus = 1.0 + frac * (k - 1.0 - 1.0);
            for (double c : {-0.9, 0.0, 0.9}) {
                double w = W(k, E_minus, c);
                CAPTURE(k, E_minus, c);
                REQUIRE(std::isfinite(w));
                REQUIRE(w >= 0.0);
            }
        }
    }
}

TEST_CASE("W is NaN when E_minus is pushed outside the physical [1,k-1] domain", "[bh][spectrum]") {
    // Precondition is documented but deliberately unguarded (matches
    // kaspectra::pgamma::x_pm's convention for unguarded discriminants):
    // p_minus/p_plus become sqrt of a negative number outside [1,k-1].
    REQUIRE(std::isnan(W(5.0, 0.5, 0.0)));   // E_minus < 1 -> p_minus NaN
    REQUIRE(std::isnan(W(5.0, 4.5, 0.0)));   // E_minus > k-1 -> p_plus NaN
    REQUIRE(std::isnan(W(2.05, 0.9, 0.0)));  // E_minus < 1 -> p_minus NaN
    REQUIRE(std::isnan(W(2.05, 1.06, 0.0))); // E_minus > k-1 -> p_plus NaN
}

TEST_CASE("W is finite at the interior symmetric point for a range of k", "[bh][spectrum]") {
    // E_minus = k/2 (equal energy split) is always inside (1, k-1) for k>2,
    // and is the point used as an anchor across all of the golden-value checks
    // above -- verify it stays well-behaved as k grows.
    for (double k : {2.01, 2.5, 3.0, 10.0, 50.0, 500.0}) {
        double w = W(k, k / 2.0, 0.0);
        CAPTURE(k);
        REQUIRE(std::isfinite(w));
        REQUIRE(w >= 0.0);
    }
}

// --- KA2008 eq.60 kinematic bound helpers ---------------------------------

TEST_CASE("E_minus_lo_bound and omega_lo_bound sit exactly at their physical floor when gamma_p==e_e",
          "[bh][spectrum][bounds]") {
    // (gamma_p^2+e_e^2)/(2 gamma_p e_e) == 1 and (gamma_p+e_e)^2/(2 gamma_p e_e) == 2
    // exactly when gamma_p==e_e (AM-GM equality case) -- matches the physical
    // floors E_minus>=1 and k(=omega)>=2 respectively, verified algebraically
    // during the draft phase, checked here numerically at several scales.
    for (double g : {1.0, 5.0, 1e6, 1e11}) {
        double E_minus_lo = bh::E_minus_lo_bound(g, g);
        double omega_lo = bh::omega_lo_bound(g, g);
        CAPTURE(g);
        REQUIRE(E_minus_lo == Catch::Approx(1.0).epsilon(1e-9));
        REQUIRE(omega_lo == Catch::Approx(2.0).epsilon(1e-9));
    }
}

TEST_CASE("E_minus_lo_bound and omega_lo_bound never fall below their physical floor", "[bh][spectrum][bounds]") {
    // AM-GM guarantees E_minus_lo_bound>=1, omega_lo_bound>=2 for ANY positive
    // (gamma_p, e_e), not just the equal case -- swept across asymmetric ratios.
    for (double g : {1e3, 1e9, 1e11}) {
        for (double ratio : {1e-4, 1e-1, 1.0, 1e1, 1e4}) {
            double e_e = g * ratio;
            double E_minus_lo = bh::E_minus_lo_bound(g, e_e);
            double omega_lo = bh::omega_lo_bound(g, e_e);
            CAPTURE(g, ratio, e_e);
            REQUIRE(E_minus_lo >= 1.0 - 1e-9);
            REQUIRE(omega_lo >= 2.0 - 1e-9);
            REQUIRE(omega_lo == Catch::Approx(E_minus_lo + 1.0).epsilon(1e-9));
        }
    }
}

TEST_CASE("eps_lo_bound matches omega_lo_bound/(2 gamma_p) self-consistently", "[bh][spectrum][bounds]") {
    // eps_lo_bound was derived (independently, during the draft phase) as
    // exactly omega_lo/(2*gamma_p) -- the eps at which the omega window
    // [omega_lo, 2*gamma_p*eps/m_e] first opens. Cross-check the two formulas
    // agree with each other directly, rather than re-deriving a separate
    // golden value.
    double m_e = constants::m_e;
    for (double gamma_p : {1e6, 1e9, 1e11}) {
        for (double E_e : {1.0, 1e4, 1e9}) {
            double e_e = E_e / m_e;
            double eps_lo = bh::eps_lo_bound(gamma_p, E_e, m_e);
            double omega_lo = bh::omega_lo_bound(gamma_p, e_e);
            CAPTURE(gamma_p, E_e);
            REQUIRE(eps_lo == Catch::Approx(omega_lo * m_e / (2.0 * gamma_p)).epsilon(1e-9));
        }
    }
}

TEST_CASE("E_p_min inverts eps_lo_bound self-consistently, and is infinite below the physical floor",
          "[bh][spectrum][bounds]") {
    double m_e = constants::m_e;
    double m_p = constants::m_p;

    // Below the floor epsilon_max*E_e <= m_e^2/4: no gamma_p (however large)
    // reaches threshold -- E_p_min must report +infinity.
    double E_e = 1e6;
    double epsilon_max_too_low = (m_e * m_e / 4.0) / E_e * 0.5;
    REQUIRE(std::isinf(bh::E_p_min(E_e, epsilon_max_too_low)));

    // Above the floor: E_p_min(E_e, epsilon_max) should be exactly the E_p at
    // which eps_lo_bound(E_p/m_p, E_e, m_e) == epsilon_max.
    for (double E_e_test : {1e2, 1e6, 1e10}) {
        for (double epsilon_max : {1e-6, 1e-9, 1e-12}) {
            double floor = (m_e * m_e / 4.0) / E_e_test;
            if (epsilon_max <= floor) continue;
            double E_p_lo = bh::E_p_min(E_e_test, epsilon_max);
            REQUIRE(std::isfinite(E_p_lo));
            double gamma_p_lo = E_p_lo / m_p;
            double eps_lo_check = bh::eps_lo_bound(gamma_p_lo, E_e_test, m_e);
            CAPTURE(E_e_test, epsilon_max, E_p_lo);
            REQUIRE(eps_lo_check == Catch::Approx(epsilon_max).epsilon(1e-6));
        }
    }
}

// --- dN_dEe / q_pair_spectrum ----------------------------------------------

TEST_CASE("dN_dEe returns exactly zero below the photon field's reach", "[bh][spectrum]") {
    BlackbodyPhotonField cmb(2.725);
    // Well below GZK-scale: eps_lo_bound for this (E_p, E_e) sits far above any
    // eps the CMB field has non-negligible density at, or above epsilon_max entirely.
    double gamma_p = 1e6; // E_p ~ 9.4e5 GeV ~ 1e15 eV, below CMB pair-production threshold
    double E_e = 1e4;
    double eps_max = 1e-6;
    double v = bh::dN_dEe(E_e, gamma_p, cmb, eps_max);
    REQUIRE(v == 0.0);

    // Directly-degenerate case: epsilon_max below the analytic eps_lo_bound.
    double gamma_p2 = 1e9;
    double E_e2 = 1e6;
    double eps_lo = bh::eps_lo_bound(gamma_p2, E_e2, constants::m_e);
    REQUIRE(bh::dN_dEe(E_e2, gamma_p2, cmb, eps_lo * 0.5) == 0.0);
}

TEST_CASE("dN_dEe is finite and non-negative at a UHECR-scale point (CMB target)",
          "[bh][spectrum][slow]") {
    // Slow (multi-second): 3-level nested integral at gamma_p~1e11 (~1e20 eV).
    // Confirmed during development to complete in a few seconds with this
    // header's tuned default tolerances (see dN_dEe's own doc comment for why
    // these differ from bh/source.hpp's tighter single-integral defaults).
    BlackbodyPhotonField cmb(2.725);
    double gamma_p = 1e11;
    double E_e = 1e9;
    double eps_max = 1e-6;

    double v = bh::dN_dEe(E_e, gamma_p, cmb, eps_max);
    CAPTURE(v);
    REQUIRE(std::isfinite(v));
    REQUIRE(v >= 0.0);
    REQUIRE(v > 0.0); // this (E_e, gamma_p) combination is well within the reachable window
}

TEST_CASE("q_pair_spectrum is finite, non-negative, and zero below the reachable window",
          "[bh][spectrum][slow]") {
    BlackbodyPhotonField cmb(2.725);
    PowerLawSpectrum J_p(1.0, 2.0);
    double E_e = 1e8;
    double eps_max = 1e-6;

    // Below the reachable window: E_p_max less than the kinematic E_p_min for this E_e.
    double E_p_min_here = bh::E_p_min(E_e, eps_max);
    REQUIRE(bh::q_pair_spectrum(E_e, J_p, cmb, E_p_min_here * 0.5, eps_max) == 0.0);

    // Just above the reachable window: a narrow E_p_max window keeps this smoke
    // test's runtime bounded (each dN_dEe evaluation inside costs real time at
    // this scale) while still exercising the outer integral's own logic.
    double E_p_max = E_p_min_here * 1.5;
    double v = bh::q_pair_spectrum(E_e, J_p, cmb, E_p_max, eps_max);
    CAPTURE(v, E_p_min_here, E_p_max);
    REQUIRE(std::isfinite(v));
    REQUIRE(v >= 0.0);
}

TEST_CASE("dN_dEe_planck returns exactly zero below the photon field's reach", "[bh][spectrum]") {
    // Same analytic eps_lo >= epsilon_max early return as the general path --
    // cheap, no integral evaluated.
    double kT = 2.725 * constants::k_boltzmann;
    REQUIRE(bh::dN_dEe_planck(1e-3, 10.0, kT, 1e-15) == 0.0);
}

TEST_CASE("dN_dEe_planck agrees with dN_dEe_general, and dN_dEe dispatches blackbody fields to it",
          "[bh][spectrum][slow]") {
    // The Planckian path is an EXACT algebraic rearrangement of eq.62 (the eps
    // integral done analytically by parts, eq.66/67, extended to finite
    // epsilon_max) -- so the two paths must agree to within the quadrature
    // tolerances themselves, not just physics-level accuracy. Warm blackbody
    // off-peak band keeps the general path's cost at ~1s/call (same config as
    // spectrum_cache_test.cpp's accuracy test).
    double gamma_p = 1e6;
    double kT = 1e-7; // GeV
    BlackbodyPhotonField warm_field(kT / constants::k_boltzmann);
    double eps_max = 20.0 * kT;
    double peak = gamma_p * constants::m_e;

    for (double frac : {0.05, 0.15}) {
        double E_e = peak * frac;
        double general = bh::dN_dEe_general(E_e, gamma_p, warm_field, eps_max);
        double planck = bh::dN_dEe_planck(E_e, gamma_p, kT, eps_max);
        CAPTURE(E_e, general, planck);
        REQUIRE(general > 0.0);
        REQUIRE(planck == Catch::Approx(general).epsilon(0.05));

        // Auto-dispatch: dN_dEe with a BlackbodyPhotonField IS the planck path.
        // Compare against a planck call fed the field's OWN kT() -- the
        // Kelvin->GeV->Kelvin construction roundtrip shifts kT by an ulp
        // relative to the raw 1e-7 used above, so only this form is bitwise
        // identical to what the dispatcher computes.
        REQUIRE(bh::dN_dEe(E_e, gamma_p, warm_field, eps_max)
                == bh::dN_dEe_planck(E_e, gamma_p, warm_field.kT(), eps_max));
    }
}

TEST_CASE("dN_dEe_planck handles the finite epsilon_max cutoff exactly", "[bh][spectrum][slow]") {
    // The by-parts boundary term makes the finite-epsilon_max form exact, so a
    // cutoff INSIDE the field's occupied range (where it actually bites -- at
    // 4kT a blackbody still has real occupancy) must reproduce the general
    // path's value, not just the epsilon_max -> infinity limit of eq.67.
    double gamma_p = 1e6;
    double kT = 1e-7;
    BlackbodyPhotonField warm_field(kT / constants::k_boltzmann);
    double E_e = gamma_p * constants::m_e * 0.1;

    double eps_max_tight = 4.0 * kT;
    double general = bh::dN_dEe_general(E_e, gamma_p, warm_field, eps_max_tight);
    double planck = bh::dN_dEe_planck(E_e, gamma_p, kT, eps_max_tight);
    CAPTURE(general, planck);
    REQUIRE(general > 0.0);
    REQUIRE(planck == Catch::Approx(general).epsilon(0.05));
}
