#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <kaspectra/bh/fits.hpp>

using kaspectra::bh::psi_over_kappa2;
using kaspectra::bh::phi_over_kappa2;

// Golden values independently recomputed (Python, same closed-form expressions
// transcribed directly from Chodorowski, Zdziarski & Sikora 1992, eq.2.4/2.5/3.14/
// 3.16/3.18), not copied from this header's own output.

TEST_CASE("psi_over_kappa2 is exactly zero below threshold", "[bh][fits]") {
    REQUIRE(psi_over_kappa2(1.5) == 0.0);
    REQUIRE(psi_over_kappa2(1.9999) == 0.0);
}

TEST_CASE("phi_over_kappa2 is exactly zero below threshold", "[bh][fits]") {
    REQUIRE(phi_over_kappa2(1.5) == 0.0);
    REQUIRE(phi_over_kappa2(1.9999) == 0.0);
}

TEST_CASE("psi_over_kappa2 is accurate at machine precision near threshold", "[bh][fits]") {
    // Golden values from a 60-digit mpmath evaluation of the printed eq.2.4.
    // psi vanishes QUARTICALLY at kappa=2 (orders 0-3 of its Taylor expansion
    // are identically zero, a4 = 1/16 exactly), so the printed form loses all
    // significant digits to cancellation by kappa-2 ~ 1e-4; the series/
    // delta-form rewrite in fits.hpp holds ~1e-12 or better everywhere.
    struct Ref { double kappa, ref; };
    const Ref refs[] = {
        {2.0 + 1e-9, 3.2724934249951868e-38},
        {2.0 + 1e-7, 3.2724917697429673e-30},
        {2.0001,     3.271936085001196e-18},
        {2.001,      3.266935223874015e-14},
        {2.5,        0.00095771683458671812},
        {3.5,        0.02468867161574065},
        {3.999,      0.049443787122196545},
    };
    for (const auto& r : refs) {
        CAPTURE(r.kappa);
        REQUIRE(psi_over_kappa2(r.kappa) == Catch::Approx(r.ref).epsilon(1e-10));
    }
    REQUIRE(psi_over_kappa2(2.0) == 0.0);

    // Continuity across the internal series/delta-form switch at kappa=2.25.
    REQUIRE(psi_over_kappa2(2.25 - 1e-12)
            == Catch::Approx(psi_over_kappa2(2.25 + 1e-12)).epsilon(1e-9));
}

TEST_CASE("psi_over_kappa2 matches independently recomputed values", "[bh][fits]") {
    REQUIRE(psi_over_kappa2(2.0) == Catch::Approx(0.0).margin(1e-10)); // eq.2.4, kappa->2+
    REQUIRE(psi_over_kappa2(3.0) == Catch::Approx(8.24180698e-03).epsilon(1e-6));
    REQUIRE(psi_over_kappa2(10.0) == Catch::Approx(5.49410383e-01).epsilon(1e-6));
    REQUIRE(psi_over_kappa2(100.0) == Catch::Approx(3.47582559e+00).epsilon(1e-6));
    REQUIRE(psi_over_kappa2(1000.0) == Catch::Approx(7.01077805e+00).epsilon(1e-6));
    REQUIRE(psi_over_kappa2(10000.0) == Catch::Approx(1.05906670e+01).epsilon(1e-6));
}

TEST_CASE("phi_over_kappa2 matches independently recomputed values", "[bh][fits]") {
    REQUIRE(phi_over_kappa2(2.0) == Catch::Approx(0.0).margin(1e-10)); // eq.3.13, kappa->2+
    REQUIRE(phi_over_kappa2(3.0) == Catch::Approx(1.49033340e-02).epsilon(1e-6));
    REQUIRE(phi_over_kappa2(10.0) == Catch::Approx(6.18333806e-01).epsilon(1e-6));
    REQUIRE(phi_over_kappa2(100.0) == Catch::Approx(1.06623147e+00).epsilon(1e-6));
    REQUIRE(phi_over_kappa2(1000.0) == Catch::Approx(4.56542351e-01).epsilon(1e-6));
    REQUIRE(phi_over_kappa2(10000.0) == Catch::Approx(1.24093878e-01).epsilon(1e-6));
}

TEST_CASE("psi_over_kappa2 is continuous across the kappa=4 seam", "[bh][fits]") {
    // eq.2.4 (kappa<4) and eq.2.5 (kappa>=4) are two independent series fits to the
    // same underlying function; the paper's own stated fractional error (<1.1e-3)
    // bounds how well they need to agree here, not exact equality.
    double below = psi_over_kappa2(3.9999);
    double above = psi_over_kappa2(4.0001);
    REQUIRE(above == Catch::Approx(below).epsilon(1e-3));
}

TEST_CASE("phi_over_kappa2 seam at kappa=25 matches the paper's own documented accuracy", "[bh][fits]") {
    // eq.3.14 (kappa<25) and eq.3.18 (kappa>=25) are independent fits; the paper
    // states relative error up to ~1.3e-3/1.5e-3 for each, so a seam mismatch on
    // that order is expected, not a bug.
    double below = phi_over_kappa2(24.9999);
    double above = phi_over_kappa2(25.0001);
    REQUIRE(above == Catch::Approx(below).epsilon(3e-3));
}

TEST_CASE("psi_over_kappa2 and phi_over_kappa2 rise monotonically from threshold in the low-kappa regime", "[bh][fits]") {
    double prev_psi = psi_over_kappa2(2.0);
    double prev_phi = phi_over_kappa2(2.0);
    for (double k : {2.5, 3.0, 3.5, 4.0}) {
        double psi = psi_over_kappa2(k);
        double phi = phi_over_kappa2(k);
        REQUIRE(psi > prev_psi);
        REQUIRE(phi > prev_phi);
        prev_psi = psi;
        prev_phi = phi;
    }
}

TEST_CASE("phi_over_kappa2 peaks and declines at high kappa, psi_over_kappa2 keeps rising", "[bh][fits]") {
    // Matches Fig.1 (psi/k2 monotonically rising to ~10 by kappa~1e4) and Fig.2
    // (phi/k2 peaking ~1.0-1.2 around kappa~25-100, then declining) in Chodorowski+1992.
    REQUIRE(psi_over_kappa2(10000.0) > psi_over_kappa2(1000.0));
    REQUIRE(psi_over_kappa2(1000.0) > psi_over_kappa2(100.0));

    REQUIRE(phi_over_kappa2(100.0) > phi_over_kappa2(1000.0));
    REQUIRE(phi_over_kappa2(1000.0) > phi_over_kappa2(10000.0));
}
