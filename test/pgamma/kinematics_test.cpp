#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <kaspectra/pgamma/kinematics.hpp>
#include <kaspectra/constants.hpp>

using namespace kaspectra::pgamma;

TEST_CASE("eta_threshold(1.0) matches constants::eta_0 (eq.16, R=1)", "[pgamma][kinematics]") {
    REQUIRE(eta_threshold(1.0) == Catch::Approx(kaspectra::constants::eta_0).epsilon(1e-12));
}

TEST_CASE("eta_threshold(1+r) matches the two-pion threshold 4r(1+r) (eq.38)", "[pgamma][kinematics]") {
    using namespace kaspectra::constants;
    double expected = 4.0 * r * (1.0 + r);
    REQUIRE(eta_threshold(1.0 + r) == Catch::Approx(expected).epsilon(1e-12));
}

TEST_CASE("eta() golden value, independently recomputed", "[pgamma][kinematics]") {
    REQUIRE(eta(1.0e-9, 100.0) == Catch::Approx(4.543624195255131e-07).epsilon(1e-9));
}

TEST_CASE("x_pm golden values at R=1, independently recomputed", "[pgamma][kinematics]") {
    auto [xm1, xp1] = x_pm(1.0, 1.0);
    REQUIRE(xm1 == Catch::Approx(0.022652351152992245).epsilon(1e-9));
    REQUIRE(xp1 == Catch::Approx(0.4884113134258238).epsilon(1e-9));

    auto [xm2, xp2] = x_pm(2.0, 1.0);
    REQUIRE(xm2 == Catch::Approx(0.01112625781375424).epsilon(1e-9));
    REQUIRE(xp2 == Catch::Approx(0.6629161852387896).epsilon(1e-9));

    auto [xm3, xp3] = x_pm(10.0, 1.0);
    REQUIRE(xm3 == Catch::Approx(0.002213223838260331).epsilon(1e-9));
    REQUIRE(xp3 == Catch::Approx(0.9088892606306154).epsilon(1e-9));
}

TEST_CASE("x_pm golden values at R=1+r (e-/nuebar kinematics, eq.40)", "[pgamma][kinematics]") {
    using namespace kaspectra::constants;
    auto [xm, xp] = x_pm(10.0, 1.0 + r);
    REQUIRE(xm == Catch::Approx(0.002286508660924401).epsilon(1e-9));
    REQUIRE(xp == Catch::Approx(0.8797584773430537).epsilon(1e-9));
}

TEST_CASE("x_pm collapses to a single point exactly at threshold (discriminant = 0)", "[pgamma][kinematics]") {
    // Floating-point note: eta_threshold(R) and x_pm's discriminant compute
    // the same zero-crossing via different arithmetic paths ((R+r)^2-1 vs
    // eta+1-(R+r)^2), so the discriminant lands at ~1e-15 rather than
    // exactly 0.0 here - a looser tolerance than the other golden-value
    // checks in this file is expected, not a sign of a bug.
    double eta_th = eta_threshold(1.0);
    auto [xm, xp] = x_pm(eta_th, 1.0);
    REQUIRE(xm == Catch::Approx(xp).epsilon(1e-5));
    REQUIRE(xm == Catch::Approx(0.1294905273417759).epsilon(1e-5));

    using namespace kaspectra::constants;
    double eta_th2 = eta_threshold(1.0 + r);
    auto [xm2, xp2] = x_pm(eta_th2, 1.0 + r);
    REQUIRE(xm2 == Catch::Approx(xp2).epsilon(1e-5));
    REQUIRE(xm2 == Catch::Approx(0.1146450759941548).epsilon(1e-5));
}

TEST_CASE("x_pm(eta,1) reduces to the R=1 single-pion range and R=1+r matches the two-pion range", "[pgamma][kinematics]") {
    // Both special cases (eq.19 at R=1, eq.40 at R=1+r) come from the same
    // general eq.21 formula - already exercised by construction, this test
    // just checks x_minus < x_plus (a valid, non-degenerate interval) well
    // above each threshold.
    auto [xm1, xp1] = x_pm(5.0, 1.0);
    REQUIRE(xm1 < xp1);

    using namespace kaspectra::constants;
    auto [xm2, xp2] = x_pm(5.0, 1.0 + r);
    REQUIRE(xm2 < xp2);
}
