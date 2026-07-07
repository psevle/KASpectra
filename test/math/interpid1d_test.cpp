#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/math/interp1d.hpp>

using kaspectra::math::Interpolator1D;
using kaspectra::math::InterpMode;

TEST_CASE("Interpolator1d rejects malformed input", "[math][interp1d]") {
    REQUIRE_THROWS_AS(Interpolator1D({1.0}, {1.0}), std::invalid_argument); // <2 points
    REQUIRE_THROWS_AS(Interpolator1D({1.0, 2.0}, {1.0}), std::invalid_argument); // mismatched sizes
    REQUIRE_THROWS_AS(Interpolator1D({-1.0, 2.0}, {1.0, 2.0}, InterpMode::LogLog), std::invalid_argument); // non-positive x in LogLog mode
    REQUIRE_NOTHROW(Interpolator1D({-1.0, 2.0}, {1.0, 2.0}, InterpMode::Linear)); // fine in linear mode
}

TEST_CASE("Interpolator1D clamps outside its domain", "[math][interp1d]") {
    Interpolator1D f({1.0, 10.0, 100.0}, {2.0, 20.0, 200.0});
    REQUIRE(f(0.001) == 2.0); // clamps to front, not extrapolated
    REQUIRE(f(1.0e6) == 200.0); // clamps to back
}

TEST_CASE("Interpolator1D passes exactly through its knows", "[math][interp1d]") {
    Interpolator1D f({1.0, 10.0, 100.0}, {2.0, 20.0, 200.0});
    REQUIRE(f(1.0) == 2.0);
    REQUIRE(f(10.0) == Catch::Approx(20.0).epsilon(1e-12));
    REQUIRE(f(100.0) == 200.0);
}

TEST_CASE("Interpolator1D LogLog preserves geometric-mean property", "[math][interp1d]") {
    Interpolator1D f({1.0, 100.0}, {2.0, 200.0}, InterpMode::LogLog);
    // sqrt(1*100)=10 should map to sqrt(2*200)=20 in genuine log-log interpolation
    REQUIRE(f(10.0) == Catch::Approx(20.0).epsilon(1e-12));
}

TEST_CASE("Interpolator1D Linear mode gives arithmetic mean at midpoint", "[math][interp1d]") {
    Interpolator1D f({0.0, 10.0}, {0.0, 100.0}, InterpMode::Linear);
    REQUIRE(f(5.0) == Catch::Approx(50.0).epsilon(1e-12));
}

TEST_CASE("Interpolator1D operator() and at() agree", "[math][interp1d]") {
    Interpolator1D f({1.0, 10.0}, {2.0, 20.0});
    REQUIRE(f(5.0) == f.at(5.0));
}