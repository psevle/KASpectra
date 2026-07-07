#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/constants.hpp>
#include <kaspectra/pp/cross_section.hpp>

using kaspectra::pp::sigma_inel;

TEST_CASE("sigma_inel is zero at and below threshold", "[pp][cross_section]") {
    REQUIRE(sigma_inel(1.0) == 0.0);
    REQUIRE(sigma_inel(kaspectra::constants::E_threshold_pp) == 0.0);
}

TEST_CASE("sigma_inel matches KAB06 eq.79 at the L=0 (1 TeV) tie point", "[pp][cross_section]") {
    // At E_p = 1 TeV, L = ln(1) = 0, so sigma = 34.3 mb * bracket^2 (bracket ~ 1)
    REQUIRE(sigma_inel(1000.0) == Catch::Approx(3.43e-26).epsilon(1e-6));
}

TEST_CASE("sigma_inel additional golden values", "[pp][cross_section]") {
    REQUIRE(sigma_inel(1.3) == Catch::Approx(1.653281e-27).epsilon(1e-5));
    REQUIRE(sigma_inel(10.0) == Catch::Approx(3.093047e-26).epsilon(1e-5));
    REQUIRE(sigma_inel(100.0) == Catch::Approx(3.129661e-26).epsilon(1e-5));
    REQUIRE(sigma_inel(1.0e5) == Catch::Approx(4.825962e-26).epsilon(1e-5));
    REQUIRE(sigma_inel(1.0e6) == Catch::Approx(5.921585e-26).epsilon(1e-5));
}

TEST_CASE("sigma_inel is monotonically increasing across decades", "[pp][cross_section]") {
    REQUIRE(sigma_inel(10.0) < sigma_inel(100.0));
    REQUIRE(sigma_inel(100.0) < sigma_inel(1000.0));
    REQUIRE(sigma_inel(1000.0) < sigma_inel(1.0e5));
}

TEST_CASE("sigma_inel is never negative", "[pp][cross_section]") {
    for (double E_p : {1.0, 1.22, 1.3, 10.0, 1000.0, 1.0e5, 1.0e6}) {
        REQUIRE(sigma_inel(E_p) >= 0.0);
    }
}
