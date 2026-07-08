#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/pgamma/phi_tables.hpp>

using namespace kaspectra::pgamma;

TEST_CASE("lookup_table_I reproduces exact tabulated rows", "[pgamma][phi_tables]") {
    auto p = lookup_table_I(1.5);
    REQUIRE(p.s == Catch::Approx(0.219));
    REQUIRE(p.delta == Catch::Approx(0.788));
    REQUIRE(p.B == Catch::Approx(1.60e-17));
}

TEST_CASE("lookup_table_I interpolates linearly between rows", "[pgamma][phi_tables]") {
    // Midpoint of the 1.2/1.3 rows (Linear mode, exact arithmetic midpoint)
    auto p = lookup_table_I(1.25);
    REQUIRE(p.s == Catch::Approx(0.144));
    REQUIRE(p.delta == Catch::Approx(0.645));
    REQUIRE(p.B == Catch::Approx(3.925e-18));
}

TEST_CASE("lookup_table_I clamps outside its tabulated range", "[pgamma][phi_tables]") {
    auto below = lookup_table_I(0.5);
    auto row_1_1 = lookup_table_I(1.1);
    REQUIRE(below.s == Catch::Approx(row_1_1.s));
    REQUIRE(below.delta == Catch::Approx(row_1_1.delta));
    REQUIRE(below.B == Catch::Approx(row_1_1.B));

    auto above = lookup_table_I(200.0);
    auto row_100 = lookup_table_I(100.0);
    REQUIRE(above.s == Catch::Approx(row_100.s));
    REQUIRE(above.B == Catch::Approx(row_100.B));
}

TEST_CASE("lookup_table_II_e_plus reproduces exact tabulated rows", "[pgamma][phi_tables]") {
    auto p = lookup_table_II_e_plus(1.5);
    REQUIRE(p.s == Catch::Approx(0.224));
    REQUIRE(p.delta == Catch::Approx(2.69));
    REQUIRE(p.B == Catch::Approx(5.48e-17));
}

TEST_CASE("lookup_table_II species columns are independently addressable", "[pgamma][phi_tables]") {
    // Below eta/eta_0 = 2.14, e+ and numubar are physically near-identical
    // (same mu+ decay chain) - table columns should match closely there,
    // but numu (different decay leg) should differ substantially.
    auto e_plus = lookup_table_II_e_plus(1.5);
    auto numubar = lookup_table_II_numubar(1.5);
    auto numu = lookup_table_II_numu(1.5);

    REQUIRE(e_plus.B == Catch::Approx(numubar.B).epsilon(0.05));
    REQUIRE(numu.B != Catch::Approx(e_plus.B).epsilon(0.05));
}

TEST_CASE("lookup_table_III_e_minus reproduces exact tabulated rows", "[pgamma][phi_tables]") {
    auto p = lookup_table_III_e_minus(5.0);
    REQUIRE(p.s == Catch::Approx(0.256));
    REQUIRE(p.delta == Catch::Approx(2.39));
    REQUIRE(p.B == Catch::Approx(1.24e-16));
}

TEST_CASE("lookup_table_III clamps below its data floor of eta/eta_0 = 3.0", "[pgamma][phi_tables]") {
    // Table III has no data below 3.0 (a domain distinct from both the
    // kinematic two-pion threshold 2.14 and the psi turn-on at 4.0) -
    // querying there returns the row-3.0 clamp, not zero.
    auto below = lookup_table_III_e_minus(2.5);
    auto row_3 = lookup_table_III_e_minus(3.0);
    REQUIRE(below.B == Catch::Approx(row_3.B));
    REQUIRE(below.s == Catch::Approx(row_3.s));
}

TEST_CASE("all table lookups return positive B across their tabulated range", "[pgamma][phi_tables]") {
    for (double x : {1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0}) {
        REQUIRE(lookup_table_I(x).B > 0.0);
        REQUIRE(lookup_table_II_e_plus(x).B > 0.0);
        REQUIRE(lookup_table_II_numubar(x).B > 0.0);
        REQUIRE(lookup_table_II_nue(x).B > 0.0);
    }
    for (double x : {3.0, 5.0, 10.0, 30.0, 100.0}) {
        REQUIRE(lookup_table_III_e_minus(x).B > 0.0);
        REQUIRE(lookup_table_III_nuebar(x).B > 0.0);
    }
}
