#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <kaspectra/pp/secondary_spectra.hpp>

using namespace kaspectra::pp;

TEST_CASE("F_gamma/F_e/F_numu match verified golden values at x=0.1, E_p=1000", "[pp][secondary_spectra]") {
    REQUIRE(F_gamma(0.1, 1000.0) == Catch::Approx(5.264899522).epsilon(1e-8));
    REQUIRE(F_e(0.1, 1000.0) == Catch::Approx(1.820619787).epsilon(1e-8));
    REQUIRE(F_numu(0.1, 1000.0) == Catch::Approx(3.069505246).epsilon(1e-8));
}

TEST_CASE("F_gamma/F_e/F_numu1 golden values across x", "[pp][secondary_spectra]") {
    REQUIRE(F_gamma(0.01, 1000.0) == Catch::Approx(176.519731684).epsilon(1e-9));
    REQUIRE(F_e(0.01, 1000.0) == Catch::Approx(88.9823797877).epsilon(1e-9));
    REQUIRE(F_numu1(0.01, 1000.0) == Catch::Approx(90.7568262701).epsilon(1e-9));

    REQUIRE(F_gamma(0.3, 1000.0) == Catch::Approx(0.317513272193).epsilon(1e-9));
    REQUIRE(F_e(0.3, 1000.0) == Catch::Approx(0.0532169200596).epsilon(1e-9));
    REQUIRE(F_numu1(0.3, 1000.0) == Catch::Approx(0.00748931012984).epsilon(1e-9));

    REQUIRE(F_gamma(0.9, 1000.0) == Catch::Approx(5.53071483997e-05).epsilon(1e-9));
}

TEST_CASE("F_numu1 has a hard kinematic cutoff at x=0.427", "[pp][secondary_spectra]") {
    REQUIRE(F_numu1(0.4, 1000.0) == Catch::Approx(1.65563665944e-05).epsilon(1e-9)); // just below cutoff: still nonzero
    REQUIRE(F_numu1(0.427, 1000.0) == 0.0); // at cutoff: zero
    REQUIRE(F_numu1(0.5, 1000.0) == 0.0);   // past cutoff: zero
}

TEST_CASE("all secondary spectra vanish outside (0,1)", "[pp][secondary_spectra]") {
    for (double x : {-1.0, 0.0, 1.0, 2.0}) {
        REQUIRE(F_gamma(x, 1000.0) == 0.0);
        REQUIRE(F_e(x, 1000.0) == 0.0);
        REQUIRE(F_numu1(x, 1000.0) == 0.0);
    }
}

TEST_CASE("F_numu is the sum of F_numu1 and F_e", "[pp][secondary_spectra]") {
    for (double x : {0.05, 0.1, 0.2, 0.4}) {
        REQUIRE(F_numu(x, 1000.0) == Catch::Approx(F_numu1(x, 1000.0) + F_e(x, 1000.0)).epsilon(1e-12));
    }
}

TEST_CASE("F_nue equals F_e per the paper's approximation", "[pp][secondary_spectra]") {
    for (double x : {0.05, 0.1, 0.2, 0.4}) {
        REQUIRE(F_nue(x, 1000.0) == F_e(x, 1000.0));
    }
}

TEST_CASE("secondary spectra are non-negative over their domain", "[pp][secondary_spectra]") {
    for (double x : {0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9}) {
        REQUIRE(F_gamma(x, 1000.0) >= 0.0);
        REQUIRE(F_e(x, 1000.0) >= 0.0);
        REQUIRE(F_numu(x, 1000.0) >= 0.0);
        REQUIRE(F_nue(x, 1000.0) >= 0.0);
    }
}
