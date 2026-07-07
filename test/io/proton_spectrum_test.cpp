#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <memory>
#include <vector>
#include <kaspectra/io/proton_spectrum.hpp>

using namespace kaspectra::io;

TEST_CASE("PowerLawSpectrum matches closed-form power law", "[io][proton_spectrum]") {
    PowerLawSpectrum J(5.0, 2.0, /*E_ref=*/1.0);
    REQUIRE(J(1.0) == Catch::Approx(5.0));
    REQUIRE(J(2.0) == Catch::Approx(5.0 / 4.0)); // (2/1)^-2 = 1/4
}

TEST_CASE("ExpCutoffPowerLawSpectrum reduces to plain power law times exp(-1) at E_cutoff", "[io][proton_spectrum]") {
    PowerLawSpectrum plain(3.0, 2.0, 1.0);
    ExpCutoffPowerLawSpectrum cutoff(3.0, 2.0, /*E_cutoff=*/10.0, 1.0);
    REQUIRE(cutoff(10.0) == Catch::Approx(plain(10.0) * std::exp(-1.0)).epsilon(1e-12));
}

TEST_CASE("BrokenPowerLawSpectrum is continuous at E_break", "[io][proton_spectrum]") {
    BrokenPowerLawSpectrum J(1.0, 2.0, 3.0, /*E_break=*/10.0, /*E_ref=*/1.0);

    // exact closed-form check: both branch formulas must agree at E_break itself
    double at = J(10.0);
    REQUIRE(at == Catch::Approx(1.0 * std::pow(10.0 / 1.0, -2.0)).epsilon(1e-12));

    // finite-difference straddle: belt-and-suspenders continuity check
    double below = J(10.0 - 1e-6);
    double above = J(10.0 + 1e-6);
    REQUIRE(at == Catch::Approx(below).epsilon(1e-4));
    REQUIRE(at == Catch::Approx(above).epsilon(1e-4));
}

TEST_CASE("TabulatedSpectrum reproduces tabulated points exactly", "[io][proton_spectrum]") {
    TabulatedSpectrum J({1.0, 10.0, 100.0}, {2.0, 0.2, 0.02});
    REQUIRE(J(10.0) == Catch::Approx(0.2));
}

TEST_CASE("ProtonSpectrum implementations are usable through the base pointer", "[io][proton_spectrum]") {
    std::vector<std::unique_ptr<ProtonSpectrum>> specs;
    specs.push_back(std::make_unique<PowerLawSpectrum>(1.0, 2.0));
    specs.push_back(std::make_unique<ExpCutoffPowerLawSpectrum>(1.0, 2.0, 100.0));
    specs.push_back(std::make_unique<BrokenPowerLawSpectrum>(1.0, 2.0, 3.0, 10.0));
    specs.push_back(std::make_unique<TabulatedSpectrum>(
        std::vector<double>{1.0, 10.0}, std::vector<double>{1.0, 0.1}));

    for (const auto& s : specs) {
        REQUIRE((*s)(5.0) >= 0.0);
    }
}
