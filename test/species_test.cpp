#include <catch2/catch_test_macros.hpp>
#include <kaspectra/species.hpp>
#include <kaspectra/pp/source.hpp>
#include <kaspectra/pgamma/source.hpp>

using kaspectra::Species;

namespace {
    constexpr Species all_species[] = {
        Species::Gamma, Species::Positron, Species::Electron,
        Species::NuMu, Species::NuMuBar, Species::NuE, Species::NuEBar
    };
}

TEST_CASE("Species enumerators are all distinct", "[species]") {
    for (std::size_t i = 0; i < 7; ++i) {
        for (std::size_t j = i + 1; j < 7; ++j) {
            REQUIRE(all_species[i] != all_species[j]);
        }
    }
}

TEST_CASE("pp::F_species dispatches every Species value without throwing", "[species][pp]") {
    double x = 0.1, E = 1000.0;
    for (Species s : all_species) {
        REQUIRE_NOTHROW(kaspectra::pp::F_species(s, x, E));
    }
}

TEST_CASE("pgamma::Phi_species dispatches every Species value without throwing", "[species][pgamma]") {
    double x = 0.5, eta = 10.0;
    for (Species s : all_species) {
        REQUIRE_NOTHROW(kaspectra::pgamma::Phi_species(s, x, eta));
    }
}

TEST_CASE("pp aliases Positron/Electron and NuMuBar/NuE to their pp-family partners", "[species][pp]") {
    double x = 0.1, E = 1000.0;
    REQUIRE(kaspectra::pp::F_species(Species::Positron, x, E) ==
            kaspectra::pp::F_species(Species::Electron, x, E));
    REQUIRE(kaspectra::pp::F_species(Species::NuMuBar, x, E) ==
            kaspectra::pp::F_species(Species::NuMu, x, E));
    REQUIRE(kaspectra::pp::F_species(Species::NuEBar, x, E) ==
            kaspectra::pp::F_species(Species::NuE, x, E));
}

TEST_CASE("pgamma does NOT alias Positron/Electron (genuinely distinct tables)", "[species][pgamma]") {
    double x = 0.5, eta = 10.0;
    REQUIRE(kaspectra::pgamma::Phi_species(Species::Positron, x, eta) !=
            kaspectra::pgamma::Phi_species(Species::Electron, x, eta));
}
