#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <memory>
#include <vector>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/constants.hpp>

using namespace kaspectra::io;

TEST_CASE("BlackbodyPhotonField matches the closed-form Planck spectrum", "[io][photon_field]") {
    using namespace kaspectra::constants;
    double T = 2.725; // CMB temperature, Kelvin
    BlackbodyPhotonField cmb(T);

    double kT = T * k_boltzmann;
    auto planck = [&](double epsilon) {
        return (epsilon * epsilon) / (pi * pi * hbar_c3 * std::expm1(epsilon / kT));
    };

    for (double eps : {0.5 * kT, kT, 2.0 * kT, 5.0 * kT}) {
        REQUIRE(cmb(eps) == Catch::Approx(planck(eps)).epsilon(1e-9));
    }
}

TEST_CASE("BlackbodyPhotonField golden values, independently recomputed", "[io][photon_field]") {
    // Recomputed directly (Python, double precision) from the same closed-form
    // formula, not copied from any prior report.
    BlackbodyPhotonField cmb(2.725);
    REQUIRE(cmb(1.1741116569475e-13) == Catch::Approx(280221554573821.97).epsilon(1e-6));
    REQUIRE(cmb(2.348223313895e-13) == Catch::Approx(423180132501827.1).epsilon(1e-6));
    REQUIRE(cmb(4.69644662779e-13) == Catch::Approx(455242665322336.94).epsilon(1e-6));
    REQUIRE(cmb(1.1741116569475e-12) == Catch::Approx(123317134011218.5).epsilon(1e-6));
}

TEST_CASE("PowerLawPhotonField matches closed-form power law with exponential cutoff", "[io][photon_field]") {
    PowerLawPhotonField f(2.0, 3.0, /*epsilon_cutoff=*/10.0, /*epsilon_ref=*/5.0);
    REQUIRE(f(1.0) == Catch::Approx(226.20935450898986).epsilon(1e-9));
    REQUIRE(f(5.0) == Catch::Approx(2.0 * std::exp(-5.0 / 10.0)).epsilon(1e-12));
}

TEST_CASE("TabulatedPhotonField reproduces tabulated points exactly", "[io][photon_field]") {
    TabulatedPhotonField f({1.0e-9, 1.0e-6, 1.0e-3}, {2.0e10, 2.0e7, 2.0e4});
    REQUIRE(f(1.0e-6) == Catch::Approx(2.0e7));
}

TEST_CASE("PhotonField implementations are usable through the base pointer", "[io][photon_field]") {
    std::vector<std::unique_ptr<PhotonField>> fields;
    fields.push_back(std::make_unique<BlackbodyPhotonField>(2.725));
    fields.push_back(std::make_unique<PowerLawPhotonField>(1.0, 2.0, 1.0e-6));
    fields.push_back(std::make_unique<TabulatedPhotonField>(
        std::vector<double>{1.0e-9, 1.0e-6}, std::vector<double>{1.0, 0.1}));

    for (const auto& f : fields) {
        REQUIRE((*f)(1.0e-7) >= 0.0);
    }
}

TEST_CASE("CompositePhotonField sums its components and is empty-safe", "[io][photon_field]") {
    CompositePhotonField empty;
    REQUIRE(empty(1.0e-13) == 0.0);
    REQUIRE(empty.size() == 0);

    BlackbodyPhotonField cmb(2.725);
    PowerLawPhotonField ir(1.0e5, 1.5, 1.0e-10);
    CompositePhotonField combo;
    combo.add(cmb);
    combo.add(ir);
    REQUIRE(combo.size() == 2);
    for (double eps : {1.0e-14, 1.0e-13, 1.0e-11, 1.0e-10}) {
        CAPTURE(eps);
        REQUIRE(combo(eps) == Catch::Approx(cmb(eps) + ir(eps)).epsilon(1e-12));
    }
}
