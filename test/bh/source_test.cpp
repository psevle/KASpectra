#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include <kaspectra/bh/source.hpp>
#include <kaspectra/io/photon_field.hpp>
#include <kaspectra/io/proton_spectrum.hpp>
#include <kaspectra/constants.hpp>

using namespace kaspectra;
using kaspectra::io::BlackbodyPhotonField;
using kaspectra::io::PowerLawSpectrum;

namespace {
    constexpr double kSecPerYear = 3.15576e7;
}

TEST_CASE("eps_min matches the kappa==2 threshold", "[bh][source]") {
    using namespace kaspectra::constants;
    double E_p = 1e10; // GeV
    double eps = bh::eps_min(E_p);
    REQUIRE(bh::kappa(eps, E_p) == Catch::Approx(2.0).epsilon(1e-9));
}

TEST_CASE("interaction_rate and energy_loss_rate are finite and non-negative across UHECR energies", "[bh][source]") {
    BlackbodyPhotonField cmb(2.725);
    double eps_max = 1e-6; // GeV, generous margin past the CMB Wien tail

    for (double E_p : {1e9, 1e10, 1e11, 1e12, 1e13}) {
        double rate = bh::interaction_rate(E_p, cmb, eps_max);
        double loss = bh::energy_loss_rate(E_p, cmb, eps_max);
        CAPTURE(E_p, rate, loss);
        REQUIRE(std::isfinite(rate));
        REQUIRE(std::isfinite(loss));
        REQUIRE(rate >= 0.0);
        REQUIRE(loss >= 0.0);
    }
}

TEST_CASE("interaction_rate and energy_loss_rate return exactly zero for a degenerate window", "[bh][source]") {
    BlackbodyPhotonField cmb(2.725);
    double E_p = 1e9; // GeV
    double eps_lo = bh::eps_min(E_p);

    // epsilon_max below the analytic threshold: window is empty by construction
    REQUIRE(bh::interaction_rate(E_p, cmb, eps_lo * 0.5) == 0.0);
    REQUIRE(bh::energy_loss_rate(E_p, cmb, eps_lo * 0.5) == 0.0);
}

TEST_CASE("interaction_rate and energy_loss_rate never hang or NaN near the threshold edge", "[bh][source]") {
    BlackbodyPhotonField cmb(2.725);
    double E_p = 1e9;
    double eps_lo = bh::eps_min(E_p);

    // epsilon_max sitting just above eps_min: a tight but non-empty window. The true
    // value here is essentially zero (phi/psi_over_kappa2 both vanish as (kappa-2)^n
    // right at threshold), so adaptive Simpson's Richardson-extrapolation correction
    // can land a few ULPs on the negative side of zero (observed: ~-7e-34, vs a
    // physically-scaled result of order 1e-25 or smaller) - allow that noise floor
    // explicitly rather than asserting strict non-negativity on a near-zero integral.
    double rate = bh::interaction_rate(E_p, cmb, eps_lo * 1.0001);
    double loss = bh::energy_loss_rate(E_p, cmb, eps_lo * 1.0001);
    REQUIRE(std::isfinite(rate));
    REQUIRE(std::isfinite(loss));
    REQUIRE(rate >= Catch::Approx(0.0).margin(1e-20));
    REQUIRE(loss >= Catch::Approx(0.0).margin(1e-20));
}

TEST_CASE("fractional energy-loss rate peaks near 2.5e19 eV for CMB, matching KA2008 Fig.11 order of magnitude", "[bh][source][validation]") {
    // Read directly off KA2008 Fig.11 (dashed e+e- pair-production curve): rises
    // from ~1e-11/yr near 1e18 eV, peaks in the 1e-10-1e-9/yr range around a few
    // times 1e19 eV, then declines toward 1e22 eV. Verified during this session by
    // plotting this exact quantity and confirming the shape/magnitude against the
    // actual figure.
    BlackbodyPhotonField cmb(2.725);
    double eps_max = 1e-6;

    auto frac_loss_per_yr = [&](double E_p_eV) {
        double E_p_GeV = E_p_eV * 1e-9;
        return bh::energy_loss_rate(E_p_GeV, cmb, eps_max) / E_p_GeV * kSecPerYear;
    };

    double at_1e18 = frac_loss_per_yr(1e18);
    double at_2p5e19 = frac_loss_per_yr(2.5e19);
    double at_1e22 = frac_loss_per_yr(1e22);

    CAPTURE(at_1e18, at_2p5e19, at_1e22);

    // Peak (near 2.5e19 eV) is well above both the low- and high-energy ends
    REQUIRE(at_2p5e19 > at_1e18);
    REQUIRE(at_2p5e19 > at_1e22);

    // Order-of-magnitude check against Fig.11's dashed curve peak
    REQUIRE(at_2p5e19 > 1e-11);
    REQUIRE(at_2p5e19 < 1e-8);
}

TEST_CASE("q_pair_rate and q_pair_energy_loss are finite and non-negative", "[bh][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    double E_p_max = 1e13;
    double eps_max = 1e-6;

    double rate = bh::q_pair_rate(J_p, cmb, E_p_max, eps_max);
    double loss = bh::q_pair_energy_loss(J_p, cmb, E_p_max, eps_max);

    CAPTURE(rate, loss);
    REQUIRE(std::isfinite(rate));
    REQUIRE(std::isfinite(loss));
    REQUIRE(rate > 0.0);
    REQUIRE(loss > 0.0);
}

TEST_CASE("q_pair_rate and q_pair_energy_loss return exactly zero when E_p_max is below the reachable window", "[bh][source]") {
    PowerLawSpectrum J_p(1.0, 2.0);
    BlackbodyPhotonField cmb(2.725);
    double eps_max = 1e-9; // ~1 eV; E_p_lo = m_p*m_e/eps_max ~ 4.8e5 GeV

    REQUIRE(bh::q_pair_rate(J_p, cmb, 1.0, eps_max) == 0.0);
    REQUIRE(bh::q_pair_energy_loss(J_p, cmb, 1.0, eps_max) == 0.0);
}
