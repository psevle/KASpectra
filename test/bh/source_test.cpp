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

TEST_CASE("interaction_rate does not hang at an extreme photon-field normalization", "[bh][source]") {
    // Regression for a real (confirmed via a standalone timed probe: >20M adaptive_simpson
    // evaluations and still climbing before this fix) near-hang. Near kappa==2,
    // psi_over_kappa2's [2,4) branch sums O(10-300)-magnitude terms that must cancel almost
    // exactly; that cancellation is not Lipschitz at the ULP level (confirmed directly: the
    // returned value fluctuates ~1e-14, sign-flipping, as kappa creeps toward 2 by successive
    // ULPs). At a photon-field normalization extreme enough to make the integrand's absolute
    // magnitude huge everywhere (including in this noisy near-threshold region), the old
    // max_depth=50 default let adaptive_simpson chase this unresolvable noise almost all the
    // way to max_depth on every affected branch. Confirmed empirically (see bh/source.hpp's
    // doc comment on interaction_rate) that the answer is bit-identical for max_depth in
    // [15,25] -- this test's real assertion is that the call returns promptly at all, not
    // its specific numeric value.
    const double gamma_p = 10.0;
    const double E_p = gamma_p * constants::m_p;
    const double kT = 1e-4; // GeV -- deliberately unphysical/hot, not a real CMB-like scale
    const double T_kelvin = kT / constants::k_boltzmann;
    const double eps_max = 20.0 * constants::m_e / (2.0 * gamma_p);

    BlackbodyPhotonField f_ph(T_kelvin);
    double rate = bh::interaction_rate(E_p, f_ph, eps_max);
    double loss = bh::energy_loss_rate(E_p, f_ph, eps_max);

    CAPTURE(rate, loss);
    REQUIRE(std::isfinite(rate));
    REQUIRE(std::isfinite(loss));
    REQUIRE(rate >= 0.0);
    REQUIRE(loss >= 0.0);
}

TEST_CASE("nucleus wrappers reduce to the proton functions at Z=A=1 and scale as Z^2 at fixed gamma",
          "[bh][source]") {
    BlackbodyPhotonField cmb(2.725);
    double E_p = 1e10, eps_max = 1e-6;

    REQUIRE(bh::nucleus_interaction_rate(E_p, 1.0, 1.0, cmb, eps_max)
            == bh::interaction_rate(E_p, cmb, eps_max));
    REQUIRE(bh::nucleus_energy_loss_rate(E_p, 1.0, 1.0, cmb, eps_max)
            == bh::energy_loss_rate(E_p, cmb, eps_max));

    // Helium at the SAME Lorentz factor as a 1e10 GeV proton (E_He = 4 E_p):
    // rate and -dE/dt are exactly Z^2 = 4x the proton's.
    REQUIRE(bh::nucleus_interaction_rate(4.0 * E_p, 2.0, 4.0, cmb, eps_max)
            == Catch::Approx(4.0 * bh::interaction_rate(E_p, cmb, eps_max)).epsilon(1e-12));
    REQUIRE(bh::nucleus_energy_loss_rate(4.0 * E_p, 2.0, 4.0, cmb, eps_max)
            == Catch::Approx(4.0 * bh::energy_loss_rate(E_p, cmb, eps_max)).epsilon(1e-12));
}

TEST_CASE("loss_timescale and interaction_length wrap the rates, infinite below threshold",
          "[bh][source]") {
    BlackbodyPhotonField cmb(2.725);
    double E_p = 1e10, eps_max = 1e-6;

    REQUIRE(bh::loss_timescale(E_p, cmb, eps_max)
            == Catch::Approx(E_p / bh::energy_loss_rate(E_p, cmb, eps_max)).epsilon(1e-12));
    REQUIRE(bh::interaction_length(E_p, cmb, eps_max)
            == Catch::Approx(constants::c_light / bh::interaction_rate(E_p, cmb, eps_max)).epsilon(1e-12));

    // Below the kappa==2 threshold both rates are exactly zero -> +infinity.
    double eps_below = constants::m_p * constants::m_e / E_p * 0.5;
    REQUIRE(std::isinf(bh::loss_timescale(E_p, cmb, eps_below)));
    REQUIRE(std::isinf(bh::interaction_length(E_p, cmb, eps_below)));
}
