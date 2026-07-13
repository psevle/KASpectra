#include <catch2/catch_test_macros.hpp>
#include <chrono>

#include <kaspectra/bh/spectrum.hpp>
#include <kaspectra/io/photon_field.hpp>

using namespace kaspectra;
using clk = std::chrono::steady_clock;

// Wall-clock ceilings on calls with known cost profiles. The bounds are
// DELIBERATELY generous (100-1000x the measured cost on a dev laptop, to
// absorb slow CI machines): the target is the hour-class failure mode, not
// percent-level drift. Two such regressions shipped undetected in this
// project because nothing timed anything: the falsely-converged 70-orders
// under-integration, and the p_minus=0 NaN grind that made a single
// dN_dEe call at the spectral peak consume 4+ CPU-hours.

namespace {
    double seconds_for(void (*f)()) {
        const auto t0 = clk::now();
        f();
        return std::chrono::duration<double>(clk::now() - t0).count();
    }
}

TEST_CASE("dN_dEe at the exact spectral peak completes in bounded time", "[perf]") {
    // Measured ~0.01 s (Planckian path). Regressed form: 4+ CPU-hours.
    const double elapsed = seconds_for([] {
        io::BlackbodyPhotonField warm_field(1e-7 / constants::k_boltzmann);
        volatile double v = bh::dN_dEe(1e6 * constants::m_e, 1e6, warm_field, 2e-6);
        (void)v;
    });
    CAPTURE(elapsed);
    REQUIRE(elapsed < 30.0);
}

TEST_CASE("dN_dEe off-peak against the CMB at UHECR scale completes in bounded time", "[perf]") {
    // Measured ~0.01-0.1 s (Planckian path); the pre-eq.67 general path
    // took ~6 s here and the config is representative of real use.
    const double elapsed = seconds_for([] {
        io::BlackbodyPhotonField cmb(2.725);
        volatile double v = bh::dN_dEe(1e11 * constants::m_e * 0.1, 1e11, cmb, 1e-6);
        (void)v;
    });
    CAPTURE(elapsed);
    REQUIRE(elapsed < 60.0);
}
