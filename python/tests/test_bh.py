import math

import numpy as np
import pytest

import kaspectra as ks
from kaspectra import bh


@pytest.fixture(scope="module")
def cmb():
    return ks.io.BlackbodyPhotonField(2.725)


@pytest.fixture(scope="module")
def J_p():
    return ks.io.PowerLawSpectrum(1.0, 2.0)


# ---- rates (single integrals: fast) ----------------------------------------


def test_rates_finite_nonnegative_across_uhecr_energies(cmb):
    for E_p in (1e9, 1e10, 1e11):
        rate = bh.interaction_rate(E_p, cmb, 1e-6)
        loss = bh.energy_loss_rate(E_p, cmb, 1e-6)
        assert math.isfinite(rate) and rate >= 0.0
        assert math.isfinite(loss) and loss >= 0.0


def test_rates_zero_below_threshold(cmb):
    # epsilon_max below the kappa==2 threshold: window empty by construction
    E_p = 1e9
    eps_lo = ks.constants.m_p * ks.constants.m_e / E_p
    assert bh.interaction_rate(E_p, cmb, eps_lo * 0.5) == 0.0
    assert bh.energy_loss_rate(E_p, cmb, eps_lo * 0.5) == 0.0


def test_population_rates_finite_positive(J_p, cmb):
    rate = bh.q_pair_rate(J_p, cmb, 1e13, 1e-6)
    loss = bh.q_pair_energy_loss(J_p, cmb, 1e13, 1e-6)
    assert math.isfinite(rate) and rate > 0.0
    assert math.isfinite(loss) and loss > 0.0


# ---- kinematic helpers (closed-form: fast) ----------------------------------


def test_E_p_min_inf_when_unreachable():
    # epsilon_max * E_e <= m_e^2/4 -> never reachable
    assert math.isinf(bh.E_p_min(1e-3, 1e-15))


def test_E_p_min_finite_and_consistent_with_window():
    E_e, eps_max = 1e5, 1e-6
    E_p_lo = bh.E_p_min(E_e, eps_max)
    assert math.isfinite(E_p_lo) and E_p_lo > 0.0
    # At gamma_p just above that minimum, E_e must fall inside E_e_window
    gamma_lo = E_p_lo / ks.constants.m_p
    w = bh.E_e_window(gamma_lo * 1.01, eps_max)
    assert w.reachable
    assert w.lo <= E_e <= w.hi


def test_E_e_window_unreachable_below_floor():
    gamma_p = 10.0
    w = bh.E_e_window(gamma_p, 0.5 * ks.constants.m_e / gamma_p)
    assert not w.reachable
    assert w.lo == 0.0 and w.hi == 0.0


def test_suggested_n_points_monotonic_and_clamped():
    assert bh.suggested_n_points(1.0, 1.0) == 4  # degenerate span clamps
    a = bh.suggested_n_points(1.0, 1e2)
    b = bh.suggested_n_points(1.0, 1e4)
    assert 4 <= a < b


# ---- dN_dEe cheap paths -----------------------------------------------------


def test_dN_dEe_zero_below_kinematic_threshold(cmb):
    # eps_lo_bound >= epsilon_max -> exact analytic 0, no integral evaluated
    assert bh.dN_dEe(1e-3, 10.0, cmb, 1e-15) == 0.0


def test_empty_table_returns_zero_everywhere():
    t = bh.DNdEeTable()
    assert not t.reachable()
    assert len(t) == 0
    assert t(1e5) == 0.0
    np.testing.assert_array_equal(t(np.array([1.0, 1e5, 1e10])), np.zeros(3))


def test_unreachable_over_gamma_p_table(cmb):
    # E_e never reachable at this epsilon_max: factory returns an empty,
    # GammaP-axis table without ever calling dN_dEe.
    t = bh.DNdEeTable.over_gamma_p(1e-3, cmb, 1e-15, gamma_p_max=1e10)
    assert not t.reachable()
    assert t.axis() == bh.CacheAxis.GammaP
    assert t.fixed_value() == 1e-3
    assert t(1e8) == 0.0


def test_unreachable_over_E_e_table(cmb):
    gamma_p = 10.0
    t = bh.DNdEeTable.over_E_e(gamma_p, cmb, 0.5 * ks.constants.m_e / gamma_p)
    assert not t.reachable()
    assert t.axis() == bh.CacheAxis.ElectronEnergy
    assert t(1.0) == 0.0


def test_cached_spectrum_rejects_wrong_axis_table(J_p, cmb):
    gamma_p = 10.0
    t = bh.DNdEeTable.over_E_e(gamma_p, cmb, 0.5 * ks.constants.m_e / gamma_p)
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(1e5, J_p, t, 1e13, 1e-6)


def test_cached_spectrum_rejects_mismatched_E_e(J_p, cmb):
    # Unreachable GammaP-axis table: built without any dN_dEe call, but still
    # carries fixed_value()==1e-3, so querying at a different E_e must raise.
    t = bh.DNdEeTable.over_gamma_p(1e-3, cmb, 1e-15, gamma_p_max=1e10)
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(2e5, J_p, t, 1e13, 1e-6)


def test_rates_array_overloads_match_scalar_loop(cmb):
    Es = np.array([1e9, 1e10, 1e11])
    np.testing.assert_allclose(
        bh.interaction_rate_array(Es, cmb, 1e-6),
        np.array([bh.interaction_rate(e, cmb, 1e-6) for e in Es]))
    np.testing.assert_allclose(
        bh.energy_loss_rate_array(Es, cmb, 1e-6),
        np.array([bh.energy_loss_rate(e, cmb, 1e-6) for e in Es]))


def test_dN_dEe_planck_zero_below_reach():
    kT = 2.725 * ks.constants.k_boltzmann
    assert bh.dN_dEe_planck(1e-3, 10.0, kT, 1e-15) == 0.0


def test_table2d_unbuilt_and_degenerate(cmb):
    t = bh.DNdEeTable2D()
    assert not t.built()
    assert t(1e5, 1e9) == 0.0
    assert not bh.DNdEeTable2D.build(cmb, 1e-6, 1e13, 0.0, 1e6).built()
    assert not bh.DNdEeTable2D.build(cmb, 1e-6, 1e13, 1e6, 1e5).built()


def test_table2d_unreachable_band_builds_cheap_and_zero(cmb):
    # epsilon_max*E_e <= m_e^2/4 across the band: every line empty, no dN_dEe
    # call made, metadata intact.
    t = bh.DNdEeTable2D.build(cmb, 1e-15, 1e13, 1.0, 10.0, n_E_e_lines=4)
    assert t.built()
    assert len(t) == 4
    assert t.epsilon_max() == 1e-15
    assert t(3.0, 1e10) == 0.0
    assert t(0.5, 1e10) == 0.0


def test_q_pair_spectrum_cached_2d_guards(J_p, cmb):
    t = bh.DNdEeTable2D.build(cmb, 1e-15, 1e13, 1.0, 10.0, n_E_e_lines=4)
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(3.0, J_p, bh.DNdEeTable2D(), 1e13, 1e-15)
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(3.0, J_p, t, 1e13, 2e-15)   # epsilon_max mismatch
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(0.1, J_p, t, 1e13, 1e-15)   # below band
    with pytest.raises(ValueError):
        bh.q_pair_spectrum_cached(3.0, J_p, t, 1e14, 1e-15)   # beyond gamma range
    assert bh.q_pair_spectrum_cached(3.0, J_p, t, 1e13, 1e-15) == 0.0


def test_nucleus_wrappers_and_observables(cmb):
    E_p = 1e10
    assert bh.nucleus_interaction_rate(E_p, 1.0, 1.0, cmb, 1e-6) == bh.interaction_rate(E_p, cmb, 1e-6)
    assert bh.nucleus_interaction_rate(4 * E_p, 2.0, 4.0, cmb, 1e-6) == pytest.approx(
        4.0 * bh.interaction_rate(E_p, cmb, 1e-6), rel=1e-12)
    assert bh.loss_timescale(E_p, cmb, 1e-6) == pytest.approx(
        E_p / bh.energy_loss_rate(E_p, cmb, 1e-6), rel=1e-12)
    assert math.isinf(bh.interaction_length(E_p, cmb, ks.constants.m_p * ks.constants.m_e / E_p * 0.5))


def test_composite_photon_field_sums(cmb):
    ir = ks.io.PowerLawPhotonField(1e5, 1.5, 1e-10)
    combo = ks.io.CompositePhotonField()
    combo.add(cmb)
    combo.add(ir)
    assert len(combo) == 2
    for eps in (1e-14, 1e-12, 1e-10):
        assert combo(eps) == pytest.approx(cmb(eps) + ir(eps), rel=1e-12)


def test_python_defined_photon_field_works_in_rates():
    kT = 2.725 * ks.constants.k_boltzmann

    class PyPlanck(ks.io.PhotonField):
        def __init__(self):
            super().__init__()

        def __call__(self, eps):
            import math as _m
            x = eps / kT
            if x > 700.0:   # C++ expm1 saturates; Python math raises OverflowError
                return 0.0
            return eps * eps / (ks.constants.pi ** 2 * ks.constants.hbar_c3 * _m.expm1(x))

    ref = ks.io.BlackbodyPhotonField(2.725)
    mine = PyPlanck()
    # Field evaluation agrees...
    assert mine(1e-13) == pytest.approx(ref(1e-13), rel=1e-12)
    # ...and it is usable from a C++ integral loop (rate is a single integral;
    # cheap even through the Python-callback overhead).
    got = bh.interaction_rate(1e10, mine, 1e-6)
    want = bh.interaction_rate(1e10, ref, 1e-6)
    assert got == pytest.approx(want, rel=1e-9)


def test_python_defined_field_in_parallel_table_build():
    # The GIL-release + trampoline combination: a Python-defined field
    # evaluated from DNdEeTable's std::async workers must neither deadlock
    # nor crash. An unreachable band keeps this instant (build makes no
    # dN_dEe calls but still exercises the guarded entry point end to end).
    class Zero(ks.io.PhotonField):
        def __call__(self, eps):
            return 0.0

    t = bh.DNdEeTable.over_gamma_p(1e-3, Zero(), 1e-15, gamma_p_max=1e10)
    assert not t.reachable()


def test_q_total_species_matches_channel_sum(cmb):
    J_p = ks.io.PowerLawSpectrum(1.0, 2.0)
    total = ks.q_total_species("gamma", 1e3, J_p, 1.0, cmb, 1e6, 1e6, 1e-9)
    from kaspectra import pp as _pp, pgamma as _pg
    parts = (_pp.q_species("gamma", 1e3, J_p, 1.0, 1e6, 1e6)
             + _pg.q_species("gamma", 1e3, J_p, cmb, 1e6, 1e-9))
    assert total == pytest.approx(parts, rel=1e-12)
    with pytest.raises(ValueError):
        ks.q_total_species("proton", 1e3, J_p, 1.0, cmb, 1e6, 1e6, 1e-9)


def test_max_depth_hits_diagnostics():
    ks.reset_max_depth_hits()
    assert ks.max_depth_hits() == 0
    # A smooth integrand certifies its tolerance: counter stays 0. (pp's
    # single integral is genuinely smooth; bh's rates legitimately report a
    # panel or two near the kappa==2 threshold where the fit's curvature
    # meets the noise-stagnation bail, so they are not a zero-hit baseline.)
    from kaspectra import pp as _pp
    _pp.q_species("gamma", 1e3, ks.io.PowerLawSpectrum(1.0, 2.0), 1.0, 1e6, 1e6)
    assert ks.max_depth_hits() == 0
    # reset clears whatever the counter holds.
    bh.interaction_rate(1e10, ks.io.BlackbodyPhotonField(2.725), 1e-6)
    ks.reset_max_depth_hits()
    assert ks.max_depth_hits() == 0


# ---- expensive parity checks (opt-in: pytest -m slow) -----------------------


@pytest.mark.slow
def test_table_matches_direct_dN_dEe():
    # Off-peak band against a warm blackbody (kT = 1e-7 GeV), mirroring
    # test/bh/spectrum_cache_test.cpp's accuracy test: calls there cost
    # ~0.5-1.5 s each instead of the tens of seconds near the spectral peak.
    gamma_p = 1e6
    kT = 1e-7  # GeV
    f_ph = ks.io.BlackbodyPhotonField(kT / ks.constants.k_boltzmann)
    eps_max = 20.0 * kT
    peak = gamma_p * ks.constants.m_e
    t = bh.DNdEeTable.over_E_e(gamma_p, f_ph, eps_max,
                               E_e_min=peak * 0.03, E_e_max=peak * 0.3, n_points=8)
    assert t.reachable()
    x = math.sqrt(t.domain_lo() * t.domain_hi())  # off-grid query point
    direct = bh.dN_dEe(x, gamma_p, f_ph, eps_max)
    assert t(x) == pytest.approx(direct, rel=0.2)
