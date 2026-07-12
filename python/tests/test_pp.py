import math

import numpy as np
import pytest

import kaspectra as ks
from kaspectra import pp


@pytest.fixture(scope="module")
def J_p():
    return ks.io.PowerLawSpectrum(1.0, 2.0)


def test_q_gamma_stitch_is_continuous(J_p):
    # Mirrors test/pp/source_test.cpp's stitch-continuity check at
    # E_delta_approx_max (100 GeV), where the delta-function approximation
    # hands over to the accurate parametrization.
    E0 = ks.constants.E_delta_approx_max
    args = (J_p, 1.0, 1e6, 1e6)
    at = pp.q_species("gamma", E0, *args)
    below = pp.q_species("gamma", E0 - 1e-3, *args)
    above = pp.q_species("gamma", E0 + 1e-3, *args)
    assert math.isclose(below, at, rel_tol=1e-4)
    assert math.isclose(above, at, rel_tol=1e-4)


def test_species_enum_and_string_agree(J_p):
    args = (1e3, J_p, 1.0, 1e6, 1e6)
    assert pp.q_species(ks.Species.Gamma, *args) == pp.q_species("gamma", *args)
    assert pp.q_species(ks.Species.Positron, *args) == pp.q_species("e+", *args)


def test_positron_equals_electron(J_p):
    # pp aliases e+ and e- to the same F_e fit (KAB2006 states e+ ~= e-).
    args = (1e3, J_p, 1.0, 1e6, 1e6)
    assert pp.q_species("e+", *args) == pp.q_species("e-", *args)


def test_unknown_species_string_raises(J_p):
    with pytest.raises(ValueError):
        pp.q_species("proton", 1e3, J_p, 1.0, 1e6, 1e6)


def test_all_species_finite_nonnegative(J_p):
    for s in ("gamma", "e+", "e-", "numu", "numubar", "nue", "nuebar"):
        q = pp.q_species(s, 1e3, J_p, 1.0, 1e6, 1e6)
        assert math.isfinite(q)
        assert q >= 0.0


def test_array_overload_matches_scalar_loop(J_p):
    Es = np.logspace(2.5, 4.5, 5)
    vec = pp.q_species("gamma", Es, J_p, 1.0, 1e6, 1e6)
    scalar = np.array([pp.q_species("gamma", e, J_p, 1.0, 1e6, 1e6) for e in Es])
    np.testing.assert_allclose(vec, scalar)


def test_n_H_scales_linearly(J_p):
    q1 = pp.q_species("gamma", 1e3, J_p, 1.0, 1e6, 1e6)
    q2 = pp.q_species("gamma", 1e3, J_p, 2.0, 1e6, 1e6)
    assert q2 == pytest.approx(2.0 * q1, rel=1e-12)
