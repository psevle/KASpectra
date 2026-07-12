import math

import pytest

import kaspectra as ks
from kaspectra import pgamma


@pytest.fixture(scope="module")
def J_p():
    return ks.io.PowerLawSpectrum(1.0, 2.0)


@pytest.fixture(scope="module")
def cmb():
    return ks.io.BlackbodyPhotonField(2.725)


def test_q_gamma_finite_nonnegative(J_p, cmb):
    q = pgamma.q_species("gamma", 1e5, J_p, cmb, 1e11, 1e-9)
    assert math.isfinite(q)
    assert q >= 0.0


def test_species_enum_and_string_agree(J_p, cmb):
    args = (1e5, J_p, cmb, 1e11, 1e-9)
    assert pgamma.q_species(ks.Species.Gamma, *args) == pgamma.q_species("gamma", *args)


def test_returns_zero_when_E_at_or_above_E_p_max(J_p, cmb):
    assert pgamma.q_species("gamma", 1e11, J_p, cmb, 1e11, 1e-9) == 0.0


def test_unknown_species_string_raises(J_p, cmb):
    with pytest.raises(ValueError):
        pgamma.q_species("mu", 1e5, J_p, cmb, 1e11, 1e-9)


@pytest.mark.slow
def test_array_overload_matches_scalar_loop(J_p, cmb):
    import numpy as np

    Es = np.array([1e5, 1e6])
    vec = pgamma.q_species("gamma", Es, J_p, cmb, 1e11, 1e-9)
    scalar = np.array([pgamma.q_species("gamma", e, J_p, cmb, 1e11, 1e-9) for e in Es])
    np.testing.assert_allclose(vec, scalar)
