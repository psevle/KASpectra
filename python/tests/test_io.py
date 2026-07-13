import math

import pytest

import kaspectra as ks


def test_power_law_spectrum_closed_form():
    J_p = ks.io.PowerLawSpectrum(2.0, 2.5)
    assert J_p(1.0) == pytest.approx(2.0)
    assert J_p(10.0) == pytest.approx(2.0 * 10.0 ** -2.5)


def test_power_law_spectrum_E_ref_default_and_override():
    assert ks.io.PowerLawSpectrum(1.0, 2.0)(100.0) == pytest.approx(1e-4)
    assert ks.io.PowerLawSpectrum(1.0, 2.0, E_ref=10.0)(100.0) == pytest.approx(1e-2)


def test_exp_cutoff_power_law_spectrum():
    J_p = ks.io.ExpCutoffPowerLawSpectrum(1.0, 2.0, 100.0)
    assert J_p(100.0) == pytest.approx(1e-4 * math.exp(-1.0))


def test_broken_power_law_is_continuous_at_break():
    J_p = ks.io.BrokenPowerLawSpectrum(1.0, 2.0, 3.0, 100.0)
    below = J_p(100.0 * (1 - 1e-9))
    above = J_p(100.0 * (1 + 1e-9))
    assert below == pytest.approx(above, rel=1e-6)


def test_tabulated_spectrum_reproduces_power_law():
    E = [1.0, 10.0, 100.0, 1000.0]
    J = [e ** -2.0 for e in E]
    J_p = ks.io.TabulatedSpectrum(E, J)
    # log-log interpolation is exact for a pure power law
    assert J_p(30.0) == pytest.approx(30.0 ** -2.0, rel=1e-12)


def test_blackbody_photon_field_planck_form():
    cmb = ks.io.BlackbodyPhotonField(2.725)
    kT = 2.725 * ks.constants.k_boltzmann
    eps = kT  # near the peak
    expected = eps * eps / (
        ks.constants.pi ** 2 * ks.constants.hbar_c3 * math.expm1(eps / kT)
    )
    assert cmb(eps) == pytest.approx(expected, rel=1e-12)


def test_power_law_photon_field():
    f_ph = ks.io.PowerLawPhotonField(1.0, 2.0, 1e-9)
    eps = 1e-9
    assert f_ph(eps) == pytest.approx((eps / 1.0) ** -2.0 * math.exp(-1.0), rel=1e-12)


def test_tabulated_photon_field_reproduces_power_law():
    eps = [1e-12, 1e-11, 1e-10, 1e-9]
    f = [e ** -1.5 for e in eps]
    f_ph = ks.io.TabulatedPhotonField(eps, f)
    assert f_ph(3e-11) == pytest.approx(3e-11 ** -1.5, rel=1e-9)


def test_abstract_bases_subclassable_but_not_callable_raw():
    # Since the trampolines landed, the bases ARE instantiable (that is what
    # makes Python subclassing work) -- but calling the pure-virtual
    # __call__ on a non-overriding instance must raise, not crash.
    class MySpec(ks.io.ProtonSpectrum):
        def __call__(self, E_p):
            return 2.0 * E_p

    assert MySpec()(3.0) == 6.0

    with pytest.raises(RuntimeError):
        ks.io.ProtonSpectrum()(1.0)
    with pytest.raises(RuntimeError):
        ks.io.PhotonField()(1.0)
