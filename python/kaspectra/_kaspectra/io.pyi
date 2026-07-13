"""
Proton spectrum and photon field input types
"""
from __future__ import annotations
__all__: list[str] = ['BlackbodyPhotonField', 'BrokenPowerLawSpectrum', 'CompositePhotonField', 'ExpCutoffPowerLawSpectrum', 'PhotonField', 'PowerLawPhotonField', 'PowerLawSpectrum', 'ProtonSpectrum', 'TabulatedPhotonField', 'TabulatedSpectrum']
class BlackbodyPhotonField(PhotonField):
    """
    Isotropic Planckian photon field at temperature T_kelvin [K] (e.g. 2.725 for the CMB).
    """
    def __init__(self, T_kelvin: float) -> None:
        ...
    def kT(self) -> float:
        """
        kT in GeV.
        """
class BrokenPowerLawSpectrum(ProtonSpectrum):
    """
    Two power laws joined continuously at E_break; norm is defined at E_ref < E_break.
    """
    def __init__(self, norm: float, index_lo: float, index_hi: float, E_break: float, E_ref: float = 1.0) -> None:
        ...
class CompositePhotonField(PhotonField):
    """
    Sum of photon fields (e.g. CMB + IR + starlight). add() stores a reference; keep_alive ties each added field's lifetime to the composite.
    """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        ...
    def add(self, field: PhotonField) -> None:
        ...
class ExpCutoffPowerLawSpectrum(ProtonSpectrum):
    """
    J_p(E_p) = norm * (E_p/E_ref)^-index * exp(-E_p/E_cutoff)
    """
    def __init__(self, norm: float, index: float, E_cutoff: float, E_ref: float = 1.0) -> None:
        ...
class PhotonField:
    """
    f_ph(epsilon): differential photon number density [GeV^-1 cm^-3]. Subclassable from Python: implement __call__(self, epsilon) -> float.
    """
    def __call__(self, epsilon: float) -> float:
        ...
    def __init__(self) -> None:
        ...
class PowerLawPhotonField(PhotonField):
    """
    f_ph(epsilon) = norm * (epsilon/epsilon_ref)^-index * exp(-epsilon/epsilon_cutoff)
    """
    def __init__(self, norm: float, index: float, epsilon_cutoff: float, epsilon_ref: float = 1.0) -> None:
        ...
class PowerLawSpectrum(ProtonSpectrum):
    """
    J_p(E_p) = norm * (E_p/E_ref)^-index
    """
    def __init__(self, norm: float, index: float, E_ref: float = 1.0) -> None:
        ...
class ProtonSpectrum:
    """
    J_p(E_p): differential proton flux/number density [GeV]. Construct one of the concrete subclasses below, use TabulatedSpectrum for arbitrary data, or subclass this in Python and implement __call__(self, E_p) -> float. No unit/normalization convention is enforced -- caller's responsibility.
    """
    def __call__(self, E_p: float) -> float:
        ...
    def __init__(self) -> None:
        ...
class TabulatedPhotonField(PhotonField):
    """
    Log-log interpolated photon field over tabulated (epsilon, f_ph) points.
    """
    def __init__(self, epsilon: list[float], f_ph: list[float]) -> None:
        ...
class TabulatedSpectrum(ProtonSpectrum):
    """
    Log-log interpolated proton spectrum over tabulated (E, J) points.
    """
    def __init__(self, E: list[float], J: list[float]) -> None:
        ...
