"""
Photomeson (pgamma) secondary spectra (KA2008)
"""
from __future__ import annotations
import kaspectra._kaspectra
import kaspectra._kaspectra.io
import numpy
import typing
__all__: list[str] = ['q_species']
@typing.overload
def q_species(species: kaspectra._kaspectra.Species, E: float, J_p: kaspectra._kaspectra.io.ProtonSpectrum, f_ph: kaspectra._kaspectra.io.PhotonField, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 50) -> float:
    """
    Population-level pgamma secondary differential spectrum q(E). NOTE: same name as kaspectra.pp.q_species but a different signature (photon field instead of n_H) -- kept in separate submodules, mirroring the two C++ namespaces, so there is no collision.
    """
@typing.overload
def q_species(species: kaspectra._kaspectra.Species, E: numpy.ndarray[numpy.float64], J_p: kaspectra._kaspectra.io.ProtonSpectrum, f_ph: kaspectra._kaspectra.io.PhotonField, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 50) -> numpy.ndarray[numpy.float64]:
    """
    Vectorized form of q_species: E may be a numpy array.
    """
