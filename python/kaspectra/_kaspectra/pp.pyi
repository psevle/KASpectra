"""
Proton-proton secondary spectra (arXiv:astro-ph/0606058)
"""
from __future__ import annotations
import kaspectra._kaspectra
import kaspectra._kaspectra.io
import numpy
import typing
__all__: list[str] = ['q_species']
@typing.overload
def q_species(species: kaspectra._kaspectra.Species, E: float, J_p: kaspectra._kaspectra.io.ProtonSpectrum, n_H: float, E_p_max: float, E_pi_max: float, K_pi: float = 0.17, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 50) -> float:
    """
    Population-level pp secondary differential spectrum q(E). Below E_delta_approx_max (100 GeV), Gamma is stitched to a delta-function approximation; other species are q_accurate extrapolated past its fitted validity range. See kaspectra/pp/source.hpp.
    """
@typing.overload
def q_species(species: kaspectra._kaspectra.Species, E: numpy.ndarray[numpy.float64], J_p: kaspectra._kaspectra.io.ProtonSpectrum, n_H: float, E_p_max: float, E_pi_max: float, K_pi: float = 0.17, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 50) -> numpy.ndarray[numpy.float64]:
    """
    Vectorized form of q_species: E may be a numpy array.
    """
