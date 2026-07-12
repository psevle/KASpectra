"""Pythonic wrapper over kaspectra._kaspectra.pp: adds species-string aliases.
See kaspectra/pp/source.hpp for the underlying physics."""
from ._kaspectra import pp as _pp
from ._species import coerce


def q_species(species, E, J_p, n_H, E_p_max, E_pi_max, **kwargs):
    """Population-level pp secondary differential spectrum q(E).

    species may be a kaspectra.Species member or an alias string
    ("gamma", "e+", "e-", "numu", "numubar", "nue", "nuebar").
    E may be a float or a numpy array (vectorized in C++, GIL released).
    """
    return _pp.q_species(coerce(species), E, J_p, n_H, E_p_max, E_pi_max, **kwargs)
