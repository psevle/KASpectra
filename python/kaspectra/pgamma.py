"""Pythonic wrapper over kaspectra._kaspectra.pgamma: adds species-string aliases.
See kaspectra/pgamma/source.hpp for the underlying physics."""
from ._kaspectra import pgamma as _pgamma
from ._species import coerce


def q_species(species, E, J_p, f_ph, E_p_max, epsilon_max, **kwargs):
    """Population-level pgamma (photomeson) secondary differential spectrum q(E).

    species may be a kaspectra.Species member or an alias string
    ("gamma", "e+", "e-", "numu", "numubar", "nue", "nuebar").
    E may be a float or a numpy array (vectorized in C++, GIL released).
    """
    return _pgamma.q_species(coerce(species), E, J_p, f_ph, E_p_max, epsilon_max, **kwargs)
