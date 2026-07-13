"""KASpectra: cosmic-ray secondary spectra (pp, pgamma, Bethe-Heitler).

Python bindings over the C++17 header-only library. Submodules mirror the
C++ namespaces: kaspectra.io (input spectra/fields), kaspectra.pp,
kaspectra.pgamma (both with species-string ergonomics), kaspectra.bh.
"""
from ._kaspectra import (
    Species,
    constants,
    io,
    bh,
    max_depth_hits,
    reset_max_depth_hits,
    _q_total_species,
)
from . import pp
from . import pgamma
from ._species import coerce as _coerce


def q_total_species(species, E, J_p, n_H, f_ph, E_p_max, E_pi_max, epsilon_max):
    """Combined pp + pgamma (+ bh for e+/e-) spectrum for one species.

    Accepts the Species enum or a string alias ("gamma", "e+", "numu", ...).
    Each channel is evaluated with its own tuned default tolerances.
    """
    return _q_total_species(_coerce(species), E, J_p, n_H, f_ph,
                            E_p_max, E_pi_max, epsilon_max)


__all__ = ["Species", "constants", "io", "pp", "pgamma", "bh",
           "q_total_species", "max_depth_hits", "reset_max_depth_hits"]
