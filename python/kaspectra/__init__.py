"""KASpectra: cosmic-ray secondary spectra (pp, pgamma, Bethe-Heitler).

Python bindings over the C++17 header-only library. Submodules mirror the
C++ namespaces: kaspectra.io (input spectra/fields), kaspectra.pp,
kaspectra.pgamma (both with species-string ergonomics), kaspectra.bh.
"""
from ._kaspectra import Species, constants, io, bh
from . import pp
from . import pgamma

__all__ = ["Species", "constants", "io", "pp", "pgamma", "bh"]
