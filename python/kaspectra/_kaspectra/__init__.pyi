"""
KASpectra: cosmic-ray secondary spectra (pp, pgamma, Bethe-Heitler). Low-level typed extension -- import kaspectra instead of this module directly.
"""
from __future__ import annotations
import typing
from . import bh
from . import constants
from . import io
from . import pgamma
from . import pp
__all__: list[str] = ['Species', 'bh', 'constants', 'io', 'max_depth_hits', 'pgamma', 'pp', 'reset_max_depth_hits']
class Species:
    """
    Secondary particle species, shared by the pp and pgamma channels.
    
    Members:
    
      Gamma
    
      Positron
    
      Electron
    
      NuMu
    
      NuMuBar
    
      NuE
    
      NuEBar
    """
    Electron: typing.ClassVar[Species]  # value = <Species.Electron: 2>
    Gamma: typing.ClassVar[Species]  # value = <Species.Gamma: 0>
    NuE: typing.ClassVar[Species]  # value = <Species.NuE: 5>
    NuEBar: typing.ClassVar[Species]  # value = <Species.NuEBar: 6>
    NuMu: typing.ClassVar[Species]  # value = <Species.NuMu: 3>
    NuMuBar: typing.ClassVar[Species]  # value = <Species.NuMuBar: 4>
    Positron: typing.ClassVar[Species]  # value = <Species.Positron: 1>
    __members__: typing.ClassVar[dict[str, Species]]  # value = {'Gamma': <Species.Gamma: 0>, 'Positron': <Species.Positron: 1>, 'Electron': <Species.Electron: 2>, 'NuMu': <Species.NuMu: 3>, 'NuMuBar': <Species.NuMuBar: 4>, 'NuE': <Species.NuE: 5>, 'NuEBar': <Species.NuEBar: 6>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
def _q_total_species(species: Species, E: float, J_p: ..., n_H: float, f_ph: ..., E_p_max: float, E_pi_max: float, epsilon_max: float) -> float:
    """
    Sum of all implemented channels for one species; each channel uses its own tuned default tolerances.
    """
def max_depth_hits() -> int:
    """
    Panels that exhausted max_depth WITHOUT meeting tolerance since the last reset -- the signature of a possibly under-converged integral.
    """
def reset_max_depth_hits() -> None:
    """
    Zero the max_depth_hits counter for this thread.
    """
