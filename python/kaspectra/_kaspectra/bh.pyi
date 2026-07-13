"""
Bethe-Heitler pair-production rates and spectra
"""
from __future__ import annotations
import kaspectra._kaspectra.io
import numpy
import typing
__all__: list[str] = ['CacheAxis', 'DNdEeTable', 'DNdEeTable2D', 'E_e_window', 'E_p_min', 'EeWindow', 'dN_dEe', 'dN_dEe_array', 'dN_dEe_fast', 'dN_dEe_general', 'dN_dEe_planck', 'energy_loss_rate', 'energy_loss_rate_array', 'interaction_length', 'interaction_rate', 'interaction_rate_array', 'loss_timescale', 'make_table_for_q_pair_spectrum', 'nucleus_energy_loss_rate', 'nucleus_interaction_rate', 'q_pair_energy_loss', 'q_pair_rate', 'q_pair_spectrum', 'q_pair_spectrum_cached', 'suggested_n_points']
class CacheAxis:
    """
    Members:
    
      GammaP
    
      ElectronEnergy
    """
    ElectronEnergy: typing.ClassVar[CacheAxis]  # value = <CacheAxis.ElectronEnergy: 1>
    GammaP: typing.ClassVar[CacheAxis]  # value = <CacheAxis.GammaP: 0>
    __members__: typing.ClassVar[dict[str, CacheAxis]]  # value = {'GammaP': <CacheAxis.GammaP: 0>, 'ElectronEnergy': <CacheAxis.ElectronEnergy: 1>}
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
class DNdEeTable:
    """
    Cached/interpolated stand-in for dN_dEe. Build once via over_gamma_p or over_E_e, then call repeatedly -- cheap log-log interpolation instead of the underlying nested integral. Immutable after construction.
    """
    @staticmethod
    def over_E_e(gamma_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, E_e_min: float = 0.0, E_e_max: float = ..., n_points: int = 0, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> DNdEeTable:
        """
        Build over E_e at fixed gamma_p (feeds a dense E_e sweep/plot).
        """
    @staticmethod
    def over_gamma_p(E_e: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, gamma_p_max: float, n_points: int = 0, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> DNdEeTable:
        """
        Build over gamma_p at fixed E_e (feeds q_pair_spectrum_cached's outer integral).
        """
    @typing.overload
    def __call__(self, x: float) -> float:
        ...
    @typing.overload
    def __call__(self, x: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]:
        ...
    def __init__(self) -> None:
        """
        Empty table: reachable()==False, __call__ always returns 0.0.
        """
    def __len__(self) -> int:
        ...
    def axis(self) -> CacheAxis:
        ...
    def domain_hi(self) -> float:
        ...
    def domain_lo(self) -> float:
        ...
    def fixed_value(self) -> float:
        ...
    def reachable(self) -> bool:
        ...
class DNdEeTable2D:
    """
    2D cache over (E_e, gamma_p): one over_gamma_p line per log-spaced E_e node, log-interpolated between adjacent lines. One build serves a whole SED sweep via q_pair_spectrum_cached.
    """
    @staticmethod
    def build(f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, E_p_max: float, E_e_min: float, E_e_max: float, n_E_e_lines: int = 0, n_points_per_line: int = 0, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> DNdEeTable2D:
        """
        Build over the E_e band [E_e_min, E_e_max] with proton budget E_p_max.
        """
    def E_e_max(self) -> float:
        ...
    def E_e_min(self) -> float:
        ...
    def __call__(self, E_e: float, gamma_p: float) -> float:
        ...
    def __init__(self) -> None:
        """
        Empty table: built()==False, __call__ always returns 0.0.
        """
    def __len__(self) -> int:
        ...
    def built(self) -> bool:
        ...
    def epsilon_max(self) -> float:
        ...
    def gamma_p_max(self) -> float:
        ...
class EeWindow:
    """
    Kinematically reachable E_e window at fixed gamma_p (see E_e_window).
    """
    @property
    def hi(self) -> float:
        ...
    @property
    def lo(self) -> float:
        ...
    @property
    def reachable(self) -> bool:
        ...
def E_e_window(gamma_p: float, epsilon_max: float, m_e: float = 0.00051099895) -> EeWindow:
    """
    Closed-form kinematically reachable E_e window at fixed gamma_p.
    """
def E_p_min(E_e: float, epsilon_max: float) -> float:
    """
    Minimum proton energy for which E_e is kinematically reachable (returns +inf if never reachable at any E_p).
    """
def dN_dEe(E_e: float, gamma_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> float:
    """
    Single-proton lab-frame e+/e- spectrum dN/dE_e [GeV^-1]. Expensive: a 3-level nested adaptive integral, up to several seconds per call at UHECR gamma_p against a CMB-like field. Prefer DNdEeTable for repeated evaluations at the same gamma_p or E_e.
    """
def dN_dEe_array(E_e: numpy.ndarray[numpy.float64], gamma_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> numpy.ndarray[numpy.float64]:
    """
    Vectorized dN_dEe over an array of E_e -- still one expensive nested integral PER element; prefer DNdEeTable for a dense sweep.
    """
def dN_dEe_fast(E_e: float, gamma_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> float:
    """
    Field-agnostic fast path: same integration-order swap as dN_dEe_planck with the eps cumulative computed numerically. Exact rearrangement of eq.62, one quadrature level cheaper.
    """
def dN_dEe_general(E_e: float, gamma_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> float:
    """
    Straightforward eq.62 triple integral -- the reference path. dN_dEe itself auto-dispatches to the faster eq.67-style forms (dN_dEe_planck for blackbody fields, dN_dEe_fast otherwise).
    """
def dN_dEe_planck(E_e: float, gamma_p: float, kT: float, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.0001, max_depth: int = 20, panels: int = 24) -> float:
    """
    Planckian-specialized dN_dEe (KA2008 eq.67 extended to finite epsilon_max): exact rearrangement of the eq.62 integral, one quadrature level cheaper. kT in GeV.
    """
def energy_loss_rate(E_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Single-proton BH energy-loss rate -dE_p/dt [GeV/s].
    """
def energy_loss_rate_array(E_p: numpy.ndarray[numpy.float64], f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> numpy.ndarray[numpy.float64]:
    """
    Vectorized energy_loss_rate over an array of E_p (GIL released).
    """
def interaction_length(E_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Mean free path c/rate [cm]; +inf below threshold. Divide by constants.cm_per_Mpc for Mpc.
    """
def interaction_rate(E_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Single-proton BH pair-production interaction rate [s^-1].
    """
def interaction_rate_array(E_p: numpy.ndarray[numpy.float64], f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> numpy.ndarray[numpy.float64]:
    """
    Vectorized interaction_rate over an array of E_p (GIL released).
    """
def loss_timescale(E_p: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    E-folding energy-loss time E_p/(dE_p/dt) [s]; +inf below threshold. Divide by constants.seconds_per_year for years.
    """
def make_table_for_q_pair_spectrum(E_e: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, E_p_max: float, n_points: int = 0) -> DNdEeTable:
    """
    Convenience: builds a DNdEeTable via over_gamma_p with gamma_p_max derived from E_p_max, matching what q_pair_spectrum_cached expects.
    """
def nucleus_energy_loss_rate(E_N: float, Z: float, A: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Z^2-scaled -dE_N/dt [GeV/s] for a nucleus (see nucleus_interaction_rate).
    """
def nucleus_interaction_rate(E_N: float, Z: float, A: float, f_ph: kaspectra._kaspectra.io.PhotonField, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Z^2-scaled interaction rate for a nucleus of energy E_N, charge Z, mass number A (evaluated at the proton energy with the same Lorentz factor).
    """
def q_pair_energy_loss(J_p: kaspectra._kaspectra.io.ProtonSpectrum, f_ph: kaspectra._kaspectra.io.PhotonField, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Population-level BH energy-loss rate, integrated over J_p.
    """
def q_pair_rate(J_p: kaspectra._kaspectra.io.ProtonSpectrum, f_ph: kaspectra._kaspectra.io.PhotonField, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-10, rel_tol: float = 1e-08, max_depth: int = 20, panels: int = 32) -> float:
    """
    Population-level BH interaction rate, integrated over J_p.
    """
def q_pair_spectrum(E_e: float, J_p: kaspectra._kaspectra.io.ProtonSpectrum, f_ph: kaspectra._kaspectra.io.PhotonField, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.001, max_depth: int = 10, panels: int = 6) -> float:
    """
    Population-level BH pair spectrum -- slow but exact (calls dN_dEe at every outer quadrature node). See q_pair_spectrum_cached for repeated calls.
    """
@typing.overload
def q_pair_spectrum_cached(E_e: float, J_p: kaspectra._kaspectra.io.ProtonSpectrum, table: DNdEeTable, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.001, max_depth: int = 10, panels: int = 6) -> float:
    """
    Like q_pair_spectrum, but consults a pre-built DNdEeTable instead of calling dN_dEe at every outer quadrature node. Raises ValueError (pybind11's default translation of std::invalid_argument) if the table doesn't match this E_e/axis, or doesn't cover E_p_max.
    """
@typing.overload
def q_pair_spectrum_cached(E_e: float, J_p: kaspectra._kaspectra.io.ProtonSpectrum, table: DNdEeTable2D, E_p_max: float, epsilon_max: float, abs_tol: float = 1e-25, rel_tol: float = 0.001, max_depth: int = 10, panels: int = 6) -> float:
    """
    SED-sweep variant: one shared DNdEeTable2D serves every E_e in its built band. Raises ValueError on epsilon_max/band/E_p_max mismatch.
    """
def suggested_n_points(x_lo: float, x_hi: float, points_per_decade: int = 12) -> int:
    """
    Grid-density heuristic for DNdEeTable build points (default 12/decade, benchmarked: ~2.5% max interpolation error near the peak).
    """
