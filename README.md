# KASpectra

C++17 header-only library (with Python bindings) for cosmic-ray secondary
spectra: gamma-rays, electrons/positrons, and neutrinos produced by
relativistic protons on gas and photon targets.

| Channel | Physics source | Species |
|---|---|---|
| `kaspectra::pp` | Kelner, Aharonian & Bugayov 2006 (astro-ph/0606058) | γ, e⁺, e⁻, ν_μ, ν̄_μ, ν_e, ν̄_e |
| `kaspectra::pgamma` | Kelner & Aharonian 2008 (arXiv:0803.0688), photomeson Φ tables | γ, e⁺, e⁻, ν_μ, ν̄_μ, ν_e, ν̄_e |
| `kaspectra::bh` | KA2008 Sec. IV + Blumenthal 1970 (eq.10) + Chodorowski et al. 1992 | e⁺/e⁻ spectrum, rates, energy loss |

Every formula was transcribed from the actual papers (not from memory), and
every module carries golden-value tests recomputed independently of the
implementation. `kaspectra::q_total_species` sums all channels per species.

## Units

GeV for energies, cm for lengths, s for time, throughout. Photon fields are
differential number densities `f_ph(ε) [GeV⁻¹ cm⁻³]`; proton spectra
`J_p(E_p)` carry whatever normalization you give them. Convenience constants
`constants::cm_per_Mpc` and `constants::seconds_per_year` convert observables.

## C++ quickstart

Header-only: add `include/` to your include path, or install and use
`find_package`:

```bash
cmake -B build && cmake --install build --prefix ~/.local
```

```cmake
find_package(KASpectra REQUIRED)
target_link_libraries(app PRIVATE KASpectra::kaspectra)
```

```cpp
#include <kaspectra/total.hpp>

kaspectra::io::PowerLawSpectrum J_p(1.0, 2.0);
kaspectra::io::BlackbodyPhotonField cmb(2.725);          // T in Kelvin

// gamma-ray emissivity at 1 TeV from pp (n_H = 1 cm^-3) + pgamma:
double q = kaspectra::q_total_species(
    kaspectra::Species::Gamma, 1e3, J_p, 1.0, cmb,
    /*E_p_max=*/1e6, /*E_pi_max=*/1e6, /*epsilon_max=*/1e-9);

// Bethe-Heitler energy-loss time of a 1e10 GeV proton on the CMB, in years:
double t_loss = kaspectra::bh::loss_timescale(1e10, cmb, 1e-6)
              / kaspectra::constants::seconds_per_year;
```

Tests and examples build by default (`KASPECTRA_BUILD_TESTS`,
`KASPECTRA_BUILD_EXAMPLES`); run `./build/test/kaspectra_tests` for the full
suite. Example CLIs in `examples/` write CSV to stdout.

## Python quickstart

```bash
pip install .        # builds the extension via scikit-build-core
```

```python
import numpy as np
import kaspectra as ks

J_p = ks.io.PowerLawSpectrum(1.0, 2.0)
cmb = ks.io.BlackbodyPhotonField(2.725)

# Species accept the enum or strings: "gamma", "e+", "e-", "numu", ...
q = ks.pp.q_species("gamma", np.logspace(2, 5, 50), J_p, 1.0, 1e6, 1e6)

# Bethe-Heitler pair SED, amortized through one 2D cache build:
table = ks.bh.DNdEeTable2D.build(cmb, 1e-6, 1e12, 1e4, 1e8)
sed = [ks.bh.q_pair_spectrum_cached(E, J_p, table, 1e12, 1e-6)
       for E in np.logspace(4, 8, 100)]
```

Custom inputs can be written directly in Python by subclassing
`ks.io.ProtonSpectrum` / `ks.io.PhotonField` and implementing
`__call__` (orders of magnitude slower per evaluation than the built-in C++
implementations, but fully supported, including from the parallel table
builds). `ks.io.CompositePhotonField` sums fields (CMB + IR + starlight).

`pytest` runs the fast suite; expensive parity tests are opt-in with
`pytest -m slow`.

## Performance notes

`bh::dN_dEe` is the expensive primitive (a nested adaptive integral over
Blumenthal's double-differential cross-section). The library keeps it fast
via, in order of preference:

1. **Planckian fast path** (`dN_dEe_planck`, auto-dispatched for blackbody
   fields): KA2008 eq.67's analytic collapse of the photon-energy integral,
   extended exactly to finite `epsilon_max`. ~250× the reference path.
2. **Field-agnostic fast path** (`dN_dEe_fast`, auto-dispatched otherwise):
   the same integration-order swap with the photon-energy cumulative
   computed numerically.
3. **Caches**: `DNdEeTable` (1D, parallel build) and `DNdEeTable2D` (one
   build per SED sweep) with log-log interpolation, ~2.5% max error at the
   default grid density (benchmarked).

`dN_dEe_general` remains available as the straightforward reference path.
If an integral may be under-converged, `kaspectra::math::max_depth_hits()`
reports panels that exhausted recursion without meeting tolerance.

## Accuracy

- Parametrizations carry the papers' own quoted accuracy (typically a few
  percent within their fitted domains).
- `bh` fit functions are evaluated in cancellation-free forms
  (machine-precision at the pair-production threshold, where the printed
  eq.2.4 loses all significant digits).
- The pγ photomeson tables clamp outside η/η₀ ∈ [1.1, 100] as published;
  the γ-ray amplitude `B` blends to the paper's stated zero at threshold.
