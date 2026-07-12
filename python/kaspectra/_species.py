from ._kaspectra import Species

_ALIASES = {
    "gamma": Species.Gamma, "g": Species.Gamma,
    "e+": Species.Positron, "positron": Species.Positron,
    "e-": Species.Electron, "electron": Species.Electron,
    "numu": Species.NuMu, "numubar": Species.NuMuBar,
    "nue": Species.NuE, "nuebar": Species.NuEBar,
}


def coerce(value):
    """Species enum member or alias string -> Species enum member."""
    if isinstance(value, Species):
        return value
    try:
        return _ALIASES[value]
    except (KeyError, TypeError) as exc:
        raise ValueError(
            f"unknown species {value!r}; use a Species member or one of {sorted(_ALIASES)}"
        ) from exc
