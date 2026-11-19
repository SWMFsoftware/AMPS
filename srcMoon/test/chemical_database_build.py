"""Authoritative shared chemistry dependencies for standalone Moon adapters.

The historical Moon CPP files remain compatibility translation units; the
physics implementation lives in src/models/chemical_database. This helper
does not discover/fallback to alternate implementations: unknown components
or absent declared source bytes fail explicitly. Build commands and manifests
must use the same returned source/header list, including this helper itself.
"""
from pathlib import Path

_ORDER = ("Runtime", "DataIO", "Photochemistry", "ElectronImpact", "ChargeExchange", "RateServices")
_COMPONENTS = {
    "photochemistry": ("DataIO", "Photochemistry"),
    "electron-impact": ("DataIO", "ElectronImpact"),
    "charge-exchange": ("DataIO", "ChargeExchange"),
    # Runtime parsers may reference every installer even when a particular
    # fixture disables it. Linking stubs would bypass real initialization.
    "runtime": _ORDER,
    "rate-services": _ORDER,
}


def _selected(components):
    if isinstance(components, str):
        components = (components,)
    selected = set()
    for component in components:
        if component not in _COMPONENTS:
            raise ValueError("unknown shared chemistry component: "+str(component))
        selected.update(_COMPONENTS[component])
    if not selected:
        raise ValueError("at least one explicit shared chemistry component is required")
    return [name for name in _ORDER if name in selected]


def shared_sources(root, components):
    """Return stable, deduplicated maintained CPP paths; no implicit fallback."""
    base = Path(root)/"src/models/chemical_database"
    return [base/(name+".cpp") for name in _selected(components)]


def shared_dependencies(root, components):
    """Return actual compiled sources/public headers plus this recipe identity."""
    base = Path(root)/"src/models/chemical_database"
    paths = []
    for name in _selected(components):
        paths.extend((base/(name+".cpp"), base/(name+".h")))
    if "Runtime" in _selected(components):
        # The shared runtime uses the maintained generic input grammar. Its
        # header is part of the actual compiler dependency, not a test parser.
        paths.append(Path(root)/"src/general/input/amps_input.h")
    paths.append(Path(root)/"srcMoon/test/chemical_database_build.py")
    return paths
