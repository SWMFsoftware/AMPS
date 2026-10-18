"""Shared paths for generated COPS products; no files are created on import.

Default output root: <AMPS>/cops_output, derived from this application's location.
COPS_OUTPUT_ROOT may override it (absolute or relative to the invoking cwd).
Relative output arguments are interpreted under that root; absolute output
arguments outside the application are honored. Output paths resolving inside
LiouvilleCOPS are rejected. Product readers prefer new outputs, then existing
cwd/application files for compatibility. Explicit ./, ../, or absolute input
paths select exactly that input instead of invoking the product search.
External raw/L4/image inputs default to
/data/vtenishe/rosetta_data_end_of_mission_cops and can be relocated with
COPS_DATA_ROOT, COPS_L2_ROOT, and COPS_SOURCES_ROOT.
"""

import os
from pathlib import Path
import sys


APP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = Path("/data/vtenishe/rosetta_data_end_of_mission_cops")


def data_root():
    """External read-only COPS data root (override with COPS_DATA_ROOT)."""
    return Path(os.environ.get("COPS_DATA_ROOT", DEFAULT_DATA_ROOT)).expanduser().resolve()


def sources_root():
    """External L4/image sources (override with COPS_SOURCES_ROOT)."""
    configured = os.environ.get("COPS_SOURCES_ROOT")
    return Path(configured).expanduser().resolve() if configured else data_root() / "sources"


def cops_l2_root():
    """External raw L2 tree (override with COPS_L2_ROOT)."""
    configured = os.environ.get("COPS_L2_ROOT")
    return Path(configured).expanduser().resolve() if configured else data_root() / "data_COPS_L2"


def output_root():
    configured = os.environ.get("COPS_OUTPUT_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    return APP_ROOT.parent.parent / "cops_output"


def output_path(value):
    path = Path(value).expanduser()
    candidate = path if path.is_absolute() else output_root() / path
    resolved = candidate.resolve()
    application = APP_ROOT.resolve()
    if resolved == application or application in resolved.parents:
        raise ValueError("generated output must be outside LiouvilleCOPS: " + str(resolved))
    return candidate


def prepare_output(value):
    """Resolve a generated filename/prefix and create its parent directory."""
    path = output_path(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)


def product_input(value, app_root=None):
    """Read generated products in the new root, with read-only legacy fallback."""
    name = os.fspath(value)
    path = Path(name).expanduser()
    if path.is_absolute() or name.startswith(("./", "../")):
        return str(path.resolve())
    primary = output_root() / path
    if primary.exists():
        return str(primary)
    for candidate in (Path.cwd() / path, (app_root or APP_ROOT) / path):
        if candidate.exists():
            print("Using existing input outside cops_output: " + str(candidate.resolve()), file=sys.stderr)
            return str(candidate.resolve())
    return str(primary)
