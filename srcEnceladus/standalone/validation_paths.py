"""Portable resolution of the external Cassini validation-data tree.

The source tree contains code and fixed model definitions. Observation tables,
archive products, and generated studies remain in a separate data tree whose
absolute location is machine dependent. Explicit paths always take precedence
over the common root; the root may be supplied on the command line or through
``ENCELADUS_VALIDATION_DATA_ROOT``.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path


VALIDATION_ROOT_ENV = "ENCELADUS_VALIDATION_DATA_ROOT"
PLUME_FIT_RELATIVE = Path("data/processed/plume_fit_v1")
CO2_FIT_RELATIVE = Path("data/processed/co2_fit_v1")
JET_CATALOG_RELATIVE = Path("data/derived/porco2014_jet_catalog.json")


def add_validation_root_argument(parser: argparse.ArgumentParser) -> None:
    """Add the common, portable data-root option to a program parser."""
    parser.add_argument(
        "--validation-data-root", type=Path,
        help=("Cassini validation-data root containing data/ and studies/. "
              f"If omitted, {VALIDATION_ROOT_ENV} is used."))


def resolve_validation_root(value: Path | str | None) -> Path | None:
    """Resolve and validate the CLI/environment validation-data root."""
    selected = value if value is not None else os.environ.get(VALIDATION_ROOT_ENV)
    if selected is None or str(selected).strip() == "":
        return None
    path = Path(selected).expanduser().resolve()
    if not path.is_dir():
        raise ValueError(f"Validation-data root is not a directory: {path}")
    return path


def resolve_input_path(explicit: Path | str | None, root: Path | None,
                       relative: Path, option_name: str) -> Path:
    """Resolve a required input, preferring its explicit path over the root."""
    path = (Path(explicit).expanduser().resolve() if explicit is not None
            else (root / relative).resolve() if root is not None else None)
    if path is None:
        raise ValueError(
            f"{option_name} is required unless --validation-data-root or "
            f"{VALIDATION_ROOT_ENV} is supplied")
    if not path.exists():
        raise ValueError(f"Required input does not exist: {path}")
    return path


def optional_root_path(explicit: Path | str | None, root: Path | None,
                       relative: Path) -> Path | None:
    """Return an explicit input or its conventional location under the root."""
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    return (root / relative).resolve() if root is not None else None
