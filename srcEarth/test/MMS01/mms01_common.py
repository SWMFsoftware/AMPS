#!/usr/bin/env python3
"""Shared contracts for the MMS1/HPCA local-product validation.

This module contains no subprocess or plotting side effects.  Preparation,
comparison, visualization, and tests import the same immutable event definition so
the QUICK and COMPLETE entry points cannot silently select different science data.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


CASE_ID = "MMS01"
CDFLIB_VERSION = "1.3.2"
HERE = Path(__file__).resolve().parent
CDF_READER_REQUIREMENTS = HERE / "requirements-cdf-reader.txt"
# Keep third-party runtime files outside the source tree while giving the normal
# test-list and direct-run commands one deterministic location to discover.  The
# directory is evidence/build state beneath test_output, not a committed dependency.
LOCAL_PYTHON_DEPENDENCIES = HERE.parents[2] / "test_output" / "MMS01" / "python_deps"
EVENT_START = "2015-11-10T13:57:54"
EVENT_STOP = "2015-11-10T14:04:24"
EVENT_DIRECTORY = "mms1_hpca_magnetosheath_20151110"
ION_FILE = "mms1_hpca_srvy_l2_ion_20151110120000_v4.1.0.cdf"
MEC_FILE = "context/mms/mms1_mec_srvy_l2_epht89d_20151110_v2.0.2.cdf"
CMAD_FILE = "HPCA_CMAD_20260430.pdf"

# These identities come from the downloaded bundle's source_manifest.json.  They
# are repeated here intentionally: changing the external manifest cannot authorize
# different bytes for a code-validation reference without a reviewed code change.
SOURCE_FILES = {
    ION_FILE: (179398777, "c80fe2b42ead85b8f839724a31eca043e6a0908a66ddf9ec43e8e5b719d24af9"),
    MEC_FILE: (2597811, "44700d52b7e2f5e09fad11ba6ad41e95d42ad2e29447949b6968ae1dc3a92544"),
    CMAD_FILE: (1736585, "7f2dd7d1cacc2d35266aa07932791d8805732c70eaac62e961f12e857fa5fcb4"),
}

# QUICK is a deterministic first/middle/last sample selected by record index before
# any reconstruction/reference residual is computed. COMPLETE evaluates every qualified
# half-spin in the frozen interval.  Neither profile changes formulas or gates.
PROFILES = {
    "QUICK": (0, 19, 38),
    "COMPLETE": tuple(range(39)),
}

EXPECTED_SWEEPS = 624
EXPECTED_DISTRIBUTIONS = 39
SWEEPS_PER_DISTRIBUTION = 16
ANODES = 16
ENERGIES = 63
REQUIRED_QUALITY = 2
REQUIRED_SCIENCE_MODE = 2

# CMAD revision 4, equation 5-13 and Table 5-15.  The document defines H+ as mass
# number one for this conversion; using the precise physical proton mass here would
# be a different, undocumented convention and measurably disagree with Level 2.
HPCA_HPLUS_MASS_EV_S2_CM2 = 1.0453453e-12

# CMAD equation 5-15, Tables 5-6 and 5-8.  Δalpha and Δbeta are the documented
# average angular-response quantities used for HPCA moments.  These values are not
# inferred from centroid spacing.
DELTA_ALPHA_RAD = 0.318
DELTA_BETA_RAD = 0.39
ANODE_CORRECTION = (
    2.09, 1.78, 1.78, 2.09, 3.13, 8.93, 8.93, 3.13,
    2.09, 1.78, 1.78, 2.09, 3.13, 8.93, 8.93, 3.13,
)

# The CMAD mass coefficient has eight significant digits.  Squaring it in equation
# 5-13 makes its last-place rounding visible at about 1e-7 in the archived product.
# A 1e-6 gate covers only documented-constant representation, not scientific model
# disagreement.  Zero-valued cells are compared exactly and never log-shifted.
MAX_RELATIVE_ERROR = 1.0e-6


def ensure_cdflib() -> Any:
    """Load the pinned CDF reader, installing it into test_output when absent.

    MMS01 is registered as a directly executable validation command. Requiring an
    activated, undocumented virtual environment would make that command fail before
    any science gate runs. We therefore use a repository-local dependency directory
    and bootstrap only the one missing pure-Python reader, pinned to an exact version.
    NumPy and Matplotlib remain ordinary declared application dependencies.

    Set ``MMS01_AUTO_INSTALL_DEPS=0`` to prohibit network/bootstrap activity. In
    that mode a missing reader fails with the exact manual installation command.
    """
    dependency_text = str(LOCAL_PYTHON_DEPENDENCIES)
    if dependency_text not in sys.path:
        # Prefer the pinned local copy over a differently versioned user install.
        sys.path.insert(0, dependency_text)
    try:
        import cdflib  # type: ignore
    except ImportError as original_error:
        auto_install = os.environ.get("MMS01_AUTO_INSTALL_DEPS", "1").strip().lower()
        if auto_install in ("0", "false", "no", "off"):
            raise RuntimeError(
                "MMS01 requires cdflib==1.3.2. Automatic installation is disabled; run: "
                f"{sys.executable} -m pip install --no-deps --target "
                f"{LOCAL_PYTHON_DEPENDENCIES} --require-hashes --requirement "
                f"{CDF_READER_REQUIREMENTS}"
            ) from original_error
        LOCAL_PYTHON_DEPENDENCIES.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
            "--no-input", "--no-deps", "--target", str(LOCAL_PYTHON_DEPENDENCIES),
            "--require-hashes", "--requirement", str(CDF_READER_REQUIREMENTS),
        ]
        print("[MMS01] installing pinned CDF reader into "
              f"{LOCAL_PYTHON_DEPENDENCIES}", flush=True)
        result = subprocess.run(command, check=False)
        if result.returncode != 0:
            raise RuntimeError(
                "MMS01 could not install cdflib==1.3.2 into its repository-local "
                f"dependency directory: {LOCAL_PYTHON_DEPENDENCIES}"
            ) from original_error
        # The failed import may have populated import-system caches before pip
        # created the package directory.
        import importlib
        importlib.invalidate_caches()
        import cdflib  # type: ignore
    version = str(getattr(cdflib, "__version__", ""))
    if version != CDFLIB_VERSION:
        raise RuntimeError(
            f"MMS01 requires cdflib=={CDFLIB_VERSION}, but Python loaded {version or 'unknown'} "
            f"from {getattr(cdflib, '__file__', 'an unknown location')}"
        )
    return cdflib


def utc_now() -> str:
    """Return a timezone-explicit timestamp for evidence manifests."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    """Hash large CDF and result files without loading them into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    """Publish JSON atomically so interrupted runs cannot leave valid-looking files."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def resolve_event_root(data_root: Path) -> Path:
    """Accept either the bundle root or the event directory, never search broadly."""
    root = data_root.expanduser().resolve()
    candidates = (root, root / EVENT_DIRECTORY)
    matches = [candidate for candidate in candidates if (candidate / ION_FILE).is_file()]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"{data_root}: expected exactly one {ION_FILE} at the supplied root or "
            f"its {EVENT_DIRECTORY} child"
        )
    return matches[0]


def verify_sources(event_root: Path) -> dict[str, dict[str, Any]]:
    """Authenticate every source used by this test before reading scientific values."""
    result: dict[str, dict[str, Any]] = {}
    for relative, (expected_size, expected_hash) in SOURCE_FILES.items():
        path = event_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"missing MMS01 source: {path}")
        actual_size = path.stat().st_size
        actual_hash = sha256_file(path)
        if actual_size != expected_size or actual_hash != expected_hash:
            raise ValueError(f"source identity mismatch: {path}")
        result[relative] = {
            "path": str(path), "bytes": actual_size, "sha256": actual_hash,
        }
    return result


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    """Write a strict, headered CSV artifact."""
    names = tuple(fieldnames)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV and reject missing or duplicate column names."""
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{path}: missing or duplicate CSV header")
        return list(reader)
