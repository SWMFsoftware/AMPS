#!/usr/bin/env python3
"""Shared fail-closed runtime helpers for executable Earth validation tests.

The scientific runners under ``srcEarth/test`` accept an ``--amps`` path so a
reviewed executable can be tested without rebuilding it.  A filesystem check is
not sufficient: one AMPS checkout can be rebuilt for another application while
leaving the same executable pathname in place.  C9/C10 experienced exactly that
failure mode when a historical absolute path began naming the SEP application.

This module deliberately checks *application identity*, not numerical results.
It does not alter a solver command, input deck, reference product, tolerance, or
pass/fail threshold.  Its other responsibility is atomic status/provenance I/O,
so an interrupted or rejected launch cannot leave an older result looking like
the result of the new invocation.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence


EARTH_HELP_BANNER = "AMPS Earth: Gridless Energetic Particle Solver"
EARTH_HELP_REQUIRED_TOKENS = (
    EARTH_HELP_BANNER,
    "-mode <3d|gridless>",
    "No solver is run.",
)


class EarthExecutableError(RuntimeError):
    """Raised when an ``--amps`` target is absent or is not the Earth driver."""


def utc_now() -> str:
    """Return a stable UTC timestamp for JSON provenance."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def sha256_file(path: Path) -> str:
    """Hash a file without loading a multi-megabyte executable into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probe_earth_executable(
    executable: Path,
    *,
    cwd: Optional[Path] = None,
    timeout_seconds: float = 30.0,
) -> Dict[str, Any]:
    """Resolve and positively identify an AMPS Earth executable.

    The probe invokes only ``-h``.  The Earth driver documents that this exits
    before model initialization, so the check does not initialize MPI, Geopack,
    SPICE, a mesh, or a physical calculation.  Requiring several Earth-specific
    help tokens prevents an executable for SEP or another AMPS application from
    reaching a costly multi-rank launch and failing after it has created a
    misleading partial output tree.
    """

    path = executable.expanduser()
    if not path.is_absolute():
        path = ((cwd or Path.cwd()) / path).resolve()
    else:
        path = path.resolve()

    if not path.is_file():
        raise EarthExecutableError("AMPS executable not found: %s" % path)
    if not os.access(str(path), os.X_OK):
        raise EarthExecutableError("AMPS executable is not executable: %s" % path)

    try:
        completed = subprocess.run(
            [str(path), "-h"],
            cwd=str((cwd or Path.cwd()).resolve()),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise EarthExecutableError(
            "AMPS executable identity probe timed out after %.1f s: %s"
            % (timeout_seconds, path)
        ) from exc
    except OSError as exc:
        raise EarthExecutableError(
            "AMPS executable identity probe could not start %s: %s" % (path, exc)
        ) from exc

    output = completed.stdout or ""
    missing = [token for token in EARTH_HELP_REQUIRED_TOKENS if token not in output]
    if completed.returncode != 0 or missing:
        details = []
        if completed.returncode != 0:
            details.append("help exit=%d" % completed.returncode)
        if missing:
            details.append("missing Earth help token(s): %s" % ", ".join(missing))
        raise EarthExecutableError(
            "AMPS executable is not the required Earth driver: %s (%s)"
            % (path, "; ".join(details))
        )

    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "identity_probe": "-h",
        "identity_banner": EARTH_HELP_BANNER,
        "identity_probe_return_code": completed.returncode,
    }


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Atomically replace a JSON status/result file in its destination directory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.tmp.%d" % (path.name, os.getpid()))
    try:
        temporary.write_text(json.dumps(dict(payload), indent=2) + "\n", encoding="utf-8")
        os.replace(str(temporary), str(path))
    finally:
        # ``os.replace`` removes the temporary path.  The cleanup matters only
        # when serialization or replacement raises midway through the write.
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def archive_previous_result(path: Path, *, run_started_utc: str) -> Optional[Path]:
    """Move a prior top-level result aside before a new physical calculation.

    The archive is recoverable and never overwrites an earlier archive.  Branch
    directories are managed by C9/C10 themselves, but their top-level result can
    otherwise survive a failed first launch and be mistaken for current output.
    """

    if not path.exists():
        return None
    token = run_started_utc.replace("-", "").replace(":", "")
    candidate = path.with_name("%s.previous-%s%s" % (path.stem, token, path.suffix))
    suffix = 1
    while candidate.exists():
        candidate = path.with_name(
            "%s.previous-%s-%d%s" % (path.stem, token, suffix, path.suffix)
        )
        suffix += 1
    os.replace(str(path), str(candidate))
    return candidate


def initial_run_status(
    *,
    schema: str,
    test_id: str,
    output_root: Path,
    argv: Sequence[str],
    execution_mode: str,
) -> Dict[str, Any]:
    """Create the fail-closed status object written before solver execution."""

    started = utc_now()
    return {
        "schema": schema,
        "test_id": test_id,
        "status": "STARTING",
        "started_utc": started,
        "finished_utc": None,
        "execution_mode": execution_mode,
        "output_root": str(output_root.resolve()),
        "argv": list(argv),
        "executable": None,
        "archived_previous_result": None,
        "result_current": False,
        "return_code": None,
        "message": "validation invocation has not completed",
    }


def finish_run_status(
    status: Dict[str, Any],
    *,
    state: str,
    return_code: int,
    message: str,
    result_current: bool,
) -> Dict[str, Any]:
    """Finalize a status object without interpreting any scientific threshold."""

    status.update(
        {
            "status": state,
            "finished_utc": utc_now(),
            "result_current": bool(result_current),
            "return_code": int(return_code),
            "message": message,
        }
    )
    return status
