"""Shared linked-application mechanics for CV02-CV12, IV01-IV06, and XM01-XM03.

This module owns evidence plumbing, not expected physics.  Every case passes a
flat, reviewed SI argument vector to the selected ``amps`` executable and owns
its independent reference formula in its own directory.  Keeping process and
provenance mechanics here makes later validation cases follow exactly the same
application-level contract without duplicating fragile subprocess checks.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence


def load_input(path: Path, case_id: str) -> Dict[str, Any]:
    """Load one versioned case input and reject cross-case substitution."""
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if (not isinstance(value, dict) or
            value.get("schema") != "srcsep-validation-case-input-v1" or
            value.get("case_id") != case_id):
        raise ValueError(f"input must use schema v1 and case_id={case_id}")
    return value


def atomic_json(path: Path, value: Dict[str, Any]) -> None:
    """Commit validation metadata atomically so partial runs are unmistakable."""
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Reload saved evidence rather than scoring unsaved in-memory values."""
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"linked model CSV is empty: {path}")
    if any(None in row for row in rows):
        raise RuntimeError(f"linked model CSV has more fields than its header: {path}")
    return rows


def write_csv(path: Path, fieldnames: Sequence[str],
              rows: Iterable[Dict[str, Any]]) -> None:
    """Write a small derived/reference CSV through a same-directory temporary."""
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _write_manifest(path: Path, case_id: str,
                    arguments: Sequence[str]) -> None:
    """Serialize strict name/value tokens for the native registry callback."""
    if not arguments or len(arguments) % 2:
        raise ValueError("native arguments must be non-empty name/value pairs")
    if any(not token or "\n" in token or "\r" in token for token in arguments):
        raise ValueError("native arguments must be non-empty single-line tokens")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"srcsep-controlled-native-args-v1:{case_id}\n")
        for token in arguments:
            stream.write(token + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def run_linked_model(*, case_id: str, arguments: Sequence[str],
                     source_root: Path, output_dir: Path, executable: Path,
                     timeout: Optional[float]) -> Dict[str, Path]:
    """Execute one native case and require its registry result to be PASS.

    The caller has already been checked by ``validation/run_case.py`` with
    ``--list-tests``.  This second check validates the selected case's actual
    report and model artifact; a process that merely returns zero cannot pass.
    """
    native = output_dir / "native"
    native.mkdir(parents=True, exist_ok=True)
    manifest = native / f"{case_id}_native_arguments.txt"
    model = native / f"{case_id}_model.csv"
    report = native / f"{case_id}_native.json"
    junit = native / f"{case_id}_native.xml"
    log = output_dir / f"{case_id}_run.log"
    _write_manifest(manifest, case_id, arguments)
    command = [
        str(executable), "--test", case_id,
        "--test-input", str(manifest),
        "--test-output-dir", str(native),
        "--test-json", str(report),
        "--test-junit", str(junit),
    ]
    # The higher-level runner prints its Python orchestration command, but the
    # linked application is the numerical system under test. Print the exact,
    # shell-escaped AMPS invocation as a separate line so a user can audit or
    # reproduce the actual model call without opening the retained log first.
    print("RUN:", shlex.join(command), flush=True)
    try:
        completed = subprocess.run(
            command, cwd=str(source_root), text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode("utf-8", errors="replace")
        log.write_text("COMMAND: " + " ".join(command) + "\n" + partial +
                       "\nTIMEOUT\n", encoding="utf-8")
        raise RuntimeError(f"linked {case_id} exceeded {timeout} seconds") from error
    log.write_text("COMMAND: " + " ".join(command) + "\n" + completed.stdout +
                   f"\nEXIT: {completed.returncode}\n", encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"linked {case_id} exited {completed.returncode}; see {log}")
    with report.open("r", encoding="utf-8") as stream:
        native_report = json.load(stream)
    results = native_report.get("results", []) if isinstance(native_report, dict) else []
    if (len(results) != 1 or str(results[0].get("id", "")).upper() != case_id or
            str(results[0].get("status", "ERROR")).upper() != "PASS"):
        raise RuntimeError(f"linked registry did not report PASS for {case_id}")
    if not model.is_file() or model.stat().st_size == 0:
        raise RuntimeError(f"linked application did not publish {model.name}")
    return {"manifest": manifest, "model": model, "report": report,
            "junit": junit, "log": log}


def metric(name: str, value: float, tolerance: float, comparison: str,
           units: str) -> Dict[str, Any]:
    """Construct the metric schema consumed by both JSON and plot runners."""
    return {"name": name, "value": float(value), "tolerance": float(tolerance),
            "comparison": comparison, "units": units}


def metrics_pass(metrics: Sequence[Dict[str, Any]]) -> bool:
    """Evaluate the deliberately small comparison vocabulary used by cases."""
    for item in metrics:
        value, tolerance = float(item["value"]), float(item["tolerance"])
        comparison = item["comparison"]
        if comparison == "<=" and not value <= tolerance:
            return False
        if comparison == ">=" and not value >= tolerance:
            return False
    return True


def finish_result(*, case_id: str, started: float, seed: int,
                  input_path: Path, executable: Path,
                  metrics: Sequence[Dict[str, Any]], artifacts: Sequence[Path],
                  message: str) -> Dict[str, Any]:
    """Return one complete registry-compatible scientific case result."""
    passed = metrics_pass(metrics)
    return {
        "id": case_id,
        "status": "PASS" if passed else "FAIL",
        "message": message if passed else message.replace("passed", "failed"),
        "elapsed_seconds": time.monotonic() - started,
        "seed": seed,
        "configuration": [
            f"input={input_path}", "execution=linked-srcsep-amps",
            f"executable={executable}", f"executable_sha256={sha256(executable)}",
            "reference=independent-python",
        ],
        "metrics": list(metrics),
        "artifacts": [str(path) for path in artifacts],
    }
