#!/usr/bin/env python3
"""Shared execution support for the srcMoon U-series verification tests.

The runner keeps four outcomes distinct.  A missing build, production
capability, or qualified dataset is SKIPPED; an infrastructure problem is
ERROR; an evaluated scientific criterion is PASS or FAIL.  In particular,
the library never upgrades source presence or a successful local probe into a
linked-runtime or observational validation result.
"""

from __future__ import annotations

import csv
import json
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any


TEST_ROOT = Path(__file__).resolve().parents[1]
SRCMOON_ROOT = TEST_ROOT.parent
REPO_ROOT = SRCMOON_ROOT.parent
DEFAULT_BUILD = REPO_ROOT / "build"
DEFAULT_OUTPUT = REPO_ROOT / "test_output" / "srcMoon" / "unit"

STATUS_EXIT = {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}


def test_directories() -> dict[str, Path]:
    """Return the frozen UXX registry discovered from per-test directories.

    Directory names supply only the stable ID.  Scientific names and execution
    policy always come from the committed acceptance contract.
    """
    result: dict[str, Path] = {}
    for path in sorted(TEST_ROOT.glob("U[0-9][0-9]_*")):
        if path.is_dir():
            result[path.name[:3]] = path
    return result


def load_acceptance(test_id: str) -> dict[str, Any]:
    """Load one machine-readable test contract without supplying defaults."""
    test_dir = test_directories()[test_id]
    with (test_dir / "reference" / "acceptance.json").open(
        encoding="utf-8"
    ) as stream:
        return json.load(stream)


def git_revision() -> str:
    """Record the source identity used by a result, or ``unknown`` on error."""
    command = ["git", "rev-parse", "HEAD"]
    completed = subprocess.run(
        command, cwd=REPO_ROOT, text=True, capture_output=True, check=False
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def _probe_compile_command(build_dir: Path, executable: Path) -> list[str]:
    """Construct the thin-probe link command against generated production code.

    Include directories point at Config.pl's generated headers, so the probe
    sees exactly the species and feature macros used by libAMPS.a.  Expected
    values remain in the probe as independent analytical/control-point oracles.
    """
    include_dirs = [
        build_dir / "pic",
        build_dir / "main",
        build_dir / "meshAMR",
        build_dir / "interface",
        build_dir / "general",
        build_dir / "models" / "exosphere",
        build_dir / "models" / "electron_impact",
        build_dir / "models" / "sputtering",
        build_dir / "models" / "dust",
        build_dir / "models" / "charge_exchange",
        build_dir / "models" / "photolytic_reactions",
        build_dir / "species",
        REPO_ROOT / "srcInterface",
        REPO_ROOT / "share" / "Library" / "src",
        REPO_ROOT,
    ]
    # Orbit-enabled Moon builds expose CSPICE types through generated headers.
    # The toolkit location is a build dependency, not validation evidence; the
    # authoritative kernel files are hashed separately by run_rotating_frame().
    spice_toolkit = Path(
        os.environ.get("MOON_SPICE_TOOLKIT_ROOT", "/home/vtenishe/SPICE/cspice")
    )
    if spice_toolkit.joinpath("include", "SpiceUsr.h").is_file():
        include_dirs.append(spice_toolkit / "include")
    compiler = os.environ.get("MPICXX", "mpicxx")
    command = [compiler, "-std=c++17", "-O2"]
    for include_dir in include_dirs:
        command.extend(["-I", str(include_dir)])
    command.extend(
        [
            str(TEST_ROOT / "common" / "production_kernel_probe.cpp"),
            str(build_dir / "libAMPS.a"),
            "-lstdc++",
            "-lmpi_cxx",
            "-o",
            str(executable),
        ]
    )
    return command


def compile_probe(build_dir: Path, output_root: Path) -> tuple[Path | None, dict[str, Any]]:
    """Build or reuse the production-kernel probe and preserve diagnostics.

    Missing generated artifacts are a capability precondition and therefore
    SKIPPED.  A compiler/linker failure with artifacts present is test
    infrastructure ERROR.  The mtime check is only a build cache; every result
    still records whether compilation occurred.
    """
    required = [
        build_dir / "libAMPS.a",
        build_dir / "main" / "Moon.h",
        build_dir / "pic" / "pic.h",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        return None, {
            "status": "SKIPPED",
            "reason": "generated AMPS build artifacts are unavailable",
            "missing": missing,
        }

    binary_dir = output_root / "_bin"
    binary_dir.mkdir(parents=True, exist_ok=True)
    executable = binary_dir / "production_kernel_probe"
    source = TEST_ROOT / "common" / "production_kernel_probe.cpp"
    library = build_dir / "libAMPS.a"

    if executable.is_file() and executable.stat().st_mtime >= max(
        source.stat().st_mtime, library.stat().st_mtime
    ):
        return executable, {"status": "PASS", "compiled": False}

    command = _probe_compile_command(build_dir, executable)
    completed = subprocess.run(
        command, cwd=REPO_ROOT, text=True, capture_output=True, check=False
    )
    detail = {
        "status": "PASS" if completed.returncode == 0 else "ERROR",
        "compiled": True,
        "command": command,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "returncode": completed.returncode,
    }
    return (executable if completed.returncode == 0 else None), detail


def run_probe_case(
    acceptance: dict[str, Any], build_dir: Path, output_root: Path
) -> dict[str, Any]:
    """Execute every probe mode and any declared configuration invariants.

    A linked kernel value alone cannot establish that the production input
    selected that kernel exactly once.  Contracts may therefore require tokens
    in the maintained input tree and in Config.pl's generated build tree.  The
    latter check is deliberately based on ``build_dir`` so an explicitly
    selected build is audited rather than the repository's default build.
    """
    executable, compile_result = compile_probe(build_dir, output_root)
    if executable is None:
        return compile_result

    outputs: list[dict[str, Any]] = []
    for mode in acceptance["probe_modes"]:
        environment = os.environ.copy()
        if mode == "lola-geometry":
            # U08 writes generated surfaces only beneath the selected test
            # output root; validation-data raw directories remain immutable.
            u08_output = output_root / "U08" / "surface_probe"
            u08_output.mkdir(parents=True, exist_ok=True)
            environment["MOON_U08_OUTPUT"] = str(u08_output)
        completed = subprocess.run(
            [str(executable), mode],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )
        outputs.append(
            {
                "mode": mode,
                "returncode": completed.returncode,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
            }
        )

    guards: list[dict[str, Any]] = []
    matching = acceptance.get("token_matching", "exact")
    if "configuration_required_tokens" in acceptance:
        guards.append(
            run_source_guard(
                {
                    "required_tokens": acceptance["configuration_required_tokens"],
                    "token_matching": matching,
                },
                REPO_ROOT,
            )
        )
    if "build_required_tokens" in acceptance:
        guards.append(
            run_source_guard(
                {
                    "required_tokens": acceptance["build_required_tokens"],
                    "token_matching": matching,
                },
                build_dir,
            )
        )

    # Infrastructure/configuration errors outrank evaluated failures.  This is
    # the same four-state ordering used by the aggregate runner and prevents a
    # missing generated definition from being reported as a physics mismatch.
    statuses = [guard["status"] for guard in guards]
    if "ERROR" in statuses:
        status = "ERROR"
    elif "FAIL" in statuses or any(item["returncode"] != 0 for item in outputs):
        status = "FAIL"
    else:
        status = "PASS"
    return {
        "status": status,
        "compile": compile_result,
        "probes": outputs,
        "configuration_guards": guards,
    }


def _sha256(path: Path) -> str:
    """Hash an external kernel without modifying the authoritative file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_rotating_frame(
    acceptance: dict[str, Any], build_dir: Path, output_root: Path
) -> dict[str, Any]:
    """Run U04 and preserve its required frame, metric, and kernel artifacts.

    Kernel presence alone is never a PASS.  Before CSPICE is called, this
    routine distinguishes an unavailable declared file (SKIPPED) from a file
    whose bytes disagree with the frozen contract (ERROR).  Only the exact
    qualified kernel set is then allowed to reach the linked production probe.
    This ordering also prevents a CSPICE abort from obscuring the scientifically
    important distinction between missing input and corrupted/substituted input.
    """
    test_output = output_root / acceptance["id"]
    test_output.mkdir(parents=True, exist_ok=True)

    kernel_root = Path(
        os.environ.get("MOON_SPICE_KERNEL_ROOT", "/home/vtenishe/SPICE/Kernels")
    )
    kernel_records: list[dict[str, Any]] = []
    missing: list[str] = []
    mismatches: list[dict[str, str]] = []
    expected_hashes = acceptance["data_hashes"]
    for relative_name in acceptance["kernel_files"]:
        path = kernel_root / relative_name
        if not path.is_file():
            missing.append(str(path))
            continue
        actual_hash = _sha256(path)
        expected_hash = expected_hashes.get(relative_name)
        kernel_records.append(
            {
                "relative_path": relative_name,
                "absolute_path": str(path),
                "bytes": path.stat().st_size,
                "sha256": actual_hash,
                "expected_sha256": expected_hash,
                "hash_matches_contract": actual_hash == expected_hash,
            }
        )
        if expected_hash is None or actual_hash != expected_hash:
            mismatches.append(
                {
                    "relative_path": relative_name,
                    "expected_sha256": expected_hash or "undeclared",
                    "actual_sha256": actual_hash,
                }
            )
    kernel_artifact = {
        "kernel_root": str(kernel_root),
        "files": kernel_records,
        "missing": missing,
        "hash_mismatches": mismatches,
    }
    with (test_output / "kernel_hashes.json").open("w", encoding="utf-8") as stream:
        json.dump(kernel_artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")

    if missing:
        return {
            "status": "SKIPPED",
            "reason": "one or more declared authoritative SPICE kernels are missing",
            "missing_kernels": missing,
        }
    if mismatches:
        return {
            "status": "ERROR",
            "reason": "one or more SPICE kernels disagree with the frozen SHA-256 contract",
            "kernel_hash_mismatches": mismatches,
        }

    result = run_probe_case(acceptance, build_dir, output_root)

    # Convert the probe's stable scalar records into a machine-readable metric
    # artifact.  Each vector component remains explicit; no norm can hide a
    # sign or component-order failure.
    metric_pattern = re.compile(
        r"^(PASS|FAIL) (\S+) actual=([+\-0-9.eE]+) expected=([+\-0-9.eE]+)$"
    )
    metrics: list[dict[str, Any]] = []
    for probe in result.get("probes", []):
        for line in probe.get("stdout", "").splitlines():
            match = metric_pattern.match(line)
            if match:
                metrics.append(
                    {
                        "status": match.group(1),
                        "name": match.group(2),
                        "actual": float(match.group(3)),
                        "expected": float(match.group(4)),
                    }
                )
    with (test_output / "term_vectors.json").open("w", encoding="utf-8") as stream:
        json.dump({"epoch_utc": acceptance["epoch_utc"], "metrics": metrics},
                  stream, indent=2, sort_keys=True)
        stream.write("\n")

    if result["status"] == "PASS" and not metrics:
        return {
            **result,
            "status": "ERROR",
            "reason": "U04 probe produced no parseable scientific metrics",
        }
    return result


def run_lola_geometry(
    acceptance: dict[str, Any], build_dir: Path, output_root: Path
) -> dict[str, Any]:
    """Run U08 locally, then apply the independent D01 readiness gate.

    Kernel PASS proves decoding/geometry behavior only.  The enclosing U08
    result remains SKIPPED until the selected external package documents its
    provenance and QA status, as required by AGENTS.md.
    """
    kernel = run_probe_case(acceptance, build_dir, output_root)
    data_root = Path(
        os.environ.get("MOON_VALIDATION_DATA_ROOT", "/data/vtenishe/moon_validation_data")
    )
    package = data_root / "lola"
    required_qualification = [
        package / "README.md",
        package / "provenance.json",
        package / "qa_report.json",
    ]
    missing = [str(path) for path in required_qualification if not path.is_file()]
    if kernel["status"] in {"ERROR", "FAIL"}:
        return kernel
    if missing:
        return {
            "status": "SKIPPED",
            "reason": (
                "LOLA production kernels passed their local checks, but D01 "
                "is not a qualified validation-data package"
            ),
            "kernel_check": kernel,
            "missing_qualification_files": missing,
            "validation_status": "NOT VALIDATED",
        }
    return kernel


def run_na_radiation_pressure(
    acceptance: dict[str, Any], build_dir: Path, output_root: Path
) -> dict[str, Any]:
    """Qualify U05 against an independently digitized publication curve.

    The C++ executable remains a thin adapter around the compiled production
    function.  Expected accelerations are loaded from a reproducibly generated
    CSV derived from Combi et al. (1997), Figure 7, rather than copied from
    ``src/species/Na.cpp``.  The reference uncertainty was frozen from raster
    calibration and line thickness before this comparison is evaluated.
    """
    invariant_result = run_probe_case(acceptance, build_dir, output_root)
    if invariant_result["status"] != "PASS":
        return invariant_result

    test_output = output_root / acceptance["id"]
    test_output.mkdir(parents=True, exist_ok=True)
    reference_dir = test_directories()[acceptance["id"]] / "reference"

    # A changed reference file without a corresponding reviewed contract is a
    # provenance failure, not a physics mismatch.  Classify it as ERROR before
    # invoking the production kernel so altered evidence cannot be scored.
    reference_hashes: list[dict[str, Any]] = []
    hash_errors: list[dict[str, str]] = []
    for relative_name, expected_hash in acceptance["data_hashes"].items():
        path = reference_dir / relative_name
        if not path.is_file():
            hash_errors.append(
                {"file": relative_name, "reason": "missing reference file"}
            )
            continue
        actual_hash = _sha256(path)
        matches = actual_hash == expected_hash
        reference_hashes.append(
            {
                "file": relative_name,
                "sha256": actual_hash,
                "expected_sha256": expected_hash,
                "matches_contract": matches,
            }
        )
        if not matches:
            hash_errors.append(
                {
                    "file": relative_name,
                    "expected_sha256": expected_hash,
                    "actual_sha256": actual_hash,
                }
            )

    if hash_errors:
        return {
            "status": "ERROR",
            "reason": "U05 reference package is missing or fails SHA-256 verification",
            "reference_hashes": reference_hashes,
            "reference_errors": hash_errors,
            "invariant_check": invariant_result,
        }

    qa_path = reference_dir / acceptance["reference_qa_file"]
    try:
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return {
            "status": "ERROR",
            "reason": f"cannot read U05 digitization QA: {error}",
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }
    if qa.get("status") != "PASS":
        return {
            "status": "ERROR",
            "reason": "U05 digitization package QA is not PASS",
            "reference_qa": qa,
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }

    reference_path = reference_dir / acceptance["reference_values_file"]
    try:
        with reference_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        references = [
            {
                "velocity_m_s": 1000.0 * float(row["velocity_km_s"]),
                "expected_cm_s2": float(row["acceleration_mean_cm_s2"]),
                "uncertainty_cm_s2": float(row["uncertainty_cm_s2"]),
            }
            for row in rows
        ]
    except (OSError, KeyError, TypeError, ValueError) as error:
        return {
            "status": "ERROR",
            "reason": f"cannot parse U05 digitized reference: {error}",
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }
    if not references:
        return {
            "status": "ERROR",
            "reason": "U05 digitized reference contains no comparison points",
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }

    executable, compile_result = compile_probe(build_dir, output_root)
    if executable is None:
        return {
            **compile_result,
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }
    command = [str(executable), "radiation-pressure-values"] + [
        f"{item['velocity_m_s']:.17g}" for item in references
    ]
    completed = subprocess.run(
        command, cwd=REPO_ROOT, text=True, capture_output=True, check=False
    )
    if completed.returncode != 0:
        return {
            "status": "ERROR" if completed.returncode == 2 else "FAIL",
            "reason": "production radiation-pressure value probe did not complete",
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }

    value_pattern = re.compile(
        r"^VALUE sodium_radiation_pressure "
        r"velocity_m_s=([+\-0-9.eE]+) acceleration_m_s2=([+\-0-9.eE]+)$"
    )
    actual_values: list[tuple[float, float]] = []
    for line in completed.stdout.splitlines():
        match = value_pattern.match(line)
        if match:
            actual_values.append((float(match.group(1)), float(match.group(2))))
    if len(actual_values) != len(references):
        return {
            "status": "ERROR",
            "reason": "production probe returned an unexpected number of U05 values",
            "expected_count": len(references),
            "actual_count": len(actual_values),
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "reference_hashes": reference_hashes,
            "invariant_check": invariant_result,
        }

    metrics: list[dict[str, Any]] = []
    for reference, (actual_velocity, actual_m_s2) in zip(
        references, actual_values
    ):
        expected_velocity = reference["velocity_m_s"]
        if abs(actual_velocity - expected_velocity) > 1.0e-9:
            return {
                "status": "ERROR",
                "reason": "production probe changed the requested U05 velocity",
                "requested_velocity_m_s": expected_velocity,
                "reported_velocity_m_s": actual_velocity,
                "reference_hashes": reference_hashes,
                "invariant_check": invariant_result,
            }
        actual_cm_s2 = 100.0 * actual_m_s2
        absolute_error = abs(actual_cm_s2 - reference["expected_cm_s2"])
        passed = absolute_error <= reference["uncertainty_cm_s2"]
        metrics.append(
            {
                "status": "PASS" if passed else "FAIL",
                "velocity_m_s": expected_velocity,
                "actual_cm_s2": actual_cm_s2,
                "expected_cm_s2": reference["expected_cm_s2"],
                "absolute_error_cm_s2": absolute_error,
                "acceptance_uncertainty_cm_s2": reference[
                    "uncertainty_cm_s2"
                ],
            }
        )

    comparison = {
        "source": acceptance["independent_reference"],
        "reference_values_file": str(reference_path),
        "reference_hashes": reference_hashes,
        "metrics": metrics,
    }
    with (test_output / "radiation_pressure_reference_comparison.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(comparison, stream, indent=2, sort_keys=True)
        stream.write("\n")

    return {
        "status": "FAIL" if any(item["status"] == "FAIL" for item in metrics) else "PASS",
        "invariant_check": invariant_result,
        "absolute_curve_check": comparison,
        "reference_qa": qa,
        "probe": {
            "command": command,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "returncode": completed.returncode,
        },
    }


def run_source_guard(
    acceptance: dict[str, Any], base_directory: Path = REPO_ROOT
) -> dict[str, Any]:
    """Check declared tokens relative to an explicit source or build root.

    Source-only guards use ``REPO_ROOT``.  Linked probes may pass their exact
    generated build directory, which makes the check compatible with the
    runner's ``--build-dir`` option and avoids mistaking a stale default build
    for the configuration under test.
    """
    checks: list[dict[str, Any]] = []
    matching = acceptance.get("token_matching", "exact")
    if matching not in {"exact", "whitespace-normalized"}:
        return {
            "status": "ERROR",
            "reason": f"unknown source-guard token_matching mode: {matching}",
        }
    for relative_name, required_tokens in acceptance["required_tokens"].items():
        path = base_directory / relative_name
        if not path.is_file():
            checks.append(
                {"file": relative_name, "status": "ERROR", "reason": "missing file"}
            )
            continue
        text = path.read_text(encoding="utf-8")
        # Callback assignments are often wrapped by the formatter.  A contract
        # may explicitly request whitespace normalization so the guard checks
        # the complete C++ token sequence without depending on line layout.
        searchable_text = " ".join(text.split()) if matching == "whitespace-normalized" else text
        for token in required_tokens:
            searchable_token = (
                " ".join(token.split())
                if matching == "whitespace-normalized"
                else token
            )
            checks.append(
                {
                    "file": relative_name,
                    "token": token,
                    "matching": matching,
                    "status": "PASS" if searchable_token in searchable_text else "FAIL",
                }
            )

    if any(item["status"] == "ERROR" for item in checks):
        status = "ERROR"
    elif any(item["status"] == "FAIL" for item in checks):
        status = "FAIL"
    else:
        status = "PASS"
    return {"status": status, "checks": checks}


def run_one(
    test_id: str,
    build_dir: Path = DEFAULT_BUILD,
    output_root: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    """Dispatch one contract and atomically describe its semantic outcome."""
    acceptance = load_acceptance(test_id)
    kind = acceptance["implementation"]

    # Dispatch labels are explicit in acceptance.json.  Tests do not infer
    # readiness from a filename, a downloaded directory, or another test.
    if kind == "linked-production-probe":
        result = run_probe_case(acceptance, build_dir, output_root)
    elif kind == "na-radiation-pressure-probe":
        result = run_na_radiation_pressure(acceptance, build_dir, output_root)
    elif kind == "rotating-frame-probe":
        result = run_rotating_frame(acceptance, build_dir, output_root)
    elif kind == "lola-geometry-probe":
        result = run_lola_geometry(acceptance, build_dir, output_root)
    elif kind == "source-guard":
        result = run_source_guard(acceptance)
    elif kind == "not-implemented":
        result = {
            "status": "SKIPPED",
            "reason": acceptance["skip_reason"],
            "missing_capability": acceptance.get("missing_capability"),
            "required_data": acceptance.get("required_data", []),
        }
    else:
        result = {"status": "ERROR", "reason": f"unknown implementation: {kind}"}

    result.update(
        {
            "id": test_id,
            "name": acceptance["name"],
            "source_revision": git_revision(),
            "test_contract": str(
                test_directories()[test_id] / "reference" / "acceptance.json"
            ),
        }
    )

    test_output = output_root / test_id
    test_output.mkdir(parents=True, exist_ok=True)
    with (test_output / "result.json").open("w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


def run_one_cli(test_id: str) -> int:
    """CLI used by each small per-ID test.py wrapper."""
    import argparse

    parser = argparse.ArgumentParser(description=f"Run srcMoon {test_id}")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run_one(test_id, args.build_dir.resolve(), args.output_dir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return STATUS_EXIT[result["status"]]


if __name__ == "__main__":
    raise SystemExit("Use srcMoon/test/run_tests.py or an individual test.py")
