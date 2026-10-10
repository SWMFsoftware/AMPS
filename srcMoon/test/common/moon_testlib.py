#!/usr/bin/env python3
"""Shared execution support for the srcMoon U-series verification tests.

The runner keeps four outcomes distinct.  A missing build, production
capability, or qualified dataset is SKIPPED; an infrastructure problem is
ERROR; an evaluated scientific criterion is PASS or FAIL.  In particular,
the library never upgrades source presence or a successful local probe into a
linked-runtime or observational validation result.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
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
    """Execute every probe mode declared by a test's frozen contract."""
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

    status = "PASS" if all(item["returncode"] == 0 for item in outputs) else "FAIL"
    return {"status": status, "compile": compile_result, "probes": outputs}


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


def run_source_guard(acceptance: dict[str, Any]) -> dict[str, Any]:
    """Check declared source tokens without calling that a runtime test."""
    checks: list[dict[str, Any]] = []
    matching = acceptance.get("token_matching", "exact")
    if matching not in {"exact", "whitespace-normalized"}:
        return {
            "status": "ERROR",
            "reason": f"unknown source-guard token_matching mode: {matching}",
        }
    for relative_name, required_tokens in acceptance["required_tokens"].items():
        path = REPO_ROOT / relative_name
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
