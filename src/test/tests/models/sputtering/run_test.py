#!/usr/bin/env python3
"""Build and run the dependency-light production sputtering-law fixture."""

from __future__ import annotations

import re
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result,
)


EXPECTED_CASES = [
    "SPUTTERING-TEOLIS-TABLE-001",
    "SPUTTERING-RUBIN-TABLE-001",
    "SPUTTERING-PRODUCT-CATALOGUE-001",
    "SPUTTERING-SPEED-WEIGHT-001",
    "SPUTTERING-WEIGHT-NORMALIZATION-001",
    "SPUTTERING-ANGLE-SELECTOR-001",
    "SPUTTERING-ANGLE-MOMENTS-001",
    "SPUTTERING-SOURCE-RATE-001",
    "SPUTTERING-PRODUCT-LEDGER-001",
    "SPUTTERING-RATE-LEDGER-CLOSURE-001",
]


def _case_results(requested: list[str], stdout: str) -> list[dict]:
    """Enforce one ordered receipt for every requested stable case ID."""
    found = re.findall(
        r"^\[([A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
        stdout,
        flags=re.MULTILINE,
    )
    observed = [identifier for identifier, _, _ in found]
    if observed == requested and len(set(observed)) == len(observed):
        return [case(identifier, status, message or "ok")
                for identifier, status, message in found]

    results = []
    for identifier in requested:
        matches = [record for record in found if record[0] == identifier]
        if len(matches) == 1:
            results.append(case(identifier, matches[0][1],
                                matches[0][2] or "ok"))
        else:
            results.append(case(
                identifier,
                "FAIL",
                "expected one exact record; observed {}".format(len(matches)),
            ))
    results[0] = case(
        results[0]["id"],
        "FAIL",
        "record contract violation: expected {} observed {}".format(
            requested, observed),
    )
    return results


def main() -> int:
    context = load_context("Test ice sputtering yield and energy laws")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    suite_dir = Path(__file__).resolve().parent
    executable = work / "test_sputtering"
    timeout = int(context["test"].get("timeout_seconds", 120))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    unknown = [identifier for identifier in requested
               if identifier not in EXPECTED_CASES]
    if not requested or unknown or len(set(requested)) != len(requested):
        raise RunnerError("invalid sputtering case selection: {}".format(
            requested))

    try:
        # Compile the real implementation with a narrow PIC dependency stub.
        # The source uses SI kg, m/s, J and eV exactly as in a configured AMPS
        # build, while avoiding MPI/mesh/particle state that cannot affect
        # these pure table and distribution laws.
        compiled = run_command(command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
            "-I", str(suite_dir / "stubs"),
            "-I", str(repo / "src" / "models" / "sputtering"),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"),
            str(suite_dir / "test_sputtering.cpp"),
            str(repo / "src" / "models" / "sputtering" / "Sputtering.cpp"),
            str(repo / "src" / "models" / "sputtering" /
                "Sputtering_Ice.cpp"),
            "-o", str(executable),
        ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context,
                "FAIL",
                "production sputtering fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested],
            )

        argv = [str(executable), "--order", "normal"]
        for identifier in requested:
            argv.extend(["--case", identifier])
        run = run_command(argv, work, timeout=timeout)
        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "sputtering-output.txt"
        transcript.write_text(run.stdout or "", encoding="utf-8")

        # This file is intentionally explicit about scope.  src/models/
        # sputtering has no angle sampler and no incident-flux/surface-area
        # source-normalization API.  Those cannot be marked tested merely
        # because yield multiplicity and emission energy passed here.
        gaps = artifacts / "sputtering-capability-gaps.txt"
        gaps.write_text(
            "SUPPORTED AND TESTED\n"
            "- Rubin and Teolis ice yield lookup, thresholds, and interpolation\n"
            "- Teolis product identities through the configured master API\n"
            "- H2O/O2 log-speed laws and relative statistical weights\n"
            "- normalization of those weights over each finite speed window\n"
            "- selectable Knudsen-cosine and uniform-hemisphere directions\n"
            "- incident flux x area x yield physical and numerical event rates\n"
            "- product identity, mass, momentum, kinetic-energy, and rate ledgers\n\n"
            "NOT YET INTEGRATED THROUGH src/models/sputtering\n"
            "- PIC particle allocation, AMR ownership, or mover insertion\n"
            "- MPI reduction of application sampling buffers\n"
            "- migration of application-local legacy angular samplers\n\n"
            "Those integration paths require separate native/application suites; "
            "this fixture does not infer them from dependency-light closure.\n",
            encoding="utf-8",
        )
        results = _case_results(requested, run.stdout or "")
        status = "PASS" if run.returncode == 0 and all(
            item["status"] == "PASS" for item in results) else "FAIL"
        return write_result(
            context,
            status,
            "{} sputtering law, angle, source-rate, and ledger cases "
            "reported".format(len(results)),
            cases=results,
            evidence=[str(transcript), str(gaps)],
        )
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
