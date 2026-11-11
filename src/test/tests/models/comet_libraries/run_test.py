#!/usr/bin/env python3
"""Run imported production-library tests without configuring native AMPS.

Each case compiles the maintained module's own independent reference tests in
the orchestrator's scratch directory. No duplicate physics implementation is
introduced here. The fluorescence provenance case also authenticates the
reviewed immutable tables; it never downloads or regenerates scientific data.
"""
from __future__ import annotations

import sys
from pathlib import Path

TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))
from lib.runner_support import (case, command, error_result, load_context,
                               main_guard, run_command, write_result)

CASES = {
    "COMET-REACTION-NETWORK-001": (
        "src/models/chemical_solver", ["test_reaction_network.cpp"]),
    "COMET-COLLISION-DATABASE-001": (
        "src/models/cometary_physics", ["test_collision_database.cpp"]),
    "COMET-RADIATIVE-COOLING-001": (
        "src/models/cometary_physics", ["test_radiative_cooling.cpp"]),
    "COMET-FLUORESCENCE-001": (
        "src/models/resonance_fluorescence",
        ["resonance_fluorescence.cpp", "tests/test_resonance_fluorescence.cpp"]),
    "COMET-REFERENCE-DATA-001": (
        "src/models/resonance_fluorescence", []),
}


def main():
    context = load_context("Validate reusable comet chemistry and radiation libraries")
    repo, work = Path(context["repo_root"]), Path(context["work_dir"])
    requested = context["test"].get("selected_cases", list(CASES))
    timeout = int(context["test"].get("timeout_seconds", 300))
    results = []
    try:
        for identifier in requested:
            print("[CASE] RUN  " + identifier, flush=True)
            directory, sources = CASES[identifier]
            module = repo / directory
            if sources:
                executable = work / identifier
                build = run_command(command(context, "cxx") + [
                    "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                    "-I", str(module), *[str(module / p) for p in sources],
                    "-o", str(executable)], repo, timeout=timeout)
                result = build if build.returncode else run_command(
                    [str(executable)], work, timeout=timeout)
                passed = result.returncode == 0
            else:
                # Both commands are read-only. A stale generated include or a
                # changed authoritative data hash must make this case fail.
                checks = [run_command(command(context, "python") + [str(module / p)] + args,
                                      work, timeout=timeout)
                          for p, args in [
                              ("tools/generate_reference_tables.py", ["--check"]),
                              ("tests/test_reference_data.py", [])]]
                passed = all(r.returncode == 0 for r in checks)
            results.append(case(identifier, "PASS" if passed else "FAIL",
                                "maintained module tests"))
        return write_result(context,
                            "PASS" if all(r["status"] == "PASS" for r in results) else "FAIL",
                            "reusable comet-library checks", cases=results)
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
