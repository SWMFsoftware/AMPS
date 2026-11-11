#!/usr/bin/env python3
"""Build and run dependency-light AMPS general-numerics cases."""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    case, command, error_result, load_context, main_guard, run_command,
    write_result)


EXPECTED_CASES = [
    "GENERAL-NUM-ARRAY-EMPTY-001",
    "GENERAL-NUM-ARRAY-INDEX-001",
    "GENERAL-NUM-ARRAY-OWNERSHIP-001",
    "GENERAL-NUM-ARRAY-RESIZE-001",
    "GENERAL-NUM-ARRAY-BOUNDS-001",
    "GENERAL-NUM-MANAGED-STACK-001",
    "GENERAL-NUM-SPECFUNC-001",
    "GENERAL-NUM-QUADRATURE-POLY-001",
    "GENERAL-NUM-QUADRATURE-CONVERGENCE-001",
    "GENERAL-NUM-QUADRATURE-DOMAIN-001",
    "PERF-CRC-STREAM-001", "PERF-CRC-SPLIT-001",
    "PERF-CRC-COPY-001", "PERF-CRC-THREAD-001",
]


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    """Mirror general headers so the test never creates repository config."""

    source = repo / "src" / "general"
    target = work / "header-mirror" / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    (work / "header-mirror" / ".general.conf").write_text("", encoding="utf-8")
    return target


def main() -> int:
    context = load_context("Test dependency-light AMPS general numerics")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_general_numerics.cpp")
    executable = work / "test_general_numerics"
    mpi_source = Path(__file__).with_name("test_array_3d_mpi.cpp")
    mpi_executable = work / "test_array_3d_mpi"
    timeout = int(context["test"].get("timeout_seconds", 180))
    expected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    try:
        include = prepare_header_mirror(repo, work)
        compiled = run_command(
            command(context, "mpi_cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-pthread",
                "-DOMPI_SKIP_MPICXX=1", "-I", str(include),
                "-I", str(TEST_ROOT / "native"),
                str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "general-numerics fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in expected])

        artifacts.mkdir(parents=True, exist_ok=True)
        observed: Dict[str, Dict] = {}
        transcript: List[str] = []
        # Process-per-case containment makes an allocator/fatal-path defect a
        # typed failure for its own ID without suppressing later evidence.
        for identifier in expected:
            result = run_command(
                [str(executable), "--case", identifier, "--order", "normal"],
                work, timeout=timeout)
            transcript.append("===== {} =====\n{}".format(
                identifier, result.stdout or ""))
            found = re.findall(
                r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                result.stdout or "", flags=re.MULTILINE)
            exact = [entry for entry in found if entry[0] == identifier]
            foreign = [entry[0] for entry in found if entry[0] != identifier]
            if len(exact) == 1 and not foreign:
                _, status, message = exact[0]
                if result.returncode != (0 if status == "PASS" else 1):
                    status = "FAIL"
                    message = "case status/exit mismatch: exit {}".format(
                        result.returncode)
                observed[identifier] = case(identifier, status, message or "ok")
            else:
                observed[identifier] = case(
                    identifier, "FAIL",
                    "expected one exact record; exact={} foreign={} exit={}".format(
                        len(exact), foreign, result.returncode))

        # Keep the public array gather contract under the existing stable
        # array-bounds ID.  The dedicated two-rank executable verifies that a
        # fixed-size array rejects an unrepresentable rank-major result before
        # entering MPI_Gather, while MPI_COMM_SELF is a value-preserving no-op.
        bounds_id = "GENERAL-NUM-ARRAY-BOUNDS-001"
        if bounds_id in expected:
            mpi_compiled = run_command(
                command(context, "mpi_cxx") + [
                    "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                    "-DOMPI_SKIP_MPICXX=1", "-I", str(include),
                    str(mpi_source), "-o", str(mpi_executable),
                ], repo, timeout=timeout)
            mpi_ok = mpi_compiled.returncode == 0
            mpi_output = mpi_compiled.stdout or ""
            mpi_exit = mpi_compiled.returncode
            if mpi_ok:
                rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
                mpi_run = run_command(command(context, "mpiexec") + [
                    rank_flag, "2", str(mpi_executable),
                ], work, timeout=timeout)
                mpi_output = mpi_run.stdout or ""
                mpi_exit = mpi_run.returncode
                mpi_ok = (mpi_run.returncode == 0 and
                          "ARRAY_3D_MPI_GATHER_PASS" in mpi_output)
            transcript.append(
                "===== ARRAY_3D TWO-RANK GATHER CONTRACT =====\n{}".format(
                    mpi_output))
            if not mpi_ok:
                observed[bounds_id] = case(
                    bounds_id, "FAIL",
                    "two-rank array_3d gather contract failed; exit {}".format(
                        mpi_exit))

        transcript_path = artifacts / "general-numerics-output.txt"
        transcript_path.write_text("\n".join(transcript), encoding="utf-8")
        cases = [observed[identifier] for identifier in expected]
        status = "PASS" if all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status, "{} general-numerics cases executed".format(len(cases)),
            cases=cases, evidence=[str(transcript_path)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
