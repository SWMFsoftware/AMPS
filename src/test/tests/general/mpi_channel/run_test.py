#!/usr/bin/env python3
"""Compile and exercise CMPI_channel with an actual two-rank MPI job."""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    case, command, error_result, load_context, main_guard, run_command, write_result)


EXPECTED_CASES = [
    "GENERAL-MPI-CHANNEL-SCALAR-001",
    "GENERAL-MPI-CHANNEL-VECTOR-001",
    "GENERAL-MPI-CHANNEL-STRUCT-001",
    "GENERAL-MPI-CHANNEL-ZERO-001",
    "GENERAL-RNG-RANK-STREAM-001",
]


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    """Provide the relative generated-config include without touching the tree."""

    source = repo / "src" / "general"
    target = work / "header-mirror" / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    # Compile rnd.cpp from the mirror so its global.h include resolves against
    # this test's empty generated configuration, not whichever application a
    # developer happened to configure in the source checkout.
    shutil.copyfile(str(source / "rnd.cpp"), str(target / "rnd.cpp"))
    (work / "header-mirror" / ".general.conf").write_text("", encoding="utf-8")
    return target


def main() -> int:
    context = load_context("Test the AMPS MPI byte-buffer channel")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_mpi_channel.cpp")
    fatal_source = Path(__file__).with_name("test_mpichannel_fatal.cpp")
    executable = work / "test_mpi_channel"
    fatal_executable = work / "test_mpichannel_fatal"
    timeout = int(context["test"].get("timeout_seconds", 120))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    try:
        include = prepare_header_mirror(repo, work)
        compiled = run_command(
            command(context, "mpi_cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-DOMPI_SKIP_MPICXX=1", "-DMPI_ON",
                "-I", str(include), str(source), str(include / "rnd.cpp"),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "MPI-channel test did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested])

        fatal_compiled = run_command(
            command(context, "mpi_cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-DOMPI_SKIP_MPICXX=1", "-I", str(include),
                str(fatal_source), "-o", str(fatal_executable),
            ], repo, timeout=timeout)
        if fatal_compiled.returncode != 0:
            return write_result(
                context, "FAIL", "MPI-channel fatal-path fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    fatal_compiled.returncode)) for identifier in requested])

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        observed = []
        transcript_chunks = []
        for identifier in requested:
            print("[CASE] RUN  {}".format(identifier), flush=True)
            fatal_ok = True
            fatal_output = []
            if identifier == "GENERAL-MPI-CHANNEL-ZERO-001":
                # Keep fatal contracts under an existing stable case ID. Each
                # probe is a standalone process, so a correct nonzero exit
                # cannot terminate or deadlock the two-rank nominal receipt.
                fatal_probes = {
                    "negative-count": "element count cannot be negative",
                    "oversized-send": "length of the data vector",
                    "missing-frame": "received an invalid frame length",
                    "count-overflow": "length of the data vector",
                }
                for probe, diagnostic in fatal_probes.items():
                    rejected = run_command(
                        [str(fatal_executable), probe], work, timeout=timeout)
                    fatal_output.append(
                        "----- fatal probe {} exit={} -----\n{}".format(
                            probe, rejected.returncode, rejected.stdout or ""))
                    if (rejected.returncode == 0 or
                            diagnostic not in (rejected.stdout or "")):
                        fatal_ok = False
                mismatched = run_command(command(context, "mpiexec") + [
                    rank_flag, "2", str(executable), "--fatal-probe",
                    "frame-payload-mismatch",
                ], work, timeout=timeout)
                fatal_output.append(
                    "----- fatal MPI probe frame-payload-mismatch exit={} -----\n{}".format(
                        mismatched.returncode, mismatched.stdout or ""))
                if (mismatched.returncode == 0 or
                        "frame header/payload length mismatch" not in
                        (mismatched.stdout or "")):
                    fatal_ok = False
            run = run_command(command(context, "mpiexec") + [
                rank_flag, "2", str(executable), "--case", identifier,
            ], work, timeout=timeout)
            transcript_chunks.append(
                "===== {} =====\n{}\n{}".format(
                    identifier, "\n".join(fatal_output), run.stdout or ""))
            matches = re.findall(
                r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                run.stdout or "", flags=re.MULTILINE)
            exact = [entry for entry in matches if entry[0] == identifier]
            foreign = [entry[0] for entry in matches if entry[0] != identifier]
            if len(exact) == 1 and not foreign and fatal_ok:
                _, status, message = exact[0]
                if run.returncode != 0 and status == "PASS":
                    status = "FAIL"
                    message = "PASS record with launcher exit {}".format(run.returncode)
                observed.append(case(identifier, status, message or "ok"))
            else:
                observed.append(case(
                    identifier, "FAIL",
                    "expected one exact record and fatal rejection; exact={} "
                    "foreign={} fatal_ok={} exit={}".format(
                        len(exact), foreign, fatal_ok, run.returncode)))

        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "mpi-channel-output.txt"
        transcript.write_text("\n".join(transcript_chunks), encoding="utf-8")
        status = "PASS" if all(item["status"] == "PASS" for item in observed) else "FAIL"
        return write_result(
            context, status, "{} two-rank channel/RNG cases executed".format(
                len(observed)), cases=observed, evidence=[str(transcript)],
            metrics={"mpi_ranks": 2})
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
