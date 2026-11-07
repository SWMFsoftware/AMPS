#!/usr/bin/env python3
"""Build and run dependency-light AMPS random/distribution cases."""

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
    "GENERAL-RNG-REPRO-001",
    "GENERAL-RNG-ENDPOINT-001",
    "GENERAL-RNG-STREAMS-001",
    "GENERAL-RNG-UNIFORM-001",
    "GENERAL-RNG-NORMAL-001",
    "GENERAL-RNG-GAMMA-001",
    "GENERAL-RNG-DIRECTION-001",
]


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    source = repo / "src" / "general"
    target = work / "header-mirror" / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    (work / "header-mirror" / ".general.conf").write_text("", encoding="utf-8")
    return target


def prepare_hybrid_header_mirror(repo: Path, work: Path) -> Path:
    """Mirror rnd.cpp with an explicit hybrid/OpenMP compile configuration."""

    source = repo / "src" / "general"
    root = work / "hybrid-header-mirror"
    target = root / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    shutil.copyfile(str(source / "rnd.cpp"), str(target / "rnd.cpp"))
    (root / ".general.conf").write_text(
        "#undef _COMPILATION_MODE_\n"
        "#define _COMPILATION_MODE_ _COMPILATION_MODE__HYBRID_\n",
        encoding="utf-8")
    return target


def main() -> int:
    context = load_context("Test dependency-light AMPS random generators")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_general_rng.cpp")
    hybrid_source = Path(__file__).with_name("test_general_rng_hybrid.cpp")
    executable = work / "test_general_rng"
    timeout = int(context["test"].get("timeout_seconds", 180))
    expected = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    try:
        include = prepare_header_mirror(repo, work)
        compiled = run_command(
            command(context, "mpi_cxx") + [
                "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                "-DOMPI_SKIP_MPICXX=1", "-I", str(include),
                "-I", str(TEST_ROOT / "native"),
                str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
                "-o", str(executable),
            ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "general-RNG fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in expected])

        hybrid_executable = None
        hybrid_compile_message = "not selected"
        if "GENERAL-RNG-STREAMS-001" in expected:
            hybrid_include = prepare_hybrid_header_mirror(repo, work)
            hybrid_executable = work / "test_general_rng_hybrid"
            hybrid_compiled = run_command(
                command(context, "mpi_cxx") + [
                    "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                    "-DOMPI_SKIP_MPICXX=1", "-fopenmp",
                    "-I", str(hybrid_include), str(hybrid_source),
                    str(hybrid_include / "rnd.cpp"), "-o",
                    str(hybrid_executable),
                ], repo, timeout=timeout)
            if hybrid_compiled.returncode != 0:
                hybrid_executable = None
                hybrid_compile_message = "hybrid compile exit {}".format(
                    hybrid_compiled.returncode)

        artifacts.mkdir(parents=True, exist_ok=True)
        observed: Dict[str, Dict] = {}
        transcript: List[str] = []
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

            if identifier == "GENERAL-RNG-STREAMS-001":
                if hybrid_executable is None:
                    observed[identifier] = case(
                        identifier, "FAIL", hybrid_compile_message)
                else:
                    hybrid = run_command(
                        [str(hybrid_executable)], work, timeout=timeout,
                        env={"OMP_NUM_THREADS": "4", "OMP_DYNAMIC": "FALSE"})
                    transcript.append("===== {} HYBRID =====\n{}".format(
                        identifier, hybrid.stdout or ""))
                    hybrid_found = re.findall(
                        r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                        hybrid.stdout or "", flags=re.MULTILINE)
                    hybrid_exact = [entry for entry in hybrid_found
                                    if entry[0] == identifier]
                    hybrid_ok = (len(hybrid_exact) == 1 and
                                 hybrid_exact[0][1] == "PASS" and
                                 hybrid.returncode == 0)
                    if (observed[identifier]["status"] != "PASS" or
                            not hybrid_ok):
                        observed[identifier] = case(
                            identifier, "FAIL",
                            "single-thread or four-worker hybrid stream "
                            "contract failed; hybrid exit={}".format(
                                hybrid.returncode))
                    else:
                        observed[identifier] = case(
                            identifier, "PASS",
                            "explicit streams and four-worker hybrid reseed pass")
        transcript_path = artifacts / "general-rng-output.txt"
        transcript_path.write_text("\n".join(transcript), encoding="utf-8")
        cases = [observed[identifier] for identifier in expected]
        status = "PASS" if all(item["status"] == "PASS" for item in cases) else "FAIL"
        return write_result(
            context, status, "{} general-RNG cases executed".format(len(cases)),
            cases=cases, evidence=[str(transcript_path)])
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
