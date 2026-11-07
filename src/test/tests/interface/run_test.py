#!/usr/bin/env python3
"""Compile and run dependency-light C/C++/Fortran interface contract tests."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


CPP_CASES = [
    "INTERFACE-GEOPACK-MATRIX-001",
    "INTERFACE-GEOPACK-CONVENTION-001",
    "INTERFACE-GEOPACK-NOSPICE-001",
]
C_CASES = ["INTERFACE-C-ABI-001"]
FORTRAN_CASES = [
    "INTERFACE-FORTRAN-LEGACY-001",
    "INTERFACE-FORTRAN-ISO-C-001",
]
PROVIDER_CASES = {
    "INTERFACE-PROVIDER-ARMS-001": (
        "interface_arms",
        "configured PIC mesh/state and a versioned ARMS data fixture"),
    "INTERFACE-PROVIDER-BATSRUS-001": (
        "interface_batsrus",
        "BATL/ModReadAmr Fortran libraries and a versioned .idl fixture"),
    "INTERFACE-PROVIDER-CCMC-001": (
        "interface_ccmc",
        "configured PIC particle state and versioned trajectory/background fixtures"),
    "INTERFACE-PROVIDER-KAMELEON-001": (
        "interface_kameleon",
        "CCMC Kameleon libraries and a versioned CDF fixture"),
    "INTERFACE-PROVIDER-SWMF-001": (
        "interface_swmf",
        "an SWMF component/coupling context and independent exchange fixture"),
}
EXPECTED_CASES = C_CASES + FORTRAN_CASES + CPP_CASES + list(PROVIDER_CASES)
CASE_PATTERN = re.compile(
    r"^\[(INTERFACE-[A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
    flags=re.MULTILINE)


def available(argv: Sequence[str]) -> bool:
    if not argv:
        return False
    executable = argv[0]
    if os.path.sep in executable:
        return Path(executable).is_file() and os.access(executable, os.X_OK)
    return shutil.which(executable) is not None


def optional_tool(context: Mapping, names: Sequence[str], environment: str,
                  fallbacks: Sequence[str], capability: str) -> Tuple[Optional[List[str]], str]:
    """Resolve an optional compiler without baking a machine path into a test.

    A declared unavailable capability wins over host auto-discovery.  This lets
    a cluster profile intentionally suppress a compiler installed only on a
    login node.  Otherwise an explicit tool vector, conventional CC/FC
    environment command, and finally common driver names are tried in order.
    """

    declaration = context.get("capabilities", {}).get(capability)
    if isinstance(declaration, dict) and not declaration.get("available", False):
        return None, "capability '{}' is declared unavailable".format(capability)

    tools = context.get("tools", {})
    candidates: List[List[str]] = []
    for name in names:
        value = tools.get(name)
        if isinstance(value, list) and value:
            candidates.append(list(value))
    value = os.environ.get(environment, "").strip()
    if value:
        candidates.append(shlex.split(value))
    candidates.extend([[name] for name in fallbacks])
    for candidate in candidates:
        if available(candidate):
            return candidate, "resolved {}".format(candidate[0])
    return None, "no {} compiler was configured or found".format(capability)


def parse_single(identifier: str, output: str, returncode: int) -> Dict:
    found = CASE_PATTERN.findall(output or "")
    exact = [entry for entry in found if entry[0] == identifier]
    foreign = [entry[0] for entry in found if entry[0] != identifier]
    if len(exact) != 1 or foreign:
        return case(identifier, "FAIL",
                    "expected one exact record; exact={} foreign={} exit={}".format(
                        len(exact), foreign, returncode))
    _, status, message = exact[0]
    if returncode != 0 and status == "PASS":
        return case(identifier, "FAIL", "PASS record with exit {}".format(returncode))
    return case(identifier, status, message or "ok")


def provider_case(identifier: str, context: Mapping) -> Dict:
    """Return typed non-evidence for provider integrations outside this suite.

    Merely finding a library is not an interface test.  These adapters need a
    configured PIC/SWMF context and independently versioned data or exchange
    fixture.  Until such a fixture is registered, both unavailable and
    nominally available providers remain visible SKIPs rather than false PASS.
    """

    capability, requirement = PROVIDER_CASES[identifier]
    declaration = context.get("capabilities", {}).get(capability, {})
    if not isinstance(declaration, dict) or not declaration.get("available", False):
        return case(identifier, "SKIP",
                    "optional capability '{}' unavailable; requires {}".format(
                        capability, requirement))
    return case(identifier, "SKIP",
                "capability '{}' is declared available, but no independent {} "
                "fixture is registered in this dependency-light suite".format(
                    capability, requirement))


def compile_cpp_contracts(context: Mapping, repo: Path, work: Path,
                          timeout: int) -> Tuple[Optional[Path], str]:
    executable = work / "test_geopack_contracts"
    source = Path(__file__).resolve().with_name("test_geopack_contracts.cpp")
    interface = repo / "src" / "interface"
    stubs = Path(__file__).resolve().with_name("geopack_stubs")
    result = run_command(
        command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-D_NO_SPICE_CALLS_",
            "-I", str(stubs), "-I", str(interface),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
            str(interface / "GeopackInterface.cpp"),
            "-o", str(executable),
        ], repo, timeout=timeout)
    if result.returncode != 0:
        return None, "Geopack contract compile exit {}".format(result.returncode)
    return executable, "compiled production Geopack interface in no-SPICE mode"


def compile_c_contract(context: Mapping, repo: Path, work: Path,
                       timeout: int) -> Tuple[Optional[Path], str]:
    compiler, reason = optional_tool(
        context, ("cc", "c"), "CC", ("cc", "gcc", "clang"), "c")
    if compiler is None:
        return None, reason
    source_dir = Path(__file__).resolve().parent
    object_file = work / "c_abi_echo.o"
    compiled = run_command(
        compiler + ["-std=c11", "-O2", "-Wall", "-Wextra", "-Wpedantic",
                    "-Werror", "-c", str(source_dir / "c_abi_echo.c"),
                    "-o", str(object_file)], repo, timeout=timeout)
    if compiled.returncode != 0:
        return None, "C echo compile exit {}".format(compiled.returncode)
    executable = work / "test_c_abi"
    linked = run_command(
        command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"),
            str(source_dir / "test_c_abi.cpp"), str(object_file),
            "-o", str(executable),
        ], repo, timeout=timeout)
    if linked.returncode != 0:
        return None, "C/C++ ABI link exit {}".format(linked.returncode)
    return executable, reason


def compile_fortran_contract(context: Mapping, repo: Path, work: Path,
                             timeout: int) -> Tuple[Optional[Path], str]:
    compiler, reason = optional_tool(
        context, ("fortran", "fc"), "FC",
        ("gfortran", "flang-new", "flang", "ifx", "ifort"), "fortran")
    if compiler is None:
        return None, reason
    source_dir = Path(__file__).resolve().parent
    object_file = work / "fortran_abi_echo.o"
    compiled = run_command(
        # Keep the optional compiler contract to widely shared Unix driver
        # switches; GNU warning flags are not portable to ifort/ifx/flang.
        compiler + ["-O2", "-c",
                    str(source_dir / "fortran_abi_echo.f90"),
                    "-o", str(object_file)], repo, timeout=timeout)
    if compiled.returncode != 0:
        return None, "Fortran echo compile exit {}".format(compiled.returncode)
    executable = work / "test_fortran_abi"
    linked = run_command(
        command(context, "cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"),
            str(source_dir / "test_fortran_abi.cpp"), str(object_file),
            "-o", str(executable),
        ], repo, timeout=timeout)
    if linked.returncode != 0:
        return None, "C++/Fortran ABI link exit {}".format(linked.returncode)
    return executable, reason


def main() -> int:
    context = load_context("Test AMPS C/C++/Fortran interface ABI contracts")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 180))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if not requested or len(set(requested)) != len(requested):
        raise RunnerError("selected_cases must be a unique nonempty list")
    unknown = [identifier for identifier in requested if identifier not in EXPECTED_CASES]
    if unknown:
        raise RunnerError("unknown interface cases: {}".format(unknown))

    work.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)
    try:
        executables: Dict[str, Optional[Path]] = {}
        reasons: Dict[str, str] = {}
        if any(identifier in CPP_CASES for identifier in requested):
            executable, reason = compile_cpp_contracts(context, repo, work, timeout)
            for identifier in CPP_CASES:
                executables[identifier], reasons[identifier] = executable, reason
        if any(identifier in C_CASES for identifier in requested):
            executable, reason = compile_c_contract(context, repo, work, timeout)
            for identifier in C_CASES:
                executables[identifier], reasons[identifier] = executable, reason
        if any(identifier in FORTRAN_CASES for identifier in requested):
            executable, reason = compile_fortran_contract(context, repo, work, timeout)
            for identifier in FORTRAN_CASES:
                executables[identifier], reasons[identifier] = executable, reason

        results: List[Dict] = []
        transcript_parts: List[str] = []
        for identifier in requested:
            if identifier in PROVIDER_CASES:
                result = provider_case(identifier, context)
                print("[CASE] SKIP {} {}".format(identifier, result["message"]),
                      flush=True)
                transcript_parts.append(
                    "===== {} =====\n[{}] SKIP {}\n".format(
                        identifier, identifier, result["message"]))
                results.append(result)
                continue
            executable = executables.get(identifier)
            reason = reasons.get(identifier, "compiler group was not prepared")
            if executable is None:
                # A missing optional C/Fortran compiler is a typed SKIP. A
                # present compiler that failed to compile/link is a test FAIL.
                status = "SKIP" if reason.startswith("no ") or "declared unavailable" in reason else "FAIL"
                print("[CASE] {} {} {}".format(status, identifier, reason),
                      flush=True)
                transcript_parts.append(
                    "===== {} =====\n[{}] {} {}\n".format(
                        identifier, identifier, status, reason))
                results.append(case(identifier, status, reason))
                continue
            run = run_command([str(executable), "--case", identifier,
                               "--order", "normal"], work, timeout=timeout)
            transcript_parts.append(
                "===== {} =====\n{}".format(identifier, run.stdout or ""))
            results.append(parse_single(identifier, run.stdout or "", run.returncode))

        transcript = artifacts / "interface-abi-output.txt"
        transcript.write_text("\n".join(transcript_parts), encoding="utf-8")
        severity = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}
        status = max(results, key=lambda item: severity[item["status"]])["status"]
        return write_result(
            context, status,
            "{} interface ABI cases reported ({} skipped)".format(
                len(results), sum(item["status"] == "SKIP" for item in results)),
            cases=results, evidence=[str(transcript)],
            metrics={
                "c_compiler_available": executables.get(C_CASES[0]) is not None
                    if any(item in requested for item in C_CASES) else None,
                "fortran_compiler_available": executables.get(FORTRAN_CASES[0]) is not None
                    if any(item in requested for item in FORTRAN_CASES) else None,
            })
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
