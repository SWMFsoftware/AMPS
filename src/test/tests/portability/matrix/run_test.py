#!/usr/bin/env python3
"""Run the bounded T12 compiler/sanitizer/backend portability matrix.

The adapter never edits repository configuration and never performs a native
AMPS build.  Compiler and sanitizer profiles all compile one production-linked
``array_1d`` fixture in the suite work directory.  Optional tools retain stable
case identities as SKIP, while a present tool that cannot compile or execute
the asserted contract is FAIL.  The strict-coverage command is intentionally
reported at face value: while the inventory is incomplete its case is FAIL,
so this suite cannot make T12 look closed merely because the failure was
anticipated.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


EXPECTED_CASES = [
    "T12-CXX-GNU-DEBUG-001",
    "T12-CXX-GNU-OPTIMIZED-001",
    "T12-CXX-CLANG-001",
    "T12-FORTRAN-ABI-001",
    "T12-SANITIZER-ASAN-001",
    "T12-SANITIZER-UBSAN-001",
    "T12-SANITIZER-TSAN-001",
    "T12-VALGRIND-LEAK-001",
    "T12-OPENMP-INVARIANCE-001",
    "T12-MPI-OPENMPI-001",
    "T12-MPI-MPICH-001",
    "T12-CUDA-COMPILE-001",
    "T12-CUDA-RUNTIME-001",
    "T12-PERFORMANCE-ARRAY-001",
    "T12-PERFORMANCE-IO-001",
    "T12-COVERAGE-STRICT-001",
]

SEVERITY = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}


def write_nested_machine_profile(context: Mapping[str, object], path: Path) -> None:
    """Mirror the resolved machine settings for a nested global-runner check.

    The public runner requires an explicit machine profile.  T12's strict
    coverage identity recursively invokes that runner, so it must not evade
    the same contract or select a hostname-specific committed file.  Only
    configuration-schema fields are copied: transient control-plane keys such
    as ``full_validation`` are execution intent and are not valid overlay
    input.  Commands are already normalized argv vectors and paths are already
    absolute, making this generated profile an exact record of the parent
    worker rather than a second source of machine assumptions.
    """

    execution = context["execution"]
    allowed_execution = (
        "build_jobs", "mpi_ranks", "mpi_rank_flag",
        "default_timeout_seconds", "keep_work", "default_campaign")
    payload = {
        "schema_version": 1,
        "paths": dict(context["paths"]),
        "tools": dict(context["tools"]),
        "execution": {
            key: execution[key] for key in allowed_execution if key in execution},
        "capabilities": dict(context.get("capabilities", {})),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def available(argv: Sequence[str]) -> bool:
    """Resolve argv[0] without invoking a shell or trusting its arguments."""

    if not argv:
        return False
    executable = os.path.expandvars(os.path.expanduser(argv[0]))
    if os.path.sep in executable:
        return Path(executable).is_file() and os.access(executable, os.X_OK)
    return shutil.which(executable) is not None


def optional(context: Mapping, tool_name: str,
             fallbacks: Sequence[str]) -> Optional[List[str]]:
    """Prefer a machine-profile tool vector, then portable executable names."""

    configured = context.get("tools", {}).get(tool_name)
    candidates: List[List[str]] = []
    if isinstance(configured, list) and configured:
        candidates.append(list(configured))
    candidates.extend([[name] for name in fallbacks])
    for candidate in candidates:
        if available(candidate):
            return candidate
    return None


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    """Mirror public general headers and an empty generated configuration.

    ``global.h`` includes ``../../.general.conf``.  A private empty file makes
    the dependency-light profile independent of whichever application is
    configured in the shared checkout, and avoids a repository write.
    """

    source = repo / "src/general"
    target = work / "header-mirror/src/general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copy2(str(item), str(target / item.name))
    (work / "header-mirror/.general.conf").write_text("", encoding="utf-8")
    return target


def compile_fixture(compiler: Sequence[str], include: Path, source: Path,
                    executable: Path, flags: Sequence[str], repo: Path,
                    timeout: int, environment: Optional[Mapping[str, str]] = None
                    ) -> Tuple[bool, str]:
    """Compile a profile and preserve a concise diagnostic for its case."""

    try:
        result = run_command(
            list(compiler) + ["-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
                              "-Werror", "-DOMPI_SKIP_MPICXX=1", "-fopenmp",
                              "-I", str(include)] + list(flags) +
            [str(source), "-o", str(executable)],
            repo, timeout=timeout, env=environment)
    except RunnerError as error:
        return False, str(error)
    detail = " ".join((result.stdout or "").split())[-500:]
    return result.returncode == 0, "compiler exit {}; {}".format(
        result.returncode, detail)


def execute(executable: Path, arguments: Sequence[str], work: Path,
            timeout: int, environment: Optional[Mapping[str, str]] = None
            ) -> Tuple[Optional[int], str]:
    """Return a typed timeout marker rather than aborting later matrix cases."""

    try:
        result = run_command([str(executable)] + list(arguments), work,
                             timeout=timeout, env=environment)
    except RunnerError as error:
        return None, str(error)
    return result.returncode, result.stdout or ""


def compile_and_run_profile(identifier: str, compiler: Sequence[str],
                            include: Path, source: Path, work: Path, repo: Path,
                            timeout: int, flags: Sequence[str],
                            environment: Optional[Mapping[str, str]] = None,
                            run_arguments: Sequence[str] = ("correctness",),
                            pass_marker: str = "CORRECTNESS status=PASS",
                            ) -> Tuple[Dict, Optional[Path], str]:
    executable = work / identifier.lower()
    ok, diagnostic = compile_fixture(
        compiler, include, source, executable, flags, repo, timeout, environment)
    if not ok:
        return case(identifier, "FAIL", "production array fixture compile failed: " + diagnostic), None, diagnostic
    returncode, output = execute(
        executable, list(run_arguments), work, timeout, environment)
    if returncode is None:
        return case(identifier, "ERROR", output), executable, output
    passed = returncode == 0 and pass_marker in output
    success_message = (
        "four-thread OpenMP production-array path passed under {}".format(
            compiler[0]) if "OPENMP" in pass_marker else
        "production array ownership/checksum under {}".format(compiler[0]))
    return (case(identifier, "PASS" if passed else "FAIL",
                 success_message if passed else
                 "profile exit={} output={}".format(
                     returncode, " ".join(output.split())[-300:])),
            executable, output)


def fortran_case(context: Mapping, repo: Path, work: Path,
                 timeout: int) -> Tuple[Dict, str]:
    identifier = "T12-FORTRAN-ABI-001"
    compiler = optional(context, "fortran", ("gfortran", "flang-new", "flang", "ifx", "ifort"))
    if compiler is None:
        return case(identifier, "SKIP", "no supported Fortran compiler configured or found"), ""
    source_dir = TEST_ROOT / "tests/interface"
    object_file = work / "fortran_abi_echo.o"
    try:
        compiled = run_command(compiler + ["-O2", "-c",
            str(source_dir / "fortran_abi_echo.f90"), "-o", str(object_file)],
            repo, timeout=timeout)
        if compiled.returncode != 0:
            return case(identifier, "FAIL", "Fortran source compile exit {}".format(
                compiled.returncode)), compiled.stdout or ""
        # Compile C++ with the configured C++ driver, then let the selected
        # Fortran driver perform the final link. The latter automatically
        # selects its own runtime (gfortran, Flang or Intel) and avoids the
        # false GNU-only ``-lgfortran`` assumption. AMPS' configured C++
        # boundary uses libstdc++, which is named explicitly when a Fortran
        # driver owns the final link line.
        main_object = work / "fortran_abi_main.o"
        test_object = work / "fortran_abi_test.o"
        cpp_flags = [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-I", str(TEST_ROOT / "native")]
        outputs = [compiled.stdout or ""]
        for source, target in (
                (TEST_ROOT / "native/test_main.cpp", main_object),
                (source_dir / "test_fortran_abi.cpp", test_object)):
            cpp = run_command(command(context, "cxx") + cpp_flags + [
                "-c", str(source), "-o", str(target)], repo, timeout=timeout)
            outputs.append(cpp.stdout or "")
            if cpp.returncode != 0:
                return case(identifier, "FAIL", "C++ ABI fixture compile exit {}".format(
                    cpp.returncode)), "\n".join(outputs)
        executable = work / "fortran_abi"
        linked = run_command(list(compiler) + [
            str(main_object), str(test_object), str(object_file),
            "-lstdc++", "-o", str(executable)], repo, timeout=timeout)
        outputs.append(linked.stdout or "")
        if linked.returncode != 0:
            return case(identifier, "FAIL", "C++/Fortran ABI link exit {}".format(
                linked.returncode)), "\n".join(outputs)
        passed = True
        for local_id in ("INTERFACE-FORTRAN-LEGACY-001", "INTERFACE-FORTRAN-ISO-C-001"):
            result = run_command([str(executable), "--case", local_id], work,
                                 timeout=timeout)
            outputs.append(result.stdout or "")
            passed = passed and result.returncode == 0 and (
                "[{}] PASS".format(local_id) in (result.stdout or ""))
        return case(identifier, "PASS" if passed else "FAIL",
                    "legacy underscore and ISO_C_BINDING ABI both passed"
                    if passed else "one or both Fortran ABI controls failed"), "\n".join(outputs)
    except RunnerError as error:
        return case(identifier, "ERROR", str(error)), str(error)


def mpi_case(identifier: str, compiler: Sequence[str], launcher: Sequence[str],
             repo: Path, work: Path, include: Path, timeout: int) -> Tuple[Dict, str]:
    executable = work / identifier.lower()
    source = TEST_ROOT / "tests/general/numerics/test_array_3d_mpi.cpp"
    try:
        built = run_command(list(compiler) + ["-std=c++17", "-O2", "-Wall", "-Wextra",
                             "-Wpedantic", "-Werror", "-DOMPI_SKIP_MPICXX=1",
                             "-I", str(include), str(source), "-o", str(executable)],
                            repo, timeout=timeout)
        if built.returncode != 0:
            return case(identifier, "FAIL", "{} compile exit {}".format(
                compiler[0], built.returncode)), built.stdout or ""
        run = run_command(list(launcher) + ["-n", "2", str(executable)], work,
                          timeout=timeout)
        output = run.stdout or ""
        restricted = ("Operation not permitted" in output or
                      "No network interfaces were found" in output or
                      "unable to create PMI port" in output)
        if run.returncode != 0 and restricted:
            return case(identifier, "SKIP",
                        "{} is installed but sandbox/scheduler socket policy prevents a two-rank launch".format(
                            launcher[0])), output
        passed = run.returncode == 0 and "ARRAY_3D_MPI_GATHER_PASS" in output
        return case(identifier, "PASS" if passed else "FAIL",
                    "two-rank production array gather rejection/identity passed"
                    if passed else "{} launcher exit {}".format(
                        launcher[0], run.returncode)), output
    except RunnerError as error:
        return case(identifier, "ERROR", str(error)), str(error)


def cuda_cases(context: Mapping, repo: Path, work: Path, timeout: int,
               run_runtime: bool) -> Tuple[Dict, Optional[Dict], str]:
    compile_id = "T12-CUDA-COMPILE-001"
    runtime_id = "T12-CUDA-RUNTIME-001"
    nvcc = optional(context, "nvcc", ("nvcc",))
    if nvcc is None:
        first = case(compile_id, "SKIP", "nvcc is unavailable")
        second = case(runtime_id, "SKIP", "CUDA runtime cannot be tested without nvcc") if run_runtime else None
        return first, second, ""
    executable = work / "cuda_constants"
    source = Path(__file__).with_name("cuda_constants.cu")
    try:
        built = run_command(list(nvcc) + ["-std=c++14", "-I", str(repo / "src/general"),
                             str(source), "-o", str(executable)], repo,
                            timeout=timeout)
    except RunnerError as error:
        return case(compile_id, "ERROR", str(error)), (
            case(runtime_id, "SKIP", "CUDA compile did not complete") if run_runtime else None), str(error)
    if built.returncode != 0:
        first = case(compile_id, "FAIL", "nvcc compile exit {}".format(
            built.returncode))
        second = case(runtime_id, "SKIP", "CUDA compile failed; runtime not attempted") if run_runtime else None
        return first, second, built.stdout or ""
    first = case(compile_id, "PASS", "production constants header compiled for CUDA")
    if not run_runtime:
        return first, None, built.stdout or ""
    try:
        run = run_command([str(executable)], work, timeout=timeout)
        if run.returncode == 77 and "CUDA_NO_DEVICE" in (run.stdout or ""):
            second = case(runtime_id, "SKIP", "CUDA compiler is present but no usable driver/device is available")
        else:
            passed = run.returncode == 0 and "CUDA_RUNTIME status=PASS" in (run.stdout or "")
            second = case(runtime_id, "PASS" if passed else "FAIL",
                          "kernel executed on an available CUDA device" if passed else
                          "CUDA runtime exit {}".format(run.returncode))
        return first, second, (built.stdout or "") + (run.stdout or "")
    except RunnerError as error:
        # Compilation is independent evidence. Preserve its PASS if device
        # launch times out or fails in the runtime control plane, and type only
        # the runtime identity as ERROR.
        return first, case(runtime_id, "ERROR", str(error)), (
            (built.stdout or "") + str(error))


def main() -> int:
    context = load_context("Run bounded T12 portability/sanitizer matrix")
    try:
        repo = Path(context["repo_root"])
        work = Path(context["work_dir"])
        artifacts = Path(context["artifact_dir"])
        timeout = int(context["test"].get("timeout_seconds", 300))
        # Selection may narrow execution, but the owner manifest must still
        # expose every identity in canonical order. Keep registry validation
        # and directory setup inside the structured receipt boundary so a
        # malformed context or unwritable work root yields per-case ERROR rows
        # rather than an unstructured exit from ``main_guard``.
        declared = context["test"].get("expected_cases")
        if declared != EXPECTED_CASES:
            raise RunnerError(
                "manifest expected_cases must exactly equal the 16-case T12 registry")
        requested = list(context["test"].get(
            "selected_cases", declared))
        if (not requested or len(set(requested)) != len(requested) or
                any(identifier not in EXPECTED_CASES for identifier in requested)):
            raise RunnerError("selected_cases is not a unique subset of the T12 registry")
        work.mkdir(parents=True, exist_ok=True)
        artifacts.mkdir(parents=True, exist_ok=True)
        include = prepare_header_mirror(repo, work)
        source = Path(__file__).with_name("portability_fixture.cpp")
        results: Dict[str, Dict] = {}
        transcripts: List[str] = []

        profile_specs = {
            "T12-CXX-GNU-DEBUG-001": (
                ["-O0", "-g", "-DAMPS_T12_REQUIRE_GNU=1"], None),
            "T12-CXX-GNU-OPTIMIZED-001": (
                ["-O3", "-DNDEBUG", "-DAMPS_T12_REQUIRE_GNU=1"], None),
            "T12-SANITIZER-ASAN-001": (["-O1", "-g", "-fsanitize=address", "-fno-omit-frame-pointer"],
                                        # LeakSanitizer cannot operate under
                                        # ptrace-based supervisors on several
                                        # nightly systems. Address bounds and
                                        # lifetime misuse remain enabled here;
                                        # the independent Valgrind case owns
                                        # the leak contract.
                                        {"ASAN_OPTIONS": "detect_leaks=0:halt_on_error=1"}),
            "T12-SANITIZER-UBSAN-001": (["-O1", "-g", "-fsanitize=undefined", "-fno-sanitize-recover=all"],
                                         {"UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"}),
            "T12-SANITIZER-TSAN-001": (["-O1", "-g", "-fsanitize=thread"],
                                        {"TSAN_OPTIONS": "halt_on_error=1"}),
        }
        # The MPI wrapper is used because production general headers include
        # mpi.h.  On the supported GNU profile it delegates to the same host
        # compiler while supplying only the MPI include/link boundary.
        mpi_compiler = command(context, "mpi_cxx")
        executables: Dict[str, Path] = {}
        for identifier, (flags, environment) in profile_specs.items():
            if identifier not in requested and not (
                    identifier == "T12-CXX-GNU-OPTIMIZED-001" and any(item in requested for item in (
                        "T12-VALGRIND-LEAK-001", "T12-OPENMP-INVARIANCE-001",
                        "T12-PERFORMANCE-ARRAY-001", "T12-PERFORMANCE-IO-001"))):
                continue
            record, executable, output = compile_and_run_profile(
                identifier, mpi_compiler, include, source, work, repo, timeout,
                flags, environment,
                ("openmp", "4") if identifier == "T12-SANITIZER-TSAN-001"
                else ("correctness",),
                "OPENMP status=PASS requested_threads=4" if identifier ==
                "T12-SANITIZER-TSAN-001" else "CORRECTNESS status=PASS")
            if identifier == "T12-SANITIZER-TSAN-001" and record["status"] == "FAIL" and (
                    "unexpected memory mapping" in output.lower() or
                    "threadsanitizer: unexpected" in output.lower() or
                    "cannot find libtsan_preinit.o" in output.lower()):
                record = case(identifier, "SKIP",
                              "ThreadSanitizer runtime/support files are unavailable or incompatible on this host")
            if identifier == "T12-SANITIZER-ASAN-001" and record["status"] == "FAIL" and (
                    "leaksanitizer does not work under ptrace" in output.lower() or
                    "leaksanitizer has encountered a fatal error" in output.lower()):
                record = case(identifier, "SKIP",
                              "LeakSanitizer execution is blocked by this host's ptrace supervision")
            results[identifier] = record
            if executable is not None: executables[identifier] = executable
            transcripts.append("===== {} =====\n{}".format(identifier, output))

        clang_id = "T12-CXX-CLANG-001"
        if clang_id in requested:
            clang = optional(context, "clang_cxx", ("clang++",))
            wrapper = optional(context, "openmpi_cxx", ("mpicxx.openmpi",))
            if clang is None or wrapper is None:
                results[clang_id] = case(clang_id, "SKIP",
                    "clang++ and an OpenMPI compiler wrapper are both required")
            else:
                record, _, output = compile_and_run_profile(
                    clang_id, wrapper, include, source, work, repo, timeout,
                    ["-O2"], {"OMPI_CXX": shlex.join(clang)})
                results[clang_id] = record
                transcripts.append("===== {} =====\n{}".format(clang_id, output))

        if "T12-FORTRAN-ABI-001" in requested:
            record, output = fortran_case(context, repo, work, timeout)
            results[record["id"]] = record
            transcripts.append("===== {} =====\n{}".format(record["id"], output))

        optimized = executables.get("T12-CXX-GNU-OPTIMIZED-001")
        valgrind_id = "T12-VALGRIND-LEAK-001"
        if valgrind_id in requested:
            valgrind = optional(context, "valgrind", ("valgrind",))
            if valgrind is None:
                results[valgrind_id] = case(valgrind_id, "SKIP", "valgrind is unavailable")
            elif optimized is None:
                results[valgrind_id] = case(valgrind_id, "FAIL", "optimized production fixture was not built")
            else:
                # Runtime libraries may intentionally retain process-lifetime
                # bookkeeping. Treat definite, indirect, and possible loss
                # as errors while retaining the complete reachable-byte
                # inventory in the transcript.
                try:
                    run = run_command(list(valgrind) + [
                        "--leak-check=full", "--show-leak-kinds=all",
                        "--errors-for-leak-kinds=definite,indirect,possible",
                        "--error-exitcode=99", str(optimized), "correctness"],
                        work, timeout=timeout)
                    passed = (run.returncode == 0 and
                              "ERROR SUMMARY: 0 errors" in (run.stdout or ""))
                    results[valgrind_id] = case(
                        valgrind_id, "PASS" if passed else "FAIL",
                        "no definite, indirect, or possible production-array leak"
                        if passed else "valgrind exit {}".format(run.returncode))
                    transcripts.append("===== {} =====\n{}".format(
                        valgrind_id, run.stdout or ""))
                except RunnerError as error:
                    # Keep a timeout local to this case so all later matrix
                    # identities are still executed and published.
                    results[valgrind_id] = case(
                        valgrind_id, "ERROR", str(error))
                    transcripts.append("===== {} =====\nERROR: {}".format(
                        valgrind_id, error))

        openmp_id = "T12-OPENMP-INVARIANCE-001"
        if openmp_id in requested:
            if optimized is None:
                results[openmp_id] = case(openmp_id, "FAIL", "optimized OpenMP fixture was not built")
            else:
                outputs = []
                checksums = []
                passed = True
                for threads in (1, 2, 4):
                    code, output = execute(optimized, ["openmp", str(threads)], work, timeout)
                    outputs.append(output)
                    match = re.search(r"checksum=([0-9.eE+-]+)", output)
                    passed = passed and code == 0 and match is not None
                    if match: checksums.append(match.group(1))
                passed = passed and len(set(checksums)) == 1
                results[openmp_id] = case(openmp_id, "PASS" if passed else "FAIL",
                    "1/2/4-thread production-array checksum is identical" if passed else
                    "thread-count checksum/execution mismatch")
                transcripts.append("===== {} =====\n{}".format(openmp_id, "".join(outputs)))

        for identifier, compiler_tool, launcher_tool, compiler_fallback, launcher_fallback in (
                ("T12-MPI-OPENMPI-001", "openmpi_cxx", "openmpi_exec",
                 "mpicxx.openmpi", "mpiexec.openmpi"),
                ("T12-MPI-MPICH-001", "mpich_cxx", "mpich_exec",
                 "mpicxx.mpich", "mpiexec.mpich")):
            if identifier in requested:
                compiler = optional(context, compiler_tool, (compiler_fallback,))
                launcher = optional(context, launcher_tool, (launcher_fallback,))
                if compiler is None or launcher is None:
                    record = case(
                        identifier, "SKIP",
                        "configured/fallback MPI implementation pair is unavailable")
                    output = ""
                else:
                    record, output = mpi_case(
                        identifier, compiler, launcher, repo, work, include, timeout)
                results[identifier] = record
                transcripts.append("===== {} =====\n{}".format(identifier, output))

        wants_cuda_compile = "T12-CUDA-COMPILE-001" in requested
        wants_cuda_runtime = "T12-CUDA-RUNTIME-001" in requested
        if wants_cuda_compile or wants_cuda_runtime:
            compile_record, runtime_record, output = cuda_cases(
                context, repo, work, timeout, wants_cuda_runtime)
            if wants_cuda_compile: results[compile_record["id"]] = compile_record
            if runtime_record is not None: results[runtime_record["id"]] = runtime_record
            transcripts.append("===== CUDA =====\n" + output)

        for identifier, mode, ceiling in (
                ("T12-PERFORMANCE-ARRAY-001", "performance-array", 30.0),
                ("T12-PERFORMANCE-IO-001", "performance-io", 15.0)):
            if identifier not in requested: continue
            if optimized is None:
                results[identifier] = case(identifier, "FAIL", "optimized performance fixture was not built")
                continue
            arguments = [mode]
            if mode == "performance-io": arguments.append(str(work / "restart-like-array.bin"))
            code, output = execute(optimized, arguments, work, timeout)
            match = re.search(r"seconds=([0-9.eE+-]+)", output)
            seconds = float(match.group(1)) if match else float("inf")
            correct = code == 0 and "status=PASS" in output
            passed = correct and seconds <= ceiling
            results[identifier] = case(identifier, "PASS" if passed else "FAIL",
                "correctness passed; elapsed {:.6g} s <= broad {:.1f} s fixture guardrail".format(
                    seconds, ceiling) if passed else
                "correct={} elapsed={} s ceiling={} s".format(correct, seconds, ceiling))
            transcripts.append("===== {} =====\n{}".format(identifier, output))

        coverage_id = "T12-COVERAGE-STRICT-001"
        if coverage_id in requested:
            nested_profile = work / "nested-machine-profile.json"
            write_nested_machine_profile(context, nested_profile)
            coverage = run_command(command(context, "python") + [
                "src/test/run_tests.py", "--config", str(nested_profile),
                "--coverage", "--require-complete-coverage"],
                repo, timeout=timeout)
            passed = coverage.returncode == 0
            results[coverage_id] = case(coverage_id, "PASS" if passed else "FAIL",
                "all inventory components are covered/exempt" if passed else
                "strict coverage remains open (exit {}); see retained inventory output".format(
                    coverage.returncode))
            transcripts.append("===== {} =====\n{}".format(
                coverage_id, coverage.stdout or ""))

        missing = [identifier for identifier in requested if identifier not in results]
        if missing:
            raise RunnerError("matrix failed to publish selected identities: {}".format(missing))
        ordered = [results[identifier] for identifier in requested]
        status = max((item["status"] for item in ordered),
                     key=lambda value: SEVERITY[value])
        transcript = artifacts / "t12-portability-matrix.txt"
        transcript.write_text("\n".join(transcripts), encoding="utf-8")
        return write_result(
            context, status,
            "bounded portability matrix executed; strict coverage and typed skips remain authoritative",
            cases=ordered, evidence=[str(transcript)],
            metrics={"matrix_cases": len(ordered),
                     "qualification": "bounded-T12-evidence-not-coverage-closure"})
    except BaseException as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
