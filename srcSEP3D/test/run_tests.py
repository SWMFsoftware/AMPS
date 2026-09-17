#!/usr/bin/env python3
"""Unified srcSEP3D test runner.

The command vocabulary intentionally follows ``srcSEP/test/run_tests.py``:
``--list``, repeatable ``--test``/``--group``/``--suite``, ``--routine``,
``--all``, and ``--output-dir`` have the same meaning.  This R0 runner joins
three evidence classes without pretending they are interchangeable:

* standalone C++ tests compile L0/L1 without AMPS or MPI;
* R0 source/ABI checks inspect the actual production manifest and AMPS pic.h;
* the strict production build runs only with a real configured AMPS checkout.

Missing production configuration is reported as SKIP, never as PASS.  That
distinction keeps a source-only development run useful while ensuring it
cannot satisfy the AMPS production-build acceptance gate.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as _xml_escape


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "test" / "stage1"


@dataclass(frozen=True)
class TestDefinition:
    test_id: str
    group: str
    name: str
    kind: str
    routine: bool = True


@dataclass
class Result:
    test_id: str
    group: str
    status: str
    message: str
    elapsed_seconds: float
    command: List[str]


# The C++ registry remains authoritative for callbacks and result semantics.
# This compact manifest exists only so --list and unknown-ID validation work
# before a binary is built, matching srcSEP's discovery-first CLI behavior.
TESTS: Tuple[TestDefinition, ...] = (
    TestDefinition("BLD01", "BLD", "Standalone binary excludes AMPS/MPI symbols", "cpp"),
    TestDefinition("HARN01", "HARN", "Empty-registry contract", "cpp"),
    TestDefinition("HARN02", "HARN", "Failure/error exit-code contract", "cpp"),
    TestDefinition("HARN03", "HARN", "Skip exit-code contract", "cpp"),
    TestDefinition("HARN04", "HARN", "JSON/JUnit writer contract", "cpp"),
    TestDefinition("LAY01", "LAY", "L0/L1 AMPS dependency exclusion", "cpp"),
    TestDefinition("LAY02", "LAY", "Layering negative control", "cpp"),
    TestDefinition("UTIL02", "UTIL", "Shared-kernel frozen record", "cpp"),
    TestDefinition("HARN02-EXITCODE", "HARN_SHELL", "Outer failure exit code", "shell", False),
    TestDefinition("HARN03-EXITCODE", "HARN_SHELL", "Outer skip exit code", "shell", False),
    TestDefinition("RUN3D01", "RUNNER", "Python runner CLI contract", "source"),
    TestDefinition("BLDL3D01", "BLDL3D", "Configured enclosing AMPS build", "source"),
    TestDefinition("BLDL3D02", "BLDL3D", "Retired production-symbol exclusion", "source"),
    TestDefinition("BLDL3D03", "BLDL3D", "AMPS mover return-code mapping", "source"),
    TestDefinition("BLDL3D04", "BLDL3D", "AMPS Pi-macro namespace hygiene", "source"),
    TestDefinition("BLDL3D05", "BLDL3D", "Source/build makefile path resolution", "source"),
)

BY_ID: Dict[str, TestDefinition] = {item.test_id: item for item in TESTS}
GROUPS: Dict[str, List[TestDefinition]] = {}
for _item in TESTS:
    GROUPS.setdefault(_item.group, []).append(_item)

SUITES: Dict[str, Tuple[str, ...]] = {
    "standalone": tuple(item.test_id for item in TESTS
                        if item.kind in ("cpp", "shell") or item.test_id == "RUN3D01"),
    "r0": ("BLDL3D01", "BLDL3D02", "BLDL3D03", "BLDL3D04", "BLDL3D05",
           "RUN3D01", "LAY01", "BLD01"),
    "production": ("BLDL3D01", "BLDL3D02", "BLDL3D03", "BLDL3D04",
                   "BLDL3D05"),
}


class RunnerError(RuntimeError):
    """Usage/setup error that must exit 2 rather than look like physics FAIL."""


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass


EPILOG = """
Examples:
  python3 test/run_tests.py --list
  python3 test/run_tests.py --routine --amps-source /path/to/AMPS
  python3 test/run_tests.py --test BLDL3D02 --output-dir test_output/r0-source
  python3 test/run_tests.py --group HARN --group BLDL3D
  python3 test/run_tests.py --suite standalone --suite production
  python3 test/run_tests.py --all --amps-source /path/to/AMPS

For a configured AMPS checkout, either place srcSEP3D in its normal
application location or provide --make-config /path/to/Makefile.conf.  A
source-only archive can run all standalone tests and BLDL3D02/03; BLDL3D01
will be recorded as SKIP until the real production configuration is present.
"""


def _utc_stamp() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run srcSEP3D standalone, R0, and production-build tests.",
        formatter_class=_HelpFormatter,
        epilog=EPILOG)
    parser.add_argument("--amps", default=os.environ.get("SEP3D_EXECUTABLE"),
                        help="linked srcSEP3D/AMPS executable (reserved for linked phases)")
    parser.add_argument("--amps-source", type=Path,
                        default=os.environ.get("AMPS_SOURCE_ROOT"),
                        help="AMPS root containing src/pic/pic.h for ABI checks")
    parser.add_argument("--make-config", type=Path,
                        default=os.environ.get("AMPS_MAKE_CONFIG"),
                        help="configured AMPS Makefile.conf for BLDL3D01")
    parser.add_argument("--list", action="store_true",
                        help="list registered tests and suites, then exit")
    parser.add_argument("--test", dest="tests", action="append", default=[],
                        metavar="ID", help="run one stable test ID; repeatable")
    parser.add_argument("--group", dest="groups", action="append", default=[],
                        metavar="GROUP", help="run one test group; repeatable")
    parser.add_argument("--routine", action="store_true",
                        help="run the bounded routine set")
    parser.add_argument("--all", action="store_true",
                        help="run every public test separately and continue after failures")
    parser.add_argument("--suite", dest="suites", action="append", default=[],
                        choices=sorted(SUITES), help="run a named suite; repeatable")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("test_output") / f"srcsep3d-tests-{_utc_stamp()}",
                        help="directory for logs and JSON/JUnit summaries")
    parser.add_argument("--sep-common-dir", type=Path,
                        help="canonical sep_common source/header directory")
    parser.add_argument("--sep-common-archive", type=Path,
                        help="path to sep_common.a")
    parser.add_argument("--cxx", default="g++", help="standalone C++ compiler")
    parser.add_argument("--rebuild", action="store_true",
                        help="force rebuilding test/stage1")
    parser.add_argument("--no-build", action="store_true",
                        help="use an existing test/stage1 binary")
    parser.add_argument("--timeout", type=float, default=300.0,
                        help="per-test command timeout in seconds")
    parser.add_argument("--verbose", action="store_true",
                        help="print complete commands and subprocess output")
    parser.add_argument("--no-color", action="store_true",
                        help="accepted for CLI parity; output is already plain text")
    return parser


def _resolve_existing(explicit: Optional[Path], candidates: Iterable[Path],
                      description: str) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise RunnerError(f"{description} does not exist: {path}")
        return path
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    rendered = "\n  ".join(str(item) for item in candidates)
    raise RunnerError(f"cannot locate {description}; checked:\n  {rendered}")


def _sep_common_paths(args: argparse.Namespace) -> Tuple[Path, Path]:
    common_dir = _resolve_existing(
        args.sep_common_dir,
        (ROOT.parent / "src" / "models" / "sep_common",
         ROOT.parent / "sep_common",
         ROOT.parent.parent / "sep_common"),
        "sep_common directory")

    archive = (args.sep_common_archive.expanduser().resolve()
               if args.sep_common_archive is not None
               else common_dir / "sep_common.a")
    if not archive.is_file():
        make = shutil.which("make")
        if make is None:
            raise RunnerError(
                f"sep_common archive is absent and make is unavailable: {archive}")
        completed = subprocess.run(
            [make, "-C", str(common_dir), "all"], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        if completed.returncode != 0:
            raise RunnerError("sep_common build failed:\n" + completed.stdout)
    if not archive.is_file():
        raise RunnerError(f"sep_common build produced no archive: {archive}")
    return common_dir, archive


def _build_standalone(args: argparse.Namespace) -> None:
    if args.no_build:
        if not BINARY.is_file():
            raise RunnerError("--no-build was requested but test/stage1 is absent")
        return

    source, archive = _sep_common_paths(args)
    inputs = [
        ROOT / "test" / "stage1.cpp",
        ROOT / "test" / "individual-test" / "test_harness.cpp",
        ROOT / "test" / "individual-test" / "test_layering.cpp",
        ROOT / "test" / "individual-test" / "test_build.cpp",
        ROOT / "test" / "individual-test" / "test_kernels.cpp",
        archive,
    ]
    newest_input = max(path.stat().st_mtime for path in inputs)
    if (BINARY.is_file() and not args.rebuild
            and BINARY.stat().st_mtime >= newest_input):
        return

    compiler = shutil.which(args.cxx)
    if compiler is None:
        raise RunnerError(f"C++ compiler not found: {args.cxx}")
    BINARY.parent.mkdir(parents=True, exist_ok=True)
    command = [
        compiler, "-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
        "-O2", f"-I{ROOT / 'core'}", f"-I{ROOT / 'background'}", f"-I{source}",
        *(str(path) for path in inputs[:-1]), str(archive), "-o", str(BINARY),
    ]
    if args.verbose:
        print("RUN:", shlex.join(command))
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               timeout=args.timeout, check=False)
    if completed.returncode != 0:
        raise RunnerError("standalone build failed:\n" + completed.stdout)
    if args.verbose and completed.stdout:
        print(completed.stdout, end="")


def _run_command(command: List[str], cwd: Path, timeout: float,
                 verbose: bool) -> Tuple[int, str, float]:
    if verbose:
        print("RUN:", shlex.join(command))
    started = time.monotonic()
    try:
        completed = subprocess.run(command, cwd=cwd, text=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   timeout=timeout, check=False)
        elapsed = time.monotonic() - started
        if verbose and completed.stdout:
            print(completed.stdout, end="")
        return completed.returncode, completed.stdout, elapsed
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + (exc.stderr or "")
        return 2, output + f"\nTimed out after {timeout:g} s", time.monotonic() - started


def _parse_cpp_result(report: Path, definition: TestDefinition,
                      command: List[str], elapsed: float, output: str) -> Result:
    if not report.is_file():
        return Result(definition.test_id, definition.group, "ERROR",
                      "C++ test wrote no JSON report. Output: " + output[-2000:],
                      elapsed, command)
    try:
        payload = json.loads(report.read_text(encoding="utf-8"))
        records = payload.get("results", [])
        record = next(item for item in records
                      if str(item.get("id", "")).upper() == definition.test_id)
    except (OSError, ValueError, StopIteration) as exc:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"cannot parse C++ report: {exc}", elapsed, command)
    return Result(definition.test_id, definition.group,
                  str(record.get("status", "ERROR")).upper(),
                  str(record.get("message", "")), elapsed, command)


def _run_cpp(definition: TestDefinition, args: argparse.Namespace,
             output_dir: Path) -> Result:
    report = output_dir / f"{definition.test_id}.json"
    command = [str(BINARY), "--test", definition.test_id,
               "--test-json", str(report)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    result = _parse_cpp_result(report, definition, command, elapsed, output)
    if code == 2 and result.status not in ("ERROR",):
        result.status = "ERROR"
        result.message += " (binary exit code 2)"
    return result


def _run_shell(definition: TestDefinition, args: argparse.Namespace) -> Result:
    script_name = {
        "HARN02-EXITCODE": "run_harn02.sh",
        "HARN03-EXITCODE": "run_harn03.sh",
    }[definition.test_id]
    command = ["sh", str(ROOT / "test" / "individual-test" / script_name)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    return Result(definition.test_id, definition.group,
                  "PASS" if code == 0 else "FAIL",
                  output.strip() or f"script exited {code}", elapsed, command)


def _production_files() -> List[Path]:
    files = [ROOT / "SEP3D.h", ROOT / "main_lib.cpp", ROOT / "main.cpp"]
    for directory in (ROOT / "core", ROOT / "background", ROOT / "amps"):
        files.extend(sorted(directory.glob("*.h")))
        files.extend(sorted(directory.glob("*.cpp")))
    return files


def _strip_cpp_comments(text: str) -> str:
    """Remove comments before checking for retired *identifiers*.

    Migration comments are allowed to explain a removed name.  The production
    gate is concerned with compiled tokens, so both // and /*...*/ comments
    are removed while quoted strings are retained (a runtime lookup by a
    retired name must still fail the check).
    """
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*?$", "", text, flags=re.MULTILINE)


def _check_retired_sources(definition: TestDefinition) -> Result:
    started = time.monotonic()
    errors: List[str] = []
    retired_file = ROOT / "SEP3D.cpp"
    if retired_file.exists():
        errors.append(
            "retired SEP3D.cpp is still present; remove this exact stale file "
            "from the installed srcSEP3D tree before rebuilding")

    makefile = (ROOT / "makefile").read_text(encoding="utf-8")
    # Strip comments before checking the active production manifest.
    active_make = "\n".join(line.split("#", 1)[0] for line in makefile.splitlines())
    if "SEP3D.cpp" in active_make or "SEP3D.o" in active_make:
        errors.append("makefile production manifest still references SEP3D.cpp")

    forbidden = (
        "Mover_Axisymmetric_SecondOrder", "TotalParticleAcceleration",
        "GlobalEnergyDistribution", "inject_particle_onto_field_line",
        "CMPI_channel", "8.760e+08", "9.445e+08",
    )
    for path in _production_files():
        text = _strip_cpp_comments(path.read_text(encoding="utf-8"))
        for token in forbidden:
            if token in text:
                errors.append(f"{path.relative_to(ROOT)} contains retired token {token}")

    main_lib = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    for call in ("PrepopulateDomain", "outputMeshDataTECPLOT", "saveMeshFile"):
        if call in main_lib:
            errors.append(f"main_lib.cpp retains prototype operation {call}")

    status = "FAIL" if errors else "PASS"
    message = "; ".join(errors) if errors else (
        "retired mover/sampler source is absent; production manifest contains only "
        "main_lib.cpp and main.cpp; wedge and prepopulation operations are absent")
    return Result(definition.test_id, definition.group, status, message,
                  time.monotonic() - started, [])


def _find_pic_header(args: argparse.Namespace) -> Optional[Path]:
    roots: List[Path] = []
    if args.amps_source is not None:
        roots.append(args.amps_source.expanduser().resolve())
    roots.extend((ROOT.parent, ROOT.parent.parent))
    for root in roots:
        for candidate in (root / "src" / "pic" / "pic.h",
                          root / "pic" / "pic.h", root / "pic.h"):
            if candidate.is_file():
                return candidate
    return None


def _macro_int(text: str, name: str) -> Optional[int]:
    match = re.search(rf"^\s*#\s*define\s+{re.escape(name)}\s+(-?\d+)\b",
                      text, re.MULTILINE)
    return int(match.group(1)) if match else None


def _check_return_codes(definition: TestDefinition,
                        args: argparse.Namespace) -> Result:
    started = time.monotonic()
    adapter = ROOT / "amps" / "amps_mover_status.h"
    adapter_text = adapter.read_text(encoding="utf-8")
    required_adapter_tokens = (
        "static_assert(_PARTICLE_DELETED_ON_THE_FACE_ == 0",
        "static_assert(_PARTICLE_LEFT_THE_DOMAIN_ == 2",
        "static_assert(_PARTICLE_MOTION_FINISHED_ == 3",
        "Core::ParticleMotionOutcome::Advanced",
        "return _PARTICLE_MOTION_FINISHED_",
        "Core::ParticleMotionOutcome::LeftDomain",
        "return _PARTICLE_LEFT_THE_DOMAIN_",
    )
    missing = [token for token in required_adapter_tokens if token not in adapter_text]
    if missing:
        return Result(definition.test_id, definition.group, "FAIL",
                      "adapter contract is incomplete: " + ", ".join(missing),
                      time.monotonic() - started, [])

    pic_header = _find_pic_header(args)
    if pic_header is None:
        return Result(definition.test_id, definition.group, "SKIP",
                      "AMPS pic.h not found; provide --amps-source to verify ABI values",
                      time.monotonic() - started, [])
    pic_text = pic_header.read_text(encoding="utf-8", errors="replace")
    expected = {
        "_PARTICLE_DELETED_ON_THE_FACE_": 0,
        "_PARTICLE_LEFT_THE_DOMAIN_": 2,
        "_PARTICLE_MOTION_FINISHED_": 3,
    }
    observed = {name: _macro_int(pic_text, name) for name in expected}
    mismatches = [f"{name}: expected {value}, found {observed[name]}"
                  for name, value in expected.items() if observed[name] != value]
    if mismatches:
        return Result(definition.test_id, definition.group, "FAIL",
                      "; ".join(mismatches), time.monotonic() - started, [])
    return Result(definition.test_id, definition.group, "PASS",
                  f"adapter mappings match {pic_header}",
                  time.monotonic() - started, [])


def _check_macro_hygiene(definition: TestDefinition,
                         args: argparse.Namespace) -> Result:
    """Compile the application constants after a synthetic AMPS Pi macro."""
    compiler = shutil.which(args.cxx)
    if compiler is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"C++ compiler not found: {args.cxx}", 0.0, [])
    source = ROOT / "test" / "compile_macro_hygiene.cpp"
    command = [compiler, "-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
               "-Werror", "-fsyntax-only", str(source)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    if code != 0:
        return Result(definition.test_id, definition.group, "FAIL",
                      "AMPS Pi macro still corrupts srcSEP3D headers:\n" + output[-3000:],
                      elapsed, command)
    return Result(definition.test_id, definition.group, "PASS",
                  "sep3d_types.h compiles after the AMPS Pi macro is defined",
                  elapsed, command)


def _check_makefile_relocation(definition: TestDefinition,
                               args: argparse.Namespace,
                               output_dir: Path) -> Result:
    """Verify path discovery from both source and AMPS/build/main locations."""
    make = shutil.which("make")
    if make is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      "make is not available", 0.0, [])

    # The temporary tree reproduces only the paths involved in discovery.  No
    # production compilation is claimed here; BLDL3D01 remains authoritative
    # for compilation against a real configured Makefile.conf.
    fixture = output_dir / "BLDL3D05-layout"
    if fixture.exists():
        shutil.rmtree(fixture)
    source_dir = fixture / "srcSEP3D"
    build_dir = fixture / "build" / "main"
    common_dir = fixture / "src" / "models" / "sep_common"
    source_dir.mkdir(parents=True)
    build_dir.mkdir(parents=True)
    common_dir.mkdir(parents=True)
    (fixture / "Makefile.conf").write_text(
        "# Empty BLDL3D05 configuration fixture.\n", encoding="utf-8")
    # The fixture's enclosing target creates valid empty archives.  It proves
    # that strict-production delegates outward and audits build/main rather
    # than attempting the invalid direct compile that lacks build/pic headers.
    (fixture / "Makefile").write_text(
        "amps:\n"
        "\t@mkdir -p build/main\n"
        "\t@ar -rc build/main/mainlib.a\n"
        "\t@ar -rc build/main/main.a\n",
        encoding="utf-8")
    shutil.copy2(ROOT / "makefile", source_dir / "makefile")
    shutil.copy2(ROOT / "makefile", build_dir / "makefile")

    expected = {
        f"AMPS_ROOT={fixture.resolve()}",
        f"AMPS_CONFIG={(fixture / 'Makefile.conf').resolve()}",
        f"SEP_COMMON_DIR={common_dir.resolve()}",
    }
    commands: List[str] = []
    elapsed = 0.0
    for location in (source_dir, build_dir):
        # Deliberately run from the fixture root, not from either makefile
        # directory.  A relative-to-CWD implementation would fail this probe.
        command = [make, "--no-print-directory", "-f",
                   str(location / "makefile"), "print-layout-paths"]
        code, output, duration = _run_command(
            command, fixture, args.timeout, args.verbose)
        commands.extend([f"cwd={fixture}", *command])
        elapsed += duration
        observed = set(output.splitlines())
        expected_for_location = expected | {
            f"SEP3D_MAKEFILE_DIR={location.resolve()}"
        }
        if code != 0:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"makefile path probe failed from {location}:\n" + output[-3000:],
                elapsed, commands)
        missing = sorted(expected_for_location - observed)
        if missing:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"makefile resolved incorrect paths from {location}; missing: " +
                ", ".join(missing), elapsed, commands)

    production_command = [
        make, "--no-print-directory", "-f", str(source_dir / "makefile"),
        "strict-production", f"AMPS_ROOT={fixture.resolve()}",
        f"AMPS_CONFIG={(fixture / 'Makefile.conf').resolve()}",
    ]
    code, output, duration = _run_command(
        production_command, output_dir, args.timeout, args.verbose)
    commands.extend([f"cwd={output_dir}", *production_command])
    elapsed += duration
    if code != 0:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "strict-production did not delegate to the enclosing AMPS target "
            "and audit build/main archives:\n" + output[-3000:],
            elapsed, commands)

    return Result(
        definition.test_id, definition.group, "PASS",
        "source srcSEP3D and copied build/main makefiles resolve the same "
        "AMPS root, Makefile.conf, and sep_common directory; production "
        "orchestration delegates to enclosing make amps",
        elapsed, commands)


def _find_make_config(args: argparse.Namespace) -> Optional[Path]:
    if args.make_config is not None:
        path = args.make_config.expanduser().resolve()
        return path if path.is_file() else None
    for candidate in (ROOT / ".." / ".." / "Makefile.conf",
                      ROOT.parent / "Makefile.conf",
                      ROOT.parent.parent / "Makefile.conf"):
        resolved = candidate.resolve()
        if resolved.is_file():
            return resolved
    return None


def _find_amps_root(args: argparse.Namespace,
                    config: Path) -> Optional[Path]:
    """Find the enclosing checkout that owns the production `make amps` target."""
    candidates: List[Path] = []
    if args.amps_source is not None:
        candidates.append(args.amps_source.expanduser().resolve())
    candidates.extend((config.parent.resolve(), ROOT.parent.resolve()))
    for candidate in candidates:
        if (candidate / "Makefile").is_file():
            return candidate
    return None


def _check_production_build(definition: TestDefinition,
                            args: argparse.Namespace) -> Result:
    config = _find_make_config(args)
    if config is None:
        return Result(definition.test_id, definition.group, "SKIP",
                      "configured AMPS Makefile.conf not found; provide --make-config",
                      0.0, [])
    amps_root = _find_amps_root(args, config)
    if amps_root is None:
        return Result(
            definition.test_id, definition.group, "ERROR",
            "enclosing AMPS Makefile not found; provide --amps-source pointing "
            "to the AMPS root", 0.0, [])
    command = ["make", "-f", str(ROOT / "makefile"), "strict-production",
               f"AMPS_ROOT={amps_root}", f"AMPS_CONFIG={config}"]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    if code == 0:
        return Result(definition.test_id, definition.group, "PASS",
                      "enclosing `make amps` completed and build/main archives "
                      "passed the retired-symbol audit",
                      elapsed, command)
    return Result(definition.test_id, definition.group, "FAIL",
                  "configured enclosing AMPS build failed:\n" + output[-4000:],
                  elapsed, command)


def _run_source(definition: TestDefinition, args: argparse.Namespace,
                output_dir: Path) -> Result:
    if definition.test_id == "BLDL3D01":
        return _check_production_build(definition, args)
    if definition.test_id == "BLDL3D02":
        return _check_retired_sources(definition)
    if definition.test_id == "BLDL3D03":
        return _check_return_codes(definition, args)
    if definition.test_id == "BLDL3D04":
        return _check_macro_hygiene(definition, args)
    if definition.test_id == "BLDL3D05":
        return _check_makefile_relocation(definition, args, output_dir)
    if definition.test_id == "RUN3D01":
        command = [sys.executable, str(ROOT / "test" / "test_python_runner.py")]
        code, output, elapsed = _run_command(
            command, ROOT, args.timeout, args.verbose)
        return Result(definition.test_id, definition.group,
                      "PASS" if code == 0 else "FAIL",
                      output.strip() or f"runner unit test exited {code}",
                      elapsed, command)
    raise RunnerError(f"no source-test implementation for {definition.test_id}")


def _select(args: argparse.Namespace) -> List[TestDefinition]:
    modes = sum((bool(args.list), bool(args.tests or args.groups),
                 bool(args.routine), bool(args.all), bool(args.suites)))
    if modes != 1:
        raise RunnerError("choose exactly one mode: --list, --test/--group, "
                          "--routine, --all, or --suite")
    if args.all and (args.tests or args.groups or args.routine or args.suites):
        raise RunnerError("--all cannot be combined with another selector")

    selected: List[TestDefinition] = []
    if args.routine:
        selected = [item for item in TESTS if item.routine]
    elif args.all:
        selected = list(TESTS)
    elif args.suites:
        for suite in args.suites:
            selected.extend(BY_ID[test_id] for test_id in SUITES[suite])
    else:
        for test_id in args.tests:
            key = test_id.upper()
            if key not in BY_ID:
                raise RunnerError(f"unknown test ID '{test_id}'; use --list")
            selected.append(BY_ID[key])
        for group in args.groups:
            key = group.upper()
            if key not in GROUPS:
                raise RunnerError(f"unknown group '{group}'; use --list")
            selected.extend(GROUPS[key])

    # Stable-ID de-duplication mirrors the srcSEP C++ registry.
    unique: Dict[str, TestDefinition] = {}
    for item in selected:
        unique.setdefault(item.test_id, item)
    return [unique[key] for key in sorted(unique)]


def _print_list() -> None:
    print(f"{'ID':<20} {'GROUP':<14} {'KIND':<10} NAME")
    print("-" * 88)
    for item in TESTS:
        print(f"{item.test_id:<20} {item.group:<14} {item.kind:<10} {item.name}")
    print("\nSuites:")
    for name, ids in sorted(SUITES.items()):
        print(f"  {name:<12} {' '.join(ids)}")


def _write_reports(results: Sequence[Result], output_dir: Path) -> None:
    counts = {status: sum(item.status == status for item in results)
              for status in ("PASS", "FAIL", "SKIP", "ERROR")}
    exit_code = 2 if counts["ERROR"] else (1 if counts["FAIL"] else 0)
    payload = {
        "schema": "srcsep3d-tests-v1",
        "generated_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "root": str(ROOT),
        "exit_code": exit_code,
        "totals": {key.lower(): value for key, value in counts.items()},
        "results": [asdict(item) for item in results],
    }
    (output_dir / "srcsep3d-tests.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    cases = []
    for item in results:
        body = ""
        if item.status == "FAIL":
            body = f'<failure message="{_xml_escape(item.message)}"/>'
        elif item.status == "ERROR":
            body = f'<error message="{_xml_escape(item.message)}"/>'
        elif item.status == "SKIP":
            body = f'<skipped message="{_xml_escape(item.message)}"/>'
        cases.append(
            f'  <testcase classname="srcSEP3D.{item.group}" '
            f'name="{item.test_id}" time="{item.elapsed_seconds:.9f}">{body}</testcase>')
    xml = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (f'<testsuite name="srcSEP3D" tests="{len(results)}" '
         f'failures="{counts["FAIL"]}" errors="{counts["ERROR"]}" '
         f'skipped="{counts["SKIP"]}">'),
        *cases, "</testsuite>", "",
    ]
    (output_dir / "srcsep3d-tests.xml").write_text("\n".join(xml), encoding="utf-8")


def _print_summary(results: Sequence[Result]) -> int:
    for item in results:
        print(f"[{item.test_id}] {item.status} ({item.elapsed_seconds:.3f}s) {item.message}")
    counts = {status: sum(item.status == status for item in results)
              for status in ("PASS", "FAIL", "SKIP", "ERROR")}
    print("\nSummary: " + " ".join(f"{key}={value}" for key, value in counts.items()))
    for status in ("FAIL", "ERROR"):
        ids = [item.test_id for item in results if item.status == status]
        print(f"{status}: " + (", ".join(ids) if ids else "none"))
    return 2 if counts["ERROR"] else (1 if counts["FAIL"] else 0)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.timeout <= 0.0:
        raise RunnerError("--timeout must be positive")
    if args.list:
        # Selection validation is still performed so contradictory modes fail
        # consistently instead of silently preferring --list.
        _select(args)
        _print_list()
        return 0

    selected = _select(args)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if any(item.kind in ("cpp", "shell") for item in selected):
        _build_standalone(args)

    results: List[Result] = []
    for definition in selected:
        if definition.kind == "cpp":
            result = _run_cpp(definition, args, output_dir)
        elif definition.kind == "shell":
            result = _run_shell(definition, args)
        else:
            result = _run_source(definition, args, output_dir)
        results.append(result)
        print(f"[{result.test_id}] {result.status}", flush=True)

    _write_reports(results, output_dir)
    return _print_summary(results)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RunnerError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
