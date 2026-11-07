#!/usr/bin/env python3
"""Build and execute the one-rank T3 AMR/geometry profile matrix.

The three dimensional profiles are compile-time alternatives, so this adapter
performs one clean repository build for every selected dimension, retains each
executable and generated configuration, and only then executes the requested
case sequence.  Retaining all executables before execution preserves the
global runner's normal/reverse/shuffle case order without rebuilding for every
case.

Every profile build follows the repository's native-build preflight exactly:
confirm the checkout and absence of another AMPS build/run, remove only the
generated ``build`` directory, select the application with ``Config.pl``, and
compile through the top-level workflow with ``MAKEFLAGS=-j16``.  This adapter
never edits ``.amps.conf`` or another Git control/configuration file directly.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


# Keep every ID literal so the independent traceability validator can prove
# registry ownership without evaluating Python string construction.
PROFILE_1D_CASES = [
    "CORE-MESH-1D-TOPOLOGY-UNIFORM-001",
    "CORE-MESH-1D-TOPOLOGY-NONUNIFORM-001",
    "CORE-MESH-1D-SEARCH-001",
    "CORE-MESH-1D-PARENT-CHILD-001",
    "CORE-MESH-1D-NEIGHBORS-001",
    "CORE-MESH-1D-CHECKSUM-001",
    "CORE-MESH-1D-ADAPT-REFINE-001",
    "CORE-MESH-1D-ADAPT-COARSEN-001",
    "CORE-MESH-1D-ADAPT-DATA-001",
    "CORE-MESH-1D-ADAPT-VOLUME-001",
    "CORE-MESH-1D-ADAPT-COUNTERS-001",
    "CORE-MESH-1D-ADAPT-LOAD-001",
    "CORE-MESH-1D-ADAPT-MAX-LEVEL-001",
    "CORE-MESH-1D-ADAPT-CYCLE-001",
]

PROFILE_2D_CASES = [
    "CORE-MESH-2D-TOPOLOGY-UNIFORM-001",
    "CORE-MESH-2D-TOPOLOGY-NONUNIFORM-001",
    "CORE-MESH-2D-SEARCH-001",
    "CORE-MESH-2D-PARENT-CHILD-001",
    "CORE-MESH-2D-NEIGHBORS-001",
    "CORE-MESH-2D-CHECKSUM-001",
    "CORE-MESH-2D-ADAPT-REFINE-001",
    "CORE-MESH-2D-ADAPT-COARSEN-001",
    "CORE-MESH-2D-ADAPT-DATA-001",
    "CORE-MESH-2D-ADAPT-VOLUME-001",
    "CORE-MESH-2D-ADAPT-COUNTERS-001",
    "CORE-MESH-2D-ADAPT-LOAD-001",
    "CORE-MESH-2D-ADAPT-MAX-LEVEL-001",
    "CORE-MESH-2D-ADAPT-CYCLE-001",
]

PROFILE_3D_MESH_CASES = [
    "CORE-MESH-3D-TOPOLOGY-UNIFORM-001",
    "CORE-MESH-3D-TOPOLOGY-NONUNIFORM-001",
    "CORE-MESH-3D-SEARCH-001",
    "CORE-MESH-3D-PARENT-CHILD-001",
    "CORE-MESH-3D-NEIGHBORS-001",
    "CORE-MESH-3D-CHECKSUM-001",
    "CORE-MESH-3D-ADAPT-REFINE-001",
    "CORE-MESH-3D-ADAPT-COARSEN-001",
    "CORE-MESH-3D-ADAPT-DATA-001",
    "CORE-MESH-3D-ADAPT-VOLUME-001",
    "CORE-MESH-3D-ADAPT-COUNTERS-001",
    "CORE-MESH-3D-ADAPT-LOAD-001",
    "CORE-MESH-3D-ADAPT-MAX-LEVEL-001",
    "CORE-MESH-3D-ADAPT-CYCLE-001",
]

INTERPOLATION_CASES = [
    "CORE-MESH-3D-INTERPOLATION-UNIFORM-001",
    "CORE-MESH-3D-INTERPOLATION-BOUNDARY-001",
    "CORE-MESH-3D-INTERPOLATION-AMR-001",
    "CORE-MESH-3D-INTERPOLATION-LEGACY-MASK-001",
    "CORE-MESH-3D-INTERPOLATION-DUPLICATE-001",
    "CORE-MESH-3D-INTERPOLATION-PERMUTATION-001",
    "CORE-MESH-3D-INTERPOLATION-CACHE-001",
    "CORE-MESH-3D-INTERPOLATION-INPUT-IMMUTABLE-001",
]

CUTCELL_CASES = [
    "MESH-CUTCELL-TRIANGLE-METRICS-001",
    "MESH-CUTCELL-RAY-001",
    "MESH-CUTCELL-TANGENT-DEGENERATE-001",
    "MESH-CUTCELL-CLOSEST-001",
    "MESH-CUTCELL-PLANE-VOLUME-001",
    "MESH-CUTCELL-SPHERE-VOLUME-001",
    "MESH-CUTCELL-SPHERE-AREA-CONVERGENCE-001",
    "MESH-CUTCELL-NASTRAN-ROUNDTRIP-001",
]

PROFILE_3D_CASES = PROFILE_3D_MESH_CASES + INTERPOLATION_CASES + CUTCELL_CASES
ALL_CASES = PROFILE_1D_CASES + PROFILE_2D_CASES + PROFILE_3D_CASES

PROFILES: Dict[str, Dict[str, Any]] = {
    "1d": {
        "application": "test/core-mesh-1d",
        "input_deck": "input/test/core-mesh-1d.input",
        "dimension": 1,
        "cases": PROFILE_1D_CASES,
    },
    "2d": {
        "application": "test/core-mesh-2d",
        "input_deck": "input/test/core-mesh-2d.input",
        "dimension": 2,
        "cases": PROFILE_2D_CASES,
    },
    "3d": {
        "application": "test/core-mesh-3d",
        "input_deck": "input/test/core-mesh-3d.input",
        "dimension": 3,
        "cases": PROFILE_3D_CASES,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(path))


def effective_macros(path: Path) -> Dict[str, str]:
    """Replay object-like define/undef directives in source order."""

    result: Dict[str, str] = {}
    text = path.read_text(encoding="utf-8", errors="replace")
    for raw in text.splitlines():
        undef = re.match(r"\s*#\s*undef\s+([A-Za-z_][A-Za-z0-9_]*)", raw)
        if undef:
            result.pop(undef.group(1), None)
            continue
        define = re.match(
            r"\s*#\s*define\s+([A-Za-z_][A-Za-z0-9_]*)\s+([^/\s]+)", raw)
        if define:
            result[define.group(1)] = define.group(2)
    return result


def validate_generated_configuration(
        repo: Path, profile: Mapping[str, Any]) -> Tuple[bool, str]:
    required = [
        repo / ".amps.conf",
        repo / "amps",
        repo / "build" / "pic" / "picGlobal.dfn",
        repo / "build" / "pic" / "pic.h",
        repo / "build" / "meshAMR" / "meshAMRdef.h",
    ]
    missing = [path.relative_to(repo).as_posix() for path in required
               if not path.is_file()]
    if missing:
        return False, "missing generated artifacts: {}".format(
            ", ".join(missing))

    selection = (repo / ".amps.conf").read_text(
        encoding="utf-8", errors="replace")
    macros = effective_macros(repo / "build" / "pic" / "picGlobal.dfn")
    application = str(profile["application"]).rsplit("/", 1)[-1]
    input_name = Path(str(profile["input_deck"])).name
    dimension = str(profile["dimension"])
    checks = {
        "APPLICATION={}".format(application):
            "APPLICATION={}".format(application) in selection,
        "InputFileAMPS={}".format(input_name):
            "InputFileAMPS={}".format(input_name) in selection,
        "SPICE=nospice": "SPICE=nospice" in selection,
        "DIM={}".format(dimension): macros.get("DIM") == dimension,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        return False, "generated profile checks failed: {}".format(
            ", ".join(failed))
    return True, "application, input, SPICE, and DIM={} verified".format(
        dimension)


def active_repo_processes(repo: Path) -> List[str]:
    """Find native AMPS/build processes rooted in this checkout."""

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active: List[str] = []
    own_pid = os.getpid()
    names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl"}
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) == own_pid:
            continue
        try:
            cwd = (entry / "cwd").resolve()
            cwd.relative_to(repo)
            raw = (entry / "cmdline").read_bytes().replace(
                b"\x00", b" ").decode("utf-8", errors="replace").strip()
        except (OSError, ValueError):
            continue
        executable = Path(raw.split()[0]).name if raw else ""
        if (executable in names or executable.startswith("amps") or
                "Config.pl" in raw or "ampsConfig.pl" in raw):
            active.append("pid={} cwd={} command={}".format(
                entry.name, cwd, raw))
    return active


def retain_profile(
        repo: Path, artifacts: Path, profile_name: str,
        profile: Mapping[str, Any]) -> Tuple[Path, Path]:
    retained = artifacts / "native-build" / profile_name
    retained.mkdir(parents=True, exist_ok=True)
    executable = retained / "amps-core-mesh-{}".format(profile_name)
    candidates = [
        repo / "amps",
        repo / ".amps.conf",
        repo / "build" / "Makefile.conf",
        repo / "build" / "pic" / "picGlobal.dfn",
        repo / "build" / "pic" / "picSpeciesMacro.dfn",
        repo / "build" / "pic" / "pic.h",
        repo / "build" / "meshAMR" / "meshAMRdef.h",
        repo / str(profile["input_deck"]),
    ]
    files: List[Dict[str, Any]] = []
    for source in candidates:
        if not source.is_file():
            continue
        target = executable if source == repo / "amps" else (
            retained / source.relative_to(repo).as_posix().replace("/", "__"))
        shutil.copy2(str(source), str(target))
        files.append({
            "repository_path": source.relative_to(repo).as_posix(),
            "retained_path": str(target),
            "size_bytes": target.stat().st_size,
            "sha256": sha256(target),
        })
    identity = retained / "native-build-identity.json"
    write_json(identity, {
        "schema_version": 1,
        "profile": profile_name,
        "dimension": int(profile["dimension"]),
        "application": str(profile["application"]),
        "input_deck": str(profile["input_deck"]),
        "configure_commands": [
            ["./Config.pl", "-application={}".format(profile["application"])],
            ["./Config.pl", "-spice-path=nospice"],
        ],
        "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
        "files": files,
    })
    if not executable.is_file():
        raise RunnerError("profile {} did not retain an executable".format(
            profile_name))
    return executable, identity


def build_profile(
        context: Mapping[str, Any], repo: Path, artifacts: Path,
        profile_name: str, profile: Mapping[str, Any], timeout: int,
        lock_stream: Any) -> Tuple[Path, Path]:
    """Perform one complete, locked native-build preflight and build."""

    actual_cwd = Path.cwd().resolve()
    if actual_cwd != repo:
        raise RunnerError("native adapter cwd {} is not repository {}".format(
            actual_cwd, repo))
    print("[preflight:{}] repository={} cwd={} status=confirmed".format(
        profile_name, repo, actual_cwd), flush=True)
    active = active_repo_processes(repo)
    if active:
        raise RunnerError("native build/run already active before {}: {}".format(
            profile_name, "; ".join(active)))
    print("[preflight:{}] active_native_build_or_test_processes=0".format(
        profile_name), flush=True)

    # The caller owns the advisory lock for the entire multi-profile matrix.
    # Reassert that fact at every preflight so the receipt distinguishes the
    # portable lock from the best-effort /proc scan above.
    fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    build = (repo / "build").resolve()
    if build != repo / "build":
        raise RunnerError("resolved build target is not repository/build")
    print("[preflight:{}] generated_removal_target={} status=validated".format(
        profile_name, build), flush=True)
    removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
    if removed.returncode != 0:
        raise RunnerError("could not remove generated build for {}".format(
            profile_name))

    application = str(profile["application"])
    configured = run_command(
        [str(repo / "Config.pl"), "-application={}".format(application)],
        repo, timeout=timeout)
    if configured.returncode != 0:
        raise RunnerError("Config.pl could not select {}".format(application))
    spice = run_command(
        [str(repo / "Config.pl"), "-spice-path=nospice"],
        repo, timeout=timeout)
    if spice.returncode != 0:
        raise RunnerError("Config.pl could not disable SPICE for {}".format(
            profile_name))

    make_command = command(context, "make") + ["amps"]
    print("[preflight:{}] build_environment=MAKEFLAGS=-j16 command={}".format(
        profile_name, " ".join(make_command)), flush=True)
    built = run_command(
        make_command, repo, timeout=timeout, env={"MAKEFLAGS": "-j16"})
    if built.returncode != 0 or not (repo / "amps").is_file():
        raise RunnerError("native {} application did not build (make exit {})".format(
            profile_name, built.returncode))
    valid, message = validate_generated_configuration(repo, profile)
    print("[preflight:{}] generated_profile={} status={}".format(
        profile_name, message, "verified" if valid else "failed"), flush=True)
    if not valid:
        raise RunnerError(message)
    return retain_profile(repo, artifacts, profile_name, profile)


def profiles_for_cases(requested: Sequence[str]) -> List[str]:
    selected: List[str] = []
    for name in ("1d", "2d", "3d"):
        owned = set(str(value) for value in PROFILES[name]["cases"])
        if any(identifier in owned for identifier in requested):
            selected.append(name)
    return selected


def profile_for_case(identifier: str) -> str:
    for name in ("1d", "2d", "3d"):
        if identifier in PROFILES[name]["cases"]:
            return name
    raise RunnerError("no compile profile owns {}".format(identifier))


def main() -> int:
    context = load_context("Run native T3 AMR topology/geometry tests")
    suite_id = str(context["test"]["id"])
    requested = list(context["test"].get("selected_cases", ALL_CASES))
    if (not requested or any(not isinstance(value, str) for value in requested)
            or len(set(requested)) != len(requested)):
        raise RunnerError("selected_cases must be a unique nonempty case list")
    unknown = [identifier for identifier in requested
               if identifier not in ALL_CASES]
    if unknown:
        raise RunnerError("suite {} does not own selected cases {}".format(
            suite_id, unknown))

    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 5400))
    if not (repo / "Config.pl").is_file() or not (repo / "src").is_dir():
        return write_result(context, "ERROR",
                            "repository preflight failed: {}".format(repo))
    if Path.cwd().resolve() != repo:
        return write_result(context, "ERROR",
                            "native adapter must run at repository root")
    artifacts.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return write_result(
                context, "ERROR",
                "another native AMPS test holds {}".format(repo / "Config.pl"))

        retained: Dict[str, Path] = {}
        evidence: List[str] = []
        for profile_name in profiles_for_cases(requested):
            executable, identity = build_profile(
                context, repo, artifacts, profile_name,
                PROFILES[profile_name], timeout, lock_stream)
            retained[profile_name] = executable
            evidence.extend([str(identity), str(executable)])

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        ranks = str(context["test"].get(
            "mpi_ranks", context["execution"].get("mpi_ranks", 1)))
        observed: Dict[str, Dict[str, Any]] = {}
        transcript_chunks: List[str] = []
        for identifier in requested:
            profile_name = profile_for_case(identifier)
            executable = retained[profile_name]
            run = run_command(
                command(context, "mpiexec") +
                [rank_flag, ranks, str(executable), "--case", identifier],
                work, timeout=timeout,
                # Geometry fixtures are repository inputs, whereas every
                # executable runs from its isolated writable work directory.
                # Pass the resolved checkout explicitly so no case depends on
                # a site's TMPDIR ancestry or the caller's launch directory.
                env={"AMPS_TEST_REPO_ROOT": str(repo)})
            transcript_chunks.append(
                "===== {} [{}] =====\n{}".format(
                    identifier, profile_name, run.stdout or ""))
            found = re.findall(
                r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                run.stdout or "", flags=re.MULTILINE)
            exact = [entry for entry in found if entry[0] == identifier]
            foreign = [entry[0] for entry in found if entry[0] != identifier]
            if len(exact) == 1 and not foreign:
                _, status, message = exact[0]
                if run.returncode != 0 and status == "PASS":
                    status = "FAIL"
                    message = "PASS record with mpiexec exit {}".format(
                        run.returncode)
                observed[identifier] = case(
                    identifier, status, message or "ok")
            else:
                observed[identifier] = case(
                    identifier, "FAIL",
                    "expected one exact record; exact={} foreign={} exit={}".format(
                        len(exact), foreign, run.returncode))

        transcript = artifacts / "native-core-mesh-output.txt"
        transcript.write_text("\n".join(transcript_chunks), encoding="utf-8")
        evidence.append(str(transcript))
        for identifier in requested:
            if identifier not in observed:
                observed[identifier] = case(
                    identifier, "FAIL", "case was not executed")
        cases = [observed[identifier] for identifier in requested]
        status = "PASS" if all(
            item["status"] == "PASS" for item in cases) else "FAIL"
        hashes = {name: sha256(path) for name, path in retained.items()}
        return write_result(
            context, status,
            "{} selected T3 cases executed from {} clean native profiles".format(
                len(cases), len(retained)),
            cases=cases, evidence=evidence,
            metrics={
                "mpi_ranks": int(ranks),
                "native_profiles": len(retained),
                "native_processes": len(requested),
                "case_order": str(context["test"].get(
                    "case_order", "normal")),
                "shuffle_seed": int(context["test"].get(
                    "shuffle_seed", 0)),
                "executable_sha256": hashes,
            })
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
