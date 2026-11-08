#!/usr/bin/env python3
"""Build and execute the collective T10 AMR/PIC MPI matrix.

Every selected case is launched at 1, 2, and 4 ranks.  The native executable
prints one assertion record per rank; this adapter rejects missing, duplicate,
foreign, or rank-count-mismatched records before it accepts the root summary.
The periodic particle case is compiled in its own profile because periodic
ownership is an AMPS compile-time contract, not a runtime test switch.
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


BASELINE_CASES = [
    "MESH-MPI-OWNER-MAP-001",
    "MESH-MPI-OWNER-GHOST-001",
    "MESH-MPI-LOAD-REPARTITION-001",
    "MESH-MPI-CUTCELL-CONSISTENCY-001",
    "PIC-MPI-HALO-REPLACE-001",
    "PIC-MPI-HALO-ADD-001",
    "PIC-MPI-SEND-MASK-001",
    "PIC-MPI-BLOCK-PACK-001",
    "PIC-MPI-PARTICLE-EXCHANGE-001",
    "PIC-MPI-PARTICLE-REPARTITION-001",
    "PIC-MPI-OUTER-LOSS-001",
    "PIC-MPI-AFFINITY-DIAGNOSTICS-001",
    "PIC-MPI-RESTART-SAME-001",
    "PIC-MPI-RESTART-REPARTITION-001",
    "PIC-MPI-RANK-INVARIANCE-001",
]
PERIODIC_CASES = ["PIC-MPI-PERIODIC-EXCHANGE-001"]
ALL_CASES = BASELINE_CASES + PERIODIC_CASES
RANK_MATRIX = (1, 2, 4)

PROFILES: Dict[str, Dict[str, Any]] = {
    "baseline": {
        "application": "test/core-mpi",
        "input_deck": "input/test/core-mpi.input",
        "cases": BASELINE_CASES,
        "periodic": False,
    },
    "periodic": {
        "application": "test/core-mpi-periodic",
        "input_deck": "input/test/core-mpi-periodic.input",
        "cases": PERIODIC_CASES,
        "periodic": True,
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


def active_repo_processes(repo: Path) -> List[str]:
    """Return native build/test processes rooted in this checkout."""

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active: List[str] = []
    own_pid = os.getpid()
    names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl",
             "mpiexec", "mpirun"}
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


def profile_for_case(identifier: str) -> str:
    for name in ("baseline", "periodic"):
        if identifier in PROFILES[name]["cases"]:
            return name
    raise RunnerError("no T10 compile profile owns {}".format(identifier))


def profiles_for_cases(requested: Sequence[str]) -> List[str]:
    return [name for name in ("baseline", "periodic")
            if any(identifier in PROFILES[name]["cases"]
                   for identifier in requested)]


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
    global_def = (repo / "build" / "pic" / "picGlobal.dfn").read_text(
        encoding="utf-8", errors="replace")
    application = str(profile["application"]).rsplit("/", 1)[-1]
    input_name = Path(str(profile["input_deck"])).name
    checks = {
        "application": "APPLICATION={}".format(application) in selection,
        "input": "InputFileAMPS={}".format(input_name) in selection,
        "spice": "SPICE=nospice" in selection,
        "dimension": re.search(r"^\s*#\s*define\s+DIM\s+3\s*$",
                               global_def, re.MULTILINE) is not None,
    }
    if bool(profile["periodic"]):
        checks["periodic"] = re.search(
            r"^\s*#\s*define\s+_PIC_BC__PERIODIC_MODE_\s+"
            r"_PIC_BC__PERIODIC_MODE_ON_\s*$",
            global_def, re.MULTILINE) is not None
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        return False, "generated profile checks failed: {}".format(
            ", ".join(failed))
    return True, "application/input/SPICE/DIM/periodic profile verified"


def retain_profile(repo: Path, artifacts: Path, name: str,
                   profile: Mapping[str, Any]) -> Tuple[Path, Path]:
    retained = artifacts / "native-build" / name
    retained.mkdir(parents=True, exist_ok=True)
    executable = retained / "amps-core-mpi-{}".format(name)
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
        "stage": "T10",
        "profile": name,
        "application": str(profile["application"]),
        "input_deck": str(profile["input_deck"]),
        "periodic": bool(profile["periodic"]),
        "configure_commands": [
            ["./Config.pl", "-application={}".format(profile["application"])],
            ["./Config.pl", "-spice-path=nospice"],
        ],
        "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
        "files": files,
    })
    if not executable.is_file():
        raise RunnerError("profile {} did not retain an executable".format(name))
    return executable, identity


def build_profile(context: Mapping[str, Any], repo: Path, artifacts: Path,
                  name: str, profile: Mapping[str, Any], timeout: int,
                  lock_stream: Any) -> Tuple[Path, Path]:
    """Run the mandatory native rebuild preflight and retain its products."""

    actual_cwd = Path.cwd().resolve()
    if actual_cwd != repo:
        raise RunnerError("native adapter cwd {} is not {}".format(
            actual_cwd, repo))
    print("[preflight:{}] repository={} cwd={} status=confirmed".format(
        name, repo, actual_cwd), flush=True)
    active = active_repo_processes(repo)
    if active:
        raise RunnerError("native build/run already active: {}".format(
            "; ".join(active)))
    print("[preflight:{}] active_native_build_or_test_processes=0".format(
        name), flush=True)
    fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    build = (repo / "build").resolve()
    if build != repo / "build":
        raise RunnerError("resolved removal target is not repository/build")
    print("[preflight:{}] generated_removal_target={} status=validated".format(
        name, build), flush=True)
    removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
    if removed.returncode != 0:
        raise RunnerError("could not remove generated build directory")
    selected = run_command(
        [str(repo / "Config.pl"),
         "-application={}".format(profile["application"])],
        repo, timeout=timeout)
    if selected.returncode != 0:
        raise RunnerError("Config.pl could not select {}".format(
            profile["application"]))
    spice = run_command([str(repo / "Config.pl"), "-spice-path=nospice"],
                        repo, timeout=timeout)
    if spice.returncode != 0:
        raise RunnerError("Config.pl could not disable SPICE")
    make_command = command(context, "make") + ["amps"]
    print("[preflight:{}] build_environment=MAKEFLAGS=-j16 command={}".format(
        name, " ".join(make_command)), flush=True)
    built = run_command(make_command, repo, timeout=timeout,
                        env={"MAKEFLAGS": "-j16"})
    if built.returncode != 0 or not (repo / "amps").is_file():
        raise RunnerError("native {} build failed (exit {})".format(
            name, built.returncode))
    valid, message = validate_generated_configuration(repo, profile)
    print("[preflight:{}] {} status={}".format(
        name, message, "verified" if valid else "failed"), flush=True)
    if not valid:
        raise RunnerError(message)
    return retain_profile(repo, artifacts, name, profile)


def main() -> int:
    context = load_context("Run collective T10 AMR/PIC MPI tests")
    suite_id = str(context["test"]["id"])
    requested = list(context["test"].get("selected_cases", ALL_CASES))
    if (not requested or any(not isinstance(value, str) for value in requested)
            or len(set(requested)) != len(requested)):
        raise RunnerError("selected_cases must be a unique nonempty list")
    unknown = [identifier for identifier in requested
               if identifier not in ALL_CASES]
    if unknown:
        raise RunnerError("suite {} does not own {}".format(suite_id, unknown))

    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 7200))
    if Path.cwd().resolve() != repo:
        return write_result(context, "ERROR",
                            "native adapter must run at repository root")
    work.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)

    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return write_result(context, "ERROR",
                                "another native AMPS test holds Config.pl")

        retained: Dict[str, Path] = {}
        evidence: List[str] = []
        for profile_name in profiles_for_cases(requested):
            executable, identity = build_profile(
                context, repo, artifacts, profile_name,
                PROFILES[profile_name], timeout, lock_stream)
            retained[profile_name] = executable
            evidence.extend([str(identity), str(executable)])

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        observed: Dict[str, Dict[str, Any]] = {}
        matrix_receipts: Dict[str, Dict[str, str]] = {}
        transcript_chunks: List[str] = []
        total_processes = 0
        local_pattern = re.compile(
            r"^\[T10\]\[case=([A-Z][A-Z0-9-]+)\]"
            r"\[rank=(\d+)/(\d+)\]\s+(PASS|FAIL)\s*(.*)$",
            re.MULTILINE)
        summary_pattern = re.compile(
            r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
            re.MULTILINE)

        for identifier in requested:
            profile_name = profile_for_case(identifier)
            executable = retained[profile_name]
            failures: List[str] = []
            matrix_receipts[identifier] = {}
            for ranks in RANK_MATRIX:
                run_dir = work / identifier.lower() / "n{}".format(ranks)
                run_dir.mkdir(parents=True, exist_ok=True)
                invocation = command(context, "mpiexec") + [
                    rank_flag, str(ranks), str(executable),
                    "--case", identifier]
                run = run_command(invocation, run_dir, timeout=timeout,
                                  env={"AMPS_TEST_REPO_ROOT": str(repo)})
                total_processes += ranks
                transcript_chunks.append(
                    "===== {} [{} n={}] exit={} =====\n{}".format(
                        identifier, profile_name, ranks, run.returncode,
                        run.stdout or ""))
                local = [entry for entry in local_pattern.findall(
                    run.stdout or "") if entry[0] == identifier]
                foreign_local = [entry[0] for entry in local_pattern.findall(
                    run.stdout or "") if entry[0] != identifier]
                local_ranks = sorted(int(entry[1]) for entry in local)
                local_sizes = {int(entry[2]) for entry in local}
                local_failures = [entry for entry in local if entry[3] != "PASS"]
                summaries = [entry for entry in summary_pattern.findall(
                    run.stdout or "") if entry[0] == identifier]
                foreign_summaries = [entry[0] for entry in summary_pattern.findall(
                    run.stdout or "") if entry[0] != identifier]
                valid = (run.returncode == 0 and not foreign_local and
                         not foreign_summaries and not local_failures and
                         local_ranks == list(range(ranks)) and
                         local_sizes == {ranks} and len(summaries) == 1 and
                         summaries[0][1] == "PASS")
                if not valid:
                    failures.append(
                        "n={}: exit={} local_ranks={} local_sizes={} "
                        "local_failures={} summaries={} foreign={}".format(
                            ranks, run.returncode, local_ranks,
                            sorted(local_sizes), len(local_failures),
                            summaries, foreign_local + foreign_summaries))
                elif valid:
                    matrix_receipts[identifier][str(ranks)] = summaries[0][2]

            if identifier == "PIC-MPI-RANK-INVARIANCE-001" and not failures:
                receipts = list(matrix_receipts[identifier].values())
                if len(set(receipts)) != 1:
                    failures.append(
                        "1/2/4-rank invariant receipts differ: {}".format(
                            matrix_receipts[identifier]))
            observed[identifier] = case(
                identifier, "FAIL" if failures else "PASS",
                "; ".join(failures) if failures else
                "1/2/4-rank collective matrix passed")

        transcript = artifacts / "native-t10-mpi-output.txt"
        transcript.write_text("\n".join(transcript_chunks), encoding="utf-8")
        receipts_path = artifacts / "t10-matrix-receipts.json"
        write_json(receipts_path, {
            "schema_version": 1,
            "rank_matrix": list(RANK_MATRIX),
            "receipts": matrix_receipts,
        })
        evidence.extend([str(transcript), str(receipts_path)])
        cases_result = [observed[identifier] for identifier in requested]
        status = "PASS" if all(item["status"] == "PASS"
                               for item in cases_result) else "FAIL"
        return write_result(
            context, status,
            "{} T10 cases executed at 1, 2, and 4 ranks".format(
                len(cases_result)),
            cases=cases_result, evidence=evidence,
            metrics={
                "rank_matrix": list(RANK_MATRIX),
                "native_profiles": len(retained),
                "native_invocations": len(requested) * len(RANK_MATRIX),
                "mpi_processes": total_processes,
                "case_order": str(context["test"].get(
                    "case_order", "normal")),
                "shuffle_seed": int(context["test"].get(
                    "shuffle_seed", 0)),
                "executable_sha256": {
                    name: sha256(path) for name, path in retained.items()},
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
