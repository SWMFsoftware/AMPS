#!/usr/bin/env python3
"""D03 configured-build, MPI, refresh, and restart release gate.

The script never substitutes a dependency-light driver for AMPS. In
development mode, an absent configured AMPS tree or site campaign manifest is
reported as a first-class SKIP. In release mode the same condition is an
INCOMPLETE gate and returns nonzero.

The site manifest contains only launch details and expected artifacts. Physics
selection remains in the production srcSEP executable: each transport case is
launched with one of the three canonical ``--particle-mover`` values. Exact
artifact checksums compare serial/MPI decompositions and uninterrupted/resumed
runs; no tolerance or result is computed by this orchestrator.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_MOVERS = ("parker", "fte-dmumu", "fte-mfp")
MOVER_OPTIONS = ("--particle-mover", "--mover", "--sep-mover")
SCHEMA = "srcsep-native-integration-v1"
EVIDENCE_SCHEMA = "srcsep-native-integration-evidence-v1"


class GateError(RuntimeError):
    """Infrastructure/manifest error distinct from a scientific case FAIL."""


def _utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def _tool_version(command: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            list(command), text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=15)
    except (OSError, subprocess.SubprocessError) as error:
        return f"unavailable: {error}"
    lines = completed.stdout.strip().splitlines()
    return lines[0] if lines else f"exit={completed.returncode}"


def _run(command: Sequence[str], cwd: Path, log: Path,
         environment: Mapping[str, str], timeout: Optional[float]) -> Tuple[int, str, float]:
    """Run one command, retain complete merged output, and echo it once.

    Native logs are intentionally not parsed from a pipe while still running;
    the short D03 cases have bounded output and subprocess.run avoids an
    additional reader thread in this standalone gate.
    """
    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command), cwd=str(cwd), env=dict(environment), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
            timeout=timeout)
        return_code = completed.returncode
        output = completed.stdout
    except subprocess.TimeoutExpired as error:
        output = (error.stdout or "")
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        output += f"\nD03 timeout after {timeout} seconds\n"
        return_code = 124
    except OSError as error:
        output = f"D03 launch failure: {error}\n"
        return_code = 127
    elapsed = time.monotonic() - started
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(output, encoding="utf-8")
    print(f"$ {shlex.join(list(command))}")
    if output:
        print(output, end="" if output.endswith("\n") else "\n")
    return return_code, output, elapsed


def validate_manifest(manifest: Mapping[str, Any]) -> List[str]:
    """Return all structural/coverage errors instead of stopping at the first.

    D03 is a release-evidence contract, so merely having executable commands is
    insufficient. The manifest must explicitly cover every mover in serial and
    two distinct MPI decompositions, SWCME refresh, exact decomposition
    equivalence, and checkpoint/resume equivalence.
    """
    if not isinstance(manifest, dict):
        return ["manifest root must be an object"]

    errors: List[str] = []
    if manifest.get("schema") != SCHEMA:
        errors.append(f"schema must be {SCHEMA}")
    if not isinstance(manifest.get("template_only", False), bool):
        errors.append("template_only must be Boolean when present")
    executable = manifest.get("executable", "amps")
    if not isinstance(executable, str) or not executable:
        errors.append("executable must be a non-empty string")
    environment = manifest.get("environment", {})
    if (not isinstance(environment, dict) or
            not all(isinstance(key, str) and isinstance(value, str)
                    for key, value in environment.items())):
        errors.append("environment must map strings to strings")
    common_args = manifest.get("common_args")
    if not isinstance(common_args, list) or not all(
            isinstance(item, str) for item in common_args):
        errors.append("common_args must be a string array")
    launcher = manifest.get("mpi_launcher", ["mpiexec", "-n", "{ranks}"])
    if not isinstance(launcher, list) or not launcher or not all(
            isinstance(item, str) for item in launcher):
        errors.append("mpi_launcher must be a non-empty string array")
    elif not any("{ranks}" in item for item in launcher):
        errors.append("mpi_launcher must contain {ranks}")

    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        return errors + ["cases must be a non-empty array"]

    ids: set[str] = set()
    ranks_by_mover: Dict[str, set[int]] = {mover: set() for mover in REQUIRED_MOVERS}
    groups_by_mover: Dict[str, Dict[str, set[int]]] = {
        mover: {} for mover in REQUIRED_MOVERS}
    restart_kinds: Dict[str, List[Mapping[str, Any]]] = {
        "uninterrupted": [], "checkpoint": [], "resume": []}
    restart_groups: Dict[str, set[str]] = {"uninterrupted": set(), "resume": set()}

    for index, raw_case in enumerate(cases):
        label = f"cases[{index}]"
        if not isinstance(raw_case, dict):
            errors.append(f"{label} must be an object")
            continue
        case_id = raw_case.get("id")
        if not isinstance(case_id, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", case_id):
            errors.append(f"{label}.id must be a filesystem-safe non-empty string")
        elif case_id in ids:
            errors.append(f"duplicate case id {case_id}")
        else:
            ids.add(case_id)
        kind = raw_case.get("kind")
        if kind not in ("transport", "uninterrupted", "checkpoint", "resume"):
            errors.append(f"{label}.kind is invalid")
        mover = raw_case.get("mover")
        if mover not in REQUIRED_MOVERS:
            errors.append(f"{label}.mover must be one of {', '.join(REQUIRED_MOVERS)}")
        ranks = raw_case.get("ranks")
        if not isinstance(ranks, int) or isinstance(ranks, bool) or ranks <= 0:
            errors.append(f"{label}.ranks must be a positive integer")
        args = raw_case.get("args", [])
        if not isinstance(args, list) or not all(isinstance(item, str) for item in args):
            errors.append(f"{label}.args must be a string array")
        expect = raw_case.get("expect", {})
        if not isinstance(expect, dict):
            errors.append(f"{label}.expect must be an object")
            expect = {}
        minimum_state_id = expect.get("minimum_source_state_id", 0)
        if (not isinstance(minimum_state_id, int) or
                isinstance(minimum_state_id, bool) or minimum_state_id < 2):
            errors.append(
                f"{label} must require minimum_source_state_id >= 2 to prove refresh")
        minimum_dispatches = expect.get("minimum_particle_dispatches", 0)
        if (not isinstance(minimum_dispatches, int) or
                isinstance(minimum_dispatches, bool) or minimum_dispatches < 1):
            errors.append(
                f"{label} must require minimum_particle_dispatches >= 1 "
                "to prove native mover execution")
        expected_exit = expect.get("exit_code", 0)
        if expected_exit != 0 or isinstance(expected_exit, bool):
            errors.append(f"{label}.expect.exit_code must be zero")
        log_regex = expect.get("log_regex", [])
        if (not isinstance(log_regex, list) or
                not all(isinstance(item, str) for item in log_regex)):
            errors.append(f"{label}.expect.log_regex must be a string array")
        case_environment = raw_case.get("environment", {})
        if (not isinstance(case_environment, dict) or
                not all(isinstance(key, str) and isinstance(value, str)
                        for key, value in case_environment.items())):
            errors.append(f"{label}.environment must map strings to strings")
        artifacts = raw_case.get("artifacts", [])
        if not isinstance(artifacts, list):
            errors.append(f"{label}.artifacts must be an array")
            artifacts = []
        groups: set[str] = set()
        for artifact_index, artifact in enumerate(artifacts):
            if not isinstance(artifact, dict) or not isinstance(artifact.get("path"), str):
                errors.append(f"{label}.artifacts[{artifact_index}] needs path")
                continue
            group = artifact.get("equivalence_group")
            if group is not None and not isinstance(group, str):
                errors.append(
                    f"{label}.artifacts[{artifact_index}].equivalence_group must be string")
            elif group:
                groups.add(group)

        if (kind == "transport" and mover in ranks_by_mover and
                isinstance(ranks, int) and not isinstance(ranks, bool) and ranks > 0):
            ranks_by_mover[mover].add(ranks)
            for group in groups:
                groups_by_mover[mover].setdefault(group, set()).add(ranks)
        if kind in restart_kinds:
            restart_kinds[kind].append(raw_case)
            if kind in restart_groups:
                restart_groups[kind].update(groups)

    for mover in REQUIRED_MOVERS:
        ranks = ranks_by_mover[mover]
        mpi_ranks = {value for value in ranks if value > 1}
        if 1 not in ranks or len(mpi_ranks) < 2:
            errors.append(
                f"{mover} needs rank 1 and at least two distinct MPI ranks; got {sorted(ranks)}")
        if ranks and not any(covered == ranks for covered in groups_by_mover[mover].values()):
            errors.append(
                f"{mover} needs one artifact equivalence_group shared by all decompositions")
    for kind, selected in restart_kinds.items():
        if not selected:
            errors.append(f"restart campaign is missing kind={kind}")
    if not restart_groups["uninterrupted"].intersection(restart_groups["resume"]):
        errors.append(
            "uninterrupted and resume cases need a shared artifact equivalence_group")
    return errors


def _expand(values: Iterable[str], replacements: Mapping[str, str]) -> List[str]:
    expanded: List[str] = []
    for value in values:
        for key, replacement in replacements.items():
            value = value.replace("{" + key + "}", replacement)
        expanded.append(value)
    return expanded


def _contains_mover_option(arguments: Iterable[str]) -> bool:
    """Recognize every production mover spelling in split or equals form."""
    return any(
        item in MOVER_OPTIONS or
        any(item.startswith(option + "=") for option in MOVER_OPTIONS)
        for item in arguments)


def _skip_or_incomplete(reason: str, release: bool, evidence_path: Path,
                        evidence: Dict[str, Any]) -> int:
    evidence.update({
        "status": "INCOMPLETE" if release else "SKIP",
        "message": reason,
        "completed_utc": _utc(),
    })
    _write_json(evidence_path, evidence)
    if release:
        print(f"D03 INCOMPLETE: {reason}", file=sys.stderr)
        return 2
    print("SRCSEP_SUITE_RESULT=SKIP")
    print(f"D03 SKIP: {reason}")
    return 0


def run_gate(args: argparse.Namespace) -> int:
    amps_source = Path(args.amps_source).expanduser().resolve()
    make_config = Path(args.make_config).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "native_integration_evidence.json"
    evidence: Dict[str, Any] = {
        "schema": EVIDENCE_SCHEMA,
        "started_utc": _utc(),
        "release_mode": bool(args.release),
        "rebuild_requested": bool(args.rebuild),
        "amps_source": str(amps_source),
        "make_config": str(make_config),
        "host": platform.platform(),
        "python": sys.version.splitlines()[0],
        "tools": {
            "make": _tool_version(["make", "--version"]),
            "compiler": _tool_version([os.environ.get("CXX", "c++"), "--version"]),
            "mpi_compiler": _tool_version([os.environ.get("MPICXX", "mpicxx"), "--version"]),
            "mpi_launcher": _tool_version([os.environ.get("MPIEXEC", "mpiexec"), "--version"]),
        },
        "build": {}, "cases": [], "artifact_groups": {},
    }
    if not amps_source.is_dir() or not (amps_source / "srcSEP" / "makefile").is_file():
        return _skip_or_incomplete(
            f"configured AMPS source tree is absent: {amps_source}",
            args.release, evidence_path, evidence)
    if not make_config.is_file():
        return _skip_or_incomplete(
            f"Makefile.conf is absent: {make_config}",
            args.release, evidence_path, evidence)
    if args.release and args.no_build:
        return _skip_or_incomplete(
            "release mode forbids --no-build", True, evidence_path, evidence)
    if args.release and not args.rebuild:
        return _skip_or_incomplete(
            "release mode requires --rebuild to prove a clean configured build",
            True, evidence_path, evidence)

    # Validate campaign coverage before starting a potentially expensive AMPS
    # build. A missing or copied-but-unedited template is a prerequisite issue,
    # not a compiler result, and must retain that distinction in evidence.
    manifest_path = Path(args.manifest).expanduser().resolve() if args.manifest else None
    if manifest_path is None or not manifest_path.is_file():
        return _skip_or_incomplete(
            "site native campaign manifest is absent; pass --native-manifest",
            args.release, evidence_path, evidence)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GateError(f"cannot read native manifest {manifest_path}: {error}")
    errors = validate_manifest(manifest)
    if errors:
        return _skip_or_incomplete(
            "invalid native manifest: " + "; ".join(errors),
            args.release, evidence_path, evidence)
    evidence["manifest"] = str(manifest_path)
    evidence["manifest_sha256"] = _sha256(manifest_path)
    if manifest.get("template_only", False):
        return _skip_or_incomplete(
            "native manifest is marked template_only; copy it, replace site "
            "placeholders, and remove template_only",
            args.release, evidence_path, evidence)

    cases_root = output_dir / "cases"
    if cases_root.exists() and any(cases_root.iterdir()):
        # Never accept an artifact left by an earlier run. Refusing a nonempty
        # directory preserves evidence without deleting user data and avoids an
        # expensive rebuild that could never yield trustworthy campaign data.
        raise GateError(
            f"native case output directory is not empty: {cases_root}; "
            "select a new --output-dir")

    environment = dict(os.environ)
    if not args.no_build:
        build_commands: List[List[str]] = []
        if args.rebuild:
            build_commands.append(["make", "-C", str(amps_source), "clean"])
        build_commands.append([
            "make", "--no-print-directory", "-f",
            str(amps_source / "srcSEP" / "makefile"), "strict-production",
            f"AMPS_ROOT={amps_source}", f"AMPS_CONFIG={make_config}"])
        for index, command in enumerate(build_commands):
            log = output_dir / f"build_{index + 1}.log"
            code, output, elapsed = _run(
                command, amps_source, log, environment, args.timeout)
            evidence["build"].setdefault("commands", []).append({
                "argv": command, "return_code": code,
                "elapsed_seconds": elapsed, "log": str(log),
                "log_sha256": _sha256(log),
            })
            if code != 0:
                evidence.update({"status": "FAIL", "message": "configured build failed",
                                 "completed_utc": _utc()})
                _write_json(evidence_path, evidence)
                return 1
        evidence["build"]["status"] = "PASS"
    else:
        evidence["build"]["status"] = "NOT_RUN"

    executable_value = manifest.get("executable", "amps")
    executable = Path(executable_value)
    if not executable.is_absolute():
        executable = amps_source / executable
    executable = executable.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        return _skip_or_incomplete(
            f"linked production executable is unavailable: {executable}",
            args.release, evidence_path, evidence)
    evidence["executable"] = str(executable)
    evidence["executable_sha256"] = _sha256(executable)
    evidence["make_config_sha256"] = _sha256(make_config)

    common_args = list(manifest.get("common_args", []))
    launcher = list(manifest.get("mpi_launcher", ["mpiexec", "-n", "{ranks}"]))
    common_environment = manifest.get("environment", {})
    if not isinstance(common_environment, dict) or not all(
            isinstance(k, str) and isinstance(v, str)
            for k, v in common_environment.items()):
        raise GateError("manifest environment must map strings to strings")

    case_directories: Dict[str, Path] = {}
    for case in manifest["cases"]:
        directory = output_dir / "cases" / case["id"]
        directory.mkdir(parents=True, exist_ok=True)
        case_directories[case["id"]] = directory

    group_records: Dict[str, List[Dict[str, str]]] = {}
    run_fingerprints: Dict[str, set[str]] = {mover: set() for mover in REQUIRED_MOVERS}
    swcme_fingerprints: set[str] = set()
    failures: List[str] = []
    for case in manifest["cases"]:
        case_id = case["id"]
        case_dir = case_directories[case_id]
        replacements = {
            "amps_source": str(amps_source), "output_dir": str(output_dir),
            "case_dir": str(case_dir), "executable": str(executable),
        }
        case_args = _expand(common_args + list(case.get("args", [])), replacements)
        # The runner, rather than a site manifest typo, owns the canonical mover
        # selector. Reject every production alias in both separated and equals
        # forms; relying on last-occurrence-wins would hide a contradictory site
        # input even though the appended canonical value ultimately prevailed.
        if _contains_mover_option(case_args):
            failures.append(
                f"{case_id}: args must not set a particle-mover option")
            continue
        case_args.extend(["--particle-mover", case["mover"]])
        ranks = case["ranks"]
        command = [str(executable)] + case_args
        if ranks > 1:
            command = [item.replace("{ranks}", str(ranks)) for item in launcher] + command
        case_environment = dict(environment)
        case_environment.update({
            key: _expand([value], replacements)[0]
            for key, value in common_environment.items()})
        raw_case_environment = case.get("environment", {})
        if not isinstance(raw_case_environment, dict):
            failures.append(f"{case_id}: environment must be an object")
            continue
        case_environment.update({
            str(key): _expand([str(value)], replacements)[0]
            for key, value in raw_case_environment.items()})
        log = case_dir / "run.log"
        code, output, elapsed = _run(
            command, case_dir, log, case_environment, args.timeout)
        expected_exit = int(case.get("expect", {}).get("exit_code", 0))
        record: Dict[str, Any] = {
            "id": case_id, "kind": case["kind"], "mover": case["mover"],
            "ranks": ranks, "argv": command, "return_code": code,
            "elapsed_seconds": elapsed, "log": str(log),
            "log_sha256": _sha256(log), "artifacts": [],
        }
        case_failures: List[str] = []
        if code != expected_exit:
            case_failures.append(f"exit {code}, expected {expected_exit}")
        if "mpi_consensus=pass" not in output:
            case_failures.append("missing mpi_consensus=pass")
        run_match = re.search(r"RunConfiguration fingerprint=([0-9a-fA-F]+)", output)
        swcme_match = re.search(
            r"SWCME configuration fingerprint=([0-9a-fA-F]+)", output)
        state_matches = re.findall(r"final_source_state_id=([0-9]+)", output)
        dispatch_matches = re.findall(
            r"completed_particle_dispatches=([0-9]+)", output)
        if not run_match:
            case_failures.append("missing RunConfiguration fingerprint")
        else:
            record["run_configuration_fingerprint"] = run_match.group(1).lower()
            if case["kind"] == "transport":
                run_fingerprints[case["mover"]].add(run_match.group(1).lower())
        if not swcme_match:
            case_failures.append("missing SWCME configuration fingerprint")
        else:
            record["swcme_configuration_fingerprint"] = swcme_match.group(1).lower()
            swcme_fingerprints.add(swcme_match.group(1).lower())
        minimum_state = int(case.get("expect", {}).get("minimum_source_state_id", 2))
        if not state_matches or int(state_matches[-1]) < minimum_state:
            case_failures.append(
                f"source state ID did not reach required {minimum_state}")
        minimum_dispatches = int(
            case.get("expect", {}).get("minimum_particle_dispatches", 1))
        if not dispatch_matches:
            case_failures.append("missing completed_particle_dispatches")
        else:
            completed_dispatches = int(dispatch_matches[-1])
            record["completed_particle_dispatches"] = completed_dispatches
            if completed_dispatches < minimum_dispatches:
                case_failures.append(
                    "completed particle dispatch count "
                    f"{completed_dispatches} is below required {minimum_dispatches}")
        for pattern in case.get("expect", {}).get("log_regex", []):
            if re.search(pattern, output, re.MULTILINE) is None:
                case_failures.append(f"missing log_regex {pattern!r}")

        for artifact in case.get("artifacts", []):
            artifact_path = Path(_expand([artifact["path"]], replacements)[0])
            if not artifact_path.is_absolute():
                artifact_path = case_dir / artifact_path
            artifact_path = artifact_path.resolve()
            if not artifact_path.is_file():
                case_failures.append(f"missing artifact {artifact_path}")
                continue
            artifact_record = {
                "path": str(artifact_path), "sha256": _sha256(artifact_path)}
            group = artifact.get("equivalence_group")
            if group:
                artifact_record["equivalence_group"] = group
                group_records.setdefault(group, []).append({
                    "case": case_id, "path": str(artifact_path),
                    "sha256": artifact_record["sha256"]})
            record["artifacts"].append(artifact_record)
        record["status"] = "PASS" if not case_failures else "FAIL"
        record["failures"] = case_failures
        evidence["cases"].append(record)
        failures.extend(f"{case_id}: {message}" for message in case_failures)

    for mover, fingerprints in run_fingerprints.items():
        if len(fingerprints) != 1:
            failures.append(
                f"{mover}: decomposition run fingerprints differ or are absent: "
                f"{sorted(fingerprints)}")
    if len(swcme_fingerprints) != 1:
        failures.append(
            "SWCME fingerprints differ across campaign: " +
            repr(sorted(swcme_fingerprints)))
    for group, records in sorted(group_records.items()):
        hashes = sorted({record["sha256"] for record in records})
        evidence["artifact_groups"][group] = {
            "status": "PASS" if len(hashes) == 1 else "FAIL",
            "sha256_values": hashes, "members": records,
        }
        if len(hashes) != 1:
            failures.append(f"artifact equivalence group {group} differs")

    evidence["completed_utc"] = _utc()
    evidence["status"] = "PASS" if not failures else "FAIL"
    evidence["message"] = (
        "clean build and native serial/MPI/refresh/restart evidence passed"
        if not failures else "; ".join(failures))
    _write_json(evidence_path, evidence)
    if failures:
        for failure in failures:
            print(f"D03 FAIL: {failure}", file=sys.stderr)
        return 1
    print(f"D03 PASS: evidence={evidence_path}")
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps-source", default=os.environ.get(
        "SRCSEP_D03_AMPS_SOURCE", str(ROOT.parent)))
    parser.add_argument("--make-config", default=os.environ.get(
        "SRCSEP_D03_MAKE_CONFIG", str(ROOT.parent / "Makefile.conf")))
    parser.add_argument("--native-manifest", dest="manifest", default=os.environ.get(
        "SRCSEP_D03_NATIVE_MANIFEST"))
    parser.add_argument("--output-dir", default=os.environ.get(
        "SRCSEP_D03_OUTPUT_DIR", str(ROOT / "test_output" / "d03-native")))
    parser.add_argument("--release", action="store_true", default=(
        os.environ.get("SRCSEP_D03_RELEASE", "").lower() in ("1", "on", "true", "yes")))
    parser.add_argument("--rebuild", action="store_true", default=(
        os.environ.get("SRCSEP_D03_REBUILD", "").lower() in ("1", "on", "true", "yes")))
    parser.add_argument("--no-build", action="store_true", default=(
        os.environ.get("SRCSEP_D03_NO_BUILD", "").lower() in ("1", "on", "true", "yes")))
    parser.add_argument("--timeout", type=float, default=(
        float(os.environ["SRCSEP_D03_TIMEOUT"])
        if os.environ.get("SRCSEP_D03_TIMEOUT") else None))
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.rebuild and args.no_build:
        raise GateError("--rebuild and --no-build are mutually exclusive")
    if args.timeout is not None and args.timeout <= 0.0:
        raise GateError("--timeout must be positive")
    return run_gate(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except GateError as error:
        print(f"D03 ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
