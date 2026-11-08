#!/usr/bin/env python3
"""Build and run the native ECSIM manufactured-operator qualification.

The legacy stencil program is useful only after its printed convergence tables
are converted into asserted case receipts.  This adapter performs the required
clean native build once, retains the generated executable/configuration, and
then runs every selected operator in a fresh one-rank process.  A zero process
exit without the exact ``[CASE] PASS`` record is a failure; diagnostic tables
alone are never accepted as numerical evidence.
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
from typing import Dict, List


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


CASES: Dict[str, List[str]] = {
    "ECSIM-STENCIL-CURL-001": ["curl_b", "--no-tecplot"],
    "ECSIM-STENCIL-GRAD-DIV-001": ["grad_div_e", "--no-tecplot"],
    "ECSIM-STENCIL-LAPLACIAN-001": ["laplacian_e_corner"],
    "ECSIM-STENCIL-CURL-CURL-001": ["curl_curl_e"],
    "ECSIM-STENCIL-MOMENTS-IDENTITIES-001": ["stencil_moments"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def active_repo_processes(repo: Path) -> List[str]:
    """Return competing native build/test processes rooted in this checkout."""

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


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(path))


def main() -> int:
    context = load_context("Run native ECSIM manufactured-operator tests")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in CASES for identifier in requested)):
        raise RunnerError("invalid ECSIM operator case selection: {}".format(
            requested))

    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 3600))
    if Path.cwd().resolve() != repo:
        return write_result(context, "ERROR",
                            "native ECSIM adapter cwd is not repository root")
    if not (repo / "Config.pl").is_file() or not (repo / "src").is_dir():
        return write_result(context, "ERROR", "repository preflight failed")
    artifacts.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return write_result(context, "ERROR",
                                "another native AMPS test holds Config.pl")

        print("[preflight] repository={} cwd={} status=confirmed".format(
            repo, Path.cwd().resolve()), flush=True)
        active = active_repo_processes(repo)
        if active:
            return write_result(
                context, "ERROR", "native build/run already active: {}".format(
                    "; ".join(active)))
        print("[preflight] active_native_build_or_test_processes=0", flush=True)

        build = (repo / "build").resolve()
        if build != repo / "build":
            raise RunnerError("resolved build target is not repository/build")
        print("[preflight] generated_removal_target={} status=validated".format(
            build), flush=True)
        removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
        if removed.returncode != 0:
            return write_result(context, "ERROR",
                                "could not remove generated build directory")

        configured = run_command(
            [str(repo / "Config.pl"), "-application=test/stencil"],
            repo, timeout=timeout)
        if configured.returncode != 0:
            return write_result(context, "FAIL",
                                "Config.pl could not select test/stencil")
        built = run_command(
            command(context, "make") + ["amps"], repo, timeout=timeout,
            env={"MAKEFLAGS": "-j16"})
        build_log = artifacts / "ecsim-operators-build.txt"
        build_log.write_text(
            (configured.stdout or "") + (built.stdout or ""), encoding="utf-8")
        if built.returncode != 0 or not (repo / "amps").is_file():
            return write_result(context, "FAIL",
                                "native ECSIM stencil application did not build",
                                evidence=[str(build_log)])

        executable = artifacts / "amps-test-stencil"
        shutil.copy2(str(repo / "amps"), str(executable))
        identity = artifacts / "ecsim-operators-build-identity.json"
        atomic_json(identity, {
            "schema_version": 1,
            "application": "test/stencil",
            "configure_command": ["./Config.pl", "-application=test/stencil"],
            "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
            "executable_sha256": sha256(executable),
            "input_deck_sha256": sha256(repo / "input/test/stencil.input"),
        })

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        results = []
        evidence = [str(build_log), str(identity), str(executable)]
        for identifier in requested:
            print("[CASE] RUN  {} manufactured periodic convergence".format(
                identifier), flush=True)
            case_work = work / identifier
            case_work.mkdir(parents=True, exist_ok=True)
            run = run_command(
                command(context, "mpiexec") +
                [rank_flag, "1", str(executable)] + CASES[identifier],
                case_work, timeout=timeout)
            transcript = artifacts / "{}.txt".format(identifier)
            transcript.write_text(run.stdout or "", encoding="utf-8")
            evidence.append(str(transcript))
            records = re.findall(
                r"^\[{}\]\s+(PASS|FAIL)\s*$".format(re.escape(identifier)),
                run.stdout or "", flags=re.MULTILINE)
            passed = run.returncode == 0 and records == ["PASS"]
            results.append(case(
                identifier, "PASS" if passed else "FAIL",
                "advertised order and independent manufactured solution accepted"
                if passed else
                "invalid receipt: exit={} records={}".format(
                    run.returncode, records)))

        passed = all(result["status"] == "PASS" for result in results)
        return write_result(
            context, "PASS" if passed else "FAIL",
            "{}/{} ECSIM manufactured-operator cases passed".format(
                sum(result["status"] == "PASS" for result in results),
                len(results)),
            cases=results, evidence=evidence)
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
