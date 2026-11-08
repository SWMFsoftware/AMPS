#!/usr/bin/env python3
"""Build and run the native one-rank spherical dust-source qualification.

The adapter owns a complete Config.pl/top-level AMPS build, so it is marked
exclusive in the manifest.  It follows the repository's native-build
preflight literally and retains enough identity and transcript evidence to
distinguish the generated DUST/individual-weight executable from every other
profile.  The application itself uses the production spherical injector,
native AMR owner lookup, particle buffer, sphere sampling buffer, and boundary
mover; its deterministic source callback removes stochastic sampling error
from the conservation oracle.
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
from typing import List


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


CASE_ID = "CORE-DUST-INJECTION-INTEGRATION-001"


def _active_repo_processes(repo: Path) -> List[str]:
    """Return native builders/runs rooted in this checkout, excluding self."""

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active: List[str] = []
    names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl"}
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
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
                any(token in raw for token in ("Config.pl", "ampsConfig.pl"))):
            active.append("pid={} cwd={} command={}".format(
                entry.name, cwd, raw))
    return active


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    context = load_context("Run native AMPS spherical dust injection test")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", [CASE_ID])))
    if requested != [CASE_ID]:
        raise RunnerError("dust injection suite owns exactly {}".format(CASE_ID))

    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 1800))
    if Path.cwd().resolve() != repo:
        return write_result(
            context, "ERROR", "native dust injection cwd is not repository root")
    if not (repo / "Config.pl").is_file() or not (repo / "src").is_dir():
        return write_result(context, "ERROR", "repository preflight failed")

    print("[preflight] repository={} cwd={} status=confirmed".format(
        repo, Path.cwd().resolve()), flush=True)
    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return write_result(
                context, "ERROR", "another native AMPS test holds Config.pl")

        active = _active_repo_processes(repo)
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
            return write_result(
                context, "ERROR", "could not remove generated build directory")

        configured = run_command(
            [str(repo / "Config.pl"),
             "-application=test/model-dust-injection"],
            repo, timeout=timeout)
        if configured.returncode != 0:
            return write_result(
                context, "FAIL",
                "Config.pl could not select test/model-dust-injection",
                cases=[case(CASE_ID, "FAIL", "configuration failed")])

        make_command = command(context, "make") + ["amps"]
        print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
            " ".join(make_command)), flush=True)
        built = run_command(
            make_command, repo, timeout=timeout, env={"MAKEFLAGS": "-j16"})
        artifacts.mkdir(parents=True, exist_ok=True)
        build_log = artifacts / "dust-injection-native-build.txt"
        build_log.write_text(
            (configured.stdout or "") + (built.stdout or ""), encoding="utf-8")
        if built.returncode != 0 or not (repo / "amps").is_file():
            return write_result(
                context, "FAIL", "native dust injection application did not build",
                cases=[case(CASE_ID, "FAIL", "make exit {}".format(
                    built.returncode))], evidence=[str(build_log)])

        retained = artifacts / "amps-test-model-dust-injection"
        shutil.copy2(str(repo / "amps"), str(retained))
        identity = artifacts / "dust-injection-native-build-identity.json"
        identity.write_text(json.dumps({
            "schema_version": 1,
            "application": "test/model-dust-injection",
            "configure_command": [
                "./Config.pl", "-application=test/model-dust-injection"],
            "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
            "executable_sha256": _sha256(retained),
            "input_deck_sha256": _sha256(
                repo / "input/test/model-dust-injection.input"),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        print("[CASE] RUN  {} production spherical source and AMR ledger".format(
            CASE_ID), flush=True)
        run = run_command(
            command(context, "mpiexec") + [rank_flag, "1", str(retained)],
            work, timeout=timeout)
        transcript = artifacts / "dust-injection-native-output.txt"
        transcript.write_text(run.stdout or "", encoding="utf-8")
        records = re.findall(
            r"^\[{}\]\s+(PASS|FAIL)\s*$".format(re.escape(CASE_ID)),
            run.stdout or "", flags=re.MULTILINE)
        exact_receipt = run.returncode == 0 and records == ["PASS"]
        message = ("native spherical dust injection closed owner and mass ledgers"
                   if exact_receipt else
                   "invalid native receipt: exit={} records={}".format(
                       run.returncode, records))
        return write_result(
            context, "PASS" if exact_receipt else "FAIL", message,
            cases=[case(CASE_ID, "PASS" if exact_receipt else "FAIL", message)],
            evidence=[str(build_log), str(identity), str(transcript)])
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
