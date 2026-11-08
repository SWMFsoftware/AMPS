#!/usr/bin/env python3
"""Build and run the configured legacy dust charging application.

This exclusive adapter owns the repository-wide generated AMPS configuration.
It follows the mandated native preflight exactly: verify the checkout/cwd and
absence of competing native processes, remove only ``build``, select
``test/model-dust`` through ``Config.pl``, and compile through top-level
``make amps`` with ``MAKEFLAGS=-j16``.  The application writes its numerical
receipt to this suite's artifact directory, never to the user's PT/plots tree.
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


CASE_ID = "CORE-DUST-CHARGING-DYNAMICS-001"


def _active_repo_processes(repo: Path) -> List[str]:
    """Find legacy build/run processes not participating in the file lock."""

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
    context = load_context("Run configured AMPS dust charging dynamics test")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", [CASE_ID])))
    if requested != [CASE_ID]:
        raise RunnerError("dust native suite owns exactly {}".format(CASE_ID))

    repo = Path(context["repo_root"]).resolve()
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 1800))
    if Path.cwd().resolve() != repo:
        return write_result(
            context, "ERROR", "native dust adapter cwd is not repository root")
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
            [str(repo / "Config.pl"), "-application=test/model-dust"],
            repo, timeout=timeout)
        if configured.returncode != 0:
            return write_result(
                context, "FAIL", "Config.pl could not select test/model-dust",
                cases=[case(CASE_ID, "FAIL", "configuration failed")])

        make_command = command(context, "make") + ["amps"]
        print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
            " ".join(make_command)), flush=True)
        built = run_command(
            make_command, repo, timeout=timeout, env={"MAKEFLAGS": "-j16"})
        artifacts.mkdir(parents=True, exist_ok=True)
        build_log = artifacts / "dust-native-build.txt"
        build_log.write_text(
            (configured.stdout or "") + (built.stdout or ""), encoding="utf-8")
        if built.returncode != 0 or not (repo / "amps").is_file():
            return write_result(
                context, "FAIL", "native dust application did not build",
                cases=[case(CASE_ID, "FAIL", "make exit {}".format(
                    built.returncode))], evidence=[str(build_log)])

        retained = artifacts / "amps-test-model-dust"
        shutil.copy2(str(repo / "amps"), str(retained))
        identity = artifacts / "dust-native-build-identity.json"
        identity.write_text(json.dumps({
            "schema_version": 1,
            "application": "test/model-dust",
            "configure_command": ["./Config.pl", "-application=test/model-dust"],
            "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
            "executable_sha256": _sha256(retained),
            "input_deck_sha256": _sha256(repo / "input/test/model-dust.input"),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        model_output = artifacts / "test_model-dust.dat"
        print("[CASE] RUN  {} configured equilibrium and transient grain charging".format(
            CASE_ID), flush=True)
        run = run_command(
            command(context, "mpiexec") + [rank_flag, "1", str(retained)],
            work, timeout=timeout,
            env={"AMPS_TEST_MODEL_DUST_OUTPUT": str(model_output)})
        transcript = artifacts / "dust-native-output.txt"
        transcript.write_text(run.stdout or "", encoding="utf-8")

        record = re.findall(
            r"^\[{}\]\s+(PASS|FAIL)\s*$".format(re.escape(CASE_ID)),
            run.stdout or "", flags=re.MULTILINE)
        numerical = model_output.read_text(
            encoding="utf-8", errors="replace") if model_output.is_file() else ""
        radii = re.findall(
            r"Charging of a grain with radius\s+([0-9.eE+-]+).*?"
            r"Equilibrium charge:\s+([0-9.eE+-]+)",
            numerical, flags=re.DOTALL)
        exact_receipt = (run.returncode == 0 and record == ["PASS"] and
                         len(radii) == 4)
        message = ("configured production charging converged for four radii"
                   if exact_receipt else
                   "invalid native receipt: exit={} records={} radii={}".format(
                       run.returncode, record, len(radii)))
        result = case(CASE_ID, "PASS" if exact_receipt else "FAIL", message)
        return write_result(
            context, "PASS" if exact_receipt else "FAIL", message,
            cases=[result], evidence=[str(build_log), str(identity),
                                     str(transcript), str(model_output)])
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
