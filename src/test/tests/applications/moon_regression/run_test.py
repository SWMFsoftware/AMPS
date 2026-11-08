#!/usr/bin/env python3
"""Run the registered Moon U01--U26 and linked I11--I20 regression campaign.

This adapter owns only global orchestration and status translation.  The Moon
individual runners retain their physical oracles and artifacts.  Because Moon
configuration rewrites the checkout-wide generated AMPS tree, the suite is
exclusive and performs the repository-mandated preflight literally: verify
the root and process state, remove only ``build``, select Moon/SPICE through
``Config.pl``, and build ``amps`` with ``MAKEFLAGS=-j16``.  The local aggregate
then runs with ``--skip-build`` so it cannot silently perform a second,
different build.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Mapping


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


LOCAL_IDS = (["U{:02d}".format(index) for index in range(1, 27)] +
             ["I{:02d}".format(index) for index in range(11, 21)])
EXPECTED_CASES = ["MOON-{}".format(identifier) for identifier in LOCAL_IDS]
STATUS_MAP = {"PASS": "PASS", "FAIL": "FAIL", "SKIPPED": "SKIP", "ERROR": "ERROR"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def active_repo_processes(repo: Path) -> List[str]:
    """Find another native AMPS build/run rooted in this checkout.

    The global scheduler serializes declared exclusive suites, while this
    Linux `/proc` check protects against a separately launched make/AMPS job.
    The advisory lock below remains the portable cross-runner guard.
    """

    proc = Path("/proc")
    if not proc.is_dir():
        return []
    active: List[str] = []
    native_names = {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl"}
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
        if (executable in native_names or executable.startswith("amps") or
                any(token in raw for token in ("Config.pl", "ampsConfig.pl"))):
            active.append("pid={} cwd={} command={}".format(entry.name, cwd, raw))
    return active


def load_summary(path: Path) -> Mapping[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("cannot read Moon aggregate summary {}: {}".format(path, error))
    if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
        raise RunnerError("Moon aggregate summary has no results array")
    return payload


def main() -> int:
    context = load_context("Run the Moon application regression campaign")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid Moon regression case selection")

    repo = Path(str(context["repo_root"])).resolve()
    artifacts = Path(str(context["artifact_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    timeout = int(context["test"].get("timeout_seconds", 14400))
    capability = context.get("capabilities", {}).get("moon-spice", {})
    available = (capability if isinstance(capability, bool) else
                 capability.get("available", False))
    if not available:
        reason = (
            "not executed: optional machine capability 'moon-spice' is unavailable; "
            "set both SPICE paths and enable it in a machine overlay")
        rows = [case(identifier, "SKIP", reason) for identifier in requested]
        for row in rows:
            print("[CASE] RUN  {} {}".format(row["id"], reason), flush=True)
        return write_result(
            context, "SKIP", reason, cases=rows,
            metrics={"selected_case_count": len(rows), "execution": "not-run"})
    spice_path_raw = str(context["paths"].get("spice_path", ""))
    spice_kernels_raw = str(context["paths"].get("spice_kernels", ""))
    if not spice_path_raw or not spice_kernels_raw:
        raise RunnerError(
            "Moon regression requires paths.spice_path and paths.spice_kernels")
    spice_path = Path(spice_path_raw)
    spice_kernels = Path(spice_kernels_raw)
    if not spice_path.is_dir() or not spice_kernels.is_dir():
        raise RunnerError("configured Moon SPICE paths do not exist")
    if Path.cwd().resolve() != repo or not (repo / "Config.pl").is_file():
        raise RunnerError("Moon native preflight is not at the AMPS repository root")

    build_transcript = artifacts / "moon-native-build.txt"
    lock_stream = (repo / "Config.pl").open("r")
    try:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise RunnerError("another native AMPS test holds Config.pl")
        active = active_repo_processes(repo)
        if active:
            raise RunnerError("native build/run already active: {}".format("; ".join(active)))
        print("[preflight] repository={} cwd={} status=confirmed".format(
            repo, Path.cwd().resolve()), flush=True)
        print("[preflight] active_native_build_or_test_processes=0", flush=True)

        build = (repo / "build").resolve()
        if build != repo / "build":
            raise RunnerError("resolved build target is not repository/build")
        print("[preflight] generated_removal_target={} status=validated".format(build), flush=True)
        removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
        if removed.returncode != 0:
            raise RunnerError("could not remove generated build directory")

        configured = run_command([
            str(repo / "Config.pl"), "-application=moon",
            "-spice-path={}".format(spice_path),
            "-spice-kernels={}".format(spice_kernels),
        ], repo, timeout=timeout)
        make_argv = command(context, "make") + ["amps"]
        print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
            " ".join(make_argv)), flush=True)
        built = run_command(
            make_argv, repo, timeout=timeout, env={"MAKEFLAGS": "-j16"})
        build_transcript.write_text(
            (removed.stdout or "") + (configured.stdout or "") + (built.stdout or ""),
            encoding="utf-8")
        if configured.returncode != 0 or built.returncode != 0 or not (repo / "amps").is_file():
            rows = [case(identifier, "FAIL", "Moon native configuration/build failed")
                    for identifier in requested]
            return write_result(
                context, "FAIL", "Moon native configuration/build failed",
                cases=rows, evidence=[str(build_transcript)])

        retained = artifacts / "amps-moon"
        shutil.copy2(str(repo / "amps"), str(retained))
        identity = artifacts / "moon-build-identity.json"
        identity.write_text(json.dumps({
            "schema_version": 1,
            "application": "moon",
            "configure_command": [
                "./Config.pl", "-application=moon",
                "-spice-path={}".format(spice_path),
                "-spice-kernels={}".format(spice_kernels),
            ],
            "build_command": ["env", "MAKEFLAGS=-j16", "make", "amps"],
            "executable_sha256": sha256(retained),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        # Python 3.8 is the control-plane floor, so strip the already-validated
        # fixed prefix without using str.removeprefix (introduced in 3.9).
        selected_local = [identifier[len("MOON-"):] for identifier in requested]
        campaign_output = artifacts / "moon-campaign"
        argv = command(context, "python") + [
            str(repo / "srcMoon" / "test" / "run_all_test.py"),
            "--skip-build", "--build-dir", str(repo / "build"),
            "--amps", str(retained), "--output-dir", str(campaign_output),
            "--campaign-mode", "regression", "--mpi-ranks",
            str(int(context["test"].get("mpi_ranks", 8))),
            "--spice-path", str(spice_path), "--spice-kernels", str(spice_kernels),
            "--case-id-prefix", "MOON-",
        ]
        for local_id in selected_local:
            argv.extend(["--only", local_id])
        local = run_command(
            argv, repo, timeout=timeout, env={"PYTHONDONTWRITEBYTECODE": "1"})
        outer_transcript = artifacts / "moon-regression.txt"
        outer_transcript.write_text(local.stdout or "", encoding="utf-8")
        summary_path = campaign_output / "run_all_summary.json"
        payload = load_summary(summary_path)
        local_rows = payload["results"]
        observed = [str(row.get("id")) for row in local_rows if isinstance(row, dict)]
        if observed != selected_local or len(observed) != len(set(observed)):
            raise RunnerError(
                "Moon selected registry mismatch: expected {} observed {}".format(
                    selected_local, observed))

        translated: Dict[str, Dict[str, str]] = {}
        for row in local_rows:
            local_id = str(row["id"])
            local_status = str(row.get("status"))
            if local_status not in STATUS_MAP:
                raise RunnerError("unsupported Moon status {!r}".format(local_status))
            identifier = "MOON-{}".format(local_id)
            message = " ".join(str(row.get("reason") or row.get("name") or local_id).split())
            translated[identifier] = case(identifier, STATUS_MAP[local_status], message)
        translated_rows = [translated[identifier] for identifier in requested]
        severity = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}
        suite_status = max(
            (row["status"] for row in translated_rows), key=severity.__getitem__)
        # run_all_test returns zero for a scientifically incomplete regression
        # containing typed SKIPPED gates, one for FAIL, and two for ERROR.
        expected_exit = 2 if suite_status == "ERROR" else (1 if suite_status == "FAIL" else 0)
        if local.returncode != expected_exit:
            raise RunnerError(
                "Moon exit {} disagrees with translated {} status".format(
                    local.returncode, suite_status))
        return write_result(
            context, suite_status,
            "imported {} exact Moon regression outcomes".format(len(translated_rows)),
            cases=translated_rows,
            evidence=[str(build_transcript), str(identity), str(outer_transcript),
                      str(summary_path)],
            metrics={"selected_case_count": len(translated_rows),
                     "local_campaign_status": payload.get("campaign_status")})
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
