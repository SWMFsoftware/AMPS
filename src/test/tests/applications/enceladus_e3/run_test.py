#!/usr/bin/env python3
"""Own the bounded 10,000-step Enceladus E3 native integration campaign.

Routine workers publish every E3 identity as SKIP unless the explicit
``enceladus-e3`` machine capability is enabled.  On an enabled host the suite
performs a clean native build and delegates all physical/numerical checks to
the application-local runner; no expected value is duplicated here.
"""

from __future__ import annotations

import fcntl
import json
import os
import shutil
import sys
from pathlib import Path
from typing import List, Mapping


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


LOCAL_CASES = ["E3R{:02d}_{}".format(index, suffix) for index, suffix in enumerate([
    "process", "configuration", "global_time_step", "completed_window",
    "exact_injection", "co2_transport_only", "particle_conservation",
    "no_plasma_output", "active_species_output", "clean_termination",
    "flyby_tecplot", "flyby_mesh_resolution", "flyby_png",
])]
EXPECTED_CASES = [
    "ENCELADUS-E3R00-PROCESS", "ENCELADUS-E3R01-CONFIGURATION",
    "ENCELADUS-E3R02-GLOBAL-TIME-STEP", "ENCELADUS-E3R03-COMPLETED-WINDOW",
    "ENCELADUS-E3R04-EXACT-INJECTION", "ENCELADUS-E3R05-CO2-TRANSPORT-ONLY",
    "ENCELADUS-E3R06-PARTICLE-CONSERVATION",
    "ENCELADUS-E3R07-NO-PLASMA-OUTPUT",
    "ENCELADUS-E3R08-ACTIVE-SPECIES-OUTPUT",
    "ENCELADUS-E3R09-CLEAN-TERMINATION", "ENCELADUS-E3R10-FLYBY-TECPLOT",
    "ENCELADUS-E3R11-FLYBY-MESH-RESOLUTION", "ENCELADUS-E3R12-FLYBY-PNG",
]


def active_repo_processes(repo: Path) -> List[str]:
    """Find checkout-local native state owners before deleting ``build``.

    The E3 suite changes generated headers and the repository-wide executable,
    so PID discovery is a safety check in addition to the advisory Config.pl
    lock.  Only commands capable of owning native generated state are matched;
    unrelated readers and documentation processes do not block the campaign.
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


def main() -> int:
    context = load_context("Run the bounded native Enceladus E3 campaign")
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    if (not requested or len(set(requested)) != len(requested) or
            any(identifier not in EXPECTED_CASES for identifier in requested)):
        raise RunnerError("invalid Enceladus E3 case selection")
    capability = context.get("capabilities", {}).get("enceladus-e3", {})
    available = (capability if isinstance(capability, bool) else
                 capability.get("available", False))
    if not available:
        # The long campaign is deliberately opt-in because a routine worker may
        # lack a permitted MPI transport or sufficient wall time.  Publishing
        # all 13 SKIP rows keeps the verification boundary visible while
        # preventing registry presence from being mistaken for a scientific
        # run.  No build or simulation side effect occurs on this path.
        reason = (
            "not executed: optional machine capability 'enceladus-e3' is unavailable; "
            "enable it only on a host approved for the 10,000-step MPI campaign")
        rows = [case(identifier, "SKIP", reason) for identifier in requested]
        for row in rows:
            print("[CASE] RUN  {} {}".format(row["id"], reason), flush=True)
        return write_result(context, "SKIP", reason, cases=rows,
                            metrics={"execution": "not-run"})

    repo = Path(str(context["repo_root"])).resolve()
    artifacts = Path(str(context["artifact_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    timeout = int(context["test"].get("timeout_seconds", 14400))
    if Path.cwd().resolve() != repo or not (repo / "Config.pl").is_file():
        raise RunnerError("Enceladus E3 preflight is not at repository root")
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
        # This literal, root-validated target is the only destructive step.
        # Configuration and compilation then follow the same production hooks
        # as an operator run; a source-only compile would not prove generated
        # species, chemistry, mover, or runtime-input integration.
        removed = run_command(["rm", "-rf", "--", "build"], repo, timeout=60)
        configured = run_command(
            [str(repo / "Config.pl"), "-application=enceladus"],
            repo, timeout=timeout)
        make_argv = command(context, "make") + ["amps"]
        print("[preflight] build_environment=MAKEFLAGS=-j16 command={}".format(
            " ".join(make_argv)), flush=True)
        built = run_command(make_argv, repo, timeout=timeout,
                            env={"MAKEFLAGS": "-j16"})
        build_log = artifacts / "enceladus-e3-build.txt"
        build_log.write_text(
            (removed.stdout or "") + (configured.stdout or "") + (built.stdout or ""),
            encoding="utf-8")
        if (removed.returncode != 0 or configured.returncode != 0 or
                built.returncode != 0 or not (repo / "amps").is_file()):
            raise RunnerError("native Enceladus E3 build failed")
        retained = artifacts / "amps-enceladus-e3"
        shutil.copy2(str(repo / "amps"), str(retained))

        for identifier in requested:
            print("[CASE] RUN  {} bounded native E3 assertion".format(
                identifier), flush=True)
        local_output = artifacts / "e3-campaign"
        rank_flag = str(context["execution"].get("mpi_rank_flag", "-n"))
        # These values reproduce the bounded E3 validation contract: 10,000
        # global iterations, ten 1,000-step sampling windows, exactly 100
        # primary macroparticles per active species per step, neutral/loss-only
        # output, and no in-domain flyby cell coarser than 20 km.  The local
        # runner derives the physical time step and checks the resulting
        # particle and species ledgers; the adapter does not invent an oracle.
        argv = command(context, "python") + [
            str(repo / "srcEnceladus" / "test" / "run_e3_simulation.py"),
            "--executable", str(retained), "--input", str(repo / "input" / "enceladus.in"),
            "--output", str(local_output), "--ranks", "1", "--launcher", "mpiexec",
            "--timeout-s", str(timeout), "--expected-iterations", "10000",
            "--expected-sample-length", "1000", "--particles-per-species-step", "100",
            "--maximum-flyby-cell-edge-m", "20000", "--expect-volume-fields", "off",
            "--rank-flag={}".format(rank_flag), "--launcher-command-json",
            json.dumps(command(context, "mpiexec")),
        ]
        completed = run_command(
            argv, repo, timeout=timeout, env={"PYTHONDONTWRITEBYTECODE": "1"})
        outer = artifacts / "enceladus-e3-runner.txt"
        outer.write_text(completed.stdout or "", encoding="utf-8")
        receipt_path = local_output / "test_receipt.json"
        try:
            payload = json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise RunnerError("cannot read E3 receipt: {}".format(error))
        checks = payload.get("checks") if isinstance(payload, dict) else None
        if not isinstance(checks, list):
            raise RunnerError("E3 receipt has no checks array")
        observed = [str(row.get("test")) for row in checks if isinstance(row, dict)]
        if observed != LOCAL_CASES or len(observed) != len(set(observed)):
            raise RunnerError("E3 registry mismatch: {}".format(observed))
        translated = {}
        for global_id, row in zip(EXPECTED_CASES, checks):
            status = str(row.get("status"))
            if status not in {"PASS", "FAIL"}:
                raise RunnerError("unsupported E3 status {!r}".format(status))
            translated[global_id] = case(
                global_id, status, str(row.get("detail", "")))
        selected_rows = [translated[identifier] for identifier in requested]
        all_pass = all(row["status"] == "PASS" for row in translated.values())
        if completed.returncode != (0 if all_pass else 1):
            raise RunnerError("E3 exit/receipt disagreement")
        status = "FAIL" if any(row["status"] == "FAIL" for row in selected_rows) else "PASS"
        evidence = [str(build_log), str(outer), str(receipt_path),
                    str(local_output / "run.log")]
        return write_result(
            context, status,
            "imported {} bounded native E3 outcomes".format(len(selected_rows)),
            cases=selected_rows, evidence=evidence,
            metrics={"local_elapsed_s": payload.get("elapsed_s")})
    except Exception as error:
        return error_result(context, error)
    finally:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        finally:
            lock_stream.close()


if __name__ == "__main__":
    main_guard(main)
