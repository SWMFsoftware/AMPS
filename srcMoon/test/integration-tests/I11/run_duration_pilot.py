#!/usr/bin/env python3
"""Run the non-qualifying 480 km I11 physical-duration pilot.

The formal I11 contract remains 500 particle steps per sample and 5000 steps
per resolution.  This separate runner removes only that stop at the coarsest
approved grid.  It lets the unchanged production-side independence,
equilibrium, and confidence-interval gates identify a candidate physical
duration.  A successful pilot is deliberately reported as ``COMPLETE`` rather
than ``PASS``: it is input to a later reviewed four-grid contract and is not
itself M1 evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import (  # noqa: E402
    REPO,
    derive_i11_input,
    descendant_pids,
    git_value,
    host_available_kib,
    manifest,
    parse_receipt,
    proc_rss_kib,
    qualification_levels,
    render_i11_stage_diagnostic,
    sha256,
)

PILOT_RESOLUTION_KM = 480
GIB_IN_KIB = 1024 * 1024


def terminate_process_group(process: subprocess.Popen[object]) -> None:
    """Stop mpiexec and every rank without leaving an orphaned AMPS process."""
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=30.0)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def latest_analysis_diagnostics(path: Path) -> list[dict[str, str]]:
    """Return only the last base-block evaluation from the append-only table."""
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return []
    last_block = rows[-1]["base_blocks"]
    return [row for row in rows if row["base_blocks"] == last_block]


def run(arguments: argparse.Namespace) -> int:
    """Execute the pilot with operational guards but no scientific step cap."""
    output = arguments.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(
            f"output directory is not empty: {output}; pilot evidence is "
            "never overwritten")
    output.mkdir(parents=True, exist_ok=True)
    arguments._output_claimed = True

    source_input = arguments.input.resolve()
    levels = qualification_levels(source_input)
    if PILOT_RESOLUTION_KM not in levels:
        raise ValueError(
            f"the reviewed pilot requires {PILOT_RESOLUTION_KM} km in the "
            f"approved sequence, got {levels}")
    level = output / f"{PILOT_RESOLUTION_KM}km"
    level.mkdir()
    derived_input = level / "application.in"
    derive_i11_input(source_input, derived_input, PILOT_RESOLUTION_KM)
    runtime = level / "runtime"
    runtime.mkdir()
    receipt = level / "receipt.txt"
    stdout_path = level / "stdout.log"
    stderr_path = level / "stderr.log"
    resource_path = level / "process_resource_metrics.csv"

    command = [
        "mpiexec", "-n", str(arguments.mpi_ranks),
        str(arguments.amps.resolve()),
        "--moon-integration-test", "I11",
        "--integration-output", str(receipt.resolve()),
        "--i11-duration-pilot",
        "-input", str(derived_input.resolve()),
    ]
    run_manifest = manifest(
        "I11", arguments.amps.resolve(), command, 731993,
        derived_input, arguments.mpi_ranks)
    run_manifest.update({
        "schema": "srcMoon-I11-duration-pilot-manifest-v1",
        "evidence_class": "pilot-only-not-I11-qualification",
        "formal_i11_status": "NOT EVALUATED",
        "pilot_resolution_km": PILOT_RESOLUTION_KM,
        "scientific_iteration_limit": None,
        "wall_budget_seconds": arguments.wall_seconds,
        "cooperative_stop_grace_seconds": arguments.stop_grace_seconds,
        "resource_guards": {
            "process_tree_rss_gib": arguments.memory_limit_gib,
            "minimum_host_available_gib": arguments.minimum_host_available_gib,
            "minimum_filesystem_free_gib": arguments.minimum_free_disk_gib,
        },
    })
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")

    # The duration is interpretable only after production storage and actual
    # neutral-Na chemistry/weighted-event accounting have passed on this exact
    # executable/input/rank layout. Preserve every prerequisite as its own
    # linked result; these PASSes cannot become formal I11 convergence.
    for runner, name in (("run_sampling_check.py", "sampling_check"),
                         ("run_budget_check.py", "budget_check")):
        preflight = [sys.executable, "-B", str(Path(__file__).with_name(runner)),
            "--amps", str(arguments.amps.resolve()), "--input", str(source_input),
            "--mpi-ranks", str(arguments.mpi_ranks), "--output-dir", str(output/name)]
        with (output/f"{name}.log").open("w") as log:
            code = subprocess.run(preflight, stdout=log, stderr=subprocess.STDOUT,
                                  check=False).returncode
        if code:
            (output/"pilot_result.json").write_text(json.dumps({
                "pilot_status": "ERROR", "formal_i11_status": "NOT EVALUATED",
                "reason": f"{name} prerequisite did not pass", "returncode": code,
                "artifacts": {"manifest": str(manifest_path), "prerequisite": str(output/name)},
            }, indent=2) + "\n")
            return 2

    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    guard_event: dict[str, object] | None = None
    wall_stop_time = None
    last_progress_minute = -1
    peak_rss_kib = 0
    minimum_available_kib: int | None = None
    minimum_disk_free_bytes: int | None = None
    started = time.monotonic()
    with stdout_path.open("w", encoding="utf-8") as stdout, \
            stderr_path.open("w", encoding="utf-8") as stderr, \
            resource_path.open("w", encoding="utf-8") as resources:
        resources.write(
            "wall_seconds,process_tree_rss_kib,host_mem_available_kib,"
            "filesystem_free_bytes\n")
        # A new session makes the launcher the process-group leader.  The
        # resource guard can therefore terminate mpiexec and every MPI rank as
        # one operation instead of leaking an expensive orphan calculation.
        process = subprocess.Popen(
            command, cwd=runtime, env=environment, stdout=stdout,
            stderr=stderr, start_new_session=True)
        try:
            while process.poll() is None:
                rss_kib = sum(
                    proc_rss_kib(pid) for pid in descendant_pids(process.pid))
                available_kib = host_available_kib()
                disk_free_bytes = shutil.disk_usage(output).free
                elapsed = time.monotonic() - started
                peak_rss_kib = max(peak_rss_kib, rss_kib)
                if available_kib >= 0:
                    minimum_available_kib = (
                        available_kib if minimum_available_kib is None else
                        min(minimum_available_kib, available_kib))
                minimum_disk_free_bytes = (
                    disk_free_bytes if minimum_disk_free_bytes is None else
                    min(minimum_disk_free_bytes, disk_free_bytes))
                resources.write(
                    f"{elapsed:.6f},{rss_kib},{available_kib},"
                    f"{disk_free_bytes}\n")
                resources.flush()
                progress_minute = int(elapsed // 60)
                if progress_minute != last_progress_minute:
                    print(f"RUNNING I11 diagnostic pilot: wall={elapsed:.0f}s "
                          f"RSS={rss_kib / 1024**2:.2f} GiB; model log={stdout_path}",
                          flush=True)
                    last_progress_minute = progress_minute

                reason = None
                if rss_kib >= arguments.memory_limit_gib * GIB_IN_KIB:
                    reason = (
                        "process-tree resident memory reached the approved "
                        f"{arguments.memory_limit_gib:g} GiB pilot guard")
                elif (available_kib >= 0 and available_kib <=
                      arguments.minimum_host_available_gib * GIB_IN_KIB):
                    reason = (
                        "host available memory reached the protective "
                        f"{arguments.minimum_host_available_gib:g} GiB floor")
                elif disk_free_bytes <= int(
                        arguments.minimum_free_disk_gib * 1024**3):
                    reason = (
                        "filesystem free space reached the protective "
                        f"{arguments.minimum_free_disk_gib:g} GiB floor")
                if reason is not None:
                    guard_event = {
                        "reason": reason,
                        "wall_seconds": elapsed,
                        "process_tree_rss_kib": rss_kib,
                        "host_mem_available_kib": available_kib,
                        "filesystem_free_bytes": disk_free_bytes,
                        "action": (
                            "terminated the complete MPI process group; increase "
                            "the affected resource allocation or review a lower-"
                            "cost approved mesh/output contract before rerunning"),
                    }
                    terminate_process_group(process)
                    break
                if elapsed >= arguments.wall_seconds and wall_stop_time is None:
                    # C++ checks this run-owned request between full production
                    # steps and writes a final ledger/receipt. The budget is
                    # operational; it never changes the formal 5000-step cap.
                    (level/"pilot_stop.request").write_text(
                        f"wall-time budget {arguments.wall_seconds:g} s reached\n")
                    wall_stop_time = time.monotonic()
                if (wall_stop_time is not None and time.monotonic()-wall_stop_time >
                        arguments.stop_grace_seconds):
                    guard_event = {"reason": "pilot did not honor cooperative stop within grace period",
                                   "wall_seconds": elapsed,
                                   "action": "terminated the owned MPI process group"}
                    terminate_process_group(process)
                    break
                time.sleep(1.0)
        except KeyboardInterrupt:
            terminate_process_group(process)
            guard_event = {
                "reason": "operator interrupt",
                "wall_seconds": time.monotonic() - started,
                "action": "terminated the complete MPI process group",
            }
        returncode = process.returncode

    common = {
        "schema": "srcMoon-I11-duration-pilot-result-v1",
        "evidence_class": "pilot-only-not-I11-qualification",
        "formal_i11_status": "NOT EVALUATED",
        "source_revision": git_value("rev-parse", "HEAD"),
        "resolution_km": PILOT_RESOLUTION_KM,
        "mpi_ranks": arguments.mpi_ranks,
        "wall_seconds": time.monotonic() - started,
        "peak_process_tree_rss_kib": peak_rss_kib,
        "minimum_host_available_kib": minimum_available_kib,
        "minimum_filesystem_free_bytes": minimum_disk_free_bytes,
        "artifacts": {
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "resource_metrics": str(resource_path),
            "manifest": str(manifest_path),
            "input_snapshot": str(derived_input),
        },
    }
    result_path = output / "pilot_result.json"
    if guard_event is not None:
        guard_path = output / "resource_guard_termination.json"
        guard_path.write_text(
            json.dumps(guard_event, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        common.update({
            "pilot_status": "ERROR", "termination_type": "RESOURCE_LIMIT",
            "reason": guard_event["reason"],
        })
        common["artifacts"]["resource_guard_termination"] = str(guard_path)
        result_path.write_text(
            json.dumps(common, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        print(f"RESOURCE_LIMIT I11 duration pilot: {guard_event['reason']}")
        return 3
    if returncode != 0 or not receipt.is_file():
        common.update({
            "pilot_status": "ERROR",
            "reason": (
                f"linked I11 pilot exited {returncode}; see {stdout_path} and "
                f"{stderr_path}"),
        })
        result_path.write_text(
            json.dumps(common, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        print(f"ERROR I11 duration pilot: linked executable exited {returncode}")
        return 2

    parsed = parse_receipt(receipt)
    required = {
        "test_id": "I11",
        "i11_execution_contract": "duration-pilot-not-qualified",
        "iteration_limit_enabled": 0.0,
        "maximum_particle_steps": 0.0,
        "particle_restart_writes_enabled": 0.0,
        "production_volume_output_suppressed": 1.0,
    }
    disagreements = {
        key: {"expected": expected, "actual": parsed.get(key)}
        for key, expected in required.items() if parsed.get(key) != expected
    }
    qualified = level / "qualified_metrics.csv"
    analysis = level / "analysis_diagnostics.csv"
    blocks = level / "block_metrics.csv"
    converged = parsed.get("convergence_criteria_met") == 1
    if disagreements or not all(
            path.is_file() for path in (analysis, blocks, level/"particle_budget.csv")) or (
                converged and not qualified.is_file()):
        common.update({
            "pilot_status": "ERROR",
            "reason": "pilot receipt/artifact contract was not satisfied",
            "receipt_disagreements": disagreements,
        })
        result_path.write_text(
            json.dumps(common, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        print("ERROR I11 duration pilot: receipt/artifact contract failed")
        return 2

    from run_budget_check import verify_ledger
    from analyze_diagnostics import produce
    budget_check = verify_ledger(level/"particle_budget.csv")
    if budget_check["status"] != "PASS" or parsed.get("particle_budget_valid") != 1:
        common.update(pilot_status="FAIL", reason="production event ledger does not close",
                      budget_check=budget_check)
        result_path.write_text(json.dumps(common, indent=2, sort_keys=True) + "\n")
        return 1
    base_blocks = int(parsed["base_blocks_executed"])
    diagnostics = produce(level, level/"diagnostics") if base_blocks else None
    common["budget_check"] = budget_check
    common["artifacts"]["particle_budget"] = str(level/"particle_budget.csv")
    common["artifacts"]["particle_budget_sha256"] = sha256(level/"particle_budget.csv")
    common["artifacts"]["particle_sources_sha256"] = sha256(level/"particle_sources.csv")
    common["independent_statistics"] = diagnostics
    if not converged:
        if parsed.get("pilot_stop_requested") != 1 or wall_stop_time is None:
            raise ValueError("pilot stopped without statistical completion or declared wall-time request")
        # A bounded diagnostic can finish normally without a candidate. Its
        # statistical evidence remains partial and cannot open any milestone.
        if base_blocks:
            diagnostic = render_i11_stage_diagnostic(output, [PILOT_RESOLUTION_KM], "PILOT")
            common["artifacts"]["stage_diagnostic"] = diagnostic["artifacts"]
        common.update(pilot_status="PARTIAL", reason="wall budget reached before unchanged statistical gates",
                      particle_steps=int(parsed["production_particle_steps_executed"]),
                      base_blocks=base_blocks, candidate_physical_duration=None)
        common["artifacts"]["receipt_sha256"] = sha256(receipt)
        result_path.write_text(json.dumps(common, indent=2, sort_keys=True) + "\n")
        print("PARTIAL I11 duration pilot: bounded diagnostic retained; formal I11 NOT EVALUATED")
        return 0

    steps = int(parsed["production_particle_steps_executed"])
    dt_s = float(parsed["particle_global_timestep_s"])
    sample_steps = int(parsed["base_block_length_steps"])
    batch_blocks = int(parsed["retained_batch_size_base_blocks"])
    batch_count = int(parsed["retained_independent_batch_count"])
    base_blocks = int(parsed["base_blocks_executed"])
    # Analyze() uses the most recent complete set of batches and discards an
    # initial remainder smaller than one batch.  These exact indices expose a
    # reviewable candidate spin-up/science window; the runner does not freeze
    # or apply it to finer grids automatically.
    discarded_blocks = base_blocks - batch_count * batch_blocks
    equilibrium_start_block = discarded_blocks + (batch_count - 9) * batch_blocks
    science_start_block = discarded_blocks + (batch_count - 5) * batch_blocks
    diagnostic = render_i11_stage_diagnostic(output, [PILOT_RESOLUTION_KM], "PILOT")
    common.update({
        "pilot_status": "COMPLETE",
        "reason": (
            "unchanged production independence, equilibrium, and confidence "
            "criteria passed at the 480 km pilot grid; investigator review is "
            "required before freezing a common-duration four-grid contract"),
        "candidate_physical_duration": {
            "particle_timestep_s": dt_s,
            "particle_steps": steps,
            "total_duration_s": steps * dt_s,
            "sample_length_steps": sample_steps,
            "sample_duration_s": sample_steps * dt_s,
            "retained_batch_size_samples": batch_blocks,
            "retained_batch_duration_s": batch_blocks * sample_steps * dt_s,
            "equilibrium_window_start_step":
                equilibrium_start_block * sample_steps,
            "equilibrium_window_start_s":
                equilibrium_start_block * sample_steps * dt_s,
            "science_window_start_step": science_start_block * sample_steps,
            "science_window_start_s": science_start_block * sample_steps * dt_s,
        },
        "latest_criterion_diagnostics": latest_analysis_diagnostics(analysis),
    })
    common["artifacts"].update({
        "receipt": str(receipt),
        "receipt_sha256": sha256(receipt),
        "block_metrics": str(blocks),
        "block_metrics_sha256": sha256(blocks),
        "analysis_diagnostics": str(analysis),
        "analysis_diagnostics_sha256": sha256(analysis),
        "qualified_metrics": str(qualified),
        "qualified_metrics_sha256": sha256(qualified),
        "stage_diagnostic": diagnostic["artifacts"],
    })
    result_path.write_text(
        json.dumps(common, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    print(
        "COMPLETE I11 duration pilot: candidate total physical duration "
        f"{steps * dt_s:.6g} s ({steps} particle steps); formal I11 was not "
        "evaluated")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--wall-seconds", type=float, default=3600.0,
                        help="request a safe post-step diagnostic stop at this wall budget (default: 3600 s)")
    parser.add_argument("--stop-grace-seconds", type=float, default=30.0)
    parser.add_argument(
        "--memory-limit-gib", type=float, default=200.0,
        help="terminate when process-tree RSS reaches this guard (default: 200)")
    parser.add_argument(
        "--minimum-host-available-gib", type=float, default=8.0,
        help="terminate before host available memory falls below this floor")
    parser.add_argument(
        "--minimum-free-disk-gib", type=float, default=100.0,
        help="terminate before the output filesystem falls below this floor")
    arguments = parser.parse_args()
    if arguments.mpi_ranks < 1:
        parser.error("--mpi-ranks must be positive")
    for name in ("memory_limit_gib", "minimum_host_available_gib",
                 "minimum_free_disk_gib", "wall_seconds", "stop_grace_seconds"):
        if not math.isfinite(getattr(arguments, name)) or getattr(arguments, name) <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    try:
        return run(arguments)
    except (OSError, ValueError, RuntimeError, KeyError) as exception:
        result = arguments.output_dir.resolve()/"pilot_result.json"
        if getattr(arguments, "_output_claimed", False) and not result.exists():
            result.write_text(json.dumps({"pilot_status": "ERROR",
                "formal_i11_status": "NOT EVALUATED", "reason": str(exception)},
                indent=2) + "\n")
        print(f"ERROR I11 duration pilot: {exception}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
