"""Run isolated I11 meshes concurrently with aggregate operational guards.

This module schedules real executables only. Scientific scoring stays in
integration_common.py. A resource or wall-clock stop preserves child artifacts
and returns an execution error; the 500/5000 particle-step contract belongs to
the linked executable and is never changed by this scheduler.
"""
from __future__ import annotations

from contextlib import ExitStack
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time


def stop(process: subprocess.Popen) -> None:
    """Terminate this launcher's private process group, escalating if needed."""
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def execute_levels(jobs: dict, output: Path, workers: int,
                   wall_seconds: float, memory_gib: float = 200.0,
                   render_previews: bool = True) -> dict:
    """Return each launch exit code while recording whole-campaign resources.

    Each job owns its input, CWD, caches and logs. The guards sum RSS over all
    live launch trees; they do not accidentally grant 200 GiB to each child.
    The eight-GiB host and 100-GiB disk floors match the approved pilot guards.
    A common deadline begins at the first launch, independent of child count.
    """
    # Import lazily so integration_common can reuse this scheduler without a
    # module-initialization cycle or a second implementation of /proc parsing.
    from integration_common import (descendant_pids, host_available_kib,
                                    proc_rss_kib, render_i11_stage_diagnostic)
    if workers not in (1, 2, 4):
        raise ValueError("parallel mesh count must be 1, 2, or 4")
    if any(not math.isfinite(value) or value <= 0
           for value in (wall_seconds, memory_gib)):
        raise ValueError("runtime and memory guards must be finite and positive")
    pending = list(jobs)
    active, codes, completed = {}, {}, []
    started = time.monotonic()
    termination = None
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        with ExitStack() as stack:
            totals = stack.enter_context((output / "campaign_resources.csv").open("w"))
            totals.write("wall_seconds,active_meshes,aggregate_rss_kib,"
                         "host_available_kib,filesystem_free_bytes\n")
            while pending or active:
                while pending and len(active) < workers:
                    key = pending.pop(0)
                    job = jobs[key]
                    directory = job["directory"]
                    stdout = stack.enter_context((directory / "stdout.log").open("w"))
                    stderr = stack.enter_context((directory / "stderr.log").open("w"))
                    resources = stack.enter_context(
                        (directory / "process_resource_metrics.csv").open("w"))
                    resources.write("wall_seconds,process_tree_rss_kib,host_mem_available_kib\n")
                    process = subprocess.Popen(job["command"], cwd=directory / "runtime",
                        env=environment, stdout=stdout, stderr=stderr,
                        start_new_session=True)
                    active[key] = (process, resources, time.monotonic())
                aggregate_rss = 0
                available = host_available_kib()
                for key, (process, resources, launched) in list(active.items()):
                    rss = sum(proc_rss_kib(pid) for pid in descendant_pids(process.pid))
                    aggregate_rss += rss
                    resources.write(f"{time.monotonic()-launched:.6f},{rss},{available}\n")
                    resources.flush()
                    code = process.poll()
                    if code is not None:
                        codes[key] = code
                        del active[key]
                        if code == 0:
                            completed.append(key)
                            # Diagnostic plots are available as soon as a mesh
                            # finishes. Formal scoring replaces this provisional
                            # status later; raw means never become a gate.
                            if render_previews:
                                render_i11_stage_diagnostic(output, completed, "PENDING")
                elapsed = time.monotonic() - started
                free = shutil.disk_usage(output).free
                totals.write(f"{elapsed:.6f},{len(active)},{aggregate_rss},{available},{free}\n")
                totals.flush()
                # A completed campaign cannot be retroactively interrupted
                # because rendering its last diagnostic crossed the deadline.
                if not pending and not active:
                    break
                if elapsed >= wall_seconds:
                    termination = f"test-runtime budget of {wall_seconds:g} s exhausted"
                elif aggregate_rss >= memory_gib * 1024**2:
                    termination = f"aggregate process-tree RSS reached {memory_gib:g} GiB"
                elif 0 <= available <= 8 * 1024**2:
                    termination = "host available memory reached 8 GiB floor"
                elif free <= 100 * 1024**3:
                    termination = "filesystem free space reached 100 GiB floor"
                if termination:
                    break
                if active:
                    time.sleep(1)
    except KeyboardInterrupt:
        termination = "operator interrupted campaign"
    finally:
        # Every early return/exception must reap launched children. Preserve
        # completed exit codes while marking interrupted/unlaunched jobs ERROR.
        for key, (process, _, _) in active.items():
            stop(process)
            codes[key] = 124
        for key in pending:
            codes[key] = 124
        if termination:
            (output / "execution_termination.json").write_text(json.dumps({
                "status": "ERROR", "reason": termination,
                "action": "stopped this campaign's MPI launches; preserved completed evidence",
                "remediation": "provide more runtime/memory/disk or reduce concurrent levels; review any science-mesh change separately",
            }, indent=2) + "\n")
    return codes
