#!/usr/bin/env python3
"""Run CRaTER transport shards concurrently, merge raw moments, and plot.

This local launcher provides the same statistical decomposition as the MPI
environment path in shieldSimCraterBenchmark.  It is convenient on a single
workstation and in CTest; clusters can launch the executable with ``mpiexec``
and call ``merge_crater_shards.py`` afterward.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
import time


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=pathlib.Path)
    parser.add_argument("--case", required=True, type=pathlib.Path)
    parser.add_argument("--workers", required=True, type=int)
    parser.add_argument("--python", default=sys.executable, type=pathlib.Path)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--status-interval-seconds", default=60.0, type=float,
        help="print an alive-worker summary at this interval (default: 60)",
    )
    args = parser.parse_args()
    args.benchmark = args.benchmark.resolve()
    args.case = args.case.resolve()
    require(args.benchmark.is_file(), f"benchmark executable is absent: {args.benchmark}")
    require(args.case.is_file(), f"case is absent: {args.case}")
    require(args.workers >= 2, "parallel launcher requires at least two workers")
    require(args.status_interval_seconds > 0.0,
            "status interval must be positive")
    case = json.loads(args.case.read_text(encoding="utf-8"))
    require(case.get("benchmark_scope") == "foreground_only_axial", "wrong case scope")
    require(args.workers <= int(case["histories"]), "workers cannot exceed histories")
    output = pathlib.Path(case["output"]["directory"])
    if not output.is_absolute():
        output = (args.case.parent / output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    token = re.sub(r"[^A-Za-z0-9_-]", "_", str(case["case_id"]))

    processes: list[tuple[int, subprocess.Popen[str], object]] = []
    started = time.monotonic()
    for index in range(args.workers):
        log_path = output / f"{token}__crater_shard_{index:04d}_of_{args.workers:04d}.log"
        log = log_path.open("w", encoding="utf-8")
        command = [
            str(args.benchmark), str(args.case), f"--shard-index={index}",
            f"--shard-count={args.workers}",
        ]
        processes.append((index, subprocess.Popen(
            command, cwd=output, text=True, stdout=log, stderr=subprocess.STDOUT
        ), log))
    # Poll the whole process set instead of blocking on shard zero.  This does
    # not affect transport or merge ordering; it only makes a multi-hour run
    # auditable from the launcher's stdout.  Shard physics output remains in
    # its dedicated log so eight Geant4 streams cannot interleave text.
    failures: list[int] = []
    finished: set[int] = set()
    next_status = started + args.status_interval_seconds
    while len(finished) != len(processes):
        for index, process, log in processes:
            if index in finished:
                continue
            return_code = process.poll()
            if return_code is None:
                continue
            log.close()
            finished.add(index)
            if return_code != 0:
                failures.append(index)
        now = time.monotonic()
        if len(finished) != len(processes) and now >= next_status:
            print(
                f"elapsed_seconds={now-started:.1f} "
                f"finished_shards={len(finished)}/{len(processes)}",
                flush=True,
            )
            next_status = now + args.status_interval_seconds
        if len(finished) != len(processes):
            time.sleep(0.25)
    require(not failures, f"CRaTER shards failed: {failures}; inspect logs in {output}")

    shards = [
        output / (f"{token}__foreground_only_axial__shard_{index:04d}_of_"
                  f"{args.workers:04d}_partial_result.json")
        for index in range(args.workers)
    ]
    require(all(path.is_file() for path in shards), "one or more shard results are missing")
    merged = output / f"{token}__foreground_only_axial_result.json"
    scripts = pathlib.Path(__file__).resolve().parent
    merge_command = [
        str(args.python), str(scripts / "merge_crater_shards.py"),
        "--output", str(merged), *(str(path) for path in shards),
    ]
    subprocess.run(merge_command, check=True, cwd=output)
    if not args.no_plot:
        prefix = output / f"{token}__foreground_only_axial_comparison"
        subprocess.run([
            str(args.python), str(scripts / "plot_crater_foreground_benchmark.py"),
            str(merged), str(prefix),
        ], check=True, cwd=output)
    elapsed = time.monotonic() - started
    print(f"merged_result={merged}")
    print(f"workers={args.workers} wall_seconds={elapsed:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
