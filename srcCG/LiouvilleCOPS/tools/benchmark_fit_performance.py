#!/usr/bin/env python3
"""Characterize Stage-A7 fitter cost on an actual static-kernel base.

The public mode launches one isolated worker for each requested enhancement
degree.  Isolation is important: Linux ``ru_maxrss`` is a process-lifetime
high-water mark, so reusing one process would make the degree-5--8 memory
figures indistinguishable.  Each worker calls the production ``fit_cops``
cache readers, ``Model``, ``Problem.residuals``, and ``Problem.jacobian``.

By default the driver also runs one deliberately truncated (``max_nfev=1``)
production fit per degree in a temporary directory.  Those runs measure the
complete A5 output package; they are output-size probes, not fitted results,
and are deleted after their sizes are recorded.  No flight-data product is
written unless ``--json`` or ``--markdown`` is explicitly requested.

The campaign estimate is a transparent capacity-planning upper scenario.  It
assumes one residual and one Jacobian call per requested optimizer evaluation;
actual trust-region call counts and convergence cannot be inferred from a
microbenchmark.  See ``A7_PERFORMANCE.md`` for interpretation.
"""

import argparse
from argparse import Namespace
import datetime as dt
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import tempfile
import time

from cops_paths import prepare_output, product_input


APP = Path(__file__).resolve().parents[1]
SCHEMA = "liouville_cops.a7_performance"
SCHEMA_VERSION = 1
DEFAULT_DEGREES = (5, 6, 7, 8)
DEFAULT_JOBS = (1, 2, 4, 8, 14)
NUMERICAL_THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
_PROCESS_STARTED = time.perf_counter()
_PARALLEL_PROBLEM = None


def sha256(path):
    """Return a streaming identity for a benchmarked input or implementation."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile_median(samples):
    """Median without importing NumPy in the orchestration process."""
    ordered = sorted(samples)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def current_rss_kib(pid=None):
    """Read current resident memory from procfs; return None off Linux."""
    status = Path("/proc/%d/status" % (pid or os.getpid()))
    try:
        for line in status.read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    except OSError:
        return None
    return None


def available_memory_mib():
    """Read MemAvailable, which is more useful than MemFree for scheduling."""
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return None


def affinity_cpus():
    """Return schedulable CPUs, honoring a scheduler/cgroup affinity mask."""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def directory_inventory(path):
    """Measure every regular file without following directory symlinks."""
    files = {}
    for candidate in sorted(Path(path).rglob("*")):
        if candidate.is_file():
            files[str(candidate.relative_to(path))] = candidate.stat().st_size
    return {"file_count": len(files), "total_bytes": sum(files.values()),
            "files": files}


def base_inventory(run):
    """Measure immutable cache storage separately from per-fit output."""
    elements = Path(run) / "elements.bin"
    kernels = directory_inventory(Path(run) / "kernels")
    table = Path(run) / "COPS_model.dat"
    return {
        "model_table_bytes": table.stat().st_size,
        "elements_bytes": elements.stat().st_size,
        "kernel_count": kernels["file_count"],
        "kernel_bytes": kernels["total_bytes"],
        "total_core_cache_bytes": (table.stat().st_size + elements.stat().st_size
                                   + kernels["total_bytes"]),
    }


def _array_bytes(objects):
    """Estimate live NumPy storage, counting aliased buffers only once.

    This is diagnostic rather than a replacement for RSS: Python, SciPy sparse
    temporaries, allocator retention, shared libraries, and memory-mapped pages
    are visible in RSS but not in this accounting.
    """
    import numpy as np

    seen_objects, buffers = set(), {}

    def visit(value):
        identity = id(value)
        if identity in seen_objects:
            return
        seen_objects.add(identity)
        if isinstance(value, np.ndarray):
            pointer = int(value.__array_interface__["data"][0])
            buffers[pointer] = max(buffers.get(pointer, 0), int(value.nbytes))
        elif isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
        elif hasattr(value, "__dict__"):
            for item in vars(value).values():
                visit(item)

    for obj in objects:
        visit(obj)
    return sum(buffers.values())


def _parallel_evaluation(repetitions):
    """One fork-pool task: production residual+Jacobian calls at the prior."""
    problem = _PARALLEL_PROBLEM
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    checksum = 0.0
    for _ in range(repetitions):
        residual = problem.residuals(problem.x0)
        jacobian = problem.jacobian(problem.x0)
        checksum += float(residual[0]) + float(jacobian[0, 0])
    return {
        "wall_seconds": time.perf_counter() - started_wall,
        "cpu_seconds": time.process_time() - started_cpu,
        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "checksum": checksum,
    }


def _apply_profile(fit_cops, run, header, rows, columns, profile_name,
                   hk_table, source_root):
    """Apply the same named-profile path used by the production fitter."""
    from cops_selection import (apply_selection_profile, load_profile,
                                profile_directory)

    _, intervals = fit_cops.run_intervals(header, source_root, run)
    saturation = fit_cops.ng_saturated(
        hk_table, [row[columns["UTC"]] for row in rows])
    profile, _ = load_profile(profile_name, profile_directory(source_root))
    selected_rows, selection = apply_selection_profile(
        rows, columns, profile, intervals, saturation)
    return selected_rows, selection


def worker(args):
    """Measure one degree in a fresh process and atomically write its JSON."""
    # Numerical libraries are imported only here.  The public driver sets all
    # thread variables in the worker environment before this import occurs.
    import numpy as np
    import scipy
    import fit_cops
    from activity_basis import evaluate_cartesian
    from cops_io import validate_fit_base

    import_finished = time.perf_counter()
    source_root = str(APP)
    run = str(Path(args.run).resolve())

    started = time.perf_counter()
    audit = validate_fit_base(run)
    cache_validation_seconds = time.perf_counter() - started
    if audit["legacy_manifest"] and not args.allow_legacy_base:
        raise RuntimeError("benchmark refuses an unmanifested base without --allow-legacy-base")
    header, rows, columns = (audit[key] for key in
                             ("headers", "rows", "columns"))

    started = time.perf_counter()
    coefficients = fit_cops.read_activity_coefficients(source_root, header)
    elements = fit_cops.read_elements(str(Path(run) / "elements.bin"))
    element_load_seconds = time.perf_counter() - started

    # Retain the full model exactly as fit_cops.py does: it is needed later to
    # write predictions for all points, not only those selected by a profile.
    started = time.perf_counter()
    full_model = fit_cops.Model(run, header, columns, rows, elements,
                                coefficients)
    full_model_load_seconds = time.perf_counter() - started
    started = time.perf_counter()
    p_ng, p_rg = full_model.evaluate(coefficients[0], full_model.T0, 1.0)
    valid = (full_model.pNG_run > 0) & (full_model.pRG_run > 0)
    consistency = {
        "points": int(np.count_nonzero(valid)),
        "max_relative_ng": float(np.max(np.abs(
            p_ng[valid] / full_model.pNG_run[valid] - 1.0))),
        "max_relative_rg": float(np.max(np.abs(
            p_rg[valid] / full_model.pRG_run[valid] - 1.0))),
    }
    consistency_seconds = time.perf_counter() - started

    selected_rows, selection = _apply_profile(
        fit_cops, run, header, rows, columns, args.profile,
        args.hk_table, source_root)
    started = time.perf_counter()
    selected_model = fit_cops.Model(
        run, header, columns, selected_rows, elements, coefficients)
    selected_model.load_kernels(selected_model.sel)
    selected_model_load_seconds = time.perf_counter() - started

    # Time the pure basis constructor independently.  Problem construction
    # below evaluates it again because that is the actual production path.
    basis_samples = []
    basis_bytes = None
    for _ in range(args.repeat):
        started = time.perf_counter()
        basis = evaluate_cartesian(elements["centroid"], 5, args.degree)
        # Reductions force allocation/materialization before stopping the clock.
        float(np.sum(basis, dtype=np.float64))
        basis_samples.append(time.perf_counter() - started)
        basis_bytes = int(basis.nbytes)
        del basis

    body_frame = "67P/C-G_CK"
    if not audit["legacy_manifest"]:
        frames = audit["manifest"].get("effective_settings", {}).get(
            "BodyFrame", [])
        if len(frames) == 1 and frames[0]:
            body_frame = frames[0]
    problem_args = Namespace(
        enhancement_degree=args.degree, fix_legacy_activity=True,
        enhancement_scale=None, degree=4, ordered_temperature=False,
        fix_temperature=True, fix_kappa=True, tbounds=(30.0, 400.0),
        body_frame=body_frame, max_nfev=1000, ftol=1e-10, xtol=1e-10,
        gtol=1e-10, jobs=1, profile=args.profile, data=None,
        apply_exclude_ng=False, exclude_ng_saturated=False)
    started = time.perf_counter()
    problem = fit_cops.Problem(
        selected_model, elements, coefficients, problem_args)
    problem.lam = args.benchmark_lambda
    problem_construction_seconds = time.perf_counter() - started
    ready_seconds = time.perf_counter() - _PROCESS_STARTED

    # Warm caches before collecting repeated operation timings.  The startup
    # timings above remain cold-process measurements.
    problem.residuals(problem.x0)
    problem.jacobian(problem.x0)
    residual_samples, jacobian_samples = [], []
    for _ in range(args.repeat):
        started = time.perf_counter()
        problem.residuals(problem.x0)
        residual_samples.append(time.perf_counter() - started)
        started = time.perf_counter()
        problem.jacobian(problem.x0)
        jacobian_samples.append(time.perf_counter() - started)

    global _PARALLEL_PROBLEM
    scaling = []
    if args.parallel_probe:
        _PARALLEL_PROBLEM = problem
        task_count = max(args.jobs)
        baseline_throughput = None
        for requested in args.jobs:
            effective = min(requested, task_count)
            started = time.perf_counter()
            if effective == 1:
                task_results = [_parallel_evaluation(args.parallel_repeats)
                                for _ in range(task_count)]
            else:
                context = multiprocessing.get_context("fork")
                with context.Pool(effective) as pool:
                    task_results = pool.map(
                        _parallel_evaluation,
                        [args.parallel_repeats] * task_count, chunksize=1)
            wall = time.perf_counter() - started
            evaluations = task_count * args.parallel_repeats
            throughput = evaluations / wall
            if baseline_throughput is None:
                baseline_throughput = throughput
            speedup = throughput / baseline_throughput
            scaling.append({
                "jobs_requested": requested,
                "jobs_effective": effective,
                "tasks": task_count,
                "evaluations_per_task": args.parallel_repeats,
                "wall_seconds": wall,
                "worker_cpu_seconds": sum(item["cpu_seconds"]
                                           for item in task_results),
                "evaluations_per_second": throughput,
                "throughput_speedup": speedup,
                "parallel_efficiency": speedup / effective,
                "maximum_worker_rss_mib": max(
                    item["max_rss_kib"] for item in task_results) / 1024.0,
                "checksum": sum(item["checksum"] for item in task_results),
            })

    result = {
        "degree": args.degree,
        "coefficient_count": int(problem.nb),
        "faces": int(len(elements["area"])),
        "all_points": int(full_model.npt),
        "selected_kernel_points": int(selected_model.npt),
        "kernel_contributions": int(len(selected_model.face)),
        "selection_counts": selection["counts"],
        "timing_seconds": {
            "python_and_module_import": import_finished - _PROCESS_STARTED,
            "cache_validation": cache_validation_seconds,
            "element_and_coefficient_load": element_load_seconds,
            "full_model_kernel_load": full_model_load_seconds,
            "base_consistency_evaluation": consistency_seconds,
            "selected_model_kernel_load": selected_model_load_seconds,
            "problem_construction_including_basis": problem_construction_seconds,
            "ready_to_evaluate": ready_seconds,
            "basis_samples": basis_samples,
            "basis_median": percentile_median(basis_samples),
            "residual_samples": residual_samples,
            "residual_median": percentile_median(residual_samples),
            "jacobian_samples": jacobian_samples,
            "jacobian_median": percentile_median(jacobian_samples),
        },
        "memory": {
            "basis_bytes": basis_bytes,
            "estimated_live_array_bytes": _array_bytes(
                (elements, full_model, selected_model, problem)),
            "current_rss_mib": ((current_rss_kib() or 0) / 1024.0),
            "peak_rss_mib": (resource.getrusage(
                resource.RUSAGE_SELF).ru_maxrss / 1024.0),
        },
        "base_consistency": consistency,
        "parallel_scaling": scaling,
        "environment": {name: os.environ.get(name)
                        for name in NUMERICAL_THREAD_VARIABLES},
        "versions": {"python": platform.python_version(),
                     "numpy": np.__version__, "scipy": scipy.__version__},
    }
    destination = Path(args.worker_result)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    os.replace(str(temporary), str(destination))
    return 0


def run_worker(script, run, degree, args, temporary, parallel_probe=False):
    """Launch one clean degree worker and return its structured measurement."""
    result_path = Path(temporary) / ("degree_%d.json" % degree)
    command = [
        sys.executable, str(script), "--worker", run,
        "--degree", str(degree), "--profile", args.profile,
        "--hk-table", args.hk_table, "--repeat", str(args.repeat),
        "--benchmark-lambda", "%.17g" % args.benchmark_lambda,
        "--worker-result", str(result_path), "--jobs",
    ] + [str(value) for value in args.jobs]
    command += ["--parallel-repeats", str(args.parallel_repeats)]
    if parallel_probe:
        command.append("--parallel-probe")
    if args.allow_legacy_base:
        command.append("--allow-legacy-base")
    environment = os.environ.copy()
    for name in NUMERICAL_THREAD_VARIABLES:
        environment[name] = str(args.threads_per_process)
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=str(APP), env=environment,
                               text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    if completed.returncode:
        raise RuntimeError("degree-%d benchmark failed:\n%s" %
                           (degree, "\n".join(completed.stdout.splitlines()[-80:])))
    result = json.loads(result_path.read_text())
    result["isolated_process_wall_seconds"] = elapsed
    result["worker_command"] = command
    result["worker_output_sha256"] = hashlib.sha256(
        completed.stdout.encode("utf-8")).hexdigest()
    return result


def run_output_probe(run, degree, args, temporary):
    """Run a one-evaluation production fit and measure its complete package."""
    output = Path(temporary) / ("output_degree_%d" % degree)
    log = Path(temporary) / ("output_degree_%d.log" % degree)
    command = [
        sys.executable, str(APP / "tools/fit_cops.py"), run,
        "--profile", args.profile, "--enhancement-degree", str(degree),
        "--fix-legacy-activity", "--fix-temperature", "--fix-kappa",
        "--lambda", "%.17g" % args.benchmark_lambda,
        "--jobs", "1", "--max-nfev", "1", "--out", str(output),
        "--hk-table", args.hk_table,
    ]
    if args.allow_legacy_base:
        command.append("--allow-legacy-base")
    environment = os.environ.copy()
    for name in NUMERICAL_THREAD_VARIABLES:
        environment[name] = str(args.threads_per_process)
    started = time.perf_counter()
    peak_kib = 0
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=str(APP), env=environment,
                                   stdout=stream, stderr=subprocess.STDOUT)
        while process.poll() is None:
            peak_kib = max(peak_kib, current_rss_kib(process.pid) or 0)
            time.sleep(0.05)
    elapsed = time.perf_counter() - started
    text = log.read_text()
    if process.returncode:
        raise RuntimeError("degree-%d output probe failed:\n%s" %
                           (degree, "\n".join(text.splitlines()[-80:])))
    inventory = directory_inventory(output)
    return {
        "method": "production fit_cops.py with max_nfev=1; temporary output deleted",
        "wall_seconds": elapsed,
        "observed_peak_rss_mib": peak_kib / 1024.0,
        "file_count": inventory["file_count"],
        "total_bytes": inventory["total_bytes"],
        "files": inventory["files"],
        "command": command,
        "stdout_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "optimizer_result_is_scientific_fit": False,
    }


def choose_jobs(scaling, cores, available_mib, worker_peak_mib,
                memory_fraction=0.70, minimum_efficiency=0.60):
    """Choose the fastest measured safe candidate under explicit bounds."""
    memory_capacity = (max(1, int(available_mib * memory_fraction /
                                  max(worker_peak_mib, 1.0)))
                       if available_mib is not None else cores)
    allowed = [row for row in scaling
               if row["jobs_effective"] <= cores
               and row["jobs_effective"] <= memory_capacity
               and row["parallel_efficiency"] >= minimum_efficiency]
    if not allowed:
        allowed = [row for row in scaling if row["jobs_effective"] == 1]
    best = max(allowed, key=lambda row: row["evaluations_per_second"])
    return {
        "jobs": best["jobs_effective"],
        "threads_per_process": 1,
        "basis": "fastest measured candidate with efficiency >= %.2f, CPU affinity, and conservative %.0f%% memory bound" %
                 (minimum_efficiency, 100.0 * memory_fraction),
        "available_memory_mib": available_mib,
        "affinity_cpus": cores,
        "conservative_memory_capacity_jobs": memory_capacity,
        "measured_parallel_efficiency": best["parallel_efficiency"],
    }


def campaign_estimate(degree_results, recommendation, profiles, lambdas,
                      max_nfev):
    """Compute an explicitly labeled planning scenario from measured costs."""
    jobs = max(1, recommendation["jobs"])
    efficiency = max(recommendation["measured_parallel_efficiency"], 1e-12)
    rows, serial_total, parallel_total, disk_total = [], 0.0, 0.0, 0
    for item in degree_results:
        timing = item["timing_seconds"]
        evaluation = timing["residual_median"] + timing["jacobian_median"]
        optimizer_cpu = profiles * lambdas * max_nfev * evaluation
        startup = profiles * timing["ready_to_evaluate"]
        serial = startup + optimizer_cpu
        parallel = startup + optimizer_cpu / (jobs * efficiency)
        output_bytes = item.get("output_probe", {}).get("total_bytes")
        disk = profiles * lambdas * output_bytes if output_bytes is not None else None
        serial_total += serial
        parallel_total += parallel
        if disk is not None:
            disk_total += disk
        rows.append({
            "degree": item["degree"], "profiles": profiles,
            "lambda_fits": profiles * lambdas, "max_nfev": max_nfev,
            "assumed_residual_plus_jacobian_seconds": evaluation,
            "serial_cpu_capacity_seconds": serial,
            "estimated_wall_seconds_at_recommended_jobs": parallel,
            "conservative_output_bytes": disk,
        })
    return {
        "scenario": {"profiles": profiles, "lambdas_per_degree": lambdas,
                     "max_nfev_per_fit": max_nfev,
                     "degrees": [item["degree"] for item in degree_results]},
        "assumption": "one residual plus one Jacobian call per max_nfev; single-fit package size multiplied by every profile/lambda result",
        "not_a_convergence_forecast": True,
        "rows": rows,
        "total_serial_cpu_capacity_seconds": serial_total,
        "total_estimated_wall_seconds_at_recommended_jobs": parallel_total,
        "total_conservative_output_bytes": (disk_total if all(
            row["conservative_output_bytes"] is not None for row in rows)
                                             else None),
    }


def bounded_output_names(profiles, degrees, campaign_id):
    """Return deterministic, non-overlapping names suitable for run_fit_profile."""
    safe = "".join(character for character in campaign_id
                   if character.isalnum() or character in "-_")
    if not safe or safe != campaign_id or len(safe) > 40:
        raise ValueError("--campaign-id must be 1..40 letters, digits, '-' or '_'")
    return ["a8_%s_l%d_%s" % (profile, degree, safe)
            for profile in profiles for degree in degrees]


def format_bytes(value):
    if value is None:
        return "not measured"
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return "%.3g %s" % (value, unit)
        value /= 1024.0


def render_markdown(report):
    """Render the machine-readable report as the required A7 resource table."""
    lines = [
        "# LiouvilleCOPS Stage-A7 performance report", "",
        "Generated: `%s`" % report["created_utc"], "",
        "Base: `%s`" % report["base"]["path"], "",
        "Profile: `%s`; numerical threads per process: `%d`." %
        (report["configuration"]["profile"],
         report["configuration"]["threads_per_process"]), "",
        "## Measured resource table", "",
        "| Degree | Coefficients | Ready/startup | Basis | Residual | Jacobian | Peak RSS | Output package |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in report["degrees"]:
        timing, memory = item["timing_seconds"], item["memory"]
        output = item.get("output_probe", {})
        lines.append(
            "| {degree} | {count} | {ready:.3f} s | {basis:.4f} s | "
            "{residual:.4f} s | {jacobian:.4f} s | {rss:.1f} MiB | {output} |".format(
                degree=item["degree"], count=item["coefficient_count"],
                ready=timing["ready_to_evaluate"], basis=timing["basis_median"],
                residual=timing["residual_median"],
                jacobian=timing["jacobian_median"], rss=memory["peak_rss_mib"],
                output=format_bytes(output.get("total_bytes"))))
    recommendation = report["recommendation"]
    lines += [
        "", "## Parallel scaling", "",
        "| Jobs | Evaluations/s | Speedup | Efficiency | Maximum worker RSS |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in report["parallel_scaling"]:
        lines.append("| {jobs} | {rate:.3f} | {speed:.3f} | {eff:.3f} | {rss:.1f} MiB |".format(
            jobs=row["jobs_effective"], rate=row["evaluations_per_second"],
            speed=row["throughput_speedup"], eff=row["parallel_efficiency"],
            rss=row["maximum_worker_rss_mib"]))
    lines += [
        "", "## Recommendation", "",
        "Use `--jobs {jobs} --threads-per-process {threads}` for a scan with at least "
        "{jobs} distinct lambda values on this measured host. For fewer lambda values, "
        "the fitter uses only that many workers.".format(
            jobs=recommendation["jobs"],
            threads=recommendation["threads_per_process"]), "",
        recommendation["basis"] + ".", "",
        "Keep `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` for "
        "multi-fit scans. A single lambda remains serial at the fitter process level.", "",
        "## Campaign planning scenario", "",
        "This is a capacity upper scenario, not a convergence or science forecast: "
        + report["campaign_estimate"]["assumption"] + ".", "",
        "Estimated total serial CPU capacity: `%.2f h`; estimated wall time at the "
        "recommended jobs: `%.2f h`; conservative output: `%s`." % (
            report["campaign_estimate"]["total_serial_cpu_capacity_seconds"] / 3600.0,
            report["campaign_estimate"]["total_estimated_wall_seconds_at_recommended_jobs"] / 3600.0,
            format_bytes(report["campaign_estimate"]["total_conservative_output_bytes"])),
        "", "Bounded output names:", "",
    ]
    lines.extend("- `%s`" % name for name in report["bounded_output_names"])
    lines += [
        "", "Every destination must be new. `run_fit_profile.py` and `fit_cops.py` "
        "refuse an existing directory rather than merge or overwrite results.", "",
        "## Limits", "",
        "Timings are host-, filesystem-cache-, library-, and load-dependent. The "
        "one-evaluation output probes are not fitted scientific results. Repeat this "
        "benchmark on the target machine immediately before a large campaign.", "",
    ]
    return "\n".join(lines)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="?",
                        default="output_fit_full_default_base",
                        help="manifested static fitting base")
    parser.add_argument("--profile", default="all_available",
                        choices=("legacy_base_mask", "all_available",
                                 "transition_cut_100m"))
    parser.add_argument("--hk-table", default="data/COPS_HK_20160930.txt")
    parser.add_argument("--degrees", type=int, nargs="+", default=DEFAULT_DEGREES)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--jobs", type=int, nargs="+", default=DEFAULT_JOBS,
                        help="candidate independent-lambda process counts")
    parser.add_argument("--threads-per-process", type=int, default=1)
    parser.add_argument("--parallel-repeats", type=int, default=1)
    parser.add_argument("--benchmark-lambda", type=float, default=1e-4)
    parser.add_argument("--skip-output-probe", action="store_true")
    parser.add_argument("--campaign-profiles", type=int, default=3)
    parser.add_argument("--campaign-lambdas", type=int, default=5)
    parser.add_argument("--campaign-max-nfev", type=int, default=1000)
    parser.add_argument("--campaign-id", default="campaign01")
    parser.add_argument("--json", help="new JSON report path")
    parser.add_argument("--markdown", help="new Markdown resource-table path")
    parser.add_argument("--allow-legacy-base", action="store_true")
    # Private worker controls are documented in code, not public operating docs.
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--degree", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--worker-result", help=argparse.SUPPRESS)
    parser.add_argument("--parallel-probe", action="store_true",
                        help=argparse.SUPPRESS)
    return parser


def validate_arguments(parser, args):
    if args.repeat < 1 or args.parallel_repeats < 1:
        parser.error("--repeat and --parallel-repeats must be positive")
    if args.threads_per_process < 1:
        parser.error("--threads-per-process must be positive")
    if any(degree not in DEFAULT_DEGREES for degree in args.degrees):
        parser.error("--degrees values must be in 5..8")
    if len(set(args.degrees)) != len(args.degrees):
        parser.error("--degrees must not contain duplicates")
    if any(value < 1 for value in args.jobs) or len(set(args.jobs)) != len(args.jobs):
        parser.error("--jobs must contain distinct positive values")
    args.jobs = sorted(args.jobs)
    if args.jobs[0] != 1:
        parser.error("--jobs must include 1 as the scaling baseline")
    if args.benchmark_lambda < 0 or not math.isfinite(args.benchmark_lambda):
        parser.error("--benchmark-lambda must be finite and nonnegative")
    if min(args.campaign_profiles, args.campaign_lambdas,
           args.campaign_max_nfev) < 1:
        parser.error("campaign counts and --campaign-max-nfev must be positive")
    try:
        bounded_output_names(
            ("legacy_base_mask", "all_available", "transition_cut_100m"),
            args.degrees, args.campaign_id)
    except ValueError as exc:
        parser.error(str(exc))


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_arguments(parser, args)
    if args.worker:
        if args.degree not in DEFAULT_DEGREES or not args.worker_result:
            parser.error("invalid private worker invocation")
        return worker(args)

    run = product_input(args.run.rstrip("/"))
    hk_table = product_input(args.hk_table)
    if not Path(run).is_dir():
        parser.error("fitting base does not exist: %s" % run)
    if not Path(hk_table).is_file():
        parser.error("housekeeping table does not exist: %s" % hk_table)
    args.hk_table = hk_table
    json_target = Path(prepare_output(args.json)) if args.json else None
    markdown_target = Path(prepare_output(args.markdown)) if args.markdown else None
    for target in (json_target, markdown_target):
        if target is not None and target.exists():
            parser.error("report already exists: %s" % target)

    print("A7 base: %s" % run, flush=True)
    print("A7 profile: %s; degrees: %s" %
          (args.profile, " ".join(str(value) for value in args.degrees)), flush=True)
    with tempfile.TemporaryDirectory(prefix="liouville_cops_a7_") as temporary:
        degree_results = []
        for degree in args.degrees:
            print("measuring degree %d ..." % degree, flush=True)
            item = run_worker(
                Path(__file__).resolve(), run, degree, args, temporary,
                parallel_probe=(degree == max(args.degrees)))
            if not args.skip_output_probe:
                print("measuring degree-%d output package ..." % degree,
                      flush=True)
                item["output_probe"] = run_output_probe(
                    run, degree, args, temporary)
            degree_results.append(item)

    scaling = next(item["parallel_scaling"] for item in degree_results
                   if item["parallel_scaling"])
    memory = available_memory_mib()
    cores = affinity_cpus()
    worker_peak = max(row["maximum_worker_rss_mib"] for row in scaling)
    recommendation = choose_jobs(scaling, cores, memory, worker_peak)
    estimates = campaign_estimate(
        degree_results, recommendation, args.campaign_profiles,
        args.campaign_lambdas, args.campaign_max_nfev)
    profiles = ("legacy_base_mask", "all_available", "transition_cut_100m")
    report = {
        "schema": SCHEMA, "schema_version": SCHEMA_VERSION,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scope": "implementation and capacity characterization; not a fitted science result",
        "base": {"path": str(Path(run).resolve()),
                 "model_table_sha256": sha256(Path(run) / "COPS_model.dat"),
                 "inventory": base_inventory(run)},
        "configuration": {
            "profile": args.profile, "degrees": args.degrees,
            "repeat": args.repeat, "candidate_jobs": args.jobs,
            "threads_per_process": args.threads_per_process,
            "benchmark_lambda": args.benchmark_lambda,
            "output_probe": not args.skip_output_probe,
        },
        "host": {"platform": platform.platform(), "affinity_cpus": cores,
                 "available_memory_mib_at_recommendation": memory},
        "implementation": {
            "benchmark_sha256": sha256(Path(__file__)),
            "fit_cops_sha256": sha256(APP / "tools/fit_cops.py"),
            "activity_basis_sha256": sha256(APP / "tools/activity_basis.py"),
        },
        "degrees": degree_results,
        "parallel_scaling": scaling,
        "recommendation": recommendation,
        "campaign_estimate": estimates,
        "bounded_output_names": bounded_output_names(
            profiles, args.degrees, args.campaign_id),
        "overwrite_policy": "all fit/report destinations must be new; existing paths are refused",
        "limitations": [
            "timings depend on host load, filesystem cache, and library versions",
            "max_nfev=1 output probes are deliberately incomplete optimizations",
            "campaign timing assumes one residual and one Jacobian per evaluation",
            "performance measurements do not validate physical adequacy or convergence",
        ],
    }
    markdown = render_markdown(report)
    print(markdown)
    if json_target is not None:
        json_target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print("JSON report: %s" % json_target)
    if markdown_target is not None:
        markdown_target.write_text(markdown + "\n")
        print("Markdown report: %s" % markdown_target)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        sys.exit(130)
