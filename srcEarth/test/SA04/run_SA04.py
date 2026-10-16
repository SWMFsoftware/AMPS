#!/usr/bin/env python3
"""Unified prepare/build/run/compare/plot orchestration for physical SA04.

The complete path creates the source-quality-eligible RBSP-B branch, records executable
and dirty-source provenance, requires complete native artifacts, evaluates the fixed
comparison contract, and hashes the report and figures. Archived RBSP-A observations
remain an audit input but cannot become a physical/comparison branch. Stage-only modes
operate on explicit existing evidence. Workflow completion and scientific
qualification remain separate records so a successfully executed disagreement cannot
be promoted to PASS.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import signal
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from compare_SA04 import compare_all, compare_convergence, validate_branch_evidence
from plot_SA04 import plot_all
from prepare_SA04 import prepare_all, add_numerical_arguments, numerical_options
from sa04_common import (
    EXPECTED_OUTPUTS, PROFILES, SPACECRAFT, atomic_json, command_text,
    resolve_data_root, sha256_file, utc_now, validate_spacecraft,
)


def repository_identity(root: Path) -> dict[str, Any]:
    """Record dirty source contents as well as HEAD; HEAD alone misses local fixes."""
    def git(*arguments: str) -> str:
        """Require successful read-only git provenance queries, retaining diagnostics."""
        process = subprocess.run(["git", *arguments], cwd=root, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f"git {' '.join(arguments)} failed: {process.stdout}")
        return process.stdout.strip()
    source_paths = [path for path in (root / "srcEarth").rglob("*")
                    if path.is_file() and path.suffix in (".cpp", ".h", ".py", ".template")]
    source_paths += [root / name for name in ("Makefile", "Makefile.local", "src/pic/pic_mover.cpp")
                     if (root / name).is_file()]
    return {"root": str(root), "head": git("rev-parse", "HEAD"),
            "status_porcelain": git("status", "--porcelain=v1").splitlines(),
            "source_sha256": {str(path.relative_to(root)): sha256_file(path)
                              for path in sorted(source_paths)}}


def reserve_campaign_output(parent: Path, profile: str) -> Path:
    """Reserve a unique evidence directory for a list-driven full campaign.

    A fixed ``--output-dir`` is appropriate for a manually named study, but it is
    intentionally non-overwritable and therefore cannot be placed directly in a
    repeatedly executable repository list.  ``--output-parent`` uses this helper to
    create an empty, exclusive child for each invocation.  The UTC time makes runs
    easy to browse; the PID and bounded collision suffix make allocation race-safe
    when two independent launchers begin during the same second.

    The returned directory is deliberately empty.  ``prepare_all`` accepts an empty
    pre-reserved root but still refuses a nonempty one, preserving the rule that no
    earlier physical evidence may be overwritten.
    """

    resolved_parent = parent.expanduser().resolve()
    resolved_parent.mkdir(parents=True, exist_ok=True)
    stamp = utc_now().replace("-", "").replace(":", "")
    stem = f"sa04_{profile.lower()}_{stamp}_pid{os.getpid()}"
    for collision in range(1000):
        suffix = "" if collision == 0 else f"_{collision:03d}"
        candidate = resolved_parent / f"{stem}{suffix}"
        try:
            candidate.mkdir()
        except FileExistsError:
            continue
        print(f"[SA04] reserved campaign output: {candidate}", flush=True)
        return candidate
    raise RuntimeError(
        f"could not reserve a unique SA04 output below {resolved_parent}"
    )


def safe_rebuild(amps_root: Path) -> dict[str, Any]:
    """Apply the mandated, narrowly checked remove-build then ``make -j`` procedure."""

    expected = Path("/data/vtenishe/GEO/AMPS").resolve()
    root = amps_root.expanduser().resolve()
    if root != expected or not (root / ".git").is_dir() or not (root / "srcEarth").is_dir():
        raise ValueError(f"refusing rebuild for unexpected repository root {root}")
    build = root / "build"
    if build.is_symlink() or build.resolve() != root / "build":
        raise ValueError(f"refusing unsafe build deletion target {build}")
    processes = subprocess.run(
        ["ps", "-eo", "pid=,comm="], text=True, stdout=subprocess.PIPE, check=True
    ).stdout.splitlines()
    # Match actual process names, not a shell/agent command merely mentioning
    # "mpiexec" or "make -j". Conservatively refuse any active launcher/build.
    active = [line.strip() for line in processes if len(line.split()) == 2
              and line.split()[1] in ("amps", "mpiexec", "mpirun", "make", "gmake")
              and int(line.split()[0]) != os.getpid()]
    if active:
        raise RuntimeError("refusing rebuild while possible AMPS/build processes are active: " + "; ".join(active))
    if build.exists():
        shutil.rmtree(build)
    log_path = root / "test_output" / "SA04" / "build_logs" / f"build_{utc_now().replace(':', '')}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    interrupted = False
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        process = subprocess.Popen(["make", "-j"], cwd=root, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert process.stdout is not None
        try:
            for line in process.stdout:
                sys.stdout.write(line); sys.stdout.flush(); log.write(line); log.flush()
            return_code = process.wait()
        except KeyboardInterrupt:
            interrupted = True
            process.terminate()
            try:
                return_code = process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill(); return_code = process.wait()
            log.write("\nRunner interrupted; child launcher terminated.\n")
    record = {"started_utc": started, "finished_utc": utc_now(),
              "command": ["make", "-j"], "return_code": return_code,
              "log": str(log_path)}
    if return_code:
        raise RuntimeError(f"mandatory rebuild failed with exit {return_code}; log={log_path}")
    return record


MODE3D_OUTPUT_PATTERNS = {
    "cutoff": "cutoff_3d_points_snapshot_*.dat",
    "spectrum": "mode3d_points_spectrum_snapshot_*.dat",
    "flux": "mode3d_points_flux_snapshot_*.dat",
    "density": "mode3d_points_density_snapshot_*.dat",
    "termination": "mode3d_termination_summary_snapshot_*.dat",
    "manifest": "standalone_mode3d_manifest_snapshot_*.json",
}


def collect_physical_outputs(branch: Path, prepared: dict[str, Any]) -> tuple[list[str], list[Path]]:
    """Return missing-product diagnostics and the complete current artifact list.

    A GRIDDED SNAPSHOT_LIST run intentionally writes one suffixed product set per
    trajectory epoch.  Requiring the exact selected-epoch count prevents a zero-exit
    partial batch from being accepted as complete.
    """

    solver = str(prepared.get("solver", "GRIDLESS")).upper()
    if solver == "GRIDLESS":
        files = [branch / name for name in EXPECTED_OUTPUTS if (branch / name).is_file()]
        missing = [name for name in EXPECTED_OUTPUTS if not (branch / name).is_file()]
        return missing, files
    if solver != "GRIDDED":
        raise ValueError(f"unsupported prepared solver {solver!r}")

    expected_count = int(prepared["trajectory"]["kept"])
    files: list[Path] = []
    missing: list[str] = []
    for product, pattern in MODE3D_OUTPUT_PATTERNS.items():
        matches = sorted(branch.glob(pattern))
        files.extend(matches)
        if len(matches) != expected_count:
            missing.append(
                f"{product}: expected {expected_count} files matching {pattern}, found {len(matches)}"
            )
    return missing, files


def terminate_process_group(process: subprocess.Popen[Any]) -> int:
    """Terminate only the launcher's dedicated process group and reap it."""

    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return process.wait()
    try:
        return process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return process.wait()


def process_identity(pid: int) -> str | None:
    """Read Linux start ticks, rejecting zombies and PID reuse while waiting.

    Merely testing ``kill(pid, 0)`` can wait forever on a zombie or accidentally
    follow an unrelated replacement process. The launcher's process group is never
    signalled by this observer; ownership of active jobs stays with their runners.
    """
    try:
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
    except FileNotFoundError:
        return None
    return None if fields[0] == "Z" else fields[19]


def wait_for_evidence(root: Path, diagnostic: bool = False) -> dict[str, Any]:
    """Wait for a recorded job to finish, never treating disappearance as success.

    Completed failure records are returned, not suppressed: they allow independent
    subsequent physical work while the failure remains available for reporting.
    A missing terminal record after the recorded runner exits is an error requiring
    inspection. A sixty-second progress interval keeps unattended waits observable.
    """
    root = root.expanduser().resolve()
    terminal = root / ("diagnostic_status.json" if diagnostic else "physical_execution_summary.json")
    active = root / ("diagnostic_status.json" if diagnostic else "active_run.json")
    def finished() -> dict[str, Any] | None:
        """Read only atomically published terminal records, not partial products."""
        if not terminal.is_file():
            return None
        record = json.loads(terminal.read_text())
        return record if record.get("status") != "RUNNING" else None
    complete = finished()
    if complete:
        return complete
    record = json.loads(active.read_text())
    pid = int(record["runner_pid"])
    identity = process_identity(pid)
    if identity is None:
        raise RuntimeError(f"{root}: recorded runner {pid} is not active and has no terminal record")
    next_message = 0.0
    while True:
        complete = finished()
        if complete:
            print(f"[WAIT COMPLETE] {root}: {complete['status']}", flush=True)
            return complete
        if process_identity(pid) != identity:
            # The atomic terminal write precedes process exit. Re-read once to
            # cover the narrow race between the two observations above.
            complete = finished()
            if complete:
                return complete
            raise RuntimeError(f"{root}: runner exited/replaced without terminal evidence")
        if time.monotonic() >= next_message:
            print(f"[WAIT] {utc_now()} runner={pid}, evidence={terminal}", flush=True)
            next_message = time.monotonic() + 60.0
        time.sleep(5.0)


def launch_branch(output: Path, spacecraft: str, amps: Path,
                  mpi_ranks: int, threads: int, identity: dict[str, Any]) -> dict[str, Any]:
    """Execute one immutable branch, preserving partial failures and complete hashes."""
    branch = output / spacecraft
    prepared_path = branch / "prepared_run.json"
    if not prepared_path.is_file():
        raise FileNotFoundError(f"prepare record missing: {prepared_path}")
    prepared = json.loads(prepared_path.read_text(encoding="utf-8"))
    for path_key, hash_key in (("input_file", "input_sha256"),
                               ("driver_file", "driver_sha256"),
                               ("spectrum_file", "spectrum_sha256"),
                               ("snapshot_list_file", "snapshot_list_sha256")):
        if prepared.get(path_key) and sha256_file(Path(prepared[path_key])) != prepared[hash_key]:
            raise ValueError(f"{branch}: prepared {path_key} changed before execution")
    if sha256_file(Path(prepared["trajectory"]["output"])) != prepared["trajectory"]["sha256"]:
        raise ValueError(f"{branch}: prepared trajectory changed before execution")
    _, existing = collect_physical_outputs(branch, prepared)
    if existing:
        raise FileExistsError(
            f"{branch} already contains physical products; select a new run directory to preserve them"
        )
    solver = str(prepared.get("solver", "GRIDLESS")).upper()
    if solver == "GRIDDED":
        mesh = prepared.get("mode3d_mesh") or {}
        command = [
            "mpiexec", "-n", str(mpi_ranks), str(amps), "-mode", "3d",
            "-i", "input/AMPS_PARAM.in", "-mode3d-field-eval", "INTERPOLATION",
            # Execute the explicit policy recorded during SA04 preparation. The
            # application default is COLLECTIVE, so spelling out EPOCH_PARALLEL here
            # keeps command-line and deck provenance aligned and prevents an unrelated
            # SNAPSHOT_LIST workflow from inheriting replicated meshes implicitly.
            "-mode3d-trajectory-snapshot-scheduler", "EPOCH_PARALLEL",
            "-mode3d-epoch-parallel-memory-fraction", "0.80",
            "-density-parallel", "OPENMP", "-density-threads", str(threads),
            "-mode3d-mpi-scheduler", "DYNAMIC", "-mode3d-mpi-dynamic-chunk", "1",
            "-mode3d-mesh-res-earth-re", str(mesh["res_earth_re"]),
            "-mode3d-mesh-res-boundary-re", str(mesh["res_boundary_re"]),
            "-mode3d-mesh-coarsening", str(mesh["coarsening"]),
            "-mode3d-mesh-exponent", str(mesh["exponent"]),
            "-mode3d-parallel-field-init",
        ]
    elif solver == "GRIDLESS":
        command = ["mpiexec", "-n", str(mpi_ranks), str(amps), "-mode", "gridless",
                   "-i", "input/AMPS_PARAM.in", "-gridless-parallel", "THREADS",
                   "-gridless-threads", str(threads), "-gridless-mpi-scheduler", "DYNAMIC"]
    else:
        raise ValueError(f"unsupported prepared solver {solver!r}")
    command += ["-mover", prepared.get("numerical_controls", {}).get("mover", "RK4")]
    environment = dict(os.environ); environment["OMP_NUM_THREADS"] = str(threads)
    started = utc_now(); log_path = branch / "amps.log"
    print(f"[START] {spacecraft} {solver}: {command_text(command)}", flush=True)
    interrupted = False
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        log.write("Command: " + command_text(command) + "\nWorking directory: " + str(branch) + "\n\n")
        log.flush()
        process = subprocess.Popen(
            command, cwd=branch, env=environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, start_new_session=True,
            bufsize=1,
        )
        atomic_json(output / "active_run.json", {
            "schema": "amps-sa04-active-run-v2", "status": "RUNNING",
            "started_utc": started, "spacecraft": spacecraft, "solver": solver,
            "runner_pid": os.getpid(), "launcher_pid": process.pid,
            "command": command, "working_directory": str(branch), "log": str(log_path),
        })
        try:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line); log.flush()
                print(line, end="", flush=True)
            return_code = process.wait()
        except KeyboardInterrupt:
            interrupted = True
            return_code = terminate_process_group(process)
            log.write("\nRunner interrupted; child launcher terminated.\n")
    missing, output_files = collect_physical_outputs(branch, prepared)
    status = {
        **prepared, "schema": "amps-sa04-physical-run-v3",
        "status": "PASS" if return_code == 0 and not missing else "ERROR",
        "started_utc": started, "finished_utc": utc_now(), "command": command,
        "working_directory": str(branch), "return_code": return_code,
        "solver": solver, "requested_mpi_ranks": mpi_ranks,
        "requested_threads_per_rank": threads,
        "effective_shared_memory_backend": "OPENMP" if solver == "GRIDDED" else "THREADS",
        "repository": identity, "executable": str(amps),
        "executable_sha256": sha256_file(amps), "missing_expected_outputs": missing,
        "interrupted": interrupted,
        "output_sha256": {path.name: sha256_file(path) for path in output_files},
    }
    atomic_json(branch / "run_status.json", status)
    atomic_json(output / "active_run.json", {
        "schema": "amps-sa04-active-run-v2", "status": status["status"],
        "started_utc": started, "finished_utc": status["finished_utc"],
        "spacecraft": spacecraft, "solver": solver, "return_code": return_code,
        "interrupted": interrupted, "command": command,
        "working_directory": str(branch), "log": str(log_path),
    })
    print(f"[{'COMPLETE' if status['status'] == 'PASS' else 'ERROR'}] "
          f"{spacecraft} {solver}: exit={return_code}, missing={len(missing)}", flush=True)
    if status["status"] != "PASS":
        reason = "interrupted; " if interrupted else ""
        raise RuntimeError(f"{spacecraft} AMPS {reason}error: exit={return_code}, missing={missing}; log={log_path}")
    return status


def diagnostic_evidence(output: Path, additions: list[Path] | None = None) -> list[dict[str, Any]]:
    """Bind explicitly relevant bounded checks to a campaign without granting coverage.

    Registration is deliberate: a timestep experiment on an older configuration
    must not silently condemn a different campaign. Once registered, a failed check
    cannot disappear merely because a later report invocation omits the CLI option.
    All comparison bytes, native diagnostic records and their recorded inputs/outputs
    are revalidated. A bounded PASS never replaces the full-population requirement.
    Use a new campaign directory when changing the qualified numerical configuration;
    this API intentionally provides no flag for discarding inconvenient evidence.
    """
    from analyze_diagnostics import completed

    path = output / "diagnostic_evidence.json"
    saved = json.loads(path.read_text()) if path.is_file() else {"records": []}
    campaign = {str(prepared.resolve()): sha256_file(prepared)
                for prepared in sorted(output.glob("*/prepared_run.json"))}
    if path.is_file() and campaign != saved["campaign_prepared_sha256"]:
        raise ValueError("registered diagnostics belong to different prepared campaign inputs")
    records = {record["path"]: record for record in saved["records"]}
    for addition in additions or []:
        resolved = addition.expanduser().resolve()
        key = str(resolved)
        digest = sha256_file(resolved)
        if key in records and records[key]["sha256"] != digest:
            raise ValueError(f"registered diagnostic evidence changed: {resolved}")
        records[key] = {"path": key, "sha256": digest}
    results = []
    for record in records.values():
        comparison_path = Path(record["path"])
        if sha256_file(comparison_path) != record["sha256"]:
            raise ValueError(f"registered diagnostic evidence changed: {comparison_path}")
        result = json.loads(comparison_path.read_text())
        if (result.get("check") not in ("exact_repeatability", "bounded_numerical_difference")
                or result.get("status") not in ("PASS", "FAIL", "INCOMPLETE")):
            raise ValueError(f"unsupported bounded diagnostic result: {comparison_path}")
        for label in ("coarse", "fine"):
            current = completed(Path(result[label]))
            if current != result["numerical_controls"][label]:
                raise ValueError(f"stale {label} diagnostic evidence: {comparison_path}")
        for artifact, expected in result["artifact_sha256"].items():
            if sha256_file(Path(artifact)) != expected:
                raise ValueError(f"changed diagnostic product: {artifact}")
        results.append({**record, "status": result["status"], "check": result["check"],
                        "maximum_fractional_change": result["maximum_fractional_change"],
                        "threshold": result["threshold"]})
    if additions:
        if not campaign:
            raise ValueError("diagnostic registration requires prepared campaign inputs")
        atomic_json(path, {"schema": "amps-sa04-registered-diagnostics-v1",
                           "updated_utc": utc_now(), "campaign_prepared_sha256": campaign,
                           "records": list(records.values())})
    return results


def qualification_status(comparison: dict[str, Any], convergence: dict[str, Any] | None,
                         diagnostics: list[dict[str, Any]] | None = None) -> str:
    """A full validation PASS requires the eligible RBSP-B branch and convergence.

    Execution-only and standalone observational comparisons have narrower meanings.
    This decision is shared by the final report and the full runner's exit status,
    preventing a successful observation score from concealing missing refinement.
    The caller must validate saved convergence provenance before using this result.
    """
    if (comparison["status"] == "FAIL" or (convergence and convergence["status"] == "FAIL")
            or any(item["status"] == "FAIL" for item in diagnostics or [])):
        return "FAIL"
    if any(item["status"] != "PASS" for item in diagnostics or []):
        return "INCOMPLETE"
    if (comparison["status"] == "PASS" and convergence and convergence["status"] == "PASS"
            and set(comparison["spacecraft"]) == set(SPACECRAFT)
            and set(convergence["evidence"]) == set(SPACECRAFT)):
        return "PASS"
    return "INCOMPLETE"


def write_report(output: Path, comparison: dict[str, Any], figures: list[str]) -> Path:
    """Separate observed agreement from qualification requiring numerical convergence."""
    convergence_path = output / "convergence_summary.json"
    convergence = json.loads(convergence_path.read_text()) if convergence_path.is_file() else None
    refinement_error_path = output / "refinement_error.json"
    refinement_error = (json.loads(refinement_error_path.read_text())
                        if refinement_error_path.is_file() else None)
    if convergence:
        # A copied/stale PASS must not qualify a newly executed branch. Revalidate
        # both saved populations and the status hashes bound into the comparison;
        # this checks the input/output bytes transitively through those statuses.
        if output.resolve() not in {Path(convergence[key]).resolve()
                                    for key in ("coarse_root", "fine_root")}:
            raise ValueError("convergence report belongs to a different campaign")
        for ship, evidence in convergence["evidence"].items():
            for label in ("coarse", "fine"):
                branch = Path(convergence[label + "_root"]) / ship
                validate_branch_evidence(branch)
                if sha256_file(branch / "run_status.json") != evidence[label + "_status_sha256"]:
                    raise ValueError(f"{ship}: stale {label} convergence evidence")
    diagnostics = diagnostic_evidence(output)
    qualification = qualification_status(comparison, convergence, diagnostics)
    if refinement_error and qualification == "PASS":
        qualification = "INCOMPLETE"
    lines = ["# SA04 validation report", "", f"Created: {utc_now()}", "",
             f"Observational status: **{comparison['status']}**", "",
             f"Overall qualification: **{qualification}**", "",
             f"Numerical convergence: **{convergence['status'] if convergence else 'NOT EVALUATED'}**", "",
             "The REPT density reference is isotropic-equivalent and derived from the same seven differential channels; it is not an independent direct density measurement.",
             "The model uses successive static T05 snapshots with no electric field, so these results do not establish dynamic electric-field transport.", "",
             "## Metrics", "",
             "| Population | Product | N | Median absolute log10 ratio | Factor two |", "|---|---|---:|---:|---:|"]
    for spacecraft, detail in comparison["spacecraft"].items():
        for product, metric in detail["products"].items():
            lines.append(f"| {spacecraft} | {product} | {metric['count']} | {metric['median_abs_log10_ratio']} | {metric['factor_two_fraction']:.4f} |")
    for product, metric in comparison["pooled"].items():
        lines.append(f"| pooled | {product} | {metric['count']} | {metric['median_abs_log10_ratio']} | {metric['factor_two_fraction']:.4f} |")
    lines += ["", "## Observation eligibility", "",
              "Eligibility is determined from immutable REPT channel/density structure before model/reference scoring.", "",
              "| Spacecraft | Modeled epochs | Comparison-eligible epochs | Ineligible modeled epochs |",
              "|---|---:|---:|---:|"]
    eligibility_details: list[str] = []
    for spacecraft, detail in comparison["spacecraft"].items():
        excluded = detail.get("observation_ineligible_modeled_epochs", [])
        modeled = detail.get("modeled_epochs", detail["matched_epochs_by_product"]["spectrum"])
        eligible = detail.get("comparison_eligible_modeled_epochs", modeled - len(excluded))
        lines.append(f"| {spacecraft} | {modeled} | {eligible} | {len(excluded)} |")
        for item in excluded:
            channels = ", ".join(item["missing_channels"]) or "none"
            eligibility_details.append(
                f"- `{spacecraft}/{item['utc']}`: missing channels `{channels}`; "
                f"missing derived density: `{item['missing_density']}`."
            )
    if eligibility_details:
        lines += ["", "Observation-ineligible modeled epochs:", ""] + eligibility_details
    failed = [gate for gate in comparison["gates"] if not gate["passed"]]
    lines += ["", "## Gate result", "", f"Failed gates: {len(failed)} of {len(comparison['gates'])}."]
    for gate in failed:
        lines.append(f"- `{gate['name']}`: {gate['value']} (required {gate['operator']} {gate['threshold']})")
    if convergence:
        lines += ["", "## Numerical convergence", "",
                  f"Maximum integrated-product fractional change: {convergence['maximum_fractional_change']:.6g} (required <= 0.02).",
                  f"Common-population product rows: {convergence['common_population_rows']}.",
                  "See `convergence_summary.json` and `.csv` for controls, provenance, every channel and epoch."]
    if refinement_error:
        lines += ["", "## Refinement evidence error", "",
                  "The requested comparison could not be evaluated; it is not a numerical PASS.",
                  refinement_error["error"]]
    if diagnostics:
        lines += ["", "## Registered bounded diagnostics", "",
                  "These checks do not supply campaign coverage. A registered FAIL remains a qualification failure even if a different refinement axis passes.", ""]
        for item in diagnostics:
            lines.append(f"- {item['check']}: **{item['status']}**, maximum change "
                         f"{item['maximum_fractional_change']:.6g}, limit {item['threshold']}; "
                         f"evidence `{item['path']}`, SHA-256 `{item['sha256']}`.")
    lines += ["", "## Figures", ""] + [f"- `{Path(path).relative_to(output)}`" for path in figures]
    lines += ["", "## Limitations", "",
              "REPT bow-tie intervals are folded as rectangular energy bins, not through a full calibrated angular/energy response.",
              "AMPS samples an instantaneous field/boundary state at each selected epoch; REPT values are five-minute means.",
              "No target-informed scaling, spacecraft-specific fit, background subtraction, or favorable exclusion was applied."]
    path = output / "SA04_validation_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    atomic_json(output / "qualification_summary.json", {
        "schema": "amps-sa04-qualification-v1", "created_utc": utc_now(),
        "status": qualification, "observational_status": comparison["status"],
        "convergence_status": convergence["status"] if convergence else "NOT_EVALUATED",
        "registered_diagnostics": diagnostics, "refinement_error": refinement_error,
        "report_sha256": sha256_file(path), "reporting_source_sha256": sha256_file(Path(__file__)),
    })
    return path


def record_workflow_completion(output: Path, comparison: dict[str, Any],
                               figures: list[str], report: Path,
                               qualification: str, profile: str) -> dict[str, Any]:
    """Record a code-validation PASS without rewriting scientific qualification.

    The repository QUICK entry needs to answer a narrower question than the SA04
    scientific qualification: did the current executable finish every eligible branch,
    produce complete comparison evidence, and render the requested figures?  A
    physically completed workflow can answer that question even when the unchanged
    observational or convergence gates report FAIL.  Keeping this in a separate
    artifact prevents an execution PASS from being mistaken for scientific PASS.

    Every artifact is rechecked here after plotting.  Comparison ``PASS`` and
    ``FAIL`` are both evaluated scientific outcomes; missing/error evidence raises
    instead of becoming workflow success.  The caller writes this record only after
    ``write_report`` has produced the independent qualification summary.
    """

    output = output.expanduser().resolve()
    if comparison.get("status") not in ("PASS", "FAIL"):
        raise ValueError("workflow completion requires an evaluated comparison status")
    physical = output / "physical_execution_summary.json"
    comparison_path = output / "comparison_summary.json"
    plot_manifest = output / "plot_manifest.json"
    qualification_path = output / "qualification_summary.json"
    required = (physical, comparison_path, plot_manifest, qualification_path, report)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "workflow completion is missing required evidence: " + ", ".join(missing)
        )
    physical_record = json.loads(physical.read_text(encoding="utf-8"))
    if physical_record.get("status") != "PASS":
        raise ValueError("workflow completion requires physical execution status PASS")
    if not figures:
        raise ValueError("workflow completion requires at least one comparison figure")
    figure_paths = [Path(path).expanduser().resolve() for path in figures]
    absent_figures = [str(path) for path in figure_paths if not path.is_file()]
    if absent_figures:
        raise FileNotFoundError(
            "workflow completion is missing generated figures: "
            + ", ".join(absent_figures)
        )

    record = {
        "schema": "amps-sa04-workflow-completion-v1",
        "created_utc": utc_now(),
        "status": "PASS",
        "profile": profile,
        "scope": "physical execution, comparison evaluation, report, and plots",
        "scientific_qualification": qualification,
        "observational_status": comparison["status"],
        "qualification_is_unchanged": True,
        "figure_count": len(figure_paths),
        "evidence_sha256": {
            str(path): sha256_file(path) for path in required
        },
        "figure_sha256": {
            str(path): sha256_file(path) for path in figure_paths
        },
    }
    atomic_json(output / "workflow_completion.json", record)
    return record


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source = result.add_mutually_exclusive_group(required=True)
    source.add_argument("--data-root", type=Path)
    source.add_argument("--data-archive", type=Path,
                        help=("Alias naming an extracted package, its model_inputs/SA04 "
                              "or model_inputs directory, or its validation_data "
                              "directory; not a tar file."))
    result.add_argument("--amps-root", type=Path, default=Path("/data/vtenishe/GEO/AMPS"))
    result.add_argument("--amps", type=Path, default=Path("/data/vtenishe/GEO/AMPS/amps"))
    destination = result.add_mutually_exclusive_group(required=True)
    destination.add_argument(
        "--output-dir", type=Path,
        help="Exact campaign directory; must be absent or empty for a full run.",
    )
    destination.add_argument(
        "--output-parent", type=Path,
        help=("Parent under which a unique UTC/PID campaign directory is reserved. "
              "Use this for repeatable test/studies list entries; full mode only."),
    )
    result.add_argument("--profile", choices=tuple(PROFILES), default="STANDARD")
    result.add_argument("--mpi-ranks", type=int, default=4)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--solver", choices=("GRIDLESS", "GRIDDED"), default="GRIDLESS")
    result.add_argument("--mode3d-mesh-res-earth-re", type=float, default=0.05)
    result.add_argument("--mode3d-mesh-res-boundary-re", type=float, default=2.0)
    result.add_argument("--mode3d-mesh-coarsening",
                        choices=("LINEAR", "LOG", "EXPONENTIAL", "POWER", "CONSTANT"),
                        default="LINEAR")
    result.add_argument("--mode3d-mesh-exponent", type=float, default=1.0)
    result.add_argument("--spacecraft", nargs="+", default=list(SPACECRAFT))
    result.add_argument("--rebuild", action="store_true")
    modes = result.add_mutually_exclusive_group()
    modes.add_argument("--prepare-only", action="store_true")
    modes.add_argument("--run-only", action="store_true", help="Run already-prepared branches without replacing them.")
    modes.add_argument("--compare-only", action="store_true")
    modes.add_argument("--plot-only", action="store_true")
    result.add_argument("--formats", nargs="+", choices=("png", "eps"), default=["png", "eps"])
    result.add_argument(
        "--exit-code-policy", choices=("QUALIFICATION", "WORKFLOW"),
        default="QUALIFICATION",
        help=("QUALIFICATION returns success only for full scientific qualification. "
              "WORKFLOW returns success after complete execution/comparison/plots "
              "while preserving the independent scientific PASS/FAIL result."),
    )
    result.add_argument("--diagnostic-evidence", action="append", type=Path, default=[],
                        help="Register a relevant bounded comparison JSON; repeatable, retained for future reports, never substitutes for campaign convergence.")
    result.add_argument("--wait-for-campaign", type=Path,
                        help="Before execution, wait for this recorded campaign runner to finish; missing terminal evidence is an error.")
    result.add_argument("--wait-for-diagnostic", action="append", type=Path, default=[],
                        help="Also wait for this recorded bounded run; repeatable, never changes or terminates it.")
    result.add_argument("--convergence-base", type=Path,
                        help="After this full run, compare against the completed coarse campaign on common epochs and refresh both reports/figures.")
    add_numerical_arguments(result)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.mpi_ranks < 1 or args.threads < 1:
        raise ValueError("--mpi-ranks and --threads must be positive")
    if args.rebuild and (args.compare_only or args.plot_only):
        raise ValueError("--rebuild is incompatible with standalone compare/plot modes")
    if ((args.wait_for_campaign or args.wait_for_diagnostic or args.convergence_base)
            and (args.prepare_only or args.compare_only or args.plot_only or args.run_only)):
        raise ValueError("waiting/refinement continuation requires the full prepare/run/compare/plot mode")
    stage_only = args.prepare_only or args.compare_only or args.plot_only or args.run_only
    if args.output_parent and stage_only:
        raise ValueError(
            "--output-parent is only valid for a full prepare/run/compare/plot "
            "campaign; stage-only modes require the exact existing --output-dir"
        )
    if args.exit_code_policy == "WORKFLOW" and stage_only:
        raise ValueError(
            "--exit-code-policy WORKFLOW is only valid for the complete workflow"
        )
    data_path = args.data_root or args.data_archive
    resolve_data_root(data_path)  # fail early before any output/build mutation
    spacecraft = validate_spacecraft(args.spacecraft)
    root = args.amps_root.expanduser().resolve(); amps = args.amps.expanduser().resolve()
    output = (args.output_dir.expanduser().resolve() if args.output_dir
              else reserve_campaign_output(args.output_parent, args.profile))
    # Wait before preparation, build or executable capture. The currently running
    # processes therefore cannot be overwritten by this sequential continuation.
    if args.wait_for_campaign or args.wait_for_diagnostic:
        # Keep the waiting record beside, not inside, the future campaign root:
        # preparation must still refuse every nonempty output directory. Exclusive
        # creation prevents two continuations from quietly targeting the same root.
        continuation = output.with_name(output.name + ".continuation.json")
        continuation.parent.mkdir(parents=True, exist_ok=True)
        record = {"status": "WAITING", "started_utc": utc_now(), "runner_pid": os.getpid(),
                  "runner_start_ticks": process_identity(os.getpid()),
                  "command": [sys.executable, *sys.argv], "output_dir": str(output),
                  "wait_for_campaign": str(args.wait_for_campaign) if args.wait_for_campaign else None,
                  "wait_for_diagnostic": [str(path) for path in args.wait_for_diagnostic]}
        with continuation.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2)
        try:
            if args.wait_for_campaign:
                wait_for_evidence(args.wait_for_campaign)
            for waiting in args.wait_for_diagnostic:
                wait_for_evidence(waiting, diagnostic=True)
        except (Exception, KeyboardInterrupt) as error:
            atomic_json(continuation, {**record, "status": "ERROR", "finished_utc": utc_now(),
                                       "error": str(error) or type(error).__name__})
            raise
        atomic_json(continuation, {**record, "status": "WAIT_COMPLETE", "finished_utc": utc_now(),
                                   "next_evidence": str(output / "active_run.json")})
    build_record = safe_rebuild(root) if args.rebuild else None
    if not (args.compare_only or args.plot_only or args.prepare_only) and not amps.is_file():
        raise FileNotFoundError(amps)
    identity = repository_identity(root)

    if args.prepare_only:
        prepare_all(data_path, output, args.profile, args.threads, spacecraft,
                    args.solver, args.mode3d_mesh_res_earth_re,
                    args.mode3d_mesh_res_boundary_re,
                    args.mode3d_mesh_coarsening, args.mode3d_mesh_exponent,
                    numerical_options(args))
        diagnostic_evidence(output, args.diagnostic_evidence)
        print(f"[PASS] preparation only: {output}"); return 0
    if args.compare_only:
        diagnostic_evidence(output, args.diagnostic_evidence)
        result = compare_all(data_path, output, spacecraft)
        write_report(output, result, [str(p) for p in sorted((output / "figures").glob("*")) if p.is_file()])
        print(f"[{result['status']}] comparison: {output / 'comparison_summary.json'}")
        return 0 if result["status"] == "PASS" else 2
    if args.plot_only:
        diagnostic_evidence(output, args.diagnostic_evidence)
        figures = plot_all(data_path, output, spacecraft, tuple(args.formats))
        write_report(output, json.loads((output / "comparison_summary.json").read_text()), figures)
        print(f"[PASS] plotting only: {len(figures)} files"); return 0

    if not args.run_only:
        prepare_all(data_path, output, args.profile, args.threads, spacecraft,
                    args.solver, args.mode3d_mesh_res_earth_re,
                    args.mode3d_mesh_res_boundary_re,
                    args.mode3d_mesh_coarsening, args.mode3d_mesh_exponent,
                    numerical_options(args))
    diagnostic_evidence(output, args.diagnostic_evidence)
    atomic_json(output / "execution_environment.json", {
        "created_utc": utc_now(), "python": sys.version, "platform": platform.platform(),
        "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "matplotlib", "pandas")},
        "repository": identity, "executable_sha256": sha256_file(amps), "build": build_record,
    })
    errors = {}
    for item in spacecraft:  # deliberately sequential, failures retain their own evidence
        try:
            launch_branch(output, item, amps, args.mpi_ranks, args.threads, identity)
        except RuntimeError as error:
            errors[item] = str(error)
            if "interrupted" in str(error):
                break  # a user interruption must not start another expensive branch
    atomic_json(output / "physical_execution_summary.json", {
        "created_utc": utc_now(), "status": "ERROR" if errors else "PASS",
        "requested_spacecraft": spacecraft, "errors": errors,
    })
    if errors:
        raise RuntimeError("physical branches incomplete; see physical_execution_summary.json: " + str(errors))
    if args.run_only:
        print(f"[PASS] physical execution only: {output}")
        return 0
    comparison = compare_all(data_path, output, spacecraft)
    if args.convergence_base:
        base = args.convergence_base.expanduser().resolve()
        # Evaluate the same saved native products at every common epoch. A FAIL
        # returns normally so figures/reports are still generated; malformed or
        # incomplete evidence remains an error and never becomes a missing PASS.
        try:
            compare_convergence(base, output, spacecraft, output / "convergence_summary")
            compare_convergence(base, output, spacecraft, base / "convergence_summary")
            base_comparison = compare_all(data_path, base, spacecraft)
            base_figures = plot_all(data_path, base, spacecraft, tuple(args.formats))
            write_report(base, base_comparison, base_figures)
        except (FileNotFoundError, KeyError, ValueError) as error:
            # A failed/incomplete coarse campaign must not discard valid fine-run
            # observational results. Preserve the exact error and still render the
            # fine report; absent convergence keeps qualification non-PASS.
            atomic_json(output / "refinement_error.json", {
                "status": "ERROR", "created_utc": utc_now(), "coarse_root": str(base),
                "fine_root": str(output), "error": str(error)})
    figures = plot_all(data_path, output, spacecraft, tuple(args.formats))
    report = write_report(output, comparison, figures)
    convergence_path = output / "convergence_summary.json"
    convergence = json.loads(convergence_path.read_text()) if convergence_path.is_file() else None
    qualification = qualification_status(comparison, convergence, diagnostic_evidence(output))
    if (output / "refinement_error.json").is_file() and qualification == "PASS":
        qualification = "INCOMPLETE"
    workflow = record_workflow_completion(
        output, comparison, figures, report, qualification, args.profile
    )
    print(
        f"[{qualification}] qualification; observational={comparison['status']}; "
        f"workflow={workflow['status']}; report={report}",
        flush=True,
    )
    # WORKFLOW is deliberately an exit-code policy, not a scientific policy. It is
    # used by repository orchestration to make failures in execution/comparison/plot
    # production actionable while the unchanged qualification remains machine-readable.
    if args.exit_code_policy == "WORKFLOW":
        return 0
    return 0 if qualification == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
