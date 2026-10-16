#!/usr/bin/env python3
"""Run the end-to-end GOES01 AMPS/GOES-16 observational validation."""

from __future__ import annotations

import argparse
import json
import os
import platform
import signal
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from compare_GOES01 import PATTERNS, compare, compare_convergence
from goes01_common import PROFILES, atomic_json, sha256_file, utc_now
from plot_GOES01 import plot
from prepare_GOES01 import prepare


def _campaign(parent: Path, profile: str) -> Path:
    """Reserve a unique evidence directory; prior campaigns are immutable."""
    parent = parent.expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    stem = (f"goes01_{profile.lower()}_" +
            datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + f"_pid{os.getpid()}")
    for suffix in range(1000):
        candidate = parent / (stem if suffix == 0 else f"{stem}_{suffix:03d}")
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError(f"cannot reserve a unique campaign below {parent}")


def _repository_identity(root: Path) -> dict[str, Any]:
    """Bind evidence to HEAD and dirty paths; HEAD alone is insufficient provenance."""
    def git(*arguments: str) -> str:
        result = subprocess.run(["git", *arguments], cwd=root, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"git {' '.join(arguments)} failed: {result.stdout}")
        return result.stdout.strip()
    return {"root": str(root), "head": git("rev-parse", "HEAD"),
            "status_porcelain": git("status", "--porcelain=v1").splitlines()}


def _terminate(process: subprocess.Popen[Any]) -> int:
    """Terminate only the launcher's private process group on user interruption."""
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


def _products(branch: Path, expected: int) -> tuple[list[str], list[Path]]:
    missing: list[str] = []
    files: list[Path] = []
    for role, pattern in PATTERNS.items():
        matches = sorted(branch.glob(pattern)); files.extend(matches)
        if len(matches) != expected:
            missing.append(f"{role}: expected {expected}, found {len(matches)} ({pattern})")
    return missing, files


def run_physics(output_root: Path, amps_root: Path, amps_argument: Path,
                mpi_ranks: int, threads: int) -> dict[str, Any]:
    """Launch the immutable GRIDDED/HC4 branch and hash every native product."""
    branch = output_root / "GOES16"
    prepared = json.loads((branch / "prepared_run.json").read_text(encoding="utf-8"))
    for path_key, hash_key in (
        ("input_file", "input_sha256"), ("driver_file", "driver_sha256"),
        ("spectrum_file", "spectrum_sha256"),
        ("snapshot_list_file", "snapshot_list_sha256"),
        ("observations_file", "observations_sha256"),
    ):
        if sha256_file(Path(prepared[path_key])) != prepared[hash_key]:
            raise ValueError(f"prepared input changed before execution: {prepared[path_key]}")
    expected = int(prepared["trajectory"]["kept"])
    _, existing = _products(branch, expected)
    if existing:
        raise FileExistsError("physical output already exists; use a new campaign directory")
    amps = amps_argument if amps_argument.is_absolute() else amps_root / amps_argument
    amps = amps.resolve()
    if not amps.is_file() or not os.access(amps, os.X_OK):
        raise FileNotFoundError(f"AMPS executable is absent or not executable: {amps}")
    mesh = prepared["mode3d_mesh"]
    command = [
        "mpiexec", "-n", str(mpi_ranks), str(amps), "-mode", "3d",
        "-i", "input/AMPS_PARAM.in", "-mode3d-field-eval", "INTERPOLATION",
        "-mode3d-trajectory-snapshot-scheduler", "EPOCH_PARALLEL",
        "-mode3d-epoch-parallel-memory-fraction", "0.80",
        "-density-parallel", "OPENMP", "-density-threads", str(threads),
        "-mode3d-mpi-scheduler", "DYNAMIC", "-mode3d-mpi-dynamic-chunk", "1",
        "-mode3d-mesh-res-earth-re", str(mesh["res_earth_re"]),
        "-mode3d-mesh-res-boundary-re", str(mesh["res_boundary_re"]),
        "-mode3d-mesh-coarsening", str(mesh["coarsening"]),
        "-mode3d-mesh-exponent", str(mesh["exponent"]),
        "-mode3d-parallel-field-init", "-mover", prepared["numerical_controls"]["mover"],
    ]
    environment = dict(os.environ); environment["OMP_NUM_THREADS"] = str(threads)
    log_path = branch / "amps.log"; started = utc_now(); interrupted = False
    atomic_json(output_root / "active_run.json", {
        "schema": "amps-goes01-active-v1", "status": "RUNNING",
        "started_utc": started, "runner_pid": os.getpid(), "command": command,
        "working_directory": str(branch), "log": str(log_path),
    })
    print("[GOES01] " + " ".join(command), flush=True)
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        log.write("Command: " + " ".join(command) + "\nWorking directory: " + str(branch) + "\n\n")
        log.flush()
        process = subprocess.Popen(command, cwd=branch, env=environment, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   start_new_session=True, bufsize=1)
        try:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line); log.flush(); print(line, end="", flush=True)
            return_code = process.wait()
        except KeyboardInterrupt:
            interrupted = True; return_code = _terminate(process)
            log.write("\nRunner interrupted; launcher process group terminated.\n")
    missing, files = _products(branch, expected)
    identity = _repository_identity(amps_root)
    status = {
        "schema": "amps-goes01-physical-run-v1", "created_utc": utc_now(),
        "started_utc": started, "status": "PASS" if return_code == 0 and not missing else "ERROR",
        "return_code": return_code, "interrupted": interrupted, "missing_outputs": missing,
        "command": command, "requested_mpi_ranks": mpi_ranks,
        "requested_threads_per_rank": threads, "repository": identity,
        "executable": str(amps), "executable_sha256": sha256_file(amps),
        "output_sha256": {path.name: sha256_file(path) for path in files},
        **{key: prepared[key] for key in ("input_sha256", "driver_sha256",
                                           "spectrum_sha256", "snapshot_list_sha256",
                                           "observations_sha256")},
    }
    atomic_json(branch / "run_status.json", status)
    atomic_json(output_root / "physical_execution_summary.json", {
        "schema": "amps-goes01-physical-summary-v1", "status": status["status"],
        "return_code": return_code, "missing_outputs": missing,
    })
    if status["status"] != "PASS":
        raise RuntimeError(f"AMPS physical branch failed: exit={return_code}, missing={missing}")
    return status


def _qualification(observational: dict[str, Any], convergence: dict[str, Any] | None) -> str:
    """Never promote observational agreement without numerical-refinement evidence."""
    if observational["status"] == "FAIL" or (convergence and convergence["status"] == "FAIL"):
        return "FAIL"
    return "PASS" if observational["status"] == "PASS" and convergence and convergence[
        "status"] == "PASS" else "INCOMPLETE"


def _report(output_root: Path, summary: dict[str, Any], convergence: dict[str, Any] | None,
            qualification: str, figures: list[Path]) -> Path:
    """Write a compact human receipt while JSON/CSV remain authoritative evidence."""
    report = output_root / "GOES01_validation_report.md"
    lines = [
        "# GOES01 validation report", "", f"- Observational status: **{summary['status']}**",
        f"- Numerical convergence: **{convergence['status'] if convergence else 'NOT EVALUATED'}**",
        f"- Overall qualification: **{qualification}**",
        f"- Profile: `{summary['profile']}`", f"- Matched epochs: {summary['matched_epochs']}",
        f"- Comparison rows: {summary['comparison_rows']}",
        f"- Maximum unresolved fraction: {summary['maximum_unresolved_fraction']:.6g}",
        "", "## Acceptance gates", "",
        "| Gate | Value | Requirement | Result |", "|---|---:|---:|---|",
    ]
    for gate in summary["gates"]:
        lines.append(f"| `{gate['name']}` | {gate['value']} | {gate['operator']} "
                     f"{gate['threshold']} | {'PASS' if gate['passed'] else 'FAIL'} |")
    lines += ["", "## Figures", ""] + [f"- `{path.relative_to(output_root)}`" for path in figures]
    lines += ["", "## Scope", "", summary["interpretation"], ""]
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True,
                        help="observational_validation_datasets root or goes16_sgps_sep2017")
    parser.add_argument("--model-data-root", type=Path, required=True,
                        help="SA01-SA04 package root or model_inputs/SA04")
    parser.add_argument("--amps-root", type=Path, default=Path.cwd())
    parser.add_argument("--amps", type=Path, default=Path("amps"))
    parser.add_argument("--output-parent", type=Path, default=Path("test_output/GOES01"))
    parser.add_argument("--output-root", type=Path,
                        help="exact evidence directory; required for RUN/COMPARE/PLOT")
    parser.add_argument("--profile", choices=tuple(PROFILES), default="QUICK")
    parser.add_argument("--stage", choices=("ALL", "PREPARE", "RUN", "COMPARE", "PLOT"),
                        default="ALL")
    parser.add_argument("--mpi-ranks", type=int, default=4)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--mover", choices=("BORIS", "HC4", "RK4", "RK6"), default="HC4")
    parser.add_argument("--trace-dt", type=float, default=0.1)
    parser.add_argument("--mode3d-mesh-res-earth-re", type=float, default=0.05)
    parser.add_argument("--mode3d-mesh-res-boundary-re", type=float, default=2.0)
    parser.add_argument("--formats", nargs="+", choices=("png", "eps"), default=("png", "eps"))
    parser.add_argument("--convergence-base", type=Path,
                        help="completed QUICK evidence root to refine with this COMPLETE run")
    parser.add_argument("--exit-code-policy", choices=("QUALIFICATION", "WORKFLOW"),
                        default="QUALIFICATION")
    args = parser.parse_args()
    if args.mpi_ranks < 1 or args.threads < 1:
        raise ValueError("MPI ranks and threads must be positive")
    if args.stage not in ("ALL", "PREPARE") and args.output_root is None:
        raise ValueError(f"--output-root is required for stage {args.stage}")
    output = (args.output_root.expanduser().resolve() if args.output_root else
              _campaign(args.output_parent, args.profile))
    if args.stage in ("ALL", "PREPARE"):
        prepare(args.data_root, args.model_data_root, output, args.profile, args.threads,
                args.mover, args.trace_dt, args.mode3d_mesh_res_earth_re,
                args.mode3d_mesh_res_boundary_re)
    if args.stage in ("ALL", "RUN"):
        run_physics(output, args.amps_root.expanduser().resolve(), args.amps,
                    args.mpi_ranks, args.threads)
    if args.stage in ("ALL", "COMPARE"):
        summary = compare(output)
    elif (output / "comparison_summary.json").is_file():
        summary = json.loads((output / "comparison_summary.json").read_text())
    else:
        summary = {"status": "NOT_COMPARED"}
    convergence = None
    if args.convergence_base and summary["status"] in ("PASS", "FAIL"):
        convergence = compare_convergence(args.convergence_base.expanduser(), output)
    elif (output / "convergence_summary.json").is_file():
        convergence = json.loads((output / "convergence_summary.json").read_text())
    qualification = (_qualification(summary, convergence)
                     if summary["status"] in ("PASS", "FAIL") else "NOT_EVALUATED")
    figures: list[Path] = []
    if args.stage in ("ALL", "PLOT"):
        figures = plot(output, tuple(args.formats))
    report = (_report(output, summary, convergence, qualification, figures)
              if summary["status"] in ("PASS", "FAIL") else None)
    if summary["status"] in ("PASS", "FAIL"):
        atomic_json(output / "qualification_summary.json", {
            "schema": "amps-goes01-qualification-v1", "created_utc": utc_now(),
            "status": qualification, "observational_status": summary["status"],
            "convergence_status": convergence["status"] if convergence else "NOT_EVALUATED",
        })
    workflow_status = "PASS"
    # A later plotting/comparison-only action must not overwrite the receipt for the
    # expensive ALL campaign. Stage-specific receipts preserve both histories.
    receipt = (output / "workflow_summary.json" if args.stage == "ALL" else
               output / f"workflow_{args.stage.lower()}_summary.json")
    atomic_json(receipt, {
        "schema": "amps-goes01-workflow-v1", "created_utc": utc_now(),
        "status": workflow_status, "stage": args.stage, "profile": args.profile,
        "qualification_status": qualification, "observational_status": summary["status"],
        "output_root": str(output),
        "report": str(report) if report else None,
        "python": sys.version, "platform": platform.platform(),
    })
    print(f"[GOES01] workflow={workflow_status} observational={summary['status']} "
          f"qualification={qualification}")
    print(f"[GOES01] evidence={output}")
    if args.exit_code_policy == "WORKFLOW":
        return 0
    return 0 if qualification in ("PASS", "NOT_EVALUATED") else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
