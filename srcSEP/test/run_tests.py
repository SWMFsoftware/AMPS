#!/usr/bin/env python3
"""Run srcSEP tests and render analytical-comparison evidence.

The C++ registry remains the sole authority for test selection and physics.
This program only orchestrates the linked executable or documented Make
targets, preserves their structured reports, and visualizes those reports.  It
never reimplements a mover or computes an "expected" value by calling the
production numerical routine.

Two plot inputs are supported:

1. A test artifact CSV containing an independent coordinate plus numerical and
   analytical columns produces a solution-overlay plot.
2. Otherwise, registry metrics produce an acceptance-comparison plot showing
   the observed error/value against the analytical reference or tolerance.

The second form is intentionally labelled metric-level evidence; it must not be
described as a pointwise solution comparison when a test reports only moments,
norms, probabilities, or convergence order.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import threading
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]

# Make target names are kept in one audited mapping.  The Python interface uses
# short stable labels while Make remains authoritative for compilers, source
# lists, sanitizers, and the exact focused executable being tested.
SOURCE_SUITES: Dict[str, str] = {
    "cli": "test-cli-unit",
    "state": "test-state-unit",
    "geometry-source": "test-geometry-source-unit",
    "mover-api": "test-mover-api-unit",
    "field-line-scope": "test-field-line-scope-unit",
    "transport-common": "test-transport-common-unit",
    "parker": "test-parker-unit",
    "fte-dmumu": "test-fte-dmumu-unit",
    "fte-mfp": "test-fte-mfp-unit",
    "coefficients": "test-coefficients-unit",
    "turbulence": "test-turbulence-core-unit",
    "reproducibility": "test-reproducibility-unit",
    "wp01-wp10": "test-wp01-wp10-unit",
    "wp11-wp20": "test-wp11-wp20-unit",
    "wp21-wp30": "test-wp21-wp30-unit",
    "wp31-wp41": "test-wp31-wp41-unit",
    "python-runner": "test-python-runner-unit",
    "acceptance": "test-acceptance-unit",
    "documentation": "test-documentation-unit",
    "controlled-analytical": "test-controlled-analytical",
    "scientific-validation": "test-scientific-validation",
    "sanitizer": "test-sanitizer",
    "stochastic-repeat": "test-stochastic-repeat",
}

# These IDs are the controlled cases whose reference is analytical or an exact
# conservation/limit invariant.  Keeping the classification explicit prevents
# a generic software contract from being advertised as an analytical solution.
ANALYTICAL_IDS = {
    "CV01", "CV02", "CV03", "CV04", "CV05", "CV06", "CV07", "CV08",
    "CV09", "CV10", "CV11", "CV12", "IV01", "IV02", "IV03", "IV04",
    "IV05", "IV06", "DXX01", "FTE01",
    "PARKER01", "TURB01", "VAL01", "CROSS02",
    *(f"PARK{i:02d}" for i in range(1, 8)),
    *(f"FTED{i:02d}" for i in range(1, 9)),
    *(f"FTEM{i:02d}" for i in range(1, 9)),
    *(f"TURB{i:02d}" for i in (2, 3, 5, 6, 7, 9, 11, 12, 13, 16, 21, 22, 23)),
}

X_COLUMNS = ("x", "time", "time_s", "s", "s_m", "mu", "radius", "radius_m",
             "energy", "energy_mev", "coordinate")
NUMERICAL_COLUMNS = ("numerical", "model", "simulated", "simulation")
ANALYTICAL_COLUMNS = ("analytic", "analytical", "exact", "reference", "expected")


# RawDescriptionHelpFormatter preserves the intentional line breaks in the
# examples below, while ArgumentDefaultsHelpFormatter keeps default executable,
# output, and plotting choices visible.  Combining them makes `--help` useful
# both as a compact reference and as a copy-and-paste quick-start guide.
class _RunnerHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                           argparse.RawDescriptionHelpFormatter):
    """Preserve annotated examples and display argument default values."""


HELP_EPILOG = r"""
Examples:
  1. Discover the tests registered in a linked srcSEP/AMPS executable:

       python3 test/run_tests.py --amps ../amps --list

     This is a pre-initialization query: it prints stable IDs and groups but
     does not start a production model run.

  2. Run one analytical test and write JSON, JUnit, PNG, and EPS evidence:

       python3 test/run_tests.py --amps ../amps --test PARK07 \
         --output-dir test_output/parker07

     Run the end-to-end CV01 validation case from its production-style input:

       python3 test/run_tests.py --amps ../amps --validation-case CV01 \
         --output-dir test_output/CV01

     CV01 invokes that exact linked application six times through its native
     --test CV01 registry entry, runs periodic/open boundaries at three
     timesteps, evaluates an independent closed-form reference, and writes both
     the standard overlay and four-panel review figure. Use --case-input to
     test one reviewed configuration variant.

     Run the remaining controlled transport cases individually or together:

       python3 test/run_tests.py --amps ../amps \
         --validation-case CV02 --validation-case CV03 \
         --validation-case CV04 --validation-case CV05 \
         --output-dir test_output/CV02-CV05

     These validate Gaussian spatial diffusion, nonuniform Ito drift,
     adiabatic cooling, and magnetic focusing. Each invokes the linked binary,
     evaluates an independent reference, and emits PNG/EPS comparisons.

     Run the advanced controlled cases with the same application contract:

       python3 test/run_tests.py --amps ../amps \
         --validation-case CV06 --validation-case CV07 \
         --validation-case CV08 --validation-case CV09 \
         --validation-case CV10 --validation-case CV11 \
         --validation-case CV12 --output-dir test_output/CV06-CV12

     CV06-CV09 are extended stochastic/distribution campaigns. CV10-CV12 are
     bounded turbulence-advection, source-history, and energy-ledger checks.

     Run all integrated manufactured cases through the linked application:

       python3 test/run_tests.py --amps ../amps \
         --validation-case IV01 --validation-case IV02 \
         --validation-case IV03 --validation-case IV04 \
         --validation-case IV05 --validation-case IV06 \
         --output-dir test_output/IV01-IV06

     These combine geometry/operators, moving grids and shocks, and nonlinear
     wave feedback. They are required nightly cases in the validation plan.

     Repeat --test to choose any collection of individual cases:

       python3 test/run_tests.py --amps ../amps \
         --test PARK01 --test FTED08 --test TURB22 \
         --output-dir test_output/selected

  3. Run one or more complete registry groups:

       python3 test/run_tests.py --amps ../amps \
         --group parker --group fte-dmumu \
         --output-dir test_output/mover-groups

     Repeated and overlapping IDs/groups are de-duplicated by the C++ registry.

  4. Run the bounded routine set used for normal development regression:

       python3 test/run_tests.py --amps ../amps --routine \
         --output-dir test_output/routine

     --routine forwards the native --all-tests policy and excludes tests marked
     extended, keeping the development gate bounded.

  5. Run every discoverable test, including extended Monte Carlo cases:

       python3 test/run_tests.py --amps ../amps --all \
         --output-dir test_output/all-tests

     --all can be substantially more expensive than --routine.  The runner
     discovers every ID with --list-tests and selects each one explicitly.

  6. Run source-only controlled tests when a linked AMPS executable is absent:

       python3 test/run_tests.py --suite controlled-analytical \
         --output-dir test_output/controlled

     Individual dependency-light suites may be repeated, for example:

       python3 test/run_tests.py --suite parker --suite fte-dmumu \
         --output-dir test_output/focused-movers

     This mode checks component kernels only. It does not replace the linked
     CV01 command in example 2.

  7. Launch the selected native tests under MPI:

       python3 test/run_tests.py --amps ../amps --group turbulence \
         --mpi-np 4 --mpiexec mpiexec \
         --output-dir test_output/turbulence-mpi

  8. Replot an existing registry report without rerunning any physics:

       python3 test/run_tests.py --from-json results.json \
         --formats png,eps --output-dir test_output/replot

     A comparison CSV produces a numerical-versus-analytical solution overlay.
     Scalar-only cases produce a clearly labelled metric/reference plot.

  9. Run tests without figures, or plot every comparable reported metric:

       python3 test/run_tests.py --amps ../amps --routine --plot none
       python3 test/run_tests.py --from-json results.json --plot all \
         --formats png --output-dir test_output/all-metric-plots

 10. Pass model-specific arguments unchanged to the linked executable:

       python3 test/run_tests.py --amps ../amps --test PARK01 \
         --output-dir test_output/parker-configured -- \
         --mover parker --turbulence-source prescribed

     The literal -- ends runner-option parsing; every following token is sent
     to srcSEP/AMPS.  Put all runner options before that separator.

Outputs:
  The selected output directory contains the command log, run manifest,
  registry JSON/JUnit evidence when available, analytical_plot_manifest.json,
  and plots/*.png plus plots/*.eps.  Test status remains authoritative: PASS,
  FAIL, and ERROR propagate through exit codes 0, 1, and 2 respectively.
"""


class RunnerError(RuntimeError):
    """A configuration, execution, report, or plotting contract failed."""


def _utc_stamp() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).strftime(
        "%Y%m%dT%H%M%SZ")


def _safe_float(value: Any) -> Optional[float]:
    """Return a finite float or None for JSON NaN/Infinity/string sentinels."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    # Reports are committed by replacement so an interrupted Python runner does
    # not leave a syntactically valid-looking partial evidence file.
    stage = path.with_name(path.name + ".stage")
    with stage.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(stage, path)


def _run_streaming(command: Sequence[str], cwd: Path, log_path: Path,
                   environment: Optional[Dict[str, str]] = None,
                   timeout_s: Optional[float] = None) -> int:
    """Execute a command while mirroring combined output to screen and log."""
    print("RUN:", shlex.join(str(item) for item in command), flush=True)
    started = _datetime.datetime.now(_datetime.timezone.utc)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{started.isoformat()}] RUN {shlex.join(command)}\n")
        log.flush()
        try:
            process = subprocess.Popen(
                list(command), cwd=str(cwd), env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1)
        except OSError as error:
            # Convert missing launchers, permission errors, and invalid working
            # environments into the runner's stable ERROR/status-2 contract
            # instead of exposing an implementation traceback to automation.
            raise RunnerError(
                f"cannot start test command {command[0]}: {error}") from error

        def copy_output() -> None:
            """Drain the pipe continuously so verbose native tests cannot block."""
            assert process.stdout is not None
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log.write(line)

        # Reading in a helper thread lets the controlling thread enforce a
        # wall-clock timeout even when a stalled test emits no newline.  A
        # simple `for line in process.stdout` followed by wait(timeout=...)
        # cannot do that because the line iterator itself may block forever.
        reader = threading.Thread(target=copy_output, daemon=True)
        reader.start()
        try:
            status = process.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
            reader.join()
            log.write(f"ERROR timeout after {timeout_s} seconds\n")
            raise RunnerError(f"test command exceeded {timeout_s} seconds")
        reader.join()
        return status


def _parse_list_output(text: str) -> List[str]:
    """Extract stable registry IDs from the documented pipe-delimited list."""
    ids: List[str] = []
    for line in text.splitlines():
        match = re.match(r"^([A-Za-z][A-Za-z0-9_-]*)\s*\|", line)
        if match:
            ids.append(match.group(1))
    if not ids:
        raise RunnerError("--list-tests returned no parseable registry IDs")
    return sorted(set(ids), key=str.casefold)


def _discover_all_ids(executable: Path, model_args: Sequence[str]) -> List[str]:
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RunnerError(f"linked srcSEP executable is not executable: {executable}")
    completed = subprocess.run(
        [str(executable), "--list-tests", *model_args], cwd=str(ROOT),
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False)
    sys.stdout.write(completed.stdout)
    if completed.returncode != 0:
        raise RunnerError(
            f"test discovery failed with exit status {completed.returncode}")
    return _parse_list_output(completed.stdout)


def _native_command(args: argparse.Namespace, report_json: Path,
                    report_junit: Path) -> List[str]:
    executable = Path(args.amps).expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RunnerError(f"linked srcSEP executable is not executable: {executable}")

    selection: List[str] = []
    if args.all:
        # The native --all-tests contract intentionally excludes extended Monte
        # Carlo cases.  Python --all means literally every discoverable test, so
        # expand the registry IDs into explicit selectors.
        for test_id in _discover_all_ids(executable, args.model_args):
            selection.extend(("--test", test_id))
    elif args.routine:
        selection.append("--all-tests")
    else:
        for test_id in args.tests:
            selection.extend(("--test", test_id))
        for group in args.groups:
            selection.extend(("--test-group", group))
    if not selection:
        raise RunnerError("select --test, --group, --routine, --all, --list, or --suite")

    command = [str(executable), *selection,
               "--test-json", str(report_json),
               "--test-junit", str(report_junit), *args.model_args]
    if args.mpi_np:
        command = [args.mpiexec, "-n", str(args.mpi_np), *command]
    return command


def _load_report(path: Path) -> Dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            report = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise RunnerError(f"cannot read test JSON {path}: {error}") from error
    if not isinstance(report, dict) or not isinstance(report.get("results"), list):
        raise RunnerError(f"test JSON has no results array: {path}")
    return report


def _merge_reports(paths: Iterable[Path], destination: Path) -> Dict[str, Any]:
    """Merge retained focused reports, de-duplicating repeated IDs."""
    by_id: Dict[str, Dict[str, Any]] = {}
    for path in sorted(paths):
        if path.resolve() == destination.resolve():
            continue
        report = _load_report(path)
        for result in report["results"]:
            if isinstance(result, dict) and result.get("id"):
                by_id[str(result["id"])] = result
    results = [by_id[key] for key in sorted(by_id, key=str.casefold)]
    totals = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0}
    status_key = {"PASS": "passed", "FAIL": "failed",
                  "SKIP": "skipped", "ERROR": "errors"}
    for result in results:
        key = status_key.get(str(result.get("status", "")).upper())
        if key:
            totals[key] += 1
    merged = {
        "schema": "srcsep-component-tests-v1",
        "exit_code": 2 if totals["errors"] else (1 if totals["failed"] else 0),
        "totals": totals,
        "results": results,
    }
    _write_json(destination, merged)
    return merged


def _column(fieldnames: Iterable[str], candidates: Sequence[str]) -> Optional[str]:
    normalized = {name.strip().casefold(): name for name in fieldnames if name}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def _read_comparison_csv(path: Path) -> Optional[Tuple[str, str, str,
                                                       List[float], List[float],
                                                       List[float], str, str]]:
    """Read a convention-based independent numerical/analytical series."""
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                return None
            x_name = _column(reader.fieldnames, X_COLUMNS)
            numerical_name = _column(reader.fieldnames, NUMERICAL_COLUMNS)
            analytical_name = _column(reader.fieldnames, ANALYTICAL_COLUMNS)
            if not x_name or not numerical_name or not analytical_name:
                return None
            x_values: List[float] = []
            numerical: List[float] = []
            analytical: List[float] = []
            quantity = "solution"
            units = ""
            for row in reader:
                x = _safe_float(row.get(x_name))
                model = _safe_float(row.get(numerical_name))
                exact = _safe_float(row.get(analytical_name))
                if x is not None and model is not None and exact is not None:
                    x_values.append(x)
                    numerical.append(model)
                    analytical.append(exact)
                    # Optional metadata columns keep the basic three-column
                    # convention backward compatible while allowing physical
                    # axes for campaign cases such as CV01.
                    if row.get("quantity"):
                        quantity = str(row["quantity"]).strip()
                    if row.get("units"):
                        units = str(row["units"]).strip()
            if len(x_values) < 2:
                return None
            return (x_name, numerical_name, analytical_name, x_values,
                    numerical, analytical, quantity, units)
    except (OSError, csv.Error):
        return None


def _comparison_operator(metric: Dict[str, Any]) -> str:
    text = str(metric.get("comparison", "")).replace(" ", "")
    for operator in ("<=", ">=", "==", "<", ">"):
        if operator in text:
            return operator
    return "reference"


def _save_figure(figure: Any, stem: Path, formats: Sequence[str]) -> List[str]:
    created: List[str] = []
    for extension in formats:
        target = stem.with_suffix("." + extension)
        figure.savefig(target, format=extension, dpi=180, bbox_inches="tight",
                       facecolor="white")
        created.append(str(target))
    return created


def _plot_series(plt: Any, test_id: str, status: str, series: Tuple[Any, ...],
                 stem: Path, formats: Sequence[str]) -> List[str]:
    (x_name, numerical_name, analytical_name, x, numerical, analytical,
     quantity, units) = series
    figure, axes = plt.subplots(
        2, 1, figsize=(8.5, 6.8), sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.25]})
    solution_axis, residual_axis = axes
    solution_axis.plot(x, analytical, color="black", linewidth=2.0,
                       label=f"Analytical ({analytical_name})")
    solution_axis.plot(x, numerical, color="#1565c0", linewidth=1.5,
                       marker="o", markersize=3.0,
                       label=f"Numerical ({numerical_name})")
    y_label = quantity + (f" [{units}]" if units else "")
    solution_axis.set_ylabel(y_label)
    solution_axis.set_title(
        f"{test_id}: numerical and analytical solution [{status}]")
    solution_axis.grid(True, color="0.85", linewidth=0.7)
    solution_axis.legend(frameon=False)

    # The residual is calculated only from values reloaded from the archived
    # CSV.  This makes the figure independently reproducible and prevents an
    # in-memory plotting path from hiding serialization or unit mistakes.
    residual = [model - exact for model, exact in zip(numerical, analytical)]
    residual_axis.axhline(0.0, color="black", linewidth=1.0)
    residual_axis.plot(x, residual, color="#c62828", linewidth=1.2,
                       marker="o", markersize=2.8)
    residual_axis.set_xlabel(x_name)
    residual_axis.set_ylabel("model - reference" +
                             (f" [{units}]" if units else ""))
    residual_axis.grid(True, color="0.88", linewidth=0.7)
    figure.tight_layout()
    outputs = _save_figure(figure, stem, formats)
    plt.close(figure)
    return outputs


def _plot_metrics(plt: Any, result: Dict[str, Any], stem: Path,
                  formats: Sequence[str]) -> List[str]:
    metrics = [item for item in result.get("metrics", [])
               if isinstance(item, dict) and
               item.get("name") != "assertion_failures" and
               _safe_float(item.get("value")) is not None and
               _safe_float(item.get("tolerance")) is not None]
    if not metrics:
        return []
    figure, axes = plt.subplots(
        len(metrics), 1, figsize=(9.0, max(3.0, 2.05 * len(metrics))),
        squeeze=False)
    for axis, metric in zip(axes[:, 0], metrics):
        observed = float(metric["value"])
        reference = float(metric["tolerance"])
        operator = _comparison_operator(metric)
        spread = max(abs(observed), abs(reference), 1.0e-30)
        lower = min(0.0, observed, reference) - 0.15 * spread
        upper = max(0.0, observed, reference) + 0.15 * spread
        if lower == upper:
            lower, upper = -1.0, 1.0
        # Green shading means the acceptance side of a one-sided analytical
        # comparison. Equality/reference metrics use a line only because the
        # registry tolerance is the exact expected value, not an uncertainty.
        if operator in ("<=", "<"):
            axis.axvspan(lower, reference, color="#dcedc8", zorder=0)
        elif operator in (">=", ">"):
            axis.axvspan(reference, upper, color="#dcedc8", zorder=0)
        axis.axvline(reference, color="black", linestyle="--", linewidth=1.3,
                     label="analytical reference / acceptance limit")
        axis.scatter([observed], [0.0], s=55, color="#1565c0", zorder=3,
                     label="numerical result")
        axis.set_xlim(lower, upper)
        axis.set_yticks([])
        units = str(metric.get("units", "")).strip()
        axis.set_xlabel(units or "reported units")
        axis.set_title(
            f"{metric.get('name', 'metric')}: {observed:.6g} {operator} {reference:.6g}",
            fontsize=10)
        axis.grid(True, axis="x", color="0.88", linewidth=0.7)
    axes[0, 0].legend(loc="best", frameon=False, fontsize=8)
    test_id = str(result.get("id", "unknown"))
    status = str(result.get("status", "UNKNOWN"))
    seed = result.get("seed")
    subtitle = "metric-level analytical/acceptance comparison"
    if seed is not None:
        subtitle += f"; seed={seed}"
    figure.suptitle(f"{test_id} [{status}]\n{subtitle}", fontsize=12)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    outputs = _save_figure(figure, stem, formats)
    plt.close(figure)
    return outputs


def generate_plots(report: Dict[str, Any], report_path: Path, output_dir: Path,
                   formats: Sequence[str], plot_mode: str) -> Dict[str, Any]:
    """Render all eligible test comparisons and return an evidence manifest."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RunnerError(
            "plotting requires matplotlib (python3 -m pip install matplotlib)") from error

    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    manifest: Dict[str, Any] = {
        "schema": "srcsep-analytical-plots-v1",
        "source_report": str(report_path),
        "source_report_sha256": _sha256(report_path),
        "formats": list(formats),
        "plots": [],
    }
    for result in report.get("results", []):
        if not isinstance(result, dict):
            continue
        test_id = str(result.get("id", "unknown"))
        artifacts = result.get("artifacts", [])
        series = None
        series_source: Optional[Path] = None
        for artifact in artifacts if isinstance(artifacts, list) else []:
            candidate = Path(str(artifact))
            if not candidate.is_absolute():
                # Native tests normally write artifacts beside their report.
                candidate = report_path.parent / candidate
            if candidate.suffix.casefold() == ".csv" and candidate.is_file():
                series = _read_comparison_csv(candidate)
                if series:
                    series_source = candidate
                    break
        eligible = plot_mode == "all" or test_id.upper() in ANALYTICAL_IDS
        if not eligible and not series:
            continue
        stem = plot_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', test_id)}_comparison"
        if series:
            outputs = _plot_series(
                plt, test_id, str(result.get("status", "UNKNOWN")),
                series, stem, formats)
            kind = "solution-series"
        else:
            outputs = _plot_metrics(plt, result, stem, formats)
            kind = "metric-acceptance"
        if outputs:
            manifest["plots"].append({
                "test_id": test_id,
                "status": result.get("status"),
                "kind": kind,
                "series_source": str(series_source) if series_source else None,
                "files": outputs,
            })
    _write_json(output_dir / "analytical_plot_manifest.json", manifest)
    return manifest


def _run_source_suites(args: argparse.Namespace, output_dir: Path,
                       log_path: Path) -> Tuple[int, Optional[Path]]:
    report_dir = output_dir / "focused_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ)
    # Focused C++ runners copy their already-validated JSON/JUnit products only
    # when this variable is set. Normal Make invocations remain artifact-free.
    environment["SRCSEP_REPORT_DIR"] = str(report_dir)
    exit_code = 0
    for suite in args.suites:
        target = SOURCE_SUITES[suite]
        status = _run_streaming(
            ["make", target], ROOT, log_path, environment, args.timeout)
        if status != 0:
            exit_code = status
            if not args.keep_going:
                break
    reports = list(report_dir.glob("*.json"))
    if not reports:
        return exit_code, None
    merged_path = output_dir / "srcsep-tests.json"
    _merge_reports(reports, merged_path)
    return exit_code, merged_path


def _run_validation_cases(args: argparse.Namespace, output_dir: Path,
                          log_path: Path) -> Tuple[int, Path, List[str]]:
    """Run application-level cases through the requested linked executable."""
    command = [sys.executable, str(ROOT / "validation" / "run_case.py"),
               "--amps", str(Path(args.amps).expanduser().resolve())]
    if args.validation_all:
        command.append("--all")
    else:
        for case_id in args.validation_cases:
            command.extend(("--case", case_id))
    if args.case_input:
        command.extend(("--input", str(args.case_input.expanduser().resolve())))
    if args.timeout is not None:
        command.extend(("--timeout", str(args.timeout)))
    command.extend(("--output-dir", str(output_dir)))
    exit_code = _run_streaming(
        command, ROOT, log_path, dict(os.environ), args.timeout)
    return exit_code, output_dir / "srcsep-tests.json", command


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run any/all srcSEP tests and plot analytical comparisons.",
        # Keep the examples readable instead of collapsing their explanations
        # into one wrapped paragraph.  Defaults are still added to each option
        # by the first formatter base class.
        formatter_class=_RunnerHelpFormatter,
        epilog=HELP_EPILOG)
    parser.add_argument("--amps", default=os.environ.get(
        "SEP_EXECUTABLE", str(ROOT.parent / "amps")),
        help="linked srcSEP/AMPS executable")
    parser.add_argument("--list", action="store_true",
                        help="list registered native tests and exit")
    parser.add_argument("--test", dest="tests", action="append", default=[],
                        metavar="ID", help="run one test; repeatable")
    parser.add_argument("--group", dest="groups", action="append", default=[],
                        help="run one registry group; repeatable")
    parser.add_argument("--routine", action="store_true",
                        help="run the bounded native --all-tests set")
    parser.add_argument("--all", action="store_true",
                        help="discover and run every test, including extended cases")
    parser.add_argument("--suite", dest="suites", action="append", default=[],
                        choices=sorted(SOURCE_SUITES),
                        help="run a dependency-light Make suite; repeatable")
    parser.add_argument("--validation-case", dest="validation_cases",
                        action="append", default=[], metavar="ID",
                        help="run one registered end-to-end validation case; repeatable")
    parser.add_argument("--validation-all", action="store_true",
                        help="run every registered end-to-end validation case")
    parser.add_argument("--case-input", type=Path,
                        help="override the input for one --validation-case")
    parser.add_argument("--from-json", type=Path,
                        help="plot an existing component-test JSON without running tests")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("test_output") / f"srcsep-tests-{_utc_stamp()}",
                        help="report, log, and figure destination")
    parser.add_argument("--plot", choices=("analytical", "all", "none"),
                        default="analytical",
                        help="which result records receive plots")
    parser.add_argument("--formats", default="png,eps",
                        help="comma-separated plot formats: png and/or eps")
    parser.add_argument("--mpi-np", type=int,
                        help="launch the native executable with this rank count")
    parser.add_argument("--mpiexec", default="mpiexec",
                        help="MPI launcher used with --mpi-np")
    parser.add_argument("--timeout", type=float,
                        help="per-command timeout in seconds")
    parser.add_argument("--keep-going", action="store_true",
                        help="continue later --suite targets after a failure")
    parser.add_argument("model_args", nargs=argparse.REMAINDER,
                        help="arguments after -- are passed unchanged to AMPS")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.model_args and args.model_args[0] == "--":
        args.model_args = args.model_args[1:]
    formats = [item.strip().casefold() for item in args.formats.split(",")
               if item.strip()]
    if not formats or any(item not in ("png", "eps") for item in formats):
        raise RunnerError("--formats accepts only png and eps")
    if args.mpi_np is not None and args.mpi_np <= 0:
        raise RunnerError("--mpi-np must be positive")

    native_selection = bool(args.tests or args.groups or args.routine or args.all)
    validation_selection = bool(args.validation_cases or args.validation_all)
    modes = sum((bool(args.suites), bool(args.from_json), native_selection,
                 validation_selection, bool(args.list)))
    if modes != 1:
        raise RunnerError(
            "choose exactly one mode: native selection, validation case, "
            "--suite, --from-json, or --list")
    if args.all and (args.routine or args.tests or args.groups):
        raise RunnerError("--all cannot be combined with other native selectors")
    if args.routine and (args.tests or args.groups):
        raise RunnerError("--routine cannot be combined with --test or --group")
    if args.validation_all and args.validation_cases:
        raise RunnerError("--validation-all cannot be combined with --validation-case")
    if args.case_input and (args.validation_all or len(args.validation_cases) != 1):
        raise RunnerError("--case-input requires exactly one --validation-case")
    if validation_selection and args.mpi_np is not None:
        raise RunnerError(
            "end-to-end validation cases currently require serial linked execution; "
            "omit --mpi-np")

    if args.list:
        executable = Path(args.amps).expanduser().resolve()
        ids = _discover_all_ids(executable, args.model_args)
        print(f"Discovered {len(ids)} registered tests.")
        return 0

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "test-run.log"
    report_path: Optional[Path] = None
    exit_code = 0
    command: Optional[List[str]] = None

    if args.from_json:
        report_path = args.from_json.expanduser().resolve()
        report = _load_report(report_path)
        exit_code = int(report.get("exit_code", 0))
    elif validation_selection:
        exit_code, report_path, command = _run_validation_cases(
            args, output_dir, log_path)
        report = _load_report(report_path) if report_path.is_file() else None
    elif args.suites:
        exit_code, report_path = _run_source_suites(args, output_dir, log_path)
        report = _load_report(report_path) if report_path else None
    else:
        report_path = output_dir / "srcsep-tests.json"
        junit_path = output_dir / "srcsep-tests.xml"
        command = _native_command(args, report_path, junit_path)
        exit_code = _run_streaming(
            command, ROOT, log_path, dict(os.environ), args.timeout)
        report = _load_report(report_path) if report_path.is_file() else None

    manifest: Dict[str, Any] = {
        "schema": "srcsep-python-test-run-v1",
        "created_utc": _datetime.datetime.now(
            _datetime.timezone.utc).isoformat(),
        "root": str(ROOT),
        "command": command,
        "source_suites": args.suites,
        "validation_cases": args.validation_cases,
        "validation_all": args.validation_all,
        "report": str(report_path) if report_path else None,
        "exit_code": exit_code,
    }
    if native_selection or validation_selection:
        executable = Path(args.amps).expanduser().resolve()
        # Hashing the linked executable makes the orchestration record useful
        # as provenance: plots can be tied to the exact production binary that
        # emitted their numerical metrics, not merely to a path that may later
        # be overwritten by another AMPS build.
        if executable.is_file():
            manifest["executable"] = str(executable)
            manifest["executable_sha256"] = _sha256(executable)
    if report_path and report_path.is_file():
        manifest["report_sha256"] = _sha256(report_path)
    _write_json(output_dir / "run_manifest.json", manifest)

    if args.plot != "none" and report is not None and report_path is not None:
        plots = generate_plots(
            report, report_path, output_dir, formats, args.plot)
        print(f"Generated {len(plots['plots'])} analytical comparison figure set(s).")
    elif args.plot != "none":
        print("No structured JSON was retained; no analytical plots were generated.")

    print(f"Results: {output_dir}")
    return exit_code


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RunnerError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
