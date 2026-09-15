"""Linked-application orchestration and scoring for XM01-XM03.

XM01 compares production stochastic characteristics with an independent
finite-volume PDE solver. XM02/XM03 compare a normalized production export
with immutable, publication-derived M-FLAMPA references. The external cases
return SKIP—not PASS—when the model export or the documented equivalence review
is absent. This distinction is part of the scientific evidence contract.
"""
from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from linked_case_common import (atomic_json, finish_result, load_input, metric,
                                read_csv, run_linked_model, sha256, write_csv)

def _xm01_arguments(case: Dict[str, Any]) -> List[str]:
    p, n = case["physics"], case["numerics"]
    return [
        "--length-m", str(p["length_m"]),
        "--speed-m-per-s", str(p["speed_m_per_s"]),
        "--duration-s", str(n["duration_s"]),
        "--time-step-s", str(n["time_step_s"]),
        "--d0-per-s", str(p["d0_per_s"]),
        "--dlnb-ds-per-m", str(p["dlnb_ds_per_m"]),
        "--divergence-per-s", str(p["divergence_per_s"]),
        "--s-bins", str(n["s_bins"]),
        "--mu-bins", str(n["mu_bins"]),
        "--particles", str(n["particles"]),
        "--campaign-seed", str(n["campaign_seed"]),
    ]

def _reference_path(case_id: str, case: Dict[str, Any]) -> Path:
    """Resolve a source-owned immutable reference without accepting traversal.

    Custom input files may live in an evidence directory. Anchoring the
    baseline at the registered case directory keeps such overrides portable
    and prevents them from silently replacing the reviewed reference.
    """
    relative = Path(str(case["reference"]["csv"]))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("reference.csv must be a case-relative path")
    case_directory = (Path(__file__).parent / case_id).resolve()
    resolved = (case_directory / relative).resolve()
    if case_directory not in resolved.parents or not resolved.is_file():
        raise ValueError(f"reference CSV is missing or outside the case: {resolved}")
    return resolved

def _formats(case: Dict[str, Any]) -> Sequence[str]:
    formats = case.get("plot", {}).get("formats", ["png", "eps"])
    if not isinstance(formats, list) or any(item not in ("png", "eps") for item in formats):
        raise ValueError("plot.formats must contain only png and/or eps")
    return formats

def _save_figure(figure, output: Path, stem: str,
                 formats: Sequence[str]) -> List[Path]:
    paths: List[Path] = []
    for extension in formats:
        path = output / f"{stem}.{extension}"
        figure.savefig(path, dpi=180 if extension == "png" else None)
        paths.append(path)
    return paths

def _xm01_plot(output: Path, model_rows: Sequence[Dict[str, str]],
               reference_rows: Sequence[Dict[str, str]],
               formats: Sequence[str]) -> List[Path]:
    """Plot spatial intensity for all five operator combinations."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = ("streaming", "scattering", "focusing", "adiabatic", "combined")
    figure, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
    for axis, name in zip(axes.flat, names):
        ref = [row for row in reference_rows if row["scenario"] == name]
        mod = [row for row in model_rows if row["scenario"] == name and
               int(row["refinement"]) == 2]
        # Intensity is repeated across pitch bins; selecting mu_bin zero gives
        # one auditable spatial sample per arclength cell.
        ref = sorted((row for row in ref if int(row["mu_bin"]) == 0),
                     key=lambda row: int(row["s_bin"]))
        mod = sorted((row for row in mod if int(row["mu_bin"]) == 0),
                     key=lambda row: int(row["s_bin"]))
        x = [0.5*(float(row["s_left_m"])+float(row["s_right_m"]))/1.0e9 for row in ref]
        axis.plot(x, [float(row["intensity"]) for row in ref], "k-",
                  label="independent FV")
        axis.plot(x, [float(row["intensity"]) for row in mod], "o", ms=3,
                  label="linked srcSEP/AMPS")
        axis.set_title(name)
        axis.set_ylabel("cell probability")
        axis.grid(alpha=0.25)
    axes.flat[-1].axis("off")
    axes.flat[0].legend()
    for axis in axes[-1, :]: axis.set_xlabel("arclength [10$^9$ m]")
    figure.suptitle("XM01 focused-transport cross-solver comparison")
    figure.tight_layout()
    paths = _save_figure(figure, output, "XM01_comparison", formats)
    plt.close(figure)
    return paths

def _rmse(values: Iterable[float]) -> float:
    values = list(values)
    return math.sqrt(sum(item*item for item in values)/len(values)) if values else math.inf

def _score_xm01(case: Dict[str, Any], model_rows: Sequence[Dict[str, str]],
                reference_rows: Sequence[Dict[str, str]]) -> List[Dict[str, Any]]:
    acceptance = case["acceptance"]
    reference = {(row["scenario"], int(row["s_bin"]), int(row["mu_bin"])): row
                 for row in reference_rows}
    errors: Dict[int, List[float]] = {1: [], 2: []}
    moment_error = 0.0
    momentum_error = 0.0
    for level in (1, 2):
        for scenario in ("streaming", "scattering", "focusing", "adiabatic", "combined"):
            rows = [row for row in model_rows if row["scenario"] == scenario and
                    int(row["refinement"]) == level]
            if len(rows) != len(reference_rows)//5:
                raise RuntimeError(f"XM01 {scenario} level {level} has an incomplete grid")
            peak = max(float(reference[(scenario, int(row["s_bin"]),
                                        int(row["mu_bin"]))]["probability"])
                       for row in rows)
            resolved = acceptance["resolved_probability_fraction_min"] * peak
            for row in rows:
                ref = reference[(scenario, int(row["s_bin"]), int(row["mu_bin"]))]
                expected = float(ref["probability"])
                observed = float(row["probability"])
                if expected >= resolved and observed > 0.0:
                    errors[level].append(math.log10(observed/expected))
            if level == 2:
                first, expected_first = float(rows[0]["anisotropy"]), \
                    float(reference[(scenario, 0, 0)]["anisotropy"])
                moment_error = max(moment_error, abs(first-expected_first))
                # Momentum is uniform for the reduced adiabatic characteristic;
                # inspect any occupied model cell and the corresponding scalar
                # reference rather than allowing empty-bin NaNs into a metric.
                occupied = next(row for row in rows if float(row["probability"]) > 0.0)
                momentum_error = max(momentum_error, abs(float(occupied["mean_log_p"]) -
                    float(reference[(scenario, 0, 0)]["mean_log_p"])))
    coarse, fine = _rmse(errors[1]), _rmse(errors[2])
    ratio = fine/coarse if coarse > 0.0 else 0.0
    return [
        metric("resolved_log10_probability_rmse", fine,
               acceptance["log10_probability_rmse_max"], "<=", "dex"),
        metric("first_angular_moment_difference", moment_error,
               acceptance["first_angular_moment_difference_max"], "<=", "dimensionless"),
        metric("log_momentum_difference", momentum_error,
               acceptance["log_momentum_difference_max"], "<=", "natural-log momentum"),
        metric("fine_to_coarse_error_ratio", ratio,
               acceptance["fine_to_coarse_error_ratio_max"], "<=", "dimensionless"),
    ]

def _external_plot(case_id: str, output: Path,
                   reference_rows: Sequence[Dict[str, str]],
                   model_rows: Optional[Sequence[Dict[str, str]]],
                   formats: Sequence[str]) -> List[Path]:
    """Render publication reference alone or overlaid with production output."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    series_key = "series" if case_id == "XM02" else "observable"
    x_key = "elapsed_hours" if case_id == "XM02" else "coordinate"
    series = sorted({row[series_key] for row in reference_rows})
    figure, axes = plt.subplots(len(series), 1, figsize=(8, max(3, 2.8*len(series))),
                                squeeze=False)
    for axis, name in zip(axes[:, 0], series):
        ref = sorted((row for row in reference_rows if row[series_key] == name),
                     key=lambda row: float(row[x_key]))
        axis.plot([float(row[x_key]) for row in ref], [float(row["intensity" if case_id == "XM02" else "value"]) for row in ref],
                  "k.-", label="digitized M-FLAMPA reference")
        if model_rows is not None:
            mod = sorted((row for row in model_rows if row[series_key] == name),
                         key=lambda row: float(row[x_key]))
            axis.plot([float(row[x_key]) for row in mod], [float(row["intensity" if case_id == "XM02" else "value"]) for row in mod],
                      "C1o-", ms=3, label="linked srcSEP/AMPS export")
        if all(float(row["intensity" if case_id == "XM02" else "value"]) > 0.0 for row in ref):
            axis.set_yscale("log")
        axis.set_title(name.replace("_", " "))
        axis.grid(alpha=0.25)
        axis.legend()
    axes[-1, 0].set_xlabel("elapsed hours" if case_id == "XM02" else "published coordinate")
    figure.suptitle(f"{case_id} publication-derived comparison")
    figure.tight_layout()
    paths = _save_figure(figure, output, f"{case_id}_comparison", formats)
    plt.close(figure)
    return paths

def _interpolate(points: Sequence[Tuple[float, float]], x: float) -> float:
    """Linearly interpolate in log-value when both neighboring values are positive."""
    ordered = sorted(points)
    if not ordered or x < ordered[0][0] or x > ordered[-1][0]:
        raise ValueError("model does not cover a reference coordinate")
    for left, right in zip(ordered[:-1], ordered[1:]):
        if left[0] <= x <= right[0]:
            if x == left[0]: return left[1]
            fraction = (x-left[0])/(right[0]-left[0])
            if left[1] > 0.0 and right[1] > 0.0:
                return 10.0**((1-fraction)*math.log10(left[1]) + fraction*math.log10(right[1]))
            return (1-fraction)*left[1] + fraction*right[1]
    return ordered[-1][1]

def _score_external(case_id: str, case: Dict[str, Any],
                    model_rows: Sequence[Dict[str, str]],
                    reference_rows: Sequence[Dict[str, str]]) -> List[Dict[str, Any]]:
    series_key = "series" if case_id == "XM02" else "observable"
    x_key = "elapsed_hours" if case_id == "XM02" else "coordinate"
    y_key = "intensity" if case_id == "XM02" else "value"
    references = sorted({row[series_key] for row in reference_rows})
    models = {row[series_key] for row in model_rows}
    coverage = len(set(references) & models)/len(references)
    a = case["acceptance"]
    log_residuals: List[float] = []
    peak_errors: List[float] = []
    onset_errors: List[float] = []
    slope_error = 0.0
    mfp_residuals: List[float] = []
    for name in references:
        if name not in models: continue
        ref = sorted((float(row[x_key]), float(row[y_key])) for row in reference_rows
                     if row[series_key] == name)
        mod = sorted((float(row[x_key]), float(row[y_key])) for row in model_rows
                     if row[series_key] == name)
        if case_id == "XM03" and name == "earth_fluence_spectral_index":
            slope_error = abs(mod[0][1]-ref[0][1]); continue
        residual_target = mfp_residuals if "mean_free_path" in name else log_residuals
        for x, expected in ref:
            observed = _interpolate(mod, x)
            if observed > 0.0 and expected > 0.0:
                residual_target.append(math.log10(observed/expected))
        if "intensity" in name or case_id == "XM02":
            peak_errors.append(abs(max(mod, key=lambda item: item[1])[0] -
                                   max(ref, key=lambda item: item[1])[0]))
            if case_id == "XM03":
                # A relative threshold avoids embedding an instrument-specific
                # absolute background in a model-to-model test. Interpolation
                # remains intentionally absent: the timing resolution is the
                # publication trace cadence and is visible in the metric.
                fraction = float(a["onset_fraction_of_peak"])
                ref_threshold = fraction*max(value for _, value in ref)
                mod_threshold = fraction*max(value for _, value in mod)
                ref_onset = next(x for x, value in ref if value >= ref_threshold)
                mod_onset = next(x for x, value in mod if value >= mod_threshold)
                onset_errors.append(abs(mod_onset-ref_onset))
    metrics = [metric("required_series_fraction", coverage,
                      a["required_series_fraction_min"], ">=", "fraction")]
    if case_id == "XM02":
        metrics += [
            metric("log10_intensity_rmse", _rmse(log_residuals),
                   a["log10_intensity_rmse_max"], "<=", "dex"),
            metric("peak_timing_error", max(peak_errors, default=math.inf),
                   a["peak_timing_error_hours_max"], "<=", "hours")]
    else:
        metrics += [
            metric("log10_intensity_rmse", _rmse(log_residuals),
                   a["log10_intensity_rmse_max"], "<=", "dex"),
            metric("onset_timing_error", max(onset_errors, default=math.inf),
                   a["onset_timing_error_hours_max"], "<=", "hours"),
            metric("peak_timing_error", max(peak_errors, default=math.inf),
                   a["peak_timing_error_hours_max"], "<=", "hours"),
            metric("fluence_spectral_index_error", slope_error,
                   a["fluence_spectral_index_error_max"], "<=", "index"),
            metric("log10_mean_free_path_rmse", _rmse(mfp_residuals),
                   a["log10_mean_free_path_rmse_max"], "<=", "dex")]
    return metrics

def _skip_result(case_id: str, started: float, case: Dict[str, Any],
                 input_path: Path, executable: Path, artifacts: Sequence[Path],
                 reason: str) -> Dict[str, Any]:
    return {"id": case_id, "status": "SKIP", "message": reason,
            "elapsed_seconds": time.monotonic()-started,
            "seed": int(case["numerics"]["campaign_seed"]),
            "configuration": [f"input={input_path}", "execution=linked-srcsep-amps",
                f"executable={executable}", f"executable_sha256={sha256(executable)}",
                "reference=digitized-published-M-FLAMPA",
                f"equivalence_reviewed={bool(case.get('equivalence_reviewed', False))}"],
            "metrics": [], "artifacts": [str(path) for path in artifacts]}

def run_cross_model_case(case_id: str, *, source_root: Path, input_path: Path,
                         output_dir: Path, executable: Path,
                         timeout: Optional[float]) -> Dict[str, Any]:
    """Run one XM case and return registry-compatible evidence."""
    started = time.monotonic()
    case = load_input(input_path, case_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = output_dir / "resolved_input.json"
    atomic_json(resolved, case)
    formats = _formats(case)

    if case_id == "XM01":
        native = run_linked_model(case_id=case_id, arguments=_xm01_arguments(case),
            source_root=source_root, output_dir=output_dir, executable=executable,
            timeout=timeout)
        reference = output_dir / "XM01_reference.csv"
        script = Path(__file__).parent / "XM01" / "reference_solution.py"
        completed = subprocess.run([sys.executable, str(script), "--input", str(resolved),
            "--output", str(reference)], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=timeout)
        if completed.returncode:
            raise RuntimeError(f"XM01 independent solver failed: {completed.stdout}")
        model_rows, reference_rows = read_csv(native["model"]), read_csv(reference)
        metrics = _score_xm01(case, model_rows, reference_rows)
        figures = _xm01_plot(output_dir, model_rows, reference_rows, formats)
        provenance = output_dir / "provenance.json"
        atomic_json(provenance, {"schema": "srcsep-validation-provenance-v1",
            "case_id": case_id, "model": "linked production focused-transport core",
            "reference": "independent conservative Python finite-volume solver",
            "sha256": {"executable": sha256(executable), "input": sha256(resolved),
                       "model": sha256(native["model"]), "reference": sha256(reference)}})
        artifacts = [resolved, native["manifest"], native["model"], native["report"],
                     native["junit"], native["log"], reference, provenance] + figures
        return finish_result(case_id=case_id, started=started,
            seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
            executable=executable, metrics=metrics, artifacts=artifacts,
            message="XM01 linked cross-solver validation passed")

    immutable = _reference_path(case_id, case)
    reference = output_dir / f"{case_id}_reference.csv"
    shutil.copyfile(immutable, reference)
    reference_rows = read_csv(reference)
    publication = Path(__file__).parent / case_id / "reference" / "provenance.json"
    artifacts: List[Path] = [resolved, reference]
    if publication.is_file(): artifacts.append(publication)
    source_text = str(case.get("model_source_csv", "")).strip()
    if not source_text:
        artifacts += _external_plot(case_id, output_dir, reference_rows, None, formats)
        return _skip_result(case_id, started, case, input_path, executable, artifacts,
            f"{case_id} reference prepared; set model_source_csv to a production srcSEP export")
    source = Path(source_text).expanduser()
    if not source.is_absolute(): source = (input_path.parent/source).resolve()
    if not source.is_file():
        artifacts += _external_plot(case_id, output_dir, reference_rows, None, formats)
        return _skip_result(case_id, started, case, input_path, executable, artifacts,
                            f"{case_id} model_source_csv does not exist: {source}")
    native = run_linked_model(case_id=case_id,
        arguments=["--source-csv", str(source), "--campaign-seed",
                   str(case["numerics"]["campaign_seed"])],
        source_root=source_root, output_dir=output_dir, executable=executable,
        timeout=timeout)
    model_rows = read_csv(native["model"])
    metrics = _score_external(case_id, case, model_rows, reference_rows)
    figures = _external_plot(case_id, output_dir, reference_rows, model_rows, formats)
    provenance = output_dir / "provenance.json"
    atomic_json(provenance, {"schema": "srcsep-validation-provenance-v1",
        "case_id": case_id, "publication_id": case["reference"]["publication_id"],
        "equivalence_reviewed": bool(case.get("equivalence_reviewed", False)),
        "sha256": {"executable": sha256(executable), "model_source": sha256(source),
                   "normalized_model": sha256(native["model"]), "reference": sha256(reference)}})
    artifacts += [native["manifest"], native["model"], native["report"], native["junit"],
                  native["log"], provenance] + figures
    if not bool(case.get("equivalence_reviewed", False)):
        skipped = _skip_result(case_id, started, case, input_path, executable, artifacts,
            f"{case_id} metrics computed, but assumption equivalence has not been reviewed")
        skipped["metrics"] = metrics
        return skipped
    return finish_result(case_id=case_id, started=started,
        seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
        executable=executable, metrics=metrics, artifacts=artifacts,
        message=f"{case_id} linked M-FLAMPA comparison passed")
