"""Linked-application orchestration and scoring for XM01-XM03.

XM01 compares production stochastic characteristics with an independent
finite-volume PDE solver. XM02 runs a linked, publication-informed controlled
transport reconstruction. XM03 compares a normalized production export with
immutable M-FLAMPA references. XM02/XM03 accept only their registered
publication inputs, so the former command-line equivalence-review input branch
is intentionally absent. This distinction is part of the evidence contract.
"""
from __future__ import annotations

import csv
import json
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

def _xm02_arguments(case: Dict[str, Any]) -> List[str]:
    """Translate the sole XM02 reconstruction into native SI arguments.

    Zhao et al. report the MFP ensemble and field-line seed radius but not the
    time-dependent shock/field-line state required to replay Figure 7.  The
    registered input therefore defines a controlled first-passage experiment:
    already-accelerated >10 MeV protons are released through a fixed causal
    source at 2.5 solar radii and transported to 1 AU with the production core.
    Every additional assumption is explicit in input.json and is serialized
    here; no command-line input or precomputed model CSV is accepted.
    """
    physics, numerics = case["physics"], case["numerics"]
    mean_free_paths = physics["far_upstream_mean_free_paths_au"]
    if not isinstance(mean_free_paths, list) or len(mean_free_paths) != 3:
        raise ValueError("XM02 requires exactly three mean-free-path values")
    if physics.get("comparison_normalization") != "unit_peak_per_series":
        raise ValueError("XM02 requires unit_peak_per_series comparison")
    fixed_choices = {
        "source_time_profile": "two-stage-exponential-release",
        "initial_pitch_distribution": "outward-isotropic-flux",
        "inner_boundary": "reflecting",
        "outer_boundary": "first-passage-observer",
    }
    for name, expected in fixed_choices.items():
        if physics.get(name) != expected:
            raise ValueError(f"XM02 {name} must be {expected}")
    return [
        "--injection-radius-solar-radii", str(physics["injection_radius_solar_radii"]),
        "--observer-radius-au", str(physics["observer_radius_au"]),
        "--particle-energy-mev", str(physics["transport_particle_energy_mev"]),
        "--mfp-0-au", str(mean_free_paths[0]),
        "--mfp-1-au", str(mean_free_paths[1]),
        "--mfp-2-au", str(mean_free_paths[2]),
        "--plasma-advection-m-per-s", str(physics["plasma_advection_m_per_s"]),
        "--dlnb-ds-per-m", str(physics["magnetic_focusing_per_m"]),
        "--source-rise-time-s", str(3600.0 * float(physics["source_rise_time_hours"])),
        "--source-decay-time-s", str(3600.0 * float(physics["source_decay_time_hours"])),
        "--particles", str(numerics["particles_per_mean_free_path"]),
        "--time-step-s", str(numerics["time_step_s"]),
        "--duration-s", str(3600.0 * float(numerics["duration_hours"])),
        "--output-cadence-s", str(3600.0 * float(numerics["output_cadence_hours"])),
        "--campaign-seed", str(numerics["campaign_seed"]),
    ]

def _reference_path(case_id: str, case: Dict[str, Any]) -> Path:
    """Resolve a source-owned immutable reference without accepting traversal.

    Other validation cases may accept a custom input file. Anchoring this
    baseline at the registered XM directory prevents any such path from
    silently replacing the reviewed reference.
    """
    relative = Path(str(case["reference"]["csv"]))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("reference.csv must be a case-relative path")
    case_directory = (Path(__file__).parent / case_id).resolve()
    resolved = (case_directory / relative).resolve()
    if case_directory not in resolved.parents or not resolved.is_file():
        raise ValueError(f"reference CSV is missing or outside the case: {resolved}")
    return resolved

def _publication_input_path(case_id: str, case: Dict[str, Any]) -> Optional[Path]:
    """Resolve and validate the publication-derived input reconstruction.

    XM02 and XM03 are comparisons with published calculations, but neither
    article supplies a complete SWMF run directory.  The reconstruction file
    therefore records every reported parameter together with assumptions and
    missing artifacts.  Keeping that record machine-readable makes the limits
    of reproducibility visible to test automation and prevents a digitized
    curve from being mistaken for a fully reproducible model setup.

    The path is constrained to the registered case directory for the same
    reason as the immutable reference CSV: a path embedded in the registered
    case may select a production result, but it must not silently replace
    source-reviewed publication evidence.
    """
    relative_text = str(case.get("publication_input_manifest", "")).strip()
    if not relative_text:
        return None
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("publication_input_manifest must be a case-relative path")
    case_directory = (Path(__file__).parent / case_id).resolve()
    resolved = (case_directory / relative).resolve()
    if case_directory not in resolved.parents or not resolved.is_file():
        raise ValueError(
            f"publication input manifest is missing or outside the case: {resolved}")

    # Validate the small set of fields on which the evidence contract relies.
    # We intentionally do not require a particular list of physics parameters:
    # future publications may expose more (or fewer) fields, while the declared
    # reconstruction status and missing-input list remain stable concepts.
    # This is deliberately a different schema from the executable case input,
    # so parse it directly rather than routing it through load_input(), whose
    # strict schema guard protects normal linked-case configurations.
    with resolved.open("r", encoding="utf-8") as stream:
        reconstruction = json.load(stream)
    if not isinstance(reconstruction, dict) or reconstruction.get("case_id") != case_id:
        raise ValueError(f"publication input manifest must use case_id={case_id}")
    if reconstruction.get("schema") != "srcsep-publication-input-reconstruction-v1":
        raise ValueError("publication input manifest has an unsupported schema")
    if reconstruction.get("reproduction_status") not in ("partial", "complete"):
        raise ValueError("publication input manifest must declare partial or complete")
    missing = reconstruction.get("missing_required_inputs")
    if not isinstance(missing, list):
        raise ValueError("publication input manifest must list missing_required_inputs")
    if reconstruction["reproduction_status"] == "complete" and missing:
        raise ValueError("a complete publication input manifest cannot list missing inputs")
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
        value_key = "intensity" if case_id == "XM02" else "value"
        reference_values = [float(row[value_key]) for row in ref]
        model_values = ([float(row[value_key]) for row in sorted(
            (row for row in model_rows if row[series_key] == name),
            key=lambda row: float(row[x_key]))] if model_rows is not None else None)
        mod = (sorted((row for row in model_rows if row[series_key] == name),
                      key=lambda row: float(row[x_key]))
               if model_rows is not None else [])
        if case_id == "XM02":
            # Zhao et al. do not provide the plasma/shock information required
            # to turn the injection coefficient into absolute pfu.  Compare
            # unit-peak shapes, retaining timing and decay sensitivity without
            # fitting a normalization to the digitized curve.
            reference_peak = max(reference_values)
            reference_values = [value / reference_peak for value in reference_values]
            if model_values is not None:
                model_peak = max(model_values)
                model_values = [value / model_peak for value in model_values]
        axis.plot([float(row[x_key]) for row in ref], reference_values,
                  "k.-", label="digitized M-FLAMPA reference")
        if model_rows is not None:
            assert model_values is not None
            model_label = ("linked controlled srcSEP reconstruction"
                           if case_id == "XM02"
                           else "linked srcSEP/AMPS export")
            axis.plot([float(row[x_key]) for row in mod], model_values,
                      "C1o-", ms=3, label=model_label)
        if all(value > 0.0 for value in reference_values):
            axis.set_yscale("log")
        axis.set_title(name.replace("_", " "))
        axis.set_ylabel("unit-peak intensity" if case_id == "XM02" else "value")
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
        if case_id == "XM02":
            # Score transport shape rather than an unknowable absolute source
            # normalization. Peak division is performed independently for the
            # two curves and does not alter their peak coordinates.
            ref_peak = max(value for _, value in ref)
            mod_peak = max(value for _, value in mod)
            if not (ref_peak > 0.0 and mod_peak > 0.0):
                continue
            ref = [(x, value / ref_peak) for x, value in ref]
            mod = [(x, value / mod_peak) for x, value in mod]
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
            metric("unit_peak_log10_intensity_rmse", _rmse(log_residuals),
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
                "input_policy=single-registered-publication-reconstruction"],
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
    if case.get("publication_input_only") is not True:
        # XM02/XM03 no longer accept an operator-selected input variant.  This
        # marker is part of their registered default input and proves that the
        # case intends to use the sole literature reconstruction shipped with
        # the source.  Refuse older/copied configurations instead of quietly
        # reviving the ambiguous equivalence-review workflow.
        raise ValueError(
            f"{case_id} input must declare publication_input_only=true")
    reference = output_dir / f"{case_id}_reference.csv"
    shutil.copyfile(immutable, reference)
    reference_rows = read_csv(reference)
    publication_input_source = _publication_input_path(case_id, case)
    publication_input: Optional[Path] = None
    if publication_input_source is not None:
        # Copy the exact reconstruction used by this invocation into the result
        # directory.  Evidence bundles then remain self-describing even when
        # detached from the source checkout that launched the validation run.
        publication_input = output_dir / f"{case_id}_publication_input.json"
        shutil.copyfile(publication_input_source, publication_input)
    publication = Path(__file__).parent / case_id / "reference" / "provenance.json"
    artifacts: List[Path] = [resolved, reference]
    if publication_input is not None: artifacts.append(publication_input)
    if publication.is_file(): artifacts.append(publication)

    if case_id == "XM02":
        # XM02 now generates its numerical solution inside the selected linked
        # application.  The sole input is the registry-owned reconstruction;
        # there is no model_source_csv and therefore no external-input branch.
        native = run_linked_model(case_id=case_id, arguments=_xm02_arguments(case),
            source_root=source_root, output_dir=output_dir, executable=executable,
            timeout=timeout)
        model_rows = read_csv(native["model"])
        metrics = _score_external(case_id, case, model_rows, reference_rows)
        figures = _external_plot(case_id, output_dir, reference_rows,
                                 model_rows, formats)
        provenance = output_dir / "provenance.json"
        atomic_json(provenance, {
            "schema": "srcsep-validation-provenance-v1",
            "case_id": case_id,
            "publication_id": case["reference"]["publication_id"],
            "model": "linked controlled first-passage reconstruction",
            "comparison_normalization": "unit_peak_per_series",
            "reproduction_status": "publication-informed-controlled-benchmark",
            "sha256": {
                "executable": sha256(executable),
                "normalized_model": sha256(native["model"]),
                "reference": sha256(reference),
                **({"publication_input": sha256(publication_input)}
                   if publication_input is not None else {})}})
        artifacts += [native["manifest"], native["model"], native["report"],
                      native["junit"], native["log"], provenance] + figures
        return finish_result(case_id=case_id, started=started,
            seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
            executable=executable, metrics=metrics, artifacts=artifacts,
            message="XM02 linked publication-informed transport comparison passed")

    source_text = str(case.get("model_source_csv", "")).strip()
    if not source_text:
        artifacts += _external_plot(case_id, output_dir, reference_rows, None, formats)
        return _skip_result(case_id, started, case, input_path, executable, artifacts,
            f"{case_id} reference prepared; registered production-result path is empty")
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
        "input_policy": "single-registered-publication-reconstruction",
        "sha256": {"executable": sha256(executable), "model_source": sha256(source),
                   "normalized_model": sha256(native["model"]), "reference": sha256(reference),
                   # A missing key means an older case definition did not
                   # provide a reconstruction; an empty string is avoided so
                   # downstream provenance readers cannot confuse it with a
                   # real digest.
                   **({"publication_input": sha256(publication_input)}
                      if publication_input is not None else {})}})
    artifacts += [native["manifest"], native["model"], native["report"], native["junit"],
                  native["log"], provenance] + figures
    return finish_result(case_id=case_id, started=started,
        seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
        executable=executable, metrics=metrics, artifacts=artifacts,
        message=f"{case_id} linked M-FLAMPA comparison passed")
