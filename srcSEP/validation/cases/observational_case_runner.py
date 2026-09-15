"""Application-level orchestration for observational cases OV01--OV05.

The module deliberately separates three things that are easy to conflate in
an event study: immutable measurements digitized from a named publication
figure, a reviewed but necessarily reduced heliospheric input, and output
created by the linked ``srcSEP/AMPS`` executable in the current invocation.
Only the linked executable advances particles.  Python resolves inputs,
preserves provenance, applies declared nuisance normalizations, and renders
the comparison evidence.
"""
from __future__ import annotations

import math
import shutil
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from linked_case_common import (atomic_json, finish_result, load_input, metric,
                                read_csv, run_linked_model, sha256, write_csv)
from cross_model_case_runner import (_formats, _publication_input_path,
    _publication_plot_label, _reference_path, _rmse, _xm03_arguments,
    _xm03_plot, _xm03_score)


def _case_file(input_path: Path, relative_text: str, purpose: str) -> Path:
    """Resolve one immutable case-owned file and reject directory traversal."""
    relative = Path(relative_text)
    case_directory = input_path.parent.resolve()
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"{purpose} must be a case-relative path")
    resolved = (case_directory / relative).resolve()
    if case_directory not in resolved.parents or not resolved.is_file():
        raise ValueError(f"{purpose} is missing or outside the case: {resolved}")
    return resolved


def _ov_arguments(case: Dict[str, Any], input_path: Path,
                  run: Dict[str, Any]) -> Tuple[List[str], Path]:
    """Translate a reviewed observer/hypothesis into native SI arguments.

    Event-wide transport coefficients live under ``physics`` and are reused
    unchanged by all runs.  The run table may vary only geometry, connection
    delay, requested channels, and source-history identity.  This structure is
    what prevents a radial comparison from silently fitting a separate mean
    free path or injection index to each spacecraft.
    """
    physics, numerics = case["physics"], case["numerics"]
    source = _case_file(input_path, str(run["source_history_csv"]),
                        "source_history_csv")
    energies = run["energies_mev"]
    if (not isinstance(energies, list) or not energies or
            any(float(value) <= 0.0 for value in energies)):
        raise ValueError("each OV run requires positive energies_mev")
    output_mode = str(run["output_mode"])
    if output_mode not in ("time-profile", "integrated-spectrum"):
        raise ValueError("OV output_mode must be time-profile or integrated-spectrum")
    arguments = [
        "--observer", str(run["name"]),
        "--output-mode", output_mode,
        "--source-radius-mode", str(physics["source_radius_mode"]),
        "--energies-mev", ",".join(str(float(value)) for value in energies),
        "--source-history-csv", str(source),
        "--injection-radius-solar-radii", str(physics["injection_radius_solar_radii"]),
        "--observer-radius-au", str(run["observer_radius_au"]),
        "--solar-wind-speed-m-per-s", str(run["solar_wind_speed_m_per_s"]),
        "--shock-speed-m-per-s", str(physics["shock_speed_m_per_s"]),
        "--solar-rotation-rate-rad-per-s", str(physics["solar_rotation_rate_rad_per_s"]),
        "--launch-offset-s", str(run.get("launch_offset_s", 0.0)),
        "--connection-delay-s", str(run.get("connection_delay_s", 0.0)),
        "--mfp-normalization-au", str(physics["mean_free_path_normalization_au"]),
        "--mfp-radial-exponent", str(physics["mean_free_path_radial_exponent"]),
        "--mfp-rigidity-exponent", str(physics["mean_free_path_rigidity_exponent"]),
        "--injection-momentum-index", str(physics["injection_momentum_index"]),
        "--particles-per-energy", str(numerics["particles_per_energy"]),
        "--time-step-s", str(numerics["time_step_s"]),
        "--duration-s", str(3600.0 * float(run["duration_hours"])),
        "--output-cadence-s", str(3600.0 * float(run["output_cadence_hours"])),
        # Derive disjoint, reproducible streams from the case seed and stable
        # run index.  Results therefore do not depend on registry ordering.
        "--campaign-seed", str(int(numerics["campaign_seed"]) + int(run["seed_offset"])),
    ]
    return arguments, source


def _interpolate(points: Sequence[Tuple[float, float]], coordinate: float) -> Optional[float]:
    """Interpolate positive traces logarithmically without inventing coverage."""
    ordered = sorted(points)
    if len(ordered) < 2 or coordinate < ordered[0][0] or coordinate > ordered[-1][0]:
        return None
    for left, right in zip(ordered[:-1], ordered[1:]):
        if left[0] <= coordinate <= right[0]:
            if coordinate == left[0]:
                return left[1]
            fraction = (coordinate - left[0]) / (right[0] - left[0])
            if left[1] > 0.0 and right[1] > 0.0:
                return 10.0 ** ((1.0 - fraction) * math.log10(left[1]) +
                                fraction * math.log10(right[1]))
            # A modeled intensity can be exactly zero before first arrival.
            # Logarithmic interpolation is undefined across that boundary, so
            # retain the interval (and therefore honest reference coverage)
            # with a linear interpolation instead of accidentally falling
            # through to the final sample of the complete time series.
            return (1.0 - fraction) * left[1] + fraction * right[1]
    return ordered[-1][1]


def _score_profiles(case: Dict[str, Any], model_rows: Sequence[Dict[str, str]],
                    reference_rows: Sequence[Dict[str, str]]) -> Tuple[
                        List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, float]]:
    """Score all trace points under the normalization declared by the case.

    A single global amplitude tests relative spacecraft, time, and energy
    behavior together.  Unit-peak-per-series is allowed only for diagnostic
    cases where the publication does not provide a cross-calibrated absolute
    product; it preserves onset, peak time, rise, and decay shape.
    """
    acceptance = case["acceptance"]
    normalization = str(case["physics"]["comparison_normalization"])
    if normalization not in ("one-global-log-amplitude", "unit-peak-per-series"):
        raise ValueError("unsupported OV comparison_normalization")
    model_by_series: Dict[str, List[Tuple[float, float]]] = {}
    for row in model_rows:
        coordinate = (float(row["elapsed_hours"])
                      if "profile" in row["series"] or
                         not row["series"].endswith("_fluence_spectrum")
                      else float(row["energy_mev"]))
        model_by_series.setdefault(row["series"], []).append(
            (coordinate, float(row["relative_intensity"])))

    raw: List[Tuple[Dict[str, str], float]] = []
    for row in reference_rows:
        coordinate = (float(row["elapsed_hours"])
                      if row["comparison_mode"] == "time-profile"
                      else float(row["energy_mev"]))
        value = _interpolate(model_by_series.get(row["series"], []), coordinate)
        if value is not None and value > 0.0 and float(row["value"]) > 0.0:
            raw.append((row, value))

    scales: Dict[str, float] = {}
    if normalization == "one-global-log-amplitude" and raw:
        offsets = [math.log10(float(row["value"]) / value) for row, value in raw]
        scales["global"] = 10.0 ** (sum(offsets) / len(offsets))
    elif normalization == "unit-peak-per-series":
        for series in sorted({row["series"] for row, _ in raw}):
            expected_peak = max(float(row["value"]) for row, _ in raw
                                if row["series"] == series)
            model_peak = max(value for row, value in raw if row["series"] == series)
            scales[series] = expected_peak / model_peak

    comparison_rows: List[Dict[str, Any]] = []
    residuals: List[float] = []
    matched = {(row["series"], row["elapsed_hours"], row["energy_mev"]): value
               for row, value in raw}
    expected_logs: List[float] = []
    model_logs: List[float] = []
    for row in reference_rows:
        key = (row["series"], row["elapsed_hours"], row["energy_mev"])
        raw_value = matched.get(key)
        scale = scales.get("global", scales.get(row["series"], 1.0))
        scaled = None if raw_value is None else scale * raw_value
        residual = (None if scaled is None else
                    math.log10(scaled / float(row["value"])))
        if residual is not None:
            residuals.append(residual)
            expected_logs.append(math.log10(float(row["value"])))
            model_logs.append(math.log10(scaled))
        comparison_rows.append({**row,
            "model_relative": "" if raw_value is None else raw_value,
            "normalization_factor": scale,
            "model_scaled": "" if scaled is None else scaled,
            "log10_model_over_reference": "" if residual is None else residual})

    correlation = -1.0
    if len(expected_logs) >= 2:
        ex_mean = sum(expected_logs) / len(expected_logs)
        mo_mean = sum(model_logs) / len(model_logs)
        numerator = sum((x-ex_mean)*(y-mo_mean)
                        for x, y in zip(expected_logs, model_logs))
        denominator = math.sqrt(sum((x-ex_mean)**2 for x in expected_logs) *
                                sum((y-mo_mean)**2 for y in model_logs))
        if denominator > 0.0:
            correlation = numerator / denominator
    coverage = len(raw) / len(reference_rows) if reference_rows else 0.0
    release_gate = str(case["validation_role"]) == "release-gate"
    metrics = [
        metric("reference_point_coverage", coverage,
               float(acceptance["reference_point_coverage_min"]), ">=", "fraction"),
        metric("log10_profile_rmse", _rmse(residuals),
               float(acceptance["log10_profile_rmse_max"]), "<=", "dex",
               gating=release_gate),
        metric("log10_profile_correlation", correlation,
               float(acceptance["log10_profile_correlation_min"]), ">=", "correlation",
               gating=release_gate),
    ]
    return metrics, comparison_rows, scales


def _plot_profiles(case_id: str, case: Dict[str, Any], output: Path,
                   model_rows: Sequence[Dict[str, str]],
                   reference_rows: Sequence[Dict[str, str]],
                   scales: Dict[str, float], formats: Sequence[str]) -> List[Path]:
    """Render publication points and linked traces with citation in the image."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = sorted({row["series"] for row in reference_rows})
    figure, axes = plt.subplots(len(series), 1,
                                figsize=(9, max(3.2, 2.7*len(series))),
                                squeeze=False)
    for axis, name in zip(axes[:, 0], series):
        reference = [row for row in reference_rows if row["series"] == name]
        mode = reference[0]["comparison_mode"]
        x_key = "elapsed_hours" if mode == "time-profile" else "energy_mev"
        reference.sort(key=lambda row: float(row[x_key]))
        model = [row for row in model_rows if row["series"] == name]
        model.sort(key=lambda row: float(row[x_key]))
        scale = scales.get("global", scales.get(name, 1.0))
        axis.plot([float(row[x_key]) for row in reference],
                  [float(row["value"]) for row in reference], "ko",
                  ms=3.5, label=f"{reference[0]['instrument']} digitized reference")
        axis.plot([float(row[x_key]) for row in model],
                  [scale*float(row["relative_intensity"]) for row in model],
                  color="#d95f02", linewidth=1.25, label="linked srcSEP/AMPS")
        if mode == "integrated-spectrum":
            axis.set_xscale("log")
            axis.set_xlabel("proton energy [MeV]")
        else:
            axis.set_xlabel("elapsed time from case epoch [h]")
        axis.set_yscale("log")
        axis.set_ylabel("intensity/fluence")
        axis.set_title(name.replace("_", " "))
        axis.grid(color="0.86", linewidth=0.6)
        axis.legend(framealpha=1.0)
    role = str(case["validation_role"]).replace("-", " ")
    figure.suptitle(f"{case_id} observational comparison ({role})\n"
                    f"{_publication_plot_label(case)}", fontsize=10)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    paths: List[Path] = []
    for extension in formats:
        path = output / f"{case_id}_observational_comparison.{extension}"
        figure.savefig(path, dpi=180 if extension == "png" else None)
        paths.append(path)
    plt.close(figure)
    return paths


def _run_ov01(*, source_root: Path, input_path: Path, output_dir: Path,
              executable: Path, timeout: Optional[float]) -> Dict[str, Any]:
    """Apply the dedicated near-Earth release gate using the XM03 kernel."""
    started = time.monotonic()
    case = load_input(input_path, "OV01")
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = output_dir / "resolved_input.json"
    atomic_json(resolved, case)
    # OV01 and XM03 intentionally point at the same reviewed measurement
    # table: OV01 applies the campaign's release criteria, while XM03 retains
    # the cross-model/event-reconstruction classification.  Resolve that
    # declared source case explicitly instead of copying 80 observations and
    # allowing the two baselines to drift independently.
    source_case = str(case["reference"].get("source_case", "OV01"))
    immutable = _reference_path(source_case, case)
    reference = output_dir / "OV01_reference.csv"
    shutil.copyfile(immutable, reference)
    reference_rows = read_csv(reference)
    arguments = _xm03_arguments(case, input_path)
    source = Path(arguments[arguments.index("--source-history-csv") + 1])
    native = run_linked_model(case_id="OV01", arguments=arguments,
        source_root=source_root, output_dir=output_dir, executable=executable,
        timeout=timeout)
    model_rows = read_csv(native["model"])
    metrics, scale, rows = _xm03_score(case, model_rows, reference_rows)
    comparison = output_dir / "OV01_observation_comparison.csv"
    write_csv(comparison, (
        "elapsed_hours", "instrument", "energy_low_mev", "energy_high_mev",
        "effective_energy_mev", "differential_intensity_pfu_per_mev", "scored",
        "model_relative", "model_scaled_pfu_per_mev",
        "log10_model_over_observation"), rows)
    figures = _xm03_plot(case, output_dir, model_rows, reference_rows, scale,
                         _formats(case), case_id="OV01")
    provenance = output_dir / "provenance.json"
    atomic_json(provenance, {"schema": "srcsep-validation-provenance-v1",
        "case_id": "OV01", "validation_role": "release-gate",
        "publication_id": case["reference"]["publication_id"],
        "comparison_normalization": {"kind": "one-global-log-amplitude",
                                     "factor": scale},
        "sha256": {"executable": sha256(executable), "input": sha256(resolved),
                   "model": sha256(native["model"]),
                   "reference": sha256(reference), "source": sha256(source)}})
    publication_input = _publication_input_path("OV01", case)
    artifacts = [resolved, reference, source, native["manifest"], native["model"],
                 native["report"], native["junit"], native["log"], comparison,
                 provenance] + figures
    if publication_input is not None:
        artifacts.append(publication_input)
    return finish_result(case_id="OV01", started=started,
        seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
        executable=executable, metrics=metrics, artifacts=artifacts,
        message="OV01 linked near-Earth release validation passed")


def run_observational_case(case_id: str, *, source_root: Path, input_path: Path,
                           output_dir: Path, executable: Path,
                           timeout: Optional[float]) -> Dict[str, Any]:
    """Run one fixed OV case and return registry-compatible evidence."""
    if case_id == "OV01":
        return _run_ov01(source_root=source_root, input_path=input_path,
                         output_dir=output_dir, executable=executable,
                         timeout=timeout)
    started = time.monotonic()
    case = load_input(input_path, case_id)
    if case.get("publication_input_only") is not True:
        raise ValueError(f"{case_id} must declare publication_input_only=true")
    if case.get("validation_role") not in ("release-gate", "diagnostic-only"):
        raise ValueError(f"{case_id} has an invalid validation_role")
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = output_dir / "resolved_input.json"
    atomic_json(resolved, case)
    immutable = _reference_path(case_id, case)
    reference = output_dir / f"{case_id}_reference.csv"
    shutil.copyfile(immutable, reference)
    reference_rows = read_csv(reference)

    model_rows: List[Dict[str, str]] = []
    artifacts: List[Path] = [resolved, reference]
    source_hashes: Dict[str, str] = {}
    for run in case["model_runs"]:
        arguments, source = _ov_arguments(case, input_path, run)
        native = run_linked_model(case_id=case_id, arguments=arguments,
            source_root=source_root, output_dir=output_dir / f"native_{run['name']}",
            executable=executable, timeout=timeout)
        model_rows.extend(read_csv(native["model"]))
        source_hashes[str(run["name"])] = sha256(source)
        artifacts.extend([source, native["manifest"], native["model"],
                          native["report"], native["junit"], native["log"]])

    combined = output_dir / f"{case_id}_linked_model.csv"
    write_csv(combined, ("observer", "series", "elapsed_hours", "energy_mev",
                         "relative_intensity", "effective_sample_count"), model_rows)
    metrics, comparison_rows, scales = _score_profiles(
        case, model_rows, reference_rows)
    comparison = output_dir / f"{case_id}_comparison.csv"
    write_csv(comparison, ("series", "observer", "instrument", "figure",
        "comparison_mode", "elapsed_hours", "energy_mev", "value",
        "digitization_relative_uncertainty", "model_relative",
        "normalization_factor", "model_scaled", "log10_model_over_reference"),
        comparison_rows)
    figures = _plot_profiles(case_id, case, output_dir, model_rows,
                             reference_rows, scales, _formats(case))
    publication_manifest = _publication_input_path(case_id, case)
    if publication_manifest is not None:
        artifacts.append(publication_manifest)
    publication_provenance = input_path.parent / "reference" / "provenance.json"
    if publication_provenance.is_file():
        artifacts.append(publication_provenance)
    provenance = output_dir / "provenance.json"
    atomic_json(provenance, {"schema": "srcsep-validation-provenance-v1",
        "case_id": case_id, "validation_role": case["validation_role"],
        "publication_id": case["reference"]["publication_id"],
        "normalization": case["physics"]["comparison_normalization"],
        "normalization_factors": scales,
        "sha256": {"executable": sha256(executable), "input": sha256(resolved),
                   "model": sha256(combined), "reference": sha256(reference),
                   "sources": source_hashes}})
    artifacts.extend([combined, comparison, provenance] + figures)
    role = "release validation" if case["validation_role"] == "release-gate" \
        else "diagnostic evidence generation"
    return finish_result(case_id=case_id, started=started,
        seed=int(case["numerics"]["campaign_seed"]), input_path=input_path,
        executable=executable, metrics=metrics, artifacts=artifacts,
        message=f"{case_id} linked observational {role} passed")
