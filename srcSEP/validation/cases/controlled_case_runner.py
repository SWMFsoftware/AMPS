"""Case-specific scoring and figures for linked validation cases CV02-CV05.

The C++ executable produces only numerical evidence.  This module invokes the
case-local independent reference program, reloads both CSVs, computes physical
acceptance metrics, and renders PNG/EPS figures from those saved artifacts.
No function here can replace the required linked ``amps --test CVxx`` stage.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from linked_case_common import (atomic_json, finish_result, load_input, metric,
                                read_csv, run_linked_model, sha256, write_csv)


def _list(values: Iterable[Any]) -> str:
    return ",".join(str(value) for value in values)


def _reference(case_id: str, input_path: Path, output_path: Path,
               timeout: Optional[float]) -> None:
    """Run the independent reference as a separate auditable Python process."""
    script = Path(__file__).resolve().parent / case_id / "reference_solution.py"
    completed = subprocess.run(
        [sys.executable, str(script), "--input", str(input_path), "--output",
         str(output_path)], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False, timeout=timeout)
    if completed.returncode:
        raise RuntimeError(f"{case_id} reference failed: {completed.stdout.strip()}")


def _plot(case_id: str, output_dir: Path, title: str,
          x: Sequence[float], model: Sequence[float],
          reference: Sequence[float], xlabel: str, ylabel: str) -> List[Path]:
    """Render identical PNG/EPS comparison figures from the solution CSV data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})
    axes[0].plot(x, reference, "k-", linewidth=2.0, label="independent reference")
    axes[0].plot(x, model, "o", markersize=4.0, label="linked srcSEP/AMPS")
    axes[0].set_ylabel(ylabel)
    axes[0].set_title(f"{case_id}: {title}")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    residual = [a - b for a, b in zip(model, reference)]
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[1].plot(x, residual, "o-", markersize=3.0)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("model-reference")
    axes[1].grid(alpha=0.25)
    figure.tight_layout()
    outputs = []
    for suffix in ("png", "eps"):
        path = output_dir / f"{case_id}_comparison.{suffix}"
        figure.savefig(path, dpi=180 if suffix == "png" else None)
        outputs.append(path)
    plt.close(figure)
    return outputs


def _order(errors: Dict[float, float]) -> float:
    """Return the endpoint refinement slope; zero errors denote exactness."""
    steps = sorted(errors, reverse=True)
    if len(steps) < 2:
        return math.nan
    coarse, fine = steps[0], steps[-1]
    if errors[fine] == 0.0:
        return 99.0
    if errors[coarse] == 0.0:
        return -99.0
    return math.log(errors[coarse] / errors[fine]) / math.log(coarse / fine)


def _provenance(case_id: str, input_path: Path, executable: Path,
                model: Path, reference: Path, output_dir: Path) -> Path:
    path = output_dir / "provenance.json"
    atomic_json(path, {
        "schema": "srcsep-validation-provenance-v1", "case_id": case_id,
        "execution": "linked-srcsep-amps", "executable": str(executable),
        "sha256": {"executable": sha256(executable), "input": sha256(input_path),
                   "model": sha256(model), "reference": sha256(reference)},
        "reference_implementation": str(Path(__file__).resolve().parent /
                                        case_id / "reference_solution.py"),
    })
    return path


def _cv02(config: Dict[str, Any], rows: List[Dict[str, str]],
          reference_rows: List[Dict[str, str]], output_dir: Path
          ) -> Tuple[List[Dict[str, Any]], List[Path]]:
    physics, domain = config["physics"], config["domain"]
    numerics, acceptance = config["numerics"], config["acceptance"]
    kappa = float(physics["kappa_parallel_m2_per_s"])
    center, final = float(domain["packet_center_m"]), float(numerics["final_time_s"])
    finest, largest = min(map(float, numerics["time_steps_s"])), max(map(int, numerics["particle_counts"]))
    seed_count = int(numerics["seed_count"])
    moments = [row for row in rows if row["row_type"] == "moment" and
               float(row["dt_s"]) == finest and int(row["particle_count"]) == largest]
    # Fit variance through the known zero-variance initial condition. Pool all
    # saved times/seeds; this is less biased than using one noisy endpoint.
    numerator = sum(float(row["time_s"]) * float(row["variance_m2"]) for row in moments)
    denominator = sum(float(row["time_s"]) ** 2 for row in moments)
    kappa_fit = 0.5 * numerator / denominator
    final_moments = [row for row in moments if float(row["time_s"]) == final]
    mean = sum(float(row["mean_m"]) for row in final_moments) / len(final_moments)
    mean_se = math.sqrt(2.0 * kappa * final / (largest * seed_count))
    mean_z = abs(mean - center) / mean_se
    escaped = max(float(row["escaped_weight"]) for row in moments)
    profile_rows = [row for row in rows if row["row_type"] == "profile" and
                    float(row["dt_s"]) == finest and int(row["particle_count"]) == largest]
    by_bin: Dict[int, List[float]] = {}
    for row in profile_rows:
        by_bin.setdefault(int(row["bin_index"]), []).append(float(row["probability"]))
    model = [sum(by_bin[i]) / len(by_bin[i]) for i in sorted(by_bin)]
    expected = [float(row["probability"]) for row in reference_rows]
    l1 = sum(abs(a - b) for a, b in zip(model, expected))
    l2 = math.sqrt(sum((a - b) ** 2 for a, b in zip(model, expected)))
    linf = max(abs(a - b) for a, b in zip(model, expected))
    mean_skew = sum(float(row["skewness"]) for row in final_moments) / len(final_moments)
    mean_kurtosis = sum(float(row["kurtosis"]) for row in final_moments) / len(final_moments)
    wrong = []
    for row in reference_rows:
        # The negative control changes only kappa and evaluates the same exact
        # bin integral locally; it never calls the production mover.
        scale = math.sqrt(4.0 * 1.5 * kappa * final)
        left, right = float(row["bin_left_m"]), float(row["bin_right_m"])
        wrong.append(0.5 * (math.erf((right - center) / scale) -
                            math.erf((left - center) / scale)))
    wrong_l1 = sum(abs(a - b) for a, b in zip(model, wrong))
    ratio = wrong_l1 / max(l1, 1.0e-30)
    metrics = [
        metric("kappa_fit_relative_error", abs(kappa_fit / kappa - 1.0),
               acceptance["kappa_relative_error_max"], "<=", "dimensionless"),
        metric("packet_mean_standard_errors", mean_z,
               acceptance["mean_standard_errors_max"], "<=", "standard-errors"),
        metric("final_profile_L1_error", l1, acceptance["profile_l1_max"], "<=", "probability"),
        metric("final_profile_L2_error", l2, acceptance["profile_l2_max"], "<=", "probability"),
        metric("final_profile_Linf_error", linf, acceptance["profile_linf_max"], "<=", "probability"),
        metric("final_skewness_absolute", abs(mean_skew), acceptance["skewness_absolute_max"], "<=", "dimensionless"),
        metric("final_kurtosis_absolute_error", abs(mean_kurtosis - 3.0), acceptance["kurtosis_absolute_error_max"], "<=", "dimensionless"),
        metric("escaped_weight", escaped, acceptance["escaped_weight_max"], "<=", "fraction"),
        metric("negative_control_error_ratio", ratio,
               acceptance["negative_control_error_ratio_min"], ">=", "dimensionless"),
    ]
    width = float(reference_rows[0]["bin_right_m"]) - float(reference_rows[0]["bin_left_m"])
    x = [(float(row["bin_left_m"]) + float(row["bin_right_m"])) * 0.5 for row in reference_rows]
    solution = output_dir / "CV02_solution.csv"
    write_csv(solution, ("coordinate", "numerical", "analytical", "residual", "units"),
              ({"coordinate": xx, "numerical": aa / width, "analytical": bb / width,
                "residual": (aa - bb) / width, "units": "probability_per_m"}
               for xx, aa, bb in zip(x, model, expected)))
    figures = _plot("CV02", output_dir, "constant spatial diffusion Green function",
                    x, [v / width for v in model], [v / width for v in expected],
                    "field-line coordinate s [m]", "probability density [m$^{-1}$]")
    return metrics, [solution] + figures


def _average_profile(rows: List[Dict[str, str]], predicate) -> List[float]:
    selected = [row for row in rows if predicate(row)]
    by_bin: Dict[int, List[float]] = {}
    for row in selected:
        by_bin.setdefault(int(row["bin_index"]), []).append(float(row["probability"]))
    if not by_bin:
        raise ValueError("requested CV03 profile is absent")
    return [sum(by_bin[i]) / len(by_bin[i]) for i in sorted(by_bin)]


def _l2(model: Sequence[float], reference: Sequence[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(model, reference)) /
                     max(sum(b * b for b in reference), 1.0e-30))


def _cv03(config: Dict[str, Any], rows: List[Dict[str, str]],
          reference_rows: List[Dict[str, str]], output_dir: Path
          ) -> Tuple[List[Dict[str, Any]], List[Path]]:
    numerics, acceptance = config["numerics"], config["acceptance"]
    physics = config["physics"]
    dts = list(map(float, numerics["time_steps_s"]))
    finest = min(dts)
    reference = [float(row["probability"]) for row in reference_rows]
    profiles = {dt: _average_profile(rows, lambda row, dt=dt:
        row["mode"] == "transient" and row["derivative_path"] == "analytic" and
        float(row["dt_s"]) == dt) for dt in dts}
    errors = {dt: _l2(profile, reference) for dt, profile in profiles.items()}
    numerical = _average_profile(rows, lambda row:
        row["mode"] == "transient" and row["derivative_path"] == "numerical" and
        float(row["dt_s"]) == finest)
    equilibrium = _average_profile(rows, lambda row: row["mode"] == "equilibrium")
    negative = _average_profile(rows, lambda row: row["mode"] == "negative-control")
    uniform = 1.0 / len(equilibrium)
    length = float(physics["line_length_m"])
    dx = length / len(equilibrium)
    kappa0 = float(physics["kappa0_m2_per_s"])
    amplitude = float(physics["sinusoidal_amplitude"])
    fluxes = []
    for index, probability in enumerate(equilibrium):
        right = (index + 1) % len(equilibrium)
        face = (index + 1.0) * dx
        kappa_face = kappa0 * (1.0 + amplitude *
                               math.sin(2.0 * math.pi * face / length))
        # Convert bin probability to density before differencing. The result
        # has units s^-1 and is the probability flux through this 1-D face.
        fluxes.append(-kappa_face * (equilibrium[right] - probability) /
                      (dx * dx))
    flux_rms = math.sqrt(sum(value * value for value in fluxes) / len(fluxes))
    closure = max(abs(sum(profile) - 1.0) for profile in
                  list(profiles.values()) + [numerical, equilibrium, negative])
    negative_ratio = _l2(negative, reference) / max(errors[finest], 1.0e-30)
    metrics = [
        metric("transient_relative_L2_error", errors[finest], acceptance["transient_l2_max"], "<=", "dimensionless"),
        metric("uniform_equilibrium_max_bin_error", max(abs(v - uniform) for v in equilibrium), acceptance["equilibrium_max_bin_error"], "<=", "probability"),
        metric("uniform_equilibrium_flux_rms", flux_rms, acceptance["equilibrium_flux_rms_per_s_max"], "<=", "s^-1"),
        metric("analytic_numeric_derivative_L2", _l2(numerical, profiles[finest]), acceptance["derivative_path_l2_max"], "<=", "dimensionless"),
        metric("refinement_order", _order(errors), acceptance["refinement_order_min"], ">=", "dimensionless"),
        metric("negative_control_error_ratio", negative_ratio, acceptance["negative_control_error_ratio_min"], ">=", "dimensionless"),
        metric("probability_closure_error", closure, acceptance["probability_closure_max"], "<=", "probability"),
    ]
    x = [(float(row["bin_left_m"]) + float(row["bin_right_m"])) * 0.5 for row in reference_rows]
    solution = output_dir / "CV03_solution.csv"
    write_csv(solution, ("coordinate", "numerical", "analytical", "residual", "units"),
              ({"coordinate": xx, "numerical": aa, "analytical": bb,
                "residual": aa - bb, "units": "bin_probability"}
               for xx, aa, bb in zip(x, profiles[finest], reference)))
    figures = _plot("CV03", output_dir, "sinusoidal diffusion transient",
                    x, profiles[finest], reference, "field-line coordinate s [m]",
                    "bin probability")
    return metrics, [solution] + figures


def _key(row: Dict[str, str]) -> Tuple[str, float, float]:
    return row["scenario"], float(row["atomic_mass_number"]), float(row["energy0_mev_per_nucleon"])


def _linear_slope(x: Sequence[float], y: Sequence[float]) -> float:
    """Least-squares slope used for the CV04 log-momentum spectrum."""
    x_mean, y_mean = sum(x) / len(x), sum(y) / len(y)
    return sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / \
        sum((a - x_mean) ** 2 for a in x)


def _cv04(config: Dict[str, Any], rows: List[Dict[str, str]],
          reference_rows: List[Dict[str, str]], output_dir: Path
          ) -> Tuple[List[Dict[str, Any]], List[Path]]:
    numerics, acceptance = config["numerics"], config["acceptance"]
    refs = {_key(row): row for row in reference_rows}
    final_by_scenario = {"constant": float(numerics["constant_final_time_s"]),
                         "spherical": float(numerics["spherical_final_time_s"])}
    finest_by_scenario = {
        "constant": min(map(float, numerics["constant_time_steps_s"])),
        "spherical": min(map(float, numerics["spherical_time_steps_s"])),
    }
    momentum_errors: Dict[str, Dict[float, float]] = {"constant": {}, "spherical": {}}
    energy_error = weight_error = 0.0
    final_rows = [row for row in rows if float(row["time_s"]) == final_by_scenario[row["scenario"]]]
    for row in final_rows:
        ref = refs[_key(row)]
        p_ref = float(ref["momentum_kg_m_per_s"])
        p_error = abs(float(row["momentum_kg_m_per_s"]) / p_ref - 1.0)
        dt = float(row["dt_s"])
        momentum_errors[row["scenario"]][dt] = max(momentum_errors[row["scenario"]].get(dt, 0.0), p_error)
        # Solution accuracy is gated on the finest declared realization; the
        # coarser energies remain in the evidence and determine convergence.
        if dt == finest_by_scenario[row["scenario"]]:
            e_ref = float(ref["kinetic_energy_mev_per_nucleon"])
            energy_error = max(energy_error, abs(float(row["kinetic_energy_mev_per_nucleon"]) / e_ref - 1.0))
        weight_error = max(weight_error,
                           abs(float(row["weight"]) / float(ref["weight"]) - 1.0))
    # Multiplicative adiabatic cooling preserves a momentum power-law slope.
    # Fit the finest spherical proton output rather than assuming weights alone
    # prove the spectrum was shifted consistently with momentum.
    spectral_rows = sorted([
        row for row in final_rows if row["scenario"] == "spherical" and
        float(row["dt_s"]) == finest_by_scenario["spherical"] and
        float(row["atomic_mass_number"]) == 1.0],
        key=lambda row: float(row["momentum_kg_m_per_s"]))
    fitted_slope = _linear_slope(
        [math.log(float(row["momentum_kg_m_per_s"])) for row in spectral_rows],
        [math.log(float(row["weight"])) for row in spectral_rows])
    slope_error = abs(fitted_slope + float(config["particles"]["momentum_spectral_index"]))
    constant_max = max(momentum_errors["constant"].values())
    spherical_finest = momentum_errors["spherical"][min(momentum_errors["spherical"])]
    metrics = [
        metric("constant_momentum_relative_error", constant_max, acceptance["constant_momentum_relative_error_max"], "<=", "dimensionless"),
        metric("spherical_momentum_relative_error", spherical_finest, acceptance["spherical_momentum_relative_error_max"], "<=", "dimensionless"),
        metric("spherical_refinement_order", _order(momentum_errors["spherical"]), acceptance["spherical_refinement_order_min"], ">=", "dimensionless"),
        metric("relativistic_energy_relative_error", energy_error, acceptance["energy_relative_error_max"], "<=", "dimensionless"),
        metric("statistical_weight_relative_error", weight_error, acceptance["weight_relative_error_max"], "<=", "dimensionless"),
        metric("momentum_spectral_slope_absolute_error", slope_error, acceptance["spectral_slope_absolute_error_max"], "<=", "dimensionless"),
    ]
    finest = min(map(float, numerics["spherical_time_steps_s"]))
    plot_rows = sorted([row for row in final_rows if row["scenario"] == "spherical" and
                        float(row["dt_s"]) == finest and float(row["atomic_mass_number"]) == 1.0],
                       key=lambda row: float(row["energy0_mev_per_nucleon"]))
    x = [float(row["energy0_mev_per_nucleon"]) for row in plot_rows]
    model = [float(row["kinetic_energy_mev_per_nucleon"]) for row in plot_rows]
    expected = [float(refs[_key(row)]["kinetic_energy_mev_per_nucleon"]) for row in plot_rows]
    solution = output_dir / "CV04_solution.csv"
    write_csv(solution, ("coordinate", "numerical", "analytical", "residual", "units"),
              ({"coordinate": xx, "numerical": aa, "analytical": bb,
                "residual": aa - bb, "units": "MeV_per_nucleon"}
               for xx, aa, bb in zip(x, model, expected)))
    figures = _plot("CV04", output_dir, "spherical-wind adiabatic energy shift",
                    x, model, expected, "initial energy [MeV/nucleon]",
                    "final energy [MeV/nucleon]")
    return metrics, [solution] + figures


def _cv05_key(row: Dict[str, str]) -> Tuple[float, int, float]:
    return (float(row["gradient_sign"]), int(row["particle_id"]),
            float(row["time_s"]))


def _cv05(config: Dict[str, Any], rows: List[Dict[str, str]],
          reference_rows: List[Dict[str, str]], output_dir: Path
          ) -> Tuple[List[Dict[str, Any]], List[Path]]:
    physics, numerics, acceptance = config["physics"], config["numerics"], config["acceptance"]
    final, path_scale = float(numerics["final_time_s"]), float(physics["speed_m_per_s"]) * float(numerics["final_time_s"])
    refs = {_cv05_key(row): row for row in reference_rows}
    dts = list(map(float, numerics["time_steps_s"]))
    errors: Dict[float, float] = {}
    mu_error = momentum_error = invariant_error = bound_excess = moment_error = 0.0
    second_moment_error = weight_closure = 0.0
    final_rows = [row for row in rows if float(row["time_s"]) == final]
    for dt in dts:
        subset = [row for row in rows if float(row["dt_s"]) == dt]
        errors[dt] = max(abs(float(row["position_m"]) - float(refs[_cv05_key(row)]["position_m"])) / path_scale for row in subset)
        mu_error = max(mu_error, *(abs(float(row["mu"]) - float(refs[_cv05_key(row)]["mu"])) for row in subset))
        bound_excess = max(bound_excess, *(max(0.0, abs(float(row["mu"])) - 1.0) for row in subset))
        p0 = {int(row["particle_id"]): float(row["momentum_kg_m_per_s"])
              for row in rows if float(row["dt_s"]) == dt and float(row["time_s"]) == 0.0}
        momentum_error = max(momentum_error, *(abs(float(row["momentum_kg_m_per_s"]) / p0[int(row["particle_id"])] - 1.0) for row in subset))
        for row in subset:
            ref = refs[_cv05_key(row)]
            inv0 = (1.0 - float(row["initial_mu"]) ** 2)
            if inv0 > 1.0e-12:
                invariant_error = max(invariant_error, abs(float(row["magnetic_moment_invariant"]) / inv0 - 1.0))
        for sign in (-1.0, 1.0):
            signed = [row for row in final_rows if float(row["dt_s"]) == dt and
                      float(row["gradient_sign"]) == sign]
            total = sum(float(row["weight"]) for row in signed)
            model_moment = sum(float(row["weight"]) * float(row["mu"]) for row in signed) / total
            expected_moment = sum(float(row["weight"]) * float(refs[_cv05_key(row)]["mu"]) for row in signed) / total
            moment_error = max(moment_error, abs(model_moment - expected_moment))
            model_second = sum(float(row["weight"]) * float(row["mu"]) ** 2 for row in signed) / total
            expected_second = sum(float(row["weight"]) * float(refs[_cv05_key(row)]["mu"]) ** 2 for row in signed) / total
            second_moment_error = max(second_moment_error,
                                      abs(model_second - expected_second))
            expected_weight = sum(float(row["weight"]) for row in reference_rows
                                  if float(row["gradient_sign"]) == sign and
                                  float(row["time_s"]) == final)
            weight_closure = max(weight_closure,
                                 abs(total / expected_weight - 1.0))
    finest = min(dts)
    metrics = [
        metric("position_error_over_path", errors[finest], acceptance["position_relative_to_path_max"], "<=", "dimensionless"),
        metric("pitch_angle_absolute_error", mu_error, acceptance["mu_absolute_error_max"], "<=", "dimensionless"),
        metric("momentum_relative_error", momentum_error, acceptance["momentum_relative_error_max"], "<=", "dimensionless"),
        metric("magnetic_moment_invariant_relative_error", invariant_error, acceptance["invariant_relative_error_max"], "<=", "dimensionless"),
        metric("angular_first_moment_absolute_error", moment_error, acceptance["angular_moment_absolute_error_max"], "<=", "dimensionless"),
        metric("angular_second_moment_absolute_error", second_moment_error, acceptance["angular_moment_absolute_error_max"], "<=", "dimensionless"),
        metric("statistical_weight_closure_relative_error", weight_closure, acceptance["weight_closure_relative_error_max"], "<=", "dimensionless"),
        metric("trajectory_refinement_order", _order(errors), acceptance["refinement_order_min"], ">=", "dimensionless"),
        metric("pitch_angle_bound_excess", bound_excess, acceptance["mu_bound_excess_max"], "<=", "dimensionless"),
    ]
    plot_rows = sorted([row for row in final_rows if float(row["dt_s"]) == finest and
                        float(row["gradient_sign"]) == -1.0], key=lambda row: float(row["initial_mu"]))
    x = [float(row["initial_mu"]) for row in plot_rows]
    model = [float(row["mu"]) for row in plot_rows]
    expected = [float(refs[_cv05_key(row)]["mu"]) for row in plot_rows]
    solution = output_dir / "CV05_solution.csv"
    write_csv(solution, ("coordinate", "numerical", "analytical", "residual", "units"),
              ({"coordinate": xx, "numerical": aa, "analytical": bb,
                "residual": aa - bb, "units": "pitch_angle_cosine"}
               for xx, aa, bb in zip(x, model, expected)))
    figures = _plot("CV05", output_dir, "magnetic focusing characteristic",
                    x, model, expected, "initial pitch-angle cosine", "final pitch-angle cosine")
    return metrics, [solution] + figures


def run_controlled_case(case_id: str, *, source_root: Path, input_path: Path,
                        output_dir: Path, executable: Path,
                        timeout: Optional[float]) -> Dict[str, Any]:
    """Run one CV02-CV05 linked model, reference, scorer, and figure stage."""
    started = time.monotonic()
    config = load_input(input_path, case_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved = output_dir / "resolved_input.json"
    atomic_json(resolved, config)
    n, p = config["numerics"], config.get("physics", {})
    if case_id == "CV02":
        d = config["domain"]
        arguments = ["--kappa-m2-per-s", str(p["kappa_parallel_m2_per_s"]), "--final-time-s", str(n["final_time_s"]), "--minimum-m", str(d["minimum_m"]), "--maximum-m", str(d["maximum_m"]), "--center-m", str(d["packet_center_m"]), "--profile-half-width-m", str(d["profile_half_width_m"]), "--campaign-seed", str(n["campaign_seed"]), "--seed-count", str(n["seed_count"]), "--profile-bins", str(n["profile_bins"]), "--sample-count", str(n["sample_count"]), "--time-steps-s", _list(n["time_steps_s"]), "--particle-counts", _list(n["particle_counts"])]
    elif case_id == "CV03":
        arguments = ["--line-length-m", str(p["line_length_m"]), "--kappa0-m2-per-s", str(p["kappa0_m2_per_s"]), "--amplitude", str(p["sinusoidal_amplitude"]), "--derivative-step-m", str(n["derivative_step_m"]), "--final-time-s", str(n["final_time_s"]), "--center-m", str(config["packet"]["center_m"]), "--particle-count", str(n["particle_count"]), "--seed-count", str(n["seed_count"]), "--campaign-seed", str(n["campaign_seed"]), "--profile-bins", str(n["profile_bins"]), "--time-steps-s", _list(n["time_steps_s"])]
    elif case_id == "CV04":
        c, particles = config["constants"], config["particles"]
        arguments = ["--light-speed-m-per-s", str(c["speed_of_light_m_per_s"]), "--proton-mass-kg", str(c["proton_mass_kg"]), "--constant-divergence-per-s", str(p["constant_divergence_per_s"]), "--constant-final-time-s", str(n["constant_final_time_s"]), "--spherical-final-time-s", str(n["spherical_final_time_s"]), "--radial-wind-speed-m-per-s", str(p["radial_wind_speed_m_per_s"]), "--initial-radius-m", str(p["initial_radius_m"]), "--constant-time-steps-s", _list(n["constant_time_steps_s"]), "--spherical-time-steps-s", _list(n["spherical_time_steps_s"]), "--energies-mev-per-nucleon", _list(particles["energies_mev_per_nucleon"]), "--atomic-mass-numbers", _list(particles["atomic_mass_numbers"]), "--spectral-index", str(particles["momentum_spectral_index"]), "--sample-count", str(n["sample_count"]), "--campaign-seed", str(n["campaign_seed"])]
    elif case_id == "CV05":
        c, particles = config["constants"], config["particles"]
        arguments = ["--speed-m-per-s", str(p["speed_m_per_s"]), "--proton-mass-kg", str(c["proton_mass_kg"]), "--light-speed-m-per-s", str(c["speed_of_light_m_per_s"]), "--focusing-length-m", str(p["focusing_length_m"]), "--initial-position-m", str(p["initial_position_m"]), "--final-time-s", str(n["final_time_s"]), "--time-steps-s", _list(n["time_steps_s"]), "--initial-mus", _list(particles["initial_mu"]), "--gradient-signs", _list(p["gradient_signs"]), "--sample-count", str(n["sample_count"]), "--campaign-seed", str(n["campaign_seed"])]
    else:
        raise ValueError(f"unsupported controlled case: {case_id}")
    native = run_linked_model(case_id=case_id, arguments=arguments,
                              source_root=source_root, output_dir=output_dir,
                              executable=executable, timeout=timeout)
    reference = output_dir / f"{case_id}_reference.csv"
    _reference(case_id, resolved, reference, timeout)
    model_rows, reference_rows = read_csv(native["model"]), read_csv(reference)
    scorer = {"CV02": _cv02, "CV03": _cv03, "CV04": _cv04, "CV05": _cv05}[case_id]
    metrics, derived = scorer(config, model_rows, reference_rows, output_dir)
    provenance = _provenance(case_id, resolved, executable, native["model"], reference, output_dir)
    artifacts = [resolved, native["manifest"], native["model"], native["report"],
                 native["junit"], native["log"], reference, provenance] + derived
    return finish_result(case_id=case_id, started=started, seed=int(n["campaign_seed"]),
                         input_path=input_path, executable=executable,
                         metrics=metrics, artifacts=artifacts,
                         message=f"{case_id} linked-application analytical validation passed")
