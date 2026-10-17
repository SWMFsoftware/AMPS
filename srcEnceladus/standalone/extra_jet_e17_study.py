#!/usr/bin/env python3
"""Search for one extra Porco-like jet that explains the E17 CO2 CA feature."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from co2_fitting_study import CO2_MASS_KG, FLYBYS, load_profile  # noqa: E402
from fitting import load_config, matrix_from_config  # noqa: E402
from joint_fitting_study import normalized_source_shape  # noqa: E402
from plume_model import (KB, Plume, radial_moment,  # noqa: E402
                         spherical_source)
from validation_paths import (CO2_FIT_RELATIVE, JET_CATALOG_RELATIVE,  # noqa: E402
                              add_validation_root_argument,
                              resolve_input_path, resolve_validation_root)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty CSV")
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def small_nonnegative_fit(matrix, values, sigma):
    """Exactly solve a small nonnegative weighted least-squares problem.

    The study has at most four columns. Enumerating every active set avoids
    optimizer tolerance changes during thousands of candidate evaluations.
    """
    matrix = np.asarray(matrix, dtype=float)
    values = np.asarray(values, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if (matrix.ndim != 2 or matrix.shape[0] != len(values)
            or sigma.shape != values.shape or matrix.shape[1] > 8
            or np.any(sigma <= 0) or not all(np.isfinite(item).all()
                                             for item in (matrix, values, sigma))):
        raise ValueError("Invalid small nonnegative fit")
    weighted = matrix / sigma[:, None]
    target = values / sigma
    scale = np.linalg.norm(weighted, axis=0)
    observable = scale > 0
    best_chi2 = float(target @ target)
    best = np.zeros(matrix.shape[1])
    indices = np.flatnonzero(observable)
    for count in range(1, len(indices) + 1):
        for active_tuple in itertools.combinations(indices, count):
            active = np.asarray(active_tuple, dtype=int)
            design = weighted[:, active] / scale[active]
            coefficient, _, _, _ = np.linalg.lstsq(design, target, rcond=None)
            tolerance = 1.0e-11 * max(1.0, float(np.max(np.abs(coefficient))))
            if np.any(coefficient < -tolerance):
                continue
            coefficient = np.maximum(coefficient, 0.0)
            residual = design @ coefficient - target
            chi2 = float(residual @ residual)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best.fill(0.0)
                best[active] = coefficient / scale[active]
    prediction = matrix @ best
    return {
        "coefficients": best,
        "predicted": prediction,
        "residual": (prediction - values) / sigma,
        "chi2": best_chi2,
        "chi2_per_sample": best_chi2 / len(values),
    }


def angular_distance_deg(latitude_a, longitude_a, latitude_b, longitude_b):
    lat_a, lat_b = np.deg2rad([latitude_a, latitude_b])
    lon_a, lon_b = np.deg2rad([-longitude_a, -longitude_b])
    cosine = (np.sin(lat_a) * np.sin(lat_b)
              + np.cos(lat_a) * np.cos(lat_b) * np.cos(lon_a - lon_b))
    return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))


def stripe_candidates(path: Path, maximum_step_deg: float):
    """Sample the modeled tiger-stripe polylines on the reference sphere."""
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("frame") != "IAU_ENCELADUS" or document.get("units") != "m":
        raise ValueError("Stripe geometry must be IAU_ENCELADUS in metres")
    candidates = []
    for stripe, points in document["stripes"].items():
        unit_points = np.asarray(points, dtype=float)
        unit_points /= np.linalg.norm(unit_points, axis=1)[:, None]
        for segment, (start, stop) in enumerate(
                zip(unit_points[:-1], unit_points[1:])):
            angle_deg = float(np.rad2deg(np.arccos(np.clip(start @ stop, -1, 1))))
            count = max(2, int(math.ceil(angle_deg / maximum_step_deg)) + 1)
            for fraction in np.linspace(0.0, 1.0, count, endpoint=False):
                point = (1.0 - fraction) * start + fraction * stop
                point /= np.linalg.norm(point)
                candidates.append({
                    "stripe": stripe,
                    "segment": segment,
                    "segment_fraction": float(fraction),
                    "latitude_deg": float(np.rad2deg(np.arcsin(point[2]))),
                    "west_longitude_deg": float(
                        (-np.rad2deg(np.arctan2(point[1], point[0]))) % 360.0),
                })
        point = unit_points[-1]
        candidates.append({
            "stripe": stripe,
            "segment": len(unit_points) - 2,
            "segment_fraction": 1.0,
            "latitude_deg": float(np.rad2deg(np.arcsin(point[2]))),
            "west_longitude_deg": float(
                (-np.rad2deg(np.arctan2(point[1], point[0]))) % 360.0),
        })
    return candidates


class CandidateEvaluator:
    def __init__(self, points, values, sigma, baseline_column,
                 background_columns, reference_radius_m, temperature_K,
                 bulk_speed_m_s, source_law, selection):
        self.points = np.asarray(points)
        self.values = np.asarray(values)
        self.sigma = np.asarray(sigma)
        self.baseline_column = np.asarray(baseline_column)
        self.background_columns = np.asarray(background_columns)
        self.reference_radius_m = float(reference_radius_m)
        self.temperature_K = float(temperature_K)
        self.bulk_speed_m_s = float(bulk_speed_m_s)
        self.source_law = source_law
        self.selection = np.asarray(selection, dtype=bool)
        self.speed_ratio = bulk_speed_m_s / math.sqrt(
            2.0 * KB * temperature_K / CO2_MASS_KG)
        if self.source_law != "tenishev2010_point":
            raise ValueError(
                "Fast source-location search currently requires the default "
                "tenishev2010_point law")

    def raw_candidate_column(self, parameters):
        latitude, longitude, tilt, azimuth = parameters
        position, normal, axis = spherical_source(
            self.reference_radius_m, latitude, longitude, tilt, azimuth)
        displacement = self.points - position
        radius = np.linalg.norm(displacement, axis=1)
        ray = displacement / radius[:, None]
        outward = ray @ normal >= 0
        return np.where(
            outward,
            radial_moment(2, ray @ axis, self.speed_ratio) / radius**2,
            0.0)

    def objective(self, parameters):
        candidate = self.raw_candidate_column(parameters)
        matrix = np.column_stack(
            (self.baseline_column, candidate, self.background_columns))
        fit = small_nonnegative_fit(
            matrix[self.selection], self.values[self.selection],
            self.sigma[self.selection])
        return fit["chi2"]


def exact_candidate_column(points, config, parameters):
    latitude, longitude, tilt, azimuth = parameters
    options = config["source_case_options"]["porco2014_98"]
    position, normal, axis = spherical_source(
        config["reference_radius_m"], latitude, longitude, tilt, azimuth)
    plume = Plume(
        tuple(position), tuple(normal), tuple(axis), options["temperature_K"],
        options["bulk_speed_m_s"], CO2_MASS_KG, config["source_law"])
    return plume.density_per_rate(points) * 1.0e28


def profile_fit(profile, matrix, labels, shape, candidate_column,
                fixed_extra_ratio=None):
    source_column = matrix @ shape
    background_indices = [
        index for index, label in enumerate(labels)
        if label.startswith("background:") and np.any(matrix[:, index] != 0)]
    backgrounds = matrix[:, background_indices]
    baseline = small_nonnegative_fit(
        np.column_stack((source_column, backgrounds)), profile["density"],
        profile["sigma"])
    augmented = small_nonnegative_fit(
        np.column_stack((source_column, candidate_column, backgrounds)),
        profile["density"], profile["sigma"])
    result = {"baseline": baseline, "augmented": augmented}
    if fixed_extra_ratio is not None:
        mixture = source_column + fixed_extra_ratio * candidate_column
        result["e17_mixture_transfer"] = small_nonnegative_fit(
            np.column_stack((mixture, backgrounds)), profile["density"],
            profile["sigma"])
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    add_validation_root_argument(parser)
    parser.add_argument(
        "--jet-catalog", type=Path,
        help="Porco catalog JSON; overrides the conventional validation-root location")
    parser.add_argument(
        "--co2-data-root", type=Path,
        help="co2_fit_v1 directory; overrides the conventional validation-root location")
    parser.add_argument("--co2-baseline-result", type=Path, required=True)
    parser.add_argument("--stripe-geometry", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--search-half-window-s", type=float, default=2.5)
    parser.add_argument("--grid-step-deg", type=float, default=2.0)
    parser.add_argument("--stripe-step-deg", type=float, default=0.025)
    parser.add_argument("--seed", type=int, default=174)
    parser.add_argument("--free-direction-seeds", type=int, default=3)
    args = parser.parse_args()
    try:
        validation_root = resolve_validation_root(args.validation_data_root)
        args.jet_catalog = resolve_input_path(
            args.jet_catalog, validation_root, JET_CATALOG_RELATIVE,
            "--jet-catalog")
        args.co2_data_root = resolve_input_path(
            args.co2_data_root, validation_root, CO2_FIT_RELATIVE,
            "--co2-data-root")
    except ValueError as exc:
        parser.error(str(exc))
    if args.output_dir.exists():
        parser.error("output directory already exists")
    if (args.search_half_window_s <= 0
            or not 0.25 <= args.grid_step_deg <= 10
            or not 0.005 <= args.stripe_step_deg <= 0.25):
        parser.error("invalid window, surface-grid step, or stripe step")
    args.output_dir.mkdir(parents=True)

    catalog = json.loads(args.jet_catalog.read_text(encoding="utf-8"))
    baseline_result = json.loads(
        args.co2_baseline_result.read_text(encoding="utf-8"))
    if baseline_result.get("source_case") != "porco2014_98":
        raise ValueError("Baseline result must be the porco2014_98 CO2 study")
    config = load_config(args.config, "porco2014_98", args.jet_catalog)
    config = copy.deepcopy(config)
    config["molecular_species"] = "CO2"
    config["molecular_mass_kg"] = CO2_MASS_KG
    options = config["source_case_options"]["porco2014_98"]
    stripe_path = (args.stripe_geometry or
                   (args.config.resolve().parent / config["stripe_geometry_file"]))
    stripe_path = stripe_path.resolve()
    stripe_grid = stripe_candidates(stripe_path, args.stripe_step_deg)

    profiles, matrices = {}, {}
    labels = None
    for flyby in FLYBYS:
        profiles[flyby] = load_profile(
            args.co2_data_root / f"{flyby.lower()}_co2_density.csv")
        matrix, current_labels = matrix_from_config(
            config, profiles[flyby]["points"])
        if labels is None:
            labels = current_labels
        elif labels != current_labels:
            raise AssertionError("Component labels changed between flybys")
        matrices[flyby] = matrix
    assert labels is not None
    combined = baseline_result["combined_shape_fits"]["shared_shape_all_co2"]
    shape = normalized_source_shape(np.asarray(combined["coefficients"]), labels)

    e17 = profiles["E17"]
    e17_matrix = matrices["E17"]
    source_column = e17_matrix @ shape
    background_indices = [
        index for index, label in enumerate(labels)
        if label.startswith("background:") and np.any(e17_matrix[:, index] != 0)]
    backgrounds = e17_matrix[:, background_indices]
    latitude_bounds = (-90.0, max(jet["latitude_deg"] for jet in catalog["jets"]))
    tilt_bounds = (0.0, max(jet["tilt_deg"] for jet in catalog["jets"]))

    radial_sensitivity = []
    primary_optimum = None
    for half_window in sorted(set((1.5, args.search_half_window_s, 4.0))):
        selection = np.abs(e17["time"]) <= half_window
        evaluator = CandidateEvaluator(
            e17["points"], e17["density"], e17["sigma"], source_column,
            backgrounds, config["reference_radius_m"], options["temperature_K"],
            options["bulk_speed_m_s"], config["source_law"], selection)
        baseline_window = small_nonnegative_fit(
            np.column_stack((source_column, backgrounds))[selection],
            e17["density"][selection], e17["sigma"][selection])
        optimum = differential_evolution(
            lambda value: evaluator.objective((value[0], value[1], 0.0, 0.0)),
            (latitude_bounds, (0.0, 360.0)), seed=args.seed,
            maxiter=60, popsize=12, tol=1.0e-8, polish=True,
            updating="immediate", workers=1)
        item = {
            "search_half_window_s": half_window,
            "samples": int(selection.sum()),
            "baseline_window_chi2": baseline_window["chi2"],
            "augmented_window_chi2": float(optimum.fun),
            "delta_chi2": baseline_window["chi2"] - float(optimum.fun),
            "latitude_deg": float(optimum.x[0]),
            "west_longitude_deg": float(optimum.x[1] % 360.0),
            "tilt_deg": 0.0,
            "azimuth_deg": 0.0,
            "optimizer_success": bool(optimum.success),
            "optimizer_message": optimum.message,
        }
        radial_sensitivity.append(item)
        if math.isclose(half_window, args.search_half_window_s):
            primary_optimum = item
    assert primary_optimum is not None
    primary_parameters = (
        primary_optimum["latitude_deg"],
        primary_optimum["west_longitude_deg"], 0.0, 0.0)

    stripe_sensitivity = []
    recommended_candidate = None
    for half_window in sorted(set((1.5, args.search_half_window_s, 4.0))):
        selection = np.abs(e17["time"]) <= half_window
        evaluator = CandidateEvaluator(
            e17["points"], e17["density"], e17["sigma"], source_column,
            backgrounds, config["reference_radius_m"], options["temperature_K"],
            options["bulk_speed_m_s"], config["source_law"], selection)
        best = min(
            stripe_grid,
            key=lambda item: evaluator.objective((
                item["latitude_deg"], item["west_longitude_deg"], 0.0, 0.0)))
        best_chi2 = evaluator.objective((
            best["latitude_deg"], best["west_longitude_deg"], 0.0, 0.0))
        baseline_fit = small_nonnegative_fit(
            np.column_stack((source_column, backgrounds))[selection],
            e17["density"][selection], e17["sigma"][selection])
        item = {
            "search_half_window_s": half_window,
            "samples": int(selection.sum()),
            "stripe": best["stripe"],
            "segment": best["segment"],
            "segment_fraction": best["segment_fraction"],
            "baseline_window_chi2": baseline_fit["chi2"],
            "augmented_window_chi2": best_chi2,
            "delta_chi2": baseline_fit["chi2"] - best_chi2,
            "latitude_deg": best["latitude_deg"],
            "west_longitude_deg": best["west_longitude_deg"],
            "tilt_deg": 0.0,
            "azimuth_deg": 0.0,
        }
        stripe_sensitivity.append(item)
        if math.isclose(half_window, args.search_half_window_s):
            recommended_candidate = item
    assert recommended_candidate is not None
    recommended_parameters = (
        recommended_candidate["latitude_deg"],
        recommended_candidate["west_longitude_deg"], 0.0, 0.0)

    primary_selection = np.abs(e17["time"]) <= args.search_half_window_s
    primary_evaluator = CandidateEvaluator(
        e17["points"], e17["density"], e17["sigma"], source_column,
        backgrounds, config["reference_radius_m"], options["temperature_K"],
        options["bulk_speed_m_s"], config["source_law"], primary_selection)
    free_direction = []
    free_bounds = (latitude_bounds, (0.0, 360.0), tilt_bounds, (0.0, 360.0))
    for offset in range(args.free_direction_seeds):
        optimum = differential_evolution(
            primary_evaluator.objective, free_bounds,
            seed=args.seed + 1009 * offset, maxiter=40, popsize=10,
            tol=1.0e-7, polish=True, updating="immediate", workers=1)
        free_direction.append({
            "seed": args.seed + 1009 * offset,
            "window_chi2": float(optimum.fun),
            "latitude_deg": float(optimum.x[0]),
            "west_longitude_deg": float(optimum.x[1] % 360.0),
            "tilt_deg": float(optimum.x[2]),
            "azimuth_deg": float(optimum.x[3] % 360.0),
            "optimizer_success": bool(optimum.success),
            "optimizer_message": optimum.message,
        })

    latitudes = np.arange(latitude_bounds[0], latitude_bounds[1] + 1.0e-9,
                          args.grid_step_deg)
    longitudes = np.arange(0.0, 360.0, args.grid_step_deg)
    baseline_window = small_nonnegative_fit(
        np.column_stack((source_column, backgrounds))[primary_selection],
        e17["density"][primary_selection], e17["sigma"][primary_selection])
    scan_rows = []
    scan_improvement = np.empty((len(latitudes), len(longitudes)))
    for lat_index, latitude in enumerate(latitudes):
        for lon_index, longitude in enumerate(longitudes):
            chi2 = primary_evaluator.objective((latitude, longitude, 0.0, 0.0))
            improvement = baseline_window["chi2"] - chi2
            scan_improvement[lat_index, lon_index] = improvement
            scan_rows.append({
                "latitude_deg": latitude,
                "west_longitude_deg": longitude,
                "window_chi2": chi2,
                "delta_chi2": improvement,
            })

    nearest = sorted(({
        "porco_id": jet["porco_id"], "id": jet["id"], "sulcus": jet["sulcus"],
        "latitude_deg": jet["latitude_deg"],
        "west_longitude_deg": jet["west_longitude_deg"],
        "angular_distance_deg": angular_distance_deg(
            recommended_parameters[0], recommended_parameters[1], jet["latitude_deg"],
            jet["west_longitude_deg"]),
    } for jet in catalog["jets"]), key=lambda item: item["angular_distance_deg"])

    candidate_columns = {
        flyby: exact_candidate_column(profiles[flyby]["points"], config,
                                      recommended_parameters)
        for flyby in FLYBYS}
    e17_fit = profile_fit(e17, e17_matrix, labels, shape,
                          candidate_columns["E17"])
    e17_coeff = e17_fit["augmented"]["coefficients"]
    if e17_coeff[0] <= 0:
        raise RuntimeError("E17 baseline rate is zero in augmented fit")
    fixed_extra_ratio = float(e17_coeff[1] / e17_coeff[0])
    profile_results = {
        flyby: profile_fit(profiles[flyby], matrices[flyby], labels, shape,
                           candidate_columns[flyby], fixed_extra_ratio)
        for flyby in FLYBYS}

    comparison_rows = []
    for flyby in FLYBYS:
        result = profile_results[flyby]
        baseline = result["baseline"]
        augmented = result["augmented"]
        transfer = result["e17_mixture_transfer"]
        comparison_rows.append({
            "flyby": flyby,
            "samples": len(profiles[flyby]["time"]),
            "baseline_chi2": baseline["chi2"],
            "extra_jet_chi2": augmented["chi2"],
            "extra_jet_delta_chi2": baseline["chi2"] - augmented["chi2"],
            "extra_jet_rate_s-1": augmented["coefficients"][1] * 1.0e28,
            "baseline_rate_s-1": augmented["coefficients"][0] * 1.0e28,
            "e17_mixture_transfer_chi2": transfer["chi2"],
            "e17_mixture_transfer_delta_chi2": baseline["chi2"] - transfer["chi2"],
        })

    write_csv(args.output_dir / "radial_window_sensitivity.csv", radial_sensitivity)
    write_csv(args.output_dir / "stripe_window_sensitivity.csv", stripe_sensitivity)
    write_csv(args.output_dir / "free_direction_solutions.csv", free_direction)
    write_csv(args.output_dir / "radial_location_scan.csv", scan_rows)
    write_csv(args.output_dir / "cross_flyby_comparison.csv", comparison_rows)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    time = e17["time"]
    axes[0].errorbar(time, e17["density"] / 1.0e12,
                     yerr=e17["sigma"] / 1.0e12, fmt=".", ms=3,
                     color="black", alpha=0.65, label="E17 mass-44 data")
    axes[0].plot(time, profile_results["E17"]["baseline"]["predicted"] / 1.0e12,
                 label="98-jet shared-shape baseline", lw=1.8)
    axes[0].plot(time, profile_results["E17"]["augmented"]["predicted"] / 1.0e12,
                 label="baseline + fitted extra stripe-constrained jet", lw=1.8)
    axes[0].axvspan(-args.search_half_window_s, args.search_half_window_s,
                    color="0.85", alpha=0.5, label="location-search window")
    axes[0].set_ylabel(r"CO$_2$ density ($10^{12}$ m$^{-3}$)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)
    for name, color in (("baseline", "C0"), ("augmented", "C1")):
        residual = ((profile_results["E17"][name]["predicted"] - e17["density"])
                    / e17["sigma"])
        axes[1].plot(time, residual, color=color, lw=1.0, label=name)
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].set_xlabel("Time from E17 closest approach (s)")
    axes[1].set_ylabel("weighted residual")
    axes[1].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.output_dir / "e17_extra_jet_fit.png", dpi=220)
    fig.savefig(args.output_dir / "e17_extra_jet_fit.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    image = ax.pcolormesh(longitudes, latitudes, scan_improvement,
                          shading="nearest", cmap="viridis")
    ax.scatter([jet["west_longitude_deg"] for jet in catalog["jets"]],
               [jet["latitude_deg"] for jet in catalog["jets"]], s=8,
               facecolors="none", edgecolors="white", linewidths=0.45,
               label="Porco catalog jets")
    ax.scatter([primary_parameters[1]], [primary_parameters[0]], marker="*",
               s=100, color="orange", edgecolors="black", linewidths=0.5,
               label="unconstrained radial optimum")
    ax.scatter([recommended_parameters[1]], [recommended_parameters[0]],
               marker="*", s=140, color="red", edgecolors="black",
               linewidths=0.5, label="stripe-constrained candidate")
    ax.set_xlim(0, 360)
    ax.set_ylim(latitude_bounds)
    ax.set_xlabel("West longitude (deg)")
    ax.set_ylabel("Planetocentric latitude (deg)")
    ax.set_title(f"E17 |t| <= {args.search_half_window_s:g} s radial-jet scan")
    fig.colorbar(image, ax=ax, label=r"window $\Delta\chi^2$")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output_dir / "radial_location_scan.png", dpi=220)
    fig.savefig(args.output_dir / "radial_location_scan.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(FLYBYS))
    ax.bar(x - 0.18, [row["extra_jet_delta_chi2"] for row in comparison_rows],
           0.36, label="independent extra-jet rate")
    ax.bar(x + 0.18,
           [row["e17_mixture_transfer_delta_chi2"] for row in comparison_rows],
           0.36, label="fixed E17 extra/baseline ratio")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xticks(x, FLYBYS)
    ax.set_ylabel(r"Improvement $\Delta\chi^2$")
    ax.set_title("Cross-flyby test of the E17-selected radial candidate")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.output_dir / "cross_flyby_improvement.png", dpi=220)
    fig.savefig(args.output_dir / "cross_flyby_improvement.pdf")
    plt.close(fig)

    position = {
        "latitude_deg": recommended_parameters[0],
        "west_longitude_deg": recommended_parameters[1],
        "tilt_deg": 0.0,
        "azimuth_deg": 0.0,
        "temperature_K": options["temperature_K"],
        "bulk_speed_m_s": options["bulk_speed_m_s"],
        "source_law": config["source_law"],
    }
    result = {
        "study": "E17 closest-approach extra-jet search and CO2 cross-flyby test",
        "source_case": "porco2014_98_plus_one_candidate",
        "search_target": {
            "flyby": "E17", "species": "CO2",
            "half_window_s": args.search_half_window_s,
            "selection": "absolute time from closest approach <= half_window_s",
        },
        "input_sha256": {
            str(args.config): sha256(args.config),
            str(args.jet_catalog): sha256(args.jet_catalog),
            str(stripe_path): sha256(stripe_path),
            str(args.co2_baseline_result): sha256(args.co2_baseline_result),
            **{str(args.co2_data_root / f"{flyby.lower()}_co2_density.csv"):
               sha256(args.co2_data_root / f"{flyby.lower()}_co2_density.csv")
               for flyby in FLYBYS},
        },
        "search_bounds": {
            "latitude_deg": latitude_bounds,
            "west_longitude_deg": (0.0, 360.0),
            "tilt_deg_free_direction_sensitivity": tilt_bounds,
            "basis": "surface-latitude and tilt envelope of the 98-jet catalog",
        },
        "primary_candidate": position,
        "unconstrained_surface_candidate": {
            "latitude_deg": primary_parameters[0],
            "west_longitude_deg": primary_parameters[1],
            "tilt_deg": 0.0,
            "azimuth_deg": 0.0,
        },
        "nearest_catalog_jets": nearest[:10],
        "radial_window_sensitivity": radial_sensitivity,
        "stripe_window_sensitivity": stripe_sensitivity,
        "free_direction_sensitivity": free_direction,
        "e17_extra_to_baseline_rate_ratio": fixed_extra_ratio,
        "cross_flyby": {
            row["flyby"]: {
                **row,
                "baseline_prediction": profile_results[row["flyby"]]["baseline"]["predicted"],
                "augmented_prediction": profile_results[row["flyby"]]["augmented"]["predicted"],
            } for row in comparison_rows
        },
        "interpretation_limits": [
            "The primary result is a radial candidate constrained to the modeled tiger-stripe polylines, not a claimed geological detection.",
            "The free-direction solutions demonstrate that one flyby track does not uniquely determine four jet-geometry parameters.",
            "The candidate was selected on E17, so its E17 chi-square improvement is not an independent validation and is subject to a look-elsewhere effect.",
            "Cross-flyby fits permit an observation-specific extra-jet rate; the fixed-mixture test separately transfers the E17 extra/baseline ratio.",
            "Counting sigmas omit common calibration and model-discrepancy errors, so delta chi-square is diagnostic rather than a formal discovery significance.",
        ],
    }
    (args.output_dir / "study_results.json").write_text(
        json.dumps(result, indent=2, default=json_default) + "\n", encoding="utf-8")

    report = [
        "# E17 extra-jet study", "",
        "The study holds the all-CO2 shared relative shape of the published 98 "
        "Porco jets fixed, adds one candidate jet, and refits only the total "
        "baseline rate, candidate rate, and empirical backgrounds.", "",
        "## Primary constrained candidate", "",
        f"The declared primary search uses a radial jet constrained to the "
        f"modeled {recommended_candidate['stripe']} polyline and the E17 |t| <= "
        f"{args.search_half_window_s:g} s closest-approach window. Its location "
        f"is latitude {position['latitude_deg']:.4f} deg and west longitude "
        f"{position['west_longitude_deg']:.4f} deg. The closest catalog source "
        f"is {nearest[0]['id']} ({nearest[0]['sulcus']}) at an angular distance "
        f"of {nearest[0]['angular_distance_deg']:.3f} deg.", "",
        "## Cross-flyby CO2 comparison", "",
        "| Flyby | Baseline chi2 | Baseline + extra chi2 | Delta chi2 | Extra rate (s^-1) | Fixed E17-mixture delta chi2 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in comparison_rows:
        report.append(
            f"| {row['flyby']} | {row['baseline_chi2']:.3f} | "
            f"{row['extra_jet_chi2']:.3f} | {row['extra_jet_delta_chi2']:.3f} | "
            f"{row['extra_jet_rate_s-1']:.6e} | "
            f"{row['e17_mixture_transfer_delta_chi2']:.3f} |")
    report += [
        "", "## Interpretation", "",
        "The radial-window and free-direction tables must be inspected before "
        "assigning geological meaning. Movement of the optimum with the time "
        "window, or separated free-direction solutions with similar chi2, means "
        "the trajectory constrains a plume signature more strongly than a unique "
        "surface vent. The E17 improvement is a development result; only the "
        "other flybys are cross-observation tests.", "",
        "See `e17_extra_jet_fit.png`, `radial_location_scan.png`, "
        "`cross_flyby_improvement.png`, and the CSV/JSON products.", "",
    ]
    (args.output_dir / "report.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir),
                      "candidate": position,
                      "nearest_catalog_jet": nearest[0]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
