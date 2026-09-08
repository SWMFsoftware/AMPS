#!/usr/bin/env python3
"""Derive storm-time cutoff morphology, lag, and hysteresis products.

Input rows are the MLT-resolved ACCESS_T50 boundaries written by
``run_morphology.py``.  All calculations are deterministic for the configured
random seed.  The script writes long-form intermediate tables so every figure
or manuscript number can be traced to specific model epochs and boundary cells.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from study_common import (
    DriverRow, default_output_root, finite_float, fit_first_harmonic,
    format_utc, interpolate_driver, load_config, parse_utc, pearson, quantile,
    read_csv, read_driver, resolve_output_path, write_csv,
)


def group_rows(rows: Sequence[Mapping[str, str]], keys: Sequence[str]):
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(row)
    return grouped


def bootstrap_correlation(x: Sequence[float], y: Sequence[float], block_length: int,
                          replicates: int, rng: random.Random) -> Tuple[Optional[float], Optional[float]]:
    """Moving-block bootstrap interval for a serially correlated correlation."""
    n = len(x)
    if n < 4 or replicates < 1:
        return None, None
    block_length = max(1, min(block_length, n))
    starts = list(range(0, n - block_length + 1))
    samples: List[float] = []
    for _ in range(replicates):
        indices: List[int] = []
        while len(indices) < n:
            start = rng.choice(starts)
            indices.extend(range(start, start + block_length))
        indices = indices[:n]
        value = pearson([x[i] for i in indices], [y[i] for i in indices])
        if value is not None and math.isfinite(value):
            samples.append(value)
    if not samples:
        return None, None
    return quantile(samples, 0.025), quantile(samples, 0.975)


def morphology_products(rows: Sequence[Mapping[str, str]], config: Mapping[str, object]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Fit the first MLT harmonic and compute analyzed-shell access fraction."""
    keys = ("epoch_utc", "altitude_km", "rigidity_gv", "hemisphere")
    harmonic_rows: List[Dict[str, object]] = []
    long_rows: List[Dict[str, object]] = []
    lat_min, lat_max = [float(value) for value in config["model"]["latitude_band_abs_deg"]]  # type: ignore[index]
    denominator = math.sin(math.radians(lat_max)) - math.sin(math.radians(lat_min))
    for key, subset in sorted(group_rows(rows, keys).items()):
        valid = []
        for row in subset:
            boundary = finite_float(row.get("boundary_aacgm_lat_deg"))
            mlt = finite_float(row.get("mlt_hour"))
            if boundary is not None and mlt is not None:
                valid.append((mlt, abs(boundary)))
        base = {
            "epoch_utc": key[0], "altitude_km": float(key[1]),
            "rigidity_gv": float(key[2]), "hemisphere": key[3],
            "n_valid_mlt": len(valid), "n_requested_mlt": len(subset),
        }
        if len(valid) >= 3:
            fit = fit_first_harmonic(
                [item[0] for item in valid], [item[1] for item in valid]
            )
            # The modeled shell stops at lat_max.  Clamp the T50 boundary to the
            # analyzed band before integrating the area poleward of the boundary.
            access_terms = []
            for _, latitude in valid:
                bounded = min(lat_max, max(lat_min, latitude))
                access_terms.append(
                    (math.sin(math.radians(lat_max)) - math.sin(math.radians(bounded)))
                    / denominator
                )
            fit["accessible_area_fraction_in_analyzed_band"] = statistics.fmean(access_terms)
            base.update(fit)
        else:
            base.update({
                "mean_latitude_deg": None, "amplitude_deg": None,
                "phase_mlt_hour": None, "fit_rms_deg": None,
                "accessible_area_fraction_in_analyzed_band": None,
            })
        harmonic_rows.append(base)

    # Boundary speed and quiet-relative erosion are calculated only after all
    # epoch fits exist, preserving the same grouping and sign convention.
    by_series = group_rows(
        [{key: str(value) if value is not None else "" for key, value in row.items()}
         for row in harmonic_rows],
        ("altitude_km", "rigidity_gv", "hemisphere"),
    )
    quiet_limit = parse_utc(config["event"]["compression_search_start_utc"])  # type: ignore[index]
    for _, series in by_series.items():
        series.sort(key=lambda row: parse_utc(row["epoch_utc"]))
        quiet_values = [finite_float(row.get("mean_latitude_deg")) for row in series
                        if parse_utc(row["epoch_utc"]) < quiet_limit]
        quiet_values = [value for value in quiet_values if value is not None]
        quiet_reference = statistics.median(quiet_values) if quiet_values else None
        for index, row in enumerate(series):
            mean_lat = finite_float(row.get("mean_latitude_deg"))
            speed = None
            if 0 < index < len(series) - 1:
                previous = series[index - 1]
                following = series[index + 1]
                previous_lat = finite_float(previous.get("mean_latitude_deg"))
                following_lat = finite_float(following.get("mean_latitude_deg"))
                hours = (parse_utc(following["epoch_utc"]) -
                         parse_utc(previous["epoch_utc"])).total_seconds() / 3600.0
                if previous_lat is not None and following_lat is not None and hours > 0:
                    speed = (following_lat - previous_lat) / hours
            row["quiet_reference_mean_latitude_deg"] = quiet_reference
            row["cutoff_erosion_deg"] = (
                None if mean_lat is None or quiet_reference is None
                else mean_lat - quiet_reference
            )
            row["boundary_speed_deg_per_hour"] = speed
            long_rows.append(dict(row))
    return harmonic_rows, long_rows


def lag_products(harmonics: Sequence[Mapping[str, object]], driver: Sequence[DriverRow],
                 config: Mapping[str, object]) -> List[Dict[str, object]]:
    """Compute driver/mean-cutoff correlations on the configured lag grid."""
    analysis = config["analysis"]  # type: ignore[index]
    min_minutes = int(round(float(analysis["lag_min_hours"]) * 60.0))
    max_minutes = int(round(float(analysis["lag_max_hours"]) * 60.0))
    step = int(analysis["lag_step_minutes"])
    replicates = int(analysis["bootstrap_replicates"])
    block_hours = float(analysis["bootstrap_block_hours"])
    rng = random.Random(int(analysis["random_seed"]))
    variables = ("pdyn_npa", "bz_nt", "symh_nt", "w1", "w2", "w3", "w4", "w5", "w6")
    output: List[Dict[str, object]] = []
    grouped = defaultdict(list)
    for row in harmonics:
        mean = finite_float(row.get("mean_latitude_deg"))
        if mean is not None:
            grouped[(row["altitude_km"], row["rigidity_gv"], row["hemisphere"])].append(row)
    for key, series in sorted(grouped.items()):
        series.sort(key=lambda row: parse_utc(str(row["epoch_utc"])))
        if len(series) > 1:
            spacings = [
                (parse_utc(str(right["epoch_utc"])) - parse_utc(str(left["epoch_utc"]))).total_seconds() / 3600.0
                for left, right in zip(series, series[1:])
            ]
            nominal_hours = statistics.median(value for value in spacings if value > 0)
        else:
            nominal_hours = block_hours
        block_length = max(1, int(round(block_hours / nominal_hours)))
        for variable in variables:
            for lag_minutes in range(min_minutes, max_minutes + 1, step):
                x: List[float] = []
                y: List[float] = []
                for row in series:
                    epoch = parse_utc(str(row["epoch_utc"]))
                    driver_epoch = epoch - timedelta(minutes=lag_minutes)
                    try:
                        sampled = interpolate_driver(driver, driver_epoch)
                    except ValueError:
                        continue
                    x.append(float(getattr(sampled, variable)))
                    y.append(float(row["mean_latitude_deg"]))
                correlation = pearson(x, y)
                ci_low, ci_high = bootstrap_correlation(
                    x, y, block_length, replicates, rng
                )
                output.append({
                    "altitude_km": key[0], "rigidity_gv": key[1],
                    "hemisphere": key[2], "driver_variable": variable,
                    "lag_minutes": lag_minutes,
                    "positive_lag_means_cutoff_follows_driver": True,
                    "n_paired_epochs": len(x), "correlation": correlation,
                    "bootstrap_ci_low": ci_low, "bootstrap_ci_high": ci_high,
                    "bootstrap_block_hours": block_hours,
                })
    return output


def hysteresis_products(boundary_rows: Sequence[Mapping[str, str]], driver: Sequence[DriverRow],
                        compression: datetime, main_phase: datetime,
                        config: Mapping[str, object]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Match main/recovery cells at similar instantaneous forcing conditions."""
    analysis = config["analysis"]  # type: ignore[index]
    symh_tol = float(analysis["hysteresis_symh_tolerance_nt"])
    pdyn_tol = float(analysis["hysteresis_pdyn_fractional_tolerance"])
    bz_tol = float(analysis["hysteresis_bz_tolerance_nt"])
    replicates = int(analysis["bootstrap_replicates"])
    rng = random.Random(int(analysis["random_seed"]) + 1)
    keys = ("altitude_km", "rigidity_gv", "hemisphere", "mlt_hour")
    pairs: List[Dict[str, object]] = []
    for key, rows in sorted(group_rows(boundary_rows, keys).items()):
        valid = []
        for row in rows:
            boundary = finite_float(row.get("boundary_aacgm_lat_deg"))
            if boundary is None:
                continue
            epoch = parse_utc(row["epoch_utc"])
            valid.append((epoch, abs(boundary), interpolate_driver(driver, epoch)))
        # Restrict the main branch to storm development.  Pre-event quiet cells
        # can share SYM-H with late recovery but are not part of the hysteresis
        # loop and would bias the matched contrast toward zero.
        main = [item for item in valid if compression <= item[0] < main_phase]
        recovery = [item for item in valid if item[0] > main_phase]
        available = set(range(len(recovery)))
        for main_epoch, main_lat, main_driver in sorted(main, reverse=True):
            candidates = []
            for index in available:
                rec_epoch, rec_lat, rec_driver = recovery[index]
                symh_difference = abs(rec_driver.symh_nt - main_driver.symh_nt)
                pdyn_fraction = abs(rec_driver.pdyn_npa - main_driver.pdyn_npa) / max(
                    1.0e-9, abs(main_driver.pdyn_npa)
                )
                bz_difference = abs(rec_driver.bz_nt - main_driver.bz_nt)
                if symh_difference <= symh_tol:
                    score = symh_difference / symh_tol + pdyn_fraction / pdyn_tol + bz_difference / bz_tol
                    candidates.append((score, index, pdyn_fraction, bz_difference))
            if not candidates:
                continue
            _, index, pdyn_fraction, bz_difference = min(candidates)
            available.remove(index)
            rec_epoch, rec_lat, rec_driver = recovery[index]
            strict = pdyn_fraction <= pdyn_tol and bz_difference <= bz_tol
            pairs.append({
                "altitude_km": float(key[0]), "rigidity_gv": float(key[1]),
                "hemisphere": key[2], "mlt_hour": float(key[3]),
                "main_epoch_utc": format_utc(main_epoch),
                "recovery_epoch_utc": format_utc(rec_epoch),
                "main_boundary_deg": main_lat, "recovery_boundary_deg": rec_lat,
                "recovery_minus_main_deg": rec_lat - main_lat,
                "main_symh_nt": main_driver.symh_nt,
                "recovery_symh_nt": rec_driver.symh_nt,
                "delta_symh_nt": rec_driver.symh_nt - main_driver.symh_nt,
                "pdyn_fractional_difference": pdyn_fraction,
                "bz_absolute_difference_nt": bz_difference,
                "strict_instantaneous_driver_match": strict,
            })

    summaries: List[Dict[str, object]] = []
    summary_groups = group_rows(
        [{key: str(value) for key, value in row.items()} for row in pairs],
        ("altitude_km", "rigidity_gv", "hemisphere"),
    )
    for key, rows in sorted(summary_groups.items()):
        for label, selected in (
            ("SYMH_ONLY", rows),
            ("STRICT", [row for row in rows if row["strict_instantaneous_driver_match"] == "True"]),
        ):
            values = [float(row["recovery_minus_main_deg"]) for row in selected]
            boot = []
            if values:
                for _ in range(replicates):
                    boot.append(statistics.median(rng.choices(values, k=len(values))))
            summaries.append({
                "altitude_km": float(key[0]), "rigidity_gv": float(key[1]),
                "hemisphere": key[2], "match_definition": label,
                "n_pairs": len(values),
                "median_recovery_minus_main_deg": statistics.median(values) if values else None,
                "bootstrap_ci_low": quantile(boot, 0.025) if boot else None,
                "bootstrap_ci_high": quantile(boot, 0.975) if boot else None,
            })
    return pairs, summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument(
        "--morphology-root", type=Path,
        default=default_output_root() / "morphology",
    )
    parser.add_argument(
        "--output-root", type=Path, default=default_output_root() / "dynamics"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root, config = load_config(args.config)
    morphology_root = resolve_output_path(args.morphology_root)
    boundary_path = morphology_root / "morphology_boundaries.csv"
    if not boundary_path.exists():
        raise SystemExit(f"missing morphology product: {boundary_path}")
    rows = read_csv(boundary_path)
    driver = read_driver(root / config["data"]["driver"])  # type: ignore[index]
    landmarks = json.loads((morphology_root / "event_landmarks.json").read_text())
    compression = parse_utc(landmarks["compression"])
    main_phase = parse_utc(landmarks["main_phase"])

    harmonics, time_series = morphology_products(rows, config)
    lags = lag_products(harmonics, driver, config)
    pairs, hysteresis = hysteresis_products(
        rows, driver, compression, main_phase, config
    )
    output = resolve_output_path(args.output_root)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "morphology_harmonics.csv", harmonics)
    write_csv(output / "cutoff_dynamics_timeseries.csv", time_series)
    write_csv(output / "lag_correlations.csv", lags)
    write_csv(output / "hysteresis_pairs.csv", pairs)
    write_csv(output / "hysteresis_summary.csv", hysteresis)
    result = {
        "n_input_boundary_rows": len(rows), "n_harmonic_rows": len(harmonics),
        "n_lag_rows": len(lags), "n_hysteresis_pairs": len(pairs),
        "n_hysteresis_summary_rows": len(hysteresis), "event_landmarks": landmarks,
    }
    (output / "dynamics_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
