#!/usr/bin/env python3
"""Run VP05: inner-to-outer multipoint ICME transit validation."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

import numpy as np

from download_data import FILE, verify
from reference_solution import arrival

THRESHOLDS = {
    "minimum_pair_count": 15,
    "maximum_production_reference_transit_s": 2.0e-8,
    "maximum_production_reference_speed_km_s": 2.0e-10,
    "maximum_observed_median_arrival_error_h": 12.0,
    "maximum_observed_median_speed_error_fraction": 0.25,
}
ALLOWED_SPACECRAFT = {"PSP", "SolarOrbiter", "STEREO-A", "Wind", "BepiColombo", "MESSENGER", "VenusExpress"}


@dataclass(frozen=True)
class Observation:
    """One catalog event boundary with usable radial speed and geometry."""

    event_id: str
    spacecraft: str
    time: datetime
    radius_au: float
    longitude_deg: float
    speed_km_s: float


@dataclass(frozen=True)
class Pair:
    """Chronological, near-radial inner/outer observation pair."""

    pair_id: str
    event_id: str
    inner: Observation
    outer: Observation
    longitude_separation_deg: float


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def circular_separation(a: float, b: float) -> float:
    difference = abs(a - b) % 360.0
    return min(difference, 360.0 - difference)


def select_pairs(path: Path) -> list[Pair]:
    """Apply only frozen geometry/completeness rules—never outcome filtering.

    The outer arrival time and speed are not read by the production driver;
    they remain held-out targets. Multiple catalog boundaries for the same
    spacecraft are retained as distinct hypotheses rather than silently
    choosing whichever produces the smallest propagation error.
    """

    groups: dict[str, list[Observation]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row.get("spacecraft") not in ALLOWED_SPACECRAFT or row.get("speed") in (None, "", "-"):
                continue
            try:
                observation = Observation(
                    row["lineupcat_id"], row["spacecraft"],
                    datetime.fromisoformat(row["event_start_time"].replace("Z", "+00:00")),
                    float(row["sc_heliodistance"]), float(row["sc_heeq_lon"]), float(row["speed"]),
                )
            except (KeyError, ValueError):
                continue
            if all(math.isfinite(value) for value in (observation.radius_au, observation.longitude_deg, observation.speed_km_s)):
                groups[observation.event_id].append(observation)
    result: list[Pair] = []
    for event_id in sorted(groups):
        candidates = sorted(groups[event_id], key=lambda value: (value.radius_au, value.time, value.spacecraft))
        for inner in candidates:
            for outer in candidates:
                separation = circular_separation(inner.longitude_deg, outer.longitude_deg)
                if (outer.radius_au - inner.radius_au < 0.08 or outer.time <= inner.time or separation > 12.0):
                    continue
                pair_id = f"P{len(result)+1:03d}_{event_id}_{inner.spacecraft}_{outer.spacecraft}"
                result.append(Pair(pair_id, event_id, inner, outer, separation))
    return result


def compile_driver(case_dir: Path, executable: Path) -> None:
    source_root = case_dir.parents[2]
    subprocess.run(["c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                    "-isystem", str(source_root), str(case_dir / "vp05_model_driver.cpp"), "-o", str(executable)], check=True)


def plot(stem: Path, rows: list[dict[str, object]]) -> None:
    """Compare transit duration and outer speed against observations and oracle."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    observed_hours = np.array([float(row["observed_transit_h"]) for row in rows])
    predicted_hours = np.array([float(row["swcme_transit_h"]) for row in rows])
    observed_speed = np.array([float(row["observed_outer_speed_km_s"]) for row in rows])
    predicted_speed = np.array([float(row["swcme_outer_speed_km_s"]) for row in rows])
    reference_hours = np.array([float(row["reference_transit_h"]) for row in rows])
    reference_speed = np.array([float(row["reference_outer_speed_km_s"]) for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.6))
    limit = max(observed_hours.max(), predicted_hours.max()) * 1.05
    axes[0].plot([0, limit], [0, limit], "k:", label="perfect observation")
    axes[0].scatter(observed_hours, predicted_hours, color="C1", label="SWCME DBM")
    axes[0].scatter(observed_hours, reference_hours, facecolors="none", edgecolors="C0", label="independent DBM")
    axes[0].set(xlabel="observed transit [h]", ylabel="predicted transit [h]")
    limit_speed = max(observed_speed.max(), predicted_speed.max()) * 1.05
    axes[1].plot([0, limit_speed], [0, limit_speed], "k:")
    axes[1].scatter(observed_speed, predicted_speed, color="C1", label="SWCME DBM")
    axes[1].scatter(observed_speed, reference_speed, facecolors="none", edgecolors="C0", label="independent DBM")
    axes[1].set(xlabel="observed outer speed [km/s]", ylabel="predicted outer speed [km/s]")
    for axis in axes:
        axis.grid(alpha=0.25); axis.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ("png", "eps"):
        fig.savefig(stem.with_suffix("." + suffix), dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    args = parser.parse_args(); case_dir = Path(__file__).resolve().parent
    data = case_dir / "data/raw" / str(FILE["name"])
    if args.download:
        subprocess.run([sys.executable, str(case_dir / "download_data.py")], check=True)
    if not data.is_file():
        raise RuntimeError(f"missing {data}; rerun with --download")
    verify(data)
    pairs = select_pairs(data)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_csv = args.output_dir / "vp05_model_input.csv"
    with input_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["pair_id", "radius0_au", "speed0_km_s", "target_au"])
        writer.writerows((pair.pair_id, pair.inner.radius_au, pair.inner.speed_km_s, pair.outer.radius_au) for pair in pairs)
    executable, model_csv = args.output_dir / "vp05_model_driver", args.output_dir / "vp05_model_output.csv"
    compile_driver(case_dir, executable)
    subprocess.run([str(executable), str(input_csv), str(model_csv)], check=True)
    production = {row["pair_id"]: row for row in csv.DictReader(model_csv.open(encoding="utf-8"))}
    if set(production) != {pair.pair_id for pair in pairs}:
        raise RuntimeError("model output identity mismatch")

    rows: list[dict[str, object]] = []
    for pair in pairs:
        model = production[pair.pair_id]
        reference_time, reference_speed = arrival(pair.inner.radius_au, pair.inner.speed_km_s, pair.outer.radius_au)
        observed_time = (pair.outer.time - pair.inner.time).total_seconds()
        rows.append({
            "pair_id": pair.pair_id, "event_id": pair.event_id,
            "inner_spacecraft": pair.inner.spacecraft, "outer_spacecraft": pair.outer.spacecraft,
            "longitude_separation_deg": pair.longitude_separation_deg,
            "inner_radius_au": pair.inner.radius_au, "outer_radius_au": pair.outer.radius_au,
            "inner_speed_km_s": pair.inner.speed_km_s, "observed_outer_speed_km_s": pair.outer.speed_km_s,
            "observed_transit_h": observed_time / 3600.0,
            "swcme_transit_h": float(model["predicted_transit_s"]) / 3600.0,
            "reference_transit_h": reference_time / 3600.0,
            "swcme_outer_speed_km_s": float(model["predicted_speed_km_s"]),
            "reference_outer_speed_km_s": reference_speed,
        })
    comparison = args.output_dir / "vp05_comparison.csv"
    with comparison.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    reference_fields = ["pair_id", "reference_transit_h", "reference_outer_speed_km_s"]
    with (args.output_dir / "vp05_reference_solution.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=reference_fields, lineterminator="\n")
        writer.writeheader(); writer.writerows({key: row[key] for key in reference_fields} for row in rows)

    reference_time_error_s = max(abs(float(row["swcme_transit_h"]) - float(row["reference_transit_h"])) * 3600.0 for row in rows)
    reference_speed_error = max(abs(float(row["swcme_outer_speed_km_s"]) - float(row["reference_outer_speed_km_s"])) for row in rows)
    arrival_errors = [abs(float(row["swcme_transit_h"]) - float(row["observed_transit_h"])) for row in rows]
    speed_errors = [abs(float(row["swcme_outer_speed_km_s"]) - float(row["observed_outer_speed_km_s"])) /
                    float(row["observed_outer_speed_km_s"]) for row in rows]
    metrics = {"pair_count": len(rows), "event_count": len({row["event_id"] for row in rows}),
               "maximum_production_reference_transit_s": reference_time_error_s,
               "maximum_production_reference_speed_km_s": reference_speed_error,
               "observed_median_arrival_error_h": statistics.median(arrival_errors),
               "observed_arrival_error_q84_h": float(np.quantile(arrival_errors, 0.84)),
               "observed_median_speed_error_fraction": statistics.median(speed_errors)}
    criteria = {"pair_count": len(rows) >= THRESHOLDS["minimum_pair_count"],
                "reference_transit": reference_time_error_s <= THRESHOLDS["maximum_production_reference_transit_s"],
                "reference_speed": reference_speed_error <= THRESHOLDS["maximum_production_reference_speed_km_s"],
                "observed_arrival": statistics.median(arrival_errors) <= THRESHOLDS["maximum_observed_median_arrival_error_h"],
                "observed_speed": statistics.median(speed_errors) <= THRESHOLDS["maximum_observed_median_speed_error_fraction"]}
    if not args.no_plots:
        plot(args.output_dir / "vp05_multipoint_comparison", rows)
    status = "PASS" if all(criteria.values()) else "FAIL"
    result = {"schema_version": 1, "validation_id": "VP05", "status": status,
              "description": "HELIO4CAST near-radial multipoint DBM propagation validation",
              "metrics": metrics, "thresholds": THRESHOLDS, "criteria": criteria,
              "data_sha256": FILE["sha256"],
              "selection": {"maximum_longitude_separation_deg": 12.0, "minimum_radial_separation_au": 0.08,
                            "required": "numeric speed and chronological outer arrival"},
              "limitations": ["LineupCAT event-start and catalog speed are ICME boundary/mean-speed proxies, not uniformly re-fitted shock onsets and shock speeds.",
                              "A one-dimensional apex DBM cannot represent flank curvature or longitudinal evolution; the 12-degree rule only limits that mismatch."]}
    (args.output_dir / "vp05_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    artifacts = [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                 for path in args.output_dir.iterdir()
                 if path.is_file() and path != executable and path.name != "vp05_artifact_manifest.json"]
    (args.output_dir / "vp05_artifact_manifest.json").write_text(json.dumps({"schema_version": 1, "artifacts": sorted(artifacts, key=lambda x: x["path"])}, indent=2) + "\n", encoding="utf-8")
    print("VP05 {}: pairs={} median arrival error={:.2f} h median speed error={:.1f}%".format(
        status, len(rows), statistics.median(arrival_errors), 100.0 * statistics.median(speed_errors)))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP05 setup/analysis error: {error}", file=sys.stderr)
        raise SystemExit(2)
