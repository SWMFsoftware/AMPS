#!/usr/bin/env python3
"""Run VP04: SWCME data-driven kinematics against CDAW height-time traces."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

from download_data import FILES, verify
from reference_solution import evaluate

SOLAR_RADIUS_KM = 695_700.0
THRESHOLDS = {
    "minimum_events": 6,
    "minimum_heldout_points": 30,
    "maximum_production_reference_radius_rs": 5.0e-12,
    "maximum_production_reference_speed_km_s": 5.0e-9,
    "maximum_fit_knot_residual_rs": 5.0e-12,
    "maximum_heldout_median_error_rs": 0.35,
    "maximum_heldout_worst_error_rs": 1.25,
}


@dataclass(frozen=True)
class Event:
    """One immutable CDAW manual leading-edge height-time sequence."""

    event_id: str
    catalog_speed_km_s: float
    times_s: list[float]
    radii_rs: list[float]


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest


def read_event(path: Path) -> Event:
    """Parse the documented YHT columns and metadata without heuristic dates."""

    timestamps: list[datetime] = []
    radii: list[float] = []
    speed = math.nan
    event_id = path.name.split(".w360h", 1)[0]
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("#SPEED:"):
            speed = float(line.split(":", 1)[1].strip())
        elif line and not line.startswith("#"):
            fields = line.split()
            if len(fields) < 3:
                raise RuntimeError(f"malformed height-time row in {path}")
            radii.append(float(fields[0]))
            timestamps.append(datetime.strptime(fields[1] + " " + fields[2], "%Y/%m/%d %H:%M:%S"))
    if len(radii) < 6 or not math.isfinite(speed):
        raise RuntimeError(f"incomplete event file {path}")
    if any(radii[index + 1] < radii[index] for index in range(len(radii) - 1)):
        raise RuntimeError(f"nonmonotone leading-edge trace {path}")
    times = [(stamp - timestamps[0]).total_seconds() for stamp in timestamps]
    if any(times[index + 1] <= times[index] for index in range(len(times) - 1)):
        raise RuntimeError(f"nonincreasing timestamps in {path}")
    return Event(event_id, speed, times, radii)


def compile_driver(case_dir: Path, executable: Path) -> None:
    source_root = case_dir.parents[2]
    subprocess.run(["c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                    "-isystem", str(source_root), str(case_dir / "vp04_model_driver.cpp"), "-o", str(executable)], check=True)


def evaluate_event(event: Event, executable: Path, work: Path) -> list[dict[str, object]]:
    """Fit alternate knots, keeping every intervening point blind to SWCME."""

    fit_indices = sorted(set([0, len(event.times_s) - 1] + list(range(0, len(event.times_s), 2))))
    fit_times = [event.times_s[index] for index in fit_indices]
    fit_radii = [event.radii_rs[index] for index in fit_indices]
    fit_path, query_path, output_path = work / "fit.csv", work / "query.csv", work / "model.csv"
    with fit_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(["time_s", "radius_rs"])
        writer.writerows(zip(fit_times, fit_radii))
    with query_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(["time_s", "observed_radius_rs"])
        writer.writerows(zip(event.times_s, event.radii_rs))
    subprocess.run([str(executable), str(fit_path), str(query_path), str(output_path)], check=True)
    produced = list(csv.DictReader(output_path.open(encoding="utf-8")))
    if len(produced) != len(event.times_s):
        raise RuntimeError(f"model output cardinality mismatch for {event.event_id}")
    rows: list[dict[str, object]] = []
    for index, model in enumerate(produced):
        radius, derivative = evaluate(fit_times, fit_radii, event.times_s[index])
        rows.append({
            "event_id": event.event_id,
            "time_s": event.times_s[index],
            "observed_radius_rs": event.radii_rs[index],
            "is_fit_knot": index in fit_indices,
            "swcme_radius_rs": float(model["swcme_radius_rs"]),
            "reference_radius_rs": radius,
            "swcme_speed_km_s": float(model["swcme_speed_km_s"]),
            "reference_speed_km_s": derivative * SOLAR_RADIUS_KM,
            "catalog_speed_km_s": event.catalog_speed_km_s,
        })
    return rows


def plot(stem: Path, events: list[Event], rows: list[dict[str, object]]) -> None:
    """Draw all six independent event traces in a compact comparison panel."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(10.0, 10.0))
    for axis, event in zip(axes.flat, events):
        selected = [row for row in rows if row["event_id"] == event.event_id]
        hours = [float(row["time_s"]) / 3600.0 for row in selected]
        axis.plot(hours, [row["reference_radius_rs"] for row in selected], "C0--", label="independent PCHIP")
        axis.plot(hours, [row["swcme_radius_rs"] for row in selected], "C1-", lw=1.8, label="SWCME")
        held = [row for row in selected if not bool(row["is_fit_knot"])]
        fit = [row for row in selected if bool(row["is_fit_knot"])]
        axis.plot([float(row["time_s"])/3600 for row in fit], [row["observed_radius_rs"] for row in fit], "ko", ms=4, label="fit knots")
        axis.plot([float(row["time_s"])/3600 for row in held], [row["observed_radius_rs"] for row in held], "rx", ms=5, label="held out")
        axis.set_title(event.event_id); axis.set_xlabel("elapsed time [h]"); axis.set_ylabel("height [solar radii]")
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7)
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
    if args.download:
        subprocess.run([sys.executable, str(case_dir / "download_data.py")], check=True)
    raw = case_dir / "data/raw"
    for name, _, size, digest in FILES:
        path = raw / name
        if not path.is_file():
            raise RuntimeError(f"missing {path}; rerun with --download")
        verify(path, size, digest)
    events = [read_event(raw / spec[0]) for spec in FILES]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    executable = args.output_dir / "vp04_model_driver"
    compile_driver(case_dir, executable)
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="vp04-") as temporary:
        for event in events:
            work = Path(temporary) / event.event_id; work.mkdir()
            rows.extend(evaluate_event(event, executable, work))
    comparison = args.output_dir / "vp04_comparison.csv"
    with comparison.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    reference_fields = ["event_id", "time_s", "reference_radius_rs", "reference_speed_km_s"]
    with (args.output_dir / "vp04_reference_solution.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=reference_fields, lineterminator="\n")
        writer.writeheader(); writer.writerows({key: row[key] for key in reference_fields} for row in rows)

    heldout = [row for row in rows if not bool(row["is_fit_knot"])]
    knots = [row for row in rows if bool(row["is_fit_knot"])]
    production_radius_error = max(abs(float(row["swcme_radius_rs"]) - float(row["reference_radius_rs"])) for row in rows)
    production_speed_error = max(abs(float(row["swcme_speed_km_s"]) - float(row["reference_speed_km_s"])) for row in rows)
    knot_error = max(abs(float(row["swcme_radius_rs"]) - float(row["observed_radius_rs"])) for row in knots)
    heldout_errors = [abs(float(row["swcme_radius_rs"]) - float(row["observed_radius_rs"])) for row in heldout]
    metrics = {"event_count": len(events), "measurement_count": len(rows), "heldout_point_count": len(heldout),
               "maximum_production_reference_radius_rs": production_radius_error,
               "maximum_production_reference_speed_km_s": production_speed_error,
               "maximum_fit_knot_residual_rs": knot_error,
               "heldout_median_absolute_error_rs": statistics.median(heldout_errors),
               "heldout_worst_absolute_error_rs": max(heldout_errors)}
    criteria = {"event_count": len(events) >= THRESHOLDS["minimum_events"],
                "heldout_count": len(heldout) >= THRESHOLDS["minimum_heldout_points"],
                "reference_radius": production_radius_error <= THRESHOLDS["maximum_production_reference_radius_rs"],
                "reference_speed": production_speed_error <= THRESHOLDS["maximum_production_reference_speed_km_s"],
                "knot_interpolation": knot_error <= THRESHOLDS["maximum_fit_knot_residual_rs"],
                "heldout_median": statistics.median(heldout_errors) <= THRESHOLDS["maximum_heldout_median_error_rs"],
                "heldout_worst": max(heldout_errors) <= THRESHOLDS["maximum_heldout_worst_error_rs"]}
    if not args.no_plots:
        plot(args.output_dir / "vp04_height_time_comparison", events, rows)
    status = "PASS" if all(criteria.values()) else "FAIL"
    result = {"schema_version": 1, "validation_id": "VP04", "status": status,
              "description": "CDAW CME leading-edge height-time validation",
              "metrics": metrics, "thresholds": THRESHOLDS, "criteria": criteria,
              "events": [event.event_id for event in events],
              "limitations": ["CDAW heights are plane-of-sky leading-edge measurements, not deprojected three-dimensional apex distances.",
                              "Alternate-point withholding tests interpolation skill inside the measured coronagraph interval; it does not validate interplanetary extrapolation."]}
    (args.output_dir / "vp04_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    artifacts = [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                 for path in args.output_dir.iterdir()
                 if path.is_file() and path != executable and path.name != "vp04_artifact_manifest.json"]
    (args.output_dir / "vp04_artifact_manifest.json").write_text(json.dumps({"schema_version": 1, "artifacts": sorted(artifacts, key=lambda x: x["path"])}, indent=2) + "\n", encoding="utf-8")
    print(f"VP04 {status}: {len(events)} events, heldout median={statistics.median(heldout_errors):.4f} Rs")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP04 setup/analysis error: {error}", file=sys.stderr)
        raise SystemExit(2)
