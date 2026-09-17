#!/usr/bin/env python3
"""Run VP02: Helios magnetic direction versus the SWCME Parker spiral."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import tarfile

import numpy as np

from download_data import FILE, verify
from reference_solution import fold_observed_angle_deg, parker_angle_deg


THRESHOLDS = {
    "minimum_daily_medians": 2000,
    "minimum_radial_bins": 10,
    "maximum_production_reference_angle_deg": 2.0e-11,
    "maximum_observed_median_absolute_error_deg": 20.0,
    "minimum_parker_polarity_fraction": 0.65,
}


@dataclass(frozen=True)
class Daily:
    """Equal-weight daily reduction after field-quality and sector filtering."""

    source: str
    radius_au: float
    latitude_deg: float
    speed_km_s: float
    angle_deg: float
    polarity_fraction: float
    samples: int


def artifact_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_daily(archive: Path) -> tuple[list[Daily], int]:
    """Stream the tarball and retain only stable, resolved Parker-sector samples.

    Status 1 is the corefit fully usable class.  The sigma-B ratio removes
    strongly varying measurement intervals, while |Br|/|B| excludes sector
    boundaries where folding by the sign of Br is ill-conditioned.  Requiring
    ten samples before a daily median prevents sparse days from having the same
    weight as a well observed day.
    """

    result: list[Daily] = []
    accepted = 0
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            if not member.isfile() or not member.name.endswith("_corefit.csv"):
                continue
            values: list[tuple[float, float, float, float, float]] = []
            binary = bundle.extractfile(member)
            if binary is None:
                continue
            reader = csv.DictReader(io.TextIOWrapper(binary, encoding="utf-8", newline=""))
            for row in reader:
                if row.get("Status", "").strip() != "1":
                    continue
                try:
                    radius = float(row["r_sun"])
                    latitude = float(row["clat"])
                    br, bt, bn = float(row["Bx"]), float(row["By"]), float(row["Bz"])
                    sigma_b, speed = float(row["sigma B"]), float(row["vp_x"])
                    bmag = math.sqrt(br * br + bt * bt + bn * bn)
                except (ValueError, KeyError):
                    continue
                if not all(math.isfinite(v) for v in (radius, latitude, br, bt, bn, sigma_b, speed, bmag)):
                    continue
                if not (0.29 <= radius <= 1.01 and 250.0 <= speed <= 850.0 and bmag > 0.0):
                    continue
                if sigma_b / bmag > 0.35 or abs(br) / bmag < 0.15:
                    continue
                angle = fold_observed_angle_deg(br, bt)
                values.append((radius, latitude, speed, angle, 1.0 if br * bt < 0.0 else 0.0))
                accepted += 1
            if len(values) >= 10:
                result.append(Daily(
                    member.name,
                    *(float(statistics.median(v[index] for v in values)) for index in range(4)),
                    float(sum(v[4] for v in values) / len(values)), len(values),
                ))
    result.sort(key=lambda value: value.source)
    return result, accepted


def make_bins(daily: list[Daily], count: int = 14) -> list[dict[str, float]]:
    """Create stable radial medians used as the model comparison coordinates."""

    edges = np.linspace(0.29, 1.01, count + 1)
    bins: list[dict[str, float]] = []
    for index in range(count):
        selected = [v for v in daily if edges[index] <= v.radius_au < edges[index + 1] or
                    (index == count - 1 and v.radius_au == edges[index + 1])]
        if len(selected) < 40:
            continue
        bins.append({
            "radius_au": float(statistics.median(v.radius_au for v in selected)),
            "latitude_deg": float(statistics.median(v.latitude_deg for v in selected)),
            "speed_km_s": float(statistics.median(v.speed_km_s for v in selected)),
            "observed_angle_deg": float(statistics.median(v.angle_deg for v in selected)),
            "observed_q16_deg": float(np.quantile([v.angle_deg for v in selected], 0.16)),
            "observed_q84_deg": float(np.quantile([v.angle_deg for v in selected], 0.84)),
            "polarity_fraction": float(sum(v.polarity_fraction * v.samples for v in selected) /
                                       sum(v.samples for v in selected)),
            "day_count": float(len(selected)),
        })
    return bins


def compile_driver(case_dir: Path, executable: Path) -> None:
    """Strictly compile the public-API adapter for standalone case execution."""

    source_root = case_dir.parents[2]
    command = ["c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
               "-isystem", str(source_root), str(case_dir / "vp02_model_driver.cpp"), "-o", str(executable)]
    subprocess.run(command, check=True)


def plot_comparison(path_stem: Path, rows: list[dict[str, float]]) -> None:
    """Write identical scientific content to the required PNG and EPS formats."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radius = np.array([row["radius_au"] for row in rows])
    observed = np.array([row["observed_angle_deg"] for row in rows])
    low = np.array([row["observed_q16_deg"] for row in rows])
    high = np.array([row["observed_q84_deg"] for row in rows])
    model = np.array([row["swcme_angle_deg"] for row in rows])
    reference = np.array([row["reference_angle_deg"] for row in rows])
    polarity = np.array([row["polarity_fraction"] for row in rows])
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
    axes[0].fill_between(radius, low, high, color="0.85", label="Helios 16–84%")
    axes[0].plot(radius, observed, "ko", label="Helios daily-bin median")
    axes[0].plot(radius, model, "C1-", lw=2, label="SWCME public API")
    axes[0].plot(radius, reference, "C0--", lw=1.5, label="independent Parker reference")
    axes[0].set_ylabel("sector-folded angle [deg]")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.25)
    axes[1].plot(radius, np.abs(observed - model), "C3o-", label="absolute angle error")
    axes[1].plot(radius, 100.0 * polarity, "C2s-", label="Parker polarity fraction [%]")
    axes[1].axhline(THRESHOLDS["maximum_observed_median_absolute_error_deg"], color="C3", ls=":")
    axes[1].set_xlabel("heliocentric radius [AU]")
    axes[1].set_ylabel("error [deg] / polarity [%]")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    for suffix in ("png", "eps"):
        fig.savefig(path_stem.with_suffix("." + suffix), dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    args = parser.parse_args()
    case_dir = Path(__file__).resolve().parent
    archive = case_dir / "data/raw/corefit.gz"
    if args.download:
        subprocess.run([sys.executable, str(case_dir / "download_data.py")], check=True)
    if not archive.is_file():
        raise RuntimeError(f"missing {archive}; rerun with --download")
    verify(archive)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    daily, accepted = read_daily(archive)
    bins = make_bins(daily)
    input_csv = args.output_dir / "vp02_model_input.csv"
    with input_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["radius_au", "latitude_deg", "speed_km_s"], lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in bins)
    executable = args.output_dir / "vp02_model_driver"
    model_csv = args.output_dir / "vp02_model_output.csv"
    compile_driver(case_dir, executable)
    subprocess.run([str(executable), str(input_csv), str(model_csv)], check=True)
    model_rows = list(csv.DictReader(model_csv.open(encoding="utf-8")))
    if len(model_rows) != len(bins):
        raise RuntimeError("model output cardinality mismatch")

    comparison: list[dict[str, float]] = []
    for observed, produced in zip(bins, model_rows):
        row = dict(observed)
        row["swcme_angle_deg"] = float(produced["parker_angle_deg"])
        row["reference_angle_deg"] = parker_angle_deg(row["radius_au"], row["latitude_deg"], row["speed_km_s"])
        comparison.append(row)
    comparison_csv = args.output_dir / "vp02_comparison.csv"
    with comparison_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(comparison)
    # Preserve the oracle as a standalone portable table as well as paired
    # columns in the comparison file. Reviewers can therefore regenerate or
    # plot the reference without executing the production driver.
    with (args.output_dir / "vp02_reference_solution.csv").open("w", encoding="utf-8", newline="") as stream:
        fields = ["radius_au", "latitude_deg", "speed_km_s", "reference_angle_deg"]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows({key: row[key] for key in fields} for row in comparison)

    production_reference = max(abs(row["swcme_angle_deg"] - row["reference_angle_deg"]) for row in comparison)
    daily_errors = [abs(value.angle_deg - parker_angle_deg(value.radius_au, value.latitude_deg, value.speed_km_s)) for value in daily]
    observed_error = float(statistics.median(daily_errors))
    polarity_fraction = float(sum(value.polarity_fraction * value.samples for value in daily) /
                              sum(value.samples for value in daily))
    metrics = {
        "accepted_samples": accepted,
        "daily_medians": len(daily),
        "radial_bins": len(bins),
        "maximum_production_reference_angle_deg": production_reference,
        "observed_median_absolute_error_deg": observed_error,
        "parker_polarity_fraction": polarity_fraction,
    }
    criteria = {
        "daily_sample_size": len(daily) >= THRESHOLDS["minimum_daily_medians"],
        "radial_coverage": len(bins) >= THRESHOLDS["minimum_radial_bins"],
        "independent_reference": production_reference <= THRESHOLDS["maximum_production_reference_angle_deg"],
        "observed_angle": observed_error <= THRESHOLDS["maximum_observed_median_absolute_error_deg"],
        "observed_polarity": polarity_fraction >= THRESHOLDS["minimum_parker_polarity_fraction"],
    }
    if not args.no_plots:
        plot_comparison(args.output_dir / "vp02_parker_comparison", comparison)
    status = "PASS" if all(criteria.values()) else "FAIL"
    result = {
        "schema_version": 1, "validation_id": "VP02", "status": status,
        "description": "Helios Parker angle and polarity comparison",
        "data_sha256": FILE["sha256"], "thresholds": THRESHOLDS,
        "metrics": metrics, "criteria": criteria,
        "limitations": ["Helios corefit is not explicitly tagged by ICME/stream context; robust daily medians describe natural background scatter.",
                        "The polarity metric tests the Parker quadrant after sector folding, not a solar-cycle sector forecast."],
    }
    result_path = args.output_dir / "vp02_result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    artifacts = [{"name": path.name, "bytes": path.stat().st_size, "sha256": artifact_hash(path)}
                 for path in sorted(args.output_dir.iterdir())
                 if path.is_file() and path != executable and path.name != "vp02_artifact_manifest.json"]
    (args.output_dir / "vp02_artifact_manifest.json").write_text(json.dumps({"schema_version": 1, "artifacts": artifacts}, indent=2) + "\n", encoding="utf-8")
    print(f"VP02 {status}: median angle error={observed_error:.3f} deg; polarity={polarity_fraction:.3f}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f"VP02 setup/analysis error: {error}", file=sys.stderr)
        raise SystemExit(2)
