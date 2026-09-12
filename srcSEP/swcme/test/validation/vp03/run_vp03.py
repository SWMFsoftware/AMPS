#!/usr/bin/env python3
"""Run VP03 against Helios background plasma and magnetic moments."""

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
from reference_solution import background_state, observed_state


THRESHOLDS = {
    "minimum_daily_medians": 2500,
    "minimum_radial_bins": 12,
    "maximum_production_reference_relative_error": 5.0e-12,
    "maximum_speed_median_factor": 1.35,
    "maximum_pressure_median_factor": 2.0,
    "maximum_sound_speed_median_factor": 1.5,
    "maximum_fast_speed_median_factor": 1.5,
}


@dataclass(frozen=True)
class Daily:
    """Daily median of a fully usable Helios proton-core interval."""

    source: str
    radius_au: float
    latitude_deg: float
    density_cm3: float
    speed_km_s: float
    pressure_pa: float
    sound_speed_km_s: float
    alfven_speed_km_s: float
    fast_speed_km_s: float
    samples: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_daily(archive: Path) -> tuple[list[Daily], int]:
    """Reduce quality-controlled RTN proton-core moments without extraction.

    Scalar temperature is `(T_parallel + 2*T_perpendicular)/3`, the trace of
    the measured gyrotropic tensor.  Characteristic speeds use measured proton
    density and |B| but the explicitly documented proton-only closure, making
    comparison with the default SWCME thermodynamics well-defined.
    """

    daily: list[Daily] = []
    accepted = 0
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle:
            if not member.isfile() or not member.name.endswith("_corefit.csv"):
                continue
            values: list[tuple[float, ...]] = []
            binary = bundle.extractfile(member)
            if binary is None:
                continue
            reader = csv.DictReader(io.TextIOWrapper(binary, encoding="utf-8", newline=""))
            for row in reader:
                if row.get("Status", "").strip() != "1":
                    continue
                try:
                    radius, latitude = float(row["r_sun"]), float(row["clat"])
                    density = float(row["n_p"])
                    vr, vt, vn = float(row["vp_x"]), float(row["vp_y"]), float(row["vp_z"])
                    br, bt, bn = float(row["Bx"]), float(row["By"]), float(row["Bz"])
                    temperature = (float(row["Tp_par"]) + 2.0 * float(row["Tp_perp"])) / 3.0
                    speed = math.sqrt(vr * vr + vt * vt + vn * vn)
                    bmag = math.sqrt(br * br + bt * bt + bn * bn)
                except (ValueError, KeyError):
                    continue
                scalars = (radius, latitude, density, speed, bmag, temperature)
                if not all(math.isfinite(value) for value in scalars):
                    continue
                if not (0.29 <= radius <= 1.01 and 0.0 < density < 1000.0 and
                        200.0 <= speed <= 900.0 and 2.0e3 <= temperature <= 2.0e6 and bmag > 0.0):
                    continue
                closure = observed_state(density, temperature, bmag)
                values.append((radius, latitude, density, speed, closure["pressure_pa"],
                               closure["sound_speed_km_s"], closure["alfven_speed_km_s"],
                               closure["fast_speed_km_s"]))
                accepted += 1
            if len(values) >= 10:
                medians = [float(statistics.median(row[index] for row in values)) for index in range(8)]
                daily.append(Daily(member.name, *medians, len(values)))
    daily.sort(key=lambda value: value.source)
    return daily, accepted


def radial_bins(daily: list[Daily], count: int = 14) -> list[dict[str, float]]:
    """Build equal-width robust summaries; each contributing day is one vote."""

    edges = np.linspace(0.29, 1.01, count + 1)
    result: list[dict[str, float]] = []
    names = ("radius_au", "latitude_deg", "density_cm3", "speed_km_s", "pressure_pa",
             "sound_speed_km_s", "alfven_speed_km_s", "fast_speed_km_s")
    for index in range(count):
        selected = [value for value in daily if edges[index] <= value.radius_au < edges[index + 1] or
                    (index == count - 1 and value.radius_au == edges[index + 1])]
        if len(selected) < 40:
            continue
        row = {"observed_" + name: float(statistics.median(getattr(value, name) for value in selected))
               for name in names[2:]}
        row.update({"radius_au": float(statistics.median(value.radius_au for value in selected)),
                    "latitude_deg": float(statistics.median(value.latitude_deg for value in selected)),
                    "day_count": float(len(selected))})
        result.append(row)
    return result


def factor(a: float, b: float) -> float:
    """Symmetric multiplicative discrepancy, invariant to choice of numerator."""

    return max(a / b, b / a)


def compile_driver(case_dir: Path, executable: Path) -> None:
    source_root = case_dir.parents[2]
    subprocess.run(["c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                    "-isystem", str(source_root), str(case_dir / "vp03_model_driver.cpp"), "-o", str(executable)], check=True)


def plot(stem: Path, rows: list[dict[str, float]]) -> None:
    """Render model, independent reference, and observational bin medians."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radius = np.array([row["radius_au"] for row in rows])
    panels = (("speed_km_s", "speed [km s$^{-1}$]"),
              ("pressure_pa", "proton pressure [Pa]"),
              ("sound_speed_km_s", "sound speed [km s$^{-1}$]"),
              ("fast_speed_km_s", "fast proxy [km s$^{-1}$]"))
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), sharex=True)
    for axis, (name, label) in zip(axes.flat, panels):
        axis.plot(radius, [row["observed_" + name] for row in rows], "ko", label="Helios")
        axis.plot(radius, [row["swcme_" + name] for row in rows], "C1-", lw=2, label="SWCME")
        axis.plot(radius, [row["reference_" + name] for row in rows], "C0--", label="independent reference")
        axis.set_ylabel(label); axis.grid(alpha=0.25)
        if name == "pressure_pa": axis.set_yscale("log")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].set_xlabel("radius [AU]"); axes[1, 1].set_xlabel("radius [AU]")
    fig.tight_layout()
    for suffix in ("png", "eps"):
        fig.savefig(stem.with_suffix("." + suffix), dpi=180)
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
    bins = radial_bins(daily)

    input_csv = args.output_dir / "vp03_model_input.csv"
    with input_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["radius_au", "latitude_deg"], lineterminator="\n")
        writer.writeheader(); writer.writerows({key: row[key] for key in writer.fieldnames} for row in bins)
    executable, model_csv = args.output_dir / "vp03_model_driver", args.output_dir / "vp03_model_output.csv"
    compile_driver(case_dir, executable)
    subprocess.run([str(executable), str(input_csv), str(model_csv)], check=True)
    model = list(csv.DictReader(model_csv.open(encoding="utf-8")))
    if len(model) != len(bins):
        raise RuntimeError("model output cardinality mismatch")

    quantities = ("density_cm3", "speed_km_s", "pressure_pa", "sound_speed_km_s", "alfven_speed_km_s", "fast_speed_km_s")
    rows: list[dict[str, float]] = []
    for observed, produced in zip(bins, model):
        row = dict(observed)
        reference = background_state(row["radius_au"], row["latitude_deg"])
        for name in quantities:
            row["swcme_" + name] = float(produced[name])
            row["reference_" + name] = reference[name]
        rows.append(row)
    comparison_csv = args.output_dir / "vp03_comparison.csv"
    with comparison_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    reference_fields = ["radius_au", "latitude_deg"] + ["reference_" + name for name in quantities]
    with (args.output_dir / "vp03_reference_solution.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=reference_fields, lineterminator="\n")
        writer.writeheader(); writer.writerows({key: row[key] for key in reference_fields} for row in rows)

    relative_errors = [abs(row["swcme_" + name] - row["reference_" + name]) /
                       max(abs(row["reference_" + name]), 1.0e-300) for row in rows for name in quantities]
    median_factors = {name: float(statistics.median(
        factor(row["observed_" + name], row["swcme_" + name]) for row in rows))
        for name in ("speed_km_s", "pressure_pa", "sound_speed_km_s", "fast_speed_km_s")}
    metrics = {"accepted_samples": accepted, "daily_medians": len(daily), "radial_bins": len(bins),
               "maximum_production_reference_relative_error": max(relative_errors),
               **{"median_factor_" + name: value for name, value in median_factors.items()}}
    criteria = {
        "daily_sample_size": len(daily) >= THRESHOLDS["minimum_daily_medians"],
        "radial_coverage": len(bins) >= THRESHOLDS["minimum_radial_bins"],
        "independent_reference": max(relative_errors) <= THRESHOLDS["maximum_production_reference_relative_error"],
        "speed": median_factors["speed_km_s"] <= THRESHOLDS["maximum_speed_median_factor"],
        "pressure": median_factors["pressure_pa"] <= THRESHOLDS["maximum_pressure_median_factor"],
        "sound_speed": median_factors["sound_speed_km_s"] <= THRESHOLDS["maximum_sound_speed_median_factor"],
        "fast_speed": median_factors["fast_speed_km_s"] <= THRESHOLDS["maximum_fast_speed_median_factor"],
    }
    if not args.no_plots:
        plot(args.output_dir / "vp03_background_comparison", rows)
    status = "PASS" if all(criteria.values()) else "FAIL"
    result = {"schema_version": 1, "validation_id": "VP03", "status": status,
              "description": "Helios speed, pressure, and characteristic-speed validation",
              "metrics": metrics, "thresholds": THRESHOLDS, "criteria": criteria,
              "data_sha256": FILE["sha256"],
              "limitations": ["The observed characteristic speeds use proton-core moments and the proton-only closure; electron and alpha pressure are outside this case.",
                              "The scalar fast proxy sqrt(cs^2+vA^2) is an upper branch, not an angle-resolved wave-mode measurement."]}
    (args.output_dir / "vp03_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    artifacts = [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                 for path in args.output_dir.iterdir()
                 if path.is_file() and path != executable and path.name != "vp03_artifact_manifest.json"]
    (args.output_dir / "vp03_artifact_manifest.json").write_text(json.dumps({"schema_version": 1, "artifacts": sorted(artifacts, key=lambda x: x["path"])}, indent=2) + "\n", encoding="utf-8")
    print("VP03 {}: speed={:.3f} pressure={:.3f} sound={:.3f} fast={:.3f}".format(
        status, median_factors["speed_km_s"], median_factors["pressure_pa"],
        median_factors["sound_speed_km_s"], median_factors["fast_speed_km_s"]))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP03 setup/analysis error: {error}", file=sys.stderr)
        raise SystemExit(2)
