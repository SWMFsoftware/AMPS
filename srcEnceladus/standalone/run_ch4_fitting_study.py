#!/usr/bin/env python3
"""Fit calibrated INMS CH4 mass-15 profiles with the standalone plume model."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

APP = Path(__file__).resolve().parent
sys.path.insert(0, str(APP))
from fitting import SOURCE_CASES, load_config, matrix_from_config, solve_amplitudes, sources_from_config  # noqa:E402
from joint_fitting_study import fit_rate_for_fixed_source_shape, normalized_source_shape, solve_combined_blocks  # noqa:E402

M_CH4 = 16.04246 * 1.66053906660e-27
FIT_FLYBYS = ("E3", "E5", "E7", "E14", "E17", "E18", "E21")
PRIMARY = ("E14", "E17", "E18")
LOW_SPEED = ("E7", "E14", "E17", "E18", "E21")
CALIBRATION_SCALE_FRACTION = 0.23
M15_SHARE_RANGE = (0.80, 1.00)
M15_SHARE_NOMINAL = 0.90


def json_default(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def read_profile(path: Path, include_flagged: bool = False) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not include_flagged:
        rows = [r for r in rows if r.get("fit_usable", "1") == "1"]
    points = np.asarray([[float(r[k]) for k in ("spacecraft_x_m", "spacecraft_y_m", "spacecraft_z_m")] for r in rows])
    return {"rows": rows, "points": points,
            "time": np.asarray([float(r["time_ca_s"]) for r in rows]),
            "density": np.asarray([float(r["density_m3"]) for r in rows]),
            "net": np.asarray([float(r["density_background_subtracted_m3"]) for r in rows]),
            "sigma": np.asarray([float(r["sigma_counting_m3"]) for r in rows])}


def source_total(vector: np.ndarray, labels: list[str]) -> float:
    return float(sum(v for l, v in zip(labels, vector) if l.startswith(("jet:", "stripe:"))) * 1e28)


def direct_result(fit: dict[str, object], labels: list[str], n: int) -> dict[str, object]:
    vector = np.asarray(fit["coefficients"])
    total = source_total(vector, labels)
    return {"total_source_rate_s": total, "chi2": float(fit["chi2"]),
            "chi2_per_sample": float(fit["chi2"] / n), "samples": n,
            "rank": int(fit["rank"]), "condition_number_scaled": float(fit["condition_number_scaled"]),
            "full_coefficients": vector.tolist()}


def combined_shape(name: str, training: tuple[str, ...], data, matrices, labels):
    combined = solve_combined_blocks([(f, matrices[f], data[f]["density"], data[f]["sigma"]) for f in training],
                                     balance_by_dataset=True)
    # Force validation now; an unobservable zero-source combined fit is invalid.
    normalized_source_shape(combined["coefficients"], labels)
    targets = {f: fit_rate_for_fixed_source_shape(matrices[f], data[f]["density"], data[f]["sigma"],
                                                   labels, combined["coefficients"])
               for f in FIT_FLYBYS}
    return {"description": name, "training_flybys": list(training),
            "combined_chi2": float(combined["chi2"]), "target_fits": targets}


def m13_check(root: Path, flyby: str) -> dict[str, float | None]:
    totals = {}
    for mass in (13, 15):
        with (root / f"{flyby.lower()}_ch4_mass{mass}.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        counts = sum(float(r["c1_counts"]) for r in rows)
        bg = sum(float(r["off_plume_background_counts_per_ip"]) for r in rows)
        totals[mass] = max(counts - bg, 0.0)
    return {"net_mass13_counts": totals[13], "net_mass15_counts": totals[15],
            "observed_net_count_ratio_m13_m15": totals[13] / totals[15] if totals[15] > 0 else None,
            "pure_ch4_calibration_ratio_m13_m15": 0.0205 / 0.723}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--validation-data-root", type=Path, required=True)
    p.add_argument("--source-case", choices=SOURCE_CASES, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--config", type=Path, default=APP / "config/tenishev2014_stage0.json")
    p.add_argument("--include-flagged", action="store_true",
                   help="sensitivity run: include cross-channel-incompatible mass-15 samples")
    args = p.parse_args()
    if args.output_dir.exists(): p.error("output directory already exists")
    args.output_dir.mkdir(parents=True)
    root = args.validation_data_root.resolve(); profiles = root / "data/processed/ch4_fit_v1"
    config = copy.deepcopy(load_config(args.config, args.source_case, validation_data_root=root))
    config["molecular_species"] = "CH4"; config["molecular_mass_kg"] = M_CH4
    _, labels = sources_from_config(config)
    data, matrices, independent = {}, {}, {}
    for flyby in FIT_FLYBYS:
        data[flyby] = read_profile(profiles / f"{flyby.lower()}_ch4_mass15.csv", args.include_flagged)
        matrix, matrix_labels = matrix_from_config(config, data[flyby]["points"])
        if matrix_labels != labels: raise RuntimeError("component label mismatch")
        matrices[flyby] = matrix
        independent[flyby] = solve_amplitudes(matrix, data[flyby]["density"], data[flyby]["sigma"])
    strategies = {
        "independent_each_observation": {"description": "Each flyby fits its own shape and total rate.",
            "training_flybys": None,
            "target_fits": {f: direct_result(independent[f], labels, len(data[f]["density"])) for f in FIT_FLYBYS}},
        "shared_shape_primary_low_impact": combined_shape("Shape fitted jointly to E14/E17/E18", PRIMARY, data, matrices, labels),
        "shared_shape_all_low_speed": combined_shape("Shape fitted jointly to E7/E14/E17/E18/E21", LOW_SPEED, data, matrices, labels),
    }
    rows = []
    for strategy, item in strategies.items():
        for flyby, fit in item["target_fits"].items():
            rate = float(fit["total_source_rate_s"])
            rows.append({"source_case": args.source_case, "strategy": strategy, "flyby": flyby,
                         "interpretation_tier": "primary" if flyby in PRIMARY else ("impact_chemistry_diagnostic" if flyby in ("E3", "E5") else "secondary"),
                         "ch4_source_rate_s-1": rate,
                         "scale_low_23pct_s-1": rate * (1-CALIBRATION_SCALE_FRACTION),
                         "scale_high_23pct_s-1": rate * (1+CALIBRATION_SCALE_FRACTION),
                         "m15_share_0p8_s-1": rate * M15_SHARE_RANGE[0] / M15_SHARE_NOMINAL,
                         "m15_share_1p0_s-1": rate * M15_SHARE_RANGE[1] / M15_SHARE_NOMINAL,
                         "chi2_per_sample": fit["chi2_per_sample"], "samples": fit["samples"]})
    with (args.output_dir / "ch4_source_rates.csv").open("w", newline="", encoding="utf-8") as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    colors = {"independent_each_observation": "#d73027", "shared_shape_primary_low_impact": "#225ea8",
              "shared_shape_all_low_speed": "#238b45"}
    fig, axes = plt.subplots(4, 2, figsize=(12, 13), constrained_layout=True)
    for ax, flyby in zip(axes.flat, FIT_FLYBYS):
        d = data[flyby]
        ax.errorbar(d["time"], d["density"], yerr=d["sigma"], fmt=".", ms=2.5, elinewidth=.3,
                    color="black", label="mass 15 CH4-equivalent (gross)")
        for name, item in strategies.items():
            fit = item["target_fits"][flyby]
            pred = matrices[flyby] @ np.asarray(fit["full_coefficients"])
            ax.plot(d["time"], pred, lw=1.25, color=colors[name], label=name.replace("_", " "))
        ax.set(title=flyby, xlabel="Time from closest approach (s)", ylabel=r"CH$_4$-equivalent density (m$^{-3}$)")
        ax.grid(alpha=.2); ax.legend(fontsize=6)
    axes.flat[-1].axis("off")
    fig.suptitle(args.source_case + ": Cassini INMS CH4 mass-15 fits")
    fig.savefig(args.output_dir / "ch4_fit_profiles.png", dpi=200); fig.savefig(args.output_dir / "ch4_fit_profiles.pdf")
    plt.close(fig)

    checks = {f: m13_check(profiles, f) for f in FIT_FLYBYS + ("E13",)}
    summary = {"study": "Cassini INMS CH4 mass-15 plume source fitting", "source_case": args.source_case,
        "molecular_mass_kg": M_CH4, "include_cross_channel_flagged_samples": args.include_flagged,
        "component_labels": labels, "strategies": strategies,
        "mass13_cross_check": checks,
        "uncertainty_model": {"pointwise": "archived Poisson C1 count error",
            "common_calibration_fraction": CALIBRATION_SCALE_FRACTION,
            "mass15_CH4_share_nominal": M15_SHARE_NOMINAL, "mass15_CH4_share_sensitivity_range": list(M15_SHARE_RANGE)},
        "interpretation_limits": [
            "E14/E17/E18 are the primary result; E3/E5 are not compositional measurements because hypervelocity grain/ice impacts alter the mass spectrum.",
            "The 90% mass-15 CH4 allocation is an accumulated-spectrum result applied to time-resolved samples; 80-100% is reported as a sensitivity bracket, not a confidence interval.",
            "Mass 13 is an aggregate cross-check only; its low counts and excess above the pure-CH4 fragmentation ratio preclude a time-resolved two-channel inversion.",
            "The 23% revised on-ram calibration uncertainty is common-mode and is not included in pointwise error bars.",
            "Individual jet coefficients are non-unique; source totals are conditional on the selected geometry, source kinematics, collisionless propagation, and nuisance-background model.",
        ]}
    (args.output_dir / "study_results.json").write_text(
        json.dumps(summary, indent=2, default=json_default, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir), "rates": rows}, indent=2)); return 0


if __name__ == "__main__": raise SystemExit(main())
