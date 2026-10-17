#!/usr/bin/env python3
"""Fit Cassini INMS CO2 profiles and transfer E3/E5 H2O/CO2 ratios."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from fitting import (SOURCE_CASES, append_extra_jet, load_config,  # noqa: E402
                     matrix_from_config,
                     solve_amplitudes, sources_from_config)
from joint_fitting_study import (fit_rate_for_fixed_source_shape,  # noqa: E402
                                 json_default, normalized_source_shape,
                                 solve_combined_blocks)
from validation_paths import (CO2_FIT_RELATIVE,  # noqa: E402
                              add_validation_root_argument,
                              resolve_input_path, resolve_validation_root)


CO2_MASS_KG = 44.0095 * 1.66053906660e-27
FLYBYS = ("E3", "E5", "E14", "E17", "E18")
COMMON_CALIBRATION_FRACTION = 0.15
H2O_SCALE_FRACTION = 0.23
FRAGMENTATION_RANGE = (2.0, 5.0)


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty table")
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_profile(path: Path) -> dict[str, object]:
    rows = read_csv(path)
    points = np.asarray([[float(row[name]) for name in
                          ("spacecraft_x_m", "spacecraft_y_m", "spacecraft_z_m")]
                         for row in rows])
    density = np.asarray([float(row["density_m3"]) for row in rows])
    sigma = np.asarray([float(row["sigma_counting_m3"]) for row in rows])
    time = np.asarray([float(row["time_ca_s"]) for row in rows])
    if (points.shape != (len(rows), 3) or np.any(density < 0)
            or np.any(sigma <= 0) or not all(np.isfinite(value).all()
                                             for value in (points, density, sigma, time))):
        raise ValueError(f"Invalid CO2 profile: {path}")
    return {"rows": rows, "points": points, "density": density,
            "sigma": sigma, "time": time}


def source_total(vector: np.ndarray, labels: list[str]) -> float:
    return float(sum(value for label, value in zip(labels, vector)
                     if label.startswith(("jet:", "stripe:"))) * 1.0e28)


def direct_fit_result(fit: dict[str, object], labels: list[str], samples: int) -> dict[str, object]:
    vector = np.asarray(fit["coefficients"])
    shape = normalized_source_shape(vector, labels)
    return {
        "total_source_rate_s": source_total(vector, labels),
        "source_shape_fractions": {
            label: float(shape[index]) for index, label in enumerate(labels)
            if label.startswith(("jet:", "stripe:"))
        },
        "full_coefficients": vector,
        "chi2": fit["chi2"],
        "chi2_per_sample": fit["chi2"] / samples,
        "rank": fit["rank"],
        "samples": samples,
        "condition_number_scaled": fit["condition_number_scaled"],
    }


def pair_h2o_co2(h2o_result: dict[str, object], co2_strategies: dict[str, object]) -> dict[str, object]:
    h2o_strategies = h2o_result["source_rate_shape_strategies"]
    pairings = {
        "independent_each_species": (
            "independent_each_observation", "independent_each_observation"),
        "shared_e3_e5_shape_each_species": (
            "shared_shape_inms_only", "shared_shape_e3_e5"),
        "h2o_e3e5_shape_co2_all_epoch_shape": (
            "shared_shape_inms_only", "shared_shape_all_co2"),
    }
    result = {}
    for name, (h2o_strategy, co2_strategy) in pairings.items():
        h2o = np.asarray([
            float(h2o_strategies[h2o_strategy]["target_fits"][event]["total_source_rate_s"])
            for event in ("E3", "E5")])
        co2 = np.asarray([
            float(co2_strategies[co2_strategy]["target_fits"][event]["total_source_rate_s"])
            for event in ("E3", "E5")])
        if np.any(co2 <= 0) or np.any(h2o <= 0):
            raise RuntimeError("Positive E3/E5 source rates are required for composition transfer")
        apparent_ratio = h2o / co2
        proportional_ratio = float(co2 @ h2o / (co2 @ co2))
        predicted = proportional_ratio * co2
        result[name] = {
            "h2o_strategy": h2o_strategy,
            "co2_strategy": co2_strategy,
            "E3_apparent_h2o_to_co2_rate_ratio": float(apparent_ratio[0]),
            "E5_apparent_h2o_to_co2_rate_ratio": float(apparent_ratio[1]),
            "proportional_through_origin_apparent_ratio": proportional_ratio,
            "two_point_proportional_residual_fraction": {
                "E3": float((predicted[0] - h2o[0]) / h2o[0]),
                "E5": float((predicted[1] - h2o[1]) / h2o[1]),
            },
            "actual_h2o_to_co2_ratio_fragmentation_bracket": [
                FRAGMENTATION_RANGE[0] * proportional_ratio,
                FRAGMENTATION_RANGE[1] * proportional_ratio,
            ],
            "nominal_log_midpoint_ratio_for_plotting_only": (
                math.sqrt(FRAGMENTATION_RANGE[0] * FRAGMENTATION_RANGE[1])
                * proportional_ratio),
            "scale_fractional_sigma_excluding_fragmentation": math.hypot(
                H2O_SCALE_FRACTION, COMMON_CALIBRATION_FRACTION),
            "correlation_status": (
                "Only two paired epochs exist. Pearson/Spearman correlation is "
                "not reported because n=2 cannot establish correlation strength."
            ),
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-case", choices=SOURCE_CASES, required=True)
    parser.add_argument("--jet-catalog", type=Path)
    parser.add_argument(
        "--extra-jet", type=Path,
        help="append one provenance-tracked development jet after resolving the base case")
    add_validation_root_argument(parser)
    parser.add_argument(
        "--co2-data-root", type=Path,
        help="co2_fit_v1 directory; overrides the conventional validation-root location")
    parser.add_argument("--h2o-study-result", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        validation_root = resolve_validation_root(args.validation_data_root)
        args.co2_data_root = resolve_input_path(
            args.co2_data_root, validation_root, CO2_FIT_RELATIVE,
            "--co2-data-root")
    except ValueError as exc:
        parser.error(str(exc))
    if args.output_dir.exists():
        parser.error("output directory already exists")
    args.output_dir.mkdir(parents=True)

    config = load_config(args.config, args.source_case, args.jet_catalog,
                         validation_root)
    config = copy.deepcopy(config)
    if args.extra_jet:
        append_extra_jet(config, args.extra_jet)
    config["molecular_species"] = "CO2"
    config["molecular_mass_kg"] = CO2_MASS_KG
    sources, labels = sources_from_config(config)
    del sources

    data, matrices, independent = {}, {}, {}
    for flyby in FLYBYS:
        path = args.co2_data_root / f"{flyby.lower()}_co2_density.csv"
        profile = load_profile(path)
        matrix, matrix_labels = matrix_from_config(config, profile["points"])
        if matrix_labels != labels:
            raise AssertionError("Component-label mismatch")
        fit = solve_amplitudes(matrix, profile["density"], profile["sigma"])
        data[flyby], matrices[flyby], independent[flyby] = profile, matrix, fit

    rate_strategies: dict[str, object] = {
        "independent_each_observation": {
            "description": "Each CO2 flyby determines its own source shape and rate.",
            "shape_training_datasets": "target flyby only",
            "target_fits": {
                flyby: direct_fit_result(independent[flyby], labels,
                                         len(data[flyby]["density"]))
                for flyby in FLYBYS},
        }
    }
    combined_definitions = {
        "shared_shape_e3_e5": ("E3 and E5", ("E3", "E5")),
        "shared_shape_all_co2": ("E3, E5, E14, E17, and E18", FLYBYS),
    }
    combined_fits = {}
    for strategy, (description, training) in combined_definitions.items():
        blocks = [(flyby, matrices[flyby], data[flyby]["density"],
                   data[flyby]["sigma"]) for flyby in training]
        combined = solve_combined_blocks(blocks, balance_by_dataset=True)
        combined_fits[strategy] = combined
        rate_strategies[strategy] = {
            "description": (
                f"One relative CO2 source shape is trained on {description}; "
                "each flyby then refits one total source rate and backgrounds."),
            "shape_training_datasets": description,
            "target_fits": {
                flyby: fit_rate_for_fixed_source_shape(
                    matrices[flyby], data[flyby]["density"], data[flyby]["sigma"],
                    labels, combined["coefficients"])
                for flyby in FLYBYS},
        }

    h2o_result = json.loads(args.h2o_study_result.read_text(encoding="utf-8"))
    source_variant = config.get("source_variant", args.source_case)
    if h2o_result.get("source_case") != args.source_case:
        raise ValueError("H2O study source case does not match CO2 source case")
    if h2o_result.get("source_variant", h2o_result.get("source_case")) != source_variant:
        raise ValueError("H2O study source variant does not match CO2 source variant")
    ratios = pair_h2o_co2(h2o_result, rate_strategies)

    rate_rows = []
    proxy_rows = []
    for strategy, result in rate_strategies.items():
        for flyby, fit in result["target_fits"].items():
            rate = float(fit["total_source_rate_s"])
            apparent = flyby in ("E3", "E5")
            rate_rows.append({
                "source_case": source_variant,
                "strategy": strategy,
                "shape_training_datasets": result["shape_training_datasets"],
                "flyby": flyby,
                "co2_source_rate_s-1": rate,
                "co2_interpretation": "apparent_CO2" if apparent else "CO2",
                "actual_co2_rate_lower_s-1": rate / FRAGMENTATION_RANGE[1] if apparent else rate,
                "actual_co2_rate_upper_s-1": rate / FRAGMENTATION_RANGE[0] if apparent else rate,
                "calibration_fractional_sigma_common": COMMON_CALIBRATION_FRACTION,
                "chi2_per_sample": fit["chi2_per_sample"],
                "samples": fit["samples"],
            })

    ratio_rows = []
    for pairing, item in ratios.items():
        ratio_rows.append({
            "source_case": source_variant,
            "pairing": pairing,
            "h2o_strategy": item["h2o_strategy"],
            "co2_strategy": item["co2_strategy"],
            "E3_apparent_h2o_to_co2": item["E3_apparent_h2o_to_co2_rate_ratio"],
            "E5_apparent_h2o_to_co2": item["E5_apparent_h2o_to_co2_rate_ratio"],
            "apparent_ratio_through_origin": item["proportional_through_origin_apparent_ratio"],
            "actual_ratio_lower_fragmentation_2x": item["actual_h2o_to_co2_ratio_fragmentation_bracket"][0],
            "actual_ratio_upper_fragmentation_5x": item["actual_h2o_to_co2_ratio_fragmentation_bracket"][1],
            "nominal_log_midpoint_plot_only": item["nominal_log_midpoint_ratio_for_plotting_only"],
        })
        co2_strategy = item["co2_strategy"]
        for flyby in ("E14", "E17", "E18"):
            co2_rate = float(rate_strategies[co2_strategy]["target_fits"][flyby]["total_source_rate_s"])
            low_ratio, high_ratio = item["actual_h2o_to_co2_ratio_fragmentation_bracket"]
            proxy_rows.append({
                "source_case": source_variant,
                "pairing": pairing,
                "co2_strategy": co2_strategy,
                "flyby": flyby,
                "co2_source_rate_s-1": co2_rate,
                "h2o_proxy_lower_s-1": low_ratio * co2_rate,
                "h2o_proxy_nominal_log_midpoint_s-1": item["nominal_log_midpoint_ratio_for_plotting_only"] * co2_rate,
                "h2o_proxy_upper_s-1": high_ratio * co2_rate,
                "proxy_limit": "2--5 E3/E5 fragmentation factor dominates; nominal is plotting convention",
            })
    write_csv(args.output_dir / "co2_source_rates.csv", rate_rows)
    write_csv(args.output_dir / "h2o_co2_source_rate_ratios.csv", ratio_rows)
    write_csv(args.output_dir / "h2o_proxy_source_rates.csv", proxy_rows)

    fig, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    axes = axes.flat
    colors = {"independent_each_observation": "#d73027",
              "shared_shape_e3_e5": "#756bb1",
              "shared_shape_all_co2": "#31a354"}
    for axis, flyby in zip(axes, FLYBYS):
        profile = data[flyby]
        axis.errorbar(profile["time"], profile["density"], yerr=profile["sigma"],
                      fmt=".", ms=2.5, elinewidth=.35, color="black",
                      label="INMS mass-44 density")
        for strategy in rate_strategies:
            fit = rate_strategies[strategy]["target_fits"][flyby]
            prediction = matrices[flyby] @ np.asarray(fit["full_coefficients"])
            axis.semilogy(profile["time"], np.maximum(prediction, 1.0),
                          color=colors[strategy], lw=1.4,
                          label=strategy.replace("_", " "))
        axis.set(title=flyby, xlabel="Time from closest approach (s)",
                 ylabel=r"CO$_2$ number density (m$^{-3}$)")
        axis.grid(True, which="both", alpha=.2)
        axis.legend(fontsize=7)
    axes[-1].axis("off")
    fig.suptitle(source_variant + ": CO2 profile fits")
    fig.savefig(args.output_dir / "co2_fit_profiles.png", dpi=190)
    fig.savefig(args.output_dir / "co2_fit_profiles.pdf")
    plt.close(fig)

    summary = {
        "study": "E3/E5/E14/E17/E18 CO2 source fitting and H2O composition transfer",
        "source_case": args.source_case,
        "source_variant": source_variant,
        "extra_jet": config.get("_extra_jet_path"),
        "extra_jet_sha256": (
            sha256(config["_extra_jet_path"])
            if config.get("_extra_jet_path") else None),
        "config": str(args.config.resolve()),
        "config_sha256": sha256(args.config),
        "co2_molecular_mass_kg": CO2_MASS_KG,
        "kinematic_assumption": (
            "CO2 uses the same source temperature and bulk velocity as H2O; "
            "the molecular mass changes the thermal speed and plume collimation."
        ),
        "data_root": str(args.co2_data_root.resolve()),
        "h2o_study_result": str(args.h2o_study_result.resolve()),
        "component_labels": labels,
        "independent_fits": independent,
        "combined_shape_fits": combined_fits,
        "source_rate_strategies": rate_strategies,
        "h2o_co2_source_rate_ratio": ratios,
        "interpretation_limits": [
            "E3/E5 are apparent CO2 because high-speed impact fragmentation can enhance mass 44 by about 2--5.",
            "The 15% INMS calibration uncertainty is common-scale and is not added independently to every point.",
            "Only two epochs have paired H2O and CO2 fits; a correlation coefficient would be uninformative.",
            "The H2O transfer proxy assumes an orbit-invariant intrinsic H2O/CO2 source ratio; its 2--5 bracket is a sensitivity range, not a confidence interval.",
            "Individual jet amplitudes remain ill-conditioned; totals are conditional on source geometry and collisionless physics.",
        ],
    }
    with (args.output_dir / "study_results.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, default=json_default, allow_nan=False)
        stream.write("\n")

    report = [
        "# CO2 fitting and H2O/CO2 transfer study", "",
        f"Source geometry: `{source_variant}`. CO2 is propagated with mass "
        f"{CO2_MASS_KG:.8e} kg at the same source temperature and bulk speed as H2O.", "",
        "## Fitted CO2 source rates", "",
        "Rates are molecules s^-1. E3/E5 values are apparent CO2; their actual "
        "CO2 brackets are the fitted rate divided by 5 and 2.", "",
        "| Strategy | Flyby | Fitted rate | Interpretation | chi-square/sample |",
        "|---|---|---:|---|---:|",
    ]
    for row in rate_rows:
        report.append(
            f"| {row['strategy']} | {row['flyby']} | {row['co2_source_rate_s-1']:.6e} | "
            f"{row['co2_interpretation']} | {row['chi2_per_sample']:.3f} |")
    report += ["", "## Paired H2O/CO2 source-rate ratios", "",
               "The proportional calibration is constrained through the origin. "
               "With only E3 and E5, it is not a correlation measurement.", "",
               "| Pairing | E3 apparent ratio | E5 apparent ratio | Combined apparent ratio | Actual ratio bracket |",
               "|---|---:|---:|---:|---:|"]
    for row in ratio_rows:
        report.append(
            f"| {row['pairing']} | {row['E3_apparent_h2o_to_co2']:.4g} | "
            f"{row['E5_apparent_h2o_to_co2']:.4g} | "
            f"{row['apparent_ratio_through_origin']:.4g} | "
            f"{row['actual_ratio_lower_fragmentation_2x']:.4g}--"
            f"{row['actual_ratio_upper_fragmentation_5x']:.4g} |")
    report += ["", "## Interpretation", ""]
    report.extend(f"- {item}" for item in summary["interpretation_limits"])
    report += ["", "See `co2_fit_profiles.png`, `co2_source_rates.csv`, "
               "`h2o_co2_source_rate_ratios.csv`, `h2o_proxy_source_rates.csv`, "
               "and `study_results.json`.", ""]
    (args.output_dir / "fitting_study.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir), "source_case": source_variant,
                      "ratios": ratios}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
