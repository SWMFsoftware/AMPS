#!/usr/bin/env python3
"""Paper-style E3/E5 comparison using the current Stage-0 fits."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fitting import load_config, matrix_from_config, solve_amplitudes  # noqa: E402


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_observations(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "time": np.asarray([float(row["time_ca_s"]) for row in rows]),
        "points": np.asarray([[float(row[name]) for name in ("x_m", "y_m", "z_m")]
                              for row in rows]),
        "density": np.asarray([float(row["density_m3"]) for row in rows]),
        "sigma": np.asarray([float(row["sigma_m3"]) for row in rows]),
    }


def calculate_curves(config_path, observations_path, result_path):
    observation = load_observations(observations_path)
    with open(result_path) as stream:
        full_result = json.load(stream)
    config = load_config(config_path)
    matrix, labels = matrix_from_config(config, observation["points"])
    if labels != full_result["component_labels"]:
        raise ValueError("fit result component order does not match current configuration")
    coefficients = np.asarray(full_result["coefficients"], dtype=float)
    full_prediction = matrix @ coefficients
    stored = np.asarray(full_result["predicted"], dtype=float)
    np.testing.assert_allclose(full_prediction, stored, rtol=2e-13, atol=1e-6)

    keep = np.asarray([not label.startswith("stripe:") for label in labels])
    no_stripe = solve_amplitudes(matrix[:, keep], observation["density"],
                                 observation["sigma"])
    no_stripe_dof = len(observation["density"]) - no_stripe["rank"]
    summary = {
        "full_model": {
            "chi2": float(full_result["chi2"]),
            "nominal_dof": int(full_result["nominal_dof"]),
            "nominal_reduced_chi2": float(full_result["nominal_reduced_chi2"]),
        },
        "no_stripe_refit": {
            "chi2": float(no_stripe["chi2"]),
            "nominal_dof": int(no_stripe_dof),
            "nominal_reduced_chi2": float(no_stripe["chi2"] / no_stripe_dof),
            "retained_components": [label for label, selected in zip(labels, keep) if selected],
        },
    }
    return observation, full_prediction, no_stripe["predicted"], summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path,
                        default=ROOT / "config/tenishev2014_stage0.json")
    parser.add_argument(
        "--data-root", type=Path, required=True,
        help="external plume_fit_v1 release directory")
    parser.add_argument("--output-stem", type=Path,
                        default=ROOT / "figures/figure3_current_stage0_fit")
    args = parser.parse_args()
    outputs = [args.output_stem.with_suffix(suffix)
               for suffix in (".png", ".pdf", ".metadata.json")]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        parser.error("refusing to overwrite: " + ", ".join(existing))

    inputs = {
        "E3": (args.data_root / "e3_h2o_sticking_corrected.csv",
               args.data_root / "e3_h2o_sticking_corrected.metadata.json",
               ROOT / "results/e3_fit.json"),
        "E5": (args.data_root / "e5_h2o_sticking_corrected.csv",
               args.data_root / "e5_h2o_sticking_corrected.metadata.json",
               ROOT / "results/e5_fit.json"),
    }
    products = {}
    provenance = {}
    for flyby, (observations, metadata, result) in inputs.items():
        products[flyby] = calculate_curves(args.config, observations, result)
        provenance[flyby] = {
            "observations": str(observations),
            "observations_sha256": sha256(observations),
            "observation_metadata": str(metadata),
            "observation_metadata_sha256": sha256(metadata),
            "fit_result": str(result),
            "fit_result_sha256": sha256(result),
            "statistics": products[flyby][3],
        }

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), sharex=True)
    for ax, flyby in zip(axes, ("E3", "E5")):
        observation, full, no_stripe, _ = products[flyby]
        ax.errorbar(observation["time"], observation["density"],
                    yerr=observation["sigma"], color="black", lw=0.8,
                    marker=".", ms=3.0, elinewidth=0.45, capsize=0,
                    label="Teolis sticking-corrected INMS")
        ax.semilogy(observation["time"], no_stripe, color="#1769e0", lw=1.8,
                    label="Current model refit without stripes")
        ax.semilogy(observation["time"], full, color="#e31a1c", lw=2.2,
                    label="Current full Stage-0 fit")
        ax.set_xlim(0, 410)
        positive = np.concatenate((observation["density"], full, no_stripe))
        ax.set_ylim(0.55 * positive.min(), 1.7 * positive.max())
        ax.set_title(flyby, fontweight="bold")
        ax.set_xlabel("Time from flyby reference (s)")
        ax.grid(True, which="both", color="0.88", lw=0.55)
    axes[0].set_ylabel(r"H$_2$O number density (m$^{-3}$)")
    handles, labels = axes[1].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    desired = ["Teolis sticking-corrected INMS",
               "Current full Stage-0 fit",
               "Current model refit without stripes"]
    fig.legend([by_label[label] for label in desired], desired,
               loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.89))
    fig.suptitle("Current Stage-0 multi-plume fits", y=0.98, fontsize=13)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.16, top=0.74,
                        wspace=0.13)
    fig.text(0.5, 0.035,
             "Fit domain is the outbound (+10 to +400 s) published sticking-corrected density; not raw PDS L1A counts.",
             ha="center", fontsize=8)
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    metadata = {
        "schema_version": 1,
        "figure": "paper-style comparison from the current Stage-0 fitting",
        "red_curve": "current 8-jet + 7-stripe + empirical-background fit",
        "blue_curve": "independent nonnegative refit with all stripe columns removed",
        "black_curve": "Teolis et al. (2010) sticking-corrected density digitized from the published trace",
        "important_boundary": "The PDS L1A V2 archive contains counts, not calibrated H2O density; no raw-count conversion is used in this figure.",
        "config": str(args.config),
        "config_sha256": sha256(args.config),
        "flybys": provenance,
    }
    with open(outputs[2], "x") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
