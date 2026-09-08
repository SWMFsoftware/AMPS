#!/usr/bin/env python3
"""Create reproducible publication figures from the study CSV products.

Every panel is derived from a machine-readable table produced elsewhere in the
package.  The plotting layer performs no boundary extraction or scientific
filtering.  Missing products cause the corresponding figure to be skipped with
an explicit message, which makes partial validation runs usable without
inventing placeholder values.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from study_common import default_output_root, load_config, read_driver, resolve_output_path


def save_figure(figure, base: Path) -> None:
    """Save both review-friendly PNG and vector PDF versions."""
    base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(base.with_suffix(".png"), dpi=220, bbox_inches="tight")
    figure.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def driver_figure(root: Path, config, output: Path) -> None:
    rows = read_driver(root / config["data"]["driver"])
    time = [row.epoch for row in rows]
    figure, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(time, [row.pdyn_npa for row in rows], color="tab:purple")
    axes[0].set_ylabel(r"$P_{dyn}$ [nPa]")
    axes[1].plot(time, [row.bz_nt for row in rows], color="tab:blue")
    axes[1].axhline(0.0, color="black", linewidth=0.7)
    axes[1].set_ylabel(r"IMF $B_z$ [nT]")
    axes[2].plot(time, [row.symh_nt for row in rows], color="tab:red")
    axes[2].axhline(0.0, color="black", linewidth=0.7)
    axes[2].set_ylabel("SYM-H [nT]")
    for index in range(1, 7):
        axes[3].plot(time, [getattr(row, f"w{index}") for row in rows],
                     linewidth=1.0, label=rf"$W_{index}$")
    axes[3].set_ylabel(r"TS05 $W_i$")
    axes[3].legend(ncol=6, fontsize=8, loc="upper right")
    axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.suptitle("December 2006 TS05 driver")
    save_figure(figure, output / "figure_driver")


def comparison_figure(path: Path, output: Path) -> None:
    if not path.exists():
        print(f"Skipping comparison figure; missing {path}")
        return
    data = pd.read_csv(path)
    data["epoch_utc"] = pd.to_datetime(data["epoch_utc"], utc=True)
    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    pamela = data[(data.dataset == "PAMELA_TABLE_S1") &
                  data.modeled_boundary_aacgm_deg.notna()].copy()
    for rigidity, group in pamela.groupby("rigidity_gv"):
        axes[0, 0].plot(group.epoch_utc, group.observed_boundary_aacgm_deg,
                        marker="o", linewidth=1, label=f"{rigidity:.3f} GV obs")
        axes[0, 0].plot(group.epoch_utc, group.modeled_boundary_aacgm_deg,
                        linestyle="--", linewidth=1)
    axes[0, 0].set_title("PAMELA: solid/markers observed, dashed AMPS")
    axes[0, 0].set_ylabel("Cutoff AACGM latitude [deg]")
    axes[0, 0].legend(fontsize=6, ncol=2)
    axes[1, 0].scatter(pamela.rigidity_gv,
                       pamela.model_minus_observation_deg, s=12, alpha=0.65)
    axes[1, 0].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 0].set_xlabel("Rigidity [GV]")
    axes[1, 0].set_ylabel("AMPS - observation [deg]")

    poes = data[(data.dataset == "NOAA_POES_METOP_SEM2") &
                data.modeled_boundary_aacgm_deg.notna()].copy()
    primary = poes[poes.used_for_primary_metrics.astype(str).str.lower().isin(("true", "1"))]
    for channel, group in primary.groupby("channel"):
        axes[0, 1].scatter(group.observed_boundary_aacgm_deg,
                           group.modeled_boundary_aacgm_deg, s=12, alpha=0.55,
                           label=channel)
    values = pd.concat([primary.observed_boundary_aacgm_deg,
                        primary.modeled_boundary_aacgm_deg]).dropna()
    if not values.empty:
        limits = (values.min() - 1, values.max() + 1)
        axes[0, 1].plot(limits, limits, color="black", linestyle="--")
        axes[0, 1].set_xlim(limits)
        axes[0, 1].set_ylim(limits)
    axes[0, 1].set_title("POES/MetOp P6/P7 paired boundaries")
    axes[0, 1].set_xlabel("Observed AACGM latitude [deg]")
    axes[0, 1].set_ylabel("AMPS AACGM latitude [deg]")
    axes[0, 1].legend()
    for channel, group in poes.groupby("channel"):
        axes[1, 1].scatter(group.mlt_hour, group.model_minus_observation_deg,
                           s=10, alpha=0.5, label=channel)
    axes[1, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 1].set_xlabel("MLT [h]")
    axes[1, 1].set_ylabel("AMPS - observation [deg]")
    axes[1, 1].legend(ncol=4, fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=0.25)
    figure.autofmt_xdate()
    figure.suptitle("Observation-equivalent AMPS comparisons")
    save_figure(figure, output / "figure_model_data_comparison")


def dynamics_figure(path: Path, output: Path) -> None:
    if not path.exists():
        print(f"Skipping dynamics figure; missing {path}")
        return
    data = pd.read_csv(path)
    data["epoch_utc"] = pd.to_datetime(data["epoch_utc"], utc=True)
    selected = [0.174013525, 0.423556372, 0.692820323, 1.131017241]
    figure, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    subset = data[(np.isclose(data.altitude_km, 850.0)) & (data.hemisphere == "N")]
    for target in selected:
        group = subset[np.isclose(subset.rigidity_gv, target)]
        if group.empty:
            continue
        label = f"{target:.3f} GV"
        axes[0].plot(group.epoch_utc, group.mean_latitude_deg, label=label)
        axes[1].plot(group.epoch_utc, group.amplitude_deg, label=label)
        axes[2].plot(group.epoch_utc,
                     group.accessible_area_fraction_in_analyzed_band, label=label)
    axes[0].set_ylabel(r"$\Lambda_0$ [deg]")
    axes[1].set_ylabel(r"$A_1$ [deg]")
    axes[2].set_ylabel("Accessible fraction")
    axes[2].set_xlabel("UTC")
    axes[0].legend(ncol=4, fontsize=8)
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.suptitle("Rigidity-dependent cutoff morphology at 850 km, north")
    save_figure(figure, output / "figure_cutoff_dynamics")


def lag_hysteresis_figure(lag_path: Path, hysteresis_path: Path, output: Path) -> None:
    if not lag_path.exists() or not hysteresis_path.exists():
        print("Skipping lag/hysteresis figure; analysis products are incomplete")
        return
    lag = pd.read_csv(lag_path)
    hys = pd.read_csv(hysteresis_path)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    target_r = 0.423556372
    selected = lag[(np.isclose(lag.altitude_km, 850.0)) &
                   (np.isclose(lag.rigidity_gv, target_r)) &
                   (lag.hemisphere == "N") &
                   (lag.driver_variable.isin(["pdyn_npa", "bz_nt", "symh_nt", "w1"]))]
    for variable, group in selected.groupby("driver_variable"):
        axes[0].plot(group.lag_minutes / 60.0, group.correlation, label=variable)
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set_xlabel("Lag [h]; positive means cutoff follows driver")
    axes[0].set_ylabel("Correlation")
    axes[0].set_title(f"Lag response at {target_r:.3f} GV")
    axes[0].legend()

    strict = hys[hys.match_definition == "STRICT"]
    for hemisphere, group in strict.groupby("hemisphere"):
        axes[1].plot(group.rigidity_gv, group.median_recovery_minus_main_deg,
                     marker="o", label=hemisphere)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_xlabel("Rigidity [GV]")
    axes[1].set_ylabel("Recovery - main boundary [deg]")
    axes[1].set_title("Matched-driver hysteresis")
    axes[1].legend()
    for axis in axes:
        axis.grid(alpha=0.25)
    save_figure(figure, output / "figure_lag_hysteresis")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument(
        "--comparison-root", type=Path,
        default=default_output_root() / "comparison",
    )
    parser.add_argument(
        "--dynamics-root", type=Path,
        default=default_output_root() / "dynamics",
    )
    parser.add_argument(
        "--output-root", type=Path, default=default_output_root() / "figures"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root, config = load_config(args.config)
    comparison, dynamics, output = map(resolve_output_path, (
        args.comparison_root, args.dynamics_root, args.output_root
    ))
    driver_figure(root, config, output)
    comparison_figure(comparison / "paired_model_observation.csv", output)
    dynamics_figure(dynamics / "cutoff_dynamics_timeseries.csv", output)
    lag_hysteresis_figure(dynamics / "lag_correlations.csv",
                          dynamics / "hysteresis_summary.csv", output)
    print(f"Figure products: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
