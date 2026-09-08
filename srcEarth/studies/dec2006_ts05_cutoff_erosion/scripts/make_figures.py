#!/usr/bin/env python3
"""Create reproducible publication figures from the study CSV products.

Every panel is derived from a machine-readable table produced elsewhere in the
package.  The plotting layer performs no boundary extraction or scientific
filtering.  Every completed figure is written as a high-resolution PNG, an EPS
vector graphic requested by common journal workflows, and a PDF vector copy.
Missing optional validation products are reported explicitly.  The dedicated
cutoff-degradation products are required when the top-level publication runner
uses ``--require-publication-products``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from study_common import default_output_root, load_config, read_driver, resolve_output_path


def save_figure(figure, base: Path) -> list[Path]:
    """Save one figure in review, journal-vector, and archival formats.

    EPS is intentionally generated directly from Matplotlib rather than by
    converting a raster.  The resulting lines and text remain vector objects
    for publication.  PDF is retained because it is convenient for internal
    review and usually preserves transparency better than PostScript.
    """

    base.parent.mkdir(parents=True, exist_ok=True)
    products = [base.with_suffix(suffix) for suffix in (".png", ".eps", ".pdf")]
    figure.savefig(products[0], dpi=300, bbox_inches="tight")
    figure.savefig(products[1], format="eps", bbox_inches="tight")
    figure.savefig(products[2], format="pdf", bbox_inches="tight")
    plt.close(figure)
    for product in products:
        print(f"Wrote figure: {product}", flush=True)
    return products


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


def _cell_edges(centers: np.ndarray, fallback_half_width: float) -> np.ndarray:
    """Convert ordered cell centers to edges for an undistorted pcolormesh."""

    centers = np.asarray(centers, dtype=float)
    if centers.size == 1:
        return np.array([centers[0] - fallback_half_width,
                         centers[0] + fallback_half_width])
    midpoints = 0.5 * (centers[1:] + centers[:-1])
    return np.concatenate((
        [centers[0] - (midpoints[0] - centers[0])],
        midpoints,
        [centers[-1] + (centers[-1] - midpoints[-1])],
    ))


def cutoff_degradation_figures(path: Path, output: Path) -> list[Path]:
    """Plot the central science result: storm-time loss of cutoff shielding.

    ``cutoff_erosion_deg`` is the fitted mean cutoff latitude minus its quiet
    reference.  Negative values therefore denote an equatorward displacement
    and reduced geomagnetic shielding.  The heat maps average north and south
    only after each hemisphere has been reduced independently; the companion
    line plot retains the hemispheres and reports the largest equatorward
    displacement at each rigidity.  Both products are calculated exclusively
    from the archived dynamics table and introduce no new filtering.
    """

    if not path.exists():
        print(f"Skipping cutoff-degradation figures; missing {path}", flush=True)
        return []
    data = pd.read_csv(path)
    required = {
        "epoch_utc", "altitude_km", "rigidity_gv", "hemisphere",
        "cutoff_erosion_deg",
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        print(f"Skipping cutoff-degradation figures; missing columns: {missing}",
              flush=True)
        return []
    data["epoch_utc"] = pd.to_datetime(data["epoch_utc"], utc=True)
    data["cutoff_erosion_deg"] = pd.to_numeric(
        data["cutoff_erosion_deg"], errors="coerce"
    )
    data = data.dropna(subset=["epoch_utc", "altitude_km", "rigidity_gv",
                               "cutoff_erosion_deg"])
    if data.empty:
        print("Skipping cutoff-degradation figures; no finite erosion rows",
              flush=True)
        return []

    # Use one symmetric color scale for both altitudes so visual differences
    # cannot be caused by separate automatic normalization.
    limit = float(np.nanmax(np.abs(data["cutoff_erosion_deg"].to_numpy())))
    limit = max(limit, 0.1)
    altitudes = sorted(data["altitude_km"].unique())
    figure, axes = plt.subplots(len(altitudes), 1,
                               figsize=(10.5, 3.5 * len(altitudes)),
                               sharex=True, squeeze=False)
    mesh = None
    for axis, altitude in zip(axes[:, 0], altitudes):
        subset = data[np.isclose(data.altitude_km, altitude)]
        # A hemispheric mean is suitable for the global erosion overview.  The
        # maximum-degradation figure below keeps N/S behavior separate.
        reduced = subset.groupby(
            ["rigidity_gv", "epoch_utc"], as_index=False
        )["cutoff_erosion_deg"].mean()
        pivot = reduced.pivot(index="rigidity_gv", columns="epoch_utc",
                              values="cutoff_erosion_deg").sort_index()
        pivot = pivot.reindex(sorted(pivot.columns), axis=1)
        x_centers = mdates.date2num(pivot.columns.to_pydatetime())
        y_centers = pivot.index.to_numpy(dtype=float)
        x_edges = _cell_edges(x_centers, 1.0 / 48.0)
        y_edges = _cell_edges(y_centers, 0.025)
        mesh = axis.pcolormesh(
            x_edges, y_edges, pivot.to_numpy(dtype=float), shading="flat",
            cmap="RdBu", norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
        )
        axis.set_ylabel("Rigidity [GV]")
        axis.set_title(f"{altitude:g} km; mean of independently fitted N/S boundaries")
        axis.grid(False)
    axes[-1, 0].set_xlabel("UTC")
    axes[-1, 0].xaxis_date()
    axes[-1, 0].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    assert mesh is not None
    colorbar = figure.colorbar(mesh, ax=axes[:, 0].tolist(), pad=0.02)
    colorbar.set_label(r"$\Delta\Lambda_c$ [deg]; negative = cutoff erosion")
    figure.suptitle("December 2006 storm-time cutoff degradation", y=0.995)
    figure.subplots_adjust(right=0.88, hspace=0.28)
    products = save_figure(figure, output / "figure_cutoff_degradation")

    # Publication summary of the peak quiet-relative equatorward motion.  Zero
    # is used when a series never moved equatorward, so the ordinate is a
    # non-negative degradation magnitude rather than a signed displacement.
    peak = data.assign(
        degradation_magnitude_deg=np.maximum(
            0.0, -data["cutoff_erosion_deg"].to_numpy(dtype=float)
        )
    ).groupby(
        ["altitude_km", "hemisphere", "rigidity_gv"], as_index=False
    )["degradation_magnitude_deg"].max()
    figure, axis = plt.subplots(figsize=(8.2, 5.0))
    for (altitude, hemisphere), group in peak.groupby(
            ["altitude_km", "hemisphere"]):
        group = group.sort_values("rigidity_gv")
        linestyle = "-" if str(hemisphere).upper() == "N" else "--"
        axis.plot(group.rigidity_gv, group.degradation_magnitude_deg,
                  marker="o", markersize=3, linewidth=1.3,
                  linestyle=linestyle,
                  label=f"{altitude:g} km, {hemisphere}")
    axis.set_xlabel("Rigidity [GV]")
    axis.set_ylabel(r"Maximum equatorward $-\Delta\Lambda_c$ [deg]")
    axis.set_title("Peak storm-time cutoff degradation relative to quiet reference")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2)
    products += save_figure(figure, output / "figure_peak_cutoff_degradation")
    return products


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
    parser.add_argument(
        "--require-publication-products", action="store_true",
        help="Fail unless both dedicated cutoff-degradation figures are written",
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
    publication_products = cutoff_degradation_figures(
        dynamics / "cutoff_dynamics_timeseries.csv", output
    )
    lag_hysteresis_figure(dynamics / "lag_correlations.csv",
                          dynamics / "hysteresis_summary.csv", output)
    expected = {
        output / f"{stem}{suffix}"
        for stem in ("figure_cutoff_degradation", "figure_peak_cutoff_degradation")
        for suffix in (".png", ".eps")
    }
    missing = sorted(str(path) for path in expected if not path.is_file())
    manifest = {
        "output_root": str(output),
        "publication_products": [str(path) for path in publication_products],
        "required_png_eps": sorted(str(path) for path in expected),
        "missing_required_products": missing,
        "passed": not missing,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Figure products: {output}", flush=True)
    if args.require_publication_products and missing:
        for path in missing:
            print(f"ERROR: required publication figure is missing: {path}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
