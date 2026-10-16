#!/usr/bin/env python3
"""Render deterministic MMS01 kernel-reconstruction/reference figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mms01_common import atomic_json, read_csv, sha256_file, utc_now


def _save(fig: plt.Figure, figures: Path, stem: str,
          formats: tuple[str, ...]) -> list[Path]:
    """Save every requested format, removing EPS-incompatible transparency."""
    paths: list[Path] = []
    for extension in formats:
        path = figures / f"{stem}.{extension}"
        changed = [(artist, artist.get_alpha()) for artist in fig.findobj()
                   if extension == "eps" and hasattr(artist, "get_alpha")
                   and artist.get_alpha() is not None and 0 < artist.get_alpha() < 1]
        try:
            for artist, _ in changed:
                artist.set_alpha(1)
            fig.savefig(path, dpi=180 if extension == "png" else None, bbox_inches="tight")
        finally:
            for artist, alpha in changed:
                artist.set_alpha(alpha)
        paths.append(path)
    plt.close(fig)
    return paths


def plot(output_root: Path, formats: tuple[str, ...]) -> dict[str, Any]:
    """Plot saved comparisons only; never recalculate candidate or reference values."""
    if not formats or any(item not in ("png", "eps", "pdf", "svg") for item in formats):
        raise ValueError("formats must be selected from png, eps, pdf, svg")
    summary = json.loads((output_root / "comparison_summary.json").read_text(encoding="utf-8"))
    if summary.get("schema") != "amps-mms01-comparison-v2":
        raise ValueError("unrecognized MMS01 comparison schema; rerun COMPARE before PLOT")
    spectrum = read_csv(output_root / "comparison" / "spectrum_comparison.csv")
    directional = read_csv(output_root / "comparison" / "directional_comparison.csv")
    positions = read_csv(output_root / "reference" / "comparison_positions.csv")
    figures = output_root / "figures"
    figures.mkdir(exist_ok=True)
    generated: list[Path] = []

    distributions = sorted({int(row["distribution_index"]) for row in spectrum})
    selected = [distributions[0], distributions[len(distributions) // 2], distributions[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7), sharex=True, sharey=True)
    for axis, distribution in zip(axes, selected):
        rows = [row for row in spectrum if int(row["distribution_index"]) == distribution]
        energy = np.asarray([float(row["energy_eV"]) for row in rows])
        reference = np.asarray([float(row["reference_angular_integral"]) for row in rows])
        reconstructed = np.asarray([
            float(row["kernel_reconstructed_angular_integral"]) for row in rows
        ])
        axis.loglog(energy, reference, "o-", ms=2.5, lw=1.0, label="HPCA L2 flux integral")
        axis.loglog(energy, reconstructed, "--", lw=1.0,
                    label="Production-kernel reconstruction from HPCA PSD")
        zero_bins = int(np.count_nonzero(reference == 0.0))
        if zero_bins:
            # Log axes cannot display physical zero.  State the count instead of
            # adding a plotting floor that could be mistaken for measured flux.
            axis.text(0.03, 0.04, f"exact zero bins: {zero_bins}",
                      transform=axis.transAxes, fontsize=8)
        axis.set_title(rows[0]["utc_center"][11:23] + " UTC")
        axis.grid(True, which="both", alpha=0.25)
    # A figure-level legend keeps the long, scientifically explicit reconstruction
    # label from hiding the first spectrum. Reserve vertical space so it also stays
    # separate from the panel timestamps and overall title in PNG and EPS.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=2, fontsize=8)
    axes[1].set_xlabel("H+ energy [eV]")
    axes[0].set_ylabel("Angular integral [cm$^{-2}$ s$^{-1}$ eV$^{-1}$]")
    fig.suptitle("MMS01 HPCA product reconstruction: first/middle/last spectra "
                 "(selection independent of residual)")
    fig.subplots_adjust(top=0.80, bottom=0.18)
    generated += _save(fig, figures, "mms01_representative_spectra", formats)

    n_time = len(distributions)
    reference_map = np.empty((63, n_time))
    reconstructed_map = np.empty_like(reference_map)
    labels = []
    energy_axis = None
    for column, distribution in enumerate(distributions):
        rows = [row for row in spectrum if int(row["distribution_index"]) == distribution]
        rows.sort(key=lambda row: int(row["energy_index"]))
        reference_map[:, column] = [float(row["reference_angular_integral"]) for row in rows]
        reconstructed_map[:, column] = [
            float(row["kernel_reconstructed_angular_integral"]) for row in rows
        ]
        if energy_axis is None:
            energy_axis = np.asarray([float(row["energy_eV"]) for row in rows])
        labels.append(rows[0]["utc_center"][11:19])
    # Build logarithmic views only on their mathematical domains.  A cell where
    # both products are exactly zero has zero comparison error but no log10 flux;
    # it remains blank in the flux maps and is represented at the documented
    # 1e-16 *display-only* floor in the error map.
    positive_reference = reference_map > 0.0
    positive_reconstructed = reconstructed_map > 0.0
    reference_log = np.full_like(reference_map, np.nan)
    reconstructed_log = np.full_like(reconstructed_map, np.nan)
    reference_log[positive_reference] = np.log10(reference_map[positive_reference])
    reconstructed_log[positive_reconstructed] = np.log10(
        reconstructed_map[positive_reconstructed])
    residual = np.full_like(reference_map, np.nan)
    residual[positive_reference] = np.abs(
        reconstructed_map[positive_reference] / reference_map[positive_reference] - 1.0)
    residual[(reference_map == 0.0) & (reconstructed_map == 0.0)] = 0.0
    residual_log = np.log10(np.maximum(residual, 1e-16))
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    values = (reference_log, reconstructed_log, residual_log)
    titles = (
        "log10 HPCA Level-2 flux angular integral (reference)",
        "log10 production-kernel reconstruction from HPCA PSD",
        "log10 reconstruction/reference relative error (display floor 1e-16)",
    )
    for axis, data, title in zip(axes, values, titles):
        image = axis.imshow(data, origin="lower", aspect="auto", interpolation="nearest")
        axis.set_title(title, loc="left", fontsize=10)
        energy_ticks = (0, 20, 40, 62)
        axis.set_yticks(energy_ticks)
        axis.set_yticklabels([f"{energy_axis[index]:.4g}" for index in energy_ticks])
        axis.set_ylabel("H+ energy [eV]")
        fig.colorbar(image, ax=axis, pad=0.01)
    tick = np.linspace(0, n_time - 1, min(7, n_time), dtype=int)
    axes[-1].set_xticks(tick)
    axes[-1].set_xticklabels([labels[index] for index in tick])
    axes[-1].set_xlabel("acquisition midpoint UTC")
    zero_spectrum_cells = int(np.count_nonzero(reference_map == 0.0))
    fig.suptitle("MMS01 HPCA product reconstruction; all 63 archived energy centroids; "
                 f"blank flux cells are {zero_spectrum_cells} exact zeroes")
    generated += _save(fig, figures, "mms01_time_energy_comparison", formats)

    # Scatter uses a deterministic stride only for rendering cost.  Metrics and CSV
    # evidence above retain every directional cell, including exact physical zeroes.
    positive = [row for row in directional
                if float(row["reference_differential_intensity"]) > 0.0]
    stride = max(1, len(positive) // 50000)
    shown = positive[::stride]
    reference = np.asarray([float(row["reference_differential_intensity"]) for row in shown])
    reconstructed = np.asarray([
        float(row["kernel_reconstructed_differential_intensity"]) for row in shown
    ])
    low, high = (min(reference.min(), reconstructed.min()),
                 max(reference.max(), reconstructed.max()))
    fig, axis = plt.subplots(figsize=(6.5, 6))
    axis.loglog(reference, reconstructed, ".", ms=1.5, alpha=0.25,
                label=f"deterministic display sample ({len(shown):,}/{len(positive):,} positive cells)")
    axis.loglog([low, high], [low, high], "k-", lw=1, label="unity")
    axis.set_xlabel("HPCA Level-2 differential intensity (reference)")
    axis.set_ylabel("Production-kernel value from HPCA PSD (reconstruction)")
    axis.grid(True, which="both", alpha=0.2)
    axis.legend(fontsize=8)
    axis.set_title("MMS01 HPCA directional-product reconstruction; "
                   f"exact zero cells: {summary['directional']['zero_reference_cells']:,}")
    generated += _save(fig, figures, "mms01_directional_scatter", formats)

    xyz = np.asarray([[float(row[name]) for name in ("x_gse_km", "y_gse_km", "z_gse_km")]
                      for row in positions])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(xyz[:, 0], xyz[:, 1], "o-")
    axes[1].plot(xyz[:, 0], xyz[:, 2], "o-")
    axes[0].set(xlabel="X GSE [km]", ylabel="Y GSE [km]")
    axes[1].set(xlabel="X GSE [km]", ylabel="Z GSE [km]")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.set_aspect("equal", adjustable="datalim")
    fig.suptitle("Context only: MMS1 locations from bracketed 30 s MEC ephemeris "
                 "(no particle comparison)")
    generated += _save(fig, figures, "mms01_mms1_locations", formats)

    manifest = {
        "schema": "amps-mms01-figures-v2", "created_utc": utc_now(),
        "profile": summary["profile"], "comparison_status": summary["status"],
        "formats": list(formats),
        "artifacts": [{"path": str(path.resolve()), "sha256": sha256_file(path)}
                      for path in generated],
    }
    atomic_json(output_root / "figure_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--formats", nargs="+", default=("png", "eps"))
    args = parser.parse_args()
    manifest = plot(args.output_root.resolve(), tuple(args.formats))
    print(f"Wrote {len(manifest['artifacts'])} MMS01 figure files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
