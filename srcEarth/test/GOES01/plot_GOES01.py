#!/usr/bin/env python3
"""Render explicitly labeled GOES01 observation-versus-AMPS evidence figures."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from goes01_common import read_csv, sha256_file


COLORS = {"WEST_MINUS_X": "#1f77b4", "EAST_PLUS_X": "#d62728"}


def _time(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def _finish(figure: plt.Figure, base: Path, formats: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    figure.tight_layout()
    for extension in formats:
        path = base.with_suffix("." + extension.lower())
        figure.savefig(path, dpi=180 if extension.lower() == "png" else None)
        paths.append(path)
    plt.close(figure)
    return paths


def plot(output_root: Path, formats: tuple[str, ...]) -> list[Path]:
    """Plot the saved comparison only; plotting can never change qualification."""
    summary = json.loads((output_root / "comparison_summary.json").read_text())
    comparison = Path(summary["comparison_table"])
    if sha256_file(comparison) != summary["comparison_table_sha256"]:
        raise ValueError("comparison table changed after qualification")
    rows = read_csv(comparison)
    prepared = json.loads((output_root / "GOES16" / "prepared_run.json").read_text())
    ledger = read_csv(Path(prepared["selection_ledger_file"]))
    figures = output_root / "figures"
    figures.mkdir(exist_ok=True)
    products: list[Path] = []
    heads = ("WEST_MINUS_X", "EAST_PLUS_X")

    # The time-series title says exactly what the two curves mean. "AMPS" is never
    # called a reconstruction of SGPS: it is the independent ERNE-boundary transport
    # result evaluated at the center of the SGPS integration interval.
    figure, axis = plt.subplots(figsize=(11.0, 5.5))
    for head in heads:
        selected = sorted((row for row in rows if row["head"] == head),
                          key=lambda row: row["comparison_midpoint_utc"])
        times = [_time(row["comparison_midpoint_utc"]) for row in selected]
        observed = np.asarray([float(row["reference"]) for row in selected])
        uncertainty = np.asarray([float(row["reference_uncertainty"]) for row in selected])
        model = np.asarray([float(row["model"]) for row in selected])
        color = COLORS[head]
        axis.errorbar(times, observed, yerr=uncertainty, fmt="o", ms=3, lw=0.8,
                      capsize=1.5, color=color, alpha=0.70,
                      label=f"GOES-16 SGPS {head} P8A observation")
        axis.plot(times, model, "-", lw=1.5, color=color,
                  label=f"AMPS transported intensity folded to {head} P8A")
    axis.set_yscale("log")
    axis.set_ylabel(r"Differential proton intensity [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$]")
    axis.set_xlabel("SGPS five-minute-bin midpoint [UTC]")
    axis.set_title("GOES01: AMPS near-Earth high-energy intensity vs GOES-16/SGPS P8A")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8, ncol=2)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    products += _finish(figure, figures / "p8a_intensity_timeseries", formats)

    figure, axis = plt.subplots(figsize=(11.0, 4.8))
    for head in heads:
        selected = sorted((row for row in rows if row["head"] == head),
                          key=lambda row: row["comparison_midpoint_utc"])
        axis.plot([_time(row["comparison_midpoint_utc"]) for row in selected],
                  [float(row["model_observed_ratio"]) for row in selected],
                  "o-", ms=3, lw=1, color=COLORS[head], label=head)
    axis.axhspan(0.5, 2.0, color="#888888", alpha=0.15, label="factor-two band")
    axis.axhline(1.0, color="black", lw=1)
    # symlog retains a physical-zero AMPS result. An ordinary log axis would hide
    # exactly the strongest access failure from the residual evidence.
    axis.set_yscale("symlog", linthresh=0.01)
    axis.set_ylabel("AMPS / GOES-16 SGPS P8A")
    axis.set_xlabel("SGPS five-minute-bin midpoint [UTC]")
    axis.set_title("GOES01 absolute-normalization residual (independent ERNE boundary)")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=8)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    products += _finish(figure, figures / "p8a_model_observed_ratio", formats)

    figure, axis = plt.subplots(figsize=(6.6, 6.2))
    all_values: list[float] = []
    for head in heads:
        selected = [row for row in rows if row["head"] == head]
        observed = np.asarray([float(row["reference"]) for row in selected])
        model = np.asarray([float(row["model"]) for row in selected])
        all_values.extend(observed.tolist() + model.tolist())
        axis.scatter(observed, model, s=22, alpha=0.70, color=COLORS[head], label=head)
    positive_values = [value for value in all_values if value > 0.0]
    if not positive_values:
        raise ValueError("scatter plot has no positive model or observation values")
    low, high = min(positive_values), max(positive_values)
    grid = np.geomspace(low, high, 200)
    axis.plot(grid, grid, "k-", lw=1.2, label="perfect agreement")
    axis.plot(grid, 2.0 * grid, "k--", lw=0.8, label="factor two")
    axis.plot(grid, 0.5 * grid, "k--", lw=0.8)
    # A symmetric-log origin keeps physical-zero model predictions visible rather
    # than dropping them from an otherwise logarithmic model-observation plot.
    linear_threshold = low / 10.0
    axis.set_xscale("symlog", linthresh=linear_threshold)
    axis.set_yscale("symlog", linthresh=linear_threshold)
    axis.set_xlabel(r"Observed SGPS P8A intensity [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$]")
    axis.set_ylabel(r"AMPS P8A-folded intensity [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$]")
    axis.set_title("GOES01 model-observation intensity comparison")
    axis.grid(True, which="both", alpha=0.25); axis.legend(fontsize=8)
    products += _finish(figure, figures / "p8a_model_observed_scatter", formats)

    figure, axis = plt.subplots(figsize=(11.0, 4.6))
    finite = [row for row in ledger if row["two_head_envelope_separation"]]
    selected = [row for row in finite if row["eligible"] == "1"]
    rejected = [row for row in finite if "two_head_inconsistent" in row["exclusion_reasons"]]
    axis.scatter([_time(row["comparison_midpoint_utc"]) for row in finite],
                 [float(row["two_head_envelope_separation"]) for row in finite], s=7, color="#aaaaaa",
                 alpha=0.45, label="all finite quality-audit candidates")
    axis.scatter([_time(row["comparison_midpoint_utc"]) for row in selected],
                 [float(row["two_head_envelope_separation"]) for row in selected], s=23, color="#2ca02c",
                 label="profile-selected comparison epochs")
    if rejected:
        axis.scatter([_time(row["comparison_midpoint_utc"]) for row in rejected],
                     [float(row["two_head_envelope_separation"]) for row in rejected], s=12, color="#d62728",
                     label="excluded: two-head inconsistency")
    axis.axhline(1.0, color="black", ls="--", lw=1,
                 label="delivered uncertainty intervals just overlap")
    # symlog retains exact zero (perfect head equality), which an ordinary log axis
    # would silently drop from this observation-selection audit.
    axis.set_yscale("symlog", linthresh=0.01)
    axis.set_ylabel(r"$|J_W-J_E| / (\sigma_W+\sigma_E)$")
    axis.set_xlabel("SGPS five-minute-bin midpoint [UTC]")
    axis.set_title("GOES01 observation-only two-head consistency selection")
    axis.grid(True, which="both", alpha=0.25); axis.legend(fontsize=8)
    axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    products += _finish(figure, figures / "p8a_two_head_consistency", formats)

    manifest = {
        "schema": "amps-goes01-figure-manifest-v1", "formats": list(formats),
        "figures": {path.name: sha256_file(path) for path in products},
    }
    from goes01_common import atomic_json
    atomic_json(output_root / "figure_manifest.json", manifest)
    return products
