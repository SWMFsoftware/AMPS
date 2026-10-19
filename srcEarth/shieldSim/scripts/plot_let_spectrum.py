#!/usr/bin/env python3
"""Render a shieldSim total step-fluence LET spectrum and its precision.

The shieldSim LET ``.dat`` product uses Tecplot-style zones.  This utility
deliberately plots only ``ZONE T=\"total ...\"``: particle-component zones are
useful diagnostics, but summing their already-rounded text values can differ
slightly from the authoritative total accumulated by the scorer.

Zero-content bins are omitted on logarithmic axes; they are never replaced by
an arbitrary epsilon.  The uncertainty band is the reported one-standard-error
Monte Carlo uncertainty divided by the LET-bin width, matching ``dPhi_dLET``.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

# This script is intended for validation and batch/HPC runs, where a display
# server is normally unavailable.  Select the backend before importing pyplot.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the total zone of a shieldSim step-fluence LET file."
    )
    parser.add_argument("datfile", type=Path, help="shieldSim *_let_spectrum.dat file")
    parser.add_argument(
        "--output-prefix",
        type=Path,
        help="output path without an extension (default: input stem)",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["pdf"],
        choices=("pdf", "png", "svg"),
        help="one or more output formats (default: pdf)",
    )
    parser.add_argument("--dpi", type=int, default=240, help="PNG resolution")
    parser.add_argument(
        "--summary-json",
        type=Path,
        help="optional pooled summary whose LET means are drawn as reference lines",
    )
    return parser.parse_args()


def read_total_zone(path: Path) -> tuple[str, np.ndarray]:
    """Return the first total zone label and its numeric rows.

    The explicit ``I=`` count prevents component-zone rows from leaking into
    the total if more zones follow it in a normal shieldSim output file.
    """

    if not path.is_file():
        raise FileNotFoundError(f"LET spectrum not found: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        match = re.match(
            r'^\s*ZONE\s+T="(total(?:\s*\|[^\"]*)?)"\s*,\s*I=(\d+)',
            line,
            flags=re.IGNORECASE,
        )
        if not match:
            continue

        zone_label = match.group(1)
        row_count = int(match.group(2))
        rows: list[list[float]] = []
        for raw in lines[index + 1 : index + 1 + row_count]:
            values = [float(value) for value in raw.split()]
            if len(values) != 9:
                raise ValueError(
                    f"expected 9 LET columns in {path}, found {len(values)}"
                )
            rows.append(values)
        if len(rows) != row_count:
            raise ValueError(
                f"total zone in {path} declares {row_count} rows; found {len(rows)}"
            )
        return zone_label, np.asarray(rows, dtype=float)

    raise ValueError(f'no ZONE T="total ..." found in {path}')


def read_let_means(path: Path | None) -> tuple[float | None, float | None]:
    """Read optional pooled track- and dose-mean LET values in keV/um."""

    if path is None:
        return None, None
    document = json.loads(path.read_text(encoding="utf-8"))
    let = document.get("quantities", {}).get("LET", {})
    return let.get("track_mean_keV_per_um"), let.get("dose_mean_keV_per_um")


def make_figure(
    zone_label: str,
    data: np.ndarray,
    track_mean: float | None,
    dose_mean: float | None,
) -> plt.Figure:
    """Build a two-panel spectrum/precision figure without altering data."""

    low, high, center = data[:, 0], data[:, 1], data[:, 2]
    fluence, fluence_se, differential = data[:, 4], data[:, 5], data[:, 6]
    widths = high - low
    differential_se = np.divide(
        fluence_se,
        widths,
        out=np.zeros_like(fluence_se),
        where=widths > 0.0,
    )

    occupied = np.isfinite(differential) & (differential > 0.0)
    if not np.any(occupied):
        raise ValueError("the total LET zone has no positive-content bins")

    relative_se = np.divide(
        fluence_se,
        fluence,
        out=np.full_like(fluence, np.nan),
        where=fluence > 0.0,
    )

    figure, (spectrum_axis, precision_axis) = plt.subplots(
        1, 2, figsize=(13.2, 5.7), constrained_layout=False
    )

    spectrum_axis.loglog(
        center[occupied],
        differential[occupied],
        color="#1565c0",
        linewidth=1.8,
        label="Total charged-particle step fluence",
    )
    # A zero or negative lower confidence edge has no logarithmic coordinate.
    # Mask only that edge; retain the measured central value and upper edge.
    lower_band = differential - differential_se
    lower_band[lower_band <= 0.0] = np.nan
    upper_band = differential + differential_se
    spectrum_axis.fill_between(
        center[occupied],
        lower_band[occupied],
        upper_band[occupied],
        color="#42a5f5",
        alpha=0.28,
        linewidth=0.0,
        label=r"Monte Carlo $\pm 1\,\mathrm{SE}$",
    )

    if track_mean is not None:
        spectrum_axis.axvline(
            track_mean,
            color="#2e7d32",
            linestyle="--",
            linewidth=1.4,
            label=f"Track mean = {track_mean:.3g} keV/um",
        )
    if dose_mean is not None:
        spectrum_axis.axvline(
            dose_mean,
            color="#c62828",
            linestyle=":",
            linewidth=1.7,
            label=f"Dose mean = {dose_mean:.3g} keV/um",
        )

    spectrum_axis.set_xlabel("Unrestricted electronic LET (keV/um)")
    spectrum_axis.set_ylabel(
        r"$d\Phi/dL$ [cm$^{-2}$ primary$^{-1}$ (keV/um)$^{-1}$]"
    )
    spectrum_axis.legend(loc="best", fontsize=8.5, frameon=False)

    precision_axis.loglog(
        center[occupied],
        relative_se[occupied],
        color="#6a1b9a",
        marker="o",
        markersize=2.6,
        linewidth=1.1,
        label="Relative standard error",
    )
    precision_axis.axhline(
        0.10,
        color="#d84315",
        linestyle="--",
        linewidth=1.3,
        label="10% requested accuracy",
    )
    precision_axis.set_xlabel("Unrestricted electronic LET (keV/um)")
    precision_axis.set_ylabel("Relative standard error")
    precision_axis.legend(loc="best", fontsize=8.5, frameon=False)

    xmin = low[np.flatnonzero(occupied)[0]]
    xmax = high[np.flatnonzero(occupied)[-1]]
    for axis in (spectrum_axis, precision_axis):
        axis.set_xlim(xmin, xmax)
        axis.grid(True, which="major", linewidth=0.5, alpha=0.4)
        axis.grid(True, which="minor", linestyle=":", linewidth=0.35, alpha=0.3)
        axis.tick_params(which="both", direction="in", top=True, right=True)

    histories = re.search(r"(\d+)\s+pooled events", zone_label, re.IGNORECASE)
    processes = re.search(r"(\d+)\s+processes", zone_label, re.IGNORECASE)
    run_label = zone_label
    if histories and processes:
        run_label = (
            f"{int(histories.group(1)):,} histories pooled from "
            f"{processes.group(1)} independent processes"
        )
    figure.suptitle("ShieldSim Pooled LET Spectrum", fontsize=15, weight="bold")
    figure.text(
        0.5,
        0.925,
        run_label,
        ha="center",
        fontsize=10,
        color="#333333",
    )
    figure.text(
        0.5,
        0.015,
        "Step-fluence LET at transported-step midpoint; source normalization is one primary per history.",
        ha="center",
        fontsize=8.5,
        color="#444444",
    )
    figure.subplots_adjust(left=0.085, right=0.98, bottom=0.13, top=0.86, wspace=0.25)
    return figure


def main() -> int:
    args = parse_args()
    zone_label, data = read_total_zone(args.datfile)
    track_mean, dose_mean = read_let_means(args.summary_json)
    figure = make_figure(zone_label, data, track_mean, dose_mean)

    prefix = args.output_prefix or args.datfile.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for output_format in dict.fromkeys(args.formats):
        output = prefix.with_suffix(f".{output_format}")
        save_options = {"bbox_inches": "tight"}
        if output_format == "png":
            save_options["dpi"] = args.dpi
        figure.savefig(output, **save_options)
        print(output)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
