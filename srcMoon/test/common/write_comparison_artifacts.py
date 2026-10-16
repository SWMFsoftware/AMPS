#!/usr/bin/env python3
"""Render an executed srcMoon comparison CSV as matching PNG and EPS plots.

This utility never creates, interpolates, or changes science values.  Its
input must be the row-level CSV emitted by an I-series scorer after an actual
model/reference comparison.  A skipped test must not call it: creating an
empty or synthetic plot would make missing evidence look like a comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys


REQUIRED_COLUMNS = (
    "comparison_id",
    "x_value",
    "x_unit",
    "reference_value",
    "model_value",
    "y_unit",
    "uncertainty_value",
    "selection",
)


def read_comparison(path: Path) -> list[dict[str, str]]:
    """Read strict row-level evidence without filling gaps or sorting rows."""
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != list(REQUIRED_COLUMNS):
            raise ValueError(
                "comparison CSV header must be exactly: "
                + ",".join(REQUIRED_COLUMNS)
            )
        rows = list(reader)
    if not rows:
        raise ValueError("comparison CSV has no evaluated rows")

    # One scientific comparison ID must retain one unit pair, but an
    # integration test may legitimately emit several independent metrics
    # (I18 has four CDFs plus a binomial fraction).  Render those groups as
    # separate panels from the same immutable CSV instead of forcing unlike
    # variables onto one axis or dropping required evidence.
    units_by_id: dict[str, tuple[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        comparison_id = row["comparison_id"]
        if not comparison_id or not row["x_unit"] or not row["y_unit"]:
            raise ValueError(
                f"line {line_number}: comparison_id and units are required")
        units = (row["x_unit"], row["y_unit"])
        previous = units_by_id.setdefault(comparison_id, units)
        if previous != units:
            raise ValueError(
                f"line {line_number}: comparison_id {comparison_id!r} mixes "
                "incompatible x/y units")
        if row["selection"] not in {
                "calibration", "hold-out", "verification",
                "training-experimental"}:
            raise ValueError(
                f"line {line_number}: selection must be calibration, hold-out, "
                "verification, or training-experimental"
            )
        for column in ("x_value", "reference_value", "model_value"):
            try:
                value = float(row[column])
            except ValueError as exception:
                raise ValueError(
                    f"line {line_number}: {column} is not numeric"
                ) from exception
            if not math.isfinite(value):
                raise ValueError(f"line {line_number}: {column} is not finite")
        # Empty uncertainty is an explicit metadata state used when the
        # authoritative archive supplies no per-observation uncertainty.  It
        # is not interpreted as zero.  A supplied value remains strict,
        # finite, and non-negative.
        if row["uncertainty_value"]:
            try:
                uncertainty = float(row["uncertainty_value"])
            except ValueError as exception:
                raise ValueError(
                    f"line {line_number}: uncertainty_value is not numeric"
                ) from exception
            if not math.isfinite(uncertainty) or uncertainty < 0.0:
                raise ValueError(
                    f"line {line_number}: uncertainty_value must be finite "
                    "and non-negative when supplied"
                )
    return rows


def render(input_csv: Path, png: Path, eps: Path, title: str) -> None:
    """Render the same immutable rows to both requested graphical formats."""
    rows = read_comparison(input_csv)

    # Import only after validating the evidence file so malformed/missing data
    # are reported independently of an optional plotting installation.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(row["comparison_id"], []).append(row)
    figure, axes = plt.subplots(
        len(groups), 1, figsize=(7.0, 4.0 * len(groups)),
        constrained_layout=True, squeeze=False)
    for axis, (comparison_id, group) in zip(axes[:, 0], groups.items()):
        x = [float(row["x_value"]) for row in group]
        reference = [float(row["reference_value"]) for row in group]
        model = [float(row["model_value"]) for row in group]
        uncertainty = [
            float(row["uncertainty_value"])
            if row["uncertainty_value"] else math.nan
            for row in group
        ]
        if any(math.isfinite(value) for value in uncertainty):
            # Matplotlib accepts NaN error bars and omits them for rows whose
            # archive uncertainty is unavailable, preserving the row itself.
            axis.errorbar(
                x, reference, yerr=uncertainty, color="black", marker="o",
                linestyle="none", capsize=2, label="independent reference")
        else:
            axis.plot(
                x, reference, color="black", marker="o", linestyle="none",
                label="independent reference (uncertainty unavailable)")
        # An eclipse training overview contains different lunar sites at
        # successive epochs. A line would imply a continuous thermal history
        # of one site that was never sampled. Keep those actual samples as
        # unconnected markers; retain historical lines for other test types.
        experimental = any(row["selection"] == "training-experimental"
                           for row in group)
        axis.plot(x, model, color="tab:blue", marker="s", label="AMPS",
                  linestyle="none" if experimental else "-")
        axis.set_xlabel(f"x [{group[0]['x_unit']}]")
        axis.set_ylabel(f"value [{group[0]['y_unit']}]")
        axis.set_title(comparison_id if len(groups)>1 else title)
        # EPS has no transparency model. Use an explicit opaque light gray so
        # the vector artifact is warning-free and matches the PNG.
        axis.grid(True, linewidth=0.4, color="#b0b0b0")
        axis.legend(framealpha=1.0)
    if len(groups)>1:
        figure.suptitle(title)

    png.parent.mkdir(parents=True, exist_ok=True)
    eps.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png, dpi=180, transparent=False)
    figure.savefig(eps, format="eps", transparent=False)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render a real srcMoon comparison CSV to PNG and EPS"
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--png", type=Path, required=True)
    parser.add_argument("--eps", type=Path, required=True)
    parser.add_argument("--title", required=True)
    arguments = parser.parse_args()
    try:
        render(arguments.input, arguments.png, arguments.eps, arguments.title)
    except (OSError, ValueError, ImportError) as exception:
        print(f"ERROR: {exception}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
