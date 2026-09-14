#!/usr/bin/env python3
"""Independent conservative finite-volume reference for CV03.

Cell probabilities satisfy ``df/dt = d/ds(kappa df/ds)`` on a periodic line.
Face diffusivities are harmonic means and explicit substeps obey the configured
stability limit.  This solver shares equations and SI inputs with srcSEP, but
no production transport code, random stream, or coefficient derivative.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def solve(config: dict) -> list[float]:
    physics, packet = config["physics"], config["packet"]
    numerics, reference = config["numerics"], config["reference"]
    length = float(physics["line_length_m"])
    kappa0 = float(physics["kappa0_m2_per_s"])
    amplitude = float(physics["sinusoidal_amplitude"])
    cells = int(reference["finite_volume_cells"])
    output_cells = int(numerics["profile_bins"])
    if cells % output_cells:
        raise ValueError("finite_volume_cells must be divisible by profile_bins")
    dx = length / cells
    centers = [(i + 0.5) * dx for i in range(cells)]
    kappa = [kappa0 * (1.0 + amplitude * math.sin(2.0 * math.pi * x / length))
             for x in centers]
    face = [2.0 * kappa[i] * kappa[(i + 1) % cells] /
            (kappa[i] + kappa[(i + 1) % cells]) for i in range(cells)]
    probability = [0.0] * cells
    packet_cell = int(float(packet["center_m"]) / dx) % cells
    probability[packet_cell] = 1.0
    final = float(numerics["final_time_s"])
    stable = float(reference["stability_factor"]) * dx * dx / max(kappa)
    steps = max(1, math.ceil(final / stable))
    dt = final / steps
    for _ in range(steps):
        updated = [0.0] * cells
        for i in range(cells):
            left = (i - 1) % cells
            right = (i + 1) % cells
            # Probability (not density) uses the same dx on every cell, so the
            # conservative face-flux update has coefficient dt/dx^2.
            updated[i] = probability[i] + dt / (dx * dx) * (
                face[i] * (probability[right] - probability[i]) -
                face[left] * (probability[i] - probability[left]))
        probability = updated
    group = cells // output_cells
    return [sum(probability[i * group:(i + 1) * group])
            for i in range(output_cells)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.input.read_text(encoding="utf-8"))
    result = solve(config)
    length = float(config["physics"]["line_length_m"])
    width = length / len(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("bin_index", "bin_left_m", "bin_right_m", "probability"))
        for index, probability in enumerate(result):
            writer.writerow((index, index * width, (index + 1) * width, probability))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
