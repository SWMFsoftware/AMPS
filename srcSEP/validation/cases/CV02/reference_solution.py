#!/usr/bin/env python3
"""Independent CV02 Gaussian Green-function reference.

The linked model reports probability in finite bins, so this reference uses
error-function bin integrals rather than point samples at bin centers.  That
removes a plotting/discretization bias from the scientific acceptance metric.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def bin_probability(left: float, right: float, center: float,
                    kappa: float, time_s: float) -> float:
    scale = math.sqrt(4.0 * kappa * time_s)
    return 0.5 * (math.erf((right - center) / scale) -
                  math.erf((left - center) / scale))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.input.read_text(encoding="utf-8"))
    physics, domain, numerics = config["physics"], config["domain"], config["numerics"]
    center = float(domain["packet_center_m"])
    half = float(domain["profile_half_width_m"])
    count = int(numerics["profile_bins"])
    final = float(numerics["final_time_s"])
    kappa = float(physics["kappa_parallel_m2_per_s"])
    width = 2.0 * half / count
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("bin_index", "bin_left_m", "bin_right_m", "probability",
                         "mean_m", "variance_m2", "skewness", "kurtosis"))
        for index in range(count):
            left = center - half + index * width
            right = left + width
            writer.writerow((index, left, right,
                             bin_probability(left, right, center, kappa, final),
                             center, 2.0 * kappa * final, 0.0, 3.0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
