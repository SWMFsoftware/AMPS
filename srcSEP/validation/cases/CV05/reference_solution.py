#!/usr/bin/env python3
"""Independent magnetic-focusing characteristics for CV05.

With ``D_mumu=0`` and constant ``g=d ln|B|/ds``, the production equation is
``dmu/dt=-(v g/2)(1-mu^2)`` and ``ds/dt=v mu``.  Integrating gives a hyperbolic
closed form for interior pitch angles; the exactly field-aligned endpoints are
handled as their limiting characteristics to avoid ``atanh(+-1)``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def characteristic(mu0: float, s0: float, speed: float, gradient: float,
                   time_s: float) -> tuple[float, float]:
    if mu0 <= -1.0:
        return s0 - speed * time_s, -1.0
    if mu0 >= 1.0:
        return s0 + speed * time_s, 1.0
    q0 = math.atanh(mu0)
    rate = -0.5 * speed * gradient
    q = q0 + rate * time_s
    mu = math.tanh(q)
    if rate == 0.0:
        return s0 + speed * mu0 * time_s, mu0
    position = s0 + speed / rate * (math.log(math.cosh(q)) -
                                     math.log(math.cosh(q0)))
    return position, mu


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.input.read_text(encoding="utf-8"))
    physics, particles, numerics = config["physics"], config["particles"], config["numerics"]
    speed = float(physics["speed_m_per_s"])
    length = float(physics["focusing_length_m"])
    s0 = float(physics["initial_position_m"])
    final = float(numerics["final_time_s"])
    sample_count = int(numerics["sample_count"])
    rows = []
    for sign in map(float, physics["gradient_signs"]):
        gradient = sign / length
        for sample_index in range(sample_count + 1):
            time_s = final * sample_index / sample_count
            for particle_id, mu0 in enumerate(map(float, particles["initial_mu"]), 1):
                position, mu = characteristic(mu0, s0, speed, gradient, time_s)
                invariant = (1.0 - mu * mu) / math.exp(gradient * position)
                rows.append((sign, particle_id, mu0, 1.0 + 0.4 * mu0, time_s,
                             position, mu, invariant))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("gradient_sign", "particle_id", "initial_mu", "weight",
                         "time_s", "position_m", "mu", "magnetic_moment_invariant"))
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
