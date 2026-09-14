#!/usr/bin/env python3
"""Independent characteristic reference for CV04 adiabatic cooling.

For constant divergence, ``p=p0 exp[-(div U)t/3]``.  For constant-speed
spherical wind, ``div U=2U/r`` and integration along ``r=r0+Ut`` gives
``p=p0 (r/r0)^(-2/3)``.  Relativistic energy conversion is evaluated here,
outside the linked application, for proton and alpha-particle mass numbers.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


EV_J = 1.602176634e-19


def momentum(energy_mev_per_nucleon: float, mass_number: float,
             proton_mass: float, light_speed: float) -> float:
    kinetic_per_nucleon = energy_mev_per_nucleon * 1.0e6 * EV_J
    kinetic = mass_number * kinetic_per_nucleon
    mass = mass_number * proton_mass
    return math.sqrt(kinetic * kinetic / (light_speed * light_speed) +
                     2.0 * mass * kinetic)


def energy_per_nucleon(p: float, mass_number: float, proton_mass: float,
                       light_speed: float) -> float:
    mass = mass_number * proton_mass
    ratio = p / (mass * light_speed)
    # Rationalized independently from the C++ implementation to avoid
    # catastrophic subtraction at the nonrelativistic end of the sweep.
    kinetic = (mass * light_speed * light_speed * ratio * ratio /
               (math.sqrt(1.0 + ratio * ratio) + 1.0))
    return kinetic / mass_number / (1.0e6 * EV_J)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.input.read_text(encoding="utf-8"))
    constants, physics = config["constants"], config["physics"]
    particles, numerics = config["particles"], config["numerics"]
    spectral_index = float(particles["momentum_spectral_index"])
    c, mp = float(constants["speed_of_light_m_per_s"]), float(constants["proton_mass_kg"])
    rows = []
    for scenario, final in (("constant", float(numerics["constant_final_time_s"])),
                            ("spherical", float(numerics["spherical_final_time_s"]))):
        for mass_number in map(float, particles["atomic_mass_numbers"]):
            for energy0 in map(float, particles["energies_mev_per_nucleon"]):
                p0 = momentum(energy0, mass_number, mp, c)
                weight = (p0 / 1.0e-19) ** (-spectral_index)
                if scenario == "constant":
                    radius = float(physics["initial_radius_m"])
                    p = p0 * math.exp(-float(physics["constant_divergence_per_s"]) * final / 3.0)
                else:
                    r0, speed = float(physics["initial_radius_m"]), float(physics["radial_wind_speed_m_per_s"])
                    radius = r0 + speed * final
                    p = p0 * (radius / r0) ** (-2.0 / 3.0)
                rows.append((scenario, mass_number, energy0, final, radius, p,
                             energy_per_nucleon(p, mass_number, mp, c), weight))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("scenario", "atomic_mass_number", "energy0_mev_per_nucleon",
                         "time_s", "radius_m", "momentum_kg_m_per_s",
                         "kinetic_energy_mev_per_nucleon", "weight"))
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
