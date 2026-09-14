#!/usr/bin/env python3
"""Independent closed-form reference for CV01 ballistic streaming.

The reference never imports or calls a srcSEP transport function.  It evaluates
the relativistic energy-speed relation and the characteristic

    s(t) = s0 + [U_parallel + mu v(E)] t

directly, then applies an analytical periodic map or the exact first crossing
of an open interval.  Keeping this code in a different language and process
from the C++ model adapter makes sign, normalization, and boundary mistakes in
the production kernel observable instead of self-consistent.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List


def load_input(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as stream:
        payload = json.load(stream)
    if not isinstance(payload, dict) or payload.get("case_id") != "CV01":
        raise ValueError("reference input must be a CV01 JSON object")
    return payload


def relativistic_state(config: Dict[str, object]) -> tuple[float, float]:
    """Return speed [m/s] and momentum [kg m/s] from kinetic energy."""
    constants = config["constants"]
    packet = config["particle_packet"]
    assert isinstance(constants, dict) and isinstance(packet, dict)
    c = float(constants["speed_of_light_m_per_s"])
    mass = float(constants["proton_mass_kg"])
    electron_volt_j = float(constants["electron_volt_j"])
    kinetic_j = float(packet["kinetic_energy_mev"]) * 1.0e6 * electron_volt_j
    gamma = 1.0 + kinetic_j / (mass * c * c)
    speed = c * math.sqrt(1.0 - 1.0 / (gamma * gamma))
    # E_total^2=(pc)^2+(mc^2)^2 gives this stable kinetic-energy form.
    momentum = math.sqrt(kinetic_j * kinetic_j / (c * c) + 2.0 * mass * kinetic_j)
    return speed, momentum


def read_particles(path: Path) -> List[Dict[str, float]]:
    particles: List[Dict[str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            particles.append({
                "particle_id": float(row["particle_id"]),
                "initial_s_m": float(row["initial_s_m"]),
                "mu": float(row["mu"]),
                "weight": float(row["weight"]),
            })
    if not particles:
        raise ValueError("initial particle CSV contains no particles")
    return particles


def _periodic(position: float, minimum: float, maximum: float) -> float:
    return minimum + (position - minimum) % (maximum - minimum)


def generate_rows(config: Dict[str, object], particles: Iterable[Dict[str, float]],
                  dt_s: float, boundary: str) -> Iterable[Dict[str, object]]:
    """Generate reference samples with exact characteristics and crossings."""
    line = config["field_line"]
    numerics = config["numerics"]
    assert isinstance(line, dict) and isinstance(numerics, dict)
    minimum = float(line["minimum_m"])
    maximum = float(line["maximum_m"])
    plasma_speed = float(line["plasma_speed_m_per_s"])
    final_time = float(numerics["final_time_s"])
    sample_interval = float(numerics["sample_interval_s"])
    speed, momentum = relativistic_state(config)
    sample_count = int(round(final_time / sample_interval))

    for sample_index in range(sample_count + 1):
        time_s = sample_index * sample_interval
        for particle in particles:
            initial = particle["initial_s_m"]
            mu = particle["mu"]
            directed_speed = plasma_speed + mu * speed
            unwrapped = initial + directed_speed * time_s
            crossing_time = math.nan
            active = True
            position = unwrapped
            if boundary == "periodic":
                position = _periodic(unwrapped, minimum, maximum)
            elif boundary == "open":
                if directed_speed > 0.0:
                    crossing_time = (maximum - initial) / directed_speed
                    boundary_position = maximum
                elif directed_speed < 0.0:
                    crossing_time = (minimum - initial) / directed_speed
                    boundary_position = minimum
                else:
                    boundary_position = initial
                if math.isfinite(crossing_time) and time_s >= crossing_time:
                    active = False
                    position = boundary_position
            else:
                raise ValueError("boundary must be periodic or open")

            yield {
                "time_s": time_s,
                "dt_s": dt_s,
                "boundary": boundary,
                "particle_id": int(particle["particle_id"]),
                "initial_mu": mu,
                "active": int(active),
                "position_m": position,
                "unwrapped_position_m": unwrapped,
                "mu": mu,
                "momentum_kg_m_per_s": momentum,
                "active_weight": particle["weight"] if active else 0.0,
                "escaped_weight": 0.0 if active else particle["weight"],
                "crossing_time_s": crossing_time if math.isfinite(crossing_time) else "",
            }


FIELDS = (
    "time_s", "dt_s", "boundary", "particle_id", "initial_mu", "active",
    "position_m", "unwrapped_position_m", "mu", "momentum_kg_m_per_s",
    "active_weight", "escaped_weight", "crossing_time_s",
)


def write_reference(config_path: Path, particles_path: Path, dt_s: float,
                    boundary: str, output_path: Path) -> None:
    config = load_input(config_path)
    particles = read_particles(particles_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(generate_rows(config, particles, dt_s, boundary))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="CV01 model-input JSON")
    parser.add_argument("--initial-particles", required=True, type=Path,
                        help="immutable packet fixture shared with the model")
    parser.add_argument("--dt-s", required=True, type=float,
                        help="model timestep recorded in the reference rows")
    parser.add_argument("--boundary", required=True, choices=("periodic", "open"))
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    write_reference(arguments.input, arguments.initial_particles,
                    arguments.dt_s, arguments.boundary, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
