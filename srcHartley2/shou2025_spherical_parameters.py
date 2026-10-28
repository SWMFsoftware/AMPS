#!/usr/bin/env python3
"""Reproduce the Shou et al. (2025) parameters used by the DIM1 case."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

from comet_shape_properties import properties, read_nastran


PAPER_HELIOCENTRIC_DISTANCE_AU = 1.06
PAPER_DIRECT_RATES = {"H2O": 2.0e27, "CO2": 3.5e26}
PAPER_NIGHT_FACTORS = {"H2O": 0.1, "CO2": 0.3}

# Shou et al. do not publish collision cross sections or a surface mean free
# path. The existing AMPS comet_sphere setup uses this constant for H2O-H2O,
# H2O-CO2, and CO2-CO2 collisions (input/comet_sphere.input).
CORE_COLLISION_CROSS_SECTION_M2 = 3.4e-19
SPECIES_MASS_AMU = {"H2O": 18.0, "CO2": 44.0}
BOLTZMANN_J_K = 1.380649e-23
ATOMIC_MASS_KG = 1.66053906660e-27

# Shou et al. use the 1.29-AU profile from Tenishev et al. (2008). The
# Hartley2 implementation labels the matching table column as 1.3 AU and uses
# 172 K on the night side (Comet.cpp, GetSurfaceTemperature, ndist=0).
TABLE_COLUMN = 1
NIGHT_TEMPERATURE_K = 172.0


def read_temperature_table(path: Path) -> list[tuple[float, ...]]:
    source = path.read_text(encoding="ascii")
    marker = "const double SurfaceTemp[90][6]={"
    try:
        table_source = source.split(marker, 1)[1].split("};", 1)[0]
    except IndexError as error:
        raise ValueError(f"{path}: SurfaceTemp[90][6] was not found") from error

    rows: list[tuple[float, ...]] = []
    number = re.compile(r"[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?")
    for line in table_source.splitlines():
        values = tuple(float(value) for value in number.findall(line))
        if values:
            if len(values) != 6:
                raise ValueError(f"{path}: malformed SurfaceTemp row: {line.strip()}")
            rows.append(values)

    if len(rows) != 90:
        raise ValueError(f"{path}: expected 90 SurfaceTemp rows, found {len(rows)}")
    return rows


def integrate_weight(lower: float, upper: float, night_factor: float) -> float:
    """Integrate max(cos(theta), a_night)*sin(theta) over one interval."""
    transition = math.acos(night_factor)
    if upper <= transition:
        return 0.5 * (math.sin(upper) ** 2 - math.sin(lower) ** 2)
    if lower >= transition:
        return night_factor * (math.cos(lower) - math.cos(upper))
    return integrate_weight(lower, transition, night_factor) + integrate_weight(
        transition, upper, night_factor
    )


def effective_temperature(
    rows: list[tuple[float, ...]], night_factor: float
) -> float:
    """Return the source-number-weighted spherical boundary temperature."""
    weighted_temperature = 0.0
    weight = 0.0

    # GetSurfaceTemperature truncates the day-side zenith angle to an integer,
    # so row i applies to the complete [i, i+1) degree interval.
    for index, row in enumerate(rows):
        lower = math.radians(index)
        upper = math.radians(index + 1)
        interval_weight = integrate_weight(lower, upper, night_factor)
        weighted_temperature += row[TABLE_COLUMN] * interval_weight
        weight += interval_weight

    night_weight = night_factor  # integral sin(theta) dtheta from pi/2 to pi
    weighted_temperature += NIGHT_TEMPERATURE_K * night_weight
    weight += night_weight
    return weighted_temperature / weight


def outgoing_surface_density(flux: float, temperature: float, mass_amu: float) -> float:
    """Density immediately above a zero-drift upward half-Maxwellian surface."""
    mass = mass_amu * ATOMIC_MASS_KG
    mean_outward_normal_speed = math.sqrt(
        2.0 * BOLTZMANN_J_K * temperature / (math.pi * mass)
    )
    return flux / mean_outward_normal_speed


def main() -> None:
    source_dir = Path(__file__).resolve().parent
    repository = source_dir.parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--temperature-table",
        type=Path,
        default=source_dir / "CG-BC.h",
    )
    parser.add_argument(
        "--shape",
        type=Path,
        default=repository / "input" / "hartley2.shape.nas.gz",
    )
    args = parser.parse_args()

    rows = read_temperature_table(args.temperature_table)
    shape = properties(*read_nastran(args.shape))
    area_m2 = float(shape["area"]) * 1.0e6  # tracked NASTRAN coordinates are km
    radius_m = float(shape["area_radius"]) * 1.0e3

    print(f"heliocentric_distance_au={PAPER_HELIOCENTRIC_DISTANCE_AU:.12g}")
    print(f"surface_area_m2={area_m2:.12f}")
    print(f"area_equivalent_radius_m={radius_m:.12f}")
    total_outgoing_density = 0.0
    for species in ("H2O", "CO2"):
        rate = PAPER_DIRECT_RATES[species]
        temperature = effective_temperature(rows, PAPER_NIGHT_FACTORS[species])
        flux = rate / area_m2
        density = outgoing_surface_density(
            flux, temperature, SPECIES_MASS_AMU[species]
        )
        total_outgoing_density += density
        print(f"{species}_direct_rate_s-1={rate:.12e}")
        print(f"{species}_uniform_flux_m-2_s-1={flux:.12e}")
        print(f"{species}_effective_temperature_k={temperature:.12f}")
        print(f"{species}_bulk_speed_m_s=0")
        print(f"{species}_outgoing_surface_density_m-3={density:.12e}")

    # Conventional hard-sphere estimate. The sqrt(2) factor is exact for one
    # equilibrium Maxwellian species; using it for this two-temperature,
    # outward-only mixture is an explicitly identified resolution estimate.
    mean_free_path = 1.0 / (
        math.sqrt(2.0)
        * total_outgoing_density
        * CORE_COLLISION_CROSS_SECTION_M2
    )
    print(f"core_collision_cross_section_m2={CORE_COLLISION_CROSS_SECTION_M2:.12e}")
    print(f"total_outgoing_surface_density_m-3={total_outgoing_density:.12e}")
    print(f"estimated_surface_mean_free_path_m={mean_free_path:.12f}")


if __name__ == "__main__":
    main()
