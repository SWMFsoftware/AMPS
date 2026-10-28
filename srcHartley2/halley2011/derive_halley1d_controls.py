#!/usr/bin/env python3
"""Derive numerical controls for the gas-only Rubin (2011) Halley case.

The paper fixes the physical source rates, temperature, parent collision
cross sections, and photochemical excess energies.  It does not fix an AMPS
mesh, statistical weight, or timestep.  This script keeps the resulting
numerical estimates reproducible and separate from values quoted directly
from the paper.

The reaction list is read from ``input/CometPhysicsDatabase.pm`` so the speed
calculation follows the network compiled into the application.  The optional
methanol-to-formaldehyde channel is excluded, matching ``halley1d.input``.
"""

from __future__ import annotations

import math
import re
from pathlib import Path


AMU_KG = 1.66053906660e-27
BOLTZMANN_J_K = 1.380649e-23
EV_J = 1.602176634e-19

NUCLEUS_RADIUS_M = 6_000.0
OUTER_RADIUS_M = 1.0e9
SURFACE_TEMPERATURE_K = 180.0
RADIAL_CELL_COUNT = 2_000
BLOCK_CELL_COUNT = 8
COURANT_NUMBER = 0.5
PRIMARY_CHARACTERISTIC_SPEED_M_S = 400.0
TARGET_SURFACE_PARTICLES_PER_SPECIES_PER_CELL = 50.0

# Direct parent rates are Rubin et al. (2011), Table 2 abundances multiplied
# by Q(H2O)=7.0e29 s^-1. Nominal mass-channel integers make the paper's Table
# 4 stoichiometry explicit. The compiled application retains the shared AMPS
# catalog masses, including its rounded SI value for atomic carbon.
PARENT_SOURCE_RATE_S = {
    "H2O": 7.00e29,
    "CO": 9.10e28,
    "CO2": 1.75e28,
    "CH3OH": 1.19e28,
    "H2CO": 1.05e28,
    "NH3": 1.05e28,
    "CH4": 4.20e27,
    "C2H6": 2.80e27,
    "C2H4": 2.10e27,
    "C2H2": 2.10e27,
    "HCN": 1.40e27,
}

MASS_AMU = {
    "H": 1,
    "H2": 2,
    "C": 12,
    "CH": 13,
    "CH2": 14,
    "N": 14,
    "CH3": 15,
    "NH": 15,
    "CH4": 16,
    "NH2": 16,
    "O": 16,
    "NH3": 17,
    "OH": 17,
    "H2O": 18,
    "C2": 24,
    "C2H": 25,
    "C2H2": 26,
    "CN": 26,
    "HCN": 27,
    "C2H4": 28,
    "CO": 28,
    "C2H5": 29,
    "HCO": 29,
    "C2H6": 30,
    "H2CO": 30,
    "CH3OH": 32,
    "CO2": 44,
}

# Only the H2O row of Rubin Table 3 is required for the stated surface
# mean-free-path estimate.  The H2O self value is evaluated below from the
# Crifo temperature law; all other entries are constant SI cross sections.
H2O_PARENT_CROSS_SECTION_M2 = {
    "CO": 3.2e-19,
    "H2CO": 3.6e-19,
    "CO2": 3.4e-19,
    "CH3OH": 3.8e-19,
    "C2H6": 3.8e-19,
    "C2H4": 3.8e-19,
    "C2H2": 3.6e-19,
    "HCN": 3.4e-19,
    "NH3": 3.6e-19,
    "CH4": 3.8e-19,
}


def outward_spatial_mean_speed(species: str) -> float:
    """Return the mean positive normal speed of the spatial distribution.

    AMPS samples a flux-weighted distribution at the boundary.  Weighting
    those samples by their residence time in a cell recovers the outgoing
    half-Maxwellian spatial density, for which Gamma=n*sqrt(2*k*T/(pi*m)).
    """

    mass_kg = MASS_AMU[species] * AMU_KG
    return math.sqrt(
        2.0 * BOLTZMANN_J_K * SURFACE_TEMPERATURE_K / (math.pi * mass_kg)
    )


def derive_surface_mean_free_path() -> tuple[float, float, float]:
    """Return total outgoing density, effective sigma, and H2O mean free path."""

    surface_area_m2 = 4.0 * math.pi * NUCLEUS_RADIUS_M**2
    density_m3 = {
        species: rate / surface_area_m2 / outward_spatial_mean_speed(species)
        for species, rate in PARENT_SOURCE_RATE_S.items()
    }
    total_density_m3 = sum(density_m3.values())

    cross_section_m2 = dict(H2O_PARENT_CROSS_SECTION_M2)
    cross_section_m2["H2O"] = 1.66e-19 * (
        SURFACE_TEMPERATURE_K / 300.0
    ) ** -0.6
    effective_cross_section_m2 = sum(
        density_m3[species] * cross_section_m2[species]
        for species in PARENT_SOURCE_RATE_S
    ) / total_density_m3

    mean_free_path_m = 1.0 / (
        math.sqrt(2.0) * total_density_m3 * effective_cross_section_m2
    )
    return total_density_m3, effective_cross_section_m2, mean_free_path_m


def geometric_expansion_ratio(first_width_m: float) -> float:
    """Solve the ideal geometric-cell sum for the requested radial domain."""

    normalized_extent = (
        OUTER_RADIUS_M - NUCLEUS_RADIUS_M
    ) / first_width_m

    def width_sum(ratio: float) -> float:
        log_ratio = math.log(ratio)
        exponent = RADIAL_CELL_COUNT * log_ratio
        # C++ ``expm1`` returns positive infinity for the deliberately broad
        # initial bracket.  Mirror that behavior instead of letting Python's
        # overflow exception interrupt the bisection.
        if exponent > math.log(float.fromhex("0x1.fffffffffffffp+1023")):
            return math.inf
        return math.expm1(exponent) / math.expm1(log_ratio)

    lower, upper = 1.0, 2.0
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        if width_sum(middle) < normalized_extent:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def reference_injection_count() -> tuple[int, dict[str, float]]:
    """Estimate the AMPS reference count for 50 parent particles/cell.

    AMPS makes the expected injection count per local parent timestep equal
    to ReferenceInjectionParticleNumber.  The mean first-cell residence time
    is dx/u_plus and dt=C*dx/Vchar, so the expected occupancy is
    Nref*Vchar/(C*u_plus).  The ceiling of the largest parent requirement
    supplies at least the target for every directly injected species.
    """

    requirements = {
        species: TARGET_SURFACE_PARTICLES_PER_SPECIES_PER_CELL
        * COURANT_NUMBER
        * outward_spatial_mean_speed(species)
        / PRIMARY_CHARACTERISTIC_SPEED_M_S
        for species in PARENT_SOURCE_RATE_S
    }
    count = math.ceil(max(requirements.values()))
    occupancy = {
        species: count
        * PRIMARY_CHARACTERISTIC_SPEED_M_S
        / (COURANT_NUMBER * outward_spatial_mean_speed(species))
        for species in PARENT_SOURCE_RATE_S
    }
    return count, occupancy


def parse_products(text: str) -> list[str]:
    """Expand Table 4 product coefficients such as ``2*H``."""

    products: list[str] = []
    for raw_term in text.split("+"):
        term = raw_term.strip()
        match = re.fullmatch(r"(?:(\d+)\*)?([A-Z0-9]+)", term)
        if match is None:
            raise ValueError(f"cannot parse reaction product {term!r}")
        products.extend([match.group(2)] * int(match.group(1) or "1"))
    return products


def load_neutral_reactions() -> list[tuple[str, str, list[str], float]]:
    """Read neutral reaction IDs, parents, products, and excess energies."""

    repository = Path(__file__).resolve().parents[2]
    database_source = (repository / "input" / "CometPhysicsDatabase.pm").read_text()
    pattern = re.compile(
        r"\['([^']+)','([^']+)','([^']+)',[^,\]]+,([0-9.eE+\-]+)\]"
    )
    reactions: list[tuple[str, str, list[str], float]] = []
    for reaction_id, parent, product_text, energy_text in pattern.findall(
        database_source
    ):
        if product_text == "NONE" or reaction_id == "CH3OH_H2CO_H2":
            continue
        reactions.append(
            (reaction_id, parent, parse_products(product_text), float(energy_text))
        )
    if len(reactions) != 46:
        raise RuntimeError(
            "expected 46 neutral Rubin reactions with the methanol-H2CO branch off, "
            f"found {len(reactions)}"
        )
    return reactions


def derive_characteristic_speeds() -> dict[str, float]:
    """Propagate conservative recoil-speed bounds through the reaction graph.

    For product i, momentum and energy conservation bound the center-of-mass
    recoil by sqrt(2*E*M_other/(m_i*M_total)).  Adding that bound to the
    parent's bound is conservative under arbitrary alignment.  Direct parents
    start at the user-supplied 400 m/s.  A species that is both injected and
    produced photochemically may therefore finish above 400 m/s.
    """

    speed = {species: 0.0 for species in MASS_AMU}
    for species in PARENT_SOURCE_RATE_S:
        speed[species] = PRIMARY_CHARACTERISTIC_SPEED_M_S

    reactions = load_neutral_reactions()
    for _ in range(len(MASS_AMU)):
        changed = False
        for reaction_id, parent, products, energy_ev in reactions:
            if speed[parent] <= 0.0:
                continue
            masses_kg = [MASS_AMU[product] * AMU_KG for product in products]
            total_mass_kg = sum(masses_kg)
            parent_mass_kg = MASS_AMU[parent] * AMU_KG
            if not math.isclose(total_mass_kg, parent_mass_kg, rel_tol=1.0e-12):
                raise RuntimeError(f"mass is not conserved in {reaction_id}")
            for product, product_mass_kg in zip(products, masses_kg):
                other_mass_kg = total_mass_kg - product_mass_kg
                recoil_m_s = math.sqrt(
                    2.0
                    * energy_ev
                    * EV_J
                    * other_mass_kg
                    / (product_mass_kg * total_mass_kg)
                )
                candidate = speed[parent] + recoil_m_s
                if candidate > speed[product] + 1.0e-10:
                    speed[product] = candidate
                    changed = True
        if not changed:
            break
    else:
        raise RuntimeError("reaction speed propagation did not converge")

    missing = [species for species, value in speed.items() if value <= 0.0]
    if missing:
        raise RuntimeError(f"no characteristic-speed path for: {', '.join(missing)}")
    return speed


def main() -> None:
    density_m3, sigma_m2, mean_free_path_m = derive_surface_mean_free_path()
    ratio = geometric_expansion_ratio(mean_free_path_m)
    outer_width_m = mean_free_path_m * ratio ** (RADIAL_CELL_COUNT - 1)
    root_cell_width_m = (OUTER_RADIUS_M - NUCLEUS_RADIUS_M) / BLOCK_CELL_COUNT
    maximum_refinement = math.ceil(
        math.log2(root_cell_width_m / mean_free_path_m)
    )
    injection_count, occupancy = reference_injection_count()
    speeds = derive_characteristic_speeds()

    print(f"surface_density_m-3={density_m3:.12e}")
    print(f"effective_H2O_cross_section_m2={sigma_m2:.12e}")
    print(f"surface_mean_free_path_m={mean_free_path_m:.12e}")
    print(f"geometric_expansion_ratio={ratio:.12e}")
    print(f"ideal_outer_cell_width_m={outer_width_m:.12e}")
    print(f"maximum_AMR_refinement_level={maximum_refinement}")
    print(f"reference_injection_particle_number={injection_count}")
    print("parent_surface_occupancy=")
    for species in PARENT_SOURCE_RATE_S:
        print(f"  {species:6s} {occupancy[species]:.6f}")
    print("characteristic_speed_m_s=")
    for species in MASS_AMU:
        # Round upward to a whole metre per second so the value written to the
        # input never falls below the calculated energy/momentum bound.
        print(f"  {species:6s} {math.ceil(speeds[species]):d}")


if __name__ == "__main__":
    main()
