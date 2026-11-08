#!/usr/bin/env python3
"""Reproduce the one runtime ``i3.in`` receipt for all six 3I/ATLAS cases.

``input/i3.conf`` fixes compiled species, dimension, block layout, maximum AMR
level, enabled core modules, collision topology, and reaction topology. This
program keeps every case-dependent calculation outside the C++ application
and writes it into named ``case`` sections in one runtime file. The result is
deterministic and network-free.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


AMU_KG = 1.66053906660e-27
BOLTZMANN_J_K = 1.380649e-23
EV_J = 1.602176634e-19
NUCLEUS_RADIUS_M = 2_500.0
DOMAIN_OUTER_RADIUS_M = 1.0e9
SURFACE_TEMPERATURE_K = 175.0
RADIAL_CELL_COUNT = 2_000
COURANT_NUMBER = 0.25
DIRECT_PARENT_SPEED_M_S = 400.0
TARGET_SURFACE_PARTICLES_PER_CELL = 800.0

MASS_AMU = {"H": 1, "H2": 2, "O": 16, "OH": 17, "H2O": 18, "CO": 28, "CO2": 44}
DIRECT_PARENTS = ("H2O", "CO", "CO2")


@dataclass(frozen=True)
class Reaction:
    reaction_id: str
    reactants: str
    products: str
    a0_s: float
    exponent: float
    excess_energy_ev: float
    energy_source: str


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _format_number(value: float) -> str:
    return f"{value:.12g}"


def _outward_spatial_mean_speed(species: str) -> float:
    """Mean positive normal speed of the residence-weighted half Maxwellian."""

    mass_kg = MASS_AMU[species] * AMU_KG
    return math.sqrt(2.0 * BOLTZMANN_J_K * SURFACE_TEMPERATURE_K / (math.pi * mass_kg))


def _parent_collision_cross_section(species0: str, species1: str) -> float:
    """Return the exact parent-pair law used by the generated input at 175 K."""

    pair = tuple(sorted((species0, species1)))
    table = {
        ("H2O", "H2O"): 1.66e-19 * (SURFACE_TEMPERATURE_K / 300.0) ** -0.6,
        ("CO", "H2O"): 3.2e-19,
        ("CO2", "H2O"): 3.4e-19,
        ("CO", "CO"): 3.2e-19,
        ("CO", "CO2"): 3.2e-19,
        ("CO2", "CO2"): 3.4e-19,
    }
    return table[pair]


def derive_surface_resolution(case: dict[str, str]) -> tuple[dict[str, float], str, float]:
    """Calculate each direct parent's surface mixture mean free path.

    Boundary injection is flux weighted, but the particles resident in the
    first cell recover the outgoing half-Maxwellian spatial distribution.
    Thus n_j=Q_j/(4*pi*R^2*<v_n>_j).  The mesh resolves the smallest of
    lambda_i=1/(sqrt(2)*sum_j n_j*sigma_ij).  The sqrt(2) convention matches
    the documented Halley surface-resolution derivation.
    """

    area = 4.0 * math.pi * NUCLEUS_RADIUS_M**2
    rates = {species: float(case[f"q_{species.lower()}_s-1"]) for species in DIRECT_PARENTS}
    density = {
        species: rates[species] / area / _outward_spatial_mean_speed(species)
        for species in DIRECT_PARENTS
    }
    mean_free_path = {
        species: 1.0
        / (
            math.sqrt(2.0)
            * sum(
                density[partner] * _parent_collision_cross_section(species, partner)
                for partner in DIRECT_PARENTS
            )
        )
        for species in DIRECT_PARENTS
    }
    limiting = min(mean_free_path, key=mean_free_path.get)
    return mean_free_path, limiting, mean_free_path[limiting]


def derive_reference_injection_count() -> int:
    """Count needed for at least 800 direct-parent particles in cell one."""

    requirements = [
        TARGET_SURFACE_PARTICLES_PER_CELL
        * COURANT_NUMBER
        * _outward_spatial_mean_speed(species)
        / DIRECT_PARENT_SPEED_M_S
        for species in DIRECT_PARENTS
    ]
    return math.ceil(max(requirements))


def _expanded_products(products: str) -> list[str]:
    if products == "NONE":
        return []
    result: list[str] = []
    for raw in products.split("+"):
        term = raw.strip()
        if "*" in term:
            coefficient_text, species = term.split("*", 1)
            result.extend([species.strip()] * int(coefficient_text))
        else:
            result.append(term)
    return result


def derive_characteristic_speeds(reactions: list[Reaction]) -> dict[str, int]:
    """Propagate conservative recoil-speed bounds through the neutral graph."""

    speed = {species: 0.0 for species in MASS_AMU}
    for species in DIRECT_PARENTS:
        speed[species] = DIRECT_PARENT_SPEED_M_S

    for _ in range(len(MASS_AMU)):
        changed = False
        for reaction in reactions:
            products = _expanded_products(reaction.products)
            if not products or speed[reaction.reactants] <= 0.0:
                continue
            masses = [MASS_AMU[product] * AMU_KG for product in products]
            total_mass = sum(masses)
            parent_mass = MASS_AMU[reaction.reactants] * AMU_KG
            if not math.isclose(total_mass, parent_mass, rel_tol=1.0e-12):
                raise ValueError(f"mass is not conserved by {reaction.reaction_id}")
            for product, product_mass in zip(products, masses):
                other_mass = total_mass - product_mass
                recoil = math.sqrt(
                    2.0
                    * reaction.excess_energy_ev
                    * EV_J
                    * other_mass
                    / (product_mass * total_mass)
                )
                candidate = speed[reaction.reactants] + recoil
                if candidate > speed[product] + 1.0e-10:
                    speed[product] = candidate
                    changed = True
        if not changed:
            break
    else:
        raise RuntimeError("characteristic-speed propagation did not converge")

    missing = [species for species, value in speed.items() if value <= 0.0]
    if missing:
        raise RuntimeError("no source/reaction path for " + ", ".join(missing))
    return {species: math.ceil(value) for species, value in speed.items()}


def _load_reactions(path: Path) -> list[Reaction]:
    reactions = []
    for row in _load_csv(path):
        reactions.append(
            Reaction(
                row["reaction_id"],
                row["reactants"],
                row["products"],
                float(row["a0_s-1"]),
                float(row["f107p_exponent"]),
                float(row["excess_energy_ev"]),
                row["energy_source"],
            )
        )
    return reactions


def _validate_observation_geometry(case: dict[str, str]) -> None:
    """Reject incomplete or internally inconsistent case geometry.

    The generator deliberately does not invent a missing observing direction.
    The reviewed CSV must contain three finite, unit-length vectors.  LOS is
    comet-to-observer, and the positive impact direction is the normalized
    projection of the antisolar radial direction into the sky plane.  These
    checks catch a column-order/sign error before a concrete AMPS input is
    written; the Perl parser and C++ runtime independently recheck the generic
    normalization and LOS/impact orthogonality constraints.
    """

    names = {
        "radial": ("er_x", "er_y", "er_z"),
        "los": ("los_x", "los_y", "los_z"),
        "impact": ("impact_x", "impact_y", "impact_z"),
    }
    vectors: dict[str, tuple[float, float, float]] = {}
    for label, columns in names.items():
        try:
            vector = tuple(float(case[column]) for column in columns)
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"case {case.get('case_id', '?')} has invalid {label} geometry"
            ) from exc
        if not all(math.isfinite(value) for value in vector):
            raise ValueError(
                f"case {case.get('case_id', '?')} has non-finite {label} geometry"
            )
        norm = math.sqrt(sum(value * value for value in vector))
        if not math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=2.0e-12):
            raise ValueError(
                f"case {case.get('case_id', '?')} {label} vector is not unit length"
            )
        vectors[label] = vector

    radial, los, impact = (
        vectors["radial"], vectors["los"], vectors["impact"]
    )
    los_impact = sum(los[index] * impact[index] for index in range(3))
    if not math.isclose(los_impact, 0.0, rel_tol=0.0, abs_tol=2.0e-12):
        raise ValueError(
            f"case {case.get('case_id', '?')} LOS and impact vectors are not orthogonal"
        )

    radial_los = sum(radial[index] * los[index] for index in range(3))
    projected = tuple(
        radial[index] - radial_los * los[index] for index in range(3)
    )
    projected_norm = math.sqrt(sum(value * value for value in projected))
    if not projected_norm > 1.0e-12:
        raise ValueError(
            f"case {case.get('case_id', '?')} has undefined projected antisolar direction"
        )
    projected = tuple(value / projected_norm for value in projected)
    mismatch = math.sqrt(
        sum((impact[index] - projected[index]) ** 2 for index in range(3))
    )
    if mismatch > 2.0e-12:
        raise ValueError(
            f"case {case.get('case_id', '?')} impact vector is not projected antisolar"
        )


def _runtime_reaction_records(reactions: list[Reaction], f107p: float) -> list[str]:
    return [
        f"    chemistry.rate.{reaction.reaction_id}_1au_s-1 = "
        f"{_format_number(reaction.a0_s * f107p**reaction.exponent)}"
        for reaction in reactions
    ]


def render_runtime_input(
    cases: list[dict[str, str]], reactions: list[Reaction]
) -> str:
    """Render common run settings plus six explicit case overlays."""

    speeds = derive_characteristic_speeds(reactions)
    lines = [
        "! Generated deterministically by srcHartley2/i3/prepare_i3_cases.py.",
        "! Compile once with input/i3.conf, then select one case explicitly:",
        "!   ./amps -i i3.in -case CASE",
        "! Do not copy case values into i3.conf; that file owns layout only.",
        "",
        "#begin application: comet",
        "  schema = 1",
        "  geometry = spherical",
        "",
        "  source.H2O.distribution = uniform",
        "  source.H2O.temperature_k = 175.0",
        "  source.H2O.bulk_speed_m_s = 0.0",
        "  source.CO.distribution = uniform",
        "  source.CO.temperature_k = 175.0",
        "  source.CO.bulk_speed_m_s = 0.0",
        "  source.CO2.distribution = uniform",
        "  source.CO2.temperature_k = 175.0",
        "  source.CO2.bulk_speed_m_s = 0.0",
        "",
        "  collision.temperature_floor_k = 175.0",
        "  collision.minimum_sample_particles = 1.0",
        "  chemistry.enabled = true",
        "  chemistry.residual_sampling = systematic",
        "  diagnostics.enabled = true",
        "  cooling.H2O.model = GOMBOSI1986_H2O",
        "",
        f"  nucleus.radius_m = {_format_number(NUCLEUS_RADIUS_M)}",
        f"  mesh.outer_radius_m = {_format_number(DOMAIN_OUTER_RADIUS_M)}",
        f"  mesh.radial_cell_count = {RADIAL_CELL_COUNT}",
        f"  simulation.reference_injection_particles = {derive_reference_injection_count()}",
        f"  simulation.courant_number = {_format_number(COURANT_NUMBER)}",
        "  simulation.inner_boundary_temperature_k = 175.0",
        "  simulation.sampling_length = 10000",
        "  simulation.weight_adaptation_enabled = true",
        "  simulation.weight_adaptation_interval = 10000",
        "  simulation.particles_per_cell.minimum = 800",
        "  simulation.particles_per_cell.maximum = 3200",
        "  simulation.chemical_only_particles_per_cell.minimum = 50",
        "  simulation.chemical_only_particles_per_cell.maximum = 200",
        "  simulation.lock_direct_source_weights = true",
        "  simulation.emergency_buffer_utilization = 0.75",
        "  simulation.allocation_low_watermark = 0.50",
        "  simulation.allocation_reserve_fraction = 0.10",
        "  simulation.max_face_transfer_multiplicity = 1",
        "  simulation.max_photochemical_products_per_reaction = 20",
        "  simulation.load_balance_interval = 50",
        "  simulation.max_iterations = 30000",
        "  simulation.convergence_tolerance = 0.05",
        "  simulation.convergence_windows = 3",
    ]
    lines.extend(
        f"  transport.characteristic_speed.{species}_m_s = {speeds[species]}"
        for species in MASS_AMU
    )
    lines.extend(
        [
            "",
            "  brightness.enabled = true",
            "  brightness.impact_parameter_grid = LOG | 2500.0 | 1.0e9 | 200",
            "  brightness.velocity_domain_policy = ABORT",
            "  brightness.emission.OH = OH | OH_A2SIGMA_X2PI_00 | "
            "OH_AX00_SA1988_LOWELL_UNQUENCHED | "
            "CORE_SCHLEICHER_AHEARN1988_LOWELL | THIN | ISOTROPIC_4PI | "
            "NONE | 0 | 0 | NONE | 308.5",
            "  brightness.emission.H = H | H_LYMAN_ALPHA | H_LYA_IKL2020 | "
            "CORE_IKL2020_CASE_SOLAR_FLUX | THIN | ISOTROPIC_4PI | "
            "LINE_CENTER | 9.5e-18 | 0.1 | "
            "DOPPLER_B8KMS_FOSC0.4164_CONSERVATIVE_DIAGNOSTIC | 0",
            "  dust.enabled = false",
        ]
    )

    for case in cases:
        _validate_observation_geometry(case)
        f107p = float(case["exposure_mean_f107p_sfu"])
        _, limiting_species, surface_mfp = derive_surface_resolution(case)
        geometry = " | ".join(
            case[name]
            for name in (
                "er_x", "er_y", "er_z", "los_x", "los_y", "los_z",
                "impact_x", "impact_y", "impact_z",
            )
        )
        lines.extend(
            [
                "",
                f"  #begin case: {case['case_id']}",
                f"    case.index = {case['case_id']}",
                f"    case.campaign_phase = {case['campaign_phase']}",
                f"    case.observation_start_utc = {case['observation_start_utc']}",
                f"    case.observation_end_utc = {case['observation_end_utc']}",
                f"    case.exposure_mean_f107p_sfu = {case['exposure_mean_f107p_sfu']}",
                f"    heliocentric_distance_au = {case['heliocentric_distance_au']}",
                f"    random_seed = {104729 + int(case['case_id'])}",
                f"    source.H2O.rate_s-1 = {case['q_h2o_s-1']}",
                f"    source.CO.rate_s-1 = {case['q_co_s-1']}",
                f"    source.CO2.rate_s-1 = {case['q_co2_s-1']}",
                "    ! Minimum mixture mean free path at the 175-K surface: "
                f"{limiting_species}.",
                f"    mesh.surface_resolution_m = {_format_number(surface_mfp)}",
            ]
        )
        lines.extend(_runtime_reaction_records(reactions, f107p))
        lines.extend(
            [
                f"    brightness.observation_geometry = {geometry}",
                "    brightness.heliocentric_radial_velocity_km_s = "
                f"{case['heliocentric_radial_velocity_km_s']}",
                "    brightness.solar_lya_flux_1au_photon_cm-2_s-1 = "
                f"{case['h_lya_photon_flux_1au_ph_cm-2_s-1']}",
                "  #end case",
            ]
        )

    lines.extend(["#end application", ""])
    return "\n".join(lines)


def generate(repository: Path, output_directory: Path) -> list[Path]:
    i3_directory = repository / "srcHartley2" / "i3"
    cases = _load_csv(i3_directory / "reference" / "i3_cases.csv")
    reactions = _load_reactions(i3_directory / "reference" / "photochemistry_channels.csv")
    output_directory.mkdir(parents=True, exist_ok=True)
    output = output_directory / "i3.in"
    output.write_text(render_runtime_input(cases, reactions), encoding="utf-8")
    return [output]


def main() -> None:
    default_repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=default_repository)
    parser.add_argument("--output-directory", type=Path)
    arguments = parser.parse_args()
    output_directory = arguments.output_directory or arguments.repository / "input"
    for output in generate(arguments.repository.resolve(), output_directory.resolve()):
        print(output)


if __name__ == "__main__":
    main()
