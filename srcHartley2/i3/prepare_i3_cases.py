#!/usr/bin/env python3
"""Reproduce the six compile-time AMPS inputs for the 3I/ATLAS campaign.

The program deliberately performs all case-dependent calculations outside
the C++ application.  AMPS inputs are compile-time configuration files, so a
separate preprocessing/build is required for each concrete case.  This tool
is deterministic and network-free: the repository CSV files contain the
reviewed observation table, solar-longitude correction, reaction parameters,
and fluorescence metadata together with their provenance.
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
BLOCK_CELL_COUNT = 8
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


def _max_refinement_level(surface_resolution_m: float) -> int:
    root_width = (DOMAIN_OUTER_RADIUS_M - NUCLEUS_RADIUS_M) / BLOCK_CELL_COUNT
    return math.ceil(math.log2(root_width / surface_resolution_m))


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


def _reaction_records(reactions: list[Reaction], f107p: float) -> str:
    records = []
    for reaction in reactions:
        rate_1au = reaction.a0_s * f107p**reaction.exponent
        records.append(
            "  ChemistryReaction="
            f"{reaction.reaction_id} | {reaction.reactants} | {reaction.products} | "
            f"HELIOCENTRIC | {_format_number(rate_1au)} | 0 | 0 | 1 | "
            f"{_format_number(reaction.excess_energy_ev)}"
        )
    return "\n".join(records)


def _speed_records(speeds: dict[str, int]) -> str:
    return "\n".join(
        f"  RadialCharacteristicSpeed={species:<3s} | {speeds[species]}"
        for species in MASS_AMU
    )


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


def _brightness_records(case: dict[str, str]) -> str:
    """Select reusable velocity-dependent models with explicit case inputs.

    The Lowell OH source is an integrated energy g factor.  The historical
    case-table photon values used 308.5 nm exactly, so that convention is now
    visible in the input instead of being hidden in a precomputed scalar.
    """

    _validate_observation_geometry(case)

    return "\n".join(
        (
            "  BrightnessMode=ON",
            "  BrightnessImpactParameterGrid=LOG | 2500.0 | 1.0e9 | 200",
            "  BrightnessObservationGeometry="
            + " | ".join(
                case[name]
                for name in (
                    "er_x", "er_y", "er_z", "los_x", "los_y", "los_z",
                    "impact_x", "impact_y", "impact_z",
                )
            ),
            "  BrightnessHeliocentricRadialVelocityKMS="
            f"{case['heliocentric_radial_velocity_km_s']}",
            "  BrightnessVelocityDomainPolicy=ABORT",
            "  BrightnessSolarLymanAlphaFlux1AUPhotonCM2S="
            f"{case['h_lya_photon_flux_1au_ph_cm-2_s-1']}",
            "  BrightnessEmission=OH | OH_A2SIGMA_X2PI_00 | "
            "OH_AX00_SA1988_LOWELL_UNQUENCHED | "
            "CORE_SCHLEICHER_AHEARN1988_LOWELL | THIN | "
            "ISOTROPIC_4PI | NONE | 0 | 0 | NONE | 308.5",
            "  BrightnessEmission=H | H_LYMAN_ALPHA | "
            "H_LYA_IKL2020 | CORE_IKL2020_CASE_SOLAR_FLUX | THIN | "
            "ISOTROPIC_4PI | LINE_CENTER | 9.5e-18 | 0.1 | "
            "DOPPLER_B8KMS_FOSC0.4164_CONSERVATIVE_DIAGNOSTIC | 0",
        )
    )


def render_case(template: str, case: dict[str, str], reactions: list[Reaction]) -> str:
    f107p = float(case["exposure_mean_f107p_sfu"])
    _, limiting_species, surface_mfp = derive_surface_resolution(case)
    speeds = derive_characteristic_speeds(reactions)
    replacements = {
        "CASE_ID": case["case_id"],
        "CAMPAIGN_PHASE": case["campaign_phase"],
        "OBSERVATION_START_UTC": case["observation_start_utc"],
        "OBSERVATION_END_UTC": case["observation_end_utc"],
        "HELIOCENTRIC_DISTANCE_AU": case["heliocentric_distance_au"],
        "Q_H2O": case["q_h2o_s-1"],
        "Q_CO": case["q_co_s-1"],
        "Q_CO2": case["q_co2_s-1"],
        "F107P_SFU": case["exposure_mean_f107p_sfu"],
        "CHEMISTRY_REACTIONS": _reaction_records(reactions, f107p),
        "RANDOM_SEED": str(104729 + int(case["case_id"])),
        "LIMITING_MFP_SPECIES": limiting_species,
        "SURFACE_MFP_M": _format_number(surface_mfp),
        "REFERENCE_INJECTION_COUNT": str(derive_reference_injection_count()),
        "MAX_REFINEMENT_LEVEL": str(_max_refinement_level(surface_mfp)),
        "CHARACTERISTIC_SPEEDS": _speed_records(speeds),
        "BRIGHTNESS_RECORDS": _brightness_records(case),
    }
    rendered = template
    for name, value in replacements.items():
        rendered = rendered.replace(f"@{name}@", value)
    unresolved = sorted({token.split("@", 1)[0] for token in rendered.split("@") [1::2]})
    if "@" in rendered:
        raise RuntimeError("unresolved input-template tokens: " + ", ".join(unresolved))
    return rendered


def generate(repository: Path, output_directory: Path) -> list[Path]:
    i3_directory = repository / "srcHartley2" / "i3"
    cases = _load_csv(i3_directory / "reference" / "i3_cases.csv")
    reactions = _load_reactions(i3_directory / "reference" / "photochemistry_channels.csv")
    template = (repository / "input" / "i3.input").read_text(encoding="utf-8")
    output_directory.mkdir(parents=True, exist_ok=True)
    outputs = []
    for case in cases:
        output = output_directory / f"i3.case{case['case_id']}.input"
        output.write_text(render_case(template, case, reactions), encoding="utf-8")
        outputs.append(output)
    return outputs


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
