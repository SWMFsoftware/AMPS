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
import os
from dataclasses import dataclass
from pathlib import Path


AMU_KG = 1.66053906660e-27
BOLTZMANN_J_K = 1.380649e-23
EV_J = 1.602176634e-19
NUCLEUS_RADIUS_M = 2_500.0
# 3e9 m keeps >=78% of an r^-2 H column at observation points inside 1e9 m
# and 62-78% at the outermost case 4-6 points (0% with the earlier 1e9 m).
DOMAIN_OUTER_RADIUS_M = 3.0e9
SURFACE_TEMPERATURE_K = 175.0
RADIAL_CELL_COUNT = 2_000
# The characteristic speeds below are conservative upper bounds on particle
# speed, so the Courant number limits the fastest particle of each species to
# 0.8 cells per step.  The first campaign (C=0.25) left H, the slowest species
# to converge, moving only ~0.14 cells per step on average.
COURANT_NUMBER = 0.8
# Upper bound on direct-parent speeds: the first campaign's output profiles
# reached bulk+3 sigma = 2.2 km/s (the earlier 400 m/s was exceeded 1.4-3.3x).
DIRECT_PARENT_SPEED_M_S = 2_500.0
TARGET_SURFACE_PARTICLES_PER_CELL = 800.0
# Must match input/i3.conf (MaxMeshRefinmentLevel, BlockCells).  The finest
# dyadic cell is the radial extent / (BlockCells * 2^level); AMPS aborts if the
# requested surface resolution is finer than that.
AMR_MAX_REFINEMENT_LEVEL = 30
AMR_BLOCK_CELLS = 8
MAX_ITERATIONS = 150_000

MASS_AMU = {"H": 1, "H2": 2, "O": 16, "OH": 17, "H2O": 18, "CO": 28, "CO2": 44}

# O I 130.4-nm brightness: three additive core-model components, each emitted
# as its own column group of pic.O.brightness.  Direct triplet pumping has
# two solar line shapes that bracket the unresolved line wings (Killen et
# al. 2022 / Gladstone 1992 versus the SOHO/SUMER atlas); Lyman-beta (Bowen)
# pumping is added to either.  The line-centre optical-depth diagnostic uses
# the 130.217-nm line (f=0.0520) with a 1 km/s Doppler parameter; like the H
# diagnostic it flags, and never attenuates, the optically thin result.
O_I_1304_TAU = "LINE_CENTER | 1.0e-17 | 0.1 | DOPPLER_B1KMS_FOSC0.0520_130.217NM_CONSERVATIVE_DIAGNOSTIC"
O_I_1304_EMISSIONS = (
    ("O_DIRECT_KVB2022", "OI_1304_TRIPLET_DIRECT_KVB2022",
     "O_I_1304_KVB2022_DIRECT_J2", "CORE_KVB2022_GLADSTONE1992_CASE_FISM2_SCALED"),
    ("O_DIRECT_SUMER2001", "OI_1304_TRIPLET_DIRECT_SUMER2001",
     "O_I_1304_SUMER2001_DIRECT_J2", "CORE_SUMER2001_ATLAS_CASE_FISM2_SCALED"),
    ("O_LYB_BOWEN", "OI_1304_TRIPLET_LYB_BOWEN_SUMER2001",
     "O_I_1304_SUMER2001_LYB_BOWEN_J2", "CORE_SUMER2001_LYB_CASE_FISM2"),
)
DIRECT_PARENTS = ("H2O", "CO", "CO2")

# OH loss variants (srcHartley2/CHEMISTRY_VARIANTS.md).  STATE_RESOLVED, the
# default, uses the chemistry.rate.OH_* keys written for every case below.
# SA1982_VH makes AMPS compute one lumped OH -> O + H rate and the OH
# ionization rate from the table file at the case's
# brightness.heliocentric_radial_velocity_km_s.  The generator always writes the
# state-resolved rates and the table path, so either variant can be selected
# later (in the input or with --chemistry-variant) without regenerating; it
# never writes the table-derived rates, so the table stays their only source.
OH_LOSS_VARIANTS = {"state_resolved": "STATE_RESOLVED", "sa1982_vh": "SA1982_VH"}
DEFAULT_OH_LIFETIME_TABLE = Path("srcHartley2/i3/reference/oh_lifetime_sa1982.csv")


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


def finest_surface_cell_m() -> float:
    """Width of the finest dyadic radial cell the compiled AMR tree allows."""

    extent = DOMAIN_OUTER_RADIUS_M - NUCLEUS_RADIUS_M
    return extent / (AMR_BLOCK_CELLS * 2**AMR_MAX_REFINEMENT_LEVEL)


def surface_resolution_request(case: dict[str, str]) -> tuple[str, float, float, bool]:
    """Requested surface cell: the mixture mean free path, floored 1% above
    the finest dyadic cell so the AMR level limit is never exceeded."""

    _, limiting, mean_free_path = derive_surface_resolution(case)
    floor = 1.01 * finest_surface_cell_m()
    return limiting, mean_free_path, max(mean_free_path, floor), mean_free_path < floor


def apply_rate_scaling(
    cases: list[dict[str, str]], scaling: list[dict[str, str]]
) -> list[dict[str, str]]:
    """Multiply the transcribed archive rates of all three direct parents by
    the reviewed per-case factor (kept in a separate file so the archive
    transcription in i3_cases.csv stays unmodified)."""

    factors = {row["case_id"]: float(row["parent_rate_factor"]) for row in scaling}
    if sorted(factors) != sorted(case["case_id"] for case in cases):
        raise ValueError("rate scaling must list every case exactly once")
    scaled = []
    for case in cases:
        factor = factors[case["case_id"]]
        if not (math.isfinite(factor) and factor > 0.0):
            raise ValueError(f"case {case['case_id']} has an invalid rate factor")
        updated = dict(case)
        for column in ("q_h2o_s-1", "q_co_s-1", "q_co2_s-1"):
            updated[column] = f"{float(case[column]) * factor:.3e}"
        scaled.append(updated)
    return scaled


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


def attach_solar_fuv_lines(
    cases: list[dict[str, str]], lines: list[dict[str, str]]
) -> list[dict[str, str]]:
    """Add the per-case O I triplet and Lyman-beta fluxes (derive_solar_fuv_lines.py)."""

    by_case = {row["case_id"]: row for row in lines}
    if set(by_case) != {case["case_id"] for case in cases}:
        raise ValueError("solar_fuv_lines.csv must list exactly the I3 cases")
    merged = []
    for case in cases:
        row = by_case[case["case_id"]]
        for key in ("earth_equivalent_start_utc", "earth_equivalent_end_utc"):
            if row[key] != case[key]:
                raise ValueError(f"case {case['case_id']}: {key} differs from i3_cases.csv")
        merged.append(
            dict(
                case,
                o_i_1304_triplet_photon_flux_1au=row[
                    "o_i_1304_triplet_photon_flux_1au_ph_cm-2_s-1"
                ],
                h_lyb_photon_flux_1au=row[
                    "h_lyb_102p4_102p8nm_photon_flux_1au_ph_cm-2_s-1"
                ],
            )
        )
    return merged


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
    cases: list[dict[str, str]],
    reactions: list[Reaction],
    oh_lifetime_table: str = "../" + DEFAULT_OH_LIFETIME_TABLE.as_posix(),
    oh_loss: str = "state_resolved",
    oh_loss_cases: dict[str, str] | None = None,
) -> str:
    """Render common run settings plus six explicit case overlays.

    oh_lifetime_table is written as given; AMPS resolves a relative path against
    the directory of the runtime input.  oh_loss selects the OH-loss variant of
    the common block and oh_loss_cases maps case IDs to per-case variants (both
    use the lower-case names of OH_LOSS_VARIANTS).
    """

    oh_loss_cases = oh_loss_cases or {}
    for variant in [oh_loss, *oh_loss_cases.values()]:
        if variant not in OH_LOSS_VARIANTS:
            raise ValueError(f"unknown OH-loss variant {variant!r}; "
                             f"known: {', '.join(OH_LOSS_VARIANTS)}")
    unknown_cases = set(oh_loss_cases) - {case["case_id"] for case in cases}
    if unknown_cases:
        raise ValueError("unknown case IDs in the per-case OH-loss selection: "
                         + ", ".join(sorted(unknown_cases)))
    speeds = derive_characteristic_speeds(reactions)
    lines = [
        "! Generated deterministically by srcHartley2/i3/prepare_i3_cases.py.",
        "! Compile once with input/i3.conf, then select one case explicitly:",
        "!   ./amps -i i3.in -case CASE",
        "! Do not copy case values into i3.conf; that file owns layout only.",
        "!",
        "! MAVEN/IUVS: the model is consistent with the MAVEN/IUVS OH (0-0), H Lyman-alpha and O I 130.4-nm",
        "! observations of all six cases with the OH chemistry implemented later, OH_LOSS = SA1982_VH, and",
        "! Q(H2O) fitted to OH. That deck, its build input (BlockCells=12) and the results are in",
        "! srcHartley2/i3/campaigns/qfit_20261123/ (README.md there; O observed for case 5 only).",
        "! This deck keeps the default STATE_RESOLVED OH loss and the source rates of",
        "! reference/i3_rate_scaling.csv; either variant can be selected for every case",
        "! (--chemistry-variant OH_LOSS=SA1982_VH or STATE_RESOLVED).",
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
        "  ! OH loss (srcHartley2/CHEMISTRY_VARIANTS.md): STATE_RESOLVED uses the",
        "  ! per-case rates of the four OH channels; SA1982_VH computes the OH rates",
        "  ! from this table at the case's heliocentric radial velocity.  A relative",
        "  ! path is resolved against the directory of this file.",
        f"  chemistry.variant.OH_LOSS = {OH_LOSS_VARIANTS[oh_loss]}",
        f"  chemistry.table.OH_TAU_SA1982.file = {oh_lifetime_table}",
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
        f"  simulation.max_iterations = {MAX_ITERATIONS}",
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
            "  brightness.impact_parameter_grid = LOG | 2500.0 | "
            f"{_format_number(DOMAIN_OUTER_RADIUS_M)} | 200",
            "  brightness.velocity_domain_policy = ABORT",
            "  brightness.emission.OH = OH | OH_A2SIGMA_X2PI_00 | "
            "OH_AX00_SA1988_LOWELL_UNQUENCHED | "
            "CORE_SCHLEICHER_AHEARN1988_LOWELL | THIN | ISOTROPIC_4PI | "
            "NONE | 0 | 0 | NONE | 308.5",
            "  brightness.emission.H = H | H_LYMAN_ALPHA | H_LYA_IKL2020 | "
            "CORE_IKL2020_CASE_SOLAR_FLUX | THIN | ISOTROPIC_4PI | "
            "LINE_CENTER | 9.5e-18 | 0.1 | "
            "DOPPLER_B8KMS_FOSC0.4164_CONSERVATIVE_DIAGNOSTIC | 0",
        ]
    )
    lines.extend(
        f"  brightness.emission.{key} = O | {transition} | {model} | "
        f"{provenance} | THIN | ISOTROPIC_4PI | {O_I_1304_TAU} | 0"
        for key, transition, model, provenance in O_I_1304_EMISSIONS
    )
    lines.append("  dust.enabled = false")

    for case in cases:
        _validate_observation_geometry(case)
        f107p = float(case["exposure_mean_f107p_sfu"])
        limiting_species, surface_mfp, surface_request, floored = (
            surface_resolution_request(case)
        )
        resolution_note = (
            f"{limiting_species} = {_format_number(surface_mfp)} m; floored at "
            f"the finest {AMR_MAX_REFINEMENT_LEVEL}-level AMR cell."
            if floored
            else f"{limiting_species}."
        )
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
                f"{resolution_note}",
                f"    mesh.surface_resolution_m = {_format_number(surface_request)}",
            ]
        )
        if case["case_id"] in oh_loss_cases:
            lines.append("    chemistry.variant.OH_LOSS = "
                         f"{OH_LOSS_VARIANTS[oh_loss_cases[case['case_id']]]}")
        lines.extend(_runtime_reaction_records(reactions, f107p))
        lines.extend(
            [
                f"    brightness.observation_geometry = {geometry}",
                "    brightness.heliocentric_radial_velocity_km_s = "
                f"{case['heliocentric_radial_velocity_km_s']}",
                "    brightness.solar_lya_flux_1au_photon_cm-2_s-1 = "
                f"{case['h_lya_photon_flux_1au_ph_cm-2_s-1']}",
                "    brightness.solar_oi1304_triplet_flux_1au_photon_cm-2_s-1 = "
                f"{case['o_i_1304_triplet_photon_flux_1au']}",
                "    brightness.solar_lyb_flux_1au_photon_cm-2_s-1 = "
                f"{case['h_lyb_photon_flux_1au']}",
                "  #end case",
            ]
        )

    lines.extend(["#end application", ""])
    return "\n".join(lines)


def generate(
    repository: Path,
    output_directory: Path,
    oh_loss: str = "state_resolved",
    oh_loss_cases: dict[str, str] | None = None,
    oh_lifetime_table: Path | None = None,
) -> list[Path]:
    i3_directory = repository / "srcHartley2" / "i3"
    table = (oh_lifetime_table or repository / DEFAULT_OH_LIFETIME_TABLE).resolve()
    if not table.is_file():
        raise FileNotFoundError(f"OH lifetime table {table} does not exist")
    cases = attach_solar_fuv_lines(
        apply_rate_scaling(
            _load_csv(i3_directory / "reference" / "i3_cases.csv"),
            _load_csv(i3_directory / "reference" / "i3_rate_scaling.csv"),
        ),
        _load_csv(i3_directory / "reference" / "solar_fuv_lines.csv"),
    )
    reactions = _load_reactions(i3_directory / "reference" / "photochemistry_channels.csv")
    output_directory.mkdir(parents=True, exist_ok=True)
    output = output_directory / "i3.in"
    # Written relative to the output directory, where AMPS resolves it.
    table_text = Path(os.path.relpath(table, output_directory.resolve())).as_posix()
    output.write_text(
        render_runtime_input(cases, reactions, table_text, oh_loss, oh_loss_cases),
        encoding="utf-8",
    )
    return [output]


def _parse_case_variant(text: str) -> tuple[str, str]:
    case_id, separator, variant = text.partition("=")
    if not separator or not case_id.strip() or variant.strip() not in OH_LOSS_VARIANTS:
        raise argparse.ArgumentTypeError(
            f"expected CASE=VARIANT with VARIANT one of {', '.join(OH_LOSS_VARIANTS)}, "
            f"got {text!r}")
    return case_id.strip(), variant.strip()


def main() -> None:
    default_repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repository", type=Path, default=default_repository)
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument(
        "--oh-loss", choices=sorted(OH_LOSS_VARIANTS), default="state_resolved",
        help="OH-loss variant of the common block (srcHartley2/CHEMISTRY_VARIANTS.md): "
             "state_resolved (default) uses the four state-resolved OH channels with "
             "the rates of photochemistry_channels.csv; sa1982_vh uses one OH -> O + H "
             "channel and OH ionization with the heliocentric-velocity-dependent "
             "lifetime table. AMPS can also select the variant at run time "
             "(--chemistry-variant OH_LOSS=SA1982_VH).")
    parser.add_argument(
        "--oh-loss-case", type=_parse_case_variant, action="append", default=[],
        metavar="CASE=VARIANT",
        help="per-case OH-loss variant, overriding --oh-loss for that case "
             "(repeatable), e.g. --oh-loss-case 3=sa1982_vh")
    parser.add_argument(
        "--oh-lifetime-table", type=Path, default=None,
        help="OH lifetime table for sa1982_vh (default: "
             f"{DEFAULT_OH_LIFETIME_TABLE}); written into i3.in relative to the "
             "output directory. Another solar-activity case is another table file "
             "with the same format.")
    arguments = parser.parse_args()
    output_directory = arguments.output_directory or arguments.repository / "input"
    for output in generate(
        arguments.repository.resolve(),
        output_directory.resolve(),
        arguments.oh_loss,
        dict(arguments.oh_loss_case),
        arguments.oh_lifetime_table,
    ):
        print(output)


if __name__ == "__main__":
    main()
