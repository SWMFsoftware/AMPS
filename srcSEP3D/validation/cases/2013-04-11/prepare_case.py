#!/usr/bin/env python3
"""Prepare OV3D01 input and checksum-owned evidence without hidden defaults.

The repository-owned parameter manifest separates publication values from
event-fit and numerical values.  This tool enforces that boundary.  It can:

* report every still-missing field in an event-fit JSON file;
* render a complete srcSEP3D input only after all required values, coordinate
  provenance, data provenance, and convergence provenance are present; and
* copy already calibrated model/reference series into an immutable Phase-V
  evidence bundle with exact SHA-256 declarations.

It deliberately does not download OMNI or particle data, choose quality flags,
fit SWCME, or invent a phase-space-density-to-birth-rate conversion.  Those
operations require scientific review and remain visible inputs to this tool.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_PARAMETERS = HERE / "parameters.json"
DEFAULT_BASE_INPUT = ROOT / "examples" / "sep3d_analytic_parker.in"
AU_M = 149597870700.0
RSUN_M = 695700000.0
SOURCE_RADIUS_M = 20.0 * RSUN_M
ELEMENTARY_CHARGE_C = 1.602176634e-19


class PreparationError(RuntimeError):
    """A requested artifact would be incomplete or physically ambiguous."""


def read_object(path: Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PreparationError(f"cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise PreparationError(f"JSON root must be an object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise PreparationError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def nested(value: Mapping[str, Any], dotted: str) -> Any:
    current: Any = value
    for component in dotted.split("."):
        if not isinstance(current, Mapping) or component not in current:
            return None
        current = current[component]
    return current


# Every field below is required because it either controls event physics,
# numerical convergence, or provenance.  Nullable values in the distributed
# template make incompleteness visible to JSON tooling and to ``status``.
REQUIRED_TEXT = (
    "coordinate_transform.name",
    "coordinate_transform.provenance",
    "ambient.averaging_start_utc",
    "ambient.averaging_end_utc",
    "ambient.source_file",
    "ambient.source_sha256",
    "swcme.thermodynamic_closure",
    "swcme.fit_provenance",
    "source.normalization_provenance",
    "mesh.convergence_provenance",
    "population_control.mode",
    "population_control.convergence_provenance",
    "turbulence.provenance",
    "run.provenance",
)

REQUIRED_POSITIVE = (
    "ambient.number_density_at_one_au_m3",
    "ambient.magnetic_field_magnitude_at_one_au_t",
    "ambient.proton_temperature_k",
    "swcme.solar_rotation_rate_rad_per_s",
    "swcme.adiabatic_index",
    "swcme.electron_temperature_k",
    "swcme.alpha_temperature_k",
    "swcme.sheath_thickness_at_one_au_m",
    "swcme.ejecta_thickness_at_one_au_m",
    "swcme.shock_smoothing_width_at_one_au_m",
    "swcme.leading_edge_smoothing_width_at_one_au_m",
    "swcme.trailing_edge_smoothing_width_at_one_au_m",
    "swcme.sheath_ramp_power",
    "swcme.leading_edge_speed_factor",
    "swcme.ejecta_density_factor",
    "swcme.ejecta_speed_factor",
    "swcme.valid_until_s",
    "source.physical_particle_rate_per_s",
    "source.maximum_energy_j",
    "source.reference_energy_mev",
    "source.samples_per_step_per_species",
    "source.relative_source_weight_per_area",
    "mesh.global_cell_size_m",
    "mesh.minimum_cell_size_m",
    "mesh.maximum_level",
    "mesh.memory_budget_bytes",
    "mesh.solar_surface_cell_size_m",
    "mesh.solar_transition_outer_radius_m",
    "mesh.refinement_tube_radius_at_reference_m",
    "mesh.refinement_tube_center_cell_size_m",
    "mesh.active_tube_radius_at_reference_m",
    "mesh.active_tube_buffer_blocks",
    "mesh.observer_collection_radius_m",
    "population_control.cadence_steps",
    "turbulence.delta_b_over_b",
    "turbulence.correlation_length_m",
    "run.time_step_s",
    "run.duration_s",
    "run.campaign_seed",
    "run.observer_cadence_s",
    "run.observer_energy_bins",
)

REQUIRED_NONNEGATIVE = (
    "swcme.alpha_to_proton_ratio",
    "swcme.drag_coefficient_per_km",
    "swcme.valid_from_s",
    "source.injection_efficiency",
    "population_control.minimum_particles_per_cell_per_species",
    "population_control.target_particles_per_cell_per_species",
    "population_control.maximum_particles_per_cell_per_species",
)

REQUIRED_INTEGER = (
    "source.samples_per_step_per_species",
    "mesh.maximum_level",
    "mesh.memory_budget_bytes",
    "mesh.active_tube_buffer_blocks",
    "population_control.minimum_particles_per_cell_per_species",
    "population_control.target_particles_per_cell_per_species",
    "population_control.maximum_particles_per_cell_per_species",
    "population_control.cadence_steps",
    "run.campaign_seed",
    "run.observer_energy_bins",
)


def finite_number(value: Any) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)))


def missing_fields(fit: Mapping[str, Any]) -> List[str]:
    missing: List[str] = []
    if fit.get("schema") != "srcsep3d-ov3d01-event-fit-v1":
        missing.append("schema=srcsep3d-ov3d01-event-fit-v1")
    if str(fit.get("case_id", "")).upper() != "OV3D01":
        missing.append("case_id=OV3D01")
    for name in REQUIRED_TEXT:
        value = nested(fit, name)
        if not isinstance(value, str) or not value.strip():
            missing.append(name)
    for name in REQUIRED_POSITIVE:
        value = nested(fit, name)
        if not finite_number(value) or float(value) <= 0.0:
            missing.append(name)
    for name in REQUIRED_NONNEGATIVE:
        value = nested(fit, name)
        if not finite_number(value) or float(value) < 0.0:
            missing.append(name)
    if nested(fit, "ambient.magnetic_polarity") not in (-1, 1):
        missing.append("ambient.magnetic_polarity")
    matrix = nested(fit, "coordinate_transform.model_from_event_frozen_carrington")
    if (not isinstance(matrix, list) or len(matrix) != 3 or
            any(not isinstance(row, list) or len(row) != 3 for row in matrix)):
        missing.append("coordinate_transform.model_from_event_frozen_carrington")
    elif any(not finite_number(item) for row in matrix for item in row):
        missing.append("coordinate_transform.model_from_event_frozen_carrington")
    return sorted(set(missing))


def determinant(matrix: Sequence[Sequence[float]]) -> float:
    return (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))


def validate_fit(fit: Mapping[str, Any]) -> None:
    missing = missing_fields(fit)
    if missing:
        raise PreparationError("event fit is incomplete: " + ", ".join(missing))

    matrix = [[float(item) for item in row] for row in nested(
        fit, "coordinate_transform.model_from_event_frozen_carrington")]
    tolerance = 1.0e-12
    for i in range(3):
        for j in range(3):
            dot = sum(matrix[k][i] * matrix[k][j] for k in range(3))
            expected = 1.0 if i == j else 0.0
            if abs(dot - expected) > tolerance:
                raise PreparationError("coordinate transform is not orthonormal")
    if abs(determinant(matrix) - 1.0) > tolerance:
        raise PreparationError("coordinate transform must be a proper rotation")
    # srcSEP3D's current analytic Parker implementation fixes the solar
    # rotation axis to model +Z.  A general tilted transform would require a
    # corresponding generalized spiral/mesh parameterization; reject it here.
    if (abs(matrix[0][2]) > tolerance or abs(matrix[1][2]) > tolerance or
            abs(matrix[2][0]) > tolerance or abs(matrix[2][1]) > tolerance or
            abs(matrix[2][2] - 1.0) > tolerance):
        raise PreparationError(
            "coordinate transform must preserve solar north/model +Z")

    if nested(fit, "ambient.magnetic_polarity") not in (-1, 1):
        raise PreparationError("ambient.magnetic_polarity must be -1 or +1")
    for name in REQUIRED_INTEGER:
        value = nested(fit, name)
        if not finite_number(value) or int(value) != float(value):
            raise PreparationError(f"{name} must be an integer")
    closure = str(nested(fit, "swcme.thermodynamic_closure"))
    if closure not in ("proton-only", "multi-species"):
        raise PreparationError(
            "swcme.thermodynamic_closure must be proton-only or multi-species")
    efficiency = float(nested(fit, "source.injection_efficiency"))
    if not (0.0 < efficiency <= 1.0):
        raise PreparationError("source.injection_efficiency must be in (0,1]")
    if float(nested(fit, "swcme.adiabatic_index")) <= 1.0:
        raise PreparationError("swcme.adiabatic_index must exceed one")
    if float(nested(fit, "swcme.valid_until_s")) <= float(
            nested(fit, "swcme.valid_from_s")):
        raise PreparationError("SWCME valid interval is not ordered")
    if float(nested(fit, "run.duration_s")) < float(
            nested(fit, "swcme.valid_until_s")):
        raise PreparationError("run duration ends before the SWCME valid interval")
    minimum_energy_j = 1.0e4 * ELEMENTARY_CHARGE_C
    maximum_energy_j = float(nested(fit, "source.maximum_energy_j"))
    reference_energy_j = float(nested(fit, "source.reference_energy_mev")) * 1.0e6 * ELEMENTARY_CHARGE_C
    if not minimum_energy_j < reference_energy_j < maximum_energy_j:
        raise PreparationError("source reference energy must lie inside its energy interval")

    samples = int(nested(fit, "source.samples_per_step_per_species"))
    # The distributed base deck uses 13 theta intervals and 24 phi points.
    # Requiring the full vertex-count upper bound guarantees that exact patch
    # allocation cannot fail merely because more surface nodes become active.
    if samples < 13 * 24:
        raise PreparationError("source samples must be at least 13*24")
    if float(nested(fit, "mesh.active_tube_radius_at_reference_m")) < float(
            nested(fit, "mesh.refinement_tube_radius_at_reference_m")):
        raise PreparationError("active tube is narrower than refinement tube")
    if float(nested(fit, "mesh.minimum_cell_size_m")) > float(
            nested(fit, "mesh.global_cell_size_m")):
        raise PreparationError("minimum cell size exceeds global cell size")

    mode = str(nested(fit, "population_control.mode"))
    low = int(nested(fit, "population_control.minimum_particles_per_cell_per_species"))
    target = int(nested(fit, "population_control.target_particles_per_cell_per_species"))
    high = int(nested(fit, "population_control.maximum_particles_per_cell_per_species"))
    if mode == "off":
        if (low, target, high) != (0, 0, 0):
            raise PreparationError("disabled population control requires zero limits")
    elif mode == "split-merge":
        if not (2 <= low <= target <= high):
            raise PreparationError("population limits must satisfy 2<=min<=target<=max")
    else:
        raise PreparationError("population_control.mode must be off or split-merge")

    declared_hash = str(nested(fit, "ambient.source_sha256"))
    if len(declared_hash) != 64 or any(c not in "0123456789abcdefABCDEF" for c in declared_hash):
        raise PreparationError("ambient.source_sha256 is not a complete SHA-256")


def published(parameters: Mapping[str, Any], name: str) -> float:
    for item in parameters.get("published_parameters", []):
        if isinstance(item, Mapping) and item.get("name") == name:
            value = item.get("value")
            if finite_number(value):
                return float(value)
    raise PreparationError(f"published parameter is absent or nonnumeric: {name}")


def unit_from_lon_lat(longitude_deg: float, latitude_deg: float) -> List[float]:
    longitude = math.radians(longitude_deg)
    latitude = math.radians(latitude_deg)
    cosine = math.cos(latitude)
    return [cosine * math.cos(longitude), cosine * math.sin(longitude),
            math.sin(latitude)]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    return [sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)]


def number(value: float) -> str:
    return format(float(value), ".17g")


def render_assignments(fit: Mapping[str, Any], parameters: Mapping[str, Any]) -> Dict[Tuple[str, str], str]:
    matrix = [[float(item) for item in row] for row in nested(
        fit, "coordinate_transform.model_from_event_frozen_carrington")]
    earth_longitude_deg = published(
        parameters, "earth_carrington_longitude_deg")
    earth_latitude_deg = published(
        parameters, "earth_carrington_latitude_deg")
    earth = matvec(matrix, unit_from_lon_lat(
        earth_longitude_deg, earth_latitude_deg))
    cme = matvec(matrix, unit_from_lon_lat(
        published(parameters, "cme_carrington_longitude_deg"),
        published(parameters, "cme_carrington_latitude_deg")))

    wind = published(parameters, "pre_event_earth_wind_speed_m_per_s")
    omega = float(nested(fit, "swcme.solar_rotation_rate_rad_per_s"))
    # The SWCME field contains the finite source-surface correction
    # B_phi/B_r=-Omega*(r-r0)*sin(theta)/V. Its field-line integral is not the
    # Archimedean Omega*(r-r0)/V angle. Derive the source longitude from the
    # *actual fitted* Omega and reviewed wind speed every time the case is
    # rendered, so a changed event fit cannot leave a stale frozen footpoint.
    parker_winding_rad = omega / wind * (
        (AU_M - SOURCE_RADIUS_M) -
        SOURCE_RADIUS_M * math.log(AU_M / SOURCE_RADIUS_M))
    parker_event_longitude_deg = (
        earth_longitude_deg + math.degrees(parker_winding_rad)) % 360.0
    parker = matvec(matrix, unit_from_lon_lat(
        parker_event_longitude_deg, earth_latitude_deg))
    parker_longitude = math.atan2(parker[1], parker[0])
    if parker_longitude < 0.0:
        parker_longitude += 2.0 * math.pi
    parker_colatitude = math.acos(max(-1.0, min(1.0, parker[2])))

    total_field = float(nested(fit, "ambient.magnetic_field_magnitude_at_one_au_t"))
    winding = omega * (AU_M - SOURCE_RADIUS_M) * math.sin(parker_colatitude) / wind
    radial_field = total_field / math.sqrt(1.0 + winding * winding)
    time_step = float(nested(fit, "run.time_step_s"))
    samples = int(nested(fit, "source.samples_per_step_per_species"))
    physical_rate = float(nested(fit, "source.physical_particle_rate_per_s"))
    efficiency = float(nested(fit, "source.injection_efficiency"))
    relative_weight = float(nested(fit, "source.relative_source_weight_per_area"))
    macro_weight = physical_rate * efficiency * relative_weight * time_step / samples
    maximum_steps = int(math.ceil(float(nested(fit, "run.duration_s")) / time_step))
    maximum_energy_j = float(nested(fit, "source.maximum_energy_j"))
    maximum_energy_mev = maximum_energy_j / (1.0e6 * ELEMENTARY_CHARGE_C)

    assignments: Dict[Tuple[str, str], str] = {
        ("run", "time_step_s"): number(time_step),
        ("run", "maximum_time_steps"): str(maximum_steps),
        ("run", "campaign_seed"): str(int(nested(fit, "run.campaign_seed"))),
        ("domain", "coordinate_frame"): str(nested(fit, "coordinate_transform.name")),
        ("parker_spiral", "start_mode"): "explicit",
        ("parker_spiral", "initial_x_m"): number(SOURCE_RADIUS_M * parker[0]),
        ("parker_spiral", "initial_y_m"): number(SOURCE_RADIUS_M * parker[1]),
        ("parker_spiral", "initial_z_m"): number(SOURCE_RADIUS_M * parker[2]),
        ("mesh", "global_cell_size_m"): number(nested(fit, "mesh.global_cell_size_m")),
        ("mesh", "minimum_cell_size_m"): number(nested(fit, "mesh.minimum_cell_size_m")),
        ("mesh", "maximum_level"): str(int(nested(fit, "mesh.maximum_level"))),
        ("mesh", "memory_budget_bytes"): str(int(nested(fit, "mesh.memory_budget_bytes"))),
        ("mesh.solar", "surface_cell_size_m"): number(nested(fit, "mesh.solar_surface_cell_size_m")),
        ("mesh.solar", "transition_outer_radius_m"): number(nested(fit, "mesh.solar_transition_outer_radius_m")),
        ("mesh.tube", "source_longitude_rad"): number(parker_longitude),
        ("mesh.tube", "source_colatitude_rad"): number(parker_colatitude),
        ("mesh.tube", "radius_at_reference_m"): number(nested(fit, "mesh.refinement_tube_radius_at_reference_m")),
        ("mesh.tube", "center_cell_size_m"): number(nested(fit, "mesh.refinement_tube_center_cell_size_m")),
        ("mesh.active_region", "mode"): "parker-tube",
        ("mesh.active_region", "radius_at_reference_m"): number(nested(fit, "mesh.active_tube_radius_at_reference_m")),
        ("mesh.active_region", "buffer_blocks"): str(int(nested(fit, "mesh.active_tube_buffer_blocks"))),
        ("background.parker", "radial_field_at_reference_t"): number(radial_field),
        ("background.parker", "solar_rotation_rate_rad_per_s"): number(omega),
        ("background.parker", "solar_wind_speed_m_per_s"): number(wind),
        ("background.parker", "magnetic_polarity"): str(int(nested(fit, "ambient.magnetic_polarity"))),
        ("background.parker", "number_density_at_one_au_m3"): number(nested(fit, "ambient.number_density_at_one_au_m3")),
        ("background.parker", "temperature_k"): number(nested(fit, "ambient.proton_temperature_k")),
        ("turbulence", "delta_b_over_b"): number(nested(fit, "turbulence.delta_b_over_b")),
        ("turbulence", "correlation_length_m"): number(nested(fit, "turbulence.correlation_length_m")),
        ("transport", "mean_free_path_model"): "radial-rigidity-power-law",
        ("transport", "mean_free_path_reference_m"): number(published(parameters, "mean_free_path_reference_m")),
        ("transport", "mean_free_path_reference_radius_m"): number(published(parameters, "mean_free_path_reference_radius_m")),
        ("transport", "mean_free_path_reference_rigidity_v"): number(published(parameters, "mean_free_path_reference_rigidity_v")),
        ("transport", "mean_free_path_radial_exponent"): number(published(parameters, "mean_free_path_radial_exponent")),
        ("transport", "mean_free_path_rigidity_exponent"): number(published(parameters, "mean_free_path_rigidity_exponent")),
        ("source", "physical_particle_rate_per_s"): number(physical_rate),
        ("source", "injection_efficiency"): number(efficiency),
        ("source", "minimum_energy_j"): number(published(parameters, "seed_injection_energy_j")),
        ("source", "maximum_energy_j"): number(maximum_energy_j),
        ("source", "spectrum_model"): "fixed-phase-space-power-law",
        ("source", "phase_space_power_index"): "5",
        ("source", "samples_per_step"): str(samples),
        ("species", "macroparticle_weight"): number(macro_weight),
        ("observer.earth", "position_x_m"): number(AU_M * earth[0]),
        ("observer.earth", "position_y_m"): number(AU_M * earth[1]),
        ("observer.earth", "position_z_m"): number(AU_M * earth[2]),
        ("observer.earth", "collection_radius_m"): number(nested(fit, "mesh.observer_collection_radius_m")),
        ("observer.earth", "cadence_s"): number(nested(fit, "run.observer_cadence_s")),
        ("observer.earth", "energy_bins"): str(int(nested(fit, "run.observer_energy_bins"))),
        ("observer.earth", "minimum_energy_j"): number(published(parameters, "seed_injection_energy_j")),
        ("observer.earth", "maximum_energy_j"): number(maximum_energy_j),
        ("population_control", "mode"): str(nested(fit, "population_control.mode")),
        ("population_control", "minimum_particles_per_cell_per_species"): str(int(nested(fit, "population_control.minimum_particles_per_cell_per_species"))),
        ("population_control", "target_particles_per_cell_per_species"): str(int(nested(fit, "population_control.target_particles_per_cell_per_species"))),
        ("population_control", "maximum_particles_per_cell_per_species"): str(int(nested(fit, "population_control.maximum_particles_per_cell_per_species"))),
        ("population_control", "cadence_steps"): str(int(nested(fit, "population_control.cadence_steps"))),
        ("swcme", "ambient.wind_speed"): number(wind / 1000.0) + " km/s",
        ("swcme", "ambient.density_1au"): number(float(nested(fit, "ambient.number_density_at_one_au_m3")) / 1.0e6) + " cm^-3",
        ("swcme", "ambient.magnetic_field_1au"): number(total_field / 1.0e-9) + " nT",
        ("swcme", "ambient.proton_temperature"): number(nested(fit, "ambient.proton_temperature_k")) + " K",
        ("swcme", "ambient.adiabatic_index"): number(nested(fit, "swcme.adiabatic_index")),
        ("swcme", "ambient.alpha_to_proton_ratio"): number(nested(fit, "swcme.alpha_to_proton_ratio")),
        ("swcme", "ambient.electron_temperature"): number(nested(fit, "swcme.electron_temperature_k")) + " K",
        ("swcme", "ambient.alpha_temperature"): number(nested(fit, "swcme.alpha_temperature_k")) + " K",
        ("swcme", "ambient.thermodynamic_closure"): str(nested(fit, "swcme.thermodynamic_closure")).replace("-", "_"),
        ("swcme", "parker.radial_polarity"): str(int(nested(fit, "ambient.magnetic_polarity"))),
        ("swcme", "parker.sin_theta"): number(math.sin(parker_colatitude)),
        ("swcme", "parker.solar_rotation_rate_rad_per_s"): number(omega),
        ("swcme", "cme.launch_speed"): number(published(parameters, "cme_speed_m_per_s") / 1000.0) + " km/s",
        ("swcme", "cme.drag_coefficient"): number(nested(fit, "swcme.drag_coefficient_per_km")) + " 1/km",
        ("swcme", "geometry.cme_direction_x"): number(cme[0]),
        ("swcme", "geometry.cme_direction_y"): number(cme[1]),
        ("swcme", "geometry.cme_direction_z"): number(cme[2]),
        ("swcme", "geometry.sheath_thickness_1au"): number(float(nested(fit, "swcme.sheath_thickness_at_one_au_m")) / AU_M) + " AU",
        ("swcme", "geometry.ejecta_thickness_1au"): number(float(nested(fit, "swcme.ejecta_thickness_at_one_au_m")) / AU_M) + " AU",
        ("swcme", "smoothing.shock_width_1au"): number(float(nested(fit, "swcme.shock_smoothing_width_at_one_au_m")) / AU_M) + " AU",
        ("swcme", "smoothing.leading_edge_width_1au"): number(float(nested(fit, "swcme.leading_edge_smoothing_width_at_one_au_m")) / AU_M) + " AU",
        ("swcme", "smoothing.trailing_edge_width_1au"): number(float(nested(fit, "swcme.trailing_edge_smoothing_width_at_one_au_m")) / AU_M) + " AU",
        ("swcme", "shock.relative_source_weight_per_area"): number(relative_weight),
        ("swcme", "sheath.ramp_power"): number(nested(fit, "swcme.sheath_ramp_power")),
        ("swcme", "sheath.leading_edge_speed_factor"): number(nested(fit, "swcme.leading_edge_speed_factor")),
        ("swcme", "ejecta.density_factor"): number(nested(fit, "swcme.ejecta_density_factor")),
        ("swcme", "ejecta.speed_factor"): number(nested(fit, "swcme.ejecta_speed_factor")),
        ("swcme", "event.valid_from"): number(nested(fit, "swcme.valid_from_s")) + " s",
        ("swcme", "event.valid_until"): number(nested(fit, "swcme.valid_until_s")) + " s",
        ("swcme", "source.energy_min"): "0.01 MeV",
        ("swcme", "source.energy_max"): number(maximum_energy_mev) + " MeV",
        ("swcme", "source.reference_energy"): number(nested(fit, "source.reference_energy_mev")) + " MeV",
        ("swcme", "source.injection_efficiency"): number(efficiency),
    }
    return assignments


def rewrite_deck(base: str, assignments: Mapping[Tuple[str, str], str]) -> str:
    output: List[str] = []
    section = ""
    seen: Dict[Tuple[str, str], int] = {key: 0 for key in assignments}
    drop_section = False
    for raw in base.splitlines():
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            drop_section = section == "observer.inner"
            if drop_section:
                continue
        if drop_section:
            continue
        if stripped and not stripped.startswith("#") and "=" in raw:
            key = raw.split("=", 1)[0].strip().lower()
            address = (section, key)
            if address in assignments:
                seen[address] += 1
                raw = f"{key} = {assignments[address]}"
        output.append(raw)
    invalid = [f"[{section}].{key} count={seen[(section, key)]}"
               for section, key in sorted(seen) if seen[(section, key)] != 1]
    if invalid:
        raise PreparationError("base input does not contain exact replacement targets: " +
                               ", ".join(invalid))
    return "\n".join(output).rstrip() + "\n"


def write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except OSError as error:
        raise PreparationError(f"cannot publish {path}: {error}") from error


def render_input(fit_path: Path, parameters_path: Path,
                 base_path: Path, output_path: Path) -> None:
    fit = read_object(fit_path)
    validate_fit(fit)
    parameters = read_object(parameters_path)
    if parameters.get("case_id") != "OV3D01":
        raise PreparationError("parameter manifest is not OV3D01")
    try:
        base = base_path.read_text(encoding="utf-8")
    except OSError as error:
        raise PreparationError(f"cannot read base input {base_path}: {error}") from error
    rendered = rewrite_deck(base, render_assignments(fit, parameters))
    header = ("# OV3D01 generated input; do not edit derived values by hand\n"
              f"# event_fit_sha256 = {sha256(fit_path)}\n"
              f"# parameter_manifest_sha256 = {sha256(parameters_path)}\n")
    write_atomic(output_path, header + rendered)


def verify_rendered_input(executable: Path, input_path: Path,
                          timeout_seconds: float) -> None:
    executable = executable.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise PreparationError(
            f"verification executable is absent or not executable: {executable}")
    try:
        completed = subprocess.run(
            [str(executable), "--input", str(input_path.resolve()), "--dry-run"],
            cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PreparationError(f"rendered-input verification failed to run: {error}") from error
    if completed.returncode != 0:
        raise PreparationError(
            "rendered input failed executable dry-run verification: " +
            completed.stdout[-4000:])
    print(completed.stdout, end="")


def validate_series(path: Path) -> None:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["coordinate_si", "value_si"]:
                raise PreparationError(
                    f"{path} header must be coordinate_si,value_si")
            previous = -math.inf
            rows = 0
            for line, row in enumerate(reader, start=2):
                try:
                    coordinate = float(row["coordinate_si"])
                    value = float(row["value_si"])
                except (TypeError, ValueError) as error:
                    raise PreparationError(
                        f"nonnumeric row {line} in {path}") from error
                if (not math.isfinite(coordinate) or coordinate <= previous or
                        not math.isfinite(value) or value <= 0.0):
                    raise PreparationError(f"invalid row {line} in {path}")
                previous = coordinate
                rows += 1
    except OSError as error:
        raise PreparationError(f"cannot read series {path}: {error}") from error
    if rows < 2:
        raise PreparationError(f"series needs at least two rows: {path}")


def build_evidence(arguments: argparse.Namespace) -> None:
    model = arguments.model_csv.resolve()
    reference = arguments.reference_csv.resolve()
    validate_series(model)
    validate_series(reference)
    directory = arguments.output_dir.resolve()
    targets = [directory / "model.csv", directory / "reference.csv",
               directory / "manifest.json"]
    if not arguments.force and any(path.exists() for path in targets):
        raise PreparationError(
            "evidence output exists; use --force only after reviewing replacement")
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(model, targets[0])
    shutil.copyfile(reference, targets[1])
    manifest = {
        "schema": "srcsep3d-scientific-evidence-v1",
        "case_id": "OV3D01",
        "coordinate_name": arguments.coordinate_name,
        "coordinate_units": arguments.coordinate_units,
        "value_name": arguments.value_name,
        "value_units": arguments.value_units,
        "model_csv": "model.csv",
        "model_sha256": sha256(targets[0]),
        "reference_csv": "reference.csv",
        "reference_sha256": sha256(targets[1]),
        "provenance": {
            "model": arguments.model_provenance,
            "reference": arguments.reference_provenance,
        },
        "exact_pairs": [],
    }
    write_atomic(targets[2], json.dumps(
        manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subcommands = result.add_subparsers(dest="command", required=True)
    status = subcommands.add_parser("status", help="list incomplete fit fields")
    status.add_argument("--fit", type=Path, required=True)

    render = subcommands.add_parser(
        "render-input", help="render a complete input after fit validation")
    render.add_argument("--fit", type=Path, required=True)
    render.add_argument("--parameters", type=Path, default=DEFAULT_PARAMETERS)
    render.add_argument("--base-input", type=Path, default=DEFAULT_BASE_INPUT)
    render.add_argument("--output", type=Path, required=True)
    render.add_argument(
        "--verify-executable", type=Path,
        help="run the rendered deck through this linked amps --dry-run")
    render.add_argument("--verify-timeout", type=float, default=300.0)

    evidence = subcommands.add_parser(
        "build-evidence", help="build a checksum-owned OV3D01 series bundle")
    evidence.add_argument("--model-csv", type=Path, required=True)
    evidence.add_argument("--reference-csv", type=Path, required=True)
    evidence.add_argument("--output-dir", type=Path, required=True)
    evidence.add_argument("--coordinate-name", required=True)
    evidence.add_argument("--coordinate-units", required=True)
    evidence.add_argument("--value-name", required=True)
    evidence.add_argument("--value-units", required=True)
    evidence.add_argument("--model-provenance", required=True)
    evidence.add_argument("--reference-provenance", required=True)
    evidence.add_argument("--force", action="store_true")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    if arguments.command == "status":
        fit = read_object(arguments.fit)
        missing = missing_fields(fit)
        if missing:
            print("OV3D01 event fit is incomplete:")
            for name in missing:
                print(f"- {name}")
            return 1
        validate_fit(fit)
        print("OV3D01 event fit is structurally complete")
        return 0
    if arguments.command == "render-input":
        render_input(arguments.fit.resolve(), arguments.parameters.resolve(),
                     arguments.base_input.resolve(), arguments.output.resolve())
        if arguments.verify_executable is not None:
            verify_rendered_input(arguments.verify_executable,
                                  arguments.output, arguments.verify_timeout)
        print(arguments.output.resolve())
        return 0
    if arguments.command == "build-evidence":
        build_evidence(arguments)
        print((arguments.output_dir.resolve() / "manifest.json"))
        return 0
    raise PreparationError(f"unknown command: {arguments.command}")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PreparationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
