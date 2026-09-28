#!/usr/bin/env python3
"""Generate a complete Stage-0 parser fixture from the schema-5 registry.

The fixture covers every key in the nine physical sections materialized by
Stage 0.  Values are deliberately an analytic-verification configuration, not
event defaults.  Tests mutate this known-valid deck one authority at a time.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "src" / "schema5_registry.inc"
OUTPUT = ROOT / "test" / "data" / "schema5_stage0.in"
SECTIONS = (
    "run", "domain", "solar_rotation", "pfss", "current_sheet",
    "closed_field_plasma", "open_closed_interface", "solar_wind", "plasma_eos",
)


def write_text_lf(path: Path, content: str) -> None:
    """Publish deterministic UTF-8/LF text on Python 3.7 and newer.

    ``Path.write_text(newline=...)`` is unavailable on Python 3.8, whereas
    ``Path.open`` supports the newline contract required by generated files.
    """

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(content)


OVERRIDES = {
    "run.intent": "analytic-verification",
    "run.transport": "ballistic-verification",
    "run.transport_frame": "inertial",
    "run.start_time_s": "0",
    "run.end_time_s": "86400",
    "run.campaign_seed_u64": "731993",
    "domain.solar_radius_m": "6.957e8",
    "domain.qualified_source_inner_radius_m": "7.02657e8",
    "domain.outer_radius_m": "2.0871e10",
    "solar_rotation.model": "rigid",
    "solar_rotation.rigid_input_rotation_rate_rad_per_s": "2.86532908457e-6",
    "solar_rotation.input_rate_convention": "sidereal",
    "solar_rotation.synodic_conversion_ephemeris_file": "none",
    "solar_rotation.differential_rotation_coefficients_file": "none",
    "pfss.outer_boundary_radius_m": "1.73925e9",
    "pfss.maximum_degree": "8",
    "pfss.coefficients_file": "analytic-dipole.coeff",
    "pfss.apodization_degree": "0",
    "current_sheet.model": "finite-shell-schatten",
    "current_sheet.interface_radius_m": "1.53054e9",
    "current_sheet.outer_radial_radius_m": "1.73925e9",
    "current_sheet.maximum_degree": "8",
    "current_sheet.maximum_negative_area_fraction": "0",
    "current_sheet.radialization_gate": "diagnostic-only",
    "current_sheet.maximum_outer_zonal_nonmonopole_power_fraction": "0",
    "current_sheet.latitude_diagnostic_radii_m": "6.957e9,1.3914e10",
    "current_sheet.latitude_minimum_unmasked_longitude_fraction": "0",
    "current_sheet.maximum_unsigned_radial_flux_rms_fraction": "0",
    "current_sheet.maximum_unsigned_radial_flux_p95_to_p05_ratio": "0",
    "closed_field_plasma.model": "isothermal-hydrostatic",
    "closed_field_plasma.reference_radius_m": "6.957e8",
    "closed_field_plasma.base_proton_number_density_m3": "1e14",
    "closed_field_plasma.base_proton_temperature_k": "1.5e6",
    "closed_field_plasma.base_electron_temperature_k": "1.5e6",
    "closed_field_plasma.closed_polytropic_index": "0",
    "closed_field_plasma.maximum_centrifugal_to_gravity_ratio": "0.1",
    "open_closed_interface.representation": "sharp-one-sided",
    "open_closed_interface.policy": "diagnostic-kinematic",
    "open_closed_interface.state_origin": "analytic-composite",
    "solar_wind.model": "flux-tube-polytropic",
    "solar_wind.scientific_role": "analytic-verification",
    "solar_wind.base_reference_radius_m": "6.957e8",
    "solar_wind.energy_closure": "polytropic",
    "solar_wind.polytropic_index": "1.05",
    "solar_wind.temperature_model": "uniform",
    "solar_wind.uniform_base_temperature_k": "1.5e6",
    "solar_wind.alpha_to_proton_ratio": "0",
    "solar_wind.momentum_residual_acceleration_floor_m_per_s2": "1e-12",
    "solar_wind.maximum_normalized_momentum_residual": "1",
    "solar_wind.maximum_absolute_momentum_residual_m_per_s2": "100",
    "plasma_eos.adiabatic_index": "1.6666666666666667",
    "plasma_eos.electron_mass_in_density": "neglect-recorded",
}


def numeric_key(key: str, expression: str) -> bool:
    if expression in {"REQUIRED_OR_ZERO", "0", "5"}:
        return True
    suffixes = (
        "_m", "_s", "_kg", "_c", "_k", "_tesla", "_pa", "_j", "_m3",
        "_rad", "_fraction", "_ratio", "_index", "_degree", "_points",
        "_level", "_steps", "_u64", "_slot", "_samples", "_cadence",
        "_order", "_power", "_tolerance", "_quantile", "_rate", "_count",
        "_number", "_weight", "_efficiency",
    )
    return key.endswith(suffixes)


def default_value(key: str, expression: str) -> str:
    alternatives = [part.strip() for part in expression.split("|")]
    literals = [part for part in alternatives
                if part not in {"REQUIRED", "REQUIRED_OR_ZERO"}]
    if literals:
        return literals[0]
    if expression == "REQUIRED_OR_ZERO":
        return "0"
    if numeric_key(key, expression):
        return "1"
    return "fixture"


def render() -> str:
    pattern = re.compile(r'^\{"([^"]+)", "([^"]+)", "([^"]+)"\},$')
    by_section: dict[str, list[tuple[str, str]]] = {section: [] for section in SECTIONS}
    for line in REGISTRY.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match and match.group(1) in by_section:
            by_section[match.group(1)].append((match.group(2), match.group(3)))
    lines = [
        "# Generated analytic-verification fixture for Stage-0 parser tests.",
        "# It is test evidence, not an event configuration or a source of defaults.",
        "",
    ]
    for section in SECTIONS:
        if not by_section[section]:
            raise RuntimeError(f"schema section {section} has no keys")
        lines.append(f"[{section}]")
        for key, expression in by_section[section]:
            dotted = f"{section}.{key}"
            value = OVERRIDES.get(dotted, default_value(key, expression))
            lines.append(f"{key} = {value}")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    generated = render()
    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != generated:
            print("STAGE0FIXTURE ERROR: generated fixture is stale", file=sys.stderr)
            return 1
        print("STAGE0FIXTURE PASS")
        return 0
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(OUTPUT, generated)
    print(f"generated {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
