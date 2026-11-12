#!/usr/bin/env python3
"""Export the fitted Enceladus source handoff consumed by the AMPS model.

The exporter intentionally does not choose a fitting strategy.  Its manifest
is the explicit, reviewable result selected by the fitting workflow; geometry
and gas kinematics come from the versioned standalone configuration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parent
# Import also works when the dependency-free exporter is loaded by file path.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from source_contract import CASES, HYBRID, SOLAR, source_ids, project_vertices
DEFAULT_CONFIG = ROOT / "config" / "tenishev2014_stage0.json"
SOURCE_IDS = (
    "jet:I_baghdad", "jet:II_damascus", "jet:III_damascus",
    "jet:IV_alexandria", "jet:V_cairo", "jet:VI_baghdad",
    "jet:VII_baghdad", "jet:VIII_cairo", "stripe:alexandria",
    "stripe:baghdad", "stripe:cairo", "stripe:damascus",
    "stripe:baghdad_extra", "stripe:cairo_extra", "stripe:damascus_extra",
)
JET_IDS = tuple(label[4:] for label in SOURCE_IDS[:8])
STRIPE_IDS = tuple(label[7:] for label in SOURCE_IDS[8:])
PORTABLE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
SPECIES_SYMBOL = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: root must be a JSON object")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_nonnegative(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label}: boolean is not a nonnegative scalar")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label}: expected a number") from error
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{label}: value must be finite and nonnegative")
    return result


def finite_number(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label}: boolean is not a physical scalar")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label}: expected a number") from error
    if not math.isfinite(result):
        raise ValueError(f"{label}: value must be finite")
    return result


def bounded_number(value: object, label: str, lower: float, upper: float,
                   *, lower_inclusive: bool = True) -> float:
    result = finite_number(value, label)
    lower_ok = result >= lower if lower_inclusive else result > lower
    if not lower_ok or result > upper:
        opening = "[" if lower_inclusive else "("
        raise ValueError(f"{label}: value must be in {opening}{lower},{upper}]")
    return result


def validate_manifest(manifest: dict) -> tuple[str, str, dict[str, dict[str, float]]]:
    allowed = {"schema_version", "source_case", "distribution_name", "mode",
               "species", "fit_provenance"}
    unknown = sorted(set(manifest) - allowed)
    if unknown:
        raise ValueError("unknown manifest keys: " + ", ".join(unknown))
    if manifest.get("schema_version") != 1:
        raise ValueError("manifest schema_version must be integer 1")
    if manifest.get("source_case") not in CASES:
        raise ValueError("unsupported source_case")
    required_ids = source_ids(manifest["source_case"])
    name = manifest.get("distribution_name")
    if not isinstance(name, str) or PORTABLE_NAME.fullmatch(name) is None:
        raise ValueError(
            "distribution_name must match [A-Za-z0-9][A-Za-z0-9_.-]*")
    mode = manifest.get("mode")
    if mode not in ("shape", "absolute"):
        raise ValueError("mode must be shape or absolute")
    species = manifest.get("species")
    if not isinstance(species, dict) or not species:
        raise ValueError("species must be a nonempty object")
    validated: dict[str, dict[str, float]] = {}
    for symbol, strengths in species.items():
        if not isinstance(symbol, str) or SPECIES_SYMBOL.fullmatch(symbol) is None:
            raise ValueError("species names must match [A-Za-z][A-Za-z0-9_]*")
        if not isinstance(strengths, dict):
            raise ValueError(f"species.{symbol} must be an object")
        if not all(isinstance(label, str) for label in strengths):
            raise ValueError(f"species.{symbol}: every source label must be a string")
        missing = sorted(set(required_ids) - set(strengths))
        extra = sorted(set(strengths) - set(required_ids))
        if missing or extra:
            raise ValueError(
                f"species.{symbol} labels differ from the canonical {len(required_ids)}; "
                f"missing={missing}, extra={extra}")
        values = {source: finite_nonnegative(strengths[source],
                                             f"species.{symbol}.{source}")
                  for source in required_ids}
        total = math.fsum(values.values())
        if not total > 0:
            raise ValueError(f"species.{symbol}: total strength must be positive")
        if mode == "shape" and not math.isclose(total, 1.0, rel_tol=0,
                                                 abs_tol=1e-12):
            raise ValueError(
                f"species.{symbol}: shape sums to {total:.17g}, not 1")
        validated[symbol] = values
    provenance = manifest.get("fit_provenance", {})
    if not isinstance(provenance, dict):
        raise ValueError("fit_provenance must be an object when present")
    return name, mode, validated


def scalar(value: float) -> str:
    return f"{value:.17g}"


def render(manifest_path: Path, config_path: Path) -> str:
    manifest = load_json(manifest_path)
    config = load_json(config_path)
    name, mode, species = validate_manifest(manifest)
    source_case = manifest["source_case"]
    if config.get("source_case") != source_case:
        raise ValueError(f"{config_path}: wrong source_case")
    if config.get("frame") != "IAU_ENCELADUS":
        raise ValueError(f"{config_path}: frame must be IAU_ENCELADUS")
    if config.get("source_law") != "tenishev2010_point":
        raise ValueError(f"{config_path}: source_law must be tenishev2010_point")
    if config.get("schema_version") != 1:
        raise ValueError(f"{config_path}: schema_version must be integer 1")
    reference_radius = bounded_number(
        config.get("reference_radius_m"),
        f"{config_path}: reference_radius_m", 0.0, math.inf,
        lower_inclusive=False)
    # Automatic fit exports carry a fully resolved geometry snapshot. Catalog
    # data are thus not needed by AMPS, and optimized kinematics are not reset
    # to the checked-in reference configuration. Legacy explicit manifests
    # retain their original raw-geometry projection path and provenance text.
    resolved = bool(config.get("_resolved_export_geometry", False))
    stripe_path = None
    stripes = {}
    if config.get("tiger_stripes"):
        stripe_path = Path(config.get("_stripe_geometry_path") or config["stripe_geometry_file"])
        if not stripe_path.is_absolute():
            stripe_path = (config_path.parent / stripe_path).resolve()
        stripe_document = load_json(stripe_path)
        stripes = stripe_document.get("stripes")
        if not isinstance(stripes, dict):
            raise ValueError(f"{stripe_path}: stripes object not found")
    provenance = json.dumps(manifest.get("fit_provenance", {}),
                            sort_keys=True, separators=(",", ":"))
    lines = [
        "! Generated by standalone/export_forward_sources.py; do not hand-normalize.",
        f"! manifest={manifest_path.resolve()}",
        f"! manifest_sha256={sha256(manifest_path)}",
        f"! geometry_config={config_path.resolve()}",
        f"! geometry_config_sha256={sha256(config_path)}",
        f"! stripe_geometry={stripe_path}",
        f"! stripe_geometry_sha256={sha256(stripe_path) if stripe_path else 'not_applicable'}",
        f"! fit_provenance={provenance}",
        f"#begin source_distribution: {name}",
        "schema_version=2",
        f"source_case={source_case}",
        "frame=IAU_ENCELADUS",
        f"mode={mode}",
        "source_law=tenishev2010_point",
        "stripe_surface_mapping=local_vertical_to_reference_sphere",
        f"stripe_reference_radius_m={scalar(reference_radius)}",
        "",
    ]

    jets = config.get("discrete_jets")
    stripe_config = config.get("tiger_stripes")
    required = source_ids(source_case)
    expected_jets = {label[4:] for label in required if label.startswith("jet:")}
    expected_stripes = {label[7:] for label in required if label.startswith("stripe:")}
    if not isinstance(jets, list) or not isinstance(stripe_config, list):
        raise ValueError("source geometry must contain jet and stripe lists")
    if not all(isinstance(item, dict) for item in jets + stripe_config):
        raise ValueError("every source entry must be an object")
    jet_ids = [item.get("id") for item in jets]
    stripe_ids = [item.get("id") for item in stripe_config]
    if len(jet_ids) != len(expected_jets) or set(jet_ids) != expected_jets:
        raise ValueError("discrete_jets must contain each canonical jet exactly once")
    if len(stripe_ids) != len(expected_stripes) or set(stripe_ids) != expected_stripes:
        raise ValueError("tiger_stripes must contain each canonical stripe exactly once")

    def add_strengths(label: str) -> None:
        for symbol in sorted(species):
            lines.append(f"strength.{symbol}={scalar(species[symbol][label])}")

    for jet in jets:
        source_id = jet["id"]
        label = "jet:" + source_id
        latitude = bounded_number(jet["latitude_deg"],
                                  f"jet {source_id} latitude_deg", -90.0, 90.0)
        longitude = bounded_number(jet["west_longitude_deg"],
                                   f"jet {source_id} west_longitude_deg",
                                   0.0, 360.0)
        if longitude == 360.0:
            raise ValueError(
                f"jet {source_id} west_longitude_deg: value must be in [0,360)")
        tilt = bounded_number(jet["tilt_deg"], f"jet {source_id} tilt_deg",
                              0.0, 90.0)
        azimuth = bounded_number(jet["azimuth_deg"],
                                 f"jet {source_id} azimuth_deg", 0.0, 360.0)
        if azimuth == 360.0:
            raise ValueError(
                f"jet {source_id} azimuth_deg: value must be in [0,360)")
        if source_case == SOLAR and (tilt != 0.0 or azimuth != 0.0):
            raise ValueError("Portyankina solar sources require local-normal axes")
        speed = finite_nonnegative(jet["bulk_speed_m_s"],
                                   f"jet {source_id} bulk_speed_m_s")
        temperature = bounded_number(jet["temperature_K"],
                                     f"jet {source_id} temperature_K", 0.0,
                                     math.inf, lower_inclusive=False)
        lines.extend([
            f"#begin source: {source_id}",
            "geometry=point_jet",
            f"latitude_deg={scalar(latitude)}",
            f"west_longitude_deg={scalar(longitude)}",
            f"tilt_deg={scalar(tilt)}",
            f"azimuth_deg={scalar(azimuth)}",
            f"bulk_speed_m_s={scalar(speed)}",
            f"temperature_K={scalar(temperature)}",
        ])
        add_strengths(label)
        lines.extend(["#end source", ""])

    for item in stripe_config:
        source_id = item["id"]
        label = "stripe:" + source_id
        geometry_id = item["geometry_id"]
        if geometry_id != source_id or geometry_id not in stripes:
            raise ValueError(f"unexpected stripe id/geometry {source_id}/{geometry_id}")
        points = item.get("points_m") if resolved else stripes[geometry_id]
        if not isinstance(points, list) or len(points) < 2:
            raise ValueError(f"stripe {geometry_id} needs at least two vertices")
        speed = finite_nonnegative(item["bulk_speed_m_s"],
                                   f"stripe {source_id} bulk_speed_m_s")
        temperature = bounded_number(item["temperature_K"],
                                     f"stripe {source_id} temperature_K", 0.0,
                                     math.inf, lower_inclusive=False)
        lines.extend([
            f"#begin source: {source_id}",
            "geometry=polyline_stripe",
            f"bulk_speed_m_s={scalar(speed)}",
            f"temperature_K={scalar(temperature)}",
        ])
        add_strengths(label)
        previous = None
        for index, point in enumerate(points):
            if not isinstance(point, list) or len(point) != 3:
                raise ValueError(f"stripe {geometry_id} vertex {index} is invalid")
            coordinates = tuple(
                finite_number(value,
                              f"stripe {geometry_id} vertex {index} coordinate")
                for value in point)
            radius = math.sqrt(math.fsum(value * value for value in coordinates))
            if radius == 0:
                raise ValueError(f"stripe {geometry_id} vertex {index} is at r=0")
            # The published Cartesian nodes encode stripe longitude/latitude
            # on a smaller legacy shape. "Vertical" for the currently
            # configured spherical Enceladus body is the local outward radial
            # direction. Preserve that direction exactly and replace only the
            # radius, placing every emitted stripe node on the declared source
            # surface. This is a geometry mapping, not an altitude offset.
            if resolved:
                if abs(radius - reference_radius) > 1.0e-6:
                    raise ValueError("resolved fitted stripe is not on the AMPS reference sphere")
                projected = coordinates
            else:
                scale = reference_radius / radius
                projected = tuple(scale * value for value in coordinates)
            if previous == projected:
                raise ValueError(
                    f"stripe {geometry_id} segment ending at {index} has zero length")
            previous = projected
            lines.extend([
                f"#begin vertex: {index:04d}",
                f"x_m={scalar(projected[0])}",
                f"y_m={scalar(projected[1])}",
                f"z_m={scalar(projected[2])}",
                "#end vertex",
            ])
        lines.extend(["#end source", ""])
    lines.extend(["#end source_distribution", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        text = render(args.manifest.resolve(), args.config.resolve())
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation prevents a rerun from silently replacing the
        # source receipt used by a campaign.
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
    except (OSError, KeyError, TypeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
