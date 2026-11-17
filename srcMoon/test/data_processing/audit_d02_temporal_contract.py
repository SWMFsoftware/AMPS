#!/usr/bin/env python3
"""Audit frozen D02 native metadata without inventing a temporal reference.

This is data/operator preflight, never observational scoring. AVGVISIB records
the fraction of sampled times with ANY part of the solar disc visible; that
is a different observable from the mean visible disc area or incident flux.
START_TIME/STOP_TIME are retained as native fields, but their role as the
simulation's averaging interval is not established by the detached labels.

The inspected public documentation does not specify all temporal/geometry
controls for these exact 201608 products. Their unresolved fields remain null
and the result SKIPPED. No native or processed canonical data are written.
The immutable hash authority is shared with the existing D02 processor.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from prepare_d02_reference import EXPECTED, integer, scalar, sha256


UNRESOLVED_FIELDS = (
    "simulation_start_epoch",
    "simulation_stop_epoch",
    "simulation_time_system",
    "sample_count_and_endpoint_convention",
    "cadence_for_exact_201608_products",
    "finite_disc_geometry_and_ephemeris_corrections",
    "observer_height_for_exact_products",
    "earth_occultation_inclusion_policy",
)


def label_value(text: str, key: str) -> str:
    """Read one native scalar/string field, rejecting duplicates or absence.

    No date is parsed into an averaging epoch and no unit is converted here.
    Quoted strings and unquoted date tokens retain their native spelling.
    DESCRIPTION spans are interpreted only to identify the declared quantity.
    """
    matches = re.findall(r"(?m)^\s*" + re.escape(key) + r"\s*=\s*([^\r\n]+)", text)
    if len(matches) != 1:
        raise ValueError("expected exactly one detached-label field: " + key)
    return matches[0].strip().strip('"')


def inspect(package: Path, expected_hashes=None) -> dict:
    """Return a fail-closed, read-only contract audit.

    The CLI always uses committed EXPECTED hashes. The optional mapping exists
    solely for tiny corruption/parser protocol fixtures in the registered
    infrastructure self-test; those files are never scientific evidence.
    Missing required data is SKIPPED, altered bytes/invalid native metadata
    ERROR. Successful native integrity alone cannot resolve temporal fields.
    """
    expected = EXPECTED if expected_hashes is None else expected_hashes
    raw = package / "raw/illumination"
    result = {
        "schema": "srcMoon-D02-native-temporal-audit-v1",
        "scope": "native-metadata-and-integrity-only; not linked or observational validation",
        "status": "SKIPPED", "dataset": "D02", "package": str(package),
        "canonical_data_modified": False, "observational_scoring": False,
        "raw_files": [], "products": [], "missing_files": [],
        "temporal_contract": {name: None for name in UNRESOLVED_FIELDS},
        "unresolved_fields": list(UNRESOLVED_FIELDS),
        "reason": "exact-product temporal/geometry contract is unresolved; D02 remains PARTIAL",
    }
    try:
        for name, frozen in expected.items():
            path = raw / name
            if not path.is_file():
                result["missing_files"].append(str(path))
                continue
            actual = sha256(path)
            result["raw_files"].append({"path": str(path), "sha256": actual,
                "expected_sha256": frozen, "bytes": path.stat().st_size})
            if actual != frozen:
                raise ValueError("frozen D02 SHA-256 mismatch: " + str(path))
        if result["missing_files"]:
            result["reason"] = "required immutable D02 files unavailable"
            return result
        for name in sorted(expected):
            if not name.endswith(".LBL"):
                continue
            path = raw / name
            text = path.read_text(encoding="ascii")
            product = label_value(text, "PRODUCT_ID")
            if product != path.stem:
                raise ValueError("native product ID does not match file: " + str(path))
            lines, samples, bits = (integer(text, key) for key in
                                   ("LINES", "LINE_SAMPLES", "SAMPLE_BITS"))
            image = path.with_suffix(".IMG")
            if bits != 16 or lines <= 0 or samples <= 0 or (
                    image.stat().st_size != lines * samples * 2):
                raise ValueError("D02 image/label byte contract mismatch: " + str(image))
            if label_value(text, "SAMPLE_TYPE") != "LSB_INTEGER":
                raise ValueError("unsupported native D02 sample representation")
            frame = label_value(text, "COORDINATE_SYSTEM_NAME")
            if frame != "MEAN EARTH/POLAR AXIS OF DE421":
                raise ValueError("unexpected frozen D02 body-fixed frame")
            if label_value(text, "MAP_PROJECTION_TYPE") != "POLAR STEREOGRAPHIC":
                raise ValueError("unexpected frozen D02 projection")
            if label_value(text, "POSITIVE_LONGITUDE_DIRECTION") != "EAST":
                raise ValueError("unexpected frozen D02 longitude convention")
            if label_value(text, "KEYWORD_LATITUDE_TYPE") != "PLANETOCENTRIC":
                raise ValueError("unexpected frozen D02 latitude convention")
            if not label_value(text, "MAP_SCALE").endswith("<m/pix>"):
                raise ValueError("native D02 map scale is not explicitly m/pix")
            if label_value(text, "UNIT") != "NONE":
                raise ValueError("native D02 classification/visibility unit is not dimensionless")
            quantity = "permanent_shadow_classification"
            if product.startswith("AVGVISIB_"):
                normalized = " ".join(text.lower().split())
                if "any fraction of the solar disc" not in normalized:
                    raise ValueError("native AVGVISIB occurrence definition is absent")
                quantity = "fraction_of_sampled_times_any_solar_disc_visible"
            result["products"].append({
                "product_id": product, "quantity": quantity,
                "lines": lines, "samples": samples, "sample_bits": bits,
                "frame": frame, "map_scale_m_per_pixel": scalar(text, "MAP_SCALE"),
                "native_scaling_factor": scalar(text, "SCALING_FACTOR"),
                "native_offset": scalar(text, "OFFSET"),
                "native_unit": label_value(text, "UNIT"),
                "native_start_time": label_value(text, "START_TIME"),
                "native_stop_time": label_value(text, "STOP_TIME"),
                "label_times_used_as_simulation_bounds": False,
            })
        result["native_integrity_status"] = "PASS"
    except (OSError, ValueError, KeyError) as error:
        result.update(status="ERROR", reason=str(error))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=Path("/data/vtenishe/moon_validation_data"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    data_root = args.data_root.resolve()
    if out == data_root or data_root in out.parents:
        parser.error("audit artifacts must be outside the canonical validation-data tree")
    if out.exists() and any(out.iterdir()):
        parser.error("output directory must be new or empty")
    out.mkdir(parents=True, exist_ok=True)
    result = inspect(data_root / "lola")
    (out / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "reason": result["reason"],
                      "result": str(out / "result.json")}, indent=2))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
