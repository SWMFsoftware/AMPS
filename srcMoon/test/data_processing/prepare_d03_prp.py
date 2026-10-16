#!/usr/bin/env python3
"""Select native Diviner PRP triangles for the independent I20 polar table.

The processor preserves native temperature values and geometry.  It declares
an approximately equal-spacing target lattice from ``--target-spacing-m``
before any AMPS result is read, then retains the geographically nearest real
PRP triangle to each target.  It does not interpolate, average, smooth, or use
model residuals in selection.

The label defines ``temp_avg`` as annual-average temperature at 2 cm depth and
``temp_max`` as annual-maximum surface temperature.  They remain separate in
the output and must not be substituted for an instantaneous surface value.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


CANONICAL_ROOT = Path("/data/vtenishe/moon_validation_data")
LUNAR_REFERENCE_RADIUS_M = 1_737_400.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def angular_distance_rad(
    latitude_a_deg: float, longitude_a_deg: float,
    latitude_b_deg: float, longitude_b_deg: float,
) -> float:
    latitude_a = math.radians(latitude_a_deg)
    latitude_b = math.radians(latitude_b_deg)
    longitude_delta = math.radians(longitude_a_deg - longitude_b_deg)
    cosine = (math.sin(latitude_a) * math.sin(latitude_b) +
              math.cos(latitude_a) * math.cos(latitude_b) *
              math.cos(longitude_delta))
    return math.acos(max(-1.0, min(1.0, cosine)))


def target_bands(spacing_m: float, north: bool) -> list[tuple[float, int]]:
    """Return latitude and longitude-count pairs for a pole-safe lattice."""
    step_deg = math.degrees(spacing_m / LUNAR_REFERENCE_RADIUS_M)
    sign = 1.0 if north else -1.0
    latitudes: list[float] = []
    latitude = 75.0 + 0.5 * step_deg
    while latitude < 90.0:
        latitudes.append(sign * latitude)
        latitude += step_deg
    latitudes.append(sign * 90.0)
    bands: list[tuple[float, int]] = []
    for latitude in latitudes:
        circumference = (2.0 * math.pi * LUNAR_REFERENCE_RADIUS_M *
                         abs(math.cos(math.radians(latitude))))
        count = 1 if circumference < spacing_m else max(
            1, int(round(circumference / spacing_m)))
        bands.append((latitude, count))
    return bands


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=CANONICAL_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-spacing-m", type=float, default=100000.0)
    arguments = parser.parse_args()
    if (not math.isfinite(arguments.target_spacing_m) or
            arguments.target_spacing_m <= 0.0):
        parser.error("--target-spacing-m must be finite and positive")

    package = arguments.data_root.resolve() / "diviner"
    raw = package / "raw" / "prp"
    output = arguments.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    processed_rows: list[dict[str, object]] = []
    raw_hashes: dict[str, str] = {}
    counts: dict[str, object] = {}

    for hemisphere, north in (("north", True), ("south", False)):
        table = raw / f"dlre_prp_{hemisphere}.tab"
        label = raw / f"dlre_prp_{hemisphere}.xml"
        if not table.is_file() or not label.is_file():
            raise FileNotFoundError(f"missing D03 PRP pair for {hemisphere}")
        label_text = label.read_text(encoding="utf-8")
        required = (
            "<name>tri_clon</name>", "<name>tri_clat</name>",
            "<name>tri_calt</name>", "<name>temp_avg</name>",
            "<name>temp_max</name>", "a depth of 2 cm below the surface",
            "annual maximum temperature (K), calculated at",
            "<lro:latitude_type>Planetocentric</lro:latitude_type>",
            "<lro:longitude_direction>Positive East</lro:longitude_direction>",
        )
        if not all(token in label_text for token in required):
            raise ValueError(f"unsupported or incomplete D03 PRP label: {label}")
        raw_hashes[str(table.relative_to(package))] = sha256(table)
        raw_hashes[str(label.relative_to(package))] = sha256(label)

        bands = target_bands(arguments.target_spacing_m, north)
        # Keys are (latitude-band index, longitude index).  Each row examines
        # only one candidate longitude per latitude band, so all 5.76 million
        # native records can be streamed with bounded memory and no spatial
        # approximation in the final great-circle comparison.
        best: dict[tuple[int, int], tuple[float, dict[str, object]]] = {}
        input_rows = invalid_rows = 0
        with table.open(newline="", encoding="ascii") as stream:
            reader = csv.reader(stream, skipinitialspace=True)
            header = [field.strip() for field in next(reader)]
            expected = [
                "tri1_x", "tri1_y", "tri1_z", "tri2_x", "tri2_y",
                "tri2_z", "tri3_x", "tri3_y", "tri3_z", "tri_clon",
                "tri_clat", "tri_calt", "temp_avg", "temp_max", "ice_depth",
            ]
            if header != expected:
                raise ValueError(f"unexpected D03 PRP columns in {table}")
            for line_number, fields in enumerate(reader, start=2):
                input_rows += 1
                if len(fields) != 15:
                    raise ValueError(f"malformed D03 PRP row {table}:{line_number}")
                longitude = float(fields[9]) % 360.0
                latitude = float(fields[10])
                altitude_km = float(fields[11])
                average_k = float(fields[12])
                maximum_k = float(fields[13])
                if not (math.isfinite(latitude) and math.isfinite(longitude) and
                        math.isfinite(altitude_km) and average_k > 0.0 and
                        maximum_k > 0.0):
                    invalid_rows += 1
                    continue
                if abs(latitude) < 75.0:
                    continue
                for band_index, (target_latitude, longitude_count) in enumerate(bands):
                    longitude_index = int(round(
                        longitude / 360.0 * longitude_count)) % longitude_count
                    target_longitude = 360.0 * longitude_index / longitude_count
                    distance = angular_distance_rad(
                        latitude, longitude, target_latitude, target_longitude)
                    key = (band_index, longitude_index)
                    previous = best.get(key)
                    if previous is None or distance < previous[0]:
                        best[key] = (distance, {
                            "point_id": f"{hemisphere}-{band_index:02d}-{longitude_index:03d}",
                            "hemisphere": hemisphere,
                            "target_latitude_deg": f"{target_latitude:.12g}",
                            "target_longitude_east_deg": f"{target_longitude:.12g}",
                            "native_latitude_deg": f"{latitude:.12g}",
                            "native_longitude_east_deg": f"{longitude:.12g}",
                            "native_altitude_km": f"{altitude_km:.12g}",
                            "annual_average_temperature_at_2cm_k": f"{average_k:.12g}",
                            "annual_maximum_surface_temperature_k": f"{maximum_k:.12g}",
                            "uncertainty_k": "",
                            "source_record": f"{table.name}:line-{line_number}",
                            "target_distance_m": f"{distance * LUNAR_REFERENCE_RADIUS_M:.12g}",
                        })

        expected_targets = sum(count for _, count in bands)
        if len(best) != expected_targets:
            raise ValueError(
                f"D03 PRP selected {len(best)} of {expected_targets} "
                f"declared {hemisphere} targets")
        ordered = [best[(band_index, longitude_index)]
                   for band_index, (_, longitude_count) in enumerate(bands)
                   for longitude_index in range(longitude_count)]
        processed_rows.extend(item[1] for item in ordered)
        counts[hemisphere] = {
            "input_rows": input_rows,
            "invalid_temperature_or_geometry_rows": invalid_rows,
            "target_rows": expected_targets,
            "maximum_target_distance_m": max(item[0] for item in ordered) *
                                         LUNAR_REFERENCE_RADIUS_M,
        }

    table_path = output / "polar_temperature_reference.csv"
    with table_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(processed_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(processed_rows)
    qa = {
        "schema": "moon-validation-d03-prp-qa-v1",
        "status": "READY",
        "native_semantics": {
            "annual_average_temperature_at_2cm_k":
                "annual average at 2 cm depth",
            "annual_maximum_surface_temperature_k":
                "annual maximum at the surface",
            "uncertainty": "not supplied per triangle; field is empty",
            "latitude": "planetocentric",
            "longitude": "positive east",
        },
        "selection": {
            "minimum_absolute_latitude_deg": 75.0,
            "target_spacing_m": arguments.target_spacing_m,
            "operation": "nearest real native triangle; no interpolation or averaging",
            "model_residuals_used": False,
        },
        "counts": counts,
        "raw_sha256": raw_hashes,
        "processed": {"path": str(table_path), "rows": len(processed_rows),
                      "sha256": sha256(table_path)},
    }
    (output / "d03_prp_qa.json").write_text(
        json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qa, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
