#!/usr/bin/env python3
"""Reproduce the U05 Combi Figure 7 acceleration reference.

The historical machine-readable digitization used to populate ``Na.cpp`` is
not available.  This script therefore converts two independently recorded
pixel-coordinate traces of the published curve into physical units.  It does
not read AMPS output or the production table, so regenerating this reference
cannot tune the evidence to the model under test.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_MEASUREMENTS = HERE / "combi_1997_figure7_pixels.json"
EXPECTED_COLUMNS = [
    "velocity_km_s",
    "acceleration_150dpi_cm_s2",
    "acceleration_300dpi_cm_s2",
    "acceleration_mean_cm_s2",
    "cross_raster_difference_cm_s2",
    "uncertainty_cm_s2",
]


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest without changing the input file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def axis_value(pixel: float, pixel_min: float, pixel_max: float,
               value_min: float, value_max: float) -> float:
    """Map one linear plot coordinate to its documented physical axis."""
    return value_min + (pixel - pixel_min) * (value_max - value_min) / (
        pixel_max - pixel_min
    )


def decode_raster(raster: dict[str, Any]) -> dict[int, float]:
    """Convert the midpoint of each traced line thickness to cm/s^2."""
    axis = raster["axis"]
    decoded: dict[int, float] = {}
    for point in raster["points"]:
        target_velocity = int(point["target_velocity_km_s"])
        measured_velocity = axis_value(
            float(point["pixel_x"]),
            float(axis["x_pixel_at_min"]),
            float(axis["x_pixel_at_max"]),
            float(axis["x_velocity_min_km_s"]),
            float(axis["x_velocity_max_km_s"]),
        )
        # A click farther than 0.15 km/s from its named sample indicates a bad
        # axis calibration or transcription and must stop reference creation.
        if abs(measured_velocity - target_velocity) > 0.15:
            raise ValueError(
                f"{raster['dpi']} dpi point {target_velocity} km/s maps to "
                f"{measured_velocity:.6f} km/s"
            )

        line_midpoint = 0.5 * (
            float(point["curve_y_min_px"]) + float(point["curve_y_max_px"])
        )
        acceleration = axis_value(
            line_midpoint,
            float(axis["y_pixel_at_min"]),
            float(axis["y_pixel_at_max"]),
            float(axis["y_acceleration_min_cm_s2"]),
            float(axis["y_acceleration_max_cm_s2"]),
        )
        decoded[target_velocity] = acceleration
    return decoded


def write_json(path: Path, value: Any) -> None:
    """Write deterministic human-readable JSON used by provenance review."""
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measurements", type=Path, default=DEFAULT_MEASUREMENTS)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument(
        "--source-pdf",
        type=Path,
        help="Optional downloaded NTRS PDF; when supplied its SHA-256 is mandatory.",
    )
    args = parser.parse_args()

    measurements_path = args.measurements.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    measurements = json.loads(measurements_path.read_text(encoding="utf-8"))

    source_pdf_verified = False
    if args.source_pdf is not None:
        source_pdf = args.source_pdf.resolve()
        actual_pdf_hash = sha256(source_pdf)
        expected_pdf_hash = measurements["source"]["pdf_sha256"]
        if actual_pdf_hash != expected_pdf_hash:
            raise ValueError(
                "source PDF SHA-256 mismatch: "
                f"expected {expected_pdf_hash}, got {actual_pdf_hash}"
            )
        source_pdf_verified = True

    rasterizations = measurements["rasterizations"]
    if [item["dpi"] for item in rasterizations] != [150, 300]:
        raise ValueError("exactly the declared 150 and 300 dpi traces are required")
    decoded = [decode_raster(item) for item in rasterizations]
    expected_velocities = measurements["selection_rule"]["selected_velocity_km_s"]
    if list(decoded[0]) != expected_velocities or list(decoded[1]) != expected_velocities:
        raise ValueError("raster point order does not match the frozen selection rule")

    uncertainty = float(measurements["acceptance_uncertainty_cm_s2"])
    cross_limit = float(measurements["maximum_cross_raster_difference_cm_s2"])
    rows: list[dict[str, float | int]] = []
    maximum_cross_difference = 0.0
    for velocity in expected_velocities:
        first = decoded[0][velocity]
        second = decoded[1][velocity]
        cross_difference = abs(first - second)
        maximum_cross_difference = max(maximum_cross_difference, cross_difference)
        if cross_difference > cross_limit:
            raise ValueError(
                f"cross-raster mismatch at {velocity} km/s: "
                f"{cross_difference:.6f} cm/s^2 > {cross_limit:.6f} cm/s^2"
            )
        rows.append(
            {
                "velocity_km_s": velocity,
                "acceleration_150dpi_cm_s2": first,
                "acceleration_300dpi_cm_s2": second,
                "acceleration_mean_cm_s2": 0.5 * (first + second),
                "cross_raster_difference_cm_s2": cross_difference,
                "uncertainty_cm_s2": uncertainty,
            }
        )

    csv_path = output_dir / "combi_1997_figure7_acceleration.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=EXPECTED_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: row[key] if key == "velocity_km_s" else f"{row[key]:.9f}"
                    for key in EXPECTED_COLUMNS
                }
            )

    provenance = {
        "dataset": "Combi 1997 sodium radiation-pressure acceleration curve",
        "source": measurements["source"],
        "source_pdf_verified_during_generation": source_pdf_verified,
        "source_pdf_local_filename": (
            args.source_pdf.name if args.source_pdf is not None else "not supplied"
        ),
        "processing_script": Path(__file__).name,
        "processing_script_sha256": sha256(Path(__file__)),
        "measurement_file": measurements_path.name,
        "measurement_file_sha256": sha256(measurements_path),
        "transformations": [
            "PDF page rasterized independently at 150 and 300 dpi",
            "linear x-axis pixel-to-km/s conversion",
            "linear y-axis pixel-to-cm/s^2 conversion",
            "midpoint of visible curve line thickness",
            "arithmetic mean of the two raster-derived accelerations",
        ],
        "selection_rule": measurements["selection_rule"],
        "interpolation": "none",
        "smoothing": "none",
        "manual_digitization": True,
        "manual_edits_after_generation": False,
        "synthetic_data_used": False,
        "output_files": {csv_path.name: sha256(csv_path)},
        "notes": [
            "The unavailable historical digitization used by Na.cpp is not claimed as provenance.",
            "This package qualifies the published curve only to its frozen graphical uncertainty.",
        ],
    }
    provenance_path = output_dir / "provenance.json"
    write_json(provenance_path, provenance)

    qa = {
        "status": "PASS",
        "point_count": len(rows),
        "selected_velocity_km_s": expected_velocities,
        "excluded_velocity_km_s": measurements["selection_rule"][
            "excluded_velocity_km_s"
        ],
        "maximum_cross_raster_difference_cm_s2": maximum_cross_difference,
        "maximum_allowed_cross_raster_difference_cm_s2": cross_limit,
        "acceptance_uncertainty_cm_s2": uncertainty,
        "velocity_monotonic": expected_velocities == sorted(expected_velocities),
        "duplicate_velocity_count": len(expected_velocities) - len(set(expected_velocities)),
        "source_pdf_sha256_verified": source_pdf_verified,
        "rejected_or_filtered_point_count": len(
            measurements["selection_rule"]["excluded_velocity_km_s"]
        ),
        "rejection_basis": measurements["selection_rule"]["exclusion_reason"],
    }
    qa_path = output_dir / "qa_report.json"
    write_json(qa_path, qa)

    checksum_paths = [
        measurements_path,
        Path(__file__),
        csv_path,
        provenance_path,
        qa_path,
    ]
    checksum_path = output_dir / "SHA256SUMS.txt"
    checksum_path.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in checksum_paths),
        encoding="utf-8",
    )

    print(f"source_pdf_verified={source_pdf_verified}")
    print(f"records_produced={len(rows)}")
    print(
        "records_rejected_or_filtered="
        f"{len(measurements['selection_rule']['excluded_velocity_km_s'])}"
    )
    print(f"maximum_cross_raster_difference_cm_s2={maximum_cross_difference:.9f}")
    print(f"qa_status={qa['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
