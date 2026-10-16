#!/usr/bin/env python3
"""Prepare a resolution-matched Diviner GCP table and independent I16 rows.

Only the label-defined bolometric brightness temperature (``tbol``) is used.
The input is streamed because the 18 immutable latitude-band tables contain
about 25 million records.  Training-grid and verification-row membership is a
coordinate-index rule fixed before AMPS is run; model residuals cannot affect
which rows are retained.

The GCP archive does not provide a per-bin uncertainty column.  The processed
schema therefore leaves ``uncertainty_k`` empty and records that absence in QA;
zero is not substituted because zero would falsely assert an exact observation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re


CANONICAL_ROOT = Path("/data/vtenishe/moon_validation_data")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def coordinate_index(value: float, first: float, step: float) -> int:
    index = int(round((value - first) / step))
    reconstructed = first + index * step
    if not math.isclose(value, reconstructed, rel_tol=0.0, abs_tol=1.0e-7):
        raise ValueError(f"coordinate {value} is not on declared native grid")
    return index


def local_time_index(value: float) -> int:
    """Return the native quarter-hour bin index from its two-decimal centre.

    The PDS table serializes exact centres 0.125, 0.375, ... with two decimal
    places as the alternating sequence 0.12, 0.38, 0.62, 0.88.  Treating the
    printed values as an exact arithmetic 0.25-hour axis would reject every
    second record.  This decoder instead validates against that documented
    rounded representation and retains the native value in the output.
    """
    index = int(round(value * 4.0 - 0.5))
    if not 0 <= index < 96:
        raise ValueError(f"local time {value} lies outside the native day")
    centre = (index + 0.5) / 4.0
    if not math.isclose(value, centre, rel_tol=0.0, abs_tol=0.0050001):
        raise ValueError(
            f"local time {value} is not a two-decimal native bin centre")
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=CANONICAL_ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--target-spacing-m", type=float, default=100000.0,
        help="maximum coarse GCP spacing matched to the current surface receipt",
    )
    args = parser.parse_args()
    if not math.isfinite(args.target_spacing_m) or args.target_spacing_m <= 0:
        parser.error("--target-spacing-m must be finite and positive")

    package = args.data_root.resolve() / "diviner"
    raw = package / "raw" / "gcp"
    output = args.output_dir.resolve() if args.output_dir else package
    output.mkdir(parents=True, exist_ok=True)
    tables = sorted(raw.glob("global_cumul_avg_cyl_*_002.tab"))
    labels = sorted(raw.glob("global_cumul_avg_cyl_*_002.xml"))
    if len(tables) != 18 or len(labels) != 18:
        raise FileNotFoundError(
            "D03 GCP requires all 18 TAB/XML latitude bands")

    # Native centres are spaced by 0.5 degree.  Floor rather than round keeps
    # the processed spacing no larger than the requested physical resolution.
    lunar_radius_m = 1737400.0
    target_degrees = math.degrees(args.target_spacing_m / lunar_radius_m)
    stride = max(2, int(math.floor(target_degrees / 0.5)))
    # An even stride supplies a native midpoint used only for independent
    # verification.  If a future requested spacing yields an odd value, move
    # down one node rather than inventing an off-grid observation.
    if stride % 2:
        stride -= 1
    midpoint = stride // 2

    table_path = output / "temperature_table.csv"
    benchmark_path = output / "benchmark_points.csv"
    table_fields = [
        "latitude_deg", "longitude_east_deg", "local_time_hour",
        "temperature_k", "uncertainty_k", "quality_flag", "product_id",
        "source_record",
    ]
    benchmark_fields = table_fields + ["selection"]
    raw_hashes: dict[str, str] = {}
    input_rows = invalid_rows = table_rows = benchmark_rows = 0

    with table_path.open("w", newline="", encoding="utf-8") as table_stream, \
         benchmark_path.open("w", newline="", encoding="utf-8") as benchmark_stream:
        table_writer = csv.DictWriter(
            table_stream, fieldnames=table_fields, lineterminator="\n")
        benchmark_writer = csv.DictWriter(
            benchmark_stream, fieldnames=benchmark_fields, lineterminator="\n")
        table_writer.writeheader()
        benchmark_writer.writeheader()

        for path in tables:
            xml = path.with_suffix(".xml")
            text = xml.read_text(encoding="utf-8")
            required = (
                "<name>tbol</name>", "<unit>K</unit>",
                "<lro:latitude_type>Planetocentric</lro:latitude_type>",
                "<lro:longitude_direction>Positive East</lro:longitude_direction>",
            )
            if not all(token in text for token in required):
                raise ValueError(f"unsupported or incomplete D03 label: {xml}")
            product_match = re.search(r"<logical_identifier>([^<]+)</logical_identifier>", text)
            if product_match is None:
                raise ValueError(f"D03 label has no logical identifier: {xml}")
            product_id = product_match.group(1)
            raw_hashes[str(path.relative_to(package))] = sha256(path)
            raw_hashes[str(xml.relative_to(package))] = sha256(xml)

            with path.open(newline="", encoding="ascii") as stream:
                reader = csv.reader(stream, skipinitialspace=True)
                header = [field.strip() for field in next(reader)]
                if header != ["clon", "clat", "ltim", "t3", "t4", "t5",
                              "t6", "t7", "t8", "t9", "tbol"]:
                    raise ValueError(f"unexpected D03 GCP columns in {path}")
                for source_line, fields in enumerate(reader, start=2):
                    input_rows += 1
                    if len(fields) != 11:
                        raise ValueError(f"malformed D03 row {path}:{source_line}")
                    longitude_native = float(fields[0])
                    latitude = float(fields[1])
                    local_time = float(fields[2])
                    temperature = float(fields[10])
                    if temperature in (-9998.0, -9999.0) or temperature <= 0.0:
                        invalid_rows += 1
                        continue
                    longitude = longitude_native % 360.0
                    lat_index = coordinate_index(latitude, -89.75, 0.5)
                    lon_index = coordinate_index(longitude_native, -179.75, 0.5)
                    time_index = local_time_index(local_time)
                    common = {
                        "latitude_deg": f"{latitude:.2f}",
                        "longitude_east_deg": f"{longitude:.2f}",
                        "local_time_hour": f"{local_time:.2f}",
                        "temperature_k": f"{temperature:.7g}",
                        "uncertainty_k": "",
                        "quality_flag": "GOOD",
                        "product_id": product_id,
                        "source_record": f"{path.name}:line-{source_line}",
                    }
                    if lat_index % stride == 0 and lon_index % stride == 0:
                        table_writer.writerow(common)
                        table_rows += 1
                    elif (lat_index % stride == midpoint and
                          lon_index % stride == midpoint):
                        benchmark_writer.writerow({**common, "selection": "hold-out"})
                        benchmark_rows += 1

    qa = {
        "schema": "moon-validation-d03-gcp-qa-v1",
        "status": "READY" if table_rows and benchmark_rows else "ERROR",
        "native_product": {
            "type": "GCP version 1.0",
            "temperature_field": "tbol",
            "temperature_unit": "K",
            "latitude": "planetocentric",
            "longitude": "positive east",
            "native_grid": "0.5 degree latitude/longitude by 0.25 hour local time",
            "uncertainty": "not supplied per bin; processed uncertainty field is empty",
        },
        "selection": {
            "target_spacing_m": args.target_spacing_m,
            "coarse_stride_native_nodes": stride,
            "coarse_spacing_deg": 0.5 * stride,
            "holdout_offset_native_nodes": midpoint,
            "holdout_local_time_stride_native_nodes": 1,
            "linked_epoch_selection": (
                "I16 selects the one native quarter-hour bin nearest each "
                "point's physical local time at simulation_start_utc"
            ),
            "membership_frozen_before_model_scoring": True,
            "model_residuals_used": False,
        },
        "counts": {
            "input_rows": input_rows,
            "invalid_or_missing_tbol_rows": invalid_rows,
            "temperature_table_rows": table_rows,
            "holdout_rows": benchmark_rows,
        },
        "raw_sha256": raw_hashes,
        "processed": {
            "temperature_table": {"path": str(table_path),
                                  "sha256": sha256(table_path)},
            "benchmark_points": {"path": str(benchmark_path),
                                 "sha256": sha256(benchmark_path)},
        },
    }
    (output / "d03_gcp_qa.json").write_text(
        json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qa, indent=2, sort_keys=True))
    return 0 if qa["status"] == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
