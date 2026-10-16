#!/usr/bin/env python3
"""Create the independent D02 polar illumination comparison table.

The four input images are immutable PDS3 products.  This processor decodes
their signed little-endian samples, applies only the scale and offset declared
by each detached label, and maps selected pixel centres through the label's
spherical polar-stereographic projection.  Selection is a deterministic
spacing rule fixed before any AMPS result is read.

The edge flag is deliberately conservative: a sampled LPSR pixel is marked
ambiguous when any of its eight native neighbours has a different binary
classification.  I15/I20 may report those rows, but may not score them as
unambiguous point classifications.  No interpolation, smoothing, gap filling,
or model-dependent selection occurs here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import struct


CANONICAL_ROOT = Path("/data/vtenishe/moon_validation_data")
EXPECTED = {
    "AVGVISIB_65N_240M_201608.IMG":
        "4d96400a8609cfb87004010cc3affbefcb2064c4434750945a93bddc6c64254c",
    "AVGVISIB_65N_240M_201608.LBL":
        "ac839d42ea50ec34314d0c68bc3553d74a2719c88acca1ac7cf3b572a1742d71",
    "AVGVISIB_65S_240M_201608.IMG":
        "4b983fb1340c7b625a36d253e881874271c330083d2fdcc7277dc2f8ac24a438",
    "AVGVISIB_65S_240M_201608.LBL":
        "10b0d6d881713d6bd4d61a30cc64da6d53934e2c9a8cc702b36a5cd758751b03",
    "LPSR_65N_240M_201608.IMG":
        "b693690d9e7d1e31ef22cc90d1fedc06b78c7bf27d0f310e6dce69fa353fe128",
    "LPSR_65N_240M_201608.LBL":
        "0aaa186b0ff5003a32c7e6db18bd742a11c8ca22d18a8fdb82b47bd719b609a0",
    "LPSR_65S_240M_201608.IMG":
        "81eac93985118aefbe9ff5cc605018ae51708ebdcd9b80c895bdc14bfc0bb1ff",
    "LPSR_65S_240M_201608.LBL":
        "e546d988e903c64a77da73d28ea8058caf10953366f6d133d137911c1d9c5a3d",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scalar(text: str, key: str) -> float:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*([-+0-9.eE]+)", text)
    if match is None:
        raise ValueError(f"detached label is missing scalar {key}")
    return float(match.group(1))


def integer(text: str, key: str) -> int:
    return int(scalar(text, key))


def read_i16(path: Path, offset: int) -> int:
    with path.open("rb") as stream:
        stream.seek(offset)
        value = stream.read(2)
    if len(value) != 2:
        raise ValueError(f"short read from {path} at byte {offset}")
    return struct.unpack("<h", value)[0]


def inverse_polar_stereographic(
    line: int, sample: int, north: bool, scale_m: float,
    line_offset: float, sample_offset: float, radius_m: float,
) -> tuple[float, float]:
    """Return planetocentric latitude and east longitude for a pixel centre.

    PDS3 line/sample numbers are one-based.  The map-projection offsets locate
    the projection origin in that one-based system; increasing image line is
    downward, hence the sign in ``y``.  The formulas are the spherical polar
    stereographic inverse for a pole-centred, zero-central-longitude map.
    """

    x = ((sample + 1) - sample_offset - 1.0) * scale_m
    y_image_up = (line_offset + 1.0 - (line + 1)) * scale_m
    rho = math.hypot(x, y_image_up)
    angular_distance = 2.0 * math.atan2(rho, 2.0 * radius_m)
    if north:
        latitude = 90.0 - math.degrees(angular_distance)
        longitude = math.degrees(math.atan2(x, -y_image_up)) % 360.0
    else:
        latitude = -90.0 + math.degrees(angular_distance)
        longitude = math.degrees(math.atan2(x, y_image_up)) % 360.0
    return latitude, longitude


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=CANONICAL_ROOT)
    parser.add_argument(
        "--target-spacing-m", type=float, default=100000.0,
        help=("deterministic comparison-point spacing; defaults to the "
              "current 100 km production-surface receipt"),
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.target_spacing_m) or args.target_spacing_m <= 0:
        parser.error("--target-spacing-m must be finite and positive")

    package = args.data_root.resolve() / "lola"
    raw = package / "raw" / "illumination"
    output = args.output_dir.resolve() if args.output_dir else package
    output.mkdir(parents=True, exist_ok=True)

    hashes: dict[str, str] = {}
    for name, expected in EXPECTED.items():
        path = raw / name
        if not path.is_file():
            raise FileNotFoundError(f"missing immutable D02 file: {path}")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"D02 SHA-256 mismatch for {path}")
        hashes[f"raw/illumination/{name}"] = actual

    rows: list[dict[str, object]] = []
    qa_hemispheres: dict[str, object] = {}
    for hemisphere, north in (("N", True), ("S", False)):
        visibility_name = f"AVGVISIB_65{hemisphere}_240M_201608"
        psr_name = f"LPSR_65{hemisphere}_240M_201608"
        label_text = (raw / f"{visibility_name}.LBL").read_text(
            encoding="ascii")
        psr_label = (raw / f"{psr_name}.LBL").read_text()
        lines = integer(label_text, "LINES")
        samples = integer(label_text, "LINE_SAMPLES")
        if (lines, samples) != (
            integer(psr_label, "LINES"), integer(psr_label, "LINE_SAMPLES")
        ):
            raise ValueError("D02 AVGVISIB/LPSR dimensions differ")
        scale_m = scalar(label_text, "MAP_SCALE")
        radius_m = scalar(label_text, "A_AXIS_RADIUS") * 1000.0
        line_offset = scalar(label_text, "LINE_PROJECTION_OFFSET")
        sample_offset = scalar(label_text, "SAMPLE_PROJECTION_OFFSET")
        visibility_scale = scalar(label_text, "SCALING_FACTOR")
        visibility_offset = scalar(label_text, "OFFSET")
        psr_scale = scalar(psr_label, "SCALING_FACTOR")
        psr_offset = scalar(psr_label, "OFFSET")
        stride = max(1, int(round(args.target_spacing_m / scale_m)))
        visibility_path = raw / f"{visibility_name}.IMG"
        psr_path = raw / f"{psr_name}.IMG"
        retained = ambiguous = 0

        # Centre the deterministic lattice in the image rather than anchoring
        # it to one edge.  Changing AMPS output cannot change these indices.
        first = stride // 2
        for line in range(first, lines, stride):
            for sample in range(first, samples, stride):
                latitude, longitude = inverse_polar_stereographic(
                    line, sample, north, scale_m, line_offset,
                    sample_offset, radius_m)
                if abs(latitude) < 65.0:
                    continue
                offset = 2 * (line * samples + sample)
                visibility_dn = read_i16(visibility_path, offset)
                psr_dn = read_i16(psr_path, offset)
                visibility = visibility_dn * visibility_scale + visibility_offset
                psr_value = psr_dn * psr_scale + psr_offset
                if not (0.0 <= visibility <= 1.0) or not (
                    math.isclose(psr_value, 0.0, abs_tol=1.0e-12) or
                    math.isclose(psr_value, 1.0, abs_tol=1.0e-12)
                ):
                    raise ValueError("D02 scaled value violates label bounds")
                psr = int(round(psr_value))
                edge = False
                for delta_line in (-1, 0, 1):
                    for delta_sample in (-1, 0, 1):
                        if delta_line == delta_sample == 0:
                            continue
                        neighbour_line = line + delta_line
                        neighbour_sample = sample + delta_sample
                        if not (0 <= neighbour_line < lines and
                                0 <= neighbour_sample < samples):
                            edge = True
                            continue
                        neighbour_offset = 2 * (
                            neighbour_line * samples + neighbour_sample)
                        neighbour = read_i16(psr_path, neighbour_offset)
                        neighbour_value = neighbour * psr_scale + psr_offset
                        if int(round(neighbour_value)) != psr:
                            edge = True
                ambiguous += int(edge)
                retained += 1
                rows.append({
                    "point_id": f"{hemisphere}-{line + 1}-{sample + 1}",
                    "hemisphere": hemisphere,
                    "native_line": line + 1,
                    "native_sample": sample + 1,
                    "latitude_deg": f"{latitude:.12f}",
                    "longitude_east_deg": f"{longitude:.12f}",
                    "average_solar_disc_visibility": f"{visibility:.12g}",
                    "permanent_shadow": psr,
                    "edge_ambiguous": int(edge),
                    "native_resolution_m": f"{scale_m:.17g}",
                    "visibility_product_id": visibility_name,
                    "psr_product_id": psr_name,
                })
        qa_hemispheres[hemisphere] = {
            "stride_pixels": stride,
            "actual_spacing_m": stride * scale_m,
            "retained_points": retained,
            "edge_ambiguous_points": ambiguous,
        }

    table = output / "illumination_reference.csv"
    fields = list(rows[0]) if rows else []
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    qa = {
        "schema": "moon-validation-d02-processed-qa-v1",
        "status": "READY" if rows else "ERROR",
        "selection": {
            "target_spacing_m": args.target_spacing_m,
            "minimum_absolute_latitude_deg": 65.0,
            "edge_rule": "exclude from unambiguous score if any 8-neighbour LPSR value differs",
            "model_residuals_used": False,
        },
        "hemispheres": qa_hemispheres,
        "raw_sha256": hashes,
        "processed": {
            "path": str(table),
            "rows": len(rows),
            "sha256": sha256(table),
        },
        "illumination_operator": {
            "reference_cycle_years": 18.6,
            "reference_time_step_hours": 1.0,
            "horizon_azimuth_spacing_deg": 0.5,
            "solar_disc": "extended source; any visible fraction defines visibility",
            "source": "Mazarico et al. 2011; NASA PGDA product 69",
        },
    }
    (output / "d02_processed_qa.json").write_text(
        json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qa, indent=2, sort_keys=True))
    return 0 if rows else 2


if __name__ == "__main__":
    raise SystemExit(main())
