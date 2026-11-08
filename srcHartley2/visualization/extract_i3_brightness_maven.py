#!/usr/bin/env python3
"""Reproduce numerical MAVEN/IUVS brightness references from supplied PNGs.

The archive contains raster plots rather than source numerical tables.  The
reviewed pixel anchors live beside the generated table in
``references/i3_brightness_MAVEN/digitization_pixels.csv``.  This script maps
those anchors through the logarithmic axes and obtains radial-bin centers and
horizontal bounds from the independently transcribed case resolution table.
It does not use the red legacy-model curve in the source panels.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
from PIL import Image


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ARCHIVE = REPOSITORY_ROOT / "Archive.zip"
DEFAULT_CASE_TABLE = REPOSITORY_ROOT / "srcHartley2/i3/reference/i3_cases.csv"
DEFAULT_OUTPUT_DIR = (
    Path(__file__).resolve().parent / "references/i3_brightness_MAVEN"
)
DEFAULT_ANCHORS = DEFAULT_OUTPUT_DIR / "digitization_pixels.csv"
OUTPUT_NAME = "maven_iuvs_brightness.csv"
METADATA_NAME = "metadata.json"

X_MIN_1E5_KM = 0.2
X_MAX_1E5_KM = 20.0
Y_MIN_KR = 0.09
Y_MAX_KR = 20.0


class ExtractionError(RuntimeError):
    """Raised when source provenance or reviewed anchors are inconsistent."""


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Rebuild the reviewed 3I MAVEN/IUVS brightness reference table."
    )
    result.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    result.add_argument("--case-table", type=Path, default=DEFAULT_CASE_TABLE)
    result.add_argument("--anchor-table", type=Path, default=DEFAULT_ANCHORS)
    result.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return result


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance_path(path: Path) -> str:
    """Use repository-relative source names when possible for portability."""

    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(resolved)


def detect_plot_box(rgb: np.ndarray) -> Tuple[int, int, int, int]:
    """Detect left, right, top, and bottom black spines of an archive panel."""

    dark = np.mean(rgb, axis=2) < 0.5
    height, width = dark.shape
    rows = np.flatnonzero(np.sum(dark, axis=1) > 0.70 * width)
    columns = np.flatnonzero(np.sum(dark, axis=0) > 0.55 * height)
    rows = rows[(rows >= 5) & (rows < height - 5)]
    columns = columns[(columns >= 5) & (columns < width - 5)]
    if rows.size < 2 or columns.size < 2:
        raise ExtractionError("could not detect the long black plot spines")
    left, right = int(columns.min()), int(columns.max())
    top, bottom = int(rows.min()), int(rows.max())
    if right - left < width / 2 or bottom - top < height / 2:
        raise ExtractionError("detected plot box is implausibly small")
    return left, right, top, bottom


def load_resolutions(path: Path) -> Dict[int, float]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    result: Dict[int, float] = {}
    for row in rows:
        case_id = int(row["case_id"])
        resolution = float(row["resolution_1e5_km"])
        if case_id in result or not math.isfinite(resolution) or resolution <= 0.0:
            raise ExtractionError(f"invalid/duplicate case resolution for case {case_id}")
        result[case_id] = resolution
    if set(result) != set(range(1, 7)):
        raise ExtractionError("case table must contain exactly cases 1--6")
    return result


def optional_float(text: str) -> float | None:
    return None if text.strip() == "" else float(text)


def load_anchors(path: Path) -> List[Mapping[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    anchors: List[Mapping[str, object]] = []
    seen = set()
    for row in rows:
        key = (int(row["case_id"]), row["species"].strip().upper(), int(row["point_index"]))
        if key in seen:
            raise ExtractionError(f"duplicate reviewed anchor {key}")
        seen.add(key)
        top = optional_float(row["y_error_top_pixel"])
        bottom = optional_float(row["y_error_bottom_pixel"])
        if (top is None) != (bottom is None):
            raise ExtractionError(f"anchor {key} has only one vertical uncertainty endpoint")
        center = float(row["center_pixel_y"])
        if top is not None and not top < center < bottom:
            raise ExtractionError(f"anchor {key} requires top < center < bottom pixels")
        anchors.append(
            {
                "case_id": key[0],
                "species": key[1],
                "point_index": key[2],
                "center_pixel_y": center,
                "y_error_top_pixel": top,
                "y_error_bottom_pixel": bottom,
            }
        )
    if not anchors:
        raise ExtractionError("reviewed anchor table is empty")
    return anchors


def member_name(species: str, case_id: int) -> str:
    return (
        f"H_lyA_case{case_id}.png"
        if species == "H"
        else f"OH_brightness_case{case_id}.png"
    )


def x_to_pixel(value: float, left: int, right: int) -> float:
    fraction = (math.log10(value) - math.log10(X_MIN_1E5_KM)) / (
        math.log10(X_MAX_1E5_KM) - math.log10(X_MIN_1E5_KM)
    )
    return left + fraction * (right - left)


def pixel_to_y(pixel: float, top: int, bottom: int) -> float:
    fraction = (pixel - top) / (bottom - top)
    logarithm = math.log10(Y_MAX_KR) - fraction * (
        math.log10(Y_MAX_KR) - math.log10(Y_MIN_KR)
    )
    return 10.0**logarithm


def validate_black_anchor(rgb: np.ndarray, x: float, y: float, label: str) -> None:
    """Require black observation ink near every reviewed center coordinate."""

    xx, yy = int(round(x)), int(round(y))
    region = rgb[max(0, yy - 4) : yy + 5, max(0, xx - 4) : xx + 5, :]
    black = (np.mean(region, axis=2) < 0.55) & (
        np.max(region, axis=2) - np.min(region, axis=2) < 0.15
    )
    if not np.any(black):
        raise ExtractionError(f"no black observation ink near reviewed anchor {label}")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    resolutions = load_resolutions(arguments.case_table)
    anchors = load_anchors(arguments.anchor_table)

    panels: Dict[Tuple[int, str], Mapping[str, object]] = {}
    with zipfile.ZipFile(arguments.archive, "r") as archive:
        for case_id in range(1, 7):
            for species in ("H", "OH"):
                member = member_name(species, case_id)
                payload = archive.read(member)
                rgb = np.asarray(Image.open(io.BytesIO(payload)).convert("RGB"), dtype=float) / 255.0
                box = detect_plot_box(rgb)
                panels[(case_id, species)] = {
                    "member": member,
                    "payload_sha256": sha256_bytes(payload),
                    "width_pixels": int(rgb.shape[1]),
                    "height_pixels": int(rgb.shape[0]),
                    "plot_box_pixels": list(box),
                    "rgb": rgb,
                }

    output_rows: List[Dict[str, object]] = []
    for anchor in anchors:
        case_id = int(anchor["case_id"])
        species = str(anchor["species"])
        point_index = int(anchor["point_index"])
        panel = panels[(case_id, species)]
        left, right, top, bottom = panel["plot_box_pixels"]  # type: ignore[misc]
        resolution = resolutions[case_id]
        distance = point_index * resolution
        distance_lower = (point_index - 0.5) * resolution
        distance_upper = (point_index + 0.5) * resolution
        center_x = x_to_pixel(distance, left, right)
        center_y = float(anchor["center_pixel_y"])
        validate_black_anchor(
            panel["rgb"],  # type: ignore[arg-type]
            center_x,
            center_y,
            f"case={case_id},species={species},point={point_index}",
        )
        top_error = anchor["y_error_top_pixel"]
        bottom_error = anchor["y_error_bottom_pixel"]
        output_rows.append(
            {
                "case_id": case_id,
                "species": species,
                "point_index": point_index,
                "distance_1e5_km": f"{distance:.12g}",
                "distance_lower_1e5_km": f"{distance_lower:.12g}",
                "distance_upper_1e5_km": f"{distance_upper:.12g}",
                "brightness_kR": f"{pixel_to_y(center_y, top, bottom):.12g}",
                "brightness_lower_kR": (
                    ""
                    if bottom_error is None
                    else f"{pixel_to_y(float(bottom_error), top, bottom):.12g}"
                ),
                "brightness_upper_kR": (
                    ""
                    if top_error is None
                    else f"{pixel_to_y(float(top_error), top, bottom):.12g}"
                ),
                "source_member": panel["member"],
                "source_pixel_x": f"{center_x:.6f}",
                "source_pixel_y": f"{center_y:.6f}",
                "source_y_error_top_pixel": "" if top_error is None else f"{float(top_error):.6f}",
                "source_y_error_bottom_pixel": "" if bottom_error is None else f"{float(bottom_error):.6f}",
            }
        )

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = arguments.output_dir / OUTPUT_NAME
    fieldnames = list(output_rows[0])
    with output_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    metadata = {
        "schema_version": 1,
        "observable": "MAVEN/IUVS radial brightness",
        "distance_unit": "1e5 km",
        "brightness_unit": "kR",
        "source_archive": provenance_path(arguments.archive),
        "source_archive_sha256": sha256_file(arguments.archive),
        "case_table": provenance_path(arguments.case_table),
        "case_table_sha256": sha256_file(arguments.case_table),
        "anchor_table": provenance_path(arguments.anchor_table),
        "anchor_table_sha256": sha256_file(arguments.anchor_table),
        "output_table": OUTPUT_NAME,
        "output_table_sha256": sha256_file(output_path),
        "axis_limits": {
            "distance_1e5_km": [X_MIN_1E5_KM, X_MAX_1E5_KM],
            "brightness_kR": [Y_MIN_KR, Y_MAX_KR],
            "scales": ["log10", "log10"],
        },
        "distance_definition": (
            "point_index times case resolution_1e5_km; bounds are +/- one half resolution"
        ),
        "vertical_uncertainty_policy": (
            "Only visibly resolved vertical bars are tabulated; blank fields mean that no "
            "separate vertical uncertainty could be resolved from the raster."
        ),
        "panels": [
            {key: value for key, value in panel.items() if key != "rgb"}
            for _, panel in sorted(panels.items())
        ],
        "row_count": len(output_rows),
    }
    metadata_path = arguments.output_dir / METADATA_NAME
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(output_rows)} observations: {output_path}")
    print(f"wrote provenance: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
