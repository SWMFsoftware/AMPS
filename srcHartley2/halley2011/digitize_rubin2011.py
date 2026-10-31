#!/usr/bin/env python3
"""Extract the published Rubin et al. (2011) neutral-coma curves.

The journal PDF stores the curves in Figs. 1--18 as colored vector paths.
This utility converts those paths to physical coordinates using the axes
printed in the paper.  It deliberately extracts only quantities that the
generic comet application calculates directly: neutral number densities in
Figs. 1--17 and the modeled H2O radial velocity in Fig. 18.  Red NMS response
curves and measured symbols are different observables and are not mixed into
the neutral-profile comparison.

The PDF uses a split logarithmic distance axis.  From the nucleus to 1000 km
the plotted interval is compressed; outside 1000 km a different scale is
used.  AXES records the vector coordinates of the tick marks, so applying one
logarithmic transform to the full plot would be an incorrect digitization.

The input hash is checked because the vector coordinates below belong to the
specific author-supplied PDF named in the repository documentation.  The CSV
is committed with the application, so routine comparison plotting does not
require access to the original PDF or to pdftocairo.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


EXPECTED_PDF_SHA256 = "c245a8a5fada2090a8fd011685818fb1711b99181d318219e1bb9e858217ce78"
SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"

# Cairo reports the green path with a small red component inherited from the
# PDF color profile.  Matching the exact SVG spelling avoids selecting black
# axes, red NMS curves, or page text.
GREEN = "0.392151%,100%,0%"
BLUE = "0%,0%,100%"
CYAN = "0%,100%,100%"
RED = "100%,0%,0%"


@dataclass(frozen=True)
class Axis:
    """Vector coordinates and physical limits for one paper panel."""

    page: int
    x_left: float
    x_break: float
    x_right: float
    y_top: float
    y_bottom: float
    x_min_km: float = 6.0
    x_break_km: float = 1000.0
    x_max_km: float = 100000.0
    y_min: float = 1.0e1
    y_max: float = 1.0e8
    quantity: str = "number_density"
    unit: str = "cm^-3"


@dataclass(frozen=True)
class CurveSpec:
    """Identify a published vector path by panel, color, and dash pattern."""

    figure: int
    species: str
    color: str
    variant: str = "nominal"
    dash_contains: Optional[str] = None


# Coordinates are taken from the major tick paths in the SVG produced by
# pdftocairo 22.02.0.  The values are PDF points after the path transform, not
# raster pixels.  For density panels, y_top and y_bottom are 10^8 and 10^1
# cm^-3.  Figure 18 uses logarithmic velocity ticks from 400 to 3000 m/s.
AXES: Dict[int, Axis] = {
    1: Axis(8, 166.184, 233.973, 436.395, 510.422, 704.480),
    2: Axis(9, 175.418, 243.207, 445.629, 74.199, 268.258),
    3: Axis(9, 174.285, 242.078, 444.496, 367.922, 561.980),
    4: Axis(10, 164.484, 232.273, 434.695, 74.199, 268.258),
    5: Axis(10, 165.617, 233.406, 435.828, 354.840, 548.895),
    6: Axis(11, 175.418, 243.207, 445.633, 74.199, 268.258),
    7: Axis(11, 174.285, 242.078, 444.500, 367.871, 561.926),
    8: Axis(12, 164.484, 232.273, 434.695, 74.199, 268.258),
    9: Axis(12, 166.184, 233.973, 436.395, 356.891, 550.949),
    10: Axis(13, 174.852, 242.645, 445.066, 74.199, 268.258),
    11: Axis(13, 174.285, 242.078, 444.500, 356.829, 550.885),
    12: Axis(14, 164.484, 232.273, 434.695, 74.199, 268.258),
    13: Axis(14, 164.484, 232.273, 434.695, 370.198, 564.256),
    14: Axis(15, 174.285, 242.078, 444.500, 74.199, 268.258),
    15: Axis(15, 176.137, 243.926, 446.348, 371.691, 565.750),
    16: Axis(16, 164.829, 232.398, 434.172, 74.184, 267.613),
    17: Axis(16, 165.047, 232.836, 435.258, 355.762, 549.816),
    18: Axis(
        17,
        177.148,
        246.484,
        453.520,
        479.891,
        678.371,
        x_max_km=1.0e6,
        y_min=400.0,
        y_max=3000.0,
        quantity="radial_velocity",
        unit="m/s",
    ),
}


CURVES: Tuple[CurveSpec, ...] = (
    CurveSpec(1, "H", GREEN),
    CurveSpec(1, "H2", BLUE),
    CurveSpec(2, "C", GREEN),
    CurveSpec(3, "CH", GREEN),
    CurveSpec(4, "CH2", GREEN),
    CurveSpec(4, "N", BLUE),
    CurveSpec(5, "CH3", GREEN),
    CurveSpec(5, "NH", BLUE),
    CurveSpec(6, "CH4", GREEN),
    CurveSpec(6, "NH2", BLUE),
    CurveSpec(6, "O", CYAN),
    CurveSpec(7, "NH3", BLUE),
    CurveSpec(7, "OH", CYAN),
    CurveSpec(8, "H2O", CYAN),
    CurveSpec(9, "C2", BLUE),
    CurveSpec(10, "C2H", BLUE),
    CurveSpec(11, "C2H2", BLUE),
    CurveSpec(11, "CN", GREEN),
    CurveSpec(12, "HCN", GREEN),
    CurveSpec(13, "CO", GREEN),
    CurveSpec(13, "C2H4", BLUE),
    CurveSpec(14, "HCO", GREEN),
    CurveSpec(14, "C2H5", BLUE),
    CurveSpec(15, "H2CO", GREEN, dash_contains="2.7167"),
    # The asterisked curve in Fig. 15 is the paper's explicit case without
    # CH3OH photodissociation to H2CO.  halley1d.input currently selects this
    # chemistry branch, so retain both variants and identify them separately.
    CurveSpec(
        15,
        "H2CO",
        GREEN,
        variant="no_ch3oh_photodissociation",
        dash_contains="0.45279,0.67918",
    ),
    CurveSpec(15, "C2H6", BLUE),
    CurveSpec(16, "CH3OH", GREEN),
    CurveSpec(17, "CO2", GREEN),
    CurveSpec(18, "H2O", RED),
)


NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="extract neutral profiles from Rubin et al. (2011) vector figures"
    )
    parser.add_argument("paper_pdf", type=Path, help="Rubin-2011-icarus.pdf")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "reference" / "rubin2011_neutral_profiles.csv",
        help="output CSV path",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="metadata JSON path (default: OUTPUT with .metadata.json suffix)",
    )
    return parser.parse_args(argv)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def transformed_path_points(node: ET.Element) -> List[Tuple[float, float]]:
    """Return M/L path vertices after applying the SVG affine transform."""

    path_data = node.attrib.get("d", "")
    # Published profile paths contain only M and L commands.  Refuse curves
    # with another drawing command instead of silently pairing control points.
    commands = set(re.findall(r"[A-Za-z]", path_data))
    if not commands.issubset({"M", "L"}):
        return []
    values = [float(value) for value in NUMBER.findall(path_data)]
    if len(values) < 4 or len(values) % 2:
        return []
    points = list(zip(values[0::2], values[1::2]))

    transform = [float(value) for value in NUMBER.findall(node.attrib.get("transform", ""))]
    if len(transform) != 6:
        return []
    a, b, c, d, e, f = transform
    return [(a * x + c * y + e, b * x + d * y + f) for x, y in points]


def extract_svg_points(svg: Path, spec: CurveSpec) -> List[Tuple[float, float]]:
    """Collect the long vector paths belonging to one published curve."""

    axis = AXES[spec.figure]
    result: List[Tuple[float, float]] = []
    root = ET.parse(svg).getroot()
    for node in root.iter(f"{SVG_NAMESPACE}path"):
        style = node.attrib.get("style", "")
        if f"stroke:rgb({spec.color})" not in style:
            continue
        if spec.dash_contains is not None and spec.dash_contains not in style:
            continue
        points = transformed_path_points(node)
        if len(points) < 2:
            continue
        xs = [point[0] for point in points]
        ys = [point[1] for point in points]
        # Legend samples are horizontal 21-point paths.  Measured Fig. 18
        # symbols are roughly 3 points wide.  Published model segments span at
        # least 20 points and vary in y, including their clipped end pieces.
        # The steep first segment of Fig. 18 spans only about eight PDF points;
        # its measured x/square symbols remain below three points.  Density
        # paths use the wider threshold because their clipped curve pieces are
        # longer and the extra separation rejects all plot symbols.
        minimum_x_span = 5.0 if spec.figure == 18 else 10.0
        if max(xs) - min(xs) < minimum_x_span or max(ys) - min(ys) < 0.1:
            continue
        if max(xs) < axis.x_left or min(xs) > axis.x_right:
            continue
        if max(ys) < axis.y_top or min(ys) > axis.y_bottom:
            continue
        result.extend(
            (x, y)
            for x, y in points
            if axis.x_left - 0.02 <= x <= axis.x_right + 0.02
            and axis.y_top - 0.02 <= y <= axis.y_bottom + 0.02
        )
    return result


def physical_radius(axis: Axis, x: float) -> float:
    """Apply the paper's two distinct logarithmic distance scales."""

    if x <= axis.x_break:
        fraction = (x - axis.x_left) / (axis.x_break - axis.x_left)
        log_radius = math.log10(axis.x_min_km) + fraction * (
            math.log10(axis.x_break_km) - math.log10(axis.x_min_km)
        )
    else:
        fraction = (x - axis.x_break) / (axis.x_right - axis.x_break)
        log_radius = math.log10(axis.x_break_km) + fraction * (
            math.log10(axis.x_max_km) - math.log10(axis.x_break_km)
        )
    return 10.0**log_radius


def physical_value(axis: Axis, y: float) -> float:
    """Convert the logarithmic density or velocity ordinate."""

    fraction = (y - axis.y_top) / (axis.y_bottom - axis.y_top)
    log_value = math.log10(axis.y_max) + fraction * (
        math.log10(axis.y_min) - math.log10(axis.y_max)
    )
    return 10.0**log_value


def collapse_vertices(
    axis: Axis, points: Iterable[Tuple[float, float]]
) -> List[Tuple[float, float]]:
    """Collapse repeated SVG x vertices to one log-mean physical sample."""

    by_x: Dict[float, List[float]] = {}
    original_x: Dict[float, float] = {}
    for x, y in points:
        key = round(x, 6)
        by_x.setdefault(key, []).append(math.log(physical_value(axis, y)))
        original_x[key] = x
    collapsed = [
        (
            physical_radius(axis, original_x[key]),
            math.exp(sum(log_values) / len(log_values)),
        )
        for key, log_values in by_x.items()
    ]
    collapsed.sort()
    return collapsed


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not args.paper_pdf.is_file():
        print(f"digitize_rubin2011.py: PDF does not exist: {args.paper_pdf}", file=sys.stderr)
        return 2
    actual_hash = sha256(args.paper_pdf)
    if actual_hash != EXPECTED_PDF_SHA256:
        print(
            "digitize_rubin2011.py: PDF hash differs from the documented Rubin paper; "
            "the stored vector coordinates cannot be applied safely\n"
            f"expected {EXPECTED_PDF_SHA256}\nactual   {actual_hash}",
            file=sys.stderr,
        )
        return 2

    output = args.output
    metadata = args.metadata or output.with_suffix(".metadata.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata.parent.mkdir(parents=True, exist_ok=True)

    extracted: Dict[Tuple[int, str, str], List[Tuple[float, float]]] = {}
    with tempfile.TemporaryDirectory(prefix="rubin2011-svg-") as temporary:
        temporary_path = Path(temporary)
        svg_by_page: Dict[int, Path] = {}
        for page in sorted({axis.page for axis in AXES.values()}):
            svg = temporary_path / f"page-{page}.svg"
            subprocess.run(
                [
                    "pdftocairo",
                    "-f",
                    str(page),
                    "-l",
                    str(page),
                    "-svg",
                    str(args.paper_pdf),
                    str(svg),
                ],
                check=True,
            )
            svg_by_page[page] = svg

        for spec in CURVES:
            axis = AXES[spec.figure]
            points = collapse_vertices(axis, extract_svg_points(svg_by_page[axis.page], spec))
            if len(points) < 20:
                print(
                    f"digitize_rubin2011.py: only {len(points)} points extracted for "
                    f"Fig. {spec.figure} {spec.species} ({spec.variant})",
                    file=sys.stderr,
                )
                return 2
            extracted[(spec.figure, spec.species, spec.variant)] = points

    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "figure",
                "pdf_page",
                "species",
                "variant",
                "quantity",
                "radius_km",
                "value",
                "unit",
            ]
        )
        for spec in CURVES:
            axis = AXES[spec.figure]
            for radius, value in extracted[(spec.figure, spec.species, spec.variant)]:
                writer.writerow(
                    [
                        spec.figure,
                        axis.page,
                        spec.species,
                        spec.variant,
                        axis.quantity,
                        f"{radius:.10e}",
                        f"{value:.10e}",
                        axis.unit,
                    ]
                )

    records = []
    for spec in CURVES:
        points = extracted[(spec.figure, spec.species, spec.variant)]
        records.append(
            {
                "figure": spec.figure,
                "pdf_page": AXES[spec.figure].page,
                "species": spec.species,
                "variant": spec.variant,
                "quantity": AXES[spec.figure].quantity,
                "samples": len(points),
                "minimum_radius_km": points[0][0],
                "maximum_radius_km": points[-1][0],
            }
        )
    metadata.write_text(
        json.dumps(
            {
                "paper": "Rubin et al., Icarus 213 (2011) 655-677",
                "doi": "10.1016/j.icarus.2011.04.006",
                "source_pdf_sha256": actual_hash,
                "method": (
                    "pdftocairo SVG vector-path extraction with the published split-log "
                    "distance axes; red NMS observables excluded"
                ),
                "axis_coordinates_pdf_points": {
                    str(figure): asdict(axis) for figure, axis in AXES.items()
                },
                "curves": records,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {sum(len(points) for points in extracted.values())} reference points to {output}")
    print(f"wrote provenance to {metadata}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
