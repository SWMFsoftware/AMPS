#!/usr/bin/env python3
"""Extract the XM03 Earth validation data from Liu et al. Figure 12.

The arXiv source distribution contains Figure 12 as a vector PDF created by
Matplotlib.  This script converts that PDF to SVG and reads the actual vector
coordinates of the plotted error bars and Earth shock-source curve.  It does
not estimate values from raster pixels and it never reads a srcSEP result.

Two products are generated:

* ``earth_observations.csv`` contains the ACE/EPAM, GOES-13/EPEAD, and
  SOHO/ERNE Earth measurements plotted in panels (a)-(c) at 4, 12, and 36 h
  after the CME flux-rope launch;
* ``earth_shock_thermal_source.csv`` contains the Earth-connected thermal
  energy-density curve plotted in panel (d).  Liu et al. state that this
  quantity is proportional to the number of particles injected at the shock,
  so XM03 uses only its time dependence as a publication-derived source.

The calibration constants below are the major tick/grid coordinates stored in
the authors' vector figure.  Keeping them explicit makes regeneration
reviewable: a changed upstream figure must fail the coordinate checks rather
than silently producing a shifted scientific baseline.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET


# Matplotlib colors embedded in the authors' Figure 12 vector PDF.  Only the
# three Earth instruments are selected; the orange STEREO-B observations are
# outside this one-field-line validation case.
EARTH_INSTRUMENT_COLORS = {
    "ACE/EPAM": "rgb(12.156677%, 55.685425%, 95.292664%)",
    "GOES-13/EPEAD": "rgb(0%, 25.097656%, 100%)",
    "SOHO/ERNE": "rgb(0%, 54.508972%, 54.508972%)",
}
EARTH_SOURCE_COLOR = "rgb(55.293274%, 82.743835%, 78.038025%)"

# Each spectrum panel has the same logarithmic axes.  Values are vector-space
# coordinates before pdftocairo applies the final PDF-to-SVG page transform.
PANEL_X_TICK_ZERO = (130.888348, 375.596970, 620.305591)  # E=1 MeV
X_PER_ENERGY_DECADE = 74.022178
Y_AT_UNIT_INTENSITY = 306.805421
Y_PER_INTENSITY_DECADE = 19.596962
# These panel labels are elapsed time after the 07:24 UTC flux-rope launch.
# They intentionally differ from the panel-(d) source's 06:00 UTC axis origin.
SNAPSHOT_HOURS = (4.0, 12.0, 36.0)

# Panel (d) spans 2013-04-11 06:00 through 2013-04-13 02:00.  Its vertical
# grid is two decades per major interval, starting at 1e5 keV m^-3.
SOURCE_X_AT_ZERO_HOURS = 58.139886
SOURCE_X_PER_HOUR = (188.182311 - 84.149154) / 8.0
SOURCE_Y_AT_1E5 = 52.923504
SOURCE_Y_PER_TWO_DECADES = 84.368654 - 52.923504


def _simple_segment(path_data: str):
    """Return an SVG M/L segment, or ``None`` for curves and compound paths."""
    match = re.fullmatch(
        r"M ([\-\d.]+) ([\-\d.]+) L ([\-\d.]+) ([\-\d.]+) ",
        path_data,
    )
    return tuple(map(float, match.groups())) if match else None


def _spectrum_energy(panel: int, x: float) -> float:
    return 10.0 ** ((x - PANEL_X_TICK_ZERO[panel]) / X_PER_ENERGY_DECADE)


def _spectrum_intensity(y: float) -> float:
    return 10.0 ** ((y - Y_AT_UNIT_INTENSITY) / Y_PER_INTENSITY_DECADE)


def _extract_observations(root: ET.Element):
    """Extract centers and horizontal energy bounds from vector error bars.

    In the original Matplotlib PDF, each measured point has a horizontal line
    whose endpoints are the instrument energy-bin limits.  The center of that
    segment is the marker coordinate.  Selecting simple, horizontal segments
    avoids confusing the measurements with smooth model curves of the same
    color.
    """
    observations = []
    for instrument, color in EARTH_INSTRUMENT_COLORS.items():
        for element in root.iter():
            if not element.tag.endswith("path") or element.get("stroke") != color:
                continue
            segment = _simple_segment(element.get("d", ""))
            if segment is None:
                continue
            x_left, y_left, x_right, y_right = segment
            if abs(y_left - y_right) > 1.0e-6 or x_right - x_left <= 1.0:
                continue
            # The three spectrum axes occupy y=224.5..389.1 in the vector
            # coordinate system.  This excludes legend samples and panel (d).
            if not 224.5 <= y_left <= 389.2:
                continue
            center = 0.5 * (x_left + x_right)
            panel = 0 if center < 290.0 else (1 if center < 535.0 else 2)
            observations.append({
                "elapsed_hours": SNAPSHOT_HOURS[panel],
                "instrument": instrument,
                "energy_low_mev": _spectrum_energy(panel, x_left),
                "energy_high_mev": _spectrum_energy(panel, x_right),
                "effective_energy_mev": _spectrum_energy(panel, center),
                "differential_intensity_pfu_per_mev": _spectrum_intensity(y_left),
            })
    observations.sort(key=lambda row: (
        row["elapsed_hours"], row["effective_energy_mev"], row["instrument"]))
    # The vector source contains 80 Earth measurements: 15/15/19 SOHO,
    # 2/6/7 ACE, and 4/4/4 GOES points across the three snapshots.
    if len(observations) != 80:
        raise RuntimeError(
            f"Figure 12 Earth observation extraction found {len(observations)} "
            "points instead of the reviewed 80")
    return observations


def _path_points(path_data: str):
    """Read the M/L vertices used by the smooth panel-(d) line."""
    tokens = re.findall(r"[ML]|[-+]?(?:\d+(?:\.\d*)?|\.\d+)", path_data)
    points = []
    index = 0
    while index < len(tokens):
        command = tokens[index]
        if command not in ("M", "L") or index + 2 >= len(tokens):
            return []
        points.append((float(tokens[index + 1]), float(tokens[index + 2])))
        index += 3
    return points


def _extract_earth_source(root: ET.Element):
    """Return the longest Earth-colored polyline in the lower source panel."""
    candidates = []
    for element in root.iter():
        if (element.tag.endswith("path") and
                element.get("stroke") == EARTH_SOURCE_COLOR):
            points = _path_points(element.get("d", ""))
            if points and max(y for _, y in points) < 183.1:
                candidates.append(points)
    if not candidates:
        raise RuntimeError("cannot find the Earth shock-source curve in Figure 12(d)")
    points = max(candidates, key=len)
    source = []
    for x, y in points:
        elapsed_hours = (x - SOURCE_X_AT_ZERO_HOURS) / SOURCE_X_PER_HOUR
        log10_density = 5.0 + 2.0 * (
            y - SOURCE_Y_AT_1E5) / SOURCE_Y_PER_TWO_DECADES
        if elapsed_hours >= 0.0 and math.isfinite(log10_density):
            source.append({
                "elapsed_hours": elapsed_hours,
                "thermal_energy_density_kev_per_m3": 10.0 ** log10_density,
            })
    if len(source) < 20:
        raise RuntimeError("Figure 12(d) Earth source curve is unexpectedly short")
    return source


def _write_csv(path: Path, fieldnames, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: f"{row[name]:.10g}" if isinstance(row[name], float)
                             else row[name] for name in fieldnames})
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-pdf", required=True, type=Path,
                        help="Fig/1304_Fig12_Spectrum_V8.pdf from arXiv:2412.07581v2")
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--earth-source", required=True, type=Path)
    arguments = parser.parse_args()
    if not arguments.figure_pdf.is_file():
        parser.error(f"Figure PDF does not exist: {arguments.figure_pdf}")

    with tempfile.TemporaryDirectory(prefix="srcsep-xm03-digitize-") as directory:
        svg = Path(directory) / "figure12.svg"
        completed = subprocess.run(
            ["pdftocairo", "-svg", str(arguments.figure_pdf), str(svg)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False)
        if completed.returncode != 0 or not svg.is_file():
            raise RuntimeError(f"pdftocairo could not convert Figure 12: {completed.stdout}")
        root = ET.parse(svg).getroot()

    observations = _extract_observations(root)
    source = _extract_earth_source(root)
    _write_csv(arguments.observations, (
        "elapsed_hours", "instrument", "energy_low_mev", "energy_high_mev",
        "effective_energy_mev", "differential_intensity_pfu_per_mev"), observations)
    _write_csv(arguments.earth_source, (
        "elapsed_hours", "thermal_energy_density_kev_per_m3"), source)
    print(f"wrote {len(observations)} Earth observations to {arguments.observations}")
    print(f"wrote {len(source)} Earth shock-source samples to {arguments.earth_source}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
