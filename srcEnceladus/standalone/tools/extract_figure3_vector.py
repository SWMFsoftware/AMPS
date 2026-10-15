#!/usr/bin/env python3
"""Recover Figure 3 curve centerlines from the Tenishev et al. (2014) PDF.

This is vector-path recovery, not scientific redigitization from a raster.  The
published strokes were converted to narrow filled polygons by the PDF writer;
the centerline is the mean of the two polygon edges at common x coordinates.
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")
RED = "fill:rgb(100%,0%,0%)"
BLUE = "fill:rgb(0%,32.156372%,100%)"

AXES = {
    "E3": {
        "x": (132.488281, 295.167969, -400.0, 400.0),
        "y_ticks": ([245.834, 216.1385, 186.451, 156.7595, 127.0725, 97.3825],
                    [9, 10, 11, 12, 13, 14]),
    },
    "E5": {
        "x": (343.160156, 506.214844, -400.0, 400.0),
        "y_ticks": ([245.7795, 214.1815, 182.592, 150.998, 119.408],
                    [10, 11, 12, 13, 14]),
    },
}


def path_points(path_data):
    values = np.asarray([float(v) for v in NUMBER.findall(path_data)], dtype=float)
    if values.size < 4 or values.size % 2:
        raise ValueError("unexpected SVG path coordinate count")
    return values.reshape(-1, 2)


def _unique_mean(x, y):
    order = np.argsort(x)
    x, y = np.asarray(x)[order], np.asarray(y)[order]
    unique, inverse = np.unique(x, return_inverse=True)
    sums = np.bincount(inverse, weights=y)
    return unique, sums / np.bincount(inverse)


def polygon_centerline(points):
    """Return a sampled centerline for a narrow, x-monotone stroke polygon."""
    split = int(np.argmax(points[:, 0]))
    if split < 1 or split >= len(points) - 2:
        raise ValueError("cannot split stroke polygon into two edges")
    xa, ya = _unique_mean(points[:split + 1, 0], points[:split + 1, 1])
    xb, yb = _unique_mean(points[split + 1:, 0], points[split + 1:, 1])
    lo, hi = max(xa.min(), xb.min()), min(xa.max(), xb.max())
    keep = (xa >= lo) & (xa <= hi)
    xa, ya = xa[keep], ya[keep]
    if xa.size < 2:
        raise ValueError("stroke polygon edges do not overlap")
    return np.column_stack((xa, 0.5 * (ya + np.interp(xa, xb, yb))))


def _style(element):
    return element.attrib.get("style", "").replace(" ", "")


def _main_colored_paths(root, color):
    candidates = []
    for element in root.iter():
        if not element.tag.endswith("path") or color not in _style(element):
            continue
        points = path_points(element.attrib["d"])
        if np.ptp(points[:, 0]) > 100.0:
            candidates.append(points)
    if len(candidates) != 2:
        raise RuntimeError("expected two main %s paths, found %d" %
                           (color, len(candidates)))
    return candidates


def _black_paths_by_clip(root):
    curves = {"E3": [], "E5": []}
    for group in root.iter():
        clip = group.attrib.get("clip-path", "")
        flyby = "E3" if clip in ("url(#clip1)", "url(#clip2)") else (
            "E5" if clip == "url(#clip3)" else None)
        if flyby is None:
            continue
        for element in group.iter():
            if (element.tag.endswith("path") and
                    "fill:rgb(0%,0%,0%)" in _style(element)):
                points = path_points(element.attrib["d"])
                if np.ptp(points[:, 0]) > 10.0:
                    curves[flyby].append(points)
    if len(curves["E3"]) != 2 or len(curves["E5"]) != 1:
        raise RuntimeError("unexpected clipped observation path count: %r" %
                           {k: len(v) for k, v in curves.items()})
    return curves


def to_physical(flyby, centerline):
    axis = AXES[flyby]
    x0, x1, t0, t1 = axis["x"]
    time = t0 + (centerline[:, 0] - x0) * (t1 - t0) / (x1 - x0)
    pixel_y, log_ticks = axis["y_ticks"]
    slope, intercept = np.polyfit(pixel_y, log_ticks, 1)
    density = np.power(10.0, slope * centerline[:, 1] + intercept)
    return time, density, slope, intercept


def recover(svg_path):
    root = ET.parse(svg_path).getroot()
    black = _black_paths_by_clip(root)
    recovered = []
    calibrations = {}
    for name, color in (("with_tiger_stripes", RED),
                        ("without_tiger_stripes", BLUE)):
        for points in _main_colored_paths(root, color):
            flyby = "E3" if np.mean(points[:, 0]) < 320 else "E5"
            time, density, slope, intercept = to_physical(
                flyby, polygon_centerline(points))
            calibrations[flyby] = {"log10_density_slope": slope,
                                   "log10_density_intercept": intercept}
            recovered.extend((flyby, name, float(t), float(n))
                             for t, n in zip(time, density))
    for flyby, polygons in black.items():
        segments = [polygon_centerline(p) for p in polygons]
        centerline = np.vstack(segments)
        centerline = centerline[np.argsort(centerline[:, 0])]
        time, density, slope, intercept = to_physical(flyby, centerline)
        calibrations[flyby] = {"log10_density_slope": slope,
                               "log10_density_intercept": intercept}
        recovered.extend((flyby, "INMS", float(t), float(n))
                         for t, n in zip(time, density))
    return sorted(recovered, key=lambda row: (row[0], row[1], row[2])), calibrations


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paper_pdf", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("output_metadata", type=Path)
    args = parser.parse_args()
    for output in (args.output_csv, args.output_metadata):
        if output.exists():
            parser.error("refusing to overwrite %s" % output)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="enceladus-fig3-") as directory:
        svg = Path(directory) / "page5.svg"
        subprocess.run(["pdftocairo", "-f", "5", "-l", "5", "-svg",
                        str(args.paper_pdf), str(svg)], check=True)
        rows, calibration = recover(svg)
    with open(args.output_csv, "x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["flyby", "series", "time_ca_s", "density_m3"])
        writer.writerows((f, s, "%.6f" % t, "%.9e" % n)
                         for f, s, t, n in rows)
    metadata = {
        "schema_version": 1,
        "source": "Tenishev et al. (2014), Figure 3",
        "doi": "10.1002/2014JE004700",
        "paper_pdf": str(args.paper_pdf),
        "paper_pdf_sha256": sha256(args.paper_pdf),
        "method": "centerlines recovered from filled vector stroke polygons on PDF page 5",
        "interpretation": {
            "INMS": "black published observational trace",
            "with_tiger_stripes": "red published 2014 model trace",
            "without_tiger_stripes": "blue published 2014 model trace"
        },
        "axis_calibration": calibration,
        "limitations": [
            "This reproduces the curves drawn in the paper; it is not a rerun of the 2014 simulation.",
            "Finite PDF stroke width and vector control-point sampling limit numerical precision."
        ]
    }
    with open(args.output_metadata, "x") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
