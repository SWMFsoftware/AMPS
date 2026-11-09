#!/usr/bin/env python3
"""Overlay particle-resolved AMPS brightness on MAVEN/IUVS references.

The default observations are the reviewed numerical raster digitization in
``visualization/references/i3_brightness_MAVEN``. Runtime plotting therefore
does not depend on the external source-image archive. The archive, reviewed
pixel anchors, member hashes, and reproducible extraction utility remain
available for provenance and regeneration.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


DEFAULT_REFERENCE_DIR = (
    Path(__file__).resolve().parent / "references/i3_brightness_MAVEN"
)
REFERENCE_TABLE_NAME = "maven_iuvs_brightness.csv"
DEFAULT_X_LIMITS = (0.2, 20.0)  # 10^5 km; common source-raster limits
DEFAULT_Y_LIMITS = (0.09, 20.0)  # kR; common source-raster limits


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Overlay an AMPS H or OH brightness profile on the extracted "
            "3I/ATLAS MAVEN/IUVS references and save a PNG."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  # Case number and species are inferred from a standard run path/table.
  %(prog)s \
    validation_runs/i3_case1_run/production/PT/plots/pic.H.brightness.s=0.out=0.dat

  # Add the outer-boundary (edge) correction from the collisionless reference
  # model; the run directory is three levels above the table.
  %(prog)s \
    RUN/PT/plots/pic.H.brightness.s=0.out=6.dat --case 4 --edge-correction

  # Use the same CSV schema from another reference-data directory.
  %(prog)s PT/plots/pic.OH.brightness.s=3.out=0.dat \
    --case 1 --species OH --reference-dir /data/references/i3_MAVEN

The Tecplot input must contain Impact Parameter [m] and Brightness ... [kR].
The reference directory must contain maven_iuvs_brightness.csv. The default
table is committed under srcHartley2/visualization/references and does not
require Archive.zip at plotting time.
""",
    )
    parser.add_argument("brightness_file", type=Path, help="AMPS brightness Tecplot file")
    parser.add_argument(
        "--case",
        type=int,
        choices=range(1, 7),
        metavar="N",
        help="3I case number, 1--6 (otherwise inferred from the input path)",
    )
    parser.add_argument(
        "--species",
        choices=("H", "OH", "h", "oh"),
        help="brightness species (otherwise inferred from Tecplot variables/name)",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help=(
            "directory containing maven_iuvs_brightness.csv "
            f"(default: {DEFAULT_REFERENCE_DIR})"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "output PNG path (default: PT/plots/i3_caseN_SPECIES_brightness_comparison.png "
            "under the current working directory)"
        ),
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=DEFAULT_X_LIMITS[0],
        help="minimum distance in units of 10^5 km (default: 0.2)",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=DEFAULT_X_LIMITS[1],
        help="maximum distance in units of 10^5 km (default: 20)",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=DEFAULT_Y_LIMITS[1],
        help="comparison maximum brightness, kR (default: 20)",
    )
    parser.add_argument(
        "--ymin",
        type=float,
        help=(
            "comparison minimum brightness, kR; by default it is expanded "
            "up to two decades below the largest in-range AMPS value"
        ),
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG resolution (default: 180)")
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="save the PNG without opening an interactive Matplotlib window",
    )
    parser.add_argument(
        "--title",
        help="custom title (default identifies the case, species, and output index)",
    )
    edge = parser.add_argument_group(
        "edge correction",
        "Overlay AMPS x N_full/N_trunc from the collisionless reference model in "
        "comet_edge_correction.py (see srcHartley2/EDGE_CORRECTION.md).",
    )
    edge.add_argument(
        "--edge-correction",
        action="store_true",
        help="also plot the brightness corrected for the outer domain boundary",
    )
    edge.add_argument(
        "--run-dir",
        type=Path,
        help="run directory holding i3.in, i3.conf and PT/plots (default: three levels above the table)",
    )
    edge.add_argument(
        "--runtime-input",
        type=Path,
        help="runtime input for the correction (default: RUN/i3.in)",
    )
    edge.add_argument(
        "--build-input",
        type=Path,
        help="compile-time input for the correction (default: RUN/i3.conf, else input/i3.conf)",
    )
    edge.add_argument(
        "--edge-parents",
        type=int,
        default=3_000_000,
        help="parent test particles per source in the reference model (default: 3000000)",
    )
    edge.add_argument(
        "--edge-seed",
        type=int,
        default=20261109,
        help="random seed of the reference model (default: 20261109)",
    )
    return parser


def parse_variables(line: str) -> List[str]:
    """Parse the quoted comma-separated Tecplot VARIABLES record."""

    _, separator, value = line.partition("=")
    if not separator:
        raise ValueError("malformed Tecplot VARIABLES record")
    return [item.strip() for item in next(csv.reader([value], skipinitialspace=True))]


def read_brightness(path: Path) -> Tuple[List[str], np.ndarray, str]:
    """Read numeric point data and the title from an AMPS Tecplot table."""

    variables: List[str] = []
    rows: List[List[float]] = []
    title = ""
    with path.open("r", encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("TITLE="):
                title = line.partition("=")[2].strip().strip('"')
                continue
            if line.startswith("VARIABLES="):
                variables = parse_variables(line)
                continue
            if not variables or line.startswith(("DATASETAUXDATA", "ZONE")):
                continue
            fields = line.split()
            if len(fields) != len(variables):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(variables)} numeric "
                    f"columns, found {len(fields)}"
                )
            try:
                rows.append([float(field) for field in fields])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-numeric data row") from error

    if not variables:
        raise ValueError(f"{path}: VARIABLES record was not found")
    if not rows:
        raise ValueError(f"{path}: no Tecplot point data were found")
    return variables, np.asarray(rows, dtype=float), title


def find_column(variables: Sequence[str], exact: str, contains: Sequence[str]) -> int:
    """Find an exact column or the unique column containing every token."""

    lowered = [value.casefold() for value in variables]
    if exact.casefold() in lowered:
        return lowered.index(exact.casefold())
    matches = [
        index
        for index, value in enumerate(lowered)
        if all(token.casefold() in value for token in contains)
    ]
    if len(matches) != 1:
        descriptions = ", ".join(repr(value) for value in variables)
        raise ValueError(
            f"could not uniquely select {exact!r}; candidate columns are {descriptions}"
        )
    return matches[0]


def infer_species(path: Path, variables: Sequence[str], requested: str | None) -> str:
    if requested:
        return requested.upper()
    joined = " ".join(variables).upper()
    if "H_LYMAN_ALPHA" in joined or re.search(r"PIC\.H\.BRIGHTNESS", path.name, re.I):
        return "H"
    if "OH_A2SIGMA_X2PI_00" in joined or re.search(r"PIC\.OH\.BRIGHTNESS", path.name, re.I):
        return "OH"
    raise ValueError("species could not be inferred; supply --species H or --species OH")


def infer_case(path: Path, requested: int | None) -> int:
    if requested is not None:
        return requested
    match = re.search(
        r"(?:^|[^a-z0-9])(?:i3[_-]?)?case[_-]?([1-6])(?:[^0-9]|$)",
        str(path),
        re.I,
    )
    if not match:
        raise ValueError("case number could not be inferred from the path; supply --case 1 ... 6")
    return int(match.group(1))


def optional_reference_float(row: Dict[str, str], name: str) -> float:
    text = row[name].strip()
    return math.nan if text == "" else float(text)


def read_reference(
    reference_dir: Path, case_number: int, species: str
) -> Tuple[Path, Dict[str, np.ndarray]]:
    """Read and strictly validate one case/species reference subset."""

    table = reference_dir / REFERENCE_TABLE_NAME
    required = {
        "case_id",
        "species",
        "point_index",
        "distance_1e5_km",
        "distance_lower_1e5_km",
        "distance_upper_1e5_km",
        "brightness_kR",
        "brightness_lower_kR",
        "brightness_upper_kR",
    }
    selected: List[Tuple[int, List[float]]] = []
    with table.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{table}: missing required columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                if int(row["case_id"]) != case_number or row["species"].strip().upper() != species:
                    continue
                selected.append(
                    (
                        int(row["point_index"]),
                        [
                            float(row["distance_1e5_km"]),
                            float(row["distance_lower_1e5_km"]),
                            float(row["distance_upper_1e5_km"]),
                            float(row["brightness_kR"]),
                            optional_reference_float(row, "brightness_lower_kR"),
                            optional_reference_float(row, "brightness_upper_kR"),
                        ],
                    )
                )
            except (KeyError, ValueError) as error:
                raise ValueError(f"{table}:{line_number}: invalid reference row") from error

    if not selected:
        raise ValueError(f"{table}: no reference rows for case {case_number}, species {species}")
    selected.sort(key=lambda item: item[0])
    indices = [item[0] for item in selected]
    if len(indices) != len(set(indices)) or indices != list(range(1, len(indices) + 1)):
        raise ValueError(f"{table}: point indices must be unique and contiguous from one")

    values = np.asarray([item[1] for item in selected], dtype=float)
    x, x_lower, x_upper, y, y_lower, y_upper = values.T
    if not (
        np.all(np.isfinite(values[:, :4]))
        and np.all(x_lower > 0.0)
        and np.all(x_lower < x)
        and np.all(x < x_upper)
        and np.all(y > 0.0)
    ):
        raise ValueError(f"{table}: non-finite/non-positive or unordered reference coordinates")
    paired_uncertainty = np.isnan(y_lower) == np.isnan(y_upper)
    valid_uncertainty = np.isnan(y_lower) | (
        (y_lower > 0.0) & (y_lower < y) & (y < y_upper)
    )
    if not np.all(paired_uncertainty & valid_uncertainty):
        raise ValueError(f"{table}: invalid or partially specified vertical uncertainty")

    return table, {
        "x": x,
        "x_lower": x_lower,
        "x_upper": x_upper,
        "y": y,
        "y_lower": y_lower,
        "y_upper": y_upper,
    }


def default_ymin(model_y: np.ndarray) -> float:
    """Keep up to two useful decades of an under-bright model visible."""

    peak = float(np.max(model_y))
    expanded = 10.0 ** (math.floor(math.log10(peak)) - 2)
    return min(DEFAULT_Y_LIMITS[0], expanded)


def output_path(requested: Path | None, case_number: int, species: str) -> Path:
    if requested is not None:
        return requested if requested.suffix.lower() == ".png" else requested.with_suffix(".png")
    return Path.cwd() / "PT" / "plots" / (
        f"i3_case{case_number}_{species}_brightness_comparison.png"
    )


def plot_edge_correction(
    axis,
    arguments: argparse.Namespace,
    case_number: int,
    species: str,
    impact_m: np.ndarray,
    brightness: np.ndarray,
    observation: Dict[str, np.ndarray],
) -> str:
    """Plot AMPS x edge factor and return a one-line summary."""

    import importlib.util

    module_path = Path(__file__).resolve().parent / "comet_edge_correction.py"
    specification = importlib.util.spec_from_file_location("comet_edge_correction", module_path)
    edge = sys.modules.get(specification.name)
    if edge is None:
        edge = importlib.util.module_from_spec(specification)
        sys.modules[specification.name] = edge  # dataclasses resolve their module
        specification.loader.exec_module(edge)

    run_dir = arguments.run_dir or edge.run_directory(arguments.brightness_file)
    factors, inputs = edge.correction_for_run(
        run_dir,
        case_number,
        edge.output_index(arguments.brightness_file),
        impact_m,
        (species,),
        arguments.edge_parents,
        arguments.edge_seed,
        arguments.runtime_input,
        arguments.build_input,
    )
    factor = factors[species]
    corrected = brightness * factor
    shown = (
        np.isfinite(corrected)
        & (corrected > 0.0)
        & (impact_m <= 2.0 / 3.0 * inputs.outer_radius_m)
        & (impact_m / 1.0e8 >= arguments.xmin)
        & (impact_m / 1.0e8 <= arguments.xmax)
    )
    axis.plot(
        impact_m[shown] / 1.0e8,
        corrected[shown],
        color="#e65100",
        linewidth=1.8,
        linestyle=(0, (5, 2)),
        label=f"AMPS x edge correction (domain edge {inputs.outer_radius_m:.1e} m)",
        zorder=5,
    )
    axis.axvline(inputs.outer_radius_m / 1.0e8, color="0.55", linewidth=0.8, linestyle="--", zorder=1)
    ok = np.isfinite(factor) & (brightness > 0.0)
    log_b, log_f = np.log(impact_m[ok]), np.log(factor[ok])
    at_points = np.exp(np.interp(np.log(observation["x"] * 1.0e8), log_b, log_f))
    return (
        f"edge correction: domain edge {inputs.outer_radius_m:.3e} m; factor at the "
        f"observation points median {np.median(at_points):.4f}, range "
        f"{at_points.min():.4f}--{at_points.max():.4f}; plotted up to 2/3 of the edge"
    )


def plot_comparison(arguments: argparse.Namespace) -> Path:
    variables, data, tecplot_title = read_brightness(arguments.brightness_file)
    species = infer_species(arguments.brightness_file, variables, arguments.species)
    case_number = infer_case(arguments.brightness_file.resolve(), arguments.case)
    reference_table, reference = read_reference(
        arguments.reference_dir, case_number, species
    )

    x_column = find_column(variables, "Impact Parameter [m]", ("impact parameter", "[m]"))
    y_column = find_column(variables, "Brightness [kR]", ("brightness", "[kr]"))
    distance = data[:, x_column] / 1.0e8  # m -> units of 10^5 km
    brightness = data[:, y_column]
    finite_positive = (
        np.isfinite(distance)
        & np.isfinite(brightness)
        & (distance > 0.0)
        & (brightness > 0.0)
    )
    in_x_range = finite_positive & (distance >= arguments.xmin) & (distance <= arguments.xmax)
    skipped = int(data.shape[0] - np.count_nonzero(finite_positive))
    if skipped:
        print(
            f"{arguments.brightness_file}: skipped {skipped} non-positive/non-finite samples",
            file=sys.stderr,
        )
    if not np.any(in_x_range):
        raise ValueError("no positive AMPS brightness samples lie inside the requested X range")
    if not (0.0 < arguments.xmin < arguments.xmax):
        raise ValueError("require 0 < --xmin < --xmax")

    x = distance[in_x_range]
    y = brightness[in_x_range]
    ymin = arguments.ymin if arguments.ymin is not None else default_ymin(y)
    if ymin <= 0.0 or ymin >= arguments.ymax:
        raise ValueError("require 0 < --ymin < --ymax")

    reference_mask = (
        (reference["x"] >= arguments.xmin) & (reference["x"] <= arguments.xmax)
    )
    if not np.any(reference_mask):
        raise ValueError("no MAVEN/IUVS reference points lie inside the requested X range")
    observation = {name: values[reference_mask] for name, values in reference.items()}

    import matplotlib

    display_available = bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    show_figure = not arguments.no_show and display_available
    if not show_figure:
        matplotlib.use("Agg")
    if not arguments.no_show and not display_available:
        print(
            "no DISPLAY or WAYLAND_DISPLAY is available; saved PNG will not be shown "
            "(use an X-forwarded/desktop terminal to open the Matplotlib window)",
            file=sys.stderr,
        )

    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(8.2, 5.7), dpi=arguments.dpi)
    y_lower_error = np.where(
        np.isnan(observation["y_lower"]),
        0.0,
        observation["y"] - observation["y_lower"],
    )
    y_upper_error = np.where(
        np.isnan(observation["y_upper"]),
        0.0,
        observation["y_upper"] - observation["y"],
    )
    axis.errorbar(
        observation["x"],
        observation["y"],
        xerr=np.vstack(
            (
                observation["x"] - observation["x_lower"],
                observation["x_upper"] - observation["x"],
            )
        ),
        yerr=np.vstack((y_lower_error, y_upper_error)),
        fmt="o-",
        color="black",
        ecolor="black",
        linewidth=1.2,
        elinewidth=1.0,
        markersize=4.0,
        capsize=0,
        label="MAVEN/IUVS (raster digitization)",
        zorder=4,
    )
    axis.plot(
        x,
        y,
        color="#1565c0",
        linewidth=2.2,
        label="AMPS particle-resolved brightness",
        zorder=5,
    )
    edge_summary = None
    if arguments.edge_correction:
        edge_summary = plot_edge_correction(
            axis, arguments, case_number, species, data[:, x_column], brightness, observation
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(arguments.xmin, arguments.xmax)
    axis.set_ylim(ymin, arguments.ymax)
    axis.set_xlabel(r"Relative distance at comet ($10^5$ km)")
    species_label = r"H Ly-$\alpha$" if species == "H" else "OH(0,0)"
    axis.set_ylabel(f"{species_label} radiance (kR)")
    axis.tick_params(which="both", direction="in", top=True, right=True)

    output_match = re.search(r"output=(\d+)", tecplot_title)
    output_index = output_match.group(1) if output_match else "?"
    default_title = f"3I/ATLAS case {case_number}: {species} brightness (AMPS output {output_index})"
    axis.set_title(arguments.title or default_title)
    if hasattr(figure.canvas.manager, "set_window_title"):
        figure.canvas.manager.set_window_title(arguments.title or default_title)
    axis.legend(loc="lower left", fontsize=8, framealpha=0.92)

    if ymin < DEFAULT_Y_LIMITS[0] < arguments.ymax:
        axis.axhline(
            DEFAULT_Y_LIMITS[0], color="0.55", linewidth=0.8, linestyle=":", zorder=2
        )
        axis.text(
            0.98,
            DEFAULT_Y_LIMITS[0],
            "lower edge of source raster",
            transform=axis.get_yaxis_transform(),
            ha="right",
            va="top",
            color="0.35",
            fontsize=7,
        )

    destination = output_path(arguments.output, case_number, species)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(destination, dpi=arguments.dpi, format="png")

    print(f"observation reference: {reference_table.resolve()}")
    print(f"observation points: {observation['x'].size}")
    print(f"AMPS profile: {arguments.brightness_file}")
    print(f"plotted AMPS samples: {x.size}; Y range: {ymin:g}--{arguments.ymax:g} kR")
    if edge_summary is not None:
        print(edge_summary)
    print(f"wrote: {destination.resolve()}")
    if show_figure:
        print("opening interactive Matplotlib window; close the window to return to the shell")
        plt.show()
    plt.close(figure)
    return destination


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        plot_comparison(arguments)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
