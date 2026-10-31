#!/usr/bin/env python3
"""Compare AMPS neutral-coma profiles with Rubin et al. (2011), Figs. 1--18.

The paper's neutral DSMC curves are digitized directly from vector paths in
the journal PDF by ``halley2011/digitize_rubin2011.py``.  Each output figure
overlays those published values with the selected AMPS sampling window.  The
red simulated/measured Giotto NMS mass-channel signals in Figs. 2--17 are a
different observable and are deliberately not presented as neutral density.

Density is converted from the AMPS SI output (m^-3) to the paper's cm^-3.
Figure 18 uses the sampled radial water velocity in m/s.  A CSV table at fixed
cometocentric radii accompanies the figures so comparisons need not be made by
reading values from raster plots.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

# This is a batch comparison utility.  Selecting Agg before importing pyplot
# makes it safe on compute nodes without X11 and keeps PNG/EPS output stable.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Reuse the validated Tecplot reader instead of duplicating assumptions about
# the AMPS column layout.  The script lives beside plot_1d.py, so this import
# works both from the source tree and through the executable shebang.
from plot_1d import PlotError, read_amps_table


@dataclass(frozen=True)
class Curve:
    """One neutral curve and the line style used for it in the paper."""

    species: str
    linestyle: str


@dataclass(frozen=True)
class PaperFigure:
    """Neutral-species subset represented in one Rubin et al. figure."""

    number: int
    curves: Tuple[Curve, ...]


# Captions in Rubin et al. (2011) identify the neutral curve styles.  Figs.
# 2--17 also contain NMS curves/points that are outside the present model and
# are deliberately absent here.
PAPER_DENSITY_FIGURES: Tuple[PaperFigure, ...] = (
    PaperFigure(1, (Curve("H", "-."), Curve("H2", ":"))),
    PaperFigure(2, (Curve("C", "-."),)),
    PaperFigure(3, (Curve("CH", "-."),)),
    PaperFigure(4, (Curve("CH2", "-."), Curve("N", ":"))),
    PaperFigure(5, (Curve("CH3", "-."), Curve("NH", ":"))),
    PaperFigure(6, (Curve("CH4", "-."), Curve("NH2", ":"), Curve("O", "--"))),
    PaperFigure(7, (Curve("NH3", ":"), Curve("OH", "--"))),
    PaperFigure(8, (Curve("H2O", "--"),)),
    PaperFigure(9, (Curve("C2", ":"),)),
    PaperFigure(10, (Curve("C2H", ":"),)),
    PaperFigure(11, (Curve("C2H2", ":"), Curve("CN", "-."))),
    PaperFigure(12, (Curve("HCN", "-."),)),
    PaperFigure(13, (Curve("CO", "-."), Curve("C2H4", ":"))),
    PaperFigure(14, (Curve("HCO", "-."), Curve("C2H5", ":"))),
    PaperFigure(15, (Curve("H2CO", "-."), Curve("C2H6", ":"))),
    PaperFigure(16, (Curve("CH3OH", "-."),)),
    PaperFigure(17, (Curve("CO2", "-."),)),
)

FILE_PATTERN = re.compile(
    r"^pic\.(?P<species>.+)\.s=(?P<species_index>\d+)\.out=(?P<output_index>\d+)\.dat$"
)
DEFAULT_RADII_KM = (10.0, 100.0, 1_000.0, 10_000.0, 100_000.0)
DEFAULT_PAPER_DATA = (
    Path(__file__).resolve().parents[1]
    / "halley2011"
    / "reference"
    / "rubin2011_neutral_profiles.csv"
)


def build_parser() -> argparse.ArgumentParser:
    """Define a reproducible command line for a completed or partial run."""

    parser = argparse.ArgumentParser(
        description=(
            "Plot AMPS neutral profiles using the species groupings and units "
            "of Rubin et al. (2011), Icarus 213, Figs. 1--18."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "example:\n"
            "  python3 srcHartley2/visualization/compare_rubin2011.py \\\n"
            "    PT/plots --run-log halley1d.run.log\n\n"
            "When --output-index is omitted, the latest index present for all "
            "27 configured species is used."
        ),
    )
    parser.add_argument("plots_dir", type=Path, help="directory containing pic.*.dat files")
    parser.add_argument(
        "--output-index",
        type=int,
        help="AMPS output-window index (default: latest complete index)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="destination (default: PLOTS_DIR/rubin2011_comparison)",
    )
    parser.add_argument(
        "--run-log",
        type=Path,
        help="optional run log from which convergence status is copied",
    )
    parser.add_argument(
        "--paper-data",
        type=Path,
        default=DEFAULT_PAPER_DATA,
        help=(
            "digitized Rubin neutral profiles (default: "
            "srcHartley2/halley2011/reference/rubin2011_neutral_profiles.csv)"
        ),
    )
    parser.add_argument(
        "--sampling-length",
        type=int,
        default=500,
        help="iterations per output window, used for labels (default: 500)",
    )
    parser.add_argument(
        "--density-range",
        type=float,
        nargs=2,
        metavar=("MIN_KM", "MAX_KM"),
        default=(10.0, 100_000.0),
        help="cometocentric range for Figs. 1--17 (default: 10 100000)",
    )
    parser.add_argument(
        "--velocity-range",
        type=float,
        nargs=2,
        metavar=("MIN_KM", "MAX_KM"),
        default=(10.0, 1_000_000.0),
        help="cometocentric range for Fig. 18 (default: 10 1000000)",
    )
    parser.add_argument(
        "--table-radii",
        type=float,
        nargs="+",
        default=DEFAULT_RADII_KM,
        metavar="KM",
        help="radii sampled in comparison_profiles.csv",
    )
    parser.add_argument("--dpi", type=int, default=180, help="PNG resolution (default: 180)")
    return parser


def validate_positive_range(name: str, values: Sequence[float]) -> Tuple[float, float]:
    """Validate an inclusive positive plotting interval."""

    lower, upper = values
    if not (math.isfinite(lower) and math.isfinite(upper) and 0.0 < lower < upper):
        raise PlotError(f"{name} requires finite positive MIN < MAX")
    return lower, upper


def discover_outputs(plots_dir: Path) -> Dict[int, Dict[str, Path]]:
    """Index AMPS files by output window and neutral-species name."""

    if not plots_dir.is_dir():
        raise PlotError(f"plot directory does not exist: {plots_dir}")

    outputs: Dict[int, Dict[str, Path]] = {}
    for path in sorted(plots_dir.iterdir()):
        match = FILE_PATTERN.match(path.name)
        if match is None:
            continue
        output_index = int(match.group("output_index"))
        species = match.group("species")
        species_table = outputs.setdefault(output_index, {})
        if species in species_table:
            raise PlotError(
                f"duplicate species {species!r} for output index {output_index}"
            )
        species_table[species] = path

    if not outputs:
        raise PlotError(f"no pic.SPECIES.s=N.out=M.dat files found in {plots_dir}")
    return outputs


def required_species() -> Tuple[str, ...]:
    """Return every species needed for the paper-ordered comparison."""

    result: List[str] = []
    for paper_figure in PAPER_DENSITY_FIGURES:
        for curve in paper_figure.curves:
            if curve.species not in result:
                result.append(curve.species)
    # Figure 18 needs H2O, already included by Figure 8.
    return tuple(result)


def select_output(
    outputs: Dict[int, Dict[str, Path]], requested_index: Optional[int]
) -> Tuple[int, Dict[str, Path]]:
    """Select a requested window or the latest one complete for all species."""

    expected = set(required_species())
    if requested_index is not None:
        if requested_index not in outputs:
            choices = ", ".join(str(index) for index in sorted(outputs))
            raise PlotError(
                f"output index {requested_index} is absent; available indices: {choices}"
            )
        missing = sorted(expected - set(outputs[requested_index]))
        if missing:
            raise PlotError(
                f"output index {requested_index} lacks species: {', '.join(missing)}"
            )
        return requested_index, outputs[requested_index]

    for output_index in sorted(outputs, reverse=True):
        if expected.issubset(outputs[output_index]):
            return output_index, outputs[output_index]
    raise PlotError("no output index contains every species required by Figs. 1--18")


def named_column(variables: Sequence[str], name: str, path: Path) -> int:
    """Return an exact AMPS column, with a useful format error if absent."""

    try:
        return variables.index(name)
    except ValueError as exc:
        raise PlotError(f"{path}: required column {name!r} is absent") from exc


def load_profiles(files: Dict[str, Path]) -> Dict[str, Dict[str, np.ndarray]]:
    """Load radius, density, and radial velocity from the selected window."""

    profiles: Dict[str, Dict[str, np.ndarray]] = {}
    for species in required_species():
        path = files[species]
        variables, values = read_amps_table(path)
        radius = values[:, named_column(variables, "X", path)] / 1000.0
        density = values[:, named_column(variables, "Number Density[1/m^3]", path)]
        radial_velocity = values[:, named_column(variables, "Vx [m/s]", path)]
        finite_radius = np.isfinite(radius) & (radius > 0.0)
        order = np.argsort(radius[finite_radius])
        profiles[species] = {
            "radius_km": radius[finite_radius][order],
            "density_cm3": density[finite_radius][order] * 1.0e-6,
            "radial_velocity_ms": radial_velocity[finite_radius][order],
        }
    return profiles


PaperKey = Tuple[int, str, str]


def load_paper_profiles(path: Path) -> Dict[PaperKey, Dict[str, np.ndarray]]:
    """Load the committed vector digitization of the Rubin DSMC curves."""

    required_columns = {
        "figure",
        "species",
        "variant",
        "quantity",
        "radius_km",
        "value",
        "unit",
    }
    rows: Dict[PaperKey, List[Tuple[float, float]]] = {}
    quantities: Dict[PaperKey, str] = {}
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
                missing = sorted(required_columns - set(reader.fieldnames or ()))
                raise PlotError(
                    f"{path}: paper-data columns are missing: {', '.join(missing)}"
                )
            for line_number, row in enumerate(reader, start=2):
                try:
                    figure = int(row["figure"])
                    species = row["species"].strip()
                    variant = row["variant"].strip()
                    quantity = row["quantity"].strip()
                    radius = float(row["radius_km"])
                    value = float(row["value"])
                except (TypeError, ValueError) as exc:
                    raise PlotError(
                        f"{path}:{line_number}: invalid digitized numeric value"
                    ) from exc
                if not species or not variant:
                    raise PlotError(
                        f"{path}:{line_number}: species and variant must be nonempty"
                    )
                if not (
                    math.isfinite(radius)
                    and radius > 0.0
                    and math.isfinite(value)
                    and value > 0.0
                ):
                    raise PlotError(
                        f"{path}:{line_number}: radius and value must be positive and finite"
                    )
                key = (figure, species, variant)
                if key in quantities and quantities[key] != quantity:
                    raise PlotError(f"{path}:{line_number}: quantity changes within {key}")
                quantities[key] = quantity
                rows.setdefault(key, []).append((radius, value))
    except OSError as exc:
        raise PlotError(f"cannot read digitized paper data {path}: {exc}") from exc

    result: Dict[PaperKey, Dict[str, np.ndarray]] = {}
    for key, samples in rows.items():
        samples.sort()
        radius = np.asarray([sample[0] for sample in samples], dtype=float)
        value = np.asarray([sample[1] for sample in samples], dtype=float)
        if len(radius) < 2 or np.any(np.diff(radius) <= 0.0):
            raise PlotError(f"{path}: {key} must contain increasing radii")
        result[key] = {
            "radius_km": radius,
            "value": value,
        }

    for paper_figure in PAPER_DENSITY_FIGURES:
        for curve in paper_figure.curves:
            key = (paper_figure.number, curve.species, "nominal")
            if key not in result:
                raise PlotError(f"{path}: digitized paper curve is missing: {key}")
            if quantities[key] != "number_density":
                raise PlotError(f"{path}: {key} is not a number-density curve")
    velocity_key = (18, "H2O", "nominal")
    if velocity_key not in result or quantities[velocity_key] != "radial_velocity":
        raise PlotError(f"{path}: digitized H2O velocity curve is missing")
    return result


def paper_variants(
    profiles: Dict[PaperKey, Dict[str, np.ndarray]], figure: int, species: str
) -> List[Tuple[str, Dict[str, np.ndarray]]]:
    """Return nominal first, followed by any explicitly published variants."""

    variants = [
        (variant, profile)
        for (candidate_figure, candidate_species, variant), profile in profiles.items()
        if candidate_figure == figure and candidate_species == species
    ]
    return sorted(variants, key=lambda item: (item[0] != "nominal", item[0]))


def density_figure(
    paper_figure: PaperFigure,
    profiles: Dict[str, Dict[str, np.ndarray]],
    published_profiles: Dict[PaperKey, Dict[str, np.ndarray]],
    output_index: int,
    iteration: int,
    radius_range: Tuple[float, float],
) -> plt.Figure:
    """Overlay AMPS and digitized Rubin neutral densities on common axes."""

    figure, axes = plt.subplots(figsize=(7.2, 5.0))
    populated = 0
    for curve_index, curve in enumerate(paper_figure.curves):
        color = f"C{curve_index % 10}"
        profile = profiles[curve.species]
        radius = profile["radius_km"]
        density = profile["density_cm3"]
        valid = (
            np.isfinite(density)
            & (density > 0.0)
            & (radius >= radius_range[0])
            & (radius <= radius_range[1])
        )
        if not np.any(valid):
            # Keep a statistically absent AMPS curve visible as an annotation
            # while still plotting the digitized paper reference below.
            axes.text(
                0.02,
                0.04 + 0.06 * curve_index,
                f"{curve.species}: no positive samples in selected window",
                transform=axes.transAxes,
                fontsize=9,
            )
        else:
            axes.plot(
                radius[valid],
                density[valid],
                color=color,
                linestyle="-",
                linewidth=1.8,
                label=f"AMPS {curve.species}",
            )
            populated += 1

        for variant, published in paper_variants(
            published_profiles, paper_figure.number, curve.species
        ):
            paper_radius = published["radius_km"]
            paper_density = published["value"]
            paper_valid = (
                np.isfinite(paper_density)
                & (paper_density > 0.0)
                & (paper_radius >= radius_range[0])
                & (paper_radius <= radius_range[1])
            )
            if not np.any(paper_valid):
                continue
            if variant == "nominal":
                label = f"Rubin DSMC {curve.species}"
                linestyle = curve.linestyle
                linewidth = 1.6
            else:
                label = f"Rubin DSMC {curve.species} (no CH3OH→H2CO)"
                linestyle = ":"
                linewidth = 1.3
            axes.plot(
                paper_radius[paper_valid],
                paper_density[paper_valid],
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label,
            )
            populated += 1

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(*radius_range)
    axes.set_xlabel("Cometocentric distance [km]")
    axes.set_ylabel(r"Neutral number density [cm$^{-3}$]")
    axes.set_title(
        f"Rubin et al. (2011) Fig. {paper_figure.number}: AMPS vs. published DSMC\n"
        f"output {output_index}, iteration {iteration}"
    )
    axes.grid(True, which="both", color="0.85")
    if populated:
        axes.legend(framealpha=1.0)
    figure.tight_layout()
    return figure


def velocity_figure(
    profiles: Dict[str, Dict[str, np.ndarray]],
    published_profiles: Dict[PaperKey, Dict[str, np.ndarray]],
    output_index: int,
    iteration: int,
    radius_range: Tuple[float, float],
) -> plt.Figure:
    """Overlay AMPS and digitized Rubin H2O radial bulk speeds."""

    profile = profiles["H2O"]
    radius = profile["radius_km"]
    velocity = profile["radial_velocity_ms"]
    valid = (
        np.isfinite(velocity)
        & (radius >= radius_range[0])
        & (radius <= radius_range[1])
    )
    if not np.any(valid):
        raise PlotError("H2O has no finite velocity samples in the requested range")

    figure, axes = plt.subplots(figsize=(7.2, 5.0))
    axes.plot(
        radius[valid],
        velocity[valid],
        "-",
        color="C0",
        linewidth=1.8,
        label="AMPS H2O",
    )
    published = published_profiles[(18, "H2O", "nominal")]
    paper_radius = published["radius_km"]
    paper_velocity = published["value"]
    paper_valid = (
        np.isfinite(paper_velocity)
        & (paper_velocity > 0.0)
        & (paper_radius >= radius_range[0])
        & (paper_radius <= radius_range[1])
    )
    axes.plot(
        paper_radius[paper_valid],
        paper_velocity[paper_valid],
        "--",
        color="C1",
        linewidth=1.6,
        label="Rubin DSMC H2O",
    )
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(*radius_range)
    axes.set_xlabel("Cometocentric distance [km]")
    axes.set_ylabel("Radial bulk velocity [m/s]")
    axes.set_title(
        "Rubin et al. (2011) Fig. 18: AMPS vs. published H2O velocity\n"
        f"output {output_index}, iteration {iteration}"
    )
    axes.grid(True, which="both", color="0.85")
    axes.legend(framealpha=1.0)
    figure.tight_layout()
    return figure


def save_figure(figure: plt.Figure, base: Path, dpi: int) -> None:
    """Write both the screen-friendly and publication-vector formats."""

    figure.savefig(base.with_suffix(".png"), dpi=dpi)
    figure.savefig(base.with_suffix(".eps"))


def positive_log_interpolate(radius: np.ndarray, values: np.ndarray, target: float) -> float:
    """Log-log interpolate a positive density profile without inventing zeros."""

    valid = np.isfinite(radius) & np.isfinite(values) & (radius > 0.0) & (values > 0.0)
    if np.count_nonzero(valid) < 2:
        return math.nan
    x = radius[valid]
    y = values[valid]
    if target < x[0] or target > x[-1]:
        return math.nan
    return float(np.exp(np.interp(np.log(target), np.log(x), np.log(y))))


def finite_log_radius_interpolate(radius: np.ndarray, values: np.ndarray, target: float) -> float:
    """Interpolate signed velocity linearly in log radius."""

    valid = np.isfinite(radius) & np.isfinite(values) & (radius > 0.0)
    if np.count_nonzero(valid) < 2:
        return math.nan
    x = radius[valid]
    y = values[valid]
    if target < x[0] or target > x[-1]:
        return math.nan
    return float(np.interp(np.log(target), np.log(x), y))


def write_comparison_table(
    output_path: Path,
    profiles: Dict[str, Dict[str, np.ndarray]],
    published_profiles: Dict[PaperKey, Dict[str, np.ndarray]],
    radii_km: Iterable[float],
) -> None:
    """Write AMPS and published profiles at common comparison radii."""

    species_order = required_species()
    species_figure = {
        curve.species: paper_figure.number
        for paper_figure in PAPER_DENSITY_FIGURES
        for curve in paper_figure.curves
    }
    fieldnames = ["radius_km"]
    for species in species_order:
        fieldnames.extend(
            [
                f"AMPS_{species}_density_cm-3",
                f"Rubin2011_{species}_density_cm-3",
            ]
        )
        if species == "H2CO":
            fieldnames.append(
                "Rubin2011_H2CO_no_ch3oh_photodissociation_density_cm-3"
            )
    fieldnames.extend(
        [
            "AMPS_H2O_radial_velocity_m_s-1",
            "Rubin2011_H2O_radial_velocity_m_s-1",
        ]
    )
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for radius_km in radii_km:
            row = {"radius_km": f"{radius_km:.9g}"}
            for species in species_order:
                profile = profiles[species]
                value = positive_log_interpolate(
                    profile["radius_km"], profile["density_cm3"], radius_km
                )
                row[f"AMPS_{species}_density_cm-3"] = (
                    "" if math.isnan(value) else f"{value:.9e}"
                )
                published = published_profiles[
                    (species_figure[species], species, "nominal")
                ]
                paper_value = positive_log_interpolate(
                    published["radius_km"], published["value"], radius_km
                )
                row[f"Rubin2011_{species}_density_cm-3"] = (
                    "" if math.isnan(paper_value) else f"{paper_value:.9e}"
                )
                if species == "H2CO":
                    alternative = published_profiles[
                        (15, "H2CO", "no_ch3oh_photodissociation")
                    ]
                    alternative_value = positive_log_interpolate(
                        alternative["radius_km"], alternative["value"], radius_km
                    )
                    row[
                        "Rubin2011_H2CO_no_ch3oh_photodissociation_density_cm-3"
                    ] = (
                        ""
                        if math.isnan(alternative_value)
                        else f"{alternative_value:.9e}"
                    )
            water = profiles["H2O"]
            velocity = finite_log_radius_interpolate(
                water["radius_km"], water["radial_velocity_ms"], radius_km
            )
            row["AMPS_H2O_radial_velocity_m_s-1"] = (
                "" if math.isnan(velocity) else f"{velocity:.9e}"
            )
            paper_water = published_profiles[(18, "H2O", "nominal")]
            paper_velocity = positive_log_interpolate(
                paper_water["radius_km"], paper_water["value"], radius_km
            )
            row["Rubin2011_H2O_radial_velocity_m_s-1"] = (
                "" if math.isnan(paper_velocity) else f"{paper_velocity:.9e}"
            )
            writer.writerow(row)


def convergence_line_for_iteration(
    run_log: Optional[Path], selected_iteration: int
) -> Optional[str]:
    """Extract the newest convergence record belonging to the selected output.

    A live run may already be sampling a later window while an earlier output
    is plotted explicitly.  Restricting the copied status to records at or
    before the selected iteration prevents a later transient/converged state
    from being attached to the wrong files.
    """

    if run_log is None:
        return None
    try:
        lines = run_log.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PlotError(f"cannot read run log {run_log}: {exc}") from exc
    matches: List[Tuple[int, str]] = []
    for line in lines:
        if "radial convergence" not in line:
            continue
        iteration_match = re.search(r"\biteration=(\d+)\b", line)
        if iteration_match is None:
            continue
        iteration = int(iteration_match.group(1))
        if iteration <= selected_iteration:
            matches.append((iteration, line.strip()))
    return max(matches, key=lambda item: item[0])[1] if matches else None


def file_sha256(path: Path) -> str:
    """Return a stable provenance digest for the digitized reference table."""

    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise PlotError(f"cannot hash digitized paper data {path}: {exc}") from exc
    return digest.hexdigest()


def write_manifest(
    output_dir: Path,
    plots_dir: Path,
    output_index: int,
    iteration: int,
    files: Dict[str, Path],
    convergence_line: Optional[str],
    density_range: Tuple[float, float],
    velocity_range: Tuple[float, float],
    paper_data: Path,
) -> None:
    """Record provenance and model-to-paper comparison limits."""

    manifest = {
        "paper": "Rubin et al., Icarus 213 (2011) 655-677",
        "paper_figures": "neutral-species subsets of Figs. 1-17 and H2O velocity in Fig. 18",
        "plots_directory": str(plots_dir.resolve()),
        "output_index": output_index,
        "nominal_iteration": iteration,
        "density_distance_range_km": list(density_range),
        "velocity_distance_range_km": list(velocity_range),
        "density_conversion": "AMPS m^-3 multiplied by 1e-6 to obtain cm^-3",
        "paper_data_file": str(paper_data.resolve()),
        "paper_data_sha256": file_sha256(paper_data),
        "paper_data_method": (
            "vector paths extracted from the Rubin et al. PDF with its "
            "piecewise logarithmic distance axes"
        ),
        "convergence_record": convergence_line,
        "source_files": {species: path.name for species, path in sorted(files.items())},
        "limitations": [
            "Figs. 2-17 also show Giotto NMS mass-channel observables; those are not neutral number densities and are excluded from this overlay.",
            "The current halley1d.input disables dust, extended sources, and the optional CH3OH-to-H2CO photodissociation branch.",
            "A selected output is scientifically steady only if the run log reports the required consecutive convergence windows.",
        ],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Load one complete AMPS window and create every comparison artifact."""

    args = build_parser().parse_args(argv)
    try:
        if args.sampling_length <= 0:
            raise PlotError("--sampling-length must be positive")
        if args.dpi <= 0:
            raise PlotError("--dpi must be positive")
        density_range = validate_positive_range("--density-range", args.density_range)
        velocity_range = validate_positive_range("--velocity-range", args.velocity_range)
        if any(not math.isfinite(radius) or radius <= 0.0 for radius in args.table_radii):
            raise PlotError("--table-radii values must be positive and finite")

        outputs = discover_outputs(args.plots_dir)
        output_index, files = select_output(outputs, args.output_index)
        output_dir = args.output_dir or args.plots_dir / "rubin2011_comparison"
        output_dir.mkdir(parents=True, exist_ok=True)
        profiles = load_profiles(files)
        published_profiles = load_paper_profiles(args.paper_data)
        iteration = (output_index + 1) * args.sampling_length

        pdf_path = output_dir / "rubin2011_neutral_profiles.pdf"
        with PdfPages(pdf_path) as pdf:
            for paper_figure in PAPER_DENSITY_FIGURES:
                figure = density_figure(
                    paper_figure,
                    profiles,
                    published_profiles,
                    output_index,
                    iteration,
                    density_range,
                )
                base = output_dir / f"rubin2011_fig{paper_figure.number:02d}_neutral"
                save_figure(figure, base, args.dpi)
                pdf.savefig(figure)
                plt.close(figure)

            figure = velocity_figure(
                profiles,
                published_profiles,
                output_index,
                iteration,
                velocity_range,
            )
            save_figure(figure, output_dir / "rubin2011_fig18_h2o_velocity", args.dpi)
            pdf.savefig(figure)
            plt.close(figure)

        write_comparison_table(
            output_dir / "comparison_profiles.csv",
            profiles,
            published_profiles,
            sorted(set(args.table_radii)),
        )
        convergence_line = convergence_line_for_iteration(args.run_log, iteration)
        write_manifest(
            output_dir,
            args.plots_dir,
            output_index,
            iteration,
            files,
            convergence_line,
            density_range,
            velocity_range,
            args.paper_data,
        )
    except PlotError as exc:
        print(f"compare_rubin2011.py: error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"compare_rubin2011.py: error: {exc}", file=sys.stderr)
        return 2

    print(f"wrote Rubin et al. comparison for output {output_index} to {output_dir}")
    print(f"combined PDF: {pdf_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
