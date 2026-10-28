#!/usr/bin/env python3
"""Create neutral-coma plots comparable to Rubin et al. (2011), Figs. 1--18.

The paper's Figs. 2--17 contain both neutral DSMC profiles and simulated or
measured Giotto NMS mass-channel signals.  The generic comet application
currently writes neutral-species properties only.  This script therefore
reproduces the paper's species groupings, distance axes, density units, and
line styles for the neutral curves; it does not manufacture an NMS response.

Density is converted from the AMPS SI output (m^-3) to the paper's cm^-3.
Figure 18 uses the sampled radial water velocity in m/s.  A CSV table at fixed
cometocentric radii accompanies the figures so comparisons need not be made by
reading values from raster plots.
"""

from __future__ import annotations

import argparse
import csv
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


def density_figure(
    paper_figure: PaperFigure,
    profiles: Dict[str, Dict[str, np.ndarray]],
    output_index: int,
    iteration: int,
    radius_range: Tuple[float, float],
) -> plt.Figure:
    """Draw the neutral subset of one density figure on paper-compatible axes."""

    figure, axes = plt.subplots(figsize=(7.2, 5.0))
    populated = 0
    for curve in paper_figure.curves:
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
            # Methanol is intentionally disabled in the current Halley input.
            # Keep the missing curve visible as an annotation instead of
            # connecting zeros on logarithmic axes or silently omitting it.
            axes.text(
                0.02,
                0.04 + 0.06 * populated,
                f"{curve.species}: no positive samples in selected window",
                transform=axes.transAxes,
                fontsize=9,
            )
            continue
        axes.plot(
            radius[valid],
            density[valid],
            color=f"C{populated}",
            linestyle=curve.linestyle,
            linewidth=1.7,
            label=curve.species,
        )
        populated += 1

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(*radius_range)
    axes.set_xlabel("Cometocentric distance [km]")
    axes.set_ylabel(r"Neutral number density [cm$^{-3}$]")
    axes.set_title(
        f"Rubin et al. (2011) Fig. {paper_figure.number}: AMPS neutrals\n"
        f"output {output_index}, iteration {iteration}"
    )
    axes.grid(True, which="both", color="0.85")
    if populated:
        axes.legend(framealpha=1.0)
    figure.tight_layout()
    return figure


def velocity_figure(
    profiles: Dict[str, Dict[str, np.ndarray]],
    output_index: int,
    iteration: int,
    radius_range: Tuple[float, float],
) -> plt.Figure:
    """Draw the modeled water radial bulk speed corresponding to paper Fig. 18."""

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
    axes.plot(radius[valid], velocity[valid], "-", color="C0", linewidth=1.7, label="H2O")
    axes.set_xscale("log")
    axes.set_xlim(*radius_range)
    axes.set_xlabel("Cometocentric distance [km]")
    axes.set_ylabel("Radial bulk velocity [m/s]")
    axes.set_title(
        "Rubin et al. (2011) Fig. 18: AMPS water velocity\n"
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
    radii_km: Iterable[float],
) -> None:
    """Write numeric profiles at common paper-reading distances."""

    species_order = required_species()
    fieldnames = ["radius_km"] + [f"{species}_density_cm-3" for species in species_order]
    fieldnames.append("H2O_radial_velocity_m_s-1")
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
                row[f"{species}_density_cm-3"] = "" if math.isnan(value) else f"{value:.9e}"
            water = profiles["H2O"]
            velocity = finite_log_radius_interpolate(
                water["radius_km"], water["radial_velocity_ms"], radius_km
            )
            row["H2O_radial_velocity_m_s-1"] = (
                "" if math.isnan(velocity) else f"{velocity:.9e}"
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


def write_manifest(
    output_dir: Path,
    plots_dir: Path,
    output_index: int,
    iteration: int,
    files: Dict[str, Path],
    convergence_line: Optional[str],
    density_range: Tuple[float, float],
    velocity_range: Tuple[float, float],
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
        "convergence_record": convergence_line,
        "source_files": {species: path.name for species, path in sorted(files.items())},
        "limitations": [
            "Figs. 2-17 in the paper also show measured and modeled Giotto NMS mass-channel signals; this output contains neutral AMPS profiles only.",
            "The current halley1d.input intentionally disables methanol, dust, and extended sources.",
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
        iteration = (output_index + 1) * args.sampling_length

        pdf_path = output_dir / "rubin2011_neutral_profiles.pdf"
        with PdfPages(pdf_path) as pdf:
            for paper_figure in PAPER_DENSITY_FIGURES:
                figure = density_figure(
                    paper_figure,
                    profiles,
                    output_index,
                    iteration,
                    density_range,
                )
                base = output_dir / f"rubin2011_fig{paper_figure.number:02d}_neutral"
                save_figure(figure, base, args.dpi)
                pdf.savefig(figure)
                plt.close(figure)

            figure = velocity_figure(
                profiles, output_index, iteration, velocity_range
            )
            save_figure(figure, output_dir / "rubin2011_fig18_h2o_velocity", args.dpi)
            pdf.savefig(figure)
            plt.close(figure)

        write_comparison_table(
            output_dir / "comparison_profiles.csv", profiles, sorted(set(args.table_radii))
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
