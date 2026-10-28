#!/usr/bin/env python3
"""Plot a physical quantity from one or more AMPS DIM1 output files.

AMPS writes radial results as whitespace-separated Tecplot tables. The first
``VARIABLES=`` record contains quoted column names and the remaining numeric
records contain one radial sample per row. This tool reads those names rather
than relying on fixed column numbers, because enabled AMPS diagnostics can
change the output layout between builds.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np


# Short names intentionally map to the exact labels written by AMPS. Keeping
# the mapping in one table makes command-line behavior and --help agree. In a
# DIM1 spherical run, X and Vx are the cometocentric radial coordinate and
# radial velocity; speed is the scalar magnitude written as |V|.
QUANTITY_ALIASES = {
    "r": "X",
    "level": "Maximum Refinment Level",
    "id": "Temp_ID",
    "rank": "Thread",
    "dens": "Number Density[1/m^3]",
    "part": "Particle Number",
    "vel": "Vx [m/s]",
    "vx": "Vx [m/s]",
    "vy": "Vy [m/s]",
    "vz": "Vz [m/s]",
    "speed": "|V| [m/s]",
    "temp": "Translational Temperature [K]",
    "tpar": "Parallel Translational Temperature [K]",
    "ttan": "Tangential Translational Temperature [K]",
    "dt": "Local Time Step",
    "weight": "log10(Local Particle Weight)",
}

# Collision columns depend on the compiled species list. These aliases cover
# every species in the current generic comet configuration and remain harmless
# when a smaller build omits one of the corresponding output columns.
COLLISION_SPECIES = (
    "H", "H2", "C", "CH", "CH2", "N", "CH3", "NH", "CH4", "NH2",
    "O", "NH3", "OH", "H2O", "C2", "C2H", "C2H2", "CN", "HCN",
    "C2H4", "CO", "C2H5", "HCO", "C2H6", "H2CO", "CH3OH", "CO2",
)
for _species in COLLISION_SPECIES:
    QUANTITY_ALIASES[f"coll-{_species.lower()}"] = (
        f"Coll Freq({_species}) [m^{{-3}} s^{{-1}}]"
    )


class PlotError(RuntimeError):
    """Report an input or selection problem without a Python traceback."""


def parse_variable_names(line: str, path: Path) -> List[str]:
    """Return quoted Tecplot variable names from a VARIABLES record."""

    names = re.findall(r'"([^"]+)"', line)
    if not names:
        raise PlotError(
            f"{path}: VARIABLES record does not contain quoted column names"
        )
    return names


def read_amps_table(path: Path) -> Tuple[List[str], np.ndarray]:
    """Read one AMPS/Tecplot table and validate every numeric row.

    TITLE, ZONE, and AUXDATA records are accepted so that the same reader also
    works with DIM1 files that include additional Tecplot metadata. A row with
    the wrong number of values is rejected instead of silently shifting a
    requested physical quantity into the wrong column.
    """

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise PlotError(f"cannot read {path}: {exc}") from exc

    variables = None
    rows = []
    metadata_prefixes = ("TITLE", "ZONE", "AUXDATA")

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        upper_line = line.upper()
        if upper_line.startswith("VARIABLES"):
            variables = parse_variable_names(line, path)
            continue
        if upper_line.startswith(metadata_prefixes):
            continue
        if variables is None:
            raise PlotError(
                f"{path}:{line_number}: numeric data precede the VARIABLES record"
            )

        fields = line.split()
        if len(fields) != len(variables):
            raise PlotError(
                f"{path}:{line_number}: found {len(fields)} values; "
                f"expected {len(variables)}"
            )
        try:
            rows.append([float(field) for field in fields])
        except ValueError as exc:
            raise PlotError(
                f"{path}:{line_number}: record is not numeric"
            ) from exc

    if variables is None:
        raise PlotError(f"{path}: VARIABLES record was not found")
    if not rows:
        raise PlotError(f"{path}: no numeric data records were found")

    return variables, np.asarray(rows, dtype=float)


def normalize_quantity_name(name: str) -> str:
    """Normalize a CLI name while retaining species text such as H2O or CO2."""

    return "".join(character.lower() for character in name if character.isalnum())


def resolve_quantity(selector: str, variables: Sequence[str], path: Path) -> int:
    """Resolve an exact name or an unambiguous abbreviated column name.

    Exact matching is case-insensitive. For convenience, punctuation and
    spaces may be omitted and a unique prefix may be used, so
    ``--quantity 'Number Density'`` selects ``Number Density[1/m^3]``. Vague
    selectors such as ``Temperature`` are rejected when several columns match.
    """

    requested_selector = selector
    selector = QUANTITY_ALIASES.get(selector.strip().casefold(), selector)
    lowered_selector = selector.casefold()
    exact = [
        index
        for index, variable in enumerate(variables)
        if variable.casefold() == lowered_selector
    ]
    if len(exact) == 1:
        return exact[0]

    normalized_selector = normalize_quantity_name(selector)
    if not normalized_selector:
        raise PlotError("quantity name is empty after removing punctuation")

    normalized_exact = [
        index
        for index, variable in enumerate(variables)
        if normalize_quantity_name(variable) == normalized_selector
    ]
    if len(normalized_exact) == 1:
        return normalized_exact[0]

    prefix_matches = [
        index
        for index, variable in enumerate(variables)
        if normalize_quantity_name(variable).startswith(normalized_selector)
    ]
    if len(prefix_matches) == 1:
        return prefix_matches[0]
    if len(prefix_matches) > 1:
        choices = ", ".join(variables[index] for index in prefix_matches)
        raise PlotError(
            f"{path}: quantity {requested_selector!r} is ambiguous; "
            f"matches: {choices}"
        )

    # Users commonly omit leading qualifiers such as "Translational" and
    # enter "Temperature". Accept that form only when it identifies one
    # column; otherwise report every candidate so no physical quantity is
    # chosen implicitly.
    substring_matches = [
        index
        for index, variable in enumerate(variables)
        if normalized_selector in normalize_quantity_name(variable)
    ]
    if len(substring_matches) == 1:
        return substring_matches[0]
    if len(substring_matches) > 1:
        choices = ", ".join(variables[index] for index in substring_matches)
        raise PlotError(
            f"{path}: quantity {requested_selector!r} is ambiguous; "
            f"matches: {choices}"
        )

    raise PlotError(
        f"{path}: quantity {requested_selector!r} was not found; "
        "use --list-quantities to inspect the available names"
    )


def list_quantities(paths: Sequence[Path]) -> None:
    """Print the available named columns for each requested file."""

    for file_index, path in enumerate(paths):
        variables, _ = read_amps_table(path)
        if file_index:
            print()
        print(f"{path}:")
        for column_index, variable in enumerate(variables, start=1):
            aliases = [
                alias
                for alias, target in QUANTITY_ALIASES.items()
                if target.casefold() == variable.casefold()
            ]
            alias_note = f"  [alias: {', '.join(aliases)}]" if aliases else ""
            print(f"  {column_index:2d}: {variable}{alias_note}")


def help_epilog_text() -> str:
    """Build the alias reference and copyable examples shown by --help."""

    regular_aliases = [
        ("r", "radius X"),
        ("level", "AMR refinement level"),
        ("id", "temporary mesh ID"),
        ("rank", "MPI rank/thread"),
        ("dens", "number density"),
        ("part", "model-particle number"),
        ("vel, vx", "radial velocity Vx"),
        ("vy", "Y velocity"),
        ("vz", "Z velocity"),
        ("speed", "speed magnitude |V|"),
        ("temp", "translational temperature"),
        ("tpar", "parallel translational temperature"),
        ("ttan", "tangential translational temperature"),
        ("dt", "local timestep"),
        ("weight", "local particle-weight column"),
    ]
    lines = ["quantity aliases:"]
    lines.extend(f"  {alias:<12} {description}" for alias, description in regular_aliases)
    collision_aliases = [f"coll-{species.lower()}" for species in COLLISION_SPECIES]
    lines.append("  coll-SPECIES collision frequency for the named species")
    lines.append("  collision aliases:")
    for start in range(0, len(collision_aliases), 7):
        lines.append("    " + ", ".join(collision_aliases[start:start + 7]))
    lines.extend(
        (
            "",
            "examples:",
            "  List quantities and their aliases:",
            "    python3 srcHartley2/visualization/plot_1d.py \\",
            "      PT/plots/pic.H2O.s=13.out=0.dat --list-quantities",
            "",
            "  Plot water density, save PNG/EPS, and open the viewer:",
            "    python3 srcHartley2/visualization/plot_1d.py \\",
            "      PT/plots/pic.H2O.s=13.out=0.dat \\",
            "      -q dens --yscale log --radius-unit km --rmax 1000 \\",
            "      --output PT/plots/h2o-density",
            "",
            "  Compare temperature from two runs:",
            "    python3 srcHartley2/visualization/plot_1d.py \\",
            "      run-a/pic.H2O.s=13.out=0.dat run-b/pic.H2O.s=13.out=4.dat \\",
            "      -q temp --yscale linear --output temperature-comparison",
            "",
            "  Save density on a batch node without opening a viewer:",
            "    python3 srcHartley2/visualization/plot_1d.py \\",
            "      PT/plots/pic.CO2.s=26.out=0.dat \\",
            "      -q dens --yscale log --rrange 6000 100000 \\",
            "      --output PT/plots/co2-density --no-show",
        )
    )
    return "\n".join(lines)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plot one physical quantity from one or more AMPS DIM1 output "
            "files. The radius axis is always logarithmic."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=help_epilog_text(),
    )
    parser.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="AMPS Tecplot output files to draw on the same axes",
    )
    parser.add_argument(
        "-q",
        "--quantity",
        help=(
            "column name or short alias to plot; exact names and unambiguous "
            "prefixes are accepted"
        ),
    )
    parser.add_argument(
        "--yscale",
        choices=("linear", "log", "symlog"),
        help="Y-axis scale; required when creating a plot",
    )
    parser.add_argument(
        "--linthresh",
        type=float,
        default=1.0,
        help="linear range around zero for --yscale symlog (default: 1.0)",
    )
    parser.add_argument(
        "--radius-unit",
        choices=("m", "km"),
        default="m",
        help="unit used for the cometocentric radius axis (default: m)",
    )
    radius_limits = parser.add_mutually_exclusive_group()
    radius_limits.add_argument(
        "--max-distance",
        "--rmax",
        type=float,
        help=(
            "maximum cometocentric distance to plot, in --radius-unit; "
            "--rmax is the short form"
        ),
    )
    radius_limits.add_argument(
        "--distance-range",
        "--rrange",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help=(
            "inclusive cometocentric distance range to plot, in "
            "--radius-unit; --rrange is the short form"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "output basename or file; PNG and EPS are always written, and a "
            "different requested extension is written as well"
        ),
    )
    parser.add_argument("--title", help="optional plot title")
    parser.add_argument(
        "--list-quantities",
        action="store_true",
        help="print available columns and exit without plotting",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="resolution of a saved raster figure (default: 150)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="save without opening a plot viewer",
    )
    return parser


def display_is_available() -> bool:
    """Return whether this process has a plausible graphical display."""

    if sys.platform.startswith("linux"):
        return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    return True


def requested_output_paths(output: Path) -> List[Path]:
    """Return the requested format plus mandatory PNG and EPS paths."""

    base = output.with_suffix("") if output.suffix else output
    paths = []
    if output.suffix:
        paths.append(output)
    paths.extend((base.with_suffix(".png"), base.with_suffix(".eps")))

    # Preserve order while removing duplicates when --output already names a
    # PNG or EPS file.
    unique_paths = []
    for path in paths:
        if path not in unique_paths:
            unique_paths.append(path)
    return unique_paths


def open_saved_plot(path: Path) -> None:
    """Open a saved PNG in the platform's desktop image viewer.

    The viewer is detached so closing the plotting script does not close the
    image window and the script does not wait for the user to dismiss it.
    Failure to find or start a viewer is reported without discarding the plots
    that were already written successfully.
    """

    try:
        if sys.platform.startswith("linux"):
            opener = shutil.which("xdg-open")
            if opener is None:
                raise OSError("xdg-open was not found")
            subprocess.Popen(
                [opener, str(path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
        elif sys.platform == "darwin":
            subprocess.Popen(
                ["open", str(path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
        elif os.name == "nt":
            os.startfile(str(path))  # type: ignore[attr-defined]
        else:
            raise OSError(f"no viewer launcher is defined for {sys.platform}")
    except OSError as exc:
        print(f"could not open plot viewer for {path}: {exc}", file=sys.stderr)
        return
    print(f"opened viewer for {path}")


def plot_files(args: argparse.Namespace) -> None:
    """Load all requested files and draw their selected profiles together."""

    if args.no_show and not args.output:
        raise PlotError("--no-show requires --output because no figure would be visible")

    show_viewer = not args.no_show
    if show_viewer and not display_is_available():
        print(
            "graphical display is unavailable; plots will be saved without "
            "opening a viewer",
            file=sys.stderr,
        )
        show_viewer = False

    # Saved plots use a non-interactive rendering backend and are displayed by
    # the desktop image viewer after both formats have been written. Without an
    # output path, retain matplotlib's interactive viewer behavior.
    if args.output or not show_viewer:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if args.linthresh <= 0.0 or not np.isfinite(args.linthresh):
        raise PlotError("--linthresh must be a positive finite value")
    if args.dpi <= 0:
        raise PlotError("--dpi must be positive")

    requested_minimum_radius = None
    requested_maximum_radius = None
    if args.max_distance is not None:
        if not np.isfinite(args.max_distance) or args.max_distance <= 0.0:
            raise PlotError("--max-distance/--rmax must be positive and finite")
        requested_maximum_radius = args.max_distance
    elif args.distance_range is not None:
        requested_minimum_radius, requested_maximum_radius = args.distance_range
        if (
            not np.isfinite(requested_minimum_radius)
            or not np.isfinite(requested_maximum_radius)
            or requested_minimum_radius <= 0.0
            or requested_maximum_radius <= requested_minimum_radius
        ):
            raise PlotError(
                "--distance-range/--rrange requires positive finite MIN and "
                "MAX with MIN < MAX"
            )

    figure, axes = plt.subplots(figsize=(8.0, 5.2))
    first_quantity_name = None
    radius_scale = 1.0 if args.radius_unit == "m" else 1000.0
    minimum_radius = np.inf
    maximum_radius = -np.inf
    plotted_file_count = 0

    for path in args.files:
        variables, values = read_amps_table(path)
        try:
            radius_index = variables.index("X")
        except ValueError as exc:
            raise PlotError(f"{path}: required radial coordinate column X is absent") from exc

        quantity_index = resolve_quantity(args.quantity, variables, path)
        quantity_name = variables[quantity_index]
        if first_quantity_name is None:
            first_quantity_name = quantity_name

        radius = values[:, radius_index] / radius_scale
        quantity = values[:, quantity_index]
        valid_radius = np.isfinite(radius) & (radius > 0.0)
        if requested_minimum_radius is not None:
            valid_radius &= radius >= requested_minimum_radius
        if requested_maximum_radius is not None:
            valid_radius &= radius <= requested_maximum_radius
        if not np.any(valid_radius):
            raise PlotError(
                f"{path}: no positive finite radial coordinates are inside "
                "the requested cometocentric distance limits"
            )
        minimum_radius = min(minimum_radius, float(np.min(radius[valid_radius])))
        maximum_radius = max(maximum_radius, float(np.max(radius[valid_radius])))
        valid_quantity = np.isfinite(quantity)

        # A logarithmic ordinate cannot represent zero or negative samples.
        # Mask them explicitly and report how many points were omitted rather
        # than relying on backend-specific matplotlib clipping behavior.
        nonpositive_count = 0
        if args.yscale == "log":
            nonpositive = valid_radius & valid_quantity & (quantity <= 0.0)
            nonpositive_count = int(np.count_nonzero(nonpositive))
            valid_quantity &= quantity > 0.0

        plottable = valid_radius & valid_quantity
        if not np.any(plottable):
            # A chemistry network commonly produces species that are absent in
            # an early sampling window. When several files are requested, an
            # all-zero log-Y profile must not prevent populated species from
            # being compared. Report and omit the empty curve; a global check
            # below still rejects a request in which every file is empty.
            print(
                f"{path}: skipped file because no plottable {quantity_name} "
                "points remain in the requested distance range",
                file=sys.stderr,
            )
            continue

        if nonpositive_count:
            print(
                f"{path}: skipped {nonpositive_count} non-positive "
                f"{quantity_name} points for logarithmic Y scale",
                file=sys.stderr,
            )

        # Retain selected radial coordinates and replace invalid ordinates with
        # NaN. Matplotlib omits those points and breaks the line at each gap;
        # deleting the rows would connect positive samples across a skipped
        # zero or negative logarithmic-Y value.
        order = np.argsort(radius[valid_radius])
        plot_radius = radius[valid_radius][order]
        plot_quantity = quantity[valid_radius][order].copy()
        selected_quantity_is_valid = valid_quantity[valid_radius][order]
        plot_quantity[~selected_quantity_is_valid] = np.nan
        axes.plot(
            plot_radius,
            plot_quantity,
            linewidth=1.5,
            label=path.name,
        )
        plotted_file_count += 1

    if plotted_file_count == 0:
        raise PlotError(
            f"none of the {len(args.files)} input files contains plottable "
            f"{args.quantity!r} data in the requested distance range"
        )

    axes.set_xscale("log")
    # Preserve the computational extent even if logarithmic-Y masking removes
    # zero-valued outer cells from a curve. Otherwise a young coma can appear
    # to fill an artificially shortened radial domain.
    plot_minimum_radius = (
        requested_minimum_radius
        if requested_minimum_radius is not None
        else minimum_radius
    )
    plot_maximum_radius = (
        requested_maximum_radius
        if requested_maximum_radius is not None
        else maximum_radius
    )
    if plot_minimum_radius < plot_maximum_radius:
        axes.set_xlim(plot_minimum_radius, plot_maximum_radius)
    if args.yscale == "symlog":
        axes.set_yscale("symlog", linthresh=args.linthresh)
    else:
        axes.set_yscale(args.yscale)
    axes.set_xlabel(f"Cometocentric radius [{args.radius_unit}]")
    axes.set_ylabel(first_quantity_name)
    if args.title:
        axes.set_title(args.title)
    # Opaque styling renders consistently in both PNG and PostScript output.
    axes.grid(True, which="both", color="0.85")
    axes.legend(title="AMPS output file", framealpha=1.0)
    figure.tight_layout()

    if args.output:
        try:
            output_paths = requested_output_paths(args.output)
            for output_path in output_paths:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                figure.savefig(output_path, dpi=args.dpi)
        except (OSError, ValueError) as exc:
            raise PlotError(f"cannot write {args.output}: {exc}") from exc
        for output_path in output_paths:
            print(f"wrote {output_path}")

    if args.output:
        if show_viewer:
            png_path = next(
                path for path in output_paths if path.suffix.lower() == ".png"
            )
            open_saved_plot(png_path)
        plt.close(figure)
    elif show_viewer:
        plt.show()
    else:
        plt.close(figure)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        if args.list_quantities:
            list_quantities(args.files)
            return 0
        if not args.quantity:
            parser.error("--quantity is required unless --list-quantities is used")
        if not args.yscale:
            parser.error("--yscale is required when creating a plot")
        plot_files(args)
    except PlotError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
