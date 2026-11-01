#!/usr/bin/env python3
"""Validate and plot AMPS Tecplot flyby-density zones."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import re


@dataclass
class Zone:
    output: int
    iterations: int
    sample_steps: int
    rows: list[list[float]]


@dataclass(frozen=True)
class ObservationSeries:
    """One independently calibrated INMS observable on the AMPS time axis."""

    species_column: str
    label: str
    times_reference_s: list[float]
    densities_m3: list[float]
    point_sigma_m3: list[float]
    zero_count: list[bool]
    caveat: str


def _finite_float(row: dict[str, str], column: str, row_number: int) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as exception:
        raise ValueError(
            f"observation row {row_number} has invalid {column!r}") from exception
    if not math.isfinite(value):
        raise ValueError(
            f"observation row {row_number} has non-finite {column!r}")
    return value


def _read_observation_rows(path: Path,
                           required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        columns = set(reader.fieldnames or [])
        missing = required-columns
        if missing:
            raise ValueError(
                f"{path} is missing observation columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} has no observation rows")
    return rows


def read_inms_h2o(path: Path,
                  closest_approach_offset_s: float) -> ObservationSeries:
    """Read the released E3 sticking-corrected INMS H2O profile.

    The release time is measured from closest approach, whereas AMPS records
    time from the declared E3 UTC reference.  If closest approach occurs at
    ``offset`` on the AMPS axis, the exact coordinate conversion is
    ``t_reference = t_closest_approach + offset``.
    """
    if not math.isfinite(closest_approach_offset_s):
        raise ValueError("closest-approach offset must be finite")
    required = {"flyby", "time_ca_s", "density_m3", "sigma_m3"}
    rows = _read_observation_rows(path, required)
    times: list[float] = []
    densities: list[float] = []
    sigmas: list[float] = []
    for row_number, row in enumerate(rows, start=2):
        if row["flyby"] != "E3":
            raise ValueError(
                f"observation row {row_number} is not an E3 sample")
        time = (_finite_float(row, "time_ca_s", row_number) +
                closest_approach_offset_s)
        density = _finite_float(row, "density_m3", row_number)
        sigma = _finite_float(row, "sigma_m3", row_number)
        if density <= 0.0 or sigma <= 0.0:
            raise ValueError(
                f"observation row {row_number} requires positive density/sigma")
        times.append(time)
        densities.append(density)
        sigmas.append(sigma)
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("H2O observation times must increase strictly")
    return ObservationSeries(
        species_column="n_H2O_m-3",
        label="INMS-derived H2O (sticking corrected)",
        times_reference_s=times,
        densities_m3=densities,
        point_sigma_m3=sigmas,
        zero_count=[False]*len(times),
        caveat=("H2O bars are independent digitization uncertainty; the 23% "
                "common density-scale uncertainty is not shown."))


def read_inms_co2(path: Path,
                  closest_approach_offset_s: float) -> ObservationSeries:
    """Read the released E3 high-speed INMS mass-44 observable.

    This product is deliberately called *apparent* CO2.  At E3 encounter
    speed, impact fragmentation of heavy grains can contribute mass-44 counts,
    so displaying the calibrated channel as composition-resolved CO2 without
    that qualifier would overstate what INMS measured.
    """
    if not math.isfinite(closest_approach_offset_s):
        raise ValueError("closest-approach offset must be finite")
    required = {"flyby", "time_ca_s", "mass_per_charge_u",
                "species_interpretation", "density_m3",
                "sigma_counting_m3", "counts_per_integration"}
    rows = _read_observation_rows(path, required)
    times: list[float] = []
    densities: list[float] = []
    sigmas: list[float] = []
    zero_count: list[bool] = []
    for row_number, row in enumerate(rows, start=2):
        if row["flyby"] != "E3":
            raise ValueError(
                f"observation row {row_number} is not an E3 sample")
        if row["species_interpretation"] != "apparent_CO2":
            raise ValueError(
                f"observation row {row_number} is not apparent_CO2")
        mass = _finite_float(row, "mass_per_charge_u", row_number)
        if mass != 44.0:
            raise ValueError(
                f"observation row {row_number} is not mass 44")
        time = (_finite_float(row, "time_ca_s", row_number) +
                closest_approach_offset_s)
        density = _finite_float(row, "density_m3", row_number)
        sigma = _finite_float(row, "sigma_counting_m3", row_number)
        counts = _finite_float(row, "counts_per_integration", row_number)
        if (density < 0.0 or sigma <= 0.0 or counts < 0.0 or
                not counts.is_integer()):
            raise ValueError(
                f"observation row {row_number} has invalid count/density/sigma")
        if (counts == 0.0) != (density == 0.0):
            raise ValueError(
                f"observation row {row_number} has inconsistent zero count")
        times.append(time)
        densities.append(density)
        sigmas.append(sigma)
        zero_count.append(counts == 0.0)
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("CO2 observation times must increase strictly")
    return ObservationSeries(
        species_column="n_CO2_m-3",
        label="INMS m/z=44 (apparent CO2)",
        times_reference_s=times,
        densities_m3=densities,
        point_sigma_m3=sigmas,
        zero_count=zero_count,
        caveat=("Apparent-CO2 bars: Poisson counting uncertainty. Common 15% "
                "calibration and 2-5 impact-fragmentation range are not shown."))


def read_tecplot(path: Path) -> tuple[list[str], list[Zone]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    variables: list[str] = []
    zones: list[Zone] = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if line.startswith("VARIABLES="):
            variables = re.findall(r'"([^"]+)"', line)
        elif line.startswith("ZONE "):
            header = re.fullmatch(
                r'ZONE T="output=(\d+), iterations=(\d+), '
                r'sample_steps=(\d+)" I=(\d+) DATAPACKING=POINT '
                r'STRANDID=1 SOLUTIONTIME=([0-9.eE+-]+)', line)
            if header is None:
                raise ValueError(f"unrecognized Tecplot zone header: {line}")
            row_count = int(header.group(4))
            rows: list[list[float]] = []
            for row_index in range(row_count):
                index += 1
                if index >= len(lines):
                    raise ValueError("Tecplot zone ends before its declared I count")
                values = [float(value) for value in lines[index].split()]
                if variables and len(values) != len(variables):
                    raise ValueError(
                        f"zone row {row_index} has {len(values)} values, "
                        f"expected {len(variables)}")
                rows.append(values)
            zones.append(Zone(int(header.group(1)), int(header.group(2)),
                              int(header.group(3)), rows))
        index += 1
    if not variables:
        raise ValueError("Tecplot VARIABLES record is missing")
    if not zones:
        raise ValueError("Tecplot file has no zones")
    return variables, zones


def plot_density(tecplot: Path, png: Path,
                 observations: list[ObservationSeries] | None = None,
                 closest_approach_offset_s: float | None = None
                 ) -> tuple[int, int]:
    # Import lazily so validation of the Tecplot artifact remains usable in
    # environments where the optional rendering dependency is unavailable.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variables, zones = read_tecplot(tecplot)
    required = {"time_from_reference_s", "in_domain"}
    if not required.issubset(variables):
        raise ValueError(f"missing required columns: {sorted(required-set(variables))}")
    species_columns = [name for name in variables
                       if name.startswith("n_") and name.endswith("_m-3") and
                       name != "n_total_primary_m-3"]
    if not species_columns:
        raise ValueError("no species density columns are present")
    time_column = variables.index("time_from_reference_s")
    domain_column = variables.index("in_domain")
    in_domain_count = sum(
        row[domain_column] == 1.0 for row in zones[-1].rows)
    covered_times = [row[time_column] for row in zones[-1].rows
                     if row[domain_column] == 1.0]
    if not covered_times:
        raise ValueError("the final Tecplot zone has no in-domain trajectory rows")
    coverage_min = min(covered_times)
    coverage_max = max(covered_times)
    observation_map = {
        series.species_column: series for series in (observations or [])}
    if len(observation_map) != len(observations or []):
        raise ValueError("only one observation series per species may be plotted")
    unknown_observations = set(observation_map)-set(species_columns)
    if unknown_observations:
        raise ValueError(
            f"observations have no modeled species column: "
            f"{sorted(unknown_observations)}")

    figure, axes = plt.subplots(len(species_columns), 1,
                                figsize=(12.0, 4.2*len(species_columns)),
                                sharex=True)
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.13, top=0.86,
                           hspace=0.065)
    if len(species_columns) == 1:
        axes = [axes]
    colors = plt.cm.viridis_r(
        [index/max(1, len(zones)-1) for index in range(len(zones))])
    for axis, column_name in zip(axes, species_columns):
        density_column = variables.index(column_name)
        any_positive = False
        for zone_index, zone in enumerate(zones):
            selected = [(row[time_column], row[density_column])
                        for row in zone.rows if row[domain_column] == 1.0]
            times = [item[0] for item in selected]
            densities = [item[1] if item[1] > 0.0 else float("nan")
                         for item in selected]
            any_positive = any_positive or any(value > 0.0
                                                for _, value in selected)
            label = (f"window {zone.output}: steps "
                     f"{zone.iterations-zone.sample_steps+1}-{zone.iterations}")
            axis.plot(times, densities, color=colors[zone_index],
                      linewidth=2.2 if zone_index == len(zones)-1 else 0.9,
                      alpha=1.0 if zone_index == len(zones)-1 else 0.55,
                      label=label)
        observation = observation_map.get(column_name)
        if observation is not None:
            selected_observations = [
                (time, density, sigma, is_zero)
                for time, density, sigma, is_zero in zip(
                    observation.times_reference_s,
                    observation.densities_m3,
                    observation.point_sigma_m3,
                    observation.zero_count)
                if coverage_min <= time <= coverage_max]
            if not selected_observations:
                raise ValueError(
                    f"{observation.label} has no samples in the AMR overlap")
            positive = [item for item in selected_observations if not item[3]]
            if positive:
                # A one-count Poisson interval reaches zero.  Zero cannot be
                # represented on a logarithmic axis, so omit only that lower
                # segment while retaining the measured point and upper bar.
                # The figure caveat calls out this rendering decision.
                lower_errors = [
                    item[2] if item[1]-item[2] > 0.0 else 0.0
                    for item in positive]
                axis.errorbar(
                    [item[0] for item in positive],
                    [item[1] for item in positive],
                    yerr=[lower_errors, [item[2] for item in positive]],
                    fmt="o", markersize=4.2, color="black", ecolor="black",
                    elinewidth=0.9, capsize=2.0, zorder=20,
                    label=(f"{observation.label}, overlap "
                           f"(N={len(selected_observations)})"))
            zero = [item for item in selected_observations if item[3]]
            if zero:
                # A zero-count sample has no positive ordinate on a logarithmic
                # density axis.  Plot the released one-count sigma as an open
                # downward triangle: a visual upper-limit marker, not a
                # fabricated positive measurement.
                axis.scatter(
                    [item[0] for item in zero],
                    [item[2] for item in zero], marker="v", s=32,
                    facecolors="none", edgecolors="black", linewidths=0.9,
                    zorder=21,
                    label="zero count (one-count upper-limit marker)")
        if any_positive:
            axis.set_yscale("log")
        species_label = {
            "n_H2O_m-3": r"n(H$_2$O) [m$^{-3}$]",
            "n_CO2_m-3": r"n(CO$_2$) [m$^{-3}$]",
        }.get(column_name, column_name.replace("_m-3", " [m$^{-3}$]"))
        axis.set_ylabel(species_label)
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7, ncol=2, loc="best")
    reference_label = "time from E3 reference UTC [s]"
    if closest_approach_offset_s is not None:
        reference_label += (f"; closest approach = "
                            f"{closest_approach_offset_s:+g} s")
    axes[-1].set_xlabel(reference_label)
    title_subject = "AMPS sampling windows"
    final_title_line = "Outside-domain rows are excluded"
    if observations:
        title_subject += " and INMS data"
        final_title_line += "; INMS points are clipped to the AMR overlap"
    figure.suptitle(
        f"Cassini E3 local density: {title_subject}\n"
        f"AMR overlap: {in_domain_count}/{len(zones[-1].rows)} trajectory "
        f"points, {coverage_min:g} to {coverage_max:g} s\n"
        f"{final_title_line}")
    if observations:
        caveat_lines = [series.caveat for series in observations]
        if any(series.species_column == "n_CO2_m-3"
               for series in observations):
            caveat_lines.append(
                "Log-axis rendering: lower Poisson bars reaching zero are "
                "omitted; open triangles mark zero counts at one-count "
                "sensitivity.")
        figure.text(
            0.01, 0.012, "\n".join(caveat_lines),
            ha="left", va="bottom", fontsize=7.5)
    png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png, dpi=180)
    plt.close(figure)
    return len(zones), in_domain_count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tecplot", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--inms-h2o", type=Path,
        help="released E3 sticking-corrected H2O observation CSV")
    parser.add_argument(
        "--inms-co2", type=Path,
        help="released E3 apparent-CO2/mass-44 observation CSV")
    parser.add_argument(
        "--closest-approach-offset-s", type=float,
        help=("E3 closest-approach time on the AMPS reference-UTC axis; "
              "required when an INMS CSV is supplied"))
    args = parser.parse_args()
    tecplot = args.tecplot.resolve()
    output = (args.output.resolve() if args.output else
              tecplot.with_suffix(".png"))
    if ((args.inms_h2o or args.inms_co2) and
            args.closest_approach_offset_s is None):
        parser.error("--closest-approach-offset-s is required with INMS data")
    observations: list[ObservationSeries] = []
    if args.inms_h2o:
        observations.append(read_inms_h2o(
            args.inms_h2o.resolve(), args.closest_approach_offset_s))
    if args.inms_co2:
        observations.append(read_inms_co2(
            args.inms_co2.resolve(), args.closest_approach_offset_s))
    zone_count, in_domain_count = plot_density(
        tecplot, output, observations,
        args.closest_approach_offset_s if observations else None)
    print(f"PASS zones={zone_count} in_domain_points={in_domain_count} "
          f"loaded_observation_rows="
          f"{sum(len(item.times_reference_s) for item in observations)} "
          f"png={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
