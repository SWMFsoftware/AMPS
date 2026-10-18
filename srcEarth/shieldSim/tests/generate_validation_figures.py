#!/usr/bin/env python3
"""Render human-review figures from numerical validation-test evidence.

The numerical tests remain authoritative: this script runs only after them and
does not compare pixels or introduce new pass/fail tolerances.  Its outputs are
review aids written to the build tree.  Observation-only panels are labeled so
that an ingestion plot cannot be mistaken for model validation.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import textwrap

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
except ImportError as error:  # pragma: no cover - depends on the build environment
    raise SystemExit(f"validation figure generation requires matplotlib: {error}")


if len(sys.argv) != 5:
    raise SystemExit(
        "usage: generate_validation_figures.py SOURCE_ROOT EVIDENCE_DIR OUTPUT_DIR SOURCE_TOOL"
    )

ROOT = pathlib.Path(sys.argv[1]).resolve()
EVIDENCE = pathlib.Path(sys.argv[2]).resolve()
OUTPUT = pathlib.Path(sys.argv[3]).resolve()
SOURCE_TOOL = pathlib.Path(sys.argv[4]).resolve()
FIXTURES = ROOT / "tests" / "reference_spectra"

COLORS = {
    "blue": "#2166ac",
    "orange": "#d6604d",
    "green": "#1b9e77",
    "purple": "#7570b3",
    "gray": "#666666",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_evidence(name: str) -> dict:
    path = EVIDENCE / name
    require(path.is_file(), f"validation evidence is missing: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    require(
        payload.get("schema_version") == "shieldSim.validation_evidence.v1",
        f"unexpected evidence schema in {path}",
    )
    comparisons = payload.get("comparisons")
    require(isinstance(comparisons, list) and comparisons, f"empty evidence in {path}")
    for item in comparisons:
        for field in ("id", "unit", "expected", "actual"):
            require(field in item, f"{path} comparison lacks {field}")
        require(
            math.isfinite(float(item["expected"])) and math.isfinite(float(item["actual"])),
            f"non-finite comparison in {path}",
        )
    return payload


def wrap_label(value: str, width: int = 25) -> str:
    return textwrap.fill(value.replace("_", " "), width=width)


def save_figure(fig, stem: str, caption: str, entries: list[dict]) -> None:
    """Write both review-friendly PNG and vector SVG forms."""

    fig.text(0.01, 0.008, textwrap.fill(caption, width=145), fontsize=8, color="#333333")
    # Reserve explicit space for both the suptitle and the scientific-scope
    # caption.  Without this, dense multi-panel figures can overlap their
    # titles even though the raster file itself is technically valid.
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.92))
    for extension in ("png", "svg"):
        path = OUTPUT / f"{stem}.{extension}"
        kwargs = {"bbox_inches": "tight"}
        if extension == "png":
            kwargs["dpi"] = 170
        fig.savefig(path, **kwargs)
        require(path.stat().st_size > 1000, f"generated figure is unexpectedly small: {path}")
        entries.append(
            {
                "caption": caption,
                "file": path.name,
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    plt.close(fig)


def comparison_figure(evidence: dict, title: str, stem: str, caption: str,
                      entries: list[dict], selected: list[dict] | None = None) -> None:
    """Plot dimensionless expected/actual agreement and signed residuals."""

    comparisons = selected if selected is not None else evidence["comparisons"]
    require(bool(comparisons), f"no comparisons selected for {stem}")
    labels = [wrap_label(str(item["id"])) for item in comparisons]
    normalized_expected: list[float] = []
    normalized_actual: list[float] = []
    residual: list[float] = []
    annotations: list[str] = []
    for item in comparisons:
        expected = float(item["expected"])
        actual = float(item["actual"])
        scale = max(abs(expected), abs(actual))
        if scale == 0.0:
            normalized_expected.append(0.0)
            normalized_actual.append(0.0)
            residual.append(0.0)
        else:
            normalized_expected.append(expected / scale)
            normalized_actual.append(actual / scale)
            residual.append((actual - expected) / scale)
        annotations.append(f"{actual:.6g} {item['unit']}")

    height = max(5.2, 0.58 * len(comparisons) + 2.6)
    fig, (agreement, difference) = plt.subplots(
        1, 2, figsize=(13.5, height), gridspec_kw={"width_ratios": [1.25, 1.0]}
    )
    positions = list(range(len(comparisons)))
    offset = 0.18
    agreement.barh(
        [position + offset for position in positions], normalized_expected,
        height=0.32, color="#bdbdbd", label="analytic/reference",
    )
    agreement.barh(
        [position - offset for position in positions], normalized_actual,
        height=0.32, color=COLORS["blue"], label="computed/test value",
    )
    agreement.axvline(1.0, color="#333333", linewidth=0.8, linestyle=":")
    agreement.set_yticks(positions)
    agreement.set_yticklabels(labels, fontsize=8)
    agreement.invert_yaxis()
    agreement.set_xlabel("value normalized by max(|reference|, |computed|)")
    agreement.legend(loc="best", fontsize=8)
    agreement.grid(axis="x", alpha=0.25)
    for position, value, annotation in zip(positions, normalized_actual, annotations):
        agreement.text(value, position - offset, f"  {annotation}", va="center", fontsize=7)

    difference.barh(positions, residual, color=COLORS["orange"])
    difference.axvline(0.0, color="#222222", linewidth=0.9)
    difference.set_yticks(positions)
    difference.set_yticklabels([])
    difference.invert_yaxis()
    difference.set_xlabel("(computed - reference) / comparison scale")
    difference.grid(axis="x", alpha=0.25)
    max_residual = max((abs(value) for value in residual), default=0.0)
    difference.set_title(f"maximum normalized residual = {max_residual:.3g}", fontsize=9)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    save_figure(fig, stem, caption, entries)


def stage4_nist_figure(evidence: dict, entries: list[dict]) -> None:
    """Compactly show every registered NIST comparison and its own limit."""

    comparisons = evidence["comparisons"]
    labels = [str(item["id"]).replace("_stopping_power", "_S")
              .replace("_mass_attenuation_coefficient", "_mu")
              .replace("_csda_range", "_R") for item in comparisons]
    ratios = [float(item["relative_difference"]) / float(item["acceptance_relative"])
              for item in comparisons]
    colors = []
    for item in comparisons:
        identifier = str(item["id"])
        colors.append(COLORS["blue"] if "stopping" in identifier
                      else COLORS["green"] if "csda" in identifier
                      else COLORS["purple"])
    fig, (values, summary) = plt.subplots(
        2, 1, figsize=(15.0, 9.0), gridspec_kw={"height_ratios": [2.4, 1.0]}
    )
    positions = list(range(len(comparisons)))
    values.scatter(positions, ratios, c=colors, s=25, edgecolors="none")
    values.axhline(1.0, color=COLORS["orange"], linestyle="--", label="registered limit")
    values.set_ylabel("|Geant4 - NIST| / row-specific acceptance")
    values.set_xticks(positions)
    values.set_xticklabels(labels, rotation=75, ha="right", fontsize=6)
    values.grid(axis="y", alpha=0.25)
    values.legend(fontsize=8)
    quantities = ["electronic_stopping_power", "total_stopping_power", "csda_range",
                  "mass_attenuation_coefficient"]
    maxima = [max(float(item["relative_difference"]) for item in comparisons
                  if str(item["id"]).endswith(quantity)) for quantity in quantities]
    limits = [max(float(item["acceptance_relative"]) for item in comparisons
                  if str(item["id"]).endswith(quantity)) for quantity in quantities]
    x = list(range(len(quantities)))
    summary.bar([value - 0.18 for value in x], maxima, 0.36, color=COLORS["blue"],
                label="maximum observed relative difference")
    summary.bar([value + 0.18 for value in x], limits, 0.36, color="#bdbdbd",
                label="largest registered row limit")
    summary.set_xticks(x)
    summary.set_xticklabels(["electronic S", "total S", "CSDA range", "photon mu/rho"])
    summary.set_ylabel("relative difference")
    summary.grid(axis="y", alpha=0.25)
    summary.legend(fontsize=8, ncol=2)
    fig.suptitle("Stage 4 Geant4 electromagnetic tables vs pinned NIST references",
                 fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage4_nist_reference_comparisons",
        "All 63 exact-node PSTAR, ASTAR, and XCOM photon mass-attenuation comparisons are shown "
        "against their preregistered row-specific limits. Photon attenuation is not interpreted "
        "as absorbed dose. Numerical assertions are authoritative; formal Stage-4 sign-off still "
        "requires independent review.",
        entries,
    )


def stage4_transport_figure(evidence: dict, entries: list[dict]) -> None:
    """Show all analytic ledger/mass conversions plus numerical TID gates."""

    comparisons = evidence["comparisons"]
    residuals = []
    for item in comparisons:
        expected = float(item["expected"])
        actual = float(item["actual"])
        scale = max(abs(expected), abs(actual), 1.0e-300)
        residuals.append((actual - expected) / scale)
    fig, (exact, convergence) = plt.subplots(2, 1, figsize=(15.0, 9.2))
    positions = list(range(len(comparisons)))
    exact.scatter(positions, residuals, color=COLORS["blue"], s=24)
    exact.axhline(0.0, color="#222222", linewidth=0.9)
    exact.set_xticks(positions)
    exact.set_xticklabels([wrap_label(str(item["id"]), 18) for item in comparisons],
                          rotation=75, ha="right", fontsize=6)
    exact.set_ylabel("(computed - analytic) / scale")
    exact.set_title("Every analytic mass/dose conversion and ledger partition closure")
    exact.grid(axis="y", alpha=0.25)

    names: list[str] = []
    values: list[float] = []
    limits: list[float] = []
    for family, points in evidence["convergence"].items():
        names.append(f"{family}\nfinal refinement")
        values.append(float(points[-1]["final_refinement_relative_change"]))
        limits.append(0.02)
    names.append("EM option 4\nsensitivity")
    values.append(float(evidence["em_option4_transport_sensitivity_relative"]))
    limits.append(float(evidence["em_option4_transport_sensitivity_relative_max"]))
    names.append("max accepted\nMC standard u")
    accepted_mc = [float(item["statistical_relative_standard_uncertainty"])
                   for item in evidence["material_species_matrix"]]
    for points in evidence["convergence"].values():
        accepted_mc.extend(float(item["statistical_relative_standard_uncertainty"])
                           for item in points)
    values.append(max(accepted_mc))
    limits.append(float(evidence["monte_carlo_relative_standard_uncertainty_max"]))
    x = list(range(len(names)))
    convergence.bar([value - 0.18 for value in x], values, 0.36, color=COLORS["green"],
                    label="observed")
    convergence.bar([value + 0.18 for value in x], limits, 0.36, color="#bdbdbd",
                    label="registered limit")
    convergence.set_xticks(x)
    convergence.set_xticklabels(names, fontsize=8)
    convergence.set_ylabel("relative change or uncertainty")
    convergence.grid(axis="y", alpha=0.25)
    convergence.legend(fontsize=8)
    fig.suptitle("Stage 4 finite-slab TID scoring and convergence",
                 fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage4_tid_transport_comparisons",
        "The upper panel includes every exact analytic slab-mass, Gy conversion, and independent "
        "particle/ancestry/process ledger-closure case. The lower panel shows quantity-specific "
        "cut, step, lateral-escape, EM-model, and Monte Carlo gates. Lateral size is tested on "
        "deposited energy because changing slab area changes mass and Gy per compact-plane primary "
        "by definition.",
        entries,
    )


def run_source_tool(command: str, options: list[str], output: pathlib.Path) -> dict:
    completed = subprocess.run(
        [str(SOURCE_TOOL), command, *options, f"--output={output}"],
        text=True,
        capture_output=True,
        check=False,
    )
    require(completed.returncode == 0, f"source tool failed: {completed.stderr}")
    return json.loads(output.read_text(encoding="utf-8"))


def ams_figure(canonical: dict, entries: list[dict]) -> None:
    info = json.loads(
        (FIXTURES / "gcr" / "ams02_20131118_hapi_info.json").read_text(encoding="utf-8")
    )
    with (FIXTURES / "gcr" / "ams02_20131118_hapi_data.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        row = next(csv.reader(stream))
    parameters = {item["name"]: item for item in info["parameters"]}
    raw = {
        "proton": {
            "x": [1000.0 * float(value) for value in parameters["ProFlux"]["bins"][0]["centers"]],
            "y": [float(value) * 1.0e-7 for value in row[2:47]],
        },
        "helium": {
            "x": [1000.0 * float(value) for value in parameters["HeFlux"]["bins"][0]["centers"]],
            "y": [float(value) * 1.0e-7 for value in row[92:132]],
        },
    }
    components = {item["component_id"]: item for item in canonical["components"]}

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex="col")
    for column, component_id in enumerate(("proton", "helium")):
        component = components[component_id]
        source = raw[component_id]
        x = [float(value) for value in component["bin_lower"]]
        y = [float(value) for value in component["value"]]
        uncertainty = [float(value) for value in component["total_uncertainty"]]
        spectrum = axes[0][column]
        spectrum.errorbar(
            x, y, yerr=uncertainty, fmt="o", markersize=3, linewidth=0.7,
            color=COLORS["blue"], label="adapter canonical values",
        )
        spectrum.scatter(
            source["x"], source["y"], marker="x", s=25,
            color=COLORS["orange"], label="published HAPI values (unit-converted)",
        )
        spectrum.set_xscale("log")
        spectrum.set_yscale("log")
        spectrum.set_ylabel("differential intensity\n[1/(cm² sr s (MeV/n))]")
        spectrum.set_title(component_id.capitalize())
        spectrum.grid(which="both", alpha=0.22)
        spectrum.legend(fontsize=8)

        residual = axes[1][column]
        intensity_residual = [
            100.0 * (actual - expected) / expected
            for actual, expected in zip(y, source["y"])
        ]
        coordinate_shift = [
            100.0 * (actual - published) / actual
            for actual, published in zip(x, source["x"])
        ]
        residual.semilogx(x, intensity_residual, "o-", markersize=3,
                          color=COLORS["blue"], label="intensity conversion residual")
        residual.semilogx(x, coordinate_shift, "s-", markersize=3,
                          color=COLORS["orange"], label="rounded HAPI edge shift")
        residual.axhline(0.0, color="#222222", linewidth=0.8)
        residual.set_xlabel("kinetic energy per nucleon [MeV/n]")
        residual.set_ylabel("relative difference [%]")
        residual.grid(which="both", alpha=0.22)
        residual.legend(fontsize=8)

    fig.suptitle("AMS-02 BR2460 observation and adapter transformation", fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage2_ams02_observation_adapter",
        "Published HAPI measurements are overlaid with the canonical adapter output. "
        "This validates ingestion, units, grids, and uncertainty retention only; the near-Earth "
        "AMS spectrum is not a complete mission environment or a transport-physics validation.",
        entries,
    )


def ams_solar_cycle_figure(canonical: dict, entries: list[dict]) -> None:
    """Compare BR2581 adapter values with all seven primary CSV tables."""

    directory = FIXTURES / "gcr" / "ams02_2025"
    interpretation = json.loads((directory / "interpretation.v1.json").read_text(encoding="utf-8"))
    raw: dict[str, list[dict[str, str]]] = {}
    for definition in interpretation["files"]:
        with (directory / definition["path"]).open(encoding="utf-8", newline="") as stream:
            raw[definition["symbol"]] = [
                row for row in csv.DictReader(stream)
                if int(row["bartels_rotation_number"]) == canonical["bartels_rotation"]
            ]
    components = {item["component_id"]: item for item in canonical["components"]}
    fig, axes = plt.subplots(2, 4, figsize=(16, 9.5))
    for axis, definition in zip(axes.flat, interpretation["files"]):
        symbol = definition["symbol"]
        name = definition["published_name"]
        component = components[symbol]
        rows = raw[symbol]
        centers = [math.sqrt(float(row["rigidity_min_GV"]) * float(row["rigidity_max_GV"]))
                   for row in rows]
        flux_name = f"{name}_flux m^-2sr^-1s^-1GV^-1"
        raw_flux = [float(row[flux_name]) for row in rows]
        values = [float(value) for value in component["value"]]
        display_uncertainty = [
            math.sqrt(stat * stat + time_dep * time_dep + time_ind * time_ind)
            for stat, time_dep, time_ind in zip(
                component["statistical_uncertainty"],
                component["time_dependent_systematic_uncertainty"],
                component["time_independent_systematic_uncertainty"],
            )
        ]
        axis.errorbar(centers, values, yerr=display_uncertainty, fmt="o", markersize=3,
                      linewidth=0.7, color=COLORS["blue"], label="adapter ± display total")
        axis.scatter(centers, raw_flux, marker="x", s=23, color=COLORS["orange"],
                     label="primary CSV")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(symbol)
        axis.set_xlabel("rigidity [GV]")
        axis.set_ylabel("intensity [1/(m² sr s GV)]")
        axis.grid(which="both", alpha=0.22)
    axes.flat[0].legend(fontsize=7)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(
        0.04, 0.82,
        "BR2581: chronological holdout\n\n"
        "Error bars show per-bin quadrature\nfor display only. The adapter retains\n"
        "statistical, time-dependent, and\ntime-independent terms separately.\n\n"
        "No model curve is plotted.",
        va="top", fontsize=10,
    )
    fig.suptitle("AMS-02 2025 solar-cycle tables: published values vs adapter",
                 fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage2_ams02_solar_cycle_observation_adapter",
        "All seven AMS-02 elemental spectra in held-out Bartels rotation 2581 are compared with "
        "the canonical adapter output. Display error bars combine the three published per-bin "
        "components only for visualization; time-independent correlations across rotations are "
        "not assumed. This is ingestion evidence, not a mission environment or transport validation.",
        entries,
    )


def ace_cris_figure(canonical: dict, entries: list[dict]) -> None:
    """Compare selected ACE/CRIS components with the pinned HAPI record."""

    with (FIXTURES / "gcr" / "ace_cris" / "ac_h3_cris_20191201_all.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        row = next(csv.reader(stream))
    symbols = ["B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
               "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
    components = {item["component_id"]: item for item in canonical["components"]}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    for axis, symbol in zip(axes.flat, ("O", "Si", "Fe", "Co")):
        component = components[symbol]
        species_index = symbols.index(symbol)
        raw_flux = [float(value) for value in row[1 + species_index * 7:1 + (species_index + 1) * 7]]
        raw_count = [float(value) for value in
                     row[1 + 24 * 7 + species_index * 7:1 + 24 * 7 + (species_index + 1) * 7]]
        lower = [float(value) for value in component["bin_lower"]]
        upper = [float(value) for value in component["bin_upper"]]
        centers = [(low + high) / 2.0 for low, high in zip(lower, upper)]
        valid = [bool(value) for value in component["valid"]]
        valid_indices = [index for index, keep in enumerate(valid) if keep]
        x = [centers[index] for index in valid_indices]
        xerr = [[centers[index] - lower[index] for index in valid_indices],
                [upper[index] - centers[index] for index in valid_indices]]
        values = [float(component["value"][index]) for index in valid_indices]
        errors = [float(component["statistical_uncertainty"][index]) for index in valid_indices]
        raw_values = [raw_flux[index] for index in valid_indices]
        axis.errorbar(x, values, xerr=xerr, yerr=errors, fmt="o", markersize=4,
                      linewidth=0.8, color=COLORS["blue"], label="adapter ± count statistic")
        axis.scatter(x, raw_values, marker="x", s=30, color=COLORS["orange"], label="HAPI value")
        invalid = [index for index, keep in enumerate(valid) if not keep]
        if invalid:
            floor = min(values) / 3.0
            axis.scatter([centers[index] for index in invalid], [floor] * len(invalid),
                         marker="v", facecolors="none", edgecolors=COLORS["purple"],
                         label="zero-count bin (invalid)")
        require(all(raw_count[index] > 0.0 for index in valid_indices),
                f"ACE/CRIS {symbol} valid mask disagrees with counts")
        axis.set_yscale("log")
        axis.set_title(symbol)
        axis.set_xlabel("kinetic energy/nucleon bin midpoint [MeV/n]")
        axis.set_ylabel("intensity [1/(cm² sr s (MeV/n))]")
        axis.grid(which="both", alpha=0.22)
        axis.legend(fontsize=7)
    fig.suptitle("ACE/CRIS daily B–Ni source: selected published values vs adapter",
                 fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage2_ace_cris_observation_adapter",
        "Selected ACE/CRIS elements compare the pinned HAPI values with canonical CDF-label "
        "energy bounds and count-derived statistical error bars. Open triangles identify "
        "zero-count bins excluded from comparisons. Systematic errors are not plotted because "
        "the processing-version mapping is unresolved. No model or transport result is shown.",
        entries,
    )


def parse_utc(value: str) -> dt.datetime:
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00"))


def goes_figure(canonical: dict, entries: list[dict]) -> None:
    csv_path = FIXTURES / "sep" / "goes16_sgps_gle72_timeseries.csv"
    with csv_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(line for line in stream if not line.startswith("#")))
    definition = json.loads(
        (FIXTURES / "sep" / "goes16_gle72_event_definition.v1.json").read_text(encoding="utf-8")
    )

    fig, axes = plt.subplots(3, 1, figsize=(14, 13))
    time_axis = axes[0]
    for sensor, color in (("0", COLORS["blue"]), ("1", COLORS["orange"])):
        selected = [row for row in rows if row["sensor_unit"] == sensor and row["channel"] == "P7"]
        times = [parse_utc(row["sample_midpoint_utc"]) for row in selected]
        values = [
            float(row["proton_flux_per_cm2_sr_s_MeV"])
            if row["proton_flux_per_cm2_sr_s_MeV"] != "NA" else math.nan
            for row in selected
        ]
        usable = [
            row["l1b_records_in_average"] == "300"
            and row["valid_1s_samples"] == "300"
            and row["dead_time_flagged_samples"] == "0"
            and row["out_of_band_flagged_samples"] == "0"
            and row["dynamic_error_flagged_samples"] == "0"
            and row["proton_flux_per_cm2_sr_s_MeV"] != "NA"
            and row["flux_uncertainty_per_cm2_sr_s_MeV"] != "NA"
            for row in selected
        ]
        time_axis.plot(times, values, color=color, linewidth=1.0, label=f"SGPS sensor {sensor}")
        rejected_times = [time for time, good in zip(times, usable) if not good]
        rejected_values = [value for value, good in zip(values, usable) if not good]
        if rejected_times:
            time_axis.scatter(rejected_times, rejected_values, marker="x", s=28,
                              color=color, label=f"sensor {sensor} rejected")
    windows = {item["role"]: item for item in definition["windows"]}
    sensitivity = windows["sensitivity"]
    primary = windows["primary"]
    time_axis.axvspan(parse_utc(sensitivity["start_utc"]), parse_utc(sensitivity["end_utc"]),
                      color="#fee090", alpha=0.35, label=">10 MeV sensitivity window")
    time_axis.axvspan(parse_utc(primary["start_utc"]), parse_utc(primary["end_utc"]),
                      color="#91bfdb", alpha=0.28, label=">100 MeV primary window")
    time_axis.axvline(parse_utc(definition["markers"][0]["time_utc"]), color=COLORS["purple"],
                      linestyle=":", linewidth=1.3, label="NMDB onset marker")
    time_axis.set_yscale("log")
    time_axis.set_ylabel("P7 differential proton flux\n[1/(cm² sr s MeV)]")
    time_axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M", tz=dt.timezone.utc))
    time_axis.grid(which="both", alpha=0.22)
    time_axis.legend(ncol=3, fontsize=8, loc="upper right")

    fluence_axis = axes[1]
    coverage_axis = axes[2]
    styles = {
        ("primary", "sgps_0"): (COLORS["blue"], "-", "o"),
        ("primary", "sgps_1"): (COLORS["orange"], "-", "o"),
        ("sensitivity", "sgps_0"): (COLORS["blue"], "--", "s"),
        ("sensitivity", "sgps_1"): (COLORS["orange"], "--", "s"),
    }
    for window in canonical["windows"]:
        role = window["role"]
        for sensor in ("sgps_0", "sgps_1"):
            selected = sorted(
                (item for item in window["channels"] if item["sensor_id"] == sensor),
                key=lambda item: item["energy_low_MeV"],
            )
            color, linestyle, marker = styles[(role, sensor)]
            energies = [math.sqrt(item["energy_low_MeV"] * item["energy_high_MeV"]) for item in selected]
            fluence = [item["partial_fluence_per_cm2_sr_MeV"] for item in selected]
            coverage = [item["coverage_fraction"] for item in selected]
            label = f"{role}, {sensor}"
            fluence_axis.plot(energies, fluence, linestyle=linestyle, marker=marker,
                              markersize=4, color=color, label=label)
            coverage_axis.plot(energies, coverage, linestyle=linestyle, marker=marker,
                               markersize=4, color=color, label=label)
    fluence_axis.set_xscale("log")
    fluence_axis.set_yscale("log")
    fluence_axis.set_ylabel("partial differential fluence\n[1/(cm² sr MeV)]")
    fluence_axis.grid(which="both", alpha=0.22)
    fluence_axis.legend(ncol=2, fontsize=8)
    coverage_axis.set_xscale("log")
    coverage_axis.set_ylim(0.0, 1.03)
    coverage_axis.set_xlabel("channel geometric-mean energy [MeV]")
    coverage_axis.set_ylabel("valid-time coverage fraction")
    coverage_axis.grid(which="both", alpha=0.22)
    coverage_axis.legend(ncol=2, fontsize=8)
    fig.suptitle("GOES-16 SGPS GLE72 observations and event integration", fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "stage2_goes16_observation_event",
        "P7 observations, approved integration windows, the independent NMDB marker, integrated "
        "channel fluences, and coverage are shown with both directions separate. Rejected samples "
        "are not interpolated. This is a geosynchronous observation/adapter validation, not an "
        "isotropic transport-boundary or radiation-proxy validation.",
        entries,
    )


def crater_figure(entries: list[dict]) -> None:
    path = FIXTURES / "observations" / "crater_gcr_let_spectra_v1_3.csv"
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    labels = rows[0][1:]
    let = [float(row[0]) for row in rows[1:]]
    values = [[float(item) for item in row[1:]] for row in rows[1:]]
    selected_indices = sorted({0, len(labels) // 3, 2 * len(labels) // 3, len(labels) - 1})

    fig, (curves, heatmap) = plt.subplots(2, 1, figsize=(13, 10))
    for index in selected_indices:
        spectrum = [row[index] if row[index] > 0.0 else math.nan for row in values]
        curves.plot(let, spectrum, linewidth=1.2, label=labels[index])
    curves.set_xscale("log")
    curves.set_yscale("log")
    curves.set_xlabel("LET [keV/µm]")
    curves.set_ylabel("published CRaTER spectral value\n[source-product units]")
    curves.grid(which="both", alpha=0.22)
    curves.legend(ncol=2, fontsize=8)
    curves.text(
        0.02, 0.04, "OBSERVATION ONLY — no shieldSim model curve at Stage 2",
        transform=curves.transAxes, fontsize=9, color="#8b0000",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )

    log_values = [
        [math.log10(value) if value > 0.0 else math.nan for value in row]
        for row in values
    ]
    let_edges = [let[0] * let[0] / math.sqrt(let[0] * let[1])]
    let_edges.extend(math.sqrt(lower * upper) for lower, upper in zip(let, let[1:]))
    let_edges.append(let[-1] * let[-1] / let_edges[-1])
    mesh = heatmap.pcolormesh(
        list(range(len(labels) + 1)), let_edges, log_values,
        shading="flat", cmap="viridis",
    )
    heatmap.set_yscale("log")
    heatmap.set_xlabel("published half-year interval index")
    heatmap.set_ylabel("LET [keV/µm]")
    colorbar = fig.colorbar(mesh, ax=heatmap, pad=0.02)
    colorbar.set_label("log10(published spectral value; source-product units)")
    fig.suptitle("CRaTER GCR LET observational reference", fontsize=14, fontweight="bold")
    save_figure(
        fig,
        "crater_observation_reference",
        "The pinned CRaTER v1.3 product is visualized for data-integrity review only. The CSV does "
        "not supply a generic slab-model comparison, and instrument geometry/response, selections, "
        "lunar environment, and preregistered validation intervals remain future gates.",
        entries,
    )


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    core = load_evidence("stage1_core.json")
    boundaries = load_evidence("stage1_case_boundaries.json")
    transport = load_evidence("stage1_transport_smoke.json")
    stage2 = load_evidence("stage2_source_adapters.json")
    published_transport = load_evidence("stage2_published_transport.json")
    stage3_core = load_evidence("stage3_core.json")
    stage3_weights = load_evidence("stage3_response_weights.json")
    stage3_transport = load_evidence("stage3_transport.json")
    stage4_core = load_evidence("stage4_core.json")
    stage4_nist = load_evidence("stage4_nist_reference.json")
    stage4_transport = load_evidence("stage4_tid_transport.json")
    stage5_core = load_evidence("stage5_core.json")
    stage5_transport = load_evidence("stage5_transport.json")
    stage6_core = load_evidence("stage6_core.json")
    stage6_reference = load_evidence("stage6_reference.json")
    stage6_transport = load_evidence("stage6_transport.json")
    stage7_core = load_evidence("stage7_core.json")
    stage7_boundaries = load_evidence("stage7_boundaries.json")
    stage7_published = load_evidence("stage7_published_qualification.json")
    entries: list[dict] = []

    comparison_figure(
        core,
        "Stage 1 core analytic contracts",
        "stage1_core_analytic",
        "Expected analytic conversions and uncertainty scaling are compared with the exact values "
        "accepted by stage1_core. Numerical assertions, not the image, define pass/fail.",
        entries,
    )
    comparison_figure(
        boundaries,
        "Stage 1 case-boundary analytic checks",
        "stage1_case_boundaries_analytic",
        "Single/split-layer equivalence, the inclusive thickness boundary, and density-dependent "
        "areal-to-geometric conversions are shown from stage1_case_boundaries evidence.",
        entries,
    )
    comparison_figure(
        transport,
        "Stage 1 transport geometry and ledger checks",
        "stage1_transport_smoke_analytic",
        "Resolved Geant4 material density, analytic slab mass, and deposited-energy ledger closure "
        "are shown for both independent absorbers. This does not validate transport physics.",
        entries,
    )
    comparison_figure(
        published_transport,
        "Downloaded AMS-02 spectrum: analytic source integral vs Geant4 input",
        "stage2_ams02_published_transport_analytic",
        "The independently reconstructed sum(J_i Delta E_i) is compared with the source "
        "normalization consumed by the full case/Geant4/proxy run. This closes the executable "
        "source boundary; it does not validate transport physics or make BR2460 a complete "
        "mission GCR environment.",
        entries,
    )
    comparison_figure(
        stage3_core,
        "Stage 3 exact geometry and response algebra",
        "stage3_core_analytic",
        "Analytic slab volume/mass and hand-calculated response folding/covariance values are "
        "compared with the deterministic Stage-3 implementation. Numerical assertions, not the "
        "image, define pass/fail; no transport-physics validation is claimed.",
        entries,
    )
    comparison_figure(
        stage3_transport,
        "Stage 3 transport equivalence and convergence",
        "stage3_transport_comparisons",
        "Five-seed single/split and direct/matrix-fold comparisons are shown with the final "
        "response-grid, cut, step, and transverse-size refinements. The numerical test enforces "
        "the preregistered uncertainty/equivalence limits. These are raw transport-engine "
        "acceptance checks pending independent review, not validation of a science quantity.",
        entries,
    )
    comparison_figure(
        stage3_weights,
        "Stage 3 response-grid interpolation algebra",
        "stage3_response_weights_analytic",
        "The production matrix assembler, explicit-weight builder, and folder are compared with "
        "an analytic tensor response affine in log-energy, mu, and periodic phi. Exact-node, "
        "extrapolation, and incomplete-grid rejection are asserted separately.",
        entries,
    )
    comparison_figure(
        stage4_core,
        "Stage 4 exact TID units and covariance",
        "stage4_core_analytic",
        "Exact MeV/kg-to-Gy, Gy-to-rad, and hand-calculated structured-covariance results are "
        "compared with the Stage-4 implementation. The image is a review aid; deterministic "
        "numeric assertions define pass/fail.",
        entries,
    )
    stage4_nist_figure(stage4_nist, entries)
    stage4_transport_figure(stage4_transport, entries)
    comparison_figure(
        stage5_core,
        "Stage 5 silicon damage-response algebra",
        "stage5_core_analytic",
        "Authenticated SR-NIEL exact nodes, log-log interpolation and analytic integration, the "
        "ASTM E722 neutron-group fold, and the 95 MeV-mb silicon neutron-equivalent convention "
        "are compared with independent closed-form values. Strict unsupported-species and "
        "out-of-domain assertions are numerical tests and are not encoded as plotted values.",
        entries,
    )
    comparison_figure(
        stage5_transport,
        "Stage 5 full-transport absorber-entry damage folds",
        "stage5_transport_comparisons",
        "Full Geant4 runs for exact-node proton, helium-4, and neutron cases are compared with "
        "independent table-fold values at the absorber entrance; neutron n_eq is checked against "
        "the ASTM E722 95 MeV-mb reference. The numerical suite separately enforces strict "
        "out-of-domain rejection and secondary-neutron sensitivity. Formal accelerator/device "
        "validation and independent review remain pending.",
        entries,
    )
    comparison_figure(
        stage6_core,
        "Stage 6 exact LET-bin and hardness mathematics",
        "stage6_core_analytic",
        "Exact constant and power-law partial-bin integrals, the analytic H100/10 fold, and a "
        "correlated ratio-of-means are compared with the deterministic Stage-6 implementation. "
        "Fieller zero/unresolved-denominator and non-clamping edge assertions remain numerical.",
        entries,
    )
    comparison_figure(
        stage6_reference,
        "Stage 6 electronic stopping references",
        "stage6_electronic_stopping_reference",
        "Geant4 unrestricted electronic mass stopping is compared at exact nodes with pinned "
        "NIST PSTAR/ASTAR data for protons and He-4 in Al, Si, and water, and with pinned GSI "
        "ATIMA v1.41 responses for C-12 in silicon. The C-12 rows are a 15% implementation "
        "screen, not an uncertainty estimate or experimental validation.",
        entries,
    )
    comparison_figure(
        stage6_transport,
        "Stage 6 full-transport LET and hardness checks",
        "stage6_transport_comparisons",
        "Thin-slab track-mean LET is compared with deterministic midpoint stopping references "
        "for proton, He-4, and C-12 cases; separate step, cut, and geometry convergence changes "
        "are shown. The weighted two-node H100/10 run is compared with its exact 2/3 solution. "
        "Numerical tests also enforce histogram closure and undefined zero denominators.",
        entries,
    )
    comparison_figure(
        stage7_core,
        "Stage 7 time, covariance, and response-quadrature algebra",
        "stage7_core_analytic",
        "Log-energy response integration, piecewise time integration with an explicit gap, "
        "full interval covariance propagation, and sampled-ensemble propagation are compared "
        "with independent analytic values. These deterministic checks validate algebra and "
        "strict boundaries, not a complete GCR or SEP environment.",
        entries,
    )
    comparison_figure(
        stage7_boundaries,
        "Stage 7 executable-boundary analytic checks",
        "stage7_boundaries_analytic",
        "The production Stage-7 CLI is compared with exact gap-preserving time integration, "
        "analytic adaptive response quadrature, and an exact Neyman history total. Separate "
        "numeric assertions reject hidden gaps, extrapolation, false full-GCR labels, and GOES "
        "transport without a boundary mapping.",
        entries,
    )
    comparison_figure(
        stage7_published,
        "Stage 7 downloaded AMS proton-component qualification",
        "stage7_published_component_comparisons",
        "Five-seed direct transport, explicit energy/angle quadrature, and response-matrix folds "
        "at all 45 published energy representatives are compared for both Al/HDPE layer orders; "
        "the nested-grid monitors are also shown. Direct and response uncertainties use the "
        "larger of within- and between-replicate standard errors for the rare TID tail. This "
        "qualifies the downloaded AMS A=1 proton component only. It is not a full GCR "
        "environment, and GOES transport remains disabled pending angular mapping.",
        entries,
    )

    # Keep transient adapter products beside the configured build artifacts.
    # This honors the repository test-output boundary and avoids hidden /tmp
    # dependencies; TemporaryDirectory still removes the intermediates after
    # the review figures and manifest have been written.
    with tempfile.TemporaryDirectory(
        prefix="shieldsim-validation-figures-", dir=OUTPUT
    ) as temporary:
        temporary_path = pathlib.Path(temporary)
        ams = run_source_tool(
            "ams02",
            [
                f"--info={FIXTURES / 'gcr' / 'ams02_20131118_hapi_info.json'}",
                f"--data={FIXTURES / 'gcr' / 'ams02_20131118_hapi_data.csv'}",
                f"--interpretation={FIXTURES / 'gcr' / 'ams02_20131118_interpretation.v1.json'}",
            ],
            temporary_path / "ams.json",
        )
        ams_cycle = run_source_tool(
            "ams02-cycle",
            [
                f"--interpretation={FIXTURES / 'gcr' / 'ams02_2025' / 'interpretation.v1.json'}",
                "--bartels-rotation=2581",
            ],
            temporary_path / "ams-cycle.json",
        )
        cris = run_source_tool(
            "ace-cris",
            [
                f"--info={FIXTURES / 'gcr' / 'ace_cris' / 'ac_h3_cris_hapi_info.json'}",
                f"--data={FIXTURES / 'gcr' / 'ace_cris' / 'ac_h3_cris_20191201_all.csv'}",
                f"--interpretation={FIXTURES / 'gcr' / 'ace_cris' / 'interpretation.v1.json'}",
            ],
            temporary_path / "cris.json",
        )
        goes = run_source_tool(
            "goes16",
            [
                f"--csv={FIXTURES / 'sep' / 'goes16_sgps_gle72_timeseries.csv'}",
                f"--definition={FIXTURES / 'sep' / 'goes16_gle72_event_definition.v1.json'}",
            ],
            temporary_path / "goes.json",
        )
        ams_figure(ams, entries)
        ams_solar_cycle_figure(ams_cycle, entries)
        ace_cris_figure(cris, entries)
        goes_figure(goes, entries)

    ams_comparisons = [
        item for item in stage2["comparisons"]
        if str(item["id"]).startswith("ams_") or str(item["id"]).startswith("ams2025_")
    ]
    comparison_figure(
        stage2,
        "AMS-02 analytic transformations and published-row comparisons",
        "stage2_ams02_analytic_comparisons",
        "The production adapter's first transformed proton/helium bin edges and differential "
        "intensities are compared with independent relativistic and primary-table Jacobian "
        "calculations; held-out BR2581 helium/oxygen values are compared with independently "
        "read primary CSV rows. Numerical tolerances remain in stage2_source_adapters.",
        entries,
        selected=ams_comparisons,
    )
    cris_comparisons = [
        item for item in stage2["comparisons"] if str(item["id"]).startswith("ace_cris_")
    ]
    comparison_figure(
        stage2,
        "ACE/CRIS observation and count-statistic comparisons",
        "stage2_ace_cris_observation_comparisons",
        "The production adapter is compared with an independently read oxygen flux and the "
        "analytic flux/sqrt(count) statistical uncertainty. Systematics remain explicitly "
        "unassessed for this CDF processing version.",
        entries,
        selected=cris_comparisons,
    )
    goes_comparisons = [
        item for item in stage2["comparisons"] if str(item["id"]).startswith("goes_")
    ]
    comparison_figure(
        stage2,
        "GOES event integration: analytic and independent-oracle comparisons",
        "stage2_goes16_analytic_comparisons",
        "Production adapter values are compared with exact overlap/linear analytic solutions and "
        "the independently summed published P7 event oracle. Numerical tests remain authoritative.",
        entries,
        selected=goes_comparisons,
    )
    crater_figure(entries)

    expected_stems = {
        "stage1_core_analytic",
        "stage1_case_boundaries_analytic",
        "stage1_transport_smoke_analytic",
        "stage2_ams02_observation_adapter",
        "stage2_ams02_published_transport_analytic",
        "stage2_ams02_solar_cycle_observation_adapter",
        "stage2_ams02_analytic_comparisons",
        "stage2_ace_cris_observation_adapter",
        "stage2_ace_cris_observation_comparisons",
        "stage2_goes16_observation_event",
        "stage2_goes16_analytic_comparisons",
        "crater_observation_reference",
        "stage3_core_analytic",
        "stage3_transport_comparisons",
        "stage3_response_weights_analytic",
        "stage4_core_analytic",
        "stage4_nist_reference_comparisons",
        "stage4_tid_transport_comparisons",
        "stage5_core_analytic",
        "stage5_transport_comparisons",
        "stage6_core_analytic",
        "stage6_electronic_stopping_reference",
        "stage6_transport_comparisons",
        "stage7_core_analytic",
        "stage7_boundaries_analytic",
        "stage7_published_component_comparisons",
    }
    require(
        {pathlib.Path(item["file"]).stem for item in entries} == expected_stems,
        "validation figure set is incomplete",
    )
    manifest = {
        "schema_version": "shieldSim.validation_figure_manifest.v1",
        "authority": "review_aid_only_numeric_tests_are_authoritative",
        "physical_validation_claim": False,
        "figures": sorted(entries, key=lambda item: item["file"]),
        "source_evidence": sorted(path.name for path in EVIDENCE.glob("*.json")),
    }
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Generated {len(expected_stems)} validation figures in PNG and SVG: {OUTPUT}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL validation figure generation: {error}", file=sys.stderr)
        raise SystemExit(1)
