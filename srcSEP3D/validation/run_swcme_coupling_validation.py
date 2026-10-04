#!/usr/bin/env python3
"""Check a native SWCME/SEP3D shock history and plot frozen observations.

This consumer does not implement CME dynamics or manufacture native evidence.
Use --amps to launch the native source-free propagation example, retain its
input/binary/log identities and plot its real installed-provider history.
An independently reviewed --event-fit is required for an absolute UTC
comparison with the downloaded observations.
--demo produces explicitly synthetic plotting fixtures, never a scientific PASS.
Only Matplotlib (and its NumPy dependency) is needed beyond the standard library.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
import sys

AU_M = 149597870700.0
SOLAR_RADIUS_M = 695700000.0
CASE_ID = "CME3D02"
SCHEMA = "srcsep3d-swcme-coupling-evidence-v1"
REFERENCE_SCHEMA = "srcsep3d-swcme-reference-v1"
DEFAULT_BUNDLE = Path(__file__).resolve().parent / "reference_data" / CASE_ID
MODEL_COLUMNS = ("time_s", "tick", "shock_radius_m", "shock_speed_m_s",
                 "shock_active", "generation", "particle_count",
                 "injected_particle_count", "mpi_radius_spread_m",
                 "mpi_clock_spread_s")
OBS_COLUMNS = ("time_utc", "radius_m", "sigma_radius_m", "quality", "role")
DEFAULT_ACCEPTANCE = dict(maximum_arrival_error_s=21600.0,
                          maximum_radius_normalized_rmse=2.0,
                          minimum_holdout_points=5,
                          minimum_holdout_coverage=0.9,
                          maximum_history_gap_s=300.0,
                          maximum_mpi_radius_spread_m=1.0,
                          maximum_mpi_clock_spread_s=1.0e-9)


class EvidenceError(ValueError):
    """Unusable/mismatched evidence, rather than a scientific disagreement."""


def require(condition, message):
    if not condition:
        raise EvidenceError(message)


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key: " + key)
            result[key] = value
        return result
    def nonfinite(value):
        raise EvidenceError("nonfinite JSON value: " + value)
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"),
                           object_pairs_hook=pairs, parse_constant=nonfinite)
    except (OSError, ValueError) as error:
        raise EvidenceError(str(error)) from error
    require(isinstance(value, dict), "JSON root must be an object")
    return value


def number(value, label, positive=False, integer=False):
    require(not isinstance(value, bool), label + " is a boolean")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise EvidenceError(label + " must be numeric") from error
    require(math.isfinite(value), label + " must be finite")
    require(not positive or value > 0, label + " must be positive")
    require(not integer or (value >= 0 and value.is_integer()),
            label + " must be a nonnegative integer")
    return int(value) if integer else value


def utc(value):
    require(isinstance(value, str), "UTC timestamp must be a string")
    try:
        stamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise EvidenceError("invalid UTC timestamp: " + value) from error
    require(stamp.tzinfo is not None and stamp.utcoffset() == timedelta(0),
            "timestamps must explicitly use UTC: " + value)
    return stamp.astimezone(timezone.utc)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_file(bundle, entry, label):
    require(isinstance(entry, dict), label + " file descriptor is absent")
    require(isinstance(entry.get("file"), str), label + " filename is absent")
    path = (bundle / entry["file"]).resolve()
    require(bundle in path.parents, label + " path escapes bundle")
    require(path.is_file(), label + " file is missing: " + str(path))
    digest = entry.get("sha256")
    require(isinstance(digest, str) and len(digest) == 64 and
            digest == sha256(path), label + " SHA256 mismatch")
    return path


def csv_rows(path, columns):
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None and
                    len(set(reader.fieldnames)) == len(reader.fieldnames),
                    "missing/duplicate CSV headers: " + str(path))
            require(set(columns).issubset(reader.fieldnames),
                    "missing CSV columns: " + str(path))
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise EvidenceError(str(error)) from error
    require(rows and all(None not in row and None not in row.values()
                         for row in rows), "empty/ragged CSV: " + str(path))
    return rows


def interpolate(rows, time_s, field):
    """Linear interpolation only inside actual native temporal coverage.

    Extrapolation could hide an early termination before 1 AU. Arrival is
    separately obtained by bracketing the front crossing, never by fitting an
    independent drag model here or shifting the modeled clock to observations.
    """
    if time_s < rows[0]["time_s"] or time_s > rows[-1]["time_s"]:
        return None
    import bisect
    i = bisect.bisect_left([row["time_s"] for row in rows], time_s)
    if i == 0 or rows[i]["time_s"] == time_s:
        return rows[i][field]
    left, right = rows[i - 1], rows[i]
    f = (time_s - left["time_s"]) / (right["time_s"] - left["time_s"])
    return left[field] + f * (right[field] - left[field])


def crossing(rows, target):
    for left, right in zip(rows, rows[1:]):
        if left["shock_radius_m"] <= target <= right["shock_radius_m"]:
            f = ((target - left["shock_radius_m"]) /
                 (right["shock_radius_m"] - left["shock_radius_m"]))
            return left["time_s"] + f * (right["time_s"] - left["time_s"])
    return None


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n"); stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)


def plot_comparison(output, manifest, model, observations, arrival_s, metrics):
    """Vector EPS and 600-dpi PNG share the same figure and unshifted UTC axis.

    Opaque symbols/error bars work in PostScript without transparency warnings.
    Calibration and withheld observations have distinct symbols. There is no
    plasma-speed overlay: bulk wind speed is not the front propagation speed.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as dates
    except ImportError as error:
        raise EvidenceError("plotting needs Matplotlib and NumPy") from error
    epoch = utc(manifest["launch_epoch_utc"])
    times = [epoch + timedelta(seconds=row["time_s"]) for row in model]
    synthetic = manifest["reference_kind"] == "synthetic"
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.linewidth": 0.8, "ps.fonttype": 42}):
        arrival_only = manifest.get("comparison_scope") == "arrival-only"
        # ``layout=`` is unavailable in the older Matplotlib shipped on some
        # supported AMPS systems.  The long-standing keyword produces the same
        # constrained layout without changing model data or scored metrics.
        fig, axes = plt.subplots(
            2 if arrival_only else 3, 1,
            figsize=(7.2, 5.3 if arrival_only else 7.5), sharex=True,
            constrained_layout=True)
        axes[0].plot(times, [r["shock_radius_m"] / AU_M for r in model],
                     color="#0072B2", label="Synthetic front" if synthetic else "SWCME front in srcSEP3D")
        axes[1].plot(times, [r["shock_speed_m_s"] / 1000 for r in model],
                     color="#0072B2", label="Modeled shock speed")
        for role, marker, label in (("calibration", "s", "Launch fit"),
                                     ("validation", "o", "Withheld shock track")):
            selected = [r for r in observations if r["role"] == role]
            if not selected:
                continue
            axes[0].errorbar([r["stamp"] for r in selected],
                            [r["radius_m"] / AU_M for r in selected],
                            yerr=[r["sigma_radius_m"] / AU_M for r in selected],
                            fmt=marker, ms=3.5, color="black", capsize=2,
                            mfc="white" if role == "calibration" else "black",
                            label=label)
        selected = [r for r in observations if r["role"] == "validation" and
                    interpolate(model, r["time_s"], "shock_radius_m") is not None]
        if not arrival_only:
            axes[2].axhline(0, color="0.5", lw=0.8)
            axes[2].errorbar([r["stamp"] for r in selected],
                [(interpolate(model, r["time_s"], "shock_radius_m") - r["radius_m"])
                 / SOLAR_RADIUS_M for r in selected],
                yerr=[r["sigma_radius_m"] / SOLAR_RADIUS_M for r in selected],
                fmt="o", ms=3.5, color="black", capsize=2)
        observed = utc(manifest["arrival"]["time_utc"])
        if arrival_only:
            axes[0].errorbar([observed], [manifest["arrival"]["heliocentric_radius_m"] / AU_M],
                xerr=timedelta(seconds=manifest["arrival"]["uncertainty_s"]),
                fmt="o", color="#D55E00", ms=5, label="Wind shock position/time")
        for ax in axes:
            ax.axvline(observed, color="#D55E00", ls="--", lw=1,
                       label="Observed spacecraft shock")
            if arrival_s is not None:
                ax.axvline(epoch + timedelta(seconds=arrival_s),
                           color="#0072B2", ls=":", lw=1,
                           label="Modeled spacecraft crossing")
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(axis="y", color="0.9", lw=0.5)
        axes[0].axhline(1, color="0.5", ls="-.", lw=0.7)
        axes[0].set_ylabel("Shock radius (AU)")
        axes[1].set_ylabel("Shock speed (km s$^{-1}$)")
        if not arrival_only:
            axes[2].set_ylabel("Radius residual ($R_\\odot$)")
        axes[-1].set_xlabel("UTC (unshifted)")
        locator = dates.AutoDateLocator(minticks=3, maxticks=8)
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(dates.ConciseDateFormatter(locator, tz=timezone.utc))
        axes[0].legend(loc="upper left", fontsize=7, framealpha=1)
        title = str(manifest["event_name"])
        if synthetic:
            title = "SYNTHETIC SOFTWARE CHECK - NOT OBSERVATIONS"
        axes[0].set_title(title)
        for label, ax in zip(("(a)", "(b)", "(c)"), axes):
            ax.text(0.98, 0.92, label, transform=ax.transAxes, ha="right")
        names = []
        for suffix in ("png", "eps"):
            path = output / ("swcme-coupling-comparison." + suffix)
            fig.savefig(path, dpi=600, facecolor="white")
            names.append(str(path))
        plt.close(fig)
    caption = ("Shock-front propagation from 20 solar radii. Squares denote "
               "launch-fit points; circles denote withheld shock-track points. "
               "Error bars show the supplied one-sigma radial uncertainties. "
               "Dashed and dotted lines identify the observed and modeled shock "
               "crossings at the spacecraft radius, which is distinct from 1 AU. "
               "No time shift or downstream refit is applied. Radius residuals "
               "are model minus observation. The speed curve is shock speed, "
               "not solar-wind bulk speed. Plasma/IMF evolution is outside this "
               "kinematic comparison. " +
               ("All curves and points are SYNTHETIC; not publication evidence."
                if synthetic else "See result.json for provenance and metrics."))
    if arrival_only:
        caption = ("Modeled shock propagation from 20 solar radii compared with the "
                   "Wind shock position/time. The observed crossing uses the actual "
                   "heliocentric spacecraft radius and the declared timing uncertainty. "
                   "The continuous curves are model outputs; no continuous "
                   "measured radial shock track is implied. The lower panel shows "
                   "modeled shock speed, not measured plasma bulk speed. Dashed and "
                   "dotted markers indicate observed and modeled spacecraft crossings "
                   "on an unshifted UTC axis. This is an arrival diagnostic, not "
                   "qualification of the full radial evolution or CME plasma/IMF. " +
                   ("All curves and points are SYNTHETIC; not publication evidence."
                    if synthetic else "See result.json for native provenance and metrics."))
    (output / "figure-caption.txt").write_text(caption + "\n", encoding="utf-8")
    return names + [str(output / "figure-caption.txt")]


def reference_bundle(bundle, output, manifest):
    """Verify downloaded references even while the native campaign is unavailable.

    A reference package must not invent executable identities, native histories,
    fitted launch clocks or shock radii to satisfy the full evidence schema.
    A later native-manifest.json supplies those owned native fields. Only arrival
    is compared until an independently identified shock-radius track exists.
    """
    require(manifest.get("case_id") == CASE_ID and manifest.get("reference_kind") == "observations",
            "reference schema/case mismatch")
    require(manifest.get("comparison_scope") == "arrival-only", "unsupported reference comparison scope")
    files = manifest.get("files")
    require(isinstance(files, list) and files, "reference file checksums are absent")
    verified = {row["file"]: checked_file(bundle, row, "reference") for row in files}
    require(len(verified) == len(files), "duplicate reference filenames")
    required = {"arrival.json", "wind_protons.csv", "wind_magnetic_1min.csv",
                "wind_magnetic_3sec_shock.csv", "wind_orbit.csv", "helcats_time_elongation.csv"}
    require(required.issubset(verified), "reference bundle lacks required datasets")
    arrival = manifest.get("arrival")
    require(isinstance(arrival, dict) and read_json(verified["arrival.json"]) == arrival,
            "arrival differs from checked arrival.json")
    require(arrival.get("spacecraft") == "Wind" and arrival.get("feature") == "shock", "not a Wind shock reference")
    utc(arrival.get("time_utc"))
    require(0.9 * AU_M < number(arrival.get("heliocentric_radius_m"), "spacecraft radius") < 1.1 * AU_M,
            "invalid reference spacecraft radius")
    number(arrival.get("uncertainty_s"), "arrival uncertainty", positive=True)
    contracts = {
        "wind_protons.csv": ("time_utc", "proton_bulk_speed_m_s", "proton_density_m3", "proton_temperature_k", "quality"),
        "wind_magnetic_1min.csv": ("time_utc", "magnitude_t", "bx_gse_t", "by_gse_t", "bz_gse_t", "quality"),
        "wind_magnetic_3sec_shock.csv": ("time_utc", "magnitude_t", "bx_gse_t", "by_gse_t", "bz_gse_t", "quality"),
        "wind_orbit.csv": ("time_utc", "x_hec_m", "y_hec_m", "z_hec_m", "heliocentric_radius_m"),
        "helcats_time_elongation.csv": ("time_utc", "event_id", "trace_id", "elongation_deg", "position_angle_deg", "feature", "role"),
    }
    actual_counts = {}
    for name, columns in contracts.items():
        rows = csv_rows(verified[name], columns)
        actual_counts[name] = len(rows)
        previous = None
        for row in rows:
            epoch = utc(row["time_utc"])
            if name != "helcats_time_elongation.csv":
                require(previous is None or epoch > previous, "nonmonotonic reference times: " + name)
            previous = epoch
            if "quality" in row:
                require(row["quality"] in ("good", "bad"), "invalid reference quality")
                if row["quality"] == "bad":
                    continue
            for key in columns:
                if key.endswith(("_m_s", "_m3", "_k", "_t", "_m", "_deg")):
                    number(row[key], "reference " + key)
    relative = manifest.get("native_manifest")
    require(isinstance(relative, str) and relative, "native manifest location is absent")
    native_path = (bundle / relative).resolve()
    require(bundle in native_path.parents, "native manifest path escapes bundle")
    if not native_path.is_file():
        return dict(id=CASE_ID, status="SKIP",
            message="reference data verified; installed native shock history and native-manifest.json are still required",
            metrics=dict(reference_counts=actual_counts,
                         observed_arrival_utc=arrival["time_utc"], spacecraft_radius_m=arrival["heliocentric_radius_m"]),
            artifacts=[str(path) for name, path in verified.items() if name.endswith((".png", ".eps"))],
            reference_ready=True, comparison_scope="arrival-only", production_release_qualified=False,
            full_plasma_coupling_validated=False, shock_radius_track_validated=False,
            reference_manifest_sha256=sha256(bundle / "manifest.json"))
    native = read_json(native_path)
    # Native inputs are independent of observed arrival. Reuse only the checked
    # spacecraft marker and never infer a launch time/speed from that marker.
    effective = dict(native, arrival=arrival, reference_kind="observations", comparison_scope="arrival-only",
                     event_name=manifest["event_name"])
    result = evaluate(bundle, output, _manifest=effective)
    result.update(reference_ready=True, reference_manifest_sha256=sha256(bundle / "manifest.json"),
                  native_manifest_sha256=sha256(native_path), shock_radius_track_validated=False)
    return result


def evaluate(bundle, output, _manifest=None):
    """Consume an immutable bundle; a missing bundle is explicitly SKIP.

    The producer declaration is provenance to audit, not proof that an arbitrary
    CSV came from AMPS. Retain its input, executable identity and execution log.
    A passing comparison does not qualify a production runtime by itself.
    """
    bundle, output = Path(bundle).resolve(), Path(output).resolve()
    result = dict(id=CASE_ID, status="SKIP", message="provide a reviewed native/observation bundle",
                  metrics={}, artifacts=[], production_release_qualified=False,
                  full_plasma_coupling_validated=False)
    if not (bundle / "manifest.json").is_file():
        return result
    manifest = read_json(bundle / "manifest.json") if _manifest is None else _manifest
    if manifest.get("schema") == REFERENCE_SCHEMA:
        return reference_bundle(bundle, output, manifest)
    require(manifest.get("schema") == SCHEMA and manifest.get("case_id") == CASE_ID,
            "unsupported evidence schema/case")
    require(manifest.get("reference_kind") in ("observations", "synthetic"),
            "reference_kind must be observations or synthetic")
    require(manifest.get("feature") == "shock-front", "ejecta tracks cannot validate shock radius")
    require(manifest.get("geometry") == "sun-centered-sphere", "unsupported front geometry")
    require(isinstance(manifest.get("event_name"), str) and manifest["event_name"], "event_name is absent")
    require(manifest.get("fit_used_holdout") is False, "downstream validation data were used to fit the model")
    epoch = utc(manifest.get("launch_epoch_utc"))
    fit_end = utc(manifest.get("calibration_end_utc"))
    require(fit_end >= epoch, "calibration cutoff precedes launch")
    synthetic = manifest["reference_kind"] == "synthetic"
    scope = manifest.get("comparison_scope", "shock-track-and-arrival")
    require(scope in ("arrival-only", "shock-track-and-arrival"), "unsupported comparison_scope")
    arrival_only = scope == "arrival-only"
    expected_producer = "synthetic-fixture" if synthetic else "srcSEP3D-native"
    require(manifest.get("producer") == expected_producer, "history must come from the declared native producer")
    require(isinstance(manifest.get("command"), list) and manifest["command"] and
            all(isinstance(x, str) for x in manifest["command"]), "command argv is absent")
    require(isinstance(manifest.get("executable_sha256"), str) and
            len(manifest["executable_sha256"]) == 64, "executable SHA256 is absent")
    number(manifest.get("mpi_ranks"), "mpi_ranks", positive=True, integer=True)
    input_path = checked_file(bundle, manifest.get("input"), "input")
    log_path = checked_file(bundle, manifest.get("execution_log"), "execution log")
    model_path = checked_file(bundle, manifest.get("model"), "native history")
    obs_path = None if arrival_only else checked_file(bundle, manifest.get("observations"), "observations")
    if not arrival_only:
        require(isinstance(manifest["observations"].get("provenance"), str) and
                manifest["observations"]["provenance"], "observation provenance is absent")
    require(isinstance(manifest.get("parameter_provenance"), str) and
            manifest["parameter_provenance"], "launch/drag/ambient provenance is absent")
    dt = number(manifest.get("time_step_s"), "time_step_s", positive=True)
    acceptance = manifest.get("acceptance")
    require(isinstance(acceptance, dict) and set(acceptance) == set(DEFAULT_ACCEPTANCE),
            "complete frozen acceptance thresholds are required")
    a = {key: number(value, key, positive=True) for key, value in acceptance.items()}
    require(a["minimum_holdout_coverage"] <= 1 and
            a["minimum_holdout_points"].is_integer(), "invalid coverage/count threshold")
    rows = csv_rows(model_path, MODEL_COLUMNS)
    model = [{key: number(row[key], key, integer=key in
              ("tick", "generation", "shock_active", "particle_count", "injected_particle_count"))
              for key in MODEL_COLUMNS} for row in rows]
    require(len(model) >= 2, "native history must contain a propagation interval")
    for previous, row in zip(model, model[1:]):
        require(row["time_s"] > previous["time_s"], "native time is not strictly increasing")
        require(row["shock_radius_m"] > previous["shock_radius_m"], "front radius is not strictly increasing")
    for row in model:
        require(row["shock_radius_m"] > 0 and row["shock_speed_m_s"] > 0 and
                row["mpi_radius_spread_m"] >= 0 and row["mpi_clock_spread_s"] >= 0,
                "invalid native geometry/spread")
        require(row["shock_active"] in (0, 1), "shock_active must be 0/1")
        require(abs(row["time_s"] - row["tick"] * dt) <= 1.e-8 * max(1, row["time_s"]),
                "native time does not agree with integer clock")
    observations = []
    rejected = 0
    for row in ([] if arrival_only else csv_rows(obs_path, OBS_COLUMNS)):
        require(row["quality"] in ("good", "bad"), "observation quality must be good/bad")
        require(row["role"] in ("calibration", "validation"), "invalid observation role")
        if row["quality"] == "bad":
            rejected += 1
            continue  # Declared instrument-invalid data never enter metrics.
        stamp = utc(row["time_utc"])
        require(stamp >= epoch, "accepted track precedes launch epoch")
        radius = number(row["radius_m"], "radius_m", positive=True)
        sigma = number(row["sigma_radius_m"], "sigma_radius_m", positive=True)
        require((row["role"] == "calibration" and stamp <= fit_end) or
                (row["role"] == "validation" and stamp > fit_end),
                "calibration and withheld observations overlap")
        require(radius >= 20 * SOLAR_RADIUS_M, "track precedes the 20-Rsun launch")
        observations.append(dict(stamp=stamp, time_s=(stamp-epoch).total_seconds(),
                                 radius_m=radius, sigma_radius_m=sigma, role=row["role"]))
    require(arrival_only or (observations and all(b["stamp"] > a_["stamp"] for a_, b in
                                zip(observations, observations[1:]))),
            "accepted observation times must be strictly increasing")
    arrival = manifest.get("arrival")
    require(isinstance(arrival, dict) and arrival.get("feature") == "shock" and
            arrival.get("role") == "validation", "withheld spacecraft shock marker is required")
    require(isinstance(arrival.get("spacecraft"), str) and arrival["spacecraft"] and
            isinstance(arrival.get("provenance"), str) and arrival["provenance"],
            "spacecraft/arrival provenance is absent")
    arrival_stamp = utc(arrival.get("time_utc"))
    require(arrival_stamp > fit_end, "arrival was part of the calibration interval")
    target = number(arrival.get("heliocentric_radius_m"), "spacecraft radius", positive=True)
    require(0.9 * AU_M < target < 1.1 * AU_M, "this case requires a near-Earth observer")
    number(arrival.get("uncertainty_s"), "arrival uncertainty", positive=True)
    observer_s = crossing(model, target)
    one_au_s = crossing(model, AU_M)
    holdout = [r for r in observations if r["role"] == "validation"]
    covered = [r for r in holdout if interpolate(model, r["time_s"], "shock_radius_m") is not None]
    residuals = [(interpolate(model, r["time_s"], "shock_radius_m") - r["radius_m"])
                 / r["sigma_radius_m"] for r in covered]
    metrics = dict(holdout_points=len(holdout), covered_holdout_points=len(covered),
                   holdout_coverage=len(covered) / len(holdout) if holdout else 0,
                   rejected_quality_points=rejected,
                   radius_normalized_rmse=math.sqrt(sum(r*r for r in residuals) / len(residuals)) if residuals else None,
                   modeled_one_au_time_s=one_au_s, modeled_spacecraft_time_s=observer_s,
                   arrival_error_s=(observer_s - (arrival_stamp-epoch).total_seconds()) if observer_s is not None else None,
                   maximum_history_gap_s=max(b["time_s"]-a_["time_s"] for a_, b in zip(model, model[1:])))
    violations = []
    def gate(condition, message):
        if not condition:
            violations.append(message)
    gate(manifest.get("source_enabled") is False, "particle source was enabled")
    gate(model[0]["tick"] == 0 and model[0]["time_s"] == 0, "history does not start at launch tick zero")
    gate(abs(model[0]["shock_radius_m"] - 20 * SOLAR_RADIUS_M) <= 1.0,
         "shock does not launch at 20 solar radii")
    gate(all(r["particle_count"] == 0 and r["injected_particle_count"] == 0 for r in model),
         "particle population/injection was nonzero")
    gate(all(r["shock_active"] == 1 for r in model), "shock inactive during recorded propagation")
    gate(all(r["generation"] > 0 for r in model), "active shock has no provider generation")
    gate(all(b["generation"] > a_["generation"] for a_, b in zip(model, model[1:])),
         "stale/non-increasing shock generation")
    gate(max(r["mpi_radius_spread_m"] for r in model) <= a["maximum_mpi_radius_spread_m"] and
         max(r["mpi_clock_spread_s"] for r in model) <= a["maximum_mpi_clock_spread_s"], "MPI ranks disagree")
    gate(one_au_s is not None, "history did not reach 1 AU")
    gate(metrics["maximum_history_gap_s"] <= a["maximum_history_gap_s"], "native history has excessive cadence/gaps")
    if not arrival_only:
        gate(len(holdout) >= a["minimum_holdout_points"], "too few withheld shock-track points")
        gate(metrics["holdout_coverage"] >= a["minimum_holdout_coverage"], "insufficient withheld temporal coverage")
        gate(metrics["radius_normalized_rmse"] is not None and
             metrics["radius_normalized_rmse"] <= a["maximum_radius_normalized_rmse"], "shock-track residual threshold missed")
    gate(metrics["arrival_error_s"] is not None and
         abs(metrics["arrival_error_s"]) <= a["maximum_arrival_error_s"], "shock-arrival threshold missed")
    output.mkdir(parents=True, exist_ok=True)
    artifacts = plot_comparison(output, manifest, model, observations, observer_s, metrics)
    result.update(status="FAIL" if violations else ("SKIP" if synthetic else "PASS"),
                  message="; ".join(violations) if violations else
                  ("synthetic plotting check; no observational qualification" if synthetic else
                   ("spacecraft arrival diagnostic threshold satisfied; radial track not validated" if arrival_only else
                    "withheld shock-track and arrival thresholds satisfied")),
                  metrics=metrics, artifacts=artifacts, violations=violations,
                  reference_kind=manifest["reference_kind"], acceptance=acceptance,
                  evidence_manifest_sha256=sha256(bundle / "manifest.json"),
                  input_sha256=sha256(input_path), execution_log_sha256=sha256(log_path),
                  model_sha256=sha256(model_path), observations_sha256=sha256(obs_path) if obs_path is not None else None,
                  comparison_scope=scope, shock_radius_track_validated=not arrival_only and not synthetic,
                  native_producer_required=True, native_exporter_implemented_here=False)
    atomic_json(output / "evidence-manifest.json", manifest)
    return result


def write_demo(bundle):
    """Constant-speed regression fixture, not a SWCME implementation or event fit."""
    bundle.mkdir(parents=True, exist_ok=True)
    epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
    speed, dt, initial = 700000.0, 60.0, 20 * SOLAR_RADIUS_M
    model = [{"time_s": t, "tick": int(t / dt), "shock_radius_m": initial + speed*t,
              "shock_speed_m_s": speed, "shock_active": 1, "generation": int(t / dt)+1,
              "particle_count": 0, "injected_particle_count": 0,
              "mpi_radius_spread_m": 0, "mpi_clock_spread_s": 0}
             for t in range(0, 240001, 60)]
    with (bundle / "history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=MODEL_COLUMNS)
        writer.writeheader(); writer.writerows(model)
    with (bundle / "observations.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=OBS_COLUMNS)
        writer.writeheader()
        for i, t in enumerate(range(0, 200001, 10000)):
            writer.writerow(dict(time_utc=(epoch+timedelta(seconds=t)).isoformat(),
                radius_m=initial+speed*t + SOLAR_RADIUS_M*0.2*math.sin(i),
                sigma_radius_m=SOLAR_RADIUS_M, quality="good",
                role="calibration" if i < 2 else "validation"))
    (bundle / "input.in").write_text("SYNTHETIC FIXTURE - NOT AN AMPS INPUT\n", encoding="utf-8")
    (bundle / "native.log").write_text("SYNTHETIC FIXTURE - NO NATIVE RUN\n", encoding="utf-8")
    manifest = dict(schema=SCHEMA, case_id=CASE_ID, event_name="Synthetic fixture",
        reference_kind="synthetic", producer="synthetic-fixture", feature="shock-front",
        geometry="sun-centered-sphere", fit_used_holdout=False, source_enabled=False,
        launch_epoch_utc=epoch.isoformat(), calibration_end_utc=(epoch+timedelta(seconds=10000)).isoformat(),
        time_step_s=dt, mpi_ranks=1, command=["synthetic-fixture"], executable_sha256="0"*64,
        parameter_provenance="synthetic constant-speed mechanics test", acceptance=DEFAULT_ACCEPTANCE,
        arrival=dict(spacecraft="synthetic observer", time_utc=(epoch+timedelta(seconds=(AU_M-initial)/speed)).isoformat(),
                     feature="shock", role="validation", heliocentric_radius_m=AU_M,
                     uncertainty_s=60, provenance="synthetic exact crossing"))
    for key, name in (("input", "input.in"), ("execution_log", "native.log"),
                      ("model", "history.csv"), ("observations", "observations.csv")):
        manifest[key] = dict(file=name, sha256=sha256(bundle / name))
    manifest["observations"]["provenance"] = "synthetic fixture, not observations"
    atomic_json(bundle / "manifest.json", manifest)


def read_native_history(path, runtime):
    """Check completed native telemetry; no dynamics are computed here."""
    require(runtime.get("schema") == "srcsep3d-native-shock-runtime-v1" and
            runtime.get("producer") == "srcSEP3D-native" and
            runtime.get("run_intent") == "shock-propagation" and
            runtime.get("source_enabled") is False, "wrong native producer/mode")
    dt = number(runtime.get("time_step_s"), "native time step", positive=True)
    model = [{key: number(row[key], key, integer=key in
             ("tick", "generation", "shock_active", "particle_count", "injected_particle_count"))
             for key in MODEL_COLUMNS} for row in csv_rows(path, MODEL_COLUMNS)]
    require(len(model) >= 2 and model[0]["tick"] == 0 and model[0]["time_s"] == 0,
            "native history lacks launch/propagation interval")
    require(len(model) == number(runtime.get("completed_steps"), "completed steps", integer=True) + 1,
            "runtime completion count differs from history")
    require(model[-1]["time_s"] == number(runtime.get("final_time_s"), "final time") and
            model[-1]["shock_radius_m"] == number(runtime.get("final_radius_m"), "final radius"),
            "runtime final state differs from history")
    for index, row in enumerate(model):
        require(row["tick"] == index and row["generation"] > 0 and row["shock_active"] in (0, 1) and
                abs(row["time_s"] - index * dt) <= 1e-8 * max(1, row["time_s"]),
                "native clock/generation is invalid or has missing ticks")
        require(row["shock_radius_m"] > 0 and row["shock_speed_m_s"] > 0 and
                row["mpi_radius_spread_m"] >= 0 and row["mpi_clock_spread_s"] >= 0,
                "invalid native geometry/spread")
        if index:
            require(row["shock_radius_m"] > model[index-1]["shock_radius_m"] and
                    row["generation"] > model[index-1]["generation"], "stale native front/generation")
    maximum=number(runtime.get("maximum_time_steps"), "maximum steps", positive=True, integer=True)
    require(runtime["completed_steps"] <= maximum, "native history exceeds its step budget")
    if runtime.get("stop_reason") == "step-budget":
        require(runtime["completed_steps"] == maximum, "native history stopped before its declared step budget")
    target = number(runtime.get("stop_shock_radius_m"), "radius stop")
    require(runtime.get("stop_reason") in ("shock-radius", "step-budget", "amps-termination"),
            "unknown native stop reason")
    if runtime["stop_reason"] == "shock-radius":
        require(target > 0 and model[-2]["shock_radius_m"] < target <= model[-1]["shock_radius_m"],
                "radius stop does not bracket the declared target")
    return model


def plot_native_control(output, model):
    """Relative clock only: an unfitted control has no observed UTC origin."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "ps.fonttype": 42}):
        fig, axes = plt.subplots(2, 1, figsize=(6.7, 5.1), sharex=True)
        hours = [row["time_s"] / 3600 for row in model]
        axes[0].plot(hours, [row["shock_radius_m"] / AU_M for row in model], color="black", lw=1.3)
        axes[0].axhline(1, color="0.5", ls="--", lw=.8, label="1 AU diagnostic")
        axes[0].set_ylabel("Front radius [AU]"); axes[0].legend(frameon=False)
        axes[1].plot(hours, [row["shock_speed_m_s"] / 1000 for row in model], color="black", lw=1.3)
        inactive = [row for row in model if row["shock_active"] == 0]
        if inactive:
            axes[1].plot([row["time_s"]/3600 for row in inactive],
                         [row["shock_speed_m_s"]/1000 for row in inactive],
                         "x", color="red", ms=3, label="Physical shock inactive")
            axes[1].legend(frameon=False)
        axes[1].set_ylabel("Front speed [km/s]"); axes[1].set_xlabel("Time since 20-Rsun launch [h]")
        axes[0].set_title("Native SWCME propagation control — event parameters unfitted", fontsize=10)
        for ax in axes: ax.grid(color="0.85", lw=.5); ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout()
        artifacts=[]
        for suffix in ("png", "eps"):
            path=output/("native-propagation-control."+suffix)
            fig.savefig(path, dpi=600); artifacts.append(str(path))
        plt.close(fig)
        return artifacts


def native_campaign(args, output):
    """Own the launcher/input/log boundary, then consume exported provider facts.

    Output filename rewrites are the only changes to the submitted input. Freeze
    both spellings and bind a reviewed event fit to the original input digest;
    no launch epoch is inferred from the withheld spacecraft arrival.
    """
    executable=args.amps.expanduser().resolve()
    input_path=args.input.expanduser().resolve()
    require(executable.is_file() and os.access(executable, os.X_OK), "AMPS executable is missing/not executable")
    require(input_path.is_file(), "native input is missing")
    ranks=number(args.ranks, "ranks", positive=True, integer=True)
    native=output/"native"; native.mkdir()
    original=output/"original-input.in"; shutil.copyfile(input_path, original)
    frozen=native/"input.in"
    text=original.read_text(encoding="utf-8")
    for key in ("initialization_mesh_tecplot_file", "initialization_parker_line_tecplot_file", "initialization_data_tecplot_file"):
        def replace(match):
            leaf=Path(match[2].strip()).name
            require(leaf not in ("", ".", ".."), "invalid initialization output filename")
            return match[1]+str(native/leaf)
        text,count=re.subn(r"(?m)^(\s*"+key+r"\s*=\s*)([^#\n]+)", replace, text)
        require(count == 1, "native input must define exactly one "+key)
    frozen.write_text(text, encoding="utf-8")
    before_input, before_executable=sha256(frozen),sha256(executable)
    base=[str(executable), "--input", str(frozen), "--output-dir", str(native)]
    print("[progress] native propagation parser/preflight", flush=True)
    preflight=subprocess.run(base+["--dry-run"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=180, check=False)
    (output/"preflight.log").write_bytes(preflight.stdout)
    require(preflight.returncode == 0, "native preflight failed; rebuild AMPS with propagation support; see preflight.log")
    resolved=dict(line.split("=",1) for line in preflight.stdout.decode("utf-8",errors="replace").splitlines() if "=" in line)
    require(resolved.get("run_intent") == "shock-propagation" and resolved.get("source_enabled") == "false",
            "native preflight does not confirm source-free propagation; rebuild the executable")
    try:
        launcher=shlex.split(args.launcher.format(ranks=ranks))
    except (KeyError, ValueError) as error:
        raise EvidenceError("invalid launcher template: "+str(error)) from error
    require(launcher and "{ranks}" in args.launcher, "launcher must contain {ranks}")
    command=launcher+base
    log=output/"native.log"
    atomic_json(output/"launch.json", dict(command=command, mpi_ranks=ranks,
                executable_sha256=before_executable, input_sha256=before_input,
                original_input_sha256=sha256(original), input_file=str(frozen), log_file=str(log)))
    print("[progress] launching "+shlex.join(command), flush=True)
    # Preserve raw bytes; lossy decoding is limited to console presentation, so
    # compiler/MPI diagnostics with non-UTF8 bytes cannot destroy the report.
    with log.open("wb") as stream:
        process=subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b""):
            stream.write(line); stream.flush()
            print(line.decode("utf-8",errors="replace"), end="", flush=True)
        code=process.wait()
    require(code == 0, "native MPI run exited "+str(code)+"; see "+str(log))
    require(sha256(frozen) == before_input and sha256(executable) == before_executable,
            "input/executable changed during the native run")
    runtime_path=native/"native-runtime.json"
    runtime=read_json(runtime_path)
    require(runtime.get("mpi_ranks") == ranks, "native MPI rank count differs from launcher request")
    require(runtime.get("application_configuration_fingerprint") == resolved.get("physics_fingerprint"),
            "native physics identity differs from parser preflight")
    require(runtime.get("provider_identity") == "canonical-swcme3d-standalone" and
            bool(runtime.get("provider_configuration_fingerprint")), "wrong/missing installed SWCME identity")
    require(runtime.get("history_file") == "shock-history.csv", "unexpected native history path")
    history=native/"shock-history.csv"
    model=read_native_history(history,runtime)
    violations=[]
    if abs(model[0]["shock_radius_m"]-20*SOLAR_RADIUS_M)>1: violations.append("front did not start at 20 solar radii")
    if any(row["particle_count"] or row["injected_particle_count"] for row in model): violations.append("particles/injection were nonzero")
    if any(row["shock_active"] != 1 for row in model): violations.append("physical shock became inactive")
    if any(row["mpi_radius_spread_m"]>1 or row["mpi_clock_spread_s"]>1e-9 for row in model): violations.append("MPI ranks disagree")
    arrival_s=crossing(model,AU_M)
    if arrival_s is None: violations.append("front did not reach 1 AU")
    artifacts=plot_native_control(output,model)+[str(history),str(runtime_path),str(log),str(output/"launch.json")]
    metrics=dict(native_rows=len(model), completed_steps=runtime["completed_steps"], modeled_one_au_time_s=arrival_s,
                 final_radius_au=model[-1]["shock_radius_m"]/AU_M, stop_reason=runtime["stop_reason"])
    result=dict(id=CASE_ID,status="FAIL" if violations else "SKIP", message="; ".join(violations) if violations else
                "native propagation checks passed; independently reviewed launch UTC/parameter fit required for observation comparison",
                native_history_ready=True, native_mechanics_pass=not violations, metrics=metrics, artifacts=artifacts,
                production_release_qualified=False, full_plasma_coupling_validated=False, shock_radius_track_validated=False)
    if violations: return result
    if args.event_fit is None:
        reference=evaluate(args.bundle,output/"reference-check")
        result["reference_ready"]=reference.get("reference_ready",False)
        result["artifacts"]+=reference.get("artifacts",[])
        return result
    fit=read_json(args.event_fit)
    require(fit.get("schema") == "srcsep3d-swcme-event-fit-v1" and fit.get("fit_used_holdout") is False,
            "event fit is unreviewed or used downstream holdout")
    require(fit.get("input_sha256") == sha256(original), "event fit is not bound to this original input")
    utc(fit.get("launch_epoch_utc")); utc(fit.get("calibration_end_utc"))
    require(isinstance(fit.get("parameter_provenance"),str) and fit["parameter_provenance"], "event fit lacks parameter provenance")
    evidence=output/"evidence"
    shutil.copytree(args.bundle.resolve(),evidence)
    # First verify all frozen downloaded products before adding our owned facts.
    reference_bundle(evidence,output/"reference-check",read_json(evidence/"manifest.json"))
    manifest=dict(schema=SCHEMA,case_id=CASE_ID,producer="srcSEP3D-native",reference_kind="observations",
                  feature="shock-front",geometry="sun-centered-sphere",fit_used_holdout=False,source_enabled=False,
                  launch_epoch_utc=fit["launch_epoch_utc"],calibration_end_utc=fit["calibration_end_utc"],
                  parameter_provenance=fit["parameter_provenance"],time_step_s=runtime["time_step_s"],
                  mpi_ranks=ranks,command=command,executable_sha256=before_executable,acceptance=DEFAULT_ACCEPTANCE)
    for key,path in (("input",frozen),("model",history),("execution_log",log)):
        destination=evidence/("native-"+path.name);shutil.copyfile(path,destination)
        manifest[key]=dict(file=destination.name,sha256=sha256(destination))
    shutil.copyfile(args.event_fit,evidence/"event-fit.json")
    manifest["event_fit"]=dict(file="event-fit.json",sha256=sha256(evidence/"event-fit.json"))
    atomic_json(evidence/"native-manifest.json",manifest)
    compared=evaluate(evidence,output)
    compared.update(native_history_ready=True,native_mechanics_pass=True,native_metrics=metrics)
    compared["artifacts"]+=artifacts
    return compared

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE,
                      help="reviewed directory containing manifest.json (default: installed validation/reference_data/CME3D02)")
    mode.add_argument("--demo", action="store_true", help="explicit synthetic figure preview only")
    parser.add_argument("--output-dir", type=Path, default=Path("test_output/swcme-coupling"))
    parser.add_argument("--amps", type=Path, help="launch rebuilt native AMPS and export source-free front history")
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parents[1]/"examples/sep3d_swcme_20rs_1au.in")
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--launcher", default="mpiexec -n {ranks}", help="launcher argv template; no shell evaluation")
    parser.add_argument("--event-fit", type=Path, help="independently reviewed event-fit JSON, bound to original input SHA256")
    parser.add_argument("--require-evidence", action="store_true", help="missing/synthetic evidence returns exit 2")
    args = parser.parse_args(argv)
    require(not (args.demo and args.amps), "--demo cannot launch native AMPS")
    require(args.event_fit is None or args.amps is not None, "--event-fit requires --amps")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Fresh output directory prevents old publication figures being mistaken
    # for a failed invocation's new artifacts. Never silently overwrite a run.
    require(not any(output.iterdir()), "output directory must be empty; choose a fresh run directory")
    print("[progress] checking native history and observation provenance", flush=True)
    if args.demo:
        write_demo(output / "synthetic-bundle")
    try:
        result = native_campaign(args, output) if args.amps else evaluate(output / "synthetic-bundle" if args.demo else args.bundle, output)
    except (EvidenceError, OSError, TypeError, KeyError, OverflowError, subprocess.TimeoutExpired) as error:
        result = dict(id=CASE_ID, status="ERROR", message=str(error), metrics={}, artifacts=[])
    atomic_json(output / "result.json", result)
    print("[" + CASE_ID + "] " + result["status"] + " " + result["message"], flush=True)
    for path in result["artifacts"]:
        print("artifact=" + path, flush=True)
    print("result_json=" + str(output / "result.json"), flush=True)
    return 2 if result["status"] == "ERROR" or (args.require_evidence and result["status"] == "SKIP") else (1 if result["status"] == "FAIL" else 0)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print("ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(2)
