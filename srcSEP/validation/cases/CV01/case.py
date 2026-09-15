"""End-to-end CV01 model, reference, scoring, and visualization pipeline."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


CASE_ID = "CV01"


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Dict[str, Any]) -> None:
    """Write auditable case metadata atomically in the evidence directory."""
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _load_and_validate(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        config = json.load(stream)
    if not isinstance(config, dict) or config.get("case_id") != CASE_ID:
        raise ValueError("CV01 input must be a JSON object with case_id=CV01")
    if config.get("schema") != "srcsep-validation-case-input-v1":
        raise ValueError("CV01 input uses an unsupported schema")

    model = config.get("model", {})
    operators = config.get("operators", {})
    line = config.get("field_line", {})
    packet = config.get("particle_packet", {})
    numerics = config.get("numerics", {})
    if model.get("mover") != "fte-dmumu":
        raise ValueError("CV01 currently requires the production fte-dmumu core")
    required_disabled = (
        "focusing", "adiabatic_momentum_change", "perpendicular_transport",
        "sources", "losses", "particle_wave_coupling",
    )
    if any(bool(operators.get(name)) for name in required_disabled):
        raise ValueError("CV01 requires every non-streaming operator to be disabled")
    if float(operators.get("pitch_angle_diffusion_per_s", math.nan)) != 0.0:
        raise ValueError("CV01 requires D_mumu=0 s^-1")
    for name in ("d_ln_abs_b_ds_per_m", "parallel_velocity_gradient_per_s",
                 "velocity_divergence_per_s", "field_aligned_strain_per_s"):
        if float(line.get(name, math.nan)) != 0.0:
            raise ValueError(f"CV01 uniform background requires {name}=0")
    if config.get("boundaries") != ["periodic", "open"]:
        raise ValueError("CV01 must run periodic and open boundary controls")
    mus = [float(value) for value in packet.get("pitch_angle_cosines", [])]
    if mus != [-1.0, -0.5, 0.5, 1.0]:
        raise ValueError("CV01 pitch-angle set is frozen at -1,-0.5,0.5,1")
    dts = [float(value) for value in numerics.get("time_steps_s", [])]
    if len(dts) != 3 or any(value <= 0.0 for value in dts):
        raise ValueError("CV01 requires exactly three positive timesteps")
    final_time = float(numerics["final_time_s"])
    sample_interval = float(numerics["sample_interval_s"])
    campaign_seed = numerics.get("campaign_seed")
    if (not isinstance(campaign_seed, int) or isinstance(campaign_seed, bool) or
            campaign_seed < 0 or campaign_seed > 2**64 - 1):
        raise ValueError("campaign_seed must be an unsigned 64-bit integer")
    for value in dts:
        if not math.isclose(final_time / value, round(final_time / value), abs_tol=1e-12):
            raise ValueError("final_time_s must be an integer multiple of every timestep")
        if not math.isclose(sample_interval / value,
                            round(sample_interval / value), abs_tol=1e-12):
            raise ValueError("sample_interval_s must be an integer multiple of every timestep")
    return config


def _relativistic_speed(config: Dict[str, Any]) -> float:
    constants = config["constants"]
    packet = config["particle_packet"]
    c = float(constants["speed_of_light_m_per_s"])
    mass = float(constants["proton_mass_kg"])
    kinetic_j = (float(packet["kinetic_energy_mev"]) * 1.0e6 *
                 float(constants["electron_volt_j"]))
    gamma = 1.0 + kinetic_j / (mass * c * c)
    return c * math.sqrt(1.0 - 1.0 / (gamma * gamma))


def _write_initial_particles(config: Dict[str, Any], path: Path) -> None:
    """Create one deterministic weighted Gaussian packet for each pitch angle.

    Positions, pitch angles, and statistical weights are written once and then
    consumed byte-for-byte by both the C++ model and Python reference.  Sharing
    this immutable input avoids hiding a model initialization discrepancy while
    keeping the two propagation algorithms independent.
    """
    packet = config["particle_packet"]
    count = int(packet["particles_per_pitch_angle"])
    if count < 5 or count % 2 == 0:
        raise ValueError("particles_per_pitch_angle must be an odd integer >=5")
    center = float(packet["center_m"])
    sigma = float(packet["sigma_m"])
    half_width = float(packet["gaussian_half_width_sigma"])
    target_weight = float(packet["total_weight_per_pitch_angle"])
    if sigma <= 0.0 or half_width <= 0.0 or target_weight <= 0.0:
        raise ValueError("packet width and weight must be positive")

    offsets = [(-half_width + 2.0 * half_width * index / (count - 1))
               for index in range(count)]
    raw_weights = [math.exp(-0.5 * offset * offset) for offset in offsets]
    normalization = target_weight / sum(raw_weights)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("particle_id", "initial_s_m", "mu", "weight"))
        particle_id = 1
        for mu in packet["pitch_angle_cosines"]:
            for offset, raw_weight in zip(offsets, raw_weights):
                writer.writerow((particle_id, center + sigma * offset,
                                 float(mu), raw_weight * normalization))
                particle_id += 1


def _run(command: Sequence[str], cwd: Path, log: Path,
         environment: Optional[Dict[str, str]] = None,
         timeout: Optional[float] = None) -> None:
    """Run one stage, append its exact command/output, and propagate failure."""
    # CV01 owns six linked realizations plus independent reference commands
    # instead of using linked_case_common.run_linked_model(). Echo every exact,
    # shell-escaped invocation here so --all offers the same audit trail for
    # this multi-realization case as it does for every other CV/IV/XM case.
    print("RUN:", shlex.join(command), flush=True)
    with log.open("a", encoding="utf-8") as stream:
        stream.write("COMMAND: " + " ".join(command) + "\n")
        try:
            completed = subprocess.run(command, cwd=str(cwd), text=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, check=False,
                                       env=environment, timeout=timeout)
        except subprocess.TimeoutExpired as error:
            # Persist partial child output before raising so timeout failures
            # retain the same forensic value as ordinary nonzero exits.
            partial = error.stdout or ""
            if isinstance(partial, bytes):
                partial = partial.decode("utf-8", errors="replace")
            stream.write(partial)
            stream.write("TIMEOUT\n")
            raise RuntimeError(
                f"command exceeded {timeout} seconds; see {log}") from error
        stream.write(completed.stdout)
        stream.write(f"EXIT: {completed.returncode}\n")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with status {completed.returncode}; see {log}")


def _write_native_arguments(path: Path, arguments: Sequence[str]) -> None:
    """Serialize one linked-model invocation in the narrow native protocol.

    The reviewed JSON remains the authoritative human-facing configuration.
    This generated file contains one argument token per line because the C++
    registry callback should not duplicate a permissive JSON parser. Newlines
    are rejected so every token has exactly one interpretation.
    """
    if not arguments or len(arguments) % 2:
        raise ValueError("native CV01 arguments must be non-empty name/value pairs")
    if any(not value or "\n" in value or "\r" in value for value in arguments):
        raise ValueError("native CV01 argument tokens must be non-empty single lines")
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("srcsep-cv01-native-args-v1\n")
        for value in arguments:
            stream.write(value + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _validate_native_report(path: Path, expected_model: Path) -> None:
    """Require the linked registry—not process completion alone—to pass CV01."""
    with path.open("r", encoding="utf-8") as stream:
        report = json.load(stream)
    results = report.get("results", []) if isinstance(report, dict) else []
    if len(results) != 1 or str(results[0].get("id", "")).upper() != CASE_ID:
        raise RuntimeError("linked application returned an unexpected CV01 registry report")
    if str(results[0].get("status", "ERROR")).upper() != "PASS":
        raise RuntimeError(
            "linked application CV01 model stage did not PASS: " +
            str(results[0].get("message", "no diagnostic")))
    if not expected_model.is_file() or expected_model.stat().st_size == 0:
        raise RuntimeError("linked application did not publish the CV01 model CSV")


def _read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _row_key(row: Dict[str, str]) -> Tuple[str, float, int, float, float]:
    # C++ emits an integral-valued double as ``60`` while Python's CSV writer
    # may emit ``60.0``.  Numeric keys keep that harmless spelling difference
    # from masquerading as missing physics evidence.
    return (row["boundary"], float(row["dt_s"]), int(row["particle_id"]),
            float(row["initial_mu"]), float(row["time_s"]))


def _float(row: Dict[str, str], name: str) -> float:
    return float(row[name])


def _group_moments(rows: Iterable[Dict[str, str]]) -> Dict[Tuple[str, float, float, float], Dict[str, float]]:
    """Calculate weighted packet moments from saved rows, never live arrays."""
    groups: Dict[Tuple[str, float, float, float], List[Dict[str, str]]] = {}
    for row in rows:
        key = (row["boundary"], float(row["dt_s"]),
               float(row["initial_mu"]), float(row["time_s"]))
        groups.setdefault(key, []).append(row)
    result: Dict[Tuple[str, float, float, float], Dict[str, float]] = {}
    for key, group in groups.items():
        active_weight = sum(_float(row, "active_weight") for row in group)
        escaped_weight = sum(_float(row, "escaped_weight") for row in group)
        coordinate_name = ("unwrapped_position_m" if key[0] == "periodic"
                           else "position_m")
        if active_weight > 0.0:
            centroid = sum(_float(row, coordinate_name) *
                           _float(row, "active_weight") for row in group) / active_weight
            variance = sum(_float(row, "active_weight") *
                           (_float(row, coordinate_name) - centroid) ** 2
                           for row in group) / active_weight
        else:
            centroid = math.nan
            variance = math.nan
        result[key] = {
            "centroid_m": centroid,
            "variance_m2": variance,
            "active_weight": active_weight,
            "escaped_weight": escaped_weight,
        }
    return result


def _comparison(config: Dict[str, Any], model_rows: List[Dict[str, str]],
                reference_rows: List[Dict[str, str]], output_dir: Path) -> Tuple[List[Dict[str, Any]], bool]:
    """Score particle trajectories, moments, boundaries, and conservation."""
    references = {_row_key(row): row for row in reference_rows}
    if len(references) != len(reference_rows):
        raise ValueError("reference contains duplicate particle/time keys")
    if len(model_rows) != len(reference_rows):
        raise ValueError("model/reference row counts differ")

    line = config["field_line"]
    length = float(line["maximum_m"]) - float(line["minimum_m"])
    final_time = float(config["numerics"]["final_time_s"])
    maximum_position_error = 0.0
    maximum_mu_error = 0.0
    maximum_momentum_relative_error = 0.0
    maximum_crossing_error = 0.0
    active_mismatches = 0
    for model in model_rows:
        reference = references.get(_row_key(model))
        if reference is None:
            raise ValueError(f"reference lacks model key {_row_key(model)}")
        active_mismatches += int(model["active"] != reference["active"])
        # After an open-boundary escape the particle no longer has a modeled
        # exterior trajectory.  Compare its exact endpoint and crossing time;
        # compare unwrapped coordinates only for the periodic control where the
        # characteristic remains represented for the full interval.
        position_errors = [
            abs(_float(model, "position_m") - _float(reference, "position_m"))]
        if model["boundary"] == "periodic":
            position_errors.append(
                abs(_float(model, "unwrapped_position_m") -
                    _float(reference, "unwrapped_position_m")))
        maximum_position_error = max(maximum_position_error, *position_errors)
        maximum_mu_error = max(
            maximum_mu_error, abs(_float(model, "mu") - _float(reference, "mu")))
        reference_momentum = _float(reference, "momentum_kg_m_per_s")
        maximum_momentum_relative_error = max(
            maximum_momentum_relative_error,
            abs(_float(model, "momentum_kg_m_per_s") - reference_momentum) /
            max(abs(reference_momentum), sys.float_info.min),
        )
        if model["crossing_time_s"] and reference["crossing_time_s"]:
            maximum_crossing_error = max(
                maximum_crossing_error,
                abs(_float(model, "crossing_time_s") -
                    _float(reference, "crossing_time_s")),
            )

    model_moments = _group_moments(model_rows)
    reference_moments = _group_moments(reference_rows)
    maximum_centroid_error = 0.0
    maximum_variance_error = 0.0
    maximum_weight_error = 0.0
    summary_path = output_dir / "CV01_packet_moments.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as stream:
        fields = ("boundary", "dt_s", "mu", "time_s", "model_centroid_m",
                  "reference_centroid_m", "centroid_error_m",
                  "model_variance_m2", "reference_variance_m2",
                  "variance_error_m2", "model_active_weight",
                  "reference_active_weight", "model_escaped_weight",
                  "reference_escaped_weight", "weight_closure_error")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for key in sorted(model_moments):
            model = model_moments[key]
            reference = reference_moments[key]
            centroid_error = (0.0 if math.isnan(model["centroid_m"]) and
                              math.isnan(reference["centroid_m"]) else
                              abs(model["centroid_m"] - reference["centroid_m"]))
            variance_error = (0.0 if math.isnan(model["variance_m2"]) and
                              math.isnan(reference["variance_m2"]) else
                              abs(model["variance_m2"] - reference["variance_m2"]))
            closure_error = abs(model["active_weight"] +
                                model["escaped_weight"] - 1.0)
            maximum_centroid_error = max(maximum_centroid_error, centroid_error)
            maximum_variance_error = max(maximum_variance_error, variance_error)
            maximum_weight_error = max(
                maximum_weight_error, closure_error,
                abs(model["active_weight"] - reference["active_weight"]),
                abs(model["escaped_weight"] - reference["escaped_weight"]),
            )
            writer.writerow({
                "boundary": key[0], "dt_s": key[1], "mu": key[2],
                "time_s": key[3], "model_centroid_m": model["centroid_m"],
                "reference_centroid_m": reference["centroid_m"],
                "centroid_error_m": centroid_error,
                "model_variance_m2": model["variance_m2"],
                "reference_variance_m2": reference["variance_m2"],
                "variance_error_m2": variance_error,
                "model_active_weight": model["active_weight"],
                "reference_active_weight": reference["active_weight"],
                "model_escaped_weight": model["escaped_weight"],
                "reference_escaped_weight": reference["escaped_weight"],
                "weight_closure_error": closure_error,
            })

    # The primary series follows one representative outward characteristic.
    # The generic Python runner recognizes the numerical/analytical columns and
    # creates an additional standard overlay plus residual panel from this file.
    finest = min(float(value) for value in config["numerics"]["time_steps_s"])
    primary_path = output_dir / "CV01_solution.csv"
    with primary_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time_s", "numerical", "analytical", "residual",
                         "quantity", "units"))
        for key in sorted(model_moments, key=lambda value: float(value[3])):
            if (key[0] == "periodic" and math.isclose(float(key[1]), finest) and
                    math.isclose(float(key[2]), 0.5)):
                model = model_moments[key]["centroid_m"]
                reference = reference_moments[key]["centroid_m"]
                writer.writerow((float(key[3]), model, reference,
                                 model - reference,
                                 "packet centroid along field line", "m"))

    # A sign-reversed characteristic is an explicit negative control.  Its
    # failure proves that the configured case can detect a flipped pitch-angle
    # or field-line orientation rather than passing any translating packet.
    speed = _relativistic_speed(config)
    plasma_speed = float(line["plasma_speed_m_per_s"])
    packet = config["particle_packet"]
    wrong_sign_centroid = (float(packet["center_m"]) +
                           (plasma_speed - speed) * final_time)
    # CSV formatting can choose either integral or decimal spellings. Locate by
    # numerical equality to keep the scorer independent of text formatting.
    correct_centroid = next(
        value["centroid_m"] for key, value in model_moments.items()
        if key[0] == "periodic" and math.isclose(float(key[1]), finest) and
        math.isclose(float(key[2]), 1.0) and
        math.isclose(float(key[3]), final_time))
    negative_control = abs(correct_centroid - wrong_sign_centroid) / length

    acceptance = config["acceptance"]
    values = {
        "trajectory_error_over_line_length": maximum_position_error / length,
        "centroid_error_over_line_length": maximum_centroid_error / length,
        "variance_error_over_line_length_squared": maximum_variance_error /
        (length * length),
        "weight_closure_relative_error": maximum_weight_error,
        "boundary_time_error_over_final_time": maximum_crossing_error / final_time,
        "pitch_angle_absolute_error": maximum_mu_error,
        "momentum_relative_error": maximum_momentum_relative_error,
        "negative_control_error_over_line_length": negative_control,
        "active_state_mismatches": float(active_mismatches),
    }
    limits = {
        "trajectory_error_over_line_length": ("<=", float(acceptance["trajectory_error_over_line_length_max"])),
        "centroid_error_over_line_length": ("<=", float(acceptance["centroid_error_over_line_length_max"])),
        "variance_error_over_line_length_squared": ("<=", float(acceptance["variance_error_over_line_length_squared_max"])),
        "weight_closure_relative_error": ("<=", float(acceptance["weight_closure_relative_error_max"])),
        "boundary_time_error_over_final_time": ("<=", float(acceptance["boundary_time_error_over_final_time_max"])),
        "pitch_angle_absolute_error": ("<=", float(acceptance["pitch_angle_absolute_error_max"])),
        "momentum_relative_error": ("<=", float(acceptance["momentum_relative_error_max"])),
        "negative_control_error_over_line_length": (">=", float(acceptance["negative_control_minimum_error_over_line_length"])),
        "active_state_mismatches": ("==", 0.0),
    }
    units = {
        "trajectory_error_over_line_length": "dimensionless",
        "centroid_error_over_line_length": "dimensionless",
        "variance_error_over_line_length_squared": "dimensionless",
        "weight_closure_relative_error": "dimensionless",
        "boundary_time_error_over_final_time": "dimensionless",
        "pitch_angle_absolute_error": "dimensionless",
        "momentum_relative_error": "dimensionless",
        "negative_control_error_over_line_length": "dimensionless",
        "active_state_mismatches": "count",
    }
    metrics: List[Dict[str, Any]] = []
    passed = True
    for name, value in values.items():
        operator, tolerance = limits[name]
        accepted = ((value <= tolerance) if operator == "<=" else
                    (value >= tolerance) if operator == ">=" else
                    math.isclose(value, tolerance, abs_tol=0.0))
        passed = passed and accepted
        metrics.append({"name": name, "value": value,
                        "tolerance": tolerance, "comparison": operator,
                        "units": units[name]})
    metrics.append({"name": "assertion_failures", "value": 0.0 if passed else 1.0,
                    "tolerance": 0.0, "comparison": "<=", "units": "count"})
    return metrics, passed


def _write_profile(model_rows: List[Dict[str, str]],
                   reference_rows: List[Dict[str, str]],
                   config: Dict[str, Any], output_dir: Path) -> Path:
    """Bin final saved particle states for the profile panel and reviewers."""
    line = config["field_line"]
    bins = int(config["plot"]["profile_bins"])
    minimum = float(line["minimum_m"])
    maximum = float(line["maximum_m"])
    width = (maximum - minimum) / bins
    finest = min(float(value) for value in config["numerics"]["time_steps_s"])
    final_time = float(config["numerics"]["final_time_s"])

    def histogram(rows: Iterable[Dict[str, str]], boundary: str,
                  mu: float) -> List[float]:
        values = [0.0] * bins
        for row in rows:
            if (row["boundary"] != boundary or
                    not math.isclose(_float(row, "dt_s"), finest) or
                    not math.isclose(_float(row, "initial_mu"), mu) or
                    not math.isclose(_float(row, "time_s"), final_time) or
                    row["active"] != "1"):
                continue
            index = min(bins - 1, max(0, int((_float(row, "position_m") - minimum) / width)))
            values[index] += _float(row, "active_weight") / width
        return values

    path = output_dir / "CV01_final_profiles.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("boundary", "mu", "bin_center_m",
                         "model_weight_density_per_m",
                         "reference_weight_density_per_m"))
        for boundary in config["boundaries"]:
            for mu in config["particle_packet"]["pitch_angle_cosines"]:
                model = histogram(model_rows, boundary, float(mu))
                reference = histogram(reference_rows, boundary, float(mu))
                for index in range(bins):
                    writer.writerow((boundary, float(mu), minimum + (index + 0.5) * width,
                                     model[index], reference[index]))
    return path


def _plot_four_panel(config: Dict[str, Any], moments_path: Path,
                     profile_path: Path, output_dir: Path) -> List[Path]:
    """Render the four diagnostics prescribed by the CV01 campaign plan."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    moments = _read_rows(moments_path)
    profiles = _read_rows(profile_path)
    dts = sorted(float(value) for value in config["numerics"]["time_steps_s"])
    finest = min(dts)
    line_length = (float(config["field_line"]["maximum_m"]) -
                   float(config["field_line"]["minimum_m"]))
    au = float(config["constants"]["astronomical_unit_m"])

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    profile_axis, centroid_axis, variance_axis, weight_axis = axes.flat
    colors = {0.5: "#1565c0", -0.5: "#c62828"}
    # Use two surviving packets in the profile panel so both curves remain
    # visible at final time.  The |mu|=1 packets deliberately escape through
    # the open boundaries and are represented in the weight-ledger panel.
    for boundary, mu in (("periodic", 0.5), ("open", -0.5)):
        selected = [row for row in profiles if row["boundary"] == boundary and
                    math.isclose(float(row["mu"]), mu)]
        x = [float(row["bin_center_m"]) / au for row in selected]
        model = [float(row["model_weight_density_per_m"]) * au for row in selected]
        reference = [float(row["reference_weight_density_per_m"]) * au for row in selected]
        label = f"{boundary}, mu={mu:g}"
        profile_axis.plot(x, reference, color=colors[mu], linewidth=2.0,
                          label=f"reference: {label}")
        profile_axis.plot(x, model, color=colors[mu], linestyle="--",
                          linewidth=1.4, label=f"model: {label}")
    profile_axis.set_title("Final packet profiles")
    profile_axis.set_xlabel("field-line position [AU]")
    profile_axis.set_ylabel("statistical-weight density [AU$^{-1}$]")
    profile_axis.legend(fontsize=7, frameon=False)

    for dt_s in dts:
        grouped: Dict[float, Tuple[float, float]] = {}
        for row in moments:
            if not math.isclose(float(row["dt_s"]), dt_s):
                continue
            time_s = float(row["time_s"])
            centroid = float(row["centroid_error_m"]) / line_length
            variance = float(row["variance_error_m2"]) / (line_length * line_length)
            old = grouped.get(time_s, (0.0, 0.0))
            grouped[time_s] = (max(old[0], centroid), max(old[1], variance))
        times = sorted(grouped)
        centroid_axis.semilogy(times, [max(grouped[t][0], 1e-18) for t in times],
                               label=f"dt={dt_s:g} s")
        variance_axis.semilogy(times, [max(grouped[t][1], 1e-30) for t in times],
                               label=f"dt={dt_s:g} s")
    centroid_axis.set_title("Maximum centroid error")
    centroid_axis.set_xlabel("time [s]")
    centroid_axis.set_ylabel("|error| / line length")
    centroid_axis.legend(fontsize=8, frameon=False)
    variance_axis.set_title("Maximum variance error")
    variance_axis.set_xlabel("time [s]")
    variance_axis.set_ylabel("|error| / line length$^2$")
    variance_axis.legend(fontsize=8, frameon=False)

    for mu in config["particle_packet"]["pitch_angle_cosines"]:
        selected = [row for row in moments if row["boundary"] == "open" and
                    math.isclose(float(row["dt_s"]), finest) and
                    math.isclose(float(row["mu"]), float(mu))]
        selected.sort(key=lambda row: float(row["time_s"]))
        weight_axis.plot([float(row["time_s"]) for row in selected],
                         [float(row["model_active_weight"]) for row in selected],
                         label=f"retained, mu={float(mu):g}")
        weight_axis.plot([float(row["time_s"]) for row in selected],
                         [float(row["model_escaped_weight"]) for row in selected],
                         linestyle="--",
                         label=f"escaped, mu={float(mu):g}")
    weight_axis.set_title("Open-boundary weight ledger")
    weight_axis.set_xlabel("time [s]")
    weight_axis.set_ylabel("fraction of initial weight")
    weight_axis.set_ylim(-0.03, 1.03)
    weight_axis.legend(fontsize=7, ncol=2, frameon=False)

    for axis in axes.flat:
        axis.grid(True, color="0.87", linewidth=0.7)
    figure.suptitle("CV01 — ballistic streaming on a uniform field line", fontsize=14)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    outputs: List[Path] = []
    for extension in config["plot"]["formats"]:
        path = output_dir / f"CV01_four_panel.{extension}"
        figure.savefig(path, format=extension, dpi=180, bbox_inches="tight",
                       facecolor="white")
        outputs.append(path)
    plt.close(figure)
    return outputs


def run_case(*, source_root: Path, input_path: Path, output_dir: Path,
             executable: Path, timeout: Optional[float]) -> Dict[str, Any]:
    """Run CV01 and return one standard component-test result object."""
    started = time.monotonic()
    config = _load_and_validate(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_input = output_dir / "resolved_input.json"
    _write_json(resolved_input, config)
    particles_path = output_dir / "initial_particles.csv"
    _write_initial_particles(config, particles_path)
    log_path = output_dir / "CV01_run.log"
    log_path.write_text("", encoding="utf-8")

    driver_source = source_root / "validation" / "cases" / CASE_ID / "cv01_model.cpp"
    reference_source = source_root / "validation" / "cases" / CASE_ID / "reference_solution.py"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"linked srcSEP/AMPS executable is unavailable: {executable}")

    constants = config["constants"]
    line = config["field_line"]
    numerics = config["numerics"]
    speed = _relativistic_speed(config)
    model_rows: List[Dict[str, str]] = []
    reference_rows: List[Dict[str, str]] = []
    commands: List[List[str]] = []
    for boundary in config["boundaries"]:
        for dt_s in numerics["time_steps_s"]:
            tag = f"{boundary}_dt_{float(dt_s):g}".replace(".", "p")
            # Each native invocation owns a separate directory because the
            # linked registry publishes a fixed, documented model filename.
            # Isolation also keeps the six native JSON/JUnit reports auditable.
            native_dir = output_dir / "native" / tag
            native_dir.mkdir(parents=True, exist_ok=True)
            native_input = native_dir / "CV01_native.args"
            model_path = native_dir / "CV01_model.csv"
            native_json = native_dir / "srcsep-native-tests.json"
            native_junit = native_dir / "srcsep-native-tests.xml"
            reference_path = output_dir / f"reference_{tag}.csv"
            native_arguments = [
                "--boundary", str(boundary),
                "--minimum-m", str(line["minimum_m"]),
                "--maximum-m", str(line["maximum_m"]),
                "--speed-m-per-s", format(speed, ".17g"),
                "--mass-kg", str(constants["proton_mass_kg"]),
                "--light-speed-m-per-s", str(constants["speed_of_light_m_per_s"]),
                "--plasma-speed-m-per-s", str(line["plasma_speed_m_per_s"]),
                "--dt-s", str(dt_s), "--final-time-s", str(numerics["final_time_s"]),
                "--sample-interval-s", str(numerics["sample_interval_s"]),
                "--campaign-seed", str(numerics["campaign_seed"]),
                "--initial-particles", str(particles_path),
            ]
            _write_native_arguments(native_input, native_arguments)
            model_command = [
                str(executable), "--test", CASE_ID,
                "--test-input", str(native_input),
                "--test-output-dir", str(native_dir),
                "--test-json", str(native_json),
                "--test-junit", str(native_junit),
            ]
            reference_command = [
                sys.executable, str(reference_source), "--input", str(resolved_input),
                "--initial-particles", str(particles_path), "--dt-s", str(dt_s),
                "--boundary", str(boundary), "--output", str(reference_path),
            ]
            _run(model_command, source_root, log_path, timeout=timeout)
            _validate_native_report(native_json, model_path)
            _run(reference_command, source_root, log_path, timeout=timeout)
            commands.extend((model_command, reference_command))
            model_rows.extend(_read_rows(model_path))
            reference_rows.extend(_read_rows(reference_path))

    metrics, passed = _comparison(config, model_rows, reference_rows, output_dir)
    moments_path = output_dir / "CV01_packet_moments.csv"
    profile_path = _write_profile(model_rows, reference_rows, config, output_dir)
    figures = _plot_four_panel(config, moments_path, profile_path, output_dir)

    provenance_path = output_dir / "provenance.json"
    _write_json(provenance_path, {
        "schema": "srcsep-validation-case-provenance-v1",
        "case_id": CASE_ID,
        "input": {"path": str(input_path), "sha256": _hash(input_path)},
        "resolved_input": {"path": str(resolved_input), "sha256": _hash(resolved_input)},
        "initial_particles": {"path": str(particles_path), "sha256": _hash(particles_path)},
        "linked_application": {"path": str(executable),
                               "sha256": _hash(executable)},
        "model_driver": {"path": str(driver_source), "sha256": _hash(driver_source),
                         "execution": "compiled-into-linked-srcsep-amps"},
        "production_sources": [
            {"path": str(source_root / "util" / "sep_transport_common.cpp"),
             "sha256": _hash(source_root / "util" / "sep_transport_common.cpp")},
            {"path": str(source_root / "util" / "sep_focused_transport_core.cpp"),
             "sha256": _hash(source_root / "util" / "sep_focused_transport_core.cpp")},
        ],
        "reference_generator": {"path": str(reference_source),
                                "sha256": _hash(reference_source),
                                "algorithm": "closed-form-characteristic-v1"},
        "commands": commands,
        "outputs": [{"path": str(path.relative_to(output_dir)),
                     "sha256": _hash(path)}
                    for path in sorted(output_dir.rglob("*"))
                    if path.is_file() and path != provenance_path],
    })

    relative = lambda path: str(Path(CASE_ID) / path.name)
    artifacts = [relative(output_dir / "CV01_solution.csv"),
                 relative(moments_path), relative(profile_path),
                 relative(resolved_input), relative(particles_path),
                 relative(log_path), relative(provenance_path)]
    artifacts.extend(str(path.relative_to(output_dir.parent))
                     for path in sorted((output_dir / "native").glob("*/CV01_model.csv")))
    artifacts.extend(str(path.relative_to(output_dir.parent))
                     for path in sorted((output_dir / "native").glob("*/srcsep-native-tests.*")))
    artifacts.extend(relative(path) for path in sorted(output_dir.glob("reference_*.csv")))
    artifacts.extend(relative(path) for path in figures)
    return {
        "id": CASE_ID,
        "status": "PASS" if passed else "FAIL",
        "message": ("linked srcSEP/AMPS focused-transport characteristics, packet "
                    "moments, boundaries, and weight ledger match the independent solution"
                    if passed else
                    "one or more ballistic-streaming acceptance metrics failed"),
        "elapsed_seconds": time.monotonic() - started,
        "seed": int(numerics["campaign_seed"]),
        "configuration": [
            "mover=fte-dmumu", "background=analytic-uniform",
            "equation_mode=full-gyrotropic", "Dmumu_s^-1=0",
            "dlnB_ds_m^-1=0", "divU_s^-1=0", "sources=false",
            "losses=false", "particle_wave_coupling=false",
            "boundaries=periodic,open", "mu=-1,-0.5,0.5,1",
            "time_steps_s=" + ",".join(str(value) for value in numerics["time_steps_s"]),
            "temporal_order=exact_for_constant_coefficients",
            f"input_sha256={_hash(input_path)}",
        ],
        "metrics": metrics,
        "artifacts": artifacts,
    }
