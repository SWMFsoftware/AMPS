#!/usr/bin/env python3
"""Authenticate inputs and prepare the model-independent GOES01 population."""

from __future__ import annotations

import bisect
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from goes01_common import (
    EVENT_END_EXCLUSIVE, EVENT_START, P8A_INDEX, PROFILES, TWO_HEAD_ENVELOPE_MAX,
    atomic_json, iso_utc, load_sgps, parse_utc, read_csv, resolve_model_inputs,
    resolve_observations, sha256_file, utc_now, write_csv,
)


HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "AMPS_PARAM.in.template"
NETCDF_EPOCH = datetime(2000, 1, 1, 12, tzinfo=timezone.utc)
HEADS = ((0, "WEST_MINUS_X"), (1, "EAST_PLUS_X"))


def _table_epochs(path: Path) -> list[datetime]:
    """Read a whitespace time-series axis and require strict chronological order."""
    epochs = [parse_utc(line.split()[0]) for line in path.read_text().splitlines()
              if line.strip() and not line.lstrip().startswith(("#", "!"))]
    if len(epochs) < 2 or any(right <= left for left, right in zip(epochs, epochs[1:])):
        raise ValueError(f"{path}: time axis is absent, duplicated, or unordered")
    return epochs


def _has_bracket(epoch: datetime, axis: list[datetime]) -> bool:
    """Require real rows on both sides; endpoint extrapolation is not validation."""
    index = bisect.bisect_right(axis, epoch)
    if index == 0 or index == len(axis):
        return False
    return ((epoch - axis[index - 1]).total_seconds() <= 900.0 and
            (axis[index] - epoch).total_seconds() <= 900.0)


def shift_boundary_to_midpoints(source: Path, destination: Path) -> None:
    """Relabel each five-minute ERNE average from bin start to bin midpoint.

    The package preprocessing code groups one-minute records with ``floor_time``;
    raw ERNE CDF metadata defines Epoch as start time. Values are therefore interval
    averages labeled by bin start, while GOES01 evaluates the same interval at its
    midpoint. Only timestamps change here. Every intensity token is copied exactly.
    """
    output = [
        "# GOES01 midpoint-labeled copy of the authenticated ERNE five-minute averages.",
        "# Each source timestamp is shifted by +150 s; intensity text is unchanged.",
    ]
    count = 0
    for raw in source.read_text(encoding="ascii").splitlines():
        if not raw.strip() or raw.lstrip().startswith(("#", "!")):
            output.append(raw)
            continue
        fields = raw.split(maxsplit=1)
        if len(fields) != 2:
            raise ValueError(f"{source}: malformed boundary row {raw!r}")
        output.append(iso_utc(parse_utc(fields[0]) + timedelta(seconds=150)) + " " + fields[1])
        count += 1
    if count < 2:
        raise ValueError(f"{source}: fewer than two boundary averages")
    destination.write_text("\n".join(output) + "\n", encoding="ascii")


def _ephemeris(path: Path) -> dict[str, dict[str, str]]:
    """Index the published half-minute positions, rejecting duplicate identities."""
    result: dict[str, dict[str, str]] = {}
    for row in read_csv(path):
        key = row["position_epoch_utc"]
        if key in result:
            raise ValueError(f"{path}: duplicate position epoch {key}")
        result[key] = row
    return result


def select_observations(netcdf: Path, ephemeris_path: Path, boundary: Path,
                        driver: Path, profile_name: str) -> tuple[list[dict[str, Any]],
                                                                  list[dict[str, Any]]]:
    """Freeze eligible epochs without consulting a modeled value.

    The two-head statistic is |Jw-Je|/(sw+se). Requiring a value no greater than one
    means the two delivered uncertainty intervals overlap. Unlike a quadrature
    denominator, this rule does not invent unavailable cross-head independence or
    covariance. It still demonstrates only consistency of two broad viewing cones,
    not a complete pitch-angle distribution.
    """
    data = load_sgps(netcdf)
    ephemeris = _ephemeris(ephemeris_path)
    # Boundary rows are start-labeled five-minute averages. The matching physical
    # value is exact at source-start + 150 s after preparation relabels the table.
    boundary_midpoints = {epoch + timedelta(seconds=150) for epoch in _table_epochs(boundary)}
    driver_axis = _table_epochs(driver)
    cadence = PROFILES[profile_name]["cadence_minutes"]
    records: list[dict[str, Any]] = []
    ledger: list[dict[str, Any]] = []

    # The NetCDF energy bounds are in keV. Keeping each head's actual P8A interval
    # prevents an undocumented nominal-energy substitution in later integration.
    lower = data["lower"][:, P8A_INDEX] / 1000.0
    upper = data["upper"][:, P8A_INDEX] / 1000.0
    expected = ((83.7, 98.5), (82.9, 99.8))
    for index, bounds in enumerate(zip(lower, upper)):
        if any(abs(float(value) - target) > 1.0e-6
               for value, target in zip(bounds, expected[index])):
            raise ValueError(f"unexpected P8A bounds for sensor {index}: {bounds}")

    for record_index, seconds in enumerate(data["time"]):
        start = NETCDF_EPOCH + timedelta(seconds=float(seconds))
        if start < EVENT_START or start >= EVENT_END_EXCLUSIVE:
            continue
        midpoint = start + timedelta(seconds=150)
        reasons: list[str] = []
        flux = data["flux"][record_index, :, P8A_INDEX]
        uncertainty = data["uncertainty"][record_index, :, P8A_INDEX]
        if int(data["records"][record_index]) != 300:
            reasons.append("not_300_l1b_records")
        if int(data["yaw"][record_index]) != 0:
            reasons.append("yaw_flip_not_zero")
        for array_name in ("valid", "dtc", "oob", "err"):
            values = data[array_name][record_index, :, P8A_INDEX]
            target = 300 if array_name == "valid" else 0
            if any(int(value) != target for value in values):
                reasons.append(f"{array_name}_quality")
        if (any(not math.isfinite(float(value)) or float(value) <= 0.0 for value in flux)
                or any(not math.isfinite(float(value)) or float(value) <= 0.0
                       for value in uncertainty)):
            reasons.append("nonpositive_or_fill_flux")
            consistency = None
        else:
            consistency = abs(float(flux[0] - flux[1])) / float(
                uncertainty[0] + uncertainty[1])
            if consistency > TWO_HEAD_ENVELOPE_MAX:
                reasons.append("two_head_inconsistent")
        ephem = ephemeris.get(iso_utc(midpoint))
        if ephem is None:
            reasons.append("missing_exact_midpoint_ephemeris")
        elif ephem["ephemeris_status"] != "good" or int(ephem["ephemeris_flag"]) != 0:
            reasons.append("ephemeris_not_good")
        if midpoint not in boundary_midpoints:
            reasons.append("missing_exact_boundary_average")
        if not _has_bracket(midpoint, driver_axis):
            reasons.append("driver_not_bracketed")
        elapsed_minutes = int((start - EVENT_START).total_seconds() // 60)
        if elapsed_minutes % cadence:
            reasons.append("profile_cadence")
        ledger.append({
            "record_start_utc": iso_utc(start), "comparison_midpoint_utc": iso_utc(midpoint),
            "two_head_envelope_separation": consistency if consistency is not None else "",
            "eligible": int(not reasons), "exclusion_reasons": ";".join(reasons),
        })
        if reasons:
            continue
        assert ephem is not None and consistency is not None
        for head_index, head_name in HEADS:
            records.append({
                "record_start_utc": iso_utc(start),
                "comparison_midpoint_utc": iso_utc(midpoint),
                "head": head_name,
                "energy_min_mev": float(lower[head_index]),
                "energy_max_mev": float(upper[head_index]),
                # NOAA differential intensity is per keV; AMPS is per MeV.
                "observed_intensity_cm2_s_sr_mev": float(flux[head_index]) * 1000.0,
                "observed_uncertainty_cm2_s_sr_mev": float(uncertainty[head_index]) * 1000.0,
                "two_head_envelope_separation": consistency,
                "geocentric_latitude_deg": float(ephem["geocentric_latitude_deg"]),
                "longitude_deg": float(ephem["longitude_deg"]),
                "altitude_km": float(ephem["spherical_altitude_re6371p2_km"]),
                "gse_local_time_hours": float(ephem["gse_local_time_hours"]),
            })
    epoch_count = len({row["comparison_midpoint_utc"] for row in records})
    if epoch_count < 12:
        raise ValueError(f"only {epoch_count} independently qualified comparison epochs")
    return records, ledger


def _render(replacements: dict[str, Any]) -> str:
    text = TEMPLATE.read_text(encoding="utf-8")
    for name, value in replacements.items():
        text = text.replace("{{" + name + "}}", str(value))
    if "{{" in text or "}}" in text:
        raise ValueError("unresolved GOES01 AMPS template token")
    return text


def prepare(observation_root: Path, model_root: Path, output_root: Path,
            profile_name: str, threads: int, mover: str, trace_dt: float,
            mesh_earth_re: float, mesh_boundary_re: float) -> dict[str, Any]:
    """Create one non-overwritable GOES16 physical branch and its audit ledger."""
    if profile_name not in PROFILES:
        raise ValueError(f"unknown profile {profile_name}")
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty campaign {output_root}")
    if mover not in ("BORIS", "HC4", "RK4", "RK6") or trace_dt <= 0.0:
        raise ValueError("invalid mover or trace timestep")
    event_root, netcdf, ephemeris = resolve_observations(observation_root)
    processed, boundary_source, driver = resolve_model_inputs(model_root)
    observations, ledger = select_observations(
        netcdf, ephemeris, boundary_source, driver, profile_name
    )
    branch = output_root / "GOES16"
    input_dir = branch / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    boundary = input_dir / "boundary_soho_erne_5min_midpoints.dat"
    shift_boundary_to_midpoints(boundary_source, boundary)
    write_csv(input_dir / "selected_observations.csv", observations)
    write_csv(input_dir / "epoch_selection_ledger.csv", ledger)

    epoch_rows: dict[str, dict[str, Any]] = {}
    for row in observations:
        epoch_rows.setdefault(row["comparison_midpoint_utc"], row)
    trajectory = input_dir / "goes16_trajectory_geo.txt"
    trajectory.write_text(
        "# UTC geocentric_latitude_deg longitude_deg altitude_km\n"
        "# Epochs are midpoints of five-minute SGPS averages; positions are exact.\n" +
        "\n".join(
            f"{utc} {row['geocentric_latitude_deg']:.12g} "
            f"{row['longitude_deg']:.12g} {row['altitude_km']:.12g}"
            for utc, row in sorted(epoch_rows.items())
        ) + "\n", encoding="ascii")
    snapshots = input_dir / "snapshot_epochs.txt"
    snapshots.write_text("\n".join(sorted(epoch_rows)) + "\n", encoding="ascii")

    profile = PROFILES[profile_name]
    spectrum_nodes = profile["ds_intervals"] + 1
    direction_nmu = max(24, math.ceil(profile["angular_samples"] / 48))
    deck = _render({
        "RUN_ID": f"GOES01_GOES16_{profile_name}",
        "FIRST_EPOCH": sorted(epoch_rows)[0].rstrip("Z"),
        "DRIVER_FILE": driver, "SPECTRUM_FILE": boundary,
        "TRAJECTORY_FILE": trajectory.resolve(), "SNAPSHOT_FILE": snapshots.resolve(),
        "DS_INTERVALS": profile["ds_intervals"], "SPECTRUM_NODES": spectrum_nodes,
        "TOTAL_PARTICLES": spectrum_nodes * profile["angular_samples"],
        "DIRECTION_NMU": direction_nmu, "THREADS": threads,
        "TRACE_DT": f"{trace_dt:.17g}", "MESH_EARTH_RE": mesh_earth_re,
        "MESH_BOUNDARY_RE": mesh_boundary_re,
    })
    deck_path = input_dir / "AMPS_PARAM.in"
    deck_path.write_text(deck, encoding="utf-8")
    record = {
        "schema": "amps-goes01-prepared-v1", "created_utc": utc_now(),
        "case_id": "GOES01", "spacecraft": "GOES16", "profile": profile_name,
        "profile_parameters": {**profile, "spectrum_nodes": spectrum_nodes,
                               "requested_directions_per_energy": profile["angular_samples"],
                               "direction_nmu": direction_nmu, "direction_nphi": 48},
        "solver": "GRIDDED", "numerical_controls": {"mover": mover, "trace_dt": trace_dt},
        "requested_threads_per_rank": threads,
        "input_file": str(deck_path.resolve()), "input_sha256": sha256_file(deck_path),
        "driver_file": str(driver), "driver_sha256": sha256_file(driver),
        "spectrum_file": str(boundary), "spectrum_sha256": sha256_file(boundary),
        "trajectory": {"output": str(trajectory.resolve()), "sha256": sha256_file(trajectory),
                       "kept": len(epoch_rows), "first_epoch": sorted(epoch_rows)[0],
                       "last_epoch": sorted(epoch_rows)[-1]},
        "snapshot_list_file": str(snapshots.resolve()),
        "snapshot_list_sha256": sha256_file(snapshots),
        "observations_file": str((input_dir / "selected_observations.csv").resolve()),
        "observations_sha256": sha256_file(input_dir / "selected_observations.csv"),
        "selection_ledger_file": str((input_dir / "epoch_selection_ledger.csv").resolve()),
        "selection_ledger_sha256": sha256_file(input_dir / "epoch_selection_ledger.csv"),
        "source": {"event_root": str(event_root), "processed_model_root": str(processed),
                   "netcdf": str(netcdf), "netcdf_sha256": sha256_file(netcdf),
                   "ephemeris": str(ephemeris), "ephemeris_sha256": sha256_file(ephemeris),
                   "boundary_start_labeled": str(boundary_source),
                   "boundary_start_labeled_sha256": sha256_file(boundary_source),
                   "boundary_time_transform": "timestamp + 150 seconds; intensities unchanged"},
        "mode3d_mesh": {"res_earth_re": mesh_earth_re,
                        "res_boundary_re": mesh_boundary_re,
                        "coarsening": "LINEAR", "exponent": 1.0},
        "science_contract": {
            "observable": "GOES-16 SGPS P8A differential proton intensity",
            "event_interval": "2017-09-11T07:30:00Z/2017-09-14T00:00:00Z",
            "comparison_time": "midpoint of each five-minute SGPS average",
            "quality": "300/300 samples, zero DQ sums, yaw flag zero, good exact ephemeris",
            "direction_policy": "both heads retained only when delivered uncertainty intervals overlap",
            "two_head_envelope_separation_max": TWO_HEAD_ENVELOPE_MAX,
            "claims_excluded": ["particle density", "directional response",
                                "SEP acceleration", "heliospheric propagation", "onset timing"],
        },
    }
    atomic_json(branch / "prepared_run.json", record)
    atomic_json(output_root / "preparation_summary.json", {
        "schema": "amps-goes01-preparation-summary-v1", "status": "PASS",
        "profile": profile_name, "modeled_epochs": len(epoch_rows),
        "comparison_rows": len(observations), "branch": str(branch.resolve()),
    })
    return record
