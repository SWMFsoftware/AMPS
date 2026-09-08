#!/usr/bin/env python3
"""Run the 15-min December 2006 AMPS cutoff-morphology experiment.

This runner extends the observation-specific C9/C10 calculations to a common
rigidity grid and two fixed altitude shells.  It deliberately imports the C10
Tecplot parser, GEO-to-AACGM conversion, and ACCESS_T50 reducer; therefore the
research product and the POES/MetOp validation cannot drift into different
definitions of access state or half transmission.

The script supports a preparation-only mode.  That mode renders every AMPS
input, writes the complete command inventory, and validates the requested epochs
without launching MPI.  This is useful on login nodes and for advance review of
the computational cost.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import shutil
import sys
from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

from study_common import (
    DriverRow, default_output_root, driver_dict, format_utc,
    interpolate_driver, load_config, parse_utc, read_driver,
    resolve_output_path, write_csv,
)


def load_c10_module(root: Path):
    """Load the bundled C10 runner as a private scientific-method library."""
    path = root / "vendor" / "C10" / "run_C10.py"
    spec = importlib.util.spec_from_file_location("dec2006_c10_core", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load C10 implementation from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def cadence_epochs(start: datetime, end: datetime, minutes: int) -> List[datetime]:
    """Return an inclusive, UTC-aligned cadence sequence."""
    if minutes < 1:
        raise ValueError("cadence must be positive")
    result: List[datetime] = []
    current = start
    step = timedelta(minutes=minutes)
    while current <= end:
        result.append(current)
        current += step
    if result[-1] != end:
        result.append(end)
    return result


def theil_sen_slope(rows: Sequence[DriverRow]) -> float:
    """Robust SYM-H slope in nT/h for the objective quiet-window test."""
    slopes: List[float] = []
    for left_index, left in enumerate(rows):
        for right in rows[left_index + 1:]:
            hours = (right.epoch - left.epoch).total_seconds() / 3600.0
            if hours > 0:
                slopes.append((right.symh_nt - left.symh_nt) / hours)
    slopes.sort()
    if not slopes:
        return float("inf")
    middle = len(slopes) // 2
    return slopes[middle] if len(slopes) % 2 else 0.5 * (
        slopes[middle - 1] + slopes[middle]
    )


def event_landmarks(driver: Sequence[DriverRow], config: Mapping[str, object]) -> Dict[str, datetime]:
    """Select quiet, compression, and main-phase times from fixed rules."""
    event = config["event"]  # type: ignore[index]
    search_start = parse_utc(event["compression_search_start_utc"])  # type: ignore[index]
    search_end = parse_utc(event["compression_search_end_utc"])  # type: ignore[index]
    arrival = [row for row in driver if search_start <= row.epoch <= search_end]
    if len(arrival) < 2:
        raise ValueError("compression search window is not covered by the driver")
    compression = max(
        zip(arrival[1:], arrival[:-1]),
        key=lambda pair: math.log(pair[0].pdyn_npa) - math.log(pair[1].pdyn_npa),
    )[0].epoch

    production_start = parse_utc(event["production_start_utc"])  # type: ignore[index]
    production_end = parse_utc(event["production_end_utc"])  # type: ignore[index]
    production = [row for row in driver if production_start <= row.epoch <= production_end]
    main_phase = min(production, key=lambda row: row.symh_nt).epoch

    window_hours = float(event["quiet_window_hours"])  # type: ignore[index]
    max_median = float(event["quiet_median_abs_symh_max_nt"])  # type: ignore[index]
    max_slope = float(event["quiet_abs_theil_sen_slope_max_nt_per_hour"])  # type: ignore[index]
    window = timedelta(hours=window_hours)
    candidates: List[Tuple[datetime, List[DriverRow]]] = []
    for end_row in driver:
        if end_row.epoch >= compression:
            break
        start_time = end_row.epoch - window
        subset = [row for row in driver if start_time <= row.epoch <= end_row.epoch]
        expected = int(round(window.total_seconds() / 300.0)) + 1
        if len(subset) == expected:
            candidates.append((end_row.epoch, subset))
    eligible = []
    for end_time, subset in candidates:
        absolute = sorted(abs(row.symh_nt) for row in subset)
        median = absolute[len(absolute) // 2]
        if median <= max_median and abs(theil_sen_slope(subset)) <= max_slope:
            eligible.append((end_time, subset))
    if not eligible:
        raise ValueError("no driver interval satisfies the configured quiet-window rule")
    quiet_rows = eligible[-1][1]
    quiet = quiet_rows[len(quiet_rows) // 2].epoch
    return {"quiet": quiet, "compression": compression, "main_phase": main_phase}


def selected_epochs(config: Mapping[str, object], profile: str,
                    landmarks: Mapping[str, datetime]) -> List[datetime]:
    """Construct SMOKE, ROUTINE, or FULL epochs and include exact landmarks."""
    profile_cfg = config["profiles"][profile]  # type: ignore[index]
    if "explicit_epochs_utc" in profile_cfg:
        return sorted({parse_utc(value) for value in profile_cfg["explicit_epochs_utc"]})
    event = config["event"]  # type: ignore[index]
    start = parse_utc(event["production_start_utc"])  # type: ignore[index]
    end = parse_utc(event["production_end_utc"])  # type: ignore[index]
    result = set(cadence_epochs(start, end, int(profile_cfg["cadence_minutes"])))
    if profile_cfg.get("include_landmarks", False):
        result.update(landmarks.values())
    if "rapid_window_half_width_minutes" in profile_cfg:
        half_width = int(profile_cfg["rapid_window_half_width_minutes"])
        rapid_cadence = int(profile_cfg["rapid_window_cadence_minutes"])
        for center in (landmarks["compression"], landmarks["main_phase"]):
            result.update(cadence_epochs(
                max(start, center - timedelta(minutes=half_width)),
                min(end, center + timedelta(minutes=half_width)),
                rapid_cadence,
            ))
    return sorted(result)


def runner_namespace(config: Mapping[str, object], args: argparse.Namespace,
                     altitude_km: float) -> SimpleNamespace:
    """Translate the study configuration to the exact controls used by C10."""
    model = config["model"]  # type: ignore[index]
    execution = config["execution"]  # type: ignore[index]
    return SimpleNamespace(
        rigidity_min_gv=min(model["rigidities_gv"]),
        rigidity_max_gv=max(model["rigidities_gv"]),
        cutoff_evaluation="DIRECT_ACCESS",
        cutoff_scan_n=120,
        rigidities_gv=list(model["rigidities_gv"]),
        access_abs_lat_min_deg=float(model["latitude_band_abs_deg"][0]),
        access_abs_lat_max_deg=float(model["latitude_band_abs_deg"][1]),
        max_trace_time=float(model["max_trace_time_s"]),
        altitude_km=altitude_km,
        shell_lon_res_deg=float(model["shell_longitude_step_deg"]),
        shell_lat_res_deg=float(model["shell_latitude_step_deg"]),
        mode3d_mesh_res_earth_re=float(model["mode3d_mesh_res_earth_re"]),
        mode3d_mesh_res_boundary_re=float(model["mode3d_mesh_res_boundary_re"]),
        mode3d_mesh_coarsening=str(model["mode3d_mesh_coarsening"]),
        mode3d_mesh_exponent=float(model["mode3d_mesh_exponent"]),
        mpirun=args.mpirun,
        np=args.np if args.np is not None else int(execution["mpi_ranks"]),
        nt=args.nt if args.nt is not None else int(execution["threads_per_rank"]),
        scheduler=str(execution["scheduler"]),
        dynamic_chunk=int(execution["dynamic_chunk"]),
        cutoff_trace_policy="ACCURATE",
        mode3d_parallel_field_init=bool(execution["parallel_field_initialization"]),
        mover=str(model["mover"]),
    )


def command_text(command: Sequence[str]) -> str:
    """Make a readable command record without shell-dependent quoting tricks."""
    import shlex
    return " ".join(shlex.quote(token) for token in command)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--driver", type=Path,
                        help="Override the configured driver for a labeled TS05 sensitivity run")
    parser.add_argument("--profile", choices=("SMOKE", "ROUTINE", "FULL"), default="SMOKE")
    parser.add_argument("--amps", type=Path, default=Path("./amps"))
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("-np", type=int, default=None)
    parser.add_argument("-nt", type=int, default=None)
    parser.add_argument(
        "--output-root", type=Path, default=default_output_root() / "morphology",
        help="Morphology output directory beneath the shared study run root",
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-run", action="store_true",
                        help="Parse existing AMPS outputs without launching the executable")
    parser.add_argument("--keep", action="store_true",
                        help="Do not overwrite an existing per-epoch directory")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root, config = load_config(args.config)
    c10 = load_c10_module(root)
    driver_path = args.driver or (root / config["data"]["driver"])  # type: ignore[index]
    if not driver_path.is_absolute():
        driver_path = (root / driver_path).resolve()
    driver = read_driver(driver_path)
    landmarks = event_landmarks(driver, config)
    epochs = selected_epochs(config, args.profile, landmarks)
    for epoch in epochs:
        interpolate_driver(driver, epoch)  # coverage and interpolation guard

    output_root = resolve_output_path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "event_landmarks.json").write_text(json.dumps(
        {name: format_utc(value) for name, value in landmarks.items()}, indent=2
    ) + "\n", encoding="utf-8")

    commands: List[Dict[str, object]] = []
    boundaries: List[Dict[str, object]] = []
    driver_samples: List[Dict[str, object]] = []
    failures: List[str] = []
    model = config["model"]  # type: ignore[index]
    mlt_bins = [3.0 * index for index in range(8)]
    n_case_slots = len(epochs) * len(model["shell_altitudes_km"])
    # With --keep, existing access products are reused and therefore are not
    # counted as launches remaining in this invocation.  This makes the counter
    # truthful when resuming a partially completed production calculation.
    if args.keep and not (args.prepare_only or args.skip_run):
        n_launches = sum(
            not (
                output_root
                / f"alt_{float(altitude):g}km/{epoch.strftime('%Y%m%dT%H%M%S')}"
                / "cutoff_3d_shells_access.dat"
            ).exists()
            for epoch in epochs
            for altitude in model["shell_altitudes_km"]
        )
    else:
        n_launches = n_case_slots
    launch_index = 0
    print(
        "Morphology execution plan: "
        f"{len(epochs)} epoch(s) x {len(model['shell_altitudes_km'])} "
        f"altitude(s) = {n_case_slots} case(s); "
        f"{n_launches} AMPS launch(es) required by this invocation",
        flush=True,
    )

    amps = args.amps.expanduser()
    if not amps.is_absolute():
        amps = (Path.cwd() / amps).resolve()
    if not (args.prepare_only or args.skip_run):
        if not amps.is_file() or not os.access(amps, os.X_OK):
            raise SystemExit(f"AMPS executable is missing or not executable: {amps}")

    for epoch in epochs:
        driver_samples.append(driver_dict(interpolate_driver(driver, epoch)))
        for altitude in model["shell_altitudes_km"]:
            altitude = float(altitude)
            controls = runner_namespace(config, args, altitude)
            tag = f"alt_{altitude:g}km/{epoch.strftime('%Y%m%dT%H%M%S')}"
            case_dir = output_root / tag
            case_dir.mkdir(parents=True, exist_ok=True)
            input_name = f"AMPS_PARAM_DEC2006_{int(round(altitude))}km.in"
            template = root / "inputs" / input_name
            local_driver = case_dir / "ts05_driver.txt"
            if not local_driver.exists():
                shutil.copy2(driver_path, local_driver)
            generated_input = case_dir / "AMPS_PARAM_C10.in"
            if not args.skip_run:
                c10.render_input(template, generated_input, epoch, local_driver, controls, "GRIDDED")
            command = c10.command_for(controls, amps, "GRIDDED", epoch)
            commands.append({
                "epoch_utc": format_utc(epoch), "altitude_km": altitude,
                "cwd": str(case_dir), "command": command,
                "command_line": command_text(command),
            })
            print(f"[{tag}] {command_text(command)}")
            if args.prepare_only:
                continue
            access_path = case_dir / "cutoff_3d_shells_access.dat"
            if not args.skip_run:
                if args.keep and access_path.exists():
                    print(f"[{tag}] keeping existing AMPS output")
                else:
                    launch_index += 1
                    # Each epoch/altitude pair is currently a separate AMPS
                    # launch.  Show completed and remaining work before MPI
                    # starts so users can estimate progress from the terminal
                    # or scheduler log without opening command_inventory.json.
                    print(
                        "Morphology AMPS progress: "
                        f"completed={launch_index - 1}/{n_launches}; "
                        f"remaining={n_launches - launch_index + 1}; "
                        f"starting={launch_index}/{n_launches}; case={tag}",
                        flush=True,
                    )
                    # Reuse C10's live tee because morphology intentionally uses
                    # the same AMPS execution and observation-operator code.  It
                    # writes the full per-case log while showing MPI, field-init,
                    # and trajectory diagnostics on screen as they occur.
                    return_code = c10.run_process(
                        command, case_dir, case_dir / "AMPS.log"
                    )
                    if return_code != 0:
                        failures.append(f"{tag}: AMPS exited with {return_code}")
                        continue
            if not access_path.exists():
                failures.append(f"{tag}: missing {access_path.name}")
                continue
            try:
                access = c10.parse_tecplot_shell_access(access_path)
                access = c10.select_common_access_band(
                    access, controls.access_abs_lat_min_deg,
                    controls.access_abs_lat_max_deg,
                )
                c10.add_aacgm_lat_mlt(access, epoch, altitude)
                estimates, profile_rows = c10.estimate_access_t50_boundaries(
                    access, controls.rigidities_gv, mlt_bins, ("N", "S"), 8,
                    0.25, 0.66, 1.0,
                )
                unresolved = sum(row.access_state == 2 for row in access)
                unresolved_fraction = unresolved / len(access) if access else 1.0
                for estimate in estimates:
                    for mlt, boundary in sorted(estimate.boundary_by_mlt.items()):
                        boundaries.append({
                            "epoch_utc": format_utc(epoch),
                            "altitude_km": altitude,
                            "rigidity_gv": estimate.rigidity_gv,
                            "hemisphere": estimate.hemisphere,
                            "mlt_hour": mlt,
                            "boundary_aacgm_lat_deg": boundary,
                            "n_valid_mlt": estimate.n_valid_mlt,
                            "n_requested_mlt": estimate.n_requested_mlt,
                            "unresolved_access_fraction": unresolved_fraction,
                            "field_model": "IGRF+TS05",
                            "observation_operator": "VERTICAL_ACCESS_T50",
                        })
                c10.write_dict_rows(case_dir / "snapshot_boundaries.csv", [
                    c10._estimate_row(estimate) for estimate in estimates
                ])
                c10.write_dict_rows(case_dir / "snapshot_t50_profiles.csv", profile_rows)
            except Exception as exc:  # keep other expensive cases usable
                failures.append(f"{tag}: postprocessing failed: {exc}")

    (output_root / "command_inventory.json").write_text(
        json.dumps(commands, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(output_root / "driver_at_model_epochs.csv", driver_samples)
    if boundaries:
        write_csv(output_root / "morphology_boundaries.csv", boundaries)
    result = {
        "study_id": config["study_id"], "profile": args.profile,
        "n_epochs": len(epochs), "n_altitudes": len(model["shell_altitudes_km"]),
        "n_cases": len(commands), "n_boundary_rows": len(boundaries),
        "prepare_only": args.prepare_only, "skip_run": args.skip_run,
        "failures": failures, "passed": not failures,
    }
    (output_root / "morphology_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
