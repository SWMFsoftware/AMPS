#!/usr/bin/env python3
"""Validate immutable inputs and prepare absolute-path SA01 AMPS branches."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

from sa01_common import (
    PROFILES, SPACECRAFT, START_UTC, atomic_json, iter_data_lines, parse_utc,
    resolve_data_root, sha256_file, utc_now, validate_data_package,
    validate_spacecraft,
)


HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "AMPS_PARAM.in.template"


def boundary_epochs(path: Path) -> set[str]:
    """Read the packaged complete-spectrum epochs; never bridge absent boundary rows."""
    epochs = {line.split()[0] for line in iter_data_lines(path)}
    if len(epochs) < 2:
        raise ValueError(f"{path}: fewer than two complete boundary epochs")
    return epochs


def prepare_trajectory(source: Path, boundary: Path, destination: Path,
                       cadence_minutes: int) -> dict[str, Any]:
    """Select by frozen cadence and exact observed-boundary availability only."""

    available = boundary_epochs(boundary)
    start = parse_utc(START_UTC)
    kept: list[str] = []
    excluded: list[dict[str, str]] = []
    seen: set[str] = set()
    for line in iter_data_lines(source):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"{source}: malformed trajectory row {line!r}")
        utc = fields[0]
        if utc in seen:
            raise ValueError(f"{source}: duplicate trajectory epoch {utc}")
        seen.add(utc)
        elapsed_minutes = int((parse_utc(utc) - start).total_seconds() // 60)
        if elapsed_minutes % cadence_minutes:
            excluded.append({"utc": utc, "reason": "profile_cadence"})
        elif utc not in available:
            excluded.append({"utc": utc, "reason": "missing_exact_boundary_spectrum"})
        else:
            kept.append(line)
    if len(kept) < 12:
        raise ValueError(f"{source}: only {len(kept)} trajectory rows survive")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        "# UTC geocentric_latitude_deg longitude_deg altitude_km\n"
        f"# Frozen profile cadence: {cadence_minutes} minutes\n"
        "# Each epoch has an exact complete SOHO/ERNE five-minute boundary row\n"
        + "\n".join(kept) + "\n", encoding="ascii",
    )
    return {
        "source": str(source.resolve()), "output": str(destination.resolve()),
        "source_rows": len(seen), "kept": len(kept), "excluded": len(excluded),
        "excluded_by_reason": {
            reason: sum(row["reason"] == reason for row in excluded)
            for reason in sorted({row["reason"] for row in excluded})
        },
        "excluded_epochs": excluded, "first_epoch": kept[0].split()[0],
        "last_epoch": kept[-1].split()[0], "sha256": sha256_file(destination),
    }


def render_deck(replacements: dict[str, Any]) -> str:
    """Render the repository template and reject any unresolved placeholder."""
    text = TEMPLATE.read_text(encoding="utf-8")
    for key, value in replacements.items():
        text = text.replace("{{" + key + "}}", str(value))
    if "{{" in text or "}}" in text:
        raise ValueError("unresolved AMPS template token")
    return text


def prepare_branch(validation_root: Path, output_root: Path, spacecraft: str,
                   profile_name: str, threads: int, solver: str = "GRIDLESS",
                   mesh_res_earth_re: float = 0.05,
                   mesh_res_boundary_re: float = 2.0,
                   mesh_coarsening: str = "LINEAR",
                   mesh_exponent: float = 1.0,
                   numerical_controls: dict[str, Any] | None = None) -> dict[str, Any]:
    """Prepare one immutable spacecraft branch for direct or mesh-backed tracing.

    GRIDDED uses Mode3D's explicit SNAPSHOT_LIST contract.  A trajectory epoch is
    therefore evaluated only against its matching frozen field grid; ordinary
    TIME_SERIES semantics would incorrectly form an epoch-by-location Cartesian
    product.  Mesh controls are always written explicitly so a gridded result never
    depends on an undocumented application default.
    """

    solver = solver.upper()
    if solver not in ("GRIDLESS", "GRIDDED"):
        raise ValueError("solver must be GRIDLESS or GRIDDED")
    if mesh_res_earth_re <= 0.0 or mesh_res_boundary_re <= 0.0 or mesh_exponent <= 0.0:
        raise ValueError("Mode3D mesh resolution controls must be positive")
    mesh_coarsening = mesh_coarsening.upper()
    if mesh_coarsening not in ("LINEAR", "LOG", "EXPONENTIAL", "POWER", "CONSTANT"):
        raise ValueError(f"unsupported Mode3D mesh coarsening {mesh_coarsening!r}")
    controls = dict(numerical_controls or {})
    mover = controls.get("mover", "RK4")
    trace_dt = float(controls.get("trace_dt", 0.1))
    if mover not in ("BORIS", "HC4", "RK4", "RK6") or not math.isfinite(trace_dt) or trace_dt <= 0:
        raise ValueError("invalid full-orbit mover or trace timestep")
    controls.update(mover=mover, trace_dt=trace_dt)
    profile = dict(PROFILES[profile_name])
    if controls.get("cadence_minutes") is not None:
        cadence = int(controls["cadence_minutes"])
        if cadence < 5 or cadence % 5:
            raise ValueError("cadence must be a positive multiple of five minutes")
        profile["cadence_minutes"] = cadence
    spectrum_nodes = profile["ds_intervals"] + 1
    # Preserve QUICK/STANDARD's historical parent grid and deterministic samples.
    # PUBLICATION requires 1536 directions, exceeding the old hard-coded 1152
    # ceiling; explicitly enlarge the parent grid instead of silently truncating.
    direction_nmu = max(24, (profile["angular_samples"] + 47) // 48)
    total_particles = spectrum_nodes * profile["angular_samples"]
    branch = output_root / spacecraft
    if branch.exists() and any(branch.iterdir()):
        raise FileExistsError(
            f"{branch} is nonempty; select a new output directory to preserve prior evidence"
        )
    input_dir = branch / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    processed = validation_root / "processed" / "SA01"
    driver = (processed / "ts05_driver_5min.txt").resolve()
    spectrum = (processed / "boundary_soho_erne_5min.dat").resolve()
    trajectory = prepare_trajectory(
        processed / f"{spacecraft.lower()}_trajectory_geo_5min.txt", spectrum,
        input_dir / "trajectory_geo.txt", profile["cadence_minutes"],
    )
    snapshot_list = input_dir / "snapshot_epochs.txt"
    selected_epochs = [line.split()[0] for line in iter_data_lines(input_dir / "trajectory_geo.txt")]
    if solver == "GRIDDED":
        snapshot_list.write_text("\n".join(selected_epochs) + "\n", encoding="ascii")

    if solver == "GRIDDED":
        mode3d_mesh_block = "\n".join((
            "#MODE3D_MESH",
            f"  MODE3D_MESH_RES_EARTH_RE {mesh_res_earth_re:.12g}",
            f"  MODE3D_MESH_RES_BOUNDARY_RE {mesh_res_boundary_re:.12g}",
            f"  MODE3D_MESH_COARSENING {mesh_coarsening}",
            f"  MODE3D_MESH_EXPONENT {mesh_exponent:.12g}",
        ))
        temporal_block = "\n".join((
            "#TEMPORAL",
            "  TEMPORAL_MODE SNAPSHOT_LIST",
            f"  SNAPSHOT_LIST_FILE {snapshot_list.resolve()}",
        ))
        # The direct THREADS backend currently schedules complete Mode3D density
        # locations.  A SNAPSHOT_LIST trajectory contains one active location per
        # epoch, so OPENMP is selected deliberately: its inner energy/direction loops
        # use the requested workers while all of them read the same immutable grid.
        # SA01 also opts into EPOCH_PARALLEL explicitly.  This keeps its previously
        # qualified private-field execution model after the application default moved
        # to memory-safe COLLECTIVE; it must not be inferred from SNAPSHOT_LIST alone.
        parallel_block = "\n".join((
            "  DENSITY_PARALLEL OPENMP",
            f"  DENSITY_THREADS {threads}",
            "  MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER EPOCH_PARALLEL",
            "  MODE3D_EPOCH_PARALLEL_MEMORY_FRACTION 0.80",
            "  MODE3D_MPI_SCHEDULER DYNAMIC",
            "  MODE3D_MPI_DYNAMIC_CHUNK 1",
        ))
    else:
        mode3d_mesh_block = "! Mode3D mesh controls are not used by GRIDLESS."
        temporal_block = "! GRIDLESS evaluates the timestamp attached to each trajectory row."
        parallel_block = "\n".join((
            "  GRIDLESS_PARALLEL THREADS",
            f"  GRIDLESS_THREADS {threads}",
            "  GRIDLESS_MPI_SCHEDULER DYNAMIC",
            "  GRIDLESS_MPI_DYNAMIC_CHUNK 0",
        ))

    deck = render_deck({
        "RUN_ID": f"SA01_{spacecraft}_{profile_name}",
        "FIRST_EPOCH": trajectory["first_epoch"][:-1]
        if trajectory["first_epoch"].endswith("Z") else trajectory["first_epoch"],
        "DRIVER_FILE": driver, "SPECTRUM_FILE": spectrum,
        "TRAJECTORY_FILE": (input_dir / "trajectory_geo.txt").resolve(),
        "CUTOFF_NODES": profile["cutoff_nodes"],
        "DS_INTERVALS": profile["ds_intervals"],
        "SPECTRUM_NODES": spectrum_nodes,
        "TOTAL_PARTICLES": total_particles, "THREADS": threads,
        "FIELD_EVAL_METHOD": "GRID_3D" if solver == "GRIDDED" else "GRIDLESS",
        "MODE3D_MESH_BLOCK": mode3d_mesh_block,
        "TEMPORAL_BLOCK": temporal_block,
        "PARALLEL_BLOCK": parallel_block,
    })
    deck_path = input_dir / "AMPS_PARAM.in"
    deck = deck.replace("#DENSITY_SPECTRUM", "#DENSITY_SPECTRUM\n"
                        f"  DS_DIRECTION_NMU {direction_nmu}\n  DS_DIRECTION_NPHI 48")
    # The mover is a CLI option for backward tracing. Record it explicitly and
    # bind it to the immutable preparation instead of relying on executable defaults.
    deck = deck.replace("  DT_TRACE 0.1", f"  DT_TRACE {trace_dt:.17g}")
    deck_path.write_text(deck, encoding="utf-8")
    record = {
        "schema": "amps-sa01-prepared-run-v2", "created_utc": utc_now(),
        "case_id": "SA01", "spacecraft": spacecraft, "profile": profile_name,
        "solver": solver,
        "numerical_controls": controls,
        "profile_parameters": {
            **profile, "spectrum_nodes": spectrum_nodes,
            "total_trajectories_per_point_cap": total_particles,
            "requested_directions_per_energy": profile["angular_samples"],
            "direction_nmu": direction_nmu, "direction_nphi": 48,
        },
        "requested_threads_per_rank": threads,
        "input_file": str(deck_path.resolve()), "input_sha256": sha256_file(deck_path),
        "driver_file": str(driver), "driver_sha256": sha256_file(driver),
        "spectrum_file": str(spectrum), "spectrum_sha256": sha256_file(spectrum),
        "trajectory": trajectory, "status": "PREPARED",
        "snapshot_list_file": str(snapshot_list.resolve()) if solver == "GRIDDED" else None,
        "snapshot_list_sha256": sha256_file(snapshot_list) if solver == "GRIDDED" else None,
        "mode3d_mesh": ({
            "field_evaluation": "INTERPOLATION",
            "res_earth_re": mesh_res_earth_re,
            "res_boundary_re": mesh_res_boundary_re,
            "coarsening": mesh_coarsening,
            "exponent": mesh_exponent,
            "parallel_field_initialization": True,
            "trajectory_snapshot_scheduler": "EPOCH_PARALLEL",
            "epoch_parallel_memory_fraction": 0.80,
        } if solver == "GRIDDED" else None),
        # Record the closure policy independently of the rendered input hash so reports
        # can distinguish the preserved 300-s primary classification from the bounded
        # convergence work used to resolve only limit-terminated directions.
        "trajectory_closure": {
            "primary_trace_limit_seconds": 300.0,
            "unresolved_tolerance": 0.01,
            "fail_on_unresolved": True,
            "unresolved_extension_passes": 3,
            "unresolved_extension_factor": 2.0,
            "maximum_extension_limit_seconds": 2400.0,
            "drift_recurrence_detection": True,
            "drift_minimum_revolutions": 3,
            "drift_profile_bins": 24,
            "drift_minimum_profile_coverage": 0.70,
            "drift_minimum_matched_bin_fraction": 0.75,
            "drift_maximum_mean_radius_change_re": 1.0,
        },
        "science_contract": {
            "field": "standalone T05 with Qin-Denton five-minute drivers",
            "electric_field": "NONE", "boundary": "isotropic and spatially uniform",
            "energy_interval_mev": [18.5, 116.0], "maximum_spectrum_gap_seconds": 900,
        },
    }
    atomic_json(branch / "prepared_run.json", record)
    return record


def prepare_all(data_path: Path, output_dir: Path, profile: str, threads: int,
                spacecraft: tuple[str, ...], solver: str = "GRIDLESS",
                mesh_res_earth_re: float = 0.05,
                mesh_res_boundary_re: float = 2.0,
                mesh_coarsening: str = "LINEAR",
                mesh_exponent: float = 1.0,
                numerical_controls: dict[str, Any] | None = None) -> dict[str, Any]:
    package_root, validation_root = resolve_data_root(data_path)
    validation = validate_data_package(validation_root)
    output = output_dir.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"{output} is nonempty; use a new --output-dir (compare/plot modes may reuse existing output)"
        )
    output.mkdir(parents=True, exist_ok=True)
    atomic_json(output / "data_validation.json", validation)
    branches = [prepare_branch(
        validation_root, output, item, profile, threads, solver,
        mesh_res_earth_re, mesh_res_boundary_re, mesh_coarsening, mesh_exponent,
        numerical_controls,
    ) for item in spacecraft]
    summary = {
        "schema": "amps-sa01-preparation-summary-v2", "created_utc": utc_now(),
        "package_root": str(package_root), "validation_data_root": str(validation_root),
        "profile": profile, "threads": threads, "solver": solver.upper(),
        "spacecraft": list(spacecraft),
        "branches": branches, "status": "PASS",
    }
    atomic_json(output / "preparation_summary.json", summary)
    return summary


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    source = result.add_mutually_exclusive_group(required=True)
    source.add_argument("--data-root", type=Path)
    source.add_argument("--data-archive", type=Path,
                        help="Alias for --data-root; this must name an extracted directory, not a tar file.")
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--profile", choices=tuple(PROFILES), default="STANDARD")
    result.add_argument("--threads", type=int, default=1)
    result.add_argument("--solver", choices=("GRIDLESS", "GRIDDED"), default="GRIDLESS")
    result.add_argument("--mode3d-mesh-res-earth-re", type=float, default=0.05)
    result.add_argument("--mode3d-mesh-res-boundary-re", type=float, default=2.0)
    result.add_argument("--mode3d-mesh-coarsening",
                        choices=("LINEAR", "LOG", "EXPONENTIAL", "POWER", "CONSTANT"),
                        default="LINEAR")
    result.add_argument("--mode3d-mesh-exponent", type=float, default=1.0)
    result.add_argument("--spacecraft", nargs="+", default=list(SPACECRAFT))
    add_numerical_arguments(result)
    return result


def add_numerical_arguments(parser: argparse.ArgumentParser) -> None:
    """Share explicit numerical controls between preparation and unified CLI."""
    parser.add_argument("--mover", choices=("RK4", "BORIS", "HC4", "RK6"), default="RK4",
                        help="Full-orbit mover, saved in preparation and passed explicitly to AMPS.")
    parser.add_argument("--trace-dt", type=float, default=0.1,
                        help="Maximum integration step in seconds; an explicit refinement control.")
    parser.add_argument("--cadence-minutes", type=int,
                        help="Override cadence to keep the same population during numerical refinement.")


def numerical_options(args: argparse.Namespace) -> dict[str, Any]:
    """Persist explicit solver choices so run-only cannot silently switch numerics."""
    return {"mover": args.mover, "trace_dt": args.trace_dt,
            "cadence_minutes": args.cadence_minutes}


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    summary = prepare_all(args.data_root or args.data_archive, args.output_dir,
                          args.profile, args.threads, validate_spacecraft(args.spacecraft),
                          args.solver, args.mode3d_mesh_res_earth_re,
                          args.mode3d_mesh_res_boundary_re,
                          args.mode3d_mesh_coarsening, args.mode3d_mesh_exponent,
                          numerical_options(args))
    print(f"[PASS] prepared {len(summary['branches'])} SA01 branch(es) in {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, KeyError, ValueError) as error:
        print(f"ERROR: {error}")
        raise SystemExit(2)
