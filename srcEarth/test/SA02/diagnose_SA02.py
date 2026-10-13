#!/usr/bin/env python3
"""Reproduce selected SA02 epochs without replacing the original failure evidence.

This diagnostic keeps the complete energy/direction grid and strict unresolved
gate. Optional mover, timestep and extension changes are recorded as numerical
experiments, never used to select observations or relax acceptance requirements.
"""
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess

from sa02_common import atomic_json, sha256_file, utc_now


def setting(deck: str, key: str) -> str:
    """Read a unique active keyword; ambiguous decks are not valid controls."""
    values = re.findall(r"^\s*" + re.escape(key) + r"\s+([^!\n]+)", deck, re.M)
    if len(values) != 1:
        raise ValueError(f"expected exactly one {key}, found {len(values)}")
    return values[0].strip()


def replace_setting(deck: str, key: str, value: object) -> str:
    """Change exactly one verified keyword without touching unrelated physics."""
    setting(deck, key)
    return re.sub(r"^(\s*" + re.escape(key) + r"\s+)[^!\n]+",
                  lambda match: match[1] + str(value), deck, flags=re.M)


def override_directions_per_energy(deck: str, requested: int,
                                   requested_n_mu: int = 0,
                                   requested_n_phi: int = 0) -> tuple[str, dict[str, int]]:
    """Change the explicit parent grid/sample count and describe the exact change.

    ``DS_MAX_PARTICLES`` is a cap for the complete energy-by-direction population,
    not a direction count.  Computing it from the deck's existing energy grid is
    essential: assigning ``requested`` directly would silently reduce the number of
    directions at every energy. By default the parent ``NMu x NPhi`` grid is retained,
    so the experiment changes only how many approximately uniformly spaced parent-grid
    directions are selected. A caller may provide both parent dimensions to continue
    refinement after exhausting the source grid; those dimensions and both old/new
    capacities are then recorded explicitly. Refusing a partial grid override,
    non-integral source population, or request larger than the selected parent grid
    keeps a malformed diagnostic from looking like valid quadrature evidence.
    """
    if requested < 1:
        raise ValueError("directions per energy must be positive")
    intervals = int(setting(deck, "DS_NINTERVALS"))
    energy_nodes = intervals + 1
    source_n_mu = int(setting(deck, "DS_DIRECTION_NMU"))
    source_n_phi = int(setting(deck, "DS_DIRECTION_NPHI"))
    if (requested_n_mu < 0 or requested_n_phi < 0 or
            bool(requested_n_mu) != bool(requested_n_phi)):
        raise ValueError("direction parent-grid override requires two positive dimensions")
    target_n_mu = requested_n_mu or source_n_mu
    target_n_phi = requested_n_phi or source_n_phi
    if target_n_mu < 1 or target_n_phi < 1:
        raise ValueError("direction parent-grid dimensions must be positive")
    source_parent_directions = source_n_mu * source_n_phi
    target_parent_directions = target_n_mu * target_n_phi
    source_particles = int(setting(deck, "DS_MAX_PARTICLES"))
    if source_particles % energy_nodes:
        raise ValueError(
            "source DS_MAX_PARTICLES is not an integral direction population "
            f"({source_particles} particles over {energy_nodes} energy nodes)"
        )
    if requested > target_parent_directions:
        raise ValueError(
            f"requested {requested} directions per energy exceeds the explicit "
            f"{target_n_mu}x{target_n_phi} parent grid ({target_parent_directions})"
        )
    source_directions = source_particles // energy_nodes
    updated_particles = requested * energy_nodes
    # Rewrite the dimensions even when they equal the source values. This keeps the
    # helper's behavior uniform and still relies on replace_setting() to reject an
    # ambiguous/missing active keyword.
    deck = replace_setting(deck, "DS_DIRECTION_NMU", target_n_mu)
    deck = replace_setting(deck, "DS_DIRECTION_NPHI", target_n_phi)
    deck = replace_setting(deck, "DS_MAX_PARTICLES", updated_particles)
    return deck, {
        "energy_nodes": energy_nodes,
        "source_direction_nmu": source_n_mu,
        "source_direction_nphi": source_n_phi,
        "source_parent_direction_count": source_parent_directions,
        "requested_direction_nmu": target_n_mu,
        "requested_direction_nphi": target_n_phi,
        "requested_parent_direction_count": target_parent_directions,
        "source_directions_per_energy": source_directions,
        "requested_directions_per_energy": requested,
        "source_max_particles": source_particles,
        "requested_max_particles": updated_particles,
    }


def launch_command(amps: Path, ranks: int, threads: int, backend: str,
                   field_init: str, mover: str | None, mpi_binding: str = "DEFAULT") -> list[str]:
    """Keep field construction and orbit-worker concurrency independently testable.

    A single density worker alone does not disable the separate field initializer:
    its caller can participate alongside an initialization worker. Explicit SERIAL
    field initialization omits the enabling flag, using the production default.
    No diagnostic path substitutes a different field or trajectory implementation.
    """
    command = ["mpiexec", "-n", str(ranks)]
    # Open MPI's small-job default core binding can make independent experiments
    # contend on the same cores even while the rest of the allocation is idle.
    # This explicit optional layout control is recorded; it never changes numerics.
    if mpi_binding == "NONE":
        command += ["--bind-to", "none"]
    command += [str(amps), "-mode", "3d",
               "-i", "input/AMPS_PARAM.in", "-mode3d-field-eval", "INTERPOLATION",
               "-density-parallel", backend, "-density-threads", str(threads)]
    if field_init == "PARALLEL":
        command.append("-mode3d-parallel-field-init")
    if mover:
        command += ["-mover", mover]
    return command


def main() -> int:
    """Create a fresh, labelled experiment and retain its full execution evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-branch", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--epochs", nargs="+", required=True)
    parser.add_argument("--amps", type=Path, default=Path("./amps"))
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--mpi-binding", choices=("DEFAULT", "NONE"), default="DEFAULT",
                        help="Optional Open MPI CPU binding override; NONE avoids small-job core contention")
    parser.add_argument("--field-init", choices=("PARALLEL", "SERIAL"), default="PARALLEL",
                        help="Vary field construction independently of trajectory workers")
    parser.add_argument("--density-backend", choices=("OPENMP", "THREADS", "SERIAL"), default="OPENMP",
                        help="Production trajectory backend; SERIAL isolates worker concurrency")
    parser.add_argument("--mover", choices=("RK4", "BORIS", "HC4", "RK6"))
    parser.add_argument("--dt", type=float)
    parser.add_argument("--directions-per-energy", type=int,
                        help="Diagnostic angular quadrature size; retains the source parent grid")
    parser.add_argument("--direction-grid", nargs=2, type=int, metavar=("NMU", "NPHI"),
                        help="Explicit replacement parent grid; requires --directions-per-energy")
    parser.add_argument("--extension-passes", type=int)
    args = parser.parse_args()
    if args.mpi_ranks < 1 or args.threads < 1 or (args.dt is not None and (not math.isfinite(args.dt) or args.dt <= 0)):
        parser.error("resource counts and timestep must be positive")
    if args.extension_passes is not None and not 0 <= args.extension_passes <= 8:
        parser.error("extension passes must be in [0,8]")
    if args.directions_per_energy is not None and args.directions_per_energy < 1:
        parser.error("directions per energy must be positive")
    if args.direction_grid is not None:
        if args.directions_per_energy is None:
            parser.error("--direction-grid requires --directions-per-energy")
        if any(value < 1 for value in args.direction_grid):
            parser.error("direction-grid dimensions must be positive")
    source = args.source_branch.resolve()
    output = args.output_dir.resolve()
    amps = args.amps.resolve()
    deck_path = source / "input/AMPS_PARAM.in"
    deck = deck_path.read_text()
    trajectory = Path(setting(deck, "TRAJ_FILE"))
    epochs = [epoch.rstrip("Z") for epoch in args.epochs]
    if len(epochs) != len(set(epochs)):
        raise ValueError("duplicate diagnostic epoch")
    rows = [row for row in trajectory.read_text().splitlines()
            if row.strip() and row.split()[0].rstrip("Z") in epochs]
    if len(rows) != len(epochs) or {row.split()[0].rstrip("Z") for row in rows} != set(epochs):
        raise ValueError("requested epoch absent or duplicated in the preserved trajectory")
    if setting(deck, "DS_UNRESOLVED_TOL") != "0.01" or setting(deck, "DS_FAIL_ON_UNRESOLVED") != "T":
        raise ValueError("diagnostic source must retain the strict SA02 gate")
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "input"
    inputs.mkdir()
    (inputs / "trajectory_geo.txt").write_text("\n".join(rows) + "\n")
    (inputs / "snapshot_epochs.txt").write_text("\n".join(epoch + "Z" for epoch in epochs) + "\n")
    deck = replace_setting(deck, "TRAJ_FILE", inputs / "trajectory_geo.txt")
    deck = replace_setting(deck, "SNAPSHOT_LIST_FILE", inputs / "snapshot_epochs.txt")
    if args.dt is not None:
        deck = replace_setting(deck, "DT_TRACE", args.dt)
    direction_override = None
    if args.directions_per_energy is not None:
        # Record enough information to reconstruct both populations.  Merely saving
        # the rewritten deck would prove the final setting but conceal which source
        # quadrature it refined, weakening the comparison's provenance.
        deck, direction_override = override_directions_per_energy(
            deck, args.directions_per_energy,
            *(args.direction_grid or (0, 0)))
    if args.extension_passes is not None:
        deck = replace_setting(deck, "DS_UNRESOLVED_EXTENSION_PASSES", args.extension_passes)
    (inputs / "AMPS_PARAM.in").write_text(deck)
    # Reuse only the immutable topology cache, never physics output or field values.
    for cache in source.glob("amr.sig=*"):
        if cache.is_file() and not cache.is_symlink():
            shutil.copy2(cache, output / cache.name)
    command = launch_command(amps, args.mpi_ranks, args.threads, args.density_backend,
                             args.field_init, args.mover, args.mpi_binding)
    # Both audits are opt-in, read-only forensics. Full field fingerprints distinguish
    # a changed numerical field from changed integration on an identical field.
    environment = dict(os.environ, OMP_NUM_THREADS=str(args.threads),
                       AMPS_DENSITY_TRACE_DIAGNOSTICS="1", AMPS_MODE3D_FIELD_AUDIT="1")
    record = {"status": "RUNNING", "started_utc": utc_now(), "command": command,
              "source_deck_sha256": sha256_file(deck_path),
              "input_sha256": sha256_file(inputs / "AMPS_PARAM.in"),
              "executable_sha256": sha256_file(amps), "epochs": epochs,
              "mover_override": args.mover, "dt_override": args.dt,
              "direction_population_override": direction_override,
              "field_initialization": args.field_init, "density_backend": args.density_backend,
              "mpi_ranks": args.mpi_ranks, "threads": args.threads,
              "mpi_binding": args.mpi_binding,
              "extension_passes_override": args.extension_passes,
              "interpretation": "bounded numerical diagnosis, not observational validation"}
    # Hash referenced bytes as well as the deck: unchanged path text cannot prove
    # that the driver, spectrum or trajectory contents stayed fixed between runs.
    record["referenced_input_sha256"] = {
        str(Path(setting(deck, key)).resolve()): sha256_file(Path(setting(deck, key)))
        for key in ("DRIVER_FILE", "SPEC_TABLE_FILE", "TRAJ_FILE", "SNAPSHOT_LIST_FILE")}
    interrupted = False
    with (output / "amps.log").open("w") as log:
        process = subprocess.Popen(command, cwd=output, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        record.update(runner_pid=os.getpid(), launcher_pid=process.pid)
        atomic_json(output / "diagnostic_status.json", record)
        print(f"[START] diagnostic {output}; launcher={process.pid}", flush=True)
        try:
            code = process.wait()
        except KeyboardInterrupt:
            # The launcher owns a dedicated process group. Never signal by name:
            # another user's AMPS job may coexist on the same host.
            interrupted = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                code = process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                code = process.wait()
    complete = all(len(list(output.glob(pattern))) == len(epochs) for pattern in (
        "mode3d_points_density_snapshot_*.dat", "mode3d_points_spectrum_snapshot_*.dat",
        "mode3d_points_flux_snapshot_*.dat", "mode3d_termination_summary_snapshot_*.dat",
        "cutoff_3d_points_snapshot_*.dat", "standalone_mode3d_manifest_snapshot_*.json"))
    record.update(status="COMPLETE" if code == 0 and complete and not interrupted else "ERROR",
                  finished_utc=utc_now(), return_code=code, complete_products=complete,
                  interrupted=interrupted, log_sha256=sha256_file(output / "amps.log"),
                  output_sha256={str(path.resolve()): sha256_file(path)
                                 for pattern in ("*_snapshot_*.dat", "*_snapshot_*.json")
                                 for path in output.glob(pattern)})
    atomic_json(output / "diagnostic_status.json", record)
    print(f"[{record['status']}] {output}; exit={code}", flush=True)
    return 0 if record["status"] == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
