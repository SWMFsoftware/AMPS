#!/usr/bin/env python3
"""srcSEP3D point sampling: native AMPS acceptance (PSN01-PSN06).

Runs the built ./amps on a copy of the reduced-front injection smoke deck
srcSEP3D/examples/application-input/amps-injection-smoke.in (source.shock.model
= reduced-shock-front; an already-formed front at 50 R_sun along +z injects
H_PLUS from step 1, about 1000 macroparticles per step; observer publication
every step) extended to two AMPS sampling intervals (200 steps) and to a
'#begin sampling: points' block with block points and a point file, on one
and on four MPI ranks, and checks the products of the AMPS-synchronized point
sampler (srcSEP3D/amps/amps_point_sampler.h).  The four points lie in four
different populated cells near the front apex (selected from a 102-point
survey run; see srcSEP3D/test/README.md):

  PSN01.rankR  the run completes and installs the sampler (4 points, 1 species)
  PSN02.rankR  one product per AMPS output: the set of indices N of
               <prefix>.point-sampling.out=N.dat equals the set of the AMPS
               mesh files pic.H_PLUS.s=0.out=N.dat (same file indexing)
  PSN03.rankR  each product covers exactly the AMPS sampling length L (read
               from the AMPS log) and the times of interval N:
               first = N L dt, last = first + (L-1) dt
  PSN04.rankR  the product is non-trivial (samples > 0) and the cumulative
               fluence of every row equals the running sum of the interval
               fluences of the earlier products
  PSN05.rankR  independent density check (each rank count): the interval-mean total
               density at each point (all channels plus the below/above-range
               tails, plasma frame) equals the mean, over the same L particle
               states, of the number density that srcSEP3D's separate
               observer publication (output/sampling.cpp, cells.csv, written
               every step) reports for the same cell, identified by the
               position-hash cell id of the cell centre in the product
  PSN06        one rank and four ranks give the same files, the same sample
               counts and the same values within a relative 1e-10 (MPI_SUM
               may change the last bits)

Each line '[ID] STATUS message' is read by tools/sep_test_orchestrator.py
(read_bracket_text).  The exit status is 1 when any check fails or errors.

Usage (from the AMPS root, after the native srcSEP3D build):
  python3 srcSEP3D/test/run_point_sampling_native.py --amps ./amps \\
      --output-dir test_output/point-sampling-native

After a failure of only one leg, the retained run of the other leg can be
re-checked instead of relaunched, e.g. '--reuse-ranks 1' (the leg's exit
status is kept in rank<R>/exit_status; the reuse is refused when ./amps is
newer than the retained run.log).

Each rank count runs in its own working directory, because AMPS writes its
own output (and the point-sampling products) to the compiled
OutputDirectory (PT/plots) relative to the working directory and clears it at
start-up.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]           # AMPS root
SMOKE = ROOT / "srcSEP3D/examples/application-input"
DECK = SMOKE / "amps-injection-smoke.in"
PREFIX = "sep3d"                                      # output.prefix of sep3d.in

SAMPLING_BLOCK = """
! ---- point sampling under test (run_point_sampling_native.py) ----
#begin sampling: points
  products = density, flux, spectrum, fluence
  point_file = points.txt
  #begin species: H_PLUS
    minimum_energy_mev = 0.01
    maximum_energy_mev = 100
    energy_intervals = 8
    energy_scale = log
  #end species
  ! Near the apex of the injecting front (50 R_sun = 3.48e10 m along +z).
  #begin point: apex
    position_x_m = 0
    position_y_m = 0
    position_z_m = 3.55e10
  #end point
  #begin point: flank_x
    position_x_m = 3.544086e9
    position_y_m = 0
    position_z_m = 3.532265e10
  #end point
#end sampling
"""
POINT_FILE = "! x y z [name] in metres (HCI)\n0 3.544086e9 3.532265e10 flank_y\n0 0 3.8e10\n"


def emit(test_id: str, status: str, message: str, results: list) -> None:
    print(f"[{test_id}] {status} {message}", flush=True)
    results.append(status)


def write_deck(work: Path, steps: int) -> Path:
    """Copy the smoke deck with absolute #include paths, the step budget and
    the sampling block appended."""
    text = DECK.read_text()
    text, budget = re.subn(r"maximum_time_steps = \d+",
                           f"maximum_time_steps = {steps}", text, count=1)
    if budget != 1:
        raise SystemExit(f"{DECK}: step budget anchor not found")
    text = re.sub(r'#include "([^"/]+)"',
                  lambda m: f'#include "{(SMOKE / m.group(1)).as_posix()}"', text)
    deck = work / "deck.in"
    deck.write_text(text + SAMPLING_BLOCK)
    (work / "points.txt").write_text(POINT_FILE)
    return deck


def parse_product(path: Path) -> dict:
    """Parse the Tecplot point-sampling product into
    {'aux': {...}, 'variables': [...], 'zones': [{'title','aux','rows'}]}."""
    product = {"aux": {}, "variables": [], "zones": []}
    zone = None
    for line in path.read_text().splitlines():
        if line.startswith("VARIABLES="):
            product["variables"] = re.findall(r'"([^"]*)"', line)
        elif line.startswith("DATASETAUXDATA "):
            key, value = re.match(r'DATASETAUXDATA (\w+)="([^"]*)"', line).groups()
            product["aux"][key] = value
        elif line.startswith("ZONE "):
            zone = {"title": re.search(r'T="([^"]*)"', line).group(1),
                    "aux": {}, "rows": []}
            product["zones"].append(zone)
        elif line.startswith("AUXDATA "):
            key, value = re.match(r'AUXDATA (\w+)="([^"]*)"', line).groups()
            zone["aux"][key] = value
        elif zone is not None and line.strip():
            zone["rows"].append([float(v) for v in line.split()])
    return product


def column(product: dict, name: str) -> int:
    for index, variable in enumerate(product["variables"]):
        if variable == name or variable.startswith(name + " "):
            return index
    raise KeyError(name)


def stable_cell_id(x: float, y: float, z: float) -> int:
    """srcSEP3D StableCellId (main_lib.cpp): FNV-1a 64 over the little-endian
    bytes of the three centre doubles."""
    value = 14695981039346656037
    for coordinate in (x, y, z):
        for byte in struct.pack("<d", coordinate):
            value ^= byte
            value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value or 1


def run_rank(args, work: Path, deck: Path, ranks: int, results: list) -> dict | None:
    """Launch (or, with --reuse-ranks, re-read) the run on `ranks` MPI ranks.

    A reused run is checked exactly like a fresh one, from its retained
    run.log, exit status file and PT/plots products; it is the caller's
    responsibility that the retained run was made by the same ./amps binary
    and deck (the runner refuses a reuse when the binary is newer than the
    retained log, see main()).
    """
    rundir = work / f"rank{ranks}"
    log = rundir / "run.log"
    status_file = rundir / "exit_status"
    if ranks in args.reuse_ranks:
        if not log.is_file() or not status_file.is_file():
            emit(f"PSN01.rank{ranks}", "ERROR",
                 f"--reuse-ranks {ranks}: no retained run.log/exit_status in {rundir}",
                 results)
            return None
        code = int(status_file.read_text().strip())
    else:
        shutil.rmtree(rundir, ignore_errors=True)
        rundir.mkdir(parents=True)
        command = [args.mpiexec, "-n", str(ranks), str(Path(args.amps).resolve()),
                   "-input", str(deck), "--output-dir", str(rundir / "publications")]
        try:
            with log.open("w") as stream:
                code = subprocess.run(command, cwd=rundir, stdout=stream,
                                      stderr=subprocess.STDOUT,
                                      timeout=args.timeout).returncode
        except subprocess.TimeoutExpired:
            emit(f"PSN01.rank{ranks}", "ERROR", f"timeout after {args.timeout} s; log {log}",
                 results)
            return None
        # Retained so that a later invocation can reuse this run (--reuse-ranks).
        status_file.write_text(f"{code}\n")
    text = log.read_text(errors="replace")
    if code != 0 or "[srcSEP3D] point sampling: 4 point(s), 1 species" not in text:
        emit(f"PSN01.rank{ranks}", "FAIL",
             f"exit {code} or sampler not installed; log {log}", results)
        return None
    emit(f"PSN01.rank{ranks}", "PASS", f"run completed; log {log}", results)
    length = re.search(r"current sample length:(\d+)", text)
    dt = re.search(r"dt_s=([0-9.eE+-]+)", text)
    if not length or not dt:
        emit(f"PSN02.rank{ranks}", "ERROR", "sample length or dt not in the log", results)
        return None
    plots = rundir / "PT/plots"
    products = {int(m.group(1)): p for p in plots.glob(f"{PREFIX}.point-sampling.out=*.dat")
                if (m := re.search(r"out=(\d+)\.dat$", p.name))}
    meshes = {int(m.group(1)) for p in plots.glob("pic.H_PLUS.s=0.out=*.dat")
              if (m := re.search(r"out=(\d+)\.dat$", p.name))}
    return {"rundir": rundir, "length": int(length.group(1)),
            "dt": float(dt.group(1)), "files": products, "meshes": meshes,
            "parsed": {n: parse_product(p) for n, p in sorted(products.items())}}


def check_rank(run: dict, ranks: int, steps: int, results: list) -> None:
    tag = f"rank{ranks}"
    expected = set(range(steps // run["length"]))
    if not run["files"] or set(run["files"]) != run["meshes"] or set(run["files"]) != expected:
        emit(f"PSN02.{tag}", "FAIL",
             f"product indices {sorted(run['files'])} vs AMPS mesh indices "
             f"{sorted(run['meshes'])}, expected {sorted(expected)}", results)
    else:
        emit(f"PSN02.{tag}", "PASS",
             f"products out={sorted(run['files'])} match the AMPS mesh files", results)

    length, dt = run["length"], run["dt"]
    problems = []
    for n, product in run["parsed"].items():
        aux = product["aux"]
        first, last = float(aux["first_sample_time_s"]), float(aux["last_sample_time_s"])
        if int(aux["sample_calls"]) != length or int(aux["amps_output_file_number"]) != n:
            problems.append(f"out={n}: sample_calls {aux['sample_calls']} != {length}")
        # The product writes times with 11 significant digits, so compare at
        # a relative 1e-9 (the interval arithmetic itself is exact to ~1e-15).
        if not math.isclose(first, n * length * dt, rel_tol=1e-9, abs_tol=1e-9 * dt) or \
           not math.isclose(last - first, (length - 1) * dt, rel_tol=1e-9):
            problems.append(f"out={n}: times {first}..{last} do not span the interval")
    emit(f"PSN03.{tag}", "FAIL" if problems else "PASS",
         "; ".join(problems) or f"every product covers L={length} calls of dt={dt} s",
         results)

    samples, mismatches = 0.0, []
    running: dict = {}
    for n, product in run["parsed"].items():
        for frame in ("plasma", "inertial"):
            s_col = column(product, f"samples_{frame}")
            fi = column(product, f"fluence_interval_{frame}")
            fc = column(product, f"fluence_{frame}")
            for zone in product["zones"]:
                for k, row in enumerate(zone["rows"]):
                    samples += row[s_col]
                    key = (zone["title"], frame, k)
                    running[key] = running.get(key, 0.0) + row[fi]
                    if not math.isclose(row[fc], running[key], rel_tol=1e-9, abs_tol=1e-300):
                        mismatches.append(f"out={n} {zone['title']} {frame} row {k}")
    if samples <= 0 or mismatches:
        emit(f"PSN04.{tag}", "FAIL",
             f"samples={samples:g}; cumulative-fluence mismatches: {mismatches[:3]}", results)
    else:
        emit(f"PSN04.{tag}", "PASS",
             f"{samples:g} macroparticle samples; cumulative fluence = running "
             "sum of interval fluences", results)


def check_density_against_publication(run: dict, ranks: int, results: list) -> None:
    """PSN05: compare with the cells.csv number density of the same cell.
    On four ranks this also guards the observer publication across the AMPS
    block rebalancing that follows each AMPS output (a stale cell index once
    dropped every migrated cell from cells.csv after the first output)."""
    snapshots = {}   # simulation time -> {cell_id: density} for species 0
    for manifest in sorted((run["rundir"] / "publications").glob(f"{PREFIX}.*/manifest.txt")):
        time = float(re.search(r"simulation_time_s=(\S+)", manifest.read_text()).group(1))
        densities = {}
        for line in (manifest.parent / "cells.csv").read_text().splitlines()[1:]:
            fields = line.split(",")
            if fields[1] == "0":
                densities[int(fields[0])] = float(fields[3])
        snapshots[time] = densities
    if not snapshots:
        emit(f"PSN05.rank{ranks}", "ERROR", "no cells.csv publications found", results)
        return
    times = sorted(snapshots)
    length, dt = run["length"], run["dt"]
    problems, compared, nonzero = [], 0, 0
    for n, product in run["parsed"].items():
        first = float(product["aux"]["first_sample_time_s"])
        n_col = column(product, "n_plasma")
        for zone in product["zones"]:
            aux = zone["aux"]
            cell = stable_cell_id(float(aux["cell_center_x_m"]),
                                  float(aux["cell_center_y_m"]),
                                  float(aux["cell_center_z_m"]))
            sampled = 1e6 * (sum(row[n_col] for row in zone["rows"]) +
                             float(aux["below_range_density_plasma_cm3"]) +
                             float(aux["above_range_density_plasma_cm3"]))
            total = 0.0
            for c in range(length):
                t = first + c * dt
                match = [s for s in times if abs(s - t) <= 1e-6 * dt]
                if match:
                    total += snapshots[match[0]].get(cell, 0.0)
                elif abs(t) > 1e-6 * dt:   # only the initial state precedes publication
                    problems.append(f"no publication at t={t}")
            expected = total / length
            compared += 1
            nonzero += expected > 0.0
            if not math.isclose(sampled, expected, rel_tol=1e-10, abs_tol=1e-300):
                problems.append(f"out={n} {zone['title']}: sampled {sampled:.12e} m^-3, "
                                f"publication mean {expected:.12e} m^-3")
    if problems or nonzero == 0:
        emit(f"PSN05.rank{ranks}", "FAIL",
             "; ".join(problems[:4]) or "every compared density is zero", results)
    else:
        emit(f"PSN05.rank{ranks}", "PASS", f"{compared} interval densities ({nonzero} "
             "non-zero) equal the cells.csv mean over the same particle states", results)


def check_rank_agreement(a: dict, b: dict, results: list) -> None:
    if set(a["parsed"]) != set(b["parsed"]):
        emit("PSN06", "FAIL", "different product files on one and four ranks", results)
        return
    worst, problems = 0.0, []
    for n in a["parsed"]:
        pa, pb = a["parsed"][n], b["parsed"][n]
        sample_cols = {i for i, v in enumerate(pa["variables"]) if v.startswith("samples_")}
        for za, zb in zip(pa["zones"], pb["zones"]):
            for ra, rb in zip(za["rows"], zb["rows"]):
                for i, (x, y) in enumerate(zip(ra, rb)):
                    if i in sample_cols and x != y:
                        problems.append(f"out={n} {za['title']}: samples {x} vs {y}")
                    scale = max(abs(x), abs(y))
                    if scale > 0:
                        worst = max(worst, abs(x - y) / scale)
    if problems or worst > 1e-10:
        emit("PSN06", "FAIL", f"max relative difference {worst:.3e}; {problems[:3]}", results)
    else:
        emit("PSN06", "PASS", f"identical samples; max relative difference {worst:.3e}",
             results)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--amps", default=str(ROOT / "amps"))
    parser.add_argument("--output-dir", default=str(ROOT / "test_output/point-sampling-native"))
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--steps", type=int, default=200,
                        help="time steps per run (two AMPS intervals of 100)")
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--reuse-ranks", type=int, nargs="*", default=[],
                        choices=(1, 4),
                        help="re-read the retained run of these rank counts "
                             "instead of relaunching it (e.g. after a failure "
                             "of only the other leg); refused when ./amps is "
                             "newer than the retained run.log")
    args = parser.parse_args()
    if not Path(args.amps).is_file():
        print(f"[PSN01.rank1] ERROR amps executable {args.amps} not found")
        return 1
    for ranks in args.reuse_ranks:
        log = Path(args.output_dir).resolve() / f"rank{ranks}" / "run.log"
        if log.is_file() and Path(args.amps).stat().st_mtime > log.stat().st_mtime:
            print(f"[PSN01.rank{ranks}] ERROR --reuse-ranks {ranks}: {args.amps} "
                  f"is newer than the retained {log}; rerun that leg")
            return 1
    work = Path(args.output_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    deck = write_deck(work, args.steps)
    results: list = []
    runs = {}
    for ranks in (1, 4):
        run = run_rank(args, work, deck, ranks, results)
        if run is None:
            continue
        check_rank(run, ranks, args.steps, results)
        runs[ranks] = run
    for ranks in (1, 4):
        if ranks in runs:
            check_density_against_publication(runs[ranks], ranks, results)
        else:
            emit(f"PSN05.rank{ranks}", "ERROR", "the run is unavailable", results)
    if 1 in runs and 4 in runs:
        check_rank_agreement(runs[1], runs[4], results)
    else:
        emit("PSN06", "ERROR", "a run is unavailable", results)
    bad = sum(status in ("FAIL", "ERROR") for status in results)
    print(f"point-sampling native summary: {len(results) - bad} passed, {bad} failed/errored")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
