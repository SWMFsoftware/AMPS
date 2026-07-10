#!/usr/bin/env python3
"""
C17 — Dipole charge-sign and velocity-reversal symmetry.

Run from the directory containing the AMPS executable:

    python srcEarth/test/C17/run_C17.py -np 4 -nt 16

C17 is an exact internal symmetry test.  In a static magnetic field with E=0,
Lorentz trajectories obey the time-reversal relation

    T_q(x0, v, R) = T_-q(x0, -v, R).

The current gridless cutoff output exposes directional cutoff maps, not a full
T(R,Omega) table.  Therefore this runner checks the equivalent directional-map
relation

    Rc_q(lon,lat) == Rc_-q(lon+180, -lat)

for every non-polar sky-map cell at every selected observation point.  It also
constructs a derived step-transmission proxy at a configurable rigidity list.
"""

from __future__ import print_function

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

TEST_ID = "C17"
TEST_NAME = "Dipole charge-sign and velocity-reversal symmetry"

RE_KM = 6371.2
DEFAULT_RIGIDITIES_GV = (0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)


class DirCell(object):
    def __init__(self, lon_deg, lat_deg, rc_gv, emin_mev):
        self.lon_deg = lon_deg
        self.lat_deg = lat_deg
        self.rc_gv = rc_gv
        self.emin_mev = emin_mev


def parse_float_list(text, what):
    vals = []
    for item in str(text).replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        try:
            vals.append(float(item))
        except ValueError:
            raise SystemExit("Invalid %s value: %r" % (what, item))
    if not vals:
        raise SystemExit("No values supplied for %s" % what)
    return vals


def parse_mover_list(text):
    allowed = {"BORIS", "RK2", "RK4", "RK6", "GC2", "GC4", "GC6"}
    movers = []
    for item in str(text).replace(";", ",").split(","):
        item = item.strip().upper()
        if not item:
            continue
        if item not in allowed:
            raise SystemExit("Unsupported mover %r; allowed values: %s" % (item, ",".join(sorted(allowed))))
        movers.append(item)
    if not movers:
        raise SystemExit("No movers requested")
    return movers


def bool_token(text):
    t = str(text).strip().upper()
    if t in ("T", "TRUE", "1", "YES", "Y"):
        return "T"
    if t in ("F", "FALSE", "0", "NO", "N"):
        return "F"
    raise SystemExit("Expected boolean token T/F, got %r" % text)


def norm_lon(lon):
    x = lon % 360.0
    if abs(x - 360.0) < 1.0e-8 or abs(x) < 1.0e-8:
        return 0.0
    return x


def norm_lat(lat):
    if abs(lat) < 1.0e-8:
        return 0.0
    return lat


def key_lon_lat(lon, lat):
    return (round(norm_lon(lon), 6), round(norm_lat(lat), 6))


def reversed_direction_key(lon, lat):
    return key_lon_lat(lon + 180.0, -lat)


def point_xyz_km(lon_deg, lat_deg, alt_km):
    r = RE_KM + alt_km
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    c = math.cos(lat)
    return (r * c * math.cos(lon), r * c * math.sin(lon), r * math.sin(lat))


def build_points(obs_lons, obs_lats, alt_km):
    points = []
    point_id = 0
    for lon in obs_lons:
        for lat in obs_lats:
            x, y, z = point_xyz_km(lon, lat, alt_km)
            points.append({
                "point_id": point_id,
                "obs_lon_deg": lon,
                "obs_lat_deg": lat,
                "obs_alt_km": alt_km,
                "x_km": x,
                "y_km": y,
                "z_km": z,
            })
            point_id += 1
    return points


def points_block(points):
    lines = []
    for p in points:
        lines.append("POINT %.12e %.12e %.12e ! km  point_id=%d lon=%g lat=%g alt=%g" %
                     (p["x_km"], p["y_km"], p["z_km"], p["point_id"],
                      p["obs_lon_deg"], p["obs_lat_deg"], p["obs_alt_km"]))
    return "\n".join(lines)


def render_input_template(template_path, output_path, run_id, species_name, charge_e,
                          mass_amu, points, args):
    text = template_path.read_text()
    repl = {
        "__RUN_ID__": run_id,
        "__CUTOFF_EMIN_MEV__": "%.10g" % args.cutoff_emin,
        "__CUTOFF_EMAX_MEV__": "%.10g" % args.cutoff_emax,
        "__CUTOFF_MAX_PARTICLES__": str(args.cutoff_max_particles),
        "__CUTOFF_NENERGY__": str(args.cutoff_nenergy),
        "__CUTOFF_UPPER_SCAN_N__": str(args.cutoff_upper_scan_n),
        "__DIRMAP_LON_RES__": "%.10g" % args.dir_lon_res,
        "__DIRMAP_LAT_RES__": "%.10g" % args.dir_lat_res,
        "__SPECIES_NAME__": species_name,
        "__CHARGE__": str(charge_e),
        "__MASS_AMU__": "%.10g" % mass_amu,
        "__POINTS_BLOCK__": points_block(points),
        "__DT_TRACE__": "%.10g" % args.dt_trace,
        "__ADAPTIVE_DT__": bool_token(args.adaptive_dt),
        "__MAX_TRACE_TIME__": "%.10g" % args.max_trace_time,
        "__MAX_TRACE_DISTANCE_RE__": "%.10g" % args.max_trace_distance,
        "__GRIDLESS_MPI_SCHEDULER__": args.scheduler,
        "__GRIDLESS_MPI_DYNAMIC_CHUNK__": str(args.dynamic_chunk),
        "__GRIDLESS_THREADS__": str(args.nt),
    }
    for key, val in repl.items():
        text = text.replace(key, val)
    output_path.write_text(text)


def run_command(cmd, cwd, log_path):
    with log_path.open("w") as log:
        log.write("Command:\n  " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
        return proc.wait()


def parse_variables(line):
    return [n.strip().lower() for n in re.findall(r'"([^"]+)"', line)]


def pick_col(variables, candidates, fallback):
    if variables:
        norm = [v.lower().replace(" ", "_").replace("-", "_") for v in variables]
        for cand in candidates:
            c = cand.lower().replace(" ", "_").replace("-", "_")
            for i, v in enumerate(norm):
                if v == c or c in v:
                    return i
    return fallback


def parse_directional_map(path):
    variables = None
    rows = []
    in_zone = False
    with path.open("r", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            upper = line.upper()
            if upper.startswith("TITLE"):
                continue
            if upper.startswith("VARIABLES"):
                variables = parse_variables(line)
                continue
            if upper.startswith("ZONE"):
                in_zone = True
                continue
            if not in_zone:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                lon_col = pick_col(variables, ["lon_deg", "longitude", "lon"], 0)
                lat_col = pick_col(variables, ["lat_deg", "latitude", "lat"], 1)
                rc_col = pick_col(variables, ["Rc_GV", "cutoff_gv", "rigidity_gv"], 2)
                emin_col = pick_col(variables, ["Emin_MeV", "energy_mev"], 3 if len(parts) > 3 else -1)
                lon = float(parts[lon_col])
                lat = float(parts[lat_col])
                rc = float(parts[rc_col])
                emin = float(parts[emin_col]) if 0 <= emin_col < len(parts) else float("nan")
                rows.append(DirCell(lon, lat, rc, emin))
            except (ValueError, IndexError):
                continue
    if not rows:
        raise RuntimeError("No directional-map rows parsed from %s" % path)
    return rows


def map_to_dict(rows):
    d = {}
    for c in rows:
        d[key_lon_lat(c.lon_deg, c.lat_deg)] = c
    return d


def find_dir_map_file(workdir, point_id):
    candidates = [
        workdir / ("cutoff_gridless_dir_map_point_%04d.dat" % point_id),
        workdir / ("cutoff_3d_dir_map_point_%04d.dat" % point_id),
    ]
    for p in candidates:
        if p.exists():
            return p
    matches = sorted(workdir.glob("*dir*map*point*%04d*.dat" % point_id))
    if matches:
        return matches[0]
    return candidates[0]


def finite_positive(x):
    return math.isfinite(x) and x > 0.0


def step_transmission_from_cutoff(rigidity_gv, rc_gv):
    if not finite_positive(rc_gv):
        return None
    return 1 if rigidity_gv >= rc_gv else 0


def compare_charge_reversal(points, pos_dir, neg_dir, rigidities_gv, rel_tol, abs_tol_gv,
                            ignore_poles=True):
    summary = []
    pair_rows = []
    messages = []
    overall_passed = True

    for p in points:
        point_id = p["point_id"]
        pos_file = find_dir_map_file(pos_dir, point_id)
        neg_file = find_dir_map_file(neg_dir, point_id)
        if not pos_file.exists() or not neg_file.exists():
            overall_passed = False
            messages.append("Missing directional-map file for point %d: %s or %s" % (point_id, pos_file, neg_file))
            continue

        pos_rows = parse_directional_map(pos_file)
        neg_rows = parse_directional_map(neg_file)
        neg_map = map_to_dict(neg_rows)

        n_pairs = 0
        n_missing_reverse = 0
        n_bad_values = 0
        n_failed_rc = 0
        n_t_mismatch = 0
        max_abs_diff = 0.0
        max_rel_diff = 0.0
        worst = None

        for c in pos_rows:
            if ignore_poles and abs(abs(c.lat_deg) - 90.0) < 1.0e-8:
                continue
            rev_key = reversed_direction_key(c.lon_deg, c.lat_deg)
            rev = neg_map.get(rev_key)
            if rev is None:
                n_missing_reverse += 1
                pair_rows.append({
                    "point_id": point_id,
                    "obs_lon_deg": p["obs_lon_deg"],
                    "obs_lat_deg": p["obs_lat_deg"],
                    "dir_lon_deg": c.lon_deg,
                    "dir_lat_deg": c.lat_deg,
                    "rev_dir_lon_deg": rev_key[0],
                    "rev_dir_lat_deg": rev_key[1],
                    "Rc_plus_GV": c.rc_gv,
                    "Rc_minus_reversed_GV": "",
                    "abs_diff_GV": "",
                    "rel_diff": "",
                    "passed_rc": False,
                    "transmission_mismatches": "",
                    "note": "missing reversed direction in negative-charge map",
                })
                continue

            n_pairs += 1
            if not finite_positive(c.rc_gv) or not finite_positive(rev.rc_gv):
                n_bad_values += 1
                passed_rc = False
                abs_diff = float("nan")
                rel_diff = float("nan")
            else:
                abs_diff = abs(c.rc_gv - rev.rc_gv)
                denom = max(abs(c.rc_gv), abs(rev.rc_gv), 1.0e-30)
                rel_diff = abs_diff / denom
                passed_rc = (abs_diff <= abs_tol_gv) or (rel_diff <= rel_tol)
                if abs_diff > max_abs_diff:
                    max_abs_diff = abs_diff
                if rel_diff > max_rel_diff:
                    max_rel_diff = rel_diff
                if worst is None or rel_diff > worst["rel_diff"]:
                    worst = {
                        "dir_lon_deg": c.lon_deg,
                        "dir_lat_deg": c.lat_deg,
                        "rev_dir_lon_deg": rev.lon_deg,
                        "rev_dir_lat_deg": rev.lat_deg,
                        "Rc_plus_GV": c.rc_gv,
                        "Rc_minus_reversed_GV": rev.rc_gv,
                        "abs_diff_GV": abs_diff,
                        "rel_diff": rel_diff,
                    }
            if not passed_rc:
                n_failed_rc += 1

            tm = 0
            for R in rigidities_gv:
                tp = step_transmission_from_cutoff(R, c.rc_gv)
                tn = step_transmission_from_cutoff(R, rev.rc_gv)
                if tp is None or tn is None or tp != tn:
                    tm += 1
            n_t_mismatch += tm

            pair_rows.append({
                "point_id": point_id,
                "obs_lon_deg": p["obs_lon_deg"],
                "obs_lat_deg": p["obs_lat_deg"],
                "dir_lon_deg": c.lon_deg,
                "dir_lat_deg": c.lat_deg,
                "rev_dir_lon_deg": rev.lon_deg,
                "rev_dir_lat_deg": rev.lat_deg,
                "Rc_plus_GV": c.rc_gv,
                "Rc_minus_reversed_GV": rev.rc_gv,
                "abs_diff_GV": abs_diff,
                "rel_diff": rel_diff,
                "passed_rc": passed_rc,
                "transmission_mismatches": tm,
                "note": "",
            })

        point_passed = (n_missing_reverse == 0 and n_bad_values == 0 and
                        n_failed_rc == 0 and n_t_mismatch == 0)
        if not point_passed:
            overall_passed = False
            messages.append(
                "point %d lon=%g lat=%g: missing=%d bad=%d failed_rc=%d T_mismatch=%d max_rel=%.3e" %
                (point_id, p["obs_lon_deg"], p["obs_lat_deg"], n_missing_reverse,
                 n_bad_values, n_failed_rc, n_t_mismatch, max_rel_diff)
            )

        item = {
            "point_id": point_id,
            "obs_lon_deg": p["obs_lon_deg"],
            "obs_lat_deg": p["obs_lat_deg"],
            "obs_alt_km": p["obs_alt_km"],
            "n_pairs_checked": n_pairs,
            "n_missing_reverse": n_missing_reverse,
            "n_bad_values": n_bad_values,
            "n_failed_rc": n_failed_rc,
            "n_transmission_mismatches": n_t_mismatch,
            "max_abs_diff_GV": max_abs_diff,
            "max_rel_diff": max_rel_diff,
            "rel_tol": rel_tol,
            "abs_tol_GV": abs_tol_gv,
            "passed": point_passed,
        }
        if worst:
            item.update({"worst_" + k: v for k, v in worst.items()})
        summary.append(item)

    return summary, pair_rows, overall_passed, messages


def write_dict_csv(rows, path):
    if not rows:
        path.write_text("")
        return
    keys = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_reference_csv(path, points, args, rigidities):
    dir_lons = []
    nlon = int(round(360.0 / args.dir_lon_res))
    for i in range(nlon):
        dir_lons.append(i * args.dir_lon_res)
    nlat = int(math.floor(180.0 / args.dir_lat_res)) + 1
    dir_lats = []
    for j in range(nlat):
        lat = -90.0 + args.dir_lat_res * j
        if lat > 90.0:
            lat = 90.0
        dir_lats.append(lat)

    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "point_id", "obs_lon_deg", "obs_lat_deg", "obs_alt_km",
            "dir_lon_deg", "dir_lat_deg", "reversed_dir_lon_deg", "reversed_dir_lat_deg",
            "expected", "rel_tol", "abs_tol_GV", "rigidity_list_GV", "note",
        ])
        for p in points:
            for lat in dir_lats:
                for lon in dir_lons:
                    if args.ignore_poles and abs(abs(lat) - 90.0) < 1.0e-8:
                        continue
                    rev = reversed_direction_key(lon, lat)
                    w.writerow([
                        p["point_id"], p["obs_lon_deg"], p["obs_lat_deg"], p["obs_alt_km"],
                        lon, lat, rev[0], rev[1],
                        "Rc_plus(dir) == Rc_minus(reversed_dir)",
                        "%.3e" % args.rel_tol, "%.3e" % args.abs_tol_gv,
                        ";".join("%.10g" % r for r in rigidities),
                        "Lorentz-force time-reversal symmetry in static dipole, E=0",
                    ])


def make_plot(pair_rows, out_path):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    vals = []
    for r in pair_rows:
        try:
            vals.append(float(r["rel_diff"]))
        except Exception:
            pass
    if not vals:
        return
    plt.figure(figsize=(7.0, 4.5))
    plt.hist(vals, bins=50)
    plt.yscale("log")
    plt.xlabel("relative difference |Rc+ - Rc-(-v)| / max(Rc)")
    plt.ylabel("number of directional cells")
    plt.title("C17 charge/time-reversal symmetry residuals")
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="C17 charge-sign / velocity-reversal symmetry test",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples:
          python srcEarth/test/C17/run_C17.py -np 4 -nt 16
          python srcEarth/test/C17/run_C17.py --movers BORIS,RK4 --dir-lon-res 60 --dir-lat-res 30
          python srcEarth/test/C17/run_C17.py --lats=-60,-30,0,30,60 --dry-run
          python srcEarth/test/C17/run_C17.py --skip-run --workdir test_output/C17_gridless
        """),
    )
    parser.add_argument("-np", type=int, default=4, help="MPI ranks passed to mpirun; default: 4")
    parser.add_argument("-nt", type=int, default=16, help="thread count passed to AMPS; default: 16")
    parser.add_argument("--movers", default="BORIS", help="comma-separated mover list; default: BORIS")
    parser.add_argument("--lats", default="-60,-30,0,30,60", help="observation latitudes, deg; default: -60,-30,0,30,60. Use --lats=-60,... if needed")
    parser.add_argument("--lons", default="0", help="observation longitudes, deg; default: 0")
    parser.add_argument("--alt", type=float, default=9000.0, help="observation altitude, km; default: 9000")
    parser.add_argument("--dir-lon-res", type=float, default=30.0, help="directional map longitude resolution, deg; default: 30")
    parser.add_argument("--dir-lat-res", type=float, default=30.0, help="directional map latitude resolution, deg; default: 30")
    parser.add_argument("--rigidities", default=",".join("%.10g" % r for r in DEFAULT_RIGIDITIES_GV), help="rigidity list for derived step-transmission check, GV")
    parser.add_argument("--cutoff-emin", type=float, default=0.1, help="CUTOFF_EMIN, MeV/n; default: 0.1")
    parser.add_argument("--cutoff-emax", type=float, default=20000.0, help="CUTOFF_EMAX, MeV/n; default: 20000")
    parser.add_argument("--cutoff-nenergy", type=int, default=50, help="CUTOFF_NENERGY; default: 50")
    parser.add_argument("--cutoff-upper-scan-n", type=int, default=80, help="CUTOFF_UPPER_SCAN_N; default: 80")
    parser.add_argument("--cutoff-max-particles", type=int, default=64, help="CUTOFF_MAX_PARTICLES for ordinary isotropic point cutoff; directional map is separate; default: 64")
    parser.add_argument("--dt-trace", type=float, default=1.0, help="DT_TRACE, seconds; default: 1")
    parser.add_argument("--adaptive-dt", default="T", help="ADAPTIVE_DT setting T/F; default: T")
    parser.add_argument("--max-trace-time", type=float, default=600.0, help="maximum trace time, seconds; default: 600")
    parser.add_argument("--max-trace-distance", type=float, default=400.0, help="maximum cumulative trace distance, Re; default: 400")
    parser.add_argument("--rel-tol", type=float, default=1.0e-6, help="relative tolerance for Rc symmetry; default: 1e-6")
    parser.add_argument("--abs-tol-gv", type=float, default=1.0e-8, help="absolute tolerance for Rc symmetry in GV; default: 1e-8")
    parser.add_argument("--include-poles", dest="ignore_poles", action="store_false", help="include +/-90 deg sky-map cells; default skips poles because longitude is degenerate there")
    parser.set_defaults(ignore_poles=True)
    parser.add_argument("--scheduler", default="DYNAMIC", choices=["DYNAMIC", "BLOCK_CYCLIC", "STATIC"], help="gridless MPI scheduler; default: DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0, help="gridless dynamic chunk; 0 means AMPS auto heuristic; default: 0")
    parser.add_argument("--amps", default="./amps", help="AMPS executable path relative to launch directory; default: ./amps")
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher executable; default: mpirun")
    parser.add_argument("--workdir", default="test_output/C17_gridless", help="base work directory; default: test_output/C17_gridless")
    parser.add_argument("--skip-run", action="store_true", help="analyze existing output without launching AMPS")
    parser.add_argument("--keep", action="store_true", help="keep existing work directory instead of deleting it")
    parser.add_argument("--dry-run", action="store_true", help="render inputs and print commands, but do not run or analyze AMPS output")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.np < 1:
        raise SystemExit("-np must be >= 1")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    if args.dir_lon_res <= 0.0 or args.dir_lat_res <= 0.0:
        raise SystemExit("directional-map resolutions must be positive")
    nlon_exact = 360.0 / args.dir_lon_res
    if abs(nlon_exact - round(nlon_exact)) > 1.0e-8:
        raise SystemExit("--dir-lon-res must divide 360 exactly for reversal mapping")
    nlat_exact = 180.0 / args.dir_lat_res
    if abs(nlat_exact - round(nlat_exact)) > 1.0e-8:
        raise SystemExit("--dir-lat-res must divide 180 exactly for reversal mapping")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")
    if args.alt <= 0.0:
        raise SystemExit("--alt must be positive")

    movers = parse_mover_list(args.movers)
    obs_lats = parse_float_list(args.lats, "--lats")
    obs_lons = parse_float_list(args.lons, "--lons")
    rigidities = parse_float_list(args.rigidities, "--rigidities")
    points = build_points(obs_lons, obs_lats, args.alt)

    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    template_input = script_dir / "AMPS_PARAM_C17_gridless.in"
    reference_csv_src = script_dir / "reference_C17_symmetry.csv"
    base_workdir = (launch_dir / args.workdir).resolve()

    if not args.skip_run:
        if base_workdir.exists() and not args.keep:
            shutil.rmtree(base_workdir)
        base_workdir.mkdir(parents=True, exist_ok=True)
    else:
        if not base_workdir.exists():
            raise SystemExit("--skip-run requested but workdir does not exist: %s" % base_workdir)

    write_reference_csv(reference_csv_src, points, args, rigidities)
    if not args.skip_run:
        shutil.copy2(reference_csv_src, base_workdir / "reference_C17_symmetry.csv")

    amps_path = Path(args.amps)
    if not amps_path.is_absolute():
        amps_path = (launch_dir / amps_path).resolve()

    run_records = []
    all_summary = []
    all_pair_rows = []
    all_messages = []
    overall_passed = True

    for mover in movers:
        mover_dir = base_workdir / mover.lower()
        pos_dir = mover_dir / "charge_plus"
        neg_dir = mover_dir / "charge_minus_reversed"

        if not args.skip_run:
            pos_dir.mkdir(parents=True, exist_ok=True)
            neg_dir.mkdir(parents=True, exist_ok=True)
            render_input_template(template_input, pos_dir / "AMPS_PARAM_C17.in",
                                  "C17_%s_charge_plus" % mover.lower(),
                                  "PROTON", +1, 1.0073, points, args)
            # Use the proton mass with charge -1.  This isolates charge-sign symmetry
            # from species mass / kinetic-energy-to-rigidity conversion differences.
            render_input_template(template_input, neg_dir / "AMPS_PARAM_C17.in",
                                  "C17_%s_charge_minus" % mover.lower(),
                                  "NEGATIVE_PROTON", -1, 1.0073, points, args)

        for label, wdir in (("charge_plus", pos_dir), ("charge_minus_reversed", neg_dir)):
            cmd = [
                args.mpirun, "-np", str(args.np), str(amps_path),
                "-mode", "gridless",
                "-i", "AMPS_PARAM_C17.in",
                "-mover", mover,
                "-gridless-mpi-scheduler", args.scheduler,
                "-gridless-mpi-dynamic-chunk", str(args.dynamic_chunk),
                "-density-parallel", "THREADS",
                "-density-threads", str(args.nt),
                "-cutoff-search", "UPPER_SCAN",
                "-cutoff-upper-scan-n", str(args.cutoff_upper_scan_n),
            ]
            log_file = wdir / ("C17_%s_%s_amps.log" % (mover.lower(), label))
            run_records.append({"mover": mover, "case": label, "workdir": str(wdir),
                                "command": cmd, "log_file": str(log_file)})
            if args.dry_run:
                print("[%s %s] %s" % (mover, label, " ".join(cmd)))
                continue
            if not args.skip_run:
                print("\nRunning C17 %s %s in %s" % (mover, label, wdir))
                print(" ".join(cmd))
                rc = run_command(cmd, cwd=wdir, log_path=log_file)
                if rc != 0:
                    overall_passed = False
                    all_messages.append("AMPS run failed for %s %s with exit code %d; see %s" %
                                        (mover, label, rc, log_file))

        if args.dry_run:
            continue

        if any("AMPS run failed for %s" % mover in msg for msg in all_messages):
            continue

        try:
            summary, pair_rows, passed, messages = compare_charge_reversal(
                points, pos_dir, neg_dir, rigidities, args.rel_tol, args.abs_tol_gv,
                ignore_poles=args.ignore_poles)
        except Exception as exc:
            overall_passed = False
            all_messages.append("C17 comparison failed for mover %s: %s" % (mover, exc))
            continue

        for r in summary:
            r["mover"] = mover
        for r in pair_rows:
            r["mover"] = mover
        all_summary.extend(summary)
        all_pair_rows.extend(pair_rows)
        if not passed:
            overall_passed = False
            all_messages.extend(["%s: %s" % (mover, m) for m in messages])

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "passed": overall_passed if not args.dry_run else None,
        "dry_run": args.dry_run,
        "points": points,
        "rigidities_GV": rigidities,
        "rel_tol": args.rel_tol,
        "abs_tol_GV": args.abs_tol_gv,
        "ignore_poles": args.ignore_poles,
        "run_records": run_records,
        "summary": all_summary,
        "messages": all_messages,
    }

    if args.dry_run:
        base_workdir.mkdir(parents=True, exist_ok=True)
        (base_workdir / "C17_result.json").write_text(json.dumps(result, indent=2))
        print("\nDry run complete.  Inputs and commands were prepared under %s" % base_workdir)
        return 0

    write_dict_csv(all_summary, base_workdir / "C17_summary.csv")
    write_dict_csv(all_pair_rows, base_workdir / "C17_pairwise_directional_residuals.csv")
    make_plot(all_pair_rows, base_workdir / "C17_residual_histogram.png")
    (base_workdir / "C17_result.json").write_text(json.dumps(result, indent=2))

    print("\nC17 summary written to:")
    print("  %s" % (base_workdir / "C17_summary.csv"))
    print("  %s" % (base_workdir / "C17_pairwise_directional_residuals.csv"))
    print("  %s" % (base_workdir / "C17_result.json"))

    if overall_passed:
        print("\nC17 PASS: charge-sign / velocity-reversal symmetry satisfied.")
        return 0
    else:
        print("\nC17 FAIL:")
        for msg in all_messages:
            print("  - " + msg)
        return 1


if __name__ == "__main__":
    sys.exit(main())
