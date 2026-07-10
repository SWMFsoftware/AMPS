#!/usr/bin/env python3
"""
C14 — Mode3D versus gridless cross-solver consistency.

Run from the directory that contains the AMPS executable:

    python srcEarth/test/C14/run_C14.py -np 4 -nt 16

C14 launches the same centered aligned-dipole vertical-cutoff shell problem
through the standalone Mode3D and gridless solver paths.  The primary reference
is cross-solver consistency.  The dipole case also has an independent analytical
vertical Størmer reference,

    Rc = R0 cos^4(lambda) / r_RE^2 .

The default grid is the validation-plan C14 shell set: altitudes 500 and
9000 km, longitudes every 30 degrees, and latitudes -60, -30, 0, +30, +60 deg.
"""

from __future__ import annotations

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
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

TEST_ID = "C14"
TEST_NAME = "Mode3D versus gridless cross-solver consistency"

RE_KM = 6371.2
STORMER_R0_GV = 0.299792458 * 0.25 * 3.12e-5 * (RE_KM * 1000.0)
DEFAULT_LONS_DEG = tuple(float(i) for i in range(0, 360, 30))
DEFAULT_LATS_DEG = (-60.0, -30.0, 0.0, 30.0, 60.0)
DEFAULT_ALTS_KM = (500.0, 9000.0)


class ShellRow:
    def __init__(self, alt_km: float, lon_deg: float, lat_deg: float,
                 rc_num_gv: float, rc_code_stormer_gv: Optional[float],
                 rel_err_code: Optional[float]):
        self.alt_km = alt_km
        self.lon_deg = lon_deg
        self.lat_deg = lat_deg
        self.rc_num_gv = rc_num_gv
        self.rc_code_stormer_gv = rc_code_stormer_gv
        self.rel_err_code = rel_err_code


def preprocess_negative_option_values(argv: Sequence[str]) -> List[str]:
    """Allow '--lats -60,-30,0' as well as '--lats=-60,-30,0'."""
    options = {"--lats", "--target-lats", "--alts", "--target-alts", "--lons", "--target-lons"}
    out: List[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in options and i + 1 < len(argv):
            b = argv[i + 1]
            if b.startswith("-") and not b.startswith("--"):
                out.append(a + "=" + b)
                i += 2
                continue
        out.append(a)
        i += 1
    return out


def parse_float_list(text: str, name: str) -> List[float]:
    vals: List[float] = []
    for item in str(text).split(','):
        item = item.strip()
        if not item:
            continue
        try:
            vals.append(float(item))
        except ValueError:
            raise SystemExit("Could not parse %s value '%s' in '%s'" % (name, item, text))
    if not vals:
        raise SystemExit("%s list is empty" % name)
    return vals


def normalize_lon(lon: float) -> float:
    x = lon % 360.0
    if abs(x - 360.0) < 1.0e-9 or abs(x) < 1.0e-9:
        return 0.0
    return x


def close_to(value: float, target: float, tol: float = 1.0e-6) -> bool:
    return abs(value - target) <= tol


def nearest_key(value: float, targets: Iterable[float], tol: float = 1.0e-6) -> Optional[float]:
    for t in targets:
        if abs(value - t) <= tol:
            return t
    return None


def nearest_lon_key(value: float, targets: Iterable[float], tol: float = 1.0e-6) -> Optional[float]:
    x = normalize_lon(value)
    for t in targets:
        y = normalize_lon(t)
        if abs(x - y) <= tol or abs(abs(x - y) - 360.0) <= tol:
            return y
    return None


def stormer_vertical_gv(lat_deg: float, alt_km: float) -> float:
    r_re = (RE_KM + alt_km) / RE_KM
    c4 = math.cos(math.radians(lat_deg)) ** 4
    return STORMER_R0_GV * c4 / (r_re * r_re)


def parse_variables(line: str) -> List[str]:
    return [n.strip().lower() for n in re.findall(r'"([^"]+)"', line)]


def normalize_name(name: str) -> str:
    return name.lower().replace(" ", "_").replace("-", "_").replace("/", "_")


def pick_column(variables: Optional[Sequence[str]], candidates: Sequence[str], fallback: int) -> int:
    if variables:
        norm = [normalize_name(v) for v in variables]
        for cand in candidates:
            c = normalize_name(cand)
            for i, v in enumerate(norm):
                if v == c or c in v:
                    return i
    return fallback


def parse_tecplot_shell_output(path: Path) -> List[ShellRow]:
    rows: List[ShellRow] = []
    current_alt: Optional[float] = None
    variables: Optional[List[str]] = None
    zone_alt_patterns = [
        re.compile(r'alt[_\s]*km\s*=\s*([0-9eE+\-.]+)', re.IGNORECASE),
        re.compile(r'alt\s*=\s*([0-9eE+\-.]+)\s*km', re.IGNORECASE),
        re.compile(r'altitude\s*=\s*([0-9eE+\-.]+)', re.IGNORECASE),
    ]

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
                for pat in zone_alt_patterns:
                    m = pat.search(line)
                    if m:
                        current_alt = float(m.group(1))
                        break
                continue
            if current_alt is None:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                lon_col = pick_column(variables, ["lon_deg", "longitude", "lon"], 0)
                lat_col = pick_column(variables, ["lat_deg", "latitude", "lat"], 1)
                rc_col = pick_column(
                    variables,
                    ["rc_num_gv", "rc_gv", "cutoff_gv", "cutoff_rigidity_gv", "rigidity_gv"],
                    5 if len(parts) > 5 else 2,
                )
                rc_ref_col = pick_column(
                    variables,
                    ["rc_vert_gv", "rc_stormer_gv", "rc_analytic_gv", "stormer_gv"],
                    6 if len(parts) > 6 else -1,
                )
                rel_col = pick_column(variables, ["rel_err", "relative_error"], 7 if len(parts) > 7 else -1)
                lon = float(parts[lon_col])
                lat = float(parts[lat_col])
                rc_num = float(parts[rc_col])
                rc_ref = float(parts[rc_ref_col]) if rc_ref_col >= 0 and rc_ref_col < len(parts) else None
                rel = float(parts[rel_col]) if rel_col >= 0 and rel_col < len(parts) else None
            except (ValueError, IndexError):
                continue
            rows.append(ShellRow(current_alt, lon, lat, rc_num, rc_ref, rel))

    if not rows:
        raise RuntimeError("No data rows were parsed from %s" % path)
    return rows


def select_target_means(rows: Sequence[ShellRow], lons: Sequence[float], lats: Sequence[float],
                        alts: Sequence[float]) -> Tuple[Dict[Tuple[float, float, float], Dict[str, float]], List[str]]:
    """Return one averaged cutoff value per requested (alt, lat, lon) target."""
    targets_lons = [normalize_lon(x) for x in lons]
    buckets: Dict[Tuple[float, float, float], List[ShellRow]] = {}
    for row in rows:
        alt_key = nearest_key(row.alt_km, alts)
        lat_key = nearest_key(row.lat_deg, lats)
        lon_key = nearest_lon_key(row.lon_deg, targets_lons)
        if alt_key is None or lat_key is None or lon_key is None:
            continue
        buckets.setdefault((alt_key, lat_key, lon_key), []).append(row)

    out: Dict[Tuple[float, float, float], Dict[str, float]] = {}
    messages: List[str] = []
    for alt in alts:
        for lat in lats:
            for lon in targets_lons:
                key = (alt, lat, lon)
                items = buckets.get(key, [])
                if not items:
                    messages.append("Missing shell output target alt=%g km, lat=%g deg, lon=%g deg" % (alt, lat, lon))
                    continue
                rc_vals = [x.rc_num_gv for x in items]
                code_refs = [x.rc_code_stormer_gv for x in items if x.rc_code_stormer_gv is not None]
                out[key] = {
                    "Rc_num_GV": sum(rc_vals) / len(rc_vals),
                    "Rc_num_min_GV": min(rc_vals),
                    "Rc_num_max_GV": max(rc_vals),
                    "n_rows": float(len(items)),
                    "Rc_code_reference_GV": (sum(code_refs) / len(code_refs)) if code_refs else float("nan"),
                }
    return out, messages


def write_reference_csv(path: Path, lons: Sequence[float], lats: Sequence[float], alts: Sequence[float],
                        cross_tol: float, cross_high_lat_tol: float,
                        analytic_mid_tol: float, analytic_high_lat_tol: float) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "alt_km", "lat_deg", "lon_deg", "Rc_stormer_GV", "cross_solver_tol",
            "analytic_stormer_tol", "check_type", "description",
        ])
        for alt in alts:
            for lat in lats:
                for lon in lons:
                    high = abs(lat) > 30.0
                    w.writerow([
                        "%.6g" % alt,
                        "%.6g" % lat,
                        "%.6g" % normalize_lon(lon),
                        "%.12e" % stormer_vertical_gv(lat, alt),
                        "%.6e" % (cross_high_lat_tol if high else cross_tol),
                        "%.6e" % (analytic_high_lat_tol if high else analytic_mid_tol),
                        "analytic_cross_solver_reference",
                        "C14 centered-dipole Mode3D/gridless/Størmer target",
                    ])


def write_csv_dicts(rows: Sequence[Dict[str, object]], path: Path) -> None:
    if not rows:
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def compare_solvers(mode3d: Dict[Tuple[float, float, float], Dict[str, float]],
                    gridless: Dict[Tuple[float, float, float], Dict[str, float]],
                    lons: Sequence[float], lats: Sequence[float], alts: Sequence[float],
                    args: argparse.Namespace) -> Tuple[List[Dict[str, object]], bool, List[str]]:
    summary: List[Dict[str, object]] = []
    messages: List[str] = []
    passed = True
    targets_lons = [normalize_lon(x) for x in lons]

    for alt in alts:
        for lat in lats:
            for lon in targets_lons:
                key = (alt, lat, lon)
                ref = stormer_vertical_gv(lat, alt)
                if key not in mode3d:
                    passed = False
                    messages.append("Mode3D missing target alt=%g lat=%g lon=%g" % key)
                    continue
                if key not in gridless:
                    passed = False
                    messages.append("Gridless missing target alt=%g lat=%g lon=%g" % key)
                    continue
                rc_3d = mode3d[key]["Rc_num_GV"]
                rc_gl = gridless[key]["Rc_num_GV"]
                high = abs(lat) > 30.0
                cross_tol = args.cross_high_lat_tol if high else args.cross_tol
                analytic_tol = args.analytic_high_lat_tol if high else args.analytic_mid_tol
                cross_rel = abs(rc_3d - rc_gl) / max(ref, 1.0e-30)
                rel_3d = (rc_3d - ref) / ref if ref > 0.0 else 0.0
                rel_gl = (rc_gl - ref) / ref if ref > 0.0 else 0.0
                row_pass = True
                if cross_rel > cross_tol:
                    passed = False
                    row_pass = False
                    messages.append(
                        "Mode3D/gridless mismatch alt=%g lat=%g lon=%g: rel_diff=%.3e > %.3e" %
                        (alt, lat, lon, cross_rel, cross_tol)
                    )
                if abs(rel_3d) > analytic_tol:
                    passed = False
                    row_pass = False
                    messages.append(
                        "Mode3D/Størmer residual alt=%g lat=%g lon=%g: abs_rel=%.3e > %.3e" %
                        (alt, lat, lon, abs(rel_3d), analytic_tol)
                    )
                if abs(rel_gl) > analytic_tol:
                    passed = False
                    row_pass = False
                    messages.append(
                        "Gridless/Størmer residual alt=%g lat=%g lon=%g: abs_rel=%.3e > %.3e" %
                        (alt, lat, lon, abs(rel_gl), analytic_tol)
                    )
                summary.append({
                    "check_type": "point_cross_solver_and_stormer",
                    "alt_km": alt,
                    "lat_deg": lat,
                    "lon_deg": lon,
                    "Rc_stormer_GV": ref,
                    "Rc_mode3d_GV": rc_3d,
                    "Rc_gridless_GV": rc_gl,
                    "mode3d_minus_gridless_over_stormer": (rc_3d - rc_gl) / max(ref, 1.0e-30),
                    "abs_cross_rel_diff": cross_rel,
                    "cross_tol": cross_tol,
                    "mode3d_rel_error_vs_stormer": rel_3d,
                    "gridless_rel_error_vs_stormer": rel_gl,
                    "analytic_tol": analytic_tol,
                    "passed": bool(row_pass),
                })

    # Longitude spread at fixed altitude/latitude.
    for solver_name, data, tol in [
        ("mode3d", mode3d, args.longitude_tol),
        ("gridless", gridless, args.longitude_tol),
    ]:
        for alt in alts:
            for lat in lats:
                vals = []
                for lon in targets_lons:
                    row = data.get((alt, lat, lon))
                    if row:
                        vals.append(row["Rc_num_GV"])
                if len(vals) < 2:
                    continue
                ref = stormer_vertical_gv(lat, alt)
                spread = (max(vals) - min(vals)) / max(ref, 1.0e-30)
                row_pass = spread <= tol
                if not row_pass:
                    passed = False
                    messages.append("%s longitude spread too large at alt=%g lat=%g: %.3e > %.3e" % (solver_name, alt, lat, spread, tol))
                summary.append({
                    "check_type": "longitude_symmetry",
                    "solver": solver_name,
                    "alt_km": alt,
                    "lat_deg": lat,
                    "lon_deg": "all",
                    "Rc_stormer_GV": ref,
                    "relative_spread": spread,
                    "tolerance": tol,
                    "passed": bool(row_pass),
                })

    # North/south symmetry at fixed altitude/longitude.
    lat_set = set(lats)
    for solver_name, data, tol in [
        ("mode3d", mode3d, args.ns_tol),
        ("gridless", gridless, args.ns_tol),
    ]:
        for alt in alts:
            for lat in [x for x in lats if x > 0.0 and -x in lat_set]:
                for lon in targets_lons:
                    north = data.get((alt, lat, lon))
                    south = data.get((alt, -lat, lon))
                    if not north or not south:
                        continue
                    ref = stormer_vertical_gv(lat, alt)
                    diff = abs(north["Rc_num_GV"] - south["Rc_num_GV"]) / max(ref, 1.0e-30)
                    row_pass = diff <= tol
                    if not row_pass:
                        passed = False
                        messages.append("%s north/south difference too large at alt=%g |lat|=%g lon=%g: %.3e > %.3e" % (solver_name, alt, lat, lon, diff, tol))
                    summary.append({
                        "check_type": "north_south_symmetry",
                        "solver": solver_name,
                        "alt_km": alt,
                        "abs_lat_deg": lat,
                        "lon_deg": lon,
                        "Rc_north_GV": north["Rc_num_GV"],
                        "Rc_south_GV": south["Rc_num_GV"],
                        "Rc_stormer_GV": ref,
                        "relative_difference": diff,
                        "tolerance": tol,
                        "passed": bool(row_pass),
                    })
    return summary, passed, messages


def render_input_template(template_path: Path, output_path: Path, args: argparse.Namespace,
                          alts: Sequence[float]) -> None:
    text = template_path.read_text()
    replacements = {
        "__CUTOFF_EMIN__": "%.12g" % args.cutoff_emin,
        "__CUTOFF_EMAX__": "%.12g" % args.cutoff_emax,
        "__MAX_TRACE_TIME__": "%.12g" % args.max_trace_time,
        "__MAX_TRACE_DISTANCE__": "%.12g" % args.max_trace_distance,
        "__DT_TRACE__": "%.12g" % args.dt_trace,
        "__SHELL_COUNT__": str(len(alts)),
        "__SHELL_ALTS_KM__": " ".join("%.12g" % x for x in alts),
        "__SHELL_RES_DEG__": "%.12g" % args.shell_res_deg,
    }
    for key, val in replacements.items():
        text = text.replace(key, val)
    text += (
        "\n! ── C14 harness run-time settings, supplied through CLI ─────────\n"
        "! C14_MODE3D_FIELD_EVAL       %s\n"
        "! C14_MPI_SCHEDULER           %s\n"
        "! C14_MPI_DYNAMIC_CHUNK       %d\n"
        "! C14_CUTOFF_SEARCH           UPPER_SCAN\n"
        "! C14_CUTOFF_SCAN_N           %d\n"
    ) % (args.mode3d_field_eval, args.scheduler, args.dynamic_chunk, args.cutoff_scan_n)
    output_path.write_text(text)


def run_command(cmd: List[str], cwd: Path, log_path: Path) -> int:
    with log_path.open("w") as log:
        log.write("Command:\n  " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
        return proc.wait()


def find_output_file(workdir: Path, mode: str, user_output: Optional[str]) -> Path:
    if user_output:
        p = Path(user_output)
        if not p.is_absolute():
            p = workdir / p
        return p
    if mode == "3d":
        candidates = ["cutoff_3d_shells_dipole_compare.dat", "cutoff_3d_shells.dat"]
    else:
        candidates = ["cutoff_gridless_shells_dipole_compare.dat", "cutoff_gridless_shells.dat"]
    for name in candidates:
        p = workdir / name
        if p.exists():
            return p
    matches = sorted(workdir.glob("cutoff_*shell*.dat"))
    if matches:
        return matches[0]
    return workdir / candidates[0]


def make_plot(summary: Sequence[Dict[str, object]], path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    rows = [r for r in summary if r.get("check_type") == "point_cross_solver_and_stormer"]
    if not rows:
        return
    x = [float(r["Rc_gridless_GV"]) for r in rows]
    y = [float(r["Rc_mode3d_GV"]) for r in rows]
    lo = min(x + y)
    hi = max(x + y)
    plt.figure(figsize=(5.2, 5.0))
    plt.plot(x, y, marker="o", linestyle="None", label="C14 points")
    plt.plot([lo, hi], [lo, hi], linestyle="--", label="1:1")
    plt.xlabel("Gridless Rc (GV)")
    plt.ylabel("Mode3D Rc (GV)")
    plt.title("C14: Mode3D/gridless cross-solver consistency")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        C14: Mode3D versus gridless cross-solver consistency

        Runs the same centered aligned-dipole vertical cutoff shell problem with
        Mode3D and gridless solvers.  The primary check is Mode3D/gridless
        agreement; the dipole case is also compared with the analytical vertical
        Størmer cutoff.
        """),
        epilog=textwrap.dedent("""
        Examples:
          python srcEarth/test/C14/run_C14.py
        W python srcEarth/test/C14/run_C14.py -np 4 -nt 16
          python srcEarth/test/C14/run_C14.py --mode3d-field-eval MESH -np 4 -nt 16
          python srcEarth/test/C14/run_C14.py --scheduler STATIC --dynamic-chunk 0
          python srcEarth/test/C14/run_C14.py --lons 0,90,180,270 --lats -60,-30,0,30,60
          python srcEarth/test/C14/run_C14.py --dry-run
          python srcEarth/test/C14/run_C14.py --skip-run --workdir test_output/C14_cross_solver
        """)
    )
    parser.add_argument("-np", type=int, default=4, help="number of MPI processes passed to mpirun; default: 4")
    parser.add_argument("-nt", type=int, default=16, help="number of threads per MPI rank; default: 16")
    parser.add_argument("--amps", default="./amps", help="AMPS executable path, relative to launch directory; default: ./amps")
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher executable; default: mpirun")
    parser.add_argument("--workdir", default=None, help="directory where test runs and outputs are written; default: test_output/C14_cross_solver")
    parser.add_argument("--mode3d-field-eval", "--field-eval", dest="mode3d_field_eval", default="ANALYTIC", choices=["ANALYTIC", "MESH", "GRID_3D"], help="Mode3D trajectory field backend; default: ANALYTIC")
    parser.add_argument("--scheduler", default="DYNAMIC", choices=["DYNAMIC", "BLOCK_CYCLIC", "STATIC"], help="MPI scheduler for both solvers; default: DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0, help="dynamic MPI chunk size; 0 means AMPS auto heuristic; default: 0")
    parser.add_argument("--cutoff-scan-n", type=int, default=100, help="number of UPPER_SCAN rigidity samples passed to AMPS; default: 100")
    parser.add_argument("--cutoff-emin", type=float, default=1.0, help="CUTOFF_EMIN in MeV/n; default: 1")
    parser.add_argument("--cutoff-emax", type=float, default=20000.0, help="CUTOFF_EMAX in MeV/n; default: 20000")
    parser.add_argument("--lons", "--target-lons", dest="lons", default=",".join(str(int(x)) for x in DEFAULT_LONS_DEG), help="comma-separated longitudes; default: 0,30,...,330")
    parser.add_argument("--lats", "--target-lats", dest="lats", default=",".join(str(int(x)) for x in DEFAULT_LATS_DEG), help="comma-separated latitudes; default: -60,-30,0,30,60")
    parser.add_argument("--alts", "--target-alts", dest="alts", default=",".join(str(int(x)) for x in DEFAULT_ALTS_KM), help="comma-separated altitudes in km; default: 500,9000")
    parser.add_argument("--shell-res-deg", type=float, default=30.0, help="SHELL_RES_DEG written to both inputs; default: 30")
    parser.add_argument("--dt-trace", type=float, default=1.0, help="DT_TRACE written into both inputs; default: 1 s")
    parser.add_argument("--max-trace-time", type=float, default=600.0, help="MAX_TRACE_TIME/CUTOFF_MAX_TRAJ_TIME; default: 600 s")
    parser.add_argument("--max-trace-distance", type=float, default=400.0, help="MAX_TRACE_DISTANCE in Re; default: 400")
    parser.add_argument("--cross-tol", type=float, default=1.0e-2, help="Mode3D/gridless relative-difference tolerance for |lat|<=30 deg; default: 1e-2")
    parser.add_argument("--cross-high-lat-tol", type=float, default=5.0e-2, help="Mode3D/gridless relative-difference tolerance for |lat|>30 deg; default: 5e-2")
    parser.add_argument("--analytic-mid-tol", type=float, default=5.0e-2, help="Størmer relative-error tolerance for |lat|<=30 deg; default: 5e-2")
    parser.add_argument("--analytic-high-lat-tol", type=float, default=2.5e-1, help="Størmer relative-error tolerance for |lat|>30 deg; default: 0.25")
    parser.add_argument("--longitude-tol", type=float, default=3.0e-2, help="longitude spread tolerance at fixed alt/lat; default: 3e-2")
    parser.add_argument("--ns-tol", type=float, default=3.0e-2, help="north/south relative-difference tolerance; default: 3e-2")
    parser.add_argument("--mode3d-output-file", default=None, help="explicit Mode3D shell output file; default: auto-detect")
    parser.add_argument("--gridless-output-file", default=None, help="explicit gridless shell output file; default: auto-detect")
    parser.add_argument("--no-cutoff-search-cli", action="store_true", help="do not pass -cutoff-search/-cutoff-upper-scan-n; useful for older CLI checkouts")
    parser.add_argument("--skip-run", action="store_true", help="do not run AMPS; analyze existing mode3d/gridless subdirectories under --workdir")
    parser.add_argument("--keep", action="store_true", help="keep an existing work directory instead of replacing it")
    parser.add_argument("--dry-run", action="store_true", help="create rendered inputs and print AMPS commands without executing them")
    return parser.parse_args(preprocess_negative_option_values(list(argv) if argv is not None else sys.argv[1:]))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.np < 1:
        raise SystemExit("-np must be >= 1")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")
    if args.cutoff_scan_n < 2:
        raise SystemExit("--cutoff-scan-n must be >= 2")
    if args.cutoff_emin <= 0.0 or args.cutoff_emax <= args.cutoff_emin:
        raise SystemExit("cutoff energy bounds must satisfy 0 < CUTOFF_EMIN < CUTOFF_EMAX")
    if args.shell_res_deg <= 0.0:
        raise SystemExit("--shell-res-deg must be > 0")
    if args.max_trace_distance < 0.0:
        raise SystemExit("--max-trace-distance must be >= 0")

    lons = parse_float_list(args.lons, "lons")
    lats = parse_float_list(args.lats, "lats")
    alts = parse_float_list(args.alts, "alts")

    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    root_workdir = Path(args.workdir) if args.workdir else Path("test_output") / "C14_cross_solver"
    if not root_workdir.is_absolute():
        root_workdir = (launch_dir / root_workdir).resolve()

    mode3d_dir = root_workdir / "mode3d"
    gridless_dir = root_workdir / "gridless"

    if not args.skip_run and not args.keep and root_workdir.exists():
        shutil.rmtree(root_workdir)
    root_workdir.mkdir(parents=True, exist_ok=True)

    generated_reference = root_workdir / "reference_C14_cross_solver_generated.csv"
    write_reference_csv(generated_reference, lons, lats, alts, args.cross_tol,
                        args.cross_high_lat_tol, args.analytic_mid_tol,
                        args.analytic_high_lat_tol)

    planned_commands: List[Dict[str, object]] = []
    for mode, case_dir, template_name in [
        ("3d", mode3d_dir, "AMPS_PARAM_C14_mode3d.in"),
        ("gridless", gridless_dir, "AMPS_PARAM_C14_gridless.in"),
    ]:
        if not args.skip_run:
            case_dir.mkdir(parents=True, exist_ok=True)
            render_input_template(script_dir / template_name, case_dir / "AMPS_PARAM_C14.in", args, alts)
            shutil.copy2(generated_reference, case_dir / "reference_C14_cross_solver.csv")
            amps_path = Path(args.amps)
            if not amps_path.is_absolute():
                amps_path = (launch_dir / amps_path).resolve()
            cmd = [args.mpirun, "-np", str(args.np), str(amps_path), "-mode", mode, "-i", "AMPS_PARAM_C14.in"]
            if not args.no_cutoff_search_cli:
                cmd += ["-cutoff-search", "UPPER_SCAN", "-cutoff-upper-scan-n", str(args.cutoff_scan_n)]
            if mode == "3d":
                cmd += [
                    "-mode3d-field-eval", args.mode3d_field_eval,
                    "-mode3d-parallel", "THREADS",
                    "-mode3d-threads", str(args.nt),
                    "-mode3d-mpi-scheduler", args.scheduler,
                    "-mode3d-mpi-dynamic-chunk", str(args.dynamic_chunk),
                ]
            else:
                cmd += [
                    "-gridless-mpi-scheduler", args.scheduler,
                    "-gridless-mpi-dynamic-chunk", str(args.dynamic_chunk),
                    "-density-parallel", "THREADS",
                    "-density-threads", str(args.nt),
                ]
            planned_commands.append({"mode": mode, "cwd": str(case_dir), "command": cmd})
            if args.dry_run:
                print("Dry run %s/%s in %s" % (TEST_ID, mode, case_dir))
                print(" ".join(cmd))
                continue
            print("Running %s/%s in %s" % (TEST_ID, mode, case_dir))
            print(" ".join(cmd))
            rc = run_command(cmd, cwd=case_dir, log_path=case_dir / ("C14_%s_amps.log" % mode))
            if rc != 0:
                print("AMPS command failed with exit code %d for %s. See %s" % (rc, mode, case_dir), file=sys.stderr)
                return rc
        else:
            if not case_dir.exists():
                raise SystemExit("Missing case directory for --skip-run: %s" % case_dir)

    if args.dry_run:
        planned_path = root_workdir / "C14_dry_run_commands.json"
        with planned_path.open("w") as f:
            json.dump(planned_commands, f, indent=2)
            f.write("\n")
        print("Dry run complete. Rendered inputs are under %s" % root_workdir)
        print("Planned-command JSON: %s" % planned_path)
        return 0

    mode3d_output = find_output_file(mode3d_dir, "3d", args.mode3d_output_file)
    gridless_output = find_output_file(gridless_dir, "gridless", args.gridless_output_file)
    if not mode3d_output.exists():
        print("Expected Mode3D shell output was not found: %s" % mode3d_output, file=sys.stderr)
        return 2
    if not gridless_output.exists():
        print("Expected gridless shell output was not found: %s" % gridless_output, file=sys.stderr)
        return 2

    mode3d_rows = parse_tecplot_shell_output(mode3d_output)
    gridless_rows = parse_tecplot_shell_output(gridless_output)
    mode3d_targets, mode3d_messages = select_target_means(mode3d_rows, lons, lats, alts)
    gridless_targets, gridless_messages = select_target_means(gridless_rows, lons, lats, alts)
    summary, passed, messages = compare_solvers(mode3d_targets, gridless_targets, lons, lats, alts, args)
    if mode3d_messages:
        passed = False
        messages.extend(["Mode3D: " + x for x in mode3d_messages])
    if gridless_messages:
        passed = False
        messages.extend(["Gridless: " + x for x in gridless_messages])

    summary_csv = root_workdir / "C14_summary.csv"
    result_json = root_workdir / "C14_result.json"
    plot_png = root_workdir / "C14_mode3d_vs_gridless.png"
    write_csv_dicts(summary, summary_csv)
    make_plot(summary, plot_png)

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "passed": bool(passed),
        "messages": messages,
        "np": args.np,
        "nt": args.nt,
        "mode3d_field_eval": args.mode3d_field_eval,
        "scheduler": args.scheduler,
        "dynamic_chunk": args.dynamic_chunk,
        "cutoff_search_cli_enabled": not args.no_cutoff_search_cli,
        "cutoff_scan_n": args.cutoff_scan_n,
        "target_lons_deg": [normalize_lon(x) for x in lons],
        "target_lats_deg": lats,
        "target_alts_km": alts,
        "workdir": str(root_workdir),
        "mode3d_output_file": str(mode3d_output),
        "gridless_output_file": str(gridless_output),
        "summary_csv": str(summary_csv),
        "reference_csv": str(generated_reference),
        "plot_png": str(plot_png) if plot_png.exists() else None,
    }
    with result_json.open("w") as f:
        json.dump(result, f, indent=2)
        f.write("\n")

    print("\nC14 summary")
    print("===========")
    print("Mode3D output:   %s" % mode3d_output)
    print("Gridless output: %s" % gridless_output)
    point_rows = [r for r in summary if r.get("check_type") == "point_cross_solver_and_stormer"]
    if point_rows:
        max_cross = max(float(r["abs_cross_rel_diff"]) for r in point_rows)
        max_3d = max(abs(float(r["mode3d_rel_error_vs_stormer"])) for r in point_rows)
        max_gl = max(abs(float(r["gridless_rel_error_vs_stormer"])) for r in point_rows)
        print("Max |Mode3D-gridless|/Stormer: %.3e" % max_cross)
        print("Max |Mode3D-Stormer|/Stormer:  %.3e" % max_3d)
        print("Max |Gridless-Stormer|/Stormer: %.3e" % max_gl)
    print("Wrote: %s" % summary_csv)
    print("Wrote: %s" % result_json)
    if plot_png.exists():
        print("Wrote: %s" % plot_png)
    if messages:
        print("\nMessages:")
        for msg in messages:
            print("  - " + msg)
    print("\nRESULT: %s" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
