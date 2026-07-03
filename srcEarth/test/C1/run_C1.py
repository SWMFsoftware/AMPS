#!/usr/bin/env python3
"""
C1 — Pure dipole vertical Stormer cutoff validation test.

The script is intentionally self-contained so it can be run from the directory
that contains the AMPS executable:

    python srcEarth/test/C1/run_C1.py -np 4 -nt 16

C1 can now be run in either standalone Mode3D or gridless mode:

    python srcEarth/test/C1/run_C1.py --mode 3d
    python srcEarth/test/C1/run_C1.py --mode gridless

For Mode3D, the script also exposes the field-evaluation backend used by the
trajectory tracer.  The common production/validation command is:

    mpirun -np 8 ./amps -mode 3d -i AMPS_PARAM_C1.in -mode3d-field-eval ANALYTIC

Using ANALYTIC bypasses mesh interpolation for the dipole field and therefore
isolates the trajectory pusher, cutoff search, and classifier.  Using MESH keeps
the original mesh-stored Mode3D path and validates interpolation effects.
"""

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
from typing import Dict, Iterable, List, Tuple

TEST_ID = "C1"
TEST_NAME = "Pure dipole vertical Stormer cutoff"

# These constants intentionally match the dipole/Størmer convention used by the
# AMPS cutoff comparison writers.  The coefficient is equivalent to
#   R0 = 0.299792458 * 0.25 * B_eq(Re) * Re
# with B_eq(Re)=3.12e-5 T and Re=6371.2 km, yielding about 14.8983 GV.
RE_KM = 6371.2
STORMER_R0_GV = 0.299792458 * 0.25 * 3.12e-5 * (RE_KM * 1000.0)
TARGET_LATS_DEG = (-60.0, -30.0, 0.0, 30.0, 60.0)
TARGET_ALTS_KM = (500.0, 9000.0)


class ShellRow(object):
    """Container for one parsed Tecplot row.

    This intentionally avoids dataclasses and postponed annotations so that the
    test harness works with older Python 3 installations commonly found on HPC
    systems, including Python 3.6.
    """
    def __init__(self, alt_km, lon_deg, lat_deg, rc_num_gv, rc_code_stormer_gv, rel_err_code):
        self.alt_km = alt_km
        self.lon_deg = lon_deg
        self.lat_deg = lat_deg
        self.rc_num_gv = rc_num_gv
        self.rc_code_stormer_gv = rc_code_stormer_gv
        self.rel_err_code = rel_err_code


def stormer_vertical_gv(lat_deg: float, alt_km: float) -> float:
    """Independent analytical C1 reference used by the test harness."""
    r_re = (RE_KM + alt_km) / RE_KM
    c4 = math.cos(math.radians(lat_deg)) ** 4
    return STORMER_R0_GV * c4 / (r_re * r_re)


def write_reference_csv(path: Path) -> None:
    """Write the analytical reference table used by this test."""
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["alt_km", "lat_deg", "Rc_stormer_GV"])
        for alt in TARGET_ALTS_KM:
            for lat in TARGET_LATS_DEG:
                w.writerow(["%.1f" % alt, "%.1f" % lat, "%.9e" % stormer_vertical_gv(lat, alt)])


def _parse_variables(line):
    """Return a list of Tecplot variable names from a VARIABLES line, if present."""
    names = re.findall(r'"([^"]+)"', line)
    return [n.strip().lower() for n in names]


def _pick_column(variables, candidates, fallback):
    """Find a variable column by name, falling back to a known historical index."""
    if variables:
        norm = []
        for v in variables:
            norm.append(v.lower().replace(" ", "_").replace("-", "_"))
        for cand in candidates:
            c = cand.lower().replace(" ", "_").replace("-", "_")
            for i, v in enumerate(norm):
                if v == c or c in v:
                    return i
    return fallback


def parse_tecplot_shell_output(path: Path) -> List[ShellRow]:
    """
    Parse Mode3D or gridless shell cutoff output.

    The newest Mode3D dipole-comparison writer emits

        lon_deg, lat_deg, x_km, y_km, z_km, Rc_num_GV, Rc_vert_GV, rel_err

    while gridless or non-comparison writers may emit only a numerical cutoff
    column.  This parser therefore uses variable names when available and falls
    back to the historical convention that the first two columns are lon/lat and
    the numerical cutoff is column 5.
    """
    rows = []
    current_alt = None
    variables = None

    # Accept common zone labels such as T="alt_km=9000" or T="alt=9000 km".
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
                variables = _parse_variables(line)
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
                lon_col = _pick_column(variables, ["lon_deg", "longitude", "lon"], 0)
                lat_col = _pick_column(variables, ["lat_deg", "latitude", "lat"], 1)
                rc_col = _pick_column(
                    variables,
                    ["rc_num_gv", "rc_gv", "cutoff_gv", "cutoff_rigidity_gv", "rigidity_gv"],
                    5 if len(parts) > 5 else 2,
                )
                rc_ref_col = _pick_column(
                    variables,
                    ["rc_vert_gv", "rc_stormer_gv", "rc_analytic_gv", "stormer_gv"],
                    6 if len(parts) > 6 else -1,
                )
                rel_col = _pick_column(variables, ["rel_err", "relative_error"], 7 if len(parts) > 7 else -1)

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


def close_to_any(value: float, targets: Iterable[float], tol: float = 1.0e-6) -> bool:
    return any(abs(value - t) <= tol for t in targets)


def summarize(rows: List[ShellRow], mode: str, mode3d_field_eval: str) -> Tuple[List[Dict[str, float]], bool, List[str]]:
    """
    Build per-altitude/per-latitude summary and apply C1 acceptance criteria.

    The acceptance criteria intentionally separate mid-latitude and high-latitude
    points.  High-latitude cutoff is more sensitive to penumbra, finite-domain
    escape classification, and mesh interpolation than the clean analytic
    Størmer expression.  Direct analytic field evaluation should be the tightest
    case; mesh-backed Mode3D is allowed a slightly looser margin.
    """
    summary = []
    messages = []
    passed = True

    direct_field = (mode == "gridless") or (mode == "3d" and mode3d_field_eval.upper() == "ANALYTIC")

    for alt in TARGET_ALTS_KM:
        for lat in TARGET_LATS_DEG:
            selected = [r for r in rows if abs(r.alt_km - alt) < 1.0e-6 and abs(r.lat_deg - lat) < 1.0e-6]
            if not selected:
                passed = False
                messages.append("Missing output rows for alt=%g km, lat=%g deg" % (alt, lat))
                continue

            ref = stormer_vertical_gv(lat, alt)
            rc_values = [r.rc_num_gv for r in selected]
            code_refs = [r.rc_code_stormer_gv for r in selected if r.rc_code_stormer_gv is not None]
            mean_rc = sum(rc_values) / len(rc_values)
            min_rc = min(rc_values)
            max_rc = max(rc_values)
            mean_code_ref = sum(code_refs) / len(code_refs) if code_refs else None
            rel_mean = (mean_rc - ref) / ref if ref > 0.0 else 0.0
            rel_max_abs = max(abs((x - ref) / ref) for x in rc_values) if ref > 0.0 else 0.0
            lon_spread = (max_rc - min_rc) / ref if ref > 0.0 else 0.0
            rel_code_ref = abs((mean_code_ref - ref) / ref) if mean_code_ref is not None and ref > 0.0 else None

            # Tight checks for the independent analytical expression written by
            # AMPS itself, when that column is present.  Gridless output may not
            # include it, so the script keeps its own independent reference.
            if rel_code_ref is not None and rel_code_ref > 5.0e-3:
                passed = False
                messages.append(
                    "AMPS analytic Rc column differs from C1 reference by %.3e at alt=%g km, lat=%g deg" %
                    (rel_code_ref, alt, lat)
                )

            if direct_field:
                mean_tol_mid = 2.0e-2
                max_tol_mid = 5.0e-2
                mean_tol_high = 1.5e-1
                max_tol_high = 2.5e-1
            else:
                mean_tol_mid = 5.0e-2
                max_tol_mid = 1.0e-1
                mean_tol_high = 2.5e-1
                max_tol_high = 3.5e-1

            if abs(lat) <= 30.0:
                mean_tol = mean_tol_mid
                max_tol = max_tol_mid
            else:
                mean_tol = mean_tol_high
                max_tol = max_tol_high

            if abs(rel_mean) > mean_tol or rel_max_abs > max_tol:
                passed = False
                messages.append(
                    "Numerical Rc outside tolerance at alt=%g km, lat=%g deg: "
                    "mean_rel=%.3e (tol %.2e), max_abs_rel=%.3e (tol %.2e)" %
                    (alt, lat, rel_mean, mean_tol, rel_max_abs, max_tol)
                )

            summary.append({
                "mode": mode,
                "mode3d_field_eval": mode3d_field_eval if mode == "3d" else "N/A",
                "alt_km": alt,
                "lat_deg": lat,
                "n_lon": float(len(selected)),
                "Rc_reference_GV": ref,
                "Rc_code_reference_mean_GV": mean_code_ref if mean_code_ref is not None else float("nan"),
                "Rc_num_mean_GV": mean_rc,
                "Rc_num_min_GV": min_rc,
                "Rc_num_max_GV": max_rc,
                "rel_err_mean": rel_mean,
                "rel_err_max_abs": rel_max_abs,
                "lon_spread_over_ref": lon_spread,
            })

    return summary, passed, messages


def write_summary_csv(summary: List[Dict[str, float]], path: Path) -> None:
    if not summary:
        return
    keys = list(summary[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in summary:
            w.writerow(row)


def make_plot(summary: List[Dict[str, float]], path: Path) -> None:
    """Generate a compact comparison plot when matplotlib is available."""
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    by_alt = {}
    for row in summary:
        by_alt.setdefault(row["alt_km"], []).append(row)

    plt.figure(figsize=(7.0, 4.5))
    for alt, items in sorted(by_alt.items()):
        items = sorted(items, key=lambda x: x["lat_deg"])
        lat = [x["lat_deg"] for x in items]
        ref = [x["Rc_reference_GV"] for x in items]
        num = [x["Rc_num_mean_GV"] for x in items]
        plt.plot(lat, ref, marker="o", linestyle="-", label="Stormer reference, alt=%g km" % alt)
        plt.plot(lat, num, marker="x", linestyle="--", label="AMPS mean, alt=%g km" % alt)
    plt.xlabel("Magnetic latitude (deg)")
    plt.ylabel("Vertical cutoff rigidity (GV)")
    plt.title("C1: pure dipole vertical Stormer cutoff")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def run_command(cmd: List[str], cwd: Path, log_path: Path) -> int:
    """Run AMPS and tee stdout/stderr into a log file."""
    with log_path.open("w") as log:
        log.write("Command:\n  " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
        return proc.wait()


def render_input_template(template_path, output_path, nt, scheduler, dynamic_chunk):
    """Copy a parser-compatible C1 input template into the run directory.

    The active input keywords are intentionally kept close to the known-working
    AMPS_PARAM_test.in layout.  Current-code controls that have caused parser
    trouble in some checkouts (scheduler, thread count, cutoff-search algorithm,
    Mode3D field-evaluation backend) are passed through the AMPS command line by
    run_C1.py.  A short commented provenance block is appended so the generated
    AMPS_PARAM_C1.in remains self-documenting without adding active parser
    keywords.
    """
    text = template_path.read_text()
    text += (
        "\n! ── C1 harness run-time settings, supplied through CLI ─────────\n"
        "! C1_NT                  %d\n"
        "! C1_MPI_SCHEDULER       %s\n"
        "! C1_MPI_DYNAMIC_CHUNK   %d\n"
        "! C1_CUTOFF_SEARCH       UPPER_SCAN\n"
        "! C1_CUTOFF_SCAN_N       80\n"
    ) % (nt, scheduler, dynamic_chunk)
    output_path.write_text(text)


def find_output_file(workdir, mode, user_output):
    if user_output:
        p = Path(user_output)
        if not p.is_absolute():
            p = workdir / p
        return p

    if mode == "3d":
        candidates = [
            "cutoff_3d_shells_dipole_compare.dat",
            "cutoff_3d_shells.dat",
        ]
    else:
        candidates = [
            "cutoff_gridless_shells_dipole_compare.dat",
            "cutoff_gridless_shells.dat",
        ]
    for name in candidates:
        p = workdir / name
        if p.exists():
            return p

    matches = sorted(workdir.glob("cutoff_*shell*.dat"))
    if matches:
        return matches[0]
    return workdir / candidates[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        C1: Pure dipole vertical Stormer cutoff

        This test runs the cutoff calculation for a centered aligned dipole and
        vertical arrival directions.  The numerical result is compared with the
        analytical Stormer vertical cutoff,

            Rc = R0 cos^4(lambda) / r_RE^2 .

        The test can be run in either standalone Mode3D or gridless mode.  In
        Mode3D, --mode3d-field-eval controls whether the trajectory field is
        obtained analytically or from the Mode3D mesh.  ANALYTIC is the cleanest
        pusher/classifier test; MESH validates the mesh-interpolation path.

        The script is expected to be launched from the directory containing the
        AMPS executable.  AMPS itself is launched with mpirun.  The default
        parallel configuration is -np 4 -nt 16.
        """),
        epilog=textwrap.dedent("""
        Examples:
          python srcEarth/test/C1/run_C1.py
        W python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval ANALYTIC -np 18 -nt 16
        W python srcEarth/test/C1/run_C1.py --mode 3d --mode3d-field-eval MESH -np 18 -nt 16 
        W python srcEarth/test/C1/run_C1.py --mode gridless -np 4 -nt 16
          python srcEarth/test/C1/run_C1.py --no-cutoff-search-cli  # fallback for older CLI
          python srcEarth/test/C1/run_C1.py --skip-run --workdir test_output/C1_3d
        """)
    )
    parser.add_argument("-np", type=int, default=4, help="number of MPI processes passed to mpirun; default: 4")
    parser.add_argument("-nt", type=int, default=16, help="number of threads per MPI rank; default: 16")
    parser.add_argument("--mode", choices=["3d", "gridless"], default="3d", help="AMPS mode to validate; default: 3d")
    parser.add_argument("--mode3d-field-eval", "--field-eval", dest="mode3d_field_eval", default="ANALYTIC", choices=["ANALYTIC", "MESH", "GRID_3D"], help="Mode3D field backend passed as -mode3d-field-eval; default: ANALYTIC")
    parser.add_argument("--scheduler", default="DYNAMIC", choices=["DYNAMIC", "BLOCK_CYCLIC", "STATIC"], help="MPI scheduler to use for the selected mode; default: DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0, help="dynamic MPI chunk size; 0 means AMPS auto heuristic; default: 0")
    parser.add_argument("--no-cutoff-search-cli", action="store_true", help="do not pass -cutoff-search/-cutoff-upper-scan-n; useful for older parser/CLI checkouts")
    parser.add_argument("--amps", default="./amps", help="AMPS executable path, relative to the launch directory; default: ./amps")
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher executable; default: mpirun")
    parser.add_argument("--workdir", default=None, help="directory where the test is run and outputs are written; default: test_output/C1_<mode>")
    parser.add_argument("--output-file", default=None, help="explicit AMPS cutoff shell output file to parse; default: auto-detect")
    parser.add_argument("--skip-run", action="store_true", help="do not run AMPS; analyze existing output in --workdir")
    parser.add_argument("--keep", action="store_true", help="keep an existing work directory instead of replacing it")
    parser.add_argument("--run-in-place", action="store_true", help="run AMPS in the launch directory rather than in --workdir; output files may be overwritten")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.np < 1:
        raise SystemExit("-np must be >= 1")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")

    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    if args.mode == "3d":
        template_input = script_dir / "AMPS_PARAM_C1_mode3d.in"
    else:
        template_input = script_dir / "AMPS_PARAM_C1_gridless.in"
    # Backward-compatible fallback for older checkouts of this test directory.
    if not template_input.exists():
        template_input = script_dir / "AMPS_PARAM_C1.in"
    reference_csv = script_dir / "reference_C1_stormer.csv"

    if args.workdir is None:
        if args.mode == "3d":
            args.workdir = "test_output/C1_3d_%s" % args.mode3d_field_eval.lower()
        else:
            args.workdir = "test_output/C1_gridless"

    workdir = (launch_dir / args.workdir).resolve()
    if not args.skip_run and not args.run_in_place:
        if workdir.exists() and not args.keep:
            shutil.rmtree(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        render_input_template(template_input, workdir / "AMPS_PARAM_C1.in", args.nt, args.scheduler, args.dynamic_chunk)
        shutil.copy2(reference_csv, workdir / "reference_C1_stormer.csv")
    elif args.run_in_place:
        workdir = launch_dir
        render_input_template(template_input, workdir / "AMPS_PARAM_C1.in", args.nt, args.scheduler, args.dynamic_chunk)
        shutil.copy2(reference_csv, workdir / "reference_C1_stormer.csv")
    else:
        if not workdir.exists():
            raise SystemExit("--skip-run requested but workdir does not exist: %s" % workdir)

    log_file = workdir / "C1_amps.log"
    result_json = workdir / "C1_result.json"
    summary_csv = workdir / "C1_summary.csv"
    plot_png = workdir / "C1_stormer_comparison.png"

    if not args.skip_run:
        amps_path = Path(args.amps)
        if not amps_path.is_absolute():
            amps_path = (launch_dir / amps_path).resolve()
        cmd = [
            args.mpirun,
            "-np", str(args.np),
            str(amps_path),
            "-mode", args.mode,
            "-i", "AMPS_PARAM_C1.in",
        ]
        if not args.no_cutoff_search_cli:
            cmd += [
                "-cutoff-search", "UPPER_SCAN",
                "-cutoff-upper-scan-n", "80",
            ]
        if args.mode == "3d":
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
        print("Running %s in %s" % (TEST_ID, workdir))
        print(" ".join(cmd))
        rc = run_command(cmd, cwd=workdir, log_path=log_file)
        if rc != 0:
            print("AMPS command failed with exit code %d. See %s" % (rc, log_file), file=sys.stderr)
            return rc

    output_file = find_output_file(workdir, args.mode, args.output_file)
    if not output_file.exists():
        print("Expected AMPS output was not found: %s" % output_file, file=sys.stderr)
        return 2

    rows = parse_tecplot_shell_output(output_file)
    summary, passed, messages = summarize(rows, args.mode, args.mode3d_field_eval)
    write_summary_csv(summary, summary_csv)
    make_plot(summary, plot_png)

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "passed": passed,
        "messages": messages,
        "np": args.np,
        "nt": args.nt,
        "mode": args.mode,
        "mode3d_field_eval": args.mode3d_field_eval if args.mode == "3d" else None,
        "scheduler": args.scheduler,
        "dynamic_chunk": args.dynamic_chunk,
        "cutoff_search_cli_enabled": not args.no_cutoff_search_cli,
        "workdir": str(workdir),
        "output_file": str(output_file),
        "summary_csv": str(summary_csv),
        "plot_png": str(plot_png) if plot_png.exists() else None,
        "stormer_R0_GV": STORMER_R0_GV,
        "summary": summary,
    }
    result_json.write_text(json.dumps(result, indent=2))

    print("\nC1 summary")
    print("==========")
    print("mode=%s np=%d nt=%d scheduler=%s" % (args.mode, args.np, args.nt, args.scheduler))
    if args.mode == "3d":
        print("mode3d_field_eval=%s" % args.mode3d_field_eval)
    for row in summary:
        print(
            "alt=%7.1f km lat=%6.1f deg Rc_num_mean=%10.5e GV Rc_ref=%10.5e GV "
            "rel_mean=% .3e max_abs_rel=% .3e" %
            (row['alt_km'], row['lat_deg'], row['Rc_num_mean_GV'], row['Rc_reference_GV'],
             row['rel_err_mean'], row['rel_err_max_abs'])
        )
    if messages:
        print("\nMessages:")
        for m in messages:
            print("  - " + m)
    print("\nWrote: %s" % summary_csv)
    print("Wrote: %s" % result_json)
    if plot_png.exists():
        print("Wrote: %s" % plot_png)
    print("\nRESULT: %s" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
