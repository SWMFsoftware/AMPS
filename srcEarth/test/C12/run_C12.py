#!/usr/bin/env python3
"""
C12 — Particle-mover cross-validation against analytical Størmer cutoff.

Run from the directory that contains the AMPS executable:

    python srcEarth/test/C12/run_C12.py -np 4 -nt 16

C12 is a gridless cutoff test.  It runs the same centered-dipole vertical-cutoff
calculation with every selected particle mover and compares the numerical cutoff
with the analytical Størmer vertical cutoff,

    Rc = R0 cos^4(lambda) / r_RE^2 .

The full-orbit movers are checked on the standard outer shell used by C1/C3
(default 9000 km) and on a high shell.  The guiding-center movers are checked on
that high shell only by default, because the guiding-center approximation is not
intended to be validated very close to the inner boundary.
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

TEST_ID = "C12"
TEST_NAME = "Particle-mover cross-validation against Størmer"

RE_KM = 6371.2
STORMER_R0_GV = 0.299792458 * 0.25 * 3.12e-5 * (RE_KM * 1000.0)

FULL_ORBIT_MOVERS = ("BORIS", "RK2", "RK4", "RK6")
GC_MOVERS = ("GC2", "GC4", "GC6")
DEFAULT_MOVERS = FULL_ORBIT_MOVERS + GC_MOVERS

# The acceptance levels follow the v8 validation-campaign intent.  They are
# intentionally applied only to the shell where each mover family is meaningful.
MOVER_TOLERANCE = {
    "BORIS": 1.0e-3,
    "RK2":   5.0e-3,
    "RK4":   1.0e-3,
    "RK6":   5.0e-4,
    "GC2":   5.0e-3,
    "GC4":   5.0e-3,
    "GC6":   5.0e-3,
}


class ShellRow(object):
    def __init__(self, alt_km, lon_deg, lat_deg, rc_num_gv, rc_code_stormer_gv, rel_err_code):
        self.alt_km = alt_km
        self.lon_deg = lon_deg
        self.lat_deg = lat_deg
        self.rc_num_gv = rc_num_gv
        self.rc_code_stormer_gv = rc_code_stormer_gv
        self.rel_err_code = rel_err_code


def stormer_vertical_gv(lat_deg, alt_km):
    r_re = (RE_KM + alt_km) / RE_KM
    c4 = math.cos(math.radians(lat_deg)) ** 4
    return STORMER_R0_GV * c4 / (r_re * r_re)


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
    movers = []
    allowed = set(DEFAULT_MOVERS + ("HYBRID",))
    for item in str(text).replace(";", ",").split(","):
        item = item.strip().upper()
        if not item:
            continue
        if item not in allowed:
            raise SystemExit("Unsupported mover %r; allowed values: %s" % (item, ",".join(sorted(allowed))))
        movers.append(item)
    if not movers:
        raise SystemExit("No movers were requested")
    return movers


def _parse_variables(line):
    return [n.strip().lower() for n in re.findall(r'"([^"]+)"', line)]


def _pick_column(variables, candidates, fallback):
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


def parse_tecplot_shell_output(path):
    rows = []
    current_alt = None
    variables = None
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


def close_to(a, b, tol=1.0e-6):
    return abs(a - b) <= tol


def render_input_template(template_path, output_path, mover, shell_alts_km, shell_res_deg,
                          cutoff_emin_mev, cutoff_emax_mev, dt_trace, adaptive_dt,
                          max_trace_time, max_trace_distance_re, target_lats):
    text = template_path.read_text()
    repl = {
        "__MOVER__": mover,
        "__SHELL_COUNT__": str(len(shell_alts_km)),
        "__SHELL_ALT_KM__": " ".join("%.10g" % x for x in shell_alts_km),
        "__SHELL_RES_DEG__": "%.10g" % shell_res_deg,
        "__CUTOFF_EMIN_MEV__": "%.10g" % cutoff_emin_mev,
        "__CUTOFF_EMAX_MEV__": "%.10g" % cutoff_emax_mev,
        "__DT_TRACE__": "%.10g" % dt_trace,
        "__ADAPTIVE_DT__": adaptive_dt,
        "__MAX_TRACE_TIME__": "%.10g" % max_trace_time,
        "__MAX_TRACE_DISTANCE_RE__": "%.10g" % max_trace_distance_re,
        "__TARGET_LATS_DEG__": ",".join("%.10g" % x for x in target_lats),
    }
    for k, v in repl.items():
        text = text.replace(k, v)
    # Replace shell count after shell alt substitution; older template versions may not have this placeholder.
    text = text.replace("SHELL_COUNT            1", "SHELL_COUNT            %d" % len(shell_alts_km))
    output_path.write_text(text)


def write_reference_csv(path, movers, full_alt_km, gc_alt_km, target_lats):
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "mover", "mover_family", "alt_km", "lat_deg", "Rc_stormer_GV",
            "required_for_pass", "rel_tol", "note",
        ])
        for mover in movers:
            family = "GC" if mover.startswith("GC") else ("HYBRID" if mover == "HYBRID" else "FULL_ORBIT")
            for alt in (full_alt_km, gc_alt_km):
                if family == "GC":
                    required = abs(alt - gc_alt_km) < 1.0e-6
                    note = "guiding-center mover is enforced only on the high-altitude shell"
                else:
                    required = True
                    note = "full-orbit mover checked on both shells"
                tol = MOVER_TOLERANCE.get(mover, 5.0e-3)
                for lat in target_lats:
                    w.writerow([
                        mover, family, "%.1f" % alt, "%.1f" % lat,
                        "%.9e" % stormer_vertical_gv(lat, alt),
                        "T" if required else "F", "%.3e" % tol, note,
                    ])


def find_output_file(workdir, user_output):
    if user_output:
        p = Path(user_output)
        if not p.is_absolute():
            p = workdir / p
        return p
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


def run_command(cmd, cwd, log_path):
    with log_path.open("w") as log:
        log.write("Command:\n  " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
        return proc.wait()


def summarize_mover(mover, rows, target_lats, checked_alts, gc_alt_km):
    family = "GC" if mover.startswith("GC") else ("HYBRID" if mover == "HYBRID" else "FULL_ORBIT")
    tol = MOVER_TOLERANCE.get(mover, 5.0e-3)
    summary = []
    messages = []
    passed = True

    for alt in checked_alts:
        for lat in target_lats:
            selected = [r for r in rows if close_to(r.alt_km, alt) and close_to(r.lat_deg, lat)]
            if not selected:
                required = not (family == "GC" and not close_to(alt, gc_alt_km))
                if required:
                    passed = False
                    messages.append("%s missing rows at alt=%g km lat=%g deg" % (mover, alt, lat))
                continue
            ref = stormer_vertical_gv(lat, alt)
            vals = [r.rc_num_gv for r in selected]
            mean_rc = sum(vals) / len(vals)
            min_rc = min(vals)
            max_rc = max(vals)
            rel_mean = (mean_rc - ref) / ref if ref > 0.0 else 0.0
            rel_max_abs = max(abs((x - ref) / ref) for x in vals) if ref > 0.0 else 0.0
            lon_spread = (max_rc - min_rc) / ref if ref > 0.0 else 0.0
            required_for_pass = not (family == "GC" and not close_to(alt, gc_alt_km))
            status = "CHECKED" if required_for_pass else "DIAGNOSTIC_ONLY"
            if required_for_pass and rel_max_abs > tol:
                passed = False
                messages.append(
                    "%s outside tolerance at alt=%g km lat=%g deg: max_abs_rel=%.3e tol=%.3e" %
                    (mover, alt, lat, rel_max_abs, tol)
                )
            summary.append({
                "mover": mover,
                "mover_family": family,
                "alt_km": alt,
                "lat_deg": lat,
                "n_lon": float(len(selected)),
                "Rc_reference_GV": ref,
                "Rc_num_mean_GV": mean_rc,
                "Rc_num_min_GV": min_rc,
                "Rc_num_max_GV": max_rc,
                "rel_err_mean": rel_mean,
                "rel_err_max_abs": rel_max_abs,
                "lon_spread_over_ref": lon_spread,
                "rel_tol": tol,
                "required_for_pass": required_for_pass,
                "status": status,
            })
    return summary, passed, messages


def write_summary_csv(summary, path):
    if not summary:
        return
    keys = list(summary[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in summary:
            w.writerow(row)


def make_plot(summary, path):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    if not summary:
        return

    alts = sorted(set(row["alt_km"] for row in summary if row["required_for_pass"]))
    for alt in alts:
        plt.figure(figsize=(7.0, 4.5))
        for mover in sorted(set(row["mover"] for row in summary)):
            items = [r for r in summary if r["mover"] == mover and close_to(r["alt_km"], alt) and r["required_for_pass"]]
            if not items:
                continue
            items.sort(key=lambda x: x["lat_deg"])
            plt.plot([x["lat_deg"] for x in items], [x["rel_err_mean"] for x in items], marker="o", label=mover)
        plt.axhline(0.0, linestyle="--", linewidth=1.0)
        plt.xlabel("Magnetic latitude (deg)")
        plt.ylabel("Mean relative error vs Størmer")
        plt.title("C12: mover residuals, alt=%g km" % alt)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=8)
        plt.tight_layout()
        out = path.with_name(path.stem + "_alt%gkm" % alt + path.suffix)
        plt.savefig(out, dpi=160)
        plt.close()


def cross_mover_checks(summary, alt_km, target_lats, mover_group, tol):
    checks = []
    passed = True
    for lat in target_lats:
        rows = [r for r in summary if close_to(r["alt_km"], alt_km) and close_to(r["lat_deg"], lat) and r["mover"] in mover_group]
        if len(rows) < 2:
            continue
        values = [r["Rc_num_mean_GV"] for r in rows]
        ref = stormer_vertical_gv(lat, alt_km)
        spread = (max(values) - min(values)) / ref if ref > 0.0 else 0.0
        ok = spread <= tol
        if not ok:
            passed = False
        checks.append({
            "alt_km": alt_km,
            "lat_deg": lat,
            "movers": ",".join(mover_group),
            "spread_over_ref": spread,
            "rel_tol": tol,
            "passed": ok,
        })
    return checks, passed


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        C12: Particle-mover cross-validation against analytical Størmer cutoff.

        The test is gridless by design because the selectable mover layer is the
        shared gridless particle-mover implementation.  The runner launches AMPS
        once per requested mover, parses the shell cutoff output, and compares
        each mover with the analytical Størmer vertical cutoff.
        """),
        epilog=textwrap.dedent("""
        Examples:
          python srcEarth/test/C12/run_C12.py -np 4 -nt 16
          python srcEarth/test/C12/run_C12.py --movers BORIS,RK4,RK6 -np 2 -nt 8
          python srcEarth/test/C12/run_C12.py --movers GC4 --gc-alt 25000 --skip-run
          python srcEarth/test/C12/run_C12.py --enforce-cross-mover
        """)
    )
    parser.add_argument("-np", type=int, default=4, help="number of MPI ranks passed to mpirun; default: 4")
    parser.add_argument("-nt", type=int, default=16, help="thread count passed through density/thread CLI settings; default: 16")
    parser.add_argument("--movers", default=",".join(DEFAULT_MOVERS), help="comma-separated mover list; default: BORIS,RK2,RK4,RK6,GC2,GC4,GC6")
    parser.add_argument("--lats", default="-60,-40,-20,0,20,40,60", help="comma-separated latitudes to check; default: -60,-40,-20,0,20,40,60. Use --lats=-60,... if the first value is negative.")
    parser.add_argument("--full-alt", type=float, default=9000.0, help="standard shell altitude for full-orbit movers, km; default: 9000")
    parser.add_argument("--gc-alt", type=float, default=25000.0, help="high shell altitude where GC movers are enforced, km; default: 25000")
    parser.add_argument("--shell-res", type=float, default=20.0, help="SHELL_RES_DEG used by AMPS; default: 20")
    parser.add_argument("--cutoff-emin", type=float, default=0.1, help="CUTOFF_EMIN in MeV/n; default: 0.1 so high-altitude 60-deg cutoffs are above Rmin")
    parser.add_argument("--cutoff-emax", type=float, default=20000.0, help="CUTOFF_EMAX in MeV/n; default: 20000")
    parser.add_argument("--dt-trace", type=float, default=1.0, help="DT_TRACE, seconds; default: 1")
    parser.add_argument("--adaptive-dt", choices=["T", "F", "t", "f", "TRUE", "FALSE", "true", "false", "1", "0"], default="T", help="ADAPTIVE_DT setting written to input; default: T")
    parser.add_argument("--max-trace-time", type=float, default=600.0, help="maximum trajectory time, seconds; default: 600")
    parser.add_argument("--max-trace-distance", type=float, default=400.0, help="maximum cumulative trace distance, Re; default: 400")
    parser.add_argument("--scheduler", default="DYNAMIC", choices=["DYNAMIC", "BLOCK_CYCLIC", "STATIC"], help="gridless MPI scheduler; default: DYNAMIC")
    parser.add_argument("--dynamic-chunk", type=int, default=0, help="gridless dynamic chunk size; 0 means AMPS auto heuristic; default: 0")
    parser.add_argument("--no-cutoff-search-cli", action="store_true", help="do not pass -cutoff-search/-cutoff-upper-scan-n; useful for older CLI checkouts")
    parser.add_argument("--enforce-cross-mover", action="store_true", help="make BORIS/RK4/RK6/GC4 high-altitude cross-mover spread a hard failure")
    parser.add_argument("--cross-mover-tol", type=float, default=5.0e-3, help="relative tolerance for optional cross-mover spread; default: 5e-3")
    parser.add_argument("--amps", default="./amps", help="AMPS executable path relative to launch directory; default: ./amps")
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher executable; default: mpirun")
    parser.add_argument("--workdir", default="test_output/C12_gridless", help="base work directory; default: test_output/C12_gridless")
    parser.add_argument("--output-file", default=None, help="explicit output file name inside each mover workdir; default: auto-detect")
    parser.add_argument("--skip-run", action="store_true", help="do not run AMPS; analyze existing mover workdirs")
    parser.add_argument("--keep", action="store_true", help="keep existing work directory instead of replacing it")
    parser.add_argument("--dry-run", action="store_true", help="prepare input files and print commands, but do not run AMPS or analyze output")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.np < 1:
        raise SystemExit("-np must be >= 1")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")
    if args.full_alt <= 0.0 or args.gc_alt <= 0.0:
        raise SystemExit("--full-alt and --gc-alt must be positive")

    movers = parse_mover_list(args.movers)
    target_lats = parse_float_list(args.lats, "--lats")
    adaptive_dt = args.adaptive_dt.upper()
    if adaptive_dt == "TRUE" or adaptive_dt == "1":
        adaptive_dt = "T"
    elif adaptive_dt == "FALSE" or adaptive_dt == "0":
        adaptive_dt = "F"

    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    template_input = script_dir / "AMPS_PARAM_C12_gridless.in"
    reference_csv_src = script_dir / "reference_C12_stormer_movers.csv"

    base_workdir = (launch_dir / args.workdir).resolve()
    if not args.skip_run:
        if base_workdir.exists() and not args.keep:
            shutil.rmtree(base_workdir)
        base_workdir.mkdir(parents=True, exist_ok=True)
    else:
        if not base_workdir.exists():
            raise SystemExit("--skip-run requested but workdir does not exist: %s" % base_workdir)

    # Always write/update the reference table in the source directory and the run directory.
    write_reference_csv(reference_csv_src, movers, args.full_alt, args.gc_alt, target_lats)
    if not args.skip_run:
        shutil.copy2(reference_csv_src, base_workdir / "reference_C12_stormer_movers.csv")

    all_summary = []
    all_messages = []
    run_records = []
    overall_passed = True

    amps_path = Path(args.amps)
    if not amps_path.is_absolute():
        amps_path = (launch_dir / amps_path).resolve()

    checked_alts = [args.full_alt]
    if not close_to(args.gc_alt, args.full_alt):
        checked_alts.append(args.gc_alt)

    for mover in movers:
        mover_dir = base_workdir / mover.lower()
        if not args.skip_run:
            mover_dir.mkdir(parents=True, exist_ok=True)
            render_input_template(
                template_input,
                mover_dir / "AMPS_PARAM_C12.in",
                mover,
                checked_alts,
                args.shell_res,
                args.cutoff_emin,
                args.cutoff_emax,
                args.dt_trace,
                adaptive_dt,
                args.max_trace_time,
                args.max_trace_distance,
                target_lats,
            )

        cmd = [
            args.mpirun,
            "-np", str(args.np),
            str(amps_path),
            "-mode", "gridless",
            "-i", "AMPS_PARAM_C12.in",
            "-mover", mover,
            "-gridless-mpi-scheduler", args.scheduler,
            "-gridless-mpi-dynamic-chunk", str(args.dynamic_chunk),
            "-density-parallel", "THREADS",
            "-density-threads", str(args.nt),
        ]
        if not args.no_cutoff_search_cli:
            cmd += ["-cutoff-search", "UPPER_SCAN", "-cutoff-upper-scan-n", "80"]

        log_file = mover_dir / ("C12_%s_amps.log" % mover.lower())
        run_records.append({"mover": mover, "workdir": str(mover_dir), "command": cmd, "log_file": str(log_file)})

        if args.dry_run:
            print("[%s] %s" % (mover, " ".join(cmd)))
            continue

        if not args.skip_run:
            print("\nRunning C12 mover %s in %s" % (mover, mover_dir))
            print(" ".join(cmd))
            rc = run_command(cmd, cwd=mover_dir, log_path=log_file)
            if rc != 0:
                overall_passed = False
                all_messages.append("AMPS run for mover %s failed with exit code %d; see %s" % (mover, rc, log_file))
                continue

        output_file = find_output_file(mover_dir, args.output_file)
        if not output_file.exists():
            overall_passed = False
            all_messages.append("Missing output for mover %s: %s" % (mover, output_file))
            continue

        try:
            rows = parse_tecplot_shell_output(output_file)
        except Exception as exc:
            overall_passed = False
            all_messages.append("Could not parse output for mover %s: %s" % (mover, exc))
            continue

        summary, passed, messages = summarize_mover(mover, rows, target_lats, checked_alts, args.gc_alt)
        all_summary.extend(summary)
        all_messages.extend(messages)
        if not passed:
            overall_passed = False

        write_summary_csv(summary, mover_dir / ("C12_%s_summary.csv" % mover.lower()))

    if args.dry_run:
        print("\nDry run complete.  Prepared commands for %d movers." % len(movers))
        return 0

    cross_checks, cross_passed = cross_mover_checks(all_summary, args.gc_alt, target_lats, [m for m in ("BORIS", "RK4", "RK6", "GC4") if m in movers], args.cross_mover_tol)
    if args.enforce_cross_mover and not cross_passed:
        overall_passed = False
        all_messages.append("High-altitude cross-mover spread exceeded tolerance; see cross_mover_checks in result JSON")

    summary_csv = base_workdir / "C12_summary.csv"
    result_json = base_workdir / "C12_result.json"
    plot_png = base_workdir / "C12_mover_residuals.png"
    write_summary_csv(all_summary, summary_csv)
    make_plot(all_summary, plot_png)

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "passed": overall_passed,
        "messages": all_messages,
        "np": args.np,
        "nt": args.nt,
        "movers": movers,
        "full_alt_km": args.full_alt,
        "gc_alt_km": args.gc_alt,
        "target_lats_deg": target_lats,
        "cutoff_emin_mev": args.cutoff_emin,
        "cutoff_emax_mev": args.cutoff_emax,
        "dt_trace": args.dt_trace,
        "adaptive_dt": adaptive_dt,
        "max_trace_time": args.max_trace_time,
        "max_trace_distance_re": args.max_trace_distance,
        "scheduler": args.scheduler,
        "dynamic_chunk": args.dynamic_chunk,
        "stormer_R0_GV": STORMER_R0_GV,
        "summary_csv": str(summary_csv),
        "reference_csv": str(base_workdir / "reference_C12_stormer_movers.csv"),
        "run_records": run_records,
        "cross_mover_checks": cross_checks,
        "cross_mover_enforced": args.enforce_cross_mover,
        "summary": all_summary,
    }
    result_json.write_text(json.dumps(result, indent=2))

    print("\nC12 summary")
    print("===========")
    print("movers=%s np=%d nt=%d scheduler=%s" % (",".join(movers), args.np, args.nt, args.scheduler))
    print("full_alt=%g km gc_alt=%g km lats=%s" % (args.full_alt, args.gc_alt, ",".join("%g" % x for x in target_lats)))
    for row in all_summary:
        if not row["required_for_pass"]:
            continue
        print(
            "%5s alt=%8.1f km lat=%6.1f deg Rc=%10.5e ref=%10.5e rel=% .3e tol=% .1e" %
            (row["mover"], row["alt_km"], row["lat_deg"], row["Rc_num_mean_GV"],
             row["Rc_reference_GV"], row["rel_err_mean"], row["rel_tol"])
        )
    if cross_checks:
        print("\nCross-mover diagnostic at alt=%g km:" % args.gc_alt)
        for c in cross_checks:
            print("  lat=%6.1f spread/ref=% .3e tol=% .1e %s" % (c["lat_deg"], c["spread_over_ref"], c["rel_tol"], "PASS" if c["passed"] else "FAIL"))
        if not args.enforce_cross_mover:
            print("  Cross-mover spread is diagnostic only; use --enforce-cross-mover to make it a hard gate.")
    if all_messages:
        print("\nMessages:")
        for msg in all_messages:
            print("  - " + msg)
    print("\nWrote: %s" % summary_csv)
    print("Wrote: %s" % result_json)
    print("RESULT: %s" % ("PASS" if overall_passed else "FAIL"))
    return 0 if overall_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
