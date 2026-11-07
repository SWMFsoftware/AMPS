#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/late_time_grid.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/late_time_grid.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 244f324f9d0adb244f744b922def240043ea3e2ceed29d7aa1d3b1694eb2559d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Writes an observation file (LiouvilleCOPS format: ObservationTime, NudeGaugeReferenceData,
#   RamGaugeReferenceData) with a uniform time grid (reference 0 = gap) plus the original observations from a given
#   time, for forward runs of the late window on a fine grid (COPS time response; user 2026-11-07).
#   --allow-missing-l2 (added 2026-11-07 for the fitting base of job 122291): keep observations without an L2 STOP time.
# USED BY: jobs/late_response.sbatch (step 1), jobs/prepare_12m_lategrid.sbatch.
# RUN: .venv/bin/python late_time_grid.py --observations <AMPS>/srcCG/RosinaMeasurements.h --hk <COPS_HK_20160930.txt>
#   --start 2016-09-30T10:32:46.282 --stop 2016-09-30T10:39:26.282 --step 1 --keep-from <UTC> [--allow-missing-l2]
#   --out late_grid.h --table late_grid_points.txt
# ==== end of header ====
"""Observation file for a forward run of the late window on a fine time grid (user 2026-11-07: compare the model with the
COPS measurement's time average, windows 10 s and 10 s + DPU 5-readout average; gauge-offset sensitivity check).

    .venv/bin/python late_time_grid.py --observations <RosinaMeasurements.h> --hk <COPS_HK_20160930.txt> \
        --start 2016-09-30T10:32:46.282 --stop 2016-09-30T10:39:26.282 --step 1 --keep-from 2016-09-30T10:33:10 \
        --out <late_grid.h> [--table <late_grid_points.txt>]

Writes a file in the format LiouvilleCOPS reads (arrays ObservationTime, NudeGaugeReferenceData, RamGaugeReferenceData;
src/Observations.cpp) with
  * the uniform grid start, start + step, ..., stop (reference values 0 = data gap: computed, never fitted), and
  * the original observations from --keep-from on, with their times and reference values unchanged (to check the new
    base against the fitting base and to give the model at the time stamps as the fits used it),
merged in time order. --table lists every point: index, UTC, kind (grid/obs), the original index for obs points, and
the L2 interval (START/STOP from the HK table) of obs points. The grid must lie inside the HK table (T_RG is
interpolated from it; LiouvilleCOPS stops outside). Standard library only; read only except --out/--table.

2026-11-07 --allow-missing-l2 (for --keep-from 2016-09-30T00:00:00, all 1150 observations; prepare_12m_lategrid.sbatch):
an observation without an HK interval ending within 1 s (no L2 bg file in the HK table: 2016-09-30T9:13:04 and 9:14:04,
NG 0, RG > 0) is kept unchanged and gets L2 "-" in --table, and is listed; without the option such an observation stops
the tool as before (default behaviour unchanged).
"""

import argparse
import re
from datetime import datetime, timedelta


def parse_times(text, name):
    body = re.search(r"%s\s*\[[^\]]*\]\s*(?:\[[^\]]*\])?\s*=\s*\{(.*?)\}" % name, text, re.S).group(1)
    return re.findall(r'"([^"]+)"', body)


def parse_numbers(text, name):
    """the values as written (strings), so kept observations are copied digit for digit"""
    body = re.search(r"%s\s*\[[^\]]*\]\s*=\s*\{(.*?)\}" % name, text, re.S).group(1)
    return [v.strip() for v in body.replace("\n", " ").split(",") if v.strip()]


def utc(s):
    d, t = s.split("T")
    h, m, sec = t.split(":")
    return datetime.fromisoformat(d) + timedelta(hours=int(h), minutes=int(m), seconds=float(sec))


def fmt(t):
    return t.strftime("%Y-%m-%dT%H:%M:%S.") + "%03d" % round(t.microsecond / 1000.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observations", required=True)
    ap.add_argument("--hk", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--stop", required=True)
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--keep-from", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--table")
    ap.add_argument("--allow-missing-l2", action="store_true",
                    help="observations without an HK interval ending within 1 s get L2 '-' in --table (2026-11-07)")
    a = ap.parse_args()
    text = re.sub(r"//[^\n]*", "", open(a.observations).read())
    times = parse_times(text, "ObservationTime")
    ng = parse_numbers(text, "NudeGaugeReferenceData")
    rg = parse_numbers(text, "RamGaugeReferenceData")
    assert len(times) == len(ng) == len(rg), (len(times), len(ng), len(rg))
    hk = []
    for line in open(a.hk):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        hk.append((utc(f[1]), utc(f[2])))
    t0, t1, keep = utc(a.start), utc(a.stop), utc(a.keep_from)
    if t0 < hk[0][0] or t1 > hk[-1][1]:
        raise SystemExit("grid %s .. %s outside the HK table %s .. %s" % (t0, t1, hk[0][0], hk[-1][1]))
    pts = []
    n = int(round((t1 - t0).total_seconds() / a.step))
    for k in range(n + 1):
        t = t0 + timedelta(seconds=k * a.step)
        pts.append((t, "grid", -1, fmt(t), "0", "0", None))
    missing = []
    for i, s in enumerate(times):
        t = utc(s)
        if t >= keep:
            iv = min(hk, key=lambda r: abs((r[1] - t).total_seconds()))
            if abs((iv[1] - t).total_seconds()) > 1.0:
                if not a.allow_missing_l2:
                    raise SystemExit("no HK interval ends within 1 s of %s" % s)
                missing.append("%d %s (NG %s, RG %s)" % (i, s, ng[i], rg[i]))
                iv = None
            pts.append((t, "obs", i, s, ng[i], rg[i], iv))
    pts.sort(key=lambda p: (p[0], p[1]))
    for p, q in zip(pts, pts[1:]):
        if p[0] == q[0]:
            raise SystemExit("duplicate time %s" % p[3])
    with open(a.out, "w") as f:
        f.write("// late-window time grid made by model_ext/late_time_grid.py from %s\n" % a.observations)
        f.write("// grid %s .. %s step %g s (reference 0 = gap) + original observations from %s\n" % (a.start, a.stop, a.step, a.keep_from))
        f.write("namespace RosinaSample {\n  const int nPoints=%d;\n" % len(pts))
        f.write("  const char ObservationTime[nPoints][100]={\n" + ",\n".join('"%s"' % p[3] for p in pts) + "};\n")
        f.write("  const double  NudeGaugeReferenceData[nPoints] ={\n" + ",\n".join(p[4] for p in pts) + "};\n")
        f.write("  const double  RamGaugeReferenceData[nPoints] ={\n" + ",\n".join(p[5] for p in pts) + "};\n}\n")
    if a.table:
        with open(a.table, "w") as f:
            f.write("# index UTC kind original_index L2_start L2_stop\n")
            for j, p in enumerate(pts):
                iv = p[6]
                f.write("%d %s %s %d %s %s\n" % (j, p[3], p[1], p[2], fmt(iv[0]) if iv else "-", fmt(iv[1]) if iv else "-"))
    print("wrote %s: %d points (%d grid, %d observations)" % (a.out, len(pts), n + 1, len(pts) - n - 1))
    if missing:
        print("%d observations without an L2 interval in the HK table (L2 '-'): %s" % (len(missing), "; ".join(missing)))


if __name__ == "__main__":
    main()
