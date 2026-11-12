#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-10, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/analysis/shadow_zero_activity_12m/shadow_classes.py (on Voyager: /data001/highresolu_vt/cops_output/analysis/shadow_zero_activity_12m/shadow_classes.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 780277ee437df88d083c2be1a7f597248b66b8e67c47d9cffa3736f7c5298eaa
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis tool. Illumination (sunlit / cast shadow / night) of each seen facet over the records of the points a
#   fit uses, by count, area and NG/RG kernel weight, split into zero-activity and active facets of a clipped fit.
#   User 2026-11-10: 'which part of the seen surface that has zero activity is in the shadow?'. CPU job 123614 (D2).
# RUN: .venv/bin/python shadow_classes.py --base <base> --points <late_grid_points.txt> --fit <fit dir> --out <new dir>
#   --workers 8   (CPU job: reads every kernel record of the base, ~170 GB; ~8 GB host memory with 8 workers)
# ==== end of header ====
"""Illumination of the seen surface with and without activity in a clipped fit (user 2026-11-10: "which part of the
seen surface that has zero activity is in the shadow?").

The forward model has no single shadow state: every kernel record (point i, facet e) stores cos(SZA) of the facet at
that observation time and the cast-shadow flag (cops_cuda/formats.py KERNEL_RECORD_DTYPE). For every facet this script
sums, over the records of the points the fit uses, the record count and the gauge weights in three illumination states:
  lit     cos > 0 and not cast-shadowed  (day-side production q_min + (q_max - q_min) cos)
  shadow  cos > 0 and cast-shadowed      (night-side production q_min, srcCG rule)
  night   cos <= 0                       (night-side production q_min)
weights: NG = Kn (density kernel), RG = KV + |KPhi| (flux kernels; KPhi can carry the ram sign), and the plain count.
Points (D2's objective, as the fit's own seen set): observation rows before the late window (2016-09-30T10:33:10) with
useNG or useRG > 0, plus all model-only 1-s grid rows (they carry the W10x5 averages of the late observations).
Activity: the fit's activity_facets.dat (columns iFacet x y z area seen legacy residual total clipped).

    shadow_classes.py --base <base> --points <late_grid_points.txt> --fit <fit dir> --out <dir> [--workers N]
Writes <out>/facet_illumination.npz (per facet sums) and <out>/summary.txt. Read only; CPU.
"""
import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import formats  # noqa: E402

LATE = "2016-09-30T10:33:10"
NFACES = None


def one(args):
    path, point, nfaces = args
    _geo, rec, _sha = formats.read_kernel(path, point, nfaces)
    f = rec["face"]
    lit = (rec["cos"] > 0) & (rec["shadow"] == 0)
    sh = (rec["cos"] > 0) & (rec["shadow"] == 1)
    nt = rec["cos"] <= 0
    wn = rec["Kn"]; wr = rec["KV"] + np.abs(rec["KPhi"])
    out = {}
    for name, m in (("lit", lit), ("shadow", sh), ("night", nt)):
        out["n_" + name] = np.bincount(f[m], minlength=nfaces).astype(np.float64)
        out["ng_" + name] = np.bincount(f[m], weights=wn[m], minlength=nfaces)
        out["rg_" + name] = np.bincount(f[m], weights=wr[m], minlength=nfaces)
    return out


def chunk(jobs):
    """Process a list of points in one worker and return ONE accumulated set of per-facet sums (job 123612 ran out of
    memory when every point's ~0.9 GB of sums was sent back to the parent)."""
    acc = None
    for job in jobs:
        r = one(job)
        if acc is None:
            acc = r
        else:
            for key in acc:
                acc[key] += r[key]
        del r
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--points", required=True, type=Path)
    ap.add_argument("--fit", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    # points: from the fit's COPS_model.dat (use flags, UTC) and the point table (kind grid)
    from cops_cuda import figures
    t = figures.read_table(a.fit / "COPS_model.dat")
    kinds = {}
    for line in a.points.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        v = line.split()
        kinds[int(v[0])] = v
    head = [l for l in a.points.read_text().splitlines() if l.startswith("#")][-1][1:].split()
    ci = {n: k for k, n in enumerate(head)}
    kind_col = ci.get("kind")
    pts = []
    late = figures.parse_utc(LATE)      # UTC strings have hours without a leading zero: compare times, not text
    for k, ip in enumerate(t["iPoint"].astype(int)):
        when = t["time"][k]
        is_grid = kind_col is not None and kinds.get(ip, [None] * len(head))[kind_col] == "grid"
        if is_grid:
            pts.append(ip)
        elif when < late and (t["useNG"][k] > 0 or t["useRG"][k] > 0):
            pts.append(ip)
    nfaces = int(formats.read_element_header(a.base / "elements.bin")[0])
    jobs = [(a.base / "kernels" / formats.kernel_file_name(p), p, nfaces) for p in pts]
    print("points used: %d (%d grid); facets %d" % (len(pts), sum(1 for p in pts if kinds.get(p, [None] * len(head))[kind_col] == "grid"), nfaces), flush=True)
    acc = None
    parts = [jobs[k::a.workers] for k in range(a.workers)]
    with mp.Pool(a.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(chunk, parts)):
            if acc is None:
                acc = r
            else:
                for key in acc:
                    acc[key] += r[key]
            del r
            print("  worker part %d / %d done" % (k + 1, len(parts)), flush=True)
    # the fit's activity and area
    cols = np.loadtxt(a.fit / "activity_facets.dat", usecols=(4, 5, 9), comments="#")
    area, seen_fit, act = cols[:, 0], cols[:, 1] > 0, cols[:, 2]
    np.savez_compressed(a.out / "facet_illumination.npz", area=area, seen_fit=seen_fit, activity=act, points=np.array(pts), **acc)
    n_tot = acc["n_lit"] + acc["n_shadow"] + acc["n_night"]
    seen = n_tot > 0
    lines = ["# illumination of the seen surface, fit %s, base %s" % (a.fit, a.base),
             "# points used %d; facets seen from them %d; the fit's own seen flag %d; both %d"
             % (len(pts), seen.sum(), seen_fit.sum(), (seen & seen_fit).sum())]
    zero = seen_fit & (act <= 0); on = seen_fit & (act > 0)
    for lab, m in (("zero activity", zero), ("active", on)):
        lines.append("## seen facets with %s: %d facets, %.4e m2" % (lab, m.sum(), area[m].sum()))
        nlit, nsh, nnt = acc["n_lit"][m], acc["n_shadow"][m], acc["n_night"][m]
        tot = nlit + nsh + nnt
        never = (nlit == 0) & (tot > 0); always = (nsh + nnt == 0) & (tot > 0); mixed = ~never & ~always & (tot > 0)
        w = area[m]
        for nm, mm in (("never sunlit while seen (shadow or night at every record)", never),
                       ("  of these: cast shadow at some record", never & (nsh > 0)),
                       ("  of these: night side only", never & (nsh == 0)),
                       ("always sunlit while seen", always), ("sunlit at some records, dark at others", mixed)):
            lines.append("  %-58s %6.2f %% of the facets, %6.2f %% of the area" % (nm, 100 * mm.sum() / m.sum(), 100 * w[mm].sum() / w.sum()))
        for g in ("n", "ng", "rg"):
            L, S, N = acc[g + "_lit"][m].sum(), acc[g + "_shadow"][m].sum(), acc[g + "_night"][m].sum()
            T = L + S + N
            lab2 = {"n": "records (equal weight)", "ng": "NG weight (Kn)", "rg": "RG weight (KV+|KPhi|)"}[g]
            lines.append("  share of %-24s lit %6.2f %%  cast shadow %6.2f %%  night %6.2f %%" % (lab2, 100 * L / T, 100 * S / T, 100 * N / T))
    (a.out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
