#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/seen_after_gap.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/seen_after_gap.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3316a98123ae36539beba2662a9776502a559ca2af9363ee6e668b88f8728d13
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: marks the facets seen by NG or RG at any point after the NG gap (from the record store); the mask for
#   refine_mesh.py --seed-mask.
# USED FOR: the 20M refinement (jobs/refine_20m.sbatch; 1,466,454 facets). GPU node (reads the store).
# RUN: .venv/bin/python <this> --run <run dir> --after 2016-09-30T10:33:10 --out <seen.npz>
# ==== end of header ====
"""Facets seen by NG or RG after the gap (user request 2026-11-04: refine the 20M mesh over the entire surface seen
after the gap by either gauge). Reads the record store (face, kn, kphi, kv, cos per point, field-major), the points
with UTC >= --after (default 2016-09-30T10:33:10, the end of the gap / grain t_end).

    .venv/bin/python <this> --run <runner dir> --out <file.npz> [--after 2016-09-30T10:33:10]
Product: <out> with seen (bool per facet), ng_weight, rg_weight (sum over the post-gap points of the
stored kernel fields |kn| and |kv|; a relative measure only), point_ids; and a JSON summary on stdout (counts, areas, and the
number of facets that carry 90/99/99.9 % of the summed NG and RG kernel weight).
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--after", default="2016-09-30T10:33:10")
    a = ap.parse_args()
    if a.out.exists():
        sys.exit("%s exists" % a.out)
    st = json.loads((a.run / "store" / "store_manifest.json").read_text())
    nf = st["nfaces"]
    # UTC of every point from the base's COPS_model.dat
    utc = {}
    names = None
    for line in (a.run / "base" / "COPS_model.dat").read_text().splitlines():
        if line.startswith("#iPoint"):
            names = line[1:].split(); continue
        if line.startswith("#") or not line.strip() or names is None:
            continue
        v = dict(zip(names, line.split())); utc[int(v["iPoint"])] = v["UTC"]
    def key(s):   # COPS_model.dat writes e.g. 2016-09-30T0:00:04 (hour not zero-padded)
        d, t = s.split("T"); h, m, sec = t.split(":"); return (d, int(h), int(m), float(sec))
    pts = [p for p in st["points"] if key(utc[p["point"]]) >= key(a.after)]
    seen = np.zeros(nf, bool); wn = np.zeros(nf); wv = np.zeros(nf)
    for p in pts:
        n = p["records"]
        buf = np.fromfile(a.run / "store" / p["store_file"], dtype=np.uint8)
        face = np.frombuffer(buf, "<i4", n, 0)
        kn = np.frombuffer(buf, "<f4", n, 4 * n)
        kv = np.frombuffer(buf, "<f4", n, 12 * n)
        seen[face] = True
        np.add.at(wn, face, np.abs(kn)); np.add.at(wv, face, np.abs(kv))
    el = a.run / "base" / "elements.bin"
    sys.path[:0] = ["/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/tools"]
    from cops_cuda import reference
    area = reference.load("fit_cops").read_elements(str(el))["area"]
    def ncarry(w):
        s = np.sort(w[w > 0])[::-1]; c = np.cumsum(s) / s.sum()
        return {str(q): int(np.searchsorted(c, q) + 1) for q in (0.9, 0.99, 0.999)}
    np.savez(a.out, seen=seen, ng_weight=wn, rg_weight=wv, point_ids=np.array([p["point"] for p in pts]))
    out = {"run": str(a.run), "after_utc": a.after, "points": len(pts),
           "first_point": pts[0]["point"], "last_point": pts[-1]["point"],
           "facets_seen": int(seen.sum()), "facets_total": int(nf), "area_seen_m2": float(area[seen].sum()),
           "area_seen_median_m2": float(np.median(area[seen])),
           "facets_carrying_NG_weight": ncarry(wn), "facets_carrying_RG_weight": ncarry(wv)}
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
