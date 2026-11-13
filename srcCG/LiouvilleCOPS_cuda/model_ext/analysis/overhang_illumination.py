#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/analysis/overhang_illumination_12m/overhang_illumination.py (on Voyager: /data001/highresolu_vt/cops_output/analysis/overhang_illumination_12m/overhang_illumination.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2426ddafaf636e56684f50a00a64f3324ccbafcf8d2146ac1194ae50e0edb3d4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis. Illumination (sunlit / night / cast shadow) of the seen facets per effective-slope bin at the map
#   times, overhangs (> 90 deg) vs flat facets (FINDINGS.md A1, 2026-11-12).
# RUN: analysis/overhang_illumination_12m/run.sbatch (copy: analysis/overhang_illumination.sbatch); job 124578.
# ==== end of header ====
"""Why do overhangs (effective slope > 90 deg) contribute so much close to the surface? (user question 2026-11-12:
"more than 90 means it is in the shadow and should not contribute much (correct?) but they are a significant part
close to the surface -- how that is possible?").
For the softplus 5-cycle final state: per post-gap observation (95) and gauge, the facet contributions of the
contribution job (source_areas_<G>.npz: faces_sorted, contrib_sorted_Pa, offsets; every facet with a non-zero
contribution) split by slope group (overhang [90,180] / 60-70 / 0-40 / other) and illumination (day / night / cast
shadow) of the facet at the map time nearest the observation (t11 = 10:33:10, t12 = 10:39:26; the rotation between
them is ~3 deg). Also the illumination of all seen overhang facets and their activity. CPU, read only.
    overhang_illumination.py <chain> <geometry dir> <out dir>
"""
import json, sys
from pathlib import Path
import numpy as np

ch, geo, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=False)
nf = np.load(ch / "c5/fit/final/nonneg_facets.npz")
binv, seen = nf["bin"].astype(int), nf["seen"]
A = np.load(ch / "c5/ft/final/finetune_state.npz")["A_final"]
area = np.fromfile(geo / "mesh_facegeom.f8").reshape(len(binv), -1)[:, -1] if False else None
group = np.full(len(binv), 3, np.int8)          # 0 overhang, 1 slope 60-70, 2 slope 0-40, 3 other
group[binv == 6] = 0; group[binv == 3] = 1; group[binv == 0] = 2
GN = ["overhang 90-180", "slope 60-70", "slope 0-40", "other (40-60, 70-90)"]
ill = {}
for t in ("t11", "t12"):
    cos = np.fromfile(geo / f"{t}.cos.f8"); sh = np.fromfile(geo / f"{t}.shadow.u1", dtype=np.uint8)
    s = np.full(len(cos), 1, np.int8)            # 1 night (cos <= 0)
    s[(cos > 0) & (sh == 0)] = 0                 # 0 day
    s[(cos > 0) & (sh != 0)] = 2                 # 2 cast shadow
    ill[t] = (s, cos)
IN = ["day", "night", "cast shadow"]
t_mid = 528503658.18234831 + (528504034.18234825 - 528503658.18234831) / 2
lines = ["# overhang illumination: softplus 5-cycle final state, post-gap observations; facet contributions of the contribution job",
         "# illumination at the nearest map time (t11 10:33:10 / t12 10:39:26)"]
summary = {}
for g in ("NG", "RG"):
    z = np.load(ch / f"contributions/source_areas_{g}.npz")
    off, fs, cs, et, utc = z["offsets"], z["faces_sorted"], z["contrib_sorted_Pa"], z["et"], z["utc"]
    tab = np.zeros((len(et), 4, 3))
    for k in range(len(et)):
        f = fs[off[k]:off[k + 1]]; c = cs[off[k]:off[k + 1]]
        s = ill["t11" if et[k] < t_mid else "t12"][0][f]
        for gi in range(4):
            for si in range(3):
                tab[k, gi, si] = c[(group[f] == gi) & (s == si)].sum()
    tot = tab.sum(axis=(1, 2))
    lines.append("## %s: share of the facet contribution per observation (first, 25th, 50th, 75th, last observation, then all summed)" % g)
    for k in (0, 23, 47, 71, 94):
        row = tab[k] / tot[k]
        lines.append("  %s  " % utc[k] + "  ".join("%s: day %.3f night %.3f shadow %.3f" % (GN[gi], *row[gi]) for gi in range(3)))
    allrow = tab.sum(0) / tot.sum()
    lines.append("  ALL      " + "  ".join("%s: day %.3f night %.3f shadow %.3f" % (GN[gi], *allrow[gi]) for gi in range(4)))
    oh = tab[:, 0, :].sum(1) / tot
    oh_dark = (tab[:, 0, 1] + tab[:, 0, 2]) / np.maximum(tab[:, 0, :].sum(1), 1e-300)
    lines.append("  overhang share of the contribution: first %.3f, max %.3f (at %s), last %.3f; of the overhang part, night + shadow: median %.3f, min %.3f, max %.3f"
                 % (oh[0], oh.max(), utc[oh.argmax()], oh[-1], np.median(oh_dark), oh_dark.min(), oh_dark.max()))
    summary[g] = {"overhang_share": oh.tolist(), "overhang_dark_fraction": oh_dark.tolist(), "utc": utc.tolist(),
                  "table": tab.tolist()}
for t in ("t11", "t12"):
    s, cos = ill[t]
    m = seen & (group == 0)
    lines.append("## seen overhang facets at %s: %d; day %.3f, night %.3f, cast shadow %.3f (fraction of facets); with A > 0: %d (day %.3f, night %.3f, shadow %.3f)"
                 % (t, m.sum(), *(np.mean(s[m] == i) for i in range(3)), (m & (A > 0)).sum(), *(np.mean(s[m & (A > 0)] == i) for i in range(3))))
    m2 = seen & (group == 2)
    lines.append("   seen 0-40 deg facets at %s: day %.3f, night %.3f, cast shadow %.3f" % (t, *(np.mean(s[m2] == i) for i in range(3))))
(out / "summary.txt").write_text("\n".join(lines) + "\n")
(out / "summary.json").write_text(json.dumps(summary))
print("\n".join(lines))
