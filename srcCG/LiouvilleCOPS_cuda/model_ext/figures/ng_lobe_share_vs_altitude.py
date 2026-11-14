#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/ng_lobe_share_vs_altitude.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/ng_lobe_share_vs_altitude.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4562e524c4ab8ee50cdc1f90e365b77ac642a00ca585b36724d054ae63e32761
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Big-lobe share and within-100-m share of the modelled NG surface signal per late observation (contribution
#   data source_areas_NG.npz), table and figure.
# RUN: srun on normal128cores (CPU)
# ==== end of header ====
"""Share of the modelled NG surface signal from each lobe per late observation (2026-11-14, panorama analysis): from the
contribution data source_areas_NG.npz (per late observation k: facets and their NG contributions, sorted; altitude from prep.npz rows by iPoint) and the lobe
labels of ng_source_regions.py (_lobe.npy). Writes <out>.json and <out>.tsv (utc, altitude, total, big-lobe share,
share within 100 m of the touchdown point) and a figure <out>.png/.eps (300 dpi)."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    for k in ("--src", "--lobe", "--frame", "--geometry", "--prep", "--out"):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    out = Path(a.out)
    if Path(str(out) + ".json").exists():
        raise SystemExit("exists")
    lobe = np.load(a.lobe)
    ptd = np.load(a.frame)["ptd"]
    fg = np.memmap(Path(a.geometry) / "mesh_facegeom.f8", dtype="<f8", mode="r").reshape(-1, 8)
    near = np.linalg.norm(np.asarray(fg[:, 0:3]) - ptd, axis=1) < 100.0
    del fg
    Z = np.load(a.src, allow_pickle=False)
    print(Z.files)
    PR = np.load(a.prep)
    altmap = {int(ip): float(al) for ip, al in zip(PR["rows_ipoint"], PR["rows_alt"])}
    off = Z["offsets"]
    fs = Z["faces_sorted"]
    cs = Z["contrib_sorted_Pa"]
    rows = []
    for k in range(len(off) - 1):
        f = fs[off[k]:off[k + 1]]
        c = cs[off[k]:off[k + 1]]
        t = float(c.sum())
        rows.append({"k": k, "utc": str(Z["utc"][k]) if "utc" in Z.files else "", "altitude_m": altmap.get(int(Z["iPoint"][k]), float("nan")), "total_Pa": t,
                     "big_lobe_share": float(c[lobe[f] == 1].sum() / t) if t > 0 else float("nan"),
                     "within_100m_share": float(c[near[f]].sum() / t) if t > 0 else float("nan")})
    json.dump({"definition": __doc__, "inputs": vars(a), "rows": rows}, open(str(out) + ".json", "x"), indent=1)
    with open(str(out) + ".tsv", "x") as fo:
        fo.write("k\tutc\taltitude_m\ttotal_Pa\tbig_lobe_share\twithin_100m_share\n")
        for r in rows:
            fo.write("%d\t%s\t%.3f\t%.6g\t%.4f\t%.4f\n" % (r["k"], r["utc"], r["altitude_m"], r["total_Pa"], r["big_lobe_share"],
                                                       r["within_100m_share"]))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    alt = np.array([r["altitude_m"] for r in rows])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(alt, [r["big_lobe_share"] for r in rows], "o-", ms=3, label="big lobe (body)")
    ax.plot(alt, [r["within_100m_share"] for r in rows], "s-", ms=3, label="within 100 m of the touchdown point")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("COPS altitude [m] (descent →)")
    ax.set_ylabel("share of the modelled NG surface signal")
    ax.set_title("95 late observations (c9 state)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for e in ("png", "eps"):
        fig.savefig(str(out) + "." + e, dpi=300)
    for r in rows[::5] + [rows[-1]]:
        print("%s %8.1f m  big lobe %.3f  near %.3f" % (r["utc"], r["altitude_m"], r["big_lobe_share"], r["within_100m_share"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
