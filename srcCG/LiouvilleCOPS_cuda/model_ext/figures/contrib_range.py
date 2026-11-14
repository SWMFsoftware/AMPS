#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/contrib_range.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/contrib_range.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 59f7c2fb9e02075747aa9cf18bdc0cc8f1fa3e4a0b5267082801a0b93ca72182
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Fixed log colour range of the contribution movies (bottom: 90 % of the surface signal at frame 0; top:
#   99.99th percentile over a frame range).
# RUN: srun on normal128cores (CPU)
# ==== end of header ====
"""Fixed log colour range of the contribution movies from frame_contrib.py output (2026-11-14): bottom = the per-area
level above which the facets carry --bottom-share of the gauge's surface signal at frame 0 (full_f00000.npz, FP64);
top = the --top-percentile percentile of the per-area values > 0 pooled over frames --top-from..--top-to (fNNNNN.npz).
Why (2026-11-14): with the literal rule (top from the last 10 frames, 1949-1958) the RG top is 7.9e-20 m^-4 s^-1: the
fitted model's RG surface flux falls from 8.6e16 m^-2 s^-1 at the last COPS observation (frame 1799) to 1.4e-18 at the
first contact (RG then sees only facets at the activity floor); the same rule on the last 10 frames with COPS data
(1790-1799) gives a usable top for both gauges.
Usage: contrib_range.py <frame_contrib dir> <out json> [--top-from 1790 --top-to 1799 --top-percentile 99.99 --bottom-share 0.9]"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdir")
    ap.add_argument("out")
    ap.add_argument("--top-from", type=int, default=1790)
    ap.add_argument("--top-to", type=int, default=1799)
    ap.add_argument("--top-percentile", type=float, default=99.99)
    ap.add_argument("--bottom-share", type=float, default=0.90)
    a = ap.parse_args()
    c = Path(a.cdir)
    if Path(a.out).exists():
        raise SystemExit("%s exists" % a.out)
    z0 = np.load(c / "full_f00000.npz")
    rng = {}
    for g, key, kf in (("NG", "vng", "ng"), ("RG", "vrg", "rg")):
        v, aa = z0[key], z0["area"]
        o = np.argsort(-v)
        cs = np.cumsum((v * aa)[o])
        j = int(np.searchsorted(cs, a.bottom_share * cs[-1]))
        bottom = float(v[o][min(j, len(o) - 1)])
        pool = np.concatenate([np.load(c / ("f%05d.npz" % k))[kf].astype(np.float64) for k in range(a.top_from, a.top_to + 1)])
        top = float(np.percentile(pool[pool > 0], a.top_percentile))
        lit = json.load(open(c / "color_range.json"))[g]["top"]
        rng[g] = {"bottom": bottom, "top": top, "decades": float(np.log10(top / bottom)),
                  "units": "m^-3 per m^2 (m^-5)" if g == "NG" else "m^-2 s^-1 per m^2 (m^-4 s^-1)",
                  "label": "NG: density per m² of surface [m⁻³ m⁻²]" if g == "NG" else "RG: flux into RG per m² [m⁻² s⁻¹ m⁻²]",
                  "bottom_rule": "per-area level above which the facets carry %g %% of the gauge's surface signal at frame 0"
                                 % (100 * a.bottom_share),
                  "top_rule": "%g th percentile of the per-area values > 0 of frames %d-%d (the last 10 frames with COPS data, "
                              "10:39:25.7-10:39:26.0 UTC)" % (a.top_percentile, a.top_from, a.top_to),
                  "top_literal_rule_frames_1949_1958": lit}
    json.dump(rng, open(a.out, "x"), indent=1)
    print(json.dumps(rng, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
