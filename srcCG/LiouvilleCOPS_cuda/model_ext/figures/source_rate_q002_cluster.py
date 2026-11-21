#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/source_rate_q002_cluster.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/source_rate_q002_cluster.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 97dec33e5d8ec8935c901d2b9f0fde28dc9ea023d3ab17e39a9441728503d79e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: cluster page of the new-law source-rate maps: the production-rate panels of every completed run, one row per
#   run, per map time (2026-11-18)
# RUN: CPU (srun) once two or more runs are done; refuses an existing cluster/
# ==== end of header ====
"""Cluster page of the new-law source-rate maps (2026-11-18; user: "prepare the surface source rate maps for completed
runs for inspection of the source rate"; lead: "a cluster page comparing the source-rate maps of all completed runs").
CPU, light: composes the existing per-run panels (pub_maps_v2_figures.py, <run>/production_rate/panels/*_a_*_global.png
and *_b_*_impact_area.png, each with its own colour classes) into one page per time, one row per completed run, with
the run's Q_H2O (total, seen) and dark share from its inspection numbers.

    source_rate_q002_cluster.py --root figures_cluster/12M/source_rate_q002 [--runs D2_q002_warm,...]

Output <root>/cluster/: 12M_source_rate_q002_cluster_<hhmmss>UTC.{png,eps} (300 dpi) + _caption.txt, cluster.json.
Refuses an existing cluster/ (move it aside as cluster.superseded_<UTC> to rebuild with more runs).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402

TIMES = (("10:33:10", "103310"), ("10:39:26", "103926"))
RUNS = Path(__file__).resolve().parent / "data" / "source_rate_q002_runs_20261118.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--runs", default=None, help="comma list (default: every run of the runs spec with panels)")
    a = ap.parse_args()
    root = Path(a.root)
    out = root / "cluster"
    if out.exists():
        raise SystemExit("%s exists (move it aside as cluster.superseded_<UTC>)" % out)
    spec = json.load(open(RUNS))
    keys = a.runs.split(",") if a.runs else list(spec["runs"])
    plt = FIG.style()
    out.mkdir(parents=True)
    from PIL import Image
    meta = {"runs": {}, "figures": []}
    for hhmmss, tag in TIMES:
        rows = []
        for k in keys:
            pdir = root / k / "production_rate" / "panels"
            pa = pdir / ("12M_templates_h2o_production_rate_%sUTC_a_h2o_production_rate_global.png" % tag)
            pb = pdir / ("12M_templates_h2o_production_rate_%sUTC_b_h2o_production_rate_impact_area.png" % tag)
            if not (pa.exists() and pb.exists()):
                continue
            nj = root / k / "inspection" / "checks" / ("inspect_%s.json" % tag)
            num = json.load(open(nj))["new"] if nj.exists() else None
            rows.append((k, pa, pb, num))
            meta["runs"].setdefault(k, {})[hhmmss] = {"panels": [str(pa), str(pb)], "numbers": str(nj) if num else None}
        if not rows:
            continue
        ims = [(np.asarray(Image.open(pa).convert("RGB")), np.asarray(Image.open(pb).convert("RGB"))) for _, pa, pb, _ in rows]
        wa = max(i[0].shape[1] for i in ims)
        wb = max(i[1].shape[1] for i in ims)
        ra = max(i[0].shape[0] / i[0].shape[1] for i in ims)
        W = FIG.DOUBLE_COLUMN_IN
        wfa = W * wa / (wa + wb)
        rh = max(wfa * ra, (W - wfa) * max(i[1].shape[0] / i[1].shape[1] for i in ims))
        H = len(rows) * (rh + 0.45) + 0.3
        fig = plt.figure(figsize=(W, H))
        gs = fig.add_gridspec(len(rows), 2, width_ratios=[wa, wb], wspace=0.02, hspace=0.45 / rh, left=0.0, right=1.0,
                              top=1.0 - 0.3 / H, bottom=0.0)
        for r, ((k, pa, pb, num), (ia, ib)) in enumerate(zip(rows, ims)):
            for c, im in enumerate((ia, ib)):
                ax = fig.add_subplot(gs[r, c])
                ax.imshow(im, interpolation="lanczos")
                ax.axis("off")
            lab = "(%s) %s" % ("abcdefghijklmnop"[r], k)
            if num:
                lab += ": Q_H₂O %.2e s⁻¹ (seen %.2e, %.0f %%), night + cast shadow %.0f %% of Q" % (
                    num["Q_H2O_s-1"], num["Q_H2O_seen_s-1"], 100 * num["seen_share_of_Q"], 100 * num["dark_share_of_Q"])
            ax0 = fig.axes[2 * r]
            ax0.apply_aspect()
            fig.text(0.0, ax0.get_position().y1 + 0.16 / fig.get_size_inches()[1], lab, fontsize=7, ha="left",
                     va="bottom")
        prefix = out / ("12M_source_rate_q002_cluster_%sUTC" % tag)
        files = FIG.save(fig, prefix)
        plt.close(fig)
        cap = ("Surface H2O production rate F_H2O (instantaneous emitted H2O molecular flux per facet, m^-2 s^-1) of the "
               "completed runs fitted with the NEW H2O law q_min = 0.02 q_max (srcCG's value, the night-side fix of "
               "2026-11-18; every fit and figure before 2026-11-18 used q_min = 0.2 q_max), at %s UTC on 2016-09-30, one "
               "row per run (%s): left the global lon/lat panel, right the impact area, exactly the panels (a) and (b) "
               "of each run's production_rate/12M_templates_h2o_production_rate_%sUTC figure (pub_maps_v2_figures.py: "
               "colour only on the facets seen during the fit, seen F_H2O = 0 light green, unseen grey; per-panel log10 "
               "classes with the 2 %% area cut, i.e. EACH PANEL HAS ITS OWN COLOUR CLASSES: compare colours only within "
               "a panel; the run's inspection figure compares it with its old-law counterpart on common classes). Row "
               "labels: Q_H2O total and of the seen facets and the night + cast-shadow share of Q_H2O at this time "
               "(<run>/inspection/checks/inspect_%s.txt). The model of each run, its fit state and caveats: the run's "
               "README.txt and captions." % (hhmmss, ", ".join(k for k, _, _, _ in rows), tag, tag))
        with open(prefix.with_name(prefix.name + "_caption.txt"), "x") as fh:
            fh.write(cap + "\n")
        meta["figures"].append({"time": hhmmss, "files": files, "runs": [k for k, _, _, _ in rows]})
    with open(out / "cluster.json", "x") as fh:
        json.dump(meta, fh, indent=1)
    print("wrote %s" % [f["files"] for f in meta["figures"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
