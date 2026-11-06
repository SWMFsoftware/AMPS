#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/ft_template.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/ft_template.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ae4b80c682a638afe5efc45b881e04b7d4c6eeaa15014892a08fffdf3fc1dc02
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: turns a fine tuning's activity correction (fine_tuning_state.npz delta_activity) into an activity template
#   (<out>/<name>.templates.npz) on the same or another mesh; the integral is conserved (nearest-facet mapping).
# USED FOR: 20M fit from the 12M fine-tuned state (template ft12m_late, grain_ft20.sbatch); the 12M fit/fine-tuning
#   cycles (template ft_c<k> from the previous cycle's fine tuning, same mesh, jobs/grain_ft_cycle*.sbatch).
# RUN: .venv/bin/python ft_template.py --finetune <fine-tuning dir> --source-restart <its fit_restart.npz>
#        --target-restart <restart of the target mesh> --out <new dir> --name <template name>
# OUT: <out>/<name>.templates.npz (+ .json: facets, area, integral, clip check). The template is then given to the fit
#   as --activity-template <npz> <name> and held at amplitude 1 with --fix c_<name>=1.
# ==== end of header ====
"""The activity correction of a fine tuning as a template on another mesh (user request 2026-11-04: "run option 2
on 20M": a 20M fit starting from the 12M fine-tuned state, the fine-tuning correction carried to 20M as a fixed
added term).

    .venv/bin/python <this> --finetune <fine-tuning dir> --source-restart <the tuned fit's fit_restart.npz>
        --target-restart <restart on the target mesh> --out <new dir> [--name ft12m_late]

* The correction is delta_activity of <finetune>/fine_tuning_state.npz (clipped H2O activity added by the fine
  tuning, per source facet; all >= 0 for the fine tunings so far). The integral is conserved (user 2026-11-04):
  each corrected source facet's delta * area goes to the target facet with the nearest centroid,
  delta_t = sum(delta_s * area_s) / area_t. (The first version carried the per-area value to the target facets
  nearest to corrected source facets; with a correction concentrated on one small facet that gave 15x the
  integral on 20M: test job 121984; that version is kept as ft_template.py.before_conserve_<UTC>.)
* Template column = delta / A_scale(target) (the fitter adds amplitude * A_scale * column to the linear activity
  before the clip); with the amplitude held at 1 it adds delta. Where the fitted linear activity is >= 0 this equals
  the fine tuning's addition to the clipped field; where it is < 0 (clipped) it adds max(lin + delta, 0) instead of
  delta (reported: the difference at the target restart's state).
* Products: <out>/<name>.templates.npz (activity_templates.write_bundle, new file), <out>/ft_template.json.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path[:0] = [str(APP), str(APP / "tools")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--finetune", required=True, type=Path)
    ap.add_argument("--source-restart", required=True, type=Path)
    ap.add_argument("--target-restart", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--name", default="ft12m_late")
    a = ap.parse_args()
    if a.out.exists():
        sys.exit("ft_template: %s exists" % a.out)
    from scipy.spatial import cKDTree
    from cops_cuda import provenance, reference, restart as rs
    templates = reference.load("activity_templates")
    templates.validate_name(a.name)
    st = np.load(a.finetune / "fine_tuning_state.npz")
    meta = json.loads(str(st["meta_json"]))
    delta = np.asarray(st["delta_activity"], float)
    src = rs.read(a.source_restart)
    if meta["source_restart"]["sha256"] != provenance.file_identity(a.source_restart)["sha256"]:
        sys.exit("ft_template: the fine tuning was not made from %s" % a.source_restart)
    if np.any(delta < 0):
        sys.exit("ft_template: negative corrections (multipliers < 1) cannot be a nonnegative template")
    tgt = rs.read(a.target_restart)
    on = np.flatnonzero(delta > 0)
    # conserve the integral (user 2026-11-04): each corrected source facet's added activity * area goes to its
    # nearest target facet centroid; delta_t = sum(delta_s * area_s) / area_t
    dist_s, idx_s = cKDTree(tgt.centroid).query(src.centroid[on])
    d_t = np.zeros(len(tgt.area))
    np.add.at(d_t, idx_s, delta[on] * src.area[on])
    d_t /= tgt.area
    scale = float(tgt.enhancement_scale)
    column = d_t / scale
    a.out.mkdir(parents=True)
    bundle = a.out / ("%s.templates.npz" % a.name)
    templates.write_bundle(bundle, [a.name], column[:, None], tgt.centroid, tgt.area, {
        "templates": {a.name: {
            "definition": "activity correction (delta_activity) of the fine tuning %s, integral carried to the nearest "
                          "target facet centroid (delta_t = sum delta_s area_s / area_t), divided by A_scale; amplitude 1 adds it" % a.finetune,
            "source_restart": provenance.file_identity(a.source_restart),
            "fine_tuning_state": provenance.file_identity(a.finetune / "fine_tuning_state.npz"),
            "target_restart": provenance.file_identity(a.target_restart), "A_scale": scale,
            "interpretation": "fixed correction carried from another mesh; not a fitted result on this mesh"}},
        "source_mesh_identity": src.mesh_identity()})
    lin = np.asarray(tgt.activity_linear, float)
    pre = np.maximum(lin + d_t, 0.0)
    post = np.maximum(lin, 0.0) + d_t
    m = d_t > 0
    sa, ta = src.area, tgt.area
    summary = {"name": a.name, "bundle": str(bundle), "A_scale_target": scale,
               "source": {"facets": int(len(on)), "area_m2": float(sa[on].sum()),
                          "area_integral": float((sa * delta).sum()), "max": float(delta.max())},
               "target": {"facets": int(m.sum()), "area_m2": float(ta[m].sum()),
                          "area_integral": float((ta * d_t).sum()), "max": float(d_t.max()),
                          "distance_m": {"median": float(np.median(dist_s)), "max": float(dist_s.max())}},
               "clip_semantics_at_target_restart": {
                   "facets_with_negative_linear_activity": int(np.count_nonzero(m & (lin < 0))),
                   "relative_area_integral_difference_pre_vs_post_clip":
                       float((ta * np.abs(pre - post))[m].sum() / max((ta * post)[m].sum(), 1e-300))}}
    (a.out / "ft_template.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1), flush=True)


if __name__ == "__main__":
    main()
