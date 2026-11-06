#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/refined_bundle.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/refined_bundle.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2fd538866a27c0e143b59df4f58f7391deb4acdc38cb594a449e57ea0bdb23df
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: template bundle for a refined mesh: every child element takes its parent's template values; keeps the source
#   bundles' gravity provenance and template definitions.
# USED FOR: jobs/refine_12m.sbatch, jobs/refine_20m.sbatch.
# RUN: .venv/bin/python <this> --base <refined run>/base --parent <mesh>/parent.npy --bundle <source bundle npz>
#        --source-restart <source fit_restart.npz> --out <new npz>   (as in jobs/refine_12m.sbatch, refine_20m.sbatch)
# ==== end of header ====
"""Template bundle for a refined mesh: every element takes its parent's template values (user 2026-11-04: "adopt
the slope from the original 12M mesh as the refinement should not change the slope significantly"; the visibility
templates are carried the same way).

    .venv/bin/python <this> --base <refined base> --parent <refine dir>/parent.npy --bundle <bundle> [--bundle ...]
        --out <new bundle .npz>
Columns of all given bundles (first occurrence of a name wins), values[child] = values[parent]; centroids and areas
from the refined base's elements.bin (the mesh identity the tools check).
The source bundles' gravity provenance (cache definition, orientation policy, integrity check) is carried over and
must be the same in every source bundle that has one; each definition keeps its source text, so the cliff/slope
provenance check of activity_templates.select_templates applies (first version dropped both: job 122002, kept as
refined_bundle.py.before_gravity_<UTC>).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path[:0] = [str(APP), str(APP / "tools")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--parent", required=True, type=Path)
    ap.add_argument("--bundle", required=True, action="append", type=Path)
    ap.add_argument("--source-restart", required=True, type=Path, help="a restart on the parent mesh (centroids)")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    from cops_cuda import provenance, reference, restart as rs
    templates = reference.load("activity_templates")
    fit_cops = reference.load("fit_cops")
    el = fit_cops.read_elements(str(a.base / "elements.bin"))
    parent = np.load(a.parent)
    if len(parent) != len(el["area"]):
        sys.exit("refined_bundle: parent map has %d entries for %d elements" % (len(parent), len(el["area"])))
    src = rs.read(a.source_restart)
    names, cols, defs, gravity = [], [], {}, None
    for b in a.bundle:
        bb = templates.load_bundle(b, src.centroid, src.area)
        g = bb["metadata"].get("gravity")
        if g is not None:
            if gravity is not None and g != gravity:
                sys.exit("refined_bundle: the source bundles carry different gravity provenance")
            gravity = g
        for k, n in enumerate(bb["names"]):
            if n in names:
                continue
            names.append(n)
            cols.append(bb["values"][parent, k])
            sdef = bb["metadata"].get("templates", {}).get(n) or {}
            defs[n] = {"definition": "values of the parent element (refined mesh); source bundle %s; source definition: %s"
                                     % (b, sdef.get("definition", "")), "source": sdef}
    extra = {} if gravity is None else {"gravity": dict(gravity, carried_from_parent_mesh=
        "evaluated on the parent mesh; each child takes its parent's value (geometry unchanged by the refinement)")}
    templates.write_bundle(a.out, names, np.column_stack(cols), el["centroid"], el["area"], {
        **extra, "templates": defs, "parent_map": provenance.file_identity(a.parent),
        "source_bundles": [provenance.file_identity(b) for b in a.bundle],
        "interpretation": "templates of the parent mesh carried to its refinement; not recomputed"})
    print("wrote %s: %s" % (a.out, names))


if __name__ == "__main__":
    main()
