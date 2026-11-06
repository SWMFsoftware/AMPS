#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/refine_mesh.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/refine_mesh.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c15019893295f42d3e41dfd5d1208d3cd79b3c14622fccf1befe315b4f736a55
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: conforming refinement (edges halved, 1 -> 4, closure 1 -> 2) of chosen elements of a NASTRAN BDF shape model,
#   geometry unchanged; selection from a fine tuning (--first-interval) or a mask (--seed-mask) + a margin.
# USED FOR: the 12M late-seen refinement (jobs/refine_12m.sbatch) and the 20M seen-after-gap refinement
#   (jobs/refine_20m.sbatch). Runs on a CPU (minutes for 12M, 47 min for 20M).
# RUN: python3 <this> --bdf <in.bdf> --centroids <restart.npz> (--first-interval <fine_tuning_prepared.npz> |
#        --seed-mask <seen.npz>) --margin <m> --out <new dir>
# OUT: <out>/<name>_refined.bdf, parent.npy, refine.json.
# ==== end of header ====
"""Refine chosen surface elements of a NASTRAN (BDF) shape model, conformingly, without changing the geometry
(user request 2026-11-04: test whether a gradient of the source rate across a surface element improves the fit;
answers: edges halved (each element -> 4), the elements the 12M fine tuning can change, i.e. first seen in its last
5 min, plus a margin).

    python3 <this> --bdf <in.bdf> --centroids <restart.npz with centroid, area> --first-interval <fine_tuning_prepared.npz>
        --margin 10 --out <new dir>

* Selection: elements with first_interval >= 0 in the fine tuning (first seen inside its window), plus every element
  whose centroid lies within --margin metres of one of their centroids.
* Or (--seed-mask <file.npz>, user 2026-11-04 for 20M: "the entire part of the surface seen after the gap by
  either NG or RG"): the elements where the boolean array `seen` of that file (model_ext/seen_after_gap.py) is
  true, plus those within --margin (0: none).
* Red refinement (midpoint subdivision, 1 -> 4) of the selection; closure: an element with two or three edges split
  becomes red too (repeated until stable); an element with one edge split is bisected (green, 1 -> 2). New vertices
  are edge midpoints, so every child lies in its parent's plane and the surface is unchanged. Orientation is kept.
* Order: the elements in input order, each refined one replaced by its children in place; new GRIDs appended.
* Checks: every edge is used by exactly two elements with opposite orientation (closed, consistently oriented),
  areas sum to the parent's.
* Products: <out>/<name>_refined.bdf, <out>/parent.npy (input element index of every output element, int64),
  <out>/refine.json (counts, areas, settings).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def read_bdf(path):
    ids, xyz, tri, header = [], [], [], []
    with open(path) as f:
        lines = iter(f)
        for line in lines:
            if line.startswith("GRID*"):
                nxt = next(lines)
                ids.append(int(line[8:24]))
                xyz.append((float(line[40:56]), float(line[56:72]), float(nxt[8:24])))
            elif line.startswith("CTRIA3"):
                tri.append((int(line[24:32]), int(line[32:40]), int(line[40:48])))
            elif line.startswith("$"):
                header.append(line.rstrip("\n"))
            elif line.startswith("ENDDATA") or not line.strip():
                continue
            elif line.startswith("GRID"):
                sys.exit("small-field GRID cards are not supported")
            else:
                sys.exit("unexpected card: %r" % line[:40])
    ids = np.asarray(ids, np.int64)
    order = np.argsort(ids)
    xyz = np.asarray(xyz, float)
    tri_ids = np.asarray(tri, np.int64)
    tri_idx = order[np.searchsorted(ids[order], tri_ids)]
    if not np.array_equal(ids[tri_idx], tri_ids):
        sys.exit("a CTRIA3 references an undefined GRID")
    return ids, xyz, tri_idx, header


def edge_keys(t):
    e = np.stack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]], axis=1)       # (nf, 3, 2) directed
    lo, hi = e.min(axis=2), e.max(axis=2)
    return lo * (1 << 32) + hi, e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdf", required=True, type=Path)
    ap.add_argument("--centroids", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--first-interval", type=Path)
    g.add_argument("--seed-mask", type=Path)
    ap.add_argument("--margin", type=float, default=10.0)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    if a.out.exists():
        sys.exit("refine_mesh: %s exists" % a.out)
    from scipy.spatial import cKDTree
    ids, xyz, tri, header = read_bdf(a.bdf)
    nf = len(tri)
    r = np.load(a.centroids)
    cen = r["centroid"]
    if len(cen) != nf:
        sys.exit("refine_mesh: %d centroids for %d elements" % (len(cen), nf))
    if a.first_interval is not None:
        seed = np.flatnonzero(np.load(a.first_interval)["first_interval"] >= 0)
    else:
        mask = np.load(a.seed_mask)["seen"]
        if len(mask) != nf:
            sys.exit("refine_mesh: seed mask has %d entries for %d elements" % (len(mask), nf))
        seed = np.flatnonzero(mask)
    if a.margin > 0:
        d, _ = cKDTree(cen[seed]).query(cen, distance_upper_bound=a.margin)
        red = np.isfinite(d)
    else:
        red = np.zeros(nf, bool); red[seed] = True
    nsel = int(red.sum())
    keys, directed = edge_keys(tri)
    # closure: an element with >= 2 split edges becomes red
    while True:
        split = np.zeros(0, np.int64) if not red.any() else np.unique(keys[red].ravel())
        nsplit = np.isin(keys, split).sum(axis=1)
        promote = (~red) & (nsplit >= 2)
        if not promote.any():
            break
        red |= promote
    green = (~red) & (nsplit == 1)
    # midpoints of the split edges: new GRIDs
    split = np.unique(keys[red].ravel())
    lo, hi = split >> 32, split & ((1 << 32) - 1)
    mid_xyz = 0.5 * (xyz[lo] + xyz[hi])
    mid_index = len(xyz) + np.arange(len(split))
    def mid(k):
        return mid_index[np.searchsorted(split, k)]
    out_tri, parent = [], []
    for i in range(nf):
        t = tri[i]
        if red[i]:
            m01, m12, m20 = mid(keys[i])
            out_tri += [(t[0], m01, m20), (m01, t[1], m12), (m20, m12, t[2]), (m01, m12, m20)]
            parent += [i] * 4
        elif green[i]:
            k = int(np.flatnonzero(np.isin(keys[i], split))[0])     # the split edge: (t[k], t[k+1])
            m = mid(keys[i][k])
            p, q, s = t[k], t[(k + 1) % 3], t[(k + 2) % 3]
            out_tri += [(p, m, s), (m, q, s)]
            parent += [i, i]
        else:
            out_tri.append(tuple(t))
            parent.append(i)
    out_tri = np.asarray(out_tri, np.int64)
    parent = np.asarray(parent, np.int64)
    allxyz = np.vstack([xyz, mid_xyz])
    # checks: closed and consistently oriented; areas preserved per parent
    k2, e2 = edge_keys(out_tri)
    k2 = k2.ravel()
    fwd = (e2[:, :, 0] < e2[:, :, 1]).ravel()
    order = np.argsort(k2, kind="stable")
    ks, fs = k2[order], fwd[order]
    uniq, start, count = np.unique(ks, return_index=True, return_counts=True)
    if not np.all(count == 2):
        sys.exit("refine_mesh: %d edges not used exactly twice" % int(np.sum(count != 2)))
    if not np.all(fs[start] != fs[start + 1]):
        sys.exit("refine_mesh: inconsistent orientation")
    def areas(xx, tt):
        return 0.5 * np.linalg.norm(np.cross(xx[tt[:, 1]] - xx[tt[:, 0]], xx[tt[:, 2]] - xx[tt[:, 0]]), axis=1)
    a_in, a_out = areas(xyz, tri), areas(allxyz, out_tri)
    per_parent = np.bincount(parent, weights=a_out, minlength=nf)
    rel = float(np.max(np.abs(per_parent - a_in) / np.maximum(a_in, 1e-300)))
    if rel > 1e-9:
        sys.exit("refine_mesh: child areas do not add up to the parent (max rel %.3g)" % rel)
    a.out.mkdir(parents=True)
    name = a.bdf.stem + "_refined"
    bdf = a.out / (name + ".bdf")
    new_ids = np.concatenate([ids, ids.max() + 1 + np.arange(len(split))])
    with open(bdf, "x") as f:
        for h in header:
            f.write(h + "\n")
        f.write("$ Refined by model_ext/refine_mesh.py (user 2026-11-04): %d elements red (1->4), %d green (1->2);\n"
                "$ new GRIDs are edge midpoints (geometry unchanged); element order = input order, refined ones in place\n"
                % (int(red.sum()), int(green.sum())))
        with open(a.bdf) as src:
            for line in src:
                if line.startswith("GRID*") or (line.startswith("*") and not line.startswith("$")):
                    f.write(line)
        for gid, (x, y, z) in zip(new_ids[len(ids):], mid_xyz):
            f.write("%-8s%16d%16d%16.9E%16.9E*%07d\n*%07d%16.9E\n" % ("GRID*", gid, 0, x, y, gid % 10000000,
                                                                    gid % 10000000, z))
        for e, (g1, g2, g3) in enumerate(new_ids[out_tri], start=1):
            f.write("%-8s%8d%8d%8d%8d%8d\n" % ("CTRIA3", e, 1, g1, g2, g3))
        f.write("ENDDATA\n")
    np.save(a.out / "parent.npy", parent)
    info = {"input_bdf": str(a.bdf), "output_bdf": str(bdf), "margin_m": a.margin,
            "seed": "first_interval >= 0 of %s" % a.first_interval if a.first_interval else "seen of %s" % a.seed_mask,
            "seed_elements": int(len(seed)), "selected_with_margin": nsel,
            "red_after_closure": int(red.sum()), "green": int(green.sum()),
            "elements_in": int(nf), "elements_out": int(len(out_tri)), "grids_in": int(len(ids)),
            "grids_out": int(len(new_ids)), "refined_area_m2": float(a_in[red | green].sum()),
            "max_rel_area_error_per_parent": rel, "closed_and_oriented": True}
    (a.out / "refine.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1), flush=True)


if __name__ == "__main__":
    main()
