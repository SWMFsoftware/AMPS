#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_mesh.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_mesh.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e433823e6e2df5272ca95119f3ddb02f05d4364d2b649b36b111ca8b9c44b3ca
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module + CLI (option C, local non-negative basis; lead 2026-11-14). Coarse triangulation of the 12M nucleus
#   surface (geodesic Voronoi / Lloyd on the facet graph, dual triangles), hat functions phi_j at every facet (barycentric
#   coordinates on the coarse triangle, partition of unity checked), P1 stiffness Lg (smoothness), facet mass matrix M.
# RUN: jobs/localnn_build.sbatch STEP mesh (frozen model_ext_frozen/localnn_20261114T081912Z).
# ==== end of header ====
"""Option C of the non-negative activity models (user 2026-11-14, relayed by the lead): a LOCAL non-negative basis.
This file builds the coarse triangulation of the 12M nucleus surface and the piecewise-linear hat functions phi_j
evaluated at every 12M facet (fit: localnn_problem.py / localnn_fit.py; design matrices: localnn_design.py).

    A_e = m_b(e) * sum_j c_j phi_j(x_e),  c_j >= 0      (phi_j >= 0, sum_j phi_j = 1 at every facet)

Local, untracked model extension: nothing in the git tree is edited and no existing model_ext file is changed.

CONSTRUCTION (all on the 12M mesh itself; deterministic; numpy/scipy only)
1. Facet graph: the 12M facets (base elements.bin order = mesh_faces.i4 of a validated nonneg_maps geometry), an edge
   between two facets that share a mesh edge, length = |centroid_1 - centroid_2| (the 12M mesh is a closed, consistently
   oriented manifold: every mesh edge has exactly two facets; checked).
2. Number of coarse nodes n = round(2 S / (sqrt(3) h^2)), S = the total facet area, h = --spacing (150 m): the node count
   of an equilateral triangulation with edge length h.
3. Seeds: an area-weighted random sample of --sample facets (numpy default_rng(--seed)), Euclidean farthest-point
   sampling of n seeds from it (first seed = the sample's first facet).
4. Lloyd iterations (--lloyd): geodesic Voronoi partition of the facets (multi-source Dijkstra on the facet graph,
   scipy.sparse.csgraph.dijkstra min_only; every facet gets the nearest seed); each seed moves to the facet of its
   region closest (Euclidean) to the region's area-weighted mean centroid. The final partition is the Voronoi
   partition of the final seeds. Coarse node j = the centroid of its seed facet (a point ON the 12M surface).
5. Coarse triangles = the dual of the partition: every mesh vertex where exactly three regions meet gives the triangle
   of these three nodes (duplicates removed); a vertex where k >= 4 regions meet gives k - 2 triangles (the cyclic
   sequence of regions around the vertex, ear clipping by the shortest diagonal). Triangles are oriented along the
   area-weighted normal of the fine facets around their dual vertex.
6. Hat functions at a facet e of region i (barycentric interpolation on the coarse mesh): for every coarse triangle
   incident to node i, the barycentric coordinates of e's centroid projected onto the triangle's plane; the triangle
   with the largest smallest coordinate is taken. If that smallest coordinate is < 0 (the centroid lies outside every
   triangle of node i), the triangles incident to the nodes of node i's triangles are tried as well. Coordinates that
   are still negative (outside every candidate) are set to 0 and the three renormalized; their number, area and size
   are reported. phi_j(e) = the coordinate of node j (three per facet). Partition of unity sum_j phi_j(e) = 1 and
   phi >= 0 hold at every facet by construction and are checked (|sum - 1| <= 4 eps).
7. Smoothness operator: the P1 finite-element stiffness matrix of the coarse mesh, Lg = sum_T area_T grad(phi)^T
   grad(phi) (the cotangent-weighted graph Laplacian; c^T Lg c = integral over the coarse surface of
   |grad sum_j c_j phi_j|^2), assembled from the per-triangle gradient operator Dg (3 rows per triangle,
   sqrt(area_T) grad lambda_k), so Lg = Dg^T Dg exactly; Lg 1 = 0 checked.
8. Mass matrix on the facets (for the start projection): M = Phi^T diag(area) Phi, Phi = the facet x node matrix of phi.

DIAGNOSTICS (grid_report.json): node count, spacing realized (coarse edge lengths), triangles, Euler characteristic,
coarse edges with 1 / 2 / >2 triangles, region pairs adjacent on the fine mesh without a coarse edge, nodes without a
triangle or without support, region areas, the facets outside every candidate triangle, the distance of facet centroids
from the planes of their triangles, and the continuity of phi across every fine-mesh edge (L1 jump
sum_j |phi_j(e1) - phi_j(e2)| = 2 - 2 sum_common min).

CLI: localnn_mesh.py --base BASE --geometry GEOM --spacing 150 --out <grid dir> [--lloyd 6] [--sample 400000]
     [--seed 20261114] [--seen-from <nonneg_facets.npz>]
Products in <grid dir> (written as <grid dir>.incomplete_<UTC>, renamed at the end; never overwritten):
  grid.npz         node_xyz, node_facet, tri, tri_area, facet_node (int32 nf x 3), facet_w (nf x 3), facet_label,
                   node_support_area, node_seen, Lg_{data,indices,indptr}, Dg_{data,indices,indptr,shape},
                   M_{data,indices,indptr}, meta_json
  grid_report.json, README.txt, coarse_mesh.vtu (ParaView: the coarse triangles with node ids), regions.npy (int32)
"""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra

HERE = Path(__file__).resolve().parent
SCHEMA = "localnn_grid_v1"
LOG = "localnn_mesh:"
EPS = np.finfo(float).eps


def say(text):
    print("%s [%s] %s" % (LOG, dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(bufsize), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_elements(path):
    """elements.bin (version 1): centroid, normal, area (the reader of nonneg_maps.read_elements, same layout)."""
    import struct
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        version, n, nb = struct.unpack("<iii", f.read(12))
        if version != 1 or nb != 25 or size != 12 + n * 34 * 8:
            raise SystemExit("unsupported element file %s" % path)
        d = np.fromfile(f, dtype="<f8", count=n * 34).reshape(n, 34)
    out = {"centroid": np.ascontiguousarray(d[:, 0:3]), "normal": np.ascontiguousarray(d[:, 3:6]),
           "area": np.ascontiguousarray(d[:, 6]), "AH2O": np.ascontiguousarray(d[:, 32]),
           "ACO2": np.ascontiguousarray(d[:, 33]), "nfaces": n}
    del d
    return out


# ------------------------------------------------------------------------------------------------ facet graph
def facet_graph(faces, nvert, centroid):
    """Symmetric CSR facet adjacency (shared mesh edge), weights = centroid distance; also the edge pairs."""
    nf = len(faces)
    a = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2]]).astype(np.int64)
    b = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0]]).astype(np.int64)
    key = np.minimum(a, b) * nvert + np.maximum(a, b)
    fid = np.tile(np.arange(nf, dtype=np.int64), 3)
    order = np.argsort(key, kind="stable")
    ks = key[order]
    same = ks[1:] == ks[:-1]
    if np.any(same[1:] & same[:-1]):
        raise SystemExit("an edge of the 12M mesh has more than two facets")
    i = np.flatnonzero(same)
    f1, f2 = fid[order[i]], fid[order[i + 1]]
    if 2 * len(i) != len(key):
        raise SystemExit("the 12M mesh has %d boundary edges (expected a closed surface)" % (len(key) - 2 * len(i)))
    w = np.sqrt(np.sum((centroid[f1] - centroid[f2]) ** 2, axis=1))
    if np.any(w <= 0):
        raise SystemExit("coincident facet centroids")
    G = sp.csr_matrix((np.concatenate([w, w]), (np.concatenate([f1, f2]), np.concatenate([f2, f1]))), shape=(nf, nf))
    return G, f1, f2


def voronoi(G, seeds):
    dist, _, src = dijkstra(G, directed=True, indices=seeds, return_predecessors=True, min_only=True)
    if np.any(src < 0):
        raise SystemExit("facets not reached from any seed (disconnected facet graph)")
    inv = np.full(G.shape[0], -1, dtype=np.int64)
    inv[seeds] = np.arange(len(seeds))
    label = inv[src]
    if np.any(label < 0):
        raise SystemExit("dijkstra returned a source that is not a seed")
    return label.astype(np.int32), dist


def farthest_points(x, n):
    d = np.full(len(x), np.inf)
    out = np.empty(n, dtype=np.int64)
    out[0] = 0
    for k in range(1, n):
        d = np.minimum(d, np.sum((x - x[out[k - 1]]) ** 2, axis=1))
        out[k] = int(np.argmax(d))
    return out


def group_argmin(label, value, n):
    order = np.lexsort((np.arange(len(value)), value, label))
    first = np.searchsorted(label[order], np.arange(n), side="left")
    if np.any(first >= len(order)) or np.any(label[order[first]] != np.arange(n)):
        raise SystemExit("an empty region")
    return order[first]


# ------------------------------------------------------------------------------------------------ dual triangles
def ear_clip(cycle, pos):
    """Triangles of a cyclic node sequence (ear clipping by the shortest diagonal)."""
    cyc = list(cycle)
    tris = []
    while len(cyc) > 3:
        best, bk = None, None
        for k in range(len(cyc)):
            a, c = cyc[k - 1], cyc[(k + 1) % len(cyc)]
            d = float(np.sum((pos[a] - pos[c]) ** 2))
            if best is None or d < best:
                best, bk = d, k
        k = bk
        tris.append((cyc[k - 1], cyc[k], cyc[(k + 1) % len(cyc)]))
        del cyc[k]
    if len(cyc) == 3:
        tris.append(tuple(cyc))
    return tris


def cyclic_labels(v, faces, inc_faces, label):
    """Region labels in cyclic order around mesh vertex v (consecutive duplicates merged)."""
    nxt = {}
    for f in inc_faces:
        t = list(faces[f])
        k = t.index(v)
        a, b = t[(k + 1) % 3], t[(k + 2) % 3]     # face (v, a, b) counter-clockwise
        nxt[a] = (b, f)
    start = next(iter(nxt))
    seq, cur = [], start
    for _ in range(len(nxt)):
        b, f = nxt[cur]
        seq.append(int(label[f]))
        cur = b
        if cur == start:
            break
    if cur != start or len(seq) != len(inc_faces):
        raise SystemExit("vertex %d: the incident facets do not form one fan" % v)
    out = [s for k, s in enumerate(seq) if s != seq[k - 1]] if len(set(seq)) > 1 else seq[:1]
    return out


def dual_triangles(faces, label, pos, normal, area):
    nv = int(faces.max()) + 1
    nf = len(faces)
    vv = faces.ravel().astype(np.int64)
    ff = np.repeat(np.arange(nf, dtype=np.int64), 3)
    ll = label[ff].astype(np.int64)
    nlab = int(label.max()) + 1
    key = vv * nlab + ll
    uk = np.unique(key)
    uv = uk // nlab
    cnt = np.bincount(uv, minlength=nv)
    tri_list = []
    tri_vertex = []
    # vertices with exactly three regions
    v3 = np.flatnonzero(cnt == 3)
    start = np.searchsorted(uv, v3)
    labs = np.stack([uk[start + k] % nlab for k in range(3)], axis=1)
    tri_list.append(np.sort(labs, axis=1))
    tri_vertex.append(v3)
    # vertices with k >= 4 regions
    vmany = np.flatnonzero(cnt >= 4)
    special = {"vertices_k_ge_4": int(len(vmany)), "k_histogram": {}, "repeated_label_vertices": 0}
    if len(vmany):
        order = np.argsort(vv, kind="stable")
        vs = vv[order]
        lo = np.searchsorted(vs, vmany, side="left")
        hi = np.searchsorted(vs, vmany, side="right")
        extra, extra_v = [], []
        for v, a, b in zip(vmany.tolist(), lo.tolist(), hi.tolist()):
            inc = ff[order[a:b]].tolist()
            cyc = cyclic_labels(v, faces, inc, label)
            k = len(cyc)
            special["k_histogram"][str(int(cnt[v]))] = special["k_histogram"].get(str(int(cnt[v])), 0) + 1
            if len(set(cyc)) != len(cyc):
                special["repeated_label_vertices"] += 1
                seen, u = set(), []
                for s in cyc:
                    if s not in seen:
                        seen.add(s)
                        u.append(s)
                cyc = u
            for t in ear_clip(cyc, pos):
                if len(set(t)) == 3:
                    extra.append(sorted(t))
                    extra_v.append(v)
        if extra:
            tri_list.append(np.array(extra, dtype=np.int64))
            tri_vertex.append(np.array(extra_v, dtype=np.int64))
    T = np.concatenate(tri_list).astype(np.int64)
    TV = np.concatenate(tri_vertex)
    tk = (T[:, 0] * nlab + T[:, 1]) * nlab + T[:, 2]
    uq, first, mult = np.unique(tk, return_index=True, return_counts=True)
    T, TV = T[first], TV[first]
    special["triple_points"] = int(len(tk))
    special["duplicate_triples_removed"] = int(len(tk) - len(uq))
    # orientation: along the area-weighted normal of the fine facets at the dual vertex
    wn = normal * area[:, None]
    vn = np.stack([np.bincount(vv, weights=wn[ff, k], minlength=nv) for k in range(3)], axis=1)
    tn = np.cross(pos[T[:, 1]] - pos[T[:, 0]], pos[T[:, 2]] - pos[T[:, 0]])
    flip = np.sum(tn * vn[TV], axis=1) < 0
    T[flip] = T[flip][:, [0, 2, 1]]
    special["flipped_for_orientation"] = int(flip.sum())
    return T.astype(np.int32), special


# ------------------------------------------------------------------------------------------------ barycentric
def bary(P, A, B, C):
    """Barycentric coordinates of the projection of P onto the plane of (A, B, C); arrays (..., 3)."""
    v0, v1, v2 = B - A, C - A, P - A
    d00 = np.sum(v0 * v0, axis=-1)
    d01 = np.sum(v0 * v1, axis=-1)
    d11 = np.sum(v1 * v1, axis=-1)
    d20 = np.sum(v2 * v0, axis=-1)
    d21 = np.sum(v2 * v1, axis=-1)
    den = d00 * d11 - d01 * d01
    beta = (d11 * d20 - d01 * d21) / den
    gamma = (d00 * d21 - d01 * d20) / den
    alpha = 1.0 - beta - gamma
    proj = A + beta[..., None] * v0 + gamma[..., None] * v1
    dist = np.sqrt(np.sum((P - proj) ** 2, axis=-1))
    return np.stack([alpha, beta, gamma], axis=-1), dist


def padded(lists, n):
    m = max(1, max((len(x) for x in lists), default=1))
    out = np.full((n, m), -1, dtype=np.int64)
    for i, x in enumerate(lists):
        out[i, :len(x)] = x
    return out


def best_triangles(xc, lab, cand, T, pos):
    """For facets with centroids xc and candidate triangle lists cand[lab] (padded with -1): best triangle (largest
    smallest coordinate), its coordinates and the plane distance."""
    n = len(xc)
    chunk = max(1 << 12, (1 << 23) // cand.shape[1])
    bt = np.full(n, -1, dtype=np.int64)
    bw = np.zeros((n, 3))
    bm = np.full(n, -np.inf)
    bd = np.zeros(n)
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        C = cand[lab[s:e]]                                   # (k, m)
        valid = C >= 0
        Cs = np.where(valid, C, 0)
        A, B, Cc = pos[T[Cs, 0]], pos[T[Cs, 1]], pos[T[Cs, 2]]
        w, d = bary(xc[s:e, None, :], A, B, Cc)
        mn = np.where(valid, w.min(axis=-1), -np.inf)
        k = np.argmax(mn, axis=1)
        r = np.arange(e - s)
        bt[s:e] = C[r, k]
        bw[s:e] = w[r, k]
        bm[s:e] = mn[r, k]
        bd[s:e] = d[r, k]
    return bt, bw, bm, bd


def stiffness(T, pos):
    """Per-triangle gradient operator Dg (3 rows per triangle: sqrt(area) * d/dx_{x,y,z} of the P1 interpolant) and
    Lg = Dg^T Dg (P1 stiffness / cotangent Laplacian)."""
    p0, p1, p2 = pos[T[:, 0]], pos[T[:, 1]], pos[T[:, 2]]
    nrm = np.cross(p1 - p0, p2 - p0)
    a2 = np.sqrt(np.sum(nrm ** 2, axis=1))
    if np.any(a2 <= 0):
        raise SystemExit("degenerate coarse triangle")
    nh = nrm / a2[:, None]
    area = 0.5 * a2
    grads = []
    for k in range(3):
        e = pos[T[:, (k + 2) % 3]] - pos[T[:, (k + 1) % 3]]           # edge opposite vertex k
        grads.append(np.cross(nh, e) / a2[:, None])                   # grad lambda_k (|.| = 1/height)
    nt = len(T)
    rows, cols, vals = [], [], []
    sq = np.sqrt(area)
    for comp in range(3):
        for k in range(3):
            rows.append(3 * np.arange(nt) + comp)
            cols.append(T[:, k])
            vals.append(sq * grads[k][:, comp])
    n = int(pos.shape[0])
    Dg = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * nt, n))
    Dg.sum_duplicates()
    Lg = (Dg.T @ Dg).tocsr()
    Lg.sum_duplicates()
    return Dg, Lg, area


def write_vtu(path, pos, T, node_ids):
    """Minimal ASCII VTK unstructured grid of the coarse triangles (ParaView)."""
    with open(path, "x") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
                '<UnstructuredGrid>\n<Piece NumberOfPoints="%d" NumberOfCells="%d">\n' % (len(pos), len(T)))
        f.write('<PointData Scalars="node"><DataArray type="Int32" Name="node" format="ascii">\n')
        f.write(" ".join(str(int(v)) for v in node_ids) + "\n</DataArray></PointData>\n")
        f.write('<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        f.write("\n".join("%.6f %.6f %.6f" % tuple(p) for p in pos) + "\n</DataArray></Points>\n")
        f.write('<Cells><DataArray type="Int32" Name="connectivity" format="ascii">\n')
        f.write("\n".join("%d %d %d" % tuple(t) for t in T) + "\n</DataArray>\n")
        f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        f.write(" ".join(str(3 * (k + 1)) for k in range(len(T))) + "\n</DataArray>\n")
        f.write('<DataArray type="UInt8" Name="types" format="ascii">\n' + " ".join(["5"] * len(T)) + "\n</DataArray>\n")
        f.write("</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")


def quant(a):
    a = np.asarray(a, float)
    if not a.size:
        return {}
    q = np.quantile(a, [0.0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0])
    return {"min": float(q[0]), "p01": float(q[1]), "p10": float(q[2]), "median": float(q[3]), "p90": float(q[4]),
            "p99": float(q[5]), "max": float(q[6]), "mean": float(np.mean(a))}


# ------------------------------------------------------------------------------------------------ main
def build(args):
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    base, geom = Path(args.base), Path(args.geometry)
    gval = json.load(open(geom / "geometry_validation.json"))
    if not gval.get("passed") or Path(gval["base"]).resolve() != base.resolve():
        raise SystemExit("geometry %s is not validated against %s" % (geom, base))
    el = read_elements(base / "elements.bin")
    nf = el["nfaces"]
    nodes = np.fromfile(geom / "mesh_nodes.f8", dtype="<f8").reshape(-1, 3)
    faces = np.fromfile(geom / "mesh_faces.i4", dtype="<i4").reshape(-1, 3)
    if len(faces) != nf or faces.min() != 0 or faces.max() != len(nodes) - 1:
        raise SystemExit("mesh_faces.i4 does not match elements.bin / mesh_nodes.f8")
    xc, area, normal = el["centroid"], el["area"], el["normal"]
    cdiff = float(np.max(np.abs((nodes[faces[:, 0]] + nodes[faces[:, 1]] + nodes[faces[:, 2]]) / 3.0 - xc)))
    if cdiff > 1e-6:
        raise SystemExit("facet centroids differ from the mesh by %.3e m" % cdiff)
    S = float(np.sum(area))
    n = int(round(2.0 * S / (math.sqrt(3.0) * args.spacing ** 2)))
    say("12M mesh: %d facets, %d vertices, area %.6e m^2; centroid check %.2e m; target %d nodes (h %.1f m)"
        % (nf, len(nodes), S, cdiff, n, args.spacing))
    G, f1, f2 = facet_graph(faces, len(nodes), xc)
    say("facet graph: %d edges (%.0f s)" % (len(f1), time.time() - t0))
    rng = np.random.default_rng(args.seed)
    sample = rng.choice(nf, size=min(args.sample, nf), replace=False, p=area / S)
    seeds = sample[farthest_points(xc[sample], n)]
    say("farthest-point seeds: %d (%.0f s)" % (len(seeds), time.time() - t0))
    lloyd = []
    for it in range(args.lloyd):
        label, dist = voronoi(G, seeds)
        ra = np.bincount(label, weights=area, minlength=n)
        mean = np.stack([np.bincount(label, weights=area * xc[:, k], minlength=n) for k in range(3)], axis=1) / ra[:, None]
        d2 = np.sum((xc - mean[label]) ** 2, axis=1)
        new = group_argmin(label, d2, n)
        move = np.sqrt(np.sum((xc[new] - xc[seeds]) ** 2, axis=1))
        lloyd.append({"iteration": it + 1, "region_area_m2": quant(ra), "seed_move_m": quant(move),
                      "max_geodesic_m": float(dist.max())})
        say("Lloyd %d: region area %.0f .. %.0f m^2 (cv %.3f), seed move median %.2f max %.2f m (%.0f s)" % (
            it + 1, ra.min(), ra.max(), ra.std() / ra.mean(), np.median(move), move.max(), time.time() - t0))
        seeds = np.sort(new)
    label, dist = voronoi(G, seeds)
    ra = np.bincount(label, weights=area, minlength=n)
    pos = xc[seeds].copy()
    say("final partition: region area %.0f .. %.0f m^2, max geodesic distance to the seed %.1f m" % (
        ra.min(), ra.max(), dist.max()))
    T, special = dual_triangles(faces, label, pos, normal, area)
    nt = len(T)
    say("dual: %d triangles (%s) (%.0f s)" % (nt, json.dumps(special), time.time() - t0))
    # coarse topology
    ekey = np.concatenate([np.minimum(T[:, a], T[:, b]).astype(np.int64) * n + np.maximum(T[:, a], T[:, b])
                           for a, b in ((0, 1), (1, 2), (2, 0))])
    ue, ecnt = np.unique(ekey, return_counts=True)
    lab1, lab2 = label[f1].astype(np.int64), label[f2].astype(np.int64)
    diff = lab1 != lab2
    adj = np.unique(np.minimum(lab1[diff], lab2[diff]) * n + np.maximum(lab1[diff], lab2[diff]))
    missing = np.setdiff1d(adj, ue)
    extra_e = np.setdiff1d(ue, adj)
    node_has_tri = np.zeros(n, dtype=bool)
    node_has_tri[T.ravel()] = True
    elen = np.sqrt(np.sum((pos[ue // n] - pos[ue % n]) ** 2, axis=1))
    topo = {"nodes": n, "triangles": nt, "edges": int(len(ue)), "euler_characteristic": int(n - len(ue) + nt),
            "edges_with_1_triangle": int(np.sum(ecnt == 1)), "edges_with_2_triangles": int(np.sum(ecnt == 2)),
            "edges_with_more_triangles": int(np.sum(ecnt > 2)),
            "adjacent_region_pairs": int(len(adj)), "adjacent_pairs_without_coarse_edge": int(len(missing)),
            "coarse_edges_between_non_adjacent_regions": int(len(extra_e)),
            "nodes_without_triangle": int(np.sum(~node_has_tri)), "edge_length_m": quant(elen)}
    say("coarse mesh: " + json.dumps({k: v for k, v in topo.items() if k != "edge_length_m"}))
    say("coarse edge length [m]: " + json.dumps(topo["edge_length_m"]))
    # node -> incident triangles, and the 2-ring candidates
    inc = [[] for _ in range(n)]
    for t, tri in enumerate(T.tolist()):
        for v in tri:
            inc[v].append(t)
    cand1 = padded(inc, n)
    ring = []
    for i in range(n):
        s = set(inc[i])
        for t in inc[i]:
            for v in T[t].tolist():
                s.update(inc[v])
        ring.append(sorted(s))
    cand2 = padded(ring, n)
    say("candidates: up to %d triangles per node, %d in the 2-ring" % (cand1.shape[1], cand2.shape[1]))
    bt, bw, bm, bd = best_triangles(xc, label, cand1, T, pos)
    out1 = bm < 0
    say("facets inside a triangle of their node: %d of %d (%.4f of the area)" % (
        nf - int(out1.sum()), nf, 1 - float(np.sum(area[out1])) / S))
    if np.any(out1):
        idx = np.flatnonzero(out1)
        bt2, bw2, bm2, bd2 = best_triangles(xc[idx], label[idx], cand2, T, pos)
        better = bm2 > bm[idx]
        j = idx[better]
        bt[j], bw[j], bm[j], bd[j] = bt2[better], bw2[better], bm2[better], bd2[better]
    outside = bm < 0
    clip = {"facets_outside_node_triangles": int(out1.sum()),
            "area_fraction_outside_node_triangles": float(np.sum(area[out1]) / S),
            "facets_outside_all_candidates": int(outside.sum()),
            "area_fraction_outside_all_candidates": float(np.sum(area[outside]) / S),
            "smallest_coordinate_outside": quant(bm[outside]) if outside.any() else {}}
    w = np.where(bw < 0, 0.0, bw)
    w = w / np.sum(w, axis=1, keepdims=True)
    facet_node = T[bt].astype(np.int32)
    ssum = np.sum(w, axis=1)
    pu = {"max_abs_sum_minus_1": float(np.max(np.abs(ssum - 1.0))), "min_weight": float(w.min()),
          "passed": bool(np.max(np.abs(ssum - 1.0)) <= 4 * EPS and w.min() >= 0.0)}
    say("partition of unity: %s; outside: %s" % (json.dumps(pu), json.dumps(clip)))
    if not pu["passed"]:
        raise SystemExit("partition of unity failed")
    # continuity across fine edges: L1 jump = 2 - 2 sum over common nodes of min(w1, w2)
    jump = np.empty(len(f1))
    for s in range(0, len(f1), 1 << 21):
        e = min(len(f1), s + (1 << 21))
        a, b = f1[s:e], f2[s:e]
        na, nb, wa, wb = facet_node[a], facet_node[b], w[a], w[b]
        common = np.zeros(e - s)
        for p in range(3):
            for q in range(3):
                common += np.where(na[:, p] == nb[:, q], np.minimum(wa[:, p], wb[:, q]), 0.0)
        jump[s:e] = 2.0 - 2.0 * common
    cont = {"L1_jump_across_fine_edges": quant(jump), "edges_jump_gt_0.25": int(np.sum(jump > 0.25)),
            "edges_jump_gt_0.5": int(np.sum(jump > 0.5)), "edges_jump_gt_1": int(np.sum(jump > 1.0)),
            "fine_edges": int(len(jump)), "typical_expected": "about 2 * facet size / h (~0.04 for 3 m facets, h 150 m)"}
    say("continuity: " + json.dumps(cont))
    # matrices
    rowsP = np.repeat(np.arange(nf, dtype=np.int64), 3)
    Phi = sp.csr_matrix((w.ravel(), (rowsP, facet_node.ravel().astype(np.int64))), shape=(nf, n))
    support = np.asarray(Phi.T @ area).ravel()
    M = (Phi.T @ sp.diags(area) @ Phi).tocsr()
    del Phi
    Dg, Lg, tri_area = stiffness(T, pos)
    l1 = float(np.max(np.abs(Lg @ np.ones(n))))
    diag = Lg.diagonal()
    say("stiffness: |Lg 1|_max %.3e (diag median %.3f), mass: %d nonzeros, support area %.0f .. %.0f m^2" % (
        l1, float(np.median(diag)), M.nnz, support.min(), support.max()))
    node_seen = None
    seen_src = None
    if args.seen_from:
        with np.load(args.seen_from, allow_pickle=False) as z:
            seen = np.asarray(z["seen"], bool)
        if len(seen) != nf:
            raise SystemExit("--seen-from has %d facets" % len(seen))
        sf = np.flatnonzero(seen)
        node_seen = np.zeros(n, dtype=bool)
        for k in range(3):
            pos_k = sf[w[sf, k] > 0]
            node_seen |= np.bincount(facet_node[pos_k, k], minlength=n) > 0
        seen_src = {"path": str(Path(args.seen_from).resolve()), "sha256": sha256_file(args.seen_from),
                    "seen_facets": int(seen.sum()), "seen_area_fraction": float(np.sum(area[seen]) / S),
                    "nodes_with_seen_support": int(node_seen.sum())}
        say("seen: %s" % json.dumps(seen_src))
    stage.mkdir(parents=True)
    report = {"schema": SCHEMA, "created_utc": utc_now(), "host": socket.gethostname(),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
              "inputs": {"base": str(base.resolve()), "geometry": str(geom.resolve()),
                         "elements_bin": str((base / "elements.bin").resolve()),
                         "mesh_nodes_sha256": sha256_file(geom / "mesh_nodes.f8"),
                         "mesh_faces_sha256": sha256_file(geom / "mesh_faces.i4"),
                         "geometry_validation_sha256": sha256_file(geom / "geometry_validation.json")},
              "settings": {"spacing_m": args.spacing, "lloyd": args.lloyd, "sample": args.sample, "seed": args.seed,
                           "node_count_rule": "round(2 S / (sqrt(3) h^2)), S = total facet area"},
              "surface_area_m2": S, "nodes": n, "triangles": nt, "lloyd_iterations": lloyd,
              "region_area_m2": quant(ra), "region_max_geodesic_m": float(dist.max()), "dual": special,
              "topology": topo, "mapping": clip, "plane_distance_m": quant(bd),
              "partition_of_unity": pu, "continuity": cont,
              "stiffness": {"max_abs_Lg_ones": l1, "diag": quant(diag), "nnz": int(Lg.nnz)},
              "mass": {"nnz": int(M.nnz)}, "node_support_area_m2": quant(support),
              "nodes_zero_support": int(np.sum(support <= 0)), "seen": seen_src,
              "code_sha256": {"localnn_mesh.py": sha256_file(__file__)}, "wall_seconds": time.time() - t0}
    meta = {k: report[k] for k in ("schema", "created_utc", "inputs", "settings", "surface_area_m2", "nodes",
                                   "triangles", "partition_of_unity", "code_sha256")}
    np.savez(stage / "grid.npz", node_xyz=pos, node_facet=seeds.astype(np.int64), tri=T, tri_area=tri_area,
             facet_node=facet_node, facet_w=w, facet_label=label, node_support_area=support,
             node_seen=(node_seen if node_seen is not None else np.zeros(0, bool)),
             Lg_data=Lg.data, Lg_indices=Lg.indices, Lg_indptr=Lg.indptr,
             Dg_data=Dg.data, Dg_indices=Dg.indices, Dg_indptr=Dg.indptr, Dg_shape=np.array(Dg.shape),
             M_data=M.data, M_indices=M.indices, M_indptr=M.indptr, meta_json=json.dumps(meta, sort_keys=True))
    write_vtu(stage / "coarse_mesh.vtu", pos, T, np.arange(n))
    with open(stage / "grid_report.json", "x") as f:
        json.dump(report, f, indent=1, sort_keys=True, default=float)
    write_readme(stage, report)
    report_sha = sha256_file(stage / "grid.npz")
    stage.rename(out)
    say("done: %s (grid.npz sha256 %s, %.0f s)" % (out, report_sha, time.time() - t0))
    return 0


def write_readme(stage, r):
    t, m = r["topology"], r["mapping"]
    L = ["Local non-negative basis (option C, user 2026-11-14): coarse triangulation of the 12M nucleus surface",
         "Written %s by localnn_mesh.py (%s) on %s (Slurm job %s)." % (r["created_utc"], r["schema"], r["host"],
                                                                     r["slurm_job_id"]),
         "",
         "Nodes %d (target spacing %.0f m: n = round(2 S / (sqrt(3) h^2)), S = %.6e m^2), triangles %d, edges %d," % (
             r["nodes"], r["settings"]["spacing_m"], r["surface_area_m2"], r["triangles"], t["edges"]),
         "Euler characteristic %d (the 12M surface has genus 4: -6); coarse edge length median %.1f m (p10 %.1f, p90 %.1f)."
         % (t["euler_characteristic"], t["edge_length_m"]["median"], t["edge_length_m"]["p10"], t["edge_length_m"]["p90"]),
         "Coarse edges with 1 / 2 / >2 triangles: %d / %d / %d; adjacent regions without a coarse edge: %d." % (
             t["edges_with_1_triangle"], t["edges_with_2_triangles"], t["edges_with_more_triangles"],
             t["adjacent_pairs_without_coarse_edge"]),
         "Facets outside every candidate triangle (coordinates clipped at 0, renormalized): %d (%.2e of the area)." % (
             m["facets_outside_all_candidates"], m["area_fraction_outside_all_candidates"]),
         "Partition of unity: max |sum_j phi_j - 1| = %.2e, min phi = %.2e (passed %s)." % (
             r["partition_of_unity"]["max_abs_sum_minus_1"], r["partition_of_unity"]["min_weight"],
             r["partition_of_unity"]["passed"]),
         "",
         "Construction: see the docstring of localnn_mesh.py (geodesic Voronoi / Lloyd on the facet graph, dual",
         "triangles, barycentric coordinates of each facet centroid projected on the plane of its coarse triangle).",
         "Files: grid.npz (arrays, see the docstring), grid_report.json (all diagnostics), coarse_mesh.vtu (ParaView).",
         "", "COMMAND", "  " + " ".join(r["command"])]
    with open(stage / "README.txt", "x") as f:
        f.write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True)
    ap.add_argument("--geometry", required=True, help="a validated nonneg_maps geometry (mesh_nodes.f8, mesh_faces.i4)")
    ap.add_argument("--spacing", type=float, required=True, help="target node spacing h [m]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--lloyd", type=int, default=6)
    ap.add_argument("--sample", type=int, default=400000)
    ap.add_argument("--seed", type=int, default=20261114)
    ap.add_argument("--seen-from", help="nonneg_facets.npz with the fit's seen flag (diagnostic node_seen)")
    return build(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
