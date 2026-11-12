#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_bumps.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_bumps.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3fa124d4b780af66458a80bd7163b725320b2b5705d5d29e8516f3943c19994b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Gaussian-bump basis of the H2O activity for the bump + slope-multiplier model (option 3a; user decisions
#   2026-11-10: about 32 centres, 26 on the seen and 6 on the unseen surface, 3-D Gaussians with the distance to the
#   nearest other centre as width, start = NNLS fit to the original map): deterministic centres (farthest point +
#   area-weighted Lloyd, snapped to facets), coverage statistics, NNLS start weights. User 2026-11-11: "the idea with
#   bumps does not work" (model dropped).
# USED BY: CPU job 123981 (run from the identical copy checks/slopebumps_20261110/) ->
#   runs/shap8_012m_lategrid_dz1.5/slope_bumps/bundle/bumps32.npz.
# ==== end of header ====
"""Gaussian-bump basis of the H2O activity for the bump + slope-multiplier model (user decisions 2026-11-10, option 3a).

    .venv/bin/python -W ignore slope_bumps.py --run <base> --seen <npz with key 'seen'> --out <new .npz>
        [--n-seen 26] [--n-unseen 6] [--width nn|median_nn] [--max-iter 200] [--threads 32]

User decisions (2026-11-10): A(e) = max(1e-5 A_orig(e), m_b(e) L(e)), L(e) = sum_k w_k phi_k(x_e), w_k >= 0, phi_k >= 0;
~32 centres, denser on the seen surface: 26 spread evenly over the seen facets (the fit's seen set) and 6 over the
unseen facets; phi_k(x) = exp(-|x - c_k|^2 / (2 s_k^2)) (3-D straight-line distance), s_k = the distance from c_k to
its nearest other centre (per centre) unless the coverage is poor (then: stop and report alternatives).

Centres (deterministic, no random numbers), separately in each set (seen / unseen facets, centroids x_e, areas a_e):
 1. farthest-point initialisation: the first centre is the facet nearest to the area-weighted mean of the set, every
    next one the facet farthest from the centres so far (ties: lowest facet index);
 2. area-weighted Lloyd k-means: each facet to its nearest centre, each centre to the area-weighted mean of its cluster,
    until no assignment changes (or --max-iter);
 3. each centre snapped to the facet of its cluster nearest to the cluster's area-weighted mean (centres on the surface).
Widths: s_k = min_{j != k} |c_k - c_j| over all centres (--width nn, the user's definition); --width median_nn (common
width = the median of those distances) is only the alternative named by the user for poor coverage, never a default.

Coverage C(e) = sum_k phi_k(x_e): min / percentiles over all facets and over the seen facets, area fractions below
1e-1 .. 1e-6, written to the metadata and printed. Start weights (user): non-negative least squares of L to the original
activity A_orig (elements.bin AH2O), area weighted, all facets: min sum_e a_e (sum_k u_k phi_k(e) - A_orig(e))^2, u >= 0,
solved exactly on the 32 x 32 normal equations (G = Phi^T diag(a) Phi = R^T R, scipy nnls(R, R^-T Phi^T diag(a) A_orig);
the same minimiser as the full problem). u are physical activities (L = Phi u); the fit's optimizer coordinates are
w = u / A_scale (A_scale = the fitter's enhancement scale, applied by slope_mult_problem). Start-map quality: area- and
emission-weighted rms of ln(L/A_orig) over facets where both > 0, and the facets / area at the floor at the start (all
multipliers 1: L <= 1e-5 A_orig).
"""

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import slope_mult as sm  # noqa: E402

CHUNK = 1 << 19
DEFINITION = "exp(-|x - c_k|^2 / (2 s_k^2)), x = facet centroid, 3-D straight-line distance [m]"


def nearest(X, C, pool=None):
    """(index of the nearest centre, squared distance) for every row of X (chunked; ties: lowest centre index)."""
    n = len(X)
    idx = np.empty(n, np.int64)
    d2 = np.empty(n)
    cc = np.sum(C * C, axis=1)

    def part(i0):
        i1 = min(i0 + CHUNK, n)
        x = X[i0:i1]
        D = np.sum(x * x, axis=1)[:, None] - 2.0 * (x @ C.T) + cc[None, :]
        j = np.argmin(D, axis=1)
        idx[i0:i1] = j
        diff = x - C[j]
        d2[i0:i1] = np.sum(diff * diff, axis=1)
    starts = range(0, n, CHUNK)
    if pool is None:
        for i0 in starts:
            part(i0)
    else:
        list(pool.map(part, starts))
    return idx, d2


def centres_for(X, a, K, max_iter, pool):
    """Deterministic area-weighted k-means with farthest-point initialisation; returns (local facet indices, record)."""
    mean = np.sum(X * a[:, None], axis=0) / np.sum(a)
    first = int(np.argmin(np.sum((X - mean) ** 2, axis=1)))
    chosen = [first]
    dmin = np.sum((X - X[first]) ** 2, axis=1)
    for _ in range(1, K):
        nxt = int(np.argmax(dmin))                      # argmax: lowest index among ties
        chosen.append(nxt)
        dmin = np.minimum(dmin, np.sum((X - X[nxt]) ** 2, axis=1))
    C = X[chosen].copy()
    assign = None
    it = 0
    for it in range(1, max_iter + 1):
        new, _ = nearest(X, C, pool)
        if assign is not None and np.array_equal(new, assign):
            break
        assign = new
        w = np.bincount(assign, weights=a, minlength=K)
        for d in range(3):
            C[:, d] = np.bincount(assign, weights=a * X[:, d], minlength=K) / np.where(w > 0, w, 1.0)
    # snap: per cluster the facet nearest to the cluster's weighted mean
    snapped = []
    for k in range(K):
        sel = np.flatnonzero(assign == k)
        if not len(sel):
            raise SystemExit("slope_bumps: empty cluster %d" % k)
        j = sel[int(np.argmin(np.sum((X[sel] - C[k]) ** 2, axis=1)))]
        snapped.append(int(j))
    w = np.bincount(assign, weights=a, minlength=K)
    rec = {"farthest_point_init_local": chosen, "lloyd_iterations": it, "converged": bool(it < max_iter),
           "cluster_area_m2": [float(v) for v in w], "cluster_facets": [int(v) for v in np.bincount(assign, minlength=K)]}
    return np.array(snapped, np.int64), rec


def basis(Xall, C, s, pool):
    """Phi (nfaces, K) float64 C-contiguous: phi_k = exp(-|x - c_k|^2 / (2 s_k^2))."""
    n, K = len(Xall), len(C)
    Phi = np.empty((n, K))
    inv = 1.0 / (2.0 * s * s)
    cc = np.sum(C * C, axis=1)

    def part(i0):
        i1 = min(i0 + CHUNK, n)
        x = Xall[i0:i1]
        D = np.sum(x * x, axis=1)[:, None] - 2.0 * (x @ C.T) + cc[None, :]
        np.maximum(D, 0.0, out=D)
        Phi[i0:i1] = np.exp(-D * inv[None, :])
    list(pool.map(part, range(0, n, CHUNK)))
    return Phi


def stats(v, area, mask=None):
    v = v if mask is None else v[mask]
    a = area if mask is None else area[mask]
    q = np.percentile(v, [0.01, 0.1, 1, 5, 25, 50, 75, 95, 99])
    out = {"min": float(v.min()), "max": float(v.max()), "facets": int(len(v)),
           "percentiles": {p: float(x) for p, x in zip(("0.01", "0.1", "1", "5", "25", "50", "75", "95", "99"), q)},
           "area_fraction_below": {t: float(a[v < float(t)].sum() / a.sum()) for t in
                                   ("1e-1", "1e-2", "1e-3", "1e-4", "1e-6")}}
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", required=True)
    p.add_argument("--seen", required=True, help="npz with a boolean 'seen' per facet (the fit's seen set)")
    p.add_argument("--out", required=True)
    p.add_argument("--n-seen", type=int, default=26)
    p.add_argument("--n-unseen", type=int, default=6)
    p.add_argument("--width", choices=("nn", "median_nn"), default="nn")
    p.add_argument("--max-iter", type=int, default=200)
    p.add_argument("--threads", type=int, default=32)
    args = p.parse_args(argv)
    from scipy.optimize import nnls
    from cops_cuda import formats, provenance, reference
    reference.bind()
    templates = reference.load("activity_templates")
    t0 = time.time()
    run = Path(args.run).resolve()
    el, _ = formats.load_elements(run / "elements.bin", lazy=True)
    X = np.ascontiguousarray(np.asarray(el["centroid"], np.float64))
    area = np.ascontiguousarray(np.asarray(el["area"], np.float64))
    A_orig = np.ascontiguousarray(np.asarray(el["AH2O"], np.float64))
    n = len(area)
    with np.load(args.seen, allow_pickle=False) as z:
        seen = np.asarray(z["seen"], bool)
    if seen.shape != (n,):
        raise SystemExit("slope_bumps: seen has %s values, the mesh %d" % (seen.shape, n))
    pool = ThreadPoolExecutor(args.threads)
    sets = (("seen", seen, args.n_seen), ("unseen", ~seen, args.n_unseen))
    centres, ids, which, recs = [], [], [], {}
    for name, mask, K in sets:
        loc = np.flatnonzero(mask)
        local, rec = centres_for(X[loc], area[loc], K, args.max_iter, pool)
        g = loc[local]
        ids += [int(v) for v in g]
        which += [name] * K
        centres.append(X[g])
        recs[name] = dict(rec, facets=int(len(loc)), area_fraction=float(area[loc].sum() / area.sum()))
        print("slope_bumps: %-6s %d centres, %d Lloyd iterations (converged %s), %d facets, area fraction %.4f" % (
            name, K, rec["lloyd_iterations"], rec["converged"], len(loc), recs[name]["area_fraction"]), flush=True)
    C = np.vstack(centres)
    K = len(C)
    D = np.sqrt(np.sum((C[:, None, :] - C[None, :, :]) ** 2, axis=2))
    np.fill_diagonal(D, np.inf)
    nn = D.min(axis=1)
    s = nn.copy() if args.width == "nn" else np.full(K, float(np.median(nn)))
    Phi = basis(X, C, s, pool)
    cov = np.empty(n)

    def csum(i0):
        i1 = min(i0 + CHUNK, n)
        cov[i0:i1] = Phi[i0:i1].sum(axis=1)
    list(pool.map(csum, range(0, n, CHUNK)))
    cov_all, cov_seen, cov_unseen = stats(cov, area), stats(cov, area, seen), stats(cov, area, ~seen)
    print("slope_bumps: widths (%s) min %.1f median %.1f max %.1f m; coverage sum_k phi_k: all min %.3e p0.1 %.3e p1 "
          "%.3e p5 %.3e median %.3e; seen min %.3e p0.1 %.3e p1 %.3e; unseen min %.3e p1 %.3e; area below 1e-2: all "
          "%.2e seen %.2e; below 1e-3: all %.2e seen %.2e" % (
              args.width, s.min(), np.median(s), s.max(), cov_all["min"], cov_all["percentiles"]["0.1"],
              cov_all["percentiles"]["1"], cov_all["percentiles"]["5"], cov_all["percentiles"]["50"],
              cov_seen["min"], cov_seen["percentiles"]["0.1"], cov_seen["percentiles"]["1"], cov_unseen["min"],
              cov_unseen["percentiles"]["1"], cov_all["area_fraction_below"]["1e-2"],
              cov_seen["area_fraction_below"]["1e-2"], cov_all["area_fraction_below"]["1e-3"],
              cov_seen["area_fraction_below"]["1e-3"]), flush=True)
    # NNLS start (area weighted, all facets), on the normal equations
    Wphi = Phi * area[:, None]
    G = Phi.T @ Wphi
    c = Wphi.T @ A_orig
    del Wphi
    evals = np.linalg.eigvalsh(G)
    Rm = np.linalg.cholesky(G).T                                   # G = R^T R
    d = np.linalg.solve(Rm.T, c)
    u, rnorm = nnls(Rm, d, maxiter=50 * K)
    L0 = sm.chunked_matvec(Phi, u, pool, args.threads)
    F = sm.FLOOR_RATIO * A_orig
    both = (L0 > 0) & (A_orig > 0)
    lr = np.log(L0[both] / A_orig[both])
    wa, we = area[both], area[both] * A_orig[both]
    at_floor = ~(L0 > F)
    obj = float(np.sum(area * (L0 - A_orig) ** 2))
    quality = {"rms_ln_L_over_A_orig_area_weighted": float(np.sqrt(np.sum(wa * lr * lr) / wa.sum())),
               "rms_ln_L_over_A_orig_emission_weighted": float(np.sqrt(np.sum(we * lr * lr) / we.sum())),
               "facets_compared": int(both.sum()),
               "facets_at_floor_start_m1": int(at_floor.sum()),
               "facets_at_floor_fraction": float(at_floor.mean()),
               "area_at_floor_fraction": float(area[at_floor].sum() / area.sum()),
               "seen_area_at_floor_fraction": float(area[at_floor & seen].sum() / area[seen].sum()),
               "total_activity_ratio": float(np.sum(area * L0) / np.sum(area * A_orig)),
               "weights_zero": int(np.sum(u == 0.0)),
               "objective_area_weighted_sq": obj,
               "relative_rms_residual": float(np.sqrt(obj / np.sum(area * A_orig ** 2)))}
    print("slope_bumps: NNLS start: %d of %d weights zero; rms ln(L/A_orig) area %.3f emission %.3f nepers; floor at "
          "the start (m = 1): %.4f of the facets, %.4f of the area, %.4f of the seen area; total activity ratio %.4f; "
          "relative rms residual %.3f; cond(G) %.2e" % (
              quality["weights_zero"], K, quality["rms_ln_L_over_A_orig_area_weighted"],
              quality["rms_ln_L_over_A_orig_emission_weighted"], quality["facets_at_floor_fraction"],
              quality["area_at_floor_fraction"], quality["seen_area_at_floor_fraction"],
              quality["total_activity_ratio"], quality["relative_rms_residual"], evals.max() / evals.min()), flush=True)
    names = ["bump_%02d" % k for k in range(K)]
    per = []
    dom, _ = nearest(X, C, pool)          # nearest centre (for the record; Gaussian dominance differs by widths)
    for k in range(K):
        per.append({"name": names[k], "set": which[k], "facet": ids[k], "centre_m": [float(v) for v in C[k]],
                    "width_m": float(s[k]), "nearest_centre_distance_m": float(nn[k]),
                    "nearest_centre_area_m2": float(area[dom == k].sum()),
                    "effective_area_m2": float(np.sum(area * Phi[:, k])), "nnls_u": float(u[k])})
    mesh = templates.mesh_identity(X, area)
    meta = {
        "created_utc": provenance.utc_now(), "base_run": str(run),
        "elements": provenance.file_identity(run / "elements.bin", "elements"),
        "seen_source": provenance.file_identity(args.seen, "seen_mask"),
        "templates": {nm: {"definition": "Gaussian bump %d of %d: %s, centre facet %d (%s set), width %.3f m" % (
            k, K, DEFINITION, ids[k], which[k], s[k]),
            "interpretation": "activity basis function phi_k of L = sum_k w_k phi_k (bump + slope-multiplier model, "
                              "slope_mult.py); its fitted amplitude is the weight w_k >= 0 in units of the fitter's "
                              "A_scale", "policy": "user 2026-11-10 option 3a"} for k, nm in enumerate(names)},
        "slope_bumps": {"names": names, "n_seen": args.n_seen, "n_unseen": args.n_unseen, "definition": DEFINITION,
                        "width_rule": args.width,
                        "method": "farthest-point initialisation + area-weighted Lloyd k-means per set, centres "
                                  "snapped to the nearest facet of their cluster; deterministic (no random seed)",
                        "sets": recs, "bumps": per, "widths_m": [float(v) for v in s],
                        "centres_m": C.tolist(), "centre_facets": ids, "centre_sets": which,
                        "coverage": {"all": cov_all, "seen": cov_seen, "unseen": cov_unseen},
                        "nnls_start": {"u": [float(v) for v in u], "rnorm_reduced": float(rnorm),
                                       "gram_eigenvalues": [float(v) for v in evals],
                                       "definition": "min sum_e area_e (Phi u - A_orig)^2, u >= 0 (all facets), on the "
                                                     "normal equations (Cholesky R, scipy nnls)",
                                       "quality": quality},
                        "bumps_sha256": sm.sha256_bytes(Phi), "nfaces": int(n), "mesh_identity": mesh,
                        "ah2o_sha256": sm.sha256_bytes(A_orig), "seen_sha256": sm.sha256_bytes(seen),
                        "seconds": time.time() - t0},
        "writer": provenance.file_identity(Path(__file__), "tool"),
        "writer_module": provenance.file_identity(Path(sm.__file__), "tool"),
    }
    templates.write_bundle(args.out, names, Phi, X, area, meta)
    Path(str(args.out)[:-4] + ".json").write_text(json.dumps(meta, indent=2) + "\n")
    print("slope_bumps: wrote %s in %.0f s" % (args.out, time.time() - t0), flush=True)
    pool.shutdown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
