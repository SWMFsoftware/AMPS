#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_start.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_start.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 234d3245131155fc37a99e0130693df349167be261aaec9ae92890daa3f9e456
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Start of the squared-harmonics model (user decision 2026-11-11: "squared: b by least squares so the
#   activity matches D2's max(0, field)"; area weighted, all facets): linear LSQ start, then deterministic
#   Levenberg-Marquardt in fixed 2^18-row chunks. Bitwise only on one CPU type (numpy's OpenBLAS is DYNAMIC_ARCH:
#   node023 vs node032 differ at 8e-15, jobs 124149-124153).
# USED BY: jobs/nonneg_start.sbatch (test 124129; production 124138 ->
#   runs/shap8_012m_lategrid_dz1.5/nonneg/bundle/squared_start.npz), jobs/nonneg_start_determinism.sbatch.
# 2026-11-14: --model sos --rank r (Gram LSQ, PSD projection, top-r; job 124994), --refine lm (not chosen).
# ==== end of header ====
"""Start of the squared-harmonics model (user decision 2026-11-11: "squared: b by least squares so the activity
matches D2's max(0, field)"; the lead: area-weighted, like the 3a NNLS start, all facets).

    nonneg_start.py --run <base> --d2 <D2 final fit_restart.npz> --ft-templates-from <D2 fit_manifest.json>
        --out <start.npz> [--threads 32] [--max-iter 300]

Problem (dimensionless, everything divided by A_s = D2's enhancement scale = the fitter's):
    minimise  sum_e area_e ( (Y b)_e^2 + S_e - t_e )^2   over b (25),
    S = sum_k T_k (D2's ft_c* templates at amplitude 1, as the model holds them), t = max(0, D2's linear field) / A_s
    (D2's activity; the linear field = activity_linear of D2's fit_restart.npz, Y a + A_s (c_cliff T_cliff + sum T_k)).
The model's harmonic coefficients are then a = A_s b (A(m = 1) = A_s ((Y b)^2 + S) = (Y a)^2 / A_s + A_s S).
Method (deterministic): b0 = the area-weighted linear least squares of Y b0 ~ sqrt(max(t - S, 0)) (the positive root:
b and -b give the same activity), then Levenberg-Marquardt on the normal equations (Gauss-Newton matrix
sum area 4 (Y b)^2 Y Y^T, damping mu diag), accepted steps only, until the relative cost decrease of an accepted step
is below 1e-13, the step below 1e-12 relative, or --max-iter. Sums in fixed 2^18-row chunks (thread-count independent).
Writes <out> (b, a, meta_json: inputs with sha256, A_s, weighting, quality, iteration history) and <out>.json.

--model sos --rank r (option A, lead 2026-11-14): the Gram matrix Q' (25 x 25 symmetric, dimensionless) of
    minimise  sum_e area_e ( y_e^T Q' y_e + S_e - t_e )^2   (linear in the 325 entries of Q'),
solved on the normal equations; the products y_i y_j span only the degree <= 8 harmonics (81 functions), so the system
is rank deficient and the MINIMUM-NORM solution is taken (eigenvalues of the normal matrix below 1e-12 x the largest
dropped). Then the PSD projection (negative eigenvalues of Q' set to 0), the top-r eigenpairs M = U_r sqrt(lambda_r), and
the lower-trapezoidal L' with L' L'^T = M M^T (QR of M^T, diagonal >= 0); the model's L = A_s L' (A_s ||L'^T y||^2 =
||L^T y||^2 / A_s). Writes <out> (L, meta_json: quality of the full LSQ Q', of the projected rank-r Q, the spectrum).
--refine lm (sos only, PROPOSAL 2026-11-14, pending the lead's OK; default none = the start above): Levenberg-Marquardt
    on the 25 + 24 + ... lower-trapezoidal entries of L' of the rank-r problem itself,
        minimise  sum_e area_e ( ||L'^T y_e||^2 + S_e - t_e )^2,
    from the rank-r PSD start above (the same LM as the squared start: Gauss-Newton matrix, damping mu diag, accepted
    steps only, the same stop rules, fixed 2^18-row chunks); the squared model is its r = 1 case. The diagonal of L'[:r,:r]
    is made >= 0 at the end (column sign; the activity is unchanged). The output records both qualities.
Every mode records the ft templates' minimum and their facet counts (> 0).
"""

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
for _p in (str(HERE), str(APP / "tools"), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CHUNK = 1 << 18


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def log(msg):
    print("nonneg_start: %s" % msg, flush=True)


def chunks(n):
    return [(c0, min(c0 + CHUNK, n)) for c0 in range(0, n, CHUNK)]


def reduce_chunks(pool, fn, n):
    parts = list(pool.map(lambda c: fn(*c), chunks(n))) if pool is not None else [fn(*c) for c in chunks(n)]
    out = parts[0]
    for p in parts[1:]:
        out = tuple(x + y for x, y in zip(out, p))
    return out


def solve(Y, S, t, w, pool, max_iter, record):
    n = Y.shape[0]
    y0 = np.sqrt(np.maximum(t - S, 0.0))

    def lin(c0, c1):
        Yc, wc = Y[c0:c1], w[c0:c1]
        return (Yc.T @ (Yc * wc[:, None]), Yc.T @ (wc * y0[c0:c1]))
    G, g = reduce_chunks(pool, lin, n)
    b = np.linalg.solve(G, g)

    def cost_of(b_):
        def part(c0, c1):
            v = Y[c0:c1] @ b_
            r = v * v + S[c0:c1] - t[c0:c1]
            return (np.array([np.sum(w[c0:c1] * r * r)]),)
        return float(reduce_chunks(pool, part, n)[0][0])

    def normal(b_):
        def part(c0, c1):
            Yc, wc = Y[c0:c1], w[c0:c1]
            v = Yc @ b_
            r = v * v + S[c0:c1] - t[c0:c1]
            jr = 2.0 * v
            return (Yc.T @ (Yc * (wc * jr * jr)[:, None]), Yc.T @ (wc * jr * r), np.array([np.sum(wc * r * r)]))
        JJ, Jr, c = reduce_chunks(pool, part, n)
        return JJ, Jr, float(c[0])

    c0 = cost_of(b)
    hist = [{"iteration": 0, "cost": c0, "mu": None, "accepted": True}]
    log("b0 (linear LSQ to sqrt(max(t - S, 0))): cost %.12e" % c0)
    mu, converged, reason = 1e-3, False, "max_iter"
    JJ, Jr, cost = normal(b)
    for it in range(1, max_iter + 1):
        d = np.diag(JJ).copy()
        d[d <= 0] = 1.0
        step = np.linalg.solve(JJ + mu * np.diag(d), -Jr)
        trial = b + step
        ct = cost_of(trial)
        acc = ct < cost
        hist.append({"iteration": it, "cost": ct, "mu": mu, "accepted": bool(acc),
                     "step_rel": float(np.linalg.norm(step) / max(np.linalg.norm(b), 1e-300))})
        if acc:
            rel = (cost - ct) / cost
            b = trial
            JJ, Jr, cost = normal(b)
            mu = max(mu / 3.0, 1e-12)
            if rel < 1e-13 or hist[-1]["step_rel"] < 1e-12:
                converged, reason = True, "relative cost decrease %.2e, relative step %.2e" % (rel, hist[-1]["step_rel"])
                break
        else:
            mu *= 4.0
            if mu > 1e16:
                converged, reason = True, "no further decrease (mu %.1e)" % mu
                break
        if it % 10 == 0:
            log("iteration %d: cost %.12e (mu %.1e)" % (it, cost, mu))
    record.update(iterations=len(hist) - 1, converged=bool(converged), stop=reason, history=hist,
                  gradient_norm=float(np.linalg.norm(Jr)), final_cost=cost, start_cost=c0)
    return b


def solve_sos(Y, S, t, w, Lp0, pool, max_iter, record):
    """LM of the rank-r problem sum_e w_e (||L'^T y_e||^2 + S_e - t_e)^2 over the lower-trapezoidal entries of L' (column
    j: rows j..24), from Lp0; the same stopping rules as solve() (the squared start, its r = 1 case)."""
    n, r = Y.shape[0], Lp0.shape[1]
    idx = [(i, j) for j in range(r) for i in range(j, 25)]
    rows = np.array([i for i, _ in idx]); cols = np.array([j for _, j in idx])

    def vec(L_):
        return L_[rows, cols].copy()

    def mat(p_):
        L_ = np.zeros((25, r)); L_[rows, cols] = p_
        return L_

    def cost_of(p_):
        L_ = mat(p_)
        def part(c0, c1):
            V = Y[c0:c1] @ L_
            res = np.sum(V * V, axis=1) + S[c0:c1] - t[c0:c1]
            return (np.array([np.sum(w[c0:c1] * res * res)]),)
        return float(reduce_chunks(pool, part, n)[0][0])

    def normal(p_):
        L_ = mat(p_)
        def part(c0, c1):
            Yc, wc = Y[c0:c1], w[c0:c1]
            V = Yc @ L_
            res = np.sum(V * V, axis=1) + S[c0:c1] - t[c0:c1]
            Jc = 2.0 * V[:, cols] * Yc[:, rows]
            return (Jc.T @ (Jc * wc[:, None]), Jc.T @ (wc * res), np.array([np.sum(wc * res * res)]))
        JJ, Jr, c = reduce_chunks(pool, part, n)
        return JJ, Jr, float(c[0])

    p = vec(Lp0)
    c0 = cost_of(p)
    hist = [{"iteration": 0, "cost": c0, "mu": None, "accepted": True}]
    log("sos LM: %d parameters; start (PSD rank-%d) cost %.12e" % (len(p), r, c0))
    mu, converged, reason = 1e-3, False, "max_iter"
    JJ, Jr, cost = normal(p)
    for it in range(1, max_iter + 1):
        d = np.diag(JJ).copy()
        d[d <= 0] = 1.0
        step = np.linalg.solve(JJ + mu * np.diag(d), -Jr)
        trial = p + step
        ct = cost_of(trial)
        acc = ct < cost
        hist.append({"iteration": it, "cost": ct, "mu": mu, "accepted": bool(acc),
                     "step_rel": float(np.linalg.norm(step) / max(np.linalg.norm(p), 1e-300))})
        if acc:
            rel = (cost - ct) / cost
            p = trial
            JJ, Jr, cost = normal(p)
            mu = max(mu / 3.0, 1e-12)
            if rel < 1e-13 or hist[-1]["step_rel"] < 1e-12:
                converged, reason = True, "relative cost decrease %.2e, relative step %.2e" % (rel, hist[-1]["step_rel"])
                break
        else:
            mu *= 4.0
            if mu > 1e16:
                converged, reason = True, "no further decrease (mu %.1e)" % mu
                break
        if it % 10 == 0:
            log("sos LM iteration %d: cost %.12e (mu %.1e)" % (it, cost, mu))
    Lp = mat(p)
    sgn = np.sign(np.diag(Lp[:r, :r]))
    sgn[sgn == 0] = 1.0
    Lp = Lp * sgn[None, :]
    record.update(iterations=len(hist) - 1, converged=bool(converged), stop=reason, history=hist,
                  gradient_norm=float(np.linalg.norm(Jr)), final_cost=cost, start_cost=c0, parameters=len(p))
    return Lp


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--d2", required=True, help="D2's final fit_restart.npz")
    ap.add_argument("--ft-templates-from", required=True, help="D2's fit_manifest.json (or fit directory)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-iter", type=int, default=300)
    ap.add_argument("--model", default="squared", choices=("squared", "sos"))
    ap.add_argument("--rank", type=int, default=None, help="sos: the rank r of L")
    ap.add_argument("--refine", default="none", choices=("none", "lm"),
                    help="sos: none (the PSD rank-r start) | lm (Levenberg-Marquardt of the rank-r problem from it)")
    args = ap.parse_args(argv)
    if args.model == "sos" and not (args.rank and 1 <= args.rank <= 25):
        ap.error("--model sos needs --rank 1..25")
    if args.refine != "none" and args.model != "sos":
        ap.error("--refine is a sos option")
    t_start = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("exists: %s" % out)
    import nonneg_problem as npb
    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    at = reference.load("activity_templates")
    run = Path(args.run).resolve()
    el = fit_cops.read_elements(str(run / "elements.bin"))
    area = np.asarray(el["area"], np.float64)
    Y = np.ascontiguousarray(el["Y"], dtype=np.float64)
    mesh = at.mesh_identity(el["centroid"], el["area"])
    with np.load(args.d2, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        field = np.asarray(z["activity_linear"], np.float64)
    scale = float(meta["enhancement_scale"])
    from cops_cuda import restart as rs
    st = rs.read(args.d2)
    if st.mesh_identity() != mesh:
        raise SystemExit("D2's restart mesh %s is not the base's %s" % (st.mesh_identity()[:16], mesh[:16]))
    if field.shape != (len(area),):
        raise SystemExit("D2's activity_linear has shape %s" % (field.shape,))
    fts, ftsrc = npb.ft_templates_from(args.ft_templates_from)
    cols = []
    for t in fts:
        bundle = at.load_bundle(t["bundle"]["path"], el["centroid"], area)
        cols.append(at.select_templates(bundle, [t["name"]]))
    T = np.ascontiguousarray(np.column_stack(cols), dtype=np.float64)
    S = T @ np.ones(T.shape[1])
    t = np.maximum(field, 0.0) / scale
    w = area / area.sum()
    pool = ThreadPoolExecutor(args.threads) if args.threads > 1 else None
    tstats = {"min": float(T.min()), "facets_positive": [int(v) for v in np.sum(T > 0, axis=0)],
              "facets_negative": [int(v) for v in np.sum(T < 0, axis=0)], "names": [x["name"] for x in fts]}
    log("ft templates %s: min %.3e, facets > 0 %s, facets < 0 %s" % (tstats["names"], tstats["min"],
                                                                     tstats["facets_positive"], tstats["facets_negative"]))
    if args.model == "sos":
        return main_sos(args, out, Y, S, t, w, area, scale, mesh, run, fts, ftsrc, tstats, pool, t_start)
    rec = {}
    log("base %s: %d facets; D2 %s (A_s %.9e); templates %s; target t = max(0, D2 field) / A_s: area with t > 0 %.4f; "
        "S min %.3e max %.3e" % (run, len(area), args.d2, scale, [x["name"] for x in fts],
                                 float(area[t > 0].sum() / area.sum()), float(S.min()), float(S.max())))
    b = solve(Y, S, t, w, pool, args.max_iter, rec)
    v = Y @ b
    model = v * v + S
    r = model - t
    q = {"relative_rms_residual": float(np.sqrt(np.sum(w * r * r) / np.sum(w * t * t))),
         "total_activity_ratio": float(np.sum(w * model) / np.sum(w * t)),
         "area_target_positive": float(np.sum(w[t > 0])), "area_target_zero": float(np.sum(w[t <= 0])),
         "mean_model_over_mean_target_where_target_zero": (float(np.sum(w[t <= 0] * model[t <= 0]) /
                                                            np.sum(w * t) * np.sum(w)) if np.any(t <= 0) else None),
         "corr_area_weighted": float(np.sum(w * (model - np.sum(w * model)) * (t - np.sum(w * t))) /
                                     np.sqrt(np.sum(w * (model - np.sum(w * model)) ** 2) *
                                             np.sum(w * (t - np.sum(w * t)) ** 2))),
         "templates_share_of_model": float(np.sum(w * S) / np.sum(w * model)),
         "iterations": rec["iterations"], "converged": rec["converged"], "stop": rec["stop"],
         "cost_start": rec["start_cost"], "cost_final": rec["final_cost"], "gradient_norm": rec["gradient_norm"]}
    smeta = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "tool": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
             "run": str(run), "elements": {"path": str(run / "elements.bin"), "sha256": sha256(run / "elements.bin")},
             "mesh_identity": mesh, "d2_restart": {"path": str(Path(args.d2).resolve()), "sha256": sha256(args.d2)},
             "A_scale": scale, "ft_templates": fts, "ft_source": ftsrc,
             "problem": "minimise sum_e area_e ((Y b)_e^2 + S_e - t_e)^2; S = sum_k T_k (amplitude 1); t = max(0, D2 "
                        "activity_linear) / A_s; all facets; a = A_s b",
             "weighting": "area (all facets), normalised by the total area (the lead 2026-11-11: like the 3a NNLS start)",
             "method": "b0 = area-weighted linear LSQ of Y b0 ~ sqrt(max(t - S, 0)); Levenberg-Marquardt on the normal "
                       "equations, fixed 2^18-row chunks",
             "quality": q, "history": rec["history"], "threads": args.threads, "seconds": time.time() - t_start,
             "ft_template_stats": tstats}
    with out.open("xb") as f:
        np.savez(f, b=b, a=scale * b, meta_json=np.array(json.dumps(smeta)))
    Path(str(out) + ".json").write_text(json.dumps({k: v for k, v in smeta.items() if k != "history"}, indent=1) + "\n")
    log("b = %s" % np.array2string(b, precision=6, max_line_width=200))
    log("quality %s; %.1f s; wrote %s" % ({k: q[k] for k in ("relative_rms_residual", "total_activity_ratio",
                                                             "corr_area_weighted", "templates_share_of_model",
                                                             "iterations", "converged", "stop")},
                                          time.time() - t_start, out))
    return 0


def main_sos(args, out, Y, S, t, w, area, scale, mesh, run, fts, ftsrc, tstats, pool, t_start):
    r = int(args.rank)
    iu, ju = np.triu_indices(25)
    fac = np.where(iu == ju, 1.0, 2.0)
    CH = 1 << 15

    def part(c0, c1):
        Yc = Y[c0:c1]
        ph = Yc[:, iu] * Yc[:, ju] * fac
        wc = w[c0:c1]
        return (ph.T @ (ph * wc[:, None]), ph.T @ (wc * (t[c0:c1] - S[c0:c1])))
    chunks_ = [(c0, min(c0 + CH, len(w))) for c0 in range(0, len(w), CH)]
    parts = list(pool.map(lambda c: part(*c), chunks_)) if pool is not None else [part(*c) for c in chunks_]
    G, g = parts[0]
    for p_ in parts[1:]:
        G, g = G + p_[0], g + p_[1]
    ev, EV = np.linalg.eigh(0.5 * (G + G.T))
    keep = ev > 1e-12 * ev.max()
    qv = EV[:, keep] @ ((EV[:, keep].T @ g) / ev[keep])
    Qf = np.zeros((25, 25))
    Qf[iu, ju] = qv
    Qf[ju, iu] = qv
    lam, U = np.linalg.eigh(Qf)
    order = np.argsort(lam)[::-1]
    lam, U = lam[order], U[:, order]
    pos = np.maximum(lam, 0.0)
    M = U[:, :r] * np.sqrt(pos[:r])[None, :]
    Qm, Rm = np.linalg.qr(M.T)               # M^T = Qm Rm, Rm (r x 25) upper-trapezoidal
    Lp = Rm.T                                # 25 x r lower-trapezoidal, Lp Lp^T = M M^T
    sgn = np.sign(np.diag(Lp[:r, :r]))
    sgn[sgn == 0] = 1.0
    Lp = Lp * sgn[None, :]
    Lp[np.triu_indices(25, 1, r)] = 0.0       # exactly lower-trapezoidal
    Lp_psd = Lp.copy()
    lmrec = None
    if args.refine == "lm":
        lmrec = {}
        Lp = solve_sos(Y, S, t, w, Lp_psd, pool, args.max_iter, lmrec)
        Lp[np.triu_indices(25, 1, r)] = 0.0
    L = scale * Lp

    def quality(model):
        res = model - t
        return {"relative_rms_residual": float(np.sqrt(np.sum(w * res * res) / np.sum(w * t * t))),
                "total_activity_ratio": float(np.sum(w * model) / np.sum(w * t)),
                "corr_area_weighted": float(np.sum(w * (model - np.sum(w * model)) * (t - np.sum(w * t))) /
                                            np.sqrt(np.sum(w * (model - np.sum(w * model)) ** 2) *
                                                    np.sum(w * (t - np.sum(w * t)) ** 2))),
                "min": float(model.min())}
    full = np.empty(len(w))
    rank_r = np.empty(len(w))
    rank_lm = np.empty(len(w)) if lmrec is not None else None
    for c0, c1 in chunks_:
        Yc = Y[c0:c1]
        full[c0:c1] = np.einsum("ij,jk,ik->i", Yc, Qf, Yc) + S[c0:c1]
        V = Yc @ Lp_psd
        rank_r[c0:c1] = np.sum(V * V, axis=1) + S[c0:c1]
        if rank_lm is not None:
            V = Yc @ Lp
            rank_lm[c0:c1] = np.sum(V * V, axis=1) + S[c0:c1]
    qf, qr = quality(full), quality(rank_r)
    ql = quality(rank_lm) if rank_lm is not None else None
    if ql is not None:
        ql.update({k: lmrec[k] for k in ("iterations", "converged", "stop", "start_cost", "final_cost",
                                         "gradient_norm", "parameters")})
    q = {"rank": r, "lsq_full_Q": qf, "rank_r_psd": qr, "refine": args.refine,
         "relative_rms_residual": (ql if ql is not None else qr)["relative_rms_residual"],
         "eigenvalues_Q": [float(v) for v in lam], "negative_eigenvalue_mass": float(-lam[lam < 0].sum()),
         "eigenvalue_kept_fraction": float(pos[:r].sum() / pos.sum()) if pos.sum() > 0 else None,
         "normal_matrix_rank": int(keep.sum()), "normal_matrix_size": int(len(ev))}
    if ql is not None:
        q["rank_r_lm"] = ql
    smeta = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "tool": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)}, "model": "sos",
             "run": str(run), "elements": {"path": str(run / "elements.bin"), "sha256": sha256(run / "elements.bin")},
             "mesh_identity": mesh, "d2_restart": {"path": str(Path(args.d2).resolve()), "sha256": sha256(args.d2)},
             "A_scale": scale, "ft_templates": fts, "ft_source": ftsrc, "ft_template_stats": tstats,
             "problem": "minimise sum_e area_e (y_e^T Q' y_e + S_e - t_e)^2 over symmetric Q' (325 entries; minimum-norm "
                        "solution: rank-deficient normal matrix); PSD projection; top-r eigenpairs; L' lower-trapezoidal "
                        "(QR); L = A_s L'" + ("; then Levenberg-Marquardt of the rank-r problem sum_e area_e "
                                              "(||L'^T y_e||^2 + S_e - t_e)^2 over L' from that start (--refine lm, "
                                              "proposal pending the lead's OK)" if lmrec is not None else ""),
             "refine": args.refine,
             "weighting": "area (all facets), normalised by the total area (as the squared start)",
             "quality": q, "threads": args.threads, "seconds": time.time() - t_start}
    with out.open("xb") as f:
        extra = {"Lprime_psd": Lp_psd} if lmrec is not None else {}
        np.savez(f, L=L, Lprime=Lp, Q_lsq=Qf, meta_json=np.array(json.dumps(
            dict(smeta, lm_history=lmrec["history"]) if lmrec is not None else smeta)), **extra)
    Path(str(out) + ".json").write_text(json.dumps(smeta, indent=1) + "\n")
    log("sos rank %d: normal matrix rank %d of %d; Q eigenvalues (top 6) %s; kept fraction %.4f, negative mass %.3e"
        % (r, q["normal_matrix_rank"], q["normal_matrix_size"], ["%.3e" % v for v in lam[:6]],
           q["eigenvalue_kept_fraction"] or float("nan"), q["negative_eigenvalue_mass"]))
    log("quality: full LSQ Q %s; rank-%d PSD %s%s; %.1f s; wrote %s" % (
        qf, r, qr, ("; rank-%d LM %s" % (r, ql)) if ql is not None else "", time.time() - t_start, out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
