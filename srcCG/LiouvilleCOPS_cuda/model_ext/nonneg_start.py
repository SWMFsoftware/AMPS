#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_start.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_start.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: aa1879e6d0def3aa0f44a9a069111d64462b5de9ec31c68b0efb37ac1e25e1a7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Start of the squared-harmonics model (user decision 2026-11-11: "squared: b by least squares so the
#   activity matches D2's max(0, field)"; area weighted, all facets): linear LSQ start, then deterministic
#   Levenberg-Marquardt in fixed 2^18-row chunks. Bitwise only on one CPU type (numpy's OpenBLAS is DYNAMIC_ARCH:
#   node023 vs node032 differ at 8e-15, jobs 124149-124153).
# USED BY: jobs/nonneg_start.sbatch (test 124129; production 124138 ->
#   runs/shap8_012m_lategrid_dz1.5/nonneg/bundle/squared_start.npz), jobs/nonneg_start_determinism.sbatch.
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--d2", required=True, help="D2's final fit_restart.npz")
    ap.add_argument("--ft-templates-from", required=True, help="D2's fit_manifest.json (or fit directory)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-iter", type=int, default=300)
    args = ap.parse_args(argv)
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
             "quality": q, "history": rec["history"], "threads": args.threads, "seconds": time.time() - t_start}
    with out.open("xb") as f:
        np.savez(f, b=b, a=scale * b, meta_json=np.array(json.dumps(smeta)))
    Path(str(out) + ".json").write_text(json.dumps({k: v for k, v in smeta.items() if k != "history"}, indent=1) + "\n")
    log("b = %s" % np.array2string(b, precision=6, max_line_width=200))
    log("quality %s; %.1f s; wrote %s" % ({k: q[k] for k in ("relative_rms_residual", "total_activity_ratio",
                                                             "corr_area_weighted", "templates_share_of_model",
                                                             "iterations", "converged", "stop")},
                                          time.time() - t_start, out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
