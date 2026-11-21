#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_ft.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_ft.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f5b0db3e7e0b654ee95c49d3e93ab9faf6a9dcb9b4c26a07a7bfff589a104d97
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: CPU fine tuning of a localnn fit (option C stage 2): the frozen nonneg_finetune.run (LOCALNN_NONNEG_CODE,
#   default model_ext_frozen/nonneg_20261113T035423Z; unchanged) with a CpuLateEngine (forward / seen /
#   sensitivities / forward_override on the late kernels) and LocalnnFitContext (A reproduced bitwise, held_prev, p0
#   = all-row surface + G_late (A_fit - A_noheld)); verify-softplus (vs the GPU fine tuning of the softplus chain,
#   optional --perturb); run (localnn_ft.json); run-nonneg (2026-11-15: the fine tuning of a nonneg chain fit, e.g.
#   B / D, with the real NonnegFitContext and the fit-state pressures of localnn_late.py --state; control vs the GPU
#   c5 fine tunings of B / D rel 6.3e-14 / 4.4e-16).
# RUN 2026-11-16: verification 131977-131979; tests 131980, 132009; every cycle of the stage-2 chains; B / D control
#   134501; fine-tuning regularization tests (frozen localnn_s2ftreg_20261116T031038Z).
# ==== end of header ====
"""Option C stage 2 (lead / user 2026-11-15): the late-facet FINE TUNING of a localnn fit on CPU, with the frozen fine
tuning of the nonneg chains run UNCHANGED (model_ext_frozen/nonneg_20261113T035423Z: nonneg_finetune.run,
positive_finetune.Rows / Objective / Jacobian / optimize / late_only_facets): only its engine and its fit context are
replaced.

The fine tuning (user decisions of the softplus chains, kept): per-facet factor exp(delta_e) on the facets seen only in the
last WINDOW_MIN (5) min (late_only), no lower bound on delta, grain A, B refitted, lambda 0, MAXNFEV 500, TRF + LSMR
(positive_finetune.optimize), validation checks 1-8; the next cycle's fit holds h_next = h + delta:
    A_e = m_b(e) (sum_j c_j phi_j(e) + A_s sum_k ft_ck(e)) exp(h_e)        (localnn_problem.py, held=)
ENGINE (CpuLateEngine, the interface of nonneg_finetune.GpuModelEngine): forward() = the fit state's instantaneous surface
pressures of every base row; seen(row mask, structural) = the precomputed answers of localnn_late.py (exactly the masks
late_only_facets asks for); sensitivities(facets) = G = dp/dA of the late facets (localnn_late.py, the frozen record
arithmetic); forward_override(facets, g) = p_fit + G (A_base g - A_fit) (exact: the pressures are linear in A).
CONTEXT (LocalnnFitContext, the attributes of nonneg_finetune.NonnegFitContext that run() reads): the localnn fit's
state (localnn_state.npz, nonneg_facets.npz with A = m (u + TS) exp(h)), D2's selection, the late layout and rows of
positive_finetune.Rows (localnn_late.build_rows); the fit's data and crossing costs (its smoothness and multiplier
penalties do not depend on delta and are not part of the fine tuning's objective, as the nonneg fine tuning has no
regularization: fit_costs regularization 0). The key "fit_restart.npz" of ctx.identity holds the identity of the
localnn_state.npz (nonneg_finetune.run records that key).

Subcommands
  verify-softplus --fit <softplus c<k>/fit/final> --gpu-ft <c<k>/ft/final> --late LK --state-name NAME --out DIR
      the softplus fine tuning of cycle k re-run on CPU (the real NonnegFitContext, p_fit of the state NAME of LK) and
      compared with the GPU run: late / parameter facets, delta, grain, costs, nfev, termination -> DIR/verify.json
  run --fit <localnn fit seg or final> --late LK --design DESIGN --base BASE --d2 D2 --point-table PT --out DIR
      [nonneg_finetune options: --facets late_only --window-min 5 --max-nfev 500 --validate ...]
  run-nonneg --fit <a nonneg chain fit segment> --late LK --state-name NAME --out DIR [nonneg_finetune options]
      (2026-11-15, the fine-tuning regularization tests of the lead / user ~21:10 CST) the fine tuning of a nonneg
      fit (softplus / squared / sos / background / logspread: the REAL NonnegFitContext of the nonneg code in
      LOCALNN_NONNEG_CODE) on CPU, p_fit = the state NAME of LK (localnn_late.py --state of the fit's
      nonneg_facets.npz A); localnn_ft.json records engine, late kernels and state
The nonneg code is model_ext_frozen/nonneg_20261113T035423Z unless LOCALNN_NONNEG_CODE names another frozen snapshot.
"""

import argparse
import datetime as dt
import json
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NONNEG_CODE = Path(os.environ.get("LOCALNN_NONNEG_CODE",
                                  "/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_20261113T035423Z"))
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for _p in (str(APP), str(APP / "tools"), str(HERE), str(NONNEG_CODE)):
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)
import nonneg_finetune as nft      # noqa: E402  (frozen, unchanged)
import positive_finetune as pf     # noqa: E402
for _m in (nft, pf):
    if Path(_m.__file__).resolve().parent != NONNEG_CODE.resolve():
        raise SystemExit("%s imported from %s, not from %s" % (_m.__name__, _m.__file__, NONNEG_CODE))
import localnn_late as llate        # noqa: E402
import localnn_problem as lp       # noqa: E402

LOG = "localnn_ft:"


def say(text):
    print("%s [%s] %s" % (LOG, dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


class CpuLateEngine:
    def __init__(self, LK, p0_ng, p0_rg, A_fit, held_prev, label):
        self.LK = LK
        self.p0_ng, self.p0_rg = np.asarray(p0_ng, float), np.asarray(p0_rg, float)
        self.A_fit, self.held_prev = np.asarray(A_fit, float), np.asarray(held_prev, float)
        nf = len(self.A_fit)
        self.pos = np.full(nf, -1, dtype=np.int64)
        self.pos[LK["late_idx"]] = np.arange(len(LK["late_idx"]))
        self.label = label
        self.calls = {"forward": 0, "seen": 0, "sensitivities": 0, "override_forward": 0}
        self.seconds = {k: 0.0 for k in self.calls}

    def _cols(self, facets):
        p = self.pos[np.asarray(facets, np.int64)]
        if np.any(p < 0):
            raise SystemExit("CpuLateEngine: a facet outside the late kernels")
        return p

    def forward(self):
        self.calls["forward"] += 1
        return self.p0_ng.copy(), self.p0_rg.copy()

    def seen(self, mask, structural=False):
        self.calls["seen"] += 1
        return self.LK["engine"].seen(mask, structural)

    def sensitivities(self, facets):
        self.calls["sensitivities"] += 1
        p = self._cols(facets)
        return np.ascontiguousarray(self.LK["G_NG"][:, p]), np.ascontiguousarray(self.LK["G_RG"][:, p])

    def forward_override(self, facets, g):
        self.calls["override_forward"] += 1
        f = np.asarray(facets, np.int64)
        p = self._cols(f)
        dA = self.A_fit[f] * (np.asarray(g, float) / np.exp(self.held_prev[f])) - self.A_fit[f]
        return self.p0_ng + self.LK["G_NG"][:, p] @ dA, self.p0_rg + self.LK["G_RG"][:, p] @ dA

    def identity(self):
        return {"engine": "localnn_ft.CpuLateEngine (CPU; late kernels of localnn_late.py)", "state": self.label,
                "late_kernels": self.LK.get("path"), "late_kernels_sha256": self.LK.get("sha256")}

    def close(self):
        pass


def load_late(path, nf):
    LK = llate.load(path, nf)
    LK["path"], LK["sha256"] = str(Path(path).resolve()), lp.sha256_file(path)
    return LK


# ------------------------------------------------------------------------------------------------ verification
def cmd_verify_softplus(args):
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    ctx = nft.NonnegFitContext(args.fit)
    LK = load_late(args.late, ctx.nf)
    p0_ng, p0_rg = LK["states"][args.state_name]
    if args.perturb:            # sensitivity test: p_fit x (1 + eps r), r = +-1 (seeded), the rounding scale of p_fit
        rng = np.random.default_rng(args.perturb_seed)
        p0_ng = p0_ng * (1.0 + args.perturb * rng.choice((-1.0, 1.0), len(p0_ng)))
        p0_rg = p0_rg * (1.0 + args.perturb * rng.choice((-1.0, 1.0), len(p0_rg)))
    eng = CpuLateEngine(LK, p0_ng, p0_rg, ctx.A_fit, ctx.held_prev, "softplus %s" % args.fit)
    gft = Path(args.gpu_ft)
    gj = json.load(open(gft / "fine_tuning.json"))
    if Path(gj["fit"]["directory"]).resolve() != Path(args.fit).resolve():
        raise SystemExit("the GPU fine tuning %s was made from %s" % (gft, gj["fit"]["directory"]))
    s = gj["settings"]
    argv = ["--fit", str(Path(args.fit).resolve()), "--out", str(out / "ft"), "--facets", s["facets"], "--window-min",
            str(s["window_min"]), "--max-nfev", str(s["max_nfev"]), "--lsmr-maxiter", str(s["lsmr_maxiter"]), "--ftol",
            repr(s["ftol"]), "--xtol", repr(s["xtol"]), "--gtol", repr(s["gtol"]), "--x-scale-grain",
            repr(s["x_scale_grain"]), "--seed", str(s["seed"]), "--save-seconds", repr(s["save_seconds"])]
    if s.get("validate"):
        argv.append("--validate")
    out.mkdir(parents=True)
    rc = nft.run(nft.parse_args(argv), engine_factory=lambda c: eng, command=[__file__, "verify-softplus"] + argv)
    with np.load(out / "ft" / "finetune_state.npz", allow_pickle=False) as z:
        c_delta, c_seen, c_late = np.asarray(z["delta"]), np.asarray(z["seen"]), np.asarray(z["late"])
        c_A, c_B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        cm = json.loads(str(z["meta_json"]))
    with np.load(gft / "finetune_state.npz", allow_pickle=False) as z:
        g_delta, g_seen, g_late = np.asarray(z["delta"]), np.asarray(z["seen"]), np.asarray(z["late"])
        g_A, g_B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        gm = json.loads(str(z["meta_json"]))
    rel = lambda a, b: abs(a / b - 1.0) if b != 0 else abs(a - b)     # noqa: E731
    res = {"created_utc": lp.utc_now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "fit": str(Path(args.fit).resolve()), "gpu_ft": str(gft.resolve()), "rc": rc,
           "perturb": args.perturb, "perturb_seed": args.perturb_seed,
           "late_equal": bool(np.array_equal(c_late, g_late)), "late_facets": int(c_late.sum()),
           "parameters_equal": bool(np.array_equal(c_seen, g_seen)), "parameters": int(c_seen.sum()),
           "delta_max_abs_diff": float(np.max(np.abs(c_delta - g_delta))),
           "delta_max_abs": float(np.max(np.abs(g_delta))),
           "grain_A": [c_A, g_A, rel(c_A, g_A)], "grain_B": [c_B, g_B],
           "costs_initial": {k: [cm["costs_initial"][k], gm["costs_initial"][k], rel(cm["costs_initial"][k],
                                                                                    gm["costs_initial"][k])]
                             for k in ("cost_total", "cost_data", "cost_crossing") if gm["costs_initial"][k] != 0},
           "costs_final": {k: [cm["costs_final"][k], gm["costs_final"][k], rel(cm["costs_final"][k], gm["costs_final"][k])]
                           for k in ("cost_total", "cost_data", "cost_crossing") if gm["costs_final"][k] != 0},
           "nfev": [cm["nfev_total"], gm["nfev_total"]], "stop": [cm["stop"], gm["stop"]],
           "validation": json.load(open(out / "ft" / "validation.json")).get("passed")
           if (out / "ft" / "validation.json").exists() else None,
           "wall_seconds": time.time() - t0}
    res["passed"] = bool(rc == 0 and res["late_equal"] and res["parameters_equal"] and
                         res["costs_initial"]["cost_total"][2] < 1e-12 and res["costs_final"]["cost_total"][2] < 1e-9)
    with open(out / "verify.json", "x") as f:
        json.dump(res, f, indent=1, default=float)
    say("verify %s: %s" % ("PASSED" if res["passed"] else "FAILED", json.dumps(res, default=float)))
    return 0 if res["passed"] else 1


# ------------------------------------------------------------------------------------------------ localnn fits
_CFG = {}


class LocalnnFitContext:
    """nonneg_finetune.NonnegFitContext's attributes for a localnn fit (see the module docstring)."""

    def __init__(self, fit_dir, pool=None):
        import localnn_mesh as lm
        cfg = _CFG
        self.fit_dir = Path(fit_dir).resolve()
        f = {k: self.fit_dir / k for k in ("localnn_state.npz", "fit_manifest.json", "nonneg_facets.npz",
                                           "COPS_model.dat", "late_window_model.dat", "fit_result.txt")}
        missing = [str(p) for p in f.values() if not p.is_file()]
        if missing:
            raise SystemExit("the localnn fit directory lacks %s" % missing)
        self.files = f
        self.identity = {"fit_restart.npz": nft.file_identity(f["localnn_state.npz"]),
                         "fit_manifest.json": nft.file_identity(f["fit_manifest.json"]),
                         "COPS_model.dat": nft.file_identity(f["COPS_model.dat"]),
                         "fit_result.txt": nft.file_identity(f["fit_result.txt"])}
        with np.load(f["localnn_state.npz"], allow_pickle=False) as z:
            self.meta = json.loads(str(z["meta_json"]))
            self.x = np.asarray(z["x"], float)
            self.names = [str(v) for v in z["names"]]
            self.c, self.m = np.asarray(z["c"], float), np.asarray(z["m"], float)
            self.kappa = float(z["kappa"])
            self.xA, self.xB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        if self.meta.get("stop") != "optimizer":
            raise SystemExit("the fit %s did not end by an optimizer condition (%s)" % (self.fit_dir, self.meta.get("stop")))
        self.fit_parameterization = self.meta["parameterization"]
        costs = self.meta["costs"]
        self.fit_costs = {"cost": float(costs["data_cost"]), "data_cost": float(costs["data_cost"]),
                          "regularization_cost": 0.0, "crossing_penalty_cost": float(costs["crossing_cost"]),
                          "lambda": 0.0, "smoothness_penalty_of_the_fit": float(costs["penalty"]),
                          "mult_penalty_of_the_fit": float(costs.get("mult_penalty", 0.0))}
        self.fit_result_txt = {"final_cost": float(costs["data_cost"]), "final_data_cost": float(costs["data_cost"])}
        self.model_cfg = {"model": "localnn", "parameterization": self.fit_parameterization,
                          "templates": self.meta["templates"], "multipliers": [float(v) for v in self.m],
                          "lambda_C": self.meta["lambda_C"],
                          "mult_penalty_f": float((self.meta.get("mult_penalty") or {}).get("f_m", 0.0))}
        self.held_record = self.meta.get("held")
        self.run = Path(self.meta["base"]).resolve()
        self.point_table = Path(self.meta["data"]["point_table"]["path"]).resolve()
        self.escale = float(self.meta["A_scale"])
        self.T = np.array(lp.USER_T0)
        with np.load(f["nonneg_facets.npz"], allow_pickle=False) as z:
            self.A_fit = np.asarray(z["A"], float)
            self.u_fit, self.TS_fit = np.asarray(z["u"], float), np.asarray(z["TS"], float)
            self.bins = np.asarray(z["bin"])
            self.seen_fit = np.asarray(z["seen"], bool)
            self.held_prev = np.asarray(z["held"], float)
        self.nf = len(self.A_fit)
        el = lm.read_elements(self.run / "elements.bin")
        self.area, self.A_orig = el["area"], el["AH2O"]
        # the activity reproduced bitwise: A = m[bin] max(0, u + TS) exp(h)
        A_chk = self.m[self.bins.astype(np.intp)] * np.maximum(self.u_fit + self.TS_fit, 0.0)
        hz = self.held_prev != 0.0
        A_chk[hz] = A_chk[hz] * np.exp(self.held_prev[hz])
        if not np.array_equal(A_chk, self.A_fit):
            raise SystemExit("the fit's activity cannot be reproduced from its facets npz (u, TS, m, held)")
        if (self.held_record is None) != (not hz.any()):
            raise SystemExit("held term and its record disagree")
        self.active_fit = self.A_fit > 0.0
        R = llate.build_rows(self.run, cfg["d2"], self.point_table, NONNEG_CODE)
        self.rows, self.layout, self.utc = R["rows"], R["layout"], R["utc"]
        self.header, self.table_rows = R["header"], R["table_rows"]
        self.s_start = np.log(np.where(self.A_fit > 0.0, self.A_fit, 1e-300))
        # the fit state's instantaneous surface pressures of every row (design route, the held term exact)
        design = lp.load_design(cfg["design"])
        sN, sR = lp.all_rows_surface(design, self.c, self.m, self.kappa, self.meta["templates"] == "on")
        LK = cfg["LK"]
        li = LK["late_idx"]
        A_noheld = self.m[self.bins[li].astype(np.intp)] * np.maximum(self.u_fit[li] + self.TS_fit[li], 0.0)
        dA = self.A_fit[li] - A_noheld
        if np.any(dA[self.held_prev[li] == 0.0] != 0.0) or np.any(self.held_prev[~LK["late"]] != 0.0):
            raise SystemExit("held facets outside the late kernels")
        self.p0_ng, self.p0_rg = sN + LK["G_NG"] @ dA, sR + LK["G_RG"] @ dA


def cmd_run(args, rest):
    import localnn_mesh as lm
    t0 = time.time()
    el_n = lm.read_elements(Path(args.base) / "elements.bin")["nfaces"]
    LK = load_late(args.late, el_n)
    _CFG.update({"d2": args.d2, "design": args.design, "LK": LK})
    holder = {}

    def ctx_factory(fit_dir, pool=None):
        c = LocalnnFitContext(fit_dir, pool)
        holder["ctx"] = c
        return c
    nft.NonnegFitContext = ctx_factory            # the context of THIS process only (the module file is unchanged)

    def engine_factory(c):
        return CpuLateEngine(LK, c.p0_ng, c.p0_rg, c.A_fit, c.held_prev, "localnn %s" % c.fit_dir)
    argv = ["--fit", args.fit, "--out", args.out] + list(rest)
    rc = nft.run(nft.parse_args(argv), engine_factory=engine_factory, command=[__file__, "run"] + argv)
    out = Path(args.out)
    if (out / "fine_tuning.json").exists():
        note = {"localnn_ft": {"engine": "CpuLateEngine", "late_kernels": {"path": LK["path"], "sha256": LK["sha256"]},
                               "nonneg_code": str(NONNEG_CODE), "code_sha256": {
                                   p: lp.sha256_file(HERE / p) for p in ("localnn_ft.py", "localnn_late.py",
                                                                         "localnn_problem.py")},
                               "seconds": time.time() - t0}}
        with open(out / "localnn_ft.json", "x") as f:
            json.dump(note, f, indent=1)
    return rc


def cmd_run_nonneg(args, rest):
    t0 = time.time()
    holder = {}

    def engine_factory(c):
        LK = load_late(args.late, c.nf)
        if args.state_name not in LK["states"]:
            raise SystemExit("late kernels %s have no state %s" % (args.late, args.state_name))
        p0_ng, p0_rg = LK["states"][args.state_name]
        holder["LK"] = LK
        return CpuLateEngine(LK, p0_ng, p0_rg, c.A_fit, c.held_prev, "nonneg %s (state %s)" % (c.fit_dir, args.state_name))
    argv = ["--fit", args.fit, "--out", args.out] + list(rest)
    rc = nft.run(nft.parse_args(argv), engine_factory=engine_factory, command=[__file__, "run-nonneg"] + argv)
    out = Path(args.out)
    if (out / "fine_tuning.json").exists() and "LK" in holder:
        LK = holder["LK"]
        note = {"localnn_ft": {"engine": "CpuLateEngine", "late_kernels": {"path": LK["path"], "sha256": LK["sha256"]},
                               "state": args.state_name, "context": "nonneg_finetune.NonnegFitContext (unchanged)",
                               "nonneg_code": str(NONNEG_CODE), "code_sha256": {
                                   p: lp.sha256_file(HERE / p) for p in ("localnn_ft.py", "localnn_late.py",
                                                                         "localnn_problem.py")},
                               "seconds": time.time() - t0}}
        with open(out / "localnn_ft.json", "x") as f:
            json.dump(note, f, indent=1)
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify-softplus")
    for k in ("fit", "gpu-ft", "late", "state-name", "out"):
        v.add_argument("--" + k, required=True)
    v.add_argument("--perturb", type=float, default=0.0, help="sensitivity test: relative perturbation of p_fit")
    v.add_argument("--perturb-seed", type=int, default=1)
    r = sub.add_parser("run")
    for k in ("fit", "late", "design", "base", "d2", "out"):
        r.add_argument("--" + k, required=True)
    rn = sub.add_parser("run-nonneg")
    for k in ("fit", "late", "state-name", "out"):
        rn.add_argument("--" + k, required=True)
    args, rest = ap.parse_known_args(argv)
    if args.cmd == "run-nonneg":
        return cmd_run_nonneg(args, rest)
    if args.cmd == "verify-softplus":
        if rest:
            raise SystemExit("unknown arguments %s" % rest)
        return cmd_verify_softplus(args)
    return cmd_run(args, rest)


if __name__ == "__main__":
    sys.exit(main())
