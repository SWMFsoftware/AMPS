#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_finetune.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_finetune.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: cd5342b391a1cf2c23c24dc3496467e8e95acb9df9d6867fc2fd9b1a532b5050
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Late-facet fine tuning of a slope-multiplier / bump fit (user decisions 2026-11-10 "add fine tuning --
#   fitting cycles as in other fits"): delta_e added to the held term, A_e = max(F, m_b L exp(h + delta)), delta >=
#   ln(F/s), grain refitted, the fit's full objective, positive_finetune's optimizer; exact host model (pressures linear
#   in A, G = dp/dA once on the GPU); validation checks 1-8, v8: checks 6 and 7 with pressure-sized steps (user
#   2026-11-10 "Fix the check step").
# USED BY: jobs/slope_finetune.sbatch (v8: GPU test 124063, 3a chain fine tunings 124064, 124073, 124082).
# ==== end of header ====
"""Late-facet fine tuning of a slope-multiplier / bump + slope-multiplier fit (user decisions 2026-11-10).

    .venv/bin/python -W ignore slope_finetune.py --fit <fit seg dir> --out <new dir> --facets late_only --window-min 5
        --devices 0,...,7 --max-nfev 500 [--l-threads 32] [--validate | --validate-only] [--wall-time slurm]
        [--resume-from <finetune_state.npz of a stopped segment>] [--segment k --nfev-before N]

Model (user 2026-11-10; the analog of the positive fits' per-element fine tuning, jobs/positive_cycles.sh):
* the tuned fit's activity A_fit(e) = max(F(e), s(e)), s = m_b(e) L(e) exp(h(e)) (h = the held term of the fit, 0 if
  none; slope_mult.py); the fine tuning gives every PARAMETER facet e an increment delta_e of the held term:
      A_e(delta) = s_e exp(delta_e)   (h -> h + delta),   delta_e >= ln(F_e / s_e)  (A >= F = 1e-5 A_orig; F = 0: none)
  i.e. dA/ddelta = A above the floor (the positive tool's exp form, applied to this model's held factor); the next
  cycle's fit holds h_next = h + delta in A = max(F, m_b L exp(h_next)), which equals the tuned A wherever the bound
  holds (A >= F).
* parameters: FACETS=late_only (user 2026-11-07/08: the facets seen exclusively in the last WINDOW_MIN minutes: a
  structural record (Kn, KPhi or KV != 0) at an objective point inside [t_end - WINDOW_MIN, t_end] and none at any
  objective point before; positive_finetune.late_only_facets, unchanged, on this model's GPU records), MINUS the late
  facets AT THE FLOOR in the tuned fit (s_e <= F_e): there dA/ddelta = 0 in this model (the fit's one-sided clip
  convention), so they are not tuned (counted and recorded). Grain A, B >= 0 refitted. Everything else held.
* objective: the fit's full objective (data rows of its selection, log residuals, W10x5 + grain for the late
  observations, the crossing penalty with the fit's weight; lambda 0 as the fit): positive_finetune.Rows /
  Objective / Jacobian / optimize (TRF + LSMR, x_scale 1, LSMR_MAXITER 50, MAXNFEV 500), unchanged.
Method (why): with the temperatures and kappa held, the surface pressures are LINEAR in A, so for the parameter facets
      p(delta) = p_fit + G (A(delta) - A_fit),   G = dp/dA (instantaneous surface model, every row; NG and RG),
exactly. p_fit is the GPU model at the fit's state (bitwise the fit's model) and G is computed once on the GPU with
the evaluator's own column code (slope_mult sensitivities: active = 1 at the facet, colbasis = its indicator); every
evaluation, J v and J^T u is then a small host product (rows x parameter facets): no GPU work per evaluation. The GPU
checks the host model against direct evaluations (validation: per-facet factor override in the GPU patch).
Threads: --l-threads N: L of the fit's state (the GPU forward) is the chunked bitwise L of the fits; the 12.4M-facet
host fields (ln A, the products) are formed in chunks on the same pool (elementwise, bitwise); OpenBLAS stays 1 thread.
Products in --out: finetune_state.npz (delta, held_prev, held_next = held_prev + delta, A_final, lower_bound, seen =
the parameters, late, grain_A_uPa, grain_B_uPa_per_100s, meta_json), finetune_result.txt, COPS_model.dat (every row
recomputed on the GPU at the final delta: late observations W10x5 + grain, grid rows instantaneous + grain),
late_window_model.dat, pressure_comparison.txt, fine_tuning.json (settings, identities, timings), validation.json.
Exit: 0 finished, 75 wall-time stop (products of the best state; the job resubmits), 1 error / failed validation.
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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP / "tools"), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import positive_finetune as pf      # Rows, Objective, Jacobian, Tracker, optimize, late_only_facets (used unchanged)
import slope_mult as sm

PARAMETERIZATION_ID = "slope_mult_late_finetune_v1"
FIT_PREFIXES = ("slope_multiplier_floor_w10x5", "bump_slope_multiplier_floor_w10x5")
POLICY = ("user 2026-11-10: fine tuning as in the positive cycles (FACETS late_only, WINDOW_MIN 5, the fit's full "
          "objective, grain refitted, lambda 0, MAXNFEV 500, no cap); for this model A_e -> max(F_e, m_b L_e exp(h_e + "
          "delta_e)) with delta_e >= ln(F_e / s_e) (the floor as the lower bound, dA/ddelta = A above the floor), the "
          "late facets at the floor in the fit are not tuned (dA/ddelta = 0 there); the next fit holds h + delta")
GRAIN_NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
CHUNK = 1 << 20
# Validation checks 6 and 7: single-facet finite differences with steps sized by the facet's pressure change (user
# 2026-11-10 "Fix the check step"; thresholds unchanged, 1e-4 and 1e-6). dA_k = TAU p_r* / |G_r*k| at the facet's
# strongest row r* (largest |G_rk| A_k / p_r), i.e. the step changes that row's pressure by the fraction TAU. The former
# steps relative to A_k (delta 1e-3; A -> 2 A) put the signal of facets with tiny A_k or tiny geometric weight at the
# rounding level of p (jobs 124013, 124015; diagnostic 124034: G confirmed per facet, error ~ eps / TAU).
# TAU6 (host, delta step h = ln(1 + dA / A)): the log residuals are curved in p, error ~ TAU6 / 2 (124034 part E:
# 1e-7 .. 3e-7 at 1e-6), rounding ~ eps / TAU6. TAU7 (GPU): the pressures are exactly linear in A, no truncation,
# error ~ eps / TAU7 (124034 part C: max 7.3e-7 at 1e-6 over all 447 facets of both fits, <= 1.5e-10 at 1e-4).
CHECK6_TAU = 1e-6
CHECK7_TAU = 1e-4


def log(message):
    print("slope_finetune: %s" % message, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def file_identity(path):
    p = Path(path)
    return {"path": str(p.resolve()), "sha256": sha256(p), "size_bytes": p.stat().st_size}


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


jsonable = pf.jsonable


def elementwise(func, n, out, pool):
    """out[i0:i1] = func(i0, i1) in chunks (on `pool` when given): elementwise numpy, bitwise any thread count."""
    starts = range(0, n, CHUNK)

    def part(i0):
        i1 = min(i0 + CHUNK, n)
        out[i0:i1] = func(i0, i1)
    if pool is None:
        for i0 in starts:
            part(i0)
    else:
        list(pool.map(part, starts))
    return out


# ------------------------------------------------------------------------------------------------ the tuned fit
class SlopeFitContext:
    """Everything the fine tuning takes from the tuned fit (CPU only)."""

    def __init__(self, fit_dir, pool=None):
        from cops_cuda import reference
        import grain_crossing as gc
        import positive_problem as pp
        reference.bind()
        fit_cops = reference.load("fit_cops")
        templates = reference.load("activity_templates")
        self.fit_dir = Path(fit_dir).resolve()
        f = {k: self.fit_dir / k for k in ("fit_restart.npz", "fit_manifest.json", "slope_mult_facets.npz",
                                           "grain_term.json", "COPS_model.dat", "late_window_model.dat",
                                           "fit_result.txt")}
        missing = [str(p) for p in f.values() if not p.is_file()]
        if missing:
            raise SystemExit("the fit directory lacks %s" % missing)
        self.files = f
        self.identity = {k: file_identity(p) for k, p in f.items() if k not in ("late_window_model.dat",
                                                                                "slope_mult_facets.npz")}
        with np.load(f["fit_restart.npz"], allow_pickle=False) as z:
            self.restart = {k: z[k] for k in z.files if k not in ("meta_json", "activity_linear", "template_columns",
                                                                  "centroid")}
            self.meta = json.loads(str(z["meta_json"]))
        ext = self.meta.get("model_extension") or {}
        self.fit_parameterization = str(ext.get("parameterization", ""))
        if not self.fit_parameterization.startswith(FIT_PREFIXES):
            raise SystemExit("%s is not a slope_mult_fit.py fit (marker %r)" % (self.fit_dir, self.fit_parameterization))
        if self.meta.get("status") != "complete":
            raise SystemExit("the fit's restart has status %r (a finished fit is required)" % self.meta.get("status"))
        self.model_cfg = dict(ext["model"])
        self.ext = ext
        names = list(self.meta["names"])
        if tuple(names[-2:]) != GRAIN_NAMES:
            raise SystemExit("the fit's last parameters are %s, not the grain" % names[-2:])
        self.names, self.x = names, np.asarray(self.restart["x"], float)
        self.xA, self.xB = float(self.x[-2]), float(self.x[-1])
        self.a = np.asarray(self.restart["a"], float)
        self.T = np.asarray(self.restart["T"], float)
        self.kappa = float(np.asarray(self.restart["kappa"]).reshape(-1)[0])
        self.b = np.asarray(self.restart["template_values"], float)
        self.lam = float(self.meta["settings"]["lambda"])
        if self.lam != 0.0:
            raise SystemExit("lambda %g: this fine tuning implements the lambda-0 objective of the slope fits" % self.lam)
        self.escale = float(self.meta["enhancement_scale"])
        # manifest
        self.manifest = json.loads(f["fit_manifest.json"].read_text())
        mext = self.manifest["model_extension"]
        self.crossing_weight = float(mext["crossing_weight"])
        rec = mext["slope_multiplier"]
        self.fit_record_config = rec["config"]
        res = self.manifest["optimizer_results"]
        self.fit_costs = {"cost": float(res[0]["cost"]), "data_cost": float(res[0]["data_cost"]),
                          "regularization_cost": float(res[0]["regularization_cost"]),
                          "crossing_penalty_cost": float(mext.get("final", {}).get("penalty_cost", float("nan"))),
                          "lambda": float(res[0]["lambda"])}
        self.fit_result_txt = {}
        for line in f["fit_result.txt"].read_text().splitlines():
            v = line.split()
            if len(v) >= 2 and v[0] in ("final_cost", "final_data_cost", "final_regularization_cost", "lambda"):
                self.fit_result_txt[v[0]] = float(v[1])
        gt = json.loads(f["grain_term.json"].read_text())
        self.grain_term = {"t_end_et": float(gt["t_end_et"]), "last_ng_before_gap_et": float(gt["last_ng_before_gap_et"])}
        # base, elements
        self.run = Path(self.fit_record_config["run"]).resolve()
        self.point_table = Path(self.fit_record_config["point_table"]).resolve()
        header, names_t, table_rows = fit_cops.read_model_file(str(self.run / "COPS_model.dat"))
        self.header, self.table_rows = header, table_rows
        self.cols = {n: i for i, n in enumerate(names_t)}
        self.el = fit_cops.read_elements(str(self.run / "elements.bin"))
        self.area = np.asarray(self.el["area"], float)
        self.A_orig = np.asarray(self.el["AH2O"], float)
        self.nf = len(self.area)
        self.coef = fit_cops.read_activity_coefficients(str(Path(fit_cops.__file__).resolve().parents[1]), header)
        # template basis (bins, bumps) from the fit's bundles, in the fit's order
        sb = ext["slope_bundle"]
        bundles = [(sb["path"], sb["sha256"], list(sm.BIN_NAMES))]
        if self.model_cfg.get("basis") == "bumps":
            bb = ext["bump_bundle"]
            bundles.append((bb["path"], bb["sha256"], list(self.model_cfg["bump_names"])))
        cols = []
        for path, digest, tnames in bundles:
            if sha256(path) != digest:
                raise SystemExit("bundle %s changed since the fit" % path)
            bundle = templates.load_bundle(path, self.el["centroid"], self.area)
            cols.append(templates.select_templates(bundle, tnames))
        self.basis = np.ascontiguousarray(np.column_stack(cols), dtype=np.float64)
        del cols
        if self.basis.shape[1] != len(self.b):
            raise SystemExit("%d template columns, %d template values" % (self.basis.shape[1], len(self.b)))
        # the fit's per-facet field
        with np.load(f["slope_mult_facets.npz"], allow_pickle=False) as z:
            self.A_fit = np.asarray(z["A"], float)
            self.L_fit = np.asarray(z["L"], float)
            self.F = np.asarray(z["floor"], float)
            self.bins = np.asarray(z["bin"])
            self.seen_fit = np.asarray(z["seen"], bool)
            self.held_prev = np.asarray(z["held"], float) if "held" in z.files else np.zeros(self.nf)
        hs = (self.model_cfg.get("held") or {}).get("sha256")
        if hs is not None and sm.sha256_bytes(self.held_prev) != hs:
            raise SystemExit("the facets npz held term differs from the fit's configured held term")
        if hs is None and np.any(self.held_prev != 0.0):
            raise SystemExit("held term in the facets npz but none in the fit's configuration")
        nb = len(sm.BIN_NAMES)
        g = np.exp(self.held_prev) if hs is not None else None
        s, A, act = sm.activity(self.L_fit, self.b[:nb], self.bins, self.F, g)
        if not np.array_equal(A, self.A_fit):
            raise SystemExit("the fit's activity cannot be reproduced from its facets npz (L, m, floor, held)")
        self.s_fit, self.active_fit = s, act
        # rows, selection, window, grain, crossing (positive_finetune.FitContext's construction)
        c = self.cols
        ipoint = np.array([int(r[c["iPoint"]]) for r in table_rows])
        if not np.array_equal(ipoint, np.asarray(self.restart["point_ids"])):
            raise SystemExit("the restart's point ids differ from the base table's rows")
        refNG = np.array([float(r[c["refNG"]]) for r in table_rows])
        refRG = np.array([float(r[c["refRG"]]) for r in table_rows])
        self.utc = [r[c["UTC"]] for r in table_rows]
        ets = gc.et_by_point(self.run)
        et = np.array([ets[int(i)] for i in ipoint])
        table = pp.read_point_table(str(self.point_table))
        self.layout = pp.LateLayout(table, pp.read_run_points(self.run), ipoint)
        window = pp.WindowOperator(self.layout, np.arange(len(ipoint)))
        self.rows = pf.Rows(ipoint, et, refNG, refRG, self.restart["use_ng"], self.restart["use_rg"], window,
                            self.crossing_weight, gc.GAP_SEARCH_S, (gc.A_UNIT, gc.B_UNIT))
        if (self.rows.et_end, self.rows.et_before) != (self.grain_term["t_end_et"],
                                                         self.grain_term["last_ng_before_gap_et"]):
            raise SystemExit("grain onset differs from the fit's grain_term.json")
        # exponent form for positive_finetune.Objective: s_start = ln A (A = 0: ln 1e-300, a diagnostic proxy only)
        self.s_start = np.empty(self.nf)
        A_ = self.A_fit
        elementwise(lambda i0, i1: np.log(np.where(A_[i0:i1] > 0.0, A_[i0:i1], 1e-300)), self.nf, self.s_start, pool)


# ------------------------------------------------------------------------------------------------ engines
class HostLinearEngine:
    """positive_finetune.Objective's engine interface on the exact linear model p = p_fit + G (A - A_fit)."""

    def __init__(self, p0_ng, p0_rg, G_ng, G_rg, A0, s_start, idx):
        self.p0_ng, self.p0_rg = np.asarray(p0_ng, float), np.asarray(p0_rg, float)
        self.G_ng, self.G_rg = np.ascontiguousarray(G_ng, float), np.ascontiguousarray(G_rg, float)
        self.A0 = np.asarray(A0, float)
        self.s0_idx = np.asarray(s_start, float)[idx]
        self.idx = np.asarray(idx, np.int64)
        self.bufs = []
        self.reuploads = 0
        self.seconds = {"forward": 0.0, "jv": 0.0, "jtu": 0.0}

    def set_parameters(self, idx):
        if not np.array_equal(np.asarray(idx), self.idx):
            raise SystemExit("engine parameters differ from the objective's")
        return {"parameters": int(len(self.idx)), "engine": "host linear model p_fit + G (A - A_fit)"}

    def A(self, s):
        return self.A0 * np.exp(np.asarray(s, float)[self.idx] - self.s0_idx)

    def forward(self, s):
        t0 = time.perf_counter()
        dA = self.A(s) - self.A0
        out = (None, self.p0_ng + self.G_ng @ dA, self.p0_rg + self.G_rg @ dA)
        self.seconds["forward"] += time.perf_counter() - t0
        return out

    def direction_buffer(self, k):
        while len(self.bufs) < k:
            self.bufs.append(np.zeros(len(self.idx)))
        return self.bufs[:k]

    def jv(self, st, k):
        t0 = time.perf_counter()
        A = self.A(st.s)
        JN, JR = [], []
        for j in range(k):
            w = A * self.bufs[j]
            JN.append(self.G_ng @ w)
            JR.append(self.G_rg @ w)
        self.seconds["jv"] += time.perf_counter() - t0
        return JN, JR

    def jtu(self, st, uN, uR):
        t0 = time.perf_counter()
        g = self.A(st.s) * (self.G_ng.T @ uN + self.G_rg.T @ uR)
        self.seconds["jtu"] += time.perf_counter() - t0
        return g


class GpuModelEngine:
    """The fit's model on the GPU (slope_mult patch): p at the fit's state, structural seen masks, dp/dA, forwards
    with a per-facet factor override. `backend` (optional): an already bound backend of the fit's rows and basis (the
    one-job chain driver); default: bind one here, exactly as the fit's MultiGpuBackend (host activity)."""

    def __init__(self, ctx, devices, backend=None, l_threads=None, l_mode="chunked"):
        from cops_cuda import backend as bk
        from cops_cuda import gpu_model, precision, reference
        reference.bind()
        fit_cops = reference.load("fit_cops")
        self.ctx = ctx
        self.timing = {}
        t0 = time.time()
        if not sm._STATE["installed"]:
            sm.install_gpu(ctx.model_cfg)
        elif sm.canonical(sm.installed_config()) != sm.canonical(ctx.model_cfg):
            raise SystemExit("the GPU patch is installed with another configuration than the fit's")
        if l_threads:
            sm.set_L_threads(int(l_threads), l_mode)
        self.state = bk.ParameterState(ctx.a, ctx.T, ctx.kappa, ctx.b)
        if backend is None:
            os.environ[precision.RECORD_PRECISION_ENVIRONMENT] = "fp32"           # the fits' --precision mixed
            gpu_model.REGISTRY.configure(devices[0], "tree", 4096, devices=devices, transport="local",
                                         activity="host")
            Cls = gpu_model.gpu_model_class()
            M = Cls(str(ctx.run), ctx.header, ctx.cols, ctx.table_rows, ctx.el, ctx.coef)
            if not np.array_equal(np.asarray(M.rows), np.arange(len(ctx.table_rows))):
                raise SystemExit("the model rows are not every row of the base")
            self.M = M
            backend = gpu_model.REGISTRY.backend_for(M, ctx.basis, ctx.escale)
            self.own_backend = True
        else:
            self.M = None
            self.own_backend = False
        if bk.array_digest(ctx.basis) != backend.enhancement_basis_sha256 or float(backend.enhancement_scale) != \
                ctx.escale:
            raise SystemExit("the backend's template basis / scale is not the fit's")
        fit_digest = self._fit_basis_digest()
        if fit_digest is not None and fit_digest != backend.enhancement_basis_sha256:
            raise SystemExit("the fit's recorded basis digest %s differs from this basis %s" % (
                fit_digest, backend.enhancement_basis_sha256))
        self.backend = backend
        self.timing["bind_s"] = time.time() - t0
        self.calls = {"forward": 0, "seen": 0, "sensitivities": 0, "override_forward": 0}
        self.seconds = {"forward": 0.0, "seen": 0.0, "sensitivities": 0.0, "override_forward": 0.0}

    def _fit_basis_digest(self):
        try:
            for be in self.ctx.manifest["compute_backend"]["registry"]["backends"]:
                if be.get("enhancement_basis_sha256"):
                    return be["enhancement_basis_sha256"]
        except Exception:
            return None
        return None

    def forward(self):
        from cops_cuda import backend as bk
        t0 = time.perf_counter()
        r = self.backend.evaluate(self.state, bk.PRESSURES_ONLY)
        self.calls["forward"] += 1
        self.seconds["forward"] += time.perf_counter() - t0
        return np.array(r.pNG), np.array(r.pRG)

    def seen(self, mask, structural=False):
        t0 = time.perf_counter()
        out = self.backend.slopemul_seen_rows(np.asarray(mask, bool), structural)
        self.calls["seen"] += 1
        self.seconds["seen"] += time.perf_counter() - t0
        return np.asarray(out, bool)

    def sensitivities(self, facets):
        t0 = time.perf_counter()
        ng, rg = self.backend.slopemul_sensitivities(self.state, np.asarray(facets, np.int64))
        self.calls["sensitivities"] += 1
        self.seconds["sensitivities"] += time.perf_counter() - t0
        return ng, rg

    def forward_override(self, facets, g_values):
        t0 = time.perf_counter()
        self.backend.slopemul_set_override(np.asarray(facets, np.int64), np.asarray(g_values, np.float64))
        try:
            out = self.forward()
        finally:
            self.backend.slopemul_set_override(None, None)
        self.calls["override_forward"] += 1
        self.seconds["override_forward"] += time.perf_counter() - t0
        return out

    def identity(self):
        try:
            return jsonable(self.backend.identity())
        except Exception as exc:
            return {"error": "%s: %s" % (exc.__class__.__name__, exc)}

    def close(self):
        if self.own_backend:
            from cops_cuda import gpu_model
            gpu_model.REGISTRY.close()


class _SeenAdapter:
    """positive_finetune.late_only_facets's engine.seen(mask, structural)."""

    def __init__(self, engine):
        self.engine = engine

    def seen(self, mask, structural=False):
        return self.engine.seen(mask, structural)


class _CtxAdapter:
    def __init__(self, ctx):
        self.utc, self.area, self.A_orig = ctx.utc, ctx.area, ctx.A_orig


# ------------------------------------------------------------------------------------------------ validation
def validate(obj, ctx, eng, z0, idx, held_g_idx, support, seed=20261110):
    """(record, passed): 1 start cost = the fit's final cost (1e-8); 2 adjoint; 3 J v vs central FD (host, random
    direction); 4 host model vs a direct GPU evaluation at a random delta (override); 5 strict support (rows before
    the window bitwise unchanged on the GPU, rows inside change); 6 single facets J vs FD (host); 7 G vs a GPU FD of
    single facets; 8 zero-sensitivity parameters."""
    rng = np.random.default_rng(seed)
    rec, ok = {}, True
    st0 = obj.evaluate(z0)
    # 1
    fit_total = ctx.fit_costs["data_cost"] + ctx.fit_costs["regularization_cost"] + ctx.fit_costs["crossing_penalty_cost"]
    r1 = {"total_start": st0.costs["total"], "data_start": st0.costs["data"], "crossing_start": st0.costs["crossing"],
          "fit_data_cost": ctx.fit_costs["data_cost"], "fit_crossing": ctx.fit_costs["crossing_penalty_cost"],
          "rel_total": pf.rel(st0.costs["total"], fit_total), "rel_data": pf.rel(st0.costs["data"],
                                                                                 ctx.fit_costs["data_cost"]),
          "data_bitwise": st0.costs["data"] == ctx.fit_costs["data_cost"]}
    r1["passed"] = r1["rel_total"] <= 1e-8 and r1["rel_data"] <= 1e-8
    rec["1_start_cost"] = r1
    ok &= r1["passed"]
    log("check 1 start cost: total %.15e (fit data + crossing %.15e, rel %.1e); data %.15e (fit %.15e, bitwise %s)" % (
        st0.costs["total"], fit_total, r1["rel_total"], st0.costs["data"], ctx.fit_costs["data_cost"],
        r1["data_bitwise"]))
    # 2 adjoint
    J = pf.Jacobian(obj, st0)
    v = rng.standard_normal(obj.n)
    u = rng.standard_normal(obj.m)
    lhs, rhs = float(u @ J.matvec(v)), float(v @ J.rmatvec(u))
    rec["2_adjoint"] = {"rel": pf.rel(lhs, rhs), "passed": pf.rel(lhs, rhs) <= 1e-12}
    ok &= rec["2_adjoint"]["passed"]
    # 3 J v vs central FD on the host objective (a random direction, step 1e-6)
    d = rng.standard_normal(obj.n)
    d[-2:] *= 0.0
    h = 1e-3          # one facet's share of a point's pressure is ~1e-8: a 1e-6 step would be below the rounding
    zp, zm = z0 + h * d, z0 - h * d
    zp = np.maximum(zp, obj.lb)
    zm = np.maximum(zm, obj.lb)
    fd = (obj.evaluate(zp).r - obj.evaluate(zm).r) / (2 * h)
    jv = J.matvec(d)
    nd = obj.rows.nd
    r3 = float(np.max(np.abs(fd[:nd] - jv[:nd])) / max(float(np.max(np.abs(fd[:nd]))), 1e-300))
    rec["3_jv_fd_random_direction"] = {"rel_data": r3, "step": h, "passed": r3 <= 1e-4}
    ok &= rec["3_jv_fd_random_direction"]["passed"]
    # grain columns
    for kk, nm in ((-2, "grain_A"), (-1, "grain_B")):
        e = np.zeros(obj.n)
        e[kk] = 1.0
        hh = 1e-6 * max(abs(z0[kk]), 1.0)
        zp, zm = z0 + hh * e, np.maximum(z0 - hh * e, obj.lb)
        fd = (obj.evaluate(zp).r - obj.evaluate(zm).r) / (zp[kk] - zm[kk])
        col = J.matvec(e)
        sc = max(float(np.max(np.abs(fd))), 1e-300)
        rec["3_" + nm] = {"rel": float(np.max(np.abs(fd - col)) / sc)}
    # 4 + 5: host model vs a direct GPU evaluation at a random delta, strict support
    dz = rng.uniform(-0.5, 0.5, len(idx))
    dz = np.maximum(dz, obj.lb_delta)
    z = z0.copy()
    z[:len(idx)] = dz
    st = obj.evaluate(z)
    g_new = held_g_idx * np.exp(dz) if held_g_idx is not None else np.exp(dz)
    # the GPU's activity for these facets: (m L) * g with g = exp(h + delta) (the next fit's form)
    hprev = ctx.held_prev[idx]
    g_next = np.exp(hprev + dz)
    pg_ng, pg_rg = eng.forward_override(idx, g_next)
    host_ng, host_rg = st.pNG_s, st.pRG_s
    rel4 = float(max(np.max(np.abs(host_ng / pg_ng - 1.0)), np.max(np.abs(host_rg / pg_rg - 1.0))))
    rec["4_host_vs_gpu_random_delta"] = {"max_rel": rel4, "passed": rel4 <= 1e-10,
                                        "note": "host p_fit + G (A - A_fit) vs the GPU model with g = exp(h + delta)"}
    ok &= rec["4_host_vs_gpu_random_delta"]["passed"]
    p0_ng, p0_rg = eng.forward()
    before, inside = support["before_rows"], support["inside_rows"]
    same_before = bool(np.array_equal(pg_ng[before], p0_ng[before]) and np.array_equal(pg_rg[before], p0_rg[before]))
    changed_in = int(np.sum((pg_ng[inside] != p0_ng[inside]) | (pg_rg[inside] != p0_rg[inside])))
    rec["5_strict_support"] = {"rows_before": int(before.sum()), "before_bitwise_unchanged": same_before,
                               "rows_inside": int(inside.sum()), "rows_inside_changed": changed_in,
                               "max_abs_ln_change_inside_NG": float(np.max(np.abs(np.log(pg_ng[inside] / p0_ng[inside])))),
                               "passed": same_before and changed_in > 0}
    ok &= rec["5_strict_support"]["passed"]
    log("check 4 host model vs GPU at a random delta: max rel %.2e; check 5 strict support: %d rows before the window "
        "bitwise unchanged %s, %d of %d rows inside changed" % (rel4, int(before.sum()), same_before, changed_in,
                                                                int(inside.sum())))
    # 6 single facets (host): J column vs a central FD in delta with a step sized by the facet's pressure change
    #   dA_k = CHECK6_TAU p_r* / |G_r*k| (strongest row over NG and RG), h_k = ln(1 + dA_k / A_k) (the equivalent step in
    #   delta), the minus side clipped at the lower bound (hm); for such steps the derivative of exp is replaced by its
    #   exact secant over the step: FD vs col (e^h - e^-hm) / (h + hm)
    eng0 = obj.engine
    A_now = eng0.A0 * np.exp(z0[:len(idx)])

    def strongest(k, kinds=("ng", "rg")):
        """(row pressure / |G| at the facet's strongest row, kind, row, share) or None (zero columns)."""
        best = None
        for kind in kinds:
            G = eng0.G_ng[:, k] if kind == "ng" else eng0.G_rg[:, k]
            p = eng0.p0_ng if kind == "ng" else eng0.p0_rg
            ok_ = (p != 0.0) & (G != 0.0)
            if not np.any(ok_):
                continue
            q = np.zeros(len(p))
            q[ok_] = np.abs(G[ok_]) * A_now[k] / np.abs(p[ok_])
            r_ = int(np.argmax(q))
            if best is None or q[r_] > best[3]:
                best = (abs(p[r_]) / abs(G[r_]), kind, r_, float(q[r_]))
        return best

    k6 = rng.choice(len(idx), size=min(8, len(idx)), replace=False)
    worst6, rec6 = 0.0, []
    for k in k6:
        e = np.zeros(obj.n)
        e[k] = 1.0
        sb = strongest(k)
        if sb is None:
            rec6.append({"facet": int(idx[k]), "skipped": "zero columns"})
            continue
        dA = CHECK6_TAU * sb[0]
        hh = math.log1p(dA / A_now[k])
        hm = min(hh, z0[k] - obj.lb[k])
        fd = (obj.evaluate(z0 + hh * e).r - obj.evaluate(np.maximum(z0 - hh * e, obj.lb)).r) / (hh + hm)
        col = J.matvec(e) * ((math.exp(hh) - math.exp(-hm)) / (hh + hm))
        sc = float(np.max(np.abs(fd[:nd])))
        err = float(np.max(np.abs(fd[:nd] - col[:nd])) / sc) if sc > 0 else 0.0
        worst6 = max(worst6, err)
        rec6.append({"facet": int(idx[k]), "row": sb[2], "kind": sb[1], "share": sb[3], "dA": dA, "h": hh, "h_minus": hm,
                     "rel": err})
    rec["6_single_facets_host"] = {"facets": int(len(k6)), "worst_rel": worst6, "tau": CHECK6_TAU,
                                   "step": "dA = tau p_r*/|G_r*k|, h = ln(1 + dA/A), secant of exp", "per_facet": rec6,
                                   "passed": worst6 <= 1e-4}
    ok &= rec["6_single_facets_host"]["passed"]
    # 7 G vs a GPU finite difference (a few facets; NG rows as before): the absolute step dA_k = CHECK7_TAU p_r* /
    #   |G_r*k| at the facet's strongest NG row through the override (the pressures are exactly linear in A)
    worst7, k7 = 0.0, rng.choice(len(idx), size=min(4, len(idx)), replace=False)
    rec7 = []
    for k in k7:
        sb = strongest(k, kinds=("ng",))
        if sb is None:
            rec7.append({"facet": int(idx[k]), "skipped": "zero NG column"})
            continue
        dA = CHECK7_TAU * sb[0]
        g_old = math.exp(hprev[k])
        g_new = g_old * (obj.engine.A0[k] + dA) / obj.engine.A0[k]
        pn, pr = eng.forward_override(idx[k:k + 1], np.array([g_new]))
        dA_eff = obj.engine.A0[k] * (g_new / g_old) - obj.engine.A0[k]
        rows_ = np.flatnonzero(np.abs(obj.engine.G_ng[:, k]) > 0)
        fdn = (pn - p0_ng) / dA_eff
        err = float(np.max(np.abs(fdn[rows_] - obj.engine.G_ng[rows_, k])) / np.max(np.abs(obj.engine.G_ng[rows_, k])))
        rows_r = np.flatnonzero(np.abs(obj.engine.G_rg[:, k]) > 0)
        err_rg = (float(np.max(np.abs((pr - p0_rg)[rows_r] / dA_eff - obj.engine.G_rg[rows_r, k])) /
                        np.max(np.abs(obj.engine.G_rg[rows_r, k]))) if len(rows_r) else None)
        worst7 = max(worst7, err)
        rec7.append({"facet": int(idx[k]), "row": sb[2], "share": sb[3], "dA": dA, "dA_over_A": dA / obj.engine.A0[k],
                     "rel": err, "rel_rg_information": err_rg})
    rec["7_G_vs_gpu_fd"] = {"facets": int(len(k7)), "worst_rel": worst7, "tau": CHECK7_TAU,
                            "step": "dA = tau p_r*/|G_r*k| (strongest NG row), override g -> g (A + dA)/A",
                            "per_facet": rec7, "passed": worst7 <= 1e-6}
    ok &= rec["7_G_vs_gpu_fd"]["passed"]
    # 8 zero sensitivity
    rec["8_zero_sensitivity"] = pf.zero_sensitivity(obj, st0)
    log("checks 2 adjoint %.1e, 3 J v vs FD %.1e (grain A %.1e B %.1e), 6 single facets %.1e, 7 G vs GPU FD %.1e, 8 "
        "zero-sensitivity parameters %d" % (rec["2_adjoint"]["rel"], r3, rec["3_grain_A"]["rel"], rec["3_grain_B"]["rel"],
                                           worst6, worst7, rec["8_zero_sensitivity"]["zero_sensitivity"]))
    rec["passed"] = bool(ok)
    return rec, bool(ok)


# ------------------------------------------------------------------------------------------------ products
def write_tables(out, ctx, obj, st_fit, st, pNG_gpu, pRG_gpu, utc):
    """COPS_model.dat (every row recomputed on the GPU at the final delta + grain; W10x5 for the late observations),
    late_window_model.dat and pressure_comparison.txt (positive_finetune.write_tables' formats)."""
    rows = obj.rows
    out = Path(out)
    g = st.g
    pNG_w, pRG_w = rows.window(pNG_gpu), rows.window(pRG_gpu)
    pNG, pRG = pNG_w + g, pRG_w
    lines = ctx.files["COPS_model.dat"].read_text().splitlines(True)
    ip2row = {int(p): i for i, p in enumerate(rows.ipoint)}
    note = ["# SlopeFineTuning (model_ext/slope_finetune.py, %s): pNG/pRG of EVERY row recomputed on the GPU with the "
            "late-facet fine tuning (A = max(F, m L exp(h + delta)), finetune_state.npz) of the fit %s; grain A %.9e Pa, "
            "B %.9e Pa/s refitted; late observations: W10x5(surface) + grain(t_obs); grid rows: instantaneous surface "
            "+ grain (see ClipRefitGridRows below; late_window_model.dat lists them); every other column and header "
            "line is the fit's\n" % (utc, ctx.fit_dir, st.z[-2] * rows.A_UNIT, st.z[-1] * rows.B_UNIT)]
    new, tok = pf.replace_columns(lines, (pNG, pRG), ip2row, note)
    with (out / "COPS_model.dat").open("x") as f:
        f.write("".join(new))
    lay = ctx.layout
    stop = {int(r): float(t) for r, t in zip(lay.late_rows, lay.stop_et)}
    lr = sorted([int(x) for x in lay.grid_rows] + [int(x) for x in lay.late_rows], key=lambda i: rows.et[i])
    with (out / "late_window_model.dat").open("x") as f:
        f.write("# late window of the fine tuning (model_ext/slope_finetune.py, %s) of %s: every 1-s grid row (model "
                "only) and every observation from %s UTC on\n" % (utc, ctx.fit_dir, lay.summary["late_start_utc"]))
        f.write("# kind grid: surface_* = instantaneous surface model at et; kind late_obs: surface_* = W10x5 of the "
                "surface model ending at the L2 STOP time (stop_et); grain_NG_Pa = grain(et) (NG only, unaveraged);\n"
                "# table_pNG/pRG = COPS_model.dat (grid: instantaneous surface + grain; late_obs: W10x5 + grain); "
                "re-average only surface_* of grid rows, never table_pNG\n")
        f.write("#iPoint UTC et kind stop_et surface_pNG_Pa surface_pRG_Pa grain_NG_Pa table_pNG_Pa table_pRG_Pa "
                "refNG_Pa refRG_Pa fitNG fitRG\n")
        for i in lr:
            ip = int(rows.ipoint[i])
            t = tok[ip]
            kind = "late_obs" if i in stop else "grid"
            f.write("%d %s %.6f %s %s %.9e %.9e %.9e %s %s %s %s %d %d\n" % (
                ip, t[4], rows.et[i], kind, "%.6f" % stop[i] if i in stop else "-", pNG_w[i], pRG_w[i],
                g[i], t[0], t[1], t[2], t[3], int(rows.useNG[i]), int(rows.useRG[i])))
    with (out / "pressure_comparison.txt").open("x") as f:
        f.write("# per row of the base: the tuned fit (delta = 0, the fit's grain; recomputed here) and the fine "
                "tuning (GPU at the final delta); kind obs | late_obs (W10x5 + grain) | grid (model only); ln_* = "
                "ln(model / reference) where the reference is > 0 (nan otherwise); pNG includes the grain\n")
        f.write("#iPoint UTC et kind useNG useRG refNG_Pa refRG_Pa pNG_fit_Pa pNG_tuned_Pa pRG_fit_Pa pRG_tuned_Pa "
                "ln_NG_fit ln_NG_tuned ln_RG_fit ln_RG_tuned\n")
        lnf = lambda p, r: math.log(p / r) if r > 0 and p > 0 else float("nan")  # noqa: E731
        for i in range(rows.n):
            kind = "grid" if rows.is_grid[i] else ("late_obs" if rows.is_late[i] else "obs")
            f.write("%d %s %.6f %s %d %d %.6e %.6e %.9e %.9e %.9e %.9e %.6f %.6f %.6f %.6f\n" % (
                rows.ipoint[i], ctx.utc[i], rows.et[i], kind, int(rows.useNG[i]), int(rows.useRG[i]),
                rows.refNG[i], rows.refRG[i], st_fit.pNG[i], pNG[i], st_fit.pRG[i], pRG[i],
                lnf(st_fit.pNG[i], rows.refNG[i]), lnf(pNG[i], rows.refNG[i]),
                lnf(st_fit.pRG[i], rows.refRG[i]), lnf(pRG[i], rows.refRG[i])))
    return pNG, pRG


def write_state(out, ctx, obj, st, idx, late, meta, pool):
    n = ctx.nf
    delta = np.zeros(n)
    delta[idx] = st.z[:len(idx)]
    lb = np.full(n, np.nan)
    lb[idx] = obj.lb_delta
    held_next = np.empty(n)
    hp = ctx.held_prev
    elementwise(lambda i0, i1: hp[i0:i1] + delta[i0:i1], n, held_next, pool)
    A_final = ctx.A_fit.copy()
    A_final[idx] = obj.engine.A(st.s)
    path = Path(out) / "finetune_state.npz"
    with path.open("xb") as f:
        np.savez_compressed(f, delta=delta, held_prev=hp, held_next=held_next, A_final=A_final, lower_bound=lb,
                            seen=obj.seen, late=late, grain_A_uPa=np.array(st.z[-2]),
                            grain_B_uPa_per_100s=np.array(st.z[-1]),
                            meta_json=np.array(json.dumps(jsonable(meta), sort_keys=True)))
    return path, held_next


def write_result_txt(out, ctx, obj, s_init, s_final, opt, rec, utc, seg):
    with (Path(out) / "finetune_result.txt").open("x") as f:
        f.write("# late-facet fine tuning (model_ext/slope_finetune.py, %s) of %s\n" % (utc, ctx.fit_dir))
        f.write("# A_e = max(F_e, m_b L_e exp(h_e + delta_e)), delta_e >= ln(F_e / s_e); %d parameter facets + grain A, "
                "B; lambda 0; costs: total = data + crossing (the fit's final_cost = data)\n" % obj.npar)
        lo = rec["seen"].get("late_only") or {}
        f.write("facets late_only  WINDOW_MIN %g  window %s .. %s (et %.3f .. %.3f); late facets %d, at the floor in "
                "the fit (not tuned) %d, parameters %d\n" % (
                    lo.get("window_min", float("nan")), lo.get("t_start_utc"), lo.get("t_end_utc"),
                    lo.get("t_start_et", float("nan")), lo.get("t_end_et", float("nan")), rec["seen"]["late_facets"],
                    rec["seen"]["late_facets_at_floor"], obj.npar))
        for key in pf.RESULT_KEYS:
            a, b = s_init.get(key), s_final.get(key)
            fmt = (lambda v: "%.10e" % v if isinstance(v, float) else str(v))  # noqa: E731
            f.write("%-46s initial %s  final %s\n" % (key, fmt(a), fmt(b)))
        f.write("fit_final_cost(fit_result.txt) %s  fit_manifest cost %.15e data %.15e crossing %.15e\n" % (
            ctx.fit_result_txt.get("final_cost"), ctx.fit_costs["cost"], ctx.fit_costs["data_cost"],
            ctx.fit_costs["crossing_penalty_cost"]))
        reason = pf.stop_reason(opt["termination"])
        f.write("stop %s\ncomplete %s\n" % (reason, "no (resubmit: the next segment continues)" if reason == "wall_time"
                                            else "yes"))
        f.write("termination %s\nmessage %s\nsegment %d\nnfev %s\nnfev_total %d\nnjev %s\noptimality %s\n"
                "wall_seconds %.1f\n" % (opt["termination"], opt["message"], seg["segment"], opt["nfev"],
                                         seg["nfev_total"], opt["njev"], opt["optimality"], opt["seconds"]))
        if "active_mask" in opt:
            f.write("active_bounds lower %d upper %d grain %s\n" % (opt["active_mask"]["lower"],
                                                                   opt["active_mask"]["upper"],
                                                                   opt["active_mask"]["grain"]))


# ------------------------------------------------------------------------------------------------ main
def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--facets", choices=("late_only",), required=True)
    ap.add_argument("--window-min", type=float, required=True)
    ap.add_argument("--devices", default="0")
    ap.add_argument("--max-nfev", type=int, default=None)
    ap.add_argument("--lsmr-maxiter", type=int, default=50)
    ap.add_argument("--ftol", type=float, default=1e-10)
    ap.add_argument("--xtol", type=float, default=1e-10)
    ap.add_argument("--gtol", type=float, default=1e-10)
    ap.add_argument("--x-scale-grain", type=float, default=1.0)
    ap.add_argument("--l-threads", type=int, default=None)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--resume-from", default=None)
    ap.add_argument("--wall-time", default=None)
    ap.add_argument("--wall-margin", type=float, default=300.0)
    ap.add_argument("--save-seconds", type=float, default=900.0)
    ap.add_argument("--segment", type=int, default=1)
    ap.add_argument("--nfev-before", type=int, default=0)
    args = ap.parse_args(argv)
    if not args.validate_only and args.max_nfev is None:
        ap.error("--max-nfev is required unless --validate-only")
    if args.window_min <= 0:
        ap.error("--window-min must be > 0")
    if args.segment > 1 and not args.resume_from:
        ap.error("--segment > 1 needs --resume-from")
    return args


def run(args, engine_factory=None, backend=None, command=None):
    """The fine tuning (stand-alone main, or in process: `backend` = an already bound backend of the fit's rows and
    basis; `engine_factory(ctx)` replaces the GPU engine, e.g. a CPU stand-in in the tests). Returns the exit code."""
    started = time.time()
    t = {}
    pool = ThreadPoolExecutor(args.l_threads) if args.l_threads else None
    t0 = time.time()
    ctx = SlopeFitContext(args.fit, pool)
    t["context_s"] = time.time() - t0
    from cops_cuda import wall_time
    deadline = wall_time.resolve_deadline(args.wall_time, None, args.wall_margin, start=started)
    out = Path(args.out).resolve()
    if out.exists() and any(p for p in out.iterdir() if not p.name.startswith("provenance-")):
        raise SystemExit("--out %s exists with products" % out)
    out.mkdir(parents=True, exist_ok=True)
    devices = [int(d) for d in args.devices.split(",") if d.strip() != ""]
    record = {"schema": "slope_finetune", "parameterization": PARAMETERIZATION_ID, "policy": POLICY,
              "created_utc": utc_now(), "command": command, "host": socket.gethostname(),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "fit": {"directory": str(ctx.fit_dir), "identities": ctx.identity,
                      "parameterization": ctx.fit_parameterization, "costs": ctx.fit_costs,
                      "fit_result_txt": ctx.fit_result_txt, "model": ctx.model_cfg, "x": ctx.x.tolist(),
                      "names": ctx.names, "run": str(ctx.run), "point_table": file_identity(ctx.point_table)},
              "settings": jsonable(vars(args)), "threads": {"l_threads": args.l_threads, "pool": bool(pool)},
              "modules": {p: file_identity(HERE / p) for p in ("slope_finetune.py", "slope_mult.py",
                                                               "positive_finetune.py", "positive_problem.py",
                                                               "grain_crossing.py")}}
    tracker = pf.Tracker(deadline, args.save_seconds)
    eng = None
    try:
        t0 = time.time()
        eng = engine_factory(ctx) if engine_factory is not None else GpuModelEngine(
            ctx, devices, backend=backend, l_threads=args.l_threads)
        t["engine_setup_s"] = time.time() - t0
        rows = ctx.rows
        own, late_used = rows.objective_points()
        t0 = time.time()
        late, late_rec, support = pf.late_only_facets(_SeenAdapter(eng), rows, _CtxAdapter(ctx), args.window_min, own,
                                                      late_used)
        t["seen_s"] = time.time() - t0
        at_floor = late & ~ctx.active_fit
        params = late & ctx.active_fit
        idx = np.flatnonzero(params)
        record["seen"] = {"late_only": late_rec, "late_facets": int(late.sum()),
                          "late_facets_at_floor": int(at_floor.sum()), "parameters": int(len(idx)),
                          "parameter_area_fraction": float(ctx.area[params].sum() / ctx.area.sum()),
                          "rule": pf.FACET_RULES["late_only"] + "; minus the late facets at the floor in the fit "
                                                                "(dA/ddelta = 0 there)"}
        log("late_only: %d facets, %d at the floor in the fit (not tuned), %d parameters" % (
            int(late.sum()), int(at_floor.sum()), len(idx)))
        if not len(idx):
            raise SystemExit("no parameter facet")
        t0 = time.time()
        p0_ng, p0_rg = eng.forward()
        t["forward_fit_state_s"] = time.time() - t0
        t0 = time.time()
        G_ng, G_rg = eng.sensitivities(idx)
        t["sensitivities_s"] = time.time() - t0
        log("p at the fit's state and G = dp/dA for %d facets (%.1f s, %.1f s)" % (len(idx), t["forward_fit_state_s"],
                                                                                    t["sensitivities_s"]))
        host = HostLinearEngine(p0_ng, p0_rg, G_ng, G_rg, ctx.A_fit[idx], ctx.s_start, idx)
        obj = pf.Objective(rows, host, ctx.s_start, ctx.s_start, ctx.area, ctx.A_orig, params, 0.0, ctx.model_cfg[
            "floor_ratio"], tracker=tracker, fit_identity={"fit_directory": str(ctx.fit_dir),
                                                          "fit_restart_sha256": ctx.identity["fit_restart.npz"]["sha256"]})
        z_fit = obj.start_vector(ctx.xA, ctx.xB)
        z0 = z_fit.copy()
        if args.resume_from:
            with np.load(args.resume_from, allow_pickle=False) as z:
                m_ = json.loads(str(z["meta_json"]))
                if Path(m_["fit"]["fit_directory"]).resolve() != ctx.fit_dir:
                    raise SystemExit("--resume-from tuned another fit")
                z0[:len(idx)] = np.asarray(z["delta"], float)[idx]
                z0[-2], z0[-1] = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
            record["resumed_from"] = file_identity(args.resume_from)
        st_fit = obj.evaluate(z_fit)
        s_init = obj.summary(obj.evaluate(z0))
        if args.validate or args.validate_only:
            t0 = time.time()
            held_g = np.exp(ctx.held_prev[idx]) if (ctx.model_cfg.get("held") or {}) else None
            vrec, passed = validate(obj, ctx, eng, z_fit, idx, held_g, support, args.seed)
            t["validation_s"] = time.time() - t0
            (out / "validation.json").write_text(json.dumps(jsonable(vrec), indent=2) + "\n")
            record["validation"] = {k: v.get("passed") if isinstance(v, dict) else v for k, v in vrec.items()}
            if not passed:
                raise SystemExit("validation failed: %s" % record["validation"])
            if args.validate_only:
                return 0
        log("start: total cost %.12e (data %.12e, crossing %.6e); %d parameters + grain" % (
            s_init["cost_total"], s_init["cost_data"], s_init["cost_crossing"], len(idx)))
        t0 = time.time()
        st, opt = pf.optimize(obj, z0, args)
        t["optimize_s"] = time.time() - t0
        s_final = obj.summary(st)
        t0 = time.time()
        g_final = np.exp(ctx.held_prev[idx] + st.z[:len(idx)])
        pNG_gpu, pRG_gpu = eng.forward_override(idx, g_final)
        t["final_gpu_forward_s"] = time.time() - t0
        relf = float(max(np.max(np.abs(st.pNG_s / pNG_gpu - 1.0)), np.max(np.abs(st.pRG_s / pRG_gpu - 1.0))))
        utc = utc_now()
        nfev_total = int(args.nfev_before) + int(opt["nfev"])
        seg = {"segment": int(args.segment), "nfev_total": nfev_total}
        reason = pf.stop_reason(opt["termination"])
        meta = {"kind": "slope_finetune_state", "parameterization": PARAMETERIZATION_ID, "created_utc": utc,
                "fit": {"fit_directory": str(ctx.fit_dir), "fit_restart_sha256": ctx.identity["fit_restart.npz"]["sha256"],
                        "parameterization": ctx.fit_parameterization, "model": ctx.model_cfg},
                "facets": "late_only", "window_min": args.window_min, "parameters": int(len(idx)),
                "costs_initial": s_init, "costs_final": s_final, "stop": reason, "complete": reason != "wall_time",
                "nfev_total": nfev_total, "segment": int(args.segment),
                "held_next_definition": "held_prev + delta: the next fit's A = max(F, m_b L exp(held_next))"}
        spath, held_next = write_state(out, ctx, obj, st, idx, late, meta, pool)
        pNG_t, pRG_t = write_tables(out, ctx, obj, st_fit, st, pNG_gpu, pRG_gpu, utc)
        write_result_txt(out, ctx, obj, s_init, s_final, opt, record, utc, seg)
        record.update(optimizer=jsonable(opt), costs_initial=s_init, costs_final=s_final,
                      final_host_vs_gpu_max_rel=relf, segment=seg, stop=reason,
                      products={p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
                      engine_seconds={"host": host.seconds,
                                      "gpu": getattr(eng, "seconds", None), "gpu_calls": getattr(eng, "calls", None)},
                      timings=t, wall_s=time.time() - started,
                      backend=eng.identity() if hasattr(eng, "identity") else None)
        log("final: total %.12e (data %.12e, crossing %.6e) after %d evaluations (%s); host vs GPU at the solution max "
            "rel %.2e; products in %s" % (s_final["cost_total"], s_final["cost_data"], s_final["cost_crossing"],
                                          opt["nfev"], opt["termination"], relf, out))
        return 75 if reason == "wall_time" else 0
    finally:
        record.setdefault("timings", t)
        record["wall_s"] = time.time() - started
        (out / "fine_tuning.json").write_text(json.dumps(jsonable(record), indent=2) + "\n") if not \
            (out / "fine_tuning.json").exists() else None
        if eng is not None and hasattr(eng, "close") and backend is None:
            eng.close()
        if pool is not None:
            pool.shutdown()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    command = [str(Path(__file__).resolve())] + argv
    args = parse_args(argv)
    if any(k.startswith("COPS_POS_") and v for k, v in os.environ.items()):
        raise SystemExit("slope_finetune: COPS_POS_* set; unset them")
    return run(args, command=command)


if __name__ == "__main__":
    sys.exit(main())
