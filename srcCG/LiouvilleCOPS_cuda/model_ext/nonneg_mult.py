# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_mult.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_mult.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 560dbac48ea6912f09f021e7763b8865017e11679421cfb289c4d98a8c904c01
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Non-negative-by-construction H2O activity with slope-bin multipliers (user decisions 2026-11-11:
#   "implement option 1 and 2 ...; use 6 cliff templates as multipliers"; D2's ft_c1..ft_c6 held at 1, no floor):
#   softplus A = m_b s ln(1 + exp(L/s)) exp(h), s = 1e-5 A_orig (hard clip where A_orig = 0), and squared A = m_b (A_s
#   (Y b)^2 + A_s sum ft_c*) exp(h); host fields and column factors by chunked bitwise threads, held term per stage on
#   the host, GPU patch of the host-activity path.
# USED BY: nonneg_fit.py, nonneg_finetune.py, chain_onejob.py --kind nonneg (snapshot
#   model_ext_frozen/nonneg_20261111T044920Z; production chains 124162 softplus, 124165 squared).
# ==== end of header ====
"""Non-negative H2O activity models in the slope-multiplier framework (user decisions 2026-11-11, handoff.md section 10,
"softplus and squared-harmonics models (D2 step back)"): model math, CPU stand-in and the GPU patch of the
host-activity path. Fit-level part: nonneg_problem.py; driver: nonneg_fit.py; squared start: nonneg_start.py; fine
tuning: nonneg_finetune.py.

Local, untracked model extension: nothing in the git tree and no existing model_ext file is edited (slope_mult.py is
imported for its helpers only: bins, chunked matvec, BLAS info, the seen kernel). The patch is installed ONLY by the
nonneg drivers (install_gpu); it refuses to coexist with the slope_mult or positive_activity patch.

User decisions (2026-11-11)
* "implement option 1 and 2 ...; use 6 cliff templates as multipliers; add fine tuning cycles at the end; submit
  fitting + tuning in a single job and use multiple thread parallelization on the CPU".
* multipliers m_b: the 7 effective-gravitational-slope bins of slope_mult (bundle slope_bins.npz), m_[0,40) = 1 held,
  six free >= 0, start 1; NO floor in either model.
* D2's six fine-tuning templates ft_c1..ft_c6 held at amplitude 1 (D2's template files); NO cliff_ge_60deg.
* softplus: A_e = m_b(e) g_e(L_e) exp(h_e), L = Y a + A_s sum_k c_k T_k (25 signed degree-<=4 harmonics),
  g_e(L) = s_e ln(1 + exp(L / s_e)), s_e = 1e-5 A_orig,e (the run's srcCG map, elements.bin AH2O, as the slope model's
  floor); where A_orig,e = 0 (0.71 % of the facets): g_e(L) = max(0, L) ("Hard clip max(0, L) there").
  Derivative: g'_e(L) = sigmoid(L / s_e) (s > 0), [L > 0] (s = 0).
* squared: A_e = m_b(e) (A_s (Y b)_e^2 + A_s sum_k c_k T_k) exp(h_e) (templates outside the square). The fitter's
  harmonic coefficients a are b A_s (b dimensionless): A_s (Y b)^2 = (Y a)^2 / A_s.
* h = the held fine-tuning term (0 in cycle 1); A_s = the fitter's enhancement scale.

Implementation (why): the whole activity A and its column factors are formed on the HOST of the parent process (the
fitter), by chunked Python threads (bitwise: every 12.4M-facet step is a row-chunked matvec or an elementwise ufunc;
set_L_threads), with the held term of the CURRENT STAGE (set_stage_held, per stage, never in the frozen worker
configuration: requirement A of model_ext/logs/onejob_stage_interface_3a.md, option (a)):
* act_lin := A (>= 0). The unchanged evaluator then forms af = max(A, 0) = A and active = A > 0 (where A = 0 every
  column factor below is 0 as well, so the original mask removes nothing). The local workers read A from the
  transport's shared activity buffer (unchanged path).
* columns: the evaluator's kernels compute sum_r w_r active(f_r) colbasis[c][f_r]; the patch builds colbasis on the
  device as D Y_j (harmonics), [bin = k] P (multipliers), E A_s T_k (templates) from three per-facet host factors
      softplus: D = m g' exp(h), E = D,          P = g(L) exp(h)
      squared:  D = m (2 u / A_s) exp(h), E = m exp(h), P = (u^2 / A_s + A_s sum c T) exp(h)    (u = Y a)
  that the parent writes into two shared-memory arrays it owns (D; E and P: names in COPS_NONNEG_SHM) before every
  Jacobian; Y^T and A_s T^T of the device facets are resident on the device (uploaded once per bound basis).
* the fine tuning's per-facet override (exp(h + delta) at its facets) is applied in the same host fields (parent
  only; set_override).
* every other step (pressures, T and kappa columns) is the unchanged kernel's with af = A.
"""

import atexit
import hashlib
import inspect
import json
import os
import socket
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import slope_mult as sm             # helpers only (bins, chunk bounds, BLAS info, the seen kernel); nothing installed

HERE = Path(__file__).resolve().parent
BIN_NAMES, BIN_EDGES_DEG, BIN_LABELS = sm.BIN_NAMES, sm.BIN_EDGES_DEG, sm.BIN_LABELS
FIXED_BIN, FIXED_VALUE = sm.FIXED_BIN, sm.FIXED_VALUE
MODELS = ("softplus", "squared")
S_RATIO = 1e-5                      # softplus width s_e = 1e-5 A_orig,e (user 2026-11-11)
TEST_S_RATIOS = (0.0,)              # test only: s = 0 everywhere = the hard clip max(0, L) (softplus -> clip check)
ENV_FROZEN = "COPS_NONNEG_FROZEN_CONFIG"
ENV_SHM = "COPS_NONNEG_SHM"
CONFIG_VERSION = 1
BYTES_PER_FACE = 296                # device per facet: Y^T 25x8, A_s T^T up to 6x8, D E P 3x8, bins 4, margin
L_ALIGN_ROWS = sm.L_ALIGN_ROWS
POLICY = {
    "softplus": ("user 2026-11-11: A = m_b s ln(1 + exp(L / s)) exp(h), L = Y a + A_s sum ft_c1..ft_c6 (held at 1), "
                 "s = 1e-5 A_orig per facet, hard clip max(0, L) where A_orig = 0; 7 slope-bin multipliers (0-40 held "
                 "at 1, six free >= 0); no floor; no cliff_ge_60deg"),
    "squared": ("user 2026-11-11: A = m_b (A_s (Y b)^2 + A_s sum ft_c1..ft_c6) exp(h) (templates outside the square, "
                "held at 1; a = A_s b); 7 slope-bin multipliers (0-40 held at 1, six free >= 0); no floor; no "
                "cliff_ge_60deg")}
WORKER_REQUESTS = {"nonneg_seen": "nonneg_seen_rows", "nonneg_sens": "nonneg_sensitivities",
                   "nonneg_detach": "nonneg_detach"}
LOG = "nonneg_mult:"

_STATE = {"installed": False, "config": None, "originals": {}, "role": None,
          "pool": None, "threads": None, "mode": None,
          "held": None, "override": None, "A_orig": None, "A_orig_sha256": None,
          "shm": {}, "timing": {"fields_calls": 0, "fields_s": 0.0, "first_s": None, "max_s": 0.0},
          "checks": {"colbasis": None}}


def say(text):
    print("%s %s" % (LOG, text), flush=True)


sha256_bytes, canonical, blas_info = sm.sha256_bytes, sm.canonical, sm.blas_info


def module_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ threads
def set_L_threads(n, mode="chunked"):
    """N Python threads for every 12.4M-facet host step of THIS process (Y a and B c row chunks, elementwise
    transforms), OpenBLAS at 1 thread (bitwise to 1 thread). Idempotent: the same (n, mode) keeps the pool; another
    shuts the old pool down first (one-job requirement B). n = None or 1: no pool. Returns the record."""
    if mode != "chunked":
        raise ValueError("nonneg_mult: only the bitwise chunked mode is provided")
    n = None if n is None else int(n)
    if n is not None and n < 1:
        raise ValueError("threads must be >= 1")
    if _STATE["threads"] == n and _STATE["mode"] == mode and (n is None or _STATE["pool"] is not None):
        return {"threads": n, "mode": mode, "reused": True}
    if _STATE["pool"] is not None:
        _STATE["pool"].shutdown(wait=True)
        _STATE["pool"] = None
    info = blas_info()
    if n is not None and n > 1:
        if info["threads"] not in (1, None):
            raise RuntimeError("nonneg_mult: chunked host steps need numpy's OpenBLAS at 1 thread (got %r)"
                               % info["threads"])
        _STATE["pool"] = ThreadPoolExecutor(max_workers=n, thread_name_prefix="nonneg")
    _STATE["threads"], _STATE["mode"] = n, mode
    return {"threads": n, "mode": mode, "reused": False, "openblas_threads": info["threads"],
            "align_rows": L_ALIGN_ROWS, "bitwise_by_construction": True}


def thread_info():
    return {"threads": _STATE["threads"], "mode": _STATE["mode"], "pool": _STATE["pool"] is not None}


def matvec(M, v):
    """M @ v (M C-contiguous float64 (n, k)); row chunks on the pool when set (bitwise: test_nonneg)."""
    v = np.ascontiguousarray(v, dtype=np.float64)
    pool, n = _STATE["pool"], _STATE["threads"]
    if pool is None or n is None or n <= 1:
        return M @ v
    return sm.chunked_matvec(M, v, pool, n)


def elementwise(fn, nrows, *outs):
    """fn(c0, c1) fills outs[..][c0:c1] for row chunks (on the pool when set). Elementwise numpy: bitwise."""
    pool, n = _STATE["pool"], _STATE["threads"]
    if pool is None or n is None or n <= 1:
        fn(0, nrows)
        return
    futures = [pool.submit(fn, c0, c1) for c0, c1 in sm.chunk_bounds(nrows, n)]
    for f in futures:
        f.result()


# ------------------------------------------------------------------------------------------------ model math
def softplus(L, s):
    """(g, g') per facet: s > 0: g = s (max(x, 0) + log1p(exp(-|x|))), x = L / s, g' = sigmoid(x) (stable form);
    s = 0: g = max(0, L), g' = [L > 0]. Float64, elementwise."""
    L = np.asarray(L, np.float64)
    s = np.asarray(s, np.float64)
    pos = s > 0.0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        x = np.where(pos, L / np.where(pos, s, 1.0), 0.0)
        ax = np.abs(x)
        e = np.exp(-ax)
        g_s = s * (np.maximum(x, 0.0) + np.log1p(e))
        sig = np.where(x >= 0.0, 1.0 / (1.0 + e), e / (1.0 + e))
    g = np.where(pos, g_s, np.maximum(L, 0.0))
    d = np.where(pos, sig, (L > 0.0).astype(np.float64))
    return g, d


def fields(model, u, TS, s, mrow, gh, scale, want_cols):
    """Host fields of one chunk: (A, D, E, P) (D, E, P None unless want_cols). u = Y a, TS = A_s sum c T, s the softplus
    width, mrow = m[bin], gh = exp(h) (or None), scale = A_s. The one definition used by the GPU patch, the CPU stand-in,
    the products and the fine tuning (same operations, same order)."""
    if model == "softplus":
        L = u + TS
        q, dq = softplus(L, s)
    elif model == "squared":
        q = (u * u) / scale + TS
        dq = None
    else:
        raise ValueError("unknown model %r" % model)
    P = q if gh is None else q * gh
    A = mrow * P
    if not want_cols:
        return A, None, None, None
    if model == "softplus":
        D = mrow * dq
        if gh is not None:
            D = D * gh
        E = D
    else:
        D = mrow * ((2.0 * u) / scale)
        E = mrow.copy()
        if gh is not None:
            D = D * gh
            E = E * gh
    return A, D, E, P


def host_fields(model, Y, a, Bft, c, scale, bins, m, s, gh, want_cols, override=None, out=None):
    """The model's per-facet fields for every facet (chunked): returns dict(A, D, E, P, u, TS, q). override:
    (facets, g_values) replacing gh at those facets (fine tuning)."""
    n = Y.shape[0]
    u = matvec(Y, a)
    if Bft is not None and Bft.shape[1]:
        TS = scale * matvec(Bft, c)
    else:
        TS = np.zeros(n)
    mrow = np.asarray(m, np.float64)[np.asarray(bins).astype(np.intp)]
    if override is not None:
        fac, gv = override
        gh = (np.ones(n) if gh is None else np.array(gh, dtype=np.float64, copy=True))
        gh[np.asarray(fac, np.int64)] = np.asarray(gv, np.float64)
    o = out if out is not None else {}
    A = np.empty(n)
    D = np.empty(n) if want_cols else None
    E = np.empty(n) if want_cols else None
    P = np.empty(n) if want_cols else None

    def fn(c0, c1):
        r = fields(model, u[c0:c1], TS[c0:c1], None if s is None else s[c0:c1], mrow[c0:c1],
                   None if gh is None else gh[c0:c1], scale, want_cols)
        A[c0:c1] = r[0]
        if want_cols:
            D[c0:c1], E[c0:c1], P[c0:c1] = r[1], r[2], r[3]
    elementwise(fn, n)
    o.update(A=A, D=D, E=E, P=P, u=u, TS=TS)
    return o


# ------------------------------------------------------------------------------------------------ configuration
def process_config(model, slope_meta, ft_templates, s_ratio=S_RATIO, slope_bundle=None):
    """The constant configuration every process shares (no held term, no stage data): model, s ratio, bins, A_orig
    digest, mesh, the ft template names and bundle digests. ft_templates: [{"name", "bundle": {"path", "sha256"}}]."""
    if model not in MODELS:
        raise ValueError("nonneg_mult: model %r not in %s" % (model, MODELS))
    smb = (slope_meta or {}).get("slope_multiplier_bins")
    if not smb:
        raise ValueError("nonneg_mult: %s is not a slope_bins.py bundle" % slope_bundle)
    if list(smb.get("names") or []) != list(BIN_NAMES) or [float(v) for v in smb["edges_deg"]] != list(BIN_EDGES_DEG):
        raise ValueError("nonneg_mult: bundle bins are not the slope model's")
    f = float(s_ratio)
    if model == "softplus" and not (f == S_RATIO or f in TEST_S_RATIOS):
        raise ValueError("nonneg_mult: softplus s ratio %r is not the decision %r (test only: %s)" % (f, S_RATIO,
                                                                                                    TEST_S_RATIOS))
    if model == "squared":
        f = None
    return {"version": CONFIG_VERSION, "model": model, "s_ratio": f,
            "test_only": model == "softplus" and f != S_RATIO,
            "bin_names": list(BIN_NAMES), "bin_edges_deg": list(BIN_EDGES_DEG), "fixed_bin": FIXED_BIN,
            "nfaces": int(smb["nfaces"]), "bins_sha256": smb["bins_sha256"],
            "bin_counts": [int(v) for v in smb["counts"]], "ah2o_sha256": smb["ah2o_sha256"],
            "mesh_identity": smb["mesh_identity"],
            "templates": [{"name": t["name"], "bundle_sha256": t["bundle"]["sha256"]} for t in ft_templates],
            "n_templates": len(ft_templates)}


def installed_config():
    return dict(_STATE["config"]) if _STATE["config"] else None


def n_bins_templates(cfg):
    return len(cfg["bin_names"]), int(cfg.get("n_templates", 0))


def set_A_orig(ah2o):
    """Register A_orig (elements.bin AH2O) in this process (digest checked against the installed config); the softplus
    width s = s_ratio A_orig is formed from it. Idempotent."""
    a = np.ascontiguousarray(np.asarray(ah2o, dtype=np.float64))
    d = sha256_bytes(a)
    cfg = _STATE["config"]
    if cfg is not None and d != cfg["ah2o_sha256"]:
        raise ValueError("nonneg_mult: A_orig digest %s differs from the configured %s" % (d[:16], cfg["ah2o_sha256"][:16]))
    if _STATE["A_orig_sha256"] != d:
        _STATE["A_orig"], _STATE["A_orig_sha256"] = a, d
        _STATE.pop("s_cache", None)
    return d


def width(cfg):
    """s = s_ratio A_orig (softplus; cached); None for the squared model."""
    if cfg["model"] != "softplus":
        return None
    if _STATE["A_orig"] is None:
        raise RuntimeError("nonneg_mult: A_orig not registered (set_A_orig)")
    key = (_STATE["A_orig_sha256"], cfg["s_ratio"])
    c = _STATE.get("s_cache")
    if c is None or c[0] != key:
        s = float(cfg["s_ratio"]) * _STATE["A_orig"]
        _STATE["s_cache"] = c = (key, s)
    return c[1]


# ------------------------------------------------------------------------------------------------ the stage's held term
def load_held(path, key, nfaces):
    with np.load(str(path), allow_pickle=False) as z:
        h = np.ascontiguousarray(np.asarray(z[key], dtype=np.float64))
    if h.shape != (int(nfaces),) or not np.all(np.isfinite(h)):
        raise ValueError("nonneg_mult: held term %s[%s] has shape %s / non-finite values" % (path, key, h.shape))
    return h


def held_record(path, key, h):
    return {"path": str(Path(path).resolve()), "key": key, "sha256": sha256_bytes(h),
            "nonzero_facets": int(np.count_nonzero(h)), "min": float(h.min()), "max": float(h.max())}


def set_stage_held(path=None, key="held_next", nfaces=None, expect_sha256=None):
    """The held term of THIS stage (host side, parent process): A = m q exp(h). None clears it (cycle 1). Every stage
    sets its own at its start (nonneg_fit.main from --held, nonneg_finetune from the tuned fit); nothing of a previous
    stage survives. Returns the record (None without a held term)."""
    if path is None:
        _STATE["held"] = None
        return None
    nf = nfaces if nfaces is not None else (_STATE["config"] or {}).get("nfaces")
    h = load_held(path, key, nf)
    rec = held_record(path, key, h)
    if expect_sha256 is not None and rec["sha256"] != expect_sha256:
        raise ValueError("nonneg_mult: held term %s changed (sha256 %s, expected %s)" % (path, rec["sha256"][:16],
                                                                                        expect_sha256[:16]))
    gh = np.empty_like(h)

    def fn(c0, c1):
        gh[c0:c1] = np.exp(h[c0:c1])
    elementwise(fn, len(h))
    _STATE["held"] = dict(rec, h=h, g=gh)
    return rec


def stage_held():
    held = _STATE["held"]
    return None if held is None else {k: v for k, v in held.items() if k not in ("h", "g")}


def held_g():
    return None if _STATE["held"] is None else _STATE["held"]["g"]


def held_h():
    return None if _STATE["held"] is None else _STATE["held"]["h"]


def set_override(facets, g_values=None):
    """Fine tuning (parent, host fields): replace g = exp(h) at the given global facets by g_values for the following
    evaluations; None clears."""
    if facets is None:
        _STATE["override"] = None
        return 0
    _STATE["override"] = (np.asarray(facets, np.int64).copy(), np.asarray(g_values, np.float64).copy())
    return len(facets)


# ------------------------------------------------------------------------------------------------ per-basis host data
def basis_parts(obj, cfg):
    """(bins int8, B_ft C-contiguous (n, K)) of the bound basis of a backend / model object (cached per digest)."""
    nb, nt = n_bins_templates(cfg)
    basis = np.asarray(obj.enhancement_basis)
    key = sha256_bytes(np.ascontiguousarray(basis[:, :nb])) if not hasattr(obj, "enhancement_basis_sha256") else \
        obj.enhancement_basis_sha256
    cache = getattr(obj, "_nonneg_parts", None)
    if cache is None or cache[0] != key:
        if basis.ndim != 2 or basis.shape[1] != nb + nt:
            raise ValueError("nonneg_mult: bound basis has shape %s; the model needs %d bins + %d templates"
                             % (basis.shape, nb, nt))
        bins = sm.bins_from_basis(basis[:, :nb], nb)
        if sha256_bytes(bins) != cfg["bins_sha256"]:
            raise ValueError("nonneg_mult: the bound basis is not the configured bin assignment")
        Bft = np.ascontiguousarray(basis[:, nb:nb + nt], dtype=np.float64)
        obj._nonneg_parts = cache = (key, bins, Bft)
    return cache[1], cache[2]


def state_fields(obj, cfg, state, want_cols):
    """host_fields for a backend / model object at a ParameterState-like state (a, b) with this process's held term
    and override."""
    nb, nt = n_bins_templates(cfg)
    bins, Bft = basis_parts(obj, cfg)
    b = np.asarray(state.b, np.float64)
    t0 = time.perf_counter()
    out = host_fields(cfg["model"], obj.Y, np.asarray(state.a, np.float64), Bft, b[nb:nb + nt],
                      float(obj.enhancement_scale), bins, b[:nb], width(cfg), held_g(), want_cols,
                      override=_STATE["override"])
    dt_ = time.perf_counter() - t0
    tm = _STATE["timing"]
    tm["fields_calls"] += 1
    tm["fields_s"] += dt_
    tm["max_s"] = max(tm["max_s"], dt_)
    if tm["first_s"] is None:
        tm["first_s"] = dt_
    return out


# ------------------------------------------------------------------------------------------------ CPU stand-in
def cpu_model_class(fit_cops, cfg):
    """fit_cops.Model with this model's activity and columns (CPU; the tests' stand-in for the GPU model). b None ->
    the reference evaluate unchanged. The held term and the override are this process's (set_stage_held /
    set_override); sensitivities(T, facets): dp/dA_e."""
    base = fit_cops.Model
    nb, nt = n_bins_templates(cfg)

    class NonnegCpuModel(base):
        def __init__(self, run, header, cols, rows, elements, coef):
            super().__init__(run, header, cols, rows, elements, coef)
            set_A_orig(elements["AH2O"])

        def _temps(self, T):
            KT = fit_cops.KT
            Tmax, T75, Tmin = T
            Tl = Tmax + KT * (Tmax - T75) * self.xlit
            use_lit = self.lit & (Tl > Tmin)
            Te = np.where(use_lit, Tl, Tmin)
            return use_lit, Te, 1.0 / np.sqrt(Te)

        def sensitivities(self, T, facets):
            _, _, u = self._temps(T)
            cNG, cRG, bm = self.cNG, self.cRG, self.bm
            Bm = np.zeros((self.nfaces, len(facets)))
            Bm[np.asarray(facets, np.int64), np.arange(len(facets))] = 1.0
            ng = cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ Bm)
            rg = cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ Bm)
            return ng, rg

        def evaluate(self, a, T, kappa, b=None, enhancement_basis=None, enhancement_scale=1.0, jac=False,
                     basis_index=None, enhancement_index=None):
            if b is None:
                return base.evaluate(self, a, T, kappa, b=None, enhancement_basis=enhancement_basis,
                                     enhancement_scale=enhancement_scale, jac=jac, basis_index=basis_index,
                                     enhancement_index=enhancement_index)
            KT = fit_cops.KT

            class _Obj:
                pass
            o = _Obj()
            o.Y, o.enhancement_basis, o.enhancement_scale = self.Y, enhancement_basis, enhancement_scale
            o.enhancement_basis_sha256 = sha256_bytes(np.ascontiguousarray(np.asarray(enhancement_basis)))

            class _St:
                pass
            st = _St()
            st.a, st.b = a, b
            f = state_fields(o, cfg, st, jac)
            A = f["A"]
            use_lit, Te, u = self._temps(T)
            af = np.maximum(A, 0.0)[self.face]
            bc = lambda w: np.bincount(self.pt, weights=w, minlength=self.npt)     # noqa: E731
            cNG, cRG, bm = self.cNG, self.cRG, self.bm
            pNG_H = cNG["H2O"] * bm["H2O"] * bc(self.wn["H2O"] * af * u)
            pRG_H = cRG["H2O"] * (bc(self.wp["H2O"] * af) + bm["H2O"] * bc(self.wv["H2O"] * af * u))
            pNG_C = cNG["CO2"] * bm["CO2"] * bc(self.wn["CO2"] * u)
            pRG_C = cRG["CO2"] * (self.PhiCO2 + bm["CO2"] * bc(self.wv["CO2"] * u))
            pNG = pNG_H + kappa * pNG_C
            pRG = pRG_H + kappa * pRG_C
            if not jac:
                return pNG, pRG
            J = {}
            active = (A > 0.0)[:, None]
            bins, Bft = basis_parts(o, cfg)

            def activity_columns(Bm):
                masked = active * Bm
                JaNG = cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ masked)
                JaRG = cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ masked)
                return JaNG, JaRG
            basis_index = [] if basis_index is None else list(basis_index)
            if basis_index:
                JaNG, JaRG = activity_columns(f["D"][:, None] * self.Y[:, basis_index])
                for k, j in enumerate(basis_index):
                    J[("a", j)] = (JaNG[:, k], JaRG[:, k])
            enhancement_index = [] if enhancement_index is None else list(enhancement_index)
            if enhancement_index:
                cols = []
                for k in enhancement_index:
                    if k < nb:
                        cols.append(np.where(bins == k, f["P"], 0.0))
                    else:
                        cols.append(f["E"] * (enhancement_scale * Bft[:, k - nb]))
                JbNG, JbRG = activity_columns(np.column_stack(cols))
                for k, j in enumerate(enhancement_index):
                    J[("b", j)] = (JbNG[:, k], JbRG[:, k])
            du = -0.5 * u / Te
            dl, dn_ = np.where(use_lit, du, 0.0), np.where(use_lit, 0.0, du)
            dlx = dl * self.xlit

            def sums(w):
                s_l, s_lx, s_n = bc(w * dl), bc(w * dlx), bc(w * dn_)
                return KT * s_lx + s_l, -KT * s_lx, s_n
            gN = cNG["H2O"] * bm["H2O"] * self.wn["H2O"] * af + kappa * cNG["CO2"] * bm["CO2"] * self.wn["CO2"]
            tN = sums(gN)
            tRH = sums(bm["H2O"] * self.wv["H2O"] * af)
            tRC = sums(bm["CO2"] * self.wv["CO2"])
            for k in range(3):
                J[("T", k)] = (tN[k], cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k])
            J[("kappa", 0)] = (pNG_C, pRG_C)
            return pNG, pRG, J

    NonnegCpuModel.nonneg_model = cfg["model"]
    return NonnegCpuModel


# ------------------------------------------------------------------------------------------------ shared memory
def _shm_names():
    raw = os.environ.get(ENV_SHM)
    return json.loads(raw) if raw else None


def _shm_arrays(nfaces, create):
    """(D, E, P) global float64 arrays in shared memory (parent: create on first use; worker: attach)."""
    shm = _STATE["shm"]
    if shm.get("arrays") is not None:
        return shm["arrays"]
    from multiprocessing import shared_memory
    names = _shm_names()
    if names is None:
        raise RuntimeError("nonneg_mult: %s not set (the parent's install_gpu sets it)" % ENV_SHM)
    segs, arrs = [], []
    for nm_ in (names["D"], names["E"], names["P"]):
        if create:
            seg = shared_memory.SharedMemory(name=nm_, create=True, size=int(nfaces) * 8)
        else:
            seg = shared_memory.SharedMemory(name=nm_)
        segs.append(seg)
        arrs.append(np.ndarray((int(nfaces),), dtype=np.float64, buffer=seg.buf))
    shm.update(segments=segs, arrays=tuple(arrs), owner=bool(create))
    return shm["arrays"]


def _shm_close(unlink=False):
    shm = _STATE["shm"]
    segs = shm.get("segments") or []
    shm["arrays"] = None
    for seg in segs:
        try:
            seg.close()
            if unlink:
                seg.unlink()
        except Exception:
            pass
    shm["segments"] = []


# ------------------------------------------------------------------------------------------------ GPU patch
_SELECT = {}


def _select_kernel(cupy):
    if "kernel" not in _SELECT:
        _SELECT["kernel"] = cupy.ElementwiseKernel("int32 b, float64 P, int32 k", "float64 out",
                                                   "out = (b == k) ? P : 0.0", "nonneg_select_bin")
    return _SELECT["kernel"]


def install_gpu(cfg):
    """Patch the evaluator classes of THIS process (idempotent for the same constant cfg; refuses another)."""
    for k, v in os.environ.items():
        if (k.startswith("COPS_POS_") or k == sm.ENV_FROZEN) and v not in ("", None):
            raise RuntimeError("nonneg_mult: %s set (another model patch of the same methods); unset it" % k)
    if sm._STATE.get("installed"):
        raise RuntimeError("nonneg_mult: the slope_mult patch is installed in this process")
    text = canonical(cfg)
    if _STATE["installed"]:
        if canonical(_STATE["config"]) != text:
            raise RuntimeError("nonneg_mult: already installed with another configuration")
        return False
    frozen = os.environ.get(ENV_FROZEN)
    if frozen is not None and frozen != text:
        raise RuntimeError("nonneg_mult: this process's configuration %s differs from the inherited frozen one %s"
                           % (text, frozen))
    role = "worker" if frozen is not None else "driver"
    _STATE["role"] = role
    from cops_cuda import evaluator, multi_gpu, transport
    from cops_cuda import backend as bk
    from cops_cuda.errors import BackendError, IdentityError, UnsupportedConfiguration
    from cops_cuda.memory import MemoryPlan, scratch_reservation
    sm._check_sources(Path(evaluator.__file__).resolve().parent)
    B, MG, WS = evaluator.CudaBackend, multi_gpu.MultiGpuBackend, transport.WorkerState
    orig = _STATE["originals"]
    for name in ("__init__", "activity_linear", "_upload_activity", "_upload_columns", "identity", "release"):
        orig["B." + name] = getattr(B, name)
    for name in ("activity_linear", "_evaluate", "identity"):
        orig["MG." + name] = getattr(MG, name)
    orig["WS.handle"] = WS.handle
    orig["pipe_worker_main"] = transport.pipe_worker_main
    init_signature = inspect.signature(B.__init__)
    nb, nt = n_bins_templates(cfg)

    def __init__(self, *args, **kwargs):
        bound = init_signature.bind(self, *args, **kwargs)
        bound.apply_defaults()
        self._nonneg = None
        self._nonneg_state = None
        bound.arguments["scratch_per_face"] = int(bound.arguments["scratch_per_face"]) + BYTES_PER_FACE
        orig["B.__init__"](*bound.args, **bound.kwargs)

    def activity_linear(self, state):
        if state.b is None:
            return orig["B.activity_linear"](self, state)
        f = state_fields(self, cfg, state, False)
        return f["A"]

    def multi_activity_linear(self, state):
        if state.b is None:
            return orig["MG.activity_linear"](self, state)
        pend = getattr(self, "_nonneg_pending", None)
        if pend is not None and pend[0] is state:
            self._nonneg_pending = None
            return pend[1]
        return state_fields(self, cfg, state, False)["A"]

    def multi_evaluate(self, state, request):
        if state.b is not None:
            if self.activity != "host":
                raise UnsupportedConfiguration("nonneg_mult: the activity is a host field; use --activity host")
            if self.transport != "local":
                raise UnsupportedConfiguration("nonneg_mult: the column factors go through local shared memory; "
                                               "transport %r is not supported" % self.transport)
            if request.keys():
                f = state_fields(self, cfg, state, True)
                D, E, P = _shm_arrays(len(f["A"]), create=True)
                D[:], E[:], P[:] = f["D"], f["E"], f["P"]
                self._nonneg_pending = (state, f["A"])
        try:
            return orig["MG._evaluate"](self, state, request)
        finally:
            self._nonneg_pending = None

    def _free(self):
        d = getattr(self, "_nonneg", None)
        if d is None:
            return
        self._nonneg = None
        res = d.pop("_reservation", None)
        d.clear()
        if res is not None:
            res.release()

    def _ensure(self):
        d = getattr(self, "_nonneg", None)
        if d is not None and d["basis_sha256"] == self.enhancement_basis_sha256:
            return d
        _free(self)
        cupy, ds = self.cupy, self.dataset
        if self.enhancement_basis is None:
            raise BackendError("nonneg_mult: no template basis bound (bins + templates are required)")
        bins, Bft = basis_parts(self, cfg)
        nf = int(ds.nfaces_device)
        bins_loc = np.ascontiguousarray(ds.gather_faces(bins).astype(np.int32))
        yt = np.ascontiguousarray(np.asarray(ds.gather_faces(self.Y), dtype=np.float64).T)
        at = np.ascontiguousarray((float(self.enhancement_scale) * np.asarray(ds.gather_faces(Bft), np.float64)).T)
        plan = MemoryPlan("nonneg_mult: per-facet arrays of the column bases")
        plan.add_array("nonneg.YT", "static", yt.shape, np.float64, "Y^T of the device facets")
        plan.add_array("nonneg.AT", "static", (max(at.shape[0], 1), nf), np.float64, "A_s T^T of the device facets")
        for name in ("D", "E", "P"):
            plan.add_array("nonneg." + name, "scratch", (nf,), np.float64)
        plan.add_array("nonneg.bins", "static", (nf,), np.int32)
        with cupy.cuda.Device(self.device_id):
            reservation = self.budget.reserve(plan)
            try:
                with self.budget.allocation_guard(plan):
                    d = {"YT": cupy.asarray(yt), "AT": cupy.asarray(at) if at.shape[0] else None,
                         "D": cupy.empty(nf, dtype=np.float64), "E": cupy.empty(nf, dtype=np.float64),
                         "P": cupy.empty(nf, dtype=np.float64), "bins": cupy.asarray(bins_loc)}
            except BaseException:
                reservation.release()
                raise
        d.update(_reservation=reservation, basis_sha256=self.enhancement_basis_sha256, bins_host=bins_loc,
                 checked=False, check=None, columns_built=0,
                 info={"device_faces": nf, "device_id": int(self.device_id), "host": socket.gethostname(),
                       "pid": os.getpid(), "bytes_reserved": int(plan.total_bytes), "role": role,
                       "process_blas": blas_info()})
        self._nonneg = d
        return d

    def _upload_activity(self, state):
        self._nonneg_state = state if state.b is not None else None
        return orig["B._upload_activity"](self, state)

    def _factors(self, state):
        """Global (D, E, P): a worker reads them from the parent's shared memory, a single process computes them."""
        if self.activity_source is not None:
            return _shm_arrays(int(self.dataset.inputs.nfaces), create=False)
        f = state_fields(self, cfg, state, True)
        return f["D"], f["E"], f["P"]

    def _upload_columns(self, request):
        state = getattr(self, "_nonneg_state", None)
        if state is None:
            return orig["B._upload_columns"](self, request)
        d = _ensure(self)
        cupy, ds = self.cupy, self.dataset
        nf = int(ds.nfaces_device)
        D, E, P = self._factors(state)
        Dl = np.ascontiguousarray(ds.gather_faces(D))
        El = np.ascontiguousarray(ds.gather_faces(E))
        Pl = np.ascontiguousarray(ds.gather_faces(P))
        cb = self.arrays["colbasis"]
        with cupy.cuda.Device(self.device_id):
            d["D"].set(Dl)
            d["E"].set(El)
            d["P"].set(Pl)
            k = 0
            for j in request.basis_index:
                cupy.multiply(d["D"], d["YT"][int(j)], out=cb[k * nf:(k + 1) * nf])
                k += 1
            for i in request.enhancement_index:
                i = int(i)
                col = cb[k * nf:(k + 1) * nf]
                if 0 <= i < nb:
                    _select_kernel(cupy)(d["bins"], d["P"], np.int32(i), col)
                elif nb <= i < nb + nt:
                    cupy.multiply(d["E"], d["AT"][i - nb], out=col)
                else:
                    raise BackendError("nonneg_mult: template column %r of %d" % (i, nb + nt))
                k += 1
        self._columns_cache = None
        d["columns_built"] += 1
        if not d["checked"] and k:                  # once per bound basis: device columns vs the host expressions
            d["checked"] = True
            host = cb[:k * nf].get().reshape(k, nf)
            worst, kk = 0.0, 0
            Yl = ds.gather_faces(self.Y) if request.basis_index else None
            bins_l = d["bins_host"]
            at_h = d["AT"].get() if d["AT"] is not None else None
            for j in request.basis_index:
                worst = max(worst, float(np.max(np.abs(host[kk] - Dl * np.asarray(Yl[:, int(j)], np.float64)))))
                kk += 1
            for i in request.enhancement_index:
                i = int(i)
                want = np.where(bins_l == i, Pl, 0.0) if i < nb else El * at_h[i - nb]
                worst = max(worst, float(np.max(np.abs(host[kk] - want))))
                kk += 1
            d["check"] = {"columns": k, "max_abs_diff": worst}
            if worst != 0.0:
                raise IdentityError("nonneg_mult: device column basis differs from the host expressions (%r)" % worst)

    # ---- fine-tuning support (nonneg_finetune.py): dp/dA sensitivities, structural seen masks, detach
    def sensitivities(self, state, facets):
        from cops_cuda import evaluator as ev
        ds, cupy, c = self.dataset, self.cupy, self.constants
        facets = np.asarray(facets, np.int64)
        nf = int(ds.nfaces_device)
        if ds.local_faces is None:
            pos, ok = facets, np.ones(len(facets), bool)
        else:
            pos = np.minimum(np.searchsorted(ds.local_faces, facets), len(ds.local_faces) - 1)
            ok = ds.local_faces[pos] == facets
        npt = len(self.rows)
        ng, rg = np.zeros((npt, len(facets))), np.zeros((npt, len(facets)))
        saved_state = getattr(self, "_nonneg_state", None)
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                self.budget, "nonneg sensitivities", {"active_saved": (nf, self.arrays["active"].dtype)}):
            orig["B._upload_activity"](self, state)
            act = self.arrays["active"]
            saved = act.copy()
            cb = self.arrays["colbasis"]
            for c0 in range(0, len(facets), ev.FUSED_MAX_COLUMNS):
                cols = list(range(c0, min(c0 + ev.FUSED_MAX_COLUMNS, len(facets))))
                cb[:len(cols) * nf] = 0.0
                for k, cidx in enumerate(cols):
                    if ok[cidx]:
                        cb[k * nf + int(pos[cidx])] = 1.0
                        act[int(pos[cidx])] = 1
                codes = list(ev.PRESSURE_CODES) + [1000 + k for k in range(len(cols))] + \
                    [2000 + k for k in range(len(cols))]
                S = dict(zip(codes, self._run_sums(codes, self._params(state))))
                for k, cidx in enumerate(cols):
                    ng[:, cidx] = c.cNG["H2O"] * c.bm["H2O"] * S[1000 + k]
                    rg[:, cidx] = c.cRG["H2O"] * S[2000 + k]
                act[...] = saved
            del saved
            self._columns_cache = None
        self._nonneg_state = saved_state
        return ng, rg

    def seen_rows(self, row_mask, structural=False):
        ds, cupy = self.dataset, self.cupy
        if ds.streamed:
            raise BackendError("nonneg_mult: seen_rows needs resident records")
        if ds.record_precision != "fp32":
            raise BackendError("nonneg_mult: seen_rows is written for fp32 records (--precision mixed)")
        mask = np.asarray(row_mask, bool)
        if mask.shape != (len(self.rows),):
            raise BackendError("nonneg_mult: row mask of %s for %d rows" % (mask.shape, len(self.rows)))
        nblocks = int(self._first_block[-1])
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                self.budget, "nonneg seen rows", {"seen": (int(ds.nfaces_device), np.uint8), "mask": (len(mask), np.uint8)}):
            d_mask = cupy.asarray(mask.astype(np.uint8))
            d_seen = cupy.zeros(int(ds.nfaces_device), dtype=cupy.uint8)
            if nblocks:
                a = ds.arrays
                sm._seen_kernel(cupy)((min(nblocks, 1 << 20),), (256,),
                                      (self.arrays["blocks"], np.int64(nblocks), a["kernel.face"], a["kernel.kn"],
                                       a["kernel.kphi"], a["kernel.kv"], d_mask, d_seen,
                                       np.int32(1 if structural else 0)))
            host = d_seen.get().astype(bool)
            del d_mask, d_seen
        return ds.scatter_faces(host, False)

    def detach(self):
        _free(self)
        _shm_close(unlink=False) if _STATE["role"] == "worker" else None
        return True

    def multi_seen_rows(self, row_mask, structural=False):
        row_mask = np.asarray(row_mask, bool)
        out = None
        for w in self.workers:
            self._send(w, ("nonneg_seen", np.ascontiguousarray(row_mask[w["positions"]]), bool(structural)),
                       "nonneg_seen")
        for w in self.workers:
            m = self._receive(w, "nonneg_seen")[1]
            out = m if out is None else (out | m)
        return out

    def multi_sensitivities(self, state, facets):
        n = len(self.rows)
        ng, rg = np.zeros((n, len(facets))), np.zeros((n, len(facets)))
        for w in self.workers:
            self._send(w, ("nonneg_sens", state, np.asarray(facets, np.int64)), "nonneg_sens")
        for w in self.workers:
            _, a_ng, a_rg = self._receive(w, "nonneg_sens")
            ng[w["positions"]], rg[w["positions"]] = a_ng, a_rg
        return ng, rg

    def multi_detach(self):
        for w in self.workers:
            self._send(w, ("nonneg_detach",), "nonneg_detach")
        return [self._receive(w, "nonneg_detach")[1] for w in self.workers]

    def ws_handle(self, message):
        kind = message[0]
        if kind in WORKER_REQUESTS:
            try:
                fn = getattr(self.backend, WORKER_REQUESTS[kind])
                if kind == "nonneg_seen":
                    return ("ok", fn(message[1], message[2]))
                if kind == "nonneg_sens":
                    ng, rg = fn(message[1], message[2])
                    return ("ok", ng, rg)
                return ("ok", fn())
            except Exception:
                import traceback as _tb
                return ("error", _tb.format_exc())
        return orig["WS.handle"](self, message)

    def identity(self):
        value = orig["B.identity"](self)
        d = getattr(self, "_nonneg", None)
        info = None if d is None else dict(d["info"], colbasis_check=d["check"], columns_built=d["columns_built"])
        value["model_extension"] = {"nonneg": {"config": installed_config(), "backend": info,
                                               "module": str(Path(__file__).resolve()),
                                               "module_sha256": module_sha256(), "role": _STATE["role"],
                                               "held": stage_held()}}
        return value

    def release(self):
        try:
            orig["B.release"](self)
        finally:
            _free(self)

    def multi_identity(self):
        value = orig["MG.identity"](self)
        value["model_extension"] = {"nonneg": {"config": installed_config(), "module_sha256": module_sha256(),
                                               "driver_process_blas": blas_info(), "threads": thread_info(),
                                               "held": stage_held(), "timing": dict(_STATE["timing"])}}
        return value

    B.__init__ = __init__
    B.activity_linear, B._upload_activity, B._upload_columns = activity_linear, _upload_activity, _upload_columns
    B.identity, B.release, B.close = identity, release, release
    B._factors = _factors
    B.nonneg_sensitivities, B.nonneg_seen_rows, B.nonneg_detach = sensitivities, seen_rows, detach
    MG.activity_linear, MG._evaluate, MG.identity = multi_activity_linear, multi_evaluate, multi_identity
    MG.nonneg_seen_rows, MG.nonneg_sensitivities, MG.nonneg_detach = multi_seen_rows, multi_sensitivities, multi_detach
    WS.handle = ws_handle
    transport.pipe_worker_main = _pipe_worker_main
    os.environ[ENV_FROZEN] = text
    if role == "driver":
        tag = "nonneg_%d_%s" % (os.getpid(), hashlib.sha256(text.encode()).hexdigest()[:8])
        os.environ[ENV_SHM] = json.dumps({"D": tag + "_D", "E": tag + "_E", "P": tag + "_P"})
        atexit.register(_shm_close, True)
    _STATE["config"] = dict(cfg)
    _STATE["installed"] = True
    say("GPU patch installed (pid %d, %s, role %s): %s" % (os.getpid(), socket.gethostname(), role, text))
    return True


def install_worker_handlers():
    """Worker processes (local spawned workers; the one-job chain's worker() hook): install the patch from the
    parent's frozen configuration (idempotent)."""
    frozen = os.environ.get(ENV_FROZEN)
    if frozen is None:
        raise RuntimeError("nonneg_mult: worker without %s" % ENV_FROZEN)
    return install_gpu(json.loads(frozen))


def _pipe_worker_main(conn, bind_kwargs, device_id, shm_name, nfaces):
    """Local spawned worker: install the patch from the parent's frozen configuration, then the original loop."""
    install_worker_handlers()
    return _STATE["originals"]["pipe_worker_main"](conn, bind_kwargs, device_id, shm_name, nfaces)
