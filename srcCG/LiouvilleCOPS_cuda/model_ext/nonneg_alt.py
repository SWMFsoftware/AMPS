# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_alt.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_alt.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 74d4c5568e2a68642e14d054341d7d11055f9392d66d1bb63ef8e82223ce97cd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Alternative non-negative activity models (options A, B, D; user 2026-11-14): sos A = m_b
#   (||L^T y||^2 / A_s + A_s sum T) exp(h), background A = m_b (A_bg A_orig + softplus) exp(h), logspread A = m_b
#   (A_orig exp(u) + A_s sum T) exp(h) with the spread penalty f C_ref Var_w(u); extra parameter block, starts gram /
#   d2 / uniform, penalty rows.
# USED BY: nonneg_problem.install_outer (nonneg_fit.py --model sos|background|logspread); snapshots
#   model_ext_frozen/nonneg_20261114T084032Z (stage-1 B, A), nonneg_20261114T085147Z (stage-1 D).
# 2026-11-14: run copy updated (--spread-u0 (fixed-level spread penalty, 26 rows), --start restart dispatch; snapshot nonneg_20261115T020314Z)
# ==== end of header ====
"""Alternative non-negative activity models (options A, B, D of the lead, user 2026-11-14: "prepare suggestions for
modification of the minimization polynomial: the current approach shuts down almost all the surface, and I need a
polynomial that is non-negative and would allow fit with more seen surface contributing"): the fit-level part that
nonneg_problem.install_outer uses for --model sos | background | logspread (softplus and squared are unchanged and do
not use this module). Model math and the GPU patch: nonneg_mult.py (fields, host_fields, the column modes).

    sos (A):        A = m_b (||L^T y_e||^2 / A_s + A_s sum_k T_k) exp(h),  L 25 x r lower-trapezoidal; column 0 of L = the
                    fitter's a00..a24; the other columns are this module's extra parameters sosL_<i>_<j> (i >= j >= 1;
                    optimizer coordinate L_ij / ascale, unbounded); r = --sos-rank (25 + 24 + 23 = 72 parameters at r = 3)
    background (B): A = m_b (A_bg A_orig + s ln(1 + exp(L / s))) exp(h), L = Y a + A_s sum_k T_k, s = 1e-5 A_orig (the
                    softplus); A_bg one extra parameter (dimensionless, relative to the srcCG map A_orig), lower bound
                    --bg-min V, start V
    logspread (D):  A = m_b (A_orig exp(u) + A_s sum_k T_k) exp(h), u = Y a / A_s (= Y b, b = a / A_s dimensionless as in
                    the squared model, so that the fitter's coordinates a / ascale are O(1)); penalty (user 2026-11-14)
                    P = f C_ref Var_w(u), Var_w(u) = sum_seen w_e (u_e - u_bar)^2 (w_e = area / seen area, u_bar the
                    weighted mean): 25 rows sqrt(2 f C_ref) S b with S^T S = C = sum_seen w_e (y_e - y_bar)(y_e - y_bar)^T,
                    so that 0.5 ||rows||^2 = P exactly, on the scale of the data cost 0.5 r.r (f = --spread-f, C_ref =
                    --spread-cref); P is not part of the data cost and is reported apart (spread_report)
                    start uniform: u = ln s (a00 only), s from ln p least squares against D2's modelled pressures
Templates (D2's ft_c1..ft_c6, held at 1): all >= 0 on every facet (min 0.0, the fits' model_check), so every model is
non-negative: sos and logspread add A_s sum T outside, background has them inside the softplus (as softplus).
Layout: the extra parameters are inserted just before the two grain parameters (the grain stays last): names = [a00..a24,
c_slope_* (free), kappa, <extra>, grain_A_uPa, grain_B_uPa_per_100s]. The base problem (fit_cops / grain_crossing /
fixed_parameters / nonneg_problem.NonnegOuter) sees the vector without the extra block; this class sets the extra
parameters in nonneg_mult (host fields) before every evaluation and appends their Jacobian columns from additional
evaluations with the column mode set (sos: the factor of column j of L; background: "raw" dA/dA_bg = m A_orig exp(h)).
"""

import json
import math
from pathlib import Path

import numpy as np

import nonneg_mult as nm

SPREAD_WEIGHTS = "area / seen area (normalised: Var_w(u) is the area-weighted variance of u over the seen facets)"


def extra_names(model_cfg):
    m = model_cfg["model"]
    if m == "sos":
        r = int(model_cfg["sos_rank"])
        return ["sosL_%02d_%d" % (i, j) for j in range(1, r) for i in range(j, 25)]
    if m == "background":
        return ["A_bg"]
    return []


def extra_from_values(model_cfg, values, ascale):
    """{"L_extra": (25, r - 1)} (sos, physical = ascale x) / {"A_bg": v} (background) / {} from the extra block."""
    m = model_cfg["model"]
    v = np.asarray(values, dtype=np.float64)
    if m == "sos":
        r = int(model_cfg["sos_rank"])
        L = np.zeros((25, r - 1))
        k = 0
        for j in range(1, r):
            for i in range(j, 25):
                L[i, j - 1] = v[k] * ascale
                k += 1
        return {"L_extra": L}
    if m == "background":
        return {"A_bg": float(v[0])}
    return {}


def extra_from_restart(model_cfg, names, x, ascale):
    en = extra_names(model_cfg)
    idx = {n: k for k, n in enumerate(names)}
    missing = [n for n in en if n not in idx]
    if missing:
        raise SystemExit("nonneg_alt: the restart lacks %s" % missing[:4])
    return extra_from_values(model_cfg, [x[idx[n]] for n in en], ascale)


def spread_matrix(Y, area, seen, chunk=1 << 18):
    """(S (25 x 25) with S^T S = C, C, ybar, seen area): C = sum_seen w (y - ybar)(y - ybar)^T, w = area / seen area."""
    seen = np.asarray(seen, bool)
    w_all = np.where(seen, np.asarray(area, float), 0.0)
    W = float(w_all.sum())
    w_all = w_all / W
    G = np.zeros((25, 25))
    m1 = np.zeros(25)
    for c0 in range(0, len(w_all), chunk):
        Yc = np.asarray(Y[c0:c0 + chunk], float)
        wc = w_all[c0:c0 + chunk]
        G += Yc.T @ (Yc * wc[:, None])
        m1 += Yc.T @ wc
    C = G - np.outer(m1, m1)
    C = 0.5 * (C + C.T)
    lam, U = np.linalg.eigh(C)
    lam = np.maximum(lam, 0.0)
    S = (np.sqrt(lam)[:, None] * U.T)
    return S, C, m1, W


def make_class(Base, fit_cops, cfg, model_cfg, record, held, npb):
    """The Problem class of an alternative model (subclass of nonneg_problem's NonnegOuter)."""
    model = model_cfg["model"]
    names_x = extra_names(model_cfg)
    ne = len(names_x)
    SF = float(cfg.get("spread_f") or 0.0)
    SC = float(cfg.get("spread_cref") or 0.0)
    X = 2.0 * SF * SC                     # rows sqrt(X) S b: 0.5 ||rows||^2 = f C_ref Var_w(u)
    U0 = cfg.get("spread_u0")             # user 2026-11-14: spread about the fixed level u0 (one more row)
    bg_min = cfg.get("bg_min")

    class AltOuter(Base):
        _nonneg_defer_start = True

        def __init__(self, M, el, coef, args):
            # the base constructors evaluate the model (grain_crossing's default grain start): the extra block at its
            # default values meanwhile (background: A_bg = bg_min; sos: L columns 1.. = 0)
            nm.set_extra(extra_from_values(model_cfg, [float(bg_min)] if model == "background" else np.zeros(ne), 1.0))
            super().__init__(M, el, coef, args)
            xb_start = None
            if cfg.get("start_restart") and cfg["start_choice"] == "d2":     # background: the base start first
                if model != "background":
                    raise SystemExit("nonneg_alt: --start d2 is not a start of the %s model" % model)
                xb_start = npb.start_d2(self, el, cfg["start_restart"], record, cfg)
            names = list(self.names)
            self._alt_nbase = len(names)
            k0 = len(names) - 2
            if tuple(names[k0:]) != tuple(npb.GRAIN_NAMES):
                raise SystemExit("nonneg_alt: the grain is not last in %s" % names[-3:])
            self._alt_k0, self._alt_ne = k0, ne
            x0e = np.zeros(ne)
            loe = np.full(ne, -np.inf)
            if model == "background":
                x0e[:] = float(bg_min)
                loe[:] = float(bg_min)
            self.names = names[:k0] + names_x + names[k0:]
            self.x0 = np.concatenate([self.x0[:k0], x0e, self.x0[k0:]])
            self.lo = np.concatenate([np.asarray(self.lo, float)[:k0], loe, np.asarray(self.lo, float)[k0:]])
            self.hi = np.concatenate([np.asarray(self.hi, float)[:k0], np.full(ne, np.inf),
                                      np.asarray(self.hi, float)[k0:]])
            self._alt_n = len(self.names)
            rec = {"model": model, "extra_parameters": names_x[:6] + (["..."] if ne > 6 else []), "n_extra": ne,
                   "layout": "a00..a24, c_slope_* (free), kappa, <extra>, grain_A_uPa, grain_B_uPa_per_100s",
                   "n_parameters": self._alt_n}
            if model == "background":
                rec["bg_min"] = float(bg_min)
            if model == "sos":
                rec["sos_rank"] = int(model_cfg["sos_rank"])
            self._alt_pen = None
            self._alt_off = None
            if model == "logspread":
                seen = M.seen_faces()
                S, C, ybar, W = spread_matrix(el["Y"], el["area"], seen)
                self._alt_S = (self.ascale / float(self.enhancement_scale)) * S     # ||S_x x[a]||^2 = Var_w(u)
                self._alt_ybar = (self.ascale / float(self.enhancement_scale)) * np.asarray(ybar, float)  # u_bar = . x[a]
                self._alt_pen = math.sqrt(X) * self._alt_S
                self._alt_off = None
                if U0 is not None:            # sum w (u - u0)^2 = Var_w(u) + (u_bar - u0)^2: row sqrt(X) (u_bar - u0)
                    self._alt_pen = np.vstack([self._alt_pen, math.sqrt(X) * self._alt_ybar[None, :]])
                    self._alt_off = np.zeros(26)
                    self._alt_off[25] = math.sqrt(X) * float(U0)
                rec["spread"] = {"f": SF, "C_ref": SC, "rows_factor_sqrt_2_f_C_ref": math.sqrt(X), "u0": U0,
                                 "centre": ("fixed level u0 (user 2026-11-14): P = f C_ref [Var_w(u) + (u_bar - u0)^2], "
                                            "26 rows" if U0 is not None else "weighted mean (Var_w(u)), 25 rows"),
                                 "weights": SPREAD_WEIGHTS, "seen_facets": int(np.sum(seen)),
                                 "seen_area_m2": W, "rows": 25, "C_trace": float(np.trace(C)),
                                 "C_eigenvalues": [float(v) for v in np.linalg.eigvalsh(C)],
                                 "definition": "user 2026-11-14: penalty P = f C_ref Var_w(u), Var_w(u) = sum_seen w (u - "
                                               "u_bar)^2, w = area / seen area, u = Y a / A_s; in the fitter's cost 0.5 "
                                               "r.r as 25 rows sqrt(2 f C_ref) (ascale / A_s) S x[a] (S^T S = C, x[a] = "
                                               "a / ascale), 0.5 ||rows||^2 = P exactly; not in the data cost"}
            record["model_check"]["alternative"] = rec
            npb.say("alternative model %s: %d extra parameters %s; %s" % (
                model, ne, names_x[:3] + (["..."] if ne > 3 else []),
                ("spread penalty f C_ref Var_w(u), f %r, C_ref %r, over %d seen facets (25 rows, C trace %.4g)" % (
                    SF, SC, rec["spread"]["seen_facets"], rec["spread"]["C_trace"]) if model == "logspread" else
                 "A_bg >= %r (start %r)" % (bg_min, bg_min) if model == "background" else "")))
            if cfg.get("start_restart"):
                ch = cfg["start_choice"]
                if ch == "chain":
                    x = npb.start_chain(self, el, cfg["start_restart"], record, cfg["grain_from"], model_cfg, held,
                                        cfg.get("slope_mult_min"))
                elif ch == "d2":                # background: D2's a, kappa, grain, multipliers 1; A_bg = bg_min
                    x = self._full(xb_start, x0e)
                    record["start_restart"]["A_bg_start"] = float(bg_min)
                elif ch == "gram":              # sos
                    x = start_gram(self, el, cfg, record, npb)
                elif ch == "uniform":           # logspread
                    x = start_uniform(self, el, cfg, record, npb, fit_cops)
                elif ch == "restart":           # continuation of a complete fit of the same model (user 2026-11-14)
                    x = npb.start_restart(self, el, cfg["start_restart"], record, model_cfg, cfg)
                else:
                    raise SystemExit("nonneg_alt: start %r is not a start of the %s model" % (ch, model))
                xs = np.asarray(x, float)
                bad = [n for k, n in enumerate(self.names) if xs[k] < self.lo[k] or xs[k] > self.hi[k]]
                if bad:
                    raise SystemExit("nonneg_alt: the start of %s is outside the bounds" % bad)
                self.restart_start = xs
                record["start_restart"]["x_start"] = [float(v) for v in xs]     # the full layout
                record["start_restart"]["names"] = list(self.names)
                if model == "logspread":
                    pen = self._pen(self._split(xs)[0])
                    sv = self._alt_S @ np.asarray(self._split(xs)[0], float)[:25]
                    mu0 = float(self._alt_ybar @ np.asarray(self._split(xs)[0], float)[:25])
                    spr = record["model_check"]["alternative"]["spread"]
                    spr["penalty_cost_at_start"] = 0.5 * float(pen @ pen)
                    spr["var_w_u_at_start"] = float(sv @ sv)
                    spr["mean_u_at_start"] = mu0
                    if U0 is not None:
                        spr["dev2_mean_u0_at_start"] = (mu0 - float(U0)) ** 2

        # --- the layout
        def _full(self, xbase, xextra):
            k0 = self._alt_k0
            xbase = np.asarray(xbase, float)
            return np.concatenate([xbase[:k0], np.asarray(xextra, float), xbase[k0:]])

        def _split(self, x):
            x = np.asarray(x, dtype=float)
            if not hasattr(self, "_alt_n"):           # during the base constructors: the base layout
                return x, None
            if len(x) == self._alt_n:
                k0, ne_ = self._alt_k0, self._alt_ne
                return np.concatenate([x[:k0], x[k0 + ne_:]]), x[k0:k0 + ne_]
            if len(x) == self._alt_nbase:
                return x, None
            raise ValueError("nonneg_alt: parameter vector of length %d (full %d, base %d)" % (
                len(x), self._alt_n, self._alt_nbase))

        def _set(self, x):
            base, ext = self._split(x)
            if ext is not None:
                nm.set_extra(extra_from_values(model_cfg, ext, self.ascale))
            return base

        def _pen(self, base):
            r = self._alt_pen @ np.asarray(base, float)[:25]
            return r if self._alt_off is None else r - self._alt_off

        def spread_report(self, x):
            """logspread: the data cost and the spread penalty P = f C_ref Var_w(u) at x, apart."""
            base = self._set(x)
            sv = self._alt_S @ np.asarray(base, float)[:25]
            pen = self._pen(base)
            rd = np.asarray(self.data_residuals(x), float)
            var = float(sv @ sv)
            mu = float(self._alt_ybar @ np.asarray(base, float)[:25])
            dev2 = (mu - float(U0)) ** 2 if U0 is not None else 0.0
            out = {"f": SF, "C_ref": SC, "var_w_u": var, "rms_spread": math.sqrt(var), "mean_u": mu, "u0": U0,
                   "dev2_mean_u0": dev2 if U0 is not None else None,
                   "penalty": 0.5 * float(pen @ pen), "penalty_formula": SF * SC * (var + dev2),
                   "penalty_variance_part": SF * SC * var, "penalty_level_part": SF * SC * dev2,
                   "data_cost": 0.5 * float(rd @ rd)}
            out["data_plus_penalty"] = out["data_cost"] + out["penalty"]
            return out

        # --- the base problem's methods on the base vector
        def unpack(self, x):
            return super().unpack(self._set(x))

        def data_residuals(self, x):
            return super().data_residuals(self._set(x))

        def residuals(self, x):
            base = self._set(x)
            r = super().residuals(base)
            if getattr(self, "_alt_pen", None) is not None:
                r = np.concatenate([r, self._pen(base)])
            return r

        def regularization(self, x):
            return super().regularization(self._split(x)[0])

        def reg_jacobian(self):
            R = super().reg_jacobian()
            if not hasattr(self, "_alt_n"):
                return R
            k0 = self._alt_k0
            return np.hstack([R[:, :k0], np.zeros((R.shape[0], self._alt_ne)), R[:, k0:]])

        def deviation(self, x):
            return super().deviation(self._split(x)[0])

        def activity_linear_and_derivative(self, x, parameter_index):
            base = self._set(x)
            if not hasattr(self, "_alt_n"):
                return super().activity_linear_and_derivative(base, parameter_index)
            k0, ne_ = self._alt_k0, self._alt_ne
            if k0 <= parameter_index < k0 + ne_:
                total, _ = super().activity_linear_and_derivative(base, 0)
                return total, None
            k = parameter_index - ne_ if parameter_index >= k0 + ne_ else parameter_index
            return super().activity_linear_and_derivative(base, k)

        def _extra_columns(self, base):
            """Jacobian columns of the extra block (residual units), from evaluations with the column mode set."""
            if not self._alt_ne:
                return np.empty((self.ndata, 0))
            a, b, T, kappa = super().unpack(base)
            g, _, _ = self.grain(base, self.et_fit)
            uN, uR = self.useNG, self.useRG
            cols = []
            if model == "sos":
                jobs = [(j, list(range(j, 25)), self.ascale) for j in range(1, int(model_cfg["sos_rank"]))]
            else:
                jobs = [("raw", [0], 1.0)]
            try:
                for mode, bidx, scale in jobs:
                    nm.set_col_mode(mode)
                    pNG, pRG, J = self.M.evaluate(a, T, kappa, b=b, enhancement_basis=self.enhancement_basis,
                                                  enhancement_scale=self.enhancement_scale, jac=True,
                                                  basis_index=bidx, enhancement_index=[])
                    pNG_ = np.maximum(pNG + g, fit_cops.FLOOR)
                    pRG_ = np.maximum(pRG, fit_cops.FLOOR)
                    for i in bidx:
                        dn, dp = J[("a", i)]
                        cols.append(np.concatenate([self.wNG * dn[uN] / pNG_[uN], self.wRG * dp[uR] / pRG_[uR]]) *
                                    scale)
            finally:
                nm.set_col_mode(0)
            return np.column_stack(cols)

        def data_jacobian(self, x):
            base = self._set(x)
            nm.set_col_mode(0)
            J = super().data_jacobian(base)
            if not hasattr(self, "_alt_n"):
                return J
            Je = self._extra_columns(base)
            k0 = self._alt_k0
            return np.hstack([J[:, :k0], Je, J[:, k0:]])

        def jacobian(self, x):
            Jd = self.data_jacobian(x)
            if self.lam > 0:
                Jd = np.vstack([Jd, self.reg_jacobian()])
            if getattr(self, "_alt_pen", None) is not None:
                Jp = np.zeros((self._alt_pen.shape[0], self._alt_n))
                Jp[:, :25] = self._alt_pen
                Jd = np.vstack([Jd, Jp])
            return Jd

    return AltOuter


def _d2_state(P, el, path, npb):
    from cops_cuda import restart as rs
    state = rs.read(path)
    npb._check_source_state(P, el, state, path)
    return state


def start_gram(P, el, cfg, record, npb):
    """--start gram (sos): L from --start-file (nonneg_start.py --model sos: the area-weighted LSQ Gram matrix Q of
    D2's activity, PSD-projected, top-r eigenvectors -> lower-trapezoidal L); multipliers 1; kappa, grain from D2."""
    state = _d2_state(P, el, cfg["start_restart"], npb)
    with np.load(cfg["start_file"], allow_pickle=False) as z:
        L = np.asarray(z["L"], float)
        smeta = json.loads(str(z["meta_json"]))
    r = int(cfg["sos_rank"])
    problems = []
    if L.shape != (25, r):
        problems.append("L has shape %s, the model rank %d" % (L.shape, r))
    if smeta.get("d2_restart", {}).get("sha256") != npb.sha256(cfg["start_restart"]):
        problems.append("the start file was computed from another D2 restart")
    if float(smeta.get("A_scale")) != float(P.enhancement_scale):
        problems.append("A_scale %r vs %r" % (smeta.get("A_scale"), P.enhancement_scale))
    if smeta.get("mesh_identity") != state.mesh_identity():
        problems.append("mesh differs")
    if np.any(np.triu(L, 1) != 0.0):
        problems.append("L is not lower-trapezoidal")
    if problems:
        raise SystemExit("nonneg_alt: --start-file %s: %s" % (cfg["start_file"], "; ".join(problems)))
    names_r = list(state.meta.get("names") or [])
    x = np.array(P.x0, dtype=float, copy=True)
    for k, name in enumerate(P.names):
        if name.startswith("a") and name[1:].isdigit():
            x[k] = L[int(name[1:]), 0] / P.ascale
        elif name.startswith("sosL_"):
            i, j = int(name[5:7]), int(name[8:])
            x[k] = L[i, j] / P.ascale
        elif name.startswith("c_slope_"):
            x[k] = npb.START_MULTIPLIER
        elif name == "kappa" or name in npb.GRAIN_NAMES:
            x[k] = float(state.x[names_r.index(name)])
        else:
            raise SystemExit("nonneg_alt: --start gram: unexpected parameter %s" % name)
    gt, chk = npb._grain_check(x, P, state)
    record["start_restart"] = {"path": str(Path(cfg["start_restart"]).resolve()), "sha256": npb.sha256(cfg["start_restart"]),
                               "start_choice": "gram", "start_definition": npb.START_CHOICES["gram"] + (
                                   "; then refined by Levenberg-Marquardt of the rank-r problem (nonneg_start.py "
                                   "--refine lm)" if smeta.get("refine") == "lm" else ""),
                               "start_file": {"path": cfg["start_file"], "sha256": npb.sha256(cfg["start_file"])},
                               "start_quality": smeta.get("quality"), "names_from_restart": ["kappa"] + list(npb.GRAIN_NAMES),
                               "multipliers_set_to": npb.START_MULTIPLIER, "grain_term_check": chk,
                               "x_start": [float(v) for v in x]}
    npb.say("start gram: L (25 x %d) from %s (quality %s), kappa and grain by name from D2, multipliers 1" % (
        r, cfg["start_file"], {k: (smeta.get("quality") or {}).get(k) for k in ("relative_rms_residual",
                                                                              "eigenvalue_kept_fraction")}))
    return x


def start_uniform(P, el, cfg, record, npb, fit_cops):
    """--start uniform (logspread, user 2026-11-14): u = ln s on every facet (a00 only; every other harmonic 0), i.e.
    A = m (s A_orig + A_s sum T) exp(h); s = the least-squares fit in ln p of the model pressures to D2's modelled
    pressures over the fit's used NG and RG rows: minimise sum_rows (ln p(s) - ln p_D2)^2 (unweighted; the fit-weighted
    value is recorded as well). p_D2 = D2's stored predictions (its fit_restart.npz arrays = its COPS_model.dat),
    matched by UTC exactly as clip_problem.start_check does; p(s) = the model as the residuals see it (W10x5 for the
    late observations, + grain for NG). At fixed T, kappa and grain the pressures are linear in s, so p(s) = P0 + s P1
    from two evaluations (s = 1, 2); a third evaluation at the result checks that. Multipliers 1; kappa and grain by
    name from D2."""
    from scipy.optimize import minimize_scalar
    cp, pp = npb.cp, npb.pp
    state = _d2_state(P, el, cfg["start_restart"], npb)
    names_r = list(state.meta.get("names") or [])
    Y0 = np.asarray(el["Y"][:, 0], float)
    y00 = float(Y0[0])
    if not (y00 > 0 and float(np.max(np.abs(Y0 - y00))) <= 1e-12 * y00):
        raise SystemExit("nonneg_alt: --start uniform: Y[:, 0] is not a positive constant (%r .. %r)" % (Y0.min(), Y0.max()))
    x = np.array(P.x0, dtype=float, copy=True)
    for k, name in enumerate(P.names):
        if name.startswith("a") and name[1:].isdigit():
            x[k] = 0.0
        elif name.startswith("c_slope_"):
            x[k] = npb.START_MULTIPLIER
        elif name == "kappa" or name in npb.GRAIN_NAMES:
            x[k] = float(state.x[names_r.index(name)])
        else:
            raise SystemExit("nonneg_alt: --start uniform: unexpected parameter %s" % name)
    ia = P.names.index("a00")
    A_s = float(P.enhancement_scale)

    def x_of(t):
        x_ = x.copy()
        x_[ia] = A_s * t / y00 / P.ascale          # u = Y00 a00 / A_s = t on every facet
        return x_
    M = P.M
    info, arrays = cp.load_source(cfg["start_restart"])
    points = pp.read_run_points(M.run)
    ip = np.asarray(M.ipoint)[np.asarray(M.rows)]
    keys = np.array([pp.time_key(points[int(i)][0]) for i in ip], dtype=np.int64)
    src_keys = np.array([pp.time_key(str(u)) for u in arrays["utc"]], dtype=np.int64)
    pos = {int(k): j for j, k in enumerate(src_keys)}
    useNG, useRG = np.asarray(P.useNG, bool), np.asarray(P.useRG, bool)
    late_key = pp.time_key(pp.LATE_START_UTC)
    target, wts, early, gauge = [], [], [], []
    for g_, use, spred, w_ in (("NG", useNG, arrays["pred_ng"], P.wNG), ("RG", useRG, arrays["pred_rg"], P.wRG)):
        k_use = keys[use]
        missing = [int(k) for k in k_use if int(k) not in pos]
        if missing:
            raise SystemExit("nonneg_alt: --start uniform: %d used %s rows have no D2 row" % (len(missing), g_))
        j = np.array([pos[int(k)] for k in k_use], dtype=np.int64)
        target.append(np.log(np.maximum(np.asarray(spred, float)[j], fit_cops.FLOOR)))
        wts.append(np.broadcast_to(np.asarray(w_, float), (int(use.sum()),)))
        early.append(k_use < late_key)
        gauge.append(np.full(int(use.sum()), g_))
    target, wts, early, gauge = (np.concatenate(v) for v in (target, wts, early, gauge))

    def model_rows(x_):
        a, b, T, kappa = P.unpack(x_)
        pNG_s, pRG_s = M.evaluate(a, T, kappa, b=b if P.nb else None, enhancement_basis=P.enhancement_basis,
                                  enhancement_scale=P.enhancement_scale)
        g = P.grain(x_, P.et_fit)[0]
        return np.concatenate([np.asarray(pNG_s + g, float)[useNG], np.asarray(pRG_s, float)[useRG]])
    p1, p2 = model_rows(x_of(0.0)), model_rows(x_of(math.log(2.0)))
    P1, P0 = p2 - p1, 2.0 * p1 - p2
    if np.any(P1 <= 0) or np.any(P0 < -1e-9 * np.abs(p1)):
        raise SystemExit("nonneg_alt: --start uniform: p(s) is not P0 + s P1 with P1 > 0, P0 >= 0 (min P1 %r, min P0/p %r)"
                         % (float(P1.min()), float((P0 / p1).min())))
    P0 = np.maximum(P0, 0.0)

    def F(t, w=None):
        d = np.log(P0 + math.exp(t) * P1) - target
        return float(np.sum(d * d)) if w is None else float(np.sum((w * d) ** 2))

    def solve(w=None):
        grid = np.arange(-30.0, 30.0 + 1e-9, 0.25)
        vals = [F(t, w) for t in grid]
        t0 = float(grid[int(np.argmin(vals))])
        r_ = minimize_scalar(lambda t: F(t, w), bounds=(t0 - 0.25, t0 + 0.25), method="bounded",
                             options={"xatol": 1e-12, "maxiter": 500})
        return float(r_.x)
    t_star, t_w = solve(), solve(wts)
    xs = x_of(t_star)
    p3 = model_rows(xs)
    lin = float(np.max(np.abs(p3 / (P0 + math.exp(t_star) * P1) - 1.0)))
    if lin > 1e-9:
        raise SystemExit("nonneg_alt: --start uniform: the pressures are not linear in s (max rel %.2e)" % lin)
    d = np.log(p3) - target

    def st(m):
        return {"rows": int(m.sum()), "rms_dln": float(np.sqrt(np.mean(d[m] ** 2))) if m.any() else None,
                "mean_dln": float(np.mean(d[m])) if m.any() else None}
    uni = {"s": math.exp(t_star), "ln_s": t_star, "objective": "sum over the used NG and RG rows of (ln p(s) - ln p_D2)^2, "
           "unweighted (p_D2 = D2's stored predictions, matched by UTC as in the start check)",
           "sum_sq_dln": F(t_star), "sum_sq_dln_at_ln_s_minus_plus_0p01": [F(t_star - 0.01), F(t_star + 0.01)],
           "rows": int(len(d)),
           "dln_stats": {"NG": st(gauge == "NG"), "RG": st(gauge == "RG"), "before_late_window": st(early),
                         "late_window": st(~early), "all": st(np.ones(len(d), bool))},
           "fit_weighted_alternative": {"s": math.exp(t_w), "ln_s": t_w,
                                        "objective": "sum of (w (ln p(s) - ln p_D2))^2 with the fit's row weights "
                                                     "(recorded only)"},
           "linearity_check_max_rel": lin, "share_P0_min_max": [float((P0 / (P0 + math.exp(t_star) * P1)).min()),
                                                                 float((P0 / (P0 + math.exp(t_star) * P1)).max())],
           "a00": A_s * t_star / y00, "Y00": y00, "d2_source": {"path": info["path"], "sha256": info["sha256"]}}
    gt, chk = npb._grain_check(xs, P, state)
    record["start_restart"] = {"path": str(Path(cfg["start_restart"]).resolve()), "sha256": npb.sha256(cfg["start_restart"]),
                               "start_choice": "uniform", "start_definition": npb.START_CHOICES["uniform"],
                               "names_from_restart": ["kappa"] + list(npb.GRAIN_NAMES),
                               "multipliers_set_to": npb.START_MULTIPLIER, "grain_term_check": chk,
                               "uniform": uni, "x_start": [float(v) for v in xs]}
    record["model_check"].setdefault("alternative", {})["uniform_start"] = uni
    npb.say("start uniform: u = ln s = %.12f (s = %.10e; fit-weighted alternative s = %.6e) on every facet, rms d ln p vs "
            "D2: NG %.4f, RG %.4f, before 10:33:10 %.4f, late %.4f (linearity %.1e); multipliers 1, kappa and grain by "
            "name from D2" % (t_star, math.exp(t_star), math.exp(t_w), uni["dln_stats"]["NG"]["rms_dln"],
                              uni["dln_stats"]["RG"]["rms_dln"], uni["dln_stats"]["before_late_window"]["rms_dln"],
                              uni["dln_stats"]["late_window"]["rms_dln"] or float("nan"), lin))
    return xs
