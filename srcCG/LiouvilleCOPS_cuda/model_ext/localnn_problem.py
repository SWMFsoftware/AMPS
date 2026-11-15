# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f98d2a06d909ca0274b25c369911c8ac941fb70889516d585e58a82697213289
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module (option C). The fit problem on CPU: D2's selection, W10x5 (positive_problem.LateLayout), grain, crossing
#   rows (grain_crossing), smoothness rows, residuals, Jacobian, start projection (NNLS), lambda rule, table compare.
#   2026-11-14 ~17:30 CST (user): optional multiplier penalty P_m = f_m C_ref sum_b (ln m_b)^2 (mult_penalty_f; rows
#   sqrt(2 f_m C_ref) ln m_b; absent / 0 = the previous code bitwise).
# USED BY: localnn_fit.py, localnn_direct.py, localnn_report.py, localnn_degeneracy.py, tests (frozen localnn_anch_20261114T090437Z for
#   stage 1 without P_m; model_ext_frozen/localnn_mp_20261115T020233Z with P_m).
# ==== end of header ====
"""Option C (local non-negative basis, user 2026-11-14): the fit problem on CPU from the precomputed design matrices
(localnn_design.py) and the coarse grid (localnn_mesh.py). The driver is localnn_fit.py.

Local, untracked model extension: nothing in the git tree is edited and no existing model_ext file is changed;
positive_problem.py (LateLayout = the W10x5 operator of D2 and of every nonneg fit) is imported unchanged.

MODEL (per facet e, slope bin b(e), coarse hat functions phi_j, three per facet):
    A_e = m_b(e) * ( sum_j c_j phi_j(e) + TS_e ),   c_j >= 0, m_[0,40) = 1 held, the six others >= --slope-mult-min
    TS_e = A_s sum_k ft_ck(e) (D2's six fine-tuning templates held at 1) with --templates on, 0 with --templates off.
    The ft templates are >= 0 everywhere (checked when the design is built), so with c >= 0 the activity is >= 0
    without a wrapper; the hard clip max(0, .) around (sum c phi + TS) is applied in the per-facet products and has no
    effect (the number of clipped facets is written; it is 0).
Optimizer coordinates: x = (c / A_s [nodes], m_1..m_6, kappa, grain_A_uPa, grain_B_uPa_per_100s); A_s = D2's
enhancement scale 2.2228299497402522e18 (the area-weighted RMS of the srcCG linear field).

EVERYTHING ELSE = D2 (clip_problem.py variant D2 / grain_crossing.py / fit_cops.Problem, read from that code):
* surface pressures per base row (1551): p = sum_b m_b (G^(b) c + Tcol^(b)) + kappa Ccol (localnn_design.py; FP32 records,
  FP64 arithmetic as the GPU fits), T held at the original law (182.09922 / 136.11149 / 133 K);
* late observations (>= 2016-09-30T10:33:10 UTC, 95): W10x5 of the 1-s grid rows (positive_problem.LateLayout, same point
  table); earlier observations instantaneous; the 401 grid rows are model-only (use 0);
* NG = surface + grain(t), grain(t) = max(0, A - B (t - t_end)) for t >= t_end (unaveraged, at the time stamp), t_end =
  the first used NG point after the largest NG gap of the last hour (10:33:10 UTC); A, B >= 0; one-sided derivatives
  as grain_crossing.GrainProblem.grain;
* data residuals (fit_cops.Problem.data_residuals): w_NG (ln max(pNG, 1e-300) - ln refNG) on the used NG rows, w_RG (...)
  on the used RG rows (use > 0 and reference > 0: D2's selection, 934 NG, 1150 RG), w = 1/sqrt(N) per gauge; data
  cost = 0.5 |r|^2 (the reported cost of every fit);
* crossing rows (grain_crossing.least_squares_with): at the rows where both gauges are used, before the gap (t <=
  last NG before the gap) MU w max(0, ln pNG - ln pRG), from t_end on MU w max(0, ln pRG - ln pNG), w = sqrt(w_NG w_RG),
  MU = 10; exact Jacobian; excluded from the reported data cost;
* smoothness (option C, the lead's request): penalty lambda_C x_c^T Lg x_c (x_c = c / A_s, Lg = the P1 stiffness /
  cotangent Laplacian of the coarse mesh, localnn_mesh.py), as residual rows sqrt(2 lambda_C) R x_c with R^T R = Lg
  (R = diag(sqrt(mu)) V^T of the eigen-decomposition of Lg, eigenvalues <= 1e-12 max dropped); reported separately;
* optimizer: scipy least_squares TRF, bounds, x_scale 'jac', exact trust-region solver, ftol = xtol = gtol = 1e-10
  (D2's), max_nfev from the command line.
Objective = 0.5 |r_data|^2 + 0.5 |q_cross|^2 + lambda_C x_c^T Lg x_c [+ P_m with --mult-penalty-f f_m > 0 (user
2026-11-14 ~17:30 CST): P_m = f_m C_ref sum_b (ln m_b)^2 over the six free multipliers, rows sqrt(2 f_m C_ref) ln m_b].

START (--start d2, the lead 2026-11-14): c = the area-weighted non-negative least-squares projection of D2's activity
max(0, D2 field) onto the basis: min_c>=0 sum_e area_e (sum_j c_j phi_j(e) - target_e)^2 (exact facet mass matrix M =
Phi^T diag(area) Phi; scipy nnls on the Cholesky factor), target = max(0, D2 field) - TS (TS = 0 with --templates off),
so that the start activity at m = 1 approximates D2's activity; multipliers 1; kappa and grain A, B from D2.

LAMBDA RULE (proposed, not yet decided): lambda_C = f * cost_data(start) / (x_c0^T Lg x_c0), so that the penalty at the
start is f times the start data cost; f = 0.1 (weak), 1 (medium), 10 (strong).
"""

import datetime as dt
import hashlib
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import nnls

import positive_problem as pp      # LateLayout (W10x5), point-table readers (unchanged)

SCHEMA = "localnn_fit_v1"
PARAMETERIZATION = "localnn_hat_slope_multiplier_w10x5_v1"
NBINS = 7
BIN_NAMES = ("slope_00_40", "slope_40_50", "slope_50_60", "slope_60_70", "slope_70_80", "slope_80_90", "slope_90_180")
FLOOR = 1e-300                       # fit_cops.FLOOR
A_UNIT = 1e-6                        # grain_crossing.A_UNIT
B_UNIT = 1e-6 / 100.0                # grain_crossing.B_UNIT
GAP_SEARCH_S = 3600.0                # grain_crossing.GAP_SEARCH_S
CROSSING_WEIGHT = 10.0               # user 2026-10-27 (D2: --crossing-weight 10)
USER_T0 = (182.09922, 136.11149, 133.0)
TOL = 1e-10                          # D2's ftol / xtol / gtol
LAMBDA_FACTORS = {"weak": 0.1, "medium": 1.0, "strong": 10.0}
C_REF_C9 = 4.268806657481027e-03     # bounded softplus c9 final data cost (= the lead's 4.2688066575e-03)
MULT_PENALTY_DEFINITION = ("P_m = f_m C_ref sum_b (ln m_b)^2 over the six free slope multipliers (user 2026-11-14 ~17:30 "
                           "CST, option C), as residual rows sqrt(2 f_m C_ref) ln m_b (0.5 r.r scale), exact Jacobian "
                           "sqrt(2 f_m C_ref) / m_b; C_ref = %.15e (bounded softplus c9 data cost)" % C_REF_C9)
LOG = "localnn:"


def say(text):
    print("%s [%s] %s" % (LOG, dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(bufsize), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


jsonable = pp.jsonable


# ------------------------------------------------------------------------------------------------ inputs
def read_base_rows(base):
    """Base COPS_model.dat: header lines, column names, rows (token lists), iPoint, UTC, et."""
    lines = (Path(base) / "COPS_model.dat").read_text().splitlines()
    names, data = None, []
    for line in lines:
        if line.startswith("#iPoint"):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip():
            continue
        data.append(line.split())
    c = {n: i for i, n in enumerate(names)}
    return {"names": names, "cols": c, "rows": data, "iPoint": np.array([int(v[c["iPoint"]]) for v in data]),
            "UTC": [v[c["UTC"]] for v in data], "et": np.array([float(v[c["et"]]) for v in data]),
            "refNG": np.array([float(v[c["refNG"]]) for v in data]), "refRG": np.array([float(v[c["refRG"]]) for v in data])}


def load_grid(gdir):
    with np.load(Path(gdir) / "grid.npz", allow_pickle=False) as z:
        g = {k: np.asarray(z[k]) for k in z.files if k != "meta_json"}
        g["meta"] = json.loads(str(z["meta_json"]))
    n = len(g["node_xyz"])
    g["n"] = n
    g["Lg"] = sp.csr_matrix((g["Lg_data"], g["Lg_indices"], g["Lg_indptr"]), shape=(n, n))
    g["M"] = sp.csr_matrix((g["M_data"], g["M_indices"], g["M_indptr"]), shape=(n, n))
    g["Dg"] = sp.csr_matrix((g["Dg_data"], g["Dg_indices"], g["Dg_indptr"]), shape=tuple(g["Dg_shape"]))
    return g


def load_design(ddir):
    with np.load(Path(ddir) / "design.npz", allow_pickle=False) as z:
        d = {k: np.asarray(z[k]) for k in z.files if k != "meta_json"}
        d["meta"] = json.loads(str(z["meta_json"]))
    return d


def facet_field(grid, c):
    """sum_j c_j phi_j(e) at every facet."""
    fn, fw = grid["facet_node"], grid["facet_w"]
    return c[fn[:, 0]] * fw[:, 0] + c[fn[:, 1]] * fw[:, 1] + c[fn[:, 2]] * fw[:, 2]


# ------------------------------------------------------------------------------------------------ data layout
class Data:
    """Rows, references, selection (D2's restart), late layout (W10x5), grain timing and crossing pairs."""

    def __init__(self, base, d2_dir, point_table, design):
        self.base = Path(base)
        b = read_base_rows(base)
        self.table = b
        with np.load(Path(d2_dir) / "fit_restart.npz", allow_pickle=False) as z:
            pid = np.asarray(z["point_ids"], np.int64)
            self.ref_ng, self.ref_rg = np.asarray(z["ref_ng"], float), np.asarray(z["ref_rg"], float)
            self.use_ng, self.use_rg = np.asarray(z["use_ng"]), np.asarray(z["use_rg"])
            self.d2_x = np.asarray(z["x"], float)
            self.d2_names = json.loads(str(z["meta_json"]))["names"]
        if not np.array_equal(pid, b["iPoint"]) or not np.array_equal(design["iPoint"], b["iPoint"]):
            raise SystemExit("localnn: D2's rows / the design rows differ from the base rows")
        self.ref_check = {"refNG_equal_base": bool(np.array_equal(self.ref_ng, b["refNG"])),
                          "refRG_equal_base": bool(np.array_equal(self.ref_rg, b["refRG"]))}
        self.ipoint = b["iPoint"]
        self.et = b["et"]
        self.utc = b["UTC"]
        self.nrow = len(self.ipoint)
        self.useNG = (self.use_ng > 0) & (self.ref_ng > 0)
        self.useRG = (self.use_rg > 0) & (self.ref_rg > 0)
        self.NNG, self.NRG = int(self.useNG.sum()), int(self.useRG.sum())
        self.wNG, self.wRG = 1.0 / np.sqrt(self.NNG), 1.0 / np.sqrt(self.NRG)
        self.lnNG, self.lnRG = np.log(self.ref_ng[self.useNG]), np.log(self.ref_rg[self.useRG])
        self.rows_ng, self.rows_rg = np.flatnonzero(self.useNG), np.flatnonzero(self.useRG)
        # grain timing (grain_crossing.GrainProblem.__init__)
        et_ng = np.sort(self.et[self.useNG])
        idx = np.nonzero(et_ng >= et_ng[-1] - GAP_SEARCH_S)[0]
        k = idx[int(np.argmax(np.diff(et_ng[idx])))]
        self.et_before, self.et_end = float(et_ng[k]), float(et_ng[k + 1])
        # late layout (positive_problem, unchanged)
        self.point_table = str(Path(point_table).resolve())
        self.layout = pp.LateLayout(pp.read_point_table(point_table), pp.read_run_points(base), self.ipoint)
        lay = self.layout
        self.late_rows, self.grid_rows = lay.late_rows, lay.grid_rows
        W = lay.W.tocsr()
        W.sort_indices()
        self.W = W
        # row operator E (rows x base rows): identity, late observation rows = W over the grid rows
        Ei = sp.identity(self.nrow, format="lil")
        for kk, r in enumerate(self.late_rows):
            Ei[int(r), int(r)] = 0.0
            lo, hi = W.indptr[kk], W.indptr[kk + 1]
            for j, w in zip(W.indices[lo:hi], W.data[lo:hi]):
                Ei[int(r), int(self.grid_rows[j])] = w
        self.E = Ei.tocsr()
        if np.any(self.useNG[self.grid_rows]) or np.any(self.useRG[self.grid_rows]):
            raise SystemExit("localnn: a grid row is used")
        # crossing pairs (grain_crossing.pairs)
        both = self.useNG & self.useRG
        iNG = np.cumsum(self.useNG) - 1
        iRG = np.cumsum(self.useRG) - 1 + self.NNG
        kb = np.nonzero(both)[0]
        etb = self.et[kb]
        sign = np.where(etb >= self.et_end, -1.0, np.where(etb <= self.et_before, 1.0, 0.0))
        keep = sign != 0.0
        kb = kb[keep]
        self.cross = {"iNG": iNG[kb], "iRG": iRG[kb], "sign": sign[keep],
                      "w": CROSSING_WEIGHT * np.sqrt(self.wNG * self.wRG)}
        self.summary = {"rows": self.nrow, "used_NG": self.NNG, "used_RG": self.NRG, "t_end_et": self.et_end,
                        "last_ng_before_gap_et": self.et_before, "late_layout": lay.summary,
                        "crossing_pairs_before_gap": int(np.sum(self.cross["sign"] > 0)),
                        "crossing_pairs_after_gap": int(np.sum(self.cross["sign"] < 0)), "references": self.ref_check,
                        "point_table": {"path": self.point_table, "sha256": sha256_file(point_table)}}

    def grain(self, A_uPa, B_uPa, et):
        """grain_crossing.GrainProblem.grain (values and one-sided derivatives)."""
        A, B = A_uPa * A_UNIT, B_uPa * B_UNIT
        s = et - self.et_end
        on = s >= 0.0
        g = np.where(on, A - B * np.maximum(s, 0.0), 0.0)
        dA = np.where(on & (g >= 0.0), A_UNIT, 0.0)
        dB = np.where(on & (g > 0.0), -B_UNIT * np.maximum(s, 0.0), 0.0)
        return np.maximum(g, 0.0), dA, dB

    def model_rows(self, sNG, sRG, A_uPa, B_uPa):
        """Per base row: the fitted model (late observations W10x5 of the surface + grain(t_obs), others instantaneous
        surface + grain(t)), the surface-only model, the grain."""
        sNG_t, sRG_t = self.E @ sNG, self.E @ sRG
        g = self.grain(A_uPa, B_uPa, self.et)[0]
        return {"pNG": sNG_t + g, "pRG": sRG_t, "surface_NG": sNG_t, "surface_RG": sRG_t, "grain": g,
                "inst_NG": sNG, "inst_RG": sRG}

    def costs(self, pNG, pRG):
        """Data residuals and cost on the used rows (fit_cops), crossing penalty (grain_crossing); pNG includes the grain."""
        r = np.concatenate([self.wNG * (np.log(np.maximum(pNG[self.useNG], FLOOR)) - self.lnNG),
                            self.wRG * (np.log(np.maximum(pRG[self.useRG], FLOOR)) - self.lnRG)])
        q, d = self.crossing(r)
        return {"r": r, "data_cost": 0.5 * float(r @ r), "crossing_cost": 0.5 * float(q @ q),
                "rms_ln_NG": float(np.sqrt(np.mean((r[:self.NNG] / self.wNG) ** 2))),
                "rms_ln_RG": float(np.sqrt(np.mean((r[self.NNG:] / self.wRG) ** 2))),
                "wrong_side": int(np.sum(d > 0))}

    def crossing(self, r):
        c = self.cross
        lnNG = r[c["iNG"]] / self.wNG + self.lnNG[c["iNG"]]
        lnRG = r[c["iRG"]] / self.wRG + self.lnRG[c["iRG"] - self.NNG]
        d = c["sign"] * (lnNG - lnRG)
        return c["w"] * np.maximum(d, 0.0), d

    def late_rms(self, r):
        late = set(int(v) for v in self.late_rows)
        mng = np.array([int(v) in late for v in self.rows_ng])
        mrg = np.array([int(v) in late for v in self.rows_rg])
        rn, rr = r[:self.NNG] / self.wNG, r[self.NNG:] / self.wRG
        return {"late_rms_ln_NG": float(np.sqrt(np.mean(rn[mng] ** 2))), "late_rms_ln_RG": float(np.sqrt(np.mean(rr[mrg] ** 2))),
                "late_points": [int(mng.sum()), int(mrg.sum())]}


# ------------------------------------------------------------------------------------------------ the fit problem
class Problem:
    def __init__(self, data, design, grid, templates, slope_mult_min, A_s, smooth_lambda=0.0, mult_penalty_f=0.0):
        self.D, self.grid = data, grid
        self.templates = bool(templates)
        self.smm = float(slope_mult_min)
        self.A_s = float(A_s)
        self.n = grid["n"]
        D = data
        # effective design rows of the used points (W10x5 applied to the late observations)
        EN, ER = D.E[D.rows_ng], D.E[D.rows_rg]
        self.GN = np.stack([EN @ design["G_NG"][:, b, :] for b in range(NBINS)])       # (7, NNG, n)
        self.GR = np.stack([ER @ design["G_RG"][:, b, :] for b in range(NBINS)])
        tN = EN @ design["T_NG"] if self.templates else np.zeros((len(D.rows_ng), NBINS))
        tR = ER @ design["T_RG"] if self.templates else np.zeros((len(D.rows_rg), NBINS))
        self.tN, self.tR = tN.T.copy(), tR.T.copy()                                    # (7, N)
        self.CN, self.CR = EN @ design["C_NG"], ER @ design["C_RG"]
        self.et_ng = D.et[D.rows_ng]
        # all-row design (products)
        self.design = design
        # parameters
        self.names = ["c%04d" % j for j in range(self.n)] + ["c_" + b for b in BIN_NAMES[1:]] + \
                     ["kappa", "grain_A_uPa", "grain_B_uPa_per_100s"]
        self.ic = slice(0, self.n)
        self.im = slice(self.n, self.n + NBINS - 1)
        self.ik, self.iA, self.iB = self.n + NBINS - 1, self.n + NBINS, self.n + NBINS + 1
        self.lo = np.concatenate([np.zeros(self.n), np.full(NBINS - 1, self.smm), [0.0, 0.0, 0.0]])
        self.hi = np.full(len(self.lo), np.inf)
        # smoothness factor
        Lg = grid["Lg"].toarray()
        mu, V = np.linalg.eigh(0.5 * (Lg + Lg.T))
        keep = mu > 1e-12 * mu.max()
        self.R = (np.sqrt(mu[keep])[:, None] * V[:, keep].T)
        self.Lg = grid["Lg"]
        self.smooth = {"eigenvalues_kept": int(keep.sum()), "eigenvalues_dropped": int((~keep).sum()),
                       "min_kept": float(mu[keep].min()), "max": float(mu.max()), "min_eig": float(mu.min()),
                       "R_reconstruction_max_abs": float(np.max(np.abs(self.R.T @ self.R - Lg)))}
        self.lam = float(smooth_lambda)
        # multiplier penalty (user 2026-11-14 ~17:30 CST; off by default: no rows, the default path unchanged)
        self.fm = float(mult_penalty_f)
        if self.fm < 0:
            raise SystemExit("--mult-penalty-f must be >= 0")
        self.cref_m = C_REF_C9
        self.wm = math.sqrt(2.0 * self.fm * self.cref_m) if self.fm > 0 else 0.0
        self.nfev = 0
        self.best = None

    # ---- unpack
    def unpack(self, x):
        c = x[self.ic] * self.A_s
        m = np.concatenate([[1.0], x[self.im]])
        return c, m, float(x[self.ik]), float(x[self.iA]), float(x[self.iB])

    def surface_used(self, c, m, kappa):
        sN = np.zeros(self.GN.shape[1])
        sR = np.zeros(self.GR.shape[1])
        uN = np.empty((NBINS, self.GN.shape[1]))
        uR = np.empty((NBINS, self.GR.shape[1]))
        for b in range(NBINS):
            uN[b] = self.GN[b] @ c + self.tN[b]
            uR[b] = self.GR[b] @ c + self.tR[b]
            sN += m[b] * uN[b]
            sR += m[b] * uR[b]
        return sN + kappa * self.CN, sR + kappa * self.CR, uN, uR

    def pressures_used(self, x):
        c, m, kappa, A, B = self.unpack(x)
        sN, sR, uN, uR = self.surface_used(c, m, kappa)
        g, dA, dB = self.D.grain(A, B, self.et_ng)
        return sN + g, sR, (uN, uR, g, dA, dB)

    def data_residuals(self, x):
        pNG, pRG, _ = self.pressures_used(x)
        D = self.D
        return np.concatenate([D.wNG * (np.log(np.maximum(pNG, FLOOR)) - D.lnNG),
                               D.wRG * (np.log(np.maximum(pRG, FLOOR)) - D.lnRG)])

    def data_jacobian(self, x):
        c, m, kappa, A, B = self.unpack(x)
        pNG, pRG, (uN, uR, g, dA, dB) = self.pressures_used(x)
        D = self.D
        fN = D.wNG / np.maximum(pNG, FLOOR)
        fR = D.wRG / np.maximum(pRG, FLOOR)
        HN = np.tensordot(m, self.GN, axes=(0, 0)) * self.A_s
        HR = np.tensordot(m, self.GR, axes=(0, 0)) * self.A_s
        JN = np.empty((len(pNG), len(x)))
        JR = np.empty((len(pRG), len(x)))
        JN[:, self.ic] = HN * fN[:, None]
        JR[:, self.ic] = HR * fR[:, None]
        JN[:, self.im] = (uN[1:].T) * fN[:, None]
        JR[:, self.im] = (uR[1:].T) * fR[:, None]
        JN[:, self.ik] = self.CN * fN
        JR[:, self.ik] = self.CR * fR
        JN[:, self.iA] = dA * fN
        JN[:, self.iB] = dB * fN
        JR[:, self.iA] = 0.0
        JR[:, self.iB] = 0.0
        return np.vstack([JN, JR])

    def smooth_residuals(self, x):
        return math.sqrt(2.0 * self.lam) * (self.R @ x[self.ic])

    def penalty(self, x):
        xc = x[self.ic]
        return self.lam * float(xc @ (self.Lg @ xc))

    def mult_residuals(self, x):
        return self.wm * np.log(x[self.im])

    def mult_penalty(self, x):
        lm = np.log(x[self.im])
        return self.fm * self.cref_m * float(lm @ lm)

    def fun(self, x):
        r = self.data_residuals(x)
        q, _ = self.D.crossing(r)
        out = np.concatenate([r, q, self.smooth_residuals(x)]) if self.lam > 0 else np.concatenate([r, q])
        if self.fm > 0:
            out = np.concatenate([out, self.mult_residuals(x)])
        self.nfev += 1
        obj = 0.5 * float(out @ out)
        if self.best is None or obj < self.best[0]:
            self.best = (obj, np.array(x, copy=True))
        if self.nfev == 1 or self.nfev % 50 == 0:
            say("  eval %d: objective %.10e data %.10e" % (self.nfev, obj, 0.5 * float(r @ r)))
        if self.deadline is not None and time.time() > self.deadline:
            raise WallTime()
        return out

    deadline = None

    def jac(self, x):
        r = self.data_residuals(x)
        J = self.data_jacobian(x)
        c = self.D.cross
        _, d = self.D.crossing(r)
        Jq = c["w"] * c["sign"][:, None] * (J[c["iNG"]] / self.D.wNG - J[c["iRG"]] / self.D.wRG)
        Jq[d <= 0.0] = 0.0
        blocks = [J, Jq]
        if self.lam > 0:
            S = np.zeros((self.R.shape[0], len(x)))
            S[:, self.ic] = math.sqrt(2.0 * self.lam) * self.R
            blocks.append(S)
        if self.fm > 0:
            Mr = np.zeros((NBINS - 1, len(x)))
            Mr[np.arange(NBINS - 1), np.arange(self.im.start, self.im.stop)] = self.wm / x[self.im]
            blocks.append(Mr)
        return np.vstack(blocks)

    def evaluate(self, x):
        """Every cost term and the per-row model of all base rows (design route)."""
        c, m, kappa, A, B = self.unpack(x)
        sNG, sRG = all_rows_surface(self.design, c, m, kappa, self.templates)
        rows = self.D.model_rows(sNG, sRG, A, B)
        cst = self.D.costs(rows["pNG"], rows["pRG"])
        cst["penalty"] = self.penalty(x)
        cst["objective"] = cst["data_cost"] + cst["crossing_cost"] + cst["penalty"]
        if self.fm > 0:
            cst["mult_penalty"] = self.mult_penalty(x)
            cst["objective"] = cst["objective"] + cst["mult_penalty"]
        cst.update(self.D.late_rms(cst["r"]))
        return cst, rows


class WallTime(Exception):
    pass


def all_rows_surface(design, c, m, kappa, templates):
    """Instantaneous surface pressures of every base row from the design matrices."""
    sN = np.einsum("b,rbj,j->r", m, design["G_NG"], c)
    sR = np.einsum("b,rbj,j->r", m, design["G_RG"], c)
    if templates:
        sN = sN + design["T_NG"] @ m
        sR = sR + design["T_RG"] @ m
    return sN + kappa * design["C_NG"], sR + kappa * design["C_RG"]


# ------------------------------------------------------------------------------------------------ start
def nnls_projection(grid, area, target, chunk=None):
    """min_{c >= 0} sum_e area_e (sum_j c_j phi_j(e) - target_e)^2 with the exact facet mass matrix."""
    fn, fw = grid["facet_node"], grid["facet_w"]
    n = grid["n"]
    rhs = np.zeros(n)
    for k in range(3):
        rhs += np.bincount(fn[:, k], weights=area * fw[:, k] * target, minlength=n)
    M = grid["M"].toarray()
    ok = np.diag(M) > 0
    L = np.linalg.cholesky(M[np.ix_(ok, ok)])
    y = sla.solve_triangular(L, rhs[ok], lower=True)
    c_ok, resid = nnls(L.T, y, maxiter=50 * int(ok.sum()))
    c = np.zeros(n)
    c[ok] = c_ok
    u = facet_field(grid, c)
    tw = float(np.sum(area * target ** 2))
    return c, {"nodes_without_support": int((~ok).sum()), "nodes_at_zero": int(np.sum(c == 0.0)),
               "relative_area_weighted_misfit": float(np.sqrt(np.sum(area * (u - target) ** 2) / tw)),
               "target_area_integral": float(np.sum(area * target)), "fit_area_integral": float(np.sum(area * u)),
               "nnls_residual": float(resid)}


def lambda_from_factor(factor, data_cost_start, xc0, Lg):
    e = float(xc0 @ (Lg @ xc0))
    return factor * data_cost_start / e, e


# ------------------------------------------------------------------------------------------------ verification
def read_table_cols(path, wanted):
    names, out = None, {}
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        c = {n: i for i, n in enumerate(names)}
        out[int(t[c["iPoint"]])] = {k: t[c[k]] for k in wanted}
    return out


def compare(model, table_strings, ipoints, fmt):
    """max relative difference and the number of rows whose printed value differs (same format)."""
    rel, diff = 0.0, 0
    worst = None
    for v, ip in zip(model, ipoints):
        s = table_strings.get(int(ip))
        if s is None:
            continue
        ref = float(s)
        if ref == 0.0 and v == 0.0:
            continue
        d = abs(v / ref - 1.0)
        if d > rel:
            rel, worst = d, int(ip)
        if (fmt % v) != s and float(fmt % v) != ref:
            diff += 1
    return {"max_rel": rel, "argmax_iPoint": worst, "printed_differences": diff}
