#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_finetune.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_finetune.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: d57e6c43dc78232678270a626b9d89b5fe0d3832d2abb2f0f1a9dc6d6bb085c2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Per-element fine tuning of a positive-activity fit (user decisions 2026-11-07/08): A_e = exp(s_e + delta_e) for the
#   parameter facets, delta_e >= ln(r A_orig,e / A_e) (min with 0; r = floor ratio, 0.002), no cap; grain A, B refitted;
#   objective = the fit's full cost with its own lambda on the total ln A (data, W10x5 late window, grain, crossing,
#   Tikhonov). --facets late_only --window-min 5 (user 2026-11-08: only the facets seen exclusively in the last 5 min;
#   strict support, earlier pressures bitwise unchanged) or all_seen. Matrix-free TRF + LSMR (J v: the unchanged fused
#   kernel with one basis column; J^T u: kernels/finetune.cu). Validation checks 1-8, zero-sensitivity handling;
#   segments with resume (exit 75 at the wall time). Writes finetune_state.npz (delta, held_prev, held_next, grain,
#   mesh_identity), finetune_result.txt, COPS_model.dat, late_window_model.dat, pressure_comparison.txt, fine_tuning.json.
#   This copy = model_ext_frozen/ft_20261108T034023Z (the lambda-0 chain); the other five chains ran
#   ft_20261108T032615Z (756a7ce4: the same code without the zero-sensitivity handling; results bitwise equal when the
#   parameter block is not all zero, jobs 122594 / 122705).
# USED BY: jobs/positive_finetune.sbatch (jobs/positive_cycles.sh).
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
"""Per-element fine tuning of a positive-activity COPS fit (model_ext/positive_fit.py; user decisions 2026-11-07).

    .venv/bin/python -W ignore <this> --fit <fit seg dir> --out <new dir> --floor-ratio R --devices 0,1,...
        [--max-nfev N] [--lsmr-maxiter K] [--wall-time slurm|SECONDS] [--validate | --validate-only]
        [--resume-from <finetune_state.npz | finetune_checkpoint.npz of the same fit>]

Job script: jobs/positive_finetune.sbatch (environment from the fit, provenance). CPU tests:
model_ext/test_positive_finetune.py. Local, untracked extension: nothing in the git tree is edited; positive_fit.py,
positive_problem.py, positive_activity.py and kernels/{pressure,activity}.cu are imported/used unchanged.

Model (user decisions 2026-11-07, "activity of individual surface element can be increased or decreased down to the
floor value")
* Parameters: every facet seen from a point the objective uses gets its own exponent increment delta_e,
      A_e = exp(s_e),   s = s_fit + held + delta,   s_fit = Y a + c_cliff T_cliff (the fit's exponent field),
  held = the tuned fit's held per-element term (activity_exp_facets.npz key "held"; zeros when absent), delta = 0 at
  the start; unseen facets unchanged. "Seen": the reference predicate of fit_cops.Model.seen_faces (Kn > 0 | KPhi > 0
  | KV > 0) on the records of the points whose model enters the objective: (a) observation rows with NG or RG used
  that are not late-window observations, (b) the 1-s grid rows with a nonzero W10x5 weight for a used late
  observation. A late observation's own records do not enter the objective (its model is W10x5 of the grid rows);
  the facets seen only from them would have an exactly zero data Jacobian -- they are counted in fine_tuning.json
  ("seen") and are NOT parameters (--seen-points objective, default; --seen-points objective+late includes them).
* Lower bound A_e >= r A_orig,e (r = --floor-ratio; must equal COPS_POS_FLOOR_RATIO of the tuned fit and of this
  process, else refused): delta_e >= ln(r A_orig,e / A_e(start)); a facet already below it at the start gets
  min(0, that) (feasible start); A_orig = 0: no lower bound. No upper bound. Grain A, B >= 0 (as in the fit).
* Objective = the fit's full objective, unchanged in form: the data rows of the fit's selection (log residuals,
  equal weight per gauge), the late observations (>= 10:33:10 UTC) compared with W10x5 of the surface model with the
  grain added at the time stamp (positive_problem.LateLayout / WindowOperator), the grain A, B refitted together
  with delta (start = the fit's values), the crossing penalty with the fit's weight (grain_crossing.py's
  definitions, replicated below and tested against it), and lambda/2 D^2 with the fit's own lambda,
      D^2 = sum_all area_e (s_e - s0_e)^2 / sum area     (s0 = the fit's prior, activity_exp_facets.npz "s0"),
  as residual rows sqrt(lambda area_e / sum area) (s_e - s0_e) over ALL facets (unseen ones are constant);
  lambda = 0: no rows. Temperatures, a, c_cliff and kappa are held at the fit's values.
  Costs: "total" = data + regularization + crossing (what the optimizer minimizes). fit_result.txt final_cost
  (= fit_manifest optimizer_results cost) is data + regularization: grain_crossing.py strips the crossing term from
  the reported cost; the fit's crossing cost is fit_manifest model_extension.final.penalty_cost.

Method (why): matrix-free Gauss-Newton trust region, scipy least_squares(method="trf", tr_solver="lsmr",
jac=LinearOperator, bounds). With ~9.2 M parameters and ~3000 data + crossing rows (+ 12.4 M diagonal Tikhonov rows
for lambda > 0) neither J nor J^T J can be formed; TRF is SciPy's bounded large-scale least-squares method (the
same trust-region algorithm the fits use), it handles the millions of lower bounds by Coleman-Li scaling and
reflected steps, and it needs only J v and J^T u. For lambda = 0 the GN system is underdetermined; LSMR started at 0
returns its minimum-norm (in TRF's scaled variables) solution, as the user asked ("the minimum-norm Gauss-Newton
steps from the start decide"); --lsmr-maxiter caps its iterations per outer iteration (tr_options maxiter).
x_scale = 1 for every parameter (delta in nepers, grain in the fit's units 1e-6 Pa and 1e-6 Pa / 100 s): 'jac' is
not available with a LinearOperator and a per-facet column norm would let insensitive facets take huge steps.
GPU side (model_ext/positive_activity.py patch, mixed precision, all GPUs given by --devices, one bind of the fit's
base rows exactly as the fit's MultiGpuBackend with host activity; per evaluation the driver writes s to the
backend's shared activity buffer and every worker forms A = exp(s) itself):
* forward pressures: the unchanged evaluator (eval_fused_m_p5) -> bitwise the fit's model for the same s; with the
  CO2 flux floor (COPS_POS_FLOOR_CO2=1 of a ratio-1e-5 fit, 2026-11-08) the evaluator of positive_activity.py
  (eval_fused_m_p5f, kappa held at the fit's value) -> again the fit's own model;
* J v: the unchanged eval_fused_m_c<k> (k <= 2 columns per pass) with colbasis = v; exp activity -> column factor
  A_f v_f, records at the H2O floor contribute 0 (the fit's own Jacobian code with a per-facet basis column);
* J^T u: the new kernel kernels/finetune.cu ft_jtu_m (one launch per point, in point order: deterministic), per-
  device partial arrays (shared memory) summed on the host in device order; CO2 does not depend on the H2O activity
  (also with the CO2 flux floor: F_CO2 = max(kappa Q_CO2 ACO2, r q_max,CO2 ACO2) depends on kappa, Q_CO2 and ACO2
  only, so J v / J^T u are unchanged; the 2026-11-08 pressure.cu adds kernels and leaves static_factors_m, the H2O
  floor and the column code, which ft_jtu_m transposes, unchanged: PRESSURE_COMPATIBLE_SHA256);
* the W10x5 operator and its transpose, the grain and crossing columns and the Tikhonov diagonal on the host.
The worker requests ft_attach, ft_eval, ft_jv, ft_jtu, ft_seen, ft_stats are added to
cops_cuda.transport.WorkerState.handle in every process that runs this module (the spawned GPU workers re-import it
as __mp_main__; in-process patch, like positive_activity.py's).

Products in --out (created; it may already exist only with provenance-* directories; the job script passes
OUT/seg<NN> and links OUT/final and OUT/<product> to the finished segment):
* finetune_state.npz: delta (this tuning's change, per facet, 0 on unseen facets), held_prev (the tuned fit's held
  term, zeros if none), held_next = held_prev + delta (the term the next cycle's fit holds), s_final, A_final,
  lower_bound (bound on delta; -inf: none; NaN: not a parameter = unseen), seen (bool, the parameters),
  grain_A_uPa, grain_B_uPa_per_100s, meta_json (fit dir, its restart SHA-256, configuration, lambda, costs);
* finetune_result.txt: initial (= the fit) and final total/data/regularization/crossing cost, rms ln NG/RG, the
  late-window rms, D, nfev/iterations, termination, facets (number, area) at the lower bound, grain;
* COPS_model.dat: the fit's table with pNG/pRG of EVERY row recomputed (late observations W10x5(surface) + grain;
  grid rows instantaneous surface + grain, flagged by the fit's PositiveFitGridRows header line and listed in
  late_window_model.dat), so figures_cluster/scripts/positive_table.py works on --out; late_window_model.dat (the
  fit's format);
* pressure_comparison.txt: per row the fit's and the tuned model, references, flags;
* fine_tuning.json: settings, identities, environment, seen statistics, optimizer record, timings, memory;
* validation.json (--validate / --validate-only): the checks (1)-(6) of the request (see validate()).
* finetune_checkpoint.npz (only with --checkpoint-minutes > 0; replaced atomically at most every M minutes while the
  optimizer improves; a resume aid for --resume-from, not a product).
Wall time and exit codes: --wall-time slurm stops at the best evaluated state --save-seconds (+ --wall-margin) before
the job's end and writes the products (finetune_result.txt "stop wall_time", "complete no"; state meta stop/complete/
nfev_total/segments). Exit 0: finished (stop "tolerance" = ftol/xtol/gtol, or "max_nfev") with all products written;
75: wall-time stop (products of the best state written; continue with --segment k+1 --nfev-before <evaluations so
far> --resume-from <this --out>/finetune_state.npz into a new --out; jobs/positive_finetune.sbatch does this when the
same command is resubmitted: OUT/seg<NN>, OUT/final and the product links); 3: a bound violated; 1: an error or a
failed validation (no fit).
"""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import resource
import socket
import sys
import time
import traceback
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP / "tools"), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

KERNEL_FILE = HERE / "kernels" / "finetune.cu"
PRESSURE_KERNEL = HERE / "kernels" / "pressure.cu"
PRESSURE_SOURCE_SHA256 = "714b9d908fee29c8135931ecd8a5209fba5ef1dd6905bec5b219cf5a0a4280f3"   # finetune.cu's source
# kernels/pressure.cu versions ft_jtu_m is the transpose of (the H2O code it copies is byte-identical in each):
PRESSURE_COMPATIBLE_SHA256 = {
    PRESSURE_SOURCE_SHA256: "the 2026-11-07 copy finetune.cu was derived from",
    "85842e95149817a73664f6e1afe543fd2dff5d884063124398f5887bdb08faf7":
        "2026-11-08 CO2 flux floor (header item 5): adds co2_flux_m and the kernels eval_fused_m_p5f/p14f (template "
        "flag CF); static_factors_m, the H2O floor and the activity/template column code are unchanged and the 41 "
        "earlier kernels compile to the same PTX (test_positive_finetune.py checks the text blocks)"}
NVRTC_OPTIONS = ("--std=c++14", "--fmad=false")        # cops_cuda.precision.REFERENCE_NVRTC_OPTIONS
PARAMETERIZATION_ID = "positive_per_element_finetune_v1"
FIT_PARAMETERIZATION_PREFIX = "positive_exp_activity"
# COPS_POS_FLOOR_CO2 (2026-11-08): the CO2 flux floor of the ratio-1e-5 fits; absent (None) in the earlier fits' records
ENV_KEYS = ("COPS_POS_ENABLE", "COPS_POS_NIGHT_RATIO", "COPS_POS_ACTIVITY", "COPS_POS_FLOOR_RATIO", "COPS_POS_FLOOR_CO2")
KMAX = 2                      # J v columns per GPU pass (shared-memory direction buffer)
BLOCK = 256                   # kernels/finetune.cu BLOCK
FLOOR = 1e-300                # fit_cops.FLOOR (checked against the reference module at run time)
GRAIN_NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
POLICY = ("user 2026-11-07: per-element fine tuning of a positive-activity fit; every facet seen from a point the "
          "objective uses gets delta_e, A = exp(s_fit + held + delta), delta = 0 at the start, unseen facets "
          "unchanged; A_e >= r A_orig,e (delta_e >= min(0, ln(r A_orig,e / A_e)); A_orig = 0: no bound), no upper "
          "bound, r = --floor-ratio = COPS_POS_FLOOR_RATIO of the fit; objective = the fit's full objective (data "
          "rows, W10x5 + grain in the late window, grain A, B >= 0 refitted, crossing penalty with the fit's weight, "
          "lambda/2 sum area (s - s0)^2 / sum area with the fit's lambda on the total ln A); temperatures, a, c, kappa "
          "held; the same COPS_POS_* environment as the fit. User 2026-11-08: MAXNFEV 500 per fine tuning, no upper "
          "cap (only the floor bound); x_scale 1 and LSMR_MAXITER 50 are the method choices of this tool (recorded); "
          "the host-side speed-up (parameter-space I/O, parallel per-device sums) keeps the results bitwise")
METHOD_CHOICES = {"x_scale": "1 for every parameter (delta in nepers, grain in 1e-6 Pa and 1e-6 Pa / 100 s)",
                  "lsmr_maxiter_default": 50, "optimizer": "scipy least_squares trf, tr_solver lsmr, LinearOperator",
                  "decided": "method choices of positive_finetune.py, recorded by the user on 2026-11-08"}


# ------------------------------------------------------------------------------------------------ small helpers
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def file_identity(path):
    p = Path(path)
    return {"path": str(p.resolve()), "sha256": sha256(p), "size_bytes": p.stat().st_size}


def array_sha256(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def log(message):
    print("positive_finetune: %s" % message, flush=True)


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        v = float(value)
        return v if math.isfinite(v) else repr(v)
    if isinstance(value, (str, int)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    return repr(value)


def rel(a, b):
    """|a - b| / max(|a|, |b|) (0 when both are 0)."""
    a, b = float(a), float(b)
    s = max(abs(a), abs(b))
    return 0.0 if s == 0.0 else abs(a - b) / s


def max_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0


# ------------------------------------------------------------------------------------------------ rows (host)
class Rows:
    """The fit's rows (every row of the base, the fit model's row order), its data selection, the W10x5 window, the
    grain and the crossing pairs, exactly as positive_problem.py / grain_crossing.py define them for the fit."""

    def __init__(self, ipoint, et, refNG, refRG, use_ng, use_rg, window, crossing_weight, gap_search_s,
                 grain_units):
        self.n = len(ipoint)
        self.ipoint = np.asarray(ipoint, dtype=np.int64)
        self.et = np.asarray(et, dtype=float)
        self.refNG, self.refRG = np.asarray(refNG, float), np.asarray(refRG, float)
        # fit_cops.Problem: useNG = (M.useNG > 0) & (M.refNG > 0) over the fit rows (= all rows here)
        self.useNG = (np.asarray(use_ng) > 0) & (self.refNG > 0)
        self.useRG = (np.asarray(use_rg) > 0) & (self.refRG > 0)
        self.NNG, self.NRG = int(self.useNG.sum()), int(self.useRG.sum())
        self.wNG, self.wRG = 1.0 / np.sqrt(self.NNG), 1.0 / np.sqrt(self.NRG)
        self.lnNG, self.lnRG = np.log(self.refNG[self.useNG]), np.log(self.refRG[self.useRG])
        self.nd = self.NNG + self.NRG
        self.op = window                                     # positive_problem.WindowOperator (or a stand-in)
        self.is_late = np.zeros(self.n, bool)
        self.is_late[np.asarray(window.late_pos)] = True
        self.is_grid = np.zeros(self.n, bool)
        self.is_grid[np.asarray(window.grid_pos)] = True
        self.A_UNIT, self.B_UNIT = grain_units
        # grain onset (grain_crossing.GrainProblem.__init__)
        et_ng = np.sort(self.et[self.useNG])
        idx = np.nonzero(et_ng >= et_ng[-1] - gap_search_s)[0]
        k = idx[int(np.argmax(np.diff(et_ng[idx])))]
        self.et_before, self.et_end = float(et_ng[k]), float(et_ng[k + 1])
        # crossing pairs (grain_crossing.install.pairs)
        both = self.useNG & self.useRG
        iNG = np.cumsum(self.useNG) - 1
        iRG = np.cumsum(self.useRG) - 1 + self.NNG
        kk = np.nonzero(both)[0]
        etk = self.et[kk]
        sign = np.where(etk >= self.et_end, -1.0, np.where(etk <= self.et_before, 1.0, 0.0))
        keep = sign != 0.0
        kk = kk[keep]
        self.c_iNG, self.c_iRG, self.c_sign = iNG[kk], iRG[kk], sign[keep]
        self.c_rows = kk
        self.c_w = crossing_weight * np.sqrt(self.wNG * self.wRG)
        self.c_lnNG = self.lnNG[self.c_iNG]
        self.c_lnRG = self.lnRG[self.c_iRG - self.NNG]
        self.npen = len(self.c_sign)
        if len(np.unique(self.c_iNG)) != self.npen or len(np.unique(self.c_iRG)) != self.npen:
            raise SystemExit("crossing pairs: duplicate residual indices")

    # grain_crossing.GrainProblem.grain (x[-2] = A in 1e-6 Pa, x[-1] = B in 1e-6 Pa / 100 s)
    def grain(self, xA, xB, et=None):
        et = self.et if et is None else et
        A, B = xA * self.A_UNIT, xB * self.B_UNIT
        s = et - self.et_end
        on = s >= 0.0
        g = np.where(on, A - B * np.maximum(s, 0.0), 0.0)
        dA = np.where(on & (g >= 0.0), self.A_UNIT, 0.0)
        dB = np.where(on & (g > 0.0), -self.B_UNIT * np.maximum(s, 0.0), 0.0)
        return np.maximum(g, 0.0), dA, dB

    # GrainProblem.data_residuals (pNG includes the grain)
    def data_residuals(self, pNG, pRG):
        return np.concatenate([self.wNG * (np.log(np.maximum(pNG[self.useNG], FLOOR)) - self.lnNG),
                               self.wRG * (np.log(np.maximum(pRG[self.useRG], FLOOR)) - self.lnRG)])

    # grain_crossing.install.penalty
    def penalty(self, r):
        lnNG = r[self.c_iNG] / self.wNG + self.c_lnNG
        lnRG = r[self.c_iRG] / self.wRG + self.c_lnRG
        d = self.c_sign * (lnNG - lnRG)
        return self.c_w * np.maximum(d, 0.0), d

    def window(self, y):
        return self.op.apply(y)

    def window_T(self, u):
        """Adjoint of WindowOperator.apply: late rows -> their grid rows (W^T), late rows' own entries 0."""
        op = self.op
        out = np.array(u, dtype=float, copy=True)
        late = out[op.late_pos].copy()
        out[op.late_pos] = 0.0
        out[op.grid_pos] += op.W.T @ late
        return out

    def objective_points(self):
        """(rows whose model enters the objective, late observation rows that are used)."""
        used = self.useNG | self.useRG
        op = self.op
        late_used = used[op.late_pos]
        W = op.W[np.flatnonzero(late_used)].tocsr()
        grid_cols = np.unique(W.indices[W.data != 0.0]) if W.nnz else np.zeros(0, np.int64)
        grid_used = np.zeros(self.n, bool)
        grid_used[np.asarray(op.grid_pos)[grid_cols]] = True
        own = used & ~self.is_late
        late_rows = np.zeros(self.n, bool)
        late_rows[np.asarray(op.late_pos)[late_used]] = True
        return own | grid_used, late_rows


# ------------------------------------------------------------------------------------------------ objective
class StopFit(Exception):
    """Raised inside the optimizer's callbacks at the wall-time limit."""


class State:
    """One evaluated parameter vector: activity exponent, pressures, residuals, costs."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


class Tracker:
    """Wall clock (stop at the best evaluated state), best state, optional rolling checkpoint, counters."""

    def __init__(self, deadline=None, save_seconds=600.0, checkpoint=None, checkpoint_seconds=0.0):
        self.deadline, self.save_seconds = deadline, float(save_seconds)
        self.checkpoint, self.checkpoint_seconds = checkpoint, float(checkpoint_seconds)
        self.last_checkpoint = time.time()
        self.best = None
        self.stopped = None
        self.history = []
        self.checkpoints = 0
        self.t0 = time.time()

    def check(self):
        if self.deadline is not None and time.time() > self.deadline - self.save_seconds:
            self.stopped = utc_now()
            raise StopFit("wall-time limit: stop at the best evaluated state (%s)" % self.stopped)

    def evaluated(self, obj, st, kind):
        self.history.append({"kind": kind, "cost": st.cost, "elapsed_s": time.time() - self.t0})
        if not math.isfinite(st.cost):
            return
        if self.best is None or st.cost < self.best.cost:
            self.best = st
            if (self.checkpoint is not None and self.checkpoint_seconds > 0
                    and time.time() - self.last_checkpoint >= self.checkpoint_seconds):
                obj.write_checkpoint(self.checkpoint, st)
                self.last_checkpoint = time.time()
                self.checkpoints += 1


class Objective:
    """The fine tuning's residual vector [data (fit selection), Tikhonov rows (lambda > 0), crossing rows] over
    z = [delta (the parameter facets, increasing facet id), grain_A_uPa, grain_B_uPa_per_100s]."""

    def __init__(self, rows, engine, s_start, s0, area, A_orig, seen, lam, floor_ratio, tracker=None,
                 fit_identity=None):
        self.rows, self.engine = rows, engine
        self.s_start = np.asarray(s_start, float)
        self.s0 = np.asarray(s0, float)
        self.area = np.asarray(area, float)
        self.A_orig = np.asarray(A_orig, float)
        self.nf = len(self.area)
        self.seen = np.asarray(seen, bool)
        self.idx = np.flatnonzero(self.seen)
        self.npar = len(self.idx)
        self.n = self.npar + 2
        self.lam = float(lam)
        self.floor_ratio = float(floor_ratio)
        self.area_total = float(np.sum(self.area))
        self.reg_c = np.sqrt(self.lam * self.area / self.area_total) if self.lam > 0 else None
        self.reg_c_idx = self.reg_c[self.idx] if self.lam > 0 else None
        self.nreg = self.nf if self.lam > 0 else 0
        self.m = rows.nd + self.nreg + rows.npen
        self.tracker = tracker or Tracker()
        self.fit_identity = fit_identity or {}
        lb, below = lower_bounds(self.s_start[self.idx], self.A_orig[self.idx], self.floor_ratio)
        self.lb_delta = lb
        self.start_below_bound = below
        self.lb = np.concatenate([lb, [0.0, 0.0]])
        self.ub = np.full(self.n, np.inf)
        self.param_info = engine.set_parameters(self.idx)
        self.recent = []
        self.counters = {"fun": 0, "evaluations": 0, "jv_columns": 0, "jv_passes": 0, "jtu": 0}
        self.seconds = {"forward": 0.0, "jv": 0.0, "jtu": 0.0, "host_fun": 0.0, "host_jv": 0.0, "host_jtu": 0.0}

    # --- evaluation
    def evaluate(self, z):
        z = np.array(z, dtype=float, copy=True)
        if z.shape != (self.n,):
            raise ValueError("parameter vector of shape %s, need (%d,)" % (z.shape, self.n))
        t0 = time.perf_counter()
        s = self.s_start.copy()
        s[self.idx] += z[:self.npar]
        t1 = time.perf_counter()
        token, pNG_s, pRG_s = self.engine.forward(s)
        t2 = time.perf_counter()
        rows = self.rows
        pNG_w, pRG_w = rows.window(pNG_s), rows.window(pRG_s)
        g, dA, dB = rows.grain(z[-2], z[-1])
        pNG, pRG = pNG_w + g, pRG_w
        rd = rows.data_residuals(pNG, pRG)
        rr = self.reg_c * (s - self.s0) if self.lam > 0 else np.zeros(0)
        q, dpen = rows.penalty(rd)
        r = np.concatenate([rd, rr, q])
        with np.errstate(over="ignore", invalid="ignore"):
            costs = {"data": 0.5 * float(rd @ rd), "regularization": 0.5 * float(rr @ rr),
                     "crossing": 0.5 * float(q @ q)}
        costs["total"] = costs["data"] + costs["regularization"] + costs["crossing"]
        st = State(z=z, s=s, token=token, pNG_s=pNG_s, pRG_s=pRG_s, pNG_w=pNG_w, pRG_w=pRG_w, pNG=pNG, pRG=pRG,
                   g=g, dA=dA, dB=dB, rd=rd, rr=rr, q=q, dpen=dpen, r=r, costs=costs, cost=0.5 * float(r @ r))
        self.counters["evaluations"] += 1
        self.seconds["forward"] += t2 - t1
        self.seconds["host_fun"] += (t1 - t0) + (time.perf_counter() - t2)
        self.recent = [x for x in self.recent if x is not st][-2:] + [st]
        return st

    def find(self, z):
        z = np.asarray(z, float)
        for st in reversed(self.recent):
            if st.z.shape == z.shape and np.array_equal(st.z, z):
                return st
        if self.tracker.best is not None and np.array_equal(self.tracker.best.z, z):
            return self.tracker.best
        return None

    def fun(self, z):
        self.tracker.check()
        st = self.evaluate(z)
        self.counters["fun"] += 1
        self.tracker.evaluated(self, st, "fun")
        return st.r

    def jac(self, z, f=None):
        self.tracker.check()
        st = self.find(z)
        if st is None:
            st = self.evaluate(z)
        return Jacobian(self, st)

    def start_vector(self, xA, xB, delta=None):
        z = np.zeros(self.n)
        if delta is not None:
            z[:self.npar] = np.asarray(delta, float)[self.idx]
        z[-2], z[-1] = float(xA), float(xB)
        return z

    # --- diagnostics
    def summary(self, st):
        rows = self.rows
        rd = st.rd
        lnNG = rd[:rows.NNG] / rows.wNG
        lnRG = rd[rows.NNG:] / rows.wRG
        lateNG = rows.et[rows.useNG] >= rows.et_end
        lateRG = rows.et[rows.useRG] >= rows.et_end
        D = math.sqrt(float(np.sum(self.area * (st.s - self.s0) ** 2)) / self.area_total)
        delta = st.z[:self.npar]
        at_lb = np.isfinite(self.lb_delta) & (delta <= self.lb_delta + 1e-9 * np.maximum(1.0, np.abs(self.lb_delta)))
        aidx = self.area[self.idx]
        pen_wrong = st.dpen > 0
        return {"cost_total": st.costs["total"], "cost_data": st.costs["data"],
                "cost_regularization": st.costs["regularization"], "cost_crossing": st.costs["crossing"],
                "rms_ln_NG": float(np.sqrt(np.mean(lnNG ** 2))), "rms_ln_RG": float(np.sqrt(np.mean(lnRG ** 2))),
                "late_window_rms_ln_NG": float(np.sqrt(np.mean(lnNG[lateNG] ** 2))) if lateNG.any() else None,
                "late_window_rms_ln_RG": float(np.sqrt(np.mean(lnRG[lateRG] ** 2))) if lateRG.any() else None,
                "late_window_points": [int(lateNG.sum()), int(lateRG.sum())],
                "D_nepers": D,
                "grain_A_uPa": float(st.z[-2]), "grain_B_uPa_per_100s": float(st.z[-1]),
                "grain_A_Pa": float(st.z[-2] * rows.A_UNIT), "grain_B_Pa_per_s": float(st.z[-1] * rows.B_UNIT),
                "delta_range": [float(delta.min()), float(delta.max())] if self.npar else [0.0, 0.0],
                "delta_rms_area_weighted_parameters": float(np.sqrt(np.sum(aidx * delta ** 2) / np.sum(aidx)))
                if self.npar else 0.0,
                "facets_at_lower_bound": int(at_lb.sum()),
                "area_fraction_at_lower_bound_of_total": float(np.sum(aidx[at_lb]) / self.area_total),
                "area_fraction_at_lower_bound_of_parameters": float(np.sum(aidx[at_lb]) / np.sum(aidx))
                if self.npar else 0.0,
                "crossing_wrong_side_pairs": int(pen_wrong.sum()),
                "total_activity_over_A_orig": float(np.sum(self.area * np.exp(st.s)) /
                                                    np.sum(self.area * self.A_orig))}

    def write_checkpoint(self, path, st):
        delta = np.zeros(self.nf)
        delta[self.idx] = st.z[:self.npar]
        meta = dict(self.fit_identity, kind="positive_finetune_checkpoint", parameterization=PARAMETERIZATION_ID,
                    created_utc=utc_now(), cost=st.cost, costs=st.costs, seen_sha256=array_sha256(self.seen),
                    nfev=self.counters["fun"])
        tmp = Path(str(path) + ".tmp.%d" % os.getpid())
        with tmp.open("wb") as f:
            np.savez_compressed(f, delta=delta, seen=self.seen, grain_A_uPa=np.array(st.z[-2]),
                                grain_B_uPa_per_100s=np.array(st.z[-1]),
                                meta_json=np.array(json.dumps(jsonable(meta), sort_keys=True)))
        os.replace(str(tmp), str(path))
        log("checkpoint %s (cost %.10e, %d evaluations)" % (path, st.cost, self.counters["fun"]))


def zero_sensitivity(obj, st, seed=20261108):
    """Parameters whose data + crossing Jacobian column is exactly 0 at a state (J^T u with a random u on the data
    and crossing rows, 0 on the Tikhonov rows: a nonzero column gives a nonzero entry with probability 1; a facet whose
    records are all at the H2O floor gives an exact 0)."""
    rng = np.random.default_rng(seed)
    u = rng.standard_normal(obj.m)
    if obj.nreg:
        u[obj.rows.nd:obj.rows.nd + obj.nreg] = 0.0
    g = Jacobian(obj, st).rmatvec(u)[:obj.npar]
    z = g == 0.0
    a = obj.area[obj.idx]
    return {"parameters": int(obj.npar), "zero_sensitivity": int(z.sum()),
            "zero_sensitivity_area_fraction_of_parameters": float(a[z].sum() / a.sum()) if obj.npar else 0.0,
            "definition": "data + crossing Jacobian column exactly 0 (all the facet's records at the H2O floor or "
                          "outside the objective)"}


def lower_bounds(s_start, A_orig, ratio):
    """(lower bound on delta, start-below-bound mask): delta >= min(0, ln(r A_orig / exp(s_start))), A_orig > 0;
    A_orig = 0: -inf (no bound)."""
    s_start = np.asarray(s_start, float)
    A_orig = np.asarray(A_orig, float)
    lb = np.full(len(s_start), -np.inf)
    pos = A_orig > 0.0
    raw = np.log(ratio * A_orig[pos]) - s_start[pos]
    lb[pos] = np.minimum(0.0, raw)
    below = np.zeros(len(s_start), bool)
    below[np.flatnonzero(pos)[raw > 0.0]] = True
    return lb, below


try:
    from scipy.sparse.linalg import LinearOperator as _LinearOperator
except ImportError:                                       # pragma: no cover
    _LinearOperator = object


class Jacobian(_LinearOperator):
    """J of the full residual at one evaluated state (matrix-free; GPU for the H2O activity block)."""

    def __init__(self, obj, st):
        super().__init__(dtype=np.float64, shape=(obj.m, obj.n))
        self.obj, self.st = obj, st
        rows = obj.rows
        self.fN = rows.wNG / np.maximum(st.pNG, FLOOR)[rows.useNG]
        self.fR = rows.wRG / np.maximum(st.pRG, FLOOR)[rows.useRG]
        self.pen_coef = rows.c_w * rows.c_sign * (st.dpen > 0.0)

    def _matvec(self, v):
        v = np.asarray(v, float).ravel()
        out = np.empty(self.obj.m)
        self._columns(v[:, None], [out])
        return out

    def _matmat(self, V):
        V = np.asarray(V, dtype=float)
        if V.ndim == 1:
            V = V[:, None]
        cols = [np.empty(self.obj.m) for _ in range(V.shape[1])]
        self._columns(V, cols)
        return np.column_stack(cols)

    def _columns(self, V, outs):
        """J V[:, j] into outs[j] (k <= KMAX directions per GPU pass)."""
        obj, st, engine = self.obj, self.st, self.obj.engine
        obj.tracker.check()
        k = V.shape[1]
        for j0 in range(0, k, KMAX):
            kc = min(KMAX, k - j0)
            t0 = time.perf_counter()
            buf = engine.direction_buffer(kc)              # parameter space (k x npar)
            for j in range(kc):
                buf[j][:] = V[:obj.npar, j0 + j]
            t1 = time.perf_counter()
            JNG, JRG = engine.jv(st, kc)
            t2 = time.perf_counter()
            for j in range(kc):
                self._assemble(JNG[j], JRG[j], V[:, j0 + j], outs[j0 + j])
            obj.seconds["jv"] += t2 - t1
            obj.seconds["host_jv"] += (t1 - t0) + (time.perf_counter() - t2)
            obj.counters["jv_passes"] += 1
            obj.counters["jv_columns"] += kc

    def _assemble(self, jng_s, jrg_s, v, out):
        obj, st, rows = self.obj, self.st, self.obj.rows
        dpNG = rows.window(jng_s) + st.dA * v[-2] + st.dB * v[-1]
        dpRG = rows.window(jrg_s)
        nd = rows.nd
        dd = out[:nd]
        dd[:rows.NNG] = dpNG[rows.useNG] * self.fN
        dd[rows.NNG:] = dpRG[rows.useRG] * self.fR
        if obj.nreg:
            reg = out[nd:nd + obj.nreg]
            reg.fill(0.0)
            reg[obj.idx] = obj.reg_c_idx * v[:obj.npar]
        if rows.npen:
            out[nd + obj.nreg:] = self.pen_coef * (dd[rows.c_iNG] / rows.wNG - dd[rows.c_iRG] / rows.wRG)
        return out

    def _rmatvec(self, u):
        obj, st, rows, engine = self.obj, self.st, self.obj.rows, self.obj.engine
        obj.tracker.check()
        t0 = time.perf_counter()
        u = np.asarray(u, dtype=float).ravel()
        ud = u[:rows.nd].copy()
        if rows.npen:
            coef = self.pen_coef * u[rows.nd + obj.nreg:]
            ud[rows.c_iNG] += coef / rows.wNG
            ud[rows.c_iRG] -= coef / rows.wRG
        uNG = np.zeros(rows.n)
        uNG[rows.useNG] = ud[:rows.NNG] * self.fN
        uRG = np.zeros(rows.n)
        uRG[rows.useRG] = ud[rows.NNG:] * self.fR
        gA, gB = float(uNG @ st.dA), float(uNG @ st.dB)
        wN, wR = rows.window_T(uNG), rows.window_T(uRG)
        t1 = time.perf_counter()
        g = engine.jtu(st, wN, wR)                         # parameter space (engine buffer)
        t2 = time.perf_counter()
        out = np.empty(obj.n)
        if obj.nreg:
            np.add(g, obj.reg_c_idx * u[rows.nd:rows.nd + obj.nreg][obj.idx], out=out[:obj.npar])
        else:
            out[:obj.npar] = g
        out[-2], out[-1] = gA, gB
        obj.seconds["jtu"] += t2 - t1
        obj.seconds["host_jtu"] += (t1 - t0) + (time.perf_counter() - t2)
        obj.counters["jtu"] += 1
        return out

    def _rmatmat(self, U):
        U = np.asarray(U, float)
        return np.column_stack([self._rmatvec(U[:, j]) for j in range(U.shape[1])])


# ------------------------------------------------------------------------------------------------ GPU engine
class GpuEngine:
    """The fit's base rows bound once on the node's GPUs (MultiGpuBackend, host activity, FP32 records) plus the fine
    tuning requests. forward(s) -> (token, pNG, pRG) of the surface model (instantaneous, every row);
    set_parameters(idx) once (the parameter facets); jv(st, k) -> (JNG, JRG) (k, rows) for the k directions written
    to direction_buffer(k) (parameter space, k x npar); jtu(st, uN, uR) -> g (parameter space, npar; a shared buffer
    valid until the next call); seen(row mask) -> bool (facets). The activity buffer holds one exponent field at a
    time; a state whose field is not the buffer's is written again (token) and the workers re-upload A = exp(s).
    Speed-up 2026-11-08 ("do the ~2x speed-up first", user): directions and J^T u travel in parameter space (no host
    scatter / gather over the 12.4 M facets; the workers map their device facets to parameters on the GPU) and the sum
    of the per-device J^T u arrays is done by the workers in parallel slices (same order, device 0..7, as the parent's
    sum before); results are bitwise those of the first implementation, kept as jv_legacy / jtu_legacy (facet space)
    for the comparison in validate() check 7."""

    def __init__(self, ctx, devices, records_per_block=4096):
        from multiprocessing import shared_memory
        from cops_cuda import backend as bk
        from cops_cuda import gpu_model, multi_gpu, precision, reference
        reference.bind()
        os.environ[precision.RECORD_PRECISION_ENVIRONMENT] = "fp32"          # --precision mixed (workers inherit)
        Model = gpu_model.gpu_model_class()
        M = Model(str(ctx.run), ctx.header, ctx.cols, ctx.table_rows, ctx.el, ctx.coef)
        if not np.array_equal(np.asarray(M.rows), np.arange(ctx.rows.n)):
            raise SystemExit("the model rows are not the base rows in order")
        if not np.array_equal(np.asarray(M.ipoint), ctx.rows.ipoint):
            raise SystemExit("the model's point ids differ from the fit's")
        self.nf, self.nrows = len(ctx.area), ctx.rows.n
        t0 = time.time()
        self.backend = multi_gpu.MultiGpuBackend(M._spec, M, None, 1.0, devices=list(devices), mode="tree",
                                                 records_per_block=records_per_block, transport="local",
                                                 activity="host")
        self.bind_s = time.time() - t0
        log("bound %s in %.1f s" % (self.backend.summary(), self.bind_s))
        self.shm = []
        self.v_shm = shared_memory.SharedMemory(create=True, size=KMAX * self.nf * 8)
        self.shm.append(self.v_shm)
        self.v = np.ndarray((KMAX, self.nf), dtype=np.float64, buffer=self.v_shm.buf)
        self.g_shm, self.g = [], []
        for _ in self.backend.workers:
            s = shared_memory.SharedMemory(create=True, size=self.nf * 8)
            self.shm.append(s)
            self.g_shm.append(s)
            self.g.append(np.ndarray((self.nf,), dtype=np.float64, buffer=s.buf))
        self.state = bk.ParameterState(ctx.a, ctx.T, ctx.kappa, None)
        self.counter = 0
        self.shm_token = None
        self.reuploads = 0
        self.worker_seconds = {"ft_eval": [], "ft_jv": [], "ft_jtu": []}
        self.attach = []
        for w, g in zip(self.backend.workers, self.g_shm):
            self._send(w, ("ft_attach", {"nfaces": self.nf, "v_name": self.v_shm.name, "g_name": g.name,
                                         "kmax": KMAX, "kernel": str(KERNEL_FILE)}), "ft_attach")
        for w in self.backend.workers:
            self.attach.append(self._receive(w, "ft_attach")[1])
        for info in self.attach:
            log("worker device %s: %d points, %d records, %d device facets, streamed %s, kernel %s, GPU free %.1f of "
                "%.1f GB" % (info["device_id"], info["points"], info["records"], info["device_faces"],
                             info["streamed"], info["kernel_sha256"][:12], info["gpu_free_bytes"] / 1e9,
                             info["gpu_total_bytes"] / 1e9))

    def _send(self, w, message, what):
        self.backend._send(w, message, what)

    def _receive(self, w, what):
        reply = self.backend._receive(w, what)
        if not (isinstance(reply, tuple) and reply and reply[0] == "ok"):
            raise RuntimeError("unexpected reply to %s: %r" % (what, reply if not isinstance(reply, tuple)
                                                              else reply[:1]))
        return reply

    def _write(self, s):
        self.backend.act[:] = s
        self.counter += 1
        self.shm_token = self.counter
        return self.counter

    def _ensure(self, st):
        if st.token != self.shm_token:
            st.token = self._write(st.s)
            self.reuploads += 1
        return st.token

    def forward(self, s):
        token = self._write(s)
        for w in self.backend.workers:
            self._send(w, ("ft_eval", token, self.state), "ft_eval")
        pNG, pRG = np.empty(self.nrows), np.empty(self.nrows)
        sec = []
        for w in self.backend.workers:
            _, ng, rg, seconds = self._receive(w, "ft_eval")
            pNG[w["positions"]], pRG[w["positions"]] = ng, rg
            sec.append(seconds)
        self.worker_seconds["ft_eval"].append(max(sec))
        return token, pNG, pRG

    def set_parameters(self, idx):
        """The parameter facets (sorted ids): parameter-space shared buffers and the workers' facet maps."""
        from multiprocessing import shared_memory
        idx = np.ascontiguousarray(idx, dtype=np.int64)
        if len(idx) == 0 or np.any(np.diff(idx) <= 0):
            raise ValueError("parameter facets must be sorted and unique")
        if getattr(self, "npar", None) is not None:
            raise RuntimeError("set_parameters is called once")
        self.npar = n = len(idx)
        mk = lambda count: shared_memory.SharedMemory(create=True, size=max(8, 8 * count))  # noqa: E731
        self.idx_shm = mk(n)
        np.ndarray((n,), dtype=np.int64, buffer=self.idx_shm.buf)[:] = idx
        self.vpar_shm = mk(KMAX * n)
        self.vpar = np.ndarray((KMAX, n), dtype=np.float64, buffer=self.vpar_shm.buf)
        self.gpar_shm = [mk(n) for _ in self.backend.workers]
        self.out_shm = mk(n)
        self.out = np.ndarray((n,), dtype=np.float64, buffer=self.out_shm.buf)
        self.shm += [self.idx_shm, self.vpar_shm, self.out_shm] + self.gpar_shm
        nw = len(self.backend.workers)
        cuts = [n * k // nw for k in range(nw + 1)]
        self.slices = [(cuts[k], cuts[k + 1]) for k in range(nw)]
        spec = {"npar": n, "idx_name": self.idx_shm.name, "vpar_name": self.vpar_shm.name,
                "gpar_names": [s.name for s in self.gpar_shm], "out_name": self.out_shm.name, "kmax": KMAX}
        for k, w in enumerate(self.backend.workers):
            self._send(w, ("ft_params", dict(spec, me=k)), "ft_params")
        self.param_info = [self._receive(w, "ft_params")[1] for w in self.backend.workers]
        for info in self.param_info:
            log("worker device %s: %d of %d device facets are parameters; facets unique within every point: %s; "
                "page-locked buffers %d%s (%.2f s)" % (
                    info["device_id"], info["param_faces"], info["device_faces"], info["unique_per_point"],
                    info["pinned_buffers"], " (%s)" % info["pin_error"] if info["pin_error"] else "", info["seconds"]))
        return self.param_info

    def direction_buffer(self, k):
        if not 1 <= k <= KMAX:
            raise ValueError("1 <= k <= %d directions per pass" % KMAX)
        return self.vpar[:k]

    def jv(self, st, k):
        token = self._ensure(st)
        for w in self.backend.workers:
            self._send(w, ("ft_jv2", token, self.state, int(k)), "ft_jv2")
        JNG, JRG = np.empty((k, self.nrows)), np.empty((k, self.nrows))
        sec = []
        for w in self.backend.workers:
            _, jn, jr, seconds = self._receive(w, "ft_jv2")
            JNG[:, w["positions"]], JRG[:, w["positions"]] = jn, jr
            sec.append(seconds)
        self.worker_seconds["ft_jv"].append(max(sec))
        return JNG, JRG

    def jtu(self, st, uN, uR):
        token = self._ensure(st)
        for w in self.backend.workers:
            pos = w["positions"]
            self._send(w, ("ft_jtu2", token, self.state, np.ascontiguousarray(uN[pos]),
                           np.ascontiguousarray(uR[pos])), "ft_jtu2")
        sec = [self._receive(w, "ft_jtu2")[1] for w in self.backend.workers]
        for w, (lo, hi) in zip(self.backend.workers, self.slices):
            self._send(w, ("ft_reduce", lo, hi), "ft_reduce")
        red = [self._receive(w, "ft_reduce")[1] for w in self.backend.workers]
        self.worker_seconds["ft_jtu"].append(max(sec))
        self.worker_seconds.setdefault("ft_reduce", []).append(max(red))
        return self.out

    # --- the first implementation (facet space; host scatter, gather and sum), kept for validate() check 7
    def direction_buffer_legacy(self, k):
        if not 1 <= k <= KMAX:
            raise ValueError("1 <= k <= %d directions per pass" % KMAX)
        return self.v[:k]

    def jv_legacy(self, st, k):
        token = self._ensure(st)
        for w in self.backend.workers:
            self._send(w, ("ft_jv", token, self.state, int(k)), "ft_jv")
        JNG, JRG = np.empty((k, self.nrows)), np.empty((k, self.nrows))
        for w in self.backend.workers:
            _, jn, jr, seconds = self._receive(w, "ft_jv")
            JNG[:, w["positions"]], JRG[:, w["positions"]] = jn, jr
        return JNG, JRG

    def jtu_legacy(self, st, uN, uR):
        token = self._ensure(st)
        for w in self.backend.workers:
            pos = w["positions"]
            self._send(w, ("ft_jtu", token, self.state, np.ascontiguousarray(uN[pos]),
                           np.ascontiguousarray(uR[pos])), "ft_jtu")
        for w in self.backend.workers:
            self._receive(w, "ft_jtu")
        g = self.g[0].copy()
        for k in range(1, len(self.g)):                   # device order: deterministic
            g += self.g[k]
        return g

    def seen(self, row_mask, structural=False):
        """Facets with a record at a selected row: fit_cops' seen predicate (Kn > 0 | KPhi > 0 | KV > 0), or with
        structural=True the old late-window fine tuning's (Kn != 0 | KPhi != 0 | KV != 0)."""
        row_mask = np.asarray(row_mask, bool)
        for w in self.backend.workers:
            self._send(w, ("ft_seen", np.ascontiguousarray(row_mask[w["positions"]]), bool(structural)), "ft_seen")
        out = np.zeros(self.nf, bool)
        for k, w in enumerate(self.backend.workers):
            self._receive(w, "ft_seen")
            out |= self.g[k] != 0.0
        return out

    def stats(self):
        for w in self.backend.workers:
            self._send(w, ("ft_stats",), "ft_stats")
        return [self._receive(w, "ft_stats")[1] for w in self.backend.workers]

    def identity(self):
        try:
            return jsonable(self.backend.identity())
        except Exception as exc:                          # recorded, never silent
            return {"error": "%s: %s" % (exc.__class__.__name__, exc)}

    def close(self):
        try:
            self.backend.close()
        finally:
            for name in ("v", "g", "vpar", "out"):
                if hasattr(self, name):
                    delattr(self, name)
            for s in self.shm:
                try:
                    s.close()
                    s.unlink()
                except Exception:
                    pass
            self.shm = []


# ------------------------------------------------------------------------------------------------ worker side
class _WorkerFineTune:
    """GPU worker state of the fine tuning (lives in the spawned process next to its CudaBackend)."""

    def __init__(self, worker, spec):
        from multiprocessing import shared_memory
        from cops_cuda.errors import UnsupportedConfiguration
        from cops_cuda.memory import MemoryPlan
        import cupy
        be = worker.backend
        if be is None:
            raise RuntimeError("no backend bound")
        if not be.mixed:
            raise UnsupportedConfiguration("positive_finetune: FP32 records (--precision mixed) are required")
        ds = be.dataset
        if ds.streamed:
            raise UnsupportedConfiguration("positive_finetune: streamed records are not supported (%d of %d "
                                           "resident); use more GPUs" % (ds.nresident, ds.nrecords))
        if be.activity_mode != "host" or be.activity_source is None:
            raise UnsupportedConfiguration("positive_finetune: host activity from the driver is required")
        if be.fused is not True or be.mode != "tree":
            raise UnsupportedConfiguration("positive_finetune: fused tree mode required")
        self.cupy, self.be, self.ds = cupy, be, ds
        self.nf = int(spec["nfaces"])
        self.kmax = int(spec["kmax"])
        if ds.inputs.nfaces != self.nf:
            raise RuntimeError("mesh size %d differs from the driver's %d" % (ds.inputs.nfaces, self.nf))
        self.nfd, self.npt = int(ds.nfaces_device), int(ds.npoints)
        self.local = ds.local_faces
        self._shm = [shared_memory.SharedMemory(name=spec["v_name"]), shared_memory.SharedMemory(name=spec["g_name"])]
        self.v = np.ndarray((self.kmax, self.nf), dtype=np.float64, buffer=self._shm[0].buf)
        self.g = np.ndarray((self.nf,), dtype=np.float64, buffer=self._shm[1].buf)
        path = Path(spec["kernel"])
        source = path.read_bytes()
        self.kernel_identity = {"path": str(path), "sha256": hashlib.sha256(source).hexdigest(),
                                "options": list(NVRTC_OPTIONS)}
        c = be.constants
        self.cNbmH = c.cNG["H2O"] * c.bm["H2O"]
        self.cRG_H = np.asarray(c.cRG["H2O"], dtype=np.float64)
        if self.cRG_H.shape != (self.npt,):
            raise RuntimeError("cRG_H2O has shape %s for %d points" % (self.cRG_H.shape, self.npt))
        self.first_block = np.asarray(be._first_block, dtype=np.int64)
        self.nblocks = int(self.first_block[-1])
        with cupy.cuda.Device(be.device_id):
            self.module = cupy.RawModule(code=source.decode("utf-8"), options=NVRTC_OPTIONS)
            self.module.compile()
            self.k_jtu = self.module.get_function("ft_jtu_m")
            self.k_seen = self.module.get_function("ft_seen_m")
            self.k_dup = self.module.get_function("ft_dup_m")
            self.k_map = self.module.get_function("ft_map")
            plan = MemoryPlan("positive_finetune worker arrays")
            plan.add_array("ft.g", "scratch", (self.nfd,), np.float64, "J^T u of the device facets")
            plan.add_array("ft.seen", "scratch", (self.nfd,), np.uint8, "seen facets")
            plan.add_array("ft.mask", "scratch", (max(1, self.npt),), np.uint8, "point mask")
            plan.add_array("ft.bad", "scratch", (1,), np.uint64, "launch check")
            self.reservation = be.budget.reserve(plan)
            with be.budget.allocation_guard(plan):
                self.d_g = cupy.zeros(self.nfd, dtype=np.float64)
                self.d_seen = cupy.zeros(self.nfd, dtype=np.uint8)
                self.d_mask = cupy.zeros(max(1, self.npt), dtype=np.uint8)
                self.d_bad = cupy.zeros(1, dtype=np.uint64)
            free, total = cupy.cuda.runtime.memGetInfo()
        self.gpu = {"free": int(free), "total": int(total)}
        self.token = None
        self.uploads = 0
        self.calls = {"ft_eval": 0, "ft_jv": 0, "ft_jtu": 0, "ft_seen": 0, "ft_jv2": 0, "ft_jtu2": 0, "ft_reduce": 0}
        self.seconds = {"ft_eval": 0.0, "ft_jv": 0.0, "ft_jtu": 0.0, "upload": 0.0, "ft_jv2": 0.0, "ft_jtu2": 0.0,
                        "ft_jtu2_kernel": 0.0, "ft_reduce": 0.0}
        self.npar = None
        self.unique = None
        self.jtu_mode = None                     # "atomic" | "rmw", chosen at the first ft_jtu2 (calibration)
        self.calibration = None

    def info(self):
        ds = self.ds
        return {"device_id": self.be.device_id, "pid": os.getpid(), "host": socket.gethostname(),
                "points": self.npt, "records": int(ds.nrecords), "device_faces": self.nfd,
                "streamed": bool(ds.streamed), "blocks": self.nblocks, "kernel_sha256": self.kernel_identity["sha256"],
                "kernel": self.kernel_identity, "gpu_free_bytes": self.gpu["free"], "gpu_total_bytes": self.gpu["total"],
                "record_precision": self.be.record_precision}

    def _activity(self, token, state):
        if token != self.token:
            t0 = time.perf_counter()
            with self.cupy.cuda.Device(self.be.device_id):
                self.be._upload_activity(state)
            self.token = token
            self.uploads += 1
            self.seconds["upload"] += time.perf_counter() - t0

    def evaluate(self, token, state):
        from cops_cuda import backend as bk
        t0 = time.perf_counter()
        result = self.be.evaluate(state, bk.PRESSURES_ONLY)        # uploads A = exp(s) from the shared buffer
        self.token = token
        dt_ = time.perf_counter() - t0
        self.calls["ft_eval"] += 1
        self.seconds["ft_eval"] += dt_
        return result.pNG, result.pRG, dt_

    def jv(self, token, state, k):
        t0 = time.perf_counter()
        self._activity(token, state)
        be, cupy, nfd = self.be, self.cupy, self.nfd
        if not 1 <= k <= min(self.kmax, 16):
            raise ValueError("k = %d" % k)
        with cupy.cuda.Device(be.device_id):
            cb = be.arrays["colbasis"]
            if cb.size < k * nfd:
                raise RuntimeError("colbasis too small")
            for j in range(k):
                col = self.v[j] if self.local is None else np.take(self.v[j], self.local)
                cb[j * nfd:(j + 1) * nfd].set(np.ascontiguousarray(col, dtype=np.float64))
            be._columns_cache = None
            codes = [1000 + j for j in range(k)] + [2000 + j for j in range(k)]
            S = be._run_fused([("eval_fused_c%d" % k, 0, codes)], codes, be._params(state))
        JNG = self.cNbmH * S[:k]                 # CudaBackend._evaluate: cNG_H * bm_H * colNG
        JRG = self.cRG_H[None, :] * S[k:]        # cRG_H[:, None] * colRG
        dt_ = time.perf_counter() - t0
        self.calls["ft_jv"] += 1
        self.seconds["ft_jv"] += dt_
        return JNG, JRG, dt_

    def _coefficients(self, uN, uR):
        uN, uR = np.asarray(uN, np.float64), np.asarray(uR, np.float64)
        if uN.shape != (self.npt,) or uR.shape != (self.npt,):
            raise ValueError("u of shape %s/%s for %d points" % (uN.shape, uR.shape, self.npt))
        return self.cNbmH * uN, self.cRG_H * uR

    def _jtu_kernel(self, state, cN, cR, use_atomic):
        """J^T u of the device facets into d_g (one launch per point, point order); returns the kernel seconds."""
        be, cupy, ds = self.be, self.cupy, self.ds
        order = ("Tmax", "T75", "Tmin", "KT", "cNbmH", "kcNbmC", "bmH", "bmC")
        params = be._params(state)
        a, d = be.arrays, ds.arrays
        fb = self.first_block
        with cupy.cuda.Device(be.device_id):
            P = cupy.asarray(np.array([params[k] for k in order], dtype=np.float64))
            self.d_g.fill(0.0)
            self.d_bad.fill(0)
            blocks = a["blocks"]
            cupy.cuda.get_current_stream().synchronize()
            t0 = time.perf_counter()
            for p in range(self.npt):
                if cN[p] == 0.0 and cR[p] == 0.0:        # adding exact zeros: skipped (same sums)
                    continue
                b0, b1 = int(fb[p]), int(fb[p + 1])
                if b1 <= b0:
                    continue
                self.k_jtu((b1 - b0,), (BLOCK,),
                           (blocks[3 * b0:3 * b1], d["kernel.face"], d["kernel.kn"], d["kernel.kphi"],
                            d["kernel.kv"], d["kernel.cos"], a["pt"], np.int64(self.npt), d["mesh.ACO2"],
                            np.float64(be.constants.COS75), a["af_face"], a["active"], np.int64(self.nfd), P,
                            np.int64(p), np.float64(cN[p]), np.float64(cR[p]), self.d_g, self.d_bad,
                            np.int32(1 if use_atomic else 0)))
            bad = int(self.d_bad.get()[0])                # synchronizes
            seconds = time.perf_counter() - t0
            del P
        if bad:
            raise RuntimeError("ft_jtu_m: %d blocks with a point other than the launched one" % bad)
        return seconds

    def jtu(self, token, state, uN, uR):
        """The first implementation (facet space; atomicAdd), kept for validate() check 7."""
        t0 = time.perf_counter()
        self._activity(token, state)
        cN, cR = self._coefficients(uN, uR)
        self._jtu_kernel(state, cN, cR, True)
        with self.cupy.cuda.Device(self.be.device_id):
            g = self.d_g.get()
        if self.local is None:
            self.g[:] = g
        else:
            self.g.fill(0.0)
            self.g[self.local] = g
        dt_ = time.perf_counter() - t0
        self.calls["ft_jtu"] += 1
        self.seconds["ft_jtu"] += dt_
        return (dt_,)

    # --- parameter space (speed-up 2026-11-08)
    def set_parameters(self, spec):
        from multiprocessing import shared_memory
        from cops_cuda.memory import MemoryPlan
        cupy, be = self.cupy, self.be
        t0 = time.perf_counter()
        n = int(spec["npar"])
        s = shared_memory.SharedMemory(name=spec["idx_name"])
        view = np.ndarray((n,), dtype=np.int64, buffer=s.buf)
        idx = view.copy()
        del view
        s.close()
        local = np.arange(self.nf, dtype=np.int64) if self.local is None else np.asarray(self.local, np.int64)
        pos = np.minimum(np.searchsorted(idx, local), n - 1)
        valid = idx[pos] == local
        lsel = np.flatnonzero(valid).astype(np.int64)          # device facet positions that are parameters
        lpar = pos[valid].astype(np.int64)                      # their parameter index
        self.npar, self.nsel = n, len(lsel)
        self._pshm = [shared_memory.SharedMemory(name=spec["vpar_name"]),
                      shared_memory.SharedMemory(name=spec["out_name"])] + \
            [shared_memory.SharedMemory(name=x) for x in spec["gpar_names"]]
        self.vpar = np.ndarray((int(spec["kmax"]), n), dtype=np.float64, buffer=self._pshm[0].buf)
        self.out = np.ndarray((n,), dtype=np.float64, buffer=self._pshm[1].buf)
        self.gpar = [np.ndarray((n,), dtype=np.float64, buffer=x.buf) for x in self._pshm[2:]]
        self.gmine = self.gpar[int(spec["me"])]
        # page-lock this worker's shared buffers (direction in, J^T u out): direct DMA instead of pageable staging
        # (copies only; the values are unchanged); a failure (e.g. a memlock limit) keeps the pageable copies
        self.pinned, self.pin_error = [], None
        with cupy.cuda.Device(be.device_id):
            for arr in (self.vpar, self.gmine):
                try:
                    cupy.cuda.runtime.hostRegister(arr.ctypes.data, arr.nbytes, 0)
                    self.pinned.append(int(arr.ctypes.data))
                except Exception as exc:                  # recorded in the stats
                    self.pin_error = "%s: %s" % (exc.__class__.__name__, exc)
        with cupy.cuda.Device(be.device_id):
            plan = MemoryPlan("positive_finetune worker parameter arrays")
            plan.add_array("ft.vpar", "scratch", (n,), np.float64, "direction (parameter space)")
            plan.add_array("ft.gpar", "scratch", (n,), np.float64, "J^T u (parameter space)")
            plan.add_array("ft.lsel", "static", (max(1, self.nsel),), np.int64, "device facet -> parameter map")
            plan.add_array("ft.lpar", "static", (max(1, self.nsel),), np.int64, "device facet -> parameter map")
            plan.add_array("ft.last", "scratch", (self.nfd,), np.int32, "duplicate check")
            plan.add_array("ft.dup", "scratch", (1,), np.uint64, "duplicate check")
            self.reservation_par = be.budget.reserve(plan)
            with be.budget.allocation_guard(plan):
                self.d_vpar = cupy.zeros(n, dtype=np.float64)
                self.d_gpar = cupy.zeros(n, dtype=np.float64)
                self.d_lsel = cupy.asarray(lsel if len(lsel) else np.zeros(1, np.int64))
                self.d_lpar = cupy.asarray(lpar if len(lpar) else np.zeros(1, np.int64))
                d_last = cupy.full(self.nfd, -1, dtype=np.int32)
                d_dup = cupy.zeros(1, dtype=np.uint64)
            # facets listed twice within one point (the plain read-modify-write J^T u needs none)
            blocks, face, fb = be.arrays["blocks"], self.ds.arrays["kernel.face"], self.first_block
            for p in range(self.npt):
                b0, b1 = int(fb[p]), int(fb[p + 1])
                if b1 > b0:
                    self.k_dup((b1 - b0,), (BLOCK,), (blocks[3 * b0:3 * b1], face, np.int64(p), d_last, d_dup))
            dup = int(d_dup.get()[0])
            del d_last, d_dup
        self.unique = dup == 0
        self.duplicates = dup
        return ({"device_id": be.device_id, "param_faces": self.nsel, "device_faces": self.nfd,
                 "unique_per_point": self.unique, "duplicates_or_bad_blocks": dup,
                 "pinned_buffers": len(self.pinned), "pin_error": self.pin_error,
                 "seconds": time.perf_counter() - t0},)

    def _map(self, n, dst_idx, dst, src_idx, src):
        if n:
            self.k_map((max(1, min(65535, -(-n // 256))),), (256,), (np.int64(n), dst_idx, dst, src_idx, src))

    def jv2(self, token, state, k):
        t0 = time.perf_counter()
        self._activity(token, state)
        be, cupy, nfd = self.be, self.cupy, self.nfd
        if not 1 <= k <= min(self.kmax, 16):
            raise ValueError("k = %d" % k)
        with cupy.cuda.Device(be.device_id):
            cb = be.arrays["colbasis"]
            for j in range(k):
                self.d_vpar.set(self.vpar[j])
                seg = cb[j * nfd:(j + 1) * nfd]
                seg.fill(0.0)                                 # device facets that are not parameters
                self._map(self.nsel, self.d_lsel, seg, self.d_lpar, self.d_vpar)
            be._columns_cache = None
            codes = [1000 + j for j in range(k)] + [2000 + j for j in range(k)]
            S = be._run_fused([("eval_fused_c%d" % k, 0, codes)], codes, be._params(state))
        JNG = self.cNbmH * S[:k]
        JRG = self.cRG_H[None, :] * S[k:]
        dt_ = time.perf_counter() - t0
        self.calls["ft_jv2"] += 1
        self.seconds["ft_jv2"] += dt_
        return JNG, JRG, dt_

    def jtu2(self, token, state, uN, uR):
        t0 = time.perf_counter()
        self._activity(token, state)
        cN, cR = self._coefficients(uN, uR)
        cupy, be = self.cupy, self.be
        if self.jtu_mode is None:
            # calibration (first call): atomicAdd vs plain read-modify-write; same additions in the same order, so
            # the results must be bitwise equal; the faster one is used from now on (rmw only if no facet repeats
            # within a point)
            ta = self._jtu_kernel(state, cN, cR, True)
            with cupy.cuda.Device(be.device_id):
                ga = self.d_g.copy()
            cal = {"atomic_s": ta, "unique_per_point": self.unique}
            mode = "atomic"
            if self.unique:
                tr = self._jtu_kernel(state, cN, cR, False)
                with cupy.cuda.Device(be.device_id):
                    same = bool(cupy.array_equal(ga, self.d_g))
                del ga
                cal.update(rmw_s=tr, bitwise_equal=same)
                if same and tr < ta:
                    mode = "rmw"
                if not same:
                    self._jtu_kernel(state, cN, cR, True)      # keep the atomic result
            else:
                del ga
            cal["mode"] = mode
            self.jtu_mode, self.calibration = mode, cal
            kt = 0.0
        else:
            kt = self._jtu_kernel(state, cN, cR, self.jtu_mode == "atomic")
        with cupy.cuda.Device(be.device_id):
            self.d_gpar.fill(0.0)
            self._map(self.nsel, self.d_lpar, self.d_gpar, self.d_lsel, self.d_g)
            self.d_gpar.get(out=self.gmine)
        dt_ = time.perf_counter() - t0
        self.calls["ft_jtu2"] += 1
        self.seconds["ft_jtu2"] += dt_
        self.seconds["ft_jtu2_kernel"] += kt
        return (dt_,)

    def reduce(self, lo, hi):
        """out[lo:hi] = gpar_0 + gpar_1 + ... in device order (the order of the first implementation's host sum)."""
        t0 = time.perf_counter()
        s = self.gpar[0][lo:hi].copy()
        for k in range(1, len(self.gpar)):
            s += self.gpar[k][lo:hi]
        self.out[lo:hi] = s
        dt_ = time.perf_counter() - t0
        self.calls["ft_reduce"] += 1
        self.seconds["ft_reduce"] += dt_
        return (dt_,)

    def seen(self, mask, structural=False):
        be, cupy, ds = self.be, self.cupy, self.ds
        mask = np.asarray(mask, bool)
        if mask.shape != (self.npt,):
            raise ValueError("mask of shape %s for %d points" % (mask.shape, self.npt))
        with cupy.cuda.Device(be.device_id):
            self.d_mask[:self.npt].set(mask.astype(np.uint8))
            self.d_seen.fill(0)
            if self.nblocks:
                d = ds.arrays
                self.k_seen((min(self.nblocks, 1 << 20),), (BLOCK,),
                            (be.arrays["blocks"], np.int64(self.nblocks), d["kernel.face"], d["kernel.kn"],
                             d["kernel.kphi"], d["kernel.kv"], self.d_mask, self.d_seen,
                             np.int32(1 if structural else 0)))
            s = self.d_seen.get()
        if self.local is None:
            self.g[:] = s
        else:
            self.g.fill(0.0)
            self.g[self.local] = s
        self.calls["ft_seen"] += 1
        return (int(s.sum()),)

    def stats(self):
        cupy = self.cupy
        with cupy.cuda.Device(self.be.device_id):
            free, total = cupy.cuda.runtime.memGetInfo()
            pool = cupy.get_default_memory_pool()
            used, held = int(pool.used_bytes()), int(pool.total_bytes())
        return ({"device_id": self.be.device_id, "pid": os.getpid(), "gpu_free_bytes": int(free),
                 "gpu_total_bytes": int(total), "gpu_used_bytes": int(total - free), "pool_used_bytes": used,
                 "pool_held_bytes": held, "worker_max_rss_gb": max_rss_gb(), "activity_uploads": self.uploads,
                 "calls": dict(self.calls), "seconds": dict(self.seconds), "jtu_mode": self.jtu_mode,
                 "jtu_calibration": self.calibration, "unique_per_point": self.unique,
                 "pinned_buffers": len(getattr(self, "pinned", [])), "pin_error": getattr(self, "pin_error", None),
                 "param_faces": getattr(self, "nsel", None), "backend_timing": jsonable(self.be.timing)},)


def _ft_dispatch(worker, message):
    kind = message[0]
    if kind == "ft_attach":
        ft = _WorkerFineTune(worker, message[1])
        worker._positive_finetune = ft
        return (ft.info(),)
    ft = getattr(worker, "_positive_finetune", None)
    if ft is None or ft.be is not worker.backend:
        raise RuntimeError("positive_finetune: ft_attach is needed first (or the backend was rebound)")
    if kind == "ft_eval":
        return ft.evaluate(message[1], message[2])
    if kind == "ft_jv":
        return ft.jv(message[1], message[2], message[3])
    if kind == "ft_jtu":
        return ft.jtu(message[1], message[2], message[3], message[4])
    if kind == "ft_seen":
        return ft.seen(message[1], bool(message[2]) if len(message) > 2 else False)
    if kind == "ft_params":
        return ft.set_parameters(message[1])
    if kind == "ft_jv2":
        return ft.jv2(message[1], message[2], message[3])
    if kind == "ft_jtu2":
        return ft.jtu2(message[1], message[2], message[3], message[4])
    if kind == "ft_reduce":
        return ft.reduce(message[1], message[2])
    if kind == "ft_stats":
        return ft.stats()
    raise ValueError("unknown request %r" % (kind,))


def install_worker_handlers():
    """Add the ft_* requests to cops_cuda.transport.WorkerState.handle in this process (idempotent)."""
    from cops_cuda import transport
    current = transport.WorkerState.handle
    if getattr(current, "_positive_finetune", False):
        return False
    original = current

    def handle(self, message):
        kind = message[0] if message else None
        if isinstance(kind, str) and kind.startswith("ft_"):
            try:
                return ("ok",) + tuple(_ft_dispatch(self, message))
            except Exception:
                return ("error", traceback.format_exc())
        return original(self, message)

    handle._positive_finetune = True
    transport.WorkerState.handle = handle
    return True


def _process_setup():
    """Driver and spawned GPU workers: the GPU patch from the inherited COPS_POS_* environment, then the handlers."""
    import positive_activity
    if os.environ.get("COPS_POS_ENABLE") == "1" and not getattr(positive_activity, "_installed_by_positive_finetune",
                                                                False):
        positive_activity.install()
        positive_activity._installed_by_positive_finetune = True
    install_worker_handlers()


if __name__ in ("__main__", "__mp_main__"):
    _process_setup()


# ------------------------------------------------------------------------------------------------ the tuned fit
class FitContext:
    """Everything the fine tuning takes from the tuned fit and its base (CPU only)."""

    def __init__(self, fit_dir):
        from cops_cuda import reference
        import grain_crossing as gc
        import positive_problem as pp
        reference.bind()
        fit_cops = reference.load("fit_cops")
        if float(fit_cops.FLOOR) != FLOOR:
            raise SystemExit("fit_cops.FLOOR %r != %r" % (fit_cops.FLOOR, FLOOR))
        self.fit_dir = Path(fit_dir).resolve()
        f = {k: self.fit_dir / k for k in ("fit_restart.npz", "fit_manifest.json", "activity_exp_facets.npz",
                                           "grain_term.json", "COPS_model.dat", "late_window_model.dat",
                                           "fit_result.txt")}
        missing = [str(p) for p in f.values() if not p.is_file()]
        if missing:
            raise SystemExit("the fit directory lacks %s" % missing)
        self.files = f
        self.identity = {k: file_identity(p) for k, p in f.items() if k != "late_window_model.dat"}
        # --- restart
        with np.load(f["fit_restart.npz"], allow_pickle=False) as z:
            self.restart = {k: z[k] for k in z.files if k != "meta_json"}
            self.meta = json.loads(str(z["meta_json"]))
        ext = self.meta.get("model_extension") or {}
        self.fit_parameterization = ext.get("parameterization", "")
        if not str(self.fit_parameterization).startswith(FIT_PARAMETERIZATION_PREFIX):
            raise SystemExit("%s is not a positive_fit.py fit (restart model_extension.parameterization %r)"
                             % (self.fit_dir, self.fit_parameterization))
        if self.meta.get("status") != "complete":
            raise SystemExit("the fit's restart has status %r (a finished fit is required)" % self.meta.get("status"))
        st = self.meta["settings"]
        if not st.get("fix_temperature") or st.get("fix_kappa") or st.get("ordered_temperature"):
            raise SystemExit("unexpected fit settings %s" % st)
        names = list(self.meta["names"])
        if tuple(names[-2:]) != GRAIN_NAMES:
            raise SystemExit("the fit's last parameters are %s, not the grain" % names[-2:])
        self.names = names
        x = np.asarray(self.restart["x"], float)
        self.x = x
        self.xA, self.xB = float(x[-2]), float(x[-1])
        self.a = np.asarray(self.restart["a"], float)
        self.T = np.asarray(self.restart["T"], float)
        self.kappa = float(np.asarray(self.restart["kappa"]).reshape(-1)[0])
        self.lam = float(st["lambda"])
        if len(np.asarray(self.restart["harmonic_values"])):
            raise SystemExit("harmonic enhancement columns are not part of the positive fits")
        # --- manifest
        self.manifest = json.loads(f["fit_manifest.json"].read_text())
        mext = self.manifest["model_extension"]
        pos = mext["positive"]
        self.positive_config = pos["config"]
        self.environment = {k: str(v) for k, v in pos["environment"].items()}
        self.environment_in_process = pos.get("environment_in_process")
        self.fit_positive_activity = pos.get("positive_activity", {})
        self.crossing_weight = float(mext["crossing_weight"])
        res = self.manifest["optimizer_results"]
        if len(res) != 1:
            raise SystemExit("the fit manifest has %d optimizer results (one lambda per fit expected)" % len(res))
        self.fit_costs = {"cost": float(res[0]["cost"]), "data_cost": float(res[0]["data_cost"]),
                          "regularization_cost": float(res[0]["regularization_cost"]),
                          "crossing_penalty_cost": float(mext.get("final", {}).get("penalty_cost", float("nan"))),
                          "lambda": float(res[0]["lambda"])}
        if self.fit_costs["lambda"] != self.lam:
            raise SystemExit("lambda: restart %r, manifest %r" % (self.lam, self.fit_costs["lambda"]))
        self.fit_result_txt = {}
        for line in f["fit_result.txt"].read_text().splitlines():
            v = line.split()
            if len(v) >= 2 and v[0] in ("final_cost", "final_data_cost", "final_regularization_cost", "lambda"):
                self.fit_result_txt[v[0]] = float(v[1])
        # --- per-facet fields of the fit
        with np.load(f["activity_exp_facets.npz"], allow_pickle=False) as z:
            self.npz_keys = list(z.files)
            self.s_npz = np.asarray(z["s"], float)
            self.s0 = np.asarray(z["s0"], float)
            self.seen_fit = np.asarray(z["seen"], bool)
            self.held_prev = np.asarray(z["held"], float) if "held" in z.files else None
            self.npz_meta = json.loads(str(z["meta_json"])) if "meta_json" in z.files else {}
        self.held_present = self.held_prev is not None
        # --- grain term of the fit
        gt = json.loads(f["grain_term.json"].read_text())
        self.grain_term = {"t_end_et": float(gt["t_end_et"]), "last_ng_before_gap_et": float(gt["last_ng_before_gap_et"]),
                           "A_Pa": float(gt["result"]["A_Pa"]), "B_Pa_per_s": float(gt["result"]["B_Pa_per_s"]),
                           "crossing_weight": float(gt.get("crossing_weight", float("nan")))}
        # --- the base
        self.run = Path(self.positive_config["run"]).resolve()
        self.point_table = Path(self.positive_config["point_table"]).resolve()
        log("fit %s: lambda %g, base %s, point table %s" % (self.fit_dir, self.lam, self.run, self.point_table))
        header, names_t, table_rows = fit_cops.read_model_file(str(self.run / "COPS_model.dat"))
        self.header, self.table_rows = header, table_rows
        self.cols = {n: i for i, n in enumerate(names_t)}
        t0 = time.time()
        self.el = fit_cops.read_elements(str(self.run / "elements.bin"))
        log("elements: %d facets (%.1f s)" % (len(self.el["area"]), time.time() - t0))
        self.area = np.asarray(self.el["area"], float)
        self.A_orig = np.asarray(self.el["AH2O"], float)
        self.nf = len(self.area)
        if not np.array_equal(self.area, np.asarray(self.restart["area"], float)):
            raise SystemExit("the restart's facet areas differ from the base's elements.bin")
        # the mesh identity the next cycle's fit compares (positive_problem: activity_templates.mesh_identity of the
        # base's centroids and areas; the restart records the same function's value)
        self.mesh_identity = reference.load("activity_templates").mesh_identity(self.el["centroid"], self.el["area"])
        if self.mesh_identity != (self.meta.get("mesh") or {}).get("identity"):
            raise SystemExit("mesh identity %s of the base differs from the fit restart's %s"
                             % (self.mesh_identity, (self.meta.get("mesh") or {}).get("identity")))
        srcdir = Path(fit_cops.__file__).resolve().parents[1]
        self.coef = fit_cops.read_activity_coefficients(str(srcdir), header)
        if not np.array_equal(np.asarray(self.coef, float), np.asarray(self.restart["coef_prior"], float)):
            raise SystemExit("the base's activity coefficients differ from the restart's coef_prior")
        for name, arr in (("s", self.s_npz), ("s0", self.s0), ("seen", self.seen_fit)):
            if arr.shape != (self.nf,):
                raise SystemExit("activity_exp_facets.npz %s has shape %s" % (name, arr.shape))
        if self.held_prev is None:
            self.held_prev = np.zeros(self.nf)
        elif self.held_prev.shape != (self.nf,) or not np.all(np.isfinite(self.held_prev)):
            raise SystemExit("activity_exp_facets.npz held has shape %s or non-finite values" % (self.held_prev.shape,))
        # --- the fit's exponent field: Y a + enhancement_scale * (T @ c), the evaluator's expression
        tc = np.ascontiguousarray(np.asarray(self.restart["template_columns"], float))
        tv = np.asarray(self.restart["template_values"], float)
        escale = float(self.meta.get("enhancement_scale", 1.0))
        s_param = self.el["Y"] @ self.a
        if tc.shape[1]:
            s_param = s_param + escale * (tc @ tv)
        act = np.asarray(self.restart["activity_linear"], float)
        self.exponent_checks = {
            "max_abs_s_param_minus_restart_activity_linear": float(np.max(np.abs(s_param - act))),
            "s_param_equals_restart_activity_linear_bitwise": bool(np.array_equal(s_param, act)),
            "max_abs_s_param_minus_npz_s": float(np.max(np.abs(s_param - self.s_npz))),
            "max_abs_s_param_plus_held_minus_npz_s": float(np.max(np.abs(s_param + self.held_prev - self.s_npz))),
            "held_present": self.held_present, "enhancement_scale": escale,
            "template_values": tv.tolist()}
        ok = min(self.exponent_checks["max_abs_s_param_minus_npz_s"],
                 self.exponent_checks["max_abs_s_param_plus_held_minus_npz_s"]) <= 1e-9
        if not ok:
            raise SystemExit("the fit's exponent field cannot be reproduced: %s" % self.exponent_checks)
        self.s_param = s_param
        self.s_start = s_param + self.held_prev if self.held_present else s_param.copy()
        # --- rows, selection, window, grain, crossing
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
        self.rows = Rows(ipoint, et, refNG, refRG, self.restart["use_ng"], self.restart["use_rg"], window,
                         self.crossing_weight, gc.GAP_SEARCH_S, (gc.A_UNIT, gc.B_UNIT))
        sel = self.meta.get("selection") or {}
        if (sel.get("used_ng") is not None and (self.rows.NNG, self.rows.NRG) != (sel["used_ng"], sel["used_rg"])):
            raise SystemExit("data selection: %d NG, %d RG; the fit used %s NG, %s RG"
                             % (self.rows.NNG, self.rows.NRG, sel.get("used_ng"), sel.get("used_rg")))
        if (self.rows.et_end, self.rows.et_before) != (self.grain_term["t_end_et"],
                                                         self.grain_term["last_ng_before_gap_et"]):
            raise SystemExit("grain onset %r / %r differs from the fit's grain_term.json %s"
                             % (self.rows.et_end, self.rows.et_before, self.grain_term))
        late0 = int(self.layout.late_rows[0])
        if self.rows.et[late0] != self.rows.et_end:
            raise SystemExit("the first late-window observation is not the grain onset t_end")
        self.used_rows = self.rows.useNG | self.rows.useRG

    def fit_tables(self):
        """The fit's per-row model from COPS_model.dat (%.6e) and late_window_model.dat (%.9e surface values)."""
        out = {"pNG": np.asarray(self.restart["pred_ng"], float), "pRG": np.asarray(self.restart["pred_rg"], float)}
        lw = {}
        for line in self.files["late_window_model.dat"].read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            lw[int(v[0])] = (v[3], float(v[5]), float(v[6]), float(v[7]))
        out["late_window"] = lw
        return out


def check_environment(ctx, floor_ratio):
    """This process's COPS_POS_* and the GPU patch configuration must be the tuned fit's; --floor-ratio must be its
    COPS_POS_FLOOR_RATIO. Also the kernel copies: pressure.cu must be a version finetune.cu is the transpose of
    (PRESSURE_COMPATIBLE_SHA256) and the one the fit ran with. COPS_POS_FLOOR_CO2 (2026-11-08) is one of the compared
    variables: a fit recorded without it (every fit before the CO2 floor, and the 0.02 / 0.002 fits) needs it unset,
    a ratio-1e-5 fit (recorded "1") needs 1; the configuration's co2_floor (absent in the earlier fits' records = off)
    must agree too."""
    import positive_activity
    proc = {k: os.environ.get(k) for k in ENV_KEYS}
    fit = {k: ctx.environment.get(k) for k in ENV_KEYS}
    if proc != fit:
        raise SystemExit("COPS_POS_* of this process %s differ from the tuned fit's %s (fit_manifest "
                         "model_extension.positive.environment); export the fit's values" % (proc, fit))
    if ctx.environment_in_process and {k: ctx.environment_in_process.get(k) for k in ENV_KEYS} != fit:
        raise SystemExit("the fit's recorded environment and environment_in_process differ: %s / %s"
                         % (fit, ctx.environment_in_process))
    if float(proc["COPS_POS_FLOOR_RATIO"]) != float(floor_ratio):
        raise SystemExit("--floor-ratio %r differs from COPS_POS_FLOOR_RATIO %r of the fit (user 2026-11-07: must "
                         "match)" % (floor_ratio, proc["COPS_POS_FLOOR_RATIO"]))
    if not positive_activity._STATE["installed"]:
        raise SystemExit("positive_activity is not installed (COPS_POS_ENABLE=1 is required)")
    cfg = positive_activity.config()
    fit_cfg = (ctx.fit_positive_activity or {}).get("config", {})
    keys = ("enabled", "modified", "night_ratio", "qmin_ratio", "activity", "floor_ratio")
    mine = {k: cfg.get(k) for k in keys}
    theirs = {k: fit_cfg.get(k) for k in keys}
    mine["co2_floor"] = bool(cfg.get("co2_floor", False))
    theirs["co2_floor"] = bool(fit_cfg.get("co2_floor", False))     # absent before 2026-11-08 = no CO2 floor
    if fit_cfg and mine != theirs:
        raise SystemExit("GPU patch configuration %s differs from the fit's %s" % (mine, theirs))
    if mine["activity"] != "exp":
        raise SystemExit("COPS_POS_ACTIVITY must be exp")
    current = {name: sha256(HERE / "kernels" / (name + ".cu")) for name in ("pressure", "activity")}
    fit_k = (ctx.fit_positive_activity or {}).get("kernel_sha256") or {}
    if current["pressure"] not in PRESSURE_COMPATIBLE_SHA256:
        raise SystemExit("kernels/pressure.cu has SHA-256 %s; kernels/finetune.cu (ft_jtu_m) is the transpose of %s: "
                         "re-derive the transpose kernel and re-validate" % (current["pressure"],
                                                                            sorted(PRESSURE_COMPATIBLE_SHA256)))
    if fit_k and fit_k != current:
        raise SystemExit("the fit ran with the kernel copies %s; the current ones are %s" % (fit_k, current))
    return {"process": proc, "fit": fit, "positive_activity_config": jsonable(cfg), "kernels": current,
            "fit_kernels": fit_k, "finetune_kernel": file_identity(KERNEL_FILE),
            "finetune_kernel_derived_from_pressure_sha256": PRESSURE_SOURCE_SHA256,
            "pressure_sha256_compatibility": PRESSURE_COMPATIBLE_SHA256[current["pressure"]],
            "co2_floor": mine["co2_floor"]}


FACET_RULES = {
    "late_only": "user 2026-11-08: only the facets seen exclusively in the last WINDOW_MIN minutes (the old late-window "
                 "fine tuning's set, ../LiouvilleCOPS/FINE_TUNING_STEP.md section 2): window [t_end - WINDOW_MIN, t_end], "
                 "t_end = the last selected data point (NG or RG used); a facet is a parameter iff it has a structural "
                 "record (Kn != 0 | KPhi != 0 | KV != 0) at an objective point inside the window and none at any "
                 "objective point before the window start; objective points = every used observation row (late "
                 "window observations included) and the 1-s grid rows with a nonzero W10x5 weight (a grid row before "
                 "the window start counts as earlier)",
    "all_seen": "user 2026-11-07: every facet seen (Kn > 0 | KPhi > 0 | KV > 0) from a point whose model enters the "
                "objective (used observation rows outside the late window and the grid rows with a W10x5 weight)"}


def late_only_facets(engine, rows, ctx, window_min, own, late_used):
    """(parameter mask, record, support info) of FACETS=late_only (see FACET_RULES); variants counted for the record."""
    used = rows.useNG | rows.useRG
    i_end = int(np.flatnonzero(used)[np.argmax(rows.et[used])])
    t_end = float(rows.et[i_end])
    t_start = t_end - 60.0 * float(window_min)
    points = (used & ~rows.is_grid) | (own & rows.is_grid)       # used observation rows + grid rows with W weight
    inside = points & (rows.et >= t_start) & (rows.et <= t_end)
    before = points & (rows.et < t_start)
    after = points & (rows.et > t_end)
    s_in = engine.seen(inside, structural=True)
    s_before = engine.seen(before, structural=True)
    s_after = engine.seen(after, structural=True)
    late = s_in & ~s_before
    # variants (recorded): objective points without the late observations' own records; every row of the base
    pts_a = own
    late_a = engine.seen(pts_a & (rows.et >= t_start) & (rows.et <= t_end), structural=True) & \
        ~engine.seen(pts_a & (rows.et < t_start), structural=True)
    late_c = engine.seen((rows.et >= t_start) & (rows.et <= t_end), structural=True) & \
        ~engine.seen(rows.et < t_start, structural=True)
    late_pos = engine.seen(inside) & ~engine.seen(before)          # fit_cops' "> 0" predicate
    area = ctx.area
    at = float(np.sum(area))
    utc = lambda i: ctx.utc[int(i)]                                 # noqa: E731
    first_in = int(np.flatnonzero(inside)[np.argmin(rows.et[inside])]) if inside.any() else None
    import positive_problem as pp
    sec_end = pp.parse_utc_seconds(utc(i_end))
    t_start_utc = (dt.datetime(2016, 1, 1) + dt.timedelta(seconds=sec_end - 60.0 * float(window_min))).isoformat(
        timespec="milliseconds")
    rec = {"window_min": float(window_min), "t_end_et": t_end, "t_start_et": t_start,
           "t_end_utc": utc(i_end), "t_start_utc": t_start_utc, "last_selected_row": i_end,
           "first_point_in_window_utc": utc(first_in) if first_in is not None else None,
           "points_inside": int(inside.sum()), "points_inside_observation": int((inside & ~rows.is_grid).sum()),
           "points_inside_grid": int((inside & rows.is_grid).sum()), "points_before": int(before.sum()),
           "points_after_window_end": int(after.sum()),
           "facets": int(late.sum()), "area_fraction": float(area[late].sum() / at), "area_m2": float(area[late].sum()),
           "facets_with_record_inside": int(s_in.sum()),
           "facets_first_seen_after_window_end_only": int((s_after & ~s_in & ~s_before).sum()),
           "variant_objective_points_without_late_observation_rows": int(late_a.sum()),
           "variant_every_row_of_the_base": int(late_c.sum()),
           "variant_positive_predicate": int(late_pos.sum()),
           "old_fine_tuning_reference": "265 facets on the old 12M base (no 1-s grid rows, gauge offset 0)",
           "A_orig_zero": int(np.sum(ctx.A_orig[late] == 0.0))}
    log("FACETS late_only: window %.1f min [%s .. %s] (t_end = last selected data point, row %d); %d objective points "
        "inside (%d observation, %d grid), %d before; %d parameter facets (area %.3e m2 = %.2e of the surface; old "
        "fine tuning: 265 on the old base); variants: without the late observation rows %d, every base row %d, "
        "'> 0' predicate %d; first seen only after the window end %d" % (
            window_min, t_start_utc, utc(i_end), i_end,
            rec["points_inside"], rec["points_inside_observation"], rec["points_inside_grid"], rec["points_before"],
            rec["facets"], rec["area_m2"], rec["area_fraction"], rec["variant_objective_points_without_late_observation_rows"],
            rec["variant_every_row_of_the_base"], rec["variant_positive_predicate"],
            rec["facets_first_seen_after_window_end_only"]))
    if not late.any():
        raise SystemExit("FACETS late_only: no facet is seen only inside the window")
    support = {"t_start": t_start, "t_end": t_end, "before_rows": rows.et < t_start,
               "inside_rows": (rows.et >= t_start) & (rows.et <= t_end)}
    return late, rec, support


# ------------------------------------------------------------------------------------------------ validation
def validate(obj, ctx, z_fit, args, timings_only=False, support=None):
    """The checks of the request (2026-11-07). Returns (record, passed)."""
    rng = np.random.default_rng(args.seed)
    rec = {"created_utc": utc_now(), "seed": args.seed}
    passed = True
    rows = obj.rows
    # (1) the fit's cost at delta = 0 and the fit's grain
    st = obj.evaluate(z_fit)
    fc = ctx.fit_costs
    one = {"fine_tuning": dict(st.costs),
           "fit": {"data_plus_regularization_cost": fc["cost"], "data_cost": fc["data_cost"],
                   "regularization_cost": fc["regularization_cost"], "crossing_penalty_cost": fc["crossing_penalty_cost"],
                   "total": fc["cost"] + fc["crossing_penalty_cost"], "fit_result_txt": ctx.fit_result_txt}}
    one["relative_difference"] = {
        "total_vs_fit_cost_plus_crossing": rel(st.costs["total"], one["fit"]["total"]),
        "data_plus_regularization_vs_fit_final_cost": rel(st.costs["data"] + st.costs["regularization"], fc["cost"]),
        "data": rel(st.costs["data"], fc["data_cost"]),
        "regularization": rel(st.costs["regularization"], fc["regularization_cost"]),
        "crossing": rel(st.costs["crossing"], fc["crossing_penalty_cost"])}
    tab = ctx.fit_tables()
    obs = ~rows.is_grid
    one["pressures_vs_fit_COPS_model_dat_max_rel"] = {
        "NG": float(np.max(np.abs(st.pNG[obs] / tab["pNG"][obs] - 1.0))),
        "RG": float(np.max(np.abs(st.pRG[obs] / tab["pRG"][obs] - 1.0))),
        "note": "the table holds 7 significant digits"}
    lw_rel = {"grid_NG": 0.0, "grid_RG": 0.0, "late_obs_NG": 0.0, "late_obs_RG": 0.0, "grain": 0.0}
    pos = {int(p): i for i, p in enumerate(rows.ipoint)}
    for ip, (kind, sng, srg, gr) in tab["late_window"].items():
        i = pos[ip]
        key = "grid" if kind == "grid" else "late_obs"
        lw_rel[key + "_NG"] = max(lw_rel[key + "_NG"], abs(st.pNG_w[i] / sng - 1.0))
        lw_rel[key + "_RG"] = max(lw_rel[key + "_RG"], abs(st.pRG_w[i] / srg - 1.0))
        lw_rel["grain"] = max(lw_rel["grain"], abs(st.g[i] - gr) / max(abs(gr), 1e-30))
    one["surface_vs_fit_late_window_model_dat_max_rel"] = dict(lw_rel, note="the file holds 10 significant digits")
    rd = one["relative_difference"]
    ok1 = (rd["total_vs_fit_cost_plus_crossing"] <= 1e-8 and rd["data"] <= 1e-8 and
           (rd["regularization"] <= 1e-8 or (fc["regularization_cost"] == 0.0 and st.costs["regularization"] == 0.0))
           and (rd["crossing"] <= 1e-8 or abs(st.costs["crossing"] - fc["crossing_penalty_cost"]) <= 1e-14))
    one["passed"] = bool(ok1)
    passed &= ok1
    rec["1_cost_at_fit"] = one
    log("check 1: total %.15e (fit %.15e, rel %.2e); data %.15e (fit %.15e, rel %.2e); reg %.15e (fit %.15e, "
        "rel %.2e); crossing %.3e (fit %.3e) -> %s" % (
            st.costs["total"], one["fit"]["total"], rd["total_vs_fit_cost_plus_crossing"], st.costs["data"],
            fc["data_cost"], rd["data"], st.costs["regularization"], fc["regularization_cost"], rd["regularization"],
            st.costs["crossing"], fc["crossing_penalty_cost"], "PASS" if ok1 else "FAIL"))
    log("check 1: pressures vs the fit's COPS_model.dat max rel NG %.2e RG %.2e; vs late_window_model.dat %s" % (
        one["pressures_vs_fit_COPS_model_dat_max_rel"]["NG"], one["pressures_vs_fit_COPS_model_dat_max_rel"]["RG"],
        {k: "%.1e" % v for k, v in lw_rel.items()}))
    # timings (6) and determinism
    tm = {}
    t0 = time.perf_counter()
    for _ in range(2):
        obj.evaluate(z_fit)
    tm["forward_s"] = (time.perf_counter() - t0) / 2
    st = obj.find(z_fit)
    J = Jacobian(obj, st)
    v1 = np.zeros(obj.n)
    v1[:obj.npar] = rng.standard_normal(obj.npar)
    t0 = time.perf_counter()
    Jv_a = J.matvec(v1)
    Jv_b = J.matvec(v1)
    tm["jv_1_column_s"] = (time.perf_counter() - t0) / 2
    V2 = np.zeros((obj.n, 2))
    V2[:obj.npar] = rng.standard_normal((obj.npar, 2))
    t0 = time.perf_counter()
    J.matmat(V2)
    tm["jv_2_columns_s"] = time.perf_counter() - t0
    u1 = rng.standard_normal(obj.m)
    t0 = time.perf_counter()
    g_a = J.rmatvec(u1)
    g_b = J.rmatvec(u1)
    tm["jtu_s"] = (time.perf_counter() - t0) / 2
    tm["engine_worker_seconds_last"] = {k: (v[-1] if v else None) for k, v in obj.engine.worker_seconds.items()}
    rec["timing_per_pass"] = tm
    rec["determinism"] = {"jv_bitwise_repeat": bool(np.array_equal(Jv_a, Jv_b)),
                          "jtu_bitwise_repeat": bool(np.array_equal(g_a, g_b))}
    log("timing per pass: forward %.2f s, J v %.2f s (2 columns %.2f s), J^T u %.2f s; repeat bitwise J v %s, "
        "J^T u %s" % (tm["forward_s"], tm["jv_1_column_s"], tm["jv_2_columns_s"], tm["jtu_s"],
                      rec["determinism"]["jv_bitwise_repeat"], rec["determinism"]["jtu_bitwise_repeat"]))
    if timings_only:
        return rec, passed
    # (2) adjoint tests
    adj = []
    eng = obj.engine
    for trial in range(3):
        # (a) the GPU H2O block alone: <JNG v, uN> + <JRG v, uR> == <v, g(uN, uR)> (parameter space)
        vp = rng.standard_normal(obj.npar)
        eng.direction_buffer(1)[0][:] = vp
        JNG, JRG = eng.jv(st, 1)
        uN, uR = rng.standard_normal(rows.n), rng.standard_normal(rows.n)
        gpar = eng.jtu(st, uN, uR)
        a_ = float(JNG[0] @ uN + JRG[0] @ uR)
        b_ = float(vp @ gpar)
        norm_ = float(np.sqrt((JNG[0] @ JNG[0]) + (JRG[0] @ JRG[0])) * np.sqrt(uN @ uN + uR @ uR))
        # (b) the full operator, data + crossing rows only (u = 0 on the Tikhonov rows) and the full residual space
        v = rng.standard_normal(obj.n)
        u = rng.standard_normal(obj.m)
        if obj.nreg:
            u_dc = u.copy()
            u_dc[rows.nd:rows.nd + obj.nreg] = 0.0
        else:
            u_dc = u
        Jv = J.matvec(v)
        JTu_dc = J.rmatvec(u_dc)
        JTu = J.rmatvec(u) if obj.nreg else JTu_dc
        zero_block = bool(not np.any(JNG[0]) and not np.any(JRG[0]) and not np.any(gpar))
        # a parameter block with zero sensitivity (every parameter record at the H2O floor): both sides must be exactly 0
        rel_n = (abs(a_ - b_) / norm_) if norm_ > 0 else (0.0 if a_ == b_ == 0.0 else float("inf"))
        res = {"gpu_block": {"Jv_u": a_, "v_JTu": b_, "rel": rel(a_, b_), "rel_to_norms": rel_n,
                             "zero_block": zero_block, "both_exactly_zero": bool(a_ == 0.0 and b_ == 0.0)},
               "operator_data_crossing": {"Jv_u": float(Jv @ u_dc), "v_JTu": float(v @ JTu_dc),
                                          "rel": rel(Jv @ u_dc, v @ JTu_dc)},
               "operator_full": {"Jv_u": float(Jv @ u), "v_JTu": float(v @ JTu), "rel": rel(Jv @ u, v @ JTu)}}
        adj.append(res)
        log("check 2 (adjoint) trial %d: GPU block %s; operator data+crossing rel %.2e; full rel %.2e" % (
            trial, ("ZERO BLOCK (J v and J^T u of the parameter facets exactly 0; <Jv,u> = %r, <v,J^T u> = %r)"
                    % (a_, b_)) if zero_block else ("rel %.2e (to norms %.2e)" % (res["gpu_block"]["rel"], rel_n)),
            res["operator_data_crossing"]["rel"], res["operator_full"]["rel"]))
    worst = max(max(r["gpu_block"]["rel"], r["operator_data_crossing"]["rel"], r["operator_full"]["rel"]) for r in adj)
    ok2 = worst <= 1e-10
    zb = [r["gpu_block"]["zero_block"] for r in adj]
    if any(zb):
        ok2 &= all(r["gpu_block"]["both_exactly_zero"] for r, z in zip(adj, zb) if z)
    rec["2_adjoint"] = {"trials": adj, "worst_rel": worst, "threshold": 1e-10, "passed": bool(ok2),
                        "zero_block_trials": int(sum(zb)),
                        "note": "zero block: the parameter facets have zero data sensitivity (all their records at "
                                "the H2O floor); then <Jv,u> and <v,J^T u> must both be exactly 0"}
    passed &= ok2
    # (3) J v against central differences of the full residual
    r0 = st.r
    g0 = J.rmatvec(r0)                                      # gradient of the cost at the fit
    blocks = {"data": np.arange(rows.nd),
              "data_late_window": np.concatenate([np.flatnonzero(rows.et[rows.useNG] >= rows.et_end),
                                                  rows.NNG + np.flatnonzero(rows.et[rows.useRG] >= rows.et_end)]),
              "regularization": np.arange(rows.nd, rows.nd + obj.nreg),
              "crossing": np.arange(rows.nd + obj.nreg, obj.m)}
    # Kinks: F = max(q A, F_min) is not differentiable where q A = F_min. With q in [n q_max, q_max] (n = the night
    # ratio) and F_min = r q_max A_orig, a facet's records can sit at the kink only for r/1 <= A/A_orig <= r/n
    # (r = floor ratio): kink-free are A_orig = 0 (F_min = 0), A > (r/n) A_orig (never at the floor) and A < r A_orig
    # (always at the floor). A central difference across a kink is off by up to half the slope jump, so a direction
    # over all facets mixes in records within h of their kink (error ~ h); the pass criterion is applied to the
    # kink-free random direction, the all-facet one is reported with an h-scan.
    night = float(os.environ.get("COPS_POS_NIGHT_RATIO", "nan"))
    lnr = np.full(obj.npar, np.inf)
    Ao = obj.A_orig[obj.idx]
    pos = Ao > 0
    lnr[pos] = st.s[obj.idx][pos] - np.log(Ao[pos])
    margin = 1e-3
    kink_free = (~pos) | (lnr > math.log(obj.floor_ratio / night) + margin) | (lnr < math.log(obj.floor_ratio) - margin)
    kf = {"night_ratio": night, "floor_ratio": obj.floor_ratio, "margin_nepers": margin,
          "kink_free_facets": int(kink_free.sum()), "parameter_facets": int(obj.npar),
          "facets_A_orig_0": int((~pos).sum()),
          "facets_never_at_floor": int((pos & (lnr > math.log(obj.floor_ratio / night) + margin)).sum()),
          "facets_always_at_floor": int((pos & (lnr < math.log(obj.floor_ratio) - margin)).sum())}
    dirs = []
    if kink_free.any():
        v = np.zeros(obj.n)
        v[:obj.npar] = np.where(kink_free, rng.standard_normal(obj.npar), 0.0)
        v[:obj.npar] /= np.max(np.abs(v[:obj.npar]))
        dirs.append(("random_kink_free_facets_maxabs1_h1e-6", v, 1e-6))
        dirs.append(("random_kink_free_facets_maxabs1_h1e-4", v, 1e-4))
    v = np.zeros(obj.n)
    v[:obj.npar] = rng.standard_normal(obj.npar)
    v[:obj.npar] /= np.max(np.abs(v[:obj.npar]))
    dirs.append(("random_all_facets_maxabs1_h1e-6", v, 1e-6))
    dirs.append(("random_all_facets_maxabs1_h1e-7", v, 1e-7))
    dirs.append(("random_all_facets_maxabs1_h1e-8", v, 1e-8))
    v = np.zeros(obj.n)
    top = np.argsort(-np.abs(g0[:obj.npar]))[:1000]
    v[top] = -np.sign(g0[top])
    dirs.append(("top1000_gradient_facets_descent", v, 1e-6))
    v = np.zeros(obj.n)
    v[-2] = 1.0
    dirs.append(("grain_A", v, 1e-6 * max(abs(z_fit[-2]), 1.0)))
    v = np.zeros(obj.n)
    v[-1] = 1.0
    dirs.append(("grain_B", v, 1e-6 * max(abs(z_fit[-1]), 1.0)))
    fd_rec = {}
    worst3 = 0.0
    eps = float(np.finfo(float).eps)
    rmax = {b: float(np.max(np.abs(r0[ix]))) if len(ix) else 0.0 for b, ix in blocks.items()}
    failed3 = []
    for name, v, h in dirs:
        Jv = J.matvec(v)
        rp = obj.evaluate(z_fit + h * v).r
        rm = obj.evaluate(z_fit - h * v).r
        fd = (rp - rm) / (2.0 * h)
        col = {"h": h}
        for b, ix in blocks.items():
            if len(ix):
                scale = float(np.max(np.abs(fd[ix])))
                err = float(np.max(np.abs(fd[ix] - Jv[ix])))
                col[b] = err / scale if scale > 0 else (0.0 if err == 0 else float("inf"))
                col[b + "_abs"] = err
                col[b + "_max_abs_fd"] = scale
                # FD rounding floor: residual values carry a few ulps, divided by 2h
                col[b + "_rounding_units"] = err / (eps * rmax[b] / h) if rmax[b] > 0 else (0.0 if err == 0 else
                                                                                             float("inf"))
        col["zero_direction"] = bool(not np.any(fd) and not np.any(Jv))
        fd_rec[name] = col
        if not name.startswith("random_all_facets"):
            for b in blocks:
                if b in col and not (col[b] <= 1e-6 or col[b + "_rounding_units"] <= 1e3):
                    failed3.append("%s/%s" % (name, b))
            worst3 = max(worst3, col.get("data", 0.0), col.get("crossing", 0.0), col.get("regularization", 0.0))
        log("check 3 (J v vs central FD) %-36s h %.0e%s  data %.2e (abs %.1e, %.0f rounding units)  late %.2e  reg %s  "
            "crossing %s" % (
            name, h, " ZERO (J v and FD exactly 0)" if col["zero_direction"] else "",
            col.get("data", float("nan")), col.get("data_abs", float("nan")),
            col.get("data_rounding_units", float("nan")), col.get("data_late_window", float("nan")),
            "%.2e" % col["regularization"] if "regularization" in col else "-",
            "%.2e" % col["crossing"] if "crossing" in col else "-"))
    ok3 = not failed3
    scan = [fd_rec["random_all_facets_maxabs1_h1e-%d" % k]["data_abs"] for k in (6, 7, 8)]
    rec["3_jv_vs_fd"] = {"directions": fd_rec, "worst_rel_kink_free_directions": worst3, "failed": failed3,
                         "passed": bool(ok3), "kinks": kf, "max_abs_residual_by_block": rmax,
                         "all_facets_abs_error_h_scan_1e-6_1e-7_1e-8": scan,
                         "pass_rule": "every row block (data, late window, Tikhonov, crossing) of every kink-free "
                                      "direction (kink-free random facets at h 1e-6 and 1e-4, top-1000 gradient "
                                      "facets, grain A, grain B): relative error <= 1e-6 OR absolute error <= 1000 "
                                      "rounding units eps max|r_block| / h (the FD accuracy limit; a random-sign "
                                      "direction gives a small J v, so its relative error is dominated by rounding); "
                                      "the all-facet random direction is reported with its h-scan (kink crossings "
                                      "give an error proportional to h, rounding one proportional to 1/h)",
                         "metric": "max|J v - FD| / max|FD| per row block; rounding units = abs error / "
                                   "(eps max|r_block| / h)"}
    log("check 3: kink-free facets %d of %d (A_orig 0: %d, never at the floor: %d, always: %d); all-facet random "
        "direction abs error at h 1e-6/1e-7/1e-8: %.2e / %.2e / %.2e" % (
            kf["kink_free_facets"], kf["parameter_facets"], kf["facets_A_orig_0"], kf["facets_never_at_floor"],
            kf["facets_always_at_floor"], scan[0], scan[1], scan[2]))
    passed &= ok3
    # (4) J^T u for single facets against the per-facet finite difference
    u = rng.standard_normal(obj.m)
    if obj.nreg:
        u[rows.nd:rows.nd + obj.nreg] = 0.0
    gu = J.rmatvec(u)
    nz = np.flatnonzero(g0[:obj.npar] != 0.0)
    order = nz[np.argsort(-np.abs(g0[nz]))]
    pick = list(order[:3])
    if len(order) > 10:
        mid = order[len(order) // 100: len(order) // 10]
        pick += list(rng.choice(mid, size=min(2, len(mid)), replace=False))
    zero = np.flatnonzero(g0[:obj.npar] == 0.0)
    zero_pick = [int(k) for k in rng.choice(zero, size=min(2, len(zero)), replace=False)] if len(zero) else []
    facets = []
    worst4 = 0.0
    zero_ok = True
    for k in zero_pick:
        # zero data sensitivity: the residual must not change at all (records at the floor do not depend on A)
        h = 1e-6
        e = np.zeros(obj.n)
        e[k] = 1.0
        stp, stm = obj.evaluate(z_fit + h * e), obj.evaluate(z_fit - h * e)
        ok_k = bool(np.array_equal(stp.rd, stm.rd) and np.array_equal(stp.q, stm.q) and gu[k] == 0.0)
        zero_ok &= ok_k
        facets.append({"parameter_index": int(k), "facet": int(obj.idx[k]), "h": h, "gradient": float(g0[k]),
                       "zero_sensitivity": True, "data_and_crossing_rows_bitwise_unchanged": ok_k,
                       "JTu": float(gu[k])})
        log("check 4 (single facet %d, zero gradient) data and crossing residuals bitwise unchanged at +-%g and "
            "J^T u = 0: %s (J^T u %r)" % (int(obj.idx[k]), h, ok_k, float(gu[k])))
    for k in pick:
        h = 1e-3
        e = np.zeros(obj.n)
        e[k] = 1.0
        stp, stm = obj.evaluate(z_fit + h * e), obj.evaluate(z_fit - h * e)
        col = (stp.r - stm.r) / (2.0 * h)
        dcost = (stp.cost - stm.cost) / (2.0 * h)
        item = {"parameter_index": int(k), "facet": int(obj.idx[k]), "h": h,
                "gradient": float(g0[k]), "fd_column_dot_r": float(col @ r0), "fd_cost": float(dcost),
                "rel_gradient_vs_fd_column": rel(g0[k], col @ r0), "rel_gradient_vs_fd_cost": rel(g0[k], dcost),
                "JTu": float(gu[k]), "fd_column_dot_u": float(col @ u), "rel_JTu_vs_fd": rel(gu[k], col @ u)}
        facets.append(item)
        worst4 = max(worst4, item["rel_gradient_vs_fd_column"], item["rel_JTu_vs_fd"])
        log("check 4 (single facet %d) gradient %.6e vs FD %.6e (rel %.2e; cost FD rel %.2e); J^T u %.6e vs FD "
            "%.6e (rel %.2e)" % (item["facet"], item["gradient"], item["fd_column_dot_r"],
                                 item["rel_gradient_vs_fd_column"], item["rel_gradient_vs_fd_cost"], item["JTu"],
                                 item["fd_column_dot_u"], item["rel_JTu_vs_fd"]))
    ok4 = worst4 <= 1e-3 and zero_ok and (len(pick) > 0 or len(zero_pick) > 0)
    rec["4_jtu_single_facets_vs_fd"] = {"facets": facets, "worst_rel": worst4, "threshold": 1e-3, "passed": bool(ok4),
                                        "parameters_with_zero_gradient": int(len(zero)),
                                        "parameters_with_nonzero_gradient": int(len(nz)),
                                        "zero_gradient_facets_checked": len(zero_pick)}
    passed &= ok4
    # (7) the speed-up (2026-11-08) against the first implementation, same random vectors
    if hasattr(eng, "jv_legacy"):
        seven = {"trials": []}
        for trial in range(2):
            vp = rng.standard_normal(obj.npar)
            vf = np.zeros(obj.nf)
            vf[obj.idx] = vp
            uN, uR = rng.standard_normal(rows.n), rng.standard_normal(rows.n)
            eng.direction_buffer(1)[0][:] = vp
            t0 = time.perf_counter()
            JN, JR = eng.jv(st, 1)
            t_new_jv = time.perf_counter() - t0
            eng.direction_buffer_legacy(1)[0][:] = vf
            t0 = time.perf_counter()
            JNl, JRl = eng.jv_legacy(st, 1)
            t_old_jv = time.perf_counter() - t0
            t0 = time.perf_counter()
            g = np.array(eng.jtu(st, uN, uR), copy=True)
            t_new_jtu = time.perf_counter() - t0
            t0 = time.perf_counter()
            gl = eng.jtu_legacy(st, uN, uR)[obj.idx]
            t_old_jtu = time.perf_counter() - t0
            mr = lambda x, y: float(np.max(np.abs(x - y)) / max(np.max(np.abs(y)), 1e-300))  # noqa: E731
            item = {"jv_NG_max_rel": mr(JN, JNl), "jv_RG_max_rel": mr(JR, JRl), "jtu_max_rel": mr(g, gl),
                    "jv_bitwise": bool(np.array_equal(JN, JNl) and np.array_equal(JR, JRl)),
                    "jtu_bitwise": bool(np.array_equal(g, gl)),
                    "seconds_new": {"jv": t_new_jv, "jtu": t_new_jtu},
                    "seconds_legacy": {"jv": t_old_jv, "jtu": t_old_jtu}}
            seven["trials"].append(item)
            log("check 7 (speed-up vs first implementation) trial %d: J v max rel NG %.1e RG %.1e (bitwise %s), "
                "J^T u max rel %.1e (bitwise %s); engine seconds J v %.3f vs %.3f, J^T u %.3f vs %.3f" % (
                    trial, item["jv_NG_max_rel"], item["jv_RG_max_rel"], item["jv_bitwise"], item["jtu_max_rel"],
                    item["jtu_bitwise"], t_new_jv, t_old_jv, t_new_jtu, t_old_jtu))
        worst7 = max(max(i["jv_NG_max_rel"], i["jv_RG_max_rel"], i["jtu_max_rel"]) for i in seven["trials"])
        seven.update(worst_rel=worst7, threshold=1e-12, passed=bool(worst7 <= 1e-12),
                     note="the same additions in the same order are expected, i.e. bitwise equality")
        rec["7_speedup_vs_first_implementation"] = seven
        passed &= seven["passed"]
    # (8) strict support (FACETS late_only): a random change of every parameter leaves the surface pressures of every
    # row before the window start bitwise unchanged and changes rows inside the window
    if support is not None:
        z = np.array(z_fit, dtype=float, copy=True)
        z[:obj.npar] = np.maximum(obj.lb_delta, rng.uniform(-1.0, 1.0, obj.npar))
        s0_, s1_ = obj.evaluate(z_fit), obj.evaluate(z)
        b, i = support["before_rows"], support["inside_rows"]
        same_ng = bool(np.array_equal(s0_.pNG_s[b], s1_.pNG_s[b]))
        same_rg = bool(np.array_equal(s0_.pRG_s[b], s1_.pRG_s[b]))
        same_tab = bool(np.array_equal(s0_.pNG[b], s1_.pNG[b]) and np.array_equal(s0_.pRG[b], s1_.pRG[b]))
        changed = int(np.sum(i & ((s0_.pNG_s != s1_.pNG_s) | (s0_.pRG_s != s1_.pRG_s))))
        big = float(np.max(np.abs(np.log(s1_.pNG_s[i] / s0_.pNG_s[i])))) if i.any() else 0.0
        # a larger increase (+5 nepers above max(bound, 0)): strict support again, and an effect inside the window
        z2 = np.array(z_fit, dtype=float, copy=True)
        z2[:obj.npar] = np.maximum(obj.lb_delta, 0.0) + 5.0
        s2_ = obj.evaluate(z2)
        same2 = bool(np.array_equal(s0_.pNG_s[b], s2_.pNG_s[b]) and np.array_equal(s0_.pRG_s[b], s2_.pRG_s[b])
                     and np.array_equal(s0_.pNG[b], s2_.pNG[b]) and np.array_equal(s0_.pRG[b], s2_.pRG[b]))
        changed2 = int(np.sum(i & ((s0_.pNG_s != s2_.pNG_s) | (s0_.pRG_s != s2_.pRG_s))))
        eight = {"rows_before_window": int(b.sum()), "surface_NG_bitwise_unchanged": same_ng,
                 "surface_RG_bitwise_unchanged": same_rg, "table_pressures_bitwise_unchanged_same_grain": same_tab,
                 "rows_inside_window": int(i.sum()), "rows_inside_changed": changed,
                 "max_abs_ln_change_NG_inside": big, "perturbation": "delta ~ U(-1, 1) nepers (clipped at the bound)",
                 "plus5_rows_before_bitwise_unchanged": same2, "plus5_rows_inside_changed": changed2,
                 "parameters_without_effect": bool(changed == 0 and changed2 == 0)}
        # pass = strict support in both perturbations; no change inside the window is flagged, not failed (a zero-
        # sensitivity parameter set, e.g. every record at the floor, is a valid outcome to report)
        eight["passed"] = bool(same_ng and same_rg and same_tab and same2)
        rec["8_strict_support"] = eight
        passed &= eight["passed"]
        log("check 8 (strict support): %d rows before the window: surface NG/RG bitwise unchanged %s/%s (table %s; "
            "+5 nepers %s); rows inside the window changed: %d of %d (max |ln change| NG %.2e), +5 nepers %d%s -> %s" % (
                eight["rows_before_window"], same_ng, same_rg, same_tab, same2, changed, eight["rows_inside_window"],
                big, changed2, "; WARNING the parameters have no effect" if eight["parameters_without_effect"] else "",
                "PASS" if eight["passed"] else "FAIL"))
    # (5, start) bounds: the start is feasible
    lbmax = float(np.max(obj.lb_delta)) if obj.npar else 0.0
    rec["5_bounds_start"] = {"max_lower_bound_delta": lbmax, "start_feasible": bool(lbmax <= 0.0),
                             "facets_start_below_floor_bound": int(obj.start_below_bound.sum()),
                             "facets_without_bound_A_orig_0": int(np.sum(~np.isfinite(obj.lb_delta)))}
    passed &= lbmax <= 0.0
    rec["passed"] = bool(passed)
    return rec, passed


# ------------------------------------------------------------------------------------------------ optimizer
def optimize(obj, z0, args):
    from scipy.optimize import least_squares
    x_scale = np.ones(obj.n)
    x_scale[-2:] = args.x_scale_grain
    t0 = time.time()
    result, termination, message = None, None, None
    try:
        result = least_squares(obj.fun, z0, jac=obj.jac, bounds=(obj.lb, obj.ub), method="trf", tr_solver="lsmr",
                               tr_options={"maxiter": int(args.lsmr_maxiter)}, x_scale=x_scale,
                               max_nfev=int(args.max_nfev), ftol=args.ftol, xtol=args.xtol, gtol=args.gtol,
                               verbose=2)
        termination = {0: "max_nfev", 1: "gtol", 2: "ftol", 3: "xtol", 4: "ftol_and_xtol"}.get(result.status,
                                                                                              str(result.status))
        message = str(result.message)
        z = np.asarray(result.x, float)
    except StopFit as exc:
        termination, message = "wall_time", str(exc)
        if obj.tracker.best is None:
            raise SystemExit("stopped before any evaluation finished")
        z = obj.tracker.best.z
    st = obj.find(z)
    if st is None:
        st = obj.evaluate(z)
    if obj.tracker.best is not None and obj.tracker.best.cost < st.cost:
        log("WARNING the final state (cost %.10e) is not the best evaluated (%.10e)" % (st.cost, obj.tracker.best.cost))
    rec = {"termination": termination, "message": message, "seconds": time.time() - t0,
           "nfev": obj.counters["fun"], "njev": int(result.njev) if result is not None else None,
           "optimality": float(result.optimality) if result is not None else None,
           "status": int(result.status) if result is not None else None,
           "settings": {"method": "trf", "tr_solver": "lsmr", "lsmr_maxiter": int(args.lsmr_maxiter),
                        "max_nfev": int(args.max_nfev), "ftol": args.ftol, "xtol": args.xtol, "gtol": args.gtol,
                        "x_scale": "1 (delta), %g (grain)" % args.x_scale_grain, "verbose": 2},
           "counters": dict(obj.counters), "seconds_by_kind": dict(obj.seconds),
           "engine_reuploads": getattr(obj.engine, "reuploads", None)}
    if result is not None:
        am = np.asarray(result.active_mask)
        rec["active_mask"] = {"lower": int(np.sum(am[:obj.npar] < 0)), "upper": int(np.sum(am[:obj.npar] > 0)),
                              "grain": am[-2:].tolist()}
    return st, rec


# ------------------------------------------------------------------------------------------------ products
def replace_columns(lines, values, ipoint_to_row, note_lines):
    """The fit's COPS_model.dat lines with pNG/pRG of every row replaced (%.6e) and note lines inserted after line 0.
    Returns (new lines, {iPoint: (pNG token, pRG token, refNG token, refRG token, UTC token, useNG, useRG)})."""
    k_hdr = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    names = lines[k_hdr][1:].split()
    c = {n: i for i, n in enumerate(names)}
    pNG, pRG = values
    out, tokens, seen_rows = [], {}, set()
    for i, line in enumerate(lines):
        if i <= k_hdr or line.startswith("#") or not line.strip():
            out.append(line)
            continue
        v = line.split()
        ip = int(v[c["iPoint"]])
        r = ipoint_to_row[ip]
        if r in seen_rows:
            raise SystemExit("COPS_model.dat lists iPoint %d twice" % ip)
        seen_rows.add(r)
        v[c["pNG"]] = "%.6e" % pNG[r]
        v[c["pRG"]] = "%.6e" % pRG[r]
        out.append(" ".join(v) + "\n")
        tokens[ip] = (v[c["pNG"]], v[c["pRG"]], v[c["refNG"]], v[c["refRG"]], v[c["UTC"]], v[c["useNG"]],
                      v[c["useRG"]])
    if len(seen_rows) != len(ipoint_to_row):
        raise SystemExit("COPS_model.dat has %d rows, the model %d" % (len(seen_rows), len(ipoint_to_row)))
    return out[:1] + list(note_lines) + out[1:], tokens


def write_tables(out, ctx, obj, st_fit, st, utc):
    rows = obj.rows
    out = Path(out)
    lines = ctx.files["COPS_model.dat"].read_text().splitlines(True)
    ip2row = {int(p): i for i, p in enumerate(rows.ipoint)}
    note = ["# PositiveFineTuning (model_ext/positive_finetune.py, %s): pNG/pRG of EVERY row recomputed with the "
            "per-element fine tuning (A = exp(s_fit + held + delta), finetune_state.npz) of the fit %s; grain A %.9e "
            "Pa, B %.9e Pa/s refitted; late observations: W10x5(surface) + grain(t_obs); grid rows: instantaneous "
            "surface + grain (see PositiveFitGridRows below; late_window_model.dat lists them); every other column "
            "and header line is the fit's\n" % (utc, ctx.fit_dir, st.z[-2] * rows.A_UNIT, st.z[-1] * rows.B_UNIT)]
    new, tok = replace_columns(lines, (st.pNG, st.pRG), ip2row, note)
    with (out / "COPS_model.dat").open("x") as f:
        f.write("".join(new))
    # late_window_model.dat (positive_problem.write_late_window_products' format)
    lay = ctx.layout
    stop = {int(r): float(t) for r, t in zip(lay.late_rows, lay.stop_et)}
    lr = sorted([int(x) for x in lay.grid_rows] + [int(x) for x in lay.late_rows], key=lambda i: rows.et[i])
    with (out / "late_window_model.dat").open("x") as f:
        f.write("# late window of the fine tuning (model_ext/positive_finetune.py, %s) of %s: every 1-s grid row (model "
                "only) and every observation from %s UTC on\n" % (utc, ctx.fit_dir, ctx.layout.summary["late_start_utc"]))
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
                ip, t[4], rows.et[i], kind, "%.6f" % stop[i] if i in stop else "-", st.pNG_w[i], st.pRG_w[i],
                st.g[i], t[0], t[1], t[2], t[3], int(rows.useNG[i]), int(rows.useRG[i])))
    # pressure_comparison.txt
    with (out / "pressure_comparison.txt").open("x") as f:
        f.write("# per row of the base: the tuned fit (delta = 0, the fit's grain; recomputed here) and the fine "
                "tuning; kind obs | late_obs (W10x5 + grain) | grid (model only); ln_* = ln(model / reference) where "
                "the reference is > 0 (nan otherwise); pNG includes the grain\n")
        f.write("#iPoint UTC et kind useNG useRG refNG_Pa refRG_Pa pNG_fit_Pa pNG_tuned_Pa pRG_fit_Pa pRG_tuned_Pa "
                "ln_NG_fit ln_NG_tuned ln_RG_fit ln_RG_tuned\n")
        for i in range(rows.n):
            kind = "grid" if rows.is_grid[i] else ("late_obs" if rows.is_late[i] else "obs")
            lnf = lambda p, r: math.log(p / r) if r > 0 and p > 0 else float("nan")  # noqa: E731
            f.write("%d %s %.6f %s %d %d %.6e %.6e %.9e %.9e %.9e %.9e %.6f %.6f %.6f %.6f\n" % (
                rows.ipoint[i], ctx.utc[i], rows.et[i], kind, int(rows.useNG[i]), int(rows.useRG[i]),
                rows.refNG[i], rows.refRG[i], st_fit.pNG[i], st.pNG[i], st_fit.pRG[i], st.pRG[i],
                lnf(st_fit.pNG[i], rows.refNG[i]), lnf(st.pNG[i], rows.refNG[i]),
                lnf(st_fit.pRG[i], rows.refRG[i]), lnf(st.pRG[i], rows.refRG[i])))


def write_state(out, ctx, obj, st, meta):
    delta = np.zeros(obj.nf)
    delta[obj.idx] = st.z[:obj.npar]
    lb = np.full(obj.nf, np.nan)
    lb[obj.idx] = obj.lb_delta
    held_next = ctx.held_prev + delta
    s_final = st.s
    path = Path(out) / "finetune_state.npz"
    with path.open("xb") as f:
        np.savez_compressed(f, delta=delta, held_prev=ctx.held_prev, held_next=held_next, s_final=s_final,
                            A_final=np.exp(s_final), lower_bound=lb, seen=obj.seen,
                            grain_A_uPa=np.array(st.z[-2]), grain_B_uPa_per_100s=np.array(st.z[-1]),
                            meta_json=np.array(json.dumps(jsonable(meta), sort_keys=True)))
    check = np.max(np.abs(s_final - (ctx.s_param + held_next)))
    return path, float(check)


def stop_reason(termination):
    """'tolerance' (ftol, xtol, gtol), 'max_nfev' or 'wall_time'."""
    if termination == "wall_time":
        return "wall_time"
    if termination == "max_nfev":
        return "max_nfev"
    if termination in ("ftol", "xtol", "gtol", "ftol_and_xtol"):
        return "tolerance"
    return "other"


RESULT_KEYS = ("cost_total", "cost_data", "cost_regularization", "cost_crossing", "rms_ln_NG", "rms_ln_RG",
               "late_window_rms_ln_NG", "late_window_rms_ln_RG", "D_nepers", "grain_A_uPa", "grain_B_uPa_per_100s",
               "grain_A_Pa", "grain_B_Pa_per_s", "facets_at_lower_bound", "area_fraction_at_lower_bound_of_total",
               "area_fraction_at_lower_bound_of_parameters", "delta_range", "delta_rms_area_weighted_parameters",
               "crossing_wrong_side_pairs", "total_activity_over_A_orig")


def write_result_txt(out, ctx, obj, s_init, s_final, opt, seen_rec, utc, start_label, seg, sens=None):
    with (Path(out) / "finetune_result.txt").open("x") as f:
        f.write("# per-element fine tuning (model_ext/positive_finetune.py, %s) of %s\n" % (utc, ctx.fit_dir))
        f.write("# A = exp(s_fit + held + delta); %d parameter facets (%.4f of the area) + grain A, B; lambda %g "
                "(the fit's); floor ratio %g; costs: total = data + regularization + crossing (the fit's "
                "final_cost = data + regularization)\n" % (obj.npar, seen_rec["parameter_area_fraction"], obj.lam,
                                                          obj.floor_ratio))
        f.write("# initial = %s\n" % start_label)
        lo = seen_rec.get("late_only") or {}
        f.write("facets %s%s\n" % (seen_rec.get("facets_rule"), "  WINDOW_MIN %g  window %s .. %s (et %.3f .. %.3f)" % (
            lo["window_min"], lo["t_start_utc"], lo["t_end_utc"], lo["t_start_et"], lo["t_end_et"]) if lo else ""))
        f.write("lambda %.6e\nparameters %d\n" % (obj.lam, obj.n))
        if sens:
            f.write("parameters_with_zero_data_sensitivity  start %d  end %d  of %d facet parameters (data + crossing "
                    "Jacobian column exactly 0: every record of the facet at the H2O floor)\n" % (
                        sens["start"]["zero_sensitivity"], sens["end"]["zero_sensitivity"], obj.npar))
        for key in RESULT_KEYS:
            a, b = s_init.get(key), s_final.get(key)
            fmt = (lambda v: "%.10e" % v if isinstance(v, float) else str(v))  # noqa: E731
            f.write("%-46s initial %s  final %s\n" % (key, fmt(a), fmt(b)))
        f.write("fit_final_cost(fit_result.txt) %s  fit_manifest cost %.15e data %.15e reg %.15e crossing %.15e\n" % (
            ctx.fit_result_txt.get("final_cost"), ctx.fit_costs["cost"], ctx.fit_costs["data_cost"],
            ctx.fit_costs["regularization_cost"], ctx.fit_costs["crossing_penalty_cost"]))
        reason = stop_reason(opt["termination"])
        f.write("stop %s   (tolerance = converged on ftol/xtol/gtol; max_nfev = the evaluation budget (MAXNFEV) used up; "
                "wall_time = stopped at the best evaluated state before the job's time limit)\n" % reason)
        f.write("complete %s\n" % ("no (resubmit the same sbatch command: the next segment continues from this state)"
                                   if reason == "wall_time" else "yes"))
        f.write("termination %s\nmessage %s\nsegment %d\nnfev %s\nnfev_total %d  (all segments; MAXNFEV counts these)\n"
                "njev %s\noptimality %s\nwall_seconds %.1f\n" % (
                    opt["termination"], opt["message"], seg["segment"], opt["nfev"], seg["nfev_total"], opt["njev"],
                    opt["optimality"], opt["seconds"]))
        f.write("jv_passes %d jv_columns %d jtu %d\n" % (opt["counters"]["jv_passes"], opt["counters"]["jv_columns"],
                                                      opt["counters"]["jtu"]))
        if "active_mask" in opt:
            f.write("active_bounds lower %d upper %d grain %s\n" % (opt["active_mask"]["lower"],
                                                                   opt["active_mask"]["upper"],
                                                                   opt["active_mask"]["grain"]))


# ------------------------------------------------------------------------------------------------ main
def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True, help="the tuned fit's segment directory (positive_fit.py products)")
    ap.add_argument("--out", required=True, help="new output directory (may hold only provenance-* directories)")
    ap.add_argument("--floor-ratio", type=float, required=True,
                    help="r of the lower bound A >= r A_orig; must equal COPS_POS_FLOOR_RATIO of the fit")
    ap.add_argument("--devices", required=True, help="GPU numbers as CUDA numbers them in the job, e.g. 0,1,...,7")
    ap.add_argument("--max-nfev", type=int, default=None, help="optimizer residual evaluations (required to fit)")
    ap.add_argument("--lsmr-maxiter", type=int, default=50, help="LSMR iterations per outer iteration (default 50)")
    ap.add_argument("--ftol", type=float, default=1e-10)
    ap.add_argument("--xtol", type=float, default=1e-10)
    ap.add_argument("--gtol", type=float, default=1e-10)
    ap.add_argument("--x-scale-grain", type=float, default=1.0, help="x_scale of grain A, B (default 1)")
    ap.add_argument("--facets", choices=("late_only", "all_seen"), required=True,
                    help="the parameter facets (required, no default: user 2026-11-08): late_only = facets with a "
                         "structural record (Kn, KPhi or KV != 0) at an objective point inside the last --window-min "
                         "minutes and none at any objective point before the window start (the old late-window fine "
                         "tuning's set); all_seen = every facet seen from a point the objective uses")
    ap.add_argument("--window-min", type=float, default=None,
                    help="late_only: window length [min], ending at the last selected data point (required)")
    ap.add_argument("--seen-points", choices=("objective", "objective+late"), default="objective",
                    help="all_seen: facets seen from the points whose model enters the objective (default), or also "
                         "from the late observations' own records")
    ap.add_argument("--validate", action="store_true", help="run the checks first (refuse to fit if one fails)")
    ap.add_argument("--validate-only", action="store_true", help="the checks only (no optimization)")
    ap.add_argument("--seed", type=int, default=20261107)
    ap.add_argument("--resume-from", default=None, help="finetune_state.npz or finetune_checkpoint.npz of a previous "
                                                        "fine tuning of the same fit (its delta and grain = start)")
    ap.add_argument("--wall-time", default=None, help="slurm | seconds | [D-]HH:MM:SS (stop at the best state)")
    ap.add_argument("--wall-margin", type=float, default=300.0)
    ap.add_argument("--save-seconds", type=float, default=900.0,
                    help="time kept for the products before the wall-time limit (default 900 s)")
    ap.add_argument("--segment", type=int, default=1,
                    help="segment number of a fine tuning split by wall time (job script; > 1 needs --resume-from)")
    ap.add_argument("--nfev-before", type=int, default=0,
                    help="evaluations of the earlier segments (job script; recorded as nfev_total)")
    ap.add_argument("--checkpoint-minutes", type=float, default=0.0,
                    help="rolling finetune_checkpoint.npz at most every M minutes while improving (0: none)")
    args = ap.parse_args(argv)
    if not args.validate_only and args.max_nfev is None:
        ap.error("--max-nfev is required unless --validate-only")
    if args.lsmr_maxiter < 1:
        ap.error("--lsmr-maxiter must be >= 1")
    if args.facets == "late_only" and not (args.window_min and args.window_min > 0):
        ap.error("--facets late_only needs --window-min > 0")
    if args.facets == "all_seen" and args.window_min is not None:
        ap.error("--window-min applies to --facets late_only only")
    if args.segment > 1 and not args.resume_from:
        ap.error("--segment > 1 continues a stopped segment: --resume-from <its finetune_state.npz> is required")
    return args


def prepare_out(out):
    out = Path(out).resolve()
    if out.exists():
        extra = [p.name for p in out.iterdir() if not p.name.startswith("provenance-")]
        if extra:
            raise SystemExit("--out %s exists with %s (a new directory is required; only provenance-* may exist)"
                             % (out, extra[:5]))
    else:
        out.mkdir(parents=True)
    return out


def load_resume(path, ctx, obj):
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        delta = np.asarray(z["delta"], float)
        seen = np.asarray(z["seen"], bool) if "seen" in z.files else None
        xA, xB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    want = ctx.identity["fit_restart.npz"]["sha256"]
    got = meta.get("fit_restart_sha256")
    if got != want:
        raise SystemExit("--resume-from %s belongs to the fit restart %s, not %s" % (path, got, want))
    if delta.shape != (obj.nf,) or np.any(delta[~obj.seen] != 0.0):
        raise SystemExit("--resume-from: delta has another shape or nonzero values on non-parameter facets")
    if seen is not None and not np.array_equal(seen, obj.seen):
        raise SystemExit("--resume-from: its parameter facets differ from this run's")
    z0 = obj.start_vector(xA, xB, delta)
    bad = int(np.sum(z0 < obj.lb))
    if bad:
        raise SystemExit("--resume-from: %d parameters below their lower bound" % bad)
    return z0, {"path": str(Path(path).resolve()), "sha256": sha256(path), "meta": meta}


def main(argv=None):
    started = time.time()
    args = parse_args(sys.argv[1:] if argv is None else argv)
    command = [str(Path(__file__).resolve())] + (sys.argv[1:] if argv is None else list(argv))
    log("command: %s" % " ".join(command))
    ctx = FitContext(args.fit)
    env = check_environment(ctx, args.floor_ratio)
    log("environment %s (= the fit's); floor ratio %g" % (env["process"], args.floor_ratio))
    from cops_cuda import wall_time
    deadline = wall_time.resolve_deadline(args.wall_time, None, args.wall_margin, start=started)
    out = prepare_out(args.out)
    tracker = Tracker(deadline, args.save_seconds,
                      out / "finetune_checkpoint.npz" if args.checkpoint_minutes > 0 else None,
                      60.0 * args.checkpoint_minutes)
    devices = [int(d) for d in args.devices.split(",") if d.strip() != ""]
    engine = None
    record = {"schema": "positive_finetune", "parameterization": PARAMETERIZATION_ID, "policy": POLICY,
              "created_utc": utc_now(), "command": command, "host": socket.gethostname(),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "fit": {
                  "directory": str(ctx.fit_dir), "identities": ctx.identity, "parameterization": ctx.fit_parameterization,
                  "lambda": ctx.lam, "costs": ctx.fit_costs, "fit_result_txt": ctx.fit_result_txt,
                  "grain_term": ctx.grain_term, "x": ctx.x.tolist(), "names": ctx.names,
                  "held_present": ctx.held_present, "exponent_checks": ctx.exponent_checks,
                  "run": str(ctx.run), "point_table": file_identity(ctx.point_table),
                  "late_window": jsonable(ctx.layout.summary)},
              "environment": env, "settings": jsonable(vars(args)),
              "modules": {p: file_identity(HERE / p) for p in ("positive_finetune.py", "positive_problem.py",
                                                               "positive_activity.py", "grain_crossing.py")},
              "engine": {"devices": devices}, "method_choices": METHOD_CHOICES}
    rec_path = out / "fine_tuning.json"
    try:
        engine = GpuEngine(ctx, devices)
        record["engine"].update(bind_s=engine.bind_s, workers=engine.attach)
        # --- parameters: the facets seen from the objective's points
        rows = ctx.rows
        own, late_used = rows.objective_points()
        t0 = time.time()
        seen_obj = engine.seen(own)
        seen_late = engine.seen(own | late_used)
        seen_all = engine.seen(np.ones(rows.n, bool))
        area_total = float(np.sum(ctx.area))
        seen = seen_obj if args.seen_points == "objective" else seen_late
        seen_rec = {
            "rule": args.seen_points, "objective_points": int(own.sum()),
            "objective_observation_rows": int((own & ~rows.is_grid).sum()),
            "objective_grid_rows": int((own & rows.is_grid).sum()), "late_observation_rows": int(late_used.sum()),
            "facets_objective": int(seen_obj.sum()), "area_fraction_objective": float(ctx.area[seen_obj].sum() / area_total),
            "facets_objective_plus_late_own": int(seen_late.sum()),
            "facets_seen_only_from_late_observations": int((seen_late & ~seen_obj).sum()),
            "area_fraction_seen_only_from_late_observations": float(ctx.area[seen_late & ~seen_obj].sum() / area_total),
            "facets_all_rows": int(seen_all.sum()),
            "all_rows_equal_fit_npz_seen": bool(np.array_equal(seen_all, ctx.seen_fit)),
            "facets_fit_npz_seen": int(ctx.seen_fit.sum()),
            "parameters": int(seen.sum()), "parameter_area_fraction": float(ctx.area[seen].sum() / area_total),
            "seconds": time.time() - t0}
        log("parameters: %d facets (%.4f of the area) seen from %d objective points (%d observation + %d grid rows); "
            "+%d facets seen only from the %d late observations' own records (not parameters with rule %s); all rows "
            "%d facets (= fit's seen: %s)" % (
                seen_rec["parameters"], seen_rec["parameter_area_fraction"], seen_rec["objective_points"],
                seen_rec["objective_observation_rows"], seen_rec["objective_grid_rows"],
                seen_rec["facets_seen_only_from_late_observations"], seen_rec["late_observation_rows"], args.seen_points,
                seen_rec["facets_all_rows"], seen_rec["all_rows_equal_fit_npz_seen"]))
        support = None
        if args.facets == "late_only":
            seen, late_rec, support = late_only_facets(engine, rows, ctx, args.window_min, own, late_used)
            seen_rec.update(late_only=late_rec, parameters=int(seen.sum()),
                            parameter_area_fraction=float(ctx.area[seen].sum() / area_total))
        seen_rec["facets_rule"] = args.facets
        record["seen"] = seen_rec
        record["facets"] = {"FACETS": args.facets, "WINDOW_MIN": args.window_min,
                            "definition": FACET_RULES[args.facets]}
        fit_identity = {"fit_directory": str(ctx.fit_dir),
                        "fit_restart_sha256": ctx.identity["fit_restart.npz"]["sha256"]}
        obj = Objective(rows, engine, ctx.s_start, ctx.s0, ctx.area, ctx.A_orig, seen, ctx.lam, args.floor_ratio,
                        tracker, fit_identity)
        record["bounds"] = {"facets_start_below_floor_bound_lb0": int(obj.start_below_bound.sum()),
                            "facets_without_bound_A_orig_0": int(np.sum(~np.isfinite(obj.lb_delta))),
                            "lower_bound_delta_range_finite": [float(np.min(obj.lb_delta[np.isfinite(obj.lb_delta)])),
                                                               float(np.max(obj.lb_delta))] if obj.npar else None,
                            "rule": "delta >= min(0, ln(r A_orig / A_start)) for A_orig > 0; none for A_orig = 0; "
                                    "grain A, B >= 0; no upper bound"}
        log("objective: %d residual rows (%d data, %d Tikhonov, %d crossing), %d parameters; %d facets start below "
            "the floor bound (bound 0), %d without bound (A_orig = 0)" % (
                obj.m, rows.nd, obj.nreg, rows.npen, obj.n, record["bounds"]["facets_start_below_floor_bound_lb0"],
                record["bounds"]["facets_without_bound_A_orig_0"]))
        record["engine"]["parameter_maps"] = obj.param_info
        z_fit = obj.start_vector(ctx.xA, ctx.xB)
        st_fit = obj.evaluate(z_fit)
        s_fit_summary = obj.summary(st_fit)
        rel_cost = rel(st_fit.costs["data"] + st_fit.costs["regularization"], ctx.fit_costs["cost"])
        log("the fit's state: total %.15e (data %.15e, reg %.15e, crossing %.3e); fit cost %.15e (rel %.2e)" % (
            st_fit.costs["total"], st_fit.costs["data"], st_fit.costs["regularization"], st_fit.costs["crossing"],
            ctx.fit_costs["cost"], rel_cost))
        record["start_check"] = {"fit_state": s_fit_summary, "relative_difference_data_plus_reg_vs_fit_cost": rel_cost}
        sens = {"start": zero_sensitivity(obj, st_fit)}
        record["sensitivity"] = sens
        log("parameters with zero data sensitivity at the start: %d of %d%s" % (
            sens["start"]["zero_sensitivity"], obj.npar, " (ALL: only the grain can change the fit)"
            if sens["start"]["zero_sensitivity"] == obj.npar else ""))
        if rel_cost > 1e-8:
            raise SystemExit("the fine tuning does not reproduce the fit's cost (rel %.2e > 1e-8)" % rel_cost)
        valid = None
        if args.validate or args.validate_only:
            valid, ok = validate(obj, ctx, z_fit, args, support=support)
            with (out / "validation.json").open("x") as f:
                json.dump(jsonable(valid), f, indent=1)
            record["validation"] = {"passed": bool(ok), "file": str(out / "validation.json")}
            log("validation %s" % ("PASSED" if ok else "FAILED"))
            if not ok:
                raise SystemExit("validation failed: see %s" % (out / "validation.json"))
        if args.validate_only:
            record["stats"] = {"workers": engine.stats(), "driver_max_rss_gb": max_rss_gb()}
            return 0
        resume = None
        z0 = z_fit
        if args.resume_from:
            z0, resume = load_resume(args.resume_from, ctx, obj)
            record["resume"] = resume
        st0 = st_fit if z0 is z_fit else obj.evaluate(z0)
        s_init = obj.summary(st0)
        log("start: total %.10e (data %.10e, reg %.10e, crossing %.3e) rms ln NG %.4f RG %.4f D %.6f" % (
            s_init["cost_total"], s_init["cost_data"], s_init["cost_regularization"], s_init["cost_crossing"],
            s_init["rms_ln_NG"], s_init["rms_ln_RG"], s_init["D_nepers"]))
        obj.counters.update({k: 0 for k in obj.counters})
        obj.seconds.update({k: 0.0 for k in obj.seconds})
        for v in engine.worker_seconds.values():
            del v[:]
        st, opt = optimize(obj, z0, args)
        s_final = obj.summary(st)
        reason = stop_reason(opt["termination"])
        prev_segments = list((resume or {}).get("meta", {}).get("segments") or []) if args.segment > 1 else []
        seg = {"segment": args.segment, "nfev": opt["nfev"], "nfev_before": args.nfev_before,
               "nfev_total": args.nfev_before + opt["nfev"], "termination": opt["termination"], "stop": reason,
               "cost": s_final["cost_total"], "out": str(out), "job": os.environ.get("SLURM_JOB_ID"),
               "finished_utc": utc_now()}
        segments = prev_segments + [seg]
        record["segment"] = seg
        record["segments"] = segments
        log("final (%s): total %.10e (data %.10e, reg %.10e, crossing %.3e) rms ln NG %.4f RG %.4f D %.6f; grain A "
            "%.6e B %.6e; %d facets at the lower bound; %d evaluations in %.0f s" % (
                opt["termination"], s_final["cost_total"], s_final["cost_data"], s_final["cost_regularization"],
                s_final["cost_crossing"], s_final["rms_ln_NG"], s_final["rms_ln_RG"], s_final["D_nepers"],
                s_final["grain_A_uPa"], s_final["grain_B_uPa_per_100s"], s_final["facets_at_lower_bound"],
                opt["nfev"], opt["seconds"]))
        bounds_ok = bool(np.all(st.z >= obj.lb) and np.all(st.z <= obj.ub))
        record["bounds"]["final_respected"] = bounds_ok
        record["bounds"]["final_min_delta_minus_lb"] = float(np.min(st.z[:obj.npar] - obj.lb_delta)) if obj.npar else None
        sens["end"] = zero_sensitivity(obj, st)
        log("parameters with zero data sensitivity at the end: %d of %d" % (sens["end"]["zero_sensitivity"], obj.npar))
        stats = engine.stats()
        record["stats"] = {"workers": stats, "driver_max_rss_gb": max_rss_gb(),
                           "engine_worker_seconds_mean": {k: (float(np.mean(v)) if v else None)
                                                          for k, v in engine.worker_seconds.items()}}
        record["optimizer"] = opt
        record["optimizer_history"] = tracker.history          # every residual evaluation (cost, seconds)
        record["initial"] = s_init
        record["final"] = s_final
        record["fit_state"] = s_fit_summary
        engine.close()
        engine = None
        utc = utc_now()
        meta = dict(fit_identity, parameterization=PARAMETERIZATION_ID, policy=POLICY, created_utc=utc,
                    mesh={"identity": ctx.mesh_identity, "faces": obj.nf}, mesh_identity=ctx.mesh_identity,
                    lam=obj.lam, floor_ratio=obj.floor_ratio, environment=env["process"],
                    facets=record["facets"], late_only=seen_rec.get("late_only"), sensitivity=sens,
                    config={"seen_points": args.seen_points, "parameters": obj.n, "parameter_facets": obj.npar,
                            "lsmr_maxiter": args.lsmr_maxiter, "max_nfev": args.max_nfev},
                    costs={"initial": {k: s_init[k] for k in RESULT_KEYS[:4]},
                           "final": {k: s_final[k] for k in RESULT_KEYS[:4]}, "fit": ctx.fit_costs},
                    termination=opt["termination"], stop=reason, complete=reason != "wall_time",
                    nfev=opt["nfev"], nfev_total=seg["nfev_total"], segment=args.segment, segments=segments,
                    lower_bound_note="lower_bound is the bound on delta; -inf: none (A_orig = 0); NaN: not a "
                                     "parameter (unseen); held_next = held_prev + delta; s_final = s_fit + held_next",
                    resume=resume)
        path, chk = write_state(out, ctx, obj, st, meta)
        record["products"] = {"finetune_state.npz": {"path": str(path), "s_final_minus_s_fit_minus_held_next_max_abs":
                                                     chk}}
        write_tables(out, ctx, obj, st_fit, st, utc)
        start_label = ("the tuned fit (delta = 0, the fit's grain)" if resume is None else
                       "the resumed state %s" % resume["path"])
        write_result_txt(out, ctx, obj, s_init, s_final, opt, seen_rec, utc, start_label, seg, sens)
        for name in ("finetune_state.npz", "finetune_result.txt", "COPS_model.dat", "late_window_model.dat",
                     "pressure_comparison.txt"):
            record["products"].setdefault(name, {}).update(file_identity(out / name))
        log("products in %s" % out)
        if not bounds_ok:
            log("ERROR the final state violates a bound")
            return 3
        if reason == "wall_time":
            log("stopped at the wall-time limit (segment %d, %d evaluations in all): resubmit the same sbatch command "
                "(same FIT, OUT, MAXNFEV): segment %d continues from %s; exit 75" % (
                    args.segment, seg["nfev_total"], args.segment + 1, out / "finetune_state.npz"))
            return 75
        log("complete: stop %s (%s), %d evaluations in all" % (reason, opt["termination"], seg["nfev_total"]))
        return 0
    finally:
        record["finished_utc"] = utc_now()
        record["wall_seconds"] = time.time() - started
        if engine is not None:
            try:
                record.setdefault("stats", {"workers": engine.stats(), "driver_max_rss_gb": max_rss_gb()})
            except Exception as exc:
                record.setdefault("stats", {"error": repr(exc)})
            engine.close()
        tmp = rec_path.with_name(rec_path.name + ".tmp")
        with tmp.open("w") as f:
            json.dump(jsonable(record), f, indent=1)
        if rec_path.exists():
            os.replace(str(tmp), str(rec_path.with_name("fine_tuning.json.duplicate_%d" % os.getpid())))
        else:
            os.replace(str(tmp), str(rec_path))


if __name__ == "__main__":
    sys.exit(main())
