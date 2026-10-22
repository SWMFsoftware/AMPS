"""Late-interval activity fine tuning on the GPU: per-group multipliers and additive floors.

A second fitting stage on a completed fit (``fit_restart.npz``). The definitions of the support
are those of ``../LiouvilleCOPS/FINE_TUNING_STEP.md`` (CPU tool ``tools/fine_tune_activity.py``);
the correction law is the user's (decision 2026-10-20):

* **Window and intervals:** the last ``m`` minutes, ``n`` equal intervals, nearest-sample time
  weights. The CPU function ``fine_tuning.intervals`` itself is called, so the edges, bins and
  weights are the CPU's.
* **First structural response:** a facet belongs to the first interval in which any of Kn, KPhi,
  KV is nonzero, provided it was never nonzero at an earlier row (``kernels/restart.cu:
  first_structural_row``). Facets first seen before the window are ineligible.
* **Dominance:** within its first interval, the weighted sums of |H_NG| and |H_RG| (gauge
  pressure per unit clipped activity at the frozen temperatures) classify a facet as
  NG-dominated (R_RG <= eps R_NG), RG-dominated (R_NG <= eps R_RG) or mixed (CPU inequalities,
  ``kernels/restart.cu: first_interval_scores``).
* **Groups:** the disjoint (interval, NG/RG/mixed) atoms (the CPU's ``--third-group mixed``).

**Correction law** (differs from the CPU tool, which only adds activity). On group g:

    A_final(i) = f_g * A0(i) + c_g * A_scale,     f_g >= 0,  c_g >= 0,

with ``A0 = max(activity_linear, 0)`` the frozen fitted activity and ``A_scale`` the restart's
enhancement scale. ``f_g`` is a multiplier (f < 1 lowers the activity, f = 1 is no change);
``c_g`` is an additive floor, the only way to change facets whose fitted activity is clipped
to zero. Facets outside every group keep A0.

**Exact reduced problem.** The pressures are linear in each facet's clipped activity, and
A_final >= 0, so no new clip occurs:

    p = p0 + sum_g (f_g - 1) * WA_g + c_g * A_scale * W1_g,
    WA_g = sum_{i in g} H(i) A0(i),   W1_g = sum_{i in g} H(i),

per row and gauge. ``WA_g`` and ``W1_g`` are computed with the validated evaluator
(``restart.StateEvaluator.point_changes``, deterministic block sums), so they include every
later row of the group's facets. Rows before an interval have exactly zero response from its
groups (asserted), so every earlier pressure is bitwise unchanged.

**Objective** (the CPU's): log residuals of the selected positive late NG and RG references,
equal total weight per gauge; optional ``lambda`` times the area-weighted squared activity
change relative to the frozen power ``sum area A0^2``. Bounded TRF, analytic Jacobian,
``x_scale='jac'``, start f = 1, c = 0. Columns with no response on the selected data are fixed
at f = 1 or c = 0 and reported.

Several GPUs: ``SplitEvaluator`` loads contiguous row blocks with balanced record counts, one
per device. Point results depend only on the point's own records, so they are bitwise those of
one device; the first-interval scores are accumulated row by row in chronological order,
carried from device to device, so their summation order is the CPU's.
"""

import numpy as np

from . import formats, reference
from .memory import scratch_reservation
from .restart import StateEvaluator, record_factor_args

SCHEMA = "liouville_cops_cuda.late_interval_fine_tuning_multiplier.v1"
NEVER, BEFORE = -2, -1
KINDS = ("ng", "rg", "mixed")
UNSEEN = np.uint32(0xFFFFFFFF)


# ------------------------------------------------------------------ evaluation on 1..N GPUs
def record_counts(run, rows, cols, nfaces):
    """Record count per table row from the kernel headers (no records read)."""
    from pathlib import Path
    counts = np.empty(len(rows), dtype=np.int64)
    for k, row in enumerate(rows):
        point = int(row[cols["iPoint"]])
        path = Path(run) / "kernels" / formats.kernel_file_name(point)
        counts[k] = formats.read_kernel_header(path, point, nfaces)[1]
    return counts


class SplitEvaluator:
    """``StateEvaluator`` over contiguous row blocks, one per device; results in table order."""

    def __init__(self, run, devices):
        from pathlib import Path
        from . import multi_gpu
        devices = [int(d) for d in devices]
        if not devices:
            raise ValueError("no device given")
        run = Path(run).resolve()
        fit_cops = reference.load("fit_cops")
        _, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
        if len(devices) == 1:
            self.blocks = [(0, len(rows))]
        else:
            nfaces = formats.read_element_header(run / "elements.bin")[0]
            counts = record_counts(run, rows, {n: i for i, n in enumerate(names)}, nfaces)
            self.blocks = multi_gpu.partition(counts, len(devices))
        self.parts = []
        try:
            for (a, b), d in zip(self.blocks, devices):
                self.parts.append(StateEvaluator(run, d, row_indices=np.arange(a, b)))
        except BaseException:
            self.release()
            raise
        first = self.parts[0]
        self.run, self.rows, self.cols = first.run, first.rows, first.cols
        self.elements, self.coef, self.header = first.elements, first.coef, first.header
        self.mesh_identity = first.mesh_identity
        self.devices = devices
        self.point_ids = np.concatenate([p.point_ids for p in self.parts])
        self.nfaces = len(self.elements["area"])

    def enhancement_scale(self):
        return self.parts[0].enhancement_scale()

    def pressures(self, restart):
        """pNG, pRG of the restart state over all table rows (binds every part)."""
        out = [p.pressures(restart) for p in self.parts]
        return np.concatenate([o[0] for o in out]), np.concatenate([o[1] for o in out])

    def point_changes(self, restart, delta_af):
        out = [p.point_changes(restart, delta_af) for p in self.parts]
        return np.concatenate([o[0] for o in out]), np.concatenate([o[1] for o in out])

    def kernel_identity(self):
        """Combined kernel-set digest of the loaded parts and its check against run_manifest.json."""
        entries = [e for p in self.parts for e in p._backend.dataset.kernel_entries]
        digest = formats.kernel_set_digest(entries)
        manifest = self.parts[0]._backend.dataset.inputs.manifest or {}
        recorded = manifest.get("outputs", {}).get("kernels", {}).get("set_sha256")
        if recorded is not None and recorded != digest:
            raise ValueError("kernel set sha256 %s differs from run_manifest.json %s" % (digest, recorded))
        return {"set_sha256": digest, "files": len(entries),
                "verified_against_run_manifest": recorded is not None,
                "records": int(sum(p._backend.dataset.nrecords for p in self.parts))}

    def first_structural_row(self):
        """Per facet: the first table row with Kn, KPhi or KV != 0; -1 if none (int64)."""
        import cupy
        from . import cuda_sources
        best = np.full(self.nfaces, -1, dtype=np.int64)
        for part in self.parts:
            be = part._backend
            ds = be.dataset
            with cupy.cuda.Device(part.device_id), scratch_reservation(
                    be.budget, "first structural row", {"first": (ds.nfaces_device, np.uint32)}):
                first = cupy.full(ds.nfaces_device, UNSEEN, dtype=cupy.uint32)
                name = "first_structural_row_f32" if getattr(be, "mixed", False) else "first_structural_row"
                fn = cuda_sources.function("restart", name, part.device_id)
                for g0, g1, table, r in be.record_chunks():         # atomicMin: any order
                    fn((g1 - g0,), (256,), (table, r["face"], r["kn"], r["kphi"], r["kv"], first))
                local = ds.scatter_faces(first.get(), UNSEEN)
                del first
            # parts are in row order: the first part that sees a facet holds its earliest row
            new = (local != UNSEEN) & (best < 0)
            best[new] = part.row_indices[local[new].astype(np.int64)]
        return best

    def max_shares(self, restart, delta_af, row_mask):
        """Per facet: max over the rows in ``row_mask`` of the relative NG or RG pressure change
        caused by ``delta_af`` of that facet alone (``StateEvaluator.facet_point_shares``)."""
        out = np.zeros(self.nfaces)
        for part in self.parts:
            flagged = part.point_ids[row_mask[part.row_indices]]
            if not len(flagged):
                continue
            sn, sr, _, _ = part.facet_point_shares(restart, delta_af, flagged)
            out = np.maximum(out, np.maximum(sn, sr))
        return out

    def first_significant_row(self, restart, delta_af, tau, row_mask=None):
        """Per facet: the first table row (among ``row_mask``) at which ``delta_af`` alone changes
        NG or RG by more than ``tau`` (relative); -1 if none (``restart.cu: first_significant_row``)."""
        import cupy
        from . import cuda_sources
        best = np.full(self.nfaces, -1, dtype=np.int64)
        Tmax, T75, Tmin = (float(v) for v in restart.T)
        for part in self.parts:
            png, prg = part.pressures(restart)
            be = part._backend
            c = be.constants
            ds = be.dataset
            nf = ds.nfaces_device
            with cupy.cuda.Device(part.device_id), scratch_reservation(
                    be.budget, "first significant row", {"first": (nf, np.uint32), "dAF": (nf, np.float64)}):
                first = cupy.full(nf, UNSEEN, dtype=cupy.uint32)
                d_daf = cupy.asarray(ds.gather_faces(np.asarray(delta_af, dtype=np.float64)))
                d_png, d_prg = cupy.asarray(png), cupy.asarray(prg)
                d_crg = cupy.asarray(np.asarray(c.cRG["H2O"], dtype=np.float64))
                use = (np.ones(len(part.row_indices), np.uint8) if row_mask is None
                       else np.asarray(row_mask[part.row_indices], dtype=np.uint8))
                d_use = cupy.asarray(use)
                fn = None
                for g0, g1, table, records in be.record_chunks():   # atomicMin: any order
                    suffix, factors = record_factor_args(be, records)
                    fn = fn or cuda_sources.function("restart", "first_significant_row" + suffix, part.device_id)
                    fn((g1 - g0,), (256,),
                       (table, d_use, records["face"]) + factors + (
                        d_daf, d_png, d_prg, d_crg, np.float64(c.cNG["H2O"] * c.bm["H2O"]), np.float64(c.bm["H2O"]),
                        np.float64(Tmax), np.float64(T75), np.float64(Tmin), np.float64(c.KT), np.float64(tau),
                        first))
                local = ds.scatter_faces(first.get(), UNSEEN)
                del first, d_daf, d_png, d_prg, d_crg, d_use
            new = (local != UNSEEN) & (best < 0)
            best[new] = part.row_indices[local[new].astype(np.int64)]
        return best

    def first_interval_scores(self, restart, first_interval, bins, weights):
        """Per facet: sum over the rows of its first interval of w |H_NG| and w |H_RG| [Pa s]."""
        import cupy
        from . import cuda_sources
        ng = np.zeros(self.nfaces)
        rg = np.zeros(self.nfaces)
        flags_all = 0
        Tmax, T75, Tmin = (float(v) for v in restart.T)
        for part in self.parts:
            be = part._backend
            c = be.constants
            local_rows = np.flatnonzero(bins[part.row_indices] >= 0)
            if not len(local_rows):
                continue
            ds = be.dataset
            nf = ds.nfaces_device
            first_block = np.asarray(be._first_block, dtype=np.int64)
            with cupy.cuda.Device(part.device_id), scratch_reservation(
                    be.budget, "first interval scores", {"first": (nf, np.int32), "ng": (nf, np.float64),
                                                         "rg": (nf, np.float64), "stamp": (nf, np.uint32)}):
                d_first = cupy.asarray(ds.gather_faces(np.asarray(first_interval, dtype=np.int32)))
                d_ng, d_rg = cupy.asarray(ds.gather_faces(ng)), cupy.asarray(ds.gather_faces(rg))
                stamp = cupy.zeros(nf, dtype=cupy.uint32)
                flags = cupy.zeros(1, dtype=cupy.uint32)
                crg = np.asarray(c.cRG["H2O"], dtype=np.float64)
                fn = None
                # Chronological: chunks in block order, rows in table order within a chunk. A facet
                # has one record per point, so its sum runs over the points in table order also when
                # a point spans two chunks.
                for g0, g1, table, records in be.record_chunks():
                    suffix, factors = record_factor_args(be, records)
                    fn = fn or cuda_sources.function("restart", "first_interval_scores" + suffix, part.device_id)
                    for p in local_rows:
                        b0, b1 = max(int(first_block[p]), g0), min(int(first_block[p + 1]), g1)
                        if b1 <= b0:
                            continue
                        g = part.row_indices[p]
                        fn((b1 - b0,), (256,),
                           (table, np.int64(b0 - g0), records["face"]) + factors + (
                            d_first, np.int32(bins[g]), np.float64(weights[g]),
                            np.float64(crg[p]), np.float64(c.cNG["H2O"] * c.bm["H2O"]), np.float64(c.bm["H2O"]),
                            np.float64(Tmax), np.float64(T75), np.float64(Tmin), np.float64(c.KT),
                            d_ng, d_rg, stamp, flags))
                ng = ds.scatter_faces(d_ng.get(), out=ng.copy())
                rg = ds.scatter_faces(d_rg.get(), out=rg.copy())
                flags_all |= int(flags.get()[0])
                del d_first, d_ng, d_rg, stamp, flags
        if flags_all & 1:
            raise ValueError("a facet occurs twice within one point's kernel records")
        if flags_all & 2:
            raise ValueError("negative gauge response in the window: review the base pressure formulation")
        return ng, rg

    def release(self):
        for p in getattr(self, "parts", []):
            p.release()


# ------------------------------------------------------------------ support and groups
def reference_change(frozen, scale):
    """Per-facet reference activity change for significance: max(A0, A_scale), i.e. doubling
    the facet (f = 2) or one additive A_scale (c = 1), whichever is larger."""
    return np.maximum(frozen, scale)


def confine_total(evaluator, restart, delta_ref, baseline, bins, epsilon):
    """--confine total: the largest per-facet threshold tau whose eligible facets, changed by
    their reference change together, keep every row before the window within epsilon.

    Candidates are the facets with a nonzero response in the window. A candidate is eligible at
    tau if its own reference change alters no pre-window row by more than tau (relative, NG or
    RG). The summed change of the eligible set is exact (pressures are linear in the clipped
    activity) and grows with tau (all reference changes are increases), so the largest
    admissible tau is found by binary search over the facets' own shares, exactly as
    ``restart_tool.py enhance --confine total``.
    Returns (tau, eligible mask, record).
    """
    before, window = bins < 0, bins >= 0
    pre = evaluator.max_shares(restart, delta_ref, before)
    win = evaluator.max_shares(restart, delta_ref, window)
    candidate = win > 0
    p0 = baseline[before]

    def total(tau):
        keep = candidate & (pre <= tau)
        if not keep.any():
            return 0.0, keep
        dn, dr = evaluator.point_changes(restart, np.where(keep, delta_ref, 0.0))
        d = np.column_stack([dn, dr])[before]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(np.abs(p0) > 0, np.abs(d) / np.abs(p0), np.where(d != 0, np.inf, 0.0))
        return float(rel.max(initial=0.0)), keep

    levels = np.unique(pre[candidate])
    lo, hi = -1, len(levels) - 1
    trace = []
    while lo < hi:
        mid = (lo + hi + 1) // 2
        value, keep = total(levels[mid])
        trace.append({"tau": float(levels[mid]), "facets": int(keep.sum()), "outside_max_rel": value})
        if value <= epsilon:
            lo = mid
        else:
            hi = mid - 1
    tau = float(levels[lo]) if lo >= 0 else -1.0
    value, keep = total(tau) if lo >= 0 else (0.0, np.zeros_like(candidate))
    record = {"epsilon": epsilon, "tau": tau, "candidates_seen_in_window": int(candidate.sum()),
              "eligible": int(keep.sum()), "reference_outside_max_rel": value,
              "reference_change": "max(A0, A_scale) per facet", "search": trace}
    return tau, keep, record


def first_intervals(first_row, bins):
    """first_interval per facet: -2 never seen, -1 first seen before the window, else bin."""
    out = np.full(len(first_row), NEVER, dtype=np.int32)
    seen = first_row >= 0
    out[seen] = bins[first_row[seen]]
    return out


def classify(first_interval, ng_score, rg_score, epsilon):
    """atom = 3 j + kind (0 NG, 1 RG, 2 mixed) or -1, with the CPU's inequalities."""
    if not np.isfinite(epsilon) or not 0 < epsilon < 1:
        raise ValueError("dominance ratio must be strictly between zero and one")
    has_signal = (ng_score > 0) | (rg_score > 0)
    kind = np.full(len(first_interval), 2, dtype=np.int8)
    kind[(ng_score > 0) & (rg_score <= epsilon * ng_score)] = 0
    kind[(rg_score > 0) & (ng_score <= epsilon * rg_score)] = 1
    return np.where((first_interval >= 0) & has_signal, 3 * first_interval + kind, -1).astype(np.int32)


def threshold_near(first_interval, ng_score, rg_score, epsilon, rtol=1e-9):
    """Facets whose dominance decision would change under a relative score change of rtol."""
    lo, hi = epsilon * (1 - rtol), epsilon * (1 + rtol)
    ratio_rg = np.where(ng_score > 0, rg_score / np.where(ng_score > 0, ng_score, 1), np.inf)
    ratio_ng = np.where(rg_score > 0, ng_score / np.where(rg_score > 0, rg_score, 1), np.inf)
    near = ((ratio_rg >= lo) & (ratio_rg <= hi)) | ((ratio_ng >= lo) & (ratio_ng <= hi))
    return np.flatnonzero(near & (first_interval >= 0))


def group_names(atom, count):
    """Nonempty atoms in (interval, kind) order with their names."""
    present = np.unique(atom[atom >= 0])
    return [(int(k), "late_%03d_%s" % (k // 3, KINDS[k % 3])) for k in present]


# ------------------------------------------------------------------ the reduced problem
def response_matrices(evaluator, restart, atom, groups, frozen):
    """WA[row, gauge, g] and W1[row, gauge, g] over all table rows (exact block sums)."""
    nrow = len(evaluator.point_ids)
    WA = np.zeros((nrow, 2, len(groups)))
    W1 = np.zeros((nrow, 2, len(groups)))
    for k, (code, _) in enumerate(groups):
        member = atom == code
        WA[:, 0, k], WA[:, 1, k] = evaluator.point_changes(restart, np.where(member, frozen, 0.0))
        W1[:, 0, k], W1[:, 1, k] = evaluator.point_changes(restart, member.astype(np.float64))
    return WA, W1


def parameters(groups, WA, W1, scale):
    """Parameter list (name, group index, 'f'|'c'), the design D[row, gauge, K], reference x."""
    names, design, xref, kinds, index = [], [], [], [], []
    for k, (_, name) in enumerate(groups):
        names += ["f_" + name, "c_" + name]
        design += [WA[:, :, k], scale * W1[:, :, k]]
        xref += [1.0, 0.0]
        kinds += ["f", "c"]
        index += [k, k]
    D = np.stack(design, axis=-1) if design else np.zeros(WA.shape[:2] + (0,))
    return names, D, np.asarray(xref), kinds, np.asarray(index, dtype=int)


def gram_matrix(groups, atom, area, frozen, scale, power, kinds, index):
    """Area-weighted squared activity change, relative to the frozen power, as z^T G z with
    z = x - xref: on group g, dA = (f-1) A0 + c A_scale."""
    K = len(kinds)
    G = np.zeros((K, K))
    for k, (code, _) in enumerate(groups):
        m = atom == code
        saa = float(np.dot(area[m], frozen[m] ** 2))
        sa1 = float(np.dot(area[m], frozen[m])) * scale
        s11 = float(np.sum(area[m])) * scale ** 2
        jf, jc = 2 * k, 2 * k + 1
        G[jf, jf], G[jf, jc], G[jc, jf], G[jc, jc] = saa, sa1, sa1, s11
    return G / power


def fit(baseline, design, xref, refs, used, late, gram, lam=0.0, max_nfev=1000, confine=None, gauges=(0, 1)):
    """Bounded joint fit of all multipliers and floors to the pressures of ``gauges`` (0 NG, 1 RG):
    log(p_model/p_data) per selected late point, each pressure on its own (never the NG/RG
    ratio), equal total weight per fitted gauge.

    ``confine = (rows, epsilon)``: the summed fitted change must keep every row in ``rows``
    (the rows before the window) within ``epsilon`` relative, both gauges. The unconstrained
    TRF solution is kept if it satisfies this; otherwise the problem is solved again with these
    linear inequalities (``trust-constr``, from the feasible start f = 1, c = 0)."""
    from scipy.optimize import least_squares
    if not np.isfinite(lam) or lam < 0 or max_nfev < 1:
        raise ValueError("lambda must be finite/nonnegative and max-nfev positive")
    fitted = np.zeros(2, bool)
    fitted[list(gauges)] = True
    selected = used & late[:, None] & (refs > 0) & fitted[None, :]
    counts = selected.sum(axis=0)
    if np.any(counts[fitted] == 0):
        raise ValueError("fine tuning requires selected positive late observations of every fitted gauge")
    row, gauge = np.nonzero(selected)
    d, p, target = design[row, gauge], baseline[row, gauge], refs[row, gauge]
    if np.any(p <= 0) or not np.isfinite(p).all():
        raise ValueError("selected baseline pressure must be finite and positive")
    weights = 1 / np.sqrt(counts[gauge])
    eig, vec = np.linalg.eigh(gram) if gram.size else (np.zeros(0), np.zeros((0, 0)))
    penalty = np.sqrt(lam) * np.sqrt(np.maximum(eig, 0))[:, None] * vec.T

    def residual(x):
        z = x - xref
        return np.r_[weights * np.log((p + d @ z) / target), penalty @ z]

    def jacobian(x):
        z = x - xref
        return np.vstack([weights[:, None] * d / (p + d @ z)[:, None], penalty])

    active = np.flatnonzero(np.any(d != 0, axis=0))
    x = xref.copy()
    if len(active):
        def expand(v):
            out = xref.copy()
            out[active] = v
            return out
        result = least_squares(lambda v: residual(expand(v)), xref[active],
                               jac=lambda v: jacobian(expand(v))[:, active],
                               bounds=(0, np.inf), x_scale="jac", max_nfev=max_nfev,
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        x[active] = result.x
        status = dict(success=bool(result.success), status=int(result.status), message=str(result.message),
                      nfev=int(result.nfev), optimality=float(result.optimality))
    else:
        status = dict(success=True, status=0, message="no identifiable late-response columns", nfev=0)
    if confine is not None and len(active):
        rows, epsilon = confine
        A = design[rows][:, :, active].reshape(-1, len(active))
        p_out = np.abs(baseline[rows]).reshape(-1)
        keep = np.any(A != 0, axis=1) & (p_out > 0)
        A, p_out = A[keep] / (epsilon * p_out[keep])[:, None], p_out[keep]

        def outside(v):
            return float(np.max(np.abs(A @ (v - xref[active])), initial=0.0)) * epsilon
        status["outside_max_rel_unconstrained"] = outside(x[active])
        status["outside_constraint_rows"] = int(len(A))
        if status["outside_max_rel_unconstrained"] > epsilon:
            from scipy.optimize import Bounds, LinearConstraint, minimize
            xa = xref[active]

            def objective(v):
                rr = residual(expand(v))
                return 0.5 * float(rr @ rr), jacobian(expand(v))[:, active].T @ rr
            # SLSQP (fast for a few tens of variables); trust-constr if it fails. Both start
            # from the feasible point f = 1, c = 0; the result is accepted only if feasible.
            ineq = {"type": "ineq", "fun": lambda v: np.r_[1.0 - A @ (v - xa), 1.0 + A @ (v - xa)],
                    "jac": lambda v: np.vstack([-A, A])}
            res = minimize(objective, xa, jac=True, method="SLSQP", constraints=[ineq],
                           bounds=[(0.0, None)] * len(active),
                           options={"maxiter": max(1000, 20 * max_nfev), "ftol": 1e-14})
            # The reduced problem is badly conditioned, so SLSQP can stop early: polish with
            # trust-constr from the better of (SLSQP result if feasible, feasible start) and keep
            # the lowest-cost feasible point.
            feasible = lambda v: outside(np.maximum(v, 0.0)) <= epsilon * (1 + 1e-6)
            cost = lambda v: objective(np.maximum(v, 0.0))[0]
            candidates = [("start", xa.copy())]
            if feasible(res.x):
                candidates.append(("SLSQP", np.maximum(res.x, 0.0)))
            start_name, start = min(candidates, key=lambda c: cost(c[1]))
            con = LinearConstraint(A, -1.0 + A @ xa, 1.0 + A @ xa)
            res2 = minimize(objective, start, jac=True, method="trust-constr", constraints=[con],
                            bounds=Bounds(np.zeros(len(active)), np.full(len(active), np.inf)),
                            options={"maxiter": max(1000, 20 * max_nfev), "gtol": 1e-12, "xtol": 1e-14,
                                     "verbose": 0})
            if feasible(res2.x):
                candidates.append(("trust-constr from " + start_name, np.maximum(res2.x, 0.0)))
            method, best = min(candidates, key=lambda c: cost(c[1]))
            x[active] = best
            status["constrained_method"] = method
            status["constrained_candidates"] = {name: cost(v) for name, v in candidates}
            res = res2 if method.startswith("trust") else res
            status.update(constrained=True, constrained_message=str(res.message),
                          constrained_success=bool(res.success), constrained_niter=int(res.nit))
        else:
            status["constrained"] = False
        status["outside_max_rel"] = outside(x[active])
    jd = jacobian(x)[:len(p)][:, active] if len(active) else np.zeros((len(p), 0))
    singular = np.linalg.svd(jd, compute_uv=False) if jd.size else np.zeros(0)
    tol = (np.finfo(float).eps * max(jd.shape) * singular[0]) if len(singular) else 0.0
    status.update(x=x.tolist(), active_columns=active.tolist(),
                  zero_response_columns=np.setdiff1d(np.arange(len(x)), active).tolist(),
                  at_lower_bound=np.flatnonzero((x == 0) & np.isin(np.arange(len(x)), active)).tolist(),
                  data_singular_values=singular.tolist(),
                  data_rank=int(np.count_nonzero(singular > tol)),
                  condition_number=float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0
                  else float("inf"),
                  selected_counts=counts.tolist(), cost=float(np.sum(residual(x) ** 2) / 2),
                  initial_cost=float(np.sum(residual(xref) ** 2) / 2))
    prediction = baseline + np.einsum("igk,k->ig", design, x - xref)
    return x, prediction, status


def final_activity(frozen, atom, groups, x, kinds, index, scale):
    """A_final per facet: f_g A0 + c_g A_scale on group g, A0 elsewhere."""
    out = frozen.copy()
    for k, (code, _) in enumerate(groups):
        m = atom == code
        f = x[(index == k) & (np.asarray(kinds) == "f")][0]
        c = x[(index == k) & (np.asarray(kinds) == "c")][0]
        out[m] = f * frozen[m] + c * scale
    return out


def cpu_intervals(utc, minutes, count, end=None):
    """The CPU tool's window: edges, per-row bins (-1 before) and nearest-sample weights."""
    return reference.load("fine_tuning").intervals(utc, minutes, count, end)
