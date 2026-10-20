"""FP64 CUDA pressure/Jacobian backend (ROADMAP C06; CUDA environment).

``CudaBackend`` implements the C02 contract (``backend.PressureBackend``) for the
primary joint families B29/N30/C30/NC31 and the labelled extensions: a free
degree-0–4 baseline, free temperatures and κ, and nonnegative template
amplitudes. The fixed-baseline harmonic modes run through the same code,
because their basis is just more enhancement columns.

Division of work
----------------
Per evaluation, on the host (NumPy, the reference's own expressions):

* ``act_lin = Y @ a (+ A_scale * (B @ b))``, ``af = max(act_lin, 0)`` and
  ``active = act_lin > 0`` for every facet. Clipping decisions are therefore
  bitwise those of ``fit_cops.Model.evaluate``. A device matvec could not
  reproduce the BLAS rounding, and the upload is only about 9 bytes per facet.
* The final assembly of pressures and derivative columns from the per-point
  sums, with the reference's constants and expressions (``Model.evaluate``).

On the device (``kernels/pressure.cu``):

* once per bind, the static per-contribution factors of ``Model.load_kernels``
  and ``Model.prepare``: q·Kn, q·VB·KV, q·KPhi for both species, xlit, the lit
  flag, and the parameter-independent PhiCO2 sums;
* per evaluation, the per-contribution terms and their per-point sums: five
  pressure sums, nine temperature sums when T columns are requested, and two
  sums per requested activity or template column. Requests are grouped in
  tiles of up to 32 sums per launch.

Modes
-----
``"tree"`` (production)
    Fixed blocks of records per point; strided per-thread folds; a fixed
    shared-memory tree; block partials added in order. Deterministic, within
    the approved budget (VALIDATION.md).
``fused`` (C08, default in tree mode)
    For the evaluator's standard requests (pressures, plus optionally the nine
    temperature sums, plus NG/RG pairs of activity/template columns), the tree
    mode uses compile-time-specialized kernels: ``eval_fused_p5``,
    ``eval_fused_p14`` and ``eval_fused_c1`` … ``eval_fused_c16``. Their accumulators
    are held in registers (ptxas: no stack frame, no spills). They use the same
    per-term expressions, per-thread fold and tree as ``eval_partial``, so the
    results are bitwise identical (``tests/test_evaluator.py``,
    ``tools/c08_fused_check.py``); the full Jacobian is 2.6× faster. Any other
    code list uses the generic ``eval_partial``. ``fused=False`` forces the
    generic path.
``"sequential"`` (verification)
    One thread per point, folding records in file order, like ``np.bincount``
    and SciPy's CSR product. It reproduces the CPU sums bit for bit and proves
    that every per-contribution expression matches the reference exactly.

Supported configuration
-----------------------
Supported: ``SZAMode=continuous`` (the campaign's); all three RG pressure
models; both heliocentric models; all derivative kinds.
``SZAMode=srccg_table`` is refused (``UnsupportedConfiguration``). Its integer-degree
SZA binning would need host-side transcendental evaluation per contribution
to stay exact; it is added only if requested.
"""

import time

import numpy as np

from . import backend as bk
from . import cuda_sources, precision, provenance, reference
from .device_store import DeviceDataset
from .errors import IdentityError, UnsupportedConfiguration
from .memory import DeviceBudget, MemoryPlan, human

MAX_SUMS = 32                 # must equal MAX_SUMS in pressure.cu
FUSED_MAX_COLUMNS = 16        # eval_fused_c1 .. eval_fused_c16 in pressure.cu (2G <= MAX_SUMS)
BLOCK = 256                   # must equal BLOCK in pressure.cu
PRESSURE_CODES = (0, 1, 2, 3, 4)
TEMPERATURE_CODES = tuple(range(5, 14))
PHI_CO2_CODE = 14
MODES = ("tree", "sequential")
STATIC_FIELDS = ("wnH", "wvH", "wpH", "wnC", "wvC", "xlit")


class HostConstants:
    """The reference's per-point and global constants, from a ``fit_cops.Model``.

    The model may be the reference itself or the header-only subclass of
    ``host_model.py``; both carry the same per-point data (proven bitwise in
    ``tests/test_host_model.py``).

    Taking them from the CPU model guarantees identity with the reference.
    ``qmin`` and ``qmax`` are recomputed with ``Model.load_kernels``' own
    expressions, because the model keeps them only broadcast per contribution;
    the tests check them against the model's ``Q`` arrays.
    """

    def __init__(self, model):
        fit_cops = reference.load("fit_cops")
        self.settings = bk.ModelSettings.from_model(model)
        if self.settings.sza_mode != "continuous":
            raise UnsupportedConfiguration("CudaBackend supports SZAMode=continuous only; got %s"
                                           % self.settings.sza_mode)
        self.KT = float(fit_cops.KT)
        self.COS75 = float(fit_cops.COS75)
        self.bm = dict(model.bm)
        self.cNG = dict(model.cNG)
        self.cRG = {s: np.asarray(v, dtype=np.float64) for s, v in model.cRG.items()}
        # For tests only; None for a header-only model (host_model.py), whose
        # PhiCO2 is computed on the device alone.
        self.PhiCO2_cpu = (None if model.PhiCO2 is None
                           else np.asarray(model.PhiCO2, dtype=np.float64))
        # Model.load_kernels: rhA = rh / AU; Q_H2O ~ rhA**-p, Q_CO2 ~ rhA**-2
        rhA = model.rh / fit_cops.AU
        p = 4.2143229 if model.helio == "shap5" else 4.8
        self.qmin = {"H2O": 0.2 * 600 / rhA ** p, "CO2": 0.1 * 600 / rhA ** 2}
        self.qmax = {"H2O": 600 / rhA ** p, "CO2": 600 / rhA ** 2}
        self.vdotb = np.asarray(model.VdotB, dtype=np.float64)


class CudaBackend(bk.PressureBackend):
    """The C06 GPU evaluator behind the C02 contract."""

    name = "cuda_fp64_v1"

    def __init__(self, model, run, enhancement_basis=None, enhancement_scale=1.0,
                 device_id=0, mode="tree", budget=None, records_per_block=4096,
                 chunk_records=None, fused=True):
        if mode not in MODES:
            raise ValueError("mode must be one of %s" % (MODES,))
        nb = 0 if enhancement_basis is None else int(np.asarray(enhancement_basis).shape[1])
        super().__init__(model.rows, bk.ModelSettings.from_model(model), nb)
        import cupy
        self.cupy = cupy
        self.mode = mode
        # C08: compile-time-specialized tree kernels (bitwise equal to eval_partial).
        self.fused = bool(fused) and mode == "tree"
        self.device_id = int(device_id)
        self.records_per_block = int(records_per_block)
        self.constants = HostConstants(model)
        self.Y = model.Y                                   # the CPU's own array (exact matvec)
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale)
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else bk.array_digest(enhancement_basis))
        self.timing = {}
        self.activity_source = None
        self.arrays = {}
        self.reservation = None
        self.released = False
        self._columns_cache = None
        point_ids = [int(p) for p in np.asarray(model.ipoint)[model.rows]]
        self.owns_budget = budget is None
        self.budget = DeviceBudget(self.device_id) if budget is None else budget
        t0 = time.perf_counter()
        kwargs = {} if chunk_records is None else {"chunk_records": chunk_records}
        self.dataset = DeviceDataset.load(run, point_ids=point_ids, device_id=self.device_id,
                                          mesh_fields=("ACO2",), budget=self.budget, **kwargs)
        try:
            self._check_dataset(model)
            self.timing["dataset_load_s"] = time.perf_counter() - t0
            self._bind()
        except BaseException:
            self.release()
            raise

    # ------------------------------------------------------------------ binding
    def _check_dataset(self, model):
        """The device rows must be the model's rows, in the model's order."""
        geometry = self.dataset.scan.geometry
        if not np.array_equal(geometry[:, 1], self.constants.vdotb):
            raise IdentityError("device V.b differs from the model's kernel headers")
        if not np.array_equal(geometry[:, 8], model.rh):
            raise IdentityError("device heliocentric distances differ from the model's")
        if self.dataset.inputs.nfaces != len(self.Y):
            raise IdentityError("mesh size differs between the model and the device dataset")
        counts = np.diff(model.indptr)
        if not np.array_equal(counts, self.dataset.scan.counts):
            raise IdentityError("per-point record counts differ from the model's kernels")

    def _bind(self):
        cupy, ds = self.cupy, self.dataset
        npt, nrec, nf = ds.npoints, ds.nrecords, ds.inputs.nfaces
        blocks = self._block_table(ds.scan.counts, ds.scan.offsets)
        self.plan = MemoryPlan("cuda evaluator")
        for name in STATIC_FIELDS:
            self.plan.add_array("eval." + name, "static", (nrec,), np.float64)
        self.plan.add_array("eval.lit", "static", (nrec,), np.uint8)
        self.plan.add_array("eval.phiC_scratch", "scratch", (nrec,), np.float64)
        self.plan.add_array("eval.af_face", "scratch", (nf,), np.float64)
        self.plan.add_array("eval.active", "scratch", (nf,), np.uint8)
        self.plan.add_array("eval.colbasis", "basis", (max(1, 25 + self.nbasis_enhancement), nf),
                            np.float64, "requested activity/template columns")
        self.plan.add_array("eval.blocks", "static", (len(blocks),), np.int64)
        self.plan.add_array("eval.first_block", "static", (npt + 1,), np.int64)
        self.plan.add_array("eval.partial", "partial_sums", (len(blocks) // 3, MAX_SUMS), np.float64)
        self.plan.add_array("eval.sums", "scratch", (MAX_SUMS, npt), np.float64)
        self.plan.add_array("eval.per_point_inputs", "scratch", (5, npt), np.float64)
        self.reservation = self.budget.reserve(self.plan)
        t0 = time.perf_counter()
        with cupy.cuda.Device(self.device_id), self.budget.allocation_guard(self.plan):
            a = self.arrays
            for name in STATIC_FIELDS:
                a[name] = cupy.empty(nrec, dtype=np.float64)
            a["lit"] = cupy.empty(nrec, dtype=np.uint8)
            phiC = cupy.empty(nrec, dtype=np.float64)
            c = self.constants
            per_point = cupy.asarray(np.stack((c.qmin["H2O"], c.qmax["H2O"], c.qmin["CO2"],
                                               c.qmax["CO2"], c.vdotb)).astype(np.float64))
            d = ds.arrays
            if nrec:
                cuda_sources.function("pressure", "precompute_static", self.device_id)(
                    (max(1, min(256, -(-int(ds.scan.counts.max()) // 256))), npt), (256,),
                    (d["kernel.offsets"], d["kernel.face"], d["kernel.kn"], d["kernel.kphi"],
                     d["kernel.kv"], d["kernel.cos"], d["kernel.shadow"],
                     per_point[0], per_point[1], per_point[2], per_point[3], per_point[4],
                     d["mesh.ACO2"], np.float64(c.COS75),
                     a["wnH"], a["wvH"], a["wpH"], a["wnC"], a["wvC"], phiC,
                     a["xlit"], a["lit"]))
            a["blocks"] = cupy.asarray(blocks)
            a["first_block"] = cupy.asarray(self._first_block)
            a["partial"] = cupy.empty((len(blocks) // 3) * MAX_SUMS, dtype=np.float64)
            a["sums"] = cupy.empty(MAX_SUMS * npt, dtype=np.float64)
            a["af_face"] = cupy.zeros(nf, dtype=np.float64)
            a["active"] = cupy.zeros(nf, dtype=np.uint8)
            a["colbasis"] = cupy.empty(max(1, 25 + self.nbasis_enhancement) * nf, dtype=np.float64)
            # The parameter-independent CO2 flux sums (Model.prepare PhiCO2).
            a["phiC"] = phiC
            self.PhiCO2 = self._run_sums([PHI_CO2_CODE], None)[0]
            del a["phiC"], phiC, per_point
            self.cupy.get_default_memory_pool().free_all_blocks()
        self.timing["bind_s"] = time.perf_counter() - t0

    def _block_table(self, counts, offsets):
        """(point, begin, end) per fixed-size block of records; fixed for a bind."""
        rows, first = [], [0]
        for p, (count, start) in enumerate(zip(counts, offsets[:-1])):
            nblocks = max(1, -(-int(count) // self.records_per_block))
            for k in range(nblocks):
                begin = int(start) + k * self.records_per_block
                rows.append((p, begin, min(begin + self.records_per_block, int(start) + int(count))))
            first.append(first[-1] + nblocks)
        self._first_block = np.asarray(first, dtype=np.int64)
        return np.asarray(rows, dtype=np.int64).reshape(-1)

    # --------------------------------------------------------------- evaluation
    def capabilities(self):
        return {"sza_modes": {"continuous"},
                "heliocentric_models": {self.settings.heliocentric_model},
                "rg_models": {"equilibrium_sphere", "calibration_srccg", "dsmc_srccg"},
                "column_kinds": {"a", "b", "T", "kappa"}, "enhancement": True}

    def _params(self, state):
        c = self.constants
        Tmax, T75, Tmin = state.T
        kappa = state.kappa
        return {"Tmax": float(Tmax), "T75": float(T75), "Tmin": float(Tmin), "KT": c.KT,
                "cNbmH": float(c.cNG["H2O"] * c.bm["H2O"]),
                "kcNbmC": float(kappa * c.cNG["CO2"] * c.bm["CO2"]),
                "bmH": float(c.bm["H2O"]), "bmC": float(c.bm["CO2"])}

    def activity_linear(self, state):
        """Host act_lin = Y a (+ A_scale B b), the expression of Model.evaluate."""
        act_lin = self.Y @ state.a
        if state.b is not None:
            act_lin = act_lin + self.enhancement_scale * (self.enhancement_basis @ state.b)
        return act_lin

    def _upload_activity(self, state):
        """act_lin, then af and the active set (the reference's clip).

        ``activity_source`` (default None) lets a multi-GPU worker take act_lin
        from its parent (``multi_gpu.py``), so every device sees the same bytes
        that one process would compute.
        """
        act_lin = (self.activity_linear(state) if self.activity_source is None
                   else self.activity_source(state))
        self.arrays["af_face"].set(np.maximum(act_lin, 0.0))
        self.arrays["active"].set((act_lin > 0.0).astype(np.uint8))

    def _upload_columns(self, request):
        """Unmasked column basis: Y[:, j] and A_scale * B[:, k] (Model.evaluate)."""
        key = (request.basis_index, request.enhancement_index)
        if self._columns_cache == key:
            return
        parts = []
        if request.basis_index:
            parts.append(self.Y[:, list(request.basis_index)])
        if request.enhancement_index:
            parts.append(self.enhancement_scale *
                         self.enhancement_basis[:, list(request.enhancement_index)])
        if parts:
            matrix = np.ascontiguousarray(np.column_stack(parts).T)
            self.arrays["colbasis"][:matrix.size].set(matrix.reshape(-1))
        self._columns_cache = key

    def _fused_plan(self, codes):
        """Launches of the specialized kernels covering ``codes``, or None.

        Applies when the request is the evaluator's standard set: the five
        pressure codes, optionally the nine temperature codes, and NG/RG pairs
        of activity/template columns (codes 1000+c and 2000+c, both present).
        Returns [(kernel, col0, [codes in the kernel's output order]), ...].
        Any other list (for example PhiCO2 alone, or profiling lists) uses the
        generic kernel.
        """
        present = set(codes)
        if len(present) != len(codes):
            return None
        pressure, temperature = set(PRESSURE_CODES), set(TEMPERATURE_CODES)
        if not pressure <= present:
            return None
        rest = present - pressure
        with_t = temperature <= rest
        rest -= temperature if with_t else set()
        ng = sorted(c - 1000 for c in rest if 1000 <= c < 2000)
        rg = sorted(c - 2000 for c in rest if 2000 <= c < 3000)
        if ng != rg or len(ng) + len(rg) != len(rest):
            return None
        plan = [("eval_fused_p14", 0, list(PRESSURE_CODES) + list(TEMPERATURE_CODES)) if with_t
                else ("eval_fused_p5", 0, list(PRESSURE_CODES))]
        # Consecutive column runs, in chunks of at most FUSED_MAX_COLUMNS.
        k = 0
        while k < len(ng):
            start = k
            while (k + 1 < len(ng) and ng[k + 1] == ng[k] + 1 and
                   k + 1 - start < FUSED_MAX_COLUMNS):
                k += 1
            cols = ng[start:k + 1]
            plan.append(("eval_fused_c%d" % len(cols), cols[0],
                         [1000 + c for c in cols] + [2000 + c for c in cols]))
            k += 1
        return plan

    def _run_fused(self, plan, codes, params):
        """Run a fused plan; returns (ncodes, npt) in the order of ``codes``."""
        cupy, a, ds = self.cupy, self.arrays, self.dataset
        npt, nf = ds.npoints, ds.inputs.nfaces
        order = ("Tmax", "T75", "Tmin", "KT", "cNbmH", "kcNbmC", "bmH", "bmC")
        P = cupy.asarray(np.array([params[k] for k in order], dtype=np.float64))
        nblocks = int(self._first_block[-1])
        row = {code: i for i, code in enumerate(codes)}
        out = np.empty((len(codes), npt), dtype=np.float64)
        for name, col0, launch_codes in plan:
            n = len(launch_codes)
            partial = a["partial"][:nblocks * n]
            cuda_sources.function("pressure", name, self.device_id)(
                (nblocks,), (BLOCK,),
                (a["blocks"], ds.arrays["kernel.face"], a["wnH"], a["wvH"], a["wpH"], a["wnC"],
                 a["wvC"], a["xlit"], a["lit"], a["af_face"], a["active"], a["colbasis"],
                 np.int64(nf), np.int32(col0), P, partial))
            sums = a["sums"][:n * npt]
            cuda_sources.function("pressure", "eval_finish", self.device_id)(
                (-(-(n * npt) // 256),), (256,),
                (a["first_block"], np.int64(npt), partial, np.int32(n), sums))
            host = sums.get().reshape(n, npt)
            for i, code in enumerate(launch_codes):
                out[row[code]] = host[i]
        return out

    def _run_sums(self, codes, params):
        """Per-point sums for the given codes, in tiles of MAX_SUMS; returns (ncodes, npt)."""
        if self.fused and params is not None:
            plan = self._fused_plan(codes)
            if plan is not None:
                return self._run_fused(plan, codes, params)
        cupy, a, ds = self.cupy, self.arrays, self.dataset
        npt, nf = ds.npoints, ds.inputs.nfaces
        out = np.empty((len(codes), npt), dtype=np.float64)
        order = ("Tmax", "T75", "Tmin", "KT", "cNbmH", "kcNbmC", "bmH", "bmC")
        values = np.zeros(len(order), dtype=np.float64)
        if params is not None:
            values[:] = [params[k] for k in order]
        P = cupy.asarray(values)            # copied bit for bit (load_params in pressure.cu)
        phi = a.get("phiC", a["wnH"])          # only read for code 14
        for start in range(0, len(codes), MAX_SUMS):
            tile = np.asarray(codes[start:start + MAX_SUMS], dtype=np.int32)
            n = len(tile)
            dev_codes = cupy.asarray(tile)
            common = (ds.arrays["kernel.face"], a["wnH"], a["wvH"], a["wpH"], a["wnC"], a["wvC"], phi,
                      a["xlit"], a["lit"], a["af_face"], a["active"], a["colbasis"],
                      np.int64(nf), dev_codes, np.int32(n), P)
            sums = a["sums"][:n * npt]
            if self.mode == "sequential":
                cuda_sources.function("pressure", "eval_sequential", self.device_id)(
                    (-(-npt // 64),), (64,),
                    (ds.arrays["kernel.offsets"], np.int64(npt)) + common + (sums,))
            else:
                nblocks = int(self._first_block[-1])
                partial = a["partial"][:nblocks * n]
                cuda_sources.function("pressure", "eval_partial", self.device_id)(
                    (nblocks,), (BLOCK,), (a["blocks"],) + common + (partial,))
                cuda_sources.function("pressure", "eval_finish", self.device_id)(
                    (-(-(n * npt) // 256),), (256,),
                    (a["first_block"], np.int64(npt), partial, np.int32(n), sums))
            out[start:start + n] = sums.get().reshape(n, npt)
        return out

    def _evaluate(self, state, request):
        t0 = time.perf_counter()
        c = self.constants
        with self.cupy.cuda.Device(self.device_id):
            self._upload_activity(state)
            codes = list(PRESSURE_CODES)
            if request.temperature:
                codes += list(TEMPERATURE_CODES)
            ncol = len(request.basis_index) + len(request.enhancement_index)
            if ncol:
                self._upload_columns(request)
                codes += [1000 + k for k in range(ncol)] + [2000 + k for k in range(ncol)]
            sums = self._run_sums(codes, self._params(state))
        S = dict(zip(codes, sums))
        kappa = state.kappa
        cNG, cRG, bm = c.cNG, c.cRG, c.bm
        # Model.evaluate, expression by expression.
        pNG_H = cNG["H2O"] * bm["H2O"] * S[0]
        pRG_H = cRG["H2O"] * (S[1] + bm["H2O"] * S[2])
        pNG_C = cNG["CO2"] * bm["CO2"] * S[3]
        pRG_C = cRG["CO2"] * (self.PhiCO2 + bm["CO2"] * S[4])
        pNG = pNG_H + kappa * pNG_C
        pRG = pRG_H + kappa * pRG_C
        columns = {}
        if ncol:
            colNG = np.column_stack([S[1000 + k] for k in range(ncol)])
            colRG = np.column_stack([S[2000 + k] for k in range(ncol)])
            JNG = cNG["H2O"] * bm["H2O"] * colNG
            JRG = cRG["H2O"][:, None] * colRG
            keys = ([("a", j) for j in request.basis_index] +
                    [("b", k) for k in request.enhancement_index])
            for k, key in enumerate(keys):
                columns[key] = (JNG[:, k], JRG[:, k])
        if request.temperature:
            def triple(w):
                s_l, s_lx, s_n = S[5 + 3 * w], S[6 + 3 * w], S[7 + 3 * w]
                return c.KT * s_lx + s_l, -c.KT * s_lx, s_n
            tN, tRH, tRC = triple(0), triple(1), triple(2)
            for k in range(3):
                columns[("T", k)] = (tN[k], cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k])
        if request.kappa:
            columns[("kappa", 0)] = (pNG_C, pRG_C)
        self.timing["last_evaluate_s"] = time.perf_counter() - t0
        return bk.EvaluationResult(pNG, pRG, columns, request)

    # ------------------------------------------------------------- bookkeeping
    def identity(self):
        return {"backend": self.name, "contract_version": bk.CONTRACT_VERSION,
                "mode": self.mode, "fused_kernels": self.fused, "precision": precision.POLICY_NAME,
                "settings": self.settings.as_dict(), "constants": bk.reference_constants(),
                "rows": len(self.rows),
                "rows_sha256": provenance.canonical_sha256(list(self.rows)),
                "enhancement_columns": self.nbasis_enhancement,
                "enhancement_scale": self.enhancement_scale,
                "enhancement_basis_sha256": self.enhancement_basis_sha256,
                "records_per_block": self.records_per_block,
                "cuda_source": cuda_sources.source_identity("pressure"),
                "dataset": self.dataset.identity() if not self.dataset.released else None,
                "memory_plan": self.plan.as_dict() if hasattr(self, "plan") else None,
                "timing_s": dict(self.timing)}

    def summary(self):
        ds = self.dataset
        return "%s (%s): %d points, %d records, evaluator plan %s" % (
            self.name, self.mode, ds.npoints, ds.nrecords, human(self.plan.total_bytes))

    def release(self):
        if self.released:
            return
        self.released = True
        self.arrays.clear()
        if self.reservation is not None:
            self.reservation.release()
            self.reservation = None
        if getattr(self, "dataset", None) is not None:
            self.dataset.release()
        if self.owns_budget and self.budget is not None:
            self.budget.close()

    close = release

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.release()
        return False


def from_problem(problem, run, **kwargs):
    """A CudaBackend bound exactly like ``backend.ReferenceBackend.from_problem``."""
    return CudaBackend(problem.M, run, problem.enhancement_basis if problem.nb else None,
                       problem.enhancement_scale, **kwargs)
