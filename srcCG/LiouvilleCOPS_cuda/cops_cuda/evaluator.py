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
  ``activity="device"`` (SCALING_PLAN 3.5, user decision Q3; the default of the
  MPI transport) forms them on the device instead, with a fixed summation order
  (``kernels/activity.cu``): deterministic and partition-independent, within a
  few ulp of the terms of the host expression, checked against it at the first
  evaluation and every ``ACTIVITY_CHECK_EVERY``-th.
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

import os
import time

import numpy as np

from . import backend as bk
from . import cuda_sources, precision, provenance, reference
from .device_store import DeviceDataset, RECORD_FIELDS_FP32
from .errors import BackendError, IdentityError, UnsupportedConfiguration
from .memory import DeviceBudget, MemoryPlan, human

COMPACT_ENVIRONMENT = "COPS_CUDA_COMPACT_FACES"
# Streaming (SCALING_PLAN P2): records beyond the device budget stay in pinned host memory and
# are copied through STREAM_SLOTS device slots each pass, overlapped with the kernels.
RESIDENT_ENVIRONMENT = "COPS_CUDA_MAX_RESIDENT_RECORDS"   # cap (tests, or to force streaming)
STREAM_CHUNK_ENVIRONMENT = "COPS_CUDA_STREAM_CHUNK_RECORDS"
STREAM_SLOTS = 3
STREAM_CHUNK_RECORDS = 1 << 25
BASIS_COLUMNS_RESERVE = 8     # template columns a later rebind may add (planned when streaming)
ACTIVITY_ENVIRONMENT = "COPS_CUDA_ACTIVITY"
ACTIVITY_MODES = ("host", "device")
ACTIVITY_CHECK_EVERY = 100    # device activity: compared with the host expression every N-th evaluation
ACTIVITY_CHECK_EPS = 64       # |device - host| <= 64 eps * sum of |terms| (summation bound ~27 eps)
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
                 chunk_records=None, fused=True, record_precision=None, compact_faces=None,
                 record_store=None, activity=None, scratch_per_face=0):
        """``record_precision`` "fp32" (mixed precision, experiment 2026-10-21): kernel records
        stored in FP32, factors recomputed per record and all arithmetic in FP64
        (``pressure.cu: eval_fused_m_*``); default ``$COPS_CUDA_RECORD_PRECISION`` or fp64.
        ``compact_faces`` (default on; ``COPS_CUDA_COMPACT_FACES=0`` turns it off): per-face
        arrays only for the faces that occur in the bound records (SCALING_PLAN 3.5); results
        are bitwise those without it. ``record_store`` (default ``$COPS_CUDA_RECORD_STORE``,
        FP32 records only): read the records from an FP32 record store of this base.
        ``activity`` "host" (default, exact reference expression) or "device" (default
        ``$COPS_CUDA_ACTIVITY``): where act_lin and the clip are formed (module docstring).
        ``scratch_per_face``: bytes per device facet that callers allocate per call beyond the
        plan (restart and fine-tuning kernels), kept free when the resident records are sized."""
        if mode not in MODES:
            raise ValueError("mode must be one of %s" % (MODES,))
        self.activity_mode = activity or os.environ.get(ACTIVITY_ENVIRONMENT) or "host"
        if self.activity_mode not in ACTIVITY_MODES:
            raise ValueError("activity must be one of %s" % (ACTIVITY_MODES,))
        self.record_precision = (precision.record_precision_default() if record_precision is None
                                 else record_precision)
        if self.record_precision not in precision.RECORD_PRECISIONS:
            raise ValueError("record precision must be one of %s" % (precision.RECORD_PRECISIONS,))
        self.mixed = self.record_precision == "fp32"
        if self.mixed and (mode != "tree" or not fused):
            raise UnsupportedConfiguration("FP32 record storage supports the fused tree mode only")
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
        self.scratch_per_face = int(scratch_per_face)
        self.activity_checks = {"checks": 0, "faces_checked": 0, "max_error_eps": 0.0,
                                "max_ulp": 0, "active_flips": 0, "evaluations": 0}
        self.arrays = {}
        self.reservation = None
        self.released = False
        self._columns_cache = None
        point_ids = [int(p) for p in np.asarray(model.ipoint)[model.rows]]
        self.owns_budget = budget is None
        self.budget = DeviceBudget(self.device_id) if budget is None else budget
        t0 = time.perf_counter()
        kwargs = {} if chunk_records is None else {"chunk_records": chunk_records}
        if compact_faces is None:
            compact_faces = os.environ.get(COMPACT_ENVIRONMENT, "1") != "0"
        if record_store is None and self.mixed:
            from .record_store import environment_store
            record_store = environment_store()
        if record_store is not None:
            kwargs["record_store"] = record_store
        resident = self._resident_records(model) if self.mixed else None
        if resident is not None:
            kwargs["resident_records"] = resident
        self.dataset = DeviceDataset.load(run, point_ids=point_ids, device_id=self.device_id,
                                          mesh_fields=("ACO2",), budget=self.budget,
                                          record_precision=self.record_precision,
                                          compact_faces=bool(compact_faces), **kwargs)
        try:
            self._check_dataset(model)
            self.timing["dataset_load_s"] = time.perf_counter() - t0
            self._bind()
        except BaseException:
            self.release()
            raise

    def rebind_basis(self, enhancement_basis, enhancement_scale=1.0):
        """Change the template basis, keeping the device records (SCALING_PLAN P1).

        Only the evaluator's own arrays (static factors, PhiCO2, activity and column scratch,
        sized by the basis) are released and bound again by the same ``_bind`` as a new
        backend, so the results are bitwise those of a backend created with this basis; the
        dataset (the records, the expensive part) stays on the device.
        """
        if self.released:
            raise BackendError("backend is released")
        self._slots = None
        self.arrays.clear()
        if self.reservation is not None:
            self.reservation.release()
            self.reservation = None
        self.cupy.get_default_memory_pool().free_all_blocks()
        self.nbasis_enhancement = 0 if enhancement_basis is None else int(np.asarray(enhancement_basis).shape[1])
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale) if enhancement_basis is not None else 1.0
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else bk.array_digest(enhancement_basis))
        self._columns_cache = None
        t0 = time.perf_counter()
        self._bind()
        self.timing["rebind_basis_s"] = time.perf_counter() - t0

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
        if self.mixed:
            return self._bind_mixed()
        cupy, ds = self.cupy, self.dataset
        npt, nrec, nf = ds.npoints, ds.nrecords, ds.nfaces_device
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
        self._plan_activity(nf)
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
            self._alloc_activity(nf)
            # The parameter-independent CO2 flux sums (Model.prepare PhiCO2).
            a["phiC"] = phiC
            self.PhiCO2 = self._run_sums([PHI_CO2_CODE], None)[0]
            del a["phiC"], phiC, per_point
            self.cupy.get_default_memory_pool().free_all_blocks()
        self.timing["bind_s"] = time.perf_counter() - t0

    def _bind_mixed(self):
        """FP32 records: no stored per-record factors; per-point constants and PhiCO2 only."""
        cupy, ds = self.cupy, self.dataset
        npt, nf = ds.npoints, ds.nfaces_device
        blocks = self._block_table(ds.scan.counts, ds.scan.offsets)
        self.plan = MemoryPlan("cuda evaluator (fp32 records)")
        self.plan.add_array("eval.per_point_constants", "static", (5, npt), np.float64)
        self.plan.add_array("eval.af_face", "scratch", (nf,), np.float64)
        self.plan.add_array("eval.active", "scratch", (nf,), np.uint8)
        self.plan.add_array("eval.colbasis", "basis", (max(1, 25 + self.nbasis_enhancement), nf),
                            np.float64, "requested activity/template columns")
        self._plan_activity(nf)
        self.plan.add_array("eval.blocks", "static", (len(blocks),), np.int64)
        self.plan.add_array("eval.first_block", "static", (npt + 1,), np.int64)
        # Streaming runs every tile of a request on a chunk while it is resident, so the
        # partial sums of all tiles are kept at once (the resident path reuses one tile's).
        nsums = (14 + 2 * (25 + self.nbasis_enhancement)) if ds.streamed else MAX_SUMS
        self.plan.add_array("eval.partial", "partial_sums", (len(blocks) // 3, nsums), np.float64)
        self.plan.add_array("eval.sums", "scratch", (nsums, npt), np.float64)
        rebased = self._bind_streaming(blocks) if ds.streamed else None
        self.reservation = self.budget.reserve(self.plan)
        t0 = time.perf_counter()
        with cupy.cuda.Device(self.device_id), self.budget.allocation_guard(self.plan):
            a = self.arrays
            c = self.constants
            a["pt"] = cupy.asarray(np.ascontiguousarray(np.stack(
                (c.qmin["H2O"], c.qmax["H2O"], c.qmin["CO2"], c.qmax["CO2"], c.vdotb)).astype(np.float64)))
            a["blocks"] = cupy.asarray(blocks)
            a["first_block"] = cupy.asarray(self._first_block)
            a["partial"] = cupy.empty((len(blocks) // 3) * nsums, dtype=np.float64)
            a["sums"] = cupy.empty(nsums * npt, dtype=np.float64)
            if rebased is not None:
                self._alloc_streaming(rebased)
            a["af_face"] = cupy.zeros(nf, dtype=np.float64)
            a["active"] = cupy.zeros(nf, dtype=np.uint8)
            a["colbasis"] = cupy.empty(max(1, 25 + self.nbasis_enhancement) * nf, dtype=np.float64)
            self._alloc_activity(nf)
            self.PhiCO2 = self._run_fused([("eval_fused_m_phic", 0, [PHI_CO2_CODE])], [PHI_CO2_CODE],
                                          {k: 0.0 for k in ("Tmax", "T75", "Tmin", "KT", "cNbmH",
                                                            "kcNbmC", "bmH", "bmC")})[0]
            self.cupy.get_default_memory_pool().free_all_blocks()
        self.timing["bind_s"] = time.perf_counter() - t0

    # ----------------------------------------------------------- device activity
    def _plan_activity(self, nf):
        if self.activity_mode == "device":
            self.plan.add_array("eval.actbasis", "basis", (self.Y.shape[1] + self.nbasis_enhancement, nf),
                                np.float64, "Y and unscaled template columns of the device facets")

    def _alloc_activity(self, nf):
        """Upload Y and the raw template columns of the device's facets, term-major, once per bind."""
        if self.activity_mode != "device":
            return
        cupy, rows = self.cupy, self.dataset.local_faces
        ny, nb = self.Y.shape[1], self.nbasis_enhancement
        target = self.arrays["actbasis"] = cupy.empty((ny + nb) * nf, dtype=np.float64)
        step = 1 << 20
        for f0 in range(0, nf, step):
            f1 = min(nf, f0 + step)
            sel = slice(f0, f1) if rows is None else rows[f0:f1]
            parts = [np.asarray(self.Y[sel], dtype=np.float64)]
            if nb:
                parts.append(np.asarray(np.asarray(self.enhancement_basis)[sel], dtype=np.float64))
            block = np.ascontiguousarray(np.column_stack(parts).T)        # (ny + nb, f1 - f0)
            for j in range(ny + nb):
                target[j * nf + f0:j * nf + f1].set(block[j])

    def _device_activity(self, state, act_out=None):
        """af and active (and act_lin into ``act_out``) on the device, kernels/activity.cu."""
        cupy, a = self.cupy, self.arrays
        nf, ny = self.dataset.nfaces_device, self.Y.shape[1]
        nb = self.nbasis_enhancement if state.b is not None else 0
        coef = np.concatenate([np.asarray(state.a, dtype=np.float64).reshape(-1),
                               np.asarray(state.b if nb else [], dtype=np.float64).reshape(-1)])
        if len(coef) != ny + nb:
            raise BackendError("activity coefficients: %d values for %d columns" % (len(coef), ny + nb))
        d_coef = cupy.asarray(coef)
        cuda_sources.function("activity", "activity_linear_clip", self.device_id)(
            (max(1, min(65535, -(-nf // 256))),), (256,),
            (a["actbasis"], np.int64(nf), np.int32(ny), np.int32(nb), d_coef,
             np.float64(self.enhancement_scale), a["af_face"], a["active"],
             a["af_face"] if act_out is None else act_out, np.int32(act_out is not None)))

    def activity_check(self, state):
        """Device act_lin against the host expression on every device facet (invariant I7).

        Error measure: |device - host| / (sum_j |Y_j a_j| + A_scale sum_k |B_k b_k|) in units of
        eps (a summation of n terms differs by at most ~n eps of that sum). Raises BackendError
        beyond ``ACTIVITY_CHECK_EPS`` or for an active-set flip outside that error bound.
        """
        from .memory import scratch_reservation
        cupy, ds = self.cupy, self.dataset
        nf = ds.nfaces_device
        with cupy.cuda.Device(self.device_id), scratch_reservation(self.budget, "activity check",
                                                                   {"act": (nf, np.float64)}):
            act = cupy.empty(nf, dtype=np.float64)
            self._device_activity(state, act)
            device = act.get()
            del act
        rows = ds.local_faces
        a_abs = np.abs(np.asarray(state.a, dtype=np.float64))
        with_b = state.b is not None and self.nbasis_enhancement
        eps = np.finfo(np.float64).eps
        ordered = lambda x: np.where(x.view(np.int64) < 0, np.int64(-2**63) - x.view(np.int64), x.view(np.int64))
        result = {"faces": int(nf), "max_error_eps": 0.0, "max_ulp": 0, "active_flips": 0}
        bad = 0
        for f0 in range(0, nf, 1 << 20):                     # host side in chunks of facets
            f1 = min(nf, f0 + (1 << 20))
            sel = slice(f0, f1) if rows is None else rows[f0:f1]
            Y = self.Y[sel]
            host, magnitude = Y @ state.a, np.abs(Y) @ a_abs
            if with_b:
                B = np.asarray(self.enhancement_basis)[sel]
                host = host + self.enhancement_scale * (B @ state.b)
                magnitude = magnitude + abs(self.enhancement_scale) * (np.abs(B) @ np.abs(state.b))
            dev = device[f0:f1]
            diff = np.abs(dev - host)
            with np.errstate(divide="ignore", invalid="ignore"):
                err = np.where(diff == 0.0, 0.0, diff / (magnitude * eps))
            flips = (dev > 0.0) != (host > 0.0)
            result["max_error_eps"] = max(result["max_error_eps"], float(err.max()))
            result["max_ulp"] = max(result["max_ulp"], int(np.abs(ordered(dev) - ordered(host)).max()))
            result["active_flips"] += int(flips.sum())
            bad += int((flips & (np.abs(host) > ACTIVITY_CHECK_EPS * eps * magnitude)).sum())
            bad += int((np.isfinite(dev) != np.isfinite(host)).sum())
        c = self.activity_checks
        c["checks"] += 1
        c["faces_checked"] += result["faces"]
        c["max_error_eps"] = max(c["max_error_eps"], result["max_error_eps"])
        c["max_ulp"] = max(c["max_ulp"], result["max_ulp"])
        c["active_flips"] += result["active_flips"]
        if result["max_error_eps"] > ACTIVITY_CHECK_EPS or bad:
            raise BackendError("device activity differs from the host expression: %s (bound %d eps; "
                               "%d flips or non-finite values outside it)" % (result, ACTIVITY_CHECK_EPS, bad))
        return result

    # ---------------------------------------------------------------- streaming
    def _block_ends(self, counts):
        """Global record index at the end of every block, in block order."""
        ends, start = [], 0
        for count in counts:
            count = int(count)
            nblocks = max(1, -(-count // self.records_per_block))
            ends.extend(start + min((k + 1) * self.records_per_block, count) for k in range(nblocks))
            start += count
        return np.asarray(ends, dtype=np.int64)

    def _stream_chunk_records(self):
        return max(self.records_per_block,
                   int(os.environ.get(STREAM_CHUNK_ENVIRONMENT, STREAM_CHUNK_RECORDS)))

    def _resident_records(self, model):
        """Records kept on the device; None when all fit (no streaming).

        From the budget left after the evaluator's per-face, block and partial-sum arrays
        (upper bounds: all mesh faces, ``BASIS_COLUMNS_RESERVE`` more template columns) and,
        when streaming, the slots. Aligned down to a block end, so blocks never straddle the
        resident/streamed boundary. ``$COPS_CUDA_MAX_RESIDENT_RECORDS`` caps it (tests).
        """
        counts = np.diff(np.asarray(model.indptr, dtype=np.int64))
        nrec, npt = int(counts.sum()), len(counts)
        ends = self._block_ends(counts)
        nblocks, nf = len(ends), len(self.Y)
        ncols = 25 + max(self.nbasis_enhancement, BASIS_COLUMNS_RESERVE)
        nsums = 14 + 2 * ncols
        common = (nf * (8 + 1 + 8 + 8 * ncols + 4) + 4 * (1 << 26) + 2 * 24 * nblocks
                  + 8 * (npt + 1) + 8 * nsums * npt + 128 * npt + (64 << 20))
        if self.activity_mode == "device":                 # eval.actbasis (Y and raw templates)
            common += 8 * nf * (self.Y.shape[1] + max(self.nbasis_enhancement, BASIS_COLUMNS_RESERVE))
        # per-call scratch beyond the plan: seen_faces (1 B), the activity check (8 B, device
        # activity) and the restart/fine-tuning kernels' arrays (``scratch_per_face``)
        common += nf * (1 + (8 if self.activity_mode == "device" else 0) + self.scratch_per_face)
        avail = int(self.budget.available_bytes() * 0.99)
        cap = os.environ.get(RESIDENT_ENVIRONMENT)
        cap = None if cap in (None, "") else int(cap)
        if cap is None and nrec * 20 + common + 8 * MAX_SUMS * nblocks <= avail:
            return None
        streaming = STREAM_SLOTS * self._stream_chunk_records() * 20 + 8 * nsums * nblocks
        capacity = max(0, (avail - common - streaming) // 20)
        if cap is not None:
            capacity = min(capacity, cap)
        if capacity >= nrec:
            return None
        k = int(np.searchsorted(ends, capacity, side="right"))
        return int(ends[k - 1]) if k else 0

    def _bind_streaming(self, blocks):
        """Chunks of the streamed blocks, their rebased block tables, the slots (plan entries)."""
        ds = self.dataset
        table = blocks.reshape(-1, 3)
        nres, chunk = ds.nresident, self._stream_chunk_records()
        g_res = int(np.searchsorted(table[:, 2], nres, side="right"))
        while g_res < len(table) and table[g_res, 1] < nres:        # never split a block
            g_res -= 1
        chunks, g = [], g_res
        while g < len(table):
            a, g1 = int(table[g, 1]), g + 1
            while g1 < len(table) and int(table[g1, 2]) - a <= chunk:
                g1 += 1
            chunks.append((g, g1, a, int(table[g1 - 1, 2])))
            g = g1
        rebased = table[g_res:].copy()
        for g0, g1, a, _ in chunks:
            rebased[g0 - g_res:g1 - g_res, 1:] -= a
        self._stream = {"g_res": g_res, "chunks": chunks,
                        "slot_records": max([b - a for _, _, a, b in chunks] + [1])}
        self.plan.add_array("eval.stream_blocks", "static", (rebased.size,), np.int64,
                            "streamed blocks rebased to their chunk")
        for k in range(STREAM_SLOTS):
            for name, _ in RECORD_FIELDS_FP32:
                self.plan.add_array("eval.slot%d.%s" % (k, name), "scratch",
                                    (self._stream["slot_records"],), ds.dtypes[name])
        return rebased.reshape(-1)

    def _alloc_streaming(self, rebased):
        cupy, ds = self.cupy, self.dataset
        a = self.arrays
        a["stream_blocks"] = cupy.asarray(rebased)
        n = self._stream["slot_records"]
        for k in range(STREAM_SLOTS):                      # in self.arrays: freed with them
            for name, _ in RECORD_FIELDS_FP32:
                a["slot%d.%s" % (k, name)] = cupy.empty(n, dtype=ds.dtypes[name])
        self._slots = [{name: a["slot%d.%s" % (k, name)] for name, _ in RECORD_FIELDS_FP32}
                       for k in range(STREAM_SLOTS)]
        self._streams = (cupy.cuda.Stream(non_blocking=True), cupy.cuda.Stream(non_blocking=True))

    def _run_fused_streamed(self, plan, codes, params):
        """``_run_fused`` with streamed records: the same kernels, block for block, so the
        per-block partials and their block-order reduction are bitwise those of the resident
        path. Resident blocks run first; each streamed chunk is copied into a slot on the copy
        stream while earlier chunks compute, and every tile of the plan runs on a chunk while
        it is on the device (chunk-major: one transfer per pass)."""
        cupy, a, ds = self.cupy, self.arrays, self.dataset
        npt, nf = ds.npoints, ds.nfaces_device
        order = ("Tmax", "T75", "Tmin", "KT", "cNbmH", "kcNbmC", "bmH", "bmC")
        P = cupy.asarray(np.array([params[k] for k in order], dtype=np.float64))
        nblocks = int(self._first_block[-1])
        st, d, host = self._stream, ds.arrays, ds.host_records
        g_res, nres = st["g_res"], ds.nresident
        tiles, offset = [], 0
        for name, col0, launch_codes in plan:
            kernel = name if name.startswith("eval_fused_m_") else name.replace("eval_fused_", "eval_fused_m_")
            n = len(launch_codes)
            tiles.append((cuda_sources.function("pressure", kernel, self.device_id), col0, n, offset, launch_codes))
            offset += nblocks * n
        if offset > a["partial"].size or sum(t[2] for t in tiles) * npt > a["sums"].size:
            raise BackendError("internal: streamed request exceeds the planned partial sums")
        compute, copy = self._streams
        common = lambda: (a["pt"], np.int64(npt), d["mesh.ACO2"], np.float64(self.constants.COS75),
                          a["af_face"], a["active"], a["colbasis"], np.int64(nf))
        cupy.cuda.Device(self.device_id).synchronize()      # activity/columns uploads are complete
        if g_res:
            for fn, col0, n, off, _ in tiles:
                fn((g_res,), (BLOCK,),
                   (a["blocks"], d["kernel.face"], d["kernel.kn"], d["kernel.kphi"], d["kernel.kv"],
                    d["kernel.cos"]) + common() + (np.int32(col0), P, a["partial"][off:off + g_res * n]),
                   stream=compute)
        free = [None] * STREAM_SLOTS
        for i, (g0, g1, r0, r1) in enumerate(st["chunks"]):
            slot = self._slots[i % STREAM_SLOTS]
            if free[i % STREAM_SLOTS] is not None:
                copy.wait_event(free[i % STREAM_SLOTS])
            if r1 > r0:
                for name, _ in RECORD_FIELDS_FP32:
                    slot[name][:r1 - r0].set(host[name][r0 - nres:r1 - nres], stream=copy)
            compute.wait_event(copy.record())
            sb = a["stream_blocks"][3 * (g0 - g_res):3 * (g1 - g_res)]
            for fn, col0, n, off, _ in tiles:
                fn((g1 - g0,), (BLOCK,),
                   (sb, slot["face"], slot["kn"], slot["kphi"], slot["kv"], slot["cos"]) + common()
                   + (np.int32(col0), P, a["partial"][off + g0 * n:off + g1 * n]), stream=compute)
            free[i % STREAM_SLOTS] = compute.record()
        row = {code: i for i, code in enumerate(codes)}
        out = np.empty((len(codes), npt), dtype=np.float64)
        finish = cuda_sources.function("pressure", "eval_finish", self.device_id)
        base = 0
        for fn, col0, n, off, launch_codes in tiles:
            sums = a["sums"][base * npt:(base + n) * npt]
            finish((-(-(n * npt) // 256),), (256,),
                   (a["first_block"], np.int64(npt), a["partial"][off:off + nblocks * n], np.int32(n), sums),
                   stream=compute)
            base += n
        compute.synchronize()
        host_sums = a["sums"][:base * npt].get().reshape(base, npt)
        base = 0
        for fn, col0, n, off, launch_codes in tiles:
            for i, code in enumerate(launch_codes):
                out[row[code]] = host_sums[base + i]
            base += n
        self.timing["streamed_chunks"] = len(st["chunks"])
        return out

    def record_chunks(self):
        """The bound records in block order, for kernels over the block table (restart, fine
        tuning): yields ``(g0, g1, table, records)`` where ``table`` is the device block table
        whose first block is global block g0 and ``records`` maps face/kn/kphi/kv/cos to device
        arrays indexed as that table indexes them. Resident records come first; each streamed
        chunk is then copied into a slot (on the current stream, before the caller's launches),
        so a caller that indexes its per-block outputs by g0 + block sees every block once, in
        order. Without streaming: one chunk with the dataset's own arrays.
        """
        ds, a, d = self.dataset, self.arrays, self.dataset.arrays
        nblocks = int(self._first_block[-1])
        names = [n for n, _ in RECORD_FIELDS_FP32] if self.mixed else ["face", "kn", "kphi", "kv", "cos"]
        resident = {n: d["kernel." + n] for n in names}
        if not ds.streamed:
            yield 0, nblocks, a["blocks"], resident
            return
        st, host, nres = self._stream, ds.host_records, ds.nresident
        g_res = st["g_res"]
        if g_res:
            yield 0, g_res, a["blocks"][:3 * g_res], resident
        slot = self._slots[0]
        for g0, g1, r0, r1 in st["chunks"]:
            if r1 > r0:
                for n in names:
                    slot[n][:r1 - r0].set(host[n][r0 - nres:r1 - nres])
            yield g0, g1, a["stream_blocks"][3 * (g0 - g_res):3 * (g1 - g_res)], slot
        self.cupy.cuda.get_current_stream().synchronize()     # the slot is free for the evaluator

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
        if self.activity_mode == "device":
            n = self.activity_checks["evaluations"]
            self.activity_checks["evaluations"] = n + 1
            if n % ACTIVITY_CHECK_EVERY == 0:
                self.activity_check(state)
            self._device_activity(state)
            return
        act_lin = (self.activity_linear(state) if self.activity_source is None
                   else self.activity_source(state))
        act_lin = self.dataset.gather_faces(act_lin)          # the device's faces (compaction)
        self.arrays["af_face"].set(np.maximum(act_lin, 0.0))
        self.arrays["active"].set((act_lin > 0.0).astype(np.uint8))

    def _upload_columns(self, request):
        """Unmasked column basis: Y[:, j] and A_scale * B[:, k] (Model.evaluate)."""
        key = (request.basis_index, request.enhancement_index)
        if self._columns_cache == key:
            return
        parts = []
        rows = self.dataset.local_faces                       # None: all faces
        pick = (lambda m, cols: m[:, cols]) if rows is None else (lambda m, cols: m[np.ix_(rows, cols)])
        if request.basis_index:
            parts.append(pick(self.Y, list(request.basis_index)))
        if request.enhancement_index:
            parts.append(self.enhancement_scale *
                         pick(self.enhancement_basis, list(request.enhancement_index)))
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
        if self.dataset.streamed:
            return self._run_fused_streamed(plan, codes, params)
        cupy, a, ds = self.cupy, self.arrays, self.dataset
        npt, nf = ds.npoints, ds.nfaces_device
        order = ("Tmax", "T75", "Tmin", "KT", "cNbmH", "kcNbmC", "bmH", "bmC")
        P = cupy.asarray(np.array([params[k] for k in order], dtype=np.float64))
        nblocks = int(self._first_block[-1])
        row = {code: i for i, code in enumerate(codes)}
        out = np.empty((len(codes), npt), dtype=np.float64)
        d = ds.arrays
        for name, col0, launch_codes in plan:
            n = len(launch_codes)
            partial = a["partial"][:nblocks * n]
            if self.mixed:
                kernel = name if name.startswith("eval_fused_m_") else name.replace("eval_fused_", "eval_fused_m_")
                cuda_sources.function("pressure", kernel, self.device_id)(
                    (nblocks,), (BLOCK,),
                    (a["blocks"], d["kernel.face"], d["kernel.kn"], d["kernel.kphi"], d["kernel.kv"],
                     d["kernel.cos"], a["pt"], np.int64(npt), d["mesh.ACO2"],
                     np.float64(self.constants.COS75), a["af_face"], a["active"], a["colbasis"],
                     np.int64(nf), np.int32(col0), P, partial))
            else:
                cuda_sources.function("pressure", name, self.device_id)(
                    (nblocks,), (BLOCK,),
                    (a["blocks"], d["kernel.face"], a["wnH"], a["wvH"], a["wpH"], a["wnC"],
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
        if self.mixed:
            raise UnsupportedConfiguration("FP32 record storage evaluates the standard fused requests "
                                           "only; got codes %s" % (list(codes),))
        cupy, a, ds = self.cupy, self.arrays, self.dataset
        npt, nf = ds.npoints, ds.nfaces_device
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
                "mode": self.mode, "fused_kernels": self.fused,
                "precision": (precision.POLICY_NAME if not self.mixed else
                              {"policy": precision.MIXED_POLICY_NAME, "record_storage": "fp32",
                               "arithmetic": "float64",
                               "record_conversion": dict(self.dataset.conversion)}),
                "settings": self.settings.as_dict(), "constants": bk.reference_constants(),
                "rows": len(self.rows),
                "rows_sha256": provenance.canonical_sha256(list(self.rows)),
                "enhancement_columns": self.nbasis_enhancement,
                "enhancement_scale": self.enhancement_scale,
                "enhancement_basis_sha256": self.enhancement_basis_sha256,
                "records_per_block": self.records_per_block,
                "activity": {"mode": self.activity_mode,
                             "checks": dict(self.activity_checks) if self.activity_mode == "device" else None,
                             "cuda_source": (cuda_sources.source_identity("activity")
                                             if self.activity_mode == "device" else None)},
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
        self._slots = None
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
