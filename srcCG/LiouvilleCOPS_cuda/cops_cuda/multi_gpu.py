"""One fit on several GPUs: row-partitioned workers behind the C02 contract (ROADMAP C07).

Design (ROADMAP C07 item 1)
---------------------------
* **One optimizer.** The SciPy trajectory stays in the host process. A
  ``MultiGpuBackend`` is an ordinary ``backend.PressureBackend``, so the facade,
  ``gpu_model.GpuModel`` and the unchanged ``fit_cops.Problem`` use it like the
  single-GPU ``CudaBackend``.
* **Persistent spawned workers**, one per device. Each owns a ``CudaBackend``
  over a contiguous block of the model's rows, chosen so that the record counts
  (which dominate the work) are balanced. A worker rebuilds its model with the
  reference parser (``host_model.header_model``) from the same table rows,
  header and coefficients, and checks that its element cache is the parent's
  (digest of ``Y``). It loads only its block's kernels to its device.
* **One activity field.** Per evaluation the parent computes
  ``act_lin = Y a (+ A_scale B b)`` once, with ``CudaBackend.activity_linear``'s
  expression, and writes it to shared memory. The workers clip it
  themselves, exactly as one process would. Every device therefore sees the
  same bytes. The protocol is synchronous: the parent writes the next field
  only after every worker has replied.
* **Assembly in row order.** Worker results are scattered back into the
  original row positions.

Exactness
---------
A point's pressures and Jacobian rows depend only on its own records: the
tree blocks are fixed per point, so ``point_set_invariance_bitwise`` holds
(``c06_full_run.py``). The kernels are deterministic across identical devices
(the C04 cross-device gravity check was bitwise). A one-device and a
multi-device evaluation of the same state are therefore expected to be
**bitwise** identical. ``tools/c07_multi_gpu_check.py`` and
``tests/test_multi_gpu.py`` check exactly that.

Failure and cancellation
------------------------
* A worker exception is sent back with its traceback and raised in the parent
  as ``BackendError``.
* A dead worker (EOF) or a reply timeout is a ``BackendError``.
* After any failure the backend is unusable, and ``close()`` stops every worker:
  first a polite request, then ``terminate`` after a grace period.
* The backend is a context manager, and ``close`` is idempotent.
* Each worker owns its device's ``DeviceBudget``, so device memory is bounded
  per device and released when the worker exits.
"""

import multiprocessing as mp
from multiprocessing import shared_memory
import os
import time
import traceback

import numpy as np

from . import backend as bk
from . import provenance
from .errors import BackendError, IdentityError

DEFAULT_TIMEOUT_S = 600.0       # bind of a 788k half-dataset takes about 70 s


def partition(counts, ndev):
    """Contiguous row blocks ``[(start, stop), ...]`` with balanced record counts.

    Cut k is placed where the running record count first reaches
    ``k * total / ndev``, so every block is within one point of the ideal.
    Every block holds at least one row.
    """
    counts = np.asarray(counts, dtype=np.int64)
    n = len(counts)
    if ndev < 1 or n < ndev:
        raise ValueError("cannot split %d rows over %d devices" % (n, ndev))
    cumulative = np.cumsum(counts)
    total = int(cumulative[-1]) if n else 0
    cuts = [0]
    for k in range(1, ndev):
        cut = int(np.searchsorted(cumulative, k * total / ndev, side="left")) + 1
        cut = min(max(cut, cuts[-1] + 1), n - (ndev - k))
        cuts.append(cut)
    cuts.append(n)
    return [(cuts[k], cuts[k + 1]) for k in range(ndev)]


def _worker(conn, spec, sub_rows, device_id, mode, records_per_block, basis, scale,
            shm_name, nfaces, y_digest):
    """Worker main: bind one CudaBackend for ``sub_rows`` and serve evaluations."""
    shm = None
    backend = None
    try:
        from . import evaluator, host_model, reference
        reference.bind()
        fit_cops = reference.load("fit_cops")
        elements = fit_cops.read_elements(os.path.join(spec["run"], "elements.bin"))
        if provenance.array_sha256(elements["Y"]) != y_digest:
            raise IdentityError("worker element cache differs from the parent's")
        model = host_model.header_model(spec["run"], spec["header"], spec["cols"],
                                        spec["rows"], elements, spec["coef"])
        model.load_kernels(np.asarray(sub_rows, dtype=int))
        backend = evaluator.CudaBackend(model, spec["run"], basis, scale, device_id=device_id,
                                        mode=mode, records_per_block=records_per_block)
        shm = shared_memory.SharedMemory(name=shm_name)
        act = np.ndarray((nfaces,), dtype=np.float64, buffer=shm.buf)
        backend.activity_source = lambda state: act.copy()
        conn.send(("ready", {"device_id": device_id, "rows": len(sub_rows),
                             "records": int(backend.dataset.nrecords),
                             "bind_s": backend.timing.get("bind_s"),
                             "dataset_load_s": backend.timing.get("dataset_load_s"),
                             "pid": os.getpid()}))
    except BaseException:
        conn.send(("error", traceback.format_exc()))
        if backend is not None:
            backend.release()
        return
    try:
        while True:
            message = conn.recv()
            if message[0] == "close":
                break
            if message[0] == "identity":
                conn.send(("ok", backend.identity()))
                continue
            if message[0] == "seen":
                from . import masks
                conn.send(("ok", masks.seen_positive(backend.dataset)))
                continue
            try:
                _, state, request = message
                t0 = time.perf_counter()
                result = backend.evaluate(state, request)
                conn.send(("ok", result.pNG, result.pRG, dict(result.columns),
                           time.perf_counter() - t0))
            except Exception:
                conn.send(("error", traceback.format_exc()))
    except (EOFError, KeyboardInterrupt):
        pass
    finally:
        backend.release()
        if shm is not None:
            del act
            shm.close()
        try:
            conn.send(("closed",))
        except (BrokenPipeError, OSError):
            pass


class MultiGpuBackend(bk.PressureBackend):
    """The C06 evaluator on several devices, behind the C02 contract."""

    name = "cuda_fp64_v1_multi"

    def __init__(self, spec, model, enhancement_basis=None, enhancement_scale=1.0,
                 devices=(0, 1), mode="tree", records_per_block=4096,
                 timeout_s=DEFAULT_TIMEOUT_S):
        nb = 0 if enhancement_basis is None else int(np.asarray(enhancement_basis).shape[1])
        super().__init__(model.rows, bk.ModelSettings.from_model(model), nb)
        self.devices = [int(d) for d in devices]
        if len(set(self.devices)) != len(self.devices):
            raise ValueError("devices must be distinct")
        self.mode, self.timeout_s = mode, float(timeout_s)
        self.records_per_block = int(records_per_block)
        self.Y = model.Y
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale)
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else bk.array_digest(enhancement_basis))
        counts = np.diff(model.indptr)
        self.blocks = partition(counts, len(self.devices))
        self.block_records = [int(counts[a:b].sum()) for a, b in self.blocks]
        self.workers, self.closed, self.broken = [], False, False
        self.timing = {"worker_s": [], "evaluate_s": []}
        nf = len(self.Y)
        self.shm = shared_memory.SharedMemory(create=True, size=nf * 8)
        self.act = np.ndarray((nf,), dtype=np.float64, buffer=self.shm.buf)
        ctx = mp.get_context("spawn")
        rows = np.asarray(model.rows, dtype=int)
        y_digest = provenance.array_sha256(self.Y)
        t0 = time.perf_counter()
        try:
            for device, (a, b) in zip(self.devices, self.blocks):
                parent, child = ctx.Pipe()
                proc = ctx.Process(target=_worker, name="cops-gpu%d" % device,
                                   args=(child, spec, rows[a:b].tolist(), device, mode,
                                         self.records_per_block, enhancement_basis,
                                         self.enhancement_scale, self.shm.name, nf, y_digest),
                                   daemon=True)
                proc.start()
                child.close()
                self.workers.append({"device": device, "proc": proc, "conn": parent,
                                     "positions": np.arange(a, b)})
            self.ready = [self._receive(w, "bind") for w in self.workers]
        except BaseException:
            self.close()
            raise
        self.timing["bind_s"] = time.perf_counter() - t0

    @classmethod
    def from_problem(cls, problem, spec, **kwargs):
        return cls(spec, problem.M, problem.enhancement_basis if problem.nb else None,
                   problem.enhancement_scale, **kwargs)

    # ------------------------------------------------------------ messaging
    def _send(self, worker, message, what):
        try:
            worker["conn"].send(message)
        except (BrokenPipeError, OSError) as exc:
            self.broken = True
            raise BackendError("GPU worker on device %d is unreachable for %s (%s; exit code %s)" %
                               (worker["device"], what, exc.__class__.__name__,
                                worker["proc"].exitcode)) from None

    def _receive(self, worker, what):
        conn, proc = worker["conn"], worker["proc"]
        deadline = time.monotonic() + self.timeout_s
        while not conn.poll(1.0):
            if not proc.is_alive():
                self.broken = True
                raise BackendError("GPU worker on device %d died during %s (exit code %s)" %
                                   (worker["device"], what, proc.exitcode))
            if time.monotonic() > deadline:
                self.broken = True
                raise BackendError("GPU worker on device %d: no reply to %s within %.0f s" %
                                   (worker["device"], what, self.timeout_s))
        try:
            reply = conn.recv()
        except EOFError:
            self.broken = True
            raise BackendError("GPU worker on device %d closed its pipe during %s" %
                               (worker["device"], what)) from None
        if reply[0] == "error":
            self.broken = True
            raise BackendError("GPU worker on device %d failed during %s:\n%s" %
                               (worker["device"], what, reply[1]))
        return reply[1] if reply[0] == "ready" else reply

    # ------------------------------------------------------------ contract
    def capabilities(self):
        return {"sza_modes": {"continuous"},
                "heliocentric_models": {self.settings.heliocentric_model},
                "rg_models": {"equilibrium_sphere", "calibration_srccg", "dsmc_srccg"},
                "column_kinds": {"a", "b", "T", "kappa"}, "enhancement": True}

    def activity_linear(self, state):
        """The expression of ``CudaBackend.activity_linear`` (Model.evaluate)."""
        act_lin = self.Y @ state.a
        if state.b is not None:
            act_lin = act_lin + self.enhancement_scale * (self.enhancement_basis @ state.b)
        return act_lin

    def _evaluate(self, state, request):
        if self.closed or self.broken:
            raise BackendError("multi-GPU backend is closed or broken")
        t0 = time.perf_counter()
        self.act[:] = self.activity_linear(state)
        for w in self.workers:
            self._send(w, ("evaluate", state, request), "evaluate")
        n = len(self.rows)
        pNG, pRG = np.empty(n), np.empty(n)
        columns = {key: (np.empty(n), np.empty(n)) for key in request.keys()}
        worker_s = []
        for w in self.workers:
            _, ng, rg, cols, seconds = self._receive(w, "evaluate")
            pos = w["positions"]
            pNG[pos], pRG[pos] = ng, rg
            for key in request.keys():
                columns[key][0][pos] = cols[key][0]
                columns[key][1][pos] = cols[key][1]
            worker_s.append(seconds)
        self.timing["worker_s"] = worker_s
        self.timing["last_evaluate_s"] = time.perf_counter() - t0
        return bk.EvaluationResult(pNG, pRG, columns, request)

    def seen_faces(self):
        """``Model.seen_faces`` over all workers: the union of their device masks."""
        seen = None
        for w in self.workers:
            self._send(w, ("seen",), "seen_faces")
        for w in self.workers:
            mask = self._receive(w, "seen_faces")[1]
            seen = mask if seen is None else (seen | mask)
        return seen

    def identity(self):
        workers = []
        if not (self.closed or self.broken):
            for w in self.workers:
                self._send(w, ("identity",), "identity")
                workers.append(self._receive(w, "identity")[1])
        return {"backend": self.name, "contract_version": bk.CONTRACT_VERSION,
                "mode": self.mode, "devices": self.devices,
                "partition": [{"device": d, "rows": [a, b], "records": r}
                              for d, (a, b), r in zip(self.devices, self.blocks,
                                                      self.block_records)],
                "rows": len(self.rows),
                "rows_sha256": provenance.canonical_sha256(list(self.rows)),
                "enhancement_columns": self.nbasis_enhancement,
                "enhancement_scale": self.enhancement_scale,
                "enhancement_basis_sha256": self.enhancement_basis_sha256,
                "records_per_block": self.records_per_block,
                "workers_ready": getattr(self, "ready", None),
                "workers": workers, "timing_s": dict(self.timing)}

    def summary(self):
        return "%s (%s) on devices %s: %d points, records %s" % (
            self.name, self.mode, self.devices, len(self.rows), self.block_records)

    # ------------------------------------------------------------ lifecycle
    @property
    def released(self):
        return self.closed

    def close(self, grace_s=30.0):
        if self.closed:
            return
        self.closed = True
        for w in self.workers:
            try:
                w["conn"].send(("close",))
            except (BrokenPipeError, OSError):
                pass
        deadline = time.monotonic() + grace_s
        for w in self.workers:
            w["proc"].join(max(0.0, deadline - time.monotonic()))
            if w["proc"].is_alive():
                w["proc"].terminate()
                w["proc"].join(5.0)
            w["conn"].close()
        if getattr(self, "shm", None) is not None:
            del self.act
            self.shm.close()
            self.shm.unlink()
            self.shm = None

    release = close

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
