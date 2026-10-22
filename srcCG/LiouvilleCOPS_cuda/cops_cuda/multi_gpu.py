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

Transports (SCALING_PLAN.md 3.3, phase P3)
------------------------------------------
``transport="local"`` (default): spawned worker processes on this node, as above.
``transport="mpi"``: one MPI rank per GPU, on any number of nodes (``cops_cuda/transport.py``).
Rank 0 drives the fit and evaluates the first block of rows on its own GPU in process; ranks
>= 1 serve the other blocks. ``activity="device"`` (the MPI default, user decision Q3) forms the
activity field on every GPU for its own facets instead of in the driver; ``"host"`` keeps the
exact driver-side expression (shared memory, or ``MPI_Bcast``).

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
import time

import numpy as np

from . import backend as bk
from . import provenance, transport as tp
from .errors import BackendError

DEFAULT_TIMEOUT_S = 600.0       # bind of a 788k half-dataset takes about 70 s
BIND_S_PER_RECORD = 1e-6        # bind allowance on top: SHAP5 3M loads ~2.3e-7 s/record per worker;
                                # SHAP8 12M (2.2e9 records per worker) exceeded the fixed 600 s


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


class MultiGpuBackend(bk.PressureBackend):
    """The C06 evaluator on several devices, behind the C02 contract."""

    name = "cuda_fp64_v1_multi"

    def __init__(self, spec, model, enhancement_basis=None, enhancement_scale=1.0,
                 devices=(0, 1), mode="tree", records_per_block=4096,
                 timeout_s=DEFAULT_TIMEOUT_S, transport="local", activity=None, comm=None):
        """``transport`` "local" (spawned workers on ``devices``) or "mpi" (one worker per rank
        of ``comm``, default COMM_WORLD; rank 0 is this process and uses its own device, the
        first of ``devices``). ``activity``: "host" (default local) or "device" (default MPI)."""
        nb = 0 if enhancement_basis is None else int(np.asarray(enhancement_basis).shape[1])
        super().__init__(model.rows, bk.ModelSettings.from_model(model), nb)
        if transport not in ("local", "mpi"):
            raise ValueError("transport must be 'local' or 'mpi'")
        self.transport = transport
        self.activity = activity or ("device" if transport == "mpi" else "host")
        if self.activity not in ("host", "device"):
            raise ValueError("activity must be 'host' or 'device'")
        if transport == "mpi":
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD if comm is None else comm
            if self.comm.Get_rank() != 0:
                raise RuntimeError("the MPI multi-GPU backend is driven by rank 0")
            nworkers = self.comm.Get_size()
            self.devices = [int(devices[0]) if devices else tp.mpi_local_device(self.comm)]
        else:
            self.comm = None
            self.devices = [int(d) for d in devices]
            if len(set(self.devices)) != len(self.devices):
                raise ValueError("devices must be distinct")
            nworkers = len(self.devices)
        self.mode, self.timeout_s = mode, float(timeout_s)
        self.records_per_block = int(records_per_block)
        self.Y = model.Y
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale)
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else bk.array_digest(enhancement_basis))
        counts = np.diff(model.indptr)
        self.blocks = partition(counts, nworkers)
        self.block_records = [int(counts[a:b].sum()) for a, b in self.blocks]
        self.workers, self.closed, self.broken = [], False, False
        self.timing = {"worker_s": [], "evaluate_s": []}
        nf = len(self.Y)
        self.shm, self.act = None, None
        if self.activity == "host":
            if transport == "local":
                self.shm = shared_memory.SharedMemory(create=True, size=nf * 8)
                self.act = np.ndarray((nf,), dtype=np.float64, buffer=self.shm.buf)
            else:
                self.act = np.empty(nf, dtype=np.float64)
        rows = np.asarray(model.rows, dtype=int)
        y_digest = provenance.array_sha256(self.Y)
        common = {"spec": spec, "mode": mode, "records_per_block": self.records_per_block,
                  "basis": enhancement_basis, "scale": self.enhancement_scale,
                  "y_digest": y_digest, "activity": self.activity}
        t0 = time.perf_counter()
        try:
            if transport == "local":
                ctx = mp.get_context("spawn")
                for device, (a, b) in zip(self.devices, self.blocks):
                    channel = tp.PipeChannel(ctx, device, dict(common, sub_rows=rows[a:b].tolist()),
                                             None if self.shm is None else self.shm.name, nf)
                    self.workers.append({"device": device, "channel": channel, "proc": channel.proc,
                                         "positions": np.arange(a, b)})
            else:
                a, b = self.blocks[0]
                local = dict(common, sub_rows=rows[a:b].tolist(),
                             elements=getattr(model, "_elements", None))
                channel = tp.InProcessChannel(self.devices[0], local, self.act)
                self.workers.append({"device": self.devices[0], "channel": channel,
                                     "positions": np.arange(a, b)})
                for rank in range(1, nworkers):
                    a, b = self.blocks[rank]
                    channel = tp.MpiChannel(self.comm, rank)
                    channel.send(("bind", dict(common, sub_rows=rows[a:b].tolist()), nf))
                    self.workers.append({"device": None, "channel": channel, "positions": np.arange(a, b)})
            self.ready = [self._receive(w, "bind", self.timeout_s + BIND_S_PER_RECORD *
                                        max(self.block_records)) for w in self.workers]
            for w, info in zip(self.workers, self.ready):
                w["device"] = info["device_id"]
            if transport == "mpi":
                self.devices = [w["device"] for w in self.workers]
        except BaseException:
            self.close()
            raise
        self.timing["bind_s"] = time.perf_counter() - t0

    @classmethod
    def from_problem(cls, problem, spec, **kwargs):
        return cls(spec, problem.M, problem.enhancement_basis if problem.nb else None,
                   problem.enhancement_scale, **kwargs)

    # ------------------------------------------------------------ messaging
    def _where(self, worker):
        return worker["channel"].label

    def _send(self, worker, message, what):
        try:
            worker["channel"].send(message)
        except (BrokenPipeError, OSError) as exc:
            self.broken = True
            raise BackendError("GPU worker %s is unreachable for %s (%s; exit code %s)" %
                               (self._where(worker), what, exc.__class__.__name__,
                                worker["channel"].exitcode)) from None

    def _receive(self, worker, what, timeout_s=None):
        channel = worker["channel"]
        timeout_s = self.timeout_s if timeout_s is None else float(timeout_s)
        deadline = time.monotonic() + timeout_s
        while not channel.poll(1.0):
            if not channel.alive():
                self.broken = True
                raise BackendError("GPU worker %s died during %s (exit code %s)" %
                                   (self._where(worker), what, channel.exitcode))
            if time.monotonic() > deadline:
                self.broken = True
                raise BackendError("GPU worker %s: no reply to %s within %.0f s" %
                                   (self._where(worker), what, timeout_s))
        try:
            reply = channel.recv()
        except EOFError:
            self.broken = True
            raise BackendError("GPU worker %s closed its pipe during %s" %
                               (self._where(worker), what)) from None
        if reply[0] == "error":
            self.broken = True
            raise BackendError("GPU worker %s failed during %s:\n%s" %
                               (self._where(worker), what, reply[1]))
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
        if self.activity == "host":
            self.act[:] = self.activity_linear(state)
        for w in self.workers:
            self._send(w, ("evaluate", state, request), "evaluate")
        if self.activity == "host" and self.transport == "mpi":
            self.comm.Bcast(self.act, root=0)
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

    def rebind_basis(self, enhancement_basis, enhancement_scale=1.0):
        """Change the template basis on every worker, keeping their records (SCALING_PLAN P1)."""
        if self.closed or self.broken:
            raise BackendError("multi-GPU backend is closed or broken")
        for w in self.workers:
            self._send(w, ("rebasis", enhancement_basis, enhancement_scale), "rebind_basis")
        self.timing["rebind_basis_s"] = [self._receive(w, "rebind_basis")[1] for w in self.workers]
        self.nbasis_enhancement = (0 if enhancement_basis is None
                                   else int(np.asarray(enhancement_basis).shape[1]))
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale) if enhancement_basis is not None else 1.0
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else bk.array_digest(enhancement_basis))

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
                "transport": self.transport, "activity": self.activity,
                "ranks": None if self.comm is None else int(self.comm.Get_size()),
                "partition": [{"device": w["device"], "worker": w["channel"].label, "rows": [a, b],
                               "records": r}
                              for w, (a, b), r in zip(self.workers, self.blocks, self.block_records)],
                "rows": len(self.rows),
                "rows_sha256": provenance.canonical_sha256(list(self.rows)),
                "enhancement_columns": self.nbasis_enhancement,
                "enhancement_scale": self.enhancement_scale,
                "enhancement_basis_sha256": self.enhancement_basis_sha256,
                "records_per_block": self.records_per_block,
                "workers_ready": getattr(self, "ready", None),
                "workers": workers, "timing_s": dict(self.timing)}

    def summary(self):
        return "%s (%s, %s transport, %s activity) on devices %s: %d points, records %s" % (
            self.name, self.mode, self.transport, self.activity, self.devices, len(self.rows),
            self.block_records)

    # ------------------------------------------------------------ lifecycle
    @property
    def released(self):
        return self.closed

    def close(self, grace_s=30.0):
        if self.closed:
            return
        self.closed = True
        deadline = time.monotonic() + grace_s
        for w in self.workers:                # all workers release at once
            w["channel"].request_close()
        for w in self.workers:
            w["channel"].shutdown(deadline)
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
