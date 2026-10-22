"""Worker transports of the multi-GPU backend (SCALING_PLAN.md 3.3, phase P3).

One fit has one optimizer (the parent process, or MPI rank 0) and one ``CudaBackend`` per GPU
over a contiguous block of rows. ``multi_gpu.MultiGpuBackend`` talks to its workers through
channels with one interface (``send``, ``poll``, ``recv``, ``alive``, ``exitcode``):

* ``PipeChannel``: a spawned local process (one node; the C07 design, unchanged behaviour);
* ``MpiChannel``: an MPI rank (mpi4py), one rank per GPU, on any number of nodes. Ranks >= 1 run
  ``mpi_worker_main`` from the start and serve binds until rank 0 sends "exit";
* ``InProcessChannel``: rank 0's own GPU. Its work runs in the driver process when its reply is
  collected, which happens after the other ranks have their requests, so it overlaps them.

Every worker runs the same handler (``WorkerState.handle``), whatever the channel:
``bind``, ``evaluate``, ``rebasis``, ``seen``, ``identity``, ``close`` (release the backend,
keep serving) and ``exit`` (leave the loop). A failure inside a request is returned with its
traceback and raised by the driver as ``BackendError``; nothing is retried.

Activity field (SCALING_PLAN 3.5)
---------------------------------
* ``"host"``: the driver forms act_lin once per evaluation with the reference expression and every
  worker reads those bytes: shared memory (pipes), ``MPI_Bcast`` (MPI), the driver's array (in
  process). Bitwise one GPU for any partition.
* ``"device"``: every worker forms the activity of its own facets on its GPU
  (``kernels/activity.cu``); only the parameters travel. Deterministic and the same in every
  partition, within a few ulp of the host expression (``CudaBackend.activity_check``). Default
  of the MPI transport (user decision Q3, 2026-10-21).

MPI specifics
-------------
* Device of a rank: its rank within the node (the launcher's local-rank variable) modulo the
  visible devices, or ``devices[local_rank % len(devices)]``. Several ranks may share a GPU
  (tests). Under Slurm with ``--gpus-per-task=1`` every rank sees one device, its own.
* Arrays (template bases) travel as raw buffers in chunks of ``ARRAY_CHUNK`` elements, below the
  2**31 element count of MPI-3 libraries.
* A dead rank ends the whole MPI job (the launcher's behaviour); a reply timeout raises
  ``BackendError`` on rank 0, and ``tools/gpu_fit.py`` then calls ``MPI_Abort``.
"""

import os
import time
import traceback

import numpy as np

from . import provenance
from .errors import IdentityError

TAG_COMMAND, TAG_REPLY = 11, 12
ARRAY_CHUNK = 1 << 27                 # elements per MPI message (1 GiB of float64)
POLL_SLEEP_S = 2e-4


# ---------------------------------------------------------------- worker side
def bind_backend(spec, sub_rows, device_id, mode, records_per_block, basis, scale, y_digest,
                 activity="host", elements=None):
    """One ``CudaBackend`` for ``sub_rows`` of the driver's model, rebuilt from its frozen inputs."""
    from . import evaluator, formats, host_model, reference
    reference.bind()
    if elements is None:
        # Memory-mapped views (the file's bytes, as the reference reader's): 9.4 GB less per
        # worker at SHAP8 20M; checked against the driver's Y below.
        elements = formats.map_elements(os.path.join(spec["run"], "elements.bin"))
        if provenance.array_sha256(elements["Y"]) != y_digest:
            raise IdentityError("worker element cache differs from the driver's")
    model = host_model.header_model(spec["run"], spec["header"], spec["cols"],
                                    spec["rows"], elements, spec["coef"])
    model.load_kernels(np.asarray(sub_rows, dtype=int))
    return evaluator.CudaBackend(model, spec["run"], basis, scale, device_id=device_id,
                                 mode=mode, records_per_block=records_per_block, activity=activity)


class WorkerState:
    """The backend of one worker and the handler of the driver's requests."""

    def __init__(self, device_id):
        self.device_id = int(device_id)
        self.backend = None
        self.activity = None          # host activity buffer (global facets), when used

    def ready_info(self):
        be, ds = self.backend, self.backend.dataset
        return {"device_id": self.device_id, "rows": len(be.rows),
                "records": int(ds.nrecords),
                "streamed_records": int(ds.nrecords - ds.nresident) if ds.streamed else 0,
                "device_faces": int(be.dataset.nfaces_device),
                "bind_s": be.timing.get("bind_s"), "dataset_load_s": be.timing.get("dataset_load_s"),
                "activity": be.activity_mode, "pid": os.getpid(), "host": os.uname().nodename}

    def bind(self, kwargs, activity_buffer=None, copy=False):
        """Bind; with host activity the backend reads ``activity_buffer`` (a private copy per
        evaluation when ``copy``: the shared memory of the pipe transport)."""
        self.release()
        self.backend = bind_backend(device_id=self.device_id, **kwargs)
        if self.backend.activity_mode == "host":
            if activity_buffer is None:
                raise ValueError("host activity needs the driver's activity buffer")
            self.activity = activity_buffer
            self.backend.activity_source = ((lambda state: self.activity.copy()) if copy
                                            else (lambda state: self.activity))
        return ("ready", self.ready_info())

    def handle(self, message):
        """Reply to one request (never raises)."""
        kind = message[0]
        try:
            if self.backend is None and kind not in ("close",):
                raise RuntimeError("no backend bound for %r" % kind)
            if kind == "evaluate":
                _, state, request = message
                t0 = time.perf_counter()
                result = self.backend.evaluate(state, request)
                return ("ok", result.pNG, result.pRG, dict(result.columns), time.perf_counter() - t0)
            if kind == "rebasis":
                self.backend.rebind_basis(message[1], message[2])
                return ("ok", self.backend.timing.get("rebind_basis_s"))
            if kind == "seen":
                from . import masks
                return ("ok", masks.seen_positive(self.backend.dataset))
            if kind == "identity":
                return ("ok", self.backend.identity())
            if kind == "close":
                self.release()
                return ("closed",)
            raise ValueError("unknown request %r" % (kind,))
        except Exception:
            return ("error", traceback.format_exc())

    def release(self):
        if self.backend is not None:
            self.backend.release()
            self.backend = None
        self.activity = None


def pipe_worker_main(conn, bind_kwargs, device_id, shm_name, nfaces):
    """Spawned worker: bind, then serve requests until "close" (the process then exits)."""
    from multiprocessing import shared_memory
    worker, shm, act = WorkerState(device_id), None, None
    try:
        if shm_name is not None:
            shm = shared_memory.SharedMemory(name=shm_name)
            act = np.ndarray((nfaces,), dtype=np.float64, buffer=shm.buf)
        conn.send(worker.bind(bind_kwargs, act, copy=True))
    except BaseException:
        conn.send(("error", traceback.format_exc()))
        worker.release()
        return
    try:
        while True:
            message = conn.recv()
            if message[0] == "close":
                break
            conn.send(worker.handle(message))
    except (EOFError, KeyboardInterrupt):
        pass
    finally:
        worker.release()
        if shm is not None:
            del act
            shm.close()
        try:
            conn.send(("closed",))
        except (BrokenPipeError, OSError):
            pass


def mpi_worker_main(devices=None, comm=None):
    """Ranks >= 1: serve rank 0's requests until "exit" (returns the number of requests)."""
    from mpi4py import MPI
    comm = MPI.COMM_WORLD if comm is None else comm
    worker = WorkerState(mpi_local_device(comm, devices))
    served = 0
    try:
        while True:
            message = comm.recv(source=0, tag=TAG_COMMAND)
            kind = message[0]
            if kind == "exit":
                break
            served += 1
            if kind == "bind":
                kwargs = dict(message[1])                     # the basis follows as raw buffers
                kwargs["basis"] = recv_array(comm, 0, TAG_COMMAND) if message[2] else None
                buffer = (np.empty(message[3], dtype=np.float64)
                          if kwargs.get("activity") == "host" else None)
                try:
                    reply = worker.bind(kwargs, buffer)
                except Exception:
                    worker.release()
                    reply = ("error", traceback.format_exc())
            elif kind == "rebasis":
                basis = recv_array(comm, 0, TAG_COMMAND) if message[1] else None
                reply = worker.handle(("rebasis", basis, message[2]))
            elif kind == "evaluate":
                if message[3]:                                # host activity: collective first
                    comm.Bcast(worker.activity, root=0)
                reply = worker.handle(message[:3])
            else:
                reply = worker.handle(message)
            comm.send(reply, dest=0, tag=TAG_REPLY)
    finally:
        worker.release()
    return served


def mpi_release_workers(comm=None):
    """Rank 0, at the end of the program: every worker rank leaves ``mpi_worker_main``."""
    from mpi4py import MPI
    comm = MPI.COMM_WORLD if comm is None else comm
    for rank in range(1, comm.Get_size()):
        comm.send(("exit",), dest=rank, tag=TAG_COMMAND)


LOCAL_RANK_VARIABLES = ("OMPI_COMM_WORLD_LOCAL_RANK", "SLURM_LOCALID", "MPI_LOCALRANKID",
                        "PMI_LOCAL_RANK", "MV2_COMM_WORLD_LOCAL_RANK")


def mpi_local_rank(comm):
    """This rank's index within its node, from the launcher (Open MPI, Slurm, MPICH, MVAPICH).

    No collective call: ranks reach this at different times (rank 0 when a backend binds, the
    workers at start-up). Without a launcher variable the global rank is used (one node)."""
    for name in LOCAL_RANK_VARIABLES:
        value = os.environ.get(name)
        if value not in (None, ""):
            return int(value)
    return int(comm.Get_rank())


def mpi_local_device(comm, devices=None):
    """This rank's device: rank within the node modulo the visible (or given) devices."""
    import cupy
    local = mpi_local_rank(comm)
    if devices:
        return int(devices[local % len(devices)])
    return local % cupy.cuda.runtime.getDeviceCount()


# ---------------------------------------------------------------- arrays over MPI
def send_array(comm, array, dest, tag):
    array = np.ascontiguousarray(array)
    comm.send((array.shape, array.dtype.str), dest=dest, tag=tag)
    flat = array.reshape(-1)
    for start in range(0, flat.size, ARRAY_CHUNK):
        comm.Send(flat[start:start + ARRAY_CHUNK], dest=dest, tag=tag)


def recv_array(comm, source, tag):
    shape, dtype = comm.recv(source=source, tag=tag)
    out = np.empty(shape, dtype=np.dtype(dtype))
    flat = out.reshape(-1)
    for start in range(0, flat.size, ARRAY_CHUNK):
        comm.Recv(flat[start:start + ARRAY_CHUNK], source=source, tag=tag)
    return out


# ---------------------------------------------------------------- driver side channels
class PipeChannel:
    """A spawned worker process."""

    kind = "pipe"

    def __init__(self, ctx, device_id, bind_kwargs, shm_name, nfaces):
        parent, child = ctx.Pipe()
        self.proc = ctx.Process(target=pipe_worker_main, name="cops-gpu%d" % device_id,
                                args=(child, bind_kwargs, device_id, shm_name, nfaces), daemon=True)
        self.proc.start()
        child.close()
        self.conn, self.label = parent, "device %d" % device_id

    def send(self, message):
        self.conn.send(message)

    def poll(self, timeout_s):
        return self.conn.poll(timeout_s)

    def recv(self):
        return self.conn.recv()

    def alive(self):
        return self.proc.is_alive()

    @property
    def exitcode(self):
        return self.proc.exitcode

    def request_close(self):
        try:
            self.conn.send(("close",))
        except (BrokenPipeError, OSError):
            pass

    def shutdown(self, deadline):
        self.proc.join(max(0.0, deadline - time.monotonic()))
        if self.proc.is_alive():
            self.proc.terminate()
            self.proc.join(5.0)
        self.conn.close()


class MpiChannel:
    """A worker rank; rank 0 is the driver."""

    kind = "mpi"

    def __init__(self, comm, rank):
        self.comm, self.rank, self.label = comm, int(rank), "rank %d" % rank
        self.exitcode = None
        self.host_activity = False
        self.outstanding = 0                  # requests sent and not yet answered

    def send(self, message):
        kind = message[0]
        self.outstanding += 1
        if kind == "evaluate":                # the driver's MPI_Bcast of act_lin follows (host)
            self.comm.send(tuple(message[:3]) + (self.host_activity,), dest=self.rank, tag=TAG_COMMAND)
        elif kind == "bind":                  # ("bind", kwargs with basis, nfaces)
            kwargs = dict(message[1])
            self.host_activity = kwargs.get("activity") == "host"
            basis = kwargs.pop("basis")
            self.comm.send(("bind", kwargs, basis is not None, message[2]), dest=self.rank, tag=TAG_COMMAND)
            if basis is not None:
                send_array(self.comm, basis, self.rank, TAG_COMMAND)
        elif kind == "rebasis":
            self.comm.send(("rebasis", message[1] is not None, message[2]), dest=self.rank, tag=TAG_COMMAND)
            if message[1] is not None:
                send_array(self.comm, message[1], self.rank, TAG_COMMAND)
        else:
            self.comm.send(message, dest=self.rank, tag=TAG_COMMAND)

    def poll(self, timeout_s):
        deadline = time.monotonic() + timeout_s
        while not self.comm.Iprobe(source=self.rank, tag=TAG_REPLY):
            if time.monotonic() > deadline:
                return False
            time.sleep(POLL_SLEEP_S)
        return True

    def recv(self):
        reply = self.comm.recv(source=self.rank, tag=TAG_REPLY)
        self.outstanding -= 1
        return reply

    def alive(self):
        return True                           # a dead rank ends the MPI job

    def request_close(self):
        """Release the rank's backend; the rank keeps serving later binds."""
        try:
            self.send(("close",))
        except Exception:
            pass

    def shutdown(self, deadline):
        """Collect every unanswered reply (after a failure: stale ones too), the last "closed"."""
        try:
            while self.outstanding > 0 and self.poll(max(1.0, deadline - time.monotonic())):
                self.recv()
        except Exception:
            pass


class InProcessChannel:
    """The driver's own GPU (MPI rank 0): requests run when their reply is collected."""

    kind = "in_process"

    def __init__(self, device_id, bind_kwargs, activity_buffer):
        self.worker = WorkerState(device_id)
        self.label = "rank 0 (device %d, in process)" % device_id
        self.exitcode = None
        self._pending = [("bind", bind_kwargs, activity_buffer)]

    def send(self, message):
        self._pending.append(message)

    def poll(self, timeout_s):
        return True

    def recv(self):
        message = self._pending.pop(0)
        if message[0] == "bind":
            try:
                return self.worker.bind(message[1], message[2])
            except Exception:
                self.worker.release()
                return ("error", traceback.format_exc())
        return self.worker.handle(message)

    def alive(self):
        return True

    def request_close(self):
        self._pending.clear()

    def shutdown(self, deadline):
        self.worker.release()
