#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/mpi_support.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/mpi_support.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 717e777df01f43c465aa379a3e0855c2b4011105b4a7212e2ff7aa225fdba7e3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (no physics). Per-rank wrapper for model_ext drivers over MPI (user decisions 2026-11-10 "one fit split
#   across both nodes", "make it also runnable on one node if only one node is available"): start / end record per rank
#   (host, rank, GPU UUID and PCI bus, affinity, peak RSS), node memory sampler, optional timing and message-size
#   instrumentation of the fitting process (observation only); runs the driver unchanged as __main__ in every rank.
# USED BY: jobs/mpi_fit.sbatch (copied into provenance-<job>/model_ext/; MPI tests 123629-123649).
# ==== end of header ====
"""Per-rank wrapper for model_ext drivers run over MPI (one fit split across nodes; 2026-11-10). No physics.

    [srun --mpi=pmix ...] .venv/bin/python -W ignore mpi_support.py <driver.py> <driver arguments ...>

Every process (every MPI rank, or the single process of a local-transport run) runs THIS file, which
  1. writes a start record  $COPS_MPI_RANK_DIR/rank<NNN>.json: host, MPI rank / size / local rank (the launcher
     variable cops_cuda/transport.py reads), the CUDA device the rank uses (transport.mpi_local_device with the
     driver's --devices, exactly as tools/gpu_fit.py chooses it), its UUID (cops_cuda/environment.device_uuid,
     driver API) and PCI bus id (driver API; no CUDA context is created), CUDA_VISIBLE_DEVICES, pid,
     oom_score_adj, CPU affinity, the SLURM_* step variables, the MPI library and mpi4py version, the driver path
     and SHA-256, the command; local transport: the same for every device of --devices (the spawned workers);
  2. on the first process of each node: a daemon thread samples /proc/meminfo every $COPS_MPI_MEM_SAMPLE_S
     seconds (default 10) into $COPS_MPI_RANK_DIR/mem_<host>.tsv (host memory per node);
  3. with COPS_MPI_INSTRUMENT=1, on the fitting process only (rank 0 / the local process): wraps
     cops_cuda.multi_gpu.MultiGpuBackend (bind, _evaluate, rebind_basis) and the channels' send/recv of
     cops_cuda.transport to time every evaluation (pressures or Jacobian), the slowest/fastest worker's own time
     (the rest is communication + host work + assembly) and the bytes of every request (pickled) and reply
     (array bytes), and the Bcast of the host activity field. Observation only: no argument, array or result
     is changed; the wrapped methods call the originals;
  4. runs the driver unchanged as __main__ (runpy; sys.argv = [driver] + arguments), so the driver's own start-up
     (argument rewriting, hook installation, gpu_fit.main) happens in EVERY process; gpu_fit.main sends ranks >= 1
     into transport.mpi_worker_main;
  5. at exit: $COPS_MPI_RANK_DIR/rank<NNN>.end.json (peak RSS of the process and of its waited-for children,
     MemAvailable, exit status) and, with instrumentation, comm_rank000.json (every call + medians) and one
     summary line "mpi_support: ..." in the log.
Environment: COPS_MPI_RANK_DIR (required: a directory, created), COPS_MPI_INSTRUMENT (0|1, default 0),
COPS_MPI_MEM_SAMPLE_S (default 10). The tracked code (tools/gpu_fit.py, cops_cuda/*) and the driver are not edited.
"""

import atexit
import datetime as dt
import hashlib
import json
import os
import pickle
import resource
import runpy
import socket
import sys
import threading
import time
from pathlib import Path

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
for _p in (str(APP),):
    if _p not in sys.path:
        sys.path.append(_p)

STATE = {"exit": None}


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def meminfo():
    out = {}
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            if key in ("MemTotal", "MemAvailable", "MemFree", "Cached", "Shmem"):
                out[key + "_kB"] = int(value.split()[0])
    except OSError:
        pass
    return out


def write_json(path, data):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, sort_keys=True, default=str))
    tmp.replace(path)


def head_args(argv):
    return argv[:argv.index("--")] if "--" in argv else list(argv)


def option(head, name):
    for k, token in enumerate(head):
        if token == name and k + 1 < len(head):
            return head[k + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def pci_bus_id(index):
    """PCI bus id of CUDA device ``index`` from the driver API (no context)."""
    import ctypes
    try:
        driver = ctypes.CDLL("libcuda.so.1")
        handle, buf = ctypes.c_int(), ctypes.create_string_buffer(64)
        if (driver.cuInit(0) == 0 and driver.cuDeviceGet(ctypes.byref(handle), int(index)) == 0
                and driver.cuDeviceGetPCIBusId(buf, 64, handle) == 0):
            return buf.value.decode()
    except OSError as exc:
        return "error: %s" % exc
    return None


def device_record(index):
    """Never fails the run: an error is recorded instead."""
    rec = {"device": int(index)}
    try:
        from cops_cuda import environment
        rec["uuid"] = environment.device_uuid(int(index))
    except Exception as exc:
        rec["uuid"] = None
        rec["uuid_error"] = "%s: %s" % (exc.__class__.__name__, exc)
    rec["pci_bus_id"] = pci_bus_id(index)
    return rec


# ------------------------------------------------------------------------------------- memory sampler (per node)
def memory_sampler(path, every):
    def run():
        with open(path, "a", buffering=1) as f:
            if f.tell() == 0:
                f.write("utc\tepoch\tMemTotal_kB\tMemAvailable_kB\tused_kB\tCached_kB\n")
            while True:
                m = meminfo()
                if m:
                    f.write("%s\t%.0f\t%d\t%d\t%d\t%d\n" % (utc_now(), time.time(), m["MemTotal_kB"],
                                                            m["MemAvailable_kB"],
                                                            m["MemTotal_kB"] - m["MemAvailable_kB"],
                                                            m.get("Cached_kB", 0)))
                time.sleep(every)
    t = threading.Thread(target=run, name="mpi_support_mem", daemon=True)
    t.start()


# ------------------------------------------------------------------------------------- instrumentation (rank 0)
def array_bytes(obj):
    """Bytes of the NumPy arrays inside a reply (tuples, lists, dicts), plus 8 per scalar."""
    import numpy as np
    if isinstance(obj, np.ndarray):
        return int(obj.nbytes)
    if isinstance(obj, (tuple, list)):
        return sum(array_bytes(o) for o in obj)
    if isinstance(obj, dict):
        return sum(array_bytes(o) for o in obj.values())
    if isinstance(obj, (int, float)):
        return 8
    return 0


def instrument(rank_dir):
    import numpy as np
    from cops_cuda import multi_gpu, transport
    log = {"binds": [], "rebinds": [], "calls": [], "bytes_by_kind": {}}
    current = {"out": 0, "in": 0}

    def count(direction, kind, nbytes):
        current[direction] += nbytes
        k = log["bytes_by_kind"].setdefault(kind, {"out": 0, "in": 0, "messages": 0})
        k[direction] += nbytes
        k["messages"] += 1

    def wrap_channel(cls):
        send0, recv0 = cls.send, cls.recv

        def send(self, message):
            kind = message[0]
            try:
                if kind == "bind":
                    kw = dict(message[1])
                    basis = kw.pop("basis", None)
                    n = len(pickle.dumps((kind, kw), protocol=pickle.HIGHEST_PROTOCOL)) + array_bytes(basis)
                elif kind == "rebasis":                     # (kind, basis, scale): the basis travels as raw buffers
                    n = len(pickle.dumps((kind, message[2]), protocol=pickle.HIGHEST_PROTOCOL)) + array_bytes(message[1])
                else:
                    n = len(pickle.dumps(tuple(message[:3]), protocol=pickle.HIGHEST_PROTOCOL))
                count("out", kind, n)
            except Exception:
                pass
            return send0(self, message)

        def recv(self):
            reply = recv0(self)
            try:
                count("in", "reply", array_bytes(reply))
            except Exception:
                pass
            return reply
        cls.send, cls.recv = send, recv

    wrap_channel(transport.MpiChannel)
    wrap_channel(transport.PipeChannel)
    B = multi_gpu.MultiGpuBackend
    init0, eval0, rebind0 = B.__init__, B._evaluate, B.rebind_basis

    def init(self, spec, model, *a, **k):
        t0 = time.perf_counter()
        init0(self, spec, model, *a, **k)
        log["binds"].append({"utc": utc_now(), "seconds": time.perf_counter() - t0, "rows": len(self.rows),
                             "transport": self.transport, "activity": self.activity, "workers": len(self.workers),
                             "block_records": list(self.block_records),
                             "workers_ready": [{k: r.get(k) for k in ("host", "device_id", "rows", "records", "bind_s",
                                                                       "dataset_load_s", "activity", "pid")}
                                               for r in getattr(self, "ready", [])]})

    def evaluate(self, state, request):
        current["out"] = current["in"] = 0
        t0 = time.perf_counter()
        result = eval0(self, state, request)
        total = time.perf_counter() - t0
        ws = list(self.timing.get("worker_s") or [])
        ncols = len(list(request.keys()))
        bcast = int(self.act.nbytes) if (self.activity == "host" and self.transport == "mpi") else 0
        shm = int(self.act.nbytes) if (self.activity == "host" and self.transport == "local") else 0
        log["calls"].append({"kind": "jacobian" if ncols else "pressures", "columns": ncols, "rows": len(self.rows),
                             "seconds": total, "worker_max_s": max(ws) if ws else None,
                             "worker_min_s": min(ws) if ws else None,
                             "rank0_worker_s": ws[0] if ws else None,
                             "overhead_s": (total - max(ws)) if ws else None,
                             "request_bytes": current["out"], "reply_bytes": current["in"],
                             "bcast_bytes_per_rank": bcast, "shared_memory_bytes": shm})
        return result

    def rebind(self, basis, scale=1.0):
        t0 = time.perf_counter()
        current["out"] = current["in"] = 0
        rebind0(self, basis, scale)
        log["rebinds"].append({"seconds": time.perf_counter() - t0, "request_bytes": current["out"],
                               "columns": None if basis is None else int(np.asarray(basis).shape[1])})
    B.__init__, B._evaluate, B.rebind_basis = init, evaluate, rebind

    # The optimizer's own evaluations: scipy.optimize.least_squares is replaced BEFORE fit_cops is imported (it does
    # "from scipy.optimize import least_squares"), so this wrapper is the innermost layer under gpu_fit's wall clock
    # and the drivers' wrappers: it sees exactly the residual / Jacobian calls of the optimizer and records the cost
    # and the SHA-256 of x, of the residual vector and of the Jacobian (bitwise comparison of runs).
    import scipy.optimize
    ls0 = scipy.optimize.least_squares
    log["least_squares"] = []

    def digest(a):
        return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

    def least_squares(fun, x0, jac="2-point", **kwargs):
        call = {"utc": utc_now(), "n_x": int(np.size(x0)), "x0_sha256": digest(np.asarray(x0, float)),
                "evaluations": []}
        log["least_squares"].append(call)

        def f(x, *a, **k):
            t0 = time.perf_counter()
            value = fun(x, *a, **k)
            r = np.asarray(value, dtype=float)
            call["evaluations"].append({"kind": "fun", "cost": float(0.5 * np.dot(r, r)),
                                        "cost_repr": repr(float(0.5 * np.dot(r, r))), "x_sha256": digest(x),
                                        "f_sha256": digest(r), "seconds": time.perf_counter() - t0})
            return value

        def j(x, *a, **k):
            t0 = time.perf_counter()
            value = jac(x, *a, **k)
            call["evaluations"].append({"kind": "jac", "x_sha256": digest(x),
                                        "J_sha256": digest(value) if isinstance(value, np.ndarray) else None,
                                        "seconds": time.perf_counter() - t0})
            return value
        result = ls0(f, x0, jac=j if callable(jac) else jac, **kwargs)
        call.update(nfev=int(getattr(result, "nfev", -1)), cost=float(getattr(result, "cost", float("nan"))),
                    x_sha256=digest(np.asarray(result.x, float)), status=int(getattr(result, "status", 0)))
        return result
    scipy.optimize.least_squares = least_squares

    def summary():
        out = {"binds": log["binds"], "rebinds": log["rebinds"], "bytes_by_kind": log["bytes_by_kind"]}
        for kind in ("pressures", "jacobian"):
            calls = [c for c in log["calls"] if c["kind"] == kind]
            later = calls[1:] if len(calls) > 1 else calls
            if not calls:
                continue

            def med(key, cs=later):
                v = [c[key] for c in cs if c[key] is not None]
                return float(np.median(v)) if v else None
            out[kind] = {"calls": len(calls), "total_s": float(sum(c["seconds"] for c in calls)),
                         "first_s": calls[0]["seconds"], "median_s_after_first": med("seconds"),
                         "median_worker_max_s": med("worker_max_s"), "median_worker_min_s": med("worker_min_s"),
                         "median_rank0_worker_s": med("rank0_worker_s"), "median_overhead_s": med("overhead_s"),
                         "median_request_bytes": med("request_bytes"), "median_reply_bytes": med("reply_bytes"),
                         "bcast_bytes_per_rank": calls[0]["bcast_bytes_per_rank"],
                         "shared_memory_bytes": calls[0]["shared_memory_bytes"]}
        return out

    def dump():
        s = summary()
        write_json(Path(rank_dir) / "comm_rank000.json", {"summary": s, "calls": log["calls"]})
        write_json(Path(rank_dir) / "optimizer_history.json", log["least_squares"])
        for c in log["least_squares"]:
            costs = [e["cost"] for e in c["evaluations"] if e["kind"] == "fun"]
            print("mpi_support: least_squares call: %d residual and %d Jacobian evaluations, costs %s" % (
                len(costs), sum(1 for e in c["evaluations"] if e["kind"] == "jac"),
                " ".join("%.12e" % v for v in costs[:12]) + (" ..." if len(costs) > 12 else "")), flush=True)
        for kind in ("pressures", "jacobian"):
            if kind in s:
                v = s[kind]
                print("mpi_support: %s: %d calls, median %.3f s after the first (first %.3f s); slowest worker %s s, "
                      "fastest %s s, rank 0 %s s, overhead %s s; request %s B, reply %s B, Bcast %d B per rank" % (
                          kind, v["calls"], v["median_s_after_first"] or 0.0, v["first_s"],
                          _f(v["median_worker_max_s"]), _f(v["median_worker_min_s"]), _f(v["median_rank0_worker_s"]),
                          _f(v["median_overhead_s"]), _f(v["median_request_bytes"], "%.0f"),
                          _f(v["median_reply_bytes"], "%.0f"), v["bcast_bytes_per_rank"]), flush=True)
        for b in s["binds"]:
            print("mpi_support: bind %.1f s (%s transport, %s activity, %d workers, %d rows)" % (
                b["seconds"], b["transport"], b["activity"], b["workers"], b["rows"]), flush=True)
    return dump


def _f(v, fmt="%.4f"):
    return "-" if v is None else fmt % v


# ------------------------------------------------------------------------------------- main
def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: mpi_support.py <driver.py> <driver arguments ...>")
    driver = str(Path(sys.argv[1]).resolve())
    args = sys.argv[2:]
    rank_dir = os.environ.get("COPS_MPI_RANK_DIR")
    if not rank_dir:
        raise SystemExit("mpi_support: set COPS_MPI_RANK_DIR (the per-rank provenance directory)")
    Path(rank_dir).mkdir(parents=True, exist_ok=True)
    head = head_args(args)
    mpi = option(head, "--transport") == "mpi"
    devices = option(head, "--devices")
    devices = [int(d) for d in devices.split(",")] if devices else None
    rec = {"utc_start": utc_now(), "host": socket.gethostname(), "pid": os.getpid(), "python": sys.executable,
           "driver": driver, "driver_sha256": sha256(driver), "wrapper": str(Path(__file__).resolve()),
           "wrapper_sha256": sha256(__file__), "argv": [driver] + args, "transport": "mpi" if mpi else "local",
           "cwd": os.getcwd(), "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
           "slurm": {k: v for k, v in os.environ.items() if k in (
               "SLURM_JOB_ID", "SLURM_STEP_ID", "SLURM_PROCID", "SLURM_LOCALID", "SLURM_NODEID", "SLURM_NTASKS",
               "SLURM_JOB_NUM_NODES", "SLURM_JOB_NODELIST", "SLURM_STEP_NODELIST", "SLURM_CPUS_PER_TASK",
               "SLURM_JOB_GPUS", "SLURM_STEP_GPUS", "SLURM_GPUS_ON_NODE")},
           "mpi_environment": {k: v for k, v in os.environ.items()
                               if k.startswith(("OMPI_MCA_", "UCX_", "PMIX_MCA_", "HCOLL_"))},
           "memory_start": meminfo()}
    try:
        rec["oom_score_adj"] = int(Path("/proc/self/oom_score_adj").read_text())
    except OSError:
        rec["oom_score_adj"] = None
    try:
        rec["cpu_affinity"] = sorted(os.sched_getaffinity(0))
    except OSError:
        rec["cpu_affinity"] = None
    first_on_node = True
    if mpi:
        from mpi4py import MPI
        import mpi4py
        from cops_cuda import transport
        comm = MPI.COMM_WORLD
        rank, size = comm.Get_rank(), comm.Get_size()
        local = transport.mpi_local_rank(comm)
        dev = transport.mpi_local_device(comm, devices)
        rec.update(rank=rank, size=size, local_rank=local,
                   local_rank_variables={k: os.environ.get(k) for k in transport.LOCAL_RANK_VARIABLES},
                   mpi_library=MPI.Get_library_version().strip(), mpi_version=list(MPI.Get_version()),
                   mpi4py=mpi4py.__version__, mpi4py_file=mpi4py.__file__,
                   role="fitter + block 0 (in process)" if rank == 0 else "worker (transport.mpi_worker_main)",
                   gpu=device_record(dev))
        first_on_node = local == 0
    else:
        rank, size = 0, 1
        rec.update(rank=0, size=1, local_rank=0, role="fitter (local transport: spawned workers)",
                   gpus=[device_record(d) for d in (devices or [0])])
    write_json(Path(rank_dir) / ("rank%03d.json" % rank), rec)
    if first_on_node:
        memory_sampler(str(Path(rank_dir) / ("mem_%s.tsv" % socket.gethostname().split(".")[0])),
                       float(os.environ.get("COPS_MPI_MEM_SAMPLE_S", "10")))
    dump = None
    if rank == 0 and os.environ.get("COPS_MPI_INSTRUMENT", "0") == "1":
        dump = instrument(rank_dir)

    def at_exit():
        end = {"utc_end": utc_now(), "rank": rank, "host": socket.gethostname(), "exit": STATE["exit"],
               "maxrss_self_kB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               "maxrss_children_kB": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
               "memory_end": meminfo()}
        try:
            write_json(Path(rank_dir) / ("rank%03d.end.json" % rank), end)
        except Exception:
            pass
        if dump is not None:
            try:
                dump()
            except Exception as exc:
                print("mpi_support: instrumentation summary failed: %s" % exc, flush=True)
    atexit.register(at_exit)
    sys.argv = [driver] + args
    try:
        runpy.run_path(driver, run_name="__main__")
        STATE["exit"] = 0
    except SystemExit as exc:
        STATE["exit"] = exc.code
        raise
    except BaseException as exc:
        STATE["exit"] = "%s: %s" % (exc.__class__.__name__, exc)
        raise


if __name__ == "__main__":
    main()
