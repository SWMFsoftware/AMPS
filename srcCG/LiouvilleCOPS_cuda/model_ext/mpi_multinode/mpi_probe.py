#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/mpi_multinode/mpi_probe.py (on Voyager: /data001/highresolu_vt/cops_output/mpi_multinode/mpi_probe.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5af0b9f488060bb09ad027d74c0b9c34adcaf61753ad700b65c9f89ebfe5ea87
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test (no physics). MPI probe: per-rank host / rank / local rank / GPU (UUID, PCI bus) records, ping-pong within
#   and across nodes, Bcast of a 12M float64 field, the transport's request/reply pattern (2026-11-10).
# RUN: srun --mpi=pmix ... .venv/bin/python mpi_probe.py OUT_DIR, from jobs/mpi_probe.sbatch (job 123623 ->
#   mpi_multinode/probe-123623/)
# ==== end of header ====
"""MPI probe for the multi-node LiouvilleCOPS_cuda fits on Voyager (test script, 2026-11-10; no physics).

    srun --mpi=pmix --nodes=2 --ntasks-per-node=8 .venv/bin/python mpi_probe.py OUT_DIR [--devices 0,...,7]

Every rank writes OUT_DIR/rank<NNN>.json: host, MPI rank / size / local rank (the launcher variable that
cops_cuda/transport.py reads), CUDA_VISIBLE_DEVICES, the device transport.mpi_local_device() chooses, its UUID
(cops_cuda/environment.device_uuid: driver API, no context) and PCI bus id, the MPI library, mpi4py.
Rank 0 then times, with host (NumPy) buffers as the transport uses them:
  * ping-pong rank 0 <-> rank 1 (same node) and rank 0 <-> first rank of the other node (if any), raw buffers
    8 B .. 96 MB (comm.Send/Recv) and pickled objects (comm.send/recv) of the size of an evaluate reply;
  * MPI_Bcast from rank 0 of 12,000,000 float64 (the host activity field of a 12M mesh, 96 MB);
  * the transport's request/reply pattern: rank 0 sends a small pickled request to every rank, every rank
    answers with a pickled tuple of arrays (rows x columns of a Jacobian block), rank 0 collects.
Results: OUT_DIR/bench.json (rank 0), printed too. Nothing is written outside OUT_DIR.
"""

import json
import os
import pickle
import socket
import sys
import time
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))


def rank_info(comm, devices):
    from mpi4py import MPI
    import mpi4py
    from cops_cuda import transport, environment
    import cupy
    info = {"host": socket.gethostname(), "rank": comm.Get_rank(), "size": comm.Get_size(),
            "local_rank": transport.mpi_local_rank(comm), "pid": os.getpid(),
            "launcher_local_rank_variables": {k: os.environ.get(k) for k in transport.LOCAL_RANK_VARIABLES},
            "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "visible_device_count": cupy.cuda.runtime.getDeviceCount(),
            "mpi_library": MPI.Get_library_version().strip(), "mpi_version": MPI.Get_version(),
            "mpi4py": mpi4py.__version__,
            "slurm": {k: v for k, v in os.environ.items() if k.startswith("SLURM_") and k in (
                "SLURM_JOB_ID", "SLURM_STEP_ID", "SLURM_PROCID", "SLURM_LOCALID", "SLURM_NODEID",
                "SLURM_NTASKS", "SLURM_JOB_NUM_NODES", "SLURM_STEP_GPUS", "SLURM_JOB_GPUS", "SLURM_CPUS_PER_TASK",
                "SLURM_TOPOLOGY_ADDR")},
            "ompi_mca": {k: v for k, v in os.environ.items() if k.startswith("OMPI_MCA_") or k.startswith("UCX_")},
            "affinity": sorted(os.sched_getaffinity(0))}
    dev = transport.mpi_local_device(comm, devices)
    info["device"] = int(dev)
    info["device_uuid"] = environment.device_uuid(dev)
    try:
        info["device_pci_bus_id"] = cupy.cuda.runtime.deviceGetPCIBusId(dev)
    except Exception as exc:                       # recorded, not fatal
        info["device_pci_bus_id"] = "error: %s" % exc
    return info


def pingpong(comm, peer, nbytes, repeat):
    """Round-trip time / 2 for a raw buffer of nbytes between rank 0 and peer (rank 0's view)."""
    rank = comm.Get_rank()
    buf = np.zeros(max(nbytes // 8, 1), dtype=np.float64)
    times = []
    for _ in range(repeat + 1):
        comm.Barrier()
        t0 = time.perf_counter()
        if rank == 0:
            comm.Send(buf, dest=peer, tag=7)
            comm.Recv(buf, source=peer, tag=8)
        elif rank == peer:
            comm.Recv(buf, source=0, tag=7)
            comm.Send(buf, dest=0, tag=8)
        times.append((time.perf_counter() - t0) / 2)
    return times[1:]


def pingpong_obj(comm, peer, obj, repeat):
    rank = comm.Get_rank()
    times = []
    for _ in range(repeat + 1):
        comm.Barrier()
        t0 = time.perf_counter()
        if rank == 0:
            comm.send(obj, dest=peer, tag=9)
            comm.recv(source=peer, tag=10)
        elif rank == peer:
            o = comm.recv(source=0, tag=9)
            comm.send(o, dest=0, tag=10)
        times.append((time.perf_counter() - t0) / 2)
    return times[1:]


def stats(times, nbytes=None):
    t = np.array(times)
    out = {"median_s": float(np.median(t)), "min_s": float(t.min()), "max_s": float(t.max()), "n": len(t)}
    if nbytes:
        out["bandwidth_GBps_median"] = nbytes / out["median_s"] / 1e9
    return out


def main():
    out = Path(sys.argv[1])
    devices = None
    if "--devices" in sys.argv:
        devices = [int(d) for d in sys.argv[sys.argv.index("--devices") + 1].split(",")]
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    info = rank_info(comm, devices)
    out.mkdir(parents=True, exist_ok=True)
    (out / ("rank%03d.json" % rank)).write_text(json.dumps(info, indent=1, sort_keys=True))
    infos = comm.gather(info, root=0)
    hosts = comm.allgather(info["host"])
    bench = {}
    if rank == 0:
        print("probe: %d ranks on %s" % (size, sorted(set(hosts))), flush=True)
        for i in infos:
            print("probe: rank %2d host %s local %d CUDA_VISIBLE_DEVICES %s device %d %s %s affinity %d cpus" % (
                i["rank"], i["host"], i["local_rank"], i["CUDA_VISIBLE_DEVICES"], i["device"], i["device_uuid"],
                i["device_pci_bus_id"], len(i["affinity"])), flush=True)
        per_node = {}
        for i in infos:
            per_node.setdefault(i["host"], []).append(i["device_uuid"])
        bench["distinct_uuid_per_node"] = {h: len(set(u)) == len(u) for h, u in per_node.items()}
        bench["mpi_library"] = infos[0]["mpi_library"]
    # peers: same node (rank 1) and the first rank on another node
    other = next((r for r, h in enumerate(hosts) if h != hosts[0]), None)
    peers = {"same_node": 1 if size > 1 and hosts[1] == hosts[0] else None, "other_node": other}
    sizes = [8, 65536, 1 << 20, 16 << 20, 96 << 20]
    for label, peer in peers.items():
        if peer is None:
            continue
        for nb in sizes:
            t = pingpong(comm, peer, nb, 5 if nb >= (16 << 20) else 20)
            if rank == 0:
                bench["pingpong_%s_%dB" % (label, nb)] = stats(t, nb)
        # an evaluate reply of one block: pNG, pRG (100 rows) and 40 Jacobian column pairs
        reply = ("ok", np.ones(100), np.ones(100), {("a", k): (np.ones(100), np.ones(100)) for k in range(40)}, 0.5)
        t = pingpong_obj(comm, peer, reply, 20)
        if rank == 0:
            bench["pingpong_obj_%s_reply_%dB" % (label, len(pickle.dumps(reply)))] = stats(t, len(pickle.dumps(reply)))
    # Bcast of the host activity field (12M float64)
    act = np.full(12_000_000, float(rank), dtype=np.float64)
    times = []
    for _ in range(6):
        comm.Barrier()
        t0 = time.perf_counter()
        comm.Bcast(act, root=0)
        comm.Barrier()
        times.append(time.perf_counter() - t0)
    ok = bool(np.all(act == 0.0))
    oks = comm.gather(ok, root=0)
    if rank == 0:
        bench["bcast_12M_float64_96MB"] = dict(stats(times[1:], act.nbytes), all_ranks_correct=all(oks))
    # request/reply round (the transport's evaluate pattern), 20 rounds
    request = ("evaluate", {"a": np.ones(25), "T": np.ones(3), "kappa": 1.0, "b": np.ones(7)}, ("a", 40), False)
    times = []
    for _ in range(21):
        comm.Barrier()
        t0 = time.perf_counter()
        if rank == 0:
            for r in range(1, size):
                comm.send(request, dest=r, tag=11)
            for r in range(1, size):
                comm.recv(source=r, tag=12)
        else:
            comm.recv(source=0, tag=11)
            comm.send(("ok", np.ones(100), np.ones(100), {("a", k): (np.ones(100), np.ones(100)) for k in range(40)},
                       0.5), dest=0, tag=12)
        times.append(time.perf_counter() - t0)
    if rank == 0:
        bench["request_reply_round_all_ranks"] = stats(times[1:])
        for k, v in bench.items():
            print("probe: %s %s" % (k, json.dumps(v)), flush=True)
        (out / "bench.json").write_text(json.dumps(bench, indent=1, sort_keys=True))
    comm.Barrier()
    return 0


if __name__ == "__main__":
    sys.exit(main())
