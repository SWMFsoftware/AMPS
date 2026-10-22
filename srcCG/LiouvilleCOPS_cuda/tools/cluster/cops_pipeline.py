#!/usr/bin/env python3
"""Fitting campaign for one surface model, stage by stage and resumable (Voyager runner; C15/P0).

    .venv/bin/python tools/cluster/cops_pipeline.py plan|run|submit|status --restart-dir DIR \\
        --shape MESH --surface-bc flat|rough [--max-nfev 2000] [--fine-tune] [options]

The five campaign settings (user request 2026-10-22):

* ``--shape``: the surface model (NASTRAN ``.bdf``, SPICE DSK ``.bds``/``.dsk``; ``--shape-format``
  for other names).
* ``--surface-bc``: the angular emission law of the fitting base, ``SurfaceBC`` of the unchanged
  ``../LiouvilleCOPS/tools/rebuild_cops_boundary.py``: ``flat`` (Lambertian) or ``rough`` (uniform over the
  hemisphere). Each law has its own base and therefore its own restart directory.
* ``--restart-dir``: every product of the campaign, the stage stamps, the logs and the fits'
  ``fit_restart.npz`` files. Running the same command again resumes at the first incomplete stage;
  ``run --restart-dir DIR`` alone repeats the stored request.
* ``--max-nfev``: the optimizer evaluation limit of each fit (default 2000, the production setting).
* ``--fine-tune`` / ``--no-fine-tune`` (default): late-interval fine tuning of each fit's restart file
  (``tools/fine_tune_gpu.py``, quick_guide §4f; defaults ``--ft-minutes 5 --ft-intervals 5 --ft-method cpu``).

Stages, in this order (``plan`` prints the exact commands):

* ``base`` (group prepare): ``rebuild_cops_boundary.py --executable native/liouville_cops_cuda``, one GPU.
  Or ``--base DIR``: a base built elsewhere (copied from another machine), accepted after an identity check
  (``rebuild.json`` complete with the same SurfaceBC and mesh, ``run_manifest.json`` process complete and
  cache verified, ``elements.bin`` header);
* ``store``: ``tools/c15_record_store.py``, mixed precision only (or ``--store DIR``: its
  ``store_manifest.json`` must name this base's ``run_manifest.json`` digest);
* ``visibility``: ``tools/c05_write_bundle.py --visibility-only`` (needed by N30);
* ``gravity``, ``cache``, ``bundle``: ``c04_gravity_equivalence.py gpu`` on all GPUs of the node,
  ``c04_write_gpu_cache.py``, ``c05_write_bundle.py`` with the cache (needed by C30 and NC31);
* ``translate:<FIT>``: ``restart_tool.py translate`` of a ``--restart-from`` file of another mesh, only
  with ``--translate-restart``;
* ``fit:<FIT>`` (group fit): ``tools/gpu_fit.py``; one node: worker processes on its GPUs; several nodes
  (or ``--transport mpi``): one MPI rank per GPU (``srun``/``mpirun``). Output ``fits/<FIT>_n<N>/segNN``;
* ``finetune:<FIT>`` (group finetune): ``tools/fine_tune_gpu.py`` on the GPUs of one node.

Fits (``--fits``, default all four, templates exactly as in the 12M/20M chains): B29 (no template), N30
(``late_ng_not_rg`` of the visibility bundle), C30 (``cliff_ge_60deg`` of the full bundle), NC31 (both,
from the full bundle). The fitter arguments are those of the production chains (``DEFAULT_FIT_ARGS``).

Hardware (user decision 2026-10-22: nothing hardwired). ``run`` detects the nodes, GPUs per node, ranks,
GPU memory per device and host memory per node: Slurm variables, ``CUDA_VISIBLE_DEVICES``, ``nvidia-smi``,
the CUDA device properties (CuPy, in a child process), ``/proc/meminfo``; ``--nodes``, ``--gpus-per-node``,
``--gpu-memory-gb`` and ``--host-memory-gb`` override them. All GPUs are assumed alike. ``plan`` and
``status`` never touch a GPU: without the overrides they report the GPU memory as unknown.

Precision ``auto`` (default; user decision 2026-10-22), decided once the record count is known (after the
base stage) and recorded with the hardware in ``pipeline.json``:

* ``fp64`` when the FP64 device records (94 B each) fit resident on all fit GPUs with the headroom;
* otherwise ``mixed``: FP32 records (20 B), FP64 arithmetic, as many records resident as fit and the rest
  streamed from pinned host memory with copies overlapped with computation (the P2 engine; each MPI rank
  streams its own rows on its own node). Streaming exists for mixed precision only.
* The pinned host buffers of the streamed records, the GPU workers and the rank-0 driver must fit the host
  memory of a node (``HOST_USABLE_FRACTION``); otherwise the run is refused with the numbers. The check
  is repeated with the fit job's hardware before every fit.

Gravity precision is a separate policy (user decision 2026-10-22): ``--gravity-precision mixed`` (default;
``mixed_near_fp64_exact_far_fp32_quad3_v1``: exact FP64 terms of the faces within 6 circumradii, FP32
quadrature for the far faces) or ``fp64`` (``fp64_reference_v1``, the exact reference; about 7 h for 3M on
two A6000). It is never derived from ``--precision``, whose ``mixed`` is another policy
(``fp32_records_fp64_arithmetic_v1``: FP32 storage, FP64 arithmetic). ``pipeline.json`` records both.
``plan`` estimates the gravity time per precision from recorded gravity runs (``gpu_manifest.json`` in the
restart directory and under ``$COPS_OUTPUT_ROOT/cuda_c04``), scaled by facet pairs and GPUs; the GPU model
of those runs is printed, and without any run the time is reported as unknown.

Wall time (user decision 2026-10-22; default no limit). ``--wall-time HH:MM:SS`` (from the start of
``run``) or ``--wall-time slurm`` (the job's end) with ``--wall-margin`` (default 300 s). Fits get the
deadline (``gpu_fit.py --deadline``): a fit stops at its best evaluated state before the deadline, writes
an interrupted restart file and exits 75; the next ``run`` continues it as segment k+1 with
``--restart-from`` that file and the remaining evaluations. A segmented fit is not identical to one
continuous run. Other stages start only when their recorded duration (an earlier attempt) still fits;
a killed stage runs again from its start. ``run`` exits 75 when it stopped for the deadline.

Resume rules. A stage is complete when its stamp (``DIR/stamps/<stage>.json``) says so, the stamp's inputs
equal the current ones and its outputs exist. A stage whose stamp is ``running``, ``failed`` or
``interrupted`` (or whose output exists without a stamp) is moved aside (``<output>.incomplete_<UTC>``,
never deleted) and run again. A completed stage whose inputs changed only because an earlier product was
made again is stale and runs again; other changed inputs are refused, as is a restart directory whose
``pipeline.json`` has another mesh, emission law, fitting input, gravity body, genus or precision request.

``submit`` writes one batch script per group into ``DIR/slurm/`` (``tools/cluster/cops_pipeline.slurm``;
SITE lines for the partition, account and modules) from ``--nodes`` and ``--gpus-per-node`` (default
``auto``: the partition's GRES from ``sinfo``) and prints the ``sbatch`` commands chained with ``afterok``.
With ``--resubmit N`` (needs ``--wall-time``) each job continues itself instead: when ``status --check``
says its group is incomplete but resumable it resubmits its own script with ``afterany`` (at most N times),
and when the group is complete it submits the next group. ``submit`` never runs ``sbatch`` itself.
"""

import argparse
import contextlib
import datetime as dt
import errno
import hashlib
import json
import os
from pathlib import Path
import shlex
import signal
import socket
import struct
import subprocess
import sys
import time

APP_ROOT = Path(__file__).resolve().parents[2]
SLURM_TEMPLATE = Path(__file__).resolve().parent / "cops_pipeline.slurm"
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))

from cops_cuda import wall_time  # noqa: E402  (standard library and NumPy only)

SCHEMA = "liouville_cops_cuda.pipeline"
SCHEMA_VERSION = 1
STAMP_SCHEMA = "liouville_cops_cuda.pipeline_stage"

# Fits and their templates: (bundle, column). As in the 12M/20M chains (scratchpad chain_020m.sh).
FIT_TEMPLATES = {
    "B29": (),
    "N30": (("visibility", "late_ng_not_rg"),),
    "C30": (("bundle", "cliff_ge_60deg"),),
    "NC31": (("bundle", "late_ng_not_rg"), ("bundle", "cliff_ge_60deg")),
}
DEFAULT_FITS = ("B29", "N30", "C30", "NC31")
DEFAULT_FIT_ARGS = ("--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 "
                    "--tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6")
OWNED_FIT_ARGS = ("--out", "--max-nfev", "--activity-template", "--restart-from")
DEFAULT_MAX_NFEV = 2000
GROUPS = ("prepare", "fit", "finetune")
FINE_TUNE_DEFAULTS = {"minutes": 5.0, "intervals": 5, "method": "cpu", "max_nfev": 1000}
THREAD_ENVIRONMENT = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                      "PYTHONUNBUFFERED": "1"}
SRUN_ARGS = "--mpi=pmix --kill-on-bad-exit=1"
SUBMIT_TIMES = {"prepare": "12:00:00", "fit": "08:00:00", "finetune": "04:00:00"}
CHUNK = 8 * 1024 * 1024
EXIT_FAILED, EXIT_REFUSED, EXIT_RESUMABLE = 1, 2, wall_time.EXIT_WALL_TIME

# Memory model of the precision decision (evaluator._resident_records; DeviceDataset layouts).
KERNEL_HEADER_BYTES, KERNEL_RECORD_BYTES = 92, 40      # CPU kernel files (point=NNNN.bin)
FP64_DEVICE_BYTES_PER_RECORD = 94      # face, Kn, KPhi, KV, cos, shadow (37) + 6 static FP64 factors (48)
#                                        + lit (1) + the bind-time phiC scratch (8)
FP32_DEVICE_BYTES_PER_RECORD = 20      # record_store.BYTES_PER_RECORD
DEVICE_BYTES_PER_FACE = {"host": 286, "device": 558}   # per-face arrays for all faces, 33 basis columns;
#                                        device activity adds Y and 8 template columns (264) and the check (8)
DEVICE_FIXED_BYTES = (4 << 26) + (64 << 20) + (512 << 20)   # block tables, slack, CUDA context
STREAM_SLOT_BYTES = 3 * (1 << 25) * FP32_DEVICE_BYTES_PER_RECORD   # STREAM_SLOTS x STREAM_CHUNK_RECORDS
PARTITION_IMBALANCE = 0.03             # record partition within 2.1% of ideal at 24 ranks on 20M (P3)
DEFAULT_HEADROOM = 0.20                # memory.DEFAULT_HEADROOM_FRACTION
PINNED_BYTES_PER_STREAMED_RECORD = 28  # measured 2026-10-22, 20M on 2 A6000: 21.2 GB pinned per worker
#                                        for 0.761e9 streamed records (15.2 GB of FP32 records)
HOST_DRIVER_BYTES_PER_FACE = 500       # measured 2026-10-22, 20M: rank-0 driver 9.8 GB
HOST_WORKER_BYTES = 1.5e9              # measured 2026-10-22: GPU worker anonymous memory 1.0-1.4 GB
HOST_USABLE_FRACTION = 0.85            # of a node's memory for the fit's processes
# Precision policies (cops_cuda/precision.py POLICY_NAME, MIXED_POLICY_NAME; cops_cuda/gravity.py MIXED_POLICY).
RECORD_POLICIES = {"fp64": "fp64_reference_v1", "mixed": "fp32_records_fp64_arithmetic_v1"}
GRAVITY_POLICIES = {"fp64": "fp64_reference_v1", "mixed": "mixed_near_fp64_exact_far_fp32_quad3_v1"}
DEFAULT_GRAVITY_PRECISION = "mixed"    # user decision 2026-10-22; fp64 is the option
GRAVITY_NEAR_FACTOR = 6.0


class PipelineError(Exception):
    """A refusal or failure the user must resolve (printed without a traceback)."""


class Interrupted(BaseException):
    """SIGTERM (the Slurm time limit, scancel) while a stage runs."""


# ----------------------------------------------------------------------------- small helpers
def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def utc_stamp():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def clock_text():
    return time.strftime("%a %H:%M:%S")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def quick_identity(path):
    """Path, size and mtime of a large product (compared instead of a digest)."""
    path = Path(path)
    st = path.stat()
    return {"path": str(path), "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns}


def file_identity(path, stored=None, count_facets=False):
    """Path, size, mtime and SHA-256 of an input file.

    The digest of ``stored`` is reused when path, size and mtime are unchanged. With ``count_facets``
    the NASTRAN ``CTRIA3`` cards are counted in the same pass, and the header comment of a mesh whose
    gravity topology admission failed (convert_ply_to_bdf.py) is noted.
    """
    path = Path(path).resolve()
    st = path.stat()
    if (stored and stored.get("path") == str(path) and stored.get("size_bytes") == st.st_size
            and stored.get("mtime_ns") == st.st_mtime_ns and stored.get("sha256")
            and (not count_facets or "facets" in stored)):
        return dict(stored)
    digest = hashlib.sha256()
    facets, tail, head = 0, b"\n", b""
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            digest.update(block)
            if count_facets:
                if len(head) < 4096:
                    head += block[:4096 - len(head)]
                joined = tail + block
                facets += joined.count(b"\nCTRIA3")
                tail = joined[-6:]          # one byte short of a card: a match is never counted twice
    identity = {"path": str(path), "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns,
                "sha256": digest.hexdigest()}
    if count_facets:
        identity["facets"] = facets
        identity["topology_admission_failed"] = b"gravity topology admission FAILED" in head
    return identity


def write_json_atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.%d.tmp" % (path.name, os.getpid()))
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(str(temporary), str(path))


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def move_aside(path):
    """Rename ``path`` to ``<path>.incomplete_<UTC>`` (never delete); returns the new path."""
    path = Path(path)
    stamp = utc_stamp()
    target = path.with_name("%s.incomplete_%s" % (path.name, stamp))
    k = 1
    while os.path.lexists(str(target)):
        k += 1
        target = path.with_name("%s.incomplete_%s_%d" % (path.name, stamp, k))
    os.rename(str(path), str(target))
    return target


def elements_faces(base):
    """Facet count from the 12-byte ``elements.bin`` header (version 1, nfaces, 25)."""
    with (Path(base) / "elements.bin").open("rb") as stream:
        header = stream.read(12)
    if len(header) != 12:
        raise PipelineError("truncated elements.bin header in %s" % base)
    version, nfaces, nbasis = struct.unpack("<iii", header)
    if version != 1 or nfaces <= 0 or nbasis != 25:
        raise PipelineError("unsupported elements.bin header in %s" % base)
    return nfaces


def kernel_records(base):
    """Kernel records of a base from the kernel file sizes (no file is read)."""
    total = 0
    for entry in os.scandir(str(Path(base) / "kernels")):
        if entry.name.startswith("point=") and entry.name.endswith(".bin"):
            total += (entry.stat().st_size - KERNEL_HEADER_BYTES) // KERNEL_RECORD_BYTES
    return total


def parse_input(path):
    """Key = Value lines with # comments; repeated keys keep every value (rebuild_cops_boundary rules)."""
    values = {}
    for line in Path(path).read_text().splitlines():
        active = line.split("#", 1)[0].strip()
        if active and "=" in active:
            key, value = (s.strip() for s in active.split("=", 1))
            values.setdefault(key, []).append(value)
    return values


def input_files_missing(fit_input, reference):
    """Data files of the base input that do not exist here (observations, SPICE kernels).

    Relative paths resolve as ``rebuild_cops_boundary.py`` resolves them: from the reference tree, and
    ``Kernel`` entries from ``KernelDirectory``. ``ShapeModel`` is replaced by ``--shape``.
    """
    values = parse_input(fit_input)

    def rooted(value, base):
        p = Path(value)
        return p if p.is_absolute() else Path(base) / p
    missing = []
    for value in values.get("ObservationFile", [])[-1:]:
        if not rooted(value, reference).is_file():
            missing.append("ObservationFile = %s" % value)
    kernel_root = rooted(values.get("KernelDirectory", ["."])[-1], reference)
    if "KernelDirectory" in values and not kernel_root.is_dir():
        missing.append("KernelDirectory = %s" % values["KernelDirectory"][-1])
    for value in values.get("Kernel", []):
        if not rooted(value, kernel_root).is_file():
            missing.append("Kernel = %s" % value)
    return missing


def shape_format(path, explicit):
    if explicit:
        return explicit
    suffix = Path(path).suffix.lower()
    if suffix in (".bdf", ".nas"):
        return "nastran"
    if suffix in (".bds", ".dsk"):
        return "dsk"
    raise PipelineError("cannot tell the format of %s; give --shape-format nastran|dsk" % path)


def restart_meta(path):
    """Mesh identity and template names of a fit_restart.npz (NumPy; meta only)."""
    import numpy as np
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    return {"mesh_identity": meta["mesh"]["identity"], "faces": meta["mesh"].get("faces"),
            "templates": list(meta.get("templates", {}).get("names", []))}


def gb(value):
    return "%.1f GB" % (value / 1e9)


# ----------------------------------------------------------------------------- hardware
def run_quiet(command, timeout=60):
    try:
        return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                              universal_newlines=True, timeout=timeout).stdout
    except (OSError, subprocess.SubprocessError):
        return None


def nvidia_smi_memory():
    """Total memory of each visible GPU (bytes), from nvidia-smi (no CUDA context), or None."""
    out = run_quiet(["nvidia-smi", "--query-gpu=memory.total", "--format=csv,noheader,nounits"])
    values = [int(v.strip()) * (1 << 20) for v in (out or "").splitlines() if v.strip().isdigit()]
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if values and visible:
        picks = [int(d) for d in visible.split(",") if d.strip().isdigit()]
        values = [values[d] for d in picks if d < len(values)] or values
    return values or None


def cuda_memory(python):
    """Total memory of each visible GPU from the CUDA device properties (CuPy, child process), or None."""
    code = ("import cupy\nr = cupy.cuda.runtime\n"
            "print(' '.join(str(r.getDeviceProperties(i)['totalGlobalMem']) for i in range(r.getDeviceCount())))")
    out = run_quiet([python, "-c", code], timeout=120)
    values = [int(v) for v in (out or "").split() if v.isdigit()]
    return values or None


def host_memory():
    """(bytes, source): Slurm's per-node memory when set, else MemTotal."""
    value = os.environ.get("SLURM_MEM_PER_NODE", "")
    if value.isdigit() and int(value) > 0:
        return int(value) * (1 << 20), "SLURM_MEM_PER_NODE"
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024, "/proc/meminfo"
    except OSError:
        pass
    return None, "unknown"


def detect_hardware(args, request, live):
    """Nodes, GPUs per node, ranks, GPU and host memory, each with its source.

    ``live`` (``run``): may ask nvidia-smi and the CUDA runtime; otherwise only options, the
    environment and /proc/meminfo are read (``plan``/``status`` touch no GPU).
    """
    hw = {}
    env_nodes = os.environ.get("SLURM_JOB_NUM_NODES") or os.environ.get("SLURM_NNODES") or ""
    if args.nodes:
        hw["nodes"], hw["nodes_source"] = args.nodes, "option"
    elif env_nodes.isdigit():
        hw["nodes"], hw["nodes_source"] = int(env_nodes), "slurm"
    else:
        hw["nodes"], hw["nodes_source"] = 1, "default"
    gpus, source = None, None
    if isinstance(args.gpus_per_node, int):
        gpus, source = args.gpus_per_node, "option"
    elif os.environ.get("SLURM_GPUS_ON_NODE", "").isdigit() and int(os.environ["SLURM_GPUS_ON_NODE"]) > 0:
        gpus, source = int(os.environ["SLURM_GPUS_ON_NODE"]), "SLURM_GPUS_ON_NODE"
    elif os.environ.get("CUDA_VISIBLE_DEVICES"):
        gpus = len([d for d in os.environ["CUDA_VISIBLE_DEVICES"].split(",") if d.strip()])
        source = "CUDA_VISIBLE_DEVICES"
    memory = None
    if args.gpu_memory_gb:
        memory, msource = int(args.gpu_memory_gb * 1e9), "option"
    elif live:
        values, msource = nvidia_smi_memory(), "nvidia-smi"
        if values is None:
            values, msource = cuda_memory(request["python"]), "CUDA device properties"
        if values:
            memory = min(values)
            if gpus is None:
                gpus, source = len(values), msource
    if gpus is None:
        gpus, source = None, "unknown (give --gpus-per-node)"
    hw["gpus_per_node"], hw["gpus_source"] = gpus, source
    hw["gpu_memory_bytes"], hw["gpu_memory_source"] = memory, (msource if memory else "unknown")
    if args.host_memory_gb:
        hw["host_memory_bytes"], hw["host_memory_source"] = int(args.host_memory_gb * 1e9), "option"
    else:
        hw["host_memory_bytes"], hw["host_memory_source"] = host_memory()
    tasks = os.environ.get("SLURM_NTASKS", "")
    if tasks.isdigit() and int(tasks) > 1 and not args.nodes:
        hw["ranks"] = int(tasks)
    else:
        hw["ranks"] = hw["nodes"] * (gpus or 1)
    cpus = os.environ.get("SLURM_CPUS_ON_NODE", "")
    hw["threads"] = args.threads or (int(cpus) if cpus.isdigit() else (os.cpu_count() or 1))
    hw["launcher"] = args.launcher or ("srun" if os.environ.get("SLURM_JOB_ID") else "mpirun")
    hw["host"] = socket.gethostname()
    return hw


def sinfo_gres(partition=None):
    """(GPUs per node, CPUs per node) of a partition from ``sinfo``; (None, None) when unknown."""
    extra = ["-p", partition] if partition else []
    out = run_quiet(["sinfo", "-h", "-o", "%G|%c"] + extra)
    counts, cpus = set(), set()
    for line in (out or "").splitlines():
        gres, _, cpu = line.partition("|")
        for item in gres.split(","):
            fields = item.split("(")[0].split(":")
            if fields[0] == "gpu" and fields[-1].isdigit():
                counts.add(int(fields[-1]))
        if cpu.strip().isdigit():
            cpus.add(int(cpu.strip()))
    if len(counts) != 1:
        return None, None
    return counts.pop(), (min(cpus) if cpus else None)


def memory_estimate(records, faces, gpus, nodes, gpu_memory, host_memory, headroom, activity):
    """Resident capacity, streaming and host memory of a fit (the precision decision's model)."""
    per_face = DEVICE_BYTES_PER_FACE[activity]
    usable = gpu_memory * (1.0 - headroom) - DEVICE_FIXED_BYTES - faces * per_face
    fp64_need = records * FP64_DEVICE_BYTES_PER_RECORD
    fp64_capacity = max(0.0, usable) * gpus * (1.0 - PARTITION_IMBALANCE)
    capacity = max(0.0, usable - STREAM_SLOT_BYTES) / FP32_DEVICE_BYTES_PER_RECORD * gpus * (1.0 - PARTITION_IMBALANCE)
    if capacity >= records:                    # all resident: no slots, no streaming
        streamed = 0
    else:
        streamed = int(records - capacity)
    ranks_per_node = max(1, gpus // max(1, nodes))
    pinned = streamed / max(1, nodes) * PINNED_BYTES_PER_STREAMED_RECORD
    host_need = pinned + ranks_per_node * HOST_WORKER_BYTES + faces * HOST_DRIVER_BYTES_PER_FACE
    return {"records": int(records), "faces": int(faces), "gpus": int(gpus), "nodes": int(nodes),
            "gpu_memory_bytes": int(gpu_memory), "headroom": headroom, "activity": activity,
            "usable_bytes_per_gpu": int(usable), "fp64_bytes": int(fp64_need),
            "fp64_capacity_bytes": int(fp64_capacity), "fp64_fits": usable > 0 and fp64_need <= fp64_capacity,
            "fp32_resident_records": int(min(records, capacity)), "fp32_streamed_records": int(streamed),
            "pinned_bytes_per_node": int(pinned), "host_bytes_rank0_node": int(host_need),
            "host_memory_bytes": host_memory,
            "host_fits": host_memory is None or host_need <= HOST_USABLE_FRACTION * host_memory,
            "gpu_fits": usable > STREAM_SLOT_BYTES or capacity >= records}


def gravity_runs(roots):
    """Recorded GPU gravity runs (c04 ``gpu_manifest.json``): policy, targets, body faces, rate per GPU."""
    import math
    runs = []
    for root in roots:
        for path in sorted(Path(root).glob("cuda_c04/equivalence_*/gpu_manifest.json")):
            try:
                m = read_json(path)
            except (OSError, ValueError):
                continue
            precision = m.get("evaluator", {}).get("precision")
            policy = precision.get("policy") if isinstance(precision, dict) else precision
            kind = {v: k for k, v in GRAVITY_POLICIES.items()}.get(policy)
            rate, targets, seconds = m.get("pairs_per_second"), m.get("targets"), m.get("seconds")
            if not (kind and rate and targets and seconds):
                continue
            devices = len(m.get("devices") or [0])
            device = (m.get("manifest") or {}).get("device") or {}
            runs.append({"precision": kind, "targets": int(targets), "seconds": float(seconds),
                         "pairs": float(rate) * float(seconds), "rate_per_gpu": float(rate) / devices,
                         "devices": devices, "gpu": device.get("name", "unknown GPU"), "path": str(path.parent)})
    for run in runs:
        run["log_targets"] = math.log(run["targets"])
    return runs


def gravity_time(runs, precision, targets, body_faces, gpus):
    """(seconds, the run it is scaled from) for ``targets`` x ``body_faces`` pairs on ``gpus`` GPUs, or None.

    Scaled from the recorded run of the same policy with the nearest target count: the mixed rate grows
    with the mesh (near-face share), so a run of similar size is the best guide.
    """
    import math
    same = [r for r in runs if r["precision"] == precision]
    if not same or not targets or not body_faces or not gpus:
        return None
    run = min(same, key=lambda r: abs(r["log_targets"] - math.log(targets)))
    return float(targets) * float(body_faces) / (run["rate_per_gpu"] * gpus), run


def duration_text(seconds):
    if seconds < 120:
        return "%.0f s" % seconds
    if seconds < 7200:
        return "%.0f min" % (seconds / 60.0)
    return "%.1f h" % (seconds / 3600.0)


# ----------------------------------------------------------------------------- request
def build_request(args, stored):
    """The normalized request: stored values, overridden by the options given now."""
    previous = (stored or {}).get("request", {})
    sid = (stored or {}).get("identity", {})

    def pick(name, default=None):
        value = getattr(args, name, None)
        if value is not None:
            return value
        return previous.get(name, default)
    shape = args.shape or (sid.get("shape") or {}).get("path")
    bc = args.surface_bc or sid.get("surface_bc")
    if not shape or not bc:
        raise PipelineError("a new restart directory needs --shape and --surface-bc")
    app = Path(pick("app_root", str(APP_ROOT))).resolve()
    reference = Path(pick("reference", os.environ.get("LIOUVILLE_COPS_REFERENCE")
                          or str(app.parent / "LiouvilleCOPS"))).resolve()
    venv_python = app / ".venv" / "bin" / "python"
    python = pick("python", str(venv_python) if venv_python.exists() else sys.executable)
    fits = pick("fits", ",".join(DEFAULT_FITS))
    fits = [f.strip() for f in (fits.split(",") if isinstance(fits, str) else fits) if f.strip()]
    fits = list(dict.fromkeys(fits))
    unknown = [f for f in fits if f not in FIT_TEMPLATES]
    if unknown or not fits:
        raise PipelineError("--fits takes a subset of %s (got %s)" % (",".join(FIT_TEMPLATES), fits))
    fit_args = pick("fit_args", DEFAULT_FIT_ARGS)
    fit_args = shlex.split(fit_args) if isinstance(fit_args, str) else list(fit_args)
    owned = [a for a in fit_args if a.split("=", 1)[0] in OWNED_FIT_ARGS]
    if owned:
        raise PipelineError("--fit-args must not contain %s (the pipeline sets them)" % ", ".join(owned))
    max_nfev = int(pick("max_nfev", DEFAULT_MAX_NFEV))
    if max_nfev < 1:
        raise PipelineError("--max-nfev must be >= 1")
    restart_from = dict(previous.get("restart_from", {}))
    if args.restart_from:
        restart_from = {}
        for item in args.restart_from:
            name, sep, path = item.partition("=")
            if sep and name in FIT_TEMPLATES:
                restart_from[name] = str(Path(path).resolve())
            elif len(fits) == 1:
                restart_from[fits[0]] = str(Path(item).resolve())
            else:
                raise PipelineError("--restart-from %s: with several fits give FIT=FILE (for example "
                                    "N30=/path/fit_restart.npz)" % item)
    stray = [f for f in restart_from if f not in fits]
    if stray:
        raise PipelineError("--restart-from names fits that are not requested: %s" % ", ".join(stray))
    fine = dict(FINE_TUNE_DEFAULTS, **previous.get("fine_tune_settings", {}))
    for key, name in (("minutes", "ft_minutes"), ("intervals", "ft_intervals"), ("method", "ft_method"),
                      ("max_nfev", "ft_max_nfev")):
        if getattr(args, name, None) is not None:
            fine[key] = getattr(args, name)
    ft_fits = args.ft_fits if args.ft_fits is not None else previous.get("ft_fits")
    ft_fits = [f for f in (ft_fits.split(",") if isinstance(ft_fits, str) else (ft_fits or fits)) if f]
    if [f for f in ft_fits if f not in fits]:
        raise PipelineError("--ft-fits must be a subset of --fits")
    wall = pick("wall_time")
    if wall not in (None, "none", "slurm"):
        wall_time.parse_duration(wall)                 # validates
    request = {
        "shape": str(Path(shape).resolve()),
        "shape_format": shape_format(shape, args.shape_format or (sid.get("shape") or {}).get("format")),
        "surface_bc": bc,
        "fit_input": str(Path(pick("fit_input", str(reference / "fit_base.in"))).resolve()),
        "gravity_shape": str(Path(pick("gravity_shape", shape)).resolve()),
        "allow_genus": int(pick("allow_genus", 0)),
        "precision": pick("precision", "auto"),
        "gravity_precision": pick("gravity_precision", DEFAULT_GRAVITY_PRECISION),
        "base": pick("base"), "store": pick("store"),
        "visibility_bundle": pick("visibility_bundle"), "bundle": pick("bundle"),
        "fits": fits, "max_nfev": max_nfev, "fit_args": fit_args, "tag": pick("tag", ""),
        "restart_from": restart_from, "translate_restart": bool(pick("translate_restart", False)),
        "fine_tune": bool(pick("fine_tune", False)), "fine_tune_settings": fine, "ft_fits": ft_fits,
        "activity": pick("activity"), "transport": pick("transport"),
        "headroom_fraction": pick("headroom_fraction"),
        "wall_time": None if wall == "none" else wall,
        "wall_margin": float(pick("wall_margin", wall_time.DEFAULT_MARGIN_SECONDS)),
        "save_seconds": float(pick("save_seconds", wall_time.DEFAULT_SAVE_SECONDS)),
        "output_root": pick("output_root", os.environ.get("COPS_OUTPUT_ROOT") or None),
        "app_root": str(app), "reference": str(reference), "python": python,
        "cpu_python": pick("cpu_python", python), "store_threads": int(pick("store_threads", 8)),
        "srun_args": pick("srun_args", SRUN_ARGS),
        "fit_nodes": previous.get("fit_nodes"),
        "slurm_gpus_per_node": previous.get("slurm_gpus_per_node"),
        "resubmit": int(pick("resubmit", 0)),
        "slurm_partition": pick("slurm_partition"), "slurm_account": pick("slurm_account"),
        "slurm_modules": list(args.slurm_module) if args.slurm_module else previous.get("slurm_modules", []),
        "slurm_times": dict(previous.get("slurm_times", SUBMIT_TIMES),
                            **{g: getattr(args, "time_" + g) for g in GROUPS if getattr(args, "time_" + g, None)}),
    }
    for key in ("base", "store", "visibility_bundle", "bundle"):
        if request[key]:
            request[key] = str(Path(request[key]).resolve())
    if request["store"] and request["precision"] == "fp64":
        raise PipelineError("--store holds FP32 records: it needs --precision mixed or auto")
    if request["resubmit"] and not request["wall_time"]:
        raise PipelineError("--resubmit needs --wall-time (slurm or HH:MM:SS): a job killed at its time limit "
                            "cannot resubmit itself")
    return request


def identity_of(request, stored_identity=None):
    """The fields that define the campaign's products (a different value needs a new restart dir)."""
    old = stored_identity or {}

    def ident(path, key, **kw):
        if not Path(path).exists():
            if (old.get(key) or {}).get("path") == str(Path(path).resolve()):
                return dict(old[key], verified=False)
            raise PipelineError("%s does not exist: %s" % (key.replace("_", " "), path))
        return file_identity(path, old.get(key), **kw)
    shape = ident(request["shape"], "shape", count_facets=request["shape_format"] == "nastran")
    shape["format"] = request["shape_format"]
    gravity = None
    if request["gravity_shape"] != request["shape"]:
        gravity = ident(request["gravity_shape"], "gravity_shape",
                        count_facets=Path(request["gravity_shape"]).suffix.lower() in (".bdf", ".nas"))
    return {"shape": shape, "surface_bc": request["surface_bc"],
            "fit_input": ident(request["fit_input"], "fit_input"), "gravity_shape": gravity,
            "allow_genus": request["allow_genus"], "precision": request["precision"],
            "gravity_precision": request["gravity_precision"],
            "external": {k: request[k] for k in ("base", "store", "visibility_bundle", "bundle")}}


def identity_differences(a, b):
    """Keys whose values differ (files compared by SHA-256 and format, not by path or mtime)."""
    def key(identity):
        out = {}
        for name, value in identity.items():
            if isinstance(value, dict) and "sha256" in value:
                out[name] = (value["sha256"], value.get("format"))
            else:
                out[name] = value
        return out
    ka, kb = key(a), key(b)
    return sorted(k for k in set(ka) | set(kb) if ka.get(k) != kb.get(k))


# ----------------------------------------------------------------------------- stages
class Stage(object):
    """One step: its output, the command that makes it and the inputs that define it."""

    def __init__(self, sid, group, output, tool=None, command=None, cwd=None, env=None, inputs=None,
                 needed=True, why="", external=False, aside=None, check=None, depends=(), pending=None,
                 code=(), upstream=("run_manifest",), fit=None):
        self.sid, self.group, self.output, self.tool = sid, group, output, tool
        self.command, self.cwd, self.env = command, cwd, env or {}
        self.inputs, self.needed, self.why, self.external = inputs or {}, needed, why, external
        self.aside = aside if aside is not None else output        # what is moved aside on a retry
        self.check = check or (lambda: (Path(output).exists(), "missing %s" % output))
        self.depends, self.pending, self.code = tuple(depends), pending, tuple(code)
        self.upstream = tuple(upstream)               # input keys naming earlier products: a change means redo
        self.fit = fit                                # fit stages: the fit name (segments)


class Pipeline(object):

    def __init__(self, request, stored, hardware):
        self.r = request
        self.stored = stored or {}
        self.hw = hardware
        self.dir = Path(request["restart_dir"])
        self.app = Path(request["app_root"])
        self.reference = Path(request["reference"])
        self.identity = None
        self.end = None                               # job end (epoch seconds) of --wall-time, or None

    # products -----------------------------------------------------------------------------
    def tool(self, name):
        return self.app / "tools" / name

    def base_dir(self):
        return Path(self.r["base"]) if self.r["base"] else self.dir / "base"

    def store_dir(self):
        return Path(self.r["store"]) if self.r["store"] else self.dir / "store"

    def visibility_bundle(self):
        if self.r["visibility_bundle"]:
            return Path(self.r["visibility_bundle"])
        return self.dir / "visibility" / "visibility_templates.npz"

    def bundle(self):
        return Path(self.r["bundle"]) if self.r["bundle"] else self.dir / "bundle" / "restart_templates.npz"

    def field_dir(self):
        found = sorted((self.dir / "gravity" / "cuda_c04").glob("equivalence_*"))
        return found[-1] if found else None

    def cache_file(self):
        return self.dir / "cache" / "admitted_gravity_cache.npz"

    def fit_name(self, fit):
        return "%s_n%d%s" % (fit, self.r["max_nfev"], ("_" + self.r["tag"]) if self.r["tag"] else "")

    def fit_dir(self, fit):
        return self.dir / "fits" / self.fit_name(fit)

    def fit_sid(self, fit):
        """Stage id of a fit: its directory name (fit, evaluation limit, tag), so other settings are a new fit."""
        return "fit:" + self.fit_name(fit)

    def segment_dir(self, fit, index):
        return self.fit_dir(fit) / ("seg%02d" % index)

    def fit_result_dir(self, fit):
        """The directory of the fit's final segment (planned: seg01)."""
        stamp = read_stamp(self.dir, self.fit_sid(fit))
        if stamp and stamp.get("segments"):
            return Path(stamp["segments"][-1]["dir"])
        return self.segment_dir(fit, 1)

    def finetune_dir(self, fit):
        f = self.r["fine_tune_settings"]
        return self.dir / "finetune" / ("%s_m%g_i%d_%s" % (self.fit_name(fit), f["minutes"], f["intervals"],
                                                            f["method"]))

    def translated(self, fit):
        return self.dir / "translate" / fit / "fit_restart.npz"

    def environment(self):
        env = dict(THREAD_ENVIRONMENT)
        env["LIOUVILLE_COPS_REFERENCE"] = str(self.reference)
        env["COPS_OUTPUT_ROOT"] = self.r["output_root"] or str(self.dir)
        return env

    # decisions ----------------------------------------------------------------------------
    def base_ready(self):
        try:
            return read_json(self.base_dir() / "rebuild.json").get("status") == "complete"
        except (OSError, ValueError):
            return False

    def counts(self):
        """(records, faces) of the base, or (None, None) before it exists."""
        if not self.base_ready():
            return None, None
        base = self.base_dir()
        try:
            records = read_json(self.store_dir() / "store_manifest.json")["records"]
        except (OSError, ValueError, KeyError):
            records = kernel_records(base)
        return int(records), elements_faces(base)

    def fit_transport(self, nodes):
        return self.r["transport"] or ("mpi" if nodes > 1 else "local")

    def fit_hardware(self):
        """(gpus, nodes) the fits run on: the submitted fit job's, else this allocation's."""
        nodes = self.r["fit_nodes"] or self.hw["nodes"]
        per_node = self.r["slurm_gpus_per_node"] or self.hw["gpus_per_node"]
        return (per_node * nodes if per_node else None), nodes

    def estimate(self):
        """The memory model for the fits on the detected hardware, or (None, reason)."""
        records, faces = self.counts()
        if records is None:
            return None, "known after the base stage"
        gpus, nodes = self.fit_hardware()
        if not gpus or not self.hw["gpu_memory_bytes"]:
            return None, "GPU count/memory unknown (run detects them; or give --gpus-per-node, --gpu-memory-gb)"
        transport = self.fit_transport(nodes)
        activity = self.r["activity"] or ("device" if transport == "mpi" else "host")
        headroom = DEFAULT_HEADROOM if self.r["headroom_fraction"] is None else float(self.r["headroom_fraction"])
        return memory_estimate(records, faces, gpus, nodes, self.hw["gpu_memory_bytes"],
                               self.hw["host_memory_bytes"], headroom, activity), ""

    def precision(self):
        """(fp64|mixed, reason) or (None, reason) while undecided."""
        requested = self.r["precision"]
        if requested in ("fp64", "mixed"):
            return requested, "requested"
        resolved = self.stored.get("resolved", {}).get("precision")
        if resolved:
            return resolved["value"], resolved["reason"]
        if self.r["store"]:
            return "mixed", "an FP32 record store is given (--store)"
        est, why = self.estimate()
        if est is None:
            return None, why
        if est["fp64_fits"]:
            return "fp64", "auto: FP64 records %s fit resident on %d GPUs (%s usable)" % (
                gb(est["fp64_bytes"]), est["gpus"], gb(est["fp64_capacity_bytes"]))
        return "mixed", "auto: FP64 records %s exceed %s on %d GPUs; FP32 records, %d resident, %d streamed" % (
            gb(est["fp64_bytes"]), gb(est["fp64_capacity_bytes"]), est["gpus"], est["fp32_resident_records"],
            est["fp32_streamed_records"])

    def memory_refusal(self, precision):
        """None, or why the fits cannot run on the detected hardware with ``precision``."""
        est, _ = self.estimate()
        if est is None or precision is None:
            return None
        if precision == "fp64" and not est["fp64_fits"]:
            return ("FP64 records (%s) do not fit resident on %d GPUs of %s (%s usable with the headroom); "
                    "FP64 has no streaming path: use --precision mixed or more GPUs"
                    % (gb(est["fp64_bytes"]), est["gpus"], gb(est["gpu_memory_bytes"]),
                       gb(est["fp64_capacity_bytes"])))
        if precision == "mixed" and not est["gpu_fits"]:
            return ("the per-face arrays and the stream slots alone exceed the GPU memory (%s per GPU usable); "
                    "a smaller --headroom-fraction or larger GPUs are needed" % gb(est["usable_bytes_per_gpu"]))
        if precision == "mixed" and not est["host_fits"]:
            return ("the host memory of a node (%s) cannot hold the fit: %d records streamed need %s of pinned "
                    "buffers per node, plus the GPU workers and the driver: %s > %.0f%% of the node. Use more "
                    "nodes (--nodes) or GPUs with more memory"
                    % (gb(est["host_memory_bytes"]), est["fp32_streamed_records"], gb(est["pinned_bytes_per_node"]),
                       gb(est["host_bytes_rank0_node"]), 100 * HOST_USABLE_FRACTION))
        return None

    def needs(self):
        sources = {src for fit in self.r["fits"] for src, _ in FIT_TEMPLATES[fit]}
        return {"visibility": "visibility" in sources and not self.r["visibility_bundle"],
                "bundle": "bundle" in sources and not self.r["bundle"]}

    def gpu_list(self):
        return ",".join(str(k) for k in range(max(1, self.hw["gpus_per_node"] or 1)))

    # stage list ---------------------------------------------------------------------------
    def stages(self):
        r, py = self.r, self.r["python"]
        precision, _ = self.precision()
        need = self.needs()
        base = self.base_dir()
        out = []
        if r["base"]:
            out.append(Stage("base", "prepare", base, external=True, why="--base", check=self.check_base))
        else:
            cmd = [r["cpu_python"], str(self.reference / "tools" / "rebuild_cops_boundary.py"),
                   "--input", r["fit_input"], "--surface-bc", r["surface_bc"], "--shape-model", r["shape"],
                   "--shape-model-format", r["shape_format"], "--out", str(base),
                   "--threads", str(self.hw["threads"]),
                   "--executable", str(self.app / "native" / "liouville_cops_cuda")]
            out.append(Stage("base", "prepare", base, tool=self.reference / "tools" / "rebuild_cops_boundary.py",
                             command=cmd, cwd=self.reference, env=dict(self.environment(), COPS_CUDA_DEVICE="0"),
                             inputs={"shape_sha256": self.identity["shape"]["sha256"],
                                     "surface_bc": r["surface_bc"],
                                     "fit_input_sha256": self.identity["fit_input"]["sha256"]},
                             check=self.check_base, code=(self.app / "native" / "liouville_cops_cuda",),
                             upstream=()))
        base_inputs = {"base": str(base), "run_manifest": self.run_manifest_digest()}
        store = self.store_dir()
        if r["store"]:
            out.append(Stage("store", "prepare", store, external=True, why="--store", check=self.check_store,
                             depends=("base",)))
        else:
            cmd = [py, "-W", "ignore", str(self.tool("c15_record_store.py")), "--run", str(base),
                   "--out", str(store), "--threads", str(r["store_threads"])]
            out.append(Stage("store", "prepare", store, tool=self.tool("c15_record_store.py"), command=cmd,
                             env=self.environment(), inputs=base_inputs, needed=precision != "fp64",
                             why=("fp64 records need no store" if precision == "fp64" else
                                  "" if precision else "only with mixed precision (decided after base)"),
                             check=self.check_store, depends=("base",),
                             pending=None if precision else "precision decision"))
        vis = self.visibility_bundle()
        if r["visibility_bundle"]:
            out.append(Stage("visibility", "prepare", vis, external=True, why="--visibility-bundle",
                             check=lambda: self.check_bundle(vis, ("late_ng_not_rg",)), depends=("base",)))
        else:
            cmd = [py, "-W", "ignore", str(self.tool("c05_write_bundle.py")), "--run", str(base),
                   "--visibility-only", "--out", str(vis), "--device", "0"]
            out.append(Stage("visibility", "prepare", vis, tool=self.tool("c05_write_bundle.py"), command=cmd,
                             env=self.environment(), inputs=base_inputs, needed=need["visibility"],
                             why="" if need["visibility"] else "no requested fit uses it", aside=vis.parent,
                             check=lambda: self.check_bundle(vis, ("late_ng_not_rg",)), depends=("base",)))
        gravity_inputs = {"gravity_shape_sha256": (self.identity["gravity_shape"] or self.identity["shape"])["sha256"],
                          "allow_genus": r["allow_genus"]}
        bundle = self.bundle()
        if r["bundle"]:
            for sid in ("gravity", "cache"):
                out.append(Stage(sid, "prepare", self.dir / sid, needed=False, why="--bundle given"))
            out.append(Stage("bundle", "prepare", bundle, external=True, why="--bundle",
                             check=lambda: self.check_bundle(bundle, ("late_ng_not_rg", "cliff_ge_60deg")),
                             depends=("base",)))
        else:
            why = "" if need["bundle"] else "no requested fit uses a cliff template"
            gpre = (["--precision", "mixed", "--near-factor", "%g" % GRAVITY_NEAR_FACTOR]
                    if r["gravity_precision"] == "mixed" else ["--precision", "fp64"])
            genus = ["--allow-genus", str(r["allow_genus"])] if r["allow_genus"] else []
            cmd = ([py, "-W", "ignore", str(self.tool("c04_gravity_equivalence.py")), "gpu", "--run", str(base),
                    "--gravity-shape", r["gravity_shape"], "--devices", self.gpu_list(),
                    "--output-root", str(self.dir / "gravity")] + gpre + genus)
            out.append(Stage("gravity", "prepare", self.dir / "gravity", tool=self.tool("c04_gravity_equivalence.py"),
                             command=cmd, env=self.environment(),
                             inputs=dict(base_inputs, gravity_precision=r["gravity_precision"], **gravity_inputs),
                             needed=need["bundle"], why=why, check=self.check_gravity, depends=("base",)))
            field = self.field_dir()
            cache = self.cache_file()
            cmd = ([py, "-W", "ignore", str(self.tool("c04_write_gpu_cache.py")),
                    "--field-dir", str(field) if field else "<gravity/cuda_c04/equivalence_UTC>",
                    "--run", str(base), "--gravity-shape", r["gravity_shape"]] + genus +
                   ["--out", str(cache), "--device", "0"])
            out.append(Stage("cache", "prepare", cache, tool=self.tool("c04_write_gpu_cache.py"), command=cmd,
                             env=self.environment(),
                             inputs=dict(base_inputs, field=str(field) if field else None, **gravity_inputs),
                             needed=need["bundle"], why=why, aside=cache.parent, depends=("gravity",),
                             pending=None if field else "gravity field", upstream=("run_manifest", "field")))
            cmd = ([py, "-W", "ignore", str(self.tool("c05_write_bundle.py")), "--run", str(base),
                    "--gravity-cache", str(cache), "--gravity-shape", r["gravity_shape"]] + genus +
                   ["--out", str(bundle), "--device", "0"])
            out.append(Stage("bundle", "prepare", bundle, tool=self.tool("c05_write_bundle.py"), command=cmd,
                             env=self.environment(), inputs=dict(base_inputs, cache=str(cache), **gravity_inputs),
                             needed=need["bundle"], why=why, aside=bundle.parent,
                             check=lambda: self.check_bundle(bundle, ("late_ng_not_rg", "cliff_ge_60deg")),
                             depends=("cache",)))
        for fit in r["fits"]:
            source = r["restart_from"].get(fit)
            if not (source and r["translate_restart"]):
                continue
            target = self.template_bundle_for(fit)
            cmd = ([py, "-W", "ignore", str(self.tool("restart_tool.py")), "translate", source, "--to-run", str(base),
                    "--save", str(self.translated(fit))] + (["--target-bundle", str(target)] if target else []) +
                   ["--device", "0", "--precision", precision or "<after base>"])
            deps = ("base",) + (("store",) if precision == "mixed" else ()) + self.bundle_stages(fit)
            out.append(Stage("translate:" + fit, "prepare", self.translated(fit), tool=self.tool("restart_tool.py"),
                             command=cmd, env=self.fit_environment(precision),
                             inputs=dict(base_inputs, restart=quick_identity(source) if Path(source).exists() else source,
                                         target_bundle=str(target) if target else None),
                             aside=self.translated(fit).parent, depends=deps,
                             pending=None if precision else "precision decision"))
        for fit in r["fits"]:
            out.append(self.fit_stage(fit, precision, base_inputs))
        for fit in r["fits"]:
            out.append(self.finetune_stage(fit, precision))
        return out

    def bundle_stages(self, fit):
        return tuple(dict.fromkeys(src for src, _ in FIT_TEMPLATES[fit]))

    def template_bundle_for(self, fit):
        sources = self.bundle_stages(fit)
        if "bundle" in sources:
            return self.bundle()
        if "visibility" in sources:
            return self.visibility_bundle()
        return None

    def template_args(self, fit):
        args = []
        for source, column in FIT_TEMPLATES[fit]:
            path = self.visibility_bundle() if source == "visibility" else self.bundle()
            args += ["--activity-template", str(path), column]
        return args

    def fit_environment(self, precision):
        env = self.environment()
        if precision == "mixed":
            env["COPS_CUDA_RECORD_STORE"] = str(self.store_dir())
        return env

    def user_restart(self, fit):
        source = self.r["restart_from"].get(fit)
        if source and self.r["translate_restart"]:
            return str(self.translated(fit))
        return source

    def fit_command(self, fit, precision, out, restart, max_nfev, segment):
        r, hw = self.r, self.hw
        transport = self.fit_transport(hw["nodes"])
        driver = [r["python"], "-W", "ignore", str(self.tool("gpu_fit.py"))]
        if transport == "mpi":
            driver += ["--transport", "mpi"]
        else:
            driver += ["--devices", self.gpu_list()]
        if r["activity"]:
            driver += ["--activity", r["activity"]]
        driver += ["--precision", precision or "<after base>"]
        if precision == "mixed":
            driver += ["--record-store", str(self.store_dir())]
        if r["headroom_fraction"] is not None:
            driver += ["--headroom-fraction", repr(float(r["headroom_fraction"]))]
        driver += ["--max-nfev", str(max_nfev), "--segment", str(segment)]
        if restart:
            driver += ["--restart-from", str(restart)]
        if self.end is not None:
            driver += ["--deadline", "%.0f" % self.end, "--wall-margin", "%g" % r["wall_margin"],
                       "--save-seconds", "%g" % r["save_seconds"]]
        driver += ["--", str(self.base_dir())] + list(r["fit_args"]) + ["--out", str(out)] + self.template_args(fit)
        if transport != "mpi":
            return driver
        if hw["launcher"] == "srun":
            return ["srun"] + shlex.split(r["srun_args"]) + ["--ntasks=%d" % hw["ranks"]] + driver
        exports = []
        for name in sorted(self.fit_environment(precision)):
            exports += ["-x", name]
        return ["mpirun", "-n", str(hw["ranks"]), "--bind-to", "none"] + exports + driver

    def next_segment(self, fit):
        """(index, restart, remaining nfev, clean segments) of the fit's next segment."""
        stamp = read_stamp(self.dir, self.fit_sid(fit))
        clean = []
        if stamp and stamp.get("status") != "complete":
            clean = [s for s in stamp.get("segments", [])
                     if s["status"] == "interrupted" and Path(s["restart"]).is_file()]
        used = sum(int(s["nfev"]) for s in clean)
        restart = clean[-1]["restart"] if clean else self.user_restart(fit)
        return len(clean) + 1, restart, self.r["max_nfev"] - used, clean

    def fit_stage(self, fit, precision, base_inputs):
        index, restart, remaining, _ = self.next_segment(fit)
        out = self.segment_dir(fit, index)
        cmd = self.fit_command(fit, precision, out, restart, max(1, remaining), index)
        deps = (("base",) + (("store",) if precision == "mixed" else ()) + self.bundle_stages(fit) +
                (("translate:" + fit,) if self.r["restart_from"].get(fit) and self.r["translate_restart"] else ()))
        source = self.r["restart_from"].get(fit)
        inputs = dict(base_inputs, precision=precision, max_nfev=self.r["max_nfev"], fit_args=self.r["fit_args"],
                      templates=[[str(self.visibility_bundle() if s == "visibility" else self.bundle()), c]
                                 for s, c in FIT_TEMPLATES[fit]],
                      restart_from=(quick_identity(source) if source and Path(source).exists() else source),
                      translated=bool(source and self.r["translate_restart"]), activity=self.r["activity"])
        final = self.fit_dir(fit)
        return Stage(self.fit_sid(fit), "fit", final, tool=self.tool("gpu_fit.py"), command=cmd,
                     env=self.fit_environment(precision), inputs=inputs, depends=deps,
                     check=lambda: self.check_fit(self.fit_result_dir(fit)),
                     pending=None if precision else "precision decision",
                     code=(self.reference / "tools" / "fit_cops.py",), fit=fit)

    def finetune_stage(self, fit, precision):
        r, f = self.r, self.r["fine_tune_settings"]
        wanted = r["fine_tune"] and fit in r["ft_fits"]
        restart = self.fit_result_dir(fit) / "fit_restart.npz"
        cmd = [r["python"], "-W", "ignore", str(self.tool("fine_tune_gpu.py")), str(self.base_dir()),
               "--restart-from", str(restart), "--minutes", "%g" % f["minutes"], "--intervals", str(f["intervals"]),
               "--method", f["method"], "--max-nfev", str(f["max_nfev"]), "--devices", self.gpu_list(),
               "--precision", precision or "<after base>"]
        if r["headroom_fraction"] is not None:
            cmd += ["--headroom-fraction", repr(float(r["headroom_fraction"]))]
        out = self.finetune_dir(fit)
        cmd += ["--out", str(out)]
        inputs = {"restart": quick_identity(restart) if restart.exists() else str(restart),
                  "settings": f, "precision": precision}
        return Stage("finetune:" + self.fit_name(fit), "finetune", out, tool=self.tool("fine_tune_gpu.py"), command=cmd,
                     env=self.fit_environment(precision), inputs=inputs, needed=wanted,
                     why="" if wanted else ("--no-fine-tune" if not r["fine_tune"] else "not in --ft-fits"),
                     depends=(self.fit_sid(fit),), check=lambda: ((out / "fine_tuning.json").is_file(),
                                                             "no fine_tuning.json in %s" % out),
                     pending=None if precision else "precision decision", upstream=("restart",))

    def run_manifest_digest(self):
        try:
            return sha256_file(self.base_dir() / "run_manifest.json")
        except OSError:
            return None

    # output checks (cheap: manifests and headers, never the kernel files) ----------------------
    def check_base(self):
        base = self.base_dir()
        try:
            manifest = read_json(base / "rebuild.json")
        except (OSError, ValueError):
            return False, "no readable rebuild.json in %s" % base
        if manifest.get("status") != "complete":
            return False, "rebuild.json status %r" % manifest.get("status")
        if manifest.get("surface_bc") != self.r["surface_bc"]:
            return False, "base was built with SurfaceBC=%s, not %s" % (manifest.get("surface_bc"),
                                                                         self.r["surface_bc"])
        recorded = manifest.get("shape_model")
        if recorded and str(Path(recorded).resolve()) != self.r["shape"]:
            if Path(recorded).is_file():
                if sha256_file(recorded) != self.identity["shape"]["sha256"]:
                    return False, "base mesh %s differs from --shape %s" % (recorded, self.r["shape"])
            elif Path(recorded).name != Path(self.r["shape"]).name:
                return False, "base mesh %s differs from --shape %s" % (recorded, self.r["shape"])
        try:
            run_manifest = read_json(base / "run_manifest.json")
        except (OSError, ValueError):
            return False, "no readable run_manifest.json in %s" % base
        states = run_manifest.get("states", {})
        if states.get("process") != "complete" or states.get("cache") != "verified":
            return False, "run_manifest.json states %s (need process complete, cache verified)" % states
        recorded_digest = (manifest.get("run_manifest") or {}).get("sha256")
        if recorded_digest and recorded_digest != self.run_manifest_digest():
            return False, "run_manifest.json differs from the one rebuild.json recorded"
        try:
            n = elements_faces(base)
        except (OSError, PipelineError) as exc:
            return False, str(exc)
        facets = self.identity["shape"].get("facets")
        note = ("; the mesh file has %d CTRIA3 cards (the solver's reader keeps %d)" % (facets, n)
                if facets and facets != n else "")
        return True, "complete, SurfaceBC=%s, %d facets%s" % (manifest.get("surface_bc"), n, note)

    def check_store(self):
        store = self.store_dir()
        try:
            manifest = read_json(store / "store_manifest.json")
        except (OSError, ValueError):
            return False, "no readable store_manifest.json in %s" % store
        digest = self.run_manifest_digest()
        recorded = (manifest.get("run_manifest") or {}).get("sha256")
        if digest is None or recorded != digest:
            return False, "store %s belongs to another base (run_manifest.json digest)" % store
        if manifest.get("record_precision") != "fp32":
            return False, "store records are %r, not fp32" % manifest.get("record_precision")
        try:
            if manifest.get("nfaces") not in (None, elements_faces(self.base_dir())):
                return False, "store has %s faces, the base %d" % (manifest.get("nfaces"),
                                                                   elements_faces(self.base_dir()))
        except (OSError, PipelineError) as exc:
            return False, str(exc)
        return True, "%s records, mesh %s" % (manifest.get("records"), str(manifest.get("mesh_identity"))[:12])

    def mesh_identity(self):
        """The base's mesh identity from a product that records it (store, then bundles), else None."""
        candidates = [self.store_dir() / "store_manifest.json"]
        candidates += [Path(path).with_suffix(".json") for path in (self.visibility_bundle(), self.bundle())]
        for path in candidates:
            try:
                meta = read_json(path)
            except (OSError, ValueError):
                continue
            if path.name == "store_manifest.json":
                if (meta.get("run_manifest") or {}).get("sha256") != self.run_manifest_digest():
                    continue
            elif str(Path(meta.get("base_run", "")).resolve()) != str(self.base_dir().resolve()):
                continue
            if meta.get("mesh_identity"):
                return meta["mesh_identity"]
        return None

    def check_bundle(self, path, columns):
        meta_path = Path(path).with_suffix(".json")
        try:
            meta = read_json(meta_path)
        except (OSError, ValueError):
            return False, "no readable %s" % meta_path.name
        if not Path(path).is_file():
            return False, "missing %s" % path
        missing = [c for c in columns if c not in meta.get("names", [])]
        if missing:
            return False, "%s lacks %s" % (Path(path).name, ", ".join(missing))
        try:
            store = read_json(self.store_dir() / "store_manifest.json")
            mesh = store.get("mesh_identity") if (store.get("run_manifest") or {}).get(
                "sha256") == self.run_manifest_digest() else None
        except (OSError, ValueError):
            mesh = None
        if mesh and meta.get("mesh_identity") != mesh:
            return False, "%s was made on another mesh" % Path(path).name
        if not mesh and str(Path(meta.get("base_run", "")).resolve()) != str(self.base_dir().resolve()):
            return False, "%s was made from base %s" % (Path(path).name, meta.get("base_run"))
        return True, "columns %s" % ", ".join(meta.get("names", []))

    def check_gravity(self):
        field = self.field_dir()
        if field is None or not (field / "gpu_manifest.json").is_file():
            return False, "no gravity field (gpu_manifest.json) under %s" % (self.dir / "gravity")
        return True, str(field)

    def check_fit(self, out):
        try:
            state = read_json(out / "fit_manifest.json").get("states", {}).get("process")
        except (OSError, ValueError):
            return False, "no readable fit_manifest.json in %s" % out
        if state != "complete":
            return False, "fit process state %r" % state
        if not (out / "fit_restart.npz").is_file():
            return False, "no fit_restart.npz in %s" % out
        return True, final_cost(out)

    # restart files given by the user -----------------------------------------------------------
    def check_restarts(self):
        """Template and mesh checks of --restart-from files (notes; refusals raise)."""
        notes = []
        for fit, path in sorted(self.r["restart_from"].items()):
            if not Path(path).is_file():
                raise PipelineError("--restart-from %s: no such file %s" % (fit, path))
            meta = restart_meta(path)
            allowed = [c for _, c in FIT_TEMPLATES[fit]]
            extra = [t for t in meta["templates"] if t not in allowed and t.replace("xl_", "", 1) not in allowed]
            if extra:
                raise PipelineError("--restart-from %s: the restart has templates %s that %s does not fit "
                                    "(gpu_fit would add them to the fit)" % (fit, extra, fit))
            mesh = self.mesh_identity()
            if mesh is None:
                notes.append("%s: restart mesh %s; checked against the base when the fit starts (gpu_fit)"
                             % (fit, meta["mesh_identity"][:12]))
            elif meta["mesh_identity"] != mesh:
                if not self.r["translate_restart"]:
                    raise PipelineError("--restart-from %s: the restart belongs to another mesh (%s, %s facets) "
                                        "than the base (%s). Translating it (restart_tool.py translate: nearest "
                                        "facet centroid, amplitudes rescaled) is a choice; give "
                                        "--translate-restart to make it" % (fit, meta["mesh_identity"][:12],
                                                                           meta["faces"], mesh[:12]))
                notes.append("%s: restart of another mesh, translated (translate:%s)" % (fit, fit))
            else:
                notes.append("%s: restart on this mesh" % fit)
                if self.r["translate_restart"]:
                    raise PipelineError("--translate-restart: the %s restart is already on this mesh; "
                                        "drop the option" % fit)
        return notes


def final_cost(out):
    try:
        for line in (Path(out) / "fit_result.txt").read_text().splitlines():
            if line.startswith("final_cost"):
                return line.strip()
    except OSError:
        pass
    return "complete"


# ----------------------------------------------------------------------------- state of stages
def stamp_path(directory, sid):
    return Path(directory) / "stamps" / (sid.replace(":", "_") + ".json")


def read_stamp(directory, sid):
    try:
        return read_json(stamp_path(directory, sid))
    except (OSError, ValueError):
        return None


def stage_state(pipe, stage):
    """(state, message): external, complete, incomplete, pending, blocked, conflict or not needed.

    A completed stage whose inputs changed only in ``stage.upstream`` keys (an earlier product was made
    again) is stale and runs again; other changed inputs are refused (prepare) or a conflict (fits).
    A fit stopped at a deadline is incomplete and continues with its next segment.
    """
    if not stage.needed:
        return "not needed", stage.why
    if stage.external:
        ok, message = stage.check()
        if not ok:
            raise PipelineError("%s (%s): %s" % (stage.sid, stage.why, message))
        return "external", message
    stamp = read_stamp(pipe.dir, stage.sid)
    if stamp is not None:
        old, new = stamp.get("inputs", {}), json.loads(json.dumps(stage.inputs))
        differ = sorted(k for k in set(old) | set(new) if old.get(k) != new.get(k))
        if differ and all(k in stage.upstream for k in differ):
            return "incomplete", "stale: %s changed since it ran" % ", ".join(differ)
        if differ:
            message = ("started with other settings (%s); use a new --restart-dir, or --tag for other fit "
                       "settings" % ", ".join(k for k in differ if k not in stage.upstream))
            if stage.group == "prepare":
                raise PipelineError("stage %s was %s" % (stage.sid, message))
            return "conflict", message
        if stamp.get("status") == "complete":
            ok, message = stage.check()
            return ("complete", message) if ok else ("incomplete", "stamp complete but %s" % message)
        if stamp.get("status") == "segment":
            segment = stamp["segments"][-1]
            return "incomplete", "segment %d stopped at the deadline (%d of %d evaluations used)" % (
                segment["index"], sum(s["nfev"] for s in stamp["segments"]), pipe.r["max_nfev"])
        return "incomplete", "previous attempt %s (%s)" % (stamp.get("status"), stamp.get("finished_utc") or
                                                           stamp.get("started_utc"))
    if os.path.lexists(str(stage.aside)):
        return "incomplete", "output exists without a stamp: %s" % stage.aside
    if stage.pending:
        return "blocked", "waits for the %s" % stage.pending
    return "pending", ""


def resumable(pipe, stage, state):
    """For ``status --check``: can another job continue this stage (versus a failure to look at)?"""
    if state in ("pending", "blocked"):
        return True
    if state != "incomplete":
        return False
    stamp = read_stamp(pipe.dir, stage.sid)
    return stamp is None or stamp.get("status") in ("segment", "deferred", "interrupted", "running")


# ----------------------------------------------------------------------------- execution
def environment_text(env):
    return " ".join("%s=%s" % (k, shlex.quote(v)) for k, v in sorted(env.items()))


def stage_cwd(stage):
    return stage.cwd or (stage.tool.parents[1] if stage.tool else Path.cwd())


def command_text(stage):
    return "(cd %s && env %s %s)" % (shlex.quote(str(stage_cwd(stage))), environment_text(stage.env),
                                     " ".join(shlex.quote(c) for c in stage.command))


@contextlib.contextmanager
def terminate_on_signal():
    """Turn SIGTERM (Slurm time limit, scancel) into ``Interrupted`` inside the block."""
    def handler(signum, frame):
        raise Interrupted("signal %d" % signum)
    previous = signal.signal(signal.SIGTERM, handler)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


def run_logged(command, cwd, env, log_path, prefix):
    """Run one stage command; its output goes to the log (timestamped) and to stdout."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("=== %s %s\n=== cwd %s\n=== %s\n" % (utc_now(), environment_text(env), cwd,
                                                       " ".join(shlex.quote(c) for c in command)))
        log.flush()
        process = subprocess.Popen(command, cwd=str(cwd), env=dict(os.environ, **env), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, universal_newlines=True, bufsize=1)
        try:
            for line in process.stdout:
                stamped = "[%s] %s" % (time.strftime("%H:%M:%S"), line)
                log.write(stamped)
                log.flush()
                sys.stdout.write("%s %s" % (prefix, stamped))
                sys.stdout.flush()
            return process.wait()
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            log.write("=== %s interrupted; process ended with %s\n" % (utc_now(), process.returncode))
            raise
        finally:
            process.stdout.close()


def code_identity(pipe, stage, tree_cache):
    from cops_cuda import provenance
    if "tree" not in tree_cache:
        tree = provenance.tree_identity(pipe.app / "cops_cuda")
        tree_cache["tree"] = {"root": tree["root"], "file_count": tree["file_count"],
                              "tree_sha256": tree["tree_sha256"]}
    code = {"runner": provenance.file_identity(Path(__file__)), "cops_cuda": tree_cache["tree"]}
    for path in ((stage.tool,) if stage.tool else ()) + stage.code:
        if Path(path).is_file():
            code[Path(path).name] = provenance.file_identity(path)
    return code


def recorded_seconds(pipe, stage):
    """Longest recorded duration of earlier attempts of a stage (killed ones are lower bounds), or None."""
    values = []
    for path in (pipe.dir / "stamps").glob(stage.sid.replace(":", "_") + ".json*"):
        try:
            seconds = read_json(path).get("seconds")
        except (OSError, ValueError):
            continue
        if seconds:
            values.append(float(seconds))
    return max(values) if values else None


def time_left(pipe):
    return None if pipe.end is None else pipe.end - pipe.r["wall_margin"] - time.time()


def can_start(pipe, stage):
    """(True, "") or (False, why) under --wall-time."""
    left = time_left(pipe)
    if left is None:
        return True, ""
    if stage.fit:
        need = pipe.r["save_seconds"]
        if left < need:
            return False, "%.0f s left, a fit segment needs more than %.0f s to save a result" % (left, need)
        return True, ""
    need = recorded_seconds(pipe, stage)
    if need is not None and need > left:
        return False, "%.0f s left; an earlier attempt took %.0f s" % (left, need)
    return True, ""


def begin_stamp(pipe, stage, tree_cache, moved, extra=None):
    stamp = {"schema": STAMP_SCHEMA, "schema_version": SCHEMA_VERSION, "stage": stage.sid, "group": stage.group,
             "status": "running", "inputs": stage.inputs, "command": stage.command, "cwd": str(stage_cwd(stage)),
             "environment": stage.env, "code": code_identity(pipe, stage, tree_cache),
             "started_utc": utc_now(), "host": socket.gethostname(), "pid": os.getpid(),
             "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "moved_aside": moved, "hardware": pipe.hw,
             "deadline_utc": wall_time.utc(pipe.end)}
    stamp.update(extra or {})
    write_json_atomic(stamp_path(pipe.dir, stage.sid), stamp)
    return stamp


def run_stage_command(pipe, stage, stamp):
    t0 = time.time()
    try:
        code = run_logged(stage.command, stage_cwd(stage), stage.env,
                          pipe.dir / "logs" / (stage.sid.replace(":", "_") + ".log"), "[%s]" % stage.sid)
    except BaseException as exc:
        stamp.update(status="interrupted", finished_utc=utc_now(), error=str(exc) or exc.__class__.__name__,
                     seconds=time.time() - t0)
        write_json_atomic(stamp_path(pipe.dir, stage.sid), stamp)
        raise
    return code, time.time() - t0


def execute(pipe, stage, tree_cache):
    """Run one stage: move a previous attempt aside, stamp, run, check, stamp. Returns
    ``complete``, ``failed`` or ``deferred`` (a fit segment stopped at the deadline)."""
    if stage.fit:
        return execute_fit(pipe, stage, tree_cache)
    stamp_file = stamp_path(pipe.dir, stage.sid)
    moved = []
    if os.path.lexists(str(stage.aside)):
        moved.append(str(move_aside(stage.aside)))
    if stamp_file.exists():
        moved.append(str(move_aside(stamp_file)))
    if stage.aside != stage.output:                   # a stage directory that holds one output file
        Path(stage.aside).mkdir(parents=True)
    else:
        Path(stage.output).parent.mkdir(parents=True, exist_ok=True)
    stamp = begin_stamp(pipe, stage, tree_cache, moved)
    print("[%s] %s: start%s" % (clock_text(), stage.sid,
                                (" (previous attempt moved to %s)" % ", ".join(moved)) if moved else ""), flush=True)
    code, seconds = run_stage_command(pipe, stage, stamp)
    ok, message = (stage.check() if code == 0 else (False, "exit status %d" % code))
    stamp.update(status="complete" if ok else "failed", exit_status=code, finished_utc=utc_now(),
                 seconds=seconds, result=message)
    if ok and Path(stage.output).exists():
        stamp["outputs"] = quick_identity(stage.output)
    write_json_atomic(stamp_file, stamp)
    print("[%s] %s: %s (%s, %.0f s)" % (clock_text(), stage.sid, stamp["status"], message, seconds), flush=True)
    return "complete" if ok else "failed"


def segment_of(out):
    """The segment record gpu_fit wrote into the fit manifest (compute_backend.wall_time.segments)."""
    try:
        manifest = read_json(Path(out) / "fit_manifest.json")
    except (OSError, ValueError):
        return None
    segments = (manifest.get("compute_backend", {}).get("wall_time") or {}).get("segments") or []
    return segments[-1] if segments else None


def execute_fit(pipe, stage, tree_cache):
    """One segment of a fit: the next after the last clean stop, else the first."""
    fit = stage.fit
    stamp_file = stamp_path(pipe.dir, stage.sid)
    old = read_stamp(pipe.dir, stage.sid)
    moved = []
    index, restart, remaining, clean = pipe.next_segment(fit)
    if old is None or old.get("status") == "complete" or (old.get("inputs") != json.loads(json.dumps(stage.inputs))):
        clean = []                                    # a new fit (or a stale one): start from its beginning
        index, restart, remaining = 1, pipe.user_restart(fit), pipe.r["max_nfev"]
        if os.path.lexists(str(pipe.fit_dir(fit))):
            moved.append(str(move_aside(pipe.fit_dir(fit))))
        if stamp_file.exists():
            moved.append(str(move_aside(stamp_file)))
    out = pipe.segment_dir(fit, index)
    if os.path.lexists(str(out)):
        moved.append(str(move_aside(out)))            # a killed or failed attempt of this segment
    out.parent.mkdir(parents=True, exist_ok=True)
    refusal = pipe.memory_refusal(pipe.precision()[0])
    if refusal:
        print("[%s] %s: refused: %s" % (clock_text(), stage.sid, refusal), flush=True)
        stamp = begin_stamp(pipe, stage, tree_cache, moved, {"segments": clean})
        stamp.update(status="failed", finished_utc=utc_now(), result="refused: " + refusal)
        write_json_atomic(stamp_file, stamp)
        return "failed"
    stage.command = pipe.fit_command(fit, pipe.precision()[0], out, restart, remaining, index)
    attempts = (old or {}).get("attempts", []) if clean else []
    stamp = begin_stamp(pipe, stage, tree_cache, moved, {"segments": clean, "attempts": attempts,
                                                        "estimate": pipe.estimate()[0]})
    print("[%s] %s: segment %d, %d evaluations left%s" % (clock_text(), stage.sid, index, remaining,
                                                          (", from %s" % restart) if restart else ""), flush=True)
    code, seconds = run_stage_command(pipe, stage, stamp)
    record = segment_of(out)
    attempt = {"index": index, "dir": str(out), "exit_status": code, "seconds": seconds,
               "job": os.environ.get("SLURM_JOB_ID"), "finished_utc": utc_now()}
    ok, message = pipe.check_fit(out) if code == 0 else (False, "exit status %d" % code)
    if code == 0 and ok and record and record["status"] == "complete":
        segment = dict(record, dir=str(out), restart=str(out / "fit_restart.npz"), seconds=seconds)
        stamp.update(status="complete", segments=clean + [segment], result=message)
        outcome = "complete"
    elif (code == wall_time.EXIT_WALL_TIME and record and record["status"] == "interrupted"
          and (out / "fit_restart.npz").is_file()):
        segment = dict(record, dir=str(out), restart=str(out / "fit_restart.npz"), seconds=seconds)
        stamp.update(status="segment", segments=clean + [segment],
                     result="stopped at the deadline: %d evaluations, cost %.6e" % (record["nfev"], record["cost"]))
        outcome = "deferred"
    elif code == wall_time.EXIT_WALL_TIME:
        stamp.update(status="deferred", result="no evaluation before the deadline (nothing saved)")
        outcome = "deferred"
    else:
        stamp.update(status="failed", result=message)
        outcome = "failed"
    attempt["status"] = stamp["status"]
    stamp.update(attempts=attempts + [attempt], exit_status=code, finished_utc=utc_now(), seconds=seconds)
    write_json_atomic(stamp_file, stamp)
    record_fit(pipe, fit, stamp)
    print("[%s] %s: %s (%s, %.0f s)" % (clock_text(), stage.sid, stamp["status"], stamp["result"], seconds),
          flush=True)
    return outcome


def record_fit(pipe, fit, stamp):
    """Mirror a fit's segments into pipeline.json (count, nfev and cost per segment)."""
    stored = read_json(pipe.dir / "pipeline.json")
    stored.setdefault("fits", {})[pipe.fit_name(fit)] = {
        "status": stamp["status"], "segments": [{k: s.get(k) for k in ("index", "status", "nfev", "cost", "dir",
                                                                       "seconds", "finished_utc")}
                                                for s in stamp.get("segments", [])],
        "nfev_used": sum(int(s.get("nfev", 0)) for s in stamp.get("segments", [])),
        "max_nfev": pipe.r["max_nfev"], "attempts": len(stamp.get("attempts", []))}
    write_json_atomic(pipe.dir / "pipeline.json", stored)
    pipe.stored = stored


@contextlib.contextmanager
def directory_lock(directory):
    """One writer per restart directory (``flock``; released when the process ends)."""
    import fcntl
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / ".pipeline.lock"
    stream = path.open("a+")
    try:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        stream.close()
        if exc.errno in (errno.EAGAIN, errno.EACCES):
            raise PipelineError("another cops_pipeline.py is working in %s (lock %s)" % (directory, path))
        raise
    try:
        stream.seek(0)
        stream.truncate()
        stream.write("%s pid %d job %s %s\n" % (socket.gethostname(), os.getpid(),
                                                os.environ.get("SLURM_JOB_ID"), utc_now()))
        stream.flush()
        yield
    finally:
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        stream.close()


# ----------------------------------------------------------------------------- commands
def load(args):
    directory = Path(args.restart_dir).resolve()
    manifest = directory / "pipeline.json"
    stored = read_json(manifest) if manifest.is_file() else None
    if stored is not None and stored.get("schema") != SCHEMA:
        raise PipelineError("%s is not a %s manifest" % (manifest, SCHEMA))
    request = build_request(args, stored)
    request["restart_dir"] = str(directory)
    if stored is not None:
        resolved = stored.get("resolved", {}).get("precision", {}).get("value")
        if request["precision"] == resolved and stored["identity"].get("precision") == "auto":
            request["precision"] = "auto"            # the explicit value the stored decision already took
    identity = identity_of(request, (stored or {}).get("identity"))
    if stored is not None:
        differ = identity_differences(stored["identity"], identity)
        if differ:
            raise PipelineError("%s was made for other campaign inputs (%s differ); use a new --restart-dir"
                                % (directory, ", ".join(differ)))
    pipe = Pipeline(request, stored, detect_hardware(args, request, live=args.command == "run"))
    pipe.identity = identity
    return pipe


def save_manifest(pipe, args):
    stored = pipe.stored or {"schema": SCHEMA, "schema_version": SCHEMA_VERSION, "created_utc": utc_now(),
                             "restart_dir": str(pipe.dir), "invocations": []}
    if (pipe.dir / "pipeline.json").is_file():
        current = read_json(pipe.dir / "pipeline.json")
        stored["fits"] = current.get("fits", stored.get("fits", {}))
    stored["identity"] = pipe.identity
    stored["request"] = pipe.r
    stored.setdefault("invocations", []).append({
        "utc": utc_now(), "command": args.command, "argv": sys.argv, "host": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "hardware": pipe.hw,
        "deadline_utc": wall_time.utc(pipe.end)})
    precision, reason = pipe.precision()
    stored["policies"] = {
        "record_precision": {"value": precision, "policy": RECORD_POLICIES.get(precision), "reason": reason},
        "gravity_precision": {"value": pipe.r["gravity_precision"],
                              "policy": GRAVITY_POLICIES[pipe.r["gravity_precision"]],
                              "reason": "requested" if args.gravity_precision else
                                        "default (user decision 2026-10-22)"}}
    if precision and pipe.r["precision"] == "auto" and not stored.get("resolved", {}).get("precision"):
        stored.setdefault("resolved", {})["precision"] = {
            "value": precision, "reason": reason, "utc": utc_now(), "hardware": pipe.hw,
            "estimate": pipe.estimate()[0]}
    write_json_atomic(pipe.dir / "pipeline.json", stored)
    pipe.stored = stored


def header(pipe):
    ident, r, hw = pipe.identity, pipe.r, pipe.hw
    precision, reason = pipe.precision()
    shape = ident["shape"]
    gpus, nodes = pipe.fit_hardware()
    lines = ["restart dir   %s%s" % (pipe.dir, "" if pipe.dir.exists() else " (new)"),
             "surface model %s (%s, %s, sha256 %s)" % (shape["path"], shape["format"],
                                                       ("%d facets" % shape["facets"]) if shape.get("facets")
                                                       else "facets from the base", shape["sha256"][:16]),
             "surface BC    %s" % r["surface_bc"],
             "fit input     %s" % r["fit_input"],
             "fits          %s, max_nfev %d%s" % (",".join(r["fits"]), r["max_nfev"],
                                                  (", tag " + r["tag"]) if r["tag"] else ""),
             "fine tuning   %s" % (("%s: --minutes %g --intervals %d --method %s" % (
                 ",".join(r["ft_fits"]), r["fine_tune_settings"]["minutes"], r["fine_tune_settings"]["intervals"],
                 r["fine_tune_settings"]["method"])) if r["fine_tune"] else "off"),
             "hardware      %d node(s) [%s] x %s GPU(s) [%s]; GPU memory %s [%s]; host memory %s [%s]; "
             "launcher %s" % (hw["nodes"], hw["nodes_source"], hw["gpus_per_node"], hw["gpus_source"],
                              gb(hw["gpu_memory_bytes"]) if hw["gpu_memory_bytes"] else "unknown",
                              hw["gpu_memory_source"],
                              gb(hw["host_memory_bytes"]) if hw["host_memory_bytes"] else "unknown",
                              hw["host_memory_source"], hw["launcher"]),
             "fit GPUs      %s on %d node(s)" % (gpus, nodes),
             "precision     records %s (%s)" % (("%s, %s" % (precision, RECORD_POLICIES[precision])) if precision
                                                  else "undecided", reason),
             "gravity       %s (%s)%s" % (r["gravity_precision"], GRAVITY_POLICIES[r["gravity_precision"]],
                                         "" if pipe.needs()["bundle"] else ", not needed by the requested fits"),
             "wall time     %s" % (("%s, margin %g s, %g s to save a fit; run fixes the deadline when it starts and "
                                    "passes it to the fits (--deadline)" % (r["wall_time"], r["wall_margin"],
                                                                            r["save_seconds"]))
                                   if r["wall_time"] else "no limit")]
    est, _ = pipe.estimate()
    if est is not None:
        lines.append("memory        %d records, %d facets; FP64 %s vs %s usable; FP32 %d resident, %d streamed; "
                     "host per node %s (%s pinned)" % (est["records"], est["faces"], gb(est["fp64_bytes"]),
                                                       gb(est["fp64_capacity_bytes"]), est["fp32_resident_records"],
                                                       est["fp32_streamed_records"], gb(est["host_bytes_rank0_node"]),
                                                       gb(est["pinned_bytes_per_node"])))
    if pipe.needs()["bundle"]:
        lines.extend(gravity_lines(pipe))
    refusal = pipe.memory_refusal(precision)
    if refusal:
        lines.append("REFUSED       %s" % refusal)
    if shape.get("topology_admission_failed") and not r["allow_genus"] and pipe.needs()["bundle"]:
        lines.append("WARNING       the mesh header says its gravity topology admission FAILED (genus > 0?): the "
                     "gravity stage refuses it unless --allow-genus N (an explicit experiment decision)")
    if shape.get("verified") is False:
        lines.append("note          surface model file absent here; identity taken from pipeline.json")
    return lines


def gravity_lines(pipe):
    """Rough gravity time per precision from recorded runs (or 'unknown')."""
    roots = [pipe.dir / "gravity"] + ([Path(pipe.r["output_root"])] if pipe.r["output_root"] else [])
    runs = gravity_runs([r for r in roots if r.is_dir()])
    records, faces = pipe.counts()
    targets = faces or pipe.identity["shape"].get("facets")
    body = (pipe.identity["gravity_shape"] or pipe.identity["shape"]).get("facets") or (
        faces if not pipe.identity["gravity_shape"] else None)
    gpus = pipe.hw["gpus_per_node"]
    lines = []
    for precision in ("mixed", "fp64"):
        mark = "*" if precision == pipe.r["gravity_precision"] else " "
        estimate = gravity_time(runs, precision, targets, body, gpus)
        if estimate is None:
            why = ("no recorded %s gravity run" % precision if targets and body and gpus else
                   "facet or GPU count unknown")
            lines.append("gravity time %s%-5s unknown (%s)" % (mark, precision, why))
        else:
            seconds, run = estimate
            lines.append("gravity time %s%-5s ~%s on %d GPU(s), scaled from %s (%d targets, %d GPU(s) of %s); "
                         "other GPU models differ" % (mark, precision, duration_text(seconds), gpus,
                                                      Path(run["path"]).name, run["targets"], run["devices"],
                                                      run["gpu"]))
    return lines


def cmd_plan(pipe, args):
    for line in header(pipe):
        print(line)
    for note in pipe.check_restarts():
        print("restart       %s" % note)
    if not pipe.r["base"] and not pipe.base_ready():
        for item in input_files_missing(pipe.r["fit_input"], pipe.reference):
            print("MISSING       base input file: %s" % item)
        executable = pipe.app / "native" / "liouville_cops_cuda"
        if not executable.is_file():
            print("MISSING       %s (make -C native ARCH=sm_80 on A100)" % executable)
    print()
    for stage in pipe.stages():
        state, message = stage_state(pipe, stage)
        print("%-16s %-8s %-11s %s%s" % (stage.sid, stage.group, state, stage.output if stage.needed else "",
                                          ("  [%s]" % message) if message else ""))
        if stage.needed and stage.command and state in ("pending", "incomplete", "blocked"):
            print("    $ " + command_text(stage))
    return 0


def cmd_status(pipe, args):
    groups = GROUPS if args.group == "all" else (args.group,)
    worst = 0
    if not args.check:
        print("restart dir %s" % pipe.dir)
    for stage in pipe.stages():
        try:
            state, message = stage_state(pipe, stage)
        except PipelineError as exc:
            state, message = "refused", str(exc)
        if stage.group in groups and stage.needed:
            if state == "refused" or (state not in ("complete", "external") and not resumable(pipe, stage, state)):
                worst = EXIT_FAILED
            elif state not in ("complete", "external") and worst == 0:
                worst = EXIT_RESUMABLE
        if args.check:
            continue
        stamp = read_stamp(pipe.dir, stage.sid) if not stage.external else None
        when = ""
        if stamp:
            when = "%s .. %s" % (stamp.get("started_utc", "")[:19], (stamp.get("finished_utc") or "")[:19])
            if stamp.get("status") == "running":
                when += " running on %s (job %s)" % (stamp.get("host"), stamp.get("slurm_job_id"))
        print("%-16s %-11s %-42s %s" % (stage.sid, state, when, message))
        for segment in (stamp or {}).get("segments", []) if stage.fit else []:
            print("%18s segment %d %-11s nfev %4d  cost %.6e  %s" % ("", segment["index"], segment["status"],
                                                                     segment["nfev"], segment["cost"],
                                                                     Path(segment["restart"])))
    if args.check:
        print({0: "complete", EXIT_RESUMABLE: "incomplete (resumable)", EXIT_FAILED: "failed"}[worst])
    return worst if args.check else 0


def cmd_run(pipe, args):
    groups = GROUPS if args.group == "all" else (args.group,)
    for note in pipe.check_restarts():
        print("restart: %s" % note)
    started = time.time()
    if pipe.r["wall_time"]:
        try:
            pipe.end = wall_time.resolve_deadline(pipe.r["wall_time"], margin=0.0, start=started)
        except ValueError as exc:
            raise PipelineError(str(exc))
        print("deadline %s (margin %g s)" % (wall_time.utc(pipe.end), pipe.r["wall_margin"]), flush=True)
    with directory_lock(pipe.dir), terminate_on_signal():
        save_manifest(pipe, args)
        tree_cache, attempted, failed, deferred = {}, set(), [], []
        while True:
            pipe.stored = read_json(pipe.dir / "pipeline.json")
            stages = pipe.stages()
            states = {stage.sid: stage_state(pipe, stage)[0] for stage in stages}
            todo = None
            for stage in stages:
                if stage.group not in groups or stage.sid in attempted:
                    continue
                if states[stage.sid] == "conflict":
                    attempted.add(stage.sid)
                    failed.append(stage.sid)
                    print("[%s] %s: not run, %s" % (clock_text(), stage.sid, stage_state(pipe, stage)[1]), flush=True)
                    continue
                if states[stage.sid] in ("pending", "incomplete", "blocked"):
                    todo = stage
                    break
            if todo is None:
                break
            attempted.add(todo.sid)
            blocking = [d for d in todo.depends if states.get(d) not in ("complete", "external", "not needed")]
            if blocking:
                print("[%s] %s: not started, needs %s" % (clock_text(), todo.sid,
                                                          ", ".join("%s (%s)" % (d, states.get(d)) for d in blocking)),
                      flush=True)
                failed.append(todo.sid)
                if todo.group == "prepare":
                    break
                continue
            if todo.pending:
                raise PipelineError("%s waits for the %s, which no completed stage provides" % (todo.sid, todo.pending))
            ok, why = can_start(pipe, todo)
            if not ok:
                print("[%s] %s: not started before the deadline: %s" % (clock_text(), todo.sid, why), flush=True)
                deferred.append(todo.sid)
                break
            outcome = execute(pipe, todo, tree_cache)
            if todo.sid == "base" and outcome == "complete":
                save_manifest(pipe, args)                # records the automatic precision decision
            if outcome == "deferred":
                deferred.append(todo.sid)
                break                                    # the deadline is near: leave the rest to the next job
            if outcome == "failed":
                failed.append(todo.sid)
                if todo.group == "prepare" or args.stop_on_failure:
                    break
    if failed:
        print("failed or not started: %s" % ", ".join(failed), flush=True)
        return EXIT_FAILED
    if deferred:
        print("stopped for the deadline at %s; run again (or the next job) to continue" % ", ".join(deferred),
              flush=True)
        return EXIT_RESUMABLE
    print("done: groups %s" % ", ".join(groups), flush=True)
    return 0


def render_script(pipe, group, path, next_path, gpus, nodes, cpus):
    r = pipe.r
    if group == "fit" and nodes > 1:
        tasks = ["#SBATCH --ntasks-per-node=%d" % gpus, "#SBATCH --gpus-per-task=1",
                 "#SBATCH --gpu-bind=closest             # each rank on the socket of its GPU",
                 ("#SBATCH --cpus-per-task=%d" % max(1, cpus // gpus)) if cpus else
                 "##SBATCH --cpus-per-task=SITE          # cores per node / GPUs per node"]
    else:
        tasks = ["#SBATCH --ntasks=1", "#SBATCH --gpus-per-node=%d" % gpus,
                 "#SBATCH --exclusive                    # all cores of the node (base threads, workers)"]
    site = ["#SBATCH --partition=%s" % r["slurm_partition"] if r["slurm_partition"] else "##SBATCH --partition=SITE",
            "#SBATCH --account=%s" % r["slurm_account"] if r["slurm_account"] else "##SBATCH --account=SITE"]
    modules = (["module load %s" % m for m in r["slurm_modules"]] if r["slurm_modules"] else
               ["# SITE: module load cuda/12.x openmpi/4.x   (the MPI mpi4py was built against)"])
    tool = shlex.quote(str(Path(__file__).resolve()))
    python = shlex.quote(r["python"])
    run = "%s -W ignore %s run --restart-dir %s --group %s" % (python, tool, shlex.quote(str(pipe.dir)), group)
    chain = ["# no resubmission (cops_pipeline.py submit --resubmit N adds it)"]
    if r["resubmit"]:
        status = "%s -W ignore %s status --restart-dir %s --group %s --check" % (python, tool,
                                                                                 shlex.quote(str(pipe.dir)), group)
        nxt = ("sbatch --export=ALL,COPS_PIPELINE_RESUBMIT_LEFT=%d %s" % (r["resubmit"], shlex.quote(str(next_path)))
               if next_path else 'echo "cops_pipeline: all groups complete"')
        chain = ["# continue (cops_pipeline.py submit --resubmit %d): resubmit this group while it is resumable,"
                 % r["resubmit"], "# then submit the next group",
                 "%s; state=$?" % status,
                 "left=${COPS_PIPELINE_RESUBMIT_LEFT:-%d}" % r["resubmit"],
                 "if [ $state -eq 0 ]; then",
                 "    %s" % nxt,
                 "elif [ $state -eq %d ] && [ $left -gt 0 ]; then" % EXIT_RESUMABLE,
                 "    sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL,COPS_PIPELINE_RESUBMIT_LEFT=$((left - 1)) %s"
                 % shlex.quote(str(path)),
                 "else",
                 '    echo "cops_pipeline: group %s not complete (status $state, resubmissions left $left)" >&2' % group,
                 "fi"]
    text = SLURM_TEMPLATE.read_text()
    values = {"@JOB_NAME@": "cops_%s" % group, "@NODES@": str(nodes if group == "fit" else 1),
              "@TASK_LINES@": "\n".join(tasks), "@TIME@": r["slurm_times"][group],
              "@OUTPUT@": str(pipe.dir / "slurm" / ("%s-%%j.out" % group)), "@SITE_LINES@": "\n".join(site),
              "@MODULE_LINES@": "\n".join(modules), "@APP@": shlex.quote(r["app_root"]),
              "@RUN@": run, "@CHAIN@": "\n".join(chain), "@GROUP@": group}
    for key, value in values.items():
        text = text.replace(key, value)
    return text


def cmd_submit(pipe, args):
    for note in pipe.check_restarts():
        print("restart: %s" % note)
    gpus = args.gpus_per_node if isinstance(args.gpus_per_node, int) else None
    cpus = None
    if gpus is None:
        gpus, cpus = sinfo_gres(pipe.r["slurm_partition"])
        if gpus is None:
            raise PipelineError("cannot read the GPUs per node of the partition (sinfo -h -o %G); give "
                                "--gpus-per-node N")
        print("GPUs per node from sinfo: %d%s" % (gpus, (", %d CPUs per node" % cpus) if cpus else ""))
    nodes = args.nodes or pipe.r.get("fit_nodes") or 1
    pipe.r["fit_nodes"], pipe.r["slurm_gpus_per_node"] = nodes, gpus
    with directory_lock(pipe.dir):
        save_manifest(pipe, args)
        stages = pipe.stages()
        states = {s.sid: stage_state(pipe, s)[0] for s in stages}
        wanted = [g for g in GROUPS if any(s.group == g and states[s.sid] in ("pending", "incomplete", "blocked")
                                           for s in stages)]
        if not wanted:
            print("every needed stage is complete; nothing to submit")
            return 0
        slurm = pipe.dir / "slurm"
        slurm.mkdir(exist_ok=True)
        paths = [slurm / ("%s.sbatch" % g) for g in wanted]
        for k, (group, path) in enumerate(zip(wanted, paths)):
            if path.exists():
                move_aside(path)
            path.write_text(render_script(pipe, group, path, paths[k + 1] if k + 1 < len(paths) else None,
                                          gpus, nodes, cpus))
            path.chmod(0o755)
    print("batch scripts in %s (%d node(s) for the fits, %d GPUs per node; set the SITE lines first)."
          % (slurm, nodes, gpus))
    if pipe.r["resubmit"]:
        print("Submit the first job; each job resubmits itself (afterany, at most %d times) while its group is "
              "resumable and then submits the next group:" % pipe.r["resubmit"])
        print("  sbatch %s" % shlex.quote(str(paths[0])))
    else:
        print("Submit with (bash):")
        previous = None
        for group, path in zip(wanted, paths):
            variable = "J_" + group.upper()
            dependency = " --dependency=afterok:$%s" % previous if previous else ""
            print("  %s=$(sbatch --parsable%s %s)" % (variable, dependency, shlex.quote(str(path))))
            previous = variable
    print("Each job runs `cops_pipeline.py run --restart-dir %s --group <group>`, which resumes at the first "
          "incomplete stage, so a failed job can be resubmitted unchanged." % pipe.dir)
    return 0


def gpus_option(value):
    if value == "auto":
        return value
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("--gpus-per-node takes a number or 'auto'")
    if number < 1:
        raise argparse.ArgumentTypeError("--gpus-per-node must be >= 1")
    return number


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="cops_pipeline.py", description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="plan: stages, their state and exact commands (writes nothing, no GPU). run: execute the "
               "incomplete stages (inside an allocation). submit: write Slurm scripts per group and print "
               "the sbatch commands (does not submit). status: progress of a restart directory "
               "(--check: exit 0 complete, 75 resumable, 1 failed). See tools/cluster/README.md.")
    p.add_argument("command", choices=("plan", "run", "submit", "status"))
    g = p.add_argument_group("campaign")
    g.add_argument("--shape", help="surface model: NASTRAN .bdf or SPICE DSK (required for a new restart dir)")
    g.add_argument("--shape-format", choices=("nastran", "dsk"), help="mesh reader (default: from the extension)")
    g.add_argument("--surface-bc", choices=("flat", "rough"),
                   help="angular emission law of the base: flat (Lambertian) or rough (uniform over the hemisphere)")
    g.add_argument("--restart-dir", required=True,
                   help="all products, stamps, logs and fit_restart.npz files; rerun to resume")
    g.add_argument("--max-nfev", type=int, help="optimizer evaluation limit per fit (default 2000)")
    g.add_argument("--fine-tune", dest="fine_tune", action="store_true", default=None,
                   help="fine-tune each fit's restart (late interval; tools/fine_tune_gpu.py)")
    g.add_argument("--no-fine-tune", dest="fine_tune", action="store_false", help="no fine tuning (default)")
    f = p.add_argument_group("fits")
    f.add_argument("--fits", help="comma-separated subset of B29,N30,C30,NC31 (default all)")
    f.add_argument("--fit-args", help="fitter arguments replacing the production set: '%s'" % DEFAULT_FIT_ARGS)
    f.add_argument("--tag", help="suffix of the fit directory names (other settings in the same restart dir)")
    f.add_argument("--restart-from", action="append",
                   help="start a fit from a restart file: FIT=FILE (repeatable), or FILE with one fit")
    f.add_argument("--translate-restart", action="store_true", default=None,
                   help="translate a --restart-from file of another mesh (restart_tool.py translate)")
    f.add_argument("--precision", choices=("fp64", "mixed", "auto"),
                   help="kernel records: fp64 (validated, resident only), mixed (FP32 records, FP64 arithmetic, "
                        "streams what does not fit) or auto (default: fp64 when it fits the GPUs)")
    f.add_argument("--gravity-precision", choices=("fp64", "mixed"),
                   help="gravity field: mixed (default, user decision 2026-10-22: FP64 near, FP32 far quadrature) "
                        "or fp64 (exact reference); independent of --precision")
    f.add_argument("--activity", choices=("host", "device"), help="gpu_fit --activity (default: the tool's)")
    f.add_argument("--transport", choices=("local", "mpi"), help="fits: default local on one node, mpi on several")
    f.add_argument("--headroom-fraction", type=float, help="device memory headroom (default: the tools', 0.20)")
    w = p.add_argument_group("wall time (default: no limit)")
    w.add_argument("--wall-time", help="HH:MM:SS from the start of run, 'slurm' (the job's end) or 'none'")
    w.add_argument("--wall-margin", type=float, help="safety margin before the deadline [s] (default 300)")
    w.add_argument("--save-seconds", type=float, help="time a fit reserves to write its results [s] (default 600)")
    t = p.add_argument_group("fine tuning")
    t.add_argument("--ft-minutes", type=float, help="window length [min] (default 5)")
    t.add_argument("--ft-intervals", type=int, help="number of intervals (default 5)")
    t.add_argument("--ft-method", choices=("cpu", "multiplier"), help="default cpu (the CPU branch's method)")
    t.add_argument("--ft-max-nfev", type=int, help="fine-tuning evaluation limit (default 1000)")
    t.add_argument("--ft-fits", help="fits to fine-tune (default all requested)")
    x = p.add_argument_group("existing products (identity-checked, not rebuilt)")
    x.add_argument("--base", help="fitting base directory (rebuild.json complete, same mesh and SurfaceBC)")
    x.add_argument("--store", help="FP32 record store of that base")
    x.add_argument("--visibility-bundle", help="visibility template bundle (.npz) of that base")
    x.add_argument("--bundle", help="full template bundle (.npz, with cliff columns) of that base")
    x.add_argument("--gravity-shape", help="gravity body (default --shape)")
    x.add_argument("--allow-genus", type=int, help="admit a gravity body of genus up to N (SHAP8: 4; experiment)")
    x.add_argument("--fit-input", help="base input (default <reference>/fit_base.in)")
    e = p.add_argument_group("hardware (detected at run time; these override)")
    e.add_argument("--nodes", type=int, help="nodes (run: of this allocation; submit: of the fit job, default 1)")
    e.add_argument("--gpus-per-node", type=gpus_option,
                   help="GPUs per node (run: detected; submit: 'auto' (default) reads the partition's GRES)")
    e.add_argument("--gpu-memory-gb", type=float, help="memory per GPU [GB] (run: detected)")
    e.add_argument("--host-memory-gb", type=float, help="memory per node [GB] (default SLURM_MEM_PER_NODE or MemTotal)")
    e.add_argument("--threads", type=int, help="base build threads (default $SLURM_CPUS_ON_NODE or all cores)")
    e.add_argument("--launcher", choices=("srun", "mpirun"), help="MPI launcher (default srun inside Slurm)")
    v = p.add_argument_group("environment")
    v.add_argument("--output-root", help="COPS_OUTPUT_ROOT for the tools (default: $COPS_OUTPUT_ROOT, else the restart dir)")
    v.add_argument("--app-root", help="LiouvilleCOPS_cuda checkout (default: this tree)")
    v.add_argument("--reference", help="LiouvilleCOPS tree (default: $LIOUVILLE_COPS_REFERENCE, else ../LiouvilleCOPS)")
    v.add_argument("--python", help="CUDA environment interpreter (default <app>/.venv/bin/python)")
    v.add_argument("--cpu-python", help="interpreter of rebuild_cops_boundary.py (default --python)")
    v.add_argument("--store-threads", type=int, help="record store reader threads (default 8)")
    v.add_argument("--srun-args", help="srun options of MPI fits (default '%s')" % SRUN_ARGS)
    v.add_argument("--group", choices=GROUPS + ("all",), default="all", help="run, status: only this group")
    v.add_argument("--stop-on-failure", action="store_true", help="run: stop at the first failed fit")
    v.add_argument("--check", action="store_true", help="status: exit 0 complete, 75 resumable, 1 failed")
    s = p.add_argument_group("submit (Slurm)")
    s.add_argument("--resubmit", type=int, help="jobs resubmit themselves while resumable, at most N times "
                                                "(needs --wall-time)")
    s.add_argument("--slurm-partition")
    s.add_argument("--slurm-account")
    s.add_argument("--slurm-module", action="append", help="module to load in the jobs (repeatable)")
    for group in GROUPS:
        s.add_argument("--time-" + group, help="time limit of the %s job (default %s)" % (group, SUBMIT_TIMES[group]))
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        pipe = load(args)
        return {"plan": cmd_plan, "run": cmd_run, "submit": cmd_submit, "status": cmd_status}[args.command](pipe, args)
    except PipelineError as exc:
        print("cops_pipeline: %s" % exc, file=sys.stderr, flush=True)
        return EXIT_REFUSED
    except (Interrupted, KeyboardInterrupt) as exc:
        print("cops_pipeline: interrupted (%s); rerun to resume" % (str(exc) or exc.__class__.__name__),
              file=sys.stderr, flush=True)
        return EXIT_RESUMABLE


if __name__ == "__main__":
    sys.exit(main())
