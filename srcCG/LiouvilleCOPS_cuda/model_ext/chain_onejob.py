#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/chain_onejob.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/chain_onejob.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 560780acebbeea6d75d59e4bee1d02ca7388d20749f1303005163c6f2b4bf9a2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. One Slurm job, one process: a whole fit / fine-tuning chain (c1 fit -> c1 ft -> ... -> c<N> ft) on ONE
#   bound multi-GPU backend (user decision 2026-11-10 "yes, implement (b)"): the records are bound once, the template
#   basis is reset before every later stage, each stage runs its split job script's logic (segments, resume, move-aside,
#   final links, product checks) in this process, every model is verified against the bound dataset identity; process
#   kinds positive | slope | nonneg (v6, 2026-11-11). This copy = v6.2 (frozen copy
#   model_ext_frozen/chain_onejob_20261111T060742Z).
# USED BY: jobs/chain_onejob.sbatch (submitted with jobs/chain_onejob.sh): nonneg GPU tests 124143, 124148 (bitwise =
#   the split chains, compares 124160, 124161) and the production chains 124162 (softplus) and 124165 (squared),
#   2026-11-11. The positive-kind tests 123990, 124091-124095 ran v3 / v5 (backups *.before_*).
# 2026-11-13: v6.4: --slope-mult-min V (kind nonneg, softplus, value 1) passed to every fit argv of every cycle and
#   resumed segment, at the end of MOPT as jobs/nonneg_fit.sbatch; frozen chain_onejob_20261113T041215Z (jobs 124664,
#   124772).
# ==== end of header ====
"""One Slurm job, one driver process: a whole fit / fine-tuning chain on ONE bound multi-GPU backend
(user decision 2026-11-10, "yes, implement (b)").

    .venv/bin/python -W ignore <job provenance>/code/chain_onejob.py --chain-root <L> --lambda 1e-4 --ratio 1e-5 \\
        --cycles 3 --source <cycle-1 start fit_restart.npz> --nfev 4000 --check 1 --ft-maxnfev 500 \\
        --ft-facets late_only --ft-window-min 5 --job-provenance <L>/provenance-<job> --job-script <sbatch> \\
        --model-ext <MODEL_EXT source> [--c1-start-set kappa=1] [--plan-only]

Job script: jobs/chain_onejob.sbatch (environment, job-level provenance, the code copy this file runs from); submit
helper: jobs/chain_onejob.sh (dry run by default). CPU tests and the split-vs-one-job comparison:
model_ext/test_chain_onejob.py. Local, untracked extension: nothing in the git tree and no other model_ext file is
edited; everything below is done by run-time patches of module attributes, from this file.

Why. The split chain (jobs/positive_cycles.sh) runs every stage as its own job; each job binds the 4.74e9 kernel
records of the 12M lategrid base on 8 A100s again (118-171 s per job, ~14 min of a 3 h chain). Here one process binds
once and runs every stage of the chain on that backend, with the products, directory layout and results of the split
chain.

Stages and layout (the split chain's, jobs/positive_cycles.sh; L = --chain-root, production
<run>/positive_r<R>/lam<L>):
    c1fit -> c1ft -> c2fit -> c2ft -> ...        L/c<k>/fit/seg<NN>, L/c<k>/fit/final -> seg<NN>
                                                 L/c<k>/ft/seg<NN>,  L/c<k>/ft/final, L/c<k>/ft/<product> -> final/...
* fit stage (PositiveFitStage): the logic of jobs/positive_fit.sbatch with CHAIN=1, ported line by line (segment
  state, START_RESTART / START_SET checks, segment and evaluations left, START / GRAIN / HOLD / CHK arguments, the
  final link); the command is the one positive_fit.sbatch would run (positive_fit.py ... -- <base> ...), executed IN
  THIS PROCESS (positive_fit.main()).
* fine-tuning stage (PositiveFineTuneStage): the logic of jobs/positive_finetune.sbatch (environment of the fit, OUT
  checks, segment state, move-aside of a segment without state, resume, product check, links), executed in this
  process (positive_finetune.main(argv)); its validation (check 1 reproduces the fit's cost, ...) stays on.
Cycle k > 1 fit: START_RESTART = c<k-1>/fit/final/fit_restart.npz, HELD = GRAIN_FROM = c<k-1>/ft/finetune_state.npz
(the same paths as positive_cycles.sh).

The bound backend (one per chain)
* The first stage that evaluates binds as usual (gpu_model.REGISTRY in a fit, positive_finetune.GpuEngine in a fine
  tuning). gpu_model.REGISTRY gets three instance attributes for the whole chain: close (no-op between stages; the
  real close runs at the end of the chain), configure (accepts only the configuration already bound) and backend_for
  (refuses a second dataset; verifies EVERY model object that reaches the backend, in every stage, against the bound
  dataset identity (dataset_identity(): base, rows and point ids, record counts, every COPS_model.dat column except
  the selection columns useNG/useRG/refNG/refRG, ModelSettings and per-point constants, coefficients, digest of Y; the
  fit model M of fit_cops.main differs from MF in the selection columns and header messages and shares the backend,
  as in the split chain, where _Registry.backend_for keys on run, rows, devices, mode only).
* Before every stage after the first bind the template basis is reset (MultiGpuBackend.rebind_basis(None): the
  evaluator arrays are rebound on the kept records, ~1.5 s), i.e. the backend is in the state of a fresh bind with
  basis None, which is how both split stages start (a fit then rebinds its template basis itself, as in the split
  chain). The registry's call counters, seconds and events start empty in every fit stage, plus one event
  "chain_onejob_reuse" (fit_manifest.json compute_backend.registry).
* Fine tuning: positive_finetune.GpuEngine is subclassed per stage: during its __init__, cops_cuda.multi_gpu.
  MultiGpuBackend is replaced by a factory that returns the bound backend (after the same verification; it binds and
  registers a new one when none is bound yet, e.g. a chain resumed at a fine tuning); its close() releases the
  worker-side fine-tuning state (request chain_ft_detach: pinned host registrations, device arrays, budget
  reservations, shared-memory handles) instead of closing the backend.
* Spawned GPU workers import this file as __mp_main__: worker_setup() does what the split stages' main modules do in
  their workers (import positive_fit: COPS_POS_* + positive_activity.install(); positive_finetune._process_setup():
  the ft_* requests) and adds the chain_ft_detach request.

Module state: every stage starts from the state of a fresh process
* Process level (once, before the baseline): the positive environment (positive_problem.set_positive_environment) and
  positive_activity.install(); every cops_cuda module is imported. Baseline snapshot: sys.modules, the attribute
  dictionaries of every cops_cuda module and of positive_activity, the dictionaries of every class defined in them,
  sys.path, sys.argv, os.environ, the working directory, numpy error/print options, warnings filters.
* After every stage (also after a failure): every module imported during the stage whose file lies in the code
  directory (model_ext copy), APP/tools or the reference tools is removed from sys.modules (the next stage imports
  positive_fit, positive_problem, grain_fit, grain_crossing, gpu_fit, fit_cops, ... afresh: new globals P,
  ACTIVE_MANIFEST, HELD, wrappers); module and class attributes changed or added in the kept modules are restored or
  deleted (e.g. the held-term patch of CudaBackend / MultiGpuBackend.activity_linear, figures.read_table,
  restart.from_fit, gpu_model._CLASS, host_model._CLASS, reference._bound_root); sys.path, argv, environment, cwd,
  numpy and warnings state are restored; matplotlib figures are closed and its rc parameters reset to the rc file's;
  gc.collect().
* Per stage the environment is the baseline plus the variables the split job script would have exported (LAMBDA,
  RATIO, CYCLE, START_RESTART, ... / FIT, OUT, MAXNFEV, FACETS, ...; only environment.txt reads them).
* Kept on purpose (process level, as in any process): the bound backend and its workers, positive_activity's patch
  and _STATE, compiled kernels, imported third-party modules.

Logs, provenance
* fds 1/2 go through a pipe read by a thread (Tee): everything (this process and the GPU workers, which inherit the
  pipe) goes to the job log (<job provenance>/chain.log via the job script) and, during a stage, also to the stage
  log <stage dir>/provenance-<job>-<stage>/fit.log | finetune.log.
* <stage dir>/provenance-<job>-<stage>/: commands.sh (exports + the stage's exact command; run as is it reruns the
  stage alone in a new process), environment.txt, segment.txt, git.txt (per stage), backend.json (bind / reuse,
  seconds), the code copies of the split stage (+ chain_onejob.py) with model_ext.sha256 / run_copy.sha256, the job's
  module_list.txt, nvidia_smi_q.txt, toolchain.txt, pip_list.txt, chain_onejob.sbatch, free/df at start and end,
  (fine tuning) fit.sha256, nvidia_smi_samples.csv, free_samples.txt.
* <L>/chain_onejob_events.jsonl: one line per stage start and end (job, stage, segment directory, exit code); read on
  a resume (a fit segment that ended with 75 but has no restart file, i.e. wall time exhausted before its first
  evaluation, is moved aside and run again; any other fit segment without restart file stops the chain with 3, the
  CHAIN=1 rule of the split chain).
* <job provenance>/chain_summary.json: per stage status, exit code, segment, wall seconds, bind / reuse seconds.

Wall time, resume, exit codes
* Every stage uses its own --wall-time slurm mechanism (fit: gpu_fit/WallClock, stops at the best evaluated state,
  restart "interrupted"; fine tuning: Tracker, stop "wall_time"), so the stage stops before the job's end and the
  chain exits 75. Before a stage is started, at least --min-fit-minutes / --min-ft-minutes (+ --bind-minutes when
  nothing is bound) must remain before the job's end, else the chain exits 75 without touching the stage.
* Resubmitting the same command: completed stages are skipped after their checks (fit: restart status complete or no
  evaluations left, final link, start restart / held identity; fine tuning: the product check of
  positive_finetune.sbatch, links, the tuned fit's restart SHA-256); the current stage continues with its next
  segment exactly as the split job scripts do.
* Exit: 0 every stage complete; 75 wall time (resubmit); otherwise the failing stage's code (2 usage, 3 inconsistent
  state, 4 product check, 1 error). Nothing is deleted; moved-aside segments are kept as seg<NN>.incomplete_<UTC>.

Stage interface (for later stage types, e.g. the slope-multiplier fit and its fine tuning, user option 3a)
* A process kind (PROCESS_KINDS) owns the process-level patches: parent_setup(args) once before the baseline,
  worker_setup() in every spawned GPU worker (env COPS_CHAIN_ONEJOB_KIND). One kind per chain (positive and slope
  patch the same evaluator methods).
* A stage type (STAGE_TYPES) implements: name, directory, log_name, code_files, summary_regex, prepare() -> Plan
  (done / rc / argv / env / provenance notes), call(plan) -> exit code (the stage driver's main, in process, fresh
  modules), finish(plan, rc) -> exit code (final link, product checks). See README_ONEJOB_STAGE_INTERFACE below and
  model_ext/logs/agent_onejob_READY.txt.
"""

import argparse
import datetime as dt
import gc
import hashlib
import importlib
import json
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import threading
import time
import traceback
import warnings
import weakref
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
REFERENCE_TOOLS = APP.parent / "LiouvilleCOPS" / "tools"
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP / "tools"), str(APP)):          # as positive_fit.py / positive_finetune.py
    if _p not in sys.path:
        sys.path.insert(0, _p)

ROOT_OUT = Path("/data001/highresolu_vt/cops_output")
RUN = ROOT_OUT / "runs" / "shap8_012m_lategrid_dz1.5"
BUNDLE = ROOT_OUT / "runs" / "shap8_012m_flat_earlier" / "bundle" / "restart_templates.npz"
PY = [str(APP / ".venv" / "bin" / "python"), "-W", "ignore"]
KIND_ENV = "COPS_CHAIN_ONEJOB_KIND"
EXIT_WALL_TIME = 75
RATIO_DECISIONS = ("0.02", "0.002", "1e-5")                    # positive_fit.sbatch / positive_cycles.sh
EVENTS_NAME = "chain_onejob_events.jsonl"
DETACH_REQUEST = "chain_ft_detach"

README_ONEJOB_STAGE_INTERFACE = """\
A stage type for chain_onejob.py (positive fit / fine tuning now; slope-multiplier fit and its fine tuning later):
 1. prepare(): the split job script's logic (state of the stage directory, segment, evaluations left, move-aside,
    links) -> Plan(done | rc | argv, env, out); no GPU, no stage-module import.
 2. call(plan): the stage's driver main(argv) -> exit code, in this process; its modules are imported afresh (they
    are purged after every stage). It must not call sys.exit() other than through SystemExit, nor os._exit.
 3. The stage driver must reach the backend only through cops_cuda.gpu_model.REGISTRY (fits: automatic) or through a
    cops_cuda.multi_gpu.MultiGpuBackend(spec, model, None, 1.0, devices=..., mode="tree", records_per_block=...,
    transport="local", activity="host") constructed from the base's rows in order (fine tunings: the chain swaps the
    class for a factory during the engine's __init__) and must close it only through its engine's close().
 4. Worker-side state of a stage: one object on the WorkerState (positive_finetune: worker._positive_finetune), every
    device array, budget reservation, cudaHostRegister'ed buffer and shared-memory handle reachable from it, released
    by a detach request (chain_ft_detach here) without touching the backend's records.
 5. finish(plan, rc): the split job script's completion logic (product check, final link) -> exit code.
 6. Wall time: --wall-time slurm in the stage, stop at the best evaluated state, exit 75, resumable segments."""


def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_stamp():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def say(message):
    print("chain_onejob: %s" % message, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def readlink_f(path):
    return Path(os.path.realpath(str(path)))


def write_text(path, text):
    Path(path).write_text(text)


def run_to_file(cmd, path, shell=False, append=False):
    """A provenance command (git, free, df, du, nvidia-smi): stdout + stderr into ``path``; never raises."""
    try:
        with open(path, "a" if append else "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, shell=shell, timeout=600, cwd=str(APP))
    except Exception as exc:                       # provenance never stops a chain
        with open(path, "a") as f:
            f.write("chain_onejob: %s failed: %s\n" % (cmd, exc))


def shell_line(argv):
    return " ".join(shlex.quote(str(a)) for a in argv)


def call_main(function):
    """Exit code of ``function()`` as the process ``sys.exit(main())`` would have it."""
    try:
        code = function()
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        print(code, file=sys.stderr, flush=True)
        return 1
    except KeyboardInterrupt:
        traceback.print_exc()
        return 130
    except BaseException:                          # as the interpreter: traceback, exit 1
        traceback.print_exc()
        sys.stderr.flush()
        return 1
    if code is None:
        return 0
    if isinstance(code, int):
        return code
    print(code, file=sys.stderr, flush=True)
    return 1


def visible_devices():
    """'0,1,...,N-1' as CUDA numbers the job's GPUs (the split job scripts' DEV)."""
    cvd = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if cvd:
        n = len([d for d in cvd.split(",") if d.strip() != ""])
    else:
        out = subprocess.run(["nvidia-smi", "-L"], stdout=subprocess.PIPE, universal_newlines=True).stdout
        n = len([line for line in out.splitlines() if line.strip()])
    return ",".join(str(k) for k in range(n))


# ================================================================================================ dataset identity
# Columns of COPS_model.dat that fit_cops' data selection rewrites (named profile, --data, saturated NG) and that the GPU
# evaluator never reads (fit_cops.Model keeps them as useNG/useRG/refNG/refRG for the residuals only). The fit model M
# (fit_cops.main, after the selection) and the all-rows model MF differ exactly there and in the header's message lines
# ("selection profile ...") -- they share one backend in the split chain too (gpu_model._Registry.backend_for keys on
# run, rows, devices, mode only). Job 123990 (v3) refused M; v4 checks what identifies the records and the evaluator.
SELECTION_COLUMNS = ("useNG", "useRG", "refNG", "refRG")


def _array_digest(a):
    a = np.ascontiguousarray(np.asarray(a))
    return hashlib.sha256(str(a.dtype).encode() + str(a.shape).encode() + a.tobytes()).hexdigest()


def dataset_identity(model, y_digest=None):
    """What a bound backend's workers and evaluator depend on: the base (run), the rows and their point ids, the
    per-point record counts, every COPS_model.dat column except SELECTION_COLUMNS (geometry, kernels, T_RG, ...), the
    header-derived ModelSettings and per-point constants (rh, VdotB, T_RG, bm, cNG, cRG), the activity coefficients
    and the element basis Y (digest). Free-text header lines are not part of it (ModelSettings carries the header's
    model keys)."""
    from cops_cuda import backend as bk
    from cops_cuda import provenance
    spec = model._spec
    cols = spec["cols"]
    keep = sorted((int(i), str(n)) for n, i in cols.items() if n not in SELECTION_COLUMNS)
    rows = [int(r) for r in model.rows]
    ident = {"run": os.path.realpath(str(spec["run"])), "rows": rows,
             "columns": [n for _, n in keep],
             "row_values": hashlib.sha256(json.dumps([[r[i] for i, _ in keep] for r in spec["rows"]]).encode()
                                          ).hexdigest(),
             "settings": bk.ModelSettings.from_model(model).as_dict(),
             "coef": _array_digest(spec["coef"]),
             "ipoint_of_rows": _array_digest(np.asarray(model.ipoint)[rows]),
             "Y": y_digest or provenance.array_sha256(model.Y)}
    if getattr(model, "indptr", None) is not None:
        ident["record_counts"] = _array_digest(np.diff(np.asarray(model.indptr)))
    for name in ("rh", "VdotB", "TRG"):
        if getattr(model, name, None) is not None:
            ident[name] = _array_digest(getattr(model, name))
    for name in ("bm", "cNG", "cRG"):
        value = getattr(model, name, None)
        if isinstance(value, dict):
            ident[name] = {str(k): _array_digest(v) for k, v in sorted(value.items())}
    return ident


def identity_differences(a, b):
    """Fields of two dataset identities that differ (a field missing on one side counts as different; every
    fit_cops.Model has them all: Model.__init__ calls load_kernels on every row)."""
    return sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))


# ================================================================================================ worker side
def _detach_finetune(worker):
    """Release positive_finetune's worker state (_WorkerFineTune) without touching the backend's records."""
    import cupy
    ft = getattr(worker, "_positive_finetune", None)
    if ft is None:
        return {"attached": False, "pid": os.getpid()}
    info = {"attached": True, "pid": os.getpid(), "device_id": int(ft.be.device_id), "errors": []}
    with cupy.cuda.Device(ft.be.device_id):
        cupy.cuda.Device(ft.be.device_id).synchronize()
        unregistered = 0
        for ptr in list(getattr(ft, "pinned", []) or []):
            try:
                cupy.cuda.runtime.hostUnregister(ptr)
                unregistered += 1
            except Exception as exc:              # recorded
                info["errors"].append("hostUnregister: %s" % exc)
        ft.pinned = []
        info["host_unregistered"] = unregistered
        for name in ("d_g", "d_seen", "d_mask", "d_bad", "d_vpar", "d_gpar", "d_lsel", "d_lpar",
                     "module", "k_jtu", "k_seen", "k_dup", "k_map"):
            if hasattr(ft, name):
                delattr(ft, name)
    released = 0
    for name in ("reservation_par", "reservation"):
        r = getattr(ft, name, None)
        if r is not None:
            r.release()
            released += 1
            setattr(ft, name, None)
    info["reservations_released"] = released
    for name in ("v", "g", "vpar", "out", "gpar", "gmine"):
        if hasattr(ft, name):
            delattr(ft, name)
    closed = 0
    for s in list(getattr(ft, "_shm", []) or []) + list(getattr(ft, "_pshm", []) or []):
        try:
            s.close()
            closed += 1
        except Exception as exc:                  # e.g. BufferError: recorded, the mapping stays until exit
            info["errors"].append("shm close %s: %s" % (getattr(s, "name", "?"), exc))
    info["shm_closed"] = closed
    ft.be._columns_cache = None
    worker._positive_finetune = None
    del ft
    gc.collect()
    with cupy.cuda.Device(worker.backend.device_id):
        free, total = cupy.cuda.runtime.memGetInfo()
        pool = cupy.get_default_memory_pool()
        info.update(gpu_free_bytes=int(free), gpu_total_bytes=int(total), pool_used_bytes=int(pool.used_bytes()),
                    pool_held_bytes=int(pool.total_bytes()))
    return info


def install_detach_handler():
    """The chain_ft_detach request in this process's cops_cuda.transport.WorkerState.handle (idempotent)."""
    from cops_cuda import transport
    current = transport.WorkerState.handle
    if getattr(current, "_chain_onejob", False):
        return False
    original = current

    def handle(self, message):
        if message and message[0] == DETACH_REQUEST:
            try:
                return ("ok", _detach_finetune(self))
            except Exception:
                return ("error", traceback.format_exc())
        return original(self, message)

    handle._chain_onejob = True
    transport.WorkerState.handle = handle
    return True


def _worker_setup_positive():
    """What the split stages' main modules do in their spawned workers: positive_fit.py's module level (COPS_POS_*
    from the inherited environment, positive_activity.install()) and positive_finetune._process_setup() (the ft_*
    requests); then the chain's detach request."""
    import positive_fit  # noqa: F401  (module level: set_positive_environment, positive_activity.install)
    import positive_finetune
    positive_finetune._process_setup()
    install_detach_handler()


def worker_setup():
    kind = os.environ.get(KIND_ENV)
    if kind is None:
        return
    if kind not in PROCESS_KINDS:
        raise RuntimeError("chain_onejob: unknown %s=%s" % (KIND_ENV, kind))
    PROCESS_KINDS[kind]["worker"]()


# ================================================================================================ process state
class Baseline:
    """Snapshot of the process state before the first stage; restore() after every stage (module docstring)."""

    KEEP_PREFIXES = ("cops_cuda",)
    KEEP_NAMES = ("positive_activity",)

    def __init__(self, purge_roots, keep_names=None):
        self.purge_roots = [os.path.realpath(str(p)) for p in purge_roots]
        if keep_names is not None:                 # the process kind's patch module (positive_activity, slope_mult)
            self.KEEP_NAMES = tuple(keep_names)
        self.captured = False

    def _kept(self, name):
        return name in self.KEEP_NAMES or any(name == p or name.startswith(p + ".") for p in self.KEEP_PREFIXES)

    def capture(self):
        import pkgutil
        import cops_cuda
        for m in pkgutil.iter_modules(cops_cuda.__path__):
            importlib.import_module("cops_cuda." + m.name)
        self.modules = set(sys.modules)
        self.kept = {n: m for n, m in sys.modules.items() if self._kept(n) and m is not None}
        self.module_dicts = {n: dict(vars(m)) for n, m in self.kept.items()}
        self.class_dicts = {}
        for n, m in self.kept.items():
            for v in list(vars(m).values()):
                if isinstance(v, type) and self._kept(getattr(v, "__module__", "") or "") and id(v) not in self.class_dicts:
                    self.class_dicts[id(v)] = (v, dict(vars(v)))
        self.sys_path = list(sys.path)
        self.argv = list(sys.argv)
        self.environ = dict(os.environ)
        self.cwd = os.getcwd()
        self.np_err = np.geterr()
        self.np_print = np.get_printoptions()
        self.warn_filters = list(warnings.filters)
        self.captured = True
        return {"modules": len(self.modules), "kept_modules": sorted(self.kept), "classes": len(self.class_dicts)}

    def _under_roots(self, module):
        f = getattr(module, "__file__", None)
        if not f:
            return False
        f = os.path.realpath(f)
        return any(f.startswith(r + os.sep) for r in self.purge_roots)

    def restore(self):
        if not self.captured:
            raise RuntimeError("Baseline.restore before capture")
        report = {"purged": [], "restored": [], "deleted": [], "new_kept_modules": []}
        for name in list(sys.modules):
            if name in self.modules:
                continue
            module = sys.modules.get(name)
            if module is None or self._under_roots(module):
                del sys.modules[name]
                report["purged"].append(name)
            elif self._kept(name):                 # a cops_cuda module first imported in a stage: fresh next time
                del sys.modules[name]
                report["new_kept_modules"].append(name)
        missing = object()
        for name, snap in self.module_dicts.items():
            module = self.kept[name]
            current = vars(module)
            for key in [k for k in current if k not in snap]:
                delattr(module, key)
                report["deleted"].append("%s.%s" % (name, key))
            for key, value in snap.items():
                if current.get(key, missing) is not value:
                    setattr(module, key, value)
                    report["restored"].append("%s.%s" % (name, key))
        for cls, snap in self.class_dicts.values():
            current = dict(vars(cls))
            for key in [k for k in current if k not in snap]:
                delattr(cls, key)
                report["deleted"].append("%s.%s.%s" % (cls.__module__, cls.__qualname__, key))
            for key, value in snap.items():
                if key in ("__dict__", "__weakref__"):
                    continue
                if current.get(key, missing) is not value:
                    setattr(cls, key, value)
                    report["restored"].append("%s.%s.%s" % (cls.__module__, cls.__qualname__, key))
        sys.path[:] = self.sys_path
        sys.argv[:] = self.argv
        os.environ.clear()
        os.environ.update(self.environ)
        os.chdir(self.cwd)
        np.seterr(**self.np_err)
        np.set_printoptions(**self.np_print)
        warnings.filters[:] = self.warn_filters
        if hasattr(warnings, "_filters_mutated"):
            warnings._filters_mutated()
        if "matplotlib.pyplot" in sys.modules:
            sys.modules["matplotlib.pyplot"].close("all")
        if "matplotlib" in sys.modules:
            sys.modules["matplotlib"].rc_file_defaults()
        gc.collect()
        return report


class Tee:
    """fds 1 and 2 -> a pipe; a thread copies everything to the job log (the original fd 1) and, while a stage runs,
    to its stage log. Stage switches are synchronized with a marker written through the pipe."""

    MARK = b"\x1e@chain_onejob_switch_%d@\x1e"

    def __init__(self):
        sys.stdout.flush()
        sys.stderr.flush()
        self.job_fd = os.dup(1)
        self.err_fd = os.dup(2)
        r, w = os.pipe()
        os.dup2(w, 1)
        os.dup2(w, 2)
        os.close(w)
        self.r = r
        self.stage = None
        self.pending = {}
        self.counter = 0
        self.events = {}
        self.ended = threading.Event()
        self.lock = threading.Lock()
        self.thread = threading.Thread(target=self._pump, name="chain_onejob_tee", daemon=True)
        self.thread.start()

    def _write_all(self, fd, data):
        while data:
            n = os.write(fd, data)
            data = data[n:]

    def _emit(self, data):
        if not data:
            return
        try:
            self._write_all(self.job_fd, data)
        except OSError:
            pass
        if self.stage is not None:
            try:
                self.stage.write(data)
                self.stage.flush()
            except OSError:
                pass

    def _pump(self):
        carry = b""
        while True:
            try:
                data = os.read(self.r, 1 << 16)
            except OSError:
                data = b""
            if not data:
                self._emit(carry)
                break
            buf = carry + data
            carry = b""
            while True:
                k = buf.find(b"\x1e@chain_onejob_switch_")
                if k < 0:
                    # keep a possible partial marker at the end
                    tail = buf.rfind(b"\x1e")
                    if tail >= 0 and len(buf) - tail < 40 and b"\x1e" not in buf[tail + 1:]:
                        self._emit(buf[:tail])
                        carry = buf[tail:]
                    else:
                        self._emit(buf)
                    break
                end = buf.find(b"@\x1e", k + 1)
                if end < 0:
                    self._emit(buf[:k])
                    carry = buf[k:]
                    break
                self._emit(buf[:k])
                try:
                    n = int(buf[k + len(b"\x1e@chain_onejob_switch_"):end])
                except ValueError:
                    n = None
                if n == 0:                         # close(): stop reading (later writers get EPIPE)
                    self._emit(buf[end + 2:])
                    self.ended.set()
                    return
                self._switch_now(n)
                buf = buf[end + 2:]

    def _switch_now(self, n):
        with self.lock:
            path = self.pending.pop(n, None) if n is not None else None
            if self.stage is not None:
                try:
                    self.stage.close()
                except OSError:
                    pass
                self.stage = None
            if path is not None:
                self.stage = open(path, "ab")
            event = self.events.pop(n, None)
        if event is not None:
            event.set()

    def switch(self, path):
        """From now on the stage log is ``path`` (None: no stage log)."""
        sys.stdout.flush()
        sys.stderr.flush()
        with self.lock:
            self.counter += 1
            n = self.counter
            self.pending[n] = None if path is None else str(path)
            event = self.events[n] = threading.Event()
        os.write(1, self.MARK % n)
        if not event.wait(60.0):
            say("WARNING the log switch to %s was not confirmed within 60 s" % path)

    def close(self):
        """Back to the original fds after the pump has written everything up to the end marker."""
        self.switch(None)
        sys.stdout.flush()
        sys.stderr.flush()
        os.write(1, self.MARK % 0)
        self.ended.wait(30.0)
        os.dup2(self.job_fd, 1)
        os.dup2(self.err_fd, 2)
        self.thread.join(5.0)
        try:
            os.close(self.r)
        except OSError:
            pass


# ================================================================================================ stages
class Plan:
    def __init__(self, **kw):
        self.done = False
        self.rc = None
        self.argv = None
        self.env = {}
        self.out = None
        self.segment = None
        self.notes = []
        self.__dict__.update(kw)


def _npz_meta(path):
    with np.load(str(path), allow_pickle=False) as z:
        return json.loads(str(z["meta_json"]))


class PositiveFitStage:
    """jobs/positive_fit.sbatch with CHAIN=1 (cycle k of jobs/positive_cycles.sh)."""

    kind = "positive_fit"
    log_name = "fit.log"
    split_script = "positive_fit.sbatch"
    code_files = ("positive_fit.py", "positive_problem.py", "positive_activity.py", "positive_activity_check.py",
                  "grain_fit.py", "grain_crossing.py", "late_response_eval.py", "late_time_grid.py",
                  "test_positive_fit.py")
    run_copy_files = ("positive_fit.py", "positive_problem.py", "positive_activity.py", "grain_fit.py",
                      "grain_crossing.py")
    summary_regex = re.compile(r"^(positive_fit|consistency|initial cost|prior cost|lambda |gpu_model: bound|gpu_fit: )")
    env_regex = re.compile(r"^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|LAMBDA|NFEV|CHECK|CHECK_HELD_SPAN|RATIO|CYCLE|"
                           r"TESTNAME|START_RESTART|START_SET|GRAIN_[A-Z_]*|HELD|CHAIN|MODEL_EXT|OMP_[A-Z_]*|"
                           r"OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|"
                           r"CUDA_VISIBLE_DEVICES)=")

    def __init__(self, chain, k):
        self.chain, self.k = chain, int(k)
        self.name = "c%dfit" % self.k
        self.directory = chain.root / ("c%d" % self.k) / "fit"
        a = chain.args
        self.start_restart = (Path(a.source) if self.k == 1
                              else chain.root / ("c%d" % (self.k - 1)) / "fit" / "final" / "fit_restart.npz")
        self.held = None if self.k == 1 else chain.root / ("c%d" % (self.k - 1)) / "ft" / "finetune_state.npz"
        self.grain_from = self.held
        self.start_set = chain.c1_start_set if self.k == 1 else ""

    # --- positive_fit.sbatch state()
    def state(self, seg=None):
        D = self.directory
        segs = sorted(int(p.name[3:]) for p in D.glob("seg[0-9][0-9]")) if D.is_dir() else []
        k = int(seg) if seg else (segs[-1] if segs else 0)
        if k == 0:
            return 0, "none", 0, "-", "-"
        rs = D / ("seg%02d" % k) / "fit_restart.npz"
        if not rs.is_file():
            return k, "norestart", 0, "-", "-"
        meta = _npz_meta(rs)
        used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
        chain = (meta.get("model_extension") or {}).get("fit_chain") or {}
        start = (chain.get("start_restart") or {}).get("path") or "-"
        sset = ",".join("%s=%r" % (e["name"], float(e["value"])) for e in chain.get("start_set") or []) or "-"
        return k, meta.get("status", "unknown"), used, start, sset

    def canonical_start_set(self):
        if not self.start_set:
            return None
        return ",".join("%s=%r" % (n.strip(), float(v)) for n, v in
                        (i.split("=") for i in self.start_set.split(":") if i))

    def final_link(self, k):
        seg = "seg%02d" % k
        final = self.directory / "final"
        if final.is_symlink():
            if os.readlink(str(final)) != seg:
                print("%s -> %s, expected %s" % (final, os.readlink(str(final)), seg), file=sys.stderr, flush=True)
                return 3
        elif final.exists():
            print("%s exists and is not a link" % final, file=sys.stderr, flush=True)
            return 3
        else:
            os.symlink(seg, str(final))
        print("%s %s complete: %s -> %s" % (utc_now(), self.label(), final, seg), flush=True)
        return 0

    def label(self):
        return "positive_r%s_lam%s_c%d" % (self.chain.args.ratio, self.chain.args.lambda_, self.k)

    def check_held_identity(self, k):
        """Completed fit k > 1: the held term recorded in its restart is the previous fine tuning's state file."""
        if self.held is None:
            return None
        meta = _npz_meta(self.directory / ("seg%02d" % k) / "fit_restart.npz")
        rec = ((meta.get("model_extension") or {}).get("held_exponent") or {})
        want = sha256_file(readlink_f(self.held)) if self.held.exists() else None
        got = rec.get("sha256") if isinstance(rec, dict) else None
        if got is None or want is None or got != want:
            return "held term of %s/seg%02d (%s) is not %s (%s)" % (self.directory, k, got, self.held, want)
        return None

    def prepare(self):
        c, a, D = self.chain, self.chain.args, self.directory
        need = [c.base / "COPS_model.dat", c.base / "elements.bin", c.base / "kernels", c.point_table, c.bundle,
                HERE / "positive_fit.py", HERE / "positive_problem.py", HERE / "positive_activity.py",
                HERE / "kernels", self.start_restart] + ([self.held] if self.held else [])
        for f in need:
            if not Path(f).exists():
                return Plan(rc=2, notes=["missing: %s" % f])
        D.mkdir(parents=True, exist_ok=True)
        LAST, LSTATE, LUSED, LSTART, LSSET = self.state()
        moved = None
        if LSTATE == "norestart":
            seg_dir = D / ("seg%02d" % LAST)
            last = c.events.last(seg_dir)
            if not (last and last.get("event") == "end" and last.get("rc") == EXIT_WALL_TIME
                    and not last.get("restart_written")):
                return Plan(rc=3, notes=["%s has no restart file (a failed or killed attempt): chain stopped (CHAIN=1 "
                                         "rule of the split chain); inspect it" % seg_dir])
            moved = D / ("seg%02d.incomplete_%s" % (LAST, utc_stamp()))
            os.rename(str(seg_dir), str(moved))           # never deleted
            LAST -= 1
            if LAST >= 1:
                LAST, LSTATE, LUSED, LSTART, LSSET = self.state(LAST)
            else:
                LSTATE, LUSED, LSTART, LSSET = "none", 0, "-", "-"
        if LSTATE != "none" and LSTART != str(readlink_f(self.start_restart)):
            return Plan(rc=3, notes=["START_RESTART=%s but %s started from %s: not the same fit"
                                     % (self.start_restart, D, LSTART)])
        sset = self.canonical_start_set()
        if sset is not None and LSTATE != "none" and LSSET != sset:
            return Plan(rc=3, notes=["START_SET %s but %s started with --start-set %s: not the same fit"
                                     % (sset, D, LSSET)])
        if LSTATE == "none":
            SEG, LEFT = 1, a.nfev
        elif LSTATE == "interrupted":
            SEG, LEFT = LAST + 1, a.nfev - LUSED
        elif LSTATE == "complete":
            print("%s %s: segment %d is complete (%d evaluations); nothing to do" % (utc_now(), self.label(), LAST,
                                                                                    LUSED), flush=True)
            return self._done(LAST)
        else:
            return Plan(rc=3, notes=["%s/seg%02d/fit_restart.npz has status '%s'; not continued" % (D, LAST, LSTATE)])
        if LEFT <= 0:
            print("%s %s: %d of %d evaluations used; nothing to do" % (utc_now(), self.label(), LUSED, a.nfev),
                  flush=True)
            return self._done(LAST)
        if SEG == 1:
            START = ["--start-restart", str(readlink_f(self.start_restart))]
            if self.start_set:
                START += ["--start-set", self.start_set.replace(":", ",")]
        else:
            START = ["--restart-from", str(D / ("seg%02d" % LAST) / "fit_restart.npz")]
        GRAIN = []
        if self.grain_from is not None:
            GRAIN += ["--grain-start-from", str(readlink_f(self.grain_from))]
        HOLD = [] if self.held is None else ["--held-exponent", str(readlink_f(self.held))]
        O = D / ("seg%02d" % SEG)
        if O.exists():
            return Plan(rc=2, notes=["exists: %s" % O])
        CHK = ["--check-jacobian"] if SEG == 1 and a.check == 1 else []
        argv = ([str(HERE / "positive_fit.py"), "--point-table", str(c.point_table), "--ratio", a.ratio] + START + GRAIN
                + HOLD + ["--crossing-weight", "10"] + CHK
                + ["--devices", c.devices, "--precision", "mixed", "--max-nfev", str(LEFT), "--segment", str(SEG),
                   "--wall-time", "slurm", "--template-start", "0", "--", str(c.base),
                   "--profile", "all_available", "--degree", "4", "--enhancement-degree", "4", "--lambda", a.lambda_,
                   "--jobs", "1", "--fix-temperature", "--tbounds", "30", "400", "--ftol", "1e-10", "--xtol", "1e-10",
                   "--gtol", "1e-10", "--check-tolerance", "inf", "--out", str(O),
                   "--activity-template", str(c.bundle), "cliff_ge_60deg"])
        env = {"LAMBDA": a.lambda_, "RATIO": a.ratio, "CYCLE": str(self.k), "NFEV": str(a.nfev), "CHECK": str(a.check),
               "CHAIN": "1", "MODEL_EXT": str(c.model_ext), "START_RESTART": str(self.start_restart)}
        if self.held is not None:
            env.update(HELD=str(self.held), GRAIN_FROM=str(self.grain_from))
        if self.start_set:
            env["START_SET"] = self.start_set
        lines = ["fit %s, directory %s; ratio %s; cycle %d; start restart %s; grain A - B - from %s; held %s"
                 % (self.label(), D, a.ratio, self.k, self.start_restart, self.grain_from or "-", self.held or "-"),
                 "CO2 flux floor: %s" % ("on (COPS_POS_FLOOR_CO2=1, ratio %s)" % a.ratio if a.ratio == "1e-5"
                                         else "off (COPS_POS_FLOOR_CO2 unset)")]
        if self.start_set:
            lines.append("start set: %s (segment 1: --start-set; the fit's chain: %s)" % (self.start_set, LSSET))
        lines.append("model_ext copied and run: %s (MODEL_EXT) -> %s (the job's code copy; this process)"
                     % (c.model_ext, HERE))
        if self.held is not None:
            lines.append("held file sha256 %s" % sha256_file(readlink_f(self.held)))
        lines.append("segment %d; previous segment %d (%s, %d evaluations used); evaluations for this segment %d of %d"
                     % (SEG, LAST, LSTATE, LUSED, LEFT, a.nfev))
        if moved:
            lines.append("moved aside: %s (the last attempt ended with exit 75 before writing a restart file)" % moved)
        return Plan(argv=argv, env=env, out=O, segment=SEG, notes=lines,
                    commands_env=self._commands_env(), summary="segment %d, %d evaluations" % (SEG, LEFT))

    def _commands_env(self):
        r = self.chain.args.ratio
        return ["export COPS_POS_ENABLE=1 COPS_POS_NIGHT_RATIO=%s COPS_POS_ACTIVITY=exp COPS_POS_FLOOR_RATIO=%s" % (r, r),
                "export COPS_POS_FLOOR_CO2=1" if r == "1e-5" else "unset COPS_POS_FLOOR_CO2"]

    def _done(self, last):
        msg = self.check_held_identity(last)
        if msg:
            return Plan(rc=3, notes=[msg])
        rc = self.final_link(last)
        return Plan(done=(rc == 0), rc=(None if rc == 0 else rc), segment=last)

    def call(self, plan):
        sys.argv = list(plan.argv)
        return call_main(lambda: importlib.import_module("positive_fit").main())   # import errors -> exit code too

    def finish(self, plan, rc):
        if rc == 0:
            _, st, _, _, _ = self.state(plan.segment)
            if st == "complete":
                rc = self.final_link(plan.segment)
        if rc == EXIT_WALL_TIME:
            print("stopped before the wall-time limit: resubmit the same command to continue with segment %d"
                  % (plan.segment + 1), flush=True)
        return rc

    def restart_written(self, plan):
        return (plan.out / "fit_restart.npz").is_file()

    def extra_provenance(self, plan, P):
        pass


class PositiveFineTuneStage:
    """jobs/positive_finetune.sbatch (cycle k of jobs/positive_cycles.sh: FIT=c<k>/fit/final, OUT=c<k>/ft)."""

    kind = "positive_finetune"
    log_name = "finetune.log"
    split_script = "positive_finetune.sbatch"
    code_files = ("positive_finetune.py", "positive_problem.py", "positive_activity.py", "grain_crossing.py",
                  "test_positive_finetune.py")
    run_copy_files = code_files
    PRODUCTS = ("finetune_state.npz", "finetune_result.txt", "COPS_model.dat", "late_window_model.dat",
                "pressure_comparison.txt", "fine_tuning.json", "validation.json")
    ENV_KEYS = ("COPS_POS_ENABLE", "COPS_POS_NIGHT_RATIO", "COPS_POS_ACTIVITY", "COPS_POS_FLOOR_RATIO",
                "COPS_POS_FLOOR_CO2")
    summary_regex = re.compile(r"^positive_finetune: (bound|parameters|FACETS|objective|the fit|check|timing|validation|"
                               r"start|final|products|checkpoint|stopped|complete|ERROR)")
    env_regex = re.compile(r"^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|FIT|OUT|MAXNFEV|FACETS|WINDOW_MIN|LSMR_MAXITER|"
                           r"FLOOR_RATIO|VALIDATE[A-Z_]*|RESUME|CHECKPOINT_MIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|"
                           r"MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=")

    def __init__(self, chain, k):
        self.chain, self.k = chain, int(k)
        self.name = "c%dft" % self.k
        self.fit_given = chain.root / ("c%d" % self.k) / "fit" / "final"
        self.directory = chain.root / ("c%d" % self.k) / "ft"

    def segstate(self, nn):
        p = self.directory / ("seg%02d" % int(nn)) / "finetune_state.npz"
        if not p.is_file():
            return "nostate", 0, "-"
        meta = _npz_meta(p)
        return meta.get("stop", "unknown"), int(meta.get("nfev_total", meta.get("nfev", 0))), \
            meta.get("fit_restart_sha256", "-")

    def fit_environment(self, fit):
        m = json.loads((fit / "fit_manifest.json").read_text())
        env = m["model_extension"]["positive"]["environment"]
        out = {}
        for k in self.ENV_KEYS:
            if k == "COPS_POS_FLOOR_CO2" and env.get(k) is None:
                continue
            v = str(env[k])
            assert v.replace(".", "", 1).replace("e-", "", 1).isalnum(), v
            out[k] = v
        return out

    def product_check(self, seg):
        p = self.directory / seg / "finetune_state.npz"
        try:
            with np.load(str(p), allow_pickle=False) as z:
                meta = json.loads(str(z["meta_json"]))
                h, d, hp = z["held_next"], z["delta"], z["held_prev"]
                A, B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        except Exception as exc:
            print("product check failed: %s (%s)" % (p, exc), file=sys.stderr, flush=True)
            return 4
        ok = (meta.get("complete") is True and meta.get("stop") in ("tolerance", "max_nfev")
              and h.shape == d.shape == hp.shape and bool(np.all(np.isfinite(h))) and np.array_equal(h, hp + d)
              and np.isfinite(A) and np.isfinite(B) and A >= 0 and B >= 0)
        print("finetune_state.npz: stop %s, complete %s, %d evaluations, %d facets, grain A %.6g uPa B %.6g uPa/100s: %s"
              % (meta.get("stop"), meta.get("complete"), meta.get("nfev_total", -1), h.size, A, B,
                 "OK" if ok else "FAILED"), flush=True)
        if not ok:
            print("product check failed: %s" % p, file=sys.stderr, flush=True)
            return 4
        return 0

    def final_link(self, nn):
        seg = "seg%02d" % int(nn)
        OUT = self.directory
        rc = self.product_check(seg)
        if rc:
            return rc
        final = OUT / "final"
        if final.is_symlink():
            if os.readlink(str(final)) != seg:
                print("%s -> %s, expected %s" % (final, os.readlink(str(final)), seg), file=sys.stderr, flush=True)
                return 3
        elif final.exists():
            print("%s exists and is not a link" % final, file=sys.stderr, flush=True)
            return 3
        else:
            os.symlink(seg, str(final))
        for f in self.PRODUCTS:
            if not (OUT / seg / f).exists():
                continue
            link = OUT / f
            if link.is_symlink():
                if os.readlink(str(link)) != "final/" + f:
                    print("%s is not final/%s" % (link, f), file=sys.stderr, flush=True)
                    return 3
            elif link.exists():
                print("%s exists and is not a link" % link, file=sys.stderr, flush=True)
                return 3
            else:
                os.symlink("final/" + f, str(link))
        print("%s fine tuning complete: %s -> %s" % (utc_now(), final, seg), flush=True)
        return 0

    def prepare(self):
        c, a = self.chain, self.chain.args
        if not self.fit_given.exists():
            return Plan(rc=2, notes=["missing: %s (the fit of cycle %d is not complete)" % (self.fit_given, self.k)])
        FIT = readlink_f(self.fit_given)
        OUT = Path(os.path.abspath(str(self.directory)))
        facetopts = (["--facets", "late_only", "--window-min", a.ft_window_min] if a.ft_facets == "late_only"
                     else ["--facets", "all_seen"])
        need = [FIT / "fit_restart.npz", FIT / "fit_manifest.json", FIT / "activity_exp_facets.npz",
                HERE / "kernels" / "finetune.cu", HERE / "kernels" / "pressure.cu", HERE / "kernels" / "activity.cu"] + \
            [HERE / f for f in ("positive_finetune.py", "positive_problem.py", "positive_activity.py",
                                "grain_crossing.py", "test_positive_finetune.py")]
        for f in need:
            if not f.exists():
                return Plan(rc=2, notes=["missing: %s" % f])
        try:
            envfit = self.fit_environment(FIT)
        except Exception as exc:
            return Plan(rc=2, notes=["cannot read the COPS_POS_* environment of %s: %s" % (FIT, exc)])
        mine = {k: os.environ.get(k) for k in self.ENV_KEYS if os.environ.get(k) is not None}
        mine.pop("COPS_POS_FROZEN_CONFIG", None)
        if mine != envfit:
            return Plan(rc=2, notes=["the fit's COPS_POS_* %s differ from this chain process's %s (one process, one "
                                     "environment): not run" % (envfit, mine)])
        floor_ratio = a.ft_floor_ratio or envfit["COPS_POS_FLOOR_RATIO"]
        FITSHA = sha256_file(FIT / "fit_restart.npz")
        LAST, LSTOP, LUSED, moved = 0, "none", 0, None
        if OUT.exists():
            if not OUT.is_dir():
                return Plan(rc=2, notes=["OUT is not a directory: %s" % OUT])
            for e in OUT.iterdir():
                b = e.name
                if (re.fullmatch(r"seg[0-9][0-9]", b) or re.fullmatch(r"seg[0-9][0-9]\.incomplete_.*", b)
                        or b.startswith("provenance-") or b in ("validate", "final") or b in self.PRODUCTS):
                    continue
                return Plan(rc=2, notes=["OUT %s holds %s: not a directory of this job script (choose a new OUT)"
                                         % (OUT, b)])
            segs = [int(p.name[3:]) for p in OUT.glob("seg[0-9][0-9]")]
            LAST = max(segs) if segs else 0
        LSHA = "-"
        if LAST >= 1:
            LSTOP, LUSED, LSHA = self.segstate(LAST)
            if LSTOP == "nostate":                          # killed or failed: move aside, run that segment again
                moved = OUT / ("seg%02d.incomplete_%s" % (LAST, utc_stamp()))
                os.rename(str(OUT / ("seg%02d" % LAST)), str(moved))
                LAST, LSTOP, LUSED = LAST - 1, "none", 0
                if LAST >= 1:
                    LSTOP, LUSED, LSHA = self.segstate(LAST)
            if LAST >= 1 and LSHA != FITSHA:
                return Plan(rc=3, notes=["%s/seg%02d tunes the fit restart %s, FIT has %s: not the same fine tuning"
                                         % (OUT, LAST, LSHA, FITSHA)])
        if LSTOP == "none":
            SEG, LEFT = 1, a.ft_maxnfev
        elif LSTOP == "wall_time":
            SEG, LEFT = LAST + 1, a.ft_maxnfev - LUSED
        elif LSTOP in ("tolerance", "max_nfev"):
            print("%s %s: segment %d finished (%s, %d evaluations); nothing to do" % (utc_now(), OUT, LAST, LSTOP,
                                                                                   LUSED), flush=True)
            rc = self.final_link(LAST)
            return Plan(done=(rc == 0), rc=(None if rc == 0 else rc), segment=LAST)
        else:
            return Plan(rc=3, notes=["%s/seg%02d/finetune_state.npz has stop '%s'; not continued" % (OUT, LAST, LSTOP)])
        if LEFT <= 0:
            return Plan(rc=3, notes=["%s: %d of %d evaluations used and the last segment stopped at the wall time; "
                                     "nothing is left to run: the budget is used up but the state is marked incomplete "
                                     "(decide: raise MAXNFEV)" % (OUT, LUSED, a.ft_maxnfev)])
        O = OUT / ("seg%02d" % SEG)
        if O.exists():
            return Plan(rc=2, notes=["exists: %s" % O])
        OUT.mkdir(parents=True, exist_ok=True)
        opts = ["--max-nfev", str(LEFT), "--segment", str(SEG), "--nfev-before", str(LUSED)]
        if SEG == 1 and a.ft_validate == 1:
            opts += ["--validate"]
        if SEG > 1:
            opts += ["--resume-from", str(OUT / ("seg%02d" % LAST) / "finetune_state.npz")]
        argv = ([str(HERE / "positive_finetune.py"), "--fit", str(FIT), "--out", str(O), "--floor-ratio",
                 floor_ratio, "--devices", c.devices, "--lsmr-maxiter", str(a.ft_lsmr_maxiter), "--wall-time", "slurm",
                 "--checkpoint-minutes", str(a.ft_checkpoint_min)] + facetopts + opts)
        env = {"FIT": str(self.fit_given), "OUT": str(OUT), "MAXNFEV": str(a.ft_maxnfev), "FACETS": a.ft_facets,
               "MODEL_EXT": str(c.model_ext), "LSMR_MAXITER": str(a.ft_lsmr_maxiter),
               "CHECKPOINT_MIN": str(a.ft_checkpoint_min), "VALIDATE": str(a.ft_validate)}
        if a.ft_facets == "late_only":
            env["WINDOW_MIN"] = a.ft_window_min
        envlines = " ".join("export %s=%s" % kv for kv in envfit.items())
        lines = ["fine tuning of %s (given as %s; restart sha256 %s) -> %s" % (FIT, self.fit_given, FITSHA, O),
                 "model_ext copied and run: %s (MODEL_EXT) -> %s (the job's code copy; this process)"
                 % (c.model_ext, HERE),
                 "env: %s ; floor ratio %s; LSMR_MAXITER %s; FACETS %s WINDOW_MIN %s"
                 % (envlines, floor_ratio, a.ft_lsmr_maxiter, a.ft_facets,
                    a.ft_window_min if a.ft_facets == "late_only" else "-"),
                 "segment %d; previous segment %d (%s, %d evaluations used); evaluations for this segment %d of %d"
                 % (SEG, LAST, LSTOP, LUSED, LEFT, a.ft_maxnfev)]
        if moved:
            lines.append("moved aside: %s" % moved)
        commands_env = ["# the tuned fit's COPS_POS_* (fit_manifest.json model_extension.positive.environment)",
                        "unset COPS_POS_FLOOR_CO2"] + ["export %s=%s" % kv for kv in envfit.items()]
        return Plan(argv=argv, env=env, out=O, segment=SEG, notes=lines, fit=FIT, commands_env=commands_env,
                    summary="segment %d, %d evaluations" % (SEG, LEFT))

    def call(self, plan):
        sys.argv = list(plan.argv)

        def run():
            module = importlib.import_module("positive_finetune")
            module._process_setup()                 # as when it is the main module (idempotent install, handlers)
            self.chain.install_finetune_engine(module, self)
            return module.main(list(plan.argv[1:]))
        return call_main(run)                       # import errors -> exit code too

    def finish(self, plan, rc):
        if rc == 0:
            st, _, _ = self.segstate(plan.segment)
            if st in ("tolerance", "max_nfev"):
                rc = self.final_link(plan.segment)
            else:
                print("segment %d exited 0 with stop '%s'" % (plan.segment, st), file=sys.stderr, flush=True)
                rc = 3
        if rc == EXIT_WALL_TIME:
            print("stopped before the wall-time limit: resubmit the same command to continue with segment %d"
                  % (plan.segment + 1), flush=True)
        return rc

    def restart_written(self, plan):
        return (plan.out / "finetune_state.npz").is_file()

    def extra_provenance(self, plan, P):
        run_to_file(["sha256sum", str(plan.fit / "fit_restart.npz"), str(plan.fit / "activity_exp_facets.npz"),
                     str(plan.fit / "fit_manifest.json")], P / "fit.sha256")


# ------------------------------------------------------------------------------------------------ slope kind (v6)
# The slope-multiplier family (user 2026-11-10/11: bumps (option 3a), softplus, squared harmonics), the split chain being
# jobs/slope_cycles.sh -> jobs/slope_mult_fit.sbatch (CHAIN=1) and jobs/slope_finetune.sbatch. Model options of the fit
# (inserted where slope_mult_fit.sbatch puts its MOPT, before --held): one entry per model.
SLOPE_MODELS = {
    "bumps": {"fit_args": lambda c: ["--bump-bundle", str(readlink_f(c.args.bumps))], "c1_start": "nnls"},
}
SLOPE_FIT_FIX = ["--fix", "Tmax=182.09922", "--fix", "T75=136.11149", "--fix", "Tmin=133.0", "--fix", "c_slope_00_40=1"]
SLOPE_ENV_UNSET = "for v in $(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)'); do unset $v; done"


class SlopeFitStage(PositiveFitStage):
    """jobs/slope_mult_fit.sbatch with CHAIN=1 and FITDIR=<chain>/c<k>/fit (cycle k of jobs/slope_cycles.sh)."""

    kind = "slope_fit"
    split_script = "slope_mult_fit.sbatch"
    code_files = ("slope_mult_fit.py", "slope_mult.py", "slope_mult_problem.py", "slope_bins.py", "slope_bumps.py",
                  "slope_finetune.py", "test_slope_mult.py", "positive_finetune.py", "clip_problem.py",
                  "clip_refit.py", "grain_fit_fixed.py", "fixed_parameters.py", "grain_fit.py", "grain_crossing.py",
                  "positive_problem.py", "late_response_eval.py", "test_clip_refit.py", "test_positive_fit.py")
    run_copy_files = code_files
    summary_regex = re.compile(r"^(slope_mult|clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |"
                               r"grain_fit: |fit data|selection profile)")
    env_regex = re.compile(r"^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|NFEV|CHECK|CHECKPARAMS|TESTNAME|TESTMODE|START|"
                           r"START_RESTART|L_THREADS|L_MODE|BUMPS|HELD|GRAIN_FROM|FITDIR|BUNDLE|CHAIN|MODEL_EXT|"
                           r"SLOPE_MODEL|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|"
                           r"SLURM_JOB_GPUS|SLURM_JOB_NODELIST|SLURM_JOB_NUM_NODES|CUDA_VISIBLE_DEVICES)=")

    def __init__(self, chain, k):
        self.chain, self.k = chain, int(k)
        self.name = "c%dfit" % self.k
        self.directory = chain.root / ("c%d" % self.k) / "fit"
        a = chain.args
        self.model = SLOPE_MODELS[a.slope_model]
        if self.k == 1:
            self.start = a.slope_start or self.model["c1_start"]
            self.start_restart = Path(a.source)
            self.held = self.grain_from = None
        else:
            prev = chain.root / ("c%d" % (self.k - 1))
            self.start = "chain"
            self.start_restart = prev / "fit" / "final" / "fit_restart.npz"
            self.held = self.grain_from = prev / "ft" / "finetune_state.npz"
        self.start_set = ""

    def label(self):
        return "slope_fit_%s_%s" % (self.directory.parent.name, self.directory.name)

    def state(self, seg=None):
        """slope_mult_fit.sbatch state(): (k, status, used, start restart path, start choice)."""
        D = self.directory
        segs = sorted(int(p.name[3:]) for p in D.glob("seg[0-9][0-9]")) if D.is_dir() else []
        k = int(seg) if seg else (segs[-1] if segs else 0)
        if k == 0:
            return 0, "none", 0, "-", "-"
        rs = D / ("seg%02d" % k) / "fit_restart.npz"
        if not rs.is_file():
            return k, "norestart", 0, "-", "-"
        meta = _npz_meta(rs)
        used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
        ext = meta.get("model_extension") or {}
        chain = ext.get("fit_chain") or {}
        return (k, meta.get("status", "unknown"), used, (chain.get("start_restart") or {}).get("path") or "-",
                ext.get("start_choice") or chain.get("start_choice") or "-")

    def check_held_identity(self, k):
        return None              # the slope restart's held record is agent 1's; the fit's --held is on every segment

    def config_argv(self):
        """The cycle-1 fit's arguments without segment-dependent parts (the process-level GPU configuration)."""
        return self._argv(1, self.chain.args.nfev, self._start_options(1, None), self.directory / "seg01")[1:]

    def _start_options(self, SEG, last_restart):
        """slope_mult_fit.sbatch SOPT: segment 1 --start S --start-restart R (+ --grain-from for chain); later
        segments --start S --restart-from <the last segment's restart>."""
        if SEG == 1:
            opts = ["--start", self.start, "--start-restart", str(readlink_f(self.start_restart))]
            if self.start == "chain":
                opts += ["--grain-from", str(readlink_f(self.grain_from))]
            return opts
        return ["--start", self.start, "--restart-from", str(last_restart)]

    def _argv(self, SEG, LEFT, SOPT, O):
        a, c = self.chain.args, self.chain
        MOPT = list(self.model["fit_args"](c)) + (["--held", str(readlink_f(self.held))] if self.held else [])
        BOPT = ["--l-threads", str(a.l_threads), "--l-mode", a.l_mode] if a.l_threads else []
        CHK, GCHK, CHKP = [], [], []
        if SEG == 1 and a.check == 1:
            CHK, GCHK = ["--check-full-jacobian"], ["--check-jacobian"]
            if a.checkparams:
                CHKP = ["--check-parameters", a.checkparams.replace(":", ",")]
        return ([str(HERE / "slope_mult_fit.py"), "--point-table", str(c.point_table), "--slope-bundle",
                 str(readlink_f(a.slope_bundle))] + MOPT + SOPT + BOPT + CHK + CHKP + SLOPE_FIT_FIX
                + ["--crossing-weight", "10"] + GCHK
                + ["--devices", c.devices, "--precision", "mixed", "--max-nfev", str(LEFT), "--segment", str(SEG),
                   "--wall-time", "slurm", "--template-start", "0", "--", str(c.base),
                   "--profile", "all_available", "--degree", "4", "--enhancement-degree", "4", "--lambda", "0",
                   "--jobs", "1", "--tbounds", "30", "400", "--ftol", "1e-10", "--xtol", "1e-10", "--gtol", "1e-10",
                   "--check-tolerance", "1e-6", "--out", str(O)])

    def prepare(self):
        c, a, D = self.chain, self.chain.args, self.directory
        START_RESTART = str(readlink_f(self.start_restart))
        need = [c.base / "COPS_model.dat", c.base / "elements.bin", c.base / "kernels", c.point_table,
                Path(START_RESTART), readlink_f(a.slope_bundle)] + \
            ([readlink_f(a.bumps)] if a.slope_model == "bumps" else []) + \
            ([self.held, self.grain_from] if self.held else []) + [HERE / f for f in self.code_files]
        for f in need:
            if not Path(f).exists():
                return Plan(rc=2, notes=["missing: %s" % f])
        D.mkdir(parents=True, exist_ok=True)
        LAST, LSTATE, LUSED, LSTART, LCHOICE = self.state()
        moved = None
        if LSTATE == "norestart":
            seg_dir = D / ("seg%02d" % LAST)
            last = c.events.last(seg_dir)
            if not (last and last.get("event") == "end" and last.get("rc") == EXIT_WALL_TIME
                    and not last.get("restart_written")):
                return Plan(rc=3, notes=["%s has no restart file (a failed attempt): chain stopped (CHAIN=1); inspect "
                                         "it" % seg_dir])
            moved = D / ("seg%02d.incomplete_%s" % (LAST, utc_stamp()))
            os.rename(str(seg_dir), str(moved))
            LAST -= 1
            if LAST >= 1:
                LAST, LSTATE, LUSED, LSTART, LCHOICE = self.state(LAST)
            else:
                LSTATE, LUSED, LSTART, LCHOICE = "none", 0, "-", "-"
        if LSTATE != "none" and LSTART != START_RESTART:
            return Plan(rc=3, notes=["START_RESTART=%s but %s started from %s: not the same fit"
                                     % (START_RESTART, D, LSTART)])
        if LSTATE != "none" and LCHOICE != self.start:
            return Plan(rc=3, notes=["START=%s but %s started with --start %s: not the same fit"
                                     % (self.start, D, LCHOICE)])
        if LSTATE == "none":
            SEG, LEFT = 1, a.nfev
        elif LSTATE == "interrupted":
            SEG, LEFT = LAST + 1, a.nfev - LUSED
        elif LSTATE == "complete":
            print("%s %s: segment %d is complete (%d evaluations); nothing to do" % (utc_now(), self.label(), LAST,
                                                                                    LUSED), flush=True)
            return self._done(LAST)
        else:
            return Plan(rc=3, notes=["%s/seg%02d/fit_restart.npz has status '%s'; not continued" % (D, LAST, LSTATE)])
        if LEFT <= 0:
            print("%s %s: %d of %d evaluations used; nothing to do" % (utc_now(), self.label(), LUSED, a.nfev),
                  flush=True)
            return self._done(LAST)
        SOPT = self._start_options(SEG, D / ("seg%02d" % LAST) / "fit_restart.npz")
        O = D / ("seg%02d" % SEG)
        if O.exists():
            return Plan(rc=2, notes=["exists: %s" % O])
        argv = self._argv(SEG, LEFT, SOPT, O)
        env = {"NFEV": str(a.nfev), "CHECK": str(a.check), "START": self.start, "START_RESTART": START_RESTART,
               "BUNDLE": str(readlink_f(a.slope_bundle)), "FITDIR": str(D), "CHAIN": "1",
               "MODEL_EXT": str(c.model_ext), "SLOPE_MODEL": a.slope_model}
        if a.checkparams:
            env["CHECKPARAMS"] = a.checkparams
        if a.l_threads:
            env.update(L_THREADS=str(a.l_threads), L_MODE=a.l_mode)
        if a.slope_model == "bumps":
            env["BUMPS"] = str(readlink_f(a.bumps))
        if self.held:
            env.update(HELD=str(readlink_f(self.held)), GRAIN_FROM=str(readlink_f(self.grain_from)))
        lines = ["fit %s, directory %s; model %s; start %s; start restart (segment 1) %s; slope bundle %s"
                 % (self.label(), D, a.slope_model, self.start, START_RESTART, readlink_f(a.slope_bundle)),
                 "model options %s; held %s; grain from (segment 1) %s"
                 % (" ".join(self.model["fit_args"](c)) or "-", self.held or "none", self.grain_from or
                    "the start restart"),
                 "L = Y a of the fitter process: %s" % ("%s, %d threads" % (a.l_mode, a.l_threads) if a.l_threads
                                                        else "one call, 1 thread"),
                 "model_ext copied and run: %s (MODEL_EXT) -> %s (the job's code copy; this process)"
                 % (c.model_ext, HERE),
                 "segment %d; previous segment %d (%s, %d evaluations used); evaluations for this segment %d of %d"
                 % (SEG, LAST, LSTATE, LUSED, LEFT, a.nfev)]
        if moved:
            lines.append("moved aside: %s (the last attempt ended with exit 75 before writing a restart file)" % moved)
        return Plan(argv=argv, env=env, out=O, segment=SEG, notes=lines, commands_env=[SLOPE_ENV_UNSET],
                    summary="segment %d, %d evaluations" % (SEG, LEFT))

    def call(self, plan):
        sys.argv = list(plan.argv)
        return call_main(lambda: importlib.import_module("slope_mult_fit").main())


class SlopeFineTuneStage(PositiveFineTuneStage):
    """jobs/slope_finetune.sbatch (cycle k of jobs/slope_cycles.sh: FIT=c<k>/fit/final, OUT=c<k>/ft)."""

    kind = "slope_finetune"
    split_script = "slope_finetune.sbatch"
    code_files = ("slope_finetune.py", "slope_mult.py", "slope_mult_problem.py", "positive_finetune.py",
                  "positive_problem.py", "grain_crossing.py", "clip_problem.py", "test_slope_mult.py")
    run_copy_files = code_files
    fit_files = ("fit_restart.npz", "fit_manifest.json", "slope_mult_facets.npz")
    summary_regex = re.compile(r"^(slope_finetune|slope_mult: (GPU patch installed .*role driver|L = Y a)|Traceback|"
                               r".*Error)")
    env_regex = re.compile(r"^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|FIT|OUT|FACETS|WINDOW_MIN|MAXNFEV|LSMR_MAXITER|"
                           r"L_THREADS|VALIDATE|VALIDATE_ONLY|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PATH|"
                           r"LD_LIBRARY_PATH|SLURM_JOB_NODELIST)=")

    def segstate(self, nn):
        p = self.directory / ("seg%02d" % int(nn)) / "finetune_state.npz"
        if not p.is_file():
            return "nostate", 0, "-"
        meta = _npz_meta(p)
        return meta.get("stop", "unknown"), int(meta.get("nfev_total", 0)), \
            (meta.get("fit") or {}).get("fit_restart_sha256", "-")

    def product_check(self, seg):
        p = self.directory / seg / "finetune_state.npz"
        try:
            with np.load(str(p), allow_pickle=False) as z:
                meta = json.loads(str(z["meta_json"]))
                h, d, hp = z["held_next"], z["delta"], z["held_prev"]
                A, B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        except Exception as exc:
            print("product check failed: %s (%s)" % (p, exc), file=sys.stderr, flush=True)
            return 4
        ok = (meta.get("complete") is True and meta.get("stop") in ("tolerance", "max_nfev")
              and h.shape == d.shape == hp.shape and bool(np.all(np.isfinite(h))) and np.array_equal(h, hp + d)
              and np.isfinite(A) and np.isfinite(B) and A >= 0 and B >= 0)
        print("finetune_state.npz: stop %s, complete %s, %d evaluations, %d parameters, grain A %.6g uPa B %.6g: %s" % (
            meta.get("stop"), meta.get("complete"), meta.get("nfev_total", -1), meta.get("parameters", -1), A, B,
            "OK" if ok else "FAILED"), flush=True)
        if not ok:
            print("product check failed: %s" % p, file=sys.stderr, flush=True)
            return 4
        return 0

    def prepare(self):
        c, a = self.chain, self.chain.args
        if not self.fit_given.exists():
            return Plan(rc=2, notes=["missing: %s (the fit of cycle %d is not complete)" % (self.fit_given, self.k)])
        FIT = readlink_f(self.fit_given)
        OUT = Path(os.path.abspath(str(self.directory)))
        for f in [FIT / n for n in self.fit_files] + [HERE / f for f in self.code_files]:
            if not f.exists():
                return Plan(rc=2, notes=["missing: %s" % f])
        FITSHA = sha256_file(FIT / "fit_restart.npz")
        LAST, LSTOP, LUSED, LSHA, moved = 0, "none", 0, "-", None
        if OUT.exists():
            for e in OUT.iterdir():
                b = e.name
                if (re.fullmatch(r"seg[0-9][0-9]", b) or re.fullmatch(r"seg[0-9][0-9]\.incomplete_.*", b)
                        or b.startswith("provenance-") or b in ("validate", "final") or b in self.PRODUCTS):
                    continue
                return Plan(rc=2, notes=["OUT %s holds %s: not a directory of this job script" % (OUT, b)])
            segs = [int(p.name[3:]) for p in OUT.glob("seg[0-9][0-9]")]
            LAST = max(segs) if segs else 0
        if LAST >= 1:
            LSTOP, LUSED, LSHA = self.segstate(LAST)
            if LSTOP == "nostate":
                moved = OUT / ("seg%02d.incomplete_%s" % (LAST, utc_stamp()))
                os.rename(str(OUT / ("seg%02d" % LAST)), str(moved))
                LAST, LSTOP, LUSED = LAST - 1, "none", 0
                if LAST >= 1:
                    LSTOP, LUSED, LSHA = self.segstate(LAST)
            if LAST >= 1 and LSHA != FITSHA:
                return Plan(rc=3, notes=["%s/seg%02d tunes the fit restart %s, FIT has %s" % (OUT, LAST, LSHA, FITSHA)])
        RES = []
        if LSTOP == "none":
            SEG, LEFT = 1, a.ft_maxnfev
        elif LSTOP == "wall_time":
            SEG, LEFT = LAST + 1, a.ft_maxnfev - LUSED
            RES = ["--resume-from", str(OUT / ("seg%02d" % LAST) / "finetune_state.npz")]
        elif LSTOP in ("tolerance", "max_nfev"):
            print("%s %s: segment %d finished (%s); nothing to do" % (utc_now(), OUT, LAST, LSTOP), flush=True)
            rc = self.final_link(LAST)
            return Plan(done=(rc == 0), rc=(None if rc == 0 else rc), segment=LAST)
        else:
            return Plan(rc=3, notes=["%s/seg%02d has stop '%s'; not continued" % (OUT, LAST, LSTOP)])
        if LEFT <= 0:
            return Plan(rc=3, notes=["%s: %d of %d evaluations used and incomplete: raise MAXNFEV" % (OUT, LUSED,
                                                                                                  a.ft_maxnfev)])
        O = OUT / ("seg%02d" % SEG)
        if O.exists():
            return Plan(rc=2, notes=["exists: %s" % O])
        OUT.mkdir(parents=True, exist_ok=True)
        VOPT = ["--validate"] if (SEG == 1 and a.ft_validate == 1) else []
        TOPT = ["--l-threads", str(a.l_threads)] if a.l_threads else []
        argv = ([str(HERE / "slope_finetune.py"), "--fit", str(FIT), "--out", str(O), "--facets", "late_only",
                 "--window-min", a.ft_window_min, "--devices", c.devices, "--lsmr-maxiter", str(a.ft_lsmr_maxiter),
                 "--wall-time", "slurm", "--segment", str(SEG), "--nfev-before", str(LUSED)] + VOPT + TOPT + RES
                + ["--max-nfev", str(LEFT)])
        env = {"FIT": str(FIT), "OUT": str(OUT), "FACETS": "late_only", "WINDOW_MIN": a.ft_window_min,
               "MAXNFEV": str(a.ft_maxnfev), "LSMR_MAXITER": str(a.ft_lsmr_maxiter), "VALIDATE": str(a.ft_validate),
               "MODEL_EXT": str(c.model_ext)}
        if a.l_threads:
            env["L_THREADS"] = str(a.l_threads)
        lines = ["fine tuning of %s -> %s segment %d (previous %d: %s, %d evaluations); evaluations %d of %d"
                 % (FIT, OUT, SEG, LAST, LSTOP, LUSED, LEFT, a.ft_maxnfev),
                 "FACETS late_only WINDOW_MIN %s; L_THREADS %s; MODEL_EXT %s -> %s (this process)"
                 % (a.ft_window_min, a.l_threads or "unset", c.model_ext, HERE)]
        if moved:
            lines.append("moved aside: %s" % moved)
        return Plan(argv=argv, env=env, out=O, segment=SEG, notes=lines, fit=FIT, commands_env=[SLOPE_ENV_UNSET],
                    summary="segment %d, %d evaluations" % (SEG, LEFT))

    def call(self, plan):
        sys.argv = list(plan.argv)
        return call_main(lambda: importlib.import_module("slope_finetune").main(list(plan.argv[1:])))

    def finish(self, plan, rc):
        if rc == 0:
            rc = self.final_link(plan.segment)        # slope_finetune.sbatch: final_link (it checks the stop) on exit 0
        if rc == EXIT_WALL_TIME:
            print("stopped before the wall-time limit: resubmit the same command to continue", flush=True)
        return rc


# ------------------------------------------------------------------------------------------------ nonneg kind (v6)
# The non-negative-by-construction models (user 2026-11-11: softplus A = m_b s ln(1+exp(L/s)), squared
# A = m_b (A_s (Y b)^2 + A_s sum ft)) of agent 1 (nonneg_mult.py: model and GPU patch = the process-kind module;
# nonneg_fit.py, nonneg_finetune.py; split chain jobs/nonneg_fit.sbatch, jobs/nonneg_finetune.sbatch, laid out as the
# slope scripts). Interface (agent 1's note of 2026-11-11, via the lead): parent: cfg =
# nonneg_fit.process_config_from_argv(fit argv), nonneg_mult.install_gpu(cfg), nonneg_mult.set_L_threads(n) (idempotent);
# worker: nonneg_mult.install_worker_handlers(); the held term is the stages' own (host only, set_stage_held at each
# stage's start); fit nonneg_fit.main() (sys.argv); fine tuning nonneg_finetune.main(argv) (binds through
# gpu_model.REGISTRY, close() sends nonneg_detach).
NONNEG_MODELS = {"softplus": {"c1_start": "d2"}, "squared": {"c1_start": "lsq"}}
D2_FINAL = RUN / "clip_refit" / "D2" / "final"


class NonnegFitStage(SlopeFitStage):
    """jobs/nonneg_fit.sbatch (agent 1, snapshot nonneg_20261111T044920Z) with CHAIN=1 and FITDIR=<chain>/c<k>/fit:
    the slope_mult_fit.sbatch logic; MOPT = --model M --ft-templates-from F [--held H] [--start-file S (START lsq,
    every segment)], placed after --slope-bundle, before the start options."""

    kind = "nonneg_fit"
    split_script = "nonneg_fit.sbatch"
    code_files = ("nonneg_fit.py", "nonneg_mult.py", "nonneg_problem.py", "nonneg_start.py", "nonneg_finetune.py",
                  "test_nonneg.py", "slope_mult.py", "slope_mult_problem.py", "test_slope_mult.py",
                  "positive_finetune.py", "clip_problem.py", "clip_refit.py", "grain_fit_fixed.py", "fixed_parameters.py",
                  "grain_fit.py", "grain_crossing.py", "positive_problem.py", "late_response_eval.py",
                  "test_clip_refit.py", "test_positive_fit.py")              # nonneg_fit.sbatch FILES
    run_copy_files = code_files
    summary_regex = re.compile(r"^(nonneg|clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |"
                               r"grain_fit: |fit data|selection profile)")
    env_regex = re.compile(r"^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|NFEV|CHECK|CHECKPARAMS|TESTNAME|TESTMODE|MODEL|START|"
                           r"START_RESTART|START_FILE|FT_FROM|L_THREADS|HELD|GRAIN_FROM|FITDIR|BUNDLE|CHAIN|MODEL_EXT|"
                           r"SLOPE_MULT_MIN|"
                           r"OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|"
                           r"SLURM_JOB_GPUS|SLURM_JOB_NODELIST|SLURM_JOB_NUM_NODES|CUDA_VISIBLE_DEVICES)=")

    def __init__(self, chain, k):
        a = chain.args
        self.chain, self.k = chain, int(k)
        self.name = "c%dfit" % self.k
        self.directory = chain.root / ("c%d" % self.k) / "fit"
        self.nonneg = NONNEG_MODELS[a.slope_model]
        self.model = {"fit_args": self._model_args, "c1_start": self.nonneg["c1_start"]}
        if self.k == 1:
            self.start = a.slope_start or self.nonneg["c1_start"]
            self.start_restart = Path(a.source)
            self.held = self.grain_from = None
        else:
            prev = chain.root / ("c%d" % (self.k - 1))
            self.start = "chain"
            self.start_restart = prev / "fit" / "final" / "fit_restart.npz"
            self.held = self.grain_from = prev / "ft" / "finetune_state.npz"
        self.start_set = ""

    def _model_args(self, c):
        """nonneg_fit.sbatch MOPT without --held (SlopeFitStage._argv appends --held after it, as MOPT does): --model,
        --ft-templates-from; --start-file and --slope-mult-min come after --held in MOPT (see _argv)."""
        return ["--model", c.args.slope_model, "--ft-templates-from", str(readlink_f(c.args.ft_templates_from))]

    def label(self):
        return "nonneg_fit_%s_%s" % (self.directory.parent.name, self.directory.name)

    def _argv(self, SEG, LEFT, SOPT, O):
        """nonneg_fit.sbatch lines 79-81: MOPT = --model M --ft-templates-from F [--held H] [--start-file S (START lsq)]
        [--slope-mult-min V (SLOPE_MULT_MIN; user 2026-11-12, softplus: 1; every cycle and segment: the restart marker
        refuses a resume or chain start with another value)], then SOPT."""
        argv = SlopeFitStage._argv(self, SEG, LEFT, SOPT, O)
        argv[0] = str(HERE / "nonneg_fit.py")
        tail = []
        if self.start == "lsq":
            tail += ["--start-file", str(readlink_f(self.chain.args.start_file))]
        if self.chain.args.slope_mult_min is not None:
            tail += ["--slope-mult-min", self.chain.args.slope_mult_min]
        k = argv.index(SOPT[0]) if SOPT else len(argv)
        argv[k:k] = tail
        return argv

    def check_held_identity(self, k):
        """A completed fit: its restart's model_extension.held (agent 1: path, key, sha256 = sha256 of the held_next
        ARRAY bytes, float64 C-contiguous, nonneg_mult.sha256_bytes) must be the previous fine tuning's state file
        (cycle k > 1) and null in cycle 1."""
        meta = _npz_meta(self.directory / ("seg%02d" % k) / "fit_restart.npz")
        rec = (meta.get("model_extension") or {}).get("held")
        if self.held is None:
            return None if not rec else "%s/seg%02d has a held term %s; cycle 1 has none" % (self.directory, k, rec)
        if not isinstance(rec, dict):
            return "%s/seg%02d has no held-term record (model_extension.held)" % (self.directory, k)
        f = readlink_f(self.held)
        with np.load(str(f), allow_pickle=False) as z:
            h = np.ascontiguousarray(np.asarray(z[rec.get("key") or "held_next"], dtype=np.float64))
        want = hashlib.sha256(np.ascontiguousarray(h).tobytes()).hexdigest()
        if rec.get("sha256") != want or os.path.realpath(str(rec.get("path", ""))) != str(f):
            return "held term of %s/seg%02d (%s) is not %s (array sha256 %s)" % (self.directory, k, rec, f, want)
        return None

    def prepare(self):
        a = self.chain.args
        missing = [f for f in ([a.ft_templates_from] + ([a.start_file] if self.start == "lsq" else []))
                   if not f or not Path(f).exists()]
        if missing:
            return Plan(rc=2, notes=["missing: %s" % missing])
        plan = SlopeFitStage.prepare(self)
        if plan.argv is not None:
            plan.env.pop("BUMPS", None)
            plan.env.pop("SLOPE_MODEL", None)
            plan.env.update(MODEL=a.slope_model, FT_FROM=str(readlink_f(a.ft_templates_from)))
            if a.slope_mult_min is not None:
                plan.env["SLOPE_MULT_MIN"] = a.slope_mult_min
            if self.start == "lsq":
                plan.env["START_FILE"] = str(readlink_f(a.start_file))
            plan.env.pop("L_MODE", None)
        return plan

    def call(self, plan):
        sys.argv = list(plan.argv)
        return call_main(lambda: importlib.import_module("nonneg_fit").main())


class NonnegFineTuneStage(SlopeFineTuneStage):
    """jobs/nonneg_finetune.sbatch (agent 1; the slope_finetune.sbatch layout)."""

    kind = "nonneg_finetune"
    split_script = "nonneg_finetune.sbatch"
    code_files = ("nonneg_finetune.py", "nonneg_mult.py", "nonneg_problem.py", "slope_mult.py", "positive_finetune.py",
                  "positive_problem.py", "grain_crossing.py", "clip_problem.py", "test_nonneg.py")   # MODFILES
    run_copy_files = code_files
    fit_files = ("fit_restart.npz", "fit_manifest.json", "nonneg_facets.npz", "grain_term.json", "COPS_model.dat",
                 "late_window_model.dat", "fit_result.txt")              # agent 1, 2026-11-11 (check all)
    summary_regex = re.compile(r"^(nonneg|Traceback|.*Error)")

    def prepare(self):
        plan = SlopeFineTuneStage.prepare(self)
        if plan.argv is not None:
            plan.argv[0] = str(HERE / "nonneg_finetune.py")
        return plan

    def call(self, plan):
        sys.argv = list(plan.argv)
        return call_main(lambda: importlib.import_module("nonneg_finetune").main(list(plan.argv[1:])))


def _parent_setup_nonneg(args, chain=None):
    """Agent 1's process level: cfg = nonneg_fit.process_config_from_argv(the cycle-1 fit's argv) (constants only, never
    a held term), nonneg_mult.install_gpu(cfg), nonneg_mult.set_L_threads(n) -- before the baseline."""
    preset = sorted(k for k in os.environ if k.startswith(("COPS_POS_", "COPS_SLOPEMUL_")))
    if preset:
        raise SystemExit("chain_onejob: %s set: the nonneg kind refuses the other patches' environment" % preset)
    import nonneg_fit
    import nonneg_mult as nm
    first = next(st for st in chain.stages if isinstance(st, NonnegFitStage))
    if first.k != 1:
        raise SystemExit("chain_onejob: the nonneg layout must start with the cycle-1 fit")
    cfg = nonneg_fit.process_config_from_argv(first.config_argv())
    nm.install_gpu(cfg)
    if args.l_threads:
        nm.set_L_threads(int(args.l_threads), args.l_mode)
    return {"nonneg_config": json.dumps(cfg, default=str, sort_keys=True)[:400], "l_threads": args.l_threads}


def _worker_setup_nonneg():
    import nonneg_mult as nm
    nm.install_worker_handlers()


def _layout_nonneg(chain):
    stages = []
    for k in range(1, chain.args.cycles + 1):
        stages += [NonnegFitStage(chain, k), NonnegFineTuneStage(chain, k)]
    return stages


STAGE_TYPES = {"positive_fit": PositiveFitStage, "positive_finetune": PositiveFineTuneStage,
               "slope_fit": SlopeFitStage, "slope_finetune": SlopeFineTuneStage,
               "nonneg_fit": NonnegFitStage, "nonneg_finetune": NonnegFineTuneStage}


def _parent_setup_positive(args, chain=None):
    """Process level, before the baseline: COPS_POS_* (positive_problem.set_positive_environment with the chain's
    ratio, as positive_fit.py does at import) and positive_activity.install()."""
    import positive_problem
    env = positive_problem.set_positive_environment(ratio=args.ratio)
    import positive_activity
    if not positive_activity.install():
        raise SystemExit("chain_onejob: positive_activity.install() did nothing (COPS_POS_ENABLE?)")
    for name in ("positive_problem",):               # a stage module: imported afresh by every stage
        sys.modules.pop(name, None)
    return {"positive_environment": env, "positive_activity": positive_activity.config()}


def _layout_positive(chain):
    stages = []
    for k in range(1, chain.args.cycles + 1):
        stages += [PositiveFitStage(chain, k), PositiveFineTuneStage(chain, k)]
    return stages


def _slope_process_config(argv):
    """The GPU configuration a slope fit with these arguments installs (slope_mult_fit.py's own steps:
    slope_mult_problem.rewrite_arguments, model_config). Agent 1 may provide slope_mult_fit.process_config(argv)."""
    import slope_mult_fit
    if hasattr(slope_mult_fit, "process_config"):
        return slope_mult_fit.process_config(list(argv))
    import slope_mult_problem as smp
    head_end = argv.index("--") if "--" in argv else len(argv)
    args_in = [a for k, a in enumerate(argv)]
    for opt in ("--l-threads", "--l-mode"):                     # removed by slope_mult_fit.main before rewriting
        if opt in args_in[:head_end]:
            i = args_in.index(opt)
            del args_in[i:i + 2]
            head_end -= 2
    cfg, _, _ = smp.rewrite_arguments(args_in)
    model_cfg, _ = smp.model_config(cfg)
    return model_cfg


def _parent_setup_slope(args, chain=None):
    """Process level, before the baseline: the slope-multiplier GPU patch with the cycle-1 fit's configuration
    (slope_mult.install_gpu; every stage's own install_gpu call must then find the same configuration: the per-cycle
    held term is set per stage, requirement A). COPS_POS_* must be unset (both patches modify the same methods)."""
    preset = sorted(k for k in os.environ if k.startswith("COPS_POS_"))
    if preset:
        raise SystemExit("chain_onejob: %s set: the slope kind refuses the positive patch's environment" % preset)
    import slope_mult as sm
    first = next(st for st in chain.stages if isinstance(st, SlopeFitStage))
    if first.k != 1:
        raise SystemExit("chain_onejob: the slope layout must start with the cycle-1 fit")
    cfg = _slope_process_config(first.config_argv())
    if cfg.get("held"):
        raise SystemExit("chain_onejob: the cycle-1 configuration has a held term")
    sm.install_gpu(cfg)
    if hasattr(sm, "install_worker_handlers"):
        sm.install_worker_handlers()
    return {"slope_config": sm.canonical(cfg)[:400], "installed": sm._STATE.get("installed")}


def _worker_setup_slope():
    """Spawned workers: slope_mult (its _pipe_worker_main installs the patch from COPS_SLOPEMUL_FROZEN_CONFIG, as in
    the split stages' workers) and, when agent 1 provides them, its worker handlers."""
    import slope_mult as sm
    if hasattr(sm, "install_worker_handlers"):
        sm.install_worker_handlers()


def _layout_slope(chain):
    stages = []
    for k in range(1, chain.args.cycles + 1):
        stages += [SlopeFitStage(chain, k), SlopeFineTuneStage(chain, k)]
    return stages


PROCESS_KINDS = {"positive": {"parent": _parent_setup_positive, "worker": _worker_setup_positive,
                              "layout": _layout_positive, "keep": ("positive_activity",)},
                 "slope": {"parent": _parent_setup_slope, "worker": _worker_setup_slope, "layout": _layout_slope,
                           "keep": ("slope_mult",)},
                 "nonneg": {"parent": _parent_setup_nonneg, "worker": _worker_setup_nonneg, "layout": _layout_nonneg,
                            "keep": ("nonneg_mult",)}}


# ================================================================================================ events
class Events:
    def __init__(self, path):
        self.path = Path(path)

    def append(self, **record):
        record = dict({"utc": utc_now(), "job": os.environ.get("SLURM_JOB_ID"), "host": socket.gethostname(),
                       "pid": os.getpid()}, **record)
        with self.path.open("a") as f:
            f.write(json.dumps(record, sort_keys=True) + "\n")

    def last(self, seg_dir):
        if not self.path.is_file():
            return None
        out = None
        for line in self.path.read_text().splitlines():
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get("segment_dir") == str(seg_dir):
                out = r
        return out


# ================================================================================================ the chain
class Chain:
    def __init__(self, args):
        self.args = args
        self.root = Path(os.path.abspath(args.chain_root))
        self.base = Path(args.run) / "base"
        self.point_table = Path(args.point_table)
        self.bundle = Path(args.bundle)
        self.model_ext = Path(args.model_ext)
        self.c1_start_set = args.c1_start_set or ""
        self.devices = args.devices or (visible_devices() if not args.plan_only else "0")
        self.events = Events(self.root / EVENTS_NAME)
        self.stages = PROCESS_KINDS[args.kind]["layout"](self)
        self.job = os.environ.get("SLURM_JOB_ID") or "nojob"
        self.jobprov = Path(args.job_provenance) if args.job_provenance else None
        self.bound = None               # {"stage", "identity" (dataset_identity), "seconds", "key", "summary"}
        self.reset_verification()       # per stage: the model objects already checked, Y digests
        self.current = None             # the running stage
        self.summary = []
        self.tee = None
        self.baseline = None
        self.current_detach = None

    # ---------------------------------------------------------------- registry and the shared backend
    def reset_verification(self):
        """Per stage (v5): the model OBJECTS already checked (a WeakSet: an id() of a freed model can be reused by a
        new object, which v4's id-keyed set then skipped -- found by CPU test 5, job 124089) and the Y digests by
        array identity (the arrays are held only until the stage ends)."""
        self._verified = weakref.WeakSet()
        self._y_cache = []

    def _is_verified(self, model):
        try:
            return model in self._verified
        except TypeError:                      # not weakly referenceable: check every time
            return False

    def _mark_verified(self, model):
        try:
            self._verified.add(model)
        except TypeError:
            pass

    def _y_digest(self, Y):
        for arr, digest in self._y_cache:
            if arr is Y:
                return digest
        return None

    def _registry(self):
        from cops_cuda import gpu_model
        return gpu_model.REGISTRY

    def install_registry(self):
        from cops_cuda import gpu_model
        R = gpu_model.REGISTRY
        cls = type(R)
        chain = self

        def configure(device_id=0, mode="tree", records_per_block=4096, devices=None, transport="local",
                      activity=None):
            if not R.backends:
                return cls.configure(R, device_id, mode, records_per_block, devices=devices, transport=transport,
                                     activity=activity)
            want = {"devices": [int(d) for d in (devices if devices else [device_id])], "mode": str(mode),
                    "records_per_block": int(records_per_block), "transport": str(transport), "activity": activity}
            have = {"devices": list(R.devices), "mode": R.mode, "records_per_block": R.records_per_block,
                    "transport": R.transport, "activity": R.activity}
            if want != have:
                raise RuntimeError("chain_onejob: stage %s configures %s; the bound backend has %s"
                                   % (chain.current and chain.current.name, want, have))
            return None

        def backend_for(model, basis, scale):
            key = (str(model.run), tuple(int(r) for r in model.rows), tuple(R.devices), R.mode)
            if R.backends and key not in R.backends:
                key = chain._rekey(model, key)
            fresh = not R.backends
            if not fresh:
                chain.verify_model(model)          # every new model object, also in the binding stage (v4)
            t0 = time.perf_counter()
            backend = cls.backend_for(R, model, basis, scale)
            if fresh:
                chain.record_bind(model, backend, time.perf_counter() - t0, key)
            return backend

        def close():
            say("REGISTRY.close() of stage %s: kept (the chain closes the backend at its end)"
                % (chain.current.name if chain.current else "-"))

        R.configure, R.backend_for, R.close = configure, backend_for, close

    def _rekey(self, model, key):
        R = self._registry()
        if len(R.backends) != 1:
            raise RuntimeError("chain_onejob: %d bound backends" % len(R.backends))
        (old,) = list(R.backends)
        same = (os.path.realpath(old[0]) == os.path.realpath(key[0]) and old[1:] == key[1:])
        if not same:
            raise RuntimeError("chain_onejob: stage %s needs rows/devices/mode %s; the bound backend has %s (one bound "
                               "dataset per chain)" % (self.current and self.current.name, key[1:][1:], old[1:][1:]))
        R.backends[key] = R.backends.pop(old)
        return key

    def record_bind(self, model, backend, seconds, key):
        stage = self.current.name if self.current else None
        ident = dataset_identity(model)
        self._y_cache.append((model.Y, ident["Y"]))
        self.bound = {"stage": stage, "identity": ident, "seconds": float(seconds), "key": key,
                      "summary": backend.summary(), "utc": utc_now()}
        self._mark_verified(model)
        say("bound in stage %s in %.1f s: %s (kept for the whole chain)" % (stage, seconds, backend.summary()))

    def verify_model(self, model):
        """Every model that reaches the bound backend (each new model object of each stage, also the second model of
        the binding stage) must have the bound dataset's identity (dataset_identity: same base, rows, records,
        evaluator inputs; the selection columns and free-text header lines may differ, as between MF and M)."""
        stage = self.current.name if self.current else None
        if self._is_verified(model):
            return
        t0 = time.perf_counter()
        y = self._y_digest(model.Y)
        ident = dataset_identity(model, y)
        if y is None:
            self._y_cache.append((model.Y, ident["Y"]))
        bnd = self.bound
        problems = identity_differences(ident, bnd["identity"])
        if problems:
            raise RuntimeError("chain_onejob: stage %s: a model with another dataset identity than the one bound in "
                               "stage %s reaches the backend (fields %s): not the same base / records" %
                               (stage, bnd["stage"], problems))
        self._mark_verified(model)
        say("stage %s: model %s verified against the bound dataset (%d identity fields equal; %.1f s)"
            % (stage, type(model).__name__, len(set(ident) & set(bnd["identity"])), time.perf_counter() - t0))

    def shared_backend(self):
        R = self._registry()
        if not R.backends:
            return None, None
        ((key, (backend, basis)),) = list(R.backends.items())
        return key, backend

    def reset_basis(self):
        """The bound backend in the state of a fresh bind with basis None (both split stages start so)."""
        key, backend = self.shared_backend()
        if backend is None:
            return None
        if backend.closed or backend.broken:
            raise RuntimeError("chain_onejob: the bound backend is closed or broken")
        t0 = time.perf_counter()
        backend.rebind_basis(None, 1.0)
        seconds = time.perf_counter() - t0
        self._registry().backends[key] = (backend, None)
        say("stage %s: template basis reset on the bound records (rebind_basis(None)) in %.2f s"
            % (self.current.name, seconds))
        return seconds

    def install_finetune_engine(self, module, stage):
        """positive_finetune.GpuEngine on the bound backend (module docstring)."""
        chain = self
        Base = module.GpuEngine

        class ChainGpuEngine(Base):
            def __init__(self, ctx, devices, records_per_block=4096):
                from cops_cuda import multi_gpu
                real = multi_gpu.MultiGpuBackend

                def factory(spec, model, enhancement_basis=None, enhancement_scale=1.0, **kw):
                    return chain.finetune_backend(real, spec, model, enhancement_basis, enhancement_scale, kw)
                multi_gpu.MultiGpuBackend = factory
                try:
                    Base.__init__(self, ctx, devices, records_per_block)
                finally:
                    multi_gpu.MultiGpuBackend = real

            def close(self):
                backend = getattr(self, "backend", None)
                if backend is None:
                    return Base.close(self)

                class _Detach:
                    def close(_self):
                        chain.detach_finetune(backend)
                self.backend = _Detach()
                try:
                    Base.close(self)
                finally:
                    self.backend = backend

        ChainGpuEngine.__name__ = ChainGpuEngine.__qualname__ = "ChainGpuEngine"
        module.GpuEngine = ChainGpuEngine

    def finetune_backend(self, real, spec, model, basis, scale, kw):
        R = self._registry()
        want = {"devices": [int(d) for d in kw.get("devices", ())], "mode": kw.get("mode", "tree"),
                "records_per_block": int(kw.get("records_per_block", 4096)), "transport": kw.get("transport", "local"),
                "activity": kw.get("activity")}
        if basis is not None or float(scale) != 1.0:
            raise RuntimeError("chain_onejob: the fine tuning binds with a basis (%s, %s)" % (type(basis), scale))
        if not R.backends:                     # a chain resumed at a fine tuning: bind here, register for the fits
            R.configure(want["devices"][0], want["mode"], want["records_per_block"], devices=want["devices"],
                        transport=want["transport"], activity=want["activity"])
            t0 = time.perf_counter()
            backend = real(spec, model, basis, scale, **kw)
            seconds = time.perf_counter() - t0
            key = (str(model.run), tuple(int(r) for r in model.rows), tuple(R.devices), R.mode)
            R.backends[key] = (backend, None)
            R.seconds["bind"] += seconds
            R.events.append({"event": "bind", "rows": len(key[1]), "devices": list(R.devices),
                             "records": int(sum(backend.block_records)), "enhancement_columns": 0,
                             "enhancement_basis_sha256": None, "seconds": seconds,
                             "by": "chain_onejob (fine-tuning stage %s)" % self.current.name})
            self.record_bind(model, backend, seconds, key)
            return backend
        have = {"devices": list(R.devices), "mode": R.mode, "records_per_block": R.records_per_block,
                "transport": R.transport, "activity": R.activity}
        if want != have:
            raise RuntimeError("chain_onejob: the fine tuning asks for %s; the bound backend has %s" % (want, have))
        key, backend = self.shared_backend()
        if (os.path.realpath(key[0]) != os.path.realpath(str(model.run))
                or list(key[1]) != [int(r) for r in model.rows]):
            raise RuntimeError("chain_onejob: the fine tuning's rows differ from the bound backend's")
        if spec is not model._spec:
            raise RuntimeError("chain_onejob: unexpected spec object")
        self.verify_model(model)
        if backend.enhancement_basis is not None:
            self.reset_basis()
        say("stage %s: the fine tuning uses the backend bound in stage %s (%.1f s bind there): %s"
            % (self.current.name, self.bound["stage"], self.bound["seconds"], backend.summary()))
        return backend

    def detach_finetune(self, backend):
        if backend.closed or backend.broken:
            say("WARNING detach: the backend is closed or broken")
            return
        for w in backend.workers:
            backend._send(w, (DETACH_REQUEST,), DETACH_REQUEST)
        infos = []
        for w in backend.workers:
            reply = backend._receive(w, DETACH_REQUEST)
            if not (isinstance(reply, tuple) and reply and reply[0] == "ok"):
                raise RuntimeError("chain_onejob: unexpected reply to %s: %r" % (DETACH_REQUEST, reply))
            infos.append(reply[1])
        errors = [e for i in infos for e in i.get("errors", [])]
        say("stage %s: fine-tuning worker state released on %d workers (pinned %s, reservations %s, shm %s)%s"
            % (self.current.name, len(infos), [i.get("host_unregistered") for i in infos],
               [i.get("reservations_released") for i in infos], [i.get("shm_closed") for i in infos],
               "; ERRORS %s" % errors if errors else ""))
        self.current_detach = infos

    # ---------------------------------------------------------------- provenance
    def stage_provenance(self, stage, plan):
        P = stage.directory / ("provenance-%s-%s" % (self.job, stage.name))
        P.mkdir(parents=True, exist_ok=False)
        jp = self.jobprov
        for f in ("chain_onejob.sbatch", "module_list.txt", "nvidia_smi_q.txt", "toolchain.txt", "pip_list.txt"):
            if jp is not None and (jp / f).is_file():
                shutil.copy2(str(jp / f), str(P / f))
        for f in tuple(stage.code_files) + tuple(getattr(stage, "optional_code_files", ())) + ("chain_onejob.py",):
            if (HERE / f).is_file():
                shutil.copy2(str(HERE / f), str(P / f))
        if (HERE / "kernels").is_dir():
            shutil.copytree(str(HERE / "kernels"), str(P / "kernels"))
        with (P / "run_copy.sha256").open("w") as f:          # the files that ran (HERE) = these copies
            for name in tuple(stage.run_copy_files) + ("chain_onejob.py",):
                if (HERE / name).is_file():
                    f.write("%s  %s\n" % (sha256_file(HERE / name), name))
            for p in sorted((HERE / "kernels").iterdir()) if (HERE / "kernels").is_dir() else []:
                if p.is_file():
                    f.write("%s  kernels/%s\n" % (sha256_file(p), p.name))
        if jp is not None and (jp / "model_ext.sha256").is_file():
            shutil.copy2(str(jp / "model_ext.sha256"), str(P / "model_ext.sha256"))
        run_to_file("git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS status "
                    "--short", P / "git.txt", shell=True)
        env = "\n".join(sorted("%s=%s" % kv for kv in os.environ.items() if stage.env_regex.match("%s=%s" % kv)))
        write_text(P / "environment.txt", env + "\n")
        commit = (P / "git.txt").read_text().splitlines()[0] if (P / "git.txt").is_file() else "?"
        exports = []
        if self.args.job_script and Path(self.args.job_script).is_file():
            exports = [l for l in Path(self.args.job_script).read_text().splitlines()
                       if l.startswith("export ") and "COPS_POS_" not in l]
        lines = ["#!/bin/bash",
                 "# job %s on %s, %s; AMPS %s; one-job chain stage %s (chain_onejob.py: one process, the backend "
                 "bound once; run as is, this reruns the stage alone in a new process)"
                 % (self.job, socket.gethostname(), utc_now(), commit, stage.name)] + exports + \
            ["export MODEL_EXT=%s" % self.model_ext] + list(plan.commands_env) + ["cd %s" % APP,
                                                                                    shell_line(PY + list(plan.argv))]
        write_text(P / "commands.sh", "\n".join(lines) + "\n")
        notes = list(plan.notes) + [
            "one-job chain: job %s, stage %s (%d of %d), driver %s (pid %d); stage log %s"
            % (self.job, stage.name, self.stages.index(stage) + 1, len(self.stages), Path(__file__).resolve(),
               os.getpid(), P / stage.log_name),
            "split-chain logic replicated: jobs/%s (CHAIN=1 for fits); module_list, nvidia_smi_q, toolchain, pip_list: "
            "the job's (taken at job start)" % stage.split_script]
        write_text(P / "segment.txt", "\n".join(notes) + "\n")
        for line in notes:
            print(line, flush=True)
        run_to_file(["free", "-g"], P / "free_at_start.txt")
        run_to_file(["df", "-h", str(self.root)], P / "df_at_start.txt")
        stage.extra_provenance(plan, P)
        return P

    # ---------------------------------------------------------------- samplers (as positive_finetune.sbatch)
    def start_samplers(self, P):
        procs = []
        try:
            f = open(P / "nvidia_smi_samples.csv", "w")
            procs.append((subprocess.Popen(["nvidia-smi", "--query-gpu=timestamp,index,memory.used,memory.total,"
                                            "utilization.gpu", "--format=csv", "-l", "20"], stdout=f,
                                           stderr=subprocess.STDOUT), f))
        except Exception as exc:
            say("nvidia-smi sampler not started: %s" % exc)
        stop = threading.Event()

        def free_loop():
            with open(P / "free_samples.txt", "w") as g:
                while not stop.is_set():
                    try:
                        out = subprocess.run(["free", "-g"], stdout=subprocess.PIPE, universal_newlines=True).stdout
                        mem = [l.split() for l in out.splitlines() if l.startswith("Mem:")][0]
                        g.write("%s used_GB %s free_GB %s available_GB %s\n" % (utc_now(), mem[2], mem[3], mem[6]))
                        g.flush()
                    except Exception:
                        pass
                    stop.wait(30.0)
        t = threading.Thread(target=free_loop, daemon=True)
        t.start()
        return procs, stop, t

    @staticmethod
    def stop_samplers(samplers):
        procs, stop, t = samplers
        stop.set()
        for p, f in procs:
            try:
                p.terminate()
                p.wait(10)
            except Exception:
                pass
            f.close()
        t.join(5.0)

    # ---------------------------------------------------------------- one stage
    def remaining_seconds(self):
        from cops_cuda import wall_time
        try:
            return wall_time.slurm_end_time() - time.time()
        except ValueError:
            return None

    def execute(self, stage, plan):
        a = self.args
        self.current = stage
        self.current_detach = None
        self.reset_verification()
        rec = {"stage": stage.name, "kind": stage.kind, "directory": str(stage.directory), "segment": plan.segment,
               "out": str(plan.out), "started_utc": utc_now(), "job": self.job}
        os.environ.update(plan.env)                  # the split job script's variables (environment.txt)
        os.chdir(str(APP))
        P = self.stage_provenance(stage, plan)
        rec["provenance"] = str(P)
        bound_before = self.bound is not None
        reset_s = self.reset_basis() if bound_before else None
        R = self._registry()
        if bound_before:
            R.calls = {"pressures": 0, "jacobian": 0, "seen_faces": 0}
            R.seconds = {"pressures": 0.0, "jacobian": 0.0, "bind": 0.0}
            R.events = [{"event": "chain_onejob_reuse", "job": self.job, "stage": stage.name,
                         "bound_in_stage": self.bound["stage"], "bind_seconds_there": self.bound["seconds"],
                         "reset_basis_seconds": reset_s,
                         "note": "one-job chain (model_ext/chain_onejob.py): records bound once per job; the template "
                                 "basis was reset to None before this stage"}]
        samplers = self.start_samplers(P) if stage.kind.endswith("finetune") else None
        self.events.append(event="start", stage=stage.name, segment_dir=str(plan.out), argv=plan.argv)
        say("stage %s: %s; log %s" % (stage.name, plan.summary, P / stage.log_name))
        self.tee.switch(P / stage.log_name)
        t0 = time.time()
        rc = 1
        try:
            rc = stage.call(plan)
        finally:
            wall = time.time() - t0
            self.tee.switch(None)
            self.reset_verification()               # nothing of the stage is held past its end
            restore = self.baseline.restore()
            if samplers is not None:
                self.stop_samplers(samplers)
        say("stage %s: exit %d after %.1f s; restored the process state (%d stage modules purged, %d attributes "
            "restored, %d deleted%s)" % (stage.name, rc, wall, len(restore["purged"]), len(restore["restored"]),
                                         len(restore["deleted"]), ", new kept modules %s" % restore["new_kept_modules"]
                                         if restore["new_kept_modules"] else ""))
        log_text = (P / stage.log_name).read_text(errors="replace") if (P / stage.log_name).is_file() else ""
        for line in log_text.splitlines():
            if stage.summary_regex.match(line):
                print(line[:400], flush=True)
        for line in log_text.splitlines()[-6:]:
            print(line[:250], flush=True)
        run_to_file(["df", "-h", str(self.root)], P / "df_at_end.txt")
        run_to_file(["du", "-sh", str(plan.out)], P / "df_at_end.txt", append=True)
        if stage.kind.endswith("finetune"):
            run_to_file(["free", "-g"], P / "free_at_end.txt")
        restart = stage.restart_written(plan)
        self.events.append(event="end", stage=stage.name, segment_dir=str(plan.out), rc=rc, restart_written=restart,
                           wall_seconds=wall)
        rc_final = stage.finish(plan, rc)
        backend = {"bound_before_stage": bound_before, "bound": {k: v for k, v in (self.bound or {}).items()
                                                                  if k in ("stage", "seconds", "summary", "utc")},
                   "reset_basis_seconds": reset_s, "registry_events": R.events if stage.kind.endswith("_fit") else None,
                   "finetune_detach": self.current_detach}
        write_text(P / "backend.json", json.dumps(backend, indent=1, default=str) + "\n")
        rec.update(rc_stage=rc, rc=rc_final, wall_seconds=wall, finished_utc=utc_now(), restart_written=restart,
                   bound_here=(not bound_before and self.bound is not None),
                   bind_seconds=(self.bound["seconds"] if (not bound_before and self.bound) else 0.0),
                   reset_basis_seconds=reset_s, purged_modules=len(restore["purged"]),
                   rss_gb=self.rss_gb())
        self.summary.append(rec)
        self.write_summary()
        self.current = None
        return rc_final

    @staticmethod
    def rss_gb():
        try:
            for line in Path("/proc/self/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1048576.0
        except OSError:
            pass
        return None

    def write_summary(self):
        if self.jobprov is None:
            return
        data = {"job": self.job, "host": socket.gethostname(), "chain_root": str(self.root),
                "args": {k: v for k, v in vars(self.args).items()}, "stages": self.summary,
                "bound": {k: v for k, v in (self.bound or {}).items() if k in ("stage", "seconds", "summary", "utc")},
                "updated_utc": utc_now()}
        tmp = self.jobprov / "chain_summary.json.tmp"
        tmp.write_text(json.dumps(data, indent=1, default=str) + "\n")
        os.replace(str(tmp), str(self.jobprov / "chain_summary.json"))

    # ---------------------------------------------------------------- the whole chain
    def plan_only(self):
        """Print every stage's state and the command it would run now (no GPU, nothing written but links/mkdir are
        NOT made: prepare() is not called for stages that depend on unfinished ones)."""
        for st in self.stages:
            if isinstance(st, PositiveFitStage):
                print("%s: %s state %s" % (st.name, st.directory, st.state()[:3]))
            else:
                last = 0
                if st.directory.is_dir():
                    segs = [int(p.name[3:]) for p in st.directory.glob("seg[0-9][0-9]")]
                    last = max(segs) if segs else 0
                print("%s: %s state %s" % (st.name, st.directory, st.segstate(last)[:2] if last else ("none", 0)))
        return 0

    def parent_setup(self):
        """The process kind's setup; afterwards every stage module it imported is removed from sys.modules (only the
        kind's patch module stays) and sys.path / the reference binding are as before, so that the baseline is the
        state of a fresh process plus the patch."""
        from cops_cuda import reference
        path0, bound0 = list(sys.path), reference._bound_root
        keep = set(PROCESS_KINDS[self.args.kind]["keep"]) | {"__main__", "__mp_main__"}
        before = set(sys.modules)
        setup = PROCESS_KINDS[self.args.kind]["parent"](self.args, self)
        roots = [os.path.realpath(str(r)) for r in (HERE, APP / "tools", REFERENCE_TOOLS)]

        def under(f):
            return bool(f) and any(os.path.realpath(f).startswith(r + os.sep) for r in roots)
        # modules the patch module refers to (import X / from X import f): they stay with it (kept and restored)
        extra = set()
        for kname in PROCESS_KINDS[self.args.kind]["keep"]:
            km = sys.modules.get(kname)
            for v in (vars(km).values() if km is not None else ()):
                mod = v if isinstance(v, type(sys)) else sys.modules.get(getattr(v, "__module__", None) or "")
                if mod is not None and mod is not km and mod.__name__ in sys.modules and \
                        under(getattr(mod, "__file__", None)):
                    extra.add(mod.__name__)
        keep |= extra
        self.keep_names = tuple(PROCESS_KINDS[self.args.kind]["keep"]) + tuple(sorted(extra))
        popped = []
        for name in list(sys.modules):
            m = sys.modules.get(name)
            f = getattr(m, "__file__", None) if m is not None else None
            if (name not in before and name not in keep and f
                    and any(os.path.realpath(f).startswith(r + os.sep) for r in roots)):
                del sys.modules[name]
                popped.append(name)
        sys.path[:] = path0
        reference._bound_root = bound0
        setup = dict(setup, stage_modules_dropped_after_setup=sorted(popped), kept_with_the_patch=sorted(extra))
        return setup

    def run(self):
        a = self.args
        self.tee = Tee()
        rc = 0
        try:
            say("job %s on %s, pid %d; chain root %s; %d stages: %s; code %s; MODEL_EXT %s"
                % (self.job, socket.gethostname(), os.getpid(), self.root, len(self.stages),
                   " ".join(s.name for s in self.stages), HERE, self.model_ext))
            os.environ[KIND_ENV] = a.kind
            setup = self.parent_setup()
            say("process setup (%s): %s" % (a.kind, json.dumps(setup, default=str)[:600]))
            self.baseline = Baseline([HERE, APP / "tools", REFERENCE_TOOLS], keep_names=self.keep_names)
            info = self.baseline.capture()
            say("baseline: %d modules, %d kept modules, %d classes" % (info["modules"], len(info["kept_modules"]),
                                                                      info["classes"]))
            self.install_registry()
            for st in self.stages:
                self.current = st
                plan = st.prepare()
                self.current = None
                if plan.done:
                    say("stage %s: complete (verified); skipped" % st.name)
                    self.summary.append({"stage": st.name, "status": "complete_before", "segment": plan.segment})
                    continue
                if plan.rc is not None:
                    for line in plan.notes:
                        print(line, file=sys.stderr, flush=True)
                    say("stage %s: not runnable (exit %d); chain stopped" % (st.name, plan.rc))
                    self.summary.append({"stage": st.name, "status": "refused", "rc": plan.rc, "notes": plan.notes})
                    rc = plan.rc
                    break
                left = self.remaining_seconds()
                need = 60.0 * (a.min_fit_minutes if st.kind.endswith("_fit") else a.min_ft_minutes)
                if self.bound is None:
                    need += 60.0 * a.bind_minutes
                if left is not None and left < need:
                    say("stage %s: %.0f s left before the job's end, %.0f s needed to start it: chain stopped before the "
                        "stage (exit %d; resubmit the same command)" % (st.name, left, need, EXIT_WALL_TIME))
                    self.summary.append({"stage": st.name, "status": "not_started_wall_time", "seconds_left": left,
                                         "seconds_needed": need})
                    rc = EXIT_WALL_TIME
                    break
                rc = self.execute(st, plan)
                if rc != 0:
                    say("stage %s: exit %d: chain stopped%s" % (st.name, rc, " (resubmit the same command to continue)"
                                                                  if rc == EXIT_WALL_TIME else ""))
                    break
            if rc == 0:
                say("chain complete: every stage of %s is complete" % self.root)
        finally:
            self.write_summary()
            self.close_backend()
            self.print_table()
            self.tee.close()
        return rc

    def close_backend(self):
        try:
            from cops_cuda import gpu_model
        except Exception:
            return
        R = gpu_model.REGISTRY
        for name in ("close", "configure", "backend_for"):
            R.__dict__.pop(name, None)
        t0 = time.perf_counter()
        n = len(R.backends)
        R.close()
        if n:
            say("closed the bound backend (%d) in %.1f s" % (n, time.perf_counter() - t0))

    def print_table(self):
        say("stage table (job %s):" % self.job)
        for r in self.summary:
            if "wall_seconds" in r:
                say("  %-6s rc %-3s segment %-3s wall %8.1f s  bind %6.1f s  reset %s  rss %.1f GB"
                    % (r["stage"], r["rc"], r["segment"], r["wall_seconds"], r.get("bind_seconds") or 0.0,
                       "%.2f s" % r["reset_basis_seconds"] if r.get("reset_basis_seconds") is not None else "-",
                       r.get("rss_gb") or float("nan")))
            else:
                say("  %-6s %s" % (r["stage"], r.get("status")))


# ================================================================================================ command line
def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", default="positive", choices=sorted(PROCESS_KINDS))
    ap.add_argument("--chain-root", required=True, help="L: the chain's directory (c<k>/fit, c<k>/ft)")
    ap.add_argument("--lambda", dest="lambda_", default=None, help="positive: one lambda, verbatim (e.g. 1e-4)")
    ap.add_argument("--ratio", default=None, choices=RATIO_DECISIONS, help="positive: the ratio decision")
    ap.add_argument("--cycles", type=int, default=3)
    ap.add_argument("--source", required=True, help="cycle-1 fit's START_RESTART (fit_restart.npz, status complete)")
    ap.add_argument("--nfev", type=int, default=4000)
    ap.add_argument("--check", type=int, default=1, choices=(0, 1))
    ap.add_argument("--c1-start-set", default="", help="NAME=VALUE[:NAME=VALUE...] for the cycle-1 fit only")
    ap.add_argument("--ft-maxnfev", type=int, required=True)
    ap.add_argument("--ft-facets", required=True, choices=("late_only", "all_seen"))
    ap.add_argument("--slope-model", default=None, help="slope: the model (SLOPE_MODELS: %s)" % ", ".join(SLOPE_MODELS))
    ap.add_argument("--slope-start", default=None, help="slope: --start of the cycle-1 fit (default: the model's)")
    ap.add_argument("--slope-bundle", default=str(RUN / "slope_mult" / "bundle" / "slope_bins.npz"))
    ap.add_argument("--bumps", default=str(RUN / "slope_bumps" / "bundle" / "bumps32.npz"), help="slope model bumps")
    ap.add_argument("--l-threads", type=int, default=None, help="slope: chunked L threads (fit and fine tuning)")
    ap.add_argument("--l-mode", default="chunked", choices=("chunked", "openblas"))
    ap.add_argument("--checkparams", default="", help="slope: CHECKPARAMS of the fit (':'-separated)")
    ap.add_argument("--ft-templates-from", default=str(RUN / "clip_refit" / "D2" / "final" / "fit_manifest.json"),
                    help="nonneg: D2's fit_manifest.json (the ft_c1..ft_c6 templates, held at 1)")
    ap.add_argument("--start-file", default=None, help="nonneg squared: the LSQ start file of nonneg_start.py")
    ap.add_argument("--slope-mult-min", default=None,
                    help="nonneg: lower bound of the free slope multipliers (nonneg_fit.py --slope-mult-min, every cycle; "
                         "user 2026-11-12: 1 for the softplus refit); absent: the fits' default (commands unchanged)")
    ap.add_argument("--ft-window-min", default=None)
    ap.add_argument("--ft-lsmr-maxiter", type=int, default=50)
    ap.add_argument("--ft-checkpoint-min", default="30")
    ap.add_argument("--ft-validate", type=int, default=1, choices=(0, 1))
    ap.add_argument("--ft-floor-ratio", default=None)
    ap.add_argument("--run", default=str(RUN))
    ap.add_argument("--point-table", default=str(RUN / "late_grid_points.txt"))
    ap.add_argument("--bundle", default=str(BUNDLE))
    ap.add_argument("--model-ext", required=True, help="the MODEL_EXT source the job's code copy was taken from")
    ap.add_argument("--job-provenance", default=None)
    ap.add_argument("--job-script", default=None)
    ap.add_argument("--devices", default=None, help="default: 0..N-1 of the job's GPUs")
    ap.add_argument("--min-fit-minutes", type=float, default=25.0,
                    help="start a fit stage only with this much time left (default 25: setup + the fit's own 15-min "
                         "stop reserve)")
    ap.add_argument("--min-ft-minutes", type=float, default=25.0,
                    help="start a fine-tuning stage only with this much time left (default 25: its own 20-min reserve)")
    ap.add_argument("--bind-minutes", type=float, default=4.0, help="added to the minimum while nothing is bound")
    ap.add_argument("--plan-only", action="store_true", help="print the stages' states (no GPU, nothing run)")
    ap.add_argument("--spawn-selftest", action="store_true", help="CPU check: the spawned workers' module setup")
    args = ap.parse_args(argv)
    if not 1 <= args.cycles <= 9:
        ap.error("--cycles must be 1..9")
    if args.kind == "positive" and (args.lambda_ is None or args.ratio is None):
        ap.error("--kind positive needs --lambda and --ratio")
    if args.kind in ("slope", "nonneg"):
        models = SLOPE_MODELS if args.kind == "slope" else NONNEG_MODELS
        if args.slope_model not in models:
            ap.error("--kind %s needs --slope-model, one of %s" % (args.kind, sorted(models)))
        if args.slope_mult_min is not None:
            if args.kind != "nonneg" or args.slope_model != "softplus":
                ap.error("--slope-mult-min: nonneg softplus only (user 2026-11-12)")
            try:
                ok = float(args.slope_mult_min) == 1.0
            except ValueError:
                ok = False
            if not ok:
                ap.error("--slope-mult-min %s: the user's decision is 1 (nonneg_problem refuses any other value)"
                         % args.slope_mult_min)
        if args.kind == "nonneg" and args.l_threads and args.l_mode != "chunked":
            ap.error("--kind nonneg: --l-mode chunked only (agent 1)")
        if args.kind == "nonneg" and args.slope_model == "squared" and not args.start_file and not args.plan_only \
                and (args.slope_start or "lsq") == "lsq":
            ap.error("--slope-model squared starts from --start-file (nonneg_start.py's LSQ start)")
        if args.ft_facets != "late_only":
            ap.error("--kind slope: the slope fine tuning has FACETS late_only only")
        if args.lambda_ is not None or args.ratio is not None or args.c1_start_set:
            ap.error("--lambda / --ratio / --c1-start-set are positive-kind options")
    if args.ft_facets == "late_only" and not args.ft_window_min:
        ap.error("--ft-facets late_only needs --ft-window-min")
    if args.ft_facets == "all_seen" and args.ft_window_min:
        ap.error("--ft-window-min applies to late_only only")
    if args.c1_start_set:
        s = args.c1_start_set.replace(",", ":")
        if re.search(r"[^A-Za-z0-9_=.:+-]", s):
            ap.error("--c1-start-set NAME=VALUE[:NAME=VALUE...]")
        args.c1_start_set = s
    return args


def check_source(args):
    """positive_cycles.sh's cycle-1 source check: status complete, the chain's lambda."""
    p = Path(args.source)
    if not p.is_file():
        return "missing cycle-1 source %s" % p
    meta = _npz_meta(p)
    status, lam = meta.get("status"), (meta.get("settings") or {}).get("lambda")
    say("cycle-1 source %s: status %s, lambda %s, final cost %s" % (p, status, lam,
                                                                   (meta.get("optimizer") or {}).get("cost")))
    if status != "complete":
        return "the cycle-1 source has status %s, not complete" % status
    if lam is None or float(lam) != float(args.lambda_):
        return "the cycle-1 source has lambda %s, not %s" % (lam, args.lambda_)
    return None


def _spawn_probe(queue):
    """--spawn-selftest: what a spawned GPU worker sees after worker_setup() (this file as __mp_main__)."""
    import positive_activity
    from cops_cuda import transport
    h = transport.WorkerState.handle
    reply = h(transport.WorkerState(0), (DETACH_REQUEST,))
    queue.put({"pid": os.getpid(), "main": sys.modules["__main__"].__file__,
               "mp_main": getattr(sys.modules.get("__mp_main__"), "__name__", None),
               "positive_activity": {k: positive_activity._STATE.get(k) for k in ("installed", "role")},
               "handler_chain_onejob": bool(getattr(h, "_chain_onejob", False)),
               "detach_reply": reply, "kind_env": os.environ.get(KIND_ENV),
               "frozen": os.environ.get("COPS_POS_FROZEN_CONFIG"), "positive_fit": "positive_fit" in sys.modules,
               "positive_finetune": "positive_finetune" in sys.modules, "argv": list(sys.argv)})


def spawn_selftest(args):
    """CPU check of the worker import path: parent setup as in run(), then one spawned process (no GPU)."""
    import multiprocessing as mp
    os.environ[KIND_ENV] = args.kind
    PROCESS_KINDS[args.kind]["parent"](args)
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    proc = ctx.Process(target=_spawn_probe, args=(q,))
    proc.start()
    info = q.get(timeout=300)
    proc.join(60)
    print(json.dumps(info, indent=1, default=str))
    ok = (info["positive_activity"] == {"installed": True, "role": "worker"} and info["handler_chain_onejob"]
          and info["detach_reply"] == ("ok", {"attached": False, "pid": info["pid"]}) and info["positive_fit"]
          and info["positive_finetune"] and info["mp_main"] == "__mp_main__" and proc.exitcode == 0)
    print("SPAWN SELFTEST %s" % ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


def check_source_slope(args):
    """slope_cycles.sh: the cycle-1 source must exist (the fit's start check compares it further)."""
    p = Path(args.source)
    if not p.is_file():
        return "missing cycle-1 source %s" % p
    meta = _npz_meta(p)
    say("cycle-1 source %s: status %s, final cost %s" % (p, meta.get("status"), (meta.get("optimizer") or {}).get("cost")))
    if meta.get("status") != "complete":
        return "the cycle-1 source has status %s, not complete" % meta.get("status")
    return None


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.spawn_selftest:
        return spawn_selftest(args)
    chain = Chain(args)
    if args.plan_only:
        return chain.plan_only()
    c1 = chain.stages[0]
    if isinstance(c1, PositiveFitStage) and c1.state()[1] == "none":
        msg = check_source(args) if args.kind == "positive" else check_source_slope(args)
        if msg:
            print("chain_onejob: %s" % msg, file=sys.stderr)
            return 2
    chain.root.mkdir(parents=True, exist_ok=True)
    return chain.run()


if __name__ == "__mp_main__":                  # a spawned GPU worker (multiprocessing re-imports the main module)
    worker_setup()

if __name__ == "__main__":
    sys.exit(main())
