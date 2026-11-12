# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_mult.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_mult.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e6127bf5b9baa5d7bd919b7755755e798237c580cb53fd5ffc3cc1c9a2f06a03
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Slope-binned H2O activity multipliers with an activity floor (user decisions 2026-11-10: A = max(1e-5
#   A_orig, m_b L), one multiplier >= 0 per slope bin, 0-40 deg held at 1): model math, CPU stand-in and the GPU patch
#   of the host-activity path (kernels unchanged; installed only by slope_mult_fit.py / slope_finetune.py); bump basis
#   and held term exp(h) (option 3a), chunked bitwise L threads, the fine tuning's worker requests. This version
#   (seen-path fix 2026-11-10 21:12 UTC) is in snapshots slopemul_20261110T211229Z (v7) and slopemul_20261110T221713Z
#   (v8).
# USED BY: slope_mult_fit.py, slope_finetune.py; helpers imported by nonneg_mult.py. v8: the 3a chain from the c1 fine
#   tuning on (124064-124086); the c1 fit 124008 ran v6 and the harmonic slope-multiplier fit 123786 an earlier version
#   (snapshot slopemul_20261110T094013Z).
# ==== end of header ====
"""Slope-binned H2O activity multipliers with an activity floor (user decisions 2026-11-10): model, CPU stand-in and the
GPU patch of the host-activity path. The fit-level part is slope_mult_problem.py, the driver slope_mult_fit.py, the
template bundle (bin indicators) slope_bins.py.

Local, untracked model extension: nothing in the git tree is edited and no existing model_ext file is changed. The
patch below is installed ONLY by slope_mult_fit.py (install_gpu); every other driver runs the unchanged evaluator.

User decisions (2026-11-10, the lead's brief)
* "define multiple cliff templates (0-40, 40-50, 50-60, 60-70, 70-80, and 80-90), associate each of these templates
  with an activity multiplier that will be multiplied by the activity from spherical harmonics to get the actual
  surface activity; each of these multipliers will be fit independently; for the overall activity, clip the activity
  not by 0 but by 1E-5".
* Answers: floor = 1e-5 x the facet's original activity A_orig (elements.bin column 32 "AH2O" of the base = the run's
  srcCG activity; A_orig = 0 on 0.7 % of the area -> floor 0 there, accepted); slopes above 90 deg form a 7th bin
  90-180; multipliers >= 0, no upper bound; the 0-40 deg multiplier is held at 1 (otherwise a common factor is
  degenerate with the harmonic coefficients).
* Bins by the effective gravitational slope (reference_slope_deg of runs/shap8_012m_flat_earlier/cache/
  admitted_gravity_cache.npz, same 12M mesh): [0,40), [40,50), [50,60), [60,70), [70,80), [80,90), [90,180].

Model, per facet e with slope bin b(e):
    L(e)     = sum_j a_j Y_j(e)                         (the 25 degree-<=4 harmonics of fit_cops, signed, linear)
    s(e)     = m_b(e) * L(e)
    A_H2O(e) = max(F(e), s(e)),  F(e) = f * A_orig(e),  f = 1e-5
    m_{slope_00_40} = 1 (held, fixed_parameters), m_b >= 0 free for the other 6 bins
The flux and pressures are the unchanged GPU model with af = A_H2O (q(SZA) A, srcCG night/shadow ratios 0.2 H2O /
0.1 CO2, no flux floor, CO2 = kappa x the run's CO2 map, no CO2 floor). Pressures are linear in af.

Derivatives (the fitter's one-sided clip convention): active(e) = s(e) > F(e) (strict);
    dA/da_j = active * m_b(e) * Y_j(e),   dA/dm_k = active * [b(e) = k] * L(e).
A facet AT the floor (s <= F, equality included) contributes 0 to every activity/multiplier column (the kink is
not differentiable; the one-sided derivative from above is used only where s > F). Consequence: where F = 0
(A_orig = 0) and m_k = 0, dA/dm_k = 0 even where L > 0 (a multiplier at its bound 0 sees no gradient from such
facets). Temperature and kappa columns are the unchanged kernel's (computed with af = A_H2O, floor included).

GPU implementation (host-activity path; kernels unchanged; why: the activity is a host field, the kernels already
compute every column as sum_r w_r active(f_r) colbasis[c][f_r] with a host-uploaded colbasis and active set)
* The fitter's template basis is the 7 bin indicators 1[b(e) = k] (slope_bins.py bundle; template amplitude c_<bin> =
  multiplier m_k, lower bound 0 as for every template). The bound basis therefore carries the bin of every facet to
  every GPU process (the bind already ships it; 7 x 12.4M doubles, as D2's 7 templates).
* act_lin: CudaBackend.activity_linear / MultiGpuBackend.activity_linear return L = Y a (the original's first line,
  bitwise) instead of Y a + A_scale B b when b is given; this host field is what the local workers read from shared
  memory (and what the MPI transport broadcasts with --activity host).
* CudaBackend._upload_activity (every process that owns a GPU): s = m[bin] * L, af = max(F, s), active = s > F on the
  device (IEEE multiply / max / compare: the same bits as the host expressions; checked against the host every
  CHECK_EVERY-th evaluation and at the first one).
* CudaBackend._upload_columns: colbasis rows m[bin] * Y_j (requested harmonics) and where(bin = k, L, 0) (requested
  multipliers), formed on the device from a resident Y^T of the device facets (no per-Jacobian host upload).
* b = None (the fitter's start-up consistency check at the run's coefficients, any caller without templates): the
  ORIGINAL methods run unchanged (max(0, act_lin)); at the run's coefficients the model with m = 1 equals that map.
* device activity (--activity device, kernels/activity.cu) cannot form this field: refused when b is given.
* enhancement_scale (A_scale) is not used by this model (multipliers are dimensionless; the columns are dA/dm).

Processes (what each needs; the driver does all of it at start-up)
* every process that evaluates (the driver; each local spawned worker; each MPI rank) runs install_gpu(cfg) before
  its first bind. The driver: slope_mult_fit.py calls install_gpu at start-up (on every MPI rank too, from its own
  command line); local spawned workers: transport.pipe_worker_main is replaced by _pipe_worker_main, which installs
  from the inherited environment variable COPS_SLOPEMUL_FROZEN_CONFIG (the driver's canonical config) before the
  original worker loop binds. A process whose own config differs from an inherited frozen one refuses to start.
* files read by each GPU process: the base's elements.bin (AH2O, via the dataset loader, digest checked against the
  config) and nothing else; the bins come with the bound basis (digest checked against the config).
* refused together with COPS_POS_* (positive_activity.py patches the same methods).

Python API: config_from_bundle_meta(...), install_gpu(cfg), installed_config(), bin_index(slope), activity(...),
cpu_model_class(fit_cops, floor_ratio), canonical(cfg).
"""

import hashlib
import inspect
import json
import math
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# ------------------------------------------------------------------------------------------------ user decisions
BIN_EDGES_DEG = (0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 180.0)
BIN_NAMES = ("slope_00_40", "slope_40_50", "slope_50_60", "slope_60_70", "slope_70_80", "slope_80_90", "slope_90_180")
BIN_LABELS = ("[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]")
FIXED_BIN = "slope_00_40"           # held at 1 (the lead, 2026-11-10: a common factor is degenerate with a)
FIXED_VALUE = 1.0
FLOOR_RATIO = 1e-5                  # F = 1e-5 * A_orig (user 2026-11-10)
TEST_FLOOR_RATIOS = (0.0,)          # GPU check (a) only: floor 0 with all multipliers 1 = the plain clipped model
LEAD_COUNTS = (8851632, 1163644, 887083, 654731, 402863, 217406, 186931)   # the lead's bin counts (12M)
CLIFF60_FACETS = 1461931            # cliff_ge_60deg of runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
POLICY = ("user 2026-11-10: H2O activity A = max(1e-5 A_orig, m_b L), L = the 25 signed degree-<=4 harmonics (not "
          "exponentiated), m_b one free multiplier >= 0 (no upper bound) per effective-gravitational-slope bin [40,50), "
          "[50,60), [60,70), [70,80), [80,90), [90,180] deg, m_[0,40) = 1 held; floor F = 1e-5 x A_orig (elements.bin "
          "AH2O; F = 0 where A_orig = 0, accepted); flux q(SZA) A unchanged (srcCG night/shadow ratios, no flux floor, "
          "CO2 = kappa x the run's CO2 map, no CO2 floor); no cliff_ge_60deg and no ft_c* templates")
ENV_FROZEN = "COPS_SLOPEMUL_FROZEN_CONFIG"
CONFIG_VERSION = 1
CHECK_EVERY = 100                   # device af/active vs the host expression: first evaluation and every 100th
BYTES_PER_FACE = 240                # device scratch per facet: Y^T (25 x 8), L, mrow, s, F (4 x 8), bins (4), bool (1)
# Tracked sources the patch is written against (AMPS master; unchanged since 47fccd5): refused if any changed.
ORIGINAL_SHA256 = {
    "evaluator.py": "2f1e9d5fadbc4111e3e5222d32f8098e51e8fd01c1e85f5266730689ec53d52e",
    "multi_gpu.py": "1c3a9c1a377784f680a1dcb7a8f58e0c58f925ced241529d8f4c89a07dd516fc",
    "transport.py": "605ed2546625b0e029d277fac0757635ff978081b24526def324d18cb0256304",
    "backend.py": "84b54395fb4f590859f6648ade5ea6e356229653e8482795a2b96ef2127a5c89",
    "gpu_model.py": "c1659f5cf91b1f65f1222c3056f6bf230f07100120fd18666d9ada2e90620a50",
}
LOG = "slope_mult:"

_STATE = {"installed": False, "config": None, "originals": {}, "role": None, "backends": 0,
          "L_timing": {"calls": 0, "seconds": 0.0, "max_s": 0.0, "first_s": None, "threads": None},
          "L_threads": None, "L_lib": None, "L_mode": None, "L_pool": None, "L_fallbacks": 0}


def say(text):
    print("%s %s" % (LOG, text), flush=True)


def sha256_bytes(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def module_sha256():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def canonical(cfg):
    return json.dumps(cfg, sort_keys=True, separators=(",", ":"))


# ------------------------------------------------------------------------------------------------ BLAS threads
# User 2026-11-10 ("can CPU calculation be parallelized with OpenMP?" -> "run the threaded test fit"): the FITTER
# process may raise the thread count of numpy's bundled OpenBLAS (numpy 1.24.4: numpy.libs/libopenblas64_*.so, ILP64,
# symbols openblas_{set,get}_num_threads64_; threadpoolctl is not installed) at run time; a departure from CLAUDE.md
# section 4 (OMP/OPENBLAS/MKL_NUM_THREADS=1) for that process only. Scope (set_L_blas_threads): only the L = Y a
# matvec (the fitter's dominant CPU step, ~0.29 s of ~0.35 s per evaluation single-threaded) runs with N threads; the
# count is restored right after, so every other BLAS call of the fitter (scipy TRF's small products, the final
# Gauss-Newton algebra) is single-threaded as before. The environment is NOT changed, so the GPU worker
# processes (spawned; they initialise OpenBLAS from OPENBLAS_NUM_THREADS=1) keep 1 thread; each backend records its
# process's count (blas_info) in its identity. scipy's own OpenBLAS (scipy.libs) is not touched.
def _numpy_openblas():
    """(path, ctypes library) of the OpenBLAS numpy has loaded in this process, or (None, None)."""
    import ctypes
    import numpy  # noqa: F401  (loads its OpenBLAS)
    path = None
    try:
        with open("/proc/self/maps") as f:
            for line in f:
                if "libopenblas64_" in line and "/numpy.libs/" in line:
                    path = line.split()[-1]
                    break
    except OSError:
        return None, None
    if path is None:
        return None, None
    lib = ctypes.CDLL(path)                      # the already loaded library (same file: same handle)
    if not hasattr(lib, "openblas_get_num_threads64_"):
        return path, None
    return path, lib


def blas_info():
    """{"library", "threads", "config", "pid", "host"} of numpy's OpenBLAS in this process (threads None if unknown)."""
    import ctypes
    path, lib = _numpy_openblas()
    out = {"library": path, "threads": None, "config": None, "pid": os.getpid(), "host": socket.gethostname(),
           "env": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}}
    if lib is not None:
        out["threads"] = int(lib.openblas_get_num_threads64_())
        try:
            lib.openblas_get_config64_.restype = ctypes.c_char_p
            out["config"] = lib.openblas_get_config64_().decode()
        except Exception:                        # recorded as None, never fatal
            pass
    return out


def set_blas_threads(n):
    """Set numpy's OpenBLAS thread count of THIS process to n (1..MAX_THREADS); returns (before, after) blas_info();
    refuses if the library cannot be found or the count read back differs."""
    n = int(n)
    if n < 1:
        raise ValueError("BLAS threads must be >= 1")
    before = blas_info()
    path, lib = _numpy_openblas()
    if lib is None:
        raise RuntimeError("slope_mult: numpy's OpenBLAS (openblas_set_num_threads64_) not found in this process")
    lib.openblas_set_num_threads64_(n)
    after = blas_info()
    if after["threads"] != n:
        raise RuntimeError("slope_mult: openblas_set_num_threads64_(%d) read back %r (library %s, config %s)"
                           % (n, after["threads"], path, after["config"]))
    return before, after


def set_L_blas_threads(n):
    """Scope of --blas-threads: numpy's OpenBLAS runs L = Y a of THIS process (the fitter: MultiGpuBackend /
    CudaBackend.activity_linear) with n threads and is set back to its previous count right after; every other BLAS
    call of the process keeps the environment's count. Checked once (set, read back, restore). Returns the record."""
    n = int(n)
    if n < 1:
        raise ValueError("BLAS threads must be >= 1")
    path, lib = _numpy_openblas()
    if lib is None:
        raise RuntimeError("slope_mult: numpy's OpenBLAS (openblas_set_num_threads64_) not found in this process")
    before = int(lib.openblas_get_num_threads64_())
    lib.openblas_set_num_threads64_(n)
    got = int(lib.openblas_get_num_threads64_())
    lib.openblas_set_num_threads64_(before)
    back = int(lib.openblas_get_num_threads64_())
    if got != n or back != before:
        raise RuntimeError("slope_mult: openblas_set_num_threads64_(%d) read back %d (restore %d -> %d)"
                           % (n, got, before, back))
    _STATE["L_threads"], _STATE["L_lib"] = n, lib
    _STATE["L_timing"]["threads"] = n
    return {"scope": "L = Y a (activity_linear) of this process; other BLAS calls keep %d thread(s)" % before,
            "threads_L": n, "threads_outside_L": before, "library": path, "check_set_read_back": got,
            "check_restored": back}


# User decision 2026-11-10 ("Bitwise chunked version"): the default way to run L = Y a with N threads is a pool of N
# Python threads, each Y[c0:c1] @ a with numpy's OpenBLAS at ONE thread (numpy releases the GIL inside the BLAS call),
# chunk starts multiples of L_ALIGN_ROWS rows, all chunks written into one output array allocated per call. Every row is
# then computed by the same single-thread kernel path as in the 1-thread call (the n mod 4 tail stays the tail of the
# last chunk), so L is bitwise the 1-thread L (checked: CPU job 123976 on the 12M arrays, 8-64 threads; test_slope_mult
# CPU unit test). OpenBLAS stays at 1 thread everywhere. The threaded-OpenBLAS mode (set_L_blas_threads) is kept only as
# the explicitly named, non-default mode "openblas" (NOT bitwise: 2-4 elements of L differ at thread-partition
# boundaries and the kernel tail).
L_ALIGN_ROWS = 4096
L_MODES = ("chunked", "openblas")


def chunk_bounds(nrows, nthreads, align=L_ALIGN_ROWS):
    """[(c0, c1), ...]: at most nthreads chunks, every start a multiple of align, the last ending at nrows. A last chunk
    shorter than align rows is merged into the previous one: a 1-row chunk would go through numpy's dot special case
    ((1, k) @ (k,) -> ddot) instead of dgemv and differ from the 1-thread call (found by test_slope_mult, n = 4097)."""
    nrows, nthreads, align = int(nrows), int(nthreads), int(align)
    if nrows <= 0:
        return []
    step = -(-nrows // max(nthreads, 1))
    step = -(-step // align) * align
    out = [[c0, min(c0 + step, nrows)] for c0 in range(0, nrows, step)]
    if len(out) > 1 and out[-1][1] - out[-1][0] < align:
        last = out.pop()
        out[-1][1] = last[1]
    return [tuple(c) for c in out]


def chunked_matvec(Y, a, pool, nthreads, align=L_ALIGN_ROWS):
    """Y @ a by row chunks on `pool` (each np.matmul(Y[c0:c1], a, out=out[c0:c1])); Y C-contiguous float64 (n, k)."""
    out = np.empty(Y.shape[0], dtype=np.result_type(Y.dtype, np.asarray(a).dtype))
    futures = [pool.submit(np.matmul, Y[c0:c1], a, out=out[c0:c1]) for c0, c1 in chunk_bounds(len(Y), nthreads, align)]
    for f in futures:
        f.result()
    return out


def set_L_threads(n, mode="chunked"):
    """--l-threads N --l-mode MODE: L = Y a of THIS process (the fitter) with N threads. mode "chunked" (default, user
    2026-11-10: bitwise): a persistent pool of N Python threads, OpenBLAS at 1 thread (checked); mode "openblas"
    (explicit, non-default, NOT bitwise): set_L_blas_threads. Returns the record."""
    n = int(n)
    if n < 1:
        raise ValueError("L threads must be >= 1")
    if mode not in L_MODES:
        raise ValueError("L mode must be one of %s" % (L_MODES,))
    if mode == "openblas":
        rec = set_L_blas_threads(n)
        _STATE["L_mode"] = "openblas"
        rec.update(mode="openblas", threads=n, bitwise_by_construction=False)
        return rec
    info = blas_info()
    if info["threads"] != 1:
        raise RuntimeError("slope_mult: chunked L needs numpy's OpenBLAS at 1 thread in this process (got %r; set "
                           "OPENBLAS_NUM_THREADS=1)" % info["threads"])
    from concurrent.futures import ThreadPoolExecutor
    _STATE["L_pool"] = ThreadPoolExecutor(max_workers=n, thread_name_prefix="slopemul_L")
    _STATE["L_threads"], _STATE["L_mode"], _STATE["L_lib"] = n, "chunked", None
    _STATE["L_timing"]["threads"] = n
    return {"mode": "chunked", "threads": n, "align_rows": L_ALIGN_ROWS, "openblas_threads": info["threads"],
            "library": info["library"], "bitwise_by_construction": True,
            "scope": "L = Y a (activity_linear) of this process: %d Python threads, each Y[chunk] @ a with OpenBLAS at 1 "
                     "thread, chunk starts multiples of %d rows, one output array per call" % (n, L_ALIGN_ROWS)}


def _matvec_L(Y, a):
    """Y @ a; with --l-threads: chunked (default mode) or threaded OpenBLAS (mode openblas, restored after the call)."""
    mode, n = _STATE["L_mode"], _STATE["L_threads"]
    if mode is None or n is None:
        return Y @ a
    if mode == "chunked":
        if not (isinstance(Y, np.ndarray) and Y.ndim == 2 and Y.flags.c_contiguous and Y.dtype == np.float64):
            _STATE["L_fallbacks"] += 1               # never expected for the fitter's Y; recorded in the timing
            return Y @ a
        return chunked_matvec(Y, np.ascontiguousarray(a, dtype=np.float64), _STATE["L_pool"], n)
    lib = _STATE["L_lib"]
    prev = int(lib.openblas_get_num_threads64_())
    lib.openblas_set_num_threads64_(n)
    try:
        return Y @ a
    finally:
        lib.openblas_set_num_threads64_(prev)


def activity_linear_timing():
    return dict(_STATE["L_timing"], mode=_STATE["L_mode"], fallbacks=_STATE["L_fallbacks"])


def _time_L(seconds):
    t = _STATE["L_timing"]
    t["calls"] += 1
    t["seconds"] += seconds
    t["max_s"] = max(t["max_s"], seconds)
    if t["first_s"] is None:
        t["first_s"] = seconds


# ------------------------------------------------------------------------------------------------ model math
def bin_index(slope):
    """Bin of each facet (int8, 0..6) for the effective gravitational slope [deg]: [0,40) -> 0, [40,50) -> 1, ...,
    [80,90) -> 5, [90,180] -> 6. Refuses non-finite slopes and slopes outside [0, 180]."""
    s = np.asarray(slope, dtype=np.float64)
    if not np.all(np.isfinite(s)):
        raise ValueError("non-finite slope values")
    if s.size and (float(s.min()) < BIN_EDGES_DEG[0] or float(s.max()) > BIN_EDGES_DEG[-1]):
        raise ValueError("slope outside [%g, %g] deg: min %r max %r" % (BIN_EDGES_DEG[0], BIN_EDGES_DEG[-1],
                                                                       float(s.min()), float(s.max())))
    inner = np.asarray(BIN_EDGES_DEG[1:-1], dtype=np.float64)
    return np.searchsorted(inner, s, side="right").astype(np.int8)


def indicator_basis(bins, nbins=len(BIN_NAMES)):
    """The template matrix (nfaces, 7): column k = 1[bin = k] (float64)."""
    bins = np.asarray(bins)
    out = np.zeros((len(bins), nbins), dtype=np.float64)
    out[np.arange(len(bins)), bins.astype(np.int64)] = 1.0
    return out


def bins_from_basis(basis, nbins=len(BIN_NAMES)):
    """Bin index (int8) from the bound indicator basis; refuses anything that is not one-hot 0/1 with nbins columns."""
    B = np.asarray(basis)
    if B.ndim != 2 or B.shape[1] != nbins:
        raise ValueError("slope_mult: the template basis has shape %s; the model needs (nfaces, %d) bin indicators"
                         % (B.shape, nbins))
    if not np.all((B == 0.0) | (B == 1.0)):
        raise ValueError("slope_mult: the template basis is not 0/1 (not the bin indicators)")
    if not np.all(np.sum(B, axis=1) == 1.0):
        raise ValueError("slope_mult: a facet is in none or several bins (the basis is not one-hot)")
    return np.argmax(B, axis=1).astype(np.int8)


def activity(L, m, bins, floor, g=None):
    """(s, A, active): s = m[bin] L (times g = exp(h) with a held term), A = max(F, s), active = s > F (the host
    expressions; the device forms the same, in the same order)."""
    mrow = np.asarray(m, dtype=np.float64)[np.asarray(bins).astype(np.intp)]
    s = mrow * np.asarray(L, dtype=np.float64)
    if g is not None:
        s = s * g
    A = np.maximum(floor, s)
    active = s > floor
    return s, A, active


def bump_L(Phi_c, w, scale, pool=None, nthreads=None):
    """L = scale * (Phi w) (Phi C-contiguous (n, K)); chunked like the fitter's L when a pool is given (bitwise)."""
    w = np.ascontiguousarray(w, dtype=np.float64)
    v = chunked_matvec(Phi_c, w, pool, nthreads) if pool is not None else Phi_c @ w
    return scale * v


# ------------------------------------------------------------------------------------------------ configuration
def config_from_bundle_meta(meta, floor_ratio, bundle_path=None):
    """The model configuration every process must share, from the slope bundle's metadata (slope_bins.py writes
    ``slope_multiplier_bins``)."""
    sm = (meta or {}).get("slope_multiplier_bins")
    if not sm:
        raise ValueError("slope_mult: %s is not a slope_bins.py bundle (no slope_multiplier_bins metadata)" % bundle_path)
    if list(sm.get("names") or []) != list(BIN_NAMES) or [float(v) for v in sm.get("edges_deg")] != list(BIN_EDGES_DEG):
        raise ValueError("slope_mult: bundle bins %s / %s are not this model's %s / %s"
                         % (sm.get("names"), sm.get("edges_deg"), BIN_NAMES, BIN_EDGES_DEG))
    f = float(floor_ratio)
    if not (f == FLOOR_RATIO or f in TEST_FLOOR_RATIOS):
        raise ValueError("slope_mult: floor ratio %r is not a decision (user 2026-11-10: %r; test only: %s)"
                         % (f, FLOOR_RATIO, TEST_FLOOR_RATIOS))
    return {"version": CONFIG_VERSION, "model": "slope_multiplier_floor", "floor_ratio": f,
            "test_only": f != FLOOR_RATIO, "bin_names": list(BIN_NAMES), "bin_edges_deg": list(BIN_EDGES_DEG),
            "fixed_bin": FIXED_BIN, "nfaces": int(sm["nfaces"]), "bins_sha256": sm["bins_sha256"],
            "bin_counts": [int(v) for v in sm["counts"]], "ah2o_sha256": sm["ah2o_sha256"],
            "mesh_identity": sm["mesh_identity"]}


def add_bumps_to_config(cfg, meta, bundle_path=None):
    """Bump + slope-multiplier model (user 2026-11-10, option 3a): L = A_scale * (Phi w), Phi = the bump bundle's
    columns (slope_bumps.py metadata "slope_bumps"); the fitter's templates are [7 bin indicators, K bumps]."""
    sb = (meta or {}).get("slope_bumps")
    if not sb:
        raise ValueError("slope_mult: %s is not a slope_bumps.py bundle" % bundle_path)
    if sb["mesh_identity"] != cfg["mesh_identity"] or int(sb["nfaces"]) != int(cfg["nfaces"]):
        raise ValueError("slope_mult: the bump bundle's mesh differs from the bins bundle's")
    out = dict(cfg)
    out.update(basis="bumps", n_bumps=len(sb["names"]), bump_names=list(sb["names"]), bumps_sha256=sb["bumps_sha256"],
               model="bump_slope_multiplier_floor")
    return out


def add_held_to_config(cfg, path, key="held_next"):
    """Held per-facet term h of the fine-tuning cycles: A = max(F, m_b L exp(h)); h from <path>[key] (float64, one per
    facet); every GPU process loads it and checks the digest."""
    h = load_held(path, key, cfg["nfaces"])
    out = dict(cfg)
    out["held"] = {"path": str(Path(path).resolve()), "key": key, "sha256": sha256_bytes(h),
                   "nonzero_facets": int(np.count_nonzero(h)), "min": float(h.min()), "max": float(h.max())}
    return out


def load_held(path, key, nfaces):
    with np.load(str(path), allow_pickle=False) as z:
        h = np.ascontiguousarray(np.asarray(z[key], dtype=np.float64))
    if h.shape != (int(nfaces),) or not np.all(np.isfinite(h)):
        raise ValueError("slope_mult: held term %s[%s] has shape %s / non-finite values" % (path, key, h.shape))
    return h


def held_factor(cfg):
    """(g = exp(h), h) of the configured held term, or (None, None); cached per process; digest checked."""
    held = (cfg or {}).get("held")
    if not held:
        return None, None
    c = _STATE.setdefault("held_cache", {})
    if held["sha256"] not in c:
        h = load_held(held["path"], held["key"], cfg["nfaces"])
        if sha256_bytes(h) != held["sha256"]:
            raise ValueError("slope_mult: held term %s changed (sha256)" % held["path"])
        c[held["sha256"]] = (np.exp(h), h)
    return c[held["sha256"]]


def n_bins_bumps(cfg):
    return len(cfg["bin_names"]), int(cfg.get("n_bumps", 0))


def installed_config():
    return dict(_STATE["config"]) if _STATE["config"] else None


# ------------------------------------------------------------------------------------------------ CPU stand-in
def cpu_model_class(fit_cops, floor_ratio, nbumps=0, held_g=None):
    """fit_cops.Model with this model's activity and columns (CPU; the tests' stand-in for the GPU model).
    b None -> the reference evaluate unchanged. Expressions are the reference Model.evaluate's with af and the masked
    column bases replaced. nbumps > 0: the bump model (L = scale (Phi w), Phi = basis[:, 7:], a held at 0, a columns
    exact zeros); held_g: the factor exp(h) (None: 1); set_override(facets, g) replaces it at facets (fine tuning);
    sensitivities(T, facets): dp/dA_e of every loaded point."""
    base = fit_cops.Model
    f = float(floor_ratio)
    nb = len(BIN_NAMES)

    class SlopeMultiplierCpuModel(base):
        def __init__(self, run, header, cols, rows, elements, coef):
            self._sm_floor = f * np.asarray(elements["AH2O"], dtype=np.float64)
            self._sm_g = None if held_g is None else np.asarray(held_g, dtype=np.float64)
            self._sm_override = None
            super().__init__(run, header, cols, rows, elements, coef)

        def set_override(self, facets, g_values):
            if facets is None:
                self._sm_override = None
                return 0
            g = np.ones(self.nfaces) if self._sm_g is None else self._sm_g.copy()
            g[np.asarray(facets, np.int64)] = np.asarray(g_values, np.float64)
            self._sm_override = g
            return len(facets)

        def _g(self):
            return self._sm_override if self._sm_override is not None else self._sm_g

        def _L(self, a, b, enhancement_basis, enhancement_scale):
            if nbumps:
                Phi = np.ascontiguousarray(np.asarray(enhancement_basis, np.float64)[:, nb:nb + nbumps])
                return enhancement_scale * (Phi @ np.asarray(b, np.float64)[nb:nb + nbumps]), Phi
            return self.Y @ a, None

        def _temps(self, T):
            KT = fit_cops.KT
            Tmax, T75, Tmin = T
            Tl = Tmax + KT * (Tmax - T75) * self.xlit
            use_lit = self.lit & (Tl > Tmin)
            Te = np.where(use_lit, Tl, Tmin)
            return use_lit, Te, 1.0 / np.sqrt(Te)

        def sensitivities(self, T, facets):
            _, _, u = self._temps(T)
            cNG, cRG, bm = self.cNG, self.cRG, self.bm
            Bm = np.zeros((self.nfaces, len(facets)))
            Bm[np.asarray(facets, np.int64), np.arange(len(facets))] = 1.0
            ng = cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ Bm)
            rg = cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ Bm)
            return ng, rg

        def evaluate(self, a, T, kappa, b=None, enhancement_basis=None, enhancement_scale=1.0, jac=False,
                     basis_index=None, enhancement_index=None):
            if b is None:
                return base.evaluate(self, a, T, kappa, b=None, enhancement_basis=enhancement_basis,
                                     enhancement_scale=enhancement_scale, jac=jac, basis_index=basis_index,
                                     enhancement_index=enhancement_index)
            KT = fit_cops.KT
            bins = bins_from_basis(np.asarray(enhancement_basis)[:, :nb])
            bvec = np.asarray(b, dtype=np.float64)
            m = bvec[:nb]
            L, Phi = self._L(a, bvec, enhancement_basis, enhancement_scale)
            g = self._g()
            s, A, active_f = activity(L, m, bins, self._sm_floor, g)
            mrow = m[bins.astype(np.intp)]
            use_lit, Te, u = self._temps(T)
            af = A[self.face]
            bc = lambda w: np.bincount(self.pt, weights=w, minlength=self.npt)     # noqa: E731
            cNG, cRG, bm = self.cNG, self.cRG, self.bm
            pNG_H = cNG["H2O"] * bm["H2O"] * bc(self.wn["H2O"] * af * u)
            pRG_H = cRG["H2O"] * (bc(self.wp["H2O"] * af) + bm["H2O"] * bc(self.wv["H2O"] * af * u))
            pNG_C = cNG["CO2"] * bm["CO2"] * bc(self.wn["CO2"] * u)
            pRG_C = cRG["CO2"] * (self.PhiCO2 + bm["CO2"] * bc(self.wv["CO2"] * u))
            pNG = pNG_H + kappa * pNG_C
            pRG = pRG_H + kappa * pRG_C
            if not jac:
                return pNG, pRG
            J = {}
            active = active_f[:, None]
            gg = np.ones(len(L)) if g is None else g

            def activity_columns(Bm):
                masked = active * Bm
                JaNG = cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ masked)
                JaRG = cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ masked)
                return JaNG, JaRG
            basis_index = [] if basis_index is None else list(basis_index)
            if basis_index:
                if nbumps:
                    for j in basis_index:
                        J[("a", j)] = (np.zeros(self.npt), np.zeros(self.npt))
                else:
                    JaNG, JaRG = activity_columns((mrow * gg)[:, None] * self.Y[:, basis_index])
                    for k, j in enumerate(basis_index):
                        J[("a", j)] = (JaNG[:, k], JaRG[:, k])
            enhancement_index = [] if enhancement_index is None else list(enhancement_index)
            if enhancement_index:
                cols = []
                for k in enhancement_index:
                    if k < nb:
                        cols.append(np.where(bins == k, L * gg, 0.0))
                    else:
                        cols.append(enhancement_scale * Phi[:, k - nb] * mrow * gg)
                JbNG, JbRG = activity_columns(np.column_stack(cols))
                for k, j in enumerate(enhancement_index):
                    J[("b", j)] = (JbNG[:, k], JbRG[:, k])
            du = -0.5 * u / Te
            dl, dn_ = np.where(use_lit, du, 0.0), np.where(use_lit, 0.0, du)
            dlx = dl * self.xlit

            def sums(w):
                s_l, s_lx, s_n = bc(w * dl), bc(w * dlx), bc(w * dn_)
                return KT * s_lx + s_l, -KT * s_lx, s_n
            gN = cNG["H2O"] * bm["H2O"] * self.wn["H2O"] * af + kappa * cNG["CO2"] * bm["CO2"] * self.wn["CO2"]
            tN = sums(gN)
            tRH = sums(bm["H2O"] * self.wv["H2O"] * af)
            tRC = sums(bm["CO2"] * self.wv["CO2"])
            for k in range(3):
                J[("T", k)] = (tN[k], cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k])
            J[("kappa", 0)] = (pNG_C, pRG_C)
            return pNG, pRG, J

    SlopeMultiplierCpuModel.slope_multiplier_floor_ratio = f
    SlopeMultiplierCpuModel.slope_multiplier_nbumps = int(nbumps)
    return SlopeMultiplierCpuModel


# ------------------------------------------------------------------------------------------------ GPU patch
def _check_sources(cops_dir):
    bad = {}
    for name, want in ORIGINAL_SHA256.items():
        got = hashlib.sha256((Path(cops_dir) / name).read_bytes()).hexdigest()
        if got != want:
            bad[name] = got
    if bad:
        raise RuntimeError("slope_mult: tracked sources changed since the patch was written (%s); the patched methods "
                           "would no longer match: review slope_mult.py" % bad)


# the worker requests of the fine tuning: request kind -> the CudaBackend method (installed by install_gpu) that answers
# it in the worker (ws_handle); test_slope_mult.py --test worker_handlers checks the installed names and the dispatch
WORKER_REQUESTS = {"slopemul_override": "slopemul_set_override", "slopemul_seen": "slopemul_seen_rows",
                   "slopemul_sens": "slopemul_sensitivities"}

_SELECT = {}
_SEEN_SRC = r"""
extern "C" __global__ void slopemul_seen(const long long* __restrict__ blocks, long long nblocks,
    const int* __restrict__ face, const float* __restrict__ kn, const float* __restrict__ kphi,
    const float* __restrict__ kv, const unsigned char* __restrict__ mask, unsigned char* __restrict__ seen,
    int structural)
{
    for (long long b = blockIdx.x; b < nblocks; b += gridDim.x) {
        const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
        if (!mask[p]) continue;
        for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
            const float a = kn[r], c = kphi[r], d = kv[r];
            const bool on = structural ? (a != 0.0f || c != 0.0f || d != 0.0f) : (a > 0.0f || c > 0.0f || d > 0.0f);
            if (on) seen[face[r]] = 1;
        }
    }
}
"""


def _seen_kernel(cupy):
    """Facets with a record at a selected point (fp32 records): fit_cops' predicate (Kn > 0 | KPhi > 0 | KV > 0) or,
    structural, the old late-window fine tuning's (!= 0); the predicate of positive_finetune's ft_seen."""
    if "seen" not in _SELECT:
        _SELECT["seen"] = cupy.RawKernel(_SEEN_SRC, "slopemul_seen")
    return _SELECT["seen"]


def _select_kernel(cupy):
    """out = (bin == k) ? L : 0.0 (an exact copy or +0.0, as np.where; no temporary arrays: the pool limit of the budget
    leaves only 64 MB of slack)."""
    if "kernel" not in _SELECT:
        _SELECT["kernel"] = cupy.ElementwiseKernel("int32 b, float64 L, int32 k", "float64 out",
                                                   "out = (b == k) ? L : 0.0", "slopemul_select_bin")
    return _SELECT["kernel"]


def _device_uuid(cupy, device_id):
    try:
        props = cupy.cuda.runtime.getDeviceProperties(int(device_id))
        raw = props.get("uuid")
        if isinstance(raw, (bytes, bytearray)):
            h = raw.hex()
            return "GPU-%s-%s-%s-%s-%s" % (h[0:8], h[8:12], h[12:16], h[16:20], h[20:32])
        return str(raw)
    except Exception as exc:                                   # recorded, never fatal
        return "unknown (%s)" % exc.__class__.__name__


def install_gpu(cfg):
    """Patch the evaluator classes of THIS process for the slope-multiplier model (idempotent for the same cfg)."""
    if any(k.startswith("COPS_POS_") and v not in ("", None) for k, v in os.environ.items()):
        raise RuntimeError("slope_mult: COPS_POS_* set (positive_activity.py patches the same methods); unset them")
    text = canonical(cfg)
    if _STATE["installed"]:
        if canonical(_STATE["config"]) != text:
            raise RuntimeError("slope_mult: already installed with another configuration")
        return False
    frozen = os.environ.get(ENV_FROZEN)
    if frozen is not None and frozen != text:
        raise RuntimeError("slope_mult: this process's configuration %s differs from the inherited frozen one %s"
                           % (text, frozen))
    _STATE["role"] = "worker" if frozen is not None else "driver"
    from cops_cuda import evaluator, multi_gpu, transport
    from cops_cuda.errors import BackendError, IdentityError, UnsupportedConfiguration
    from cops_cuda.memory import MemoryPlan
    _check_sources(Path(evaluator.__file__).resolve().parent)
    B, MG = evaluator.CudaBackend, multi_gpu.MultiGpuBackend
    orig = _STATE["originals"]
    for name in ("__init__", "activity_linear", "_upload_activity", "_upload_columns", "identity", "release"):
        orig["B." + name] = getattr(B, name)
    for name in ("activity_linear", "_evaluate", "identity"):
        orig["MG." + name] = getattr(MG, name)
    orig["pipe_worker_main"] = transport.pipe_worker_main
    init_signature = inspect.signature(B.__init__)
    nbins, nbumps = n_bins_bumps(cfg)
    bumps = cfg.get("basis") == "bumps"
    from cops_cuda import backend as bk
    WS = transport.WorkerState
    orig["B._evaluate"], orig["WS.handle"] = B._evaluate, WS.handle

    def __init__(self, *args, **kwargs):
        bound = init_signature.bind(self, *args, **kwargs)
        bound.apply_defaults()
        self._slopemul = None
        bound.arguments["scratch_per_face"] = int(bound.arguments["scratch_per_face"]) + BYTES_PER_FACE
        orig["B.__init__"](*bound.args, **bound.kwargs)

    def _phi_c(obj):
        """The contiguous bump block of the bound basis (cached per backend object and basis digest)."""
        key = obj.enhancement_basis_sha256
        cache = getattr(obj, "_slopemul_phi", None)
        if cache is None or cache[0] != key:
            phi = np.ascontiguousarray(np.asarray(obj.enhancement_basis, dtype=np.float64)[:, nbins:nbins + nbumps])
            if sha256_bytes(phi) != cfg["bumps_sha256"]:
                raise IdentityError("slope_mult: the bound bump columns are not the configured bump basis")
            obj._slopemul_phi = cache = (key, phi)
        return cache[1]

    def _L_of(obj, state):
        if not bumps:
            return _matvec_L(obj.Y, state.a)            # L: the original's first expression (Y @ a)
        w = np.ascontiguousarray(np.asarray(state.b, dtype=np.float64)[nbins:nbins + nbumps])
        return obj.enhancement_scale * _matvec_L(_phi_c(obj), w)

    def activity_linear(self, state):
        if state.b is None:
            return orig["B.activity_linear"](self, state)
        t0 = time.perf_counter()
        L = _L_of(self, state)
        _time_L(time.perf_counter() - t0)
        return L

    def multi_activity_linear(self, state):
        if state.b is None:
            return orig["MG.activity_linear"](self, state)
        t0 = time.perf_counter()
        L = _L_of(self, state)
        _time_L(time.perf_counter() - t0)
        return L

    def _strip(state, request):
        """Bump mode: L does not depend on a (held at 0, refused otherwise): the a columns are exact zeros and are not
        computed on the GPU. Returns (request for the GPU, the stripped harmonic indices)."""
        if not bumps or state.b is None:
            return request, ()
        # the bump model's L does not use a (the fitter holds a00..a24 at 0; its default-start evaluations inside
        # grain_crossing pass the prior a, which the model ignores)
        if not request.basis_index:
            return request, ()
        return bk.ColumnRequest((), request.enhancement_index, temperature=request.temperature,
                                kappa=request.kappa), tuple(request.basis_index)

    def _zeros(result, stripped, request):
        if not stripped:
            return result
        n = len(result.pNG)
        cols = dict(result.columns)
        for j in stripped:
            cols[("a", j)] = (np.zeros(n), np.zeros(n))
        return bk.EvaluationResult(result.pNG, result.pRG, cols, request)

    def multi_evaluate(self, state, request):
        if state.b is not None and self.activity != "host":
            raise UnsupportedConfiguration("slope_mult: the slope-multiplier activity is a host field; use "
                                           "--activity host (got %r)" % self.activity)
        req, stripped = _strip(state, request)
        return _zeros(orig["MG._evaluate"](self, state, req), stripped, request)

    def cuda_evaluate(self, state, request):
        req, stripped = _strip(state, request)
        return _zeros(orig["B._evaluate"](self, state, req), stripped, request)

    def _free(self):
        sm = getattr(self, "_slopemul", None)
        if sm is None:
            return
        self._slopemul = None
        for k in ("YT", "L", "Lg", "mrow", "s", "F", "bins", "ab", "g"):
            sm.pop(k, None)
        res = sm.pop("_reservation", None)
        if res is not None:
            res.release()

    def _ensure(self):
        sm = getattr(self, "_slopemul", None)
        if sm is not None and sm["basis_sha256"] == self.enhancement_basis_sha256:
            return sm
        _free(self)
        cupy, ds = self.cupy, self.dataset
        if self.enhancement_basis is None:
            raise BackendError("slope_mult: no template basis bound (the bin indicators are required)")
        basis = np.asarray(self.enhancement_basis)
        if basis.shape[1] != nbins + nbumps:
            raise IdentityError("slope_mult: bound basis has %d columns, the model %d" % (basis.shape[1], nbins + nbumps))
        bins = bins_from_basis(basis[:, :nbins], nbins)
        if len(bins) != int(cfg["nfaces"]) or sha256_bytes(bins) != cfg["bins_sha256"]:
            raise IdentityError("slope_mult: the bound basis is not the configured bin assignment (sha256 %s vs %s)"
                                % (sha256_bytes(bins), cfg["bins_sha256"]))
        ah2o = np.ascontiguousarray(np.asarray(ds.inputs.elements["AH2O"], dtype=np.float64))
        if sha256_bytes(ah2o) != cfg["ah2o_sha256"]:
            raise IdentityError("slope_mult: AH2O of the dataset's elements.bin differs from the configured A_orig")
        if not np.all(np.isfinite(ah2o)) or np.any(ah2o < 0.0):
            raise IdentityError("slope_mult: A_orig (AH2O) has negative or non-finite values")
        floor = float(cfg["floor_ratio"]) * ah2o
        nf = int(ds.nfaces_device)
        bins_loc = np.ascontiguousarray(ds.gather_faces(bins).astype(np.int32))
        floor_loc = np.ascontiguousarray(ds.gather_faces(floor))
        if bumps:     # columns m g (A_scale phi_k): the scaled bump block of the device facets, transposed
            yt = np.ascontiguousarray((self.enhancement_scale * ds.gather_faces(basis[:, nbins:nbins + nbumps])).T)
        else:
            yt = np.ascontiguousarray(np.asarray(ds.gather_faces(self.Y), dtype=np.float64).T)
        ny = yt.shape[0]
        g_glob, _ = held_factor(cfg)
        g_loc = np.ascontiguousarray(ds.gather_faces(g_glob)) if g_glob is not None else None
        plan = MemoryPlan("slope_mult: per-facet arrays of the slope-multiplier activity")
        plan.add_array("slopemul.YT", "static", (ny, nf), np.float64,
                       "A_scale phi_k^T (bumps) or Y^T (harmonics) of the device facets")
        for name in ("L", "Lg", "mrow", "s", "F", "g"):
            plan.add_array("slopemul." + name, "scratch", (nf,), np.float64)
        plan.add_array("slopemul.bins", "static", (nf,), np.int32)
        plan.add_array("slopemul.active_bool", "scratch", (nf,), np.bool_)
        with cupy.cuda.Device(self.device_id):
            reservation = self.budget.reserve(plan)
            try:
                with self.budget.allocation_guard(plan):
                    d = {"YT": cupy.asarray(yt), "L": cupy.empty(nf, dtype=np.float64),
                         "mrow": cupy.empty(nf, dtype=np.float64), "s": cupy.empty(nf, dtype=np.float64),
                         "F": cupy.asarray(floor_loc), "bins": cupy.asarray(bins_loc),
                         "ab": cupy.empty(nf, dtype=np.bool_)}
                    d["g"] = cupy.asarray(g_loc) if g_loc is not None else None
                    d["Lg"] = cupy.empty(nf, dtype=np.float64) if g_loc is not None else None
            except BaseException:
                reservation.release()
                raise
        local_counts = np.bincount(bins_loc, minlength=nbins)
        sm = dict(d, _reservation=reservation, basis_sha256=self.enhancement_basis_sha256, nbins=nbins,
                  nbumps=nbumps, bins_host=bins_loc, floor_host=floor_loc, g_host=g_loc, yt_host=yt if bumps else None,
                  armed=False, evaluations=0, checks=0, check_max_diff=0.0, check_flips=0, colbasis_checked=False,
                  colbasis_check=None, ny=ny, override=None,
                  info={"device_faces": nf, "device_id": int(self.device_id), "host": socket.gethostname(),
                        "pid": os.getpid(), "device_uuid": _device_uuid(cupy, self.device_id),
                        "device_bin_counts": [int(v) for v in local_counts],
                        "device_faces_at_zero_floor": int(np.sum(floor_loc == 0.0)),
                        "bins_sha256": cfg["bins_sha256"], "ah2o_sha256": cfg["ah2o_sha256"],
                        "floor_ratio": cfg["floor_ratio"], "bytes_reserved": int(plan.total_bytes),
                        "process_blas": blas_info(), "basis": cfg.get("basis", "harmonics"),
                        "held": (cfg.get("held") or {}).get("sha256"),
                        "held_nonzero_device_facets": int(np.sum(g_loc != 1.0)) if g_loc is not None else 0})
        self._slopemul = sm
        _STATE["backends"] += 1
        return sm

    def _g_now(sm):
        """The per-facet factor in force (held exp(h), with a fine-tuning override applied), device and host."""
        if sm["override"] is not None:
            return sm["override"]["g"], sm["override"]["g_host"]
        return sm["g"], sm["g_host"]

    def _upload_activity(self, state):
        sm = getattr(self, "_slopemul", None)
        if state.b is None:
            if sm is not None:
                sm["armed"] = False
            return orig["B._upload_activity"](self, state)
        if self.activity_mode == "device":
            raise UnsupportedConfiguration("slope_mult: the slope-multiplier activity needs the host activity "
                                           "(activity='device' forms Y a + A_scale B b on the device)")
        sm = _ensure(self)
        bvec = np.asarray(state.b, dtype=np.float64).reshape(-1)
        if len(bvec) != nbins + nbumps or not np.all(np.isfinite(bvec)):
            raise BackendError("slope_mult: %d template amplitudes (finite required) for %d bins + %d bumps"
                               % (len(bvec), nbins, nbumps))
        m = np.ascontiguousarray(bvec[:nbins])
        L = self.activity_linear(state) if self.activity_source is None else self.activity_source(state)
        L = np.ascontiguousarray(self.dataset.gather_faces(L), dtype=np.float64)
        cupy = self.cupy
        sm["L"].set(L)
        m_d = cupy.asarray(m)
        cupy.take(m_d, sm["bins"], out=sm["mrow"])
        cupy.multiply(sm["mrow"], sm["L"], out=sm["s"])
        g_d, g_h = _g_now(sm)
        if g_d is not None:
            cupy.multiply(sm["s"], g_d, out=sm["s"])              # s = (m L) g
            if sm["Lg"] is None:
                sm["Lg"] = cupy.empty_like(sm["L"])
            cupy.multiply(sm["L"], g_d, out=sm["Lg"])
        cupy.maximum(sm["F"], sm["s"], out=self.arrays["af_face"])
        cupy.greater(sm["s"], sm["F"], out=sm["ab"])
        self.arrays["active"][...] = sm["ab"]
        sm["armed"] = True
        n = sm["evaluations"]
        sm["evaluations"] = n + 1
        if n % CHECK_EVERY == 0:                     # the device field vs the host expressions, bitwise
            s_h, A_h, act_h = activity(L, m, sm["bins_host"], sm["floor_host"], g_h)
            af_d = self.arrays["af_face"].get()
            act_d = self.arrays["active"].get()
            diff = float(np.max(np.abs(af_d - A_h))) if len(A_h) else 0.0
            flips = int(np.sum(act_d.astype(bool) != act_h))
            sm["checks"] += 1
            sm["check_max_diff"] = max(sm["check_max_diff"], diff)
            sm["check_flips"] += flips
            if diff != 0.0 or flips:
                raise IdentityError("slope_mult: device activity differs from the host expression (max |diff| %r, "
                                    "%d active flips)" % (diff, flips))

    def _upload_columns(self, request):
        sm = getattr(self, "_slopemul", None)
        if sm is None or not sm["armed"]:
            return orig["B._upload_columns"](self, request)
        cupy = self.cupy
        nf = int(self.dataset.nfaces_device)
        cb = self.arrays["colbasis"]
        g_d, g_h = _g_now(sm)
        Lsel = sm["Lg"] if g_d is not None else sm["L"]
        k = 0
        if bumps and request.basis_index:
            raise BackendError("slope_mult: bump model: harmonic columns requested (they are exact zeros)")
        for j in request.basis_index:
            col = cb[k * nf:(k + 1) * nf]
            cupy.multiply(sm["YT"][int(j)], sm["mrow"], out=col)
            if g_d is not None:
                cupy.multiply(col, g_d, out=col)
            k += 1
        for i in request.enhancement_index:
            i = int(i)
            col = cb[k * nf:(k + 1) * nf]
            if 0 <= i < nbins:
                _select_kernel(cupy)(sm["bins"], Lsel, np.int32(i), col)        # dA/dm_i = [bin = i] L g
            elif bumps and nbins <= i < nbins + nbumps:
                cupy.multiply(sm["YT"][i - nbins], sm["mrow"], out=col)        # dA/dw = (A_scale phi) m g
                if g_d is not None:
                    cupy.multiply(col, g_d, out=col)
            else:
                raise BackendError("slope_mult: template column %r of %d" % (i, nbins + nbumps))
            k += 1
        self._columns_cache = None                   # the original's cache key no longer describes colbasis
        if not sm["colbasis_checked"] and k:         # once: the device columns vs the host expressions, bitwise
            sm["colbasis_checked"] = True
            host = cb[:k * nf].get().reshape(k, nf)
            L = sm["L"].get()
            mrow = np.asarray(cupy.asnumpy(sm["mrow"]))
            Lh = L * g_h if g_h is not None else L
            worst, kk = 0.0, 0
            if request.basis_index:
                Y = self.dataset.gather_faces(self.Y)
            for j in request.basis_index:
                want = np.asarray(Y[:, int(j)], np.float64) * mrow
                if g_h is not None:
                    want = want * g_h
                worst = max(worst, float(np.max(np.abs(host[kk] - want))))
                kk += 1
            for i in request.enhancement_index:
                i = int(i)
                if i < nbins:
                    want = np.where(sm["bins_host"] == i, Lh, 0.0)
                else:
                    want = sm["yt_host"][i - nbins] * mrow
                    if g_h is not None:
                        want = want * g_h
                worst = max(worst, float(np.max(np.abs(host[kk] - want))))
                kk += 1
            sm["colbasis_check"] = {"columns": k, "max_abs_diff": worst}
            if worst != 0.0:
                raise IdentityError("slope_mult: device column basis differs from the host expressions (%r)" % worst)

    # ---- fine-tuning support (slope_finetune.py): a per-facet factor override and dp/dA sensitivities
    def set_override(self, facets, g_values):
        """Replace the factor g (exp(h); 1 without a held term) at the given GLOBAL facets for the following
        evaluations (facets None: clear). Facets not on this device are ignored."""
        sm = _ensure(self)
        if facets is None:
            sm["override"] = None
            return 0
        ds, cupy = self.dataset, self.cupy
        facets = np.asarray(facets, np.int64)
        vals = np.asarray(g_values, np.float64)
        base_h = sm["g_host"] if sm["g_host"] is not None else np.ones(int(ds.nfaces_device))
        gh = np.array(base_h, dtype=np.float64, copy=True)
        if ds.local_faces is None:
            pos, ok = facets, np.ones(len(facets), bool)
        else:
            pos = np.searchsorted(ds.local_faces, facets)
            pos = np.minimum(pos, len(ds.local_faces) - 1)
            ok = ds.local_faces[pos] == facets
        gh[pos[ok]] = vals[ok]
        with cupy.cuda.Device(self.device_id):
            sm["override"] = {"g": cupy.asarray(gh), "g_host": gh, "facets": int(ok.sum())}
        return int(ok.sum())

    def sensitivities(self, state, facets):
        """dpNG/dA_e and dpRG/dA_e (rows of this backend x facets; temperature law of `state`; activity-independent:
        pressures are linear in A), for the given GLOBAL facets, with the evaluator's own column code (active = 1 at
        those facets, colbasis = their indicators, at most 16 per pass)."""
        sm = _ensure(self)
        from cops_cuda import evaluator as ev
        ds, cupy, c = self.dataset, self.cupy, self.constants
        facets = np.asarray(facets, np.int64)
        nf = int(ds.nfaces_device)
        if ds.local_faces is None:
            pos, ok = facets, np.ones(len(facets), bool)
        else:
            pos = np.minimum(np.searchsorted(ds.local_faces, facets), len(ds.local_faces) - 1)
            ok = ds.local_faces[pos] == facets
        npt = len(self.rows)
        ng, rg = np.zeros((npt, len(facets))), np.zeros((npt, len(facets)))
        from cops_cuda.memory import scratch_reservation
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                self.budget, "slope_mult sensitivities", {"active_saved": (nf, self.arrays["active"].dtype)}):
            self._upload_activity(state)
            act = self.arrays["active"]
            saved = act.copy()
            cb = self.arrays["colbasis"]
            for c0 in range(0, len(facets), ev.FUSED_MAX_COLUMNS):
                cols = list(range(c0, min(c0 + ev.FUSED_MAX_COLUMNS, len(facets))))
                cb[:len(cols) * nf] = 0.0
                for k, cidx in enumerate(cols):
                    if ok[cidx]:
                        cb[k * nf + int(pos[cidx])] = 1.0
                        act[int(pos[cidx])] = 1
                codes = list(ev.PRESSURE_CODES) + [1000 + k for k in range(len(cols))] + \
                    [2000 + k for k in range(len(cols))]
                S = dict(zip(codes, self._run_sums(codes, self._params(state))))
                for k, cidx in enumerate(cols):
                    ng[:, cidx] = c.cNG["H2O"] * c.bm["H2O"] * S[1000 + k]
                    rg[:, cidx] = c.cRG["H2O"] * S[2000 + k]
                act[...] = saved
            del saved
            self._columns_cache = None
        return ng, rg

    def seen_rows(self, row_mask, structural=False):
        """Global facet mask: facets with a record (predicate above) at the rows of this backend selected by row_mask
        (local row order)."""
        ds, cupy = self.dataset, self.cupy
        if ds.streamed:
            raise BackendError("slope_mult: seen_rows needs resident records")
        if ds.record_precision != "fp32":
            raise BackendError("slope_mult: seen_rows is written for fp32 records (--precision mixed)")
        mask = np.asarray(row_mask, bool)
        if mask.shape != (len(self.rows),):
            raise BackendError("slope_mult: row mask of %s for %d rows" % (mask.shape, len(self.rows)))
        nblocks = int(self._first_block[-1])
        from cops_cuda.memory import scratch_reservation
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                self.budget, "slope_mult seen rows", {"seen": (int(ds.nfaces_device), np.uint8),
                                                      "mask": (len(mask), np.uint8)}):
            d_mask = cupy.asarray(mask.astype(np.uint8))
            d_seen = cupy.zeros(int(ds.nfaces_device), dtype=cupy.uint8)
            if nblocks:
                d = ds.arrays
                _seen_kernel(cupy)((min(nblocks, 1 << 20),), (256,),
                                   (self.arrays["blocks"], np.int64(nblocks), d["kernel.face"], d["kernel.kn"],
                                    d["kernel.kphi"], d["kernel.kv"], d_mask, d_seen, np.int32(1 if structural else 0)))
            host = d_seen.get().astype(bool)
            del d_mask, d_seen
        return ds.scatter_faces(host, False)

    def multi_seen_rows(self, row_mask, structural=False):
        row_mask = np.asarray(row_mask, bool)
        out = None
        for w in self.workers:
            self._send(w, ("slopemul_seen", np.ascontiguousarray(row_mask[w["positions"]]), bool(structural)),
                       "slopemul_seen")
        for w in self.workers:
            m = self._receive(w, "slopemul_seen")[1]
            out = m if out is None else (out | m)
        return out

    def multi_set_override(self, facets, g_values):
        for w in self.workers:
            self._send(w, ("slopemul_override", facets, g_values), "slopemul_override")
        return [self._receive(w, "slopemul_override")[1] for w in self.workers]

    def multi_sensitivities(self, state, facets):
        n = len(self.rows)
        ng, rg = np.zeros((n, len(facets))), np.zeros((n, len(facets)))
        for w in self.workers:
            self._send(w, ("slopemul_sens", state, np.asarray(facets, np.int64)), "slopemul_sens")
        for w in self.workers:
            _, a_ng, a_rg = self._receive(w, "slopemul_sens")
            ng[w["positions"]], rg[w["positions"]] = a_ng, a_rg
        return ng, rg

    def ws_handle(self, message):
        kind = message[0]
        # the worker requests of the fine tuning (never sent by a fit); each calls the method under the name this
        # patch installs on CudaBackend (WORKER_REQUESTS; 2026-11-10 fix: "slopemul_seen" called the undefined
        # seen_rows, job 123985)
        if kind == "slopemul_override":
            try:
                return ("ok", self.backend.slopemul_set_override(message[1], message[2]))
            except Exception:
                import traceback as _tb
                return ("error", _tb.format_exc())
        if kind == "slopemul_seen":
            try:
                return ("ok", self.backend.slopemul_seen_rows(message[1], message[2]))
            except Exception:
                import traceback as _tb
                return ("error", _tb.format_exc())
        if kind == "slopemul_sens":
            try:
                ng, rg = self.backend.slopemul_sensitivities(message[1], message[2])
                return ("ok", ng, rg)
            except Exception:
                import traceback as _tb
                return ("error", _tb.format_exc())
        return orig["WS.handle"](self, message)

    def identity(self):
        value = orig["B.identity"](self)
        sm = getattr(self, "_slopemul", None)
        info = None
        if sm is not None:
            info = dict(sm["info"], evaluations=sm["evaluations"], host_checks=sm["checks"],
                        host_check_max_abs_diff=sm["check_max_diff"], host_check_active_flips=sm["check_flips"],
                        colbasis_check=sm["colbasis_check"])
        value["model_extension"] = {"slope_multiplier": {"config": installed_config(), "backend": info,
                                                         "module": str(Path(__file__).resolve()),
                                                         "module_sha256": module_sha256(), "role": _STATE["role"]}}
        return value

    def release(self):
        try:
            orig["B.release"](self)
        finally:
            _free(self)

    def multi_identity(self):
        value = orig["MG.identity"](self)
        value["model_extension"] = {"slope_multiplier": {"config": installed_config(), "module_sha256": module_sha256(),
                                                         "driver_process_blas": blas_info(),
                                                         "activity_linear_timing": activity_linear_timing()}}
        return value

    B.__init__ = __init__
    B.activity_linear, B._upload_activity, B._upload_columns = activity_linear, _upload_activity, _upload_columns
    B.identity, B.release, B.close = identity, release, release      # close is an alias of release in the class
    B.slopemul_set_override, B.slopemul_sensitivities = set_override, sensitivities
    B.set_override, B.sensitivities = set_override, sensitivities
    B.slopemul_seen_rows = seen_rows
    MG.slopemul_seen_rows = multi_seen_rows
    MG.activity_linear, MG._evaluate, MG.identity = multi_activity_linear, multi_evaluate, multi_identity
    MG.slopemul_set_override, MG.slopemul_sensitivities = multi_set_override, multi_sensitivities
    if bumps:
        B._evaluate = cuda_evaluate
    WS.handle = ws_handle
    transport.pipe_worker_main = _pipe_worker_main
    os.environ[ENV_FROZEN] = text
    _STATE["config"] = dict(cfg)
    _STATE["installed"] = True
    say("GPU patch installed (pid %d, %s, role %s): %s" % (os.getpid(), socket.gethostname(), _STATE["role"], text))
    return True


def _pipe_worker_main(conn, bind_kwargs, device_id, shm_name, nfaces):
    """Local spawned worker: install the patch from the parent's frozen configuration, then the original loop."""
    frozen = os.environ.get(ENV_FROZEN)
    if frozen is None:
        raise RuntimeError("slope_mult: worker without %s" % ENV_FROZEN)
    install_gpu(json.loads(frozen))
    return _STATE["originals"]["pipe_worker_main"](conn, bind_kwargs, device_id, shm_name, nfaces)


def mpi_config_check(cfg, comm):
    """Every MPI rank: the same canonical configuration and module (allgather; collective, call on all ranks at the same
    point, before gpu_fit.main). Returns the per-rank records."""
    import socket as _s
    mine = {"config": canonical(cfg), "module_sha256": module_sha256(), "host": _s.gethostname(), "pid": os.getpid(),
            "rank": int(comm.Get_rank())}
    allr = comm.allgather(mine)
    bad = [r["rank"] for r in allr if r["config"] != mine["config"] or r["module_sha256"] != mine["module_sha256"]]
    if bad:
        raise RuntimeError("slope_mult: ranks %s have another configuration or module than rank %d" % (bad, mine["rank"]))
    return [{k: r[k] for k in ("rank", "host", "pid")} for r in allr]
