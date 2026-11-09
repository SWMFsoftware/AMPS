# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_activity.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_activity.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 90d57cab3c29df18fcfcb83b469fbcb79327a82954866a95b70e4e712297c8d8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Positive H2O activity on the GPU (user decisions 2026-11-07): installs patched kernels (kernels/pressure.cu,
#   kernels/activity.cu) into the GpuModel workers. Env: COPS_POS_ENABLE=1; COPS_POS_NIGHT_RATIO (night/shadow ratio,
#   H2O and CO2); COPS_POS_ACTIVITY=clip|exp (exp: A = exp(Y a + c T), Jacobian A * basis); COPS_POS_FLOOR_RATIO (H2O flux
#   floor F = max(q A, ratio q_max A_orig); CO2 none); COPS_POS_FROZEN_CONFIG. Refuses FP64 records and restart.cu (no
#   fine tuning) when active. With everything off the kernels are bitwise identical to the unchanged ones.
# USED BY: positive_problem.py / positive_fit.py (jobs/positive_fit.sbatch), positive_activity_check.py.
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
"""Positive H2O activity, night/shadow production ratio, per-record H2O and (optional) CO2 flux floors in the GPU evaluator.

User decisions 2026-11-07 (CO2 floor: 2026-11-08). Local, untracked model extension of the unchanged CUDA evaluator
(``srcCG/LiouvilleCOPS_cuda/cops_cuda``): nothing in the git tree is edited. ``install()`` patches the
evaluator classes and the CUDA source registry of the CURRENT process at run time; it is configured ONLY
through environment variables, so the spawned multi-GPU workers (which inherit the environment) configure
themselves identically (they run ``install()`` before they bind, see "Processes" below).

Environment interface
---------------------
COPS_POS_ENABLE=1          required for any effect; unset/0: ``install()`` does nothing (and refuses if any other
                           COPS_POS_* variable is set, so a forgotten ENABLE never silently runs the original model).
COPS_POS_NIGHT_RATIO=r     night/shadow production ratio, 0 <= r <= 1: qmin = r * 600 / rhA**p for H2O (p = 4.8, or
                           4.2143229 for HeliocentricModel shap5) and qmin = r * 600 / rhA**2 for CO2 (both species;
                           user: 0.02 "as in srcCG"). UNSET = unchanged per species (fit_cops Model.load_kernels:
                           0.2 for H2O, 0.1 for CO2). The expression form is the reference's, so r = 0.2 reproduces
                           the H2O qmin bitwise.
COPS_POS_ACTIVITY=clip|exp clip (default): af = max(act_lin, 0), d af/dc = [act_lin > 0] * basis (unchanged).
                           exp: af = exp(act_lin) for every facet (positive by construction, no clip, no reference
                           constant), d af/dc = af * basis. act_lin = Y a + A_scale * (B b) is unchanged (templates
                           enter the exponent). Host activity (CudaBackend._upload_activity) and device activity
                           (kernels/activity.cu, activity="device") are both handled.
COPS_POS_FLOOR_RATIO=f     f >= 0 (default 0 = no floor): per kernel record r of point p on facet j the H2O flux is
                           F = max(q_H2O(r) * af(j), F_min(r)), F_min(r) = qmax_H(p) * (f * A_orig(j)), A_orig = the
                           base's original H2O activity map (elements.bin column 32 "AH2O"; the run map). F replaces
                           q_H2O * af in every pressure sum (NG: Kn F u; RG: KPhi F and VdotB KV F u) and in the
                           temperature-derivative weights, with the record's own temperature law (same u). A record
                           at the floor (q_H2O * af < F_min) contributes 0 to the activity/template columns. Without
                           COPS_POS_FLOOR_CO2 CO2 has no floor ("floor form (a)": constant minimum emission at all
                           times).
COPS_POS_FLOOR_CO2=1       (2026-11-08, user: "use the same floor for CO2"; unset/0 = off, the behaviour before) the
                           same floor form for CO2 with the same ratio f = COPS_POS_FLOOR_RATIO (required > 0): per
                           kernel record r of point p on facet j the CO2 flux is
                               F_CO2 = max(kappa * qc(r), F_min,CO2(r)),  qc = Q_CO2(r) * A_CO2,run(j),
                               F_min,CO2 = qmax_C(p) * (f * A_CO2,run(j)),
                           Q_CO2 = the record's CO2 production per unit activity (night/shadow ratio as above),
                           qmax_C = 600 / rhA**2, A_CO2,run = the run's CO2 map (elements.bin column 33 "ACO2" = the
                           CO2 activity at kappa = 1, the analog of A_orig for H2O; checked >= 0 and finite). So CO2
                           emission never drops below f times the run's CO2 emission at maximum production, even for
                           kappa -> 0. A record is at the CO2 floor when kappa * qc < F_min,CO2. Pressures:
                           p = p_H2O + kappa * p_CO2(records not at the CO2 floor) + p_CO2,floor (F_min,CO2 in place
                           of kappa * qc, with the record's own temperature law u, in NG Kn, RG KPhi and VdotB KV);
                           the kappa column is d p / d kappa = p_CO2(records not at the CO2 floor) (records at the
                           floor contribute 0); the temperature columns get the floor records' terms without kappa.
                           Kernels eval_fused_m_p5f / eval_fused_m_p14f of the pressure.cu copy (header item 5);
                           the activity/template columns (H2O only) are unchanged.
COPS_POS_FROZEN_CONFIG     written by ``install()`` (canonical JSON of the configuration); a worker whose own
                           environment gives a different configuration refuses to start. "co2_floor" appears in it
                           only when the CO2 floor is on, so the canonical form of every earlier configuration is
                           unchanged.

Python API: ``install()`` (idempotent; True when patched), ``config()`` (dict of the active settings, kernel copies and
their SHA-256, for manifests), ``kernel_sha256()`` ({module: sha256 of the copy in use}), ``modifications_active()``.

What is patched (only after install() with COPS_POS_ENABLE=1)
------------------------------------------------------------
* ``cuda_sources.source_path``: modules "pressure" and "activity" are read from ``model_ext/kernels/`` (copies of the
  tracked sources with the changes documented in their headers; same kernel names and signatures). Every process
  that compiles them prints the path and SHA-256 to stderr and lists them in ``CudaBackend.identity()``.
* ``cuda_sources.module``: refuses the module "restart" (kernels/restart.cu: fine tuning, restart_tool facet shares /
  point changes) while a modification is active: those kernels implement the original clipped, floor-free model.
* ``evaluator.HostConstants``: qmin from COPS_POS_NIGHT_RATIO (when set).
* ``evaluator.CudaBackend``:
    - ``__init__``: refuses FP64 records (``record_precision`` != "fp32") when a modification is active (the
      floor and the exp derivative exist only in the mixed path, eval_fused_m_*); reserves 8 B/facet more
      (``scratch_per_face``) for the floor block when the resident records are sized (streaming).
    - ``_check_dataset`` (once per dataset, before the first bind): the device array ``mesh.ACO2`` (per device
      facet, face compaction respected) becomes [ACO2 | f * AH2O] (2 nfaces values, its own memory reservation;
      the first block is the uploaded ACO2, checked against the loader's digest). Kernels index the second block
      with ``nfaces``; every other reader of ``aco2[f]`` is unaffected. With COPS_POS_FLOOR_CO2=1 two values
      follow: [ACO2 | f * AH2O | kappa, f] (2 nfaces + 2); kappa (NaN at upload) is written before every
      evaluation, so ``DeviceDataset.verify_roundtrip`` (tools/tests only) would report that array as changed.
    - ``_evaluate`` (only with COPS_POS_FLOOR_CO2=1): ``CudaBackend._evaluate`` with the CO2 floor: writes kappa,
      runs eval_fused_m_p5f / p14f in place of eval_fused_m_p5 / p14 (the column kernels unchanged; with streamed
      records and too many sums for the planned partial buffer the head kernel and the columns run as two passes,
      the same per-code results) and assembles the pressures and columns as stated above.
    - ``_upload_activity`` / ``_device_activity``: exp activity (af = exp(act_lin), active = 2).
    - ``identity``: adds ``model_extension.positive_activity`` (this config and the compiled modules).
    - ``release`` (and its alias ``close``): returns the floor block's reservation.
* ``multi_gpu.MultiGpuBackend.identity``: adds the parent's config.
* ``transport.pipe_worker_main`` (local transport): the spawned worker's target becomes ``_pipe_worker_main``, which
  runs ``install()`` in the worker (from the inherited environment, checked against COPS_POS_FROZEN_CONFIG) before
  the original worker loop binds its backend.

With COPS_POS_ENABLE=1 and every other setting at its default, the patched code is bitwise the original evaluator
(pressures and every Jacobian column, 1 and 2 GPUs; positive_activity_check.py). Without COPS_POS_FLOOR_CO2 every
configuration evaluates exactly as before the CO2 floor was added (same host code, same kernels; the 41 kernels
that existed compile to the same PTX; positive_activity_check.py compares 0.02 and 0.002 with the 2026-11-07 module
and kernels bitwise).

Processes
---------
* One process / one GPU: call ``install()`` before the first backend binds (before ``gpu_fit.main()``).
* Local multi-GPU (``transport="local"``, the default): the workers are patched by ``_pipe_worker_main``. The module
  must be importable as ``positive_activity`` in the workers: keep ``model_ext`` on ``sys.path`` of the parent
  (spawn passes ``sys.path`` to the children).
* MPI transport: every rank must run ``install()`` itself before ``mpi_worker_main`` / the fit (the driver does it at
  start-up on every rank); rank 0's own GPU is in process.

Not patched (the driver's responsibility, documented here)
---------------------------------------------------------
* ``fit_cops`` (CPU) is unchanged: the CPU Model.evaluate, ``--cpu-comparison``, and the start-up consistency check of
  ``fit_cops.main`` (MF.evaluate vs the LiouvilleCOPS pressures of COPS_model.dat) will necessarily FAIL when a
  modification is active (they describe the original model). That check is handled by the driver; this module does
  not work around it.
* Host-side diagnostics of fit_cops / gpu_fit that apply the clip to act_lin (activity_facets.dat, activity
  diagnostics, the Tikhonov deviation D, ``fit_restart.npz`` ``activity_linear`` interpreted as clipped by the restart
  tools): in exp mode act_lin is the exponent; such products describe the clip model.
* ``Problem.enhancement_scale`` (A_scale ~ 2.2e18 on 12M) multiplies the template columns inside act_lin; in exp mode
  that scale is in the exponent (a template amplitude of 1 would give exp(2e18)). The scale is passed by the caller.
"""

import hashlib
import inspect
import json
import math
import os
import socket
import time
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
KERNEL_DIRECTORY = HERE / "kernels"
PATCHED_MODULES = ("pressure", "activity")
REFUSED_MODULES = ("restart",)
# The tracked sources the copies and wrappers were made from (AMPS master 7ea4885; unchanged since 47fccd5).
# install() refuses to patch a tree in which any of them changed: the copies/wrappers would no longer match.
ORIGINAL_SHA256 = {
    "kernels/pressure.cu": "dc4b6f4ceeba96bcd2df335c6dff607d6e545c2c363b8bf897836652010798f2",
    "kernels/activity.cu": "61dc45410249f76d42c7b27fa773a08840f6ae7b759ae911c9335322de913048",
    "evaluator.py": "2f1e9d5fadbc4111e3e5222d32f8098e51e8fd01c1e85f5266730689ec53d52e",
    "transport.py": "605ed2546625b0e029d277fac0757635ff978081b24526def324d18cb0256304",
    "multi_gpu.py": "1c3a9c1a377784f680a1dcb7a8f58e0c58f925ced241529d8f4c89a07dd516fc",
    "cuda_sources.py": "16d67d1a536a153891bd28df43c09eccf0db1239c2ee4c212e0827f326312719",
    "device_store.py": "35ae928704e5717a3ceb01e7bad102ea9c1cd1961a71e47ff1fa0a4b07e2fbe3",
}
ENV_ENABLE = "COPS_POS_ENABLE"
ENV_NIGHT = "COPS_POS_NIGHT_RATIO"
ENV_ACTIVITY = "COPS_POS_ACTIVITY"
ENV_FLOOR = "COPS_POS_FLOOR_RATIO"
ENV_FLOOR_CO2 = "COPS_POS_FLOOR_CO2"
ENV_FROZEN = "COPS_POS_FROZEN_CONFIG"
ENV_VARIABLES = (ENV_ENABLE, ENV_NIGHT, ENV_ACTIVITY, ENV_FLOOR, ENV_FLOOR_CO2)
ACTIVITY_MODES = ("clip", "exp")
REFERENCE_QMIN_RATIO = {"H2O": 0.2, "CO2": 0.1}      # fit_cops Model.load_kernels
FLOOR_BYTES_PER_FACE = 8
EXP_ACTIVE = 2                                       # active[] marker of the exp activity (pressure.cu copy)
POLICY = ("user 2026-11-07: night/shadow production ratio qmin = r qmax for H2O and CO2 (r = 0.02 'as in srcCG'); "
          "positive activity af = exp(Y a + A_scale B b) with dA/dc = A basis; per-record H2O flux floor "
          "F = max(q_H2O A, f qmax_H A_orig), A_orig = the base's AH2O map, at all times ('floor form (a)'), "
          "records at the floor do not depend on the activity, CO2 without floor")
POLICY_CO2_FLOOR = ("; user 2026-11-08 'use the same floor for CO2' (COPS_POS_FLOOR_CO2=1, the lead's reading): "
                    "per-record CO2 flux F_CO2 = max(kappa Q_CO2 A_CO2,run, f qmax_C A_CO2,run), A_CO2,run = the "
                    "base's ACO2 map (the CO2 activity at kappa = 1), same ratio f as the H2O floor; records at the "
                    "CO2 floor do not depend on kappa")
# CO2-floor sum codes of the pressure.cu copy (header item 5), after evaluator.PHI_CO2_CODE = 14
CO2_FLOOR_CODES = (15, 16, 17, 18)       # PhiCO2 of the active records; NG, RG (VdotB), RG (Phi) sums at the floor
CO2_FLOOR_T_CODES = tuple(range(19, 25))  # NG (19..21) and RG (22..24) temperature sums of the records at the floor
CO2_FLOOR_KERNELS = {"eval_fused_p5": "eval_fused_m_p5f", "eval_fused_p14": "eval_fused_m_p14f"}
CO2_FLOOR_SLOTS = 2                      # mesh.ACO2[2 nfaces] = kappa, [2 nfaces + 1] = floor ratio

_STATE = {"installed": False, "config": None, "originals": {}, "compiled": [], "role": None}


# ----------------------------------------------------------------------------------------------- configuration
def _ratio(value, name, upper=None):
    if value is None or value.strip() == "":
        return None
    try:
        r = float(value)
    except ValueError:
        raise ValueError("%s=%r is not a number" % (name, value)) from None
    if not math.isfinite(r) or r < 0.0 or (upper is not None and r > upper):
        raise ValueError("%s=%r must be finite and in [0, %s]" % (name, value, "inf" if upper is None else upper))
    return r


def _read_environment(environ=None):
    env = os.environ if environ is None else environ
    raw = {k: env.get(k) for k in ENV_VARIABLES}
    enable = (raw[ENV_ENABLE] or "").strip()
    if enable not in ("", "0", "1"):
        raise ValueError("%s=%r must be 1 (enable) or 0/unset" % (ENV_ENABLE, raw[ENV_ENABLE]))
    if enable != "1":
        others = {k: v for k, v in raw.items() if k != ENV_ENABLE and v not in (None, "")}
        if others:
            raise ValueError("positive_activity: %s set but %s is not 1; set %s=1 or unset them"
                             % (others, ENV_ENABLE, ENV_ENABLE))
        return {"enabled": False, "modified": False, "environment": raw}
    night = _ratio(raw[ENV_NIGHT], ENV_NIGHT, upper=1.0)
    activity = (raw[ENV_ACTIVITY] or "clip").strip().lower()
    if activity not in ACTIVITY_MODES:
        raise ValueError("%s=%r must be one of %s" % (ENV_ACTIVITY, raw[ENV_ACTIVITY], ACTIVITY_MODES))
    floor = _ratio(raw[ENV_FLOOR], ENV_FLOOR)
    floor = 0.0 if floor is None else floor
    co2 = (raw[ENV_FLOOR_CO2] or "").strip()
    if co2 not in ("", "0", "1"):
        raise ValueError("%s=%r must be 1 (CO2 flux floor with the ratio %s) or 0/unset"
                         % (ENV_FLOOR_CO2, raw[ENV_FLOOR_CO2], ENV_FLOOR))
    co2_floor = co2 == "1"
    if co2_floor and not floor > 0.0:
        raise ValueError("%s=1 uses the ratio %s, which is %r: set %s > 0 (the CO2 floor has no ratio of its own)"
                         % (ENV_FLOOR_CO2, ENV_FLOOR, raw[ENV_FLOOR], ENV_FLOOR))
    qmin_ratio = dict(REFERENCE_QMIN_RATIO) if night is None else {"H2O": night, "CO2": night}
    modified = night is not None or activity != "clip" or floor != 0.0 or co2_floor
    return {"enabled": True, "modified": modified, "night_ratio": night, "qmin_ratio": qmin_ratio,
            "activity": activity, "floor_ratio": floor, "co2_floor": co2_floor,
            "co2_floor_ratio": floor if co2_floor else None, "environment": raw}


def _canonical(cfg):
    """Canonical JSON of the configuration. "co2_floor" is included only when it is on: the canonical form of every
    configuration without the CO2 floor is the one of the 2026-11-07 module (COPS_POS_FROZEN_CONFIG, manifests)."""
    keys = ("enabled", "modified", "night_ratio", "qmin_ratio", "activity", "floor_ratio")
    if cfg.get("co2_floor"):
        keys += ("co2_floor",)
    return json.dumps({k: cfg.get(k) for k in keys}, sort_keys=True)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kernel_sha256():
    """{module name: SHA-256 of the kernel copy in model_ext/kernels}."""
    return {name: _sha256(KERNEL_DIRECTORY / (name + ".cu")) for name in PATCHED_MODULES}


def modifications_active():
    cfg = _STATE["config"] if _STATE["installed"] else _read_environment()
    return bool(cfg.get("modified"))


def config():
    """The active settings for manifests: the configuration, the kernel copies, this module's identity."""
    cfg = dict(_STATE["config"] if _STATE["installed"] else _read_environment())
    cfg.update(installed=bool(_STATE["installed"]), role=_STATE["role"],
               policy=POLICY + (POLICY_CO2_FLOOR if cfg.get("co2_floor") else ""),
               module={"path": str(Path(__file__).resolve()), "sha256": _sha256(__file__)},
               kernels={name: {"path": str(KERNEL_DIRECTORY / (name + ".cu")), "sha256": sha,
                               "original_sha256": ORIGINAL_SHA256["kernels/%s.cu" % name]}
                        for name, sha in kernel_sha256().items()},
               definitions={
                   "night": "qmin = r * 600 / rhA**p (H2O), r * 600 / rhA**2 (CO2); unset: 0.2 (H2O), 0.1 (CO2)",
                   "activity": {"clip": "af = max(act_lin, 0), dA/dc = [act_lin > 0] basis",
                                "exp": "af = exp(act_lin), dA/dc = af basis"},
                   "floor": "F = max(q_H2O af, F_min), F_min = qmax_H(p) * (floor_ratio * AH2O(face)); "
                            "records at the floor: dF/dc = 0; CO2 no floor" + (
                                " unless co2_floor (below)" if cfg.get("co2_floor") else ""),
                   "co2_floor": ("F_CO2 = max(kappa * qc, F_min,C), qc = Q_CO2 * ACO2(face) (per unit kappa), "
                                 "F_min,C = qmax_C(p) * (floor_ratio * ACO2(face)), ACO2 = the base's CO2 map (= the "
                                 "CO2 activity at kappa = 1); at the floor when kappa * qc < F_min,C; p = p_H2O + "
                                 "kappa p_CO2(not at the floor) + p_CO2(F_min,C); d p / d kappa = p_CO2(not at the "
                                 "floor); kernels eval_fused_m_p5f/p14f" if cfg.get("co2_floor") else
                                 "off (COPS_POS_FLOOR_CO2 unset or 0)")},
               patched=["cuda_sources.source_path (pressure, activity)", "cuda_sources.module (refuses restart "
                        "when modified)", "evaluator.HostConstants", "CudaBackend.__init__/_check_dataset/"
                        "_upload_activity/_device_activity/identity/release", "MultiGpuBackend.identity",
                        "transport.pipe_worker_main"] + (["CudaBackend._evaluate (CO2 floor)"]
                                                         if cfg.get("co2_floor") else []),
               refused=(["FP64 records (CudaBackend record_precision != fp32)", "kernels/restart.cu "
                         "(fine tuning, restart_tool shares / point changes)"] if cfg.get("modified") else []),
               compiled=[dict(c) for c in _STATE["compiled"]])
    return cfg


# ----------------------------------------------------------------------------------------------------- install
def _log(message):
    print("positive_activity[pid %d %s %s]: %s" % (os.getpid(), socket.gethostname(), _STATE["role"] or "-",
                                                   message), file=sys.stderr, flush=True)


def _check_sources(cops_cuda_dir):
    bad = {}
    for rel, sha in ORIGINAL_SHA256.items():
        got = _sha256(cops_cuda_dir / rel)
        if got != sha:
            bad[rel] = got
    if bad:
        raise RuntimeError("positive_activity: the tracked cops_cuda sources differ from the version the kernel copies "
                           "and wrappers were made from (%s); re-make the copies before using this extension" % bad)
    for name in PATCHED_MODULES:
        text = (KERNEL_DIRECTORY / (name + ".cu")).read_text()
        if ORIGINAL_SHA256["kernels/%s.cu" % name] not in text:
            raise RuntimeError("positive_activity: %s/%s.cu does not name its original's SHA-256" % (KERNEL_DIRECTORY, name))


def install():
    """Patch the GPU evaluator of this process when COPS_POS_ENABLE=1 (idempotent). Returns True when patched."""
    if _STATE["installed"]:
        if _canonical(_read_environment()) != _canonical(_STATE["config"]):
            raise RuntimeError("positive_activity: COPS_POS_* changed after install(); the patched configuration is %s"
                               % _canonical(_STATE["config"]))
        return True
    cfg = _read_environment()
    if not cfg["enabled"]:
        return False
    frozen = os.environ.get(ENV_FROZEN)
    if frozen is not None and frozen != _canonical(cfg):
        raise RuntimeError("positive_activity: this process's COPS_POS_* give %s, but the parent installed %s"
                           % (_canonical(cfg), frozen))
    _STATE["role"] = "worker" if frozen is not None else "driver"

    from cops_cuda import cuda_sources, evaluator, multi_gpu, precision, transport
    from cops_cuda.errors import BackendError, IdentityError, UnsupportedConfiguration
    from cops_cuda.memory import MemoryPlan
    _check_sources(Path(cuda_sources.__file__).resolve().parent)
    compiled = [key for key in cuda_sources._modules if key[0] in PATCHED_MODULES + REFUSED_MODULES]
    if compiled:
        raise RuntimeError("positive_activity: install() must run before any evaluator kernel is compiled; "
                           "already compiled: %s" % compiled)
    orig = _STATE["originals"]

    # --- CUDA sources ------------------------------------------------------------------------------------------
    orig["source_path"] = cuda_sources.source_path

    def source_path(name):
        if name in PATCHED_MODULES:
            orig["source_path"](name)                      # the original's name check (and that it exists)
            path = KERNEL_DIRECTORY / (name + ".cu")
            if not path.is_file():
                raise FileNotFoundError("positive_activity kernel copy not found: %s" % path)
            return path
        return orig["source_path"](name)

    orig["module"] = cuda_sources.module
    default_options = inspect.signature(cuda_sources.module).parameters["options"].default

    def module(name, device_id=0, options=default_options):
        if name in REFUSED_MODULES and cfg["modified"]:
            raise UnsupportedConfiguration(
                "positive_activity: kernels/%s.cu implements the original clipped, floor-free model; it is refused "
                "while COPS_POS_* modify the model (%s)" % (name, _canonical(cfg)))
        key = (name, int(device_id), precision.check_nvrtc_options(options))
        new = key not in cuda_sources._modules
        entry = orig["module"](name, device_id, options)
        if new and name in PATCHED_MODULES:
            identity = dict(entry[1])
            identity.update(pid=os.getpid(), host=socket.gethostname())
            _STATE["compiled"].append(identity)
            _log("compiled %s on device %d from %s (sha256 %s)" % (name, int(device_id), identity["path"],
                                                                   identity["sha256"]))
        return entry

    # --- per-point constants (night ratio) -------------------------------------------------------------------
    BaseConstants = evaluator.HostConstants

    class PositiveHostConstants(BaseConstants):
        """HostConstants with qmin = r * 600 / rhA**p (both species) when COPS_POS_NIGHT_RATIO is set."""

        def __init__(self, model):
            super().__init__(model)
            if cfg["night_ratio"] is not None:
                from cops_cuda import reference
                fit_cops = reference.load("fit_cops")
                rhA = model.rh / fit_cops.AU
                p = 4.2143229 if model.helio == "shap5" else 4.8     # Model.load_kernels
                r = cfg["night_ratio"]
                self.qmin = {"H2O": r * 600 / rhA ** p, "CO2": r * 600 / rhA ** 2}

    # --- the backend ----------------------------------------------------------------------------------------
    B = evaluator.CudaBackend
    orig["init"], orig["check_dataset"] = B.__init__, B._check_dataset
    orig["upload_activity"], orig["device_activity"] = B._upload_activity, B._device_activity
    orig["identity"], orig["release"], orig["close"] = B.identity, B.release, B.close
    orig["multi_identity"] = multi_gpu.MultiGpuBackend.identity
    init_signature = inspect.signature(B.__init__)

    def __init__(self, *args, **kwargs):
        bound = init_signature.bind(self, *args, **kwargs)
        bound.apply_defaults()
        rp = bound.arguments["record_precision"]
        rp = precision.record_precision_default() if rp is None else rp
        if cfg["modified"] and rp != "fp32":
            raise UnsupportedConfiguration(
                "positive_activity: the modified model (%s) exists only in the mixed-precision kernels "
                "(eval_fused_m_*); FP64 records are refused (use --precision mixed / "
                "COPS_CUDA_RECORD_PRECISION=fp32)" % _canonical(cfg))
        self._positive_activity = None
        bound.arguments["scratch_per_face"] = int(bound.arguments["scratch_per_face"]) + FLOOR_BYTES_PER_FACE
        orig["init"](*bound.args, **bound.kwargs)

    def _check_dataset(self, model):
        orig["check_dataset"](self, model)
        _extend_aco2(self, cfg, MemoryPlan, BackendError, IdentityError)

    def _upload_activity(self, state):
        if cfg["activity"] != "exp" or self.activity_mode == "device":
            return orig["upload_activity"](self, state)      # device: calls the patched _device_activity
        act_lin = (self.activity_linear(state) if self.activity_source is None
                   else self.activity_source(state))
        act_lin = self.dataset.gather_faces(act_lin)
        self.arrays["af_face"].set(np.exp(act_lin))
        self.arrays["active"].set(np.full(len(act_lin), EXP_ACTIVE, dtype=np.uint8))

    def _device_activity(self, state, act_out=None):
        if cfg["activity"] != "exp":
            return orig["device_activity"](self, state, act_out)
        # CudaBackend._device_activity with the kernel activity_linear_exp (same arguments)
        cupy, a = self.cupy, self.arrays
        nf, ny = self.dataset.nfaces_device, self.Y.shape[1]
        nb = self.nbasis_enhancement if state.b is not None else 0
        coef = np.concatenate([np.asarray(state.a, dtype=np.float64).reshape(-1),
                               np.asarray(state.b if nb else [], dtype=np.float64).reshape(-1)])
        if len(coef) != ny + nb:
            raise BackendError("activity coefficients: %d values for %d columns" % (len(coef), ny + nb))
        d_coef = cupy.asarray(coef)
        cuda_sources.function("activity", "activity_linear_exp", self.device_id)(
            (max(1, min(65535, -(-nf // 256))),), (256,),
            (a["actbasis"], np.int64(nf), np.int32(ny), np.int32(nb), d_coef,
             np.float64(self.enhancement_scale), a["af_face"], a["active"],
             a["af_face"] if act_out is None else act_out, np.int32(act_out is not None)))

    orig["evaluate"] = B._evaluate
    from cops_cuda import backend as bk_module

    def _evaluate(self, state, request):
        return _evaluate_co2_floor(self, state, request, evaluator, bk_module, BackendError)

    def identity(self):
        value = orig["identity"](self)
        dataset = {k: v for k, v in (getattr(self, "_positive_activity", None) or {}).items()
                   if not k.startswith("_")}
        value["model_extension"] = {"positive_activity": dict(config(), dataset=dataset)}
        return value

    def release(self):
        try:
            orig["release"](self)
        finally:
            info = getattr(self, "_positive_activity", None)
            if info and info.get("_reservation") is not None:
                info.pop("_reservation").release()

    def multi_identity(self):
        value = orig["multi_identity"](self)
        value["model_extension"] = {"positive_activity": config()}
        return value

    # --- spawned workers -------------------------------------------------------------------------------------
    orig["pipe_worker_main"] = transport.pipe_worker_main

    # apply (all at once, after every original was captured)
    cuda_sources.source_path = source_path
    cuda_sources.module = module
    evaluator.HostConstants = PositiveHostConstants
    B.__init__, B._check_dataset = __init__, _check_dataset
    B._upload_activity, B._device_activity = _upload_activity, _device_activity
    if cfg["co2_floor"]:                                            # without it: the class's own _evaluate
        B._evaluate = _evaluate
    B.identity, B.release, B.close = identity, release, release     # close is an alias of release in the class
    multi_gpu.MultiGpuBackend.identity = multi_identity
    transport.pipe_worker_main = _pipe_worker_main
    os.environ[ENV_FROZEN] = _canonical(cfg)
    _STATE["config"] = cfg
    _STATE["installed"] = True
    _log("installed: %s; kernels %s" % (_canonical(cfg), kernel_sha256()))
    return True


def _extend_aco2(backend, cfg, MemoryPlan, BackendError, IdentityError):
    """mesh.ACO2 of the backend's dataset -> [ACO2 | floor_ratio * AH2O] (device facets), once per dataset; with the
    CO2 floor [ACO2 | floor_ratio * AH2O | kappa (NaN until the first evaluation), floor_ratio]."""
    ds, cupy = backend.dataset, backend.cupy
    nf = int(ds.nfaces_device)
    old = ds.arrays.get("mesh.ACO2")
    if old is None or tuple(old.shape) != (nf,):
        raise BackendError("positive_activity: expected mesh.ACO2 of %d device facets, found %s"
                           % (nf, None if old is None else tuple(old.shape)))
    aco2 = np.ascontiguousarray(ds.gather_faces(ds.inputs.elements["ACO2"]), dtype=np.float64)
    if hashlib.sha256(aco2.tobytes()).hexdigest() != ds.host_digests.get("mesh.ACO2"):
        raise IdentityError("positive_activity: ACO2 of the device facets differs from the loader's upload")
    a_orig = np.ascontiguousarray(ds.gather_faces(ds.inputs.elements["AH2O"]), dtype=np.float64)
    if not np.all(np.isfinite(a_orig)) or np.any(a_orig < 0.0):
        raise IdentityError("positive_activity: AH2O of the base has negative or non-finite values")
    co2 = bool(cfg.get("co2_floor"))
    if co2 and (not np.all(np.isfinite(aco2)) or np.any(aco2 < 0.0)):
        raise IdentityError("positive_activity: ACO2 of the base has negative or non-finite values; the CO2 floor "
                            "max(kappa q ACO2, f qmax ACO2) is not defined for them (a user decision)")
    floor_block = cfg["floor_ratio"] * a_orig
    tail = np.array([np.nan, cfg["floor_ratio"]], dtype=np.float64) if co2 else np.zeros(0, dtype=np.float64)
    host = np.concatenate([aco2, floor_block, tail])
    plan = MemoryPlan("positive_activity: H2O floor block of mesh.ACO2")
    plan.add_array("mesh.ACO2.floor_block", "static", (nf + len(tail),), np.float64,
                   "floor_ratio * AH2O of the device facets (second block of mesh.ACO2)"
                   + ("; then kappa and floor_ratio (CO2 floor)" if co2 else ""))
    with cupy.cuda.Device(backend.device_id):
        del ds.arrays["mesh.ACO2"], old
        cupy.get_default_memory_pool().free_all_blocks()
        reservation = backend.budget.reserve(plan)
        try:
            with backend.budget.allocation_guard(plan):
                device = cupy.empty(2 * nf + len(tail), dtype=np.float64)
                device.set(host)
        except BaseException:
            reservation.release()
            raise
    ds.arrays["mesh.ACO2"] = device
    ds.host_digests["mesh.ACO2"] = hashlib.sha256(host.tobytes()).hexdigest()
    backend._positive_activity = {
        "_reservation": reservation, "device_faces": nf, "floor_ratio": cfg["floor_ratio"],
        "mesh_aco2_layout": "[ACO2 (nfaces) | floor_ratio * AH2O (nfaces)]",
        "aco2_sha256": hashlib.sha256(aco2.tobytes()).hexdigest(),
        "ah2o_sha256": hashlib.sha256(a_orig.tobytes()).hexdigest(),
        "floor_block_sha256": hashlib.sha256(floor_block.tobytes()).hexdigest()}
    if co2:
        backend._positive_activity.update(
            mesh_aco2_layout="[ACO2 (nfaces) | floor_ratio * AH2O (nfaces) | kappa (written per evaluation), "
                             "floor_ratio]",
            co2_floor=True, co2_floor_ratio=cfg["floor_ratio"], aco2_facets_zero=int(np.sum(aco2 == 0.0)),
            co2_floor_evaluations=0, co2_floor_first=None)


def _evaluate_co2_floor(self, state, request, evaluator, bk, BackendError):
    """CudaBackend._evaluate (cops_cuda/evaluator.py, the version pinned in ORIGINAL_SHA256) with the CO2 flux floor
    (COPS_POS_FLOOR_CO2=1; pressure.cu copy header item 5). Requests, uploads and the H2O assembly are the original's,
    expression by expression; the CO2 part is p_C,active * kappa + p_C,floor and the kappa column uses the active
    records only."""
    t0 = time.perf_counter()
    c, ds, a = self.constants, self.dataset, self.arrays
    if not self.mixed:
        raise BackendError("positive_activity: the CO2 floor exists only in the mixed-precision kernels")
    kappa = float(state.kappa)
    with self.cupy.cuda.Device(self.device_id):
        self._upload_activity(state)
        nf = int(ds.nfaces_device)
        d_aco2 = ds.arrays["mesh.ACO2"]
        if d_aco2.size != 2 * nf + CO2_FLOOR_SLOTS:
            raise BackendError("positive_activity: mesh.ACO2 has %d values, the CO2 floor needs %d"
                               % (d_aco2.size, 2 * nf + CO2_FLOOR_SLOTS))
        d_aco2[2 * nf:2 * nf + 1].set(np.array([kappa], dtype=np.float64))   # read by eval_fused_m_p5f / p14f
        codes = list(evaluator.PRESSURE_CODES)
        if request.temperature:
            codes += list(evaluator.TEMPERATURE_CODES)
        ncol = len(request.basis_index) + len(request.enhancement_index)
        if ncol:
            self._upload_columns(request)
            codes += [1000 + k for k in range(ncol)] + [2000 + k for k in range(ncol)]
        plan = self._fused_plan(codes)
        if not plan or plan[0][0] not in CO2_FLOOR_KERNELS:
            raise BackendError("positive_activity: unexpected fused plan %s for the CO2 floor" % plan)
        extra = list(CO2_FLOOR_CODES) + (list(CO2_FLOOR_T_CODES) if request.temperature else [])
        name0, col0, codes0 = plan[0]
        plan[0] = (CO2_FLOOR_KERNELS[name0], col0, list(codes0) + extra)
        params = self._params(state)
        groups = [plan]
        if ds.streamed:          # all tiles of one pass share the planned partial buffer (14 + 2 x columns sums)
            total = sum(len(e[2]) for e in plan)
            if (total * int(self._first_block[-1]) > a["partial"].size or total * ds.npoints > a["sums"].size):
                groups = [plan[:1], plan[1:]]
        S = {}
        for group in groups:
            if group:
                gcodes = [k for e in group for k in e[2]]
                S.update(zip(gcodes, self._run_fused(group, gcodes, params)))
    cNG, cRG, bm = c.cNG, c.cRG, c.bm
    # Model.evaluate, expression by expression; CO2: active records x kappa + records at the CO2 floor
    pNG_H = cNG["H2O"] * bm["H2O"] * S[0]
    pRG_H = cRG["H2O"] * (S[1] + bm["H2O"] * S[2])
    pNG_C = cNG["CO2"] * bm["CO2"] * S[3]
    pRG_C = cRG["CO2"] * (S[15] + bm["CO2"] * S[4])
    pNG_F = cNG["CO2"] * bm["CO2"] * S[16]
    pRG_F = cRG["CO2"] * (S[18] + bm["CO2"] * S[17])
    pNG = pNG_H + kappa * pNG_C + pNG_F
    pRG = pRG_H + kappa * pRG_C + pRG_F
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
        def triple(c0):
            s_l, s_lx, s_n = S[c0], S[c0 + 1], S[c0 + 2]
            return c.KT * s_lx + s_l, -c.KT * s_lx, s_n
        tN, tRH, tRC, tNF, tRF = triple(5), triple(8), triple(11), triple(19), triple(22)
        for k in range(3):
            columns[("T", k)] = (tN[k] + cNG["CO2"] * bm["CO2"] * tNF[k],
                                 cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k] + cRG["CO2"] * tRF[k])
    if request.kappa:
        columns[("kappa", 0)] = (pNG_C, pRG_C)          # records at the CO2 floor contribute 0
    info = getattr(self, "_positive_activity", None)
    if info is not None:
        info["co2_floor_evaluations"] = info.get("co2_floor_evaluations", 0) + 1
        if info.get("co2_floor_first") is None:
            with np.errstate(divide="ignore", invalid="ignore"):
                share = np.where(kappa * pNG_C + pNG_F > 0, pNG_F / (kappa * pNG_C + pNG_F), 0.0)
            info["co2_floor_first"] = {"kappa": kappa, "points": int(len(pNG)),
                                       "ng_co2_floor_share_max": float(np.max(share)) if len(share) else None,
                                       "ng_co2_floor_share_median": float(np.median(share)) if len(share) else None,
                                       "ng_co2_over_total_median": float(np.median(
                                           (kappa * pNG_C + pNG_F) / np.where(pNG > 0, pNG, np.inf)))
                                       if len(pNG) else None}
            _log("CO2 floor, first evaluation on device %d: %s" % (self.device_id, info["co2_floor_first"]))
    self.timing["last_evaluate_s"] = time.perf_counter() - t0
    return bk.EvaluationResult(pNG, pRG, columns, request)


def _pipe_worker_main(conn, bind_kwargs, device_id, shm_name, nfaces):
    """Target of the spawned local multi-GPU worker: install() from the inherited environment, then the original loop."""
    try:
        if not install():
            raise RuntimeError("positive_activity: the parent installed the extension but this worker's environment "
                               "has %s != 1" % ENV_ENABLE)
        target = _STATE["originals"]["pipe_worker_main"]
    except BaseException:
        conn.send(("error", "positive_activity install() failed in the worker:\n" + traceback.format_exc()))
        return
    return target(conn, bind_kwargs, device_id, shm_name, nfaces)
