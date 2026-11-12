#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_mult_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_mult_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0bc323dae028e9be1d5234c85ca2288a9cb653f903904fafbac942b15e1f7edf
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Slope-multiplier H2O activity fit (user decisions 2026-11-10): rewrites the arguments (the 7 bin
#   templates, c_slope_00_40 = 1 held, --template-start 0), installs slope_mult's GPU patch in every process (local
#   workers or MPI ranks) and slope_mult_problem's hooks around grain_fit_fixed.py -> grain_fit.py -> gpu_fit.py ->
#   fit_cops.py; --start required; --l-threads N (chunked bitwise L = Y a, user 2026-11-10 "Bitwise chunked version");
#   --bump-bundle (option 3a).
# USED BY: jobs/slope_mult_fit.sbatch; jobs/mpi_fit.sbatch with mpi_multinode/specs/slope_mult.sh. This version
#   (snapshots slopemul_20261110T191049Z ... v8): chunked test 123979, 3a chain fits 124008, 124069, 124078; the
#   harmonic fit 123786 ran an earlier version.
# ==== end of header ====
"""Slope-multiplier H2O activity fit (user decisions 2026-11-10): A = max(1e-5 A_orig, m_b L), one multiplier per
effective-gravitational-slope bin, on the lategrid base with D2's treatment of the late window (W10x5).

    .venv/bin/python -W ignore slope_mult_fit.py --point-table <run>/late_grid_points.txt --slope-bundle <bins.npz> \\
        (--start srccg|d2 --start-restart <D2 fit_restart.npz> | gpu_fit --restart-from <seg<k>/fit_restart.npz>) \\
        [--check-full-jacobian] [--check-parameters a,b,...] \\
        --fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0 \\
        --crossing-weight 10 [--check-jacobian] <gpu_fit.py options> -- <fit_cops arguments without --activity-template>

The model (slope_mult.py), the fit-level hooks (slope_mult_problem.py) and the bin bundle (slope_bins.py) are documented
in their headers. This driver, in EVERY process that runs it (the one local process, or every MPI rank started with
the same command line):
 1. rewrites the arguments (slope_mult_problem.rewrite_arguments: inserted --template-start 0, --fix c_slope_00_40=1,
    the 7 --activity-template columns of --slope-bundle in bin order; with --transport mpi: --activity host);
 2. builds the model configuration from the bundle's metadata (floor ratio, bins digest, A_orig digest) and installs
    the GPU patch (slope_mult.install_gpu) BEFORE any backend binds; local spawned workers install it themselves from
    the environment variable COPS_SLOPEMUL_FROZEN_CONFIG the patch sets; with --transport mpi every rank installs from
    its own command line and the configurations are compared across ranks (slope_mult.mpi_config_check, allgather);
 3. MPI ranks >= 1: run the unchanged grain_fit_fixed.py chain, which reaches gpu_fit.main and serves rank 0
    (transport.mpi_worker_main); nothing else;
 4. rank 0 / the local process: installs slope_mult_problem's hooks around the unchanged grain_fit_fixed.py ->
    grain_fit.py -> gpu_fit.py -> fit_cops.py chain (as clip_refit.py does with clip_problem).
--start (required for a new fit, no default; user decision 2026-11-10): srccg = the run's srcCG harmonic coefficients
(the fitter's prior x0), every free multiplier 1, kappa 1, grain A and B from --start-restart (D2) by name (production);
d2 = D2's a, kappa and grain by name, every free multiplier 1 (the first version). Recorded in the manifest
(model_extension.slope_multiplier.start_choice, start_restart, start_report), the restart marker, the start check line
and SLOPE_MULTIPLIER_MODEL.txt; a resume refuses a --start that differs from the fit's.
--l-threads N [--l-mode chunked|openblas] (optional; user 2026-11-10 "Bitwise chunked version"): the FITTER process (rank
0 / the local process) computes L = Y a (its main CPU work) with N threads: mode chunked (default) = N Python threads,
each Y[chunk] @ a with OpenBLAS at 1 thread on 4096-row-aligned chunks, bitwise the 1-thread L; mode openblas (explicit,
non-default, not bitwise) = threaded OpenBLAS for that call. Absent: unchanged (one call, 1 thread). --blas-threads (the
earlier test option) is refused. Recorded in fit_manifest.json model_extension.slope_multiplier.l_step, the start-check
line, the per-call L timing (registry) and the backend identities (worker processes: OpenBLAS 1 thread).
Refused: COPS_POS_* in the environment (positive_activity.py patches the same methods).
fit_manifest.json["model_extension"]["command"] is this driver's full command; ["slope_multiplier"] its records.
"""

import hashlib
import os
import sys
from pathlib import Path

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP / "tools"), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    command = [str(Path(__file__).resolve())] + sys.argv[1:]
    preset = sorted(k for k in os.environ if k.startswith("COPS_POS_"))
    if preset:
        raise SystemExit("slope_mult: %s set: the positive-activity patch modifies the same evaluator methods; unset "
                         "them" % ", ".join(preset))
    import slope_mult as sm
    import slope_mult_problem as smp
    args_in = list(sys.argv[1:])
    head_end = args_in.index("--") if "--" in args_in else len(args_in)
    if "--blas-threads" in args_in[:head_end]:
        raise SystemExit("slope_mult: --blas-threads (threaded OpenBLAS, the 2026-11-10 test) was replaced (user "
                         "2026-11-10: 'Bitwise chunked version'): use --l-threads N (chunked, bitwise) or, explicitly, "
                         "--l-threads N --l-mode openblas (not bitwise)")
    l_threads, l_mode = None, None                   # --l-threads N [--l-mode chunked|openblas]; absent = unchanged
    for opt in ("--l-threads", "--l-mode"):
        if opt in args_in[:head_end]:
            i = args_in.index(opt)
            if i + 1 >= len(args_in):
                raise SystemExit("slope_mult: %s needs a value" % opt)
            if opt == "--l-threads":
                l_threads = int(args_in[i + 1])
            else:
                l_mode = args_in[i + 1]
            del args_in[i:i + 2]
            head_end -= 2
    if l_mode is not None and l_threads is None:
        raise SystemExit("slope_mult: --l-mode needs --l-threads N")
    if l_threads is not None:
        l_mode = l_mode or "chunked"
        if l_threads < 1 or l_mode not in sm.L_MODES:
            raise SystemExit("slope_mult: --l-threads must be >= 1 and --l-mode one of %s" % (sm.L_MODES,))
    cfg, argv, notes = smp.rewrite_arguments(args_in)
    cfg["l_threads"], cfg["l_mode"] = l_threads, l_mode
    model_cfg, bundle_meta = smp.model_config(cfg)
    comm, rank, ranks = None, 0, None
    if cfg["transport"] == "mpi":
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    sm.install_gpu(model_cfg)
    if comm is not None:
        ranks = sm.mpi_config_check(model_cfg, comm)
    import fixed_parameters
    import grain_crossing
    import grain_fit_fixed
    if rank != 0:                                   # an MPI worker rank: the chain reaches gpu_fit.main and serves
        sys.argv = [str(HERE / "grain_fit_fixed.py")] + argv
        return grain_fit_fixed.main()
    for n in notes:
        smp.say(n)
    l_step = {"mode": l_mode, "threads": l_threads, "policy": None}
    if l_threads is not None:
        l_step["setup"] = sm.set_L_threads(l_threads, l_mode)
        l_step["policy"] = (
            "user 2026-11-10 ('Bitwise chunked version'): L = Y a of the FITTER process with %d Python threads, each "
            "Y[chunk] @ a with numpy's OpenBLAS at 1 thread (chunks aligned to %d rows), bitwise the 1-thread L; OpenBLAS "
            "stays at 1 thread everywhere (CLAUDE.md section 4); GPU worker processes unchanged" % (l_threads,
                                                                                                   sm.L_ALIGN_ROWS)
            if l_mode == "chunked" else
            "explicit non-default mode openblas: numpy's OpenBLAS runs L = Y a with %d threads (restored after every "
            "call; departure from CLAUDE.md section 4 for this step; NOT bitwise)" % l_threads)
    l_step["fitter_process_blas"] = sm.blas_info()
    smp.say("L = Y a of the fitter process: %s (OpenBLAS %s thread(s) outside the step; env OPENBLAS_NUM_THREADS=%s)" % (
        "single call, 1 thread (unchanged)" if l_threads is None else
        "%s, %d threads%s" % (l_mode, l_threads, " (%d-row chunks)" % sm.L_ALIGN_ROWS if l_mode == "chunked" else ""),
        l_step["fitter_process_blas"]["threads"], os.environ.get("OPENBLAS_NUM_THREADS")))
    import clip_problem as cp
    import positive_problem
    resume = source = arrays = None
    if cfg["restart_from"]:
        resume = smp.check_restart_marker(cfg["restart_from"], cfg, model_cfg)
        source, arrays = cp.load_source(cfg["restart_from"])
        source["role"] = "restart_from"
    else:
        source, arrays = cp.load_source(cfg["start_restart"])
        source["role"] = "start_restart"
    sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
           "parameterization": smp.PARAMETERIZATION_TEST if model_cfg["test_only"] else smp.PARAMETERIZATION,
           "description": smp.DESCRIPTION, "policy": sm.POLICY, "model": model_cfg,
           "modules": {name: {"path": str(Path(m.__file__).resolve()), "sha256": sha256(m.__file__)}
                       for name, m in (("slope_mult", sm), ("slope_mult_problem", smp), ("clip_problem", cp),
                                       ("positive_problem", positive_problem), ("grain_fit_fixed", grain_fit_fixed),
                                       ("fixed_parameters", fixed_parameters), ("grain_crossing", grain_crossing))},
           "slope_bundle": {"path": cfg["slope_bundle"], "sha256": sha256(cfg["slope_bundle"]),
                            "slope_multiplier_bins": bundle_meta.get("slope_multiplier_bins")},
           "config": dict(cfg), "inserted_or_replaced_arguments": notes, "grain_fit_fixed_argv": list(argv),
           "resume": resume, "mpi_ranks": ranks, "gpu_patch_environment": sm.ENV_FROZEN, "l_step": l_step}
    sub["start_choice"] = cfg["start_choice"]
    sub["fit_chain"] = smp.fit_chain(source, resume, cfg["start_choice"])
    smp.install_hooks(grain_crossing, fixed_parameters, cfg, model_cfg, command, sub, source, arrays)
    sys.argv = [str(HERE / "grain_fit_fixed.py")] + argv
    return grain_fit_fixed.main()


if __name__ == "__main__":
    sys.exit(main())
