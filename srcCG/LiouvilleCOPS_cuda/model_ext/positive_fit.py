#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5ea0230d74027b24f83da5a2c6f9371f2281ebff8147cae6bdc4a29996f669a8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Positive-activity fit (user decisions 2026-11-07); options before "--" (--point-table, --grain-start-uPa,
#   --crossing-weight, --check-jacobian, --devices, --precision mixed, --max-nfev, --segment, --wall-time,
#   --restart-from, --template-start) and the gpu_fit.py/fit_cops.py arguments after it.
# USED BY: jobs/positive_fit.sbatch (jobs 122307, 122310-122312, 122319, 122320), jobs/positive_fit_test.sbatch (122303).
# RUN: see the commands.sh of each fit's provenance (figures_cluster/12M/model_vs_observations_positive_lambda_scan/HOW-TO-RERUN.md).
# ==== end of header ====
"""grain_fit.py with the positive H2O activity A = exp(s), the start from the original map, the log-space Tikhonov term
and the COPS W10x5 time response in the late window (model_ext/positive_problem.py; GPU side model_ext/
positive_activity.py; user decisions 2026-11-07).

    .venv/bin/python <this> --point-table <run>/late_grid_points.txt [--grain-start-uPa A] [--check-jacobian]
        [--check-parameters a00,a02,...] <grain_fit.py / gpu_fit.py options> -- <fit_cops arguments>

Example (8 GPUs, the 122250 command adapted; no record store on this base: --precision mixed converts the kernels at
bind, ~12 min; one --lambda per process):

    R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
    B=/data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz
    .venv/bin/python -W ignore /data001/highresolu_vt/cops_output/model_ext/positive_fit.py \\
        --point-table $R/late_grid_points.txt --grain-start-uPa 6.96024299190473 --crossing-weight 10 \\
        --check-jacobian --devices 0,1,2,3,4,5,6,7 --precision mixed --max-nfev 4000 --wall-time slurm \\
        --template-start 0 -- $R/base --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 \\
        --jobs 1 --fix-temperature --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance inf \\
        --out $R/positive/<name>/fit --activity-template $B cliff_ge_60deg

Driver options (before '--'):
* --point-table: late_time_grid.py --table of the base (index UTC kind original_index L2_start L2_stop); grid rows and
  the L2 STOP times of the late observations (required);
* --grain-start-uPa: grain A start in 1e-6 Pa (default 6.96024299190473 = cycle 3's fine-tuned A, user 2026-11-07);
  grain B starts at 0;
* --check-jacobian: before fitting, analytic vs central-difference Jacobian columns of the full objective (data with
  grain and W10x5, crossing rows, Tikhonov rows if lambda > 0) at the start for --check-parameters (default a00, a02,
  a07, c_cliff_ge_60deg, kappa, grain_A_uPa); printed, recorded in fit_manifest.json["model_extension"]
  ["positive"]["jacobian_check"]; the fit continues.
Inserted when absent (stated in the record): --template-start 0, --fix-temperature, --check-tolerance inf (the
fitter's consistency check with the run's linear coefficients cannot pass with exp activity; replaced by the GPU
patch's checks and --check-jacobian). Required: exactly one --activity-template <bundle> cliff_ge_60deg. Refused:
options that contradict the user's model (positive_problem.rewrite_arguments). --restart-from (gpu_fit) is accepted
only for a restart written by this driver (resume after a wall-time stop, jobs/positive_fit.sbatch).

The unchanged grain_fit.py runs (grain + crossing, grain_crossing.py); this driver wraps grain_crossing.install:
positive_problem.install_core first (selection of the grid rows, exponent Problem, activity writers, Jacobian check;
inside the crossing wrapper), then the grain problem, then positive_problem.install_outer (grain start, late-start
check, activity_exp_facets.npz / ACTIVITY_PARAMETERIZATION.txt, restart marker). fit_manifest.json
["model_extension"]["command"] is this driver's full command; ["positive"] holds the start projection, the window
operator, the regularization, the GPU patch configuration (positive_activity.config(), SHA-256s) and the checks.

Module level (also in spawned GPU worker processes, which re-import this main module): the COPS_POS_* environment of
the GPU patch is set and positive_activity.install() is called once per process, before anything binds a backend.
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

import positive_problem  # noqa: E402

POSITIVE_ENVIRONMENT = positive_problem.set_positive_environment()     # before any GPU binding (user 2026-11-07)

import positive_activity  # noqa: E402

if not getattr(positive_activity, "_installed_by_positive_fit", False):   # one install per process
    positive_activity.install()
    positive_activity._installed_by_positive_fit = True


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    command = [str(Path(__file__).resolve())] + sys.argv[1:]
    cfg, argv, notes = positive_problem.rewrite_arguments(sys.argv[1:])
    for n in notes:
        print("positive_fit: %s" % n, flush=True)
    restart = positive_problem.check_restart_marker(cfg["restart_from"]) if cfg["restart_from"] else None
    import grain_crossing
    import grain_fit
    original = grain_crossing.install

    def install(fit_cops, weight, record):
        sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
               "module": str(Path(positive_problem.__file__).resolve()),
               "module_sha256": sha256(positive_problem.__file__),
               "positive_activity": positive_problem.describe_positive_activity(positive_activity),
               "environment": dict(POSITIVE_ENVIRONMENT),
               "environment_in_process": {k: os.environ.get(k) for k, _ in positive_problem.POSITIVE_ENV},
               "inserted_or_replaced_arguments": notes, "grain_fit_argv": list(argv),
               "restart": restart, "config": {k: v for k, v in cfg.items()},
               "consistency_check": {
                   "fitter_check_tolerance": "inf",
                   "reason": "fit_cops evaluates MF at the run's linear srcCG coefficients; with A = exp(s) these "
                             "overflow, so the printed 'consistency with LiouvilleCOPS' numbers and base_cache."
                             "python_cpp_max_relative_difference of this manifest are meaningless",
                   "replaced_by": ["positive_activity.py checks (GPU patch)", "positive_fit.py --check-jacobian"]}}
        record["positive"] = sub
        positive_problem.install_core(fit_cops, cfg, sub)
        result = original(fit_cops, weight, record)
        positive_problem.install_outer(fit_cops, cfg, sub)
        record["grain_fit_command"] = record.get("command")
        record["command"] = command
        return result
    grain_crossing.install = install
    sys.argv = [str(HERE / "grain_fit.py")] + argv
    return grain_fit.main()


if __name__ == "__main__":
    sys.exit(main())
