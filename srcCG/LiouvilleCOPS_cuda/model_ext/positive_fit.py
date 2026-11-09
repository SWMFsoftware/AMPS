#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8c1139f5c3d822e28259b7d78ce251e0c585681991d0e0b6754a459b8a91d85f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Positive-activity fit (user decisions 2026-11-07/08); options before "--" (--point-table, --grain-start-uPa,
#   --grain-start-from, --crossing-weight, --check-jacobian, --devices, --precision mixed, --max-nfev, --segment,
#   --wall-time, --restart-from (resume), --start-restart (start from another fit), --ratio 0.02|0.002, --held-exponent
#   (per-facet exponent term held, D on the total ln A), --check-held-span, --template-start) and the gpu_fit.py
#   arguments after it. 2026-11-08: ratio, start restart with start check, held term.
# USED BY: jobs/positive_fit.sbatch (jobs/positive_cycles.sh).
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
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
* --ratio R: night/shadow production ratio = H2O flux-floor ratio (COPS_POS_NIGHT_RATIO = COPS_POS_FLOOR_RATIO = R;
  user decisions: 0.02 and, since 2026-11-07 evening, 0.002; positive_problem.RATIO_DECISIONS). Default: the preset
  environment (jobs/positive_fit.sbatch RATIO), else 0.02. Recorded in the manifest, the products and the restart marker;
  a resume (--restart-from) must have its restart's ratio;
* --start-restart R.npz: start the optimizer from another positive fit's fit_restart.npz (a, c, kappa, grain by name;
  exact optimizer coordinates when the ascale is the same; the prior s0 and D's reference are unchanged; the source's
  temperature law, mesh, template and parameter layout must be this fit's). A new fit, not a resume: its wall-time
  segments start at 1. Before the first optimizer step the START CHECK compares the cost at the start with the source's
  final cost: same ratio, lambda and held term -> must reproduce it to 1e-8 relative (else the fit stops); otherwise the
  difference and its causes are printed and recorded (model_extension.positive.start_check). Excludes --restart-from;
* --start-set NAME=VALUE[,NAME=VALUE...] (2026-11-08; user "kappa = 1 where stuck": e.g. --start-set kappa=1): after
  --start-restart (required) the optimizer start of each named parameter of the fit is VALUE in its optimizer coordinate
  (kappa: kappa itself; positive_problem.START_SET_DEFINITION); lo <= VALUE <= hi; the prior is unchanged; not with
  --grain-start-* for the grain names. Recorded in the manifest (model_extension.positive.start_set, config,
  fit_chain), the restart marker (model_extension.start_set and fit_chain.start_set, kept by resumes),
  ACTIVITY_PARAMETERIZATION.txt ("Optimizer start: ... then --start-set ...") and the start check, which reports the
  start as NOT the source's state. Absent: every behaviour as before;
* --grain-start-uPa A, --grain-start-B B: grain start (A in 1e-6 Pa, B in 1e-6 Pa per 100 s). Precedence: these, else
  the start restart's grain, else A = 6.96024299190473 (cycle 3's fine-tuned A, user 2026-11-07), B = 0;
* --grain-start-from F.npz: the grain start from a fine tuning's state file (keys grain_A_uPa, grain_B_uPa_per_100s;
  path and sha256 recorded); excludes --grain-start-uPa/-B;
* --held-exponent H.npz: a per-facet term held fixed in the exponent (user 2026-11-07, cycles 2 and 3: the accumulated
  fine-tuning change): s = Y a + c T + h. Key "held_next" (the fine tuning's file) or "held" (tests); float64, one
  value per facet in mesh order, nepers. The GPU sees the full s (act_lin + h on the host: positive_problem.
  install_held_activity; host activity of the local transport only), the Jacobian is unchanged in form. Tikhonov:
  D^2 = |B (x - x0) + h|^2_area = |Lt (x - x0 + h_B)|^2 + |h_perp|^2_area (h = B h_B + h_perp), i.e. D measures the
  TOTAL ln A (fitted + held) from s0; the regularization rows are Lt (x - x0 + h_B) and the constant
  sqrt(lambda) |h_perp|_area, so the cost is 1/2 (data + lambda D^2) with the total D. Products: activity_exp_facets.npz
  "held" (s includes h), ACTIVITY_PARAMETERIZATION.txt (file, sha256, D contribution), restart marker (path, sha256;
  activity_linear includes h). A resume (--restart-from) must have the same held file (sha256);
* --check-held-span (diagnostic, tests): with a held term in the span of the fitted columns, data residuals/Jacobian
  with h at the start x == without h at x + h_B (evaluated on the GPUs both ways); recorded as held_span_check;
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



def _ratio_argument(argv):
    """--ratio R before '--' (the driver's options), or None. Spawned GPU workers get the parent's sys.argv or none;
    they inherit the environment the parent set, which set_positive_environment then checks."""
    head = argv[:argv.index("--")] if "--" in argv else argv
    for k, token in enumerate(head):
        if token == "--ratio" and k + 1 < len(head):
            return head[k + 1]
    return None


# before any GPU binding (user 2026-11-07); one ratio for night/shadow and the H2O floor
POSITIVE_ENVIRONMENT = positive_problem.set_positive_environment(ratio=_ratio_argument(sys.argv[1:]))

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
    held_record = {}
    if cfg["held_exponent"]:                         # before anything binds a backend
        run = Path(cfg["run"])
        if not (run / "elements.bin").is_file():
            from cops_cuda import reference
            run = Path(reference.load("cops_paths").product_input(str(run)))
        positive_problem.load_held(cfg["held_exponent"], positive_problem.elements_faces(run))
        positive_problem.install_held_activity(held_record)
        if positive_problem.HELD["info"]["sha256"] != cfg["held_identity"]["sha256"]:
            raise SystemExit("positive_fit: the held file changed while it was read")
    restart = positive_problem.check_restart_marker(cfg["restart_from"], cfg) if cfg["restart_from"] else None
    start_restart = (positive_problem.check_restart_marker(cfg["start_restart"], cfg) if cfg["start_restart"]
                     else None)
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
               "restart": restart, "start_restart_marker": start_restart, "config": {k: v for k, v in cfg.items()},
               "held_gpu_patch": held_record.get("held_gpu_patch"),
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
