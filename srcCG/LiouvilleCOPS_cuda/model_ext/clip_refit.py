#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/clip_refit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/clip_refit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 347862d2b78986ea2b81bffc833815bfa51faf2158d874d6dc12aefd30d1fed4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Diagnostic refit of the earlier clipped model (user 2026-11-08: "refit the earlier clipped model (original
#   temperature, with its templates) on the new base, with W10x5 and the +1.5 m offset") around the unchanged
#   grain_fit_fixed.py -> grain_fit.py -> gpu_fit.py -> fit_cops.py; --variant D1 (observation time stamps) | D2 (+ W10x5
#   of positive_problem.py for the observations from 10:33:10); --start-restart (the c6 fit's final state by name).
# USED BY: jobs/clip_refit.sbatch (D0 122880, D1 122865, D2 122866), jobs/clip_refit_bitwise.sbatch (122854).
# ==== end of header ====
"""Refit of the earlier clipped-linear model c6 on the lategrid base, variant D1 (instantaneous) or D2 (W10x5)
(user decision 2026-11-08; model_ext/clip_problem.py has the model, the grid-row treatment and the checks).

    .venv/bin/python <this> --variant D1|D2 --point-table <run>/late_grid_points.txt
        [--start-restart <c6 fit_restart.npz>] [--check-full-jacobian] [--check-parameters a00,a02,...]
        <grain_fit_fixed.py options: --fix ...> <grain_fit.py options: --crossing-weight 10 --check-jacobian>
        <gpu_fit.py options> -- <fit_cops arguments>

The c6 command (job 122221) with the new base (8 GPUs; no record store on this base: --precision mixed converts the
kernels at bind; one --lambda per process):

    R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
    E=/data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier
    C6=$E/grain/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6
    .venv/bin/python -W ignore /data001/highresolu_vt/cops_output/model_ext/clip_refit.py --variant D2 \\
        --point-table $R/late_grid_points.txt --start-restart $C6/fit_restart.npz --check-full-jacobian \\
        --fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0 --fix c_ft_c1=1 ... --fix c_ft_c6=1 \\
        --crossing-weight 10 --check-jacobian --devices 0,1,2,3,4,5,6,7 --precision mixed --max-nfev 4000 \\
        --segment 1 --wall-time slurm --template-start 0 -- $R/base --profile all_available --degree 4 \\
        --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \\
        --check-tolerance 1e-6 --out $R/clip_refit/D2/seg01 --activity-template $E/bundle/restart_templates.npz \\
        cliff_ge_60deg --activity-template $E/grain/cycles/c1/template_ft_c1/ft_c1.templates.npz ft_c1 ...

(jobs/clip_refit.sbatch writes the full command, every argument resolved, into provenance-<job>/commands.sh.)

Driver options (before '--'):
* --variant D1|D2 (required): D1 = the c6 path, model at the observation time stamps; D2 = D1 + W10x5 for the
  observations from 2016-09-30T10:33:10 UTC exactly as positive_problem.py (its own operator functions);
* --point-table: late_time_grid.py --table of the base (grid rows, L2 STOP times); required for D2 and on a base with
  grid rows (the named profile audits the 1150-row mission dataset; the grid rows get use 0, positive_problem's
  selection wrapper). Without it (a base without grid rows, e.g. the old base in the bitwise test) the driver adds
  nothing to the grain_fit_fixed path except the start restart, the start check, the optional full Jacobian check and
  the records/marker;
* --start-restart R.npz: a NEW fit from R's final state (by name, exact optimizer coordinates; mesh, names, template
  columns, held values, scales checked; prior unchanged); excludes gpu_fit --restart-from (= resume of this fit,
  jobs/clip_refit.sbatch segments; marker checked; the start check must reproduce the restart's cost to 1e-8);
* --check-full-jacobian, --check-parameters: analytic vs central-difference Jacobian columns of the full objective at the
  optimizer start (data incl. the late window, crossing rows); printed, recorded in fit_manifest.json
  ["model_extension"]["clip_refit"]["jacobian_check"]; the fit continues.
Inserted when absent: --template-start 0 (c6's value). Refused: COPS_POS_* in the environment (the positive-activity
GPU patch; this model has the srcCG night/shadow ratios and no floor), fitter options that contradict the c6 path.

The unchanged grain_fit_fixed.py runs (-> grain_fit.py -> gpu_fit.py -> fit_cops.py); this driver only installs
clip_problem.install_hooks: install_core just before grain_crossing.install (selection of the grid rows, product
filters, the W10x5 core for D2, the start check and the Jacobian check inside the crossing wrapper) and install_outer
just after fixed_parameters.install (the start restart by name, grain-window records, late_window_model.dat,
CLIP_REFIT.txt, the restart marker). fit_manifest.json["model_extension"]["command"] is this driver's full command.
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
        raise SystemExit("clip_refit: %s set: this model runs the unchanged GPU physics (srcCG night/shadow ratios, no "
                         "floor); unset them" % ", ".join(preset))
    import clip_problem as cp
    cfg, argv, notes = cp.rewrite_arguments(sys.argv[1:])
    for n in notes:
        cp.say(n)
    resume = source = arrays = None
    if cfg["restart_from"]:
        resume = cp.check_restart_marker(cfg["restart_from"], cfg)
        source, arrays = cp.load_source(cfg["restart_from"])
        source["role"] = "restart_from"
    elif cfg["start_restart"]:
        source, arrays = cp.load_source(cfg["start_restart"])
        source["role"] = "start_restart"
    import fixed_parameters
    import grain_crossing
    import grain_fit_fixed
    import positive_problem
    sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
           "module": str(Path(cp.__file__).resolve()), "module_sha256": sha256(cp.__file__),
           "uses": {"positive_problem": {"module": str(Path(positive_problem.__file__).resolve()),
                                         "sha256": sha256(positive_problem.__file__),
                                         "functions": ["read_point_table", "grid_keys", "read_run_points", "LateLayout",
                                                       "w10x5_operator", "WindowOperator", "install_window",
                                                       "install_selection", "install_figure_filter", "time_key"]},
                    "grain_fit_fixed": {"module": str(Path(grain_fit_fixed.__file__).resolve()),
                                        "sha256": sha256(grain_fit_fixed.__file__)},
                    "fixed_parameters": {"module": str(Path(fixed_parameters.__file__).resolve()),
                                         "sha256": sha256(fixed_parameters.__file__)},
                    "grain_crossing": {"module": str(Path(grain_crossing.__file__).resolve()),
                                       "sha256": sha256(grain_crossing.__file__)}},
           "config": dict(cfg), "inserted_or_replaced_arguments": notes, "grain_fit_fixed_argv": list(argv),
           "resume": resume, "environment_COPS_POS": None}
    sub["fit_chain"] = cp.fit_chain(cfg, source if cfg["start_restart"] else None, resume)
    cp.install_hooks(grain_crossing, fixed_parameters, cfg, command, sub, source, arrays)
    sys.argv = [str(HERE / "grain_fit_fixed.py")] + argv
    return grain_fit_fixed.main()


if __name__ == "__main__":
    sys.exit(main())
