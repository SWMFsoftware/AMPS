#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_fit.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_fit.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 274c89bf7020186d0f6089c7ecfb8ff93bf7308e7cc40906f60767535de15ef9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver. Fit of the softplus / squared models (user decisions 2026-11-11): argument rewriting, model
#   configuration and GPU patch, the stage's held term, --l-threads, nonneg_problem's hooks around grain_fit_fixed.py ->
#   grain_fit.py -> gpu_fit.py -> fit_cops.py; main() restores every patched attribute, so it can run repeatedly in one
#   process (one-job chain).
# USED BY: jobs/nonneg_fit.sbatch (tests 124130-124133, 124139, 124141, 124144, 124146) and chain_onejob.py --kind
#   nonneg (production 124162 softplus, 124165 squared).
# ==== end of header ====
"""Fit driver of the non-negative activity models (softplus / squared harmonics with slope-bin multipliers; user
decisions 2026-11-11). Usage (jobs/nonneg_fit.sbatch builds it):

    nonneg_fit.py --model softplus|squared --point-table T --slope-bundle B --ft-templates-from <D2 fit_manifest.json>
        --start d2|lsq|chain --start-restart R [--start-file S] [--held F --grain-from F] [--l-threads N [--l-mode chunked]]
        [--check-full-jacobian] [--check-parameters ...] [--stop-after-checks] <grain_fit_fixed / gpu_fit options>
        -- <fit_cops arguments>

What it does, in this order (every process the same: no per-stage state survives a call):
 1. nonneg_problem.rewrite_arguments: driver options, inserted options (the 7 bin templates and D2's ft_c* templates,
    --fix c_slope_00_40=1 and c_ft_*=1, --template-start 0), refusals;
 2. the constant model configuration (nonneg_problem.model_config) and nonneg_mult.install_gpu (idempotent; the local
    spawned workers install the same from COPS_NONNEG_FROZEN_CONFIG);
 3. the stage's held term on the host (nonneg_mult.set_stage_held(--held) or None): requirement A of the one-job chain;
 4. --l-threads: nonneg_mult.set_L_threads (idempotent; every 12.4M-facet host step);
 5. nonneg_problem.install_hooks around the unchanged grain_fit_fixed.py -> grain_fit.py -> gpu_fit.py path;
 6. after the fit (also on errors and SystemExit): the attributes the hooks changed in fit_cops, cops_cuda.restart,
    grain_crossing and fixed_parameters are restored, sys.argv too, and the override is cleared, so main() can be
    called again in the same process (one-job chain).
main(argv=None) -> exit code (0 done, 75 wall-time stop, other: failure); errors raise SystemExit as the tools do.
process_config_from_argv(argv) -> the constant configuration (for the one-job chain's process-kind parent).
"""

import hashlib
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
for _p in (str(HERE), str(APP / "tools"), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _split_threads(args_in):
    """Remove --l-threads N [--l-mode chunked] from the driver part; returns (args, threads, mode)."""
    args_in = list(args_in)
    head_end = args_in.index("--") if "--" in args_in else len(args_in)
    threads, mode = None, None
    for opt in ("--l-threads", "--l-mode"):
        if opt in args_in[:head_end]:
            i = args_in.index(opt)
            if i + 1 >= len(args_in):
                raise SystemExit("nonneg: %s needs a value" % opt)
            if opt == "--l-threads":
                threads = int(args_in[i + 1])
            else:
                mode = args_in[i + 1]
            del args_in[i:i + 2]
            head_end -= 2
    if mode is not None and mode != "chunked":
        raise SystemExit("nonneg: --l-mode %r: only the bitwise chunked mode is provided" % mode)
    if mode is not None and threads is None:
        raise SystemExit("nonneg: --l-mode needs --l-threads N")
    return args_in, threads, mode or ("chunked" if threads else None)


def process_config_from_argv(argv):
    """The constant configuration (nonneg_mult.process_config) a fit command would install."""
    import nonneg_problem as npb
    args_in, _, _ = _split_threads(argv)
    cfg, _, _ = npb.rewrite_arguments(args_in)
    return npb.model_config(cfg)[0]


def _snapshot(modules):
    return [(m, dict(vars(m))) for m in modules]


def _restore(snaps):
    for m, saved in snaps:
        cur = vars(m)
        for k in list(cur):
            if k not in saved:
                try:
                    delattr(m, k)
                except AttributeError:
                    pass
        for k, v in saved.items():
            if cur.get(k, object()) is not v:
                setattr(m, k, v)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    command = [str(Path(__file__).resolve())] + argv
    preset = sorted(k for k in os.environ if k.startswith("COPS_POS_") or k == "COPS_SLOPEMUL_FROZEN_CONFIG")
    if preset:
        raise SystemExit("nonneg: %s set: another model patch modifies the same evaluator methods; unset them"
                         % ", ".join(preset))
    import nonneg_mult as nm
    import nonneg_problem as npb
    args_in, threads, mode = _split_threads(argv)
    cfg, gargv, notes = npb.rewrite_arguments(args_in)
    model_cfg, bundle_meta = npb.model_config(cfg)
    nm.install_gpu(model_cfg)
    held = npb.held_of(cfg, model_cfg)
    nm.set_stage_held(cfg["held"], "held_next", model_cfg["nfaces"], expect_sha256=(held or {}).get("sha256"))
    nm.set_override(None)
    th = nm.set_L_threads(threads, "chunked") if threads else nm.set_L_threads(None, "chunked")
    from cops_cuda import reference
    from cops_cuda import restart as rs_mod
    reference.bind()
    fit_cops = reference.load("fit_cops")
    import fixed_parameters
    import grain_crossing
    import grain_fit_fixed
    import clip_problem as cp
    import positive_problem
    snaps = _snapshot([fit_cops, rs_mod, grain_crossing, fixed_parameters])
    argv_saved = list(sys.argv)
    try:
        for n in notes:
            npb.say(n)
        npb.say("held term of this stage: %s; host threads %s" % (nm.stage_held(), th))
        resume = source = arrays = None
        if cfg["restart_from"]:
            resume = npb.check_restart_marker(cfg["restart_from"], cfg, model_cfg, held)
            source, arrays = cp.load_source(cfg["restart_from"])
            source["role"] = "restart_from"
        else:
            source, arrays = cp.load_source(cfg["start_restart"])
            source["role"] = "start_restart"
        sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
               "parameterization": npb.parameterization_of(model_cfg), "description": npb.DESCRIPTION[cfg["model"]],
               "policy": nm.POLICY[cfg["model"]], "model": model_cfg, "held": held,
               "modules": {name: {"path": str(Path(m.__file__).resolve()), "sha256": sha256(m.__file__)}
                           for name, m in (("nonneg_mult", nm), ("nonneg_problem", npb), ("slope_mult", nm.sm),
                                           ("clip_problem", cp), ("positive_problem", positive_problem),
                                           ("grain_fit_fixed", grain_fit_fixed), ("fixed_parameters", fixed_parameters),
                                           ("grain_crossing", grain_crossing))},
               "slope_bundle": {"path": cfg["slope_bundle"], "sha256": sha256(cfg["slope_bundle"]),
                                "slope_multiplier_bins": bundle_meta.get("slope_multiplier_bins")},
               "config": dict(cfg), "inserted_or_replaced_arguments": notes, "grain_fit_fixed_argv": list(gargv),
               "resume": resume, "gpu_patch_environment": nm.ENV_FROZEN, "threads": th}
        sub["start_choice"] = cfg["start_choice"]
        sub["fit_chain"] = npb.fit_chain(source, resume, cfg["start_choice"])
        npb.install_hooks(grain_crossing, fixed_parameters, cfg, model_cfg, command, sub, source, arrays, held)
        sys.argv = [str(HERE / "grain_fit_fixed.py")] + gargv
        rc = grain_fit_fixed.main()
        return 0 if rc is None else rc
    finally:
        sys.argv = argv_saved
        _restore(snaps)
        grain_crossing.install = getattr(grain_crossing, "_nonneg_original_install", grain_crossing.install)
        fixed_parameters.install = getattr(fixed_parameters, "_nonneg_original_install", fixed_parameters.install)
        nm.set_override(None)


if __name__ == "__main__":
    sys.exit(main())
