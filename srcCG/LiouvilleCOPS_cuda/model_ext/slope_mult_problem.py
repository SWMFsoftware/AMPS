# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_mult_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_mult_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ba22ec7cb9f7d1f5ab7362942beb47096a95ac01d1ccf504550d9dc113c75f4b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Fit-level hooks of the slope-multiplier model around the unchanged grain_fit_fixed.py path: D2's
#   treatment through clip_problem.install_core (grid rows, W10x5 from 10:33:10 UTC, start check), the start (--start
#   srccg | d2 | nnls | chain; user 2026-11-10: srccg for the harmonic fit, nnls for option 3a), model checks,
#   diagnostics, products and restart marker.
# USED BY: slope_mult_fit.py. This version (snapshots v6 slopemul_20261110T201212Z ... v8): the 3a chain 124008-124086;
#   the harmonic fit 123786 ran an earlier version (slopemul_20261110T094013Z).
# ==== end of header ====
"""Fit-level part of the slope-multiplier model (user decisions 2026-11-10); the model and the GPU patch are
model_ext/slope_mult.py, the driver slope_mult_fit.py, the bin bundle slope_bins.py.

Local, untracked extension: fit_cops.py, gpu_fit.py, grain_fit.py, grain_crossing.py, grain_fit_fixed.py,
fixed_parameters.py, positive_problem.py and clip_problem.py are NOT edited. slope_mult_fit.py runs the unchanged
grain_fit_fixed.py path and installs this module around it, exactly where clip_refit.py installs clip_problem:
* install_core, just before grain_crossing.install: clip_problem.install_core with the D2 configuration (variant D2:
  grid-row selection, product filters, the W10x5 core of positive_problem for the observations from 10:33:10 UTC,
  grain unaveraged at the time stamp, the start check and the optional full Jacobian check) -- the D2 treatment
  verbatim, by calling clip_problem's own function; plus the per-evaluation log and the test-only stop;
* install_outer, just after fixed_parameters.install: the start (--start, required for a new fit, START_CHOICES:
  srccg = the run's srcCG harmonic coefficients (the fitter's prior x0), every free multiplier 1, kappa 1, grain A and B
  from D2 by name (production, user 2026-11-10); d2 = D2's a, kappa, grain by name, every free multiplier 1), the start
  report (start cost and the share of facets / area / seen area at the floor at the start), the model checks (template names/columns = the bin indicators of the configured bins; held values), the
  grain-window record (D2's check that the first post-gap NG point is the first late observation), this model's
  activity diagnostics, the products and the restart marker.

Everything else is D2's command (jobs/slope_mult_fit.sbatch): lategrid base (1150 observations + 401 model-only 1-s grid
rows, gauge +1.5 m), original temperature law held (Tmax 182.09922, T75 136.11149, Tmin 133 K, --fix), kappa free,
grain + crossing weight 10, lambda 0, --tbounds 30 400, ftol/xtol/gtol 1e-10, max_nfev 4000; NO cliff_ge_60deg and NO
ft_c1..ft_c6 templates: the 7 slope-bin indicators (slope_bins.py) are the fitter's templates, their amplitudes the
multipliers (c_slope_00_40 held at 1).
"""

import datetime as dt
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

import clip_problem as cp          # D2's core (grid rows, W10x5, start and Jacobian checks), used unchanged
import positive_problem as pp      # point-table helpers (via clip_problem), used unchanged
import slope_mult as sm

PARAMETERIZATION = "slope_multiplier_floor_w10x5_v1"
PARAMETERIZATION_TEST = "slope_multiplier_floor_w10x5_TEST_v1"   # floor ratio 0 (GPU check (a) only)
DESCRIPTION = ("H2O activity A = max(1e-5 A_orig, m_b L) with one multiplier per effective-gravitational-slope bin "
               "(m_[0,40) = 1 held), L = the 25 signed degree-<=4 harmonics; lategrid base (+1.5 m gauge offset); W10x5 "
               "for the observations from %s UTC (positive_problem's operator, as D2), grain unaveraged at the time "
               "stamp" % pp.LATE_START_UTC)
USER_T0 = (182.09922, 136.11149, 133.0)
GRAIN_NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
CHECK_PARAMETERS = tuple("c_" + n for n in sm.BIN_NAMES if n != sm.FIXED_BIN) + (
    "a00", "a02", "a07", "a12", "kappa", "grain_A_uPa")
START_MULTIPLIER = 1.0             # the lead's brief: all multipliers start at 1
# Optimizer start (--start, required for a new fit; user decision 2026-11-10 on the READY file's open question 1):
START_CHOICES = {
    "srccg": ("user 2026-11-10 (production): the ORIGINAL srcCG harmonic coefficients (the run's map = the fitter's prior "
              "x0, a_j = coef0_j / ascale), every free multiplier 1 (0-40 deg held at 1), kappa 1, grain A and B from "
              "D2's final state (--start-restart, by name); no other D2 parameter"),
    "d2": ("the lead's brief of 2026-11-10 (first version): D2's final a00..a24, kappa and grain by name (--start-restart), "
           "every free multiplier 1"),
    "nnls": ("user 2026-11-10 (bump model, option 3a): bump weights = the non-negative least-squares fit of L to the "
             "original srcCG activity A_orig (area weighted, all facets; slope_bumps.py, bundle metadata "
             "nnls_start.u / A_scale), every free multiplier 1, kappa 1, grain A and B from D2's final state "
             "(--start-restart, by name); a00..a24 held at 0 (not part of the bump model)"),
    "chain": ("user 2026-11-10 (fit / fine-tuning cycles): every parameter of the previous cycle's fit of this model "
              "(--start-restart = its fit_restart.npz, by name, exact optimizer coordinates), the grain A and B from "
              "the previous fine tuning (--grain-from), the accumulated fine-tuning term held (--held)")}
BUMP_START_CHOICES = ("nnls", "chain")
HARMONIC_START_CHOICES = ("srccg", "d2", "chain")
PARAMETERIZATION_BUMPS = "bump_slope_multiplier_floor_w10x5_v1"
PARAMETERIZATION_BUMPS_TEST = "bump_slope_multiplier_floor_w10x5_TEST_v1"
GRAIN_TERM_RTOL = 1e-6             # D2's grain_term.txt prints A with 7 significant digits
D2_RESTART = ("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/final/fit_restart.npz")
LOG = "slope_mult:"


def say(text):
    print("%s %s" % (LOG, text), flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


jsonable = pp.jsonable
_pop, _flag, _value, fixed_values = cp._pop, cp._flag, cp._value, cp.fixed_values


# ------------------------------------------------------------------------------------------------ arguments
def rewrite_arguments(argv):
    """Split the slope_mult_fit.py command into (config, grain_fit_fixed.py argv, notes).

    Driver options (before '--'): --point-table PATH (required), --slope-bundle NPZ (required: slope_bins.py bundle;
    the driver adds its 7 columns as --activity-template in bin order), --start-restart R.npz (a new fit: D2's final
    state, by name) | gpu_fit --restart-from (resume of this fit), --floor-ratio F (1e-5, the decision; 0 only with
    --test-floor-zero), --check-full-jacobian, --check-parameters a,b,..., --stop-after-checks (test: stop after the
    start/Jacobian checks, nothing written but the manifest), --eval-log (per-evaluation log, default on).
    Inserted when absent: gpu_fit --template-start 0 (the start comes from --start-restart), grain_fit_fixed
    --fix c_slope_00_40=1; with --transport mpi: --activity host. Refused: see the messages."""
    argv = list(argv)
    if "--" not in argv:
        raise SystemExit("slope_mult: need '--' before the fit_cops arguments")
    split = argv.index("--")
    head, tail = argv[:split], argv[split + 1:]
    notes = []
    cfg = {"point_table": _pop(head, "--point-table"), "slope_bundle": _pop(head, "--slope-bundle"),
           "start_restart": _pop(head, "--start-restart"), "floor_ratio": float(_pop(head, "--floor-ratio",
                                                                                     repr(sm.FLOOR_RATIO))),
           "test_floor_zero": _flag(head, "--test-floor-zero"),
           "check_full_jacobian": _flag(head, "--check-full-jacobian"),
           "stop_after_checks": _flag(head, "--stop-after-checks"),
           "check_parameters": list(CHECK_PARAMETERS), "variant": "D2", "start_choice": _pop(head, "--start"),
           "bump_bundle": _pop(head, "--bump-bundle"), "held": _pop(head, "--held"),
           "grain_from": _pop(head, "--grain-from")}
    for key, opt in (("bump_bundle", "--bump-bundle"), ("held", "--held"), ("grain_from", "--grain-from")):
        if cfg[key] is not None:
            cfg[key] = str(Path(cfg[key]).resolve())
            if not Path(cfg[key]).is_file():
                raise SystemExit("slope_mult: %s not found: %s" % (opt, cfg[key]))
    bumps = cfg["bump_bundle"] is not None
    names = _pop(head, "--check-parameters")
    if names:
        cfg["check_parameters"] = [n.strip() for n in names.split(",") if n.strip()]
    for key, opt in (("point_table", "--point-table"), ("slope_bundle", "--slope-bundle")):
        if cfg[key] is None:
            raise SystemExit("slope_mult: %s is required" % opt)
        cfg[key] = str(Path(cfg[key]).resolve())
        if not Path(cfg[key]).is_file():
            raise SystemExit("slope_mult: %s not found: %s" % (opt, cfg[key]))
    f = cfg["floor_ratio"]
    if f != sm.FLOOR_RATIO:
        if not (f in sm.TEST_FLOOR_RATIOS and cfg["test_floor_zero"]):
            raise SystemExit("slope_mult: --floor-ratio %r is not a decision (user 2026-11-10: %r); 0 is allowed only "
                             "with --test-floor-zero (GPU check (a))" % (f, sm.FLOOR_RATIO))
    elif cfg["test_floor_zero"]:
        raise SystemExit("slope_mult: --test-floor-zero needs --floor-ratio 0")
    if cfg["start_restart"] is not None:
        cfg["start_restart"] = str(Path(cfg["start_restart"]).resolve())
        if not Path(cfg["start_restart"]).is_file():
            raise SystemExit("slope_mult: --start-restart not found: %s" % cfg["start_restart"])
    cfg["restart_from"] = _value(head, "--restart-from")
    if cfg["restart_from"] is not None:
        cfg["restart_from"] = str(Path(cfg["restart_from"]).resolve())
        if cfg["start_restart"] is not None:
            raise SystemExit("slope_mult: --start-restart (new fit) and gpu_fit --restart-from (resume) exclude each other")
    if cfg["start_restart"] is None and cfg["restart_from"] is None:
        raise SystemExit("slope_mult: give --start srccg|d2 --start-restart <D2 fit_restart.npz> (new fit) or gpu_fit "
                         "--restart-from (resume)")
    if cfg["start_choice"] is not None and cfg["start_choice"] not in START_CHOICES:
        raise SystemExit("slope_mult: --start %r: choose one of %s" % (cfg["start_choice"], sorted(START_CHOICES)))
    if cfg["start_restart"] is not None and cfg["start_choice"] is None:
        raise SystemExit("slope_mult: a new fit needs --start srccg (production, user 2026-11-10: the run's harmonic "
                         "coefficients, multipliers 1, kappa 1, D2's grain) or --start d2 (D2's a, kappa, grain); "
                         "there is no default, so an older command line cannot silently change its meaning")
    allowed = BUMP_START_CHOICES if bumps else HARMONIC_START_CHOICES
    if cfg["start_choice"] is not None and cfg["start_choice"] not in allowed:
        raise SystemExit("slope_mult: --start %s is not a start of the %s model (%s)" % (
            cfg["start_choice"], "bump" if bumps else "harmonic", list(allowed)))
    if cfg["start_choice"] == "chain" and cfg["start_restart"] is not None:
        if cfg["held"] is None or cfg["grain_from"] is None:
            raise SystemExit("slope_mult: --start chain needs --held <previous finetune_state.npz> and --grain-from "
                             "<previous finetune_state.npz>")
    elif cfg["grain_from"] is not None and cfg["restart_from"] is None:
        raise SystemExit("slope_mult: --grain-from belongs to --start chain")
    if cfg["held"] is not None and cfg["restart_from"] is None and cfg["start_choice"] != "chain":
        raise SystemExit("slope_mult: --held belongs to --start chain (or a resume of such a fit)")
    ts = _value(head, "--template-start")
    if ts is None:
        head += ["--template-start", "0"]
        notes.append("inserted gpu_fit --template-start 0 (the start comes from --start-restart / --restart-from)")
    elif float(ts) != 0.0:
        raise SystemExit("slope_mult: --template-start must be 0 (a nonzero value replaces the start state)")
    fixed = fixed_values(head)
    key = "c_" + sm.FIXED_BIN
    if key not in fixed:
        head += ["--fix", "%s=%r" % (key, sm.FIXED_VALUE)]
        notes.append("inserted grain_fit_fixed --fix %s=%r (the lead: the 0-40 deg multiplier is held at 1)"
                     % (key, sm.FIXED_VALUE))
    elif fixed[key] != sm.FIXED_VALUE:
        raise SystemExit("slope_mult: %s must be held at %r" % (key, sm.FIXED_VALUE))
    cfg["fixed"] = fixed = fixed_values(head)
    held = tuple(fixed.get(k) for k in ("Tmax", "T75", "Tmin"))
    if held != USER_T0:
        raise SystemExit("slope_mult: the temperatures must be held at Tmax/T75/Tmin = %s with --fix (got %s)"
                         % (USER_T0, held))
    anames = ["a%02d" % j for j in range(25)]
    if bumps:     # the bump model does not use the 25 harmonics: held at 0 (fixed_parameters), inserted when absent
        bad_a = [k for k in anames if k in fixed and fixed[k] != 0.0]
        if bad_a:
            raise SystemExit("slope_mult: bump model: the harmonic coefficients are held at 0; refused %s" % bad_a)
        missing = [k for k in anames if k not in fixed]
        for k in missing:
            head += ["--fix", "%s=0" % k]
        if missing:
            notes.append("inserted grain_fit_fixed --fix a00=0 ... a24=0 (%d; the bump model does not use the "
                         "harmonics: L = A_scale Phi w)" % len(missing))
        cfg["fixed"] = fixed = fixed_values(head)
    other = sorted(k for k in fixed if k not in ("Tmax", "T75", "Tmin", key) and not (bumps and k in anames))
    if other:
        raise SystemExit("slope_mult: only Tmax, T75, Tmin and %s are held in this model; refused --fix %s"
                         % (key, other))
    transport = _value(head, "--transport") or "local"
    activity = _value(head, "--activity")
    if activity == "device":
        raise SystemExit("slope_mult: --activity device cannot form this model's field; use --activity host")
    if transport == "mpi" and activity is None:
        head += ["--activity", "host"]
        notes.append("inserted gpu_fit --activity host (MPI transport: the host field is broadcast)")
    cfg["transport"] = transport
    for bad, why in (("--restart-from", "give --restart-from to gpu_fit (before '--')"),
                     ("--activity-template", "the driver adds the 7 slope-bin columns of --slope-bundle"),
                     ("--fix-temperature", "the temperatures are held with --fix (fixed_parameters), as D2"),
                     ("--ordered-temperature", "the temperatures are held with --fix (fixed_parameters), as D2"),
                     ("--fix-kappa", "kappa is free (the lead's brief)"),
                     ("--check-jacobian", "use grain_fit --check-jacobian or --check-full-jacobian (before '--')"),
                     ("--data-only", "not a fit"),
                     ("--fix-legacy-activity", "the 25 harmonic coefficients are free"),
                     ("--enhancement-scale", "not used by this model")):
        if bad in tail or any(t.startswith(bad + "=") for t in tail):
            raise SystemExit("slope_mult: %s is not allowed after '--': %s" % (bad, why))
    for opt, want in (("--lambda", "0"), ("--degree", "4"), ("--enhancement-degree", "4")):
        v = _value(tail, opt)
        if v is None or float(v) != float(want):
            raise SystemExit("slope_mult: %s %s is required after '--' (got %r)" % (opt, want, v))
    run = next((t for t in tail if not t.startswith("-")), None)
    if run is None:
        raise SystemExit("slope_mult: no run directory after '--'")
    cfg["run"] = str(Path(run).resolve())
    for name in sm.BIN_NAMES:
        tail += ["--activity-template", cfg["slope_bundle"], name]
    notes.append("added --activity-template %s <name> for the 7 bins in order %s" % (cfg["slope_bundle"],
                                                                                   list(sm.BIN_NAMES)))
    if bumps:
        from cops_cuda import figures as fg
        bnames = list((fg.bundle_metadata(cfg["bump_bundle"]).get("slope_bumps") or {}).get("names") or [])
        if not bnames:
            raise SystemExit("slope_mult: %s is not a slope_bumps.py bundle" % cfg["bump_bundle"])
        for name in bnames:
            tail += ["--activity-template", cfg["bump_bundle"], name]
        notes.append("added --activity-template %s <name> for the %d bumps in order (%s .. %s)" % (
            cfg["bump_bundle"], len(bnames), bnames[0], bnames[-1]))
        cfg["bump_names"] = bnames
        cfg["check_parameters"] = (list(cfg["check_parameters"]) if names else
                                   ["c_" + n for n in sm.BIN_NAMES if n != sm.FIXED_BIN] + ["c_" + n for n in bnames] +
                                   ["kappa", "grain_A_uPa"])
    return cfg, head + ["--"] + tail, notes



def model_config(cfg):
    """The GPU model configuration (slope_mult.config_from_bundle_meta) from the bundle's metadata; with
    --bump-bundle the bump basis, with --held the held per-facet term (keys only present in those modes)."""
    from cops_cuda import figures as fg
    meta = fg.bundle_metadata(cfg["slope_bundle"])
    out = sm.config_from_bundle_meta(meta, cfg["floor_ratio"], cfg["slope_bundle"])
    if cfg.get("bump_bundle"):
        bmeta = fg.bundle_metadata(cfg["bump_bundle"])
        out = sm.add_bumps_to_config(out, bmeta, cfg["bump_bundle"])
        meta = dict(meta, slope_bumps=bmeta.get("slope_bumps"))
    if cfg.get("held"):
        out = sm.add_held_to_config(out, cfg["held"], "held_next")
    return out, meta


def parameterization_of(model_cfg):
    if model_cfg.get("basis") == "bumps":
        return PARAMETERIZATION_BUMPS_TEST if model_cfg["test_only"] else PARAMETERIZATION_BUMPS
    return PARAMETERIZATION_TEST if model_cfg["test_only"] else PARAMETERIZATION


# ------------------------------------------------------------------------------------------------ restarts
def check_restart_marker(path, cfg, model_cfg):
    """--restart-from (resume): the restart must be this model's, same floor, bins and point table."""
    meta = cp.restart_meta(path)
    ext = meta.get("model_extension") or {}
    want = parameterization_of(model_cfg)
    problems = []
    for k in ("basis", "bumps_sha256", "n_bumps"):
        if (ext.get("model") or {}).get(k) != model_cfg.get(k):
            problems.append("model %s %r vs %r" % (k, (ext.get("model") or {}).get(k), model_cfg.get(k)))
    if ((ext.get("model") or {}).get("held") or {}).get("sha256") != (model_cfg.get("held") or {}).get("sha256"):
        problems.append("held term differs (restart %r, this command %r)" % (
            ((ext.get("model") or {}).get("held") or {}).get("sha256"), (model_cfg.get("held") or {}).get("sha256")))
    if ext.get("parameterization") != want:
        problems.append("marker %r, this fit %r" % (ext.get("parameterization"), want))
    if (ext.get("model") or {}).get("floor_ratio") != model_cfg["floor_ratio"]:
        problems.append("floor ratio %r vs %r" % ((ext.get("model") or {}).get("floor_ratio"), model_cfg["floor_ratio"]))
    if (ext.get("model") or {}).get("bins_sha256") != model_cfg["bins_sha256"]:
        problems.append("bins differ")
    if (ext.get("point_table") or {}).get("sha256") != sha256(cfg["point_table"]):
        problems.append("point table differs")
    marked = ext.get("start_choice") or ((ext.get("fit_chain") or {}).get("start_choice"))
    if cfg.get("start_choice") is not None and marked != cfg["start_choice"]:
        problems.append("the fit started with --start %r, this command gives --start %r" % (marked, cfg["start_choice"]))
    if problems:
        raise SystemExit("slope_mult: --restart-from %s is not a segment of this fit: %s" % (path, "; ".join(problems)))
    return {"path": str(Path(path).resolve()), "sha256": sha256(path), "parameterization": want,
            "status": meta.get("status"), "segments": meta.get("segments"), "fit_directory": meta.get("fit_directory"),
            "fit_chain": ext.get("fit_chain"), "cost": float((meta.get("optimizer") or {}).get("cost", float("nan")))}


def fit_chain(source, resume, start_choice=None):
    if source is not None and source.get("role") == "start_restart":
        keep = ("path", "sha256", "fit_directory", "cost", "nfev", "status", "segments")
        return {"start": "%s: %s" % (start_choice, START_CHOICES[start_choice]), "start_choice": start_choice,
                "start_restart": {k: source.get(k) for k in keep}}
    if resume is not None:
        prev = dict(resume.get("fit_chain") or {})
        prev["resumed_from"] = list(prev.get("resumed_from") or []) + [
            {k: resume.get(k) for k in ("path", "sha256", "fit_directory", "cost", "status")}]
        return prev
    return {"start": "?"}


def grain_term_values(fit_dir):
    """A_Pa and B_Pa_per_s of a fit's grain_term.txt (None when absent)."""
    path = Path(fit_dir or "/nonexistent") / "grain_term.txt"
    if not path.is_file():
        return None
    out = {"path": str(path), "sha256": sha256(path)}
    for line in path.read_text().splitlines():
        v = line.split()
        if len(v) == 2 and v[0] in ("A_Pa", "B_Pa_per_s"):
            out[v[0]] = float(v[1])
    return out


def start_from_restart(P, el, path, record, choice):
    """Optimizer start, by name from the --start-restart file (D2), for the start choice (START_CHOICES):
    d2:    a00..a24 (exact optimizer coordinates, same ascale), kappa, grain A and B from the restart;
    srccg: a00..a24 = the fitter's prior x0 (the run's srcCG coefficients / ascale), kappa 1, grain A and B from the
           restart (by name; checked against its grain_term.txt to the printed precision);
    both:  every free multiplier START_MULTIPLIER. Mesh, run coefficients, ascale and held temperatures checked."""
    if choice not in START_CHOICES:
        raise SystemExit("slope_mult: unknown start choice %r" % choice)
    from cops_cuda import reference
    from cops_cuda import restart as rs
    state = rs.read(path)
    problems = []
    mesh = reference.load("activity_templates").mesh_identity(el["centroid"], el["area"])
    if mesh != state.mesh_identity():
        problems.append("mesh %s vs the base's %s" % (state.mesh_identity()[:16], mesh[:16]))
    if float(state.meta.get("ascale")) != float(P.ascale):
        problems.append("ascale %r vs %r" % (state.meta.get("ascale"), P.ascale))
    if not np.array_equal(np.asarray(state.coef_prior, float), np.asarray(P.coef, float)):
        problems.append("the run's activity coefficients differ")
    if [float(v) for v in state.T] != list(USER_T0):
        problems.append("the source's temperatures %s are not the held %s" % (list(state.T), USER_T0))
    if problems:
        raise SystemExit("slope_mult: --start-restart %s does not match: %s" % (path, "; ".join(problems)))
    x, report = rs.start_vector(P, state)
    free_m = ["c_" + n for n in sm.BIN_NAMES if n != sm.FIXED_BIN]
    if sorted(report["kept_problem_start"]) != sorted(free_m) or report["clipped"]:
        raise SystemExit("slope_mult: start from %s: names not in the restart %s (expected exactly the free "
                         "multipliers %s), clipped %s" % (path, report["kept_problem_start"], free_m, report["clipped"]))
    names_r = list(state.meta.get("names") or [])
    for k, name in enumerate(P.names):
        if name in free_m:
            x[k] = START_MULTIPLIER
        elif name not in names_r or float(x[k]) != float(state.x[names_r.index(name)]):
            raise SystemExit("slope_mult: start %s = %r is not the source's optimizer coordinate" % (name, x[k]))
    extra = {}
    if choice == "srccg":
        prior = np.asarray(P.x0, float)
        for k, name in enumerate(P.names):
            if name.startswith("a") and name[1:].isdigit():
                x[k] = prior[k]                       # the run's coefficient / ascale (the fitter's own prior)
            elif name == "kappa":
                x[k] = 1.0
        if float(prior[P.names.index("kappa")]) != 1.0:
            raise SystemExit("slope_mult: the fitter's prior kappa is %r, not 1" % prior[P.names.index("kappa")])
        a_s = P.unpack(x)[0]
        coef0 = np.asarray(P.coef[0], float)
        extra["a_vs_run_coefficients_max_rel_diff"] = float(np.max(np.abs(a_s - coef0)) / np.max(np.abs(coef0)))
        extra["a_equals_prior_x0"] = bool(np.array_equal(x[:25], prior[:25]))
        gt = grain_term_values(state.meta.get("fit_directory"))
        extra["source_grain_term"] = gt
        if gt is not None:
            xa = float(x[P.names.index("grain_A_uPa")]) * 1e-6
            xb = float(x[P.names.index("grain_B_uPa_per_100s")]) * (1e-6 / 100.0)
            okA = abs(xa / gt["A_Pa"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("A_Pa") else xa == 0.0
            okB = (abs(xb / gt["B_Pa_per_s"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("B_Pa_per_s") else xb == 0.0)
            extra["grain_term_check"] = {"A_Pa_restart": xa, "B_Pa_per_s_restart": xb, "A_ok": bool(okA),
                                         "B_ok": bool(okB), "rtol": GRAIN_TERM_RTOL}
            if not (okA and okB):
                raise SystemExit("slope_mult: the restart's grain (A %r Pa, B %r Pa/s) disagrees with %s" % (
                    xa, xb, gt["path"]))
    a, b, T, kappa = P.unpack(x)
    held = {"a": bool(np.array_equal(a, np.asarray(state.a, float))),
            "T": [float(v) for v in T] == [float(v) for v in state.T],
            "kappa": float(kappa) == float(np.asarray(state.kappa).reshape(-1)[0]),
            "multipliers": [float(v) for v in b]}
    if choice == "d2" and not (held["a"] and held["T"] and held["kappa"]):
        raise SystemExit("slope_mult: the physical start state differs from the source's: %s" % held)
    if choice == "srccg" and not held["T"]:
        raise SystemExit("slope_mult: the held temperatures differ from the source's: %s" % held)
    from_restart = ([n for n in P.names if n not in free_m] if choice == "d2" else list(GRAIN_NAMES))
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_vector": jsonable(report),
           "start_choice": choice, "start_definition": START_CHOICES[choice],
           "names_from_restart": from_restart, "multipliers_set_to": START_MULTIPLIER,
           "kappa_start": float(kappa), "grain_start": {n: float(x[P.names.index(n)]) for n in GRAIN_NAMES},
           "physical_state": held, "source_names": names_r, "source_x": [float(v) for v in state.x],
           "source_template_names": list(state.template_names),
           "source_template_values": [float(v) for v in state.template_values],
           "x_start": [float(v) for v in x], "prior_x0_unchanged": True,
           "note": "D2's own templates (cliff_ge_60deg, ft_c1..ft_c6) are not part of this model; their activity is "
                   "absent at the start", **extra}
    record["start_restart"] = out
    if choice == "d2":
        say("start d2 from %s: a00..a24, kappa, grain by name (exact optimizer coordinates), multipliers %s = %r; held "
            "T %s; D2's templates %s (values %s) are NOT in this model" % (
                path, free_m, START_MULTIPLIER, list(T), list(state.template_names),
                ["%.4g" % v for v in state.template_values]))
    else:
        say("start srccg: a00..a24 = the run's srcCG coefficients (the fitter's prior x0; a vs coef0 max rel. diff "
            "%.1e), multipliers %s = %r, kappa 1, grain A %.9g uPa, B %.9g uPa/100s by name from %s (grain_term.txt "
            "check %s); held T %s" % (extra["a_vs_run_coefficients_max_rel_diff"], free_m, START_MULTIPLIER,
                                      out["grain_start"]["grain_A_uPa"], out["grain_start"]["grain_B_uPa_per_100s"],
                                      path, extra.get("grain_term_check"), list(T)))
    del state
    return x


def _check_source_state(P, el, state, path):
    from cops_cuda import reference
    problems = []
    mesh = reference.load("activity_templates").mesh_identity(el["centroid"], el["area"])
    if mesh != state.mesh_identity():
        problems.append("mesh %s vs the base's %s" % (state.mesh_identity()[:16], mesh[:16]))
    if float(state.meta.get("ascale")) != float(P.ascale):
        problems.append("ascale %r vs %r" % (state.meta.get("ascale"), P.ascale))
    if not np.array_equal(np.asarray(state.coef_prior, float), np.asarray(P.coef, float)):
        problems.append("the run's activity coefficients differ")
    if [float(v) for v in state.T] != list(USER_T0):
        problems.append("the source's temperatures %s are not the held %s" % (list(state.T), USER_T0))
    if problems:
        raise SystemExit("slope_mult: --start-restart %s does not match: %s" % (path, "; ".join(problems)))


def start_nnls(P, el, path, record, bump_meta):
    """Bump model, --start nnls: bump weights w_k = u_k / A_scale (slope_bumps.py NNLS of L to A_orig), every free
    multiplier 1, kappa 1, grain A and B by name from --start-restart (D2; checked against its grain_term.txt)."""
    from cops_cuda import restart as rs
    state = rs.read(path)
    _check_source_state(P, el, state, path)
    u = np.asarray(bump_meta["nnls_start"]["u"], float)
    bnames = ["c_" + n for n in bump_meta["names"]]
    if len(u) != len(bnames):
        raise SystemExit("slope_mult: %d NNLS weights for %d bumps" % (len(u), len(bnames)))
    names_r = list(state.meta.get("names") or [])
    x = np.array(P.x0, dtype=float, copy=True)
    scale = float(P.enhancement_scale)
    for k, name in enumerate(P.names):
        if name in bnames:
            x[k] = u[bnames.index(name)] / scale
        elif name.startswith("c_slope_"):
            x[k] = START_MULTIPLIER
        elif name == "kappa":
            x[k] = 1.0
        elif name in GRAIN_NAMES:
            x[k] = float(state.x[names_r.index(name)])
        else:
            raise SystemExit("slope_mult: --start nnls: unexpected parameter %s" % name)
    if float(P.x0[P.names.index("kappa")]) != 1.0:
        raise SystemExit("slope_mult: the fitter's prior kappa is not 1")
    gt = grain_term_values(state.meta.get("fit_directory"))
    chk = None
    if gt is not None:
        xa = float(x[P.names.index("grain_A_uPa")]) * 1e-6
        xb = float(x[P.names.index("grain_B_uPa_per_100s")]) * (1e-6 / 100.0)
        okA = abs(xa / gt["A_Pa"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("A_Pa") else xa == 0.0
        okB = abs(xb / gt["B_Pa_per_s"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("B_Pa_per_s") else xb == 0.0
        chk = {"A_Pa_restart": xa, "B_Pa_per_s_restart": xb, "A_ok": bool(okA), "B_ok": bool(okB),
               "rtol": GRAIN_TERM_RTOL}
        if not (okA and okB):
            raise SystemExit("slope_mult: the restart's grain disagrees with %s" % gt["path"])
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_choice": "nnls",
           "start_definition": START_CHOICES["nnls"], "names_from_restart": list(GRAIN_NAMES),
           "multipliers_set_to": START_MULTIPLIER, "kappa_start": 1.0, "A_scale": scale,
           "nnls_u": [float(v) for v in u], "weights_start": [float(v) for v in u / scale],
           "nnls_quality": bump_meta["nnls_start"].get("quality"),
           "grain_start": {n: float(x[P.names.index(n)]) for n in GRAIN_NAMES}, "source_grain_term": gt,
           "grain_term_check": chk, "x_start": [float(v) for v in x], "prior_x0_unchanged": True}
    record["start_restart"] = out
    say("start nnls: %d bump weights from the NNLS start (%d zero; w = u / A_scale, A_scale %.6e), multipliers 1, "
        "kappa 1, grain A %.9g uPa, B %.9g uPa/100s by name from %s (grain_term.txt check %s)" % (
            len(u), int(np.sum(u == 0.0)), scale, out["grain_start"]["grain_A_uPa"],
            out["grain_start"]["grain_B_uPa_per_100s"], path, chk))
    return x


def start_chain(P, el, path, record, grain_from, model_cfg):
    """--start chain (fit / fine-tuning cycles): every parameter by name from the previous cycle's fit of THIS model
    (exact optimizer coordinates), grain A and B from the previous fine tuning (--grain-from)."""
    from cops_cuda import restart as rs
    state = rs.read(path)
    _check_source_state(P, el, state, path)
    ext = state.meta.get("model_extension") or {}
    want = parameterization_of(dict(model_cfg, test_only=model_cfg["test_only"]))
    m_src = ext.get("model") or {}
    problems = []
    if ext.get("parameterization") != want:
        problems.append("marker %r, this model %r" % (ext.get("parameterization"), want))
    for k in ("bins_sha256", "floor_ratio", "basis", "bumps_sha256"):
        if m_src.get(k) != model_cfg.get(k):
            problems.append("%s %r vs %r" % (k, m_src.get(k), model_cfg.get(k)))
    names_r = list(state.meta.get("names") or [])
    if names_r != list(P.names):
        problems.append("parameters %s vs %s" % (names_r, list(P.names)))
    if state.meta.get("status") != "complete":
        problems.append("status %r (the previous fit must be complete)" % state.meta.get("status"))
    if problems:
        raise SystemExit("slope_mult: --start chain from %s: %s" % (path, "; ".join(problems)))
    x = np.asarray(state.x, float).copy()
    with np.load(grain_from, allow_pickle=False) as z:
        gA, gB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        gmeta = json.loads(str(z["meta_json"])) if "meta_json" in z.files else {}
    fit_of_ft = (gmeta.get("fit") or {}).get("fit_directory")
    if fit_of_ft is not None and Path(fit_of_ft).resolve() != Path(state.meta.get("fit_directory")).resolve():
        raise SystemExit("slope_mult: --grain-from %s tuned the fit %s, not %s" % (grain_from, fit_of_ft,
                                                                                   state.meta.get("fit_directory")))
    x[P.names.index("grain_A_uPa")], x[P.names.index("grain_B_uPa_per_100s")] = gA, gB
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_choice": "chain",
           "start_definition": START_CHOICES["chain"], "names_from_restart": [n for n in P.names if n not in GRAIN_NAMES],
           "grain_from": {"path": str(Path(grain_from).resolve()), "sha256": sha256(grain_from),
                          "grain_A_uPa": gA, "grain_B_uPa_per_100s": gB, "fine_tuned_fit": fit_of_ft},
           "held": model_cfg.get("held"), "source_cost": float((state.meta.get("optimizer") or {}).get("cost")),
           "x_start": [float(v) for v in x], "prior_x0_unchanged": True}
    record["start_restart"] = out
    say("start chain: every parameter by name from %s (cost %.9e), grain A %.9g uPa B %.9g uPa/100s from %s, held "
        "term %s (%s nonzero facets)" % (path, out["source_cost"], gA, gB, grain_from,
                                         (model_cfg.get("held") or {}).get("sha256", "-")[:16],
                                         (model_cfg.get("held") or {}).get("nonzero_facets")))
    return x


def floor_share(P, x, el, seen=None):
    """Share of the facets / area (/ seen area) at the floor (s <= F) at x, overall and per bin."""
    st = facet_state(P, x, el)
    area = np.asarray(el["area"], float)
    at = ~st["active"]
    out = {"facets": float(np.mean(at)), "area": float(area[at].sum() / area.sum()),
           "facets_at_floor": int(at.sum()), "facets_total": int(len(at)), "per_bin": {}}
    if seen is not None:
        seen = np.asarray(seen, bool)
        out["seen_area"] = float(area[at & seen].sum() / area[seen].sum())
    for k, name in enumerate(sm.BIN_NAMES):
        sel = st["bins"] == k
        out["per_bin"][name] = {"facets": float(np.mean(at[sel])), "area": float(area[sel & at].sum() / area[sel].sum())}
    return out


def install_start_report(fit_cops, record):
    """Inner wrapper (installed before clip_problem's start and Jacobian checks): after those checks, print the start
    cost with the start choice and the share of facets / area / seen area at the floor at the start."""
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"]:
            state["done"] = True
            P = fit_cops.P
            st = (record.get("clip_core") or {}).get("start_check") or {}
            el = getattr(P, "_slope_el", None)
            seen = P.M.seen_faces() if el is not None else None
            fs = floor_share(P, np.asarray(x0, float), el, seen) if el is not None else None
            record["start_report"] = {"start_choice": (record.get("start_restart") or {}).get("start_choice"),
                                      "cost": (st.get("start") or {}).get("cost"),
                                      "data_cost": (st.get("start") or {}).get("data_cost"), "floor_share": fs}
            ls_ = record.get("l_step") or {}
            record["start_report"]["l_step"] = {"mode": ls_.get("mode"), "threads": ls_.get("threads")}
            say("start check (start %s): cost at the start %.12e (data %.12e); at the floor at the start: %.4f of the "
                "facets (%d of %d), %.4f of the area, %.4f of the seen area; per bin (area) %s; L = Y a of the fitter: "
                "mode %s, threads %s" % (
                    record["start_report"]["start_choice"], st["start"]["cost"], st["start"]["data_cost"],
                    fs["facets"], fs["facets_at_floor"], fs["facets_total"], fs["area"], fs.get("seen_area", float("nan")),
                    {k: round(v["area"], 4) for k, v in fs["per_bin"].items()},
                    ls_.get("mode") or "single call", ls_.get("threads") or 1))
            fit_cops.save_active_manifest()
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


# ------------------------------------------------------------------------------------------------ evaluation log
def install_eval_log(fit_cops, record):
    """Wrap the fun/jac that least_squares passes down: per call kind, seconds and the data cost (fun). The full log
    goes to <out>/evaluations.txt; record["evaluations"] keeps counts and timing statistics."""
    original = fit_cops.least_squares
    stats = {"fun": [], "jac": []}
    state = {"file": None, "n": 0}

    def write(kind, seconds, cost):
        if state["file"] is None and fit_cops.ACTIVE_MANIFEST_PATH:
            path = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent / "evaluations.txt"
            state["file"] = path.open("a")
            state["file"].write("# slope_mult evaluation log: index kind seconds data_cost(fun; 0.5 r.r of the fitter's "
                                "residuals, no crossing rows) utc\n")
        state["n"] += 1
        if state["file"] is not None:
            state["file"].write("%d %s %.6f %s %s\n" % (state["n"], kind, seconds,
                                                        "%.15e" % cost if cost is not None else "-", utc_now()))
            state["file"].flush()

    def least_squares(fun, x0, jac="2-point", **kwargs):
        def fun_logged(x, *a, **k):
            t0 = time.perf_counter()
            r = fun(x, *a, **k)
            dt_ = time.perf_counter() - t0
            c = 0.5 * float(np.dot(r, r))
            stats["fun"].append(dt_)
            write("fun", dt_, c)
            return r
        jac_logged = jac
        if callable(jac):
            def jac_logged(x, *a, **k):
                t0 = time.perf_counter()
                J = jac(x, *a, **k)
                dt_ = time.perf_counter() - t0
                stats["jac"].append(dt_)
                write("jac", dt_, None)
                return J
        try:
            return original(fun_logged, x0, jac=jac_logged, **kwargs)
        finally:
            record["evaluations"] = {k: {"calls": len(v), "seconds_total": float(np.sum(v)) if v else 0.0,
                                         "seconds_median": float(np.median(v)) if v else None,
                                         "seconds_first": float(v[0]) if v else None}
                                     for k, v in stats.items()}
    fit_cops.least_squares = least_squares


def install_stop_after_checks(fit_cops, record):
    """Test option: innermost wrapper; the start check and the Jacobian check (installed on top) run first."""
    original = fit_cops.least_squares

    def least_squares(fun, x0, jac="2-point", **kwargs):
        record["stopped_after_checks"] = True
        fit_cops.save_active_manifest()
        say("stopped after the start and Jacobian checks (--stop-after-checks, test option): no optimizer step, no "
            "products")
        raise SystemExit(0)
    fit_cops.least_squares = least_squares
    return original


# ------------------------------------------------------------------------------------------------ diagnostics
def model_L(model_cfg, a, b, el, basis, scale):
    """L of the model (host): harmonics Y a, or bumps scale * (Phi w) with Phi the contiguous bump block of the
    fitter's template basis (the GPU patch's expression)."""
    nb, nk = sm.n_bins_bumps(model_cfg)
    if model_cfg.get("basis") == "bumps":
        Phi = np.ascontiguousarray(np.asarray(basis, np.float64)[:, nb:nb + nk])
        return sm.bump_L(Phi, np.asarray(b, float)[nb:nb + nk], float(scale)), Phi
    return el["Y"] @ a, None


def facet_state(P, x, el, with_phi=False):
    """Per facet: L, bins, F, g (held), s, A, active (the host expressions of slope_mult.activity)."""
    a, b, T, kappa = P.unpack(x)
    cfg = P.slope_model
    nb, nk = sm.n_bins_bumps(cfg)
    bins = sm.bins_from_basis(np.asarray(P.enhancement_basis)[:, :nb])
    L, Phi = model_L(cfg, a, b, el, P.enhancement_basis, P.enhancement_scale)
    F = cfg["floor_ratio"] * np.asarray(el["AH2O"], dtype=np.float64)
    g, h = sm.held_factor(cfg)
    s, A, active = sm.activity(L, np.asarray(b, float)[:nb], bins, F, g)
    out = {"a": a, "m": np.asarray(b, float)[:nb], "w": np.asarray(b, float)[nb:nb + nk], "T": T, "kappa": kappa,
           "bins": bins, "L": L, "F": F, "s": s, "A": A, "active": active, "g": g, "h": h}
    if with_phi:
        out["Phi"] = Phi
    return out


def install_activity_diagnostics(fit_cops, floor_ratio, model_cfg=None):
    """fit_cops.activity_diagnostics in this model's semantics (same keys; fit_result.txt and the manifest):
    act = A = max(F, m[bin] L); active_area = area with s > F; clipped_area_fraction = area at the floor (s <= F);
    enhancement_relative_rms_D / enhancement_area_weighted_mean = rms / mean of A - A(m = 1) (the activity the
    multipliers add), the rms relative to the rms run map (as the fitter's D normalization)."""
    original = fit_cops.activity_diagnostics

    def activity_diagnostics(a, b, el, seen, coef, enhancement_basis, enhancement_scale):
        if not len(b):
            return original(a, b, el, seen, coef, enhancement_basis, enhancement_scale)
        cfg = model_cfg or {"bin_names": list(sm.BIN_NAMES)}
        nb, nk = sm.n_bins_bumps(cfg)
        bins = sm.bins_from_basis(np.asarray(enhancement_basis)[:, :nb])
        L, _ = model_L(cfg, a, b, el, enhancement_basis, enhancement_scale)
        F = float(floor_ratio) * np.asarray(el["AH2O"], dtype=np.float64)
        g, _ = sm.held_factor(cfg)
        s, act, active = sm.activity(L, np.asarray(b, float)[:nb], bins, F, g)
        _, act1, _ = sm.activity(L, np.ones(nb), bins, F, g)
        act0 = np.maximum(el["Y"] @ coef[0], 0.0)
        area = el["area"]
        tot, tot0 = np.sum(area * act), np.sum(area * act0)
        delta = act - act1
        n0 = np.sum(area * (el["Y"] @ coef[0]) ** 2)
        return {"active_area": np.sum(area[active]) / np.sum(area),
                "activity_ratio": tot / tot0,
                "activity_ratio_seen": np.sum((area * act)[seen]) / np.sum((area * act0)[seen]),
                "unseen_activity_fraction": np.sum((area * act)[~seen]) / tot if tot > 0 else np.nan,
                "unseen_activity_fraction_prior": np.sum((area * act0)[~seen]) / tot0,
                "unseen_area_fraction": np.sum(area[~seen]) / np.sum(area),
                "enhancement_relative_rms_D": np.sqrt(np.sum(area * delta ** 2) / n0),
                "enhancement_area_weighted_mean": np.sum(area * delta) / np.sum(area),
                "clipped_area_fraction": np.sum(area[~active]) / np.sum(area)}
    fit_cops.activity_diagnostics = activity_diagnostics


# ------------------------------------------------------------------------------------------------ products
def write_facets(out, P, r, el, seen, record):
    st = facet_state(P, r["x"], el, with_phi=True)
    path = Path(out) / "slope_mult_facets.npz"
    meta = {"created_utc": utc_now(), "parameterization": P.slope_parameterization, "model": P.slope_model,
            "bin_names": list(sm.BIN_NAMES), "bin_edges_deg": list(sm.BIN_EDGES_DEG),
            "multipliers": [float(v) for v in st["m"]], "kappa": float(st["kappa"]), "T": [float(v) for v in st["T"]],
            "definitions": {"L": "sum_j a_j Y_j (signed linear harmonic field, fit_cops' 25 degree-<=4 Y)",
                            "A": "max(F, m[bin] L): the fitted H2O activity of the model",
                            "floor": "F = floor_ratio * A_orig (elements.bin AH2O)",
                            "bin": "slope bin index 0..6 (bin_names)",
                            "active": "m[bin] L > F (above the floor; at the floor the activity columns are 0)",
                            "seen": "facet contributes to a fit point (fit_cops Model.seen_faces of the fit model)"},
            "a": [float(v) for v in st["a"]], "weights": [float(v) for v in st["w"]],
            "A_scale": float(P.enhancement_scale)}
    if P.slope_model.get("basis") == "bumps":
        meta["definitions"]["L"] = "A_scale * sum_k w_k phi_k (Gaussian bumps of the bump bundle; a held at 0)"
    extra = {}
    if st["h"] is not None:
        extra["held"] = st["h"]
        meta["definitions"]["A"] = "max(F, m[bin] L exp(held)): the fitted H2O activity with the held fine-tuning term"
        meta["definitions"]["held"] = "the held per-facet term h (accumulated fine tuning; A = max(F, m L exp(h)))"
    with path.open("xb") as f:
        np.savez_compressed(f, L=st["L"], A=st["A"], floor=st["F"], bin=st["bins"], seen=np.asarray(seen, bool),
                            active=st["active"], meta_json=np.array(json.dumps(meta)), **extra)
    return path, st


BUMP_TABLE_HEADER = """\
# slope_bumps_table.txt -- written by model_ext/slope_mult_problem.py ({utc})
# {pid}: bump weights of L = A_scale sum_k w_k phi_k, phi_k = exp(-|x - c_k|^2 / (2 s_k^2)); A_scale = {scale:.9e}
# weight = the fitted amplitude w_k (optimizer coordinate, >= 0); u = A_scale w_k (physical activity per unit phi);
# sigma = the fitter's local Gauss-Newton sigma of w_k (fit_result.txt; not a two-sided uncertainty at a bound);
# nearest_area = area of the facets whose nearest centre is c_k; effective_area = sum_e area_e phi_k(e);
# production_share = sum over facets above the floor of area m_b g A_scale w_k phi_k / sum_e area_e A_e (the share of
# the activity x area carried by bump k; the floor and q(SZA) are not attributed)
#bump set facet x_m y_m z_m width_m weight u sigma nearest_area_m2 effective_area_m2 production_share
"""


def write_bump_table(out, P, r, el, st, bump_meta):
    area = np.asarray(el["area"], float)
    Phi, w = st["Phi"], st["w"]
    scale = float(P.enhancement_scale)
    sig = {n: float(v) for n, v in zip(P.names, r["sig"])}
    mrow = st["m"][st["bins"].astype(np.intp)]
    gg = st["g"] if st["g"] is not None else 1.0
    fac = np.where(st["active"], area * mrow * gg * scale, 0.0)
    ptot = float(np.sum(area * st["A"]))
    rows = []
    bumps = bump_meta["bumps"]
    path = Path(out) / "slope_bumps_table.txt"
    with path.open("x") as f:
        f.write(BUMP_TABLE_HEADER.format(utc=utc_now(), pid=P.slope_parameterization, scale=scale))
        for k, bm in enumerate(bumps):
            share = float(np.sum(fac * Phi[:, k]) * w[k] / ptot) if ptot > 0 else float("nan")
            row = {"name": bm["name"], "set": bm["set"], "facet": bm["facet"], "centre_m": bm["centre_m"],
                   "width_m": bm["width_m"], "weight": float(w[k]), "u": scale * float(w[k]),
                   "sigma": sig.get("c_" + bm["name"]), "nearest_centre_area_m2": bm["nearest_centre_area_m2"],
                   "effective_area_m2": bm["effective_area_m2"], "production_share": share}
            rows.append(row)
            f.write("%s %s %d %.3f %.3f %.3f %.3f %.9e %.9e %s %.6e %.6e %.6e\n" % (
                row["name"], row["set"], row["facet"], row["centre_m"][0], row["centre_m"][1], row["centre_m"][2],
                row["width_m"], row["weight"], row["u"], "-" if row["sigma"] is None else "%.3e" % row["sigma"],
                row["nearest_centre_area_m2"], row["effective_area_m2"], row["production_share"]))
        f.write("# sum of the bump production shares %.6f (the rest: facets at the floor)\n" % float(
            np.nansum([r_["production_share"] for r_ in rows])))
    return path, rows


BIN_TABLE_HEADER = """\
# slope_mult_bins.txt -- written by model_ext/slope_mult_problem.py ({utc})
# {pid}: {description}
# floor ratio f = {floor}: F = f * A_orig; at the floor = m_b L <= F (equality included); seen = facet contributes to a
# fit point (fit_cops Model.seen_faces); production share = sum_bin area*A / sum_all area*A (activity x area, i.e. the
# H2O production at equal insolation; the time-dependent q(SZA) is not applied); *_seen: the same over seen facets
# only. sigma: the fitter's local Gauss-Newton sigma (fit_result.txt), not a two-sided uncertainty at a bound.
#bin interval_deg multiplier sigma facets area_fraction seen_area_fraction_of_bin facets_at_floor_frac area_at_floor_frac seen_area_at_floor_frac production_share production_share_seen a_orig_zero_area_frac
"""


def write_bin_table(out, P, r, el, seen, st):
    area = np.asarray(el["area"], float)
    seen = np.asarray(seen, bool)
    tot = float(np.sum(area))
    prod = area * st["A"]
    ptot, pseen = float(np.sum(prod)), float(np.sum(prod[seen]))
    at = ~st["active"]
    zero = st["F"] == 0.0
    sig = {n: float(v) for n, v in zip(P.names, r["sig"])}
    rows = []
    for k, name in enumerate(sm.BIN_NAMES):
        sel = st["bins"] == k
        ss = sel & seen
        a_bin, a_seen = float(area[sel].sum()), float(area[ss].sum())
        row = {"bin": name, "interval_deg": sm.BIN_LABELS[k], "multiplier": float(st["m"][k]),
               "sigma": sig.get("c_" + name), "facets": int(sel.sum()), "area_fraction": a_bin / tot,
               "seen_area_fraction_of_bin": a_seen / a_bin if a_bin > 0 else float("nan"),
               "facets_at_floor_frac": float(np.sum(sel & at)) / max(int(sel.sum()), 1),
               "area_at_floor_frac": float(area[sel & at].sum()) / a_bin if a_bin > 0 else float("nan"),
               "seen_area_at_floor_frac": float(area[ss & at].sum()) / a_seen if a_seen > 0 else float("nan"),
               "production_share": float(prod[sel].sum()) / ptot if ptot > 0 else float("nan"),
               "production_share_seen": float(prod[ss].sum()) / pseen if pseen > 0 else float("nan"),
               "a_orig_zero_area_frac": float(area[sel & zero].sum()) / a_bin if a_bin > 0 else float("nan")}
        rows.append(row)
    total = {"bin": "all", "interval_deg": "[0,180]", "multiplier": float("nan"), "sigma": None,
             "facets": int(len(area)), "area_fraction": 1.0, "seen_area_fraction_of_bin": float(area[seen].sum()) / tot,
             "facets_at_floor_frac": float(np.sum(at)) / len(area), "area_at_floor_frac": float(area[at].sum()) / tot,
             "seen_area_at_floor_frac": float(area[seen & at].sum()) / max(float(area[seen].sum()), 1e-300),
             "production_share": 1.0, "production_share_seen": 1.0,
             "a_orig_zero_area_frac": float(area[zero].sum()) / tot}
    path = Path(out) / "slope_mult_bins.txt"
    with path.open("x") as f:
        f.write(BIN_TABLE_HEADER.format(utc=utc_now(), pid=P.slope_parameterization, description=DESCRIPTION,
                                        floor=P.slope_model["floor_ratio"]))
        for row in rows + [total]:
            f.write("%s %s %s %s %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n" % (
                row["bin"], row["interval_deg"],
                "%.9e" % row["multiplier"] if row["bin"] != "all" else "-",
                ("fixed" if row["bin"] == sm.FIXED_BIN else "-" if row["sigma"] is None else "%.3e" % row["sigma"]),
                row["facets"], row["area_fraction"], row["seen_area_fraction_of_bin"], row["facets_at_floor_frac"],
                row["area_at_floor_frac"], row["seen_area_at_floor_frac"], row["production_share"],
                row["production_share_seen"], row["a_orig_zero_area_frac"]))
        f.write("# total production sum(area*A) = %.9e (seen %.9e); run map sum(area*A_orig) = %.9e\n" % (
            ptot, pseen, float(np.sum(area * np.asarray(el["AH2O"], float)))))
    return path, rows + [total]


TEXT = """\
# SLOPE_MULTIPLIER_MODEL.txt -- written by model_ext/slope_mult_problem.py ({utc})
# {pid}: {description}
# {policy}

H2O activity per facet e (slope bin b(e) of {bundle}):
    L(e) = sum_j a_j Y_j(e) (25 signed degree-<=4 harmonics of fit_cops), s(e) = m_b(e) L(e),
    A(e) = max(F(e), s(e)),  F(e) = {floor} * A_orig(e) (elements.bin AH2O of the base; F = 0 where A_orig = 0).
  Multipliers: {mult}; {fixed} held at 1. Bounds m >= 0, no upper bound.
  Derivative convention: a facet at the floor (s <= F) contributes 0 to the a and m columns (one-sided, as the fitter's
  clip); temperature and kappa columns are the unchanged kernel's with af = A.
  Flux and pressures: the unchanged GPU model with af = A (srcCG night/shadow ratios 0.2 H2O / 0.1 CO2, no flux floor,
  CO2 = kappa x the run's CO2 map, no CO2 floor). Held: {held}.
Base: {run} ({nrows} rows; {ngrid} model-only 1-s grid rows, reason {reason}, use 0). Fit data: {nng} NG, {nrg} RG.
Late window (>= {late} UTC, {nlate} observations): W10x5 of the surface model (mean of five 10-s averages ending at T,
  T-2, ..., T-8 s; T = L2 STOP), grain added at the observation time stamp (D2's treatment, clip_problem.install_core).
Optimizer start: {start}

VALID in this directory:
  fit_result.txt: cost, residuals, optimizer status, the parameter table (a00..a24 = the coefficients of L,
    c_slope_* = the MULTIPLIERS m_b with their local Gauss-Newton sigma, kappa, grain) and the "activity diagnostics"
    lines, computed in THIS model's semantics (slope_mult_problem.install_activity_diagnostics: act = A; active_area =
    area above the floor; clipped_area_fraction = area at the floor; enhancement_* = A - A(m = 1)). Its fixed header
    text ("templates ... one clip after legacy linear + enhancement", "A_scale", "Tikhonov") describes the fitter's
    generic additive layout, NOT this model; deviation_from_prior_D is the additive-model D (lambda = 0: not in the
    cost) and is NOT meaningful here.
  slope_mult_facets.npz (L, A, floor, bin, active, seen per facet), slope_mult_bins.txt (per-bin table),
  SLOPE_MULTIPLIER_MODEL.txt (this file), evaluations.txt (per-evaluation log), grain_term.*,
  start_check_model.npz, fit_manifest.json (model_extension: every argument; slope_multiplier: checks, records).
  COPS_model.dat: observation rows = this model (late observations W10x5(surface) + grain(t_obs), earlier ones
    instantaneous); GRID ROWS ({ngrid}, header line ClipRefitGridRows = D2's grid-row treatment): model-only,
    pNG = INSTANTANEOUS surface model + grain(t), pRG = instantaneous surface model; NOT observations: drop them for
    plots (figures_cluster/scripts/positive_table.py does, with late_window_model.dat); never re-average pNG.
  late_window_model.dat (kind grid | late_obs; same columns as D2's / a positive fit's), pressure_residuals.*,
  pressure_comparison_*.png, figures/ (observation rows only).
  fit_restart.npz: x, names, a, template_values (= the multipliers), kappa, T, grain, predictions: VALID for a resume of
    THIS fit (gpu_fit --restart-from with slope_mult_fit.py; the marker meta model_extension is checked).
NOT VALID (written by the unchanged fitter for its additive clipped-linear layout; do not use):
  activity_linear in fit_restart.npz (= Y a + A_scale B m, not this model's field: use slope_mult_facets.npz),
  ActivityFit.txt as a replay input (its coefficients are those of L, but the C++ applies max(0, linear)),
  ActivityEnhancementFit.txt, ActivityEnhancementFacets.dat (additive A_scale x multiplier x indicator fields),
  activity_facets.dat, activity_maps.png, activity_template_power.dat, activity_degree_power.dat,
  activity_region_diagnostics.dat (all max(0, Y a + A_scale B m)),
  cops_fit.in (LiouvilleCOPS C++ replay: additive templates, clip at 0, no floor, instantaneous model -- it does not
    reproduce this fit), the fine-tuning tools (grain_finetune.py, fine_tune_gpu.py, restart_tool.py facet shares,
    kernels/restart.cu: the original clipped floor-free additive model, instantaneous), any tool that rebuilds the
    activity from fit_restart.npz's activity_linear or template columns as additive templates.
"""


def write_text(out, P, cfg, record):
    lay = P.late_layout
    chain = record.get("fit_chain") or {}
    start = chain.get("start", "?")
    if chain.get("start_restart"):
        s = chain["start_restart"]
        start = ("--start %s = %s; source %s (sha256 %s; cost %.9e, nfev %s); D2's templates are not in this model"
                 % (chain.get("start_choice"), START_CHOICES.get(chain.get("start_choice"), "?"), s["path"], s["sha256"],
                    s["cost"], s["nfev"]))
    if chain.get("resumed_from"):
        start = "resumed from %s (wall-time segments); the fit started with: %s" % (chain["resumed_from"][-1]["path"],
                                                                                   start)
    sr = record.get("start_report") or {}
    if sr.get("cost") is not None and sr.get("floor_share"):
        fs = sr["floor_share"]
        start += ("\n  Cost at the start of this segment %.12e; at the floor at the start: %.4f of the facets, %.4f of the "
                  "area, %.4f of the seen area" % (sr["cost"], fs["facets"], fs["area"], fs.get("seen_area", float("nan"))))
    text = TEXT.format(
        utc=utc_now(), pid=P.slope_parameterization, description=DESCRIPTION, policy=sm.POLICY,
        bundle=cfg["slope_bundle"], floor=P.slope_model["floor_ratio"],
        mult=", ".join("%s %s" % (n, sm.BIN_LABELS[k]) for k, n in enumerate(sm.BIN_NAMES)), fixed=sm.FIXED_BIN,
        held=", ".join("%s=%r" % kv for kv in sorted(cfg["fixed"].items())), run=str(P.M.run),
        nrows=len(lay.ipoint) if lay is not None else "?", ngrid=len(lay.grid_rows) if lay is not None else 0,
        reason=pp.GRID_REASON, nng=int(P.NNG), nrg=int(P.NRG), late=pp.LATE_START_UTC,
        nlate=len(lay.late_rows) if lay is not None else "?", start=start)
    mc = P.slope_model
    variant = []
    if mc.get("basis") == "bumps":
        variant.append(
            "BUMP BASIS (user 2026-11-10, option 3a; supersedes the L line above): L(e) = A_scale sum_k w_k phi_k(x_e),\n"
            "  phi_k = exp(-|x - c_k|^2 / (2 s_k^2)) (%d Gaussian bumps of %s: centres, widths and coverage in its\n"
            "  metadata and in slope_bumps_table.txt), w_k >= 0 = the template amplitudes c_bump_*, A_scale = %.9e;\n"
            "  a00..a24 are held at 0 and do not enter the model (their lines in fit_result.txt are 0; their GPU columns\n"
            "  are exact zeros). Derivatives: dA/dw_k = m_b g A_scale phi_k above the floor, 0 at the floor." % (
                sm.n_bins_bumps(mc)[1], cfg.get("bump_bundle"), float(P.enhancement_scale)))
    if mc.get("held"):
        variant.append(
            "HELD TERM (fine-tuning cycles): A(e) = max(F(e), m_b L(e) exp(h(e))), h = %s[%s] (sha256 %s; %d nonzero\n"
            "  facets, %.4g .. %.4g); held fixed in this fit (every column is multiplied by g = exp(h))." % (
                mc["held"]["path"], mc["held"]["key"], mc["held"]["sha256"], mc["held"]["nonzero_facets"],
                mc["held"]["min"], mc["held"]["max"]))
    if variant:
        variant.append("Fine tuning of this model: model_ext/slope_finetune.py (not the positive/clip tools).")
        text = text.replace("\nBase: ", "\n" + "\n".join(variant) + "\nBase: ", 1)
    path = Path(out) / "SLOPE_MULTIPLIER_MODEL.txt"
    with path.open("x") as f:
        f.write(text)
    return path


def add_table_note(out, P):
    table = Path(out) / "COPS_model.dat"
    lines = table.read_text().splitlines(True)
    k = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    mc = P.slope_model
    note = ("# SlopeMultiplierModel (model_ext/slope_mult_problem.py, %s): pNG/pRG of the observation rows = A = max(%g "
            "A_orig, m_b L%s), L = %s, with the fitted slope-bin multipliers (SLOPE_MULTIPLIER_MODEL.txt); grid rows as "
            "stated in ClipRefitGridRows\n" % (P.slope_parameterization, mc["floor_ratio"],
                                               " exp(h)" if mc.get("held") else "",
                                               "A_scale sum_k w_k phi_k (bumps)" if mc.get("basis") == "bumps" else
                                               "sum_j a_j Y_j"))
    tmp = table.with_name(table.name + ".tmp.%d" % os.getpid())
    tmp.write_text("".join(lines[:k] + [note] + lines[k:]))
    os.replace(str(tmp), str(table))


_MARKER = {}          # the current call's (record, cfg, model_cfg): a second in-process call replaces it (no stale data)


def install_restart_marker(record, cfg, model_cfg):
    from cops_cuda import restart as rs
    _MARKER.update(record=record, cfg=cfg, model_cfg=model_cfg)
    if getattr(rs.from_fit, "_slope_marker", False):
        return
    original = rs.from_fit

    def from_fit(problem, result, fit_dir, manifest, elements):
        state = original(problem, result, fit_dir, manifest, elements)
        record, cfg, model_cfg = _MARKER["record"], _MARKER["cfg"], _MARKER["model_cfg"]
        lw = (record.get("clip_core") or {}).get("late_window") or {}
        bumps = model_cfg.get("basis") == "bumps"
        held = model_cfg.get("held")
        act = ("A = max(floor_ratio * A_orig, m[bin] * L%s), L = %s; template_values = [the 7 multipliers m%s]" % (
            " * exp(held)" if held else "", "A_scale * (Phi w) (bumps; a held at 0)" if bumps else "Y a",
            ", the bump weights w" if bumps else ""))
        state.meta["model_extension"] = {
            "parameterization": problem.slope_parameterization, "description": DESCRIPTION, "policy": sm.POLICY,
            "model": dict(model_cfg), "point_table": (record.get("clip_core") or {}).get("point_table"),
            "slope_bundle": {"path": cfg["slope_bundle"], "sha256": sha256(cfg["slope_bundle"])},
            "bump_bundle": ({"path": cfg["bump_bundle"], "sha256": sha256(cfg["bump_bundle"])}
                            if cfg.get("bump_bundle") else None),
            "fit_chain": record.get("fit_chain"), "fixed": cfg.get("fixed"),
            "start_choice": (record.get("fit_chain") or {}).get("start_choice"),
            "late_window": {k: lw.get(k) for k in ("late_start_utc", "late_observations", "grid_rows",
                                                   "window_installed")},
            "activity": act,
            "warning": "activity_linear in this file is Y a + A_scale B b (the fitter's additive layout), NOT this "
                       "model's activity (slope_mult_facets.npz has A); tools that rebuild the activity from "
                       "activity_linear / template columns, the positive/clip fine tuning and cops_fit.in do not "
                       "evaluate this model (its fine tuning is slope_finetune.py)"}
        return state
    from_fit._slope_marker = True
    rs.from_fit = from_fit


# ------------------------------------------------------------------------------------------------ installation
def install_core(fit_cops, cfg, record, source=None, arrays=None):
    """Before grain_crossing.install: the evaluation log, the test stop, then D2's core (clip_problem.install_core)."""
    core = {}
    record["clip_core"] = core
    record["clip_core_note"] = ("clip_problem.install_core with variant D2 (grid rows, W10x5, start check, Jacobian "
                                "check): its 'variant'/'parameterization'/'description'/'policy' labels name D2's "
                                "late-window treatment, not this model")
    if cfg.get("stop_after_checks"):
        install_stop_after_checks(fit_cops, record)
    install_start_report(fit_cops, record)
    core_cfg = {"variant": "D2", "point_table": cfg["point_table"], "check_full_jacobian": cfg["check_full_jacobian"],
                "check_parameters": cfg["check_parameters"]}
    cp.install_core(fit_cops, core_cfg, core, source, arrays)
    return core


def install_outer(fit_cops, cfg, model_cfg, record, source=None):
    """After fixed_parameters.install: model checks, start, grain-window record, diagnostics, products, marker."""
    base = fit_cops.Problem
    pid = parameterization_of(model_cfg)
    bumps = model_cfg.get("basis") == "bumps"
    nb, nk = sm.n_bins_bumps(model_cfg)
    bump_meta = None
    if bumps:
        from cops_cuda import figures as fg
        bump_meta = fg.bundle_metadata(cfg["bump_bundle"])["slope_bumps"]

    class SlopeOuter(base):
        def __init__(self, M, el, coef, args):
            super().__init__(M, el, coef, args)
            self.slope_parameterization = pid
            self.slope_model = dict(model_cfg)
            names = list(self.names)
            free_m = ["c_" + n for n in sm.BIN_NAMES if n != sm.FIXED_BIN]
            tnames = list(sm.BIN_NAMES) + (list(model_cfg["bump_names"]) if bumps else [])
            if bumps:
                expected = free_m + ["c_" + n for n in model_cfg["bump_names"]] + ["kappa"] + list(GRAIN_NAMES)
            else:
                expected = ["a%02d" % j for j in range(25)] + free_m + ["kappa"] + list(GRAIN_NAMES)
            if names != expected:
                raise SystemExit("slope_mult: parameter layout %s, expected %s" % (names, expected))
            if list(self.template_names) != tnames or self.nh:
                raise SystemExit("slope_mult: templates %s (harmonic modes %d); expected exactly %s"
                                 % (self.template_names, self.nh, tnames))
            basis = np.asarray(self.enhancement_basis)
            bins = sm.bins_from_basis(basis[:, :nb])
            if sm.sha256_bytes(bins) != model_cfg["bins_sha256"] or len(bins) != model_cfg["nfaces"]:
                raise SystemExit("slope_mult: the fitter's template basis is not the configured bin assignment")
            if bumps and sm.sha256_bytes(np.ascontiguousarray(basis[:, nb:nb + nk])) != model_cfg["bumps_sha256"]:
                raise SystemExit("slope_mult: the fitter's bump columns are not the configured bump basis")
            if bumps:
                a0 = self.unpack(self.x0)[0]
                if np.any(a0 != 0.0):
                    raise SystemExit("slope_mult: bump model: a00..a24 must be held at 0 (fixed_parameters)")
            ah2o = np.ascontiguousarray(np.asarray(el["AH2O"], dtype=np.float64))
            if sm.sha256_bytes(ah2o) != model_cfg["ah2o_sha256"]:
                raise SystemExit("slope_mult: AH2O of the base differs from the bundle's A_orig digest")
            a, b, T, kappa = self.unpack(self.x0)
            if [float(v) for v in T] != list(USER_T0):
                raise SystemExit("slope_mult: held temperatures %s, expected %s" % (list(T), USER_T0))
            if float(b[sm.BIN_NAMES.index(sm.FIXED_BIN)]) != sm.FIXED_VALUE:
                raise SystemExit("slope_mult: %s is not held at %r" % (sm.FIXED_BIN, sm.FIXED_VALUE))
            counts = [int(v) for v in np.bincount(bins, minlength=len(sm.BIN_NAMES))]
            record["model_check"] = {"bins_sha256": model_cfg["bins_sha256"], "bin_counts": counts,
                                     "ah2o_sha256": model_cfg["ah2o_sha256"], "layout": names,
                                     "held_T": [float(v) for v in T], "fixed_bin": sm.FIXED_BIN,
                                     "floor_ratio": model_cfg["floor_ratio"], "basis": model_cfg.get("basis",
                                                                                                    "harmonics"),
                                     "bumps_sha256": model_cfg.get("bumps_sha256"), "n_bumps": nk,
                                     "A_scale": float(self.enhancement_scale), "held": model_cfg.get("held")}
            # grain window (clip_problem.ClipOuter's record and D2's check)
            et_u = pp.read_run_points(M.run)
            ip = {i: int(p) for i, p in enumerate(np.asarray(M.ipoint))}
            i_end = int(np.nonzero(self.et_fit == self.et_end)[0][0])
            row_end = int(np.asarray(M.rows)[i_end])
            both = self.useNG & self.useRG
            et_b = self.et_fit[both]
            gw = {"t_end_et": float(self.et_end), "t_end_utc": et_u[ip[row_end]][0],
                  "last_ng_before_gap_et": float(self.et_before), "gap_seconds": float(self.et_end - self.et_before),
                  "pairs_before_gap": int(np.sum(et_b <= self.et_before)),
                  "pairs_after_gap": int(np.sum(et_b >= self.et_end)), "NNG": int(self.NNG), "NRG": int(self.NRG)}
            record["grain_window"] = gw
            lay = self.late_layout
            if lay is None:
                raise SystemExit("slope_mult: no late layout (the point table is required)")
            if row_end != int(lay.late_rows[0]):
                raise SystemExit("slope_mult: the first post-gap NG point (row %d) is not the first late-window "
                                 "observation (row %d, %s): the window start needs a user decision"
                                 % (row_end, int(lay.late_rows[0]), pp.LATE_START_UTC))
            say("model check: bins %s (sha256 %s...), floor ratio %r, A_orig digest ok, held T %s, %s = 1; grain window "
                "t_end %s, gap %.1f s, crossing pairs %d / %d; %d NG, %d RG points" % (
                    counts, model_cfg["bins_sha256"][:16], model_cfg["floor_ratio"], list(T), sm.FIXED_BIN,
                    gw["t_end_utc"], gw["gap_seconds"], gw["pairs_before_gap"], gw["pairs_after_gap"], gw["NNG"],
                    gw["NRG"]))
            self._slope_el = el
            if cfg.get("start_restart"):
                ch = cfg["start_choice"]
                if ch == "nnls":
                    self.restart_start = start_nnls(self, el, cfg["start_restart"], record, bump_meta)
                elif ch == "chain":
                    self.restart_start = start_chain(self, el, cfg["start_restart"], record, cfg["grain_from"],
                                                     model_cfg)
                else:
                    self.restart_start = start_from_restart(self, el, cfg["start_restart"], record, ch)

    fit_cops.Problem = SlopeOuter
    install_activity_diagnostics(fit_cops, model_cfg["floor_ratio"], model_cfg)
    original_write = fit_cops.write_fit

    def write_fit(out, r, *rest):
        record["final_x"] = [float(v) for v in r["x"]]
        diag = original_write(out, r, *rest)
        P = fit_cops.P
        el, seen = rest[6], rest[7]
        cp.write_late_window_products(out, fit_cops, r, record["clip_core"], "D2")
        add_table_note(out, P)
        fpath, st = write_facets(out, P, r, el, seen, record)
        bpath, rows = write_bin_table(out, P, r, el, seen, st)
        paths = [fpath, bpath]
        if bumps:
            kpath, krows = write_bump_table(out, P, r, el, st, bump_meta)
            paths.append(kpath)
            record["bumps_final"] = krows
        del st
        tpath = write_text(out, P, cfg, record)
        paths.append(tpath)
        record.setdefault("products", {}).update({p.name: sha256(p) for p in paths})
        record["bins_final"] = rows
        say("products: %s; multipliers %s" % (", ".join(p.name for p in paths),
                                              ["%.6g" % v["multiplier"] for v in rows[:-1]]))
        return diag
    fit_cops.write_fit = write_fit
    install_restart_marker(record, cfg, model_cfg)
    return SlopeOuter


def install_hooks(grain_crossing, fixed_parameters, cfg, model_cfg, command, sub, source=None, arrays=None):
    """As clip_problem.install_hooks: install_core just before grain_crossing.install, install_outer just after
    fixed_parameters.install; the manifest's model_extension.command becomes ``command``."""
    original_gc = grain_crossing.install
    original_fp = fixed_parameters.install
    main_record = {}

    def gc_install(fit_cops, weight, record):
        main_record["record"] = record
        record["slope_multiplier"] = sub
        install_core(fit_cops, cfg, sub, source, arrays)
        result = original_gc(fit_cops, weight, record)
        install_eval_log(fit_cops, sub)               # outside the crossing wrapper: fun/jac = data residuals
        return result

    def fp_install(fit_cops, fixed, fp_record):
        result = original_fp(fit_cops, fixed, fp_record)
        install_outer(fit_cops, cfg, model_cfg, sub, source)
        record = main_record["record"]
        record["grain_fit_fixed_command"] = record.get("command")
        record["command"] = list(command)
        return result
    grain_crossing.install = gc_install
    fixed_parameters.install = fp_install
    return original_gc, original_fp
