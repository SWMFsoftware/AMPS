# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5c5416d1f7c16a40bcceacd79c696d14e1cb41d7e0a8134def16f62b1c7338bd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Fit-level hooks of the softplus / squared models around the unchanged grain_fit_fixed.py path: D2's core
#   (clip_problem.install_core: grid rows, W10x5, start check, full Jacobian check), model checks (the 7 bin templates +
#   D2's ft_c1..ft_c6 held at 1, no cliff_ge_60deg), the start (--start d2 | lsq | chain; user 2026-11-11: both models
#   from D2), diagnostics, products, restart marker.
# USED BY: nonneg_fit.py (snapshot nonneg_20261111T044920Z).
# 2026-11-13: --slope-mult-min 1 (user 2026-11-12, softplus only): lower bound 1 instead of the fitter's 0 on the six
#   free slope-bin multipliers (start 1, no upper bound), set on the Problem's bounds in the existing hook; recorded in
#   the manifest, restart marker, NONNEG_MODEL.txt, nonneg_bins.txt; a resume / chain start with another value is
#   refused. Snapshot nonneg_20261113T035423Z (production chain 124664/124772).
# ==== end of header ====
"""Fit-level part of the non-negative activity models (softplus / squared harmonics with slope-bin multipliers; user
decisions 2026-11-11, handoff.md section 10); the model and the GPU patch are nonneg_mult.py, the driver nonneg_fit.py.

Local, untracked extension: fit_cops.py, gpu_fit.py, grain_fit.py, grain_crossing.py, grain_fit_fixed.py,
fixed_parameters.py, positive_problem.py and clip_problem.py are NOT edited. nonneg_fit.py runs the unchanged
grain_fit_fixed.py path and installs this module around it, where clip_refit.py installs clip_problem:
* install_core, just before grain_crossing.install: clip_problem.install_core with the D2 configuration (grid rows,
  W10x5 for the observations from 10:33:10 UTC, grain at the time stamp, start check, optional full Jacobian check),
  the per-evaluation log and the test-only stop;
* install_outer, just after fixed_parameters.install: model checks (templates = the 7 bin indicators + D2's ft_c*
  templates, held at 1), the start (--start d2 | lsq | chain), the start report, activity diagnostics in this model's
  semantics, the products and the restart marker.
Everything else is D2's command (jobs/nonneg_fit.sbatch): lategrid base, original T law held (--fix), kappa free, grain +
crossing weight 10, lambda 0, --tbounds 30 400, ftol/xtol/gtol 1e-10, max_nfev 4000; NO cliff_ge_60deg.
"""

import datetime as dt
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

import clip_problem as cp          # D2's core, used unchanged
import positive_problem as pp      # point-table helpers, used unchanged
import nonneg_mult as nm

PARAMETERIZATION = {"softplus": "nonneg_softplus_slope_multiplier_w10x5_v1",
                    "squared": "nonneg_squared_slope_multiplier_w10x5_v1"}
PARAMETERIZATION_TEST = {"softplus": "nonneg_softplus_slope_multiplier_w10x5_TEST_s0_v1"}
DESCRIPTION = {
    "softplus": ("H2O activity A = m_b s ln(1 + exp(L / s)) exp(h), L = Y a + A_s sum_k T_k (D2's ft_c* templates held at "
                 "1), s = 1e-5 A_orig per facet (hard clip max(0, L) where A_orig = 0), one multiplier per "
                 "effective-gravitational-slope bin (m_[0,40) = 1 held); W10x5 for the observations from %s UTC (D2), "
                 "grain at the time stamp" % pp.LATE_START_UTC),
    "squared": ("H2O activity A = m_b (A_s (Y b)^2 + A_s sum_k T_k) exp(h) (a = A_s b; D2's ft_c* templates held at 1, "
                "outside the square), one multiplier per effective-gravitational-slope bin (m_[0,40) = 1 held); W10x5 "
                "for the observations from %s UTC (D2), grain at the time stamp" % pp.LATE_START_UTC)}
USER_T0 = (182.09922, 136.11149, 133.0)
GRAIN_NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
CHECK_PARAMETERS = tuple("c_" + n for n in nm.BIN_NAMES if n != nm.FIXED_BIN) + (
    "a00", "a02", "a07", "a12", "a24", "kappa", "grain_A_uPa")
START_MULTIPLIER = 1.0
START_CHOICES = {
    "d2": ("user 2026-11-11 (softplus): D2's final harmonic coefficients a00..a24 (exact optimizer coordinates), kappa "
           "and grain A, B by name (--start-restart = D2's fit_restart.npz), every free multiplier 1; D2's ft_c* "
           "templates held at 1 as in D2; D2's cliff_ge_60deg is not part of the model"),
    "lsq": ("user 2026-11-11 (squared): b = the area-weighted least-squares fit of A(m = 1) = A_s ((Y b)^2 + sum_k T_k) "
            "to D2's activity max(0, D2 field) over all facets (nonneg_start.py, --start-file), a = A_s b; kappa and "
            "grain A, B by name from D2 (--start-restart); every free multiplier 1"),
    "chain": ("user 2026-11-11 (fit / fine-tuning cycles): every parameter of the previous cycle's fit of this model "
              "(--start-restart = its fit_restart.npz, by name, exact optimizer coordinates), the grain A and B from "
              "the previous fine tuning (--grain-from), the accumulated fine-tuning term held (--held)")}
MODEL_STARTS = {"softplus": ("d2", "chain"), "squared": ("lsq", "chain")}
GRAIN_TERM_RTOL = 1e-6
SLOPE_MULT_MIN_CHOICES = (1.0,)    # user 2026-11-12 (softplus): "do not allow the multiplier factor decrease below 1"
LOG = "nonneg:"
FT_PREFIX = "ft_c"


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


# ------------------------------------------------------------------------------------------------ templates of D2
def ft_templates_from(path):
    """D2's ft_c* templates (name, bundle path + sha256) from its fit_manifest.json (a file or the fit directory), in
    D2's order; every bundle digest checked; cliff_ge_60deg and any other template are left out (user 2026-11-11)."""
    p = Path(path)
    if p.is_dir():
        p = p / "fit_manifest.json"
    m = json.loads(p.read_text())
    out = []
    for t in (m.get("activity_enhancement") or {}).get("activity_templates") or []:
        if not str(t.get("name", "")).startswith(FT_PREFIX):
            continue
        b = t["bundle"]
        if not Path(b["path"]).is_file() or sha256(b["path"]) != b["sha256"]:
            raise SystemExit("nonneg: template bundle %s of %s is missing or changed" % (b["path"], t["name"]))
        out.append({"name": t["name"], "bundle": {"path": str(Path(b["path"]).resolve()), "sha256": b["sha256"]}})
    if not out:
        raise SystemExit("nonneg: %s lists no %s* templates" % (p, FT_PREFIX))
    vals = dict(zip((m.get("parameters") or {}).get("template_names") or [],
                    (m.get("parameters") or {}).get("template_values") or []))
    return out, {"manifest": {"path": str(p.resolve()), "sha256": sha256(p)}, "values_in_source": vals}


# ------------------------------------------------------------------------------------------------ arguments
def rewrite_arguments(argv):
    """Split the nonneg_fit.py command into (config, grain_fit_fixed.py argv, notes).

    Driver options (before '--'): --model softplus|squared (required), --point-table PATH, --slope-bundle NPZ (the 7
    bin columns are added as --activity-template in bin order), --ft-templates-from <D2 fit dir | fit_manifest.json>
    (its ft_c* templates are added in D2's order and held at 1), --start d2|lsq|chain (required for a new fit),
    --start-restart R.npz, --start-file S.npz (lsq), --held F --grain-from F (chain), --s-ratio 0 --test-clip
    (softplus test only: the hard clip everywhere), --slope-mult-min 1 (softplus, user 2026-11-12: lower bound 1 of the
    six free multipliers instead of the fitter's 0; absent: unchanged), --check-full-jacobian, --check-parameters
    a,b,.., --stop-after-checks. Inserted when absent: gpu_fit --template-start 0, --fix c_slope_00_40=1, --fix c_<ft>=1."""
    argv = list(argv)
    if "--" not in argv:
        raise SystemExit("nonneg: need '--' before the fit_cops arguments")
    split = argv.index("--")
    head, tail = argv[:split], argv[split + 1:]
    notes = []
    cfg = {"model": _pop(head, "--model"), "point_table": _pop(head, "--point-table"),
           "slope_bundle": _pop(head, "--slope-bundle"), "ft_templates_from": _pop(head, "--ft-templates-from"),
           "start_restart": _pop(head, "--start-restart"), "start_file": _pop(head, "--start-file"),
           "start_choice": _pop(head, "--start"), "held": _pop(head, "--held"), "grain_from": _pop(head, "--grain-from"),
           "s_ratio": _pop(head, "--s-ratio"), "test_clip": _flag(head, "--test-clip"),
           "slope_mult_min": _pop(head, "--slope-mult-min"),
           "check_full_jacobian": _flag(head, "--check-full-jacobian"),
           "stop_after_checks": _flag(head, "--stop-after-checks"), "check_parameters": list(CHECK_PARAMETERS),
           "variant": "D2"}
    if cfg["model"] not in nm.MODELS:
        raise SystemExit("nonneg: --model %r: choose one of %s" % (cfg["model"], nm.MODELS))
    if cfg["slope_mult_min"] is None:
        del cfg["slope_mult_min"]               # absent: today's fit, its records unchanged
    else:
        v = float(cfg["slope_mult_min"])
        if cfg["model"] != "softplus":
            raise SystemExit("nonneg: --slope-mult-min is a decision for the softplus model only (user 2026-11-12)")
        if v not in SLOPE_MULT_MIN_CHOICES:
            raise SystemExit("nonneg: --slope-mult-min %r is not a decision (user 2026-11-12: %s)" % (
                v, SLOPE_MULT_MIN_CHOICES))
        cfg["slope_mult_min"] = v
    names = _pop(head, "--check-parameters")
    if names:
        cfg["check_parameters"] = [n.strip() for n in names.split(",") if n.strip()]
    for key, opt in (("point_table", "--point-table"), ("slope_bundle", "--slope-bundle"),
                     ("ft_templates_from", "--ft-templates-from")):
        if cfg[key] is None:
            raise SystemExit("nonneg: %s is required" % opt)
        cfg[key] = str(Path(cfg[key]).resolve())
        if not Path(cfg[key]).exists():
            raise SystemExit("nonneg: %s not found: %s" % (opt, cfg[key]))
    for key, opt in (("start_restart", "--start-restart"), ("start_file", "--start-file"), ("held", "--held"),
                     ("grain_from", "--grain-from")):
        if cfg[key] is not None:
            cfg[key] = str(Path(cfg[key]).resolve())
            if not Path(cfg[key]).is_file():
                raise SystemExit("nonneg: %s not found: %s" % (opt, cfg[key]))
    if cfg["s_ratio"] is None:
        cfg["s_ratio"] = nm.S_RATIO
    else:
        cfg["s_ratio"] = float(cfg["s_ratio"])
        if cfg["model"] != "softplus":
            raise SystemExit("nonneg: --s-ratio belongs to the softplus model")
    if cfg["model"] == "softplus" and cfg["s_ratio"] != nm.S_RATIO:
        if not (cfg["s_ratio"] in nm.TEST_S_RATIOS and cfg["test_clip"]):
            raise SystemExit("nonneg: --s-ratio %r is not the decision %r; 0 only with --test-clip" % (
                cfg["s_ratio"], nm.S_RATIO))
    elif cfg["test_clip"]:
        raise SystemExit("nonneg: --test-clip needs --model softplus --s-ratio 0")
    cfg["restart_from"] = _value(head, "--restart-from")
    if cfg["restart_from"] is not None:
        cfg["restart_from"] = str(Path(cfg["restart_from"]).resolve())
        if cfg["start_restart"] is not None:
            raise SystemExit("nonneg: --start-restart (new fit) and gpu_fit --restart-from (resume) exclude each other")
    if cfg["start_restart"] is None and cfg["restart_from"] is None:
        raise SystemExit("nonneg: give --start ... --start-restart <R> (new fit) or gpu_fit --restart-from (resume)")
    if cfg["start_choice"] is not None and cfg["start_choice"] not in MODEL_STARTS[cfg["model"]]:
        raise SystemExit("nonneg: --start %r is not a start of the %s model %s" % (
            cfg["start_choice"], cfg["model"], MODEL_STARTS[cfg["model"]]))
    if cfg["start_restart"] is not None and cfg["start_choice"] is None:
        raise SystemExit("nonneg: a new fit needs --start (%s)" % ", ".join(MODEL_STARTS[cfg["model"]]))
    if cfg["start_choice"] == "lsq" and cfg["start_file"] is None:
        raise SystemExit("nonneg: --start lsq needs --start-file <nonneg_start.py output>")
    if cfg["start_file"] is not None and cfg["start_choice"] != "lsq" and cfg["restart_from"] is None:
        raise SystemExit("nonneg: --start-file belongs to --start lsq")
    if cfg["start_choice"] == "chain" and cfg["start_restart"] is not None:
        if cfg["held"] is None or cfg["grain_from"] is None:
            raise SystemExit("nonneg: --start chain needs --held and --grain-from (the previous finetune_state.npz)")
    elif cfg["grain_from"] is not None and cfg["restart_from"] is None:
        raise SystemExit("nonneg: --grain-from belongs to --start chain")
    if cfg["held"] is not None and cfg["restart_from"] is None and cfg["start_choice"] != "chain":
        raise SystemExit("nonneg: --held belongs to --start chain (or a resume of such a fit)")
    ts = _value(head, "--template-start")
    if ts is None:
        head += ["--template-start", "0"]
        notes.append("inserted gpu_fit --template-start 0 (the start comes from --start-restart / --restart-from)")
    elif float(ts) != 0.0:
        raise SystemExit("nonneg: --template-start must be 0")
    fts, ft_src = ft_templates_from(cfg["ft_templates_from"])
    cfg["ft_templates"], cfg["ft_source"] = fts, ft_src
    fixed = fixed_values(head)
    want_fixed = {"c_" + nm.FIXED_BIN: nm.FIXED_VALUE}
    want_fixed.update({"c_" + t["name"]: 1.0 for t in fts})
    for key, v in want_fixed.items():
        if key not in fixed:
            head += ["--fix", "%s=%r" % (key, v)]
            notes.append("inserted grain_fit_fixed --fix %s=%r" % (key, v))
        elif fixed[key] != v:
            raise SystemExit("nonneg: %s must be held at %r (user 2026-11-11)" % (key, v))
    cfg["fixed"] = fixed = fixed_values(head)
    if tuple(fixed.get(k) for k in ("Tmax", "T75", "Tmin")) != USER_T0:
        raise SystemExit("nonneg: the temperatures must be held at %s with --fix" % (USER_T0,))
    other = sorted(k for k in fixed if k not in ("Tmax", "T75", "Tmin") and k not in want_fixed)
    if other:
        raise SystemExit("nonneg: only Tmax, T75, Tmin, c_%s and the ft templates are held; refused --fix %s" % (
            nm.FIXED_BIN, other))
    transport = _value(head, "--transport") or "local"
    if transport != "local":
        raise SystemExit("nonneg: the column factors go through local shared memory: --transport local only")
    if _value(head, "--activity") == "device":
        raise SystemExit("nonneg: --activity device cannot form this model's field")
    cfg["transport"] = transport
    for bad, why in (("--restart-from", "give --restart-from to gpu_fit (before '--')"),
                     ("--activity-template", "the driver adds the bin and ft template columns"),
                     ("--fix-temperature", "the temperatures are held with --fix"),
                     ("--ordered-temperature", "the temperatures are held with --fix"),
                     ("--fix-kappa", "kappa is free"), ("--check-jacobian", "use --check-full-jacobian"),
                     ("--data-only", "not a fit"), ("--fix-legacy-activity", "the 25 harmonic coefficients are free"),
                     ("--enhancement-scale", "A_s is the fitter's own scale")):
        if bad in tail or any(t.startswith(bad + "=") for t in tail):
            raise SystemExit("nonneg: %s is not allowed after '--': %s" % (bad, why))
    for opt, want in (("--lambda", "0"), ("--degree", "4"), ("--enhancement-degree", "4")):
        v = _value(tail, opt)
        if v is None or float(v) != float(want):
            raise SystemExit("nonneg: %s %s is required after '--' (got %r)" % (opt, want, v))
    run = next((t for t in tail if not t.startswith("-")), None)
    if run is None:
        raise SystemExit("nonneg: no run directory after '--'")
    cfg["run"] = str(Path(run).resolve())
    for name in nm.BIN_NAMES:
        tail += ["--activity-template", cfg["slope_bundle"], name]
    for t in fts:
        tail += ["--activity-template", t["bundle"]["path"], t["name"]]
    notes.append("added --activity-template for the 7 bins (%s) and the %d ft templates %s (held at 1)" % (
        cfg["slope_bundle"], len(fts), [t["name"] for t in fts]))
    return cfg, head + ["--"] + tail, notes


def model_config(cfg):
    """The constant process configuration (nonneg_mult.process_config): no held term."""
    from cops_cuda import figures as fg
    meta = fg.bundle_metadata(cfg["slope_bundle"])
    return nm.process_config(cfg["model"], meta, cfg["ft_templates"], cfg["s_ratio"], cfg["slope_bundle"]), meta


def held_of(cfg, model_cfg):
    """The stage's held record ({path, key, sha256, ...} or None) from --held."""
    if not cfg.get("held"):
        return None
    h = nm.load_held(cfg["held"], "held_next", model_cfg["nfaces"])
    return nm.held_record(cfg["held"], "held_next", h)


def parameterization_of(model_cfg):
    if model_cfg.get("test_only"):
        return PARAMETERIZATION_TEST[model_cfg["model"]]
    return PARAMETERIZATION[model_cfg["model"]]


# ------------------------------------------------------------------------------------------------ restarts
def check_restart_marker(path, cfg, model_cfg, held):
    meta = cp.restart_meta(path)
    ext = meta.get("model_extension") or {}
    want = parameterization_of(model_cfg)
    problems = []
    if ext.get("parameterization") != want:
        problems.append("marker %r, this fit %r" % (ext.get("parameterization"), want))
    src = ext.get("model") or {}
    for k in ("model", "s_ratio", "bins_sha256", "templates", "ah2o_sha256"):
        if src.get(k) != model_cfg.get(k):
            problems.append("model %s %r vs %r" % (k, src.get(k), model_cfg.get(k)))
    if (ext.get("held") or {}).get("sha256") != (held or {}).get("sha256"):
        problems.append("held term differs (restart %r, this command %r)" % ((ext.get("held") or {}).get("sha256"),
                                                                           (held or {}).get("sha256")))
    if (ext.get("point_table") or {}).get("sha256") != sha256(cfg["point_table"]):
        problems.append("point table differs")
    marked = ext.get("start_choice")
    if cfg.get("start_choice") is not None and marked != cfg["start_choice"]:
        problems.append("the fit started with --start %r, this command gives %r" % (marked, cfg["start_choice"]))
    if ext.get("slope_mult_min") != cfg.get("slope_mult_min"):
        problems.append("multiplier lower bound (--slope-mult-min) %r in the restart, %r in this command" % (
            ext.get("slope_mult_min"), cfg.get("slope_mult_min")))
    if problems:
        raise SystemExit("nonneg: --restart-from %s is not a segment of this fit: %s" % (path, "; ".join(problems)))
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
    path = Path(fit_dir or "/nonexistent") / "grain_term.txt"
    if not path.is_file():
        return None
    out = {"path": str(path), "sha256": sha256(path)}
    for line in path.read_text().splitlines():
        v = line.split()
        if len(v) == 2 and v[0] in ("A_Pa", "B_Pa_per_s"):
            out[v[0]] = float(v[1])
    return out


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
        raise SystemExit("nonneg: --start-restart %s does not match: %s" % (path, "; ".join(problems)))


def _grain_check(x, P, state):
    gt = grain_term_values(state.meta.get("fit_directory"))
    if gt is None:
        return None, None
    xa = float(x[P.names.index("grain_A_uPa")]) * 1e-6
    xb = float(x[P.names.index("grain_B_uPa_per_100s")]) * (1e-6 / 100.0)
    okA = abs(xa / gt["A_Pa"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("A_Pa") else xa == 0.0
    okB = abs(xb / gt["B_Pa_per_s"] - 1.0) <= GRAIN_TERM_RTOL if gt.get("B_Pa_per_s") else xb == 0.0
    chk = {"A_Pa_restart": xa, "B_Pa_per_s_restart": xb, "A_ok": bool(okA), "B_ok": bool(okB),
           "rtol": GRAIN_TERM_RTOL}
    if not (okA and okB):
        raise SystemExit("nonneg: the restart's grain (A %r Pa, B %r Pa/s) disagrees with %s" % (xa, xb, gt["path"]))
    return gt, chk


def start_d2(P, el, path, record, cfg):
    """--start d2 (softplus): a00..a24, kappa, grain by name from D2 (exact optimizer coordinates), multipliers 1."""
    from cops_cuda import restart as rs
    state = rs.read(path)
    _check_source_state(P, el, state, path)
    x, report = rs.start_vector(P, state)
    free_m = ["c_" + n for n in nm.BIN_NAMES if n != nm.FIXED_BIN]
    if sorted(report["kept_problem_start"]) != sorted(free_m) or report["clipped"]:
        raise SystemExit("nonneg: start from %s: names not in the restart %s (expected the free multipliers %s), "
                         "clipped %s" % (path, report["kept_problem_start"], free_m, report["clipped"]))
    names_r = list(state.meta.get("names") or [])
    for k, name in enumerate(P.names):
        if name in free_m:
            x[k] = START_MULTIPLIER
        elif name not in names_r or float(x[k]) != float(state.x[names_r.index(name)]):
            raise SystemExit("nonneg: start %s = %r is not the source's optimizer coordinate" % (name, x[k]))
    a, b, T, kappa = P.unpack(x)
    nb = len(nm.BIN_NAMES)
    src_t = dict(zip(state.template_names, [float(v) for v in state.template_values]))
    ft_src = {t["name"]: src_t.get(t["name"]) for t in cfg["ft_templates"]}
    if any(v != 1.0 for v in ft_src.values()):
        raise SystemExit("nonneg: D2's ft template values %s are not 1 (the model holds them at 1)" % ft_src)
    held = {"a": bool(np.array_equal(a, np.asarray(state.a, float))),
            "T": [float(v) for v in T] == [float(v) for v in state.T],
            "kappa": float(kappa) == float(np.asarray(state.kappa).reshape(-1)[0]),
            "multipliers": [float(v) for v in b[:nb]], "ft_values": [float(v) for v in b[nb:]]}
    if not (held["a"] and held["T"] and held["kappa"]):
        raise SystemExit("nonneg: the physical start state differs from the source's: %s" % held)
    gt, chk = _grain_check(x, P, state)
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_vector": jsonable(report),
           "start_choice": "d2", "start_definition": START_CHOICES["d2"],
           "names_from_restart": [n for n in P.names if n not in free_m], "multipliers_set_to": START_MULTIPLIER,
           "kappa_start": float(kappa), "grain_start": {n: float(x[P.names.index(n)]) for n in GRAIN_NAMES},
           "physical_state": held, "source_template_names": list(state.template_names),
           "source_template_values": [float(v) for v in state.template_values], "ft_values_in_source": ft_src,
           "source_grain_term": gt, "grain_term_check": chk, "x_start": [float(v) for v in x]}
    record["start_restart"] = out
    say("start d2 from %s: a00..a24, kappa %.9g, grain A %.9g uPa B %.9g uPa/100s by name (exact optimizer "
        "coordinates), multipliers %s = 1; D2's ft templates %s held at 1; D2's cliff_ge_60deg (value %s) is not in "
        "this model" % (path, float(kappa), out["grain_start"]["grain_A_uPa"], out["grain_start"]["grain_B_uPa_per_100s"],
                        free_m, list(ft_src), src_t.get("cliff_ge_60deg")))
    return x


def start_lsq(P, el, path, record, cfg):
    """--start lsq (squared): a = A_s b from --start-file (nonneg_start.py), kappa and grain by name from D2,
    multipliers 1."""
    from cops_cuda import restart as rs
    state = rs.read(path)
    _check_source_state(P, el, state, path)
    with np.load(cfg["start_file"], allow_pickle=False) as z:
        b_start = np.asarray(z["b"], float)
        smeta = json.loads(str(z["meta_json"]))
    scale = float(P.enhancement_scale)
    problems = []
    if smeta.get("d2_restart", {}).get("sha256") != sha256(path):
        problems.append("the start file was computed from %s, not %s" % (smeta.get("d2_restart"), path))
    if float(smeta.get("A_scale")) != scale:
        problems.append("A_scale %r vs the fitter's %r" % (smeta.get("A_scale"), scale))
    if smeta.get("mesh_identity") != state.mesh_identity():
        problems.append("mesh differs")
    if [t["name"] for t in smeta.get("ft_templates", [])] != [t["name"] for t in cfg["ft_templates"]]:
        problems.append("ft templates %s vs %s" % ([t["name"] for t in smeta.get("ft_templates", [])],
                                                   [t["name"] for t in cfg["ft_templates"]]))
    if len(b_start) != 25:
        problems.append("%d coefficients" % len(b_start))
    if problems:
        raise SystemExit("nonneg: --start-file %s: %s" % (cfg["start_file"], "; ".join(problems)))
    names_r = list(state.meta.get("names") or [])
    x = np.array(P.x0, dtype=float, copy=True)
    for k, name in enumerate(P.names):
        if name.startswith("a") and name[1:].isdigit():
            x[k] = scale * b_start[int(name[1:])] / P.ascale
        elif name.startswith("c_slope_"):
            x[k] = START_MULTIPLIER
        elif name == "kappa" or name in GRAIN_NAMES:
            x[k] = float(state.x[names_r.index(name)])
        else:
            raise SystemExit("nonneg: --start lsq: unexpected parameter %s" % name)
    gt, chk = _grain_check(x, P, state)
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_choice": "lsq",
           "start_definition": START_CHOICES["lsq"], "start_file": {"path": cfg["start_file"],
                                                                 "sha256": sha256(cfg["start_file"])},
           "b": [float(v) for v in b_start], "A_scale": scale, "start_quality": smeta.get("quality"),
           "names_from_restart": ["kappa"] + list(GRAIN_NAMES), "multipliers_set_to": START_MULTIPLIER,
           "kappa_start": float(x[P.names.index("kappa")]),
           "grain_start": {n: float(x[P.names.index(n)]) for n in GRAIN_NAMES}, "source_grain_term": gt,
           "grain_term_check": chk, "x_start": [float(v) for v in x]}
    record["start_restart"] = out
    say("start lsq: a = A_s b from %s (area-weighted LSQ to D2's activity; quality %s), kappa %.9g and grain A %.9g "
        "uPa B %.9g uPa/100s by name from %s, multipliers 1" % (
            cfg["start_file"], {k: smeta.get("quality", {}).get(k) for k in ("relative_rms_residual", "iterations",
                                                                             "converged")},
            out["kappa_start"], out["grain_start"]["grain_A_uPa"], out["grain_start"]["grain_B_uPa_per_100s"], path))
    return x


def start_chain(P, el, path, record, grain_from, model_cfg, held, slope_mult_min=None):
    from cops_cuda import restart as rs
    state = rs.read(path)
    _check_source_state(P, el, state, path)
    ext = state.meta.get("model_extension") or {}
    m_src = ext.get("model") or {}
    problems = []
    if ext.get("parameterization") != parameterization_of(model_cfg):
        problems.append("marker %r, this model %r" % (ext.get("parameterization"), parameterization_of(model_cfg)))
    for k in ("model", "s_ratio", "bins_sha256", "templates"):
        if m_src.get(k) != model_cfg.get(k):
            problems.append("%s %r vs %r" % (k, m_src.get(k), model_cfg.get(k)))
    names_r = list(state.meta.get("names") or [])
    if names_r != list(P.names):
        problems.append("parameters %s vs %s" % (names_r, list(P.names)))
    if state.meta.get("status") != "complete":
        problems.append("status %r (the previous fit must be complete)" % state.meta.get("status"))
    if ext.get("slope_mult_min") != slope_mult_min:
        problems.append("multiplier lower bound (--slope-mult-min) %r in the previous fit, %r in this command" % (
            ext.get("slope_mult_min"), slope_mult_min))
    if problems:
        raise SystemExit("nonneg: --start chain from %s: %s" % (path, "; ".join(problems)))
    x = np.asarray(state.x, float).copy()
    with np.load(grain_from, allow_pickle=False) as z:
        gA, gB = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
        gmeta = json.loads(str(z["meta_json"])) if "meta_json" in z.files else {}
    fit_of_ft = (gmeta.get("fit") or {}).get("fit_directory")
    if fit_of_ft is not None and Path(fit_of_ft).resolve() != Path(state.meta.get("fit_directory")).resolve():
        raise SystemExit("nonneg: --grain-from %s tuned the fit %s, not %s" % (grain_from, fit_of_ft,
                                                                              state.meta.get("fit_directory")))
    x[P.names.index("grain_A_uPa")], x[P.names.index("grain_B_uPa_per_100s")] = gA, gB
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_choice": "chain",
           "start_definition": START_CHOICES["chain"], "names_from_restart": [n for n in P.names if n not in GRAIN_NAMES],
           "grain_from": {"path": str(Path(grain_from).resolve()), "sha256": sha256(grain_from),
                          "grain_A_uPa": gA, "grain_B_uPa_per_100s": gB, "fine_tuned_fit": fit_of_ft},
           "held": held, "source_cost": float((state.meta.get("optimizer") or {}).get("cost")),
           "x_start": [float(v) for v in x]}
    record["start_restart"] = out
    say("start chain: every parameter by name from %s (cost %.9e), grain A %.9g uPa B %.9g uPa/100s from %s, held "
        "term %s (%s nonzero facets)" % (path, out["source_cost"], gA, gB, grain_from,
                                         (held or {}).get("sha256", "-")[:16], (held or {}).get("nonzero_facets")))
    return x


# ------------------------------------------------------------------------------------------------ host state, reports
def facet_state(P, x, el):
    """Per facet: u = Y a, TS = A_s sum c T, q-related fields, A, bins, m (nonneg_mult.host_fields, this stage's held
    term; no override)."""
    a, b, T, kappa = P.unpack(x)
    cfg = P.nonneg_model
    nb, nt = nm.n_bins_templates(cfg)
    bins = nm.sm.bins_from_basis(np.asarray(P.enhancement_basis)[:, :nb])
    Bft = np.ascontiguousarray(np.asarray(P.enhancement_basis, np.float64)[:, nb:nb + nt])
    f = nm.host_fields(cfg["model"], el["Y"], np.asarray(a, float), Bft, np.asarray(b, float)[nb:nb + nt],
                       float(P.enhancement_scale), bins, np.asarray(b, float)[:nb], nm.width(cfg), nm.held_g(), True)
    f.update(a=a, m=np.asarray(b, float)[:nb], c=np.asarray(b, float)[nb:nb + nt], T=T, kappa=kappa, bins=bins,
             h=nm.held_h())
    return f


def zero_share(P, x, el, seen=None):
    st = facet_state(P, x, el)
    area = np.asarray(el["area"], float)
    z = st["A"] <= 0.0
    out = {"facets": float(np.mean(z)), "area": float(area[z].sum() / area.sum()), "facets_zero": int(z.sum()),
           "facets_total": int(len(z)), "per_bin": {}}
    if seen is not None:
        seen = np.asarray(seen, bool)
        out["seen_area"] = float(area[z & seen].sum() / area[seen].sum())
    for k, name in enumerate(nm.BIN_NAMES):
        sel = st["bins"] == k
        out["per_bin"][name] = {"facets": float(np.mean(z[sel])), "area": float(area[sel & z].sum() / area[sel].sum())}
    if P.nonneg_model["model"] == "softplus":
        s = nm.width(P.nonneg_model)
        out["hard_clip_facets"] = int(np.sum(s == 0.0))
        out["hard_clip_area"] = float(area[s == 0.0].sum() / area.sum())
    return out


def install_start_report(fit_cops, record):
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"]:
            state["done"] = True
            P = fit_cops.P
            st = (record.get("clip_core") or {}).get("start_check") or {}
            el = getattr(P, "_nonneg_el", None)
            seen = P.M.seen_faces() if el is not None else None
            zs = zero_share(P, np.asarray(x0, float), el, seen) if el is not None else None
            record["start_report"] = {"start_choice": (record.get("start_restart") or {}).get("start_choice"),
                                      "cost": (st.get("start") or {}).get("cost"),
                                      "data_cost": (st.get("start") or {}).get("data_cost"), "zero_share": zs,
                                      "threads": nm.thread_info()}
            say("start check (start %s): cost at the start %.12e (data %.12e); A = 0 at the start: %.4f of the facets "
                "(%d of %d), %.4f of the area, %.4f of the seen area; per bin (area) %s; host threads %s" % (
                    record["start_report"]["start_choice"], st["start"]["cost"], st["start"]["data_cost"],
                    zs["facets"], zs["facets_zero"], zs["facets_total"], zs["area"], zs.get("seen_area", float("nan")),
                    {k: round(v["area"], 4) for k, v in zs["per_bin"].items()}, nm.thread_info()))
            fit_cops.save_active_manifest()
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


def install_eval_log(fit_cops, record):
    original = fit_cops.least_squares
    stats = {"fun": [], "jac": []}
    state = {"file": None, "n": 0}

    def write(kind, seconds, cost):
        if state["file"] is None and fit_cops.ACTIVE_MANIFEST_PATH:
            path = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent / "evaluations.txt"
            state["file"] = path.open("a")
            state["file"].write("# nonneg evaluation log: index kind seconds data_cost(fun; 0.5 r.r of the fitter's "
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
            stats["fun"].append(dt_)
            write("fun", dt_, 0.5 * float(np.dot(r, r)))
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
            if state["file"] is not None:
                state["file"].close()
                state["file"] = None
            record["evaluations"] = {k: {"calls": len(v), "seconds_total": float(np.sum(v)) if v else 0.0,
                                         "seconds_median": float(np.median(v)) if v else None,
                                         "seconds_first": float(v[0]) if v else None}
                                     for k, v in stats.items()}
    fit_cops.least_squares = least_squares


def install_stop_after_checks(fit_cops, record):
    original = fit_cops.least_squares

    def least_squares(fun, x0, jac="2-point", **kwargs):
        record["stopped_after_checks"] = True
        fit_cops.save_active_manifest()
        say("stopped after the start and Jacobian checks (--stop-after-checks, test option): no optimizer step, no "
            "products")
        raise SystemExit(0)
    fit_cops.least_squares = least_squares
    return original


def install_activity_diagnostics(fit_cops, model_cfg):
    """fit_cops.activity_diagnostics in this model's semantics: act = A; active_area = area with A > 0;
    clipped_area_fraction = area with A = 0; enhancement_* = A - A(m = 1)."""
    original = fit_cops.activity_diagnostics

    def activity_diagnostics(a, b, el, seen, coef, enhancement_basis, enhancement_scale):
        if not len(b):
            return original(a, b, el, seen, coef, enhancement_basis, enhancement_scale)
        nb, nt = nm.n_bins_templates(model_cfg)
        bins = nm.sm.bins_from_basis(np.asarray(enhancement_basis)[:, :nb])
        Bft = np.ascontiguousarray(np.asarray(enhancement_basis, np.float64)[:, nb:nb + nt])
        b = np.asarray(b, float)
        f = nm.host_fields(model_cfg["model"], el["Y"], np.asarray(a, float), Bft, b[nb:nb + nt],
                           float(enhancement_scale), bins, b[:nb], nm.width(model_cfg), nm.held_g(), False)
        f1 = nm.host_fields(model_cfg["model"], el["Y"], np.asarray(a, float), Bft, b[nb:nb + nt],
                            float(enhancement_scale), bins, np.ones(nb), nm.width(model_cfg), nm.held_g(), False)
        act, act1 = f["A"], f1["A"]
        active = act > 0.0
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
    st = facet_state(P, r["x"], el)
    path = Path(out) / "nonneg_facets.npz"
    model = P.nonneg_model["model"]
    meta = {"created_utc": utc_now(), "parameterization": P.nonneg_parameterization, "model": P.nonneg_model,
            "held": nm.stage_held(), "bin_names": list(nm.BIN_NAMES), "multipliers": [float(v) for v in st["m"]],
            "template_values": [float(v) for v in st["c"]], "kappa": float(st["kappa"]), "T": [float(v) for v in st["T"]],
            "A_scale": float(P.enhancement_scale), "a": [float(v) for v in st["a"]],
            "definitions": {"u": "sum_j a_j Y_j (fit_cops' 25 degree-<=4 Y)",
                            "TS": "A_s sum_k c_k T_k (D2's ft templates, c_k = 1)",
                            "A": ("m[bin] g(u + TS) exp(h), g = s ln(1 + exp(L/s)), s = s_ratio A_orig, max(0, L) where "
                                  "s = 0" if model == "softplus" else "m[bin] (u^2 / A_s + TS) exp(h)"),
                            "P": "A / m[bin] = q exp(h) (the multiplier column base)",
                            "bin": "slope bin 0..6", "seen": "facet contributes to a fit point (Model.seen_faces)",
                            "held": "the held fine-tuning term h of this fit (absent: 0)"}}
    extra = {}
    if st["h"] is not None:
        extra["held"] = st["h"]
    with path.open("xb") as f:
        np.savez_compressed(f, u=st["u"], TS=st["TS"], A=st["A"], P=st["P"], bin=st["bins"],
                            seen=np.asarray(seen, bool), meta_json=np.array(json.dumps(meta)), **extra)
    return path, st


BIN_TABLE_HEADER = """\
# nonneg_bins.txt -- written by model_ext/nonneg_problem.py ({utc})
# {pid}: {description}
# zero = A = 0 (no floor in this model); seen = facet contributes to a fit point; production share = sum_bin area*A /
# sum_all area*A (activity x area; q(SZA) not applied); *_seen over seen facets. sigma: the fitter's local GN sigma.
#bin interval_deg multiplier sigma facets area_fraction seen_area_fraction_of_bin zero_facets_frac zero_area_frac seen_zero_area_frac production_share production_share_seen
"""


def write_bin_table(out, P, r, el, seen, st):
    area = np.asarray(el["area"], float)
    seen = np.asarray(seen, bool)
    tot = float(np.sum(area))
    prod = area * st["A"]
    ptot, pseen = float(np.sum(prod)), float(np.sum(prod[seen]))
    z = st["A"] <= 0.0
    sig = {n: float(v) for n, v in zip(P.names, r["sig"])}
    rows = []
    for k, name in enumerate(nm.BIN_NAMES):
        sel = st["bins"] == k
        ss = sel & seen
        a_bin, a_seen = float(area[sel].sum()), float(area[ss].sum())
        rows.append({"bin": name, "interval_deg": nm.BIN_LABELS[k], "multiplier": float(st["m"][k]),
                     "sigma": sig.get("c_" + name), "facets": int(sel.sum()), "area_fraction": a_bin / tot,
                     "seen_area_fraction_of_bin": a_seen / a_bin if a_bin > 0 else float("nan"),
                     "zero_facets_frac": float(np.sum(sel & z)) / max(int(sel.sum()), 1),
                     "zero_area_frac": float(area[sel & z].sum()) / a_bin if a_bin > 0 else float("nan"),
                     "seen_zero_area_frac": float(area[ss & z].sum()) / a_seen if a_seen > 0 else float("nan"),
                     "production_share": float(prod[sel].sum()) / ptot if ptot > 0 else float("nan"),
                     "production_share_seen": float(prod[ss].sum()) / pseen if pseen > 0 else float("nan")})
    rows.append({"bin": "all", "interval_deg": "[0,180]", "multiplier": float("nan"), "sigma": None,
                 "facets": int(len(area)), "area_fraction": 1.0,
                 "seen_area_fraction_of_bin": float(area[seen].sum()) / tot,
                 "zero_facets_frac": float(np.sum(z)) / len(area), "zero_area_frac": float(area[z].sum()) / tot,
                 "seen_zero_area_frac": float(area[seen & z].sum()) / max(float(area[seen].sum()), 1e-300),
                 "production_share": 1.0, "production_share_seen": 1.0})
    path = Path(out) / "nonneg_bins.txt"
    with path.open("x") as f:
        f.write(BIN_TABLE_HEADER.format(utc=utc_now(), pid=P.nonneg_parameterization,
                                        description=DESCRIPTION[P.nonneg_model["model"]]))
        for row in rows:
            f.write("%s %s %s %s %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n" % (
                row["bin"], row["interval_deg"], "%.9e" % row["multiplier"] if row["bin"] != "all" else "-",
                ("fixed" if row["bin"] == nm.FIXED_BIN else "-" if row["sigma"] is None else "%.3e" % row["sigma"]),
                row["facets"], row["area_fraction"], row["seen_area_fraction_of_bin"], row["zero_facets_frac"],
                row["zero_area_frac"], row["seen_zero_area_frac"], row["production_share"],
                row["production_share_seen"]))
        f.write("# total production sum(area*A) = %.9e (seen %.9e)\n" % (ptot, pseen))
        lb = (cfg_of(P) or {}).get("slope_mult_min")
        if lb is not None:
            atb = at_bound(rows[:-1], lb)
            f.write("# multiplier lower bound %r (user 2026-11-12; --slope-mult-min) on the six free multipliers, no "
                    "upper bound; m - lb per bin %s; at the lower bound (m - lb <= %g max(1, lb); scipy TRF keeps its "
                    "iterates strictly inside, ~1e-10 above a bound it rests on): %s; sigma is the local unconstrained "
                    "Gauss-Newton value (not a two-sided uncertainty at the bound)\n" % (
                        lb, {r_["bin"]: "%.3g" % (r_["multiplier"] - lb) for r_ in rows[:-1] if r_["bin"] != nm.FIXED_BIN},
                        AT_BOUND_RTOL, atb if atb else "none"))
    return path, rows


def cfg_of(P):
    return getattr(P, "nonneg_cfg", None)


AT_BOUND_RTOL = 1e-8                 # reporting only: a multiplier within 1e-8 max(1, lb) of its lower bound is "at" it


def at_bound(rows, lb):
    return [r_["bin"] for r_ in rows if r_["bin"] != nm.FIXED_BIN and
            r_["multiplier"] - lb <= AT_BOUND_RTOL * max(1.0, abs(lb))]


TEXT = """\
# NONNEG_MODEL.txt -- written by model_ext/nonneg_problem.py ({utc})
# {pid}: {description}
# {policy}

Model {model}: {definition}
  Multipliers: 7 slope bins of {bundle} ({fixed} held at 1, six free >= 0); templates {templates} held at 1 (from {ftsrc}).
  Held fine-tuning term: {held}.
  Column factors: {columns}
Base: {run} ({nrows} rows; {ngrid} model-only grid rows). Fit data: {nng} NG, {nrg} RG. Late window >= {late} UTC
  ({nlate} observations): W10x5 + grain at the time stamp (D2's treatment, clip_problem.install_core).
Optimizer start: {start}

VALID: fit_result.txt (cost, parameters: a00..a24 = the harmonic coefficients of u = Y a{aunits}, c_slope_* = the
  multipliers), nonneg_facets.npz (u, TS, A, P, bin, seen, held), nonneg_bins.txt, NONNEG_MODEL.txt, evaluations.txt,
  grain_term.*, start_check_model.npz, fit_manifest.json (model_extension.nonneg), COPS_model.dat (observation rows =
  this model; grid rows as stated in ClipRefitGridRows), late_window_model.dat, fit_restart.npz for a resume of THIS
  fit (nonneg_fit.py; marker checked) and as the chain start of the next cycle.
NOT VALID (the unchanged fitter's additive clipped-linear layout): activity_linear in fit_restart.npz, ActivityFit.txt as
  a replay input, ActivityEnhancementFit.txt, ActivityEnhancementFacets.dat, activity_facets.dat, activity_maps.png,
  activity_*_power.dat, activity_region_diagnostics.dat, cops_fit.in, the positive/clip/slope fine-tuning tools (this
  model's fine tuning is nonneg_finetune.py).
"""


def write_text(out, P, cfg, record):
    lay = P.late_layout
    chain = record.get("fit_chain") or {}
    start = chain.get("start", "?")
    if chain.get("resumed_from"):
        start = "resumed from %s; the fit started with: %s" % (chain["resumed_from"][-1]["path"], start)
    sr = record.get("start_report") or {}
    if sr.get("cost") is not None:
        start += "\n  Cost at the start of this segment %.12e" % sr["cost"]
    model = P.nonneg_model["model"]
    definition = {"softplus": "A = m_b s ln(1 + exp(L/s)) exp(h), L = Y a + A_s sum_k T_k, s = %r A_orig (max(0, L) "
                              "where A_orig = 0)" % P.nonneg_model["s_ratio"],
                  "squared": "A = m_b (A_s (Y b)^2 + A_s sum_k T_k) exp(h), a = A_s b"}[model]
    columns = {"softplus": "dA/da_j = m g'(L) exp(h) Y_j, dA/dm_k = [bin = k] g(L) exp(h), g' = sigmoid(L/s)",
               "squared": "dA/da_j = m (2 u / A_s) exp(h) Y_j, dA/dm_k = [bin = k] (u^2/A_s + TS) exp(h)"}[model]
    held = nm.stage_held()
    lbv = cfg.get("slope_mult_min")
    text = TEXT.format(
        utc=utc_now(), pid=P.nonneg_parameterization, description=DESCRIPTION[model], policy=nm.POLICY[model],
        model=model, definition=definition, bundle=cfg["slope_bundle"], fixed=nm.FIXED_BIN,
        templates=[t["name"] for t in cfg["ft_templates"]], ftsrc=cfg["ft_templates_from"],
        held=("none" if held is None else "%s[%s] sha256 %s (%d nonzero facets, %.4g .. %.4g)" % (
            held["path"], held["key"], held["sha256"], held["nonzero_facets"], held["min"], held["max"])),
        columns=columns, run=str(P.M.run), nrows=len(lay.ipoint) if lay is not None else "?",
        ngrid=len(lay.grid_rows) if lay is not None else 0, nng=int(P.NNG), nrg=int(P.NRG), late=pp.LATE_START_UTC,
        nlate=len(lay.late_rows) if lay is not None else "?", start=start,
        aunits=" (= A_s b in the squared model)" if model == "squared" else "")
    if lbv is not None:
        x = np.asarray(record.get("final_x") or [], float)
        names = list(P.names)
        atb = [n for k, n in enumerate(names) if n.startswith("c_slope_") and len(x) and
               x[k] - lbv <= AT_BOUND_RTOL * max(1.0, abs(lbv))]
        text = text.replace("\nBase: ", "\nMULTIPLIER LOWER BOUND (user 2026-11-12, --slope-mult-min %r): the six free "
                            "multipliers have lower bound %r (not the fitter's 0), no upper bound, start 1; at the bound "
                            "in this fit (within 1e-8; scipy TRF keeps its iterates ~1e-10 inside): %s. fit_result.txt's header line 'template coefficients have lower bound zero' "
                            "does not hold for c_slope_* (the bounds are in fit_manifest.json lower_bounds); sigma is "
                            "the local unconstrained Gauss-Newton value.\nBase: " % (lbv, lbv, atb if atb else "none"),
                            1)
    path = Path(out) / "NONNEG_MODEL.txt"
    with path.open("x") as f:
        f.write(text)
    return path


def add_table_note(out, P):
    table = Path(out) / "COPS_model.dat"
    lines = table.read_text().splitlines(True)
    k = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    note = ("# NonnegModel (model_ext/nonneg_problem.py, %s): pNG/pRG of the observation rows = the %s model "
            "(NONNEG_MODEL.txt)%s; grid rows as stated in ClipRefitGridRows\n" % (
                P.nonneg_parameterization, P.nonneg_model["model"], " with the held term" if nm.stage_held() else ""))
    tmp = table.with_name(table.name + ".tmp.%d" % os.getpid())
    tmp.write_text("".join(lines[:k] + [note] + lines[k:]))
    os.replace(str(tmp), str(table))


_MARKER = {}


def install_restart_marker(record, cfg, model_cfg):
    """The restart's model_extension marker (the current call's record: a second in-process call replaces it)."""
    from cops_cuda import restart as rs
    _MARKER.update(record=record, cfg=cfg, model_cfg=model_cfg)
    if getattr(rs.from_fit, "_nonneg_marker", False):
        return
    original = rs.from_fit

    def from_fit(problem, result, fit_dir, manifest, elements):
        state = original(problem, result, fit_dir, manifest, elements)
        record, cfg, model_cfg = _MARKER["record"], _MARKER["cfg"], _MARKER["model_cfg"]
        lw = (record.get("clip_core") or {}).get("late_window") or {}
        state.meta["model_extension"] = {
            "parameterization": problem.nonneg_parameterization, "description": DESCRIPTION[model_cfg["model"]],
            "policy": nm.POLICY[model_cfg["model"]], "model": dict(model_cfg), "held": nm.stage_held(),
            "point_table": (record.get("clip_core") or {}).get("point_table"),
            "slope_bundle": {"path": cfg["slope_bundle"], "sha256": sha256(cfg["slope_bundle"])},
            "ft_templates": cfg["ft_templates"], "ft_source": cfg["ft_source"],
            "fit_chain": record.get("fit_chain"), "fixed": cfg.get("fixed"),
            "start_choice": (record.get("fit_chain") or {}).get("start_choice"),
            "late_window": {k: lw.get(k) for k in ("late_start_utc", "late_observations", "grid_rows",
                                                   "window_installed")},
            "warning": "activity_linear in this file is Y a + A_s B b (the fitter's additive layout), NOT this "
                       "model's activity (nonneg_facets.npz has A); its fine tuning is nonneg_finetune.py"}
        if cfg.get("slope_mult_min") is not None:
            state.meta["model_extension"]["slope_mult_min"] = cfg["slope_mult_min"]
        return state
    from_fit._nonneg_marker = True
    rs.from_fit = from_fit


# ------------------------------------------------------------------------------------------------ installation
def install_core(fit_cops, cfg, record, source=None, arrays=None):
    core = {}
    record["clip_core"] = core
    record["clip_core_note"] = ("clip_problem.install_core with variant D2 (grid rows, W10x5, start check, Jacobian "
                                "check): its labels name D2's late-window treatment, not this model")
    if cfg.get("stop_after_checks"):
        install_stop_after_checks(fit_cops, record)
    install_start_report(fit_cops, record)
    core_cfg = {"variant": "D2", "point_table": cfg["point_table"], "check_full_jacobian": cfg["check_full_jacobian"],
                "check_parameters": cfg["check_parameters"]}
    cp.install_core(fit_cops, core_cfg, core, source, arrays)
    return core


def install_outer(fit_cops, cfg, model_cfg, record, source=None, held=None):
    base = fit_cops.Problem
    pid = parameterization_of(model_cfg)
    nb, nt = nm.n_bins_templates(model_cfg)
    ftn = [t["name"] for t in cfg["ft_templates"]]

    class NonnegOuter(base):
        def __init__(self, M, el, coef, args):
            nm.set_A_orig(el["AH2O"])
            super().__init__(M, el, coef, args)
            self.nonneg_parameterization = pid
            self.nonneg_model = dict(model_cfg)
            self.nonneg_cfg = cfg
            names = list(self.names)
            free_m = ["c_" + n for n in nm.BIN_NAMES if n != nm.FIXED_BIN]
            expected = ["a%02d" % j for j in range(25)] + free_m + ["kappa"] + list(GRAIN_NAMES)
            if names != expected:
                raise SystemExit("nonneg: parameter layout %s, expected %s" % (names, expected))
            if list(self.template_names) != list(nm.BIN_NAMES) + ftn or self.nh:
                raise SystemExit("nonneg: templates %s (harmonic modes %d); expected %s" % (
                    self.template_names, self.nh, list(nm.BIN_NAMES) + ftn))
            basis = np.asarray(self.enhancement_basis)
            bins = nm.sm.bins_from_basis(basis[:, :nb])
            if nm.sha256_bytes(bins) != model_cfg["bins_sha256"] or len(bins) != model_cfg["nfaces"]:
                raise SystemExit("nonneg: the fitter's template basis is not the configured bin assignment")
            ah2o = np.ascontiguousarray(np.asarray(el["AH2O"], dtype=np.float64))
            if nm.sha256_bytes(ah2o) != model_cfg["ah2o_sha256"]:
                raise SystemExit("nonneg: AH2O of the base differs from the bundle's A_orig digest")
            a, b, T, kappa = self.unpack(self.x0)
            if [float(v) for v in T] != list(USER_T0):
                raise SystemExit("nonneg: held temperatures %s, expected %s" % (list(T), USER_T0))
            if float(b[nm.BIN_NAMES.index(nm.FIXED_BIN)]) != nm.FIXED_VALUE or np.any(np.asarray(b)[nb:] != 1.0):
                raise SystemExit("nonneg: %s and the ft templates must be held at 1 (got %s)" % (nm.FIXED_BIN, list(b)))
            Bft = basis[:, nb:nb + nt]
            counts = [int(v) for v in np.bincount(bins, minlength=len(nm.BIN_NAMES))]
            record["model_check"] = {"bins_sha256": model_cfg["bins_sha256"], "bin_counts": counts,
                                     "ah2o_sha256": model_cfg["ah2o_sha256"], "layout": names,
                                     "held_T": [float(v) for v in T], "fixed_bin": nm.FIXED_BIN,
                                     "model": model_cfg["model"], "s_ratio": model_cfg["s_ratio"],
                                     "templates": ftn, "template_columns_sha256": nm.sha256_bytes(
                                         np.ascontiguousarray(Bft)),
                                     "template_min": float(np.min(Bft)) if Bft.size else None,
                                     "A_scale": float(self.enhancement_scale), "held": held,
                                     "a_orig_zero_facets": int(np.sum(ah2o == 0.0)),
                                     "a_orig_zero_area": float(np.sum(np.asarray(el["area"])[ah2o == 0.0]) /
                                                               np.sum(el["area"]))}
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
                raise SystemExit("nonneg: no late layout (the point table is required)")
            if row_end != int(lay.late_rows[0]):
                raise SystemExit("nonneg: the first post-gap NG point (row %d) is not the first late-window observation "
                                 "(row %d): the window start needs a user decision" % (row_end, int(lay.late_rows[0])))
            say("model check: %s, s ratio %r, bins %s (sha256 %s...), templates %s held at 1 (min %s), A_orig digest "
                "ok (A_orig = 0 on %d facets, %.4f of the area), held T %s, %s = 1, A_s %.9e; held term %s; grain "
                "window t_end %s, gap %.1f s; %d NG, %d RG points" % (
                    model_cfg["model"], model_cfg["s_ratio"], counts, model_cfg["bins_sha256"][:16], ftn,
                    record["model_check"]["template_min"], record["model_check"]["a_orig_zero_facets"],
                    record["model_check"]["a_orig_zero_area"], list(T), nm.FIXED_BIN, float(self.enhancement_scale),
                    (held or {}).get("sha256", "none")[:16], gw["t_end_utc"], gw["gap_seconds"], gw["NNG"], gw["NRG"]))
            if record["model_check"]["template_min"] is not None and record["model_check"]["template_min"] < 0.0:
                raise SystemExit("nonneg: a template column is negative (the model needs T_k >= 0)")
            self._nonneg_el = el
            vmin = cfg.get("slope_mult_min")
            if vmin is not None:               # user 2026-11-12: the six free multipliers >= vmin (fitter: >= 0)
                lo = np.array(self.lo, dtype=float, copy=True)
                changed = []
                for k, name in enumerate(self.names):
                    if name in free_m:
                        if not (lo[k] == 0.0 and np.isinf(self.hi[k])):
                            raise SystemExit("nonneg: %s has bounds [%r, %r], expected the fitter's [0, inf]" % (
                                name, lo[k], self.hi[k]))
                        lo[k] = vmin
                        changed.append(name)
                self.lo = lo
                record["model_check"]["multiplier_lower_bound"] = {
                    "value": vmin, "parameters": changed, "upper_bound": "inf",
                    "policy": "user 2026-11-12: do not allow the multiplier factor decrease below 1; start all of "
                              "them with 1; they can grow (fitter default lower bound 0 replaced here)",
                    "scipy_trf_start_on_bound": "least_squares(method='trf') shifts a start on a bound to "
                                                "lb + 1e-10 max(1, |lb|) (make_strictly_feasible) before its first "
                                                "evaluation; the start check evaluates the start itself"}
                say("multiplier lower bound %r on %s (user 2026-11-12; the fitter's default 0 replaced), no upper "
                    "bound" % (vmin, changed))
            if cfg.get("start_restart"):
                ch = cfg["start_choice"]
                if ch == "d2":
                    self.restart_start = start_d2(self, el, cfg["start_restart"], record, cfg)
                elif ch == "lsq":
                    self.restart_start = start_lsq(self, el, cfg["start_restart"], record, cfg)
                else:
                    self.restart_start = start_chain(self, el, cfg["start_restart"], record, cfg["grain_from"],
                                                     model_cfg, held, vmin)
            if vmin is not None and getattr(self, "restart_start", None) is not None:
                xs = np.asarray(self.restart_start, float)
                bad = [n for k, n in enumerate(self.names) if xs[k] < self.lo[k] or xs[k] > self.hi[k]]
                if bad:
                    raise SystemExit("nonneg: the start of %s is outside the bounds" % bad)

    fit_cops.Problem = NonnegOuter
    install_activity_diagnostics(fit_cops, model_cfg)
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
        del st
        tpath = write_text(out, P, cfg, record)
        record.setdefault("products", {}).update({p.name: sha256(p) for p in (fpath, bpath, tpath)})
        record["bins_final"] = rows
        say("products: %s, %s, %s; multipliers %s" % (fpath.name, bpath.name, tpath.name,
                                                       ["%.6g" % v["multiplier"] for v in rows[:-1]]))
        return diag
    fit_cops.write_fit = write_fit
    install_restart_marker(record, cfg, model_cfg)
    return NonnegOuter


def install_hooks(grain_crossing, fixed_parameters, cfg, model_cfg, command, sub, source=None, arrays=None, held=None):
    """As clip_problem.install_hooks: install_core just before grain_crossing.install, install_outer just after
    fixed_parameters.install. Re-entrant: the wrappers always wrap the modules' ORIGINAL install functions (kept on
    the modules at the first call), so a second in-process fit does not stack hooks."""
    original_gc = getattr(grain_crossing, "_nonneg_original_install", None) or grain_crossing.install
    original_fp = getattr(fixed_parameters, "_nonneg_original_install", None) or fixed_parameters.install
    grain_crossing._nonneg_original_install = original_gc
    fixed_parameters._nonneg_original_install = original_fp
    main_record = {}

    def gc_install(fit_cops, weight, record):
        main_record["record"] = record
        record["nonneg"] = sub
        install_core(fit_cops, cfg, sub, source, arrays)
        result = original_gc(fit_cops, weight, record)
        install_eval_log(fit_cops, sub)
        return result

    def fp_install(fit_cops, fixed, fp_record):
        result = original_fp(fit_cops, fixed, fp_record)
        install_outer(fit_cops, cfg, model_cfg, sub, source, held)
        record = main_record["record"]
        record["grain_fit_fixed_command"] = record.get("command")
        record["command"] = list(command)
        return result
    grain_crossing.install = gc_install
    fixed_parameters.install = fp_install
    return original_gc, original_fp
