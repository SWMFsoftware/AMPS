# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/clip_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/clip_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 49371152ffc17bffdfb9ea15ed0d9824851ef5982cba391afaebbe0960bf6e43
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module of clip_refit.py: start restart, start check, grid-row handling, W10x5 hooks (positive_problem's own
#   functions), products (CLIP_REFIT.txt, late_window_model.dat). 2026-11-08.
# ==== end of header ====
"""Refit of the earlier CLIPPED-LINEAR model (cycle c6, job 122221) on the lategrid base: variants D1 and D2
(user decision 2026-11-08: "refit the earlier clipped model (original temperature, with its templates) on the new base,
with W10x5 and the +1.5 m offset", to separate the loss of fit quality caused by the more realistic gauge treatment
from the loss caused by the positive-activity parameterization). Fit-level part; the driver is clip_refit.py.

Local, untracked extension: fit_cops.py, gpu_fit.py, grain_fit.py, grain_crossing.py, grain_fit_fixed.py,
fixed_parameters.py and positive_problem.py are NOT edited; clip_refit.py runs the unchanged grain_fit_fixed.py path and
installs this module around it (install_core before grain_crossing.install, install_outer after fixed_parameters.install).

The model (= c6, runs/shap8_012m_flat_earlier/grain/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6)
* H2O activity max(0, Y a + A_scale (c_cliff T_cliff + sum_k 1 * ft_ck)) (the fitter's clipped linear model; 25
  degree-<=4 coefficients free, template cliff_ge_60deg free >= 0, ft_c1..ft_c6 held at 1 with fixed_parameters);
  temperatures Tmax 182.09922 / T75 136.11149 / Tmin 133 K held; kappa free; dust grain in NG + crossing weight 10
  (grain_crossing.py); lambda 0; the srcCG night/shadow ratios of the GPU model (no COPS_POS_* environment: refused by
  the driver); no flux floor. Everything comes from the command line of the job (jobs/clip_refit.sbatch = c6's command
  with the new base); this module checks that a start restart has the same held values.
* New base: runs/shap8_012m_lategrid_dz1.5/base (gauge offset +1.5 m along +Z, 1150 observations + 401 model-only 1-s
  grid rows 10:32:46-10:39:26 UTC with reference 0 and use 0; point table late_grid_points.txt).

Variants
* D1 (instantaneous): the model at the observation time stamps, as in c6. The fit model M has the fitter's own rows
  (the observation rows the profile selects; the grid rows have use 0 and are NOT loaded into M); the all-rows model
  MF (1551 rows, its own GPU backend) writes COPS_model.dat (grid rows: instantaneous surface + grain, model-only).
* D2 (W10x5): D1 plus the COPS time response for the late observations (>= 2016-09-30T10:33:10 UTC), EXACTLY as
  model_ext/positive_problem.py does it, with its own functions (LateLayout, w10x5_operator, WindowOperator,
  install_window; not copies): W10x5 = mean of five 10-s running averages of the 1-s grid model ending at T, T-2, ...,
  T-8 s, T = the observation's L2 STOP time (point table); the grid rows join M's rows (M = MF rows, one backend);
  M.evaluate and MF.evaluate replace pNG, pRG and every Jacobian column at the late observations by W @ (grid rows);
  the dust grain is added afterwards at the observation time stamp, unaveraged (grain_crossing.py adds it to the
  evaluated surface model). Earlier observations keep the instantaneous model.

Grid rows (both variants; checked at run time and recorded in fit_manifest.json model_extension.clip_refit)
* selection: the named profile (all_available) audits the 1150-row mission dataset, so it is applied to the observation
  rows only and the grid rows get use 0, reason "late_grid_point" (positive_problem.install_selection, unchanged);
* cost: fit_cops.Problem uses a row only if use > 0 AND reference > 0; the grid rows have both 0, so they are not in
  the data residuals, the gauge weights (1/sqrt(N): N = 934 NG, 1150 RG as in c6) or the Jacobian rows;
* crossing pairs (grain_crossing.pairs): rows where NG and RG are both used -> no grid row;
* grain term: t_end (first post-gap NG point) and the grain start search use the used NG rows only; the grain enters
  only the residuals of used NG rows; in the written COPS_model.dat the grid rows get surface + grain(t) (as positive
  fits: model-only, flagged);
* products: COPS_model.dat keeps the grid rows with a header note (ClipRefitGridRows), late_window_model.dat lists the
  late window row by row (kind grid | late_obs, same columns as a positive fit's; positive_table.py removes the grid
  rows for plots), pressure_residuals.* / pressure_comparison_* and the gpu_fit figures have the observation rows only.

Start (both variants): --start-restart <c6 fit_restart.npz>: the optimizer starts from c6's final x by NAME
(cops_cuda.restart.start_vector; same ascale -> exact optimizer coordinates; a, c_cliff, kappa, grain A and B); the
prior P.x0 (the fitter's default) is unchanged, as in c6. Refused unless the mesh, the parameter names, the template
names AND columns, the held temperatures and held template amplitudes, the enhancement scale and the run coefficients
are those of the source (the source's base differs: that is the point). Not gpu_fit --restart-from: that is a resume
(it carries the source's 6 segments into this fit's restart and segment count). Before the first optimizer step the
START CHECK prints the cost at the start (data, crossing penalty, D) next to the source's final cost (c6: 4.781872e-03
on the old base) and the change of the model predictions (by UTC; all / before / from 10:33:10, NG and RG) against the
source's stored predictions, and checks that the data (references, use flags) are the source's. start_check_model.npz
(fit rows at the start: surface model, grain, flags) allows D1 vs D2 row comparisons.
Resume after a wall-time stop: gpu_fit --restart-from <this fit's seg<k>/fit_restart.npz>; its marker (meta
model_extension) must have this variant's parameterization and point table.
"""

import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

import positive_problem as pp      # W10x5, selection wrapper, figure filter, point-table helpers (used unchanged)

VARIANTS = ("D1", "D2")
PARAMETERIZATION = {"D1": "clip_refit_c6_D1_instantaneous_v1", "D2": "clip_refit_c6_D2_w10x5_v1"}
DESCRIPTION = {
    "D1": "clipped linear c6 model on the lategrid base (+1.5 m gauge offset); model at the observation time stamps",
    "D2": "clipped linear c6 model on the lategrid base (+1.5 m gauge offset); W10x5 for the observations from "
          "%s UTC (positive_problem's operator), grain unaveraged at the time stamp" % pp.LATE_START_UTC}
POLICY = ("user 2026-11-08: refit the earlier clipped model (c6: original temperature held, cliff_ge_60deg free >= 0, "
          "ft_c1..ft_c6 held at 1, kappa free, grain + crossing 10, lambda 0, srcCG night/shadow ratios, no floor) on the "
          "new base (gauge +1.5 m), D1 instantaneous, D2 with W10x5 as positive_problem; start = c6's final a, c_cliff, "
          "kappa, grain (by name)")
CHECK_PARAMETERS = ("a00", "a02", "a07", "c_cliff_ge_60deg", "kappa", "grain_A_uPa")
GRAIN_NAMES = ("grain_A_uPa", "grain_B_uPa_per_100s")
RESUME_RTOL = 1e-8      # a resume (gpu_fit --restart-from this fit's segment) must reproduce the restart's final cost
TNAMES = ("Tmax", "T75", "Tmin")
LOG = "clip_refit:"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


jsonable = pp.jsonable


def say(text):
    print("%s %s" % (LOG, text), flush=True)


# ------------------------------------------------------------------------------------------------ arguments
def _pop(head, name, default=None):
    if name in head:
        i = head.index(name)
        if i + 1 >= len(head):
            raise SystemExit("clip_refit: %s needs a value" % name)
        value = head[i + 1]
        del head[i:i + 2]
        return value
    return default


def _flag(head, name):
    if name in head:
        head.remove(name)
        return True
    return False


def _value(tokens, name):
    for k, token in enumerate(tokens):
        if token == name and k + 1 < len(tokens):
            return tokens[k + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def fixed_values(head):
    """{name: value} of the grain_fit_fixed --fix options (left in place)."""
    out = {}
    for k, token in enumerate(head):
        if token == "--fix" and k + 1 < len(head):
            name, value = head[k + 1].split("=")
            out[name] = float(value)
    return out


def rewrite_arguments(argv):
    """Split the clip_refit.py command into (config, grain_fit_fixed.py argv, notes).

    Driver options (before '--'): --variant D1|D2 (required), --point-table PATH (the base's late_time_grid.py table;
    required for D2 and on a base with grid rows), --start-restart R.npz (a new fit from another fit's final state, by
    name; excludes gpu_fit --restart-from = resume), --check-full-jacobian (analytic vs central-difference columns of
    the full objective at the optimizer start: data incl. the late window, crossing rows), --check-parameters a,b,...
    Everything else goes to grain_fit_fixed.py (--fix ...), grain_fit.py (--crossing-weight, --check-jacobian) and
    gpu_fit.py (--devices, --precision, --max-nfev, --wall-time, --restart-from, --segment, ...).
    Inserted when absent: gpu_fit --template-start 0 (c6's value; a nonzero value would replace the start's template
    amplitudes). Refused: fitter options that contradict the c6 path (see the messages)."""
    argv = list(argv)
    if "--" not in argv:
        raise SystemExit("clip_refit: need '--' before the fit_cops arguments")
    split = argv.index("--")
    head, tail = argv[:split], argv[split + 1:]
    notes = []
    cfg = {"variant": _pop(head, "--variant"), "point_table": _pop(head, "--point-table"),
           "start_restart": _pop(head, "--start-restart"),
           "check_full_jacobian": _flag(head, "--check-full-jacobian"),
           "check_parameters": list(CHECK_PARAMETERS)}
    names = _pop(head, "--check-parameters")
    if names:
        cfg["check_parameters"] = [n.strip() for n in names.split(",") if n.strip()]
    if cfg["variant"] not in VARIANTS:
        raise SystemExit("clip_refit: --variant D1 or D2 is required (got %r)" % cfg["variant"])
    if cfg["point_table"] is not None:
        cfg["point_table"] = str(Path(cfg["point_table"]).resolve())
        if not Path(cfg["point_table"]).is_file():
            raise SystemExit("clip_refit: point table not found: %s" % cfg["point_table"])
    elif cfg["variant"] == "D2":
        raise SystemExit("clip_refit: D2 (W10x5) needs --point-table <run>/late_grid_points.txt")
    if cfg["start_restart"] is not None:
        cfg["start_restart"] = str(Path(cfg["start_restart"]).resolve())
        if not Path(cfg["start_restart"]).is_file():
            raise SystemExit("clip_refit: --start-restart file not found: %s" % cfg["start_restart"])
    cfg["restart_from"] = _value(head, "--restart-from")
    if cfg["restart_from"] is not None:
        cfg["restart_from"] = str(Path(cfg["restart_from"]).resolve())
        if cfg["start_restart"] is not None:
            raise SystemExit("clip_refit: --start-restart (a new fit from another fit's state) and gpu_fit "
                             "--restart-from (resume of this fit's previous segment) exclude each other")
    ts = _value(head, "--template-start")
    if ts is None:
        head += ["--template-start", "0"]
        notes.append("inserted gpu_fit --template-start 0 (as c6; the start comes from the start restart)")
    elif float(ts) != 0.0:
        raise SystemExit("clip_refit: --template-start must be 0 (c6's value): a nonzero value would replace the "
                         "start restart's template amplitudes")
    cfg["fixed"] = fixed_values(head)
    for bad, why in (("--restart-from", "give --restart-from to gpu_fit (before '--')"),
                     ("--fix-temperature", "the temperatures are held with --fix (fixed_parameters), as in c6"),
                     ("--ordered-temperature", "the temperatures are held with --fix (fixed_parameters), as in c6"),
                     ("--check-jacobian", "the fitter's check stops after checking; use grain_fit --check-jacobian "
                                          "or clip_refit --check-full-jacobian (before '--')"),
                     ("--data-only", "not a fit"),
                     ("--fix-legacy-activity", "c6 fits all 25 coefficients")):
        if bad in tail or any(t.startswith(bad + "=") for t in tail):
            raise SystemExit("clip_refit: %s is not allowed after '--': %s" % (bad, why))
    cfg["templates"] = [tail[k + 2] for k, t in enumerate(tail) if t == "--activity-template" and k + 2 < len(tail)]
    run = next((t for t in tail if not t.startswith("-")), None)
    if run is None:
        raise SystemExit("clip_refit: no run directory after '--'")
    cfg["run"] = run
    return cfg, head + ["--"] + tail, notes


# ------------------------------------------------------------------------------------------------ restarts
def restart_meta(path):
    with np.load(str(path), allow_pickle=False) as z:
        return json.loads(str(z["meta_json"]))


def load_source(path):
    """What the start check needs from the start restart: identity, final cost, segments, the source fit's manifest
    numbers (data cost, crossing penalty, D, t_end) and its stored per-row arrays (UTC, references, flags,
    predictions)."""
    meta = restart_meta(path)
    ext = meta.get("model_extension") or {}
    if str(ext.get("parameterization", "")).startswith("positive"):
        raise SystemExit("clip_refit: %s is a positive (exp) fit's restart, not a clipped-linear one" % path)
    with np.load(str(path), allow_pickle=False) as z:
        arrays = {k: np.array(z[k]) for k in ("point_ids", "utc", "ref_ng", "ref_rg", "use_ng", "use_rg", "pred_ng",
                                               "pred_rg", "x")}
    opt = meta.get("optimizer") or {}
    info = {"path": str(Path(path).resolve()), "sha256": sha256(path), "cost": float(opt["cost"]),
            "nfev": opt.get("nfev"), "message": opt.get("message"), "status": meta.get("status"),
            "lambda": float((meta.get("settings") or {}).get("lambda", 0.0)),
            "segments": len(meta.get("segments") or []), "fit_directory": meta.get("fit_directory"),
            "base": meta.get("base"), "marker": ext or None}
    manifest = Path(meta.get("fit_directory") or "/nonexistent") / "fit_manifest.json"
    try:
        m = json.loads(manifest.read_text())
        o = [v for v in m["optimizer_results"] if float(v["lambda"]) == info["lambda"]][0]
        me = m.get("model_extension") or {}
        info.update(manifest=str(manifest), data_cost=float(o["data_cost"]), D=float(o["deviation_from_prior"]),
                    penalty_cost=float(me["final"]["penalty_cost"]), t_end_et=me.get("t_end_et"),
                    gap_seconds=me.get("gap_seconds"), pairs_before_gap=me["final"].get("pairs_before_gap"),
                    pairs_after_gap=me["final"].get("pairs_after_gap"))
    except Exception as exc:                          # recorded and printed, never silent
        info["manifest_error"] = "%s: %s" % (exc.__class__.__name__, exc)
    return info, arrays


def check_restart_marker(path, cfg):
    """--restart-from (resume): the restart must be this driver's, with this variant and point table."""
    meta = restart_meta(path)
    ext = meta.get("model_extension") or {}
    want = PARAMETERIZATION[cfg["variant"]]
    if ext.get("parameterization") != want:
        raise SystemExit("clip_refit: --restart-from %s has the marker %r, this fit is %r: a resume must continue the "
                         "same problem" % (path, ext.get("parameterization"), want))
    table = (ext.get("point_table") or {}).get("sha256")
    now = sha256(cfg["point_table"]) if cfg.get("point_table") else None
    if table != now:
        raise SystemExit("clip_refit: --restart-from %s was written with the point table %s, this process has %s"
                         % (path, table, now))
    return {"path": str(Path(path).resolve()), "sha256": sha256(path), "parameterization": want,
            "status": meta.get("status"), "segments": meta.get("segments"), "fit_directory": meta.get("fit_directory"),
            "fit_chain": ext.get("fit_chain"), "cost": float((meta.get("optimizer") or {}).get("cost", float("nan")))}


def fit_chain(cfg, source, resume):
    if source is not None:
        keep = ("path", "sha256", "fit_directory", "cost", "nfev", "status", "segments")
        return {"start": "another fit's final state (clip_refit.py --start-restart)",
                "start_restart": {k: source.get(k) for k in keep}}
    if resume is not None:
        prev = dict(resume.get("fit_chain") or {})
        prev["resumed_from"] = list(prev.get("resumed_from") or []) + [
            {k: resume.get(k) for k in ("path", "sha256", "fit_directory", "cost", "status")}]
        return prev
    return {"start": "the fitter's default start (no restart)"}


def start_from_restart(P, el, path, record):
    """Optimizer start of P (the reduced problem of fixed_parameters) from the start restart, by name; every
    identity of the model except the base is checked. Returns x."""
    from cops_cuda import reference
    from cops_cuda import restart as rs
    state = rs.read(path)
    problems = []
    if list(state.meta.get("names") or []) != list(P.names):
        problems.append("parameters %s vs this fit's %s" % (state.meta.get("names"), list(P.names)))
    mesh = reference.load("activity_templates").mesh_identity(el["centroid"], el["area"])
    if mesh != state.mesh_identity():
        problems.append("mesh %s vs the base's %s" % (state.mesh_identity()[:16], mesh[:16]))
    cols = np.asarray(P.enhancement_basis[:, P.nh:], float)
    if state.template_names != list(P.template_names):
        problems.append("templates %s vs %s" % (state.template_names, list(P.template_names)))
    elif not np.array_equal(state.template_columns, cols):
        problems.append("template columns differ (same names)")
    if P.nh or len(state.harmonic_modes):
        problems.append("harmonic residual modes present (%d / %d): not the c6 layout" % (P.nh, len(state.harmonic_modes)))
    if float(state.meta.get("enhancement_scale")) != float(P.enhancement_scale):
        problems.append("enhancement scale %r vs %r" % (state.meta.get("enhancement_scale"), P.enhancement_scale))
    if float(state.meta.get("ascale")) != float(P.ascale):
        problems.append("ascale %r vs %r" % (state.meta.get("ascale"), P.ascale))
    if not np.array_equal(np.asarray(state.coef_prior, float), np.asarray(P.coef, float)):
        problems.append("the run's activity coefficients differ")
    if problems:
        raise SystemExit("clip_refit: --start-restart %s does not match this model: %s" % (path, "; ".join(problems)))
    x, report = rs.start_vector(P, state)
    if report["clipped"] or report["kept_problem_start"]:
        raise SystemExit("clip_refit: --start-restart %s: clipped %s, not in the restart %s (the start must be the "
                         "source's state)" % (path, report["clipped"], report["kept_problem_start"]))
    a, b, T, kappa = P.unpack(x)
    held = {"a": bool(np.array_equal(a, np.asarray(state.a, float))),
            "template_values": bool(np.array_equal(b[P.nh:], np.asarray(state.template_values, float))),
            "T": bool(np.array_equal(np.asarray(T, float), np.asarray(state.T, float))),
            "kappa": float(kappa) == float(np.asarray(state.kappa).reshape(-1)[0])}
    if not all(held.values()):
        raise SystemExit("clip_refit: the physical state at the start differs from the source's (%s): the held values "
                         "(--fix) are not the source's (T %s vs %s; templates %s vs %s)" % (
                             held, list(T), list(state.T), list(b[P.nh:]), list(state.template_values)))
    act = el["Y"] @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    dact = float(np.max(np.abs(act - state.activity_linear)) / max(float(np.max(np.abs(state.activity_linear))), 1e-300))
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_vector": jsonable(report),
           "x_equals_source_x": bool(np.array_equal(x, np.asarray(state.x, float))),
           "physical_state_equal": held, "activity_linear_max_rel_diff": dact,
           "source_base": state.meta.get("base"), "source_fit_directory": state.meta.get("fit_directory"),
           "source_x": [float(v) for v in state.x], "prior_x0_unchanged": True,
           "held_T": [float(v) for v in T], "held_template_values": [float(v) for v in b[P.nh:]]}
    record["start_restart"] = out
    say("start restart %s: x by name (%d parameters, exact optimizer coordinates %s), held T %s and template values %s "
        "= the source's, linear activity max rel. diff %.1e; the prior x0 is the fitter's default (unchanged)" % (
            path, len(P.names), out["x_equals_source_x"], list(T), list(b[P.nh:]), dact))
    del state
    return x


# ------------------------------------------------------------------------------------------------ checks
def install_jacobian_check(fit_cops, names, record):
    """Inside the crossing wrapper (installed before grain_crossing.install): analytic vs central-difference Jacobian
    columns of the full objective (data with grain and, for D2, the window; crossing rows) at the optimizer start,
    once (positive_problem.jacobian_check's definition, this module's log prefix)."""
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"] and callable(jac):
            state["done"] = True
            record["jacobian_check"] = jacobian_check(fun, jac, x0, fit_cops.P, names)
            fit_cops.save_active_manifest()
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


def jacobian_check(fun, jac, x0, P, names):
    x = np.array(x0, dtype=float)
    r0 = np.asarray(fun(x), float)
    J = np.asarray(jac(x), float)
    nd = int(P.ndata)
    nreg = len(P.regularization(x)) if P.lam > 0 else 0
    late = np.concatenate([P.et_fit[P.useNG] >= P.et_end, P.et_fit[P.useRG] >= P.et_end])
    blocks = {"data": np.arange(nd), "data_late_window": np.nonzero(late)[0],
              "regularization": np.arange(nd, nd + nreg), "crossing": np.arange(nd + nreg, len(r0))}
    out = {"lambda": float(P.lam), "rows": {k: int(len(v)) for k, v in blocks.items()},
           "metric": "max|analytic - central FD| / max|FD| per block (abs: max|analytic - FD|); step 1e-6 max(|x|,1)",
           "columns": {}}
    for name in names:
        if name not in P.names:
            out["columns"][name] = "not a parameter"
            continue
        k = P.names.index(name)
        h = 1e-6 * max(abs(x[k]), 1.0)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        fd = (np.asarray(fun(xp), float) - np.asarray(fun(xm), float)) / (2.0 * h)
        col = {"h": h}
        for b, idx in blocks.items():
            if len(idx):
                scale = float(np.max(np.abs(fd[idx])))
                err = float(np.max(np.abs(fd[idx] - J[idx, k])))
                col[b] = err / scale if scale > 0 else (0.0 if err == 0 else float("inf"))
                col[b + "_abs"] = err
        out["columns"][name] = col
        say("Jacobian check d/d%-18s h %.2e  data %.2e  late %.2e  crossing %s" % (
            name, h, col.get("data", float("nan")), col.get("data_late_window", float("nan")),
            "%.2e" % col["crossing"] if "crossing" in col else "-"))
    worst = [v["data"] for v in out["columns"].values() if isinstance(v, dict) and "data" in v]
    out["worst_data"] = float(max(worst)) if worst else None
    return out


def install_start_check(fit_cops, cfg, record, source, arrays):
    """At the first least_squares call, inside the crossing wrapper (fun = data + crossing rows): the start check."""
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"]:
            state["done"] = True
            record["start_check"] = res = start_check(fun, x0, fit_cops, cfg, source, arrays)
            fit_cops.save_active_manifest()
            if res.get("failed"):
                raise SystemExit("clip_refit: start check failed: %s" % res["failed"])
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


def _stats(d):
    if not len(d):
        return {"points": 0}
    k = int(np.argmax(np.abs(d)))
    return {"points": int(len(d)), "max_abs_dln": float(abs(d[k])), "rms_dln": float(np.sqrt(np.mean(d * d))),
            "mean_dln": float(np.mean(d)), "at": int(k)}


def start_check(fun, x0, fit_cops, cfg, source, arrays):
    P, M = fit_cops.P, fit_cops.P.M
    x = np.array(x0, dtype=float)
    r = np.asarray(fun(x), float)
    nd = int(P.ndata)
    nreg = len(P.regularization(x)) if P.lam > 0 else 0
    rd, rr, rc = r[:nd], r[nd:nd + nreg], r[nd + nreg:]
    data, reg, pen = (0.5 * float(v @ v) for v in (rd, rr, rc))
    cost, D = data + reg, float(P.deviation(x))
    # the model at the start on the fit rows (surface; NG + grain), as the residuals see it
    a, b, T, kappa = P.unpack(x)
    pNG_s, pRG_s = M.evaluate(a, T, kappa, b=b if P.nb else None, enhancement_basis=P.enhancement_basis,
                              enhancement_scale=P.enhancement_scale)
    g = P.grain(x, P.et_fit)[0]
    pNG = np.maximum(pNG_s + g, fit_cops.FLOOR)
    pRG = np.maximum(pRG_s, fit_cops.FLOOR)
    res_ng = P.wNG * (np.log(pNG[P.useNG]) - P.lnNG)
    res_rg = P.wRG * (np.log(pRG[P.useRG]) - P.lnRG)
    self_consistent = float(np.max(np.abs(np.concatenate([res_ng, res_rg]) - rd)))
    points = pp.read_run_points(M.run)
    ip = np.asarray(M.ipoint)[np.asarray(M.rows)]
    utc = [points[int(i)][0] for i in ip]
    keys = np.array([pp.time_key(u) for u in utc], dtype=np.int64)
    late_key = pp.time_key(pp.LATE_START_UTC)
    out_dir = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent
    npz = out_dir / "start_check_model.npz"
    with npz.open("xb") as f:
        np.savez_compressed(f, ipoint=ip, rows=np.asarray(M.rows), time_key_ms=keys, et=np.asarray(P.et_fit),
                            surface_pNG=np.asarray(pNG_s, float), surface_pRG=np.asarray(pRG_s, float), grain=g,
                            useNG=np.asarray(P.useNG), useRG=np.asarray(P.useRG), x=x,
                            meta_json=np.array(json.dumps({"variant": cfg["variant"], "names": list(P.names),
                                                           "created_utc": utc_now(), "run": str(M.run)})))
    out = {"variant": cfg["variant"],
           "start": {"cost": cost, "data_cost": data, "regularization_cost": reg, "crossing_penalty_cost": pen, "D": D,
                     "NNG": int(P.NNG), "NRG": int(P.NRG)},
           "residuals_vs_direct_model_max_abs": self_consistent, "model_file": str(npz)}
    if source is not None:
        out["source"] = {k: v for k, v in source.items() if k not in ("marker", "base")}
        out["cost_rel_diff_vs_source_final"] = cost / source["cost"] - 1.0
        if "data_cost" in source:
            out["data_cost_rel_diff_vs_source"] = data / source["data_cost"] - 1.0
            out["D_diff_vs_source"] = D - source["D"]
        src_keys = np.array([pp.time_key(str(u)) for u in arrays["utc"]], dtype=np.int64)
        pos = {int(k): j for j, k in enumerate(src_keys)}
        comp = {}
        for gauge, use, model, ln_ref, sref, suse, spred in (
                ("NG", P.useNG, pNG, P.lnNG, arrays["ref_ng"], arrays["use_ng"], arrays["pred_ng"]),
                ("RG", P.useRG, pRG, P.lnRG, arrays["ref_rg"], arrays["use_rg"], arrays["pred_rg"])):
            k_use = keys[use]
            missing = [int(k) for k in k_use if int(k) not in pos]
            if missing:
                comp[gauge] = {"error": "%d used rows have no source row (e.g. time key %d)" % (len(missing), missing[0])}
                continue
            j = np.array([pos[int(k)] for k in k_use], dtype=np.int64)
            dl = np.log(model[use]) - np.log(np.maximum(spred[j], fit_cops.FLOOR))
            dref = ln_ref - np.log(np.maximum(sref[j], fit_cops.FLOOR))
            src_used = (np.asarray(suse) > 0) & (np.asarray(sref) > 0)
            early = k_use < late_key
            comp[gauge] = {"all": _stats(dl), "before_late_window": _stats(dl[early]),
                           "late_window": _stats(dl[~early]),
                           "data_ln_ref_max_abs_diff": float(np.max(np.abs(dref))) if len(dref) else 0.0,
                           "used_rows_equal_source_used_rows": bool(np.array_equal(np.sort(k_use),
                                                                                   np.sort(src_keys[src_used])))}
        out["predictions_vs_source"] = dict(comp, note="source = its COPS_model.dat (%.6e): |d ln p| ~ 5e-7 is "
                                                       "rounding; d ln p = ln(model at the start) - ln(source "
                                                       "prediction), matched by UTC")
    say("start check (variant %s): cost at the start %.12e (data %.12e, crossing penalty %.3e, D %.9f; %d NG, %d RG "
        "points; residuals vs direct model %.1e)" % (cfg["variant"], cost, data, pen, D, P.NNG, P.NRG,
                                                     self_consistent))
    if source is not None:
        say("start check: source %s final cost %.12e (data %s, crossing penalty %s, D %s) on its base %s; relative "
            "difference of the start cost %+.4e" % (
                source["path"], source["cost"], "%.12e" % source["data_cost"] if "data_cost" in source else "?",
                "%.3e" % source["penalty_cost"] if "penalty_cost" in source else "?",
                "%.9f" % source["D"] if "D" in source else "?",
                ((source.get("base") or {}).get("base_model_table") or {}).get("path"), out["cost_rel_diff_vs_source_final"]))
        pv = out["predictions_vs_source"]
        for gauge in ("NG", "RG"):
            v = pv[gauge]
            if "error" in v:
                say("WARNING start check %s: %s" % (gauge, v["error"]))
                continue
            say("start check: %s model at the start vs the source's predictions: all %d max |d ln p| %.3e rms %.3e "
                "mean %+.3e; before %s: max %.3e rms %.3e mean %+.3e; from %s (%d): max %.3e rms %.3e mean %+.3e; "
                "data ln ref max |diff| %.1e, used rows = source's %s" % (
                    gauge, v["all"]["points"], v["all"]["max_abs_dln"], v["all"]["rms_dln"], v["all"]["mean_dln"],
                    pp.LATE_START_UTC[11:], v["before_late_window"].get("max_abs_dln", 0.0),
                    v["before_late_window"].get("rms_dln", 0.0), v["before_late_window"].get("mean_dln", 0.0),
                    pp.LATE_START_UTC[11:], v["late_window"]["points"], v["late_window"].get("max_abs_dln", 0.0),
                    v["late_window"].get("rms_dln", 0.0), v["late_window"].get("mean_dln", 0.0),
                    v["data_ln_ref_max_abs_diff"], v["used_rows_equal_source_used_rows"]))
        if "manifest_error" in source:
            say("WARNING start check: the source's fit_manifest.json could not be read: %s" % source["manifest_error"])
        if source.get("role") == "restart_from":       # a resume continues the same problem: it must reproduce
            ok = abs(out["cost_rel_diff_vs_source_final"]) <= RESUME_RTOL
            out["resume_reproduction"] = {"tolerance": RESUME_RTOL, "ok": bool(ok)}
            say("start check: resume of the same problem: reproduction %s (|rel| %.1e, tolerance %.0e)" % (
                "CONFIRMED" if ok else "FAILED", abs(out["cost_rel_diff_vs_source_final"]), RESUME_RTOL))
            if not ok:
                out["failed"] = "the resume does not reproduce the restart's final cost"
    return out


# ------------------------------------------------------------------------------------------------ products
def install_pressure_filter(fit_cops, gkeys, record):
    """The A5 pressure products (pressure_residuals.dat/png, pressure_comparison_*.png) for the observation rows only;
    the activity products are the fitter's (clipped linear: valid for this model). As the pressure part of
    positive_problem.install_activity_writers."""
    original_a5 = fit_cops.write_stage_a5_report

    def write_stage_a5_report(*args, **kwargs):
        original_wpp = fit_cops.write_pressure_products

        def write_pressure_products(directory, identity, label, point_ids, utc_text, altitude, predictions,
                                    references, used, saturated, reasons, transition_intervals=()):
            idx = [i for i, u in enumerate(utc_text) if pp.time_key(u) not in gkeys]
            if len(idx) == len(utc_text):
                return original_wpp(directory, identity, label, point_ids, utc_text, altitude, predictions,
                                    references, used, saturated, reasons, transition_intervals)
            sub = lambda v: np.asarray(v)[idx]          # noqa: E731
            record["pressure_products_rows"] = {"observation_rows": len(idx),
                                                "grid_rows_dropped": len(utc_text) - len(idx)}
            return original_wpp(directory, identity, label, [point_ids[i] for i in idx], [utc_text[i] for i in idx],
                                sub(altitude), {g: sub(v) for g, v in predictions.items()},
                                {g: sub(v) for g, v in references.items()}, {g: sub(v) for g, v in used.items()},
                                sub(saturated), {g: [reasons[g][i] for i in idx] for g in reasons},
                                transition_intervals)
        fit_cops.write_pressure_products = write_pressure_products
        try:
            return original_a5(*args, **kwargs)
        finally:
            fit_cops.write_pressure_products = original_wpp
    fit_cops.write_stage_a5_report = write_stage_a5_report


TEXT = """\
# CLIP_REFIT.txt -- written by model_ext/clip_problem.py ({utc})
# {pid}: {description}
# {policy}

H2O activity per facet: max(0, Y a + A_scale (c_cliff T_cliff + sum_k c_k T_k)) -- the fitter's clipped linear model;
templates {templates}; held (fixed_parameters --fix): {held}.
Base: {run} ({nrows} rows; {ngrid} model-only 1-s grid rows, reason {reason}, use 0). Fit data: {nng} NG, {nrg} RG.
Late window (>= {late} UTC, {nlate} observations): {late_model}
Optimizer start: {start}

VALID in this directory: the fitter's products (fit_result.txt, ActivityFit.txt, ActivityEnhancementFit.txt,
  ActivityEnhancementFacets.dat, activity_*.dat, activity_maps.png: clipped linear, as in c6), grain_term.*,
  fit_restart.npz (activity_linear = the linear field before the clip; meta model_extension = this marker),
  fit_manifest.json (model_extension: every argument; clip_refit: checks and records), start_check_model.npz.
  COPS_model.dat: observation rows = this model ({table_model}); GRID ROWS ({ngrid}, header line ClipRefitGridRows):
  model-only, pNG = INSTANTANEOUS surface model + grain(t), pRG = instantaneous surface model; NOT observations: drop
  them for plots and comparisons (figures_cluster/scripts/positive_table.py does, with late_window_model.dat), never
  re-average pNG (it includes the grain).
  late_window_model.dat: the late window row by row (kind grid | late_obs; surface-only model, grain, the table's values,
  data, fit flags; same columns as a positive fit's).
  pressure_residuals.*, pressure_comparison_*.png, figures/: observation rows only.
{replay}
"""


def write_text(out, P, cfg, record, fixed):
    lay = P.late_layout
    chain = record.get("fit_chain") or {}
    if chain.get("start_restart"):
        s = chain["start_restart"]
        start = "the final state of %s (sha256 %s; cost %.9e on its base, nfev %s), by name" % (
            s["path"], s["sha256"], s["cost"], s["nfev"])
    elif chain.get("resumed_from"):
        start = "resumed from %s (wall-time segments)" % chain["resumed_from"][-1]["path"]
    else:
        start = chain.get("start", "?")
    d2 = cfg["variant"] == "D2"
    replay = ("cops_fit.in (LiouvilleCOPS C++ replay) gives the INSTANTANEOUS model: for D2 it does not reproduce the late "
              "observations of this fit (W10x5); the C++ has no grain either.\nFine-tuning tools (grain_finetune.py, "
              "fine_tune_gpu.py) evaluate the instantaneous model: for D2 not the fitted problem." if d2 else
              "cops_fit.in replays the surface model (instantaneous, as fitted; the C++ has no grain).")
    text = TEXT.format(
        utc=utc_now(), pid=PARAMETERIZATION[cfg["variant"]], description=DESCRIPTION[cfg["variant"]], policy=POLICY,
        templates=", ".join(P.template_names) or "-", held=", ".join("%s=%r" % kv for kv in sorted(fixed.items())),
        run=str(P.M.run), nrows=len(lay.ipoint) if lay is not None else "?",
        ngrid=len(lay.grid_rows) if lay is not None else 0, reason=pp.GRID_REASON, nng=int(P.NNG), nrg=int(P.NRG),
        late=pp.LATE_START_UTC, nlate=len(lay.late_rows) if lay is not None else "?",
        late_model=("W10x5 of the surface model (mean of five 10-s averages ending at T, T-2, ..., T-8 s; T = L2 STOP), "
                    "grain added at the observation time stamp" if d2 else
                    "the instantaneous model at the observation time stamp (no time response), as c6"),
        table_model=("late observations W10x5(surface) + grain(t_obs), earlier ones instantaneous" if d2 else
                     "instantaneous surface model + grain"),
        start=start, replay=replay)
    path = Path(out) / "CLIP_REFIT.txt"
    with path.open("x") as f:
        f.write(text)
    return path


def write_late_window_products(out, fit_cops, r, record, variant):
    """late_window_model.dat (surface-only model, grain, kinds; positive_problem's columns) and the header notes of
    COPS_model.dat on the grid rows and the late observations."""
    out = Path(out)
    P, MF = fit_cops.P, fit_cops.MF
    lay = P.late_layout
    a, b, T, kappa = P.unpack(r["x"])
    pNG_s, pRG_s = MF.evaluate(a, T, kappa, b=b if P.nb else None, enhancement_basis=P.enhancement_basis,
                               enhancement_scale=P.enhancement_scale)   # surface model (D2: late observations W10x5)
    pos = {int(row): k for k, row in enumerate(np.asarray(MF.rows))}
    grain = P.grain(r["x"], P.et_all)[0]
    fit_pos = {int(row): k for k, row in enumerate(np.asarray(P.M.rows))}
    stop = {int(row): float(t) for row, t in zip(lay.late_rows, lay.stop_et)}
    table = out / "COPS_model.dat"
    lines = table.read_text().splitlines(True)
    k_hdr = next(i for i, line in enumerate(lines) if line.startswith("#iPoint"))
    names = lines[k_hdr][1:].split()
    c = {n: i for i, n in enumerate(names)}
    data = [line.split() for line in lines[k_hdr + 1:] if line.strip() and not line.startswith("#")]
    by_point = {int(v[c["iPoint"]]): v for v in data}
    rows = sorted([int(x) for x in lay.grid_rows] + [int(x) for x in lay.late_rows], key=lambda i: lay.et[i])
    late_text = ("W10x5 of the surface model ending at the L2 STOP time (stop_et)" if variant == "D2" else
                 "the INSTANTANEOUS surface model at et (variant D1: no time response; stop_et for reference)")
    worst = 0.0
    with (out / "late_window_model.dat").open("x") as f:
        f.write("# late window of this fit (model_ext/clip_problem.py, variant %s, %s): every 1-s grid row (model only) "
                "and every observation from %s UTC on\n" % (variant, utc_now(), pp.LATE_START_UTC))
        f.write("# kind grid: surface_* = instantaneous surface model at et; kind late_obs: surface_* = %s; "
                "grain_NG_Pa = grain(et) (NG only, unaveraged);\n# table_pNG/pRG = COPS_model.dat (grid: instantaneous "
                "surface + grain; late_obs: surface_* + grain); re-average only surface_* of grid rows, never table_pNG\n"
                % late_text)
        f.write("#iPoint UTC et kind stop_et surface_pNG_Pa surface_pRG_Pa grain_NG_Pa table_pNG_Pa table_pRG_Pa "
                "refNG_Pa refRG_Pa fitNG fitRG\n")
        for i in rows:
            k = pos[i]
            v = by_point[int(lay.ipoint[i])]
            kind = "late_obs" if i in stop else "grid"
            fk = fit_pos.get(i)
            fitng = int(P.useNG[fk]) if fk is not None else 0
            fitrg = int(P.useRG[fk]) if fk is not None else 0
            tng, trg = float(v[c["pNG"]]), float(v[c["pRG"]])
            worst = max(worst, abs(tng / (pNG_s[k] + grain[i]) - 1.0), abs(trg / pRG_s[k] - 1.0))
            f.write("%d %s %.6f %s %s %.9e %.9e %.9e %s %s %s %s %d %d\n" % (
                int(lay.ipoint[i]), lay.utc[i], lay.et[i], kind, "%.6f" % stop[i] if i in stop else "-",
                pNG_s[k], pRG_s[k], grain[i], v[c["pNG"]], v[c["pRG"]], v[c["refNG"]], v[c["refRG"]], fitng, fitrg))
    note = ["# ClipRefitGridRows (model_ext/clip_problem.py, variant %s): the %d rows of kind grid (1-s late grid of %s; "
            "reference 0, use 0, reason %s) are MODEL-ONLY points: pNG = instantaneous surface model + grain(t), "
            "pRG = instantaneous surface model. They are not observations: drop them in plots and comparisons and do "
            "not re-average them (pNG includes the grain); late_window_model.dat lists them (kind grid) with the "
            "surface-only model.\n" % (variant, len(lay.grid_rows), record.get("point_table", {}).get("path"),
                                       pp.GRID_REASON),
            "# ClipRefitLateObservations: the %d observations from %s UTC on have %s.\n" % (
                len(lay.late_rows), pp.LATE_START_UTC,
                "pNG = W10x5(surface) + grain(t_obs), pRG = W10x5(surface) (L2 STOP times in late_window_model.dat); "
                "earlier rows: instantaneous model" if variant == "D2" else
                "the instantaneous model (+ grain for NG), as every earlier row (variant D1: no time response)")]
    tmp = table.with_name(table.name + ".tmp.%d" % os.getpid())
    tmp.write_text("".join(lines[:k_hdr] + note + lines[k_hdr:]))
    os.replace(str(tmp), str(table))
    record["late_window_products"] = {
        "late_window_model.dat": {"sha256": sha256(out / "late_window_model.dat"), "rows": len(rows)},
        "COPS_model.dat_header_note": "ClipRefitGridRows / ClipRefitLateObservations lines before #iPoint",
        "table_vs_surface_plus_grain_max_rel": float(worst)}
    say("wrote %s (%d rows; COPS_model.dat = surface + grain to %.1e); grid rows flagged in COPS_model.dat"
        % (out / "late_window_model.dat", len(rows), worst))


def install_restart_marker(record, cfg):
    from cops_cuda import restart as rs
    if getattr(rs.from_fit, "_clip_marker", False):
        return
    original = rs.from_fit

    def from_fit(problem, result, fit_dir, manifest, elements):
        state = original(problem, result, fit_dir, manifest, elements)
        lw = record.get("late_window") or {}
        state.meta["model_extension"] = {
            "parameterization": PARAMETERIZATION[cfg["variant"]], "variant": cfg["variant"],
            "description": DESCRIPTION[cfg["variant"]], "policy": POLICY,
            "point_table": record.get("point_table"), "fit_chain": record.get("fit_chain"),
            "fixed": cfg.get("fixed"),
            "late_window": {k: lw.get(k) for k in ("late_start_utc", "late_observations", "grid_rows",
                                                   "window_installed")},
            "activity": "clipped linear: max(0, activity_linear) (the fitter's semantics)",
            "warning": ("D2: the late observations are W10x5 of the surface model; tools that evaluate the instantaneous "
                        "model (fine tuning, cops_fit.in) do not evaluate this fit's problem" if cfg["variant"] == "D2"
                        else None)}
        return state
    from_fit._clip_marker = True
    rs.from_fit = from_fit


# ------------------------------------------------------------------------------------------------ installation
def install_core(fit_cops, cfg, record, source=None, arrays=None):
    """Before grain_crossing.install: grid-row selection, product filters, the window core (D2), checks."""
    variant = cfg["variant"]
    record.update(variant=variant, parameterization=PARAMETERIZATION[variant], description=DESCRIPTION[variant],
                  policy=POLICY)
    table = gkeys = None
    if cfg.get("point_table"):
        table = pp.read_point_table(cfg["point_table"])
        gkeys = pp.grid_keys(table)
        record["point_table"] = {"path": cfg["point_table"], "sha256": sha256(cfg["point_table"]), "rows": len(table),
                                 "grid_rows": len(gkeys)}
        pp.install_selection(fit_cops, gkeys, record)          # record["selection"] when the base has grid rows
        install_pressure_filter(fit_cops, gkeys, record)
        pp.install_figure_filter(gkeys, record)
    base = fit_cops.Problem

    class ClipCore(base):
        def __init__(self, M, el, coef, args):
            self.late_layout = None
            lw = None
            if table is not None:
                points = pp.read_run_points(M.run)
                layout = pp.LateLayout(table, points, M.ipoint)
                grid = layout.grid_rows
                if np.any(M.refNG[grid] != 0.0) or np.any(M.refRG[grid] != 0.0):
                    raise SystemExit("clip_refit: grid rows with reference data")
                if np.any(M.useNG[grid] != 0) or np.any(M.useRG[grid] != 0):
                    raise SystemExit("clip_refit: grid rows with use flags after the selection")
                lw = dict(layout.summary, window_installed=False)
                if variant == "D2":
                    # as positive_problem.PositiveCore: grid rows join the fit rows; W10x5 on M and on MF
                    rows = np.union1d(np.asarray(M.rows, dtype=np.int64), grid)
                    if len(rows) != len(M.rows) or np.any(rows != np.asarray(M.rows)):
                        M.load_kernels(rows)
                    M.sel = rows
                    op_fit = pp.install_window(M, layout)
                    MF = getattr(fit_cops, "MF", None)
                    op_all = pp.install_window(MF, layout) if MF is not None and MF is not M else None
                    shared = MF is not None and len(MF.rows) == len(M.rows) and bool(np.all(MF.rows == M.rows))
                    lw.update(window_installed=True, late_rows_in_fit_model=int(len(op_fit.late_pos)),
                              installed_on_all_rows_model=op_all is not None, fit_rows_equal_all_rows=shared)
                    if not shared:
                        say("WARNING the fit rows (%d) differ from all rows (%d): a second GPU backend is bound"
                            % (len(M.rows), len(MF.rows) if MF is not None else -1))
                else:
                    in_fit = np.intersect1d(np.asarray(M.rows), grid)
                    if len(in_fit):
                        raise SystemExit("clip_refit: D1: %d grid rows in the fit model" % len(in_fit))
                    MF = getattr(fit_cops, "MF", None)
                    lw.update(grid_rows_in_fit_model=0, fit_model_rows=int(len(M.rows)),
                              all_rows_model_rows=int(len(MF.rows)) if MF is not None else None,
                              note="D1: the fit model M has the fitter's own rows (use > 0; no grid row); the all-rows "
                                   "model MF (its own GPU backend) writes the grid rows")
                self.late_layout = layout
            elif variant == "D2":
                raise SystemExit("clip_refit: D2 needs the point table")
            super().__init__(M, el, coef, args)
            if self.late_layout is not None:
                pos = {int(r_): k for k, r_ in enumerate(np.asarray(M.rows))}
                gpos = [pos[int(r_)] for r_ in self.late_layout.grid_rows if int(r_) in pos]
                used = int(np.sum(self.useNG[gpos]) + np.sum(self.useRG[gpos])) if gpos else 0
                if used:
                    raise SystemExit("clip_refit: %d grid-row data entries in the fit" % used)
                lw.update(fit_rows=int(len(M.rows)), grid_rows_in_fit_rows=len(gpos), grid_data_entries_used=used,
                          NNG=int(self.NNG), NRG=int(self.NRG))
                record["late_window"] = lw

    fit_cops.Problem = ClipCore
    if cfg.get("check_full_jacobian"):
        install_jacobian_check(fit_cops, cfg["check_parameters"], record)
    install_start_check(fit_cops, cfg, record, source, arrays)    # outside the Jacobian check: runs first
    return ClipCore


def install_outer(fit_cops, cfg, record, source=None):
    """After fixed_parameters.install: start restart (by name), late-window/grain checks, products, restart marker."""
    base = fit_cops.Problem
    variant = cfg["variant"]

    class ClipOuter(base):
        def __init__(self, M, el, coef, args):
            super().__init__(M, el, coef, args)
            names = list(self.names)
            if tuple(names[-2:]) != GRAIN_NAMES:
                raise SystemExit("clip_refit: unexpected parameter layout (grain last): %s" % names[-4:])
            et_u = pp.read_run_points(M.run)
            ip = {i: int(p) for i, p in enumerate(np.asarray(M.ipoint))}
            i_end = int(np.nonzero(self.et_fit == self.et_end)[0][0])
            row_end = int(np.asarray(M.rows)[i_end])
            # crossing pairs (grain_crossing.pairs' definition), for the record
            both = self.useNG & self.useRG
            et_b = self.et_fit[both]
            gw = {"t_end_et": float(self.et_end), "t_end_utc": et_u[ip[row_end]][0],
                  "last_ng_before_gap_et": float(self.et_before), "gap_seconds": float(self.et_end - self.et_before),
                  "pairs_before_gap": int(np.sum(et_b <= self.et_before)),
                  "pairs_after_gap": int(np.sum(et_b >= self.et_end)), "NNG": int(self.NNG), "NRG": int(self.NRG)}
            if source is not None:
                for k in ("t_end_et", "gap_seconds", "pairs_before_gap", "pairs_after_gap"):
                    if source.get(k) is not None:
                        gw["source_" + k] = source[k]
            record["grain_window"] = gw
            lay = self.late_layout
            if lay is not None:
                equal = row_end == int(lay.late_rows[0])
                record.setdefault("late_window", {})["grain_t_end_row"] = row_end
                record["late_window"]["grain_t_end_equals_late_start"] = bool(equal)
                if variant == "D2" and not equal:
                    raise SystemExit("clip_refit: the first post-gap NG point (row %d) is not the first late-window "
                                     "observation (row %d, %s): the window start needs a user decision"
                                     % (row_end, int(lay.late_rows[0]), pp.LATE_START_UTC))
            say("grain window: t_end %s (et %.3f), gap %.1f s, crossing pairs %d before / %d after the gap; %d NG, %d RG "
                "points%s" % (gw["t_end_utc"], gw["t_end_et"], gw["gap_seconds"], gw["pairs_before_gap"],
                             gw["pairs_after_gap"], gw["NNG"], gw["NRG"],
                             "; source: t_end et %s, gap %s s, pairs %s / %s" % (
                                 source.get("t_end_et"), source.get("gap_seconds"), source.get("pairs_before_gap"),
                                 source.get("pairs_after_gap")) if source is not None else ""))
            if cfg.get("start_restart"):
                self.restart_start = start_from_restart(self, el, cfg["start_restart"], record)

    fit_cops.Problem = ClipOuter
    original_write = fit_cops.write_fit

    def write_fit(out, r, *rest):
        diag = original_write(out, r, *rest)
        P = fit_cops.P
        if P.late_layout is not None:
            write_late_window_products(out, fit_cops, r, record, variant)
        path = write_text(out, P, cfg, record, cfg.get("fixed") or {})
        record.setdefault("products", {})["CLIP_REFIT.txt"] = sha256(path)
        return diag
    fit_cops.write_fit = write_fit
    install_restart_marker(record, cfg)
    return ClipOuter


def install_hooks(grain_crossing, fixed_parameters, cfg, command, sub, source=None, arrays=None):
    """Hook this module into the unchanged grain_fit_fixed.py path (clip_refit.py and the CPU tests): install_core
    runs just before grain_crossing.install, install_outer just after fixed_parameters.install (the last step of
    grain_fit_fixed's install); the manifest's model_extension.command becomes ``command`` (grain_fit_fixed's own
    command is kept as grain_fit_fixed_command)."""
    original_gc = grain_crossing.install
    original_fp = fixed_parameters.install
    main_record = {}

    def gc_install(fit_cops, weight, record):
        main_record["record"] = record
        record["clip_refit"] = sub
        install_core(fit_cops, cfg, sub, source, arrays)
        return original_gc(fit_cops, weight, record)

    def fp_install(fit_cops, fixed, fp_record):
        result = original_fp(fit_cops, fixed, fp_record)
        install_outer(fit_cops, cfg, sub, source)
        record = main_record["record"]
        record["grain_fit_fixed_command"] = record.get("command")
        record["command"] = list(command)
        return result
    grain_crossing.install = gc_install
    fixed_parameters.install = fp_install
    return original_gc, original_fp
