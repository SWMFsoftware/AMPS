# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f6d8d792ccc8c46c9570909b7146f6b9d9dae940d62a4f84f8cea74f757f00b5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Fit driver of the positive-activity fits (user decisions 2026-11-07): wraps fit_cops.Problem / gpu_fit.py; activity
#   exp(Y a + c T_cliff), c >= 0; start and prior = emission-weighted (w = area A_orig) LS projection of ln A_orig;
#   area-weighted log Tikhonov D; late observations (>= 10:33:10 UTC) compared with W10x5 of the surface model (+ grain
#   at the time stamp); model-only 1-s grid rows (use 0); inserts --check-tolerance inf, --template-start 0,
#   --fix-temperature; writes activity_exp_facets.npz, ACTIVITY_PARAMETERIZATION.txt, late_window_model.dat.
# USED BY: positive_fit.py.
# ==== end of header ====
"""Positive H2O activity A = exp(s), start from the original map, log-space Tikhonov term and the COPS time response
in the late window (user decisions 2026-11-07). Fit-level part; the GPU side is model_ext/positive_activity.py.

Local, untracked extension of the unchanged fitter (``fit_cops``) and CUDA driver (``tools/gpu_fit.py``), installed by
``positive_fit.py`` around ``grain_crossing.install`` (like grain_fit_weighted.py). Nothing in the git tree is edited.

Model (user 2026-11-07)
* Activity positive by construction: A_e = exp(s_e), s = Y a + 1.0 * c_cliff * T_cliff (exponent field: the 25
  degree-<=4 real spherical harmonics Y of the fitter plus the template cliff_ge_60deg of
  runs/shap8_012m_flat_earlier/bundle/restart_templates.npz inside the exponent, i.e. a multiplicative factor
  exp(c_cliff) on cliffs). The fine-tuning templates ft_c1..ft_c6 of earlier fits are not used. The GPU patch
  (positive_activity.py, env COPS_POS_*) forms af = exp(act_lin) from act_lin = Y a + enhancement_scale * B c exactly as
  the fitter passes them, so this Problem passes enhancement_scale = 1.0 and exponent coefficients. The same patch
  applies the night/shadow ratio 0.02 and the per-record flux floor 0.02 q_max A_orig (COPS_POS_NIGHT_RATIO,
  COPS_POS_FLOOR_RATIO); this module does not reimplement them.
* Optimizer coordinates: a_j = ascale * x_j with ascale = max_j |a_start_j| (O(1) start; replaces the fitter's
  max |coef[0]|, which is the scale of the linear srcCG coefficients); c_cliff = x (amplitude, lower bound 0 as for
  every template, CLAUDE.md; with exp it is ln of the cliff factor).
* Start and prior (user 2026-11-07: "the ORIGINAL map"; weighting decided 2026-11-07 after the 12M check: EMISSION-
  WEIGHTED): least-squares projection of ln(A_orig) onto [Y (25 columns), T_cliff],
      min sum_e w_e (s_e - ln A_orig,e)^2 over the facets with A_orig > 0,  w_e = area_e * A_orig,e
  (A_orig = the run's H2O activity, elements.bin column 32 = max(0, Y coef[0]); 0.7 % of the 12M area has A_orig = 0
  and is excluded; the weight is the original emission, so low-activity facets near the zero lines of the linear map,
  where ln A_orig diverges, do not dominate). If the projected cliff coefficient is negative it is set to 0 and Y alone
  is re-projected (stated in the record). The first version (area weights w_e = area_e, backup
  positive_problem.py.before_emission_weighting_20261107T073050Z) gave 41 % of the area / 32 % of the emission off by
  more than x2 on 12M. Reported: area- and emission-weighted rms of the ln A residual over the A_orig > 0 facets, the
  area and emission fractions where exp(s0) differs from A_orig by more than x2, the area fraction with A_orig = 0,
  total activity of exp(s0) relative to A_orig. Recorded as "start_projection_weighting": "emission_weighted_log_ls".
* Temperatures: the original analytical law held (Tmax 182.09922, T75 136.11149, Tmin 133 K = the run's M.T0) with the
  fitter's own --fix-temperature (checked against these values). kappa free >= 0, start 1 (fitter default).
* Regularization (lambda scanned later, one --lambda per process): log space,
      D^2 = sum_e area_e (s_e - s0_e)^2 / sum_e area_e   over all facets (s0 = the start/prior exponent field),
  i.e. the rms deviation of ln A from the prior s0 = the emission-weighted projection of the original map (area
  weights in the penalty, emission weights in the projection: D is therefore NOT, up to a constant, the area-weighted
  deviation of ln A from ln A_orig). Implemented as the fitter's Lt with
  G = B^T diag(area) B / sum(area), B = [Y * ascale, T_cliff * 1.0] (columns = derivatives with respect to the optimizer
  coordinates), eigen-square-root like the reference (Lt^T Lt = G); prior = the start x0; reg_slice = the activity
  coordinates. D is the rms deviation of ln A in nepers.
* COPS time response (user 2026-11-07: averaging window "10 s + 5 readings", W10x5): every observation from
  2016-09-30T10:33:10 UTC on (the late window = the first post-gap NG point) is compared with the mean of five 10-s
  running averages of the model ending at T, T-2, T-4, T-6, T-8 s, T = the observation's L2 STOP time (point table of
  late_time_grid.py: L2_stop). The model at the 1-s grid points of the base (kind "grid", reference values 0, use flags
  0) is interpolated linearly in time (0.05-s sub-samples) and integrated with the trapezoidal rule, exactly the
  definitions of late_response_eval.py (window_mean/responses, copied below), as a fixed sparse matrix
  W (n_late_obs x n_grid). The grid rows are added to the fit model's rows; the model's evaluate is wrapped (a
  per-instance subclass, see install_window) so that pNG, pRG and every Jacobian column have their late-observation
  entries replaced by W @ (grid entries). Observations before 10:33:10 keep the instantaneous model. The same averaging
  is installed on the all-rows model MF, so COPS_model.dat reports the averaged surface model (+ grain) for the late
  observations and the instantaneous model elsewhere. The dust grain (grain_crossing.py) is added afterwards at the
  observation's time stamp, unaveraged: the averaging applies to the surface model only.
  Grid rows in the written products (review 2026-11-07): COPS_model.dat keeps them (instantaneous surface + grain,
  model-only) with a header note (PositiveFitGridRows) and late_window_model.dat lists the late window row by row
  (kind grid | late_obs, surface-only model, grain); the A5 pressure products (pressure_residuals.*,
  pressure_comparison_*) and the gpu_fit figures (cops_cuda.figures.read_table wrapper) contain the observation rows
  only. Downstream scripts must not re-average the table's grid rows (their pNG includes the grain).
* Template: exactly one, cliff_ge_60deg (user 2026-11-07; refused otherwise).
* Grain + crossing (grain_crossing.py) as in all grain fits (NG = surface + grain(t), crossing weight 10). Grain start
  A = 6.96024299190473e-06 Pa (cycle 3's fine-tuned value, grain/cycles/c3/finetune/.../grain_finetune.json "A_Pa"),
  B = 0 (user 2026-11-07); grain_crossing.py's own start (the observed jump) is recorded and replaced.
* Selection: the base has the observation rows plus the grid rows. The named profile (all_available) audits the
  1150-row mission dataset, so the profile is applied to the observation rows only (its audited counts unchanged) and
  the grid rows get use = 0 with the reason code "late_grid_point" (apply_selection_profile wrapper).
* Outputs: the fitter's activity writers assume max(0, linear). Patched (content in exp semantics, same file layout):
  fit_cops.activity_diagnostics (the "activity diagnostics" lines of fit_result.txt and the manifest; its fixed header
  text still names the clip) and the A5 activity products (activity_facets.dat: legacy_linear = Y coef[0] of the run,
  residual_linear = A - legacy_linear, total_linear = clipped_activity = A = exp(s); activity_template_power.dat:
  the activity added by the template factor, A - exp(Y a); activity_degree_power.dat: no rows, A has no linear degree
  decomposition; activity_region_diagnostics.dat, activity_maps.png). Not replayable and NOT to be used:
  ActivityFit.txt (exponent coefficients), ActivityEnhancementFit.txt, ActivityEnhancementFacets.dat, cops_fit.in
  (LiouvilleCOPS C++ applies max(0, linear)), fit_restart.npz's activity_linear (= s; meta "model_extension" marks
  it). Added: activity_exp_facets.npz (s, A, s0, F_min/q_max = floor ratio x A_orig, seen) and
  ACTIVITY_PARAMETERIZATION.txt.
* The fitter's start-up consistency check (MF.evaluate at the run's linear coefficients vs the LiouvilleCOPS
  pressures) cannot pass with this physics (exp of the linear srcCG coefficients overflows); positive_fit.py inserts
  --check-tolerance inf (inf, not 1e300: an inf difference would still fail "> 1e300"; nan never compares) and the
  manifest records that it is replaced by the extension's validation (the GPU patch's checks and --check-jacobian).
"""

import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import scipy.sparse as sp

# ------------------------------------------------------------------------------------------------ user decisions
USER_T0 = (182.09922, 136.11149, 133.0)            # K, the run's analytical law (held, user 2026-11-07)
LATE_START_UTC = "2016-09-30T10:33:10"             # first post-gap NG point; W10x5 from here on (user 2026-11-07)
WINDOW_WIDTH_S = 10.0                              # electrometer 10-s running average
WINDOW_LAGS_S = (0.0, 2.0, 4.0, 6.0, 8.0)          # DPU mean of the last 5 readouts (2-s readout)
WINDOW_DT_S = 0.05                                 # sub-sample step of late_response_eval.window_mean
DEFAULT_GRAIN_START_UPA = 6.96024299190473         # 6.96024299190473e-06 Pa: cycle 3 fine-tuned A_Pa (user 2026-11-07)
POSITIVE_ENV = (("COPS_POS_ENABLE", "1"), ("COPS_POS_NIGHT_RATIO", "0.02"), ("COPS_POS_ACTIVITY", "exp"),
                ("COPS_POS_FLOOR_RATIO", "0.02"))  # agent B's patch (user 2026-11-07)
PARAMETERIZATION_ID = "positive_exp_activity_log_tikhonov_w10x5_emisproj_v2"   # v2: emission-weighted start/prior
PROJECTION_WEIGHTING = "emission_weighted_log_ls"
PROJECTION_WEIGHT_DEFINITION = ("w_e = area_e * A_orig,e (A_orig > 0; normalized by max A_orig, which does not change "
                                "the solution); min sum_e w_e (s_e - ln A_orig,e)^2")
CHECK_PARAMETERS = ("a00", "a02", "a07", "c_cliff_ge_60deg", "kappa", "grain_A_uPa")
GRID_REASON = "late_grid_point"
REQUIRED_TEMPLATE = "cliff_ge_60deg"                # the user's model (2026-11-07): exactly this template
POLICY = ("user 2026-11-07: A = exp(Y a + c_cliff T_cliff) (positive by construction; templates inside the exponent; "
          "ft_c* not used); start and prior = emission-weighted (w = area A_orig) LS projection of ln A_orig on "
          "[Y, T_cliff] over A_orig > 0 (negative "
          "cliff coefficient -> 0, Y re-projected); original analytical temperature law held (fit_cops "
          "--fix-temperature); kappa free >= 0 from 1; night ratio 0.02 and flux floor 0.02 q_max A_orig on the GPU "
          "(positive_activity.py); D^2 = sum area (s - s0)^2 / sum area (log space, prior = start); late window "
          "(>= 10:33:10 UTC) compared with W10x5 of the surface model (mean of five 10-s averages ending at "
          "T, T-2, ..., T-8 s, T = L2 STOP), grain added at the observation time; grain + crossing 10, grain start "
          "A = 6.96024299190473e-06 Pa, B = 0")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def jsonable(value):
    """Plain JSON types (the fitter's manifest writer has no default=)."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (str, int)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    return repr(value)


# ------------------------------------------------------------------------------------------------ environment (B)
def set_positive_environment(environ=None):
    """Set the COPS_POS_* variables of the GPU patch (before anything binds a backend; spawned workers inherit them).
    A different preset value is refused: the values are the user's decisions of 2026-11-07."""
    environ = os.environ if environ is None else environ
    for key, value in POSITIVE_ENV:
        old = environ.get(key)
        if old is not None and old != value:
            raise SystemExit("positive_fit: %s=%s is preset; the user's decision (2026-11-07) is %s" % (key, old, value))
        environ[key] = value
    return dict(POSITIVE_ENV)


def describe_positive_activity(module):
    """Provenance of the GPU patch module: file, SHA-256, config() (settings, kernel copies and their SHA-256),
    kernel_sha256(), modifications_active() and the SHA-256 of every file in model_ext/kernels/."""
    path = Path(module.__file__).resolve()
    out = {"module": str(path), "module_sha256": sha256(path)}
    try:
        out["config"] = jsonable(module.config())
    except Exception as exc:                      # recorded, never silent
        out["config_error"] = "%s: %s" % (exc.__class__.__name__, exc)
    for name in ("kernel_sha256", "modifications_active"):
        f = getattr(module, name, None)
        if callable(f):
            try:
                out[name] = jsonable(f())
            except Exception as exc:
                out[name + "_error"] = "%s: %s" % (exc.__class__.__name__, exc)
    kernels = path.parent / "kernels"
    if kernels.is_dir():
        out["kernel_files_sha256"] = {str(p): sha256(p) for p in sorted(kernels.iterdir()) if p.is_file()}
    return out


# ------------------------------------------------------------------------------------------------ arguments
def _pop(head, name, convert=str, default=None):
    if name in head:
        i = head.index(name)
        if i + 1 >= len(head):
            raise SystemExit("positive_fit: %s needs a value" % name)
        value = convert(head[i + 1])
        del head[i:i + 2]
        return value
    return default


def _flag(head, name):
    if name in head:
        head.remove(name)
        return True
    return False


def _value(tail, name):
    """Value of a single-valued fitter option (or None)."""
    for k, token in enumerate(tail):
        if token == name and k + 1 < len(tail):
            return tail[k + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def rewrite_arguments(argv):
    """Split the positive_fit.py command into (driver config, grain_fit.py argv, notes).

    Driver options (before '--'): --point-table PATH (required), --grain-start-uPa A (default 6.96024299190473),
    --check-jacobian, --check-parameters a00,a02,... . Everything else goes to grain_fit.py / gpu_fit.py
    (--crossing-weight, --devices, --precision, --max-nfev, --wall-time, ...).
    Inserted when absent: gpu_fit --template-start 0 (the start is the projection; the default 1.0 would replace the
    projected cliff coefficient), fitter --fix-temperature; fitter --check-tolerance is set to inf. Refused: options
    that contradict the user's model (see the messages)."""
    argv = list(argv)
    if "--" not in argv:
        raise SystemExit("positive_fit: need '--' before the fit_cops arguments")
    split = argv.index("--")
    head, tail = argv[:split], argv[split + 1:]
    notes = []
    cfg = {"point_table": _pop(head, "--point-table"),
           "grain_start_uPa": _pop(head, "--grain-start-uPa", float, DEFAULT_GRAIN_START_UPA),
           "check_jacobian": _flag(head, "--check-jacobian"),
           "check_parameters": list(CHECK_PARAMETERS)}
    names = _pop(head, "--check-parameters")
    if names:
        cfg["check_parameters"] = [n.strip() for n in names.split(",") if n.strip()]
    if not cfg["point_table"]:
        raise SystemExit("positive_fit: --point-table <run>/late_grid_points.txt is required")
    cfg["point_table"] = str(Path(cfg["point_table"]).resolve())
    if not Path(cfg["point_table"]).is_file():
        raise SystemExit("positive_fit: point table not found: %s" % cfg["point_table"])
    if not cfg["grain_start_uPa"] >= 0.0:
        raise SystemExit("positive_fit: --grain-start-uPa must be >= 0")
    # gpu_fit.py options
    ts = _value(head, "--template-start")
    if ts is None:
        head += ["--template-start", "0"]
        notes.append("inserted gpu_fit --template-start 0 (start = the projection of the original map)")
    elif float(ts) != 0.0:
        raise SystemExit("positive_fit: --template-start must be 0: the template amplitude starts at the projected "
                         "value (user 2026-11-07); a nonzero value would replace it")
    cfg["restart_from"] = _value(head, "--restart-from")
    # fitter options
    for bad, why in (("--ordered-temperature", "the temperature law is held (--fix-temperature)"),
                     ("--fix-legacy-activity", "all 25 exponent coefficients are fitted"),
                     ("--fix-kappa", "kappa is free (user 2026-11-07)"),
                     ("--restart-from", "give --restart-from to gpu_fit (before '--')"),
                     ("--check-jacobian", "the fitter's own check assumes the linear clip; use positive_fit "
                                          "--check-jacobian (before '--')"),
                     ("--enhancement-scale", "the exponent enters with enhancement scale 1"),
                     ("--data-only", "not a fit")):
        if bad in tail or any(t.startswith(bad + "=") for t in tail):
            raise SystemExit("positive_fit: %s is not allowed: %s" % (bad, why))
    for name, need in (("--degree", "4"), ("--enhancement-degree", "4")):
        v = _value(tail, name)
        if v is not None and int(v) != int(need):
            raise SystemExit("positive_fit: %s must be %s (exponent field = the 25 degree-<=4 harmonics + template)"
                             % (name, need))
    templates = [tail[k + 2] for k, t in enumerate(tail) if t == "--activity-template" and k + 2 < len(tail)]
    if templates != [REQUIRED_TEMPLATE]:
        raise SystemExit("positive_fit: exactly one --activity-template <bundle> %s is required (the user's model, "
                         "2026-11-07); got %s" % (REQUIRED_TEMPLATE, templates))
    cfg["templates"] = templates
    if "--fix-temperature" not in tail:
        tail.append("--fix-temperature")
        notes.append("inserted fit_cops --fix-temperature (original temperature law held)")
    k = None
    for i, t in enumerate(tail):
        if t == "--check-tolerance":
            k = i
    if k is not None:
        if tail[k + 1] != "inf":
            notes.append("replaced fit_cops --check-tolerance %s by inf" % tail[k + 1])
        tail[k + 1] = "inf"
    else:
        tail += ["--check-tolerance", "inf"]
        notes.append("inserted fit_cops --check-tolerance inf (consistency check replaced, see the record)")
    run = next((t for t in tail if not t.startswith("-")), None)
    if run is None:
        raise SystemExit("positive_fit: no run directory after '--'")
    cfg["run"] = run
    return cfg, head + ["--"] + tail, notes


def check_restart_marker(path):
    """A --restart-from file must come from this parameterization (a linear fit's restart would be misread)."""
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    ext = meta.get("model_extension") or {}
    if ext.get("parameterization") != PARAMETERIZATION_ID:
        raise SystemExit("positive_fit: %s is not a restart of a positive_fit.py fit (meta model_extension."
                         "parameterization %r, need %r)" % (path, ext.get("parameterization"), PARAMETERIZATION_ID))
    return {"path": str(Path(path).resolve()), "parameterization": ext.get("parameterization")}


# ------------------------------------------------------------------------------------------------ tables
def parse_utc_seconds(value):
    """Seconds after 2016-01-01 for 'YYYY-MM-DDTh:mm:ss[.sss]' (hours without leading zero allowed; as
    cops_selection.parse_utc)."""
    date, clock = value.strip().rstrip("Z").split("T")
    hour, minute, second = clock.split(":")
    midnight = dt.datetime.strptime(date, "%Y-%m-%d")
    return (midnight - dt.datetime(2016, 1, 1)).total_seconds() + 3600 * int(hour) + 60 * int(minute) + float(second)


def time_key(value):
    return int(round(parse_utc_seconds(value) * 1000.0))       # ms


def read_point_table(path):
    """late_time_grid.py --table: index UTC kind original_index L2_start L2_stop ('-': none)."""
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) != 6 or f[2] not in ("grid", "obs"):
            raise SystemExit("point table %s: bad line %r" % (path, line.rstrip()))
        out.append({"index": int(f[0]), "utc": f[1], "kind": f[2], "original_index": int(f[3]),
                    "l2_start": f[4], "l2_stop": f[5], "key": time_key(f[1])})
    keys = [e["key"] for e in out]
    if len(set(keys)) != len(keys):
        raise SystemExit("point table %s has duplicate times" % path)
    return out


def grid_keys(table):
    return {e["key"] for e in table if e["kind"] == "grid"}


def read_run_points(run):
    """{iPoint: (UTC, et, refNG, refRG)} from <run>/COPS_model.dat."""
    names, out = None, {}
    with open(Path(run) / "COPS_model.dat", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#iPoint"):
                names = line[1:].split()
                continue
            if line.startswith("#") or not line.strip() or names is None:
                continue
            v = line.split()
            c = {n: i for i, n in enumerate(names)}
            out[int(v[c["iPoint"]])] = (v[c["UTC"]], float(v[c["et"]]), float(v[c["refNG"]]), float(v[c["refRG"]]))
    return out


# ------------------------------------------------------------------------------------------------ W10x5
# Copies of late_response_eval.py (2026-11-07) window_mean / responses: the definition the operator reproduces.
def window_mean(tg, yg, t_end, width, dt=WINDOW_DT_S):
    s = np.arange(t_end - width, t_end + 0.5 * dt, dt)
    s[-1] = t_end
    if s[0] < tg[0] - 1e-6 or s[-1] > tg[-1] + 1e-6:
        return np.nan
    return np.trapz(np.interp(s, tg, yg), s) / (s[-1] - s[0])


def responses(tg, yg, T):
    w10 = lambda e: window_mean(tg, yg, e, 10.0)      # noqa: E731
    return {"W4": window_mean(tg, yg, T, 4.0), "W10": w10(T), "W10x5": np.mean([w10(T - 2.0 * k) for k in range(5)])}


def window_weights(tg, t_end, width, dt=WINDOW_DT_S):
    """(grid indices, weights) with sum_j w_j y_j == window_mean(tg, y, t_end, width) for every y."""
    s = np.arange(t_end - width, t_end + 0.5 * dt, dt)
    s[-1] = t_end
    if s[0] < tg[0] - 1e-6 or s[-1] > tg[-1] + 1e-6:
        raise SystemExit("W10x5: the window [%.3f, %.3f] (et) is not covered by the grid [%.3f, %.3f]"
                         % (s[0], s[-1], tg[0], tg[-1]))
    sc = np.clip(s, tg[0], tg[-1])                      # np.interp clamps at the ends
    j = np.clip(np.searchsorted(tg, sc, side="right") - 1, 0, len(tg) - 2)
    f = (sc - tg[j]) / (tg[j + 1] - tg[j])
    d = np.diff(s)
    c = np.zeros(len(s))
    c[:-1] += 0.5 * d                                   # np.trapz
    c[1:] += 0.5 * d
    c /= s[-1] - s[0]
    return np.concatenate([j, j + 1]), np.concatenate([c * (1.0 - f), c * f])


def w10x5_operator(tg, stops, width=WINDOW_WIDTH_S, lags=WINDOW_LAGS_S):
    """Sparse W (len(stops) x len(tg)): row k = W10x5 ending at stops[k] (mean of the windows ending at T - lag)."""
    tg = np.asarray(tg, float)
    if len(tg) < 2 or not np.all(np.diff(tg) > 0):
        raise SystemExit("W10x5: the grid times must increase")
    rows, cols, vals = [], [], []
    for k, T in enumerate(stops):
        for lag in lags:
            j, w = window_weights(tg, float(T) - lag, width)
            rows.append(np.full(len(j), k))
            cols.append(j)
            vals.append(w / len(lags))
    W = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(stops), len(tg))).tocsr()
    W.sum_duplicates()
    return W


class LateLayout:
    """The late-window rows of a run and the W10x5 operator over its grid rows (run row indices)."""

    def __init__(self, table, run_points, ipoint, late_start_utc=LATE_START_UTC):
        self.table_path = None
        by_key = {e["key"]: e for e in table}
        t_late = time_key(late_start_utc)
        utc = [run_points[int(i)][0] for i in ipoint]
        et = np.array([run_points[int(i)][1] for i in ipoint])
        keys = [time_key(u) for u in utc]
        if len(set(keys)) != len(keys):
            raise SystemExit("late window: the run has duplicate times")
        matched = set()
        grid, late, stop = [], [], []
        index_equal = True
        for row, (key, u) in enumerate(zip(keys, utc)):
            e = by_key.get(key)
            if e is not None:
                matched.add(key)
                index_equal &= e["index"] == row
            if e is not None and e["kind"] == "grid":
                grid.append(row)
                continue
            if key >= t_late:
                if e is None or e["l2_stop"] == "-":
                    raise SystemExit("late window: observation %s (row %d) has no L2_stop in the point table" % (u, row))
                late.append(row)
                stop.append(et[row] + (parse_utc_seconds(e["l2_stop"]) - parse_utc_seconds(u)))
        missing = [e["utc"] for e in table if e["key"] not in matched]
        if missing:
            raise SystemExit("late window: %d points of the table are not rows of the run (e.g. %s)"
                             % (len(missing), missing[0]))
        if not grid or not late:
            raise SystemExit("late window: %d grid rows, %d late observations" % (len(grid), len(late)))
        grid = np.array(grid, dtype=np.int64)
        order = np.argsort(et[grid], kind="stable")
        self.grid_rows = grid[order]
        self.late_rows = np.array(late, dtype=np.int64)
        self.stop_et = np.array(stop)
        self.grid_et = et[self.grid_rows]
        self.late_et = et[self.late_rows]
        self.W = w10x5_operator(self.grid_et, self.stop_et)
        self.utc = list(utc)
        self.ipoint = np.asarray(ipoint, dtype=np.int64)
        self.et = et
        self.summary = {
            "late_start_utc": late_start_utc, "late_observations": int(len(self.late_rows)),
            "first_late_utc": utc[int(self.late_rows[0])], "last_late_utc": utc[int(self.late_rows[-1])],
            "grid_rows": int(len(self.grid_rows)), "grid_first_utc": utc[int(self.grid_rows[0])],
            "grid_last_utc": utc[int(self.grid_rows[-1])],
            "grid_step_s": [float(np.min(np.diff(self.grid_et))), float(np.max(np.diff(self.grid_et)))],
            "stop_minus_time_s": [float(np.min(self.stop_et - self.late_et)), float(np.max(self.stop_et - self.late_et))],
            "table_rows": len(table), "run_rows": len(ipoint), "table_index_equals_row": bool(index_equal),
            "operator": {"shape": list(self.W.shape), "nnz": int(self.W.nnz),
                         "row_sum_range": [float(self.W.sum(axis=1).min()), float(self.W.sum(axis=1).max())],
                         "definition": "W10x5 = mean_k W10(T - 2k), k = 0..4; W10(e) = trapz(interp(s, t_grid, y), s) "
                                       "/ 10 s, s = e - 10 s .. e step 0.05 s (late_response_eval.window_mean)"}}

    def operator_for(self, rows):
        """(late positions, grid positions, W restricted to the late rows present) within a model's rows."""
        pos = {int(r): k for k, r in enumerate(np.asarray(rows))}
        missing = [int(r) for r in self.grid_rows if int(r) not in pos]
        if missing:
            raise SystemExit("late window: %d grid rows are not loaded (e.g. row %d)" % (len(missing), missing[0]))
        present = [k for k, r in enumerate(self.late_rows) if int(r) in pos]
        late_pos = np.array([pos[int(self.late_rows[k])] for k in present], dtype=np.int64)
        grid_pos = np.array([pos[int(r)] for r in self.grid_rows], dtype=np.int64)
        return late_pos, grid_pos, self.W[present]


class WindowOperator:
    """y -> y with y[late] = W @ y[grid] (copy)."""

    def __init__(self, layout, rows):
        self.late_pos, self.grid_pos, self.W = layout.operator_for(rows)
        self.nrows = len(rows)
        self.calls = 0

    def apply(self, y):
        y = np.array(y, dtype=float, copy=True)
        if y.shape != (self.nrows,):
            raise ValueError("late window: array of %s for %d rows" % (y.shape, self.nrows))
        y[self.late_pos] = self.W @ y[self.grid_pos]
        return y


_AVERAGED = {}


def averaged_class(base):
    """Subclass of a Model class whose evaluate applies the instance's ``_late_window`` to pNG, pRG and every
    Jacobian column. Swapped in per instance (``__class__``), so it survives the grain writer's temporary
    ``MF.evaluate`` instance attribute and its ``del`` (grain_crossing.write_fit)."""
    if getattr(base, "_late_averaged", False):
        return base
    if base not in _AVERAGED:
        class LateAveragedModel(base):
            _late_averaged = True

            def evaluate(self, *args, **kwargs):
                out = super().evaluate(*args, **kwargs)
                op = self._late_window
                op.calls += 1
                pNG, pRG = op.apply(out[0]), op.apply(out[1])
                if len(out) == 2:
                    return pNG, pRG
                J = {key: (op.apply(dn), op.apply(dr)) for key, (dn, dr) in out[2].items()}
                return pNG, pRG, J

        LateAveragedModel.__name__ = LateAveragedModel.__qualname__ = "LateAveraged" + base.__name__
        _AVERAGED[base] = LateAveragedModel
    return _AVERAGED[base]


def install_window(model, layout):
    model.__class__ = averaged_class(type(model))
    model._late_window = WindowOperator(layout, model.rows)
    return model._late_window


# ------------------------------------------------------------------------------------------------ start, Tikhonov
def gram_matrices(Y, area, A_orig, T, chunk=1 << 20):
    """G_all = B^T diag(area) B (Tikhonov), G_pos and r_pos of the start projection: weights w = area * A_orig / max
    A_orig over A_orig > 0 (emission weighting, user 2026-11-07), rhs ln A_orig; B = [Y, T]."""
    n, ny = Y.shape
    k = T.shape[1]
    m = ny + k
    G_all, G_pos, r_pos = np.zeros((m, m)), np.zeros((m, m)), np.zeros(m)
    amax = float(np.max(A_orig))
    for i in range(0, n, chunk):
        sl = slice(i, min(i + chunk, n))
        B = np.empty((sl.stop - sl.start, m))
        B[:, :ny] = Y[sl]
        B[:, ny:] = T[sl]
        w = np.asarray(area[sl], float)
        a = np.asarray(A_orig[sl], float)
        G_all += B.T @ (w[:, None] * B)
        wp = np.where(a > 0.0, w * (a / amax), 0.0)
        G_pos += B.T @ (wp[:, None] * B)
        r_pos += B.T @ (wp * np.log(np.where(a > 0.0, a, 1.0)))
    return G_all, G_pos, r_pos


def exponent_field(Y, a, T=None, c=None, chunk=1 << 20):
    s = np.empty(len(Y))
    for i in range(0, len(Y), chunk):
        sl = slice(i, min(i + chunk, len(Y)))
        s[sl] = Y[sl] @ a
        if T is not None and T.shape[1]:
            s[sl] += T[sl] @ c
    return s


def projection_stats(s0, area, A_orig):
    """How well exp(s0) represents A_orig: area- and emission-weighted rms of ln(exp(s0)/A_orig) over A_orig > 0, the
    area (of the total) and the activity/emission (area x A_orig) fractions with a factor > 2, the A_orig = 0 area,
    total activity ratio."""
    pos = A_orig > 0.0
    res = np.where(pos, s0 - np.log(np.where(pos, A_orig, 1.0)), 0.0)
    off = pos & (np.abs(res) > math.log(2.0))
    total, pos_area, act = float(np.sum(area)), float(np.sum(area[pos])), float(np.sum(area * A_orig))
    return {"facets": int(len(area)), "area_fraction_A_orig_zero": float(np.sum(area[~pos]) / total),
            "rms_ln_residual_over_A_orig_gt_0": float(np.sqrt(np.sum(area * res ** 2) / pos_area)),
            "rms_ln_residual_emission_weighted": float(np.sqrt(np.sum(area * A_orig * res ** 2) / act)),
            "max_abs_ln_residual": float(np.max(np.abs(res))),
            "area_fraction_factor_gt_2_of_total": float(np.sum(area[off]) / total),
            "area_fraction_factor_gt_2_of_A_orig_gt_0": float(np.sum(area[off]) / pos_area),
            "activity_fraction_factor_gt_2": float(np.sum((area * A_orig)[off]) / act),
            "total_activity_exp_s0_over_A_orig": float(np.sum(area * np.exp(s0)) / act)}


def start_projection(Y, area, A_orig, T, grams=None):
    """Emission-weighted (w = area A_orig) LS projection of ln A_orig (A_orig > 0) onto [Y, T]; a negative template
    coefficient is set to 0 and Y alone re-projected (user 2026-11-07; one template at most)."""
    ny, k = Y.shape[1], T.shape[1]
    if k > 1:
        raise SystemExit("start projection: the rule for more than one template is a user decision")
    G_all, G_pos, r_pos = grams if grams is not None else gram_matrices(Y, area, A_orig, T)
    theta = np.linalg.solve(G_pos, r_pos)
    projected_c = theta[ny:].copy()
    clamped = bool(k and theta[ny] < 0.0)
    if clamped:
        theta = np.concatenate([np.linalg.solve(G_pos[:ny, :ny], r_pos[:ny]), [0.0]])
    a, c = theta[:ny], theta[ny:]
    s0 = exponent_field(Y, a, T, c)
    stats = dict(projection_stats(s0, area, A_orig),
                 weighting=PROJECTION_WEIGHTING, weight_definition=PROJECTION_WEIGHT_DEFINITION,
                 condition_number_gram_pos=float(np.linalg.cond(G_pos)),
                 projected_template_coefficients=[float(v) for v in projected_c],
                 template_coefficient_clamped_to_0_and_Y_reprojected=clamped)
    ascale = float(np.max(np.abs(a)))
    return {"a": a, "c": c, "s0": s0, "ascale": ascale, "stats": stats, "G_all": G_all}


def log_tikhonov(G_all, area_total, ascale, ny):
    """Lt with Lt^T Lt = G = S G_all S / sum(area), S = diag(ascale (Y columns), 1 (templates))."""
    m = G_all.shape[0]
    S = np.concatenate([np.full(ny, ascale), np.ones(m - ny)])
    G = (S[:, None] * G_all * S[None, :]) / area_total
    G = 0.5 * (G + G.T)
    values, vectors = np.linalg.eigh(G)
    tol = np.finfo(float).eps * m * max(1.0, float(np.max(np.abs(values))))
    if np.min(values) < -tol:
        raise ValueError("log Tikhonov Gram matrix is not positive semidefinite")
    Lt = np.sqrt(np.maximum(values, 0.0))[:, None] * vectors.T
    check = float(np.max(np.abs(Lt.T @ Lt - G)) / np.max(np.abs(G)))
    return Lt, G, {"eigenvalue_range": [float(values.min()), float(values.max())], "LtT_Lt_minus_G_rel": check}


# ------------------------------------------------------------------------------------------------ installation
def install_selection(fit_cops, gkeys, record):
    """The named profile audits the observation rows; the grid rows get use 0, reason late_grid_point."""
    original = fit_cops.apply_selection_profile

    def apply_selection_profile(rows, columns, profile, run_intervals, saturated_ng):
        is_grid = [time_key(r[columns["UTC"]]) in gkeys for r in rows]
        if not any(is_grid):
            return original(rows, columns, profile, run_intervals, saturated_ng)
        obs = [i for i, g in enumerate(is_grid) if not g]
        grid = [i for i, g in enumerate(is_grid) if g]
        for i in grid:
            r = rows[i]
            if float(r[columns["refNG"]]) != 0.0 or float(r[columns["refRG"]]) != 0.0:
                raise SystemExit("late grid row %s has reference data" % r[columns["UTC"]])
        res, report = original([rows[i] for i in obs], columns, profile, run_intervals,
                               [saturated_ng[i] for i in obs])
        result = [None] * len(rows)
        reasons = {g: [None] * len(rows) for g in report["reasons"]}
        for j, i in enumerate(obs):
            result[i] = res[j]
            for g in reasons:
                reasons[g][i] = report["reasons"][g][j]
        for i in grid:
            r = list(rows[i])
            r[columns["useNG"]] = r[columns["useRG"]] = "0"
            result[i] = r
            for g in reasons:
                reasons[g][i] = GRID_REASON
        new = dict(report)
        new["reasons"] = reasons
        new["reason_counts"] = {g: dict(sorted(dict(report["reason_counts"][g], **{GRID_REASON: len(grid)}).items()))
                                for g in reasons}
        new["reason_sha256_observation_rows"] = report["reason_sha256"]
        new["reason_sha256"] = fit_cops.canonical_sha256(
            [[str(row[columns["UTC"]]), reasons["NG"][i], reasons["RG"][i]] for i, row in enumerate(result)])
        new["counts"] = dict(report["counts"], late_grid_rows=len(grid), base_rows=len(rows))
        new["late_grid_rows_note"] = ("model_ext/positive_problem.py: the profile was applied to the %d observation "
                                      "rows (its audited counts); the %d late grid rows (reference 0) have use 0 "
                                      "and reason %s" % (len(obs), len(grid), GRID_REASON))
        record["selection"] = {"observation_rows": len(obs), "grid_rows": len(grid),
                               "profile_counts_observation_rows": jsonable(report["counts"]),
                               "reason_sha256_observation_rows": report["reason_sha256"],
                               "reason_sha256_all_rows": new["reason_sha256"]}
        print("positive_fit: profile %s applied to the %d observation rows; %d grid rows not used (%s)"
              % (profile.get("name"), len(obs), len(grid), GRID_REASON), flush=True)
        return result, new

    fit_cops.apply_selection_profile = apply_selection_profile


def install_core(fit_cops, cfg, record):
    """Before grain_crossing.install: selection wrapper, the exponent Problem, activity writers, Jacobian check."""
    table = read_point_table(cfg["point_table"])
    record.update(policy=POLICY, parameterization=PARAMETERIZATION_ID,
                  point_table={"path": cfg["point_table"], "sha256": sha256(cfg["point_table"]), "rows": len(table)},
                  user_temperature_K=list(USER_T0), late_start_utc=LATE_START_UTC)
    gkeys = grid_keys(table)
    install_selection(fit_cops, gkeys, record)
    base = fit_cops.Problem

    class PositiveCore(base):
        def __init__(self, M, el, coef, args):
            if args.degree != 4 or args.enhancement_degree != 4 or args.fix_legacy_activity:
                raise SystemExit("positive_problem: needs --degree 4 --enhancement-degree 4, free legacy coefficients")
            if not args.fix_temperature or args.fix_kappa:
                raise SystemExit("positive_problem: needs --fix-temperature and a free kappa (user 2026-11-07)")
            if not np.allclose(np.asarray(M.T0, float), USER_T0, rtol=0.0, atol=1e-9):
                raise SystemExit("positive_problem: the run's temperature law %s differs from the user's %s"
                                 % (list(M.T0), list(USER_T0)))
            if list(getattr(args, "activity_template_names", [])) != [REQUIRED_TEMPLATE]:
                raise SystemExit("positive_problem: exactly the template %s is required (user 2026-11-07), got %s"
                                 % (REQUIRED_TEMPLATE, list(getattr(args, "activity_template_names", []))))
            # late window: grid rows join the fit rows; W10x5 on the fit model and on the all-rows model MF
            points = read_run_points(M.run)
            layout = LateLayout(table, points, M.ipoint)
            rows = np.union1d(np.asarray(M.rows, dtype=np.int64), layout.grid_rows)
            if len(rows) != len(M.rows) or np.any(rows != np.asarray(M.rows)):
                M.load_kernels(rows)
            M.sel = rows
            if np.any(M.refNG[layout.grid_rows] != 0.0) or np.any(M.refRG[layout.grid_rows] != 0.0):
                raise SystemExit("positive_problem: grid rows with reference data")
            op_fit = install_window(M, layout)
            MF = getattr(fit_cops, "MF", None)
            op_all = install_window(MF, layout) if MF is not None and MF is not M else None
            shared = MF is not None and len(MF.rows) == len(M.rows) and bool(np.all(MF.rows == M.rows))
            record["late_window"] = dict(layout.summary, fit_rows=int(len(M.rows)),
                                         fit_rows_equal_all_rows=shared,
                                         late_rows_in_fit_model=int(len(op_fit.late_pos)),
                                         installed_on_all_rows_model=op_all is not None)
            if not shared:
                print("positive_fit: WARNING the fit rows (%d) differ from all rows (%d): a second GPU backend "
                      "is bound for them" % (len(M.rows), len(MF.rows) if MF is not None else -1), flush=True)
            self.late_layout = layout
            super().__init__(M, el, coef, args)
            if len(self.fit_j) != 25 or self.nh != 0 or self.nt > 1:
                raise SystemExit("positive_problem: unexpected layout (fit_j %d, nh %d, nt %d)"
                                 % (len(self.fit_j), self.nh, self.nt))
            # exponent parameterization, start = emission-weighted projection of ln A_orig
            T = np.ascontiguousarray(self.enhancement_basis[:, self.nh:])
            grams = gram_matrices(el["Y"], el["area"], el["AH2O"], T)
            proj = start_projection(el["Y"], el["area"], el["AH2O"], T, grams)
            self.ascale = proj["ascale"]
            self.enhancement_scale = 1.0
            x0 = np.array(self.x0, float)
            x0[:self.na] = proj["a"][self.fit_j] / self.ascale
            x0[self.na:self.na + self.nt] = proj["c"]
            if not (np.all(x0 >= self.lo) and np.all(x0 <= self.hi)):
                raise SystemExit("positive_problem: the projected start is outside the bounds")
            self.x0 = x0
            self.s0 = proj["s0"]
            self.area_total = float(np.sum(el["area"]))
            Lt, _, reg = log_tikhonov(grams[0], self.area_total, self.ascale, len(self.fit_j))
            self.Lt = Lt
            self.reg_slice = slice(0, self.na + self.nb)
            names = list(self.names)
            record["start_projection"] = dict(
                proj["stats"], ascale=self.ascale,
                a=[float(v) for v in proj["a"]], c={("c_" + n): float(v) for n, v in zip(self.template_names, proj["c"])},
                x0_scaled=[float(v) for v in x0], names=names)
            record["regularization"] = dict(
                reg, identifier="positive_log_activity_deviation_from_start_v1",
                formula="D^2 = sum(area*(s - s0)^2)/sum(area) over all facets; s = Y a + c T (exponent), s0 = start",
                prior="the start x0 = the emission-weighted (area A_orig) LS projection of ln A_orig; D = rms "
                      "deviation of ln A from it (area-weighted penalty)", reg_slice=[0, self.na + self.nb],
                note="fit_manifest objective.regularization (and N0 = prior.legacy_linear_area_power) describe the "
                     "reference definition and do not apply; P.N0 is left at the legacy value")
            st = proj["stats"]
            record["start_projection_weighting"] = PROJECTION_WEIGHTING
            print("positive_fit: start = emission-weighted projection of ln A_orig on [Y, %s]: rms ln residual %.4f "
                  "(area-weighted, A_orig > 0), %.4f (emission-weighted), emission with factor > 2: %.4f, "
                  "area with factor > 2: %.4f of total, A_orig = 0 on %.4f of the area, total exp(s0)/A_orig %.4f; "
                  "c = %s%s; ascale %.6g" % (
                      ", ".join(self.template_names) or "-", st["rms_ln_residual_over_A_orig_gt_0"],
                      st["rms_ln_residual_emission_weighted"], st["activity_fraction_factor_gt_2"],
                      st["area_fraction_factor_gt_2_of_total"], st["area_fraction_A_orig_zero"],
                      st["total_activity_exp_s0_over_A_orig"], [float(v) for v in proj["c"]],
                      " (projected %s < 0: set to 0, Y re-projected)" % st["projected_template_coefficients"]
                      if st["template_coefficient_clamped_to_0_and_Y_reprojected"] else "", self.ascale), flush=True)

    fit_cops.Problem = PositiveCore
    install_activity_writers(fit_cops, record, gkeys)
    if cfg.get("check_jacobian"):
        install_jacobian_check(fit_cops, cfg["check_parameters"], record)
    return PositiveCore


def install_outer(fit_cops, cfg, record):
    """After grain_crossing.install: grain start (user 2026-11-07), late-start check, exp products, restart marker."""
    base = fit_cops.Problem

    class PositiveGrainProblem(base):
        def __init__(self, M, el, coef, args):
            super().__init__(M, el, coef, args)
            names = list(self.names)
            if "grain_A_uPa" not in names:
                raise SystemExit("positive_problem: the grain parameters are missing")
            k = names.index("grain_A_uPa")
            default = float(self.x0[k])
            self.x0[k] = cfg["grain_start_uPa"]
            self.x0[k + 1] = 0.0
            record["grain_start"] = {"grain_A_uPa": float(cfg["grain_start_uPa"]), "grain_B_uPa_per_100s": 0.0,
                                     "replaced_grain_crossing_default_A_uPa": default,
                                     "source": "user 2026-11-07: cycle 3 fine-tuned A_Pa (grain/cycles/c3/finetune/"
                                               "C30_n4000_grain_cross10_ftc3_m5_i5_cpu/grain_finetune.json), B = 0"}
            # the late window starts at the first post-gap NG point (user 2026-11-07)
            layout = self.late_layout
            i_end = int(np.nonzero(self.et_fit == self.et_end)[0][0])
            row_end = int(np.asarray(M.rows)[i_end])
            if row_end != int(layout.late_rows[0]):
                raise SystemExit("positive_problem: the first post-gap NG point (row %d) is not the first late-window "
                                 "observation (row %d, %s): the window start needs a user decision"
                                 % (row_end, int(layout.late_rows[0]), LATE_START_UTC))
            record["late_window"]["grain_t_end_row"] = row_end
            record["late_window"]["grain_t_end_equals_late_start"] = True

    fit_cops.Problem = PositiveGrainProblem
    original_write = fit_cops.write_fit

    def write_fit(out, r, *rest):
        diag = original_write(out, r, *rest)
        write_exp_products(out, fit_cops.P, r, rest[6], rest[7], record)
        write_late_window_products(out, fit_cops, r, record)
        return diag
    fit_cops.write_fit = write_fit
    install_restart_marker(record)
    install_figure_filter(grid_keys(read_point_table(cfg["point_table"])), record)
    return PositiveGrainProblem


def install_figure_filter(gkeys, record):
    """gpu_fit.write_figures plots cops_cuda.figures.read_table(base) and (fit) over all rows: drop the late grid rows
    (model-only points) from both tables, so the figures show observation rows only."""
    from cops_cuda import figures
    if getattr(figures.read_table, "_positive_grid_filter", False):
        return
    original = figures.read_table

    def read_table(path):
        out = original(path)
        keep = np.array([time_key(u) not in gkeys for u in out["utc"]], dtype=bool)
        if keep.all():
            return out
        n = len(keep)
        new = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape[:1] == (n,) else v) for k, v in out.items()}
        new["utc"] = [u for u, k in zip(out["utc"], keep) if k]
        record.setdefault("figures_grid_rows_dropped", {})[str(out["path"])] = int(n - keep.sum())
        return new
    read_table._positive_grid_filter = True
    figures.read_table = read_table


# ------------------------------------------------------------------------------------------------ writers
def exp_state(P, el, a, b):
    """The exponent field s = Y a + enhancement_scale * B b (enhancement_scale = 1 here)."""
    s = exponent_field(el["Y"], np.asarray(a, float))
    for j in range(P.nb):
        s += P.enhancement_scale * float(b[j]) * P.enhancement_basis[:, j]
    return s


def install_activity_writers(fit_cops, record, gkeys=frozenset()):
    """activity_diagnostics and the A5 activity products in exp semantics (see the module docstring); the A5 pressure
    products (pressure_residuals.dat/png, pressure_comparison_*.png) for the observation rows only (no grid rows)."""

    def activity_diagnostics(a, b, el, seen, coef, enhancement_basis, enhancement_scale):
        s = el["Y"] @ a
        if len(b):
            s = s + enhancement_scale * (enhancement_basis @ b)
        act = np.exp(s)
        act0 = np.maximum(el["Y"] @ coef[0], 0.0)
        area = el["area"]
        tot, tot0 = np.sum(area * act), np.sum(area * act0)
        d = {"active_area": np.sum(area[act > 0]) / np.sum(area),
             "activity_ratio": tot / tot0,
             "activity_ratio_seen": np.sum((area * act)[seen]) / np.sum((area * act0)[seen]),
             "unseen_activity_fraction": np.sum((area * act)[~seen]) / tot if tot > 0 else np.nan,
             "unseen_activity_fraction_prior": np.sum((area * act0)[~seen]) / tot0,
             "unseen_area_fraction": np.sum(area[~seen]) / np.sum(area)}
        if len(b):
            delta = act - np.exp(el["Y"] @ a)          # activity added by the template factor(s)
            n0 = np.sum(area * (el["Y"] @ coef[0]) ** 2)
            d.update({"enhancement_relative_rms_D": np.sqrt(np.sum(area * delta ** 2) / n0),
                      "enhancement_area_weighted_mean": np.sum(area * delta) / np.sum(area),
                      "clipped_area_fraction": 0.0})
        return d

    fit_cops.activity_diagnostics = activity_diagnostics
    original_a5 = fit_cops.write_stage_a5_report

    def write_stage_a5_report(out, P, r, rows, cols, el, seen, a, b, *rest):
        original_wap = fit_cops.write_activity_products

        def write_activity_products(directory, identity, label, centroid, area, seen_, baseline, residual,
                                    total_linear, clipped, degree_components, template_components=()):
            s = exp_state(P, el, a, b)
            A = np.exp(s)
            tc = []
            for j, name in enumerate(P.template_names):
                c = float(b[P.nh + j])
                tc.append((name, c, A - np.exp(s - c * P.enhancement_scale * P.enhancement_basis[:, P.nh + j])))
            return original_wap(directory, identity, label, centroid, area, seen_, baseline, A - baseline, A, A,
                                [], template_components=tc)
        original_wpp = fit_cops.write_pressure_products

        def write_pressure_products(directory, identity, label, point_ids, utc_text, altitude, predictions,
                                    references, used, saturated, reasons, transition_intervals=()):
            idx = [i for i, u in enumerate(utc_text) if time_key(u) not in gkeys]
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
        fit_cops.write_activity_products = write_activity_products
        fit_cops.write_pressure_products = write_pressure_products
        try:
            return original_a5(out, P, r, rows, cols, el, seen, a, b, *rest)
        finally:
            fit_cops.write_activity_products = original_wap
            fit_cops.write_pressure_products = original_wpp

    fit_cops.write_stage_a5_report = write_stage_a5_report
    record["activity_writers"] = {
        "patched": ["fit_cops.activity_diagnostics (exp)", "fit_cops.write_activity_products via write_stage_a5_report "
                    "(legacy_linear = Y coef[0], residual = A - legacy, total = clipped = A = exp(s); template "
                    "component = A - exp(s - c T); no degree rows)"],
        "not_replayable_do_not_use": ["ActivityFit.txt", "ActivityEnhancementFit.txt", "ActivityEnhancementFacets.dat",
                                      "cops_fit.in", "fit_restart.npz activity_linear (= exponent s)"]}


PARAMETERIZATION_TEXT = """\
# ACTIVITY_PARAMETERIZATION.txt -- written by model_ext/positive_problem.py ({utc})
# Fit with the positive H2O activity parameterization {pid} (user decisions 2026-11-07).

H2O activity per facet:  A = exp(s),  s = sum_j a_j Y_j + c * T_cliff   (Y: the 25 degree-<=4 real spherical harmonics
of the fitter; T_cliff: template {templates} of {bundle}; c = ln of the cliff factor, >= 0).
Production per kernel record (GPU, model_ext/positive_activity.py): night/shadow ratio {night} (H2O and CO2; reference
0.2 and 0.1); H2O flux floor F >= {floor} q_max A_orig (A_orig = the run's activity, elements.bin column 32; records at
the floor do not depend on the activity parameters; CO2 has no floor). Temperature law held: Tmax {T0[0]},
T75 {T0[1]}, Tmin {T0[2]} K. Late window (>= {late} UTC): model = W10x5 of the surface model (+ grain at the time stamp).
Tikhonov: D^2 = sum area (s - s0)^2 / sum area (nepers; s0 = the start and prior = the emission-weighted, w = area A_orig,
LS projection of ln A_orig on the same basis). lambda = {lam:g}.

VALID in this directory:
  activity_exp_facets.npz       s, A = exp(s), s0, F_min_over_qmax = {floor} A_orig, seen  (per facet, mesh order)
  activity_facets.dat           legacy_linear = Y coef[0] of the run; residual_linear = A - legacy_linear;
                                total_linear = clipped_activity = A = exp(s)  (patched writer)
  activity_template_power.dat   template coefficient c and the activity added by the factor, A - exp(Y a)
  activity_region_diagnostics.dat, activity_maps.png   from the fields above ("clipped" = A; nothing is clipped)
  fit_result.txt                parameters (aNN = exponent coefficients, x ascale), costs, D; its "activity
                                diagnostics" lines are computed with A = exp(s) although the header line names the clip
  COPS_model.dat                observation rows: NG/RG of this model; the late observations (>= {late} UTC):
                                pNG = W10x5(surface) + grain(t_obs), pRG = W10x5(surface); earlier: instantaneous
                                (+ grain after t_end). GRID ROWS ({ngrid}, kind "grid" of the point table, reason
                                {reason}, use 0, listed in late_window_model.dat): model-only 1-s points, pNG =
                                INSTANTANEOUS surface + grain(t), pRG = instantaneous surface; NOT observations: drop
                                them in plots/comparisons and do NOT re-average them (pNG includes the grain); a header
                                line of COPS_model.dat states this
  late_window_model.dat         the late window row by row: kind (grid | late_obs), et, L2 STOP et, surface-only model
                                (grid: instantaneous; late_obs: W10x5), grain, the table's pNG/pRG, data, fit use
  pressure_residuals.*, pressure_comparison_*.png, figures/   observation rows only (grid rows removed)
  grain_term.*, fit_manifest.json ("model_extension": everything needed to rerun)
NOT VALID -- do not use (written by the reference fitter for the linear max(0, Y a + ...) model; LiouvilleCOPS C++
cannot replay an exp activity):
  ActivityFit.txt (the H2O row holds exponent coefficients), ActivityEnhancementFit.txt, ActivityEnhancementFacets.dat,
  cops_fit.in (would run the C++ model with these numbers as linear coefficients), activity_degree_power.dat (empty),
  fit_restart.npz "activity_linear" (= s, not a linear activity): tools that use max(activity_linear, 0)
  (fine_tune_gpu.py, fine_tuning.py, restart_tool.py, grain_finetune.py, late_response_eval.py) are not valid for this
  fit. The restart file can continue this fit with positive_fit.py --restart-from (marker meta "model_extension").
"""


def write_exp_products(out, P, r, el, seen, record):
    out = Path(out)
    a, b, T, kappa = P.unpack(r["x"])
    s = exp_state(P, el, a, b)
    A = np.exp(s)
    floor = float(os.environ.get("COPS_POS_FLOOR_RATIO", "nan"))
    night = os.environ.get("COPS_POS_NIGHT_RATIO", "?")
    meta = {"parameterization": PARAMETERIZATION_ID, "created_utc": utc_now(), "lambda": float(r["lam"]),
            "names": list(P.names), "x": [float(v) for v in r["x"]], "ascale": float(P.ascale),
            "floor_ratio": floor, "night_ratio": night, "templates": list(P.template_names),
            "run": str(getattr(P.args, "run", ""))}
    path = out / "activity_exp_facets.npz"
    with path.open("xb") as f:
        np.savez_compressed(f, s=s, A=A, s0=np.asarray(P.s0, float),
                            F_min_over_qmax=floor * np.asarray(el["AH2O"], float),
                            seen=np.asarray(seen, bool), meta_json=np.array(json.dumps(meta, sort_keys=True)))
    bundle = [src["bundle"]["path"] for src in getattr(P, "template_sources", [])]
    text = PARAMETERIZATION_TEXT.format(utc=utc_now(), pid=PARAMETERIZATION_ID,
                                        templates=", ".join(P.template_names) or "(none)",
                                        bundle=", ".join(bundle) or "-", night=night, floor=floor, T0=USER_T0,
                                        late=LATE_START_UTC, lam=float(r["lam"]),
                                        ngrid=len(P.late_layout.grid_rows), reason=GRID_REASON)
    tpath = out / "ACTIVITY_PARAMETERIZATION.txt"
    with tpath.open("x") as f:
        f.write(text)
    D = float(np.sqrt(np.sum(el["area"] * (s - P.s0) ** 2) / P.area_total))
    record["products"] = {p.name: {"sha256": sha256(p), "size_bytes": p.stat().st_size} for p in (path, tpath)}
    record["result_summary"] = {"D_direct_nepers": D, "D_fitter": float(r["D"]),
                                "total_activity_over_A_orig": float(np.sum(el["area"] * A) / np.sum(el["area"] * el["AH2O"])),
                                "A_range": [float(A.min()), float(A.max())], "kappa": float(kappa)}
    print("positive_fit: wrote %s, %s (D direct %.6f nepers, fitter %.6f)" % (path, tpath, D, r["D"]), flush=True)


def write_late_window_products(out, fit_cops, r, record):
    """late_window_model.dat (surface-only model, grain, kinds) and the header note of COPS_model.dat on the grid rows."""
    out = Path(out)
    P, MF = fit_cops.P, fit_cops.MF
    lay = P.late_layout
    a, b, T, kappa = P.unpack(r["x"])
    pNG_s, pRG_s = MF.evaluate(a, T, kappa, b=b if P.nb else None, enhancement_basis=P.enhancement_basis,
                               enhancement_scale=P.enhancement_scale)      # surface model, late observations W10x5
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
    worst = 0.0
    with (out / "late_window_model.dat").open("x") as f:
        f.write("# late window of this fit (model_ext/positive_problem.py, %s): every 1-s grid row (model only) and every "
                "observation from %s UTC on\n" % (utc_now(), LATE_START_UTC))
        f.write("# kind grid: surface_* = instantaneous surface model at et; kind late_obs: surface_* = W10x5 of the "
                "surface model ending at the L2 STOP time (stop_et); grain_NG_Pa = grain(et) (NG only, unaveraged);\n"
                "# table_pNG/pRG = COPS_model.dat (grid: instantaneous surface + grain; late_obs: W10x5 + grain); "
                "re-average only surface_* of grid rows, never table_pNG\n")
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
    note = ["# PositiveFitGridRows (model_ext/positive_problem.py): the %d rows of kind grid (1-s late grid of %s; "
            "reference 0, use 0, reason %s) are MODEL-ONLY points: pNG = instantaneous surface model + grain(t), "
            "pRG = instantaneous surface model. They are not observations: drop them in plots and comparisons and do "
            "not re-average them (pNG includes the grain); late_window_model.dat lists them (kind grid) with the "
            "surface-only model.\n" % (len(lay.grid_rows), record.get("point_table", {}).get("path"), GRID_REASON),
            "# PositiveFitLateObservations: the %d observations from %s UTC on have pNG = W10x5(surface) + grain(t_obs), "
            "pRG = W10x5(surface) (L2 STOP times in late_window_model.dat); earlier rows: instantaneous model.\n"
            % (len(lay.late_rows), LATE_START_UTC)]
    tmp = table.with_name(table.name + ".tmp.%d" % os.getpid())
    tmp.write_text("".join(lines[:k_hdr] + note + lines[k_hdr:]))
    os.replace(str(tmp), str(table))
    record["late_window_products"] = {
        "late_window_model.dat": {"sha256": sha256(out / "late_window_model.dat"), "rows": len(rows)},
        "COPS_model.dat_header_note": "PositiveFitGridRows / PositiveFitLateObservations lines before #iPoint",
        "table_vs_surface_plus_grain_max_rel": float(worst)}
    print("positive_fit: wrote %s (%d rows; COPS_model.dat = surface + grain to %.1e); grid rows flagged in "
          "COPS_model.dat" % (out / "late_window_model.dat", len(rows), worst), flush=True)


def install_restart_marker(record):
    from cops_cuda import restart as rs
    if getattr(rs.from_fit, "_positive_marker", False):
        return
    original = rs.from_fit

    def from_fit(problem, result, fit_dir, manifest, elements):
        state = original(problem, result, fit_dir, manifest, elements)
        state.meta["model_extension"] = {
            "parameterization": PARAMETERIZATION_ID, "policy": POLICY,
            "activity": "A = exp(activity_linear); activity_linear = Y a + c T (exponent)",
            "warning": "not a linear activity: tools that use max(activity_linear, 0) do not apply",
            "ascale": float(problem.ascale), "positive_activity": record.get("positive_activity"),
            "start_projection_weighting": PROJECTION_WEIGHTING,
            "late_window": {k: record.get("late_window", {}).get(k) for k in ("late_start_utc", "late_observations")}}
        return state
    from_fit._positive_marker = True
    rs.from_fit = from_fit


# ------------------------------------------------------------------------------------------------ Jacobian check
def install_jacobian_check(fit_cops, names, record):
    """A least_squares wrapper installed BEFORE grain_crossing.install, i.e. inside the crossing wrapper: it sees the
    full objective (data residuals with grain and window averaging, Tikhonov rows for lambda > 0, crossing rows) and
    compares analytic Jacobian columns with central differences at the optimizer start, once. Its 2 x len(names)
    residual evaluations pass the wall-clock wrapper (counted there; a probe could be its 'best state' only if the
    deadline came before the first improving step)."""
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
    nreg = int(P.Lt.shape[0]) if P.lam > 0 else 0
    late = np.concatenate([P.et_fit[P.useNG] >= P.et_end, P.et_fit[P.useRG] >= P.et_end])
    blocks = {"data": np.arange(nd), "data_late_window": np.nonzero(late)[0],
              "regularization": np.arange(nd, nd + nreg), "crossing": np.arange(nd + nreg, len(r0))}
    out = {"lambda": float(P.lam), "rows": {k: int(len(v)) for k, v in blocks.items()},
           "metric": "max|analytic - central FD| / max|FD| per block (abs: max|analytic - FD|)", "columns": {}}
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
        print("positive_fit: Jacobian check d/d%-18s h %.2e  data %.2e  late %.2e  crossing %s" % (
            name, h, col.get("data", float("nan")), col.get("data_late_window", float("nan")),
            "%.2e" % col["crossing"] if "crossing" in col else "-"), flush=True)
    worst = [v["data"] for v in out["columns"].values() if isinstance(v, dict) and "data" in v]
    out["worst_data"] = float(max(worst)) if worst else None
    return out
