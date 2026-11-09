# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/positive_problem.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/positive_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b234281c22bebe49dc4078b52b2b1dbd2bb37f1146cae5cf942dfb2a741570f9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Fit driver of the positive-activity fits (user decisions 2026-11-07/08): exp(Y a + c T_cliff + held), c >= 0; start and prior
#   = emission-weighted projection of ln A_orig; area-weighted log Tikhonov on the total ln A (held: Lt (x - x0 + h_B) plus
#   the constant |h_perp|); W10x5 late window; grid rows; ratio, start restart and start check; products and restart
#   marker (ratio, held file sha256).
# USED BY: positive_fit.py, positive_finetune.py.
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
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
  COPS_POS_FLOOR_RATIO); this module does not reimplement them. The ratio is one of RATIO_DECISIONS (0.02, 0.002,
  1e-5); the decision of 1e-5 (user 2026-11-08) also includes the CO2 flux floor with the same ratio
  (COPS_POS_FLOOR_CO2=1: F_CO2 = max(kappa Q_CO2 A_CO2,run, r q_max,CO2 A_CO2,run), A_CO2,run = the run's ACO2 map;
  CO2_FLOOR_DECISIONS; recorded in the manifest ("ratio", "environment"), the restart marker ("co2_floor"; a resume
  must have the same), the start check (a CO2-floor change is a change of the problem), ACTIVITY_PARAMETERIZATION.txt
  and activity_exp_facets.npz (meta co2_floor, F_min_CO2_over_qmax_CO2)).
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
import textwrap
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
                ("COPS_POS_FLOOR_RATIO", "0.02"),  # agent B's patch; the ratio entries are the DEFAULT (0.02)
                ("COPS_POS_FLOOR_CO2", None))      # 2026-11-08: set from the ratio's decision (CO2_FLOOR_DECISIONS)
ENV_NIGHT, ENV_FLOOR, ENV_FLOOR_CO2 = "COPS_POS_NIGHT_RATIO", "COPS_POS_FLOOR_RATIO", "COPS_POS_FLOOR_CO2"
DEFAULT_RATIO = "0.02"
# The night/shadow production ratio and the H2O flux-floor ratio are the user's decisions; one value R sets both
# (COPS_POS_NIGHT_RATIO = COPS_POS_FLOOR_RATIO = R). Any other value is refused until the user decides it.
RATIO_DECISIONS = {
    0.02: "user 2026-11-07 (first decision): night/shadow ratio 0.02 'as in srcCG' and H2O flux floor 0.02 q_max A_orig",
    0.002: "user 2026-11-07 (second decision): night/shadow ratio and H2O floor ratio both 0.002 (CO2 no floor)",
    1e-5: "user 2026-11-08 (third decision): 'decrease the floor to 1.0E-5 and also use that for nightside and shadows'; night/shadow ratio and H2O floor ratio both 1e-5; then 'use the same floor for CO2': CO2 flux floor with the same ratio (the lead's reading: F_CO2 = max(q kappa A_CO2,run, 1e-5 q_max,CO2 A_CO2,run), A_CO2,run = the run's CO2 map ACO2 = the CO2 activity at kappa = 1; COPS_POS_FLOOR_CO2=1)"}
# The CO2 flux floor (COPS_POS_FLOOR_CO2=1, positive_activity.py) is part of a ratio's decision: on for 1e-5 only (user
# 2026-11-08); the earlier decisions have no CO2 floor (the variable must be unset for them).
CO2_FLOOR_DECISIONS = {0.02: False, 0.002: False, 1e-5: True}
CO2_FLOOR_TEXT = ("CO2 flux floor F_CO2 >= {floor} q_max,CO2 A_CO2,run (F_CO2 = max(kappa Q_CO2 A_CO2,run, {floor} "
                  "q_max,CO2 A_CO2,run), A_CO2,run = the run's CO2 map ACO2 (elements.bin column 33) = the CO2 activity "
                  "at kappa = 1; user 2026-11-08 'use the same floor for CO2', the lead's reading; records at the CO2 "
                  "floor do not depend on kappa)")
START_REPRODUCTION_RTOL = 1e-8     # a start that must reproduce its source (same ratio, lambda, held term): cost rel.
# --start-set NAME=VALUE (2026-11-08): explicit optimizer-start values applied after --start-restart (and the grain start)
START_SET_DEFINITION = ("--start-set NAME=VALUE[,NAME=VALUE...]: after --start-restart the optimizer start of each named "
                        "parameter of this fit is VALUE, in the fitter's optimizer coordinate x of that name (kappa: kappa "
                        "itself; aNN: a_NN / ascale; c_cliff_ge_60deg: ln of the cliff factor; grain_A_uPa: 1e-6 Pa; "
                        "grain_B_uPa_per_100s: 1e-6 Pa per 100 s); every other parameter keeps the start restart's value; "
                        "the prior x0 is unchanged; lo <= VALUE <= hi is required")
START_SET_DECISIONS = {"kappa": "user 2026-11-08 'kappa = 1 where stuck': the ratio-1e-5 cycle-1 fits of lambda 1e-3, 1e-2 "
                                "and 1e-1 start kappa at 1 (the fitter's default, the run's CO2 map); every other parameter "
                                "from the ratio-0.002 cycle-1 fit (--start-restart); lambda 0, 1e-5, 1e-4 keep the "
                                "restart's kappa"}
PARAMETERIZATION_ID = "positive_exp_activity_log_tikhonov_w10x5_emisproj_v2"   # v2: emission-weighted start/prior
PROJECTION_WEIGHTING = "emission_weighted_log_ls"
PROJECTION_WEIGHT_DEFINITION = ("w_e = area_e * A_orig,e (A_orig > 0; normalized by max A_orig, which does not change "
                                "the solution); min sum_e w_e (s_e - ln A_orig,e)^2")
CHECK_PARAMETERS = ("a00", "a02", "a07", "c_cliff_ge_60deg", "kappa", "grain_A_uPa")
GRID_REASON = "late_grid_point"
REQUIRED_TEMPLATE = "cliff_ge_60deg"                # the user's model (2026-11-07): exactly this template
POLICY_TEMPLATE = (
    "user 2026-11-07: A = exp(Y a + c_cliff T_cliff) (positive by construction; templates inside the exponent; "
    "ft_c* not used); start and prior = emission-weighted (w = area A_orig) LS projection of ln A_orig on "
    "[Y, T_cliff] over A_orig > 0 (negative "
    "cliff coefficient -> 0, Y re-projected); original analytical temperature law held (fit_cops "
    "--fix-temperature); kappa free >= 0 from 1; night ratio {ratio} and flux floor {ratio} q_max A_orig on the GPU "
    "(positive_activity.py); D^2 = sum area (s - s0)^2 / sum area (log space, prior = start); late window "
    "(>= 10:33:10 UTC) compared with W10x5 of the surface model (mean of five 10-s averages ending at "
    "T, T-2, ..., T-8 s, T = L2 STOP), grain added at the observation time; grain + crossing 10, grain start "
    "A = 6.96024299190473e-06 Pa, B = 0 (unless a start restart / fine tuning gives it)")
POLICY = POLICY_TEMPLATE.format(ratio="0.02")      # the first decision's text (kept for importers)


def current_ratio(environ=None):
    """The ratio R of this process (COPS_POS_NIGHT_RATIO = COPS_POS_FLOOR_RATIO = R), as a float."""
    env = os.environ if environ is None else environ
    return float(env.get(ENV_NIGHT, DEFAULT_RATIO))


def current_co2_floor(environ=None):
    """True when this process has the CO2 flux floor (COPS_POS_FLOOR_CO2=1; set by set_positive_environment)."""
    env = os.environ if environ is None else environ
    return (env.get(ENV_FLOOR_CO2) or "").strip() == "1"


def current_policy():
    text = POLICY_TEMPLATE.format(ratio=repr(current_ratio()))
    if current_co2_floor():
        text += "; " + CO2_FLOOR_TEXT.format(floor=repr(current_ratio()))
    return text


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
def _ratio_value(text, what):
    try:
        r = float(text)
    except (TypeError, ValueError):
        raise SystemExit("positive_fit: %s=%r is not a number" % (what, text)) from None
    if r not in RATIO_DECISIONS:
        raise SystemExit("positive_fit: %s=%s is not a ratio the user decided (decided: %s); a new value is a user "
                         "decision (add it to positive_problem.RATIO_DECISIONS with its date)"
                         % (what, text, ", ".join("%g" % v for v in sorted(RATIO_DECISIONS))))
    return r


def parse_start_set(text):
    """--start-set NAME=VALUE[,NAME=VALUE...] -> [[name, value], ...] (finite floats, unique names; checked against the
    fit's parameters and bounds when the Problem exists)."""
    out = []
    for item in str(text).split(","):
        item = item.strip()
        if not item:
            continue
        if item.count("=") != 1:
            raise SystemExit("positive_fit: --start-set %r: need NAME=VALUE[,NAME=VALUE...]" % text)
        name, value = (s.strip() for s in item.split("="))
        try:
            v = float(value)
        except ValueError:
            raise SystemExit("positive_fit: --start-set %s=%r: the value is not a number" % (name, value)) from None
        if not name or not math.isfinite(v):
            raise SystemExit("positive_fit: --start-set %r: need a parameter name and a finite value" % item)
        if name in [n for n, _ in out]:
            raise SystemExit("positive_fit: --start-set names %s twice" % name)
        out.append([name, v])
    if not out:
        raise SystemExit("positive_fit: --start-set %r sets nothing" % text)
    return out


def set_positive_environment(environ=None, ratio=None):
    """Set the COPS_POS_* variables of the GPU patch (before anything binds a backend; spawned workers inherit them).

    The night/shadow ratio and the H2O floor ratio are one value R (user decisions, RATIO_DECISIONS): from ``ratio``
    (positive_fit.py --ratio), else from a preset COPS_POS_NIGHT_RATIO / COPS_POS_FLOOR_RATIO (jobs/positive_fit.sbatch
    RATIO; spawned workers inherit the parent's), else the default 0.02 (the behaviour before 2026-11-07 evening).
    Refused: an undecided R, night != floor, a preset value different from ``ratio``, other preset COPS_POS_* values
    different from the user's (ENABLE 1, ACTIVITY exp).
    CO2 flux floor (2026-11-08): COPS_POS_FLOOR_CO2 follows the ratio's decision (CO2_FLOOR_DECISIONS): R = 1e-5 sets
    it to 1 (a preset value must be 1); for 0.02 and 0.002 it must be unset (any preset value, also 0, is refused, so
    the recorded environment of those fits is the one of the earlier ones). Returned (and recorded as the fit's
    "environment") only when set."""
    environ = os.environ if environ is None else environ
    preset = {k: environ.get(k) for k in (ENV_NIGHT, ENV_FLOOR) if environ.get(k) is not None}
    values = {k: _ratio_value(v, k) for k, v in preset.items()}
    if len(set(values.values())) > 1:
        raise SystemExit("positive_fit: %s differ; the user's decisions set both to one value" % preset)
    if ratio is not None:
        r = _ratio_value(ratio, "--ratio")
        if values and set(values.values()) != {r}:
            raise SystemExit("positive_fit: --ratio %s but the environment has %s" % (ratio, preset))
        text = str(ratio).strip()
    elif preset:
        text = next(iter(preset.values())).strip()
    else:
        text = DEFAULT_RATIO
    r = _ratio_value(text, "ratio")
    co2 = CO2_FLOOR_DECISIONS[r]
    old_co2 = environ.get(ENV_FLOOR_CO2)
    if co2 and old_co2 is not None and old_co2 != "1":
        raise SystemExit("positive_fit: %s=%s is preset; the decision of ratio %s includes the CO2 floor (%s=1, user "
                         "2026-11-08)" % (ENV_FLOOR_CO2, old_co2, text, ENV_FLOOR_CO2))
    if not co2 and old_co2 is not None:
        raise SystemExit("positive_fit: %s=%s is preset, but the decision of ratio %s has no CO2 floor (only 1e-5 has "
                         "one, user 2026-11-08): unset it" % (ENV_FLOOR_CO2, old_co2, text))
    out = {}
    for key, value in POSITIVE_ENV:
        if key in (ENV_NIGHT, ENV_FLOOR):
            old = environ.get(key)
            if old is None:
                environ[key] = text
            out[key] = environ[key]
            continue
        if key == ENV_FLOOR_CO2:
            if co2:
                environ[key] = "1"
                out[key] = "1"
            continue
        old = environ.get(key)
        if old is not None and old != value:
            raise SystemExit("positive_fit: %s=%s is preset; the user's decision (2026-11-07) is %s" % (key, old, value))
        environ[key] = value
        out[key] = value
    out["ratio_decision"] = RATIO_DECISIONS[r]
    return out


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

    Driver options (before '--'): --point-table PATH (required), --ratio R (night/shadow = floor ratio; also via the
    environment, see set_positive_environment), --start-restart PATH (start from another positive fit's
    fit_restart.npz; the prior is unchanged), --grain-start-uPa A, --grain-start-B B (grain B start in 1e-6 Pa per 100 s),
    --check-jacobian, --check-parameters a00,a02,..., --start-set NAME=VALUE[,...] (2026-11-08: needs --start-restart;
    START_SET_DEFINITION; absent: nothing changes). Grain start precedence: explicit --grain-start-uPa / -B, else
    the start restart's grain, else A = 6.96024299190473, B = 0. Everything else goes to grain_fit.py / gpu_fit.py
    (--crossing-weight, --devices, --precision, --max-nfev, --wall-time, --restart-from, ...).
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
           "ratio": _pop(head, "--ratio"),
           "start_restart": _pop(head, "--start-restart"),
           "grain_start_uPa": _pop(head, "--grain-start-uPa", float),
           "grain_start_B": _pop(head, "--grain-start-B", float),
           "grain_start_from": _pop(head, "--grain-start-from"),
           "held_exponent": _pop(head, "--held-exponent"),
           "start_set": _pop(head, "--start-set"),
           "check_jacobian": _flag(head, "--check-jacobian"),
           "check_held_span": _flag(head, "--check-held-span"),
           "check_parameters": list(CHECK_PARAMETERS)}
    names = _pop(head, "--check-parameters")
    if names:
        cfg["check_parameters"] = [n.strip() for n in names.split(",") if n.strip()]
    if not cfg["point_table"]:
        raise SystemExit("positive_fit: --point-table <run>/late_grid_points.txt is required")
    cfg["point_table"] = str(Path(cfg["point_table"]).resolve())
    if not Path(cfg["point_table"]).is_file():
        raise SystemExit("positive_fit: point table not found: %s" % cfg["point_table"])
    for key, opt in (("grain_start_uPa", "--grain-start-uPa"), ("grain_start_B", "--grain-start-B")):
        if cfg[key] is not None and not (math.isfinite(cfg[key]) and cfg[key] >= 0.0):
            raise SystemExit("positive_fit: %s must be finite and >= 0 (lower bound of the grain parameters)" % opt)
    if cfg["grain_start_from"] is not None:
        if cfg["grain_start_uPa"] is not None or cfg["grain_start_B"] is not None:
            raise SystemExit("positive_fit: --grain-start-from and --grain-start-uPa/-B exclude each other")
        A, B, ident = read_grain_start(cfg["grain_start_from"])
        cfg["grain_start_from"] = ident["path"]
        cfg["grain_start_uPa"], cfg["grain_start_B"] = A, B
        cfg["grain_start_uPa_source"] = cfg["grain_start_B_source"] = "%s (sha256 %s)" % (ident["path"], ident["sha256"])
        cfg["grain_start_from_identity"] = ident
    if cfg["held_exponent"] is not None:
        path = Path(cfg["held_exponent"]).resolve()
        if not path.is_file():
            raise SystemExit("positive_fit: --held-exponent file not found: %s" % path)
        cfg["held_exponent"] = str(path)
        cfg["held_identity"] = {"path": str(path), "sha256": sha256(path)}
        for opt, bad in (("--activity", "device"), ("--transport", "mpi")):
            v = _value(head, opt)
            if v is not None and (v == bad or (opt == "--transport" and v != "local")):
                raise SystemExit("positive_fit: --held-exponent needs the host activity of the local transport (the "
                                 "held term is added to act_lin on the host; %s %s computes it elsewhere)" % (opt, v))
    if cfg["check_held_span"] and cfg["held_exponent"] is None:
        raise SystemExit("positive_fit: --check-held-span needs --held-exponent")
    if cfg["ratio"] is not None:
        _ratio_value(cfg["ratio"], "--ratio")
    if cfg["start_restart"] is not None:
        cfg["start_restart"] = str(Path(cfg["start_restart"]).resolve())
        if not Path(cfg["start_restart"]).is_file():
            raise SystemExit("positive_fit: --start-restart file not found: %s" % cfg["start_restart"])
    if cfg["start_set"] is None:
        del cfg["start_set"]                         # absent option: the configuration of the earlier fits
    else:
        if cfg["start_restart"] is None:
            raise SystemExit("positive_fit: --start-set needs --start-restart (it overrides values of that start; a "
                             "resume continues its own state)")
        cfg["start_set"] = parse_start_set(cfg["start_set"])
        for name, _ in cfg["start_set"]:
            if name in ("grain_A_uPa", "grain_B_uPa_per_100s") and (cfg["grain_start_uPa"] is not None or
                                                                   cfg["grain_start_B"] is not None):
                raise SystemExit("positive_fit: --start-set %s and --grain-start-uPa/-B/-from both set the grain start"
                                 % name)
        notes.append("--start-set %s applied after --start-restart" % ",".join("%s=%r" % (n, v) for n, v in
                                                                             cfg["start_set"]))
    # gpu_fit.py options
    ts = _value(head, "--template-start")
    if ts is None:
        head += ["--template-start", "0"]
        notes.append("inserted gpu_fit --template-start 0 (start = the projection of the original map)")
    elif float(ts) != 0.0:
        raise SystemExit("positive_fit: --template-start must be 0: the template amplitude starts at the projected "
                         "value (user 2026-11-07); a nonzero value would replace it")
    cfg["restart_from"] = _value(head, "--restart-from")
    if cfg["restart_from"] is not None and cfg["start_restart"] is not None:
        raise SystemExit("positive_fit: --start-restart (a new fit from another fit's state) and gpu_fit --restart-from "
                         "(resume of this fit's previous segment) exclude each other")
    if cfg["restart_from"] is not None:
        cfg["restart_from"] = str(Path(cfg["restart_from"]).resolve())
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


def marker_ratios(ext):
    """(night, floor) ratios of a restart marker: the explicit fields (since 2026-11-07 evening), else the GPU patch
    configuration recorded in the marker (earlier restarts)."""
    if ext.get("night_ratio") is not None:
        return float(ext["night_ratio"]), float(ext["floor_ratio"])
    cfg = ((ext.get("positive_activity") or {}).get("config") or {})
    if cfg.get("night_ratio") is None:
        return None, None
    return float(cfg["night_ratio"]), float(cfg.get("floor_ratio"))


def marker_co2_floor(ext):
    """The CO2 flux floor of a restart marker: the explicit field (since 2026-11-08), else False (the earlier
    restarts were written before the CO2 floor existed)."""
    return bool(ext.get("co2_floor", False))


def marker_held(ext):
    """The held exponent term of a restart marker: None (no held term) or {path, sha256, ...}."""
    h = ext.get("held_exponent")
    return h if h else None


def check_restart_marker(path, cfg=None):
    """A restart (--restart-from, --start-restart) must come from this parameterization (a linear fit's restart would
    be misread). For --restart-from (resume of this fit) the ratio and the held term must be this process's."""
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
    ext = meta.get("model_extension") or {}
    if ext.get("parameterization") != PARAMETERIZATION_ID:
        raise SystemExit("positive_fit: %s is not a restart of a positive_fit.py fit (meta model_extension."
                         "parameterization %r, need %r)" % (path, ext.get("parameterization"), PARAMETERIZATION_ID))
    night, floor = marker_ratios(ext)
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "parameterization": ext.get("parameterization"),
           "night_ratio": night, "floor_ratio": floor, "co2_floor": marker_co2_floor(ext),
           "held_exponent": marker_held(ext),
           "lambda": (meta.get("settings") or {}).get("lambda"), "status": meta.get("status"),
           "segments": meta.get("segments"), "fit_directory": meta.get("fit_directory"),
           "fit_chain": ext.get("fit_chain")}
    if cfg is not None and cfg.get("restart_from") and str(Path(path).resolve()) == cfg["restart_from"]:
        r = current_ratio()
        if night != r or floor != r:
            raise SystemExit("positive_fit: --restart-from %s was written with night/floor ratios %s/%s; this process "
                             "has %s: a resume must continue the same problem" % (path, night, floor, r))
        if marker_co2_floor(ext) != current_co2_floor():
            raise SystemExit("positive_fit: --restart-from %s was written with the CO2 floor %s; this process has %s: a "
                             "resume must continue the same problem" % (path, marker_co2_floor(ext),
                                                                        current_co2_floor()))
        held = held_identity(cfg)
        if (marker_held(ext) or {}).get("sha256") != (held or {}).get("sha256"):
            raise SystemExit("positive_fit: --restart-from %s has the held term %s; this process has %s: a resume "
                             "must continue the same problem" % (path, marker_held(ext), held))
    return out


def held_identity(cfg):
    """{path, sha256} of the held exponent file of this process (None: no held term)."""
    return cfg.get("held_identity")


def read_grain_start(path):
    """(A [1e-6 Pa], B [1e-6 Pa per 100 s], identity) from a fine tuning's state file (keys grain_A_uPa,
    grain_B_uPa_per_100s)."""
    path = Path(path).resolve()
    if not path.is_file():
        raise SystemExit("positive_fit: --grain-start-from file not found: %s" % path)
    with np.load(str(path), allow_pickle=False) as z:
        missing = [k for k in ("grain_A_uPa", "grain_B_uPa_per_100s") if k not in z.files]
        if missing:
            raise SystemExit("positive_fit: --grain-start-from %s lacks %s" % (path, missing))
        A, B = (float(np.asarray(z[k]).reshape(-1)[0]) if np.asarray(z[k]).size == 1 else float("nan")
                for k in ("grain_A_uPa", "grain_B_uPa_per_100s"))
    for name, v in (("grain_A_uPa", A), ("grain_B_uPa_per_100s", B)):
        if not (math.isfinite(v) and v >= 0.0):
            raise SystemExit("positive_fit: --grain-start-from %s: %s = %r (need one finite value >= 0)" % (path, name, v))
    return A, B, {"path": str(path), "sha256": sha256(path), "grain_A_uPa": A, "grain_B_uPa_per_100s": B}


# ------------------------------------------------------------------------------------------------ held exponent
# The per-facet term h held fixed in the exponent (user 2026-11-07: cycles 2 and 3 hold the accumulated fine-tuning
# change): s = Y a + c T + h. Loaded once per process by positive_fit.py (load_held + install_held_activity, before any
# backend binds); the Problem, the writers and the restart marker read HELD.
HELD = {"h": None, "info": None, "calls": 0}
HELD_KEYS = ("held_next", "held")       # held_next: written by the fine tuning; held: tests


def elements_faces(run):
    """Number of facets of <run>/elements.bin (header int32 [version, nfaces, ny])."""
    header = np.fromfile(str(Path(run) / "elements.bin"), dtype="<i4", count=3)
    return int(header[1])


def load_held(path, nfaces):
    """h (float64, one value per facet in mesh order, nepers) from an .npz with exactly one of the keys held_next /
    held. Refused: another dtype, shape or a non-finite value. A mesh identity stored in the file (key mesh_identity, or
    meta_json mesh_identity / mesh.identity) is compared with the base's when the Problem is built."""
    path = Path(path).resolve()
    with np.load(str(path), allow_pickle=False) as z:
        keys = [k for k in HELD_KEYS if k in z.files]
        if len(keys) != 1:
            raise SystemExit("positive_fit: --held-exponent %s must have exactly one of the keys %s (has %s)"
                             % (path, HELD_KEYS, z.files))
        h = np.array(z[keys[0]])
        meta = {}
        if "meta_json" in z.files:
            try:
                meta = json.loads(str(z["meta_json"]))
            except ValueError:
                meta = {"unreadable_meta_json": True}
        mesh = str(z["mesh_identity"]) if "mesh_identity" in z.files else (
            meta.get("mesh_identity") or (meta.get("mesh") or {}).get("identity") if isinstance(meta, dict) else None)
    if h.dtype != np.float64 or h.shape != (nfaces,):
        raise SystemExit("positive_fit: --held-exponent %s key %s is %s %s; need float64 (%d,) in mesh order"
                         % (path, keys[0], h.dtype, h.shape, nfaces))
    if not np.all(np.isfinite(h)):
        raise SystemExit("positive_fit: --held-exponent %s has non-finite values" % path)
    h.setflags(write=False)
    info = {"path": str(path), "sha256": sha256(path), "key": keys[0], "facets": int(nfaces),
            "array_sha256": hashlib.sha256(h.tobytes()).hexdigest(), "mesh_identity_in_file": mesh,
            "min": float(h.min()), "max": float(h.max()), "nonzero_facets": int(np.count_nonzero(h)),
            "file_meta": jsonable({k: v for k, v in meta.items()
                                   if not isinstance(v, (list, dict)) or len(json.dumps(v)) < 4000})
            if isinstance(meta, dict) else None}
    HELD["h"], HELD["info"] = h, info
    return h, info


def install_held_activity(record):
    """GPU activity with the held term: act_lin = Y a + A_scale B b + h on the host (CudaBackend.activity_linear and
    MultiGpuBackend.activity_linear of THIS process; the local multi-GPU workers read the parent's act_lin from shared
    memory, so they see h without a patch of their own). The device activity (MPI transport or --activity device)
    computes act_lin on the GPU without h: refused at the first evaluation. The Jacobian columns (d af / d coef = af
    basis, af = exp(act_lin)) need no change: h is a constant. Called after positive_activity.install()."""
    from cops_cuda import evaluator, multi_gpu
    from cops_cuda.errors import UnsupportedConfiguration
    h = HELD["h"]
    if h is None:
        return False
    B, M = evaluator.CudaBackend, multi_gpu.MultiGpuBackend
    if getattr(B.activity_linear, "_positive_held", False):
        return True
    orig_cuda, orig_multi = B.activity_linear, M.activity_linear
    orig_upload, orig_evaluate = B._upload_activity, M._evaluate

    def with_held(act):
        h_now = HELD["h"]                  # None only inside check_held_span (the shifted comparison)
        if h_now is None:
            return act
        if act.shape != h_now.shape:
            raise UnsupportedConfiguration("positive_problem: act_lin has %s values, the held term %s"
                                           % (act.shape, h_now.shape))
        HELD["calls"] += 1
        return act + h_now

    def activity_linear_cuda(self, state):
        return with_held(orig_cuda(self, state))

    def activity_linear_multi(self, state):
        return with_held(orig_multi(self, state))

    def _upload_activity(self, state):
        if self.activity_mode == "device":
            raise UnsupportedConfiguration("positive_problem: the held exponent term needs the host activity")
        return orig_upload(self, state)

    def _evaluate(self, state, request):
        if self.activity != "host":
            raise UnsupportedConfiguration("positive_problem: the held exponent term needs the host activity "
                                           "(multi-GPU activity %r)" % self.activity)
        return orig_evaluate(self, state, request)
    for f in (activity_linear_cuda, activity_linear_multi, _upload_activity, _evaluate):
        f._positive_held = True
    B.activity_linear, M.activity_linear = activity_linear_cuda, activity_linear_multi
    B._upload_activity, M._evaluate = _upload_activity, _evaluate
    record["held_gpu_patch"] = {
        "patched": ["evaluator.CudaBackend.activity_linear (+h)", "multi_gpu.MultiGpuBackend.activity_linear (+h)",
                    "CudaBackend._upload_activity (refuses device activity)",
                    "MultiGpuBackend._evaluate (refuses device activity)"],
        "process": os.getpid(), "workers": "read the parent's act_lin (shared memory; host activity)"}
    return True


def held_projection(Y, T, area, h, G_all, ascale, chunk=1 << 20):
    """h = B h_B + h_perp, area-weighted LS projection on the fitted columns B = [Y, T] (physical coordinates; h_B is
    returned in optimizer coordinates, i.e. divided by S = (ascale (Y), 1 (templates))). |h_perp|^2_area (the constant
    part of D^2), |h|^2_area, and the orthogonality residual max|B^T diag(area) h_perp| / max|B^T diag(area) h|."""
    n, ny = Y.shape
    m = G_all.shape[0]
    rhs = np.zeros(m)
    for i in range(0, n, chunk):
        sl = slice(i, min(i + chunk, n))
        Bm = np.column_stack((Y[sl], T[sl]))
        rhs += Bm.T @ (area[sl] * h[sl])
    z = np.linalg.solve(G_all, rhs)
    perp2, h2, ortho = 0.0, 0.0, np.zeros(m)
    for i in range(0, n, chunk):
        sl = slice(i, min(i + chunk, n))
        Bm = np.column_stack((Y[sl], T[sl]))
        hp = h[sl] - Bm @ z
        perp2 += float(np.sum(area[sl] * hp * hp))
        h2 += float(np.sum(area[sl] * h[sl] * h[sl]))
        ortho += Bm.T @ (area[sl] * hp)
    total = float(np.sum(area))
    S = np.concatenate([np.full(ny, ascale), np.ones(m - ny)])
    scale = float(np.max(np.abs(rhs)))
    return {"hB": z / S, "z_physical": z, "perp2": perp2 / total, "h2": h2 / total,
            "orthogonality_rel": float(np.max(np.abs(ortho)) / scale) if scale > 0 else 0.0}


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
    record.update(policy=current_policy(), parameterization=PARAMETERIZATION_ID,
                  ratio={"night_ratio": current_ratio(), "floor_ratio": float(os.environ.get(ENV_FLOOR, "nan")),
                         "decision": RATIO_DECISIONS.get(current_ratio()), "co2_floor": current_co2_floor(),
                         "co2_floor_ratio": (float(os.environ.get(ENV_FLOOR, "nan")) if current_co2_floor() else None),
                         "co2_floor_definition": (CO2_FLOOR_TEXT.format(floor=repr(current_ratio()))
                                                  if current_co2_floor() else "none (CO2 has no floor)")},
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
            self.held = HELD["h"]
            self.held_hB, self.held_perp2, self.held_row = None, 0.0, False
            if self.held is not None:
                self._install_held(el, T, grams[0], record)
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

        # ---- held exponent term h (user 2026-11-07): s = B x + h, D^2 = |B (x - x0) + h|^2_area
        #      = |Lt (x - x0 + h_B)|^2 + |h_perp|^2_area  (h = B h_B + h_perp, h_perp area-orthogonal to B)
        def _install_held(self, el, T, G_all, record):
            from cops_cuda import reference
            h, info = self.held, dict(HELD["info"])
            if h.shape != (len(el["area"]),):
                raise SystemExit("positive_problem: the held term has %s values for %d facets" % (h.shape, len(el["area"])))
            mesh = reference.load("activity_templates").mesh_identity(el["centroid"], el["area"])
            if info.get("mesh_identity_in_file") not in (None, mesh):
                raise SystemExit("positive_problem: the held term belongs to the mesh %s, the base is %s"
                                 % (info["mesh_identity_in_file"], mesh))
            info["mesh_identity_verified"] = info.get("mesh_identity_in_file") == mesh
            pr = held_projection(el["Y"], T, el["area"], h, G_all, self.ascale)
            self.held_hB = pr["hB"]
            self.held_perp2 = pr["perp2"]
            self.held_row = pr["perp2"] > 0.0
            info.update(rms_area_nepers=math.sqrt(pr["h2"]), D_perp_nepers=math.sqrt(pr["perp2"]),
                        D_perp_squared=pr["perp2"], h_B_optimizer=[float(v) for v in pr["hB"]],
                        h_B_physical=[float(v) for v in pr["z_physical"]], orthogonality_rel=pr["orthogonality_rel"],
                        regularization_rows="Lt (x - x0 + h_B)" + (" and the constant sqrt(lambda) |h_perp|_area"
                                                                   if self.held_row else ""),
                        D_at_prior_x0_nepers=math.sqrt(pr["h2"]))
            record["held_exponent"] = info
            record["held_exponent_marker"] = {k: info[k] for k in ("path", "sha256", "key", "array_sha256",
                                                                   "D_perp_nepers", "rms_area_nepers")}
            print("positive_fit: held exponent term %s (key %s, sha256 %s): rms_area %.6f nepers, %d nonzero facets; "
                  "h = B h_B + h_perp with |h_perp|_area = %.6f nepers (constant part of D), orthogonality %.1e; "
                  "s = Y a + c T + h" % (info["path"], info["key"], info["sha256"][:16], info["rms_area_nepers"],
                                         info["nonzero_facets"], info["D_perp_nepers"], pr["orthogonality_rel"]),
                  flush=True)

        def regularization(self, x):
            if self.held is None:
                return super().regularization(x)
            rs = self.reg_slice
            r = np.sqrt(self.lam) * (self.Lt @ (x[rs] - self.x0[rs] + self.held_hB))
            if self.held_row:
                r = np.concatenate([r, [np.sqrt(self.lam) * math.sqrt(self.held_perp2)]])
            return r

        def reg_jacobian(self):
            R = super().reg_jacobian()
            if self.held is not None and self.held_row:
                R = np.vstack([R, np.zeros((1, R.shape[1]))])
            return R

        def deviation(self, x):
            if self.held is None:
                return super().deviation(x)
            rs = self.reg_slice
            v = self.Lt @ (x[rs] - self.x0[rs] + self.held_hB)
            return float(np.sqrt(float(v @ v) + self.held_perp2))

    fit_cops.Problem = PositiveCore
    install_activity_writers(fit_cops, record, gkeys)
    if cfg.get("check_jacobian"):
        install_jacobian_check(fit_cops, cfg["check_parameters"], record)
    if cfg.get("check_held_span"):
        install_held_span_check(fit_cops, record)
    source = cfg.get("start_restart") or cfg.get("restart_from")
    if source:                                       # outside the Jacobian check: runs first
        role = "start_restart" if cfg.get("start_restart") else "restart_from"
        info, arrays = load_source(source, role)
        record["fit_chain"] = fit_chain(info, role)
        install_start_check(fit_cops, cfg, record, info, arrays)
    else:
        record["fit_chain"] = {"start": "the projection of the original map (no restart)"}
    return PositiveCore


def install_held_span_check(fit_cops, record):
    """Diagnostic (--check-held-span, tests): with a held term in the span of the fitted columns (h = B z), the data
    residuals and the data Jacobian at the optimizer start x must equal those WITHOUT the held term at x + h_B (the GPU
    evaluated both ways: act_lin + h on the host vs the shifted coefficients), and D(x) must equal the D of the shifted
    state. Once, at the first least_squares call; printed and recorded (model_extension.positive.held_span_check)."""
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"]:
            state["done"] = True
            P = fit_cops.P
            out = {"perp_over_h_rms": math.sqrt(P.held_perp2) / max(math.sqrt(sum(
                np.asarray(P.held, float) ** 2) / len(P.held)), 1e-300)}
            x = np.array(x0, float)
            rH, JH, DH = P.data_residuals(x), P.data_jacobian(x), P.deviation(x)
            calls = HELD["calls"]
            xs = x.copy()
            xs[P.reg_slice] += P.held_hB
            saved = HELD["h"], P.held
            HELD["h"], P.held = None, None
            try:
                rN, JN, DN = P.data_residuals(xs), P.data_jacobian(xs), P.deviation(xs)
            finally:
                HELD["h"], P.held = saved
            scale = float(np.max(np.abs(rN)))
            out.update(residual_max_abs_diff=float(np.max(np.abs(rH - rN))), residual_scale=scale,
                       jacobian_max_rel_diff=float(np.max(np.abs(JH - JN)) / np.max(np.abs(JN))),
                       D_held=float(DH), D_shifted=float(DN), held_calls_with=int(calls),
                       held_calls_without=int(HELD["calls"] - calls),
                       h_B=[float(v) for v in P.held_hB])
            record["held_span_check"] = out
            print("positive_fit: held span check (h in span(B)): data residuals with h at x vs without h at x + h_B: "
                  "max |diff| %.3e (scale %.3e); data Jacobian max rel diff %.3e; D %.12f vs %.12f; |h_perp|/rms(h) "
                  "%.1e" % (out["residual_max_abs_diff"], scale, out["jacobian_max_rel_diff"], DH, DN,
                            out["perp_over_h_rms"]), flush=True)
            fit_cops.save_active_manifest()
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


# ------------------------------------------------------------------------------------------------ start restart
def load_source(path, role):
    """What the start check needs from the restart that supplies the start (role start_restart: another fit's state;
    restart_from: this fit's previous segment): marker, final cost, lambda, the source fit's manifest numbers (data
    cost, regularization cost, D, crossing penalty), stored predictions (COPS_model.dat, 7 significant digits)."""
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        arrays = {k: np.array(z[k]) for k in ("point_ids", "pred_ng", "pred_rg", "x")}
    ext = meta.get("model_extension") or {}
    if ext.get("parameterization") != PARAMETERIZATION_ID:
        raise SystemExit("positive_fit: %s is not a restart of a positive_fit.py fit (parameterization %r)"
                         % (path, ext.get("parameterization")))
    night, floor = marker_ratios(ext)
    opt = meta.get("optimizer") or {}
    info = {"role": role, "path": str(Path(path).resolve()), "sha256": sha256(path),
            "cost": float(opt["cost"]), "lambda": float((meta.get("settings") or {}).get("lambda")),
            "nfev": opt.get("nfev"), "message": opt.get("message"), "status": meta.get("status"),
            "segments": meta.get("segments"), "night_ratio": night, "floor_ratio": floor,
            "co2_floor": marker_co2_floor(ext),
            "held_exponent": marker_held(ext), "fit_directory": meta.get("fit_directory"),
            "previous_fit_chain": ext.get("fit_chain")}
    manifest = Path(meta.get("fit_directory") or "/nonexistent") / "fit_manifest.json"
    try:
        m = json.loads(manifest.read_text())
        o = [v for v in m["optimizer_results"] if float(v["lambda"]) == info["lambda"]][0]
        info.update(manifest=str(manifest), data_cost=float(o["data_cost"]),
                    regularization_cost=float(o["regularization_cost"]), D=float(o["deviation_from_prior"]),
                    penalty_cost=float(m["model_extension"]["final"]["penalty_cost"]),
                    manifest_cost=float(o["cost"]))
    except Exception as exc:                          # recorded and printed, never silent
        info["manifest_error"] = "%s: %s" % (exc.__class__.__name__, exc)
    return info, arrays


def fit_chain(info, role):
    keep = ("path", "sha256", "fit_directory", "lambda", "cost", "night_ratio", "floor_ratio", "co2_floor",
            "held_exponent", "status", "nfev")
    entry = {k: info.get(k) for k in keep}
    if role == "start_restart":
        return {"start": "another fit's restart (positive_fit.py --start-restart)", "start_restart": entry}
    prev = dict(info.get("previous_fit_chain") or {})
    prev["resumed_from"] = list(prev.get("resumed_from") or []) + [entry]
    return prev


def start_from_restart(P, el, path, record):
    """Optimizer start of P from another positive fit's restart (cops_cuda.restart.start_vector: by name; with the same
    ascale the exact optimizer coordinates). Refused: another parameterization, parameter layout, mesh, template or
    temperature law. The prior P.x0 is not changed."""
    from cops_cuda import reference
    from cops_cuda import restart as rs
    state = rs.read(path)
    ext = state.meta.get("model_extension") or {}
    if ext.get("parameterization") != PARAMETERIZATION_ID:
        raise SystemExit("positive_fit: --start-restart %s is not a positive fit (%r)" % (path, ext.get("parameterization")))
    if list(state.meta.get("names") or []) != list(P.names):
        raise SystemExit("positive_fit: --start-restart %s has the parameters %s, this fit %s"
                         % (path, state.meta.get("names"), list(P.names)))
    templates = reference.load("activity_templates")
    if templates.mesh_identity(el["centroid"], el["area"]) != state.mesh_identity():
        raise SystemExit("positive_fit: --start-restart %s belongs to another mesh" % path)
    cols = np.asarray(P.enhancement_basis[:, P.nh:], float)
    if state.template_names != list(P.template_names) or not np.array_equal(state.template_columns, cols):
        raise SystemExit("positive_fit: --start-restart %s has other template columns (%s)" % (path, state.template_names))
    if not np.allclose(np.asarray(state.T, float), USER_T0, rtol=0.0, atol=1e-9):
        raise SystemExit("positive_fit: --start-restart %s has the temperature law %s, the user's is %s"
                         % (path, list(state.T), list(USER_T0)))
    x, report = rs.start_vector(P, state)
    same_scale = float(state.meta.get("ascale", np.nan)) == float(P.ascale)
    out = {"path": str(Path(path).resolve()), "sha256": sha256(path), "start_vector": jsonable(report),
           "same_ascale": same_scale, "x_equals_source_x": bool(np.array_equal(x, np.asarray(state.x, float))),
           "source_x": [float(v) for v in state.x], "prior_x0_unchanged": True}
    if report["clipped"] or report["kept_problem_start"]:
        print("positive_fit: WARNING start restart: clipped %s, kept the problem's start for %s"
              % (report["clipped"], report["kept_problem_start"]), flush=True)
    del state
    record["start_restart"] = out
    return x


def install_start_check(fit_cops, cfg, record, source, arrays):
    """At the first least_squares call (inside the crossing wrapper: fun = data + Tikhonov + crossing rows): the cost
    at the optimizer start against the source restart's final cost. When the source is the same problem (same ratio,
    lambda and held term) the start must reproduce it (|rel| <= START_REPRODUCTION_RTOL; otherwise the fit stops);
    otherwise the differences are reported with their causes (printed, manifest model_extension.positive.start_check)."""
    original = fit_cops.least_squares
    state = {"done": False}

    def least_squares(fun, x0, jac="2-point", **kwargs):
        if not state["done"]:
            state["done"] = True
            res = start_check(fun, x0, fit_cops, cfg, source, arrays)
            record["start_check"] = res
            fit_cops.save_active_manifest()
            if res.get("failed"):
                raise SystemExit("positive_fit: start check failed: %s" % res["failed"])
        return original(fun, x0, jac=jac, **kwargs)
    fit_cops.least_squares = least_squares


def start_check(fun, x0, fit_cops, cfg, source, arrays):
    P = fit_cops.P
    x = np.array(x0, dtype=float)
    r = np.asarray(fun(x), float)
    nd = int(P.ndata)
    nreg = len(P.regularization(x)) if P.lam > 0 else 0
    rd, rr, rc = r[:nd], r[nd:nd + nreg], r[nd + nreg:]
    data, reg, pen = (0.5 * float(v @ v) for v in (rd, rr, rc))
    cost, D = data + reg, float(P.deviation(x))
    held_now = (held_identity(cfg) or {}).get("sha256")
    held_src = (source.get("held_exponent") or {}).get("sha256")
    ratio = current_ratio()
    floor = float(os.environ.get(ENV_FLOOR, "nan"))
    co2 = current_co2_floor()
    co2_src = bool(source.get("co2_floor", False))
    same = {"ratio": source["night_ratio"] == ratio and source["floor_ratio"] == floor and co2_src == co2,
            "lambda": source["lambda"] == float(P.lam), "held_exponent": held_now == held_src,
            "x": len(arrays["x"]) == len(x) and bool(np.array_equal(x, arrays["x"]))}
    sset = cfg.get("start_set") or []
    if sset:                                          # --start-set (2026-11-08): never the source's state
        same["start_set"] = False
    expect = all(same.values())
    out = {"source": {k: v for k, v in source.items() if k != "previous_fit_chain"},
           "start": {"cost": cost, "data_cost": data, "regularization_cost": reg, "D": D, "crossing_penalty_cost": pen,
                     "night_ratio": ratio, "floor_ratio": floor, "co2_floor": co2, "lambda": float(P.lam),
                     "held_sha256": held_now},
           "same": same, "reproduction_expected": expect,
           "cost_rel_diff": cost / source["cost"] - 1.0,
           "x_start_equals_source_x": same["x"], "tolerance": START_REPRODUCTION_RTOL}
    if "data_cost" in source:
        out["data_cost_rel_diff"] = data / source["data_cost"] - 1.0
        out["D_diff"] = D - source["D"]
        out["crossing_penalty_diff"] = pen - source["penalty_cost"]
    # model at the start vs the source's stored predictions (used rows; ln p = r / w + ln ref)
    try:
        ip = np.asarray(P.M.ipoint)[np.asarray(P.M.rows)]
        pos = {int(p): k for k, p in enumerate(arrays["point_ids"])}

        def dln(use, w, ln_ref, res, stored):
            idx = np.array([pos[int(p)] for p in ip[use]], dtype=np.int64)
            d = (res / w + ln_ref) - np.log(np.maximum(stored[idx], fit_cops.FLOOR))
            k = int(np.argmax(np.abs(d)))
            return {"points": int(len(d)), "max_abs_dln": float(abs(d[k])), "at_point": int(ip[use][k]),
                    "rms_dln": float(np.sqrt(np.mean(d * d))), "mean_dln": float(np.mean(d))}
        out["predictions_vs_source"] = {"NG": dln(P.useNG, P.wNG, P.lnNG, rd[:P.NNG], arrays["pred_ng"]),
                                        "RG": dln(P.useRG, P.wRG, P.lnRG, rd[P.NNG:], arrays["pred_rg"]),
                                        "note": "source = its COPS_model.dat (%.6e): |d ln p| ~ 5e-7 is rounding"}
    except Exception as exc:
        out["predictions_vs_source"] = {"error": "%s: %s" % (exc.__class__.__name__, exc)}
    changes = []
    if not same["ratio"]:
        changes.append("ratio %s/%s -> %s/%s" % (source["night_ratio"], source["floor_ratio"], ratio, floor)
                       + ("; CO2 floor %s -> %s" % ("on" if co2_src else "off", "on" if co2 else "off")
                          if co2_src != co2 else ""))
    if not same["lambda"]:
        changes.append("lambda %g -> %g" % (source["lambda"], P.lam))
    if not same["held_exponent"]:
        changes.append("held term %s -> %s" % ((held_src or "none")[:12], (held_now or "none")[:12]))
    if sset:
        names = list(P.names)
        out["start_set"] = [{"name": n, "value": float(v), "start_value": float(x[names.index(n)]),
                             "source_value": (float(arrays["x"][names.index(n)])
                                              if len(arrays["x"]) == len(x) else None)} for n, v in sset]
        changes.append("--start-set %s (the source's %s)" % (
            ", ".join("%s=%r" % (e["name"], e["value"]) for e in out["start_set"]),
            ", ".join("%s=%r" % (e["name"], e["source_value"]) for e in out["start_set"])))
    if not same["x"]:
        diff = ([P.names[k] for k in np.nonzero(x != arrays["x"])[0]] if len(arrays["x"]) == len(x)
                else ["parameter layout"])
        changes.append("start x differs from the source's x in %s" % ", ".join(diff))
    out["changes"] = changes
    if expect and not abs(out["cost_rel_diff"]) <= START_REPRODUCTION_RTOL:
        out["failed"] = ("the start reproduces the source's problem (same ratio, lambda, held term) but its cost %.12e "
                         "differs from the source's final cost %.12e by %.3e (relative; tolerance %.0e)"
                         % (cost, source["cost"], out["cost_rel_diff"], START_REPRODUCTION_RTOL))
    print("positive_fit: start check (%s %s): cost at the start %.12e (data %.12e, regularization %.6e, D %.9f; "
          "crossing penalty %.3e); source final cost %.12e (data %s, D %s, penalty %s); relative difference %+.3e%s" % (
              source["role"], source["path"], cost, data, reg, D, pen, source["cost"],
              "%.12e" % source["data_cost"] if "data_cost" in source else "?",
              "%.9f" % source["D"] if "D" in source else "?",
              "%.3e" % source["penalty_cost"] if "penalty_cost" in source else "?", out["cost_rel_diff"],
              " (data %+.3e)" % out["data_cost_rel_diff"] if "data_cost_rel_diff" in out else ""), flush=True)
    pv = out["predictions_vs_source"]
    if "NG" in pv:
        print("positive_fit: start check: model at the start vs the source's stored predictions: NG max |d ln p| %.3e "
              "(rms %.3e, mean %+.3e), RG max %.3e (rms %.3e, mean %+.3e)" % (
                  pv["NG"]["max_abs_dln"], pv["NG"]["rms_dln"], pv["NG"]["mean_dln"], pv["RG"]["max_abs_dln"],
                  pv["RG"]["rms_dln"], pv["RG"]["mean_dln"]), flush=True)
    else:
        print("positive_fit: WARNING start check: prediction comparison unavailable: %s" % pv.get("error"), flush=True)
    if "manifest_error" in source:
        print("positive_fit: WARNING start check: the source's fit_manifest.json could not be read: %s"
              % source["manifest_error"], flush=True)
    if sset:
        print("positive_fit: start check: the start is NOT the source's state: --start-set %s"
              % "; ".join("%s = %r (source %r)" % (e["name"], e["value"], e["source_value"]) for e in out["start_set"]),
              flush=True)
    if expect:
        print("positive_fit: start check: same problem and state as the source (ratio, lambda, held term, x): "
              "reproduction %s "
              "(|rel| %.1e, tolerance %.0e)" % ("CONFIRMED" if "failed" not in out else "FAILED",
                                               abs(out["cost_rel_diff"]), START_REPRODUCTION_RTOL), flush=True)
    else:
        print("positive_fit: start check: NOT the source's problem (%s): the start cost differs from the source's "
              "final cost because of these changes (expected, recorded)" % "; ".join(changes), flush=True)
    return out


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
            if names[k + 1] != "grain_B_uPa_per_100s":
                raise SystemExit("positive_problem: unexpected grain parameter layout %s" % names[k:])
            default = float(self.x0[k])
            start = start_from_restart(self, el, cfg["start_restart"], record) if cfg.get("start_restart") else None
            # grain start: explicit (fine tuning) > the start restart's > A = cycle 3's fine-tuned value, B = 0
            values, sources = [], []
            for j, key, fallback in ((k, "grain_start_uPa", DEFAULT_GRAIN_START_UPA), (k + 1, "grain_start_B", 0.0)):
                if cfg.get(key) is not None:
                    values.append(float(cfg[key]))
                    sources.append(cfg.get(key + "_source") or "explicit (--%s)" % key.replace("_", "-"))
                elif start is not None:
                    values.append(float(start[j]))
                    sources.append("the start restart %s" % cfg["start_restart"])
                else:
                    values.append(float(fallback))
                    sources.append("user 2026-11-07: cycle 3 fine-tuned A_Pa (grain/cycles/c3/finetune/C30_n4000_"
                                   "grain_cross10_ftc3_m5_i5_cpu/grain_finetune.json)" if j == k else
                                   "user 2026-11-07: B = 0")
            self.x0[k], self.x0[k + 1] = values
            if start is not None:
                start[k], start[k + 1] = values
                if cfg.get("start_set"):                   # 2026-11-08: explicit start values after the restart's
                    apply_start_set(self, start, cfg["start_set"], record)
                self.restart_start = start                 # fit_cops.fit_one's start; the prior x0 is unchanged
            record["grain_start"] = {"grain_A_uPa": values[0], "grain_B_uPa_per_100s": values[1],
                                     "source_A": sources[0], "source_B": sources[1],
                                     "replaced_grain_crossing_default_A_uPa": default,
                                     "note": "with gpu_fit --restart-from (resume) the previous segment's values "
                                             "replace this start (gpu_fit install_template_start)"}
            print("positive_fit: grain start A %.9g uPa (%s), B %.9g uPa/100 s (%s)%s" % (
                values[0], sources[0], values[1], sources[1],
                "; optimizer start = the start restart's x (prior unchanged)" if start is not None else ""), flush=True)
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


def apply_start_set(P, start, assignments, record):
    """--start-set: start[name] = value for each assignment (optimizer coordinates, START_SET_DEFINITION); refused for
    names that are not parameters of P or values outside [lo, hi]. Recorded in record["start_set"] (manifest) and in
    record["fit_chain"]["start_set"] (restart marker, resumes, ACTIVITY_PARAMETERIZATION.txt)."""
    names = list(P.names)
    applied = []
    for name, value in assignments:
        if name not in names:
            raise SystemExit("positive_fit: --start-set %s: not a parameter of this fit (%s)" % (name, ", ".join(names)))
        j = names.index(name)
        lo, hi = float(P.lo[j]), float(P.hi[j])
        if not lo <= float(value) <= hi:
            raise SystemExit("positive_fit: --start-set %s=%r is outside its bounds [%r, %r]" % (name, value, lo, hi))
        applied.append({"name": name, "index": j, "value": float(value), "start_restart_value": float(start[j]),
                        "bounds": [lo, hi], "decision": START_SET_DECISIONS.get(name, "no recorded user decision for "
                                                                                     "this name")})
        start[j] = float(value)
    record["start_set"] = {"assignments": applied, "definition": START_SET_DEFINITION,
                           "applied_after": "--start-restart and the grain start; the prior x0 is unchanged"}
    if isinstance(record.get("fit_chain"), dict):
        record["fit_chain"]["start_set"] = applied
    print("positive_fit: --start-set: %s; the optimizer start is NOT the start restart's state in these parameters"
          % ", ".join("%s = %r (the start restart's %r)" % (a["name"], a["value"], a["start_restart_value"])
                      for a in applied), flush=True)
    return applied


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
    """The exponent field s = Y a + enhancement_scale * B b (+ the held term h; enhancement_scale = 1 here)."""
    s = exponent_field(el["Y"], np.asarray(a, float))
    for j in range(P.nb):
        s += P.enhancement_scale * float(b[j]) * P.enhancement_basis[:, j]
    h = getattr(P, "held", None)
    if h is not None:
        s += h
    return s


def install_activity_writers(fit_cops, record, gkeys=frozenset()):
    """activity_diagnostics and the A5 activity products in exp semantics (see the module docstring); the A5 pressure
    products (pressure_residuals.dat/png, pressure_comparison_*.png) for the observation rows only (no grid rows)."""

    def activity_diagnostics(a, b, el, seen, coef, enhancement_basis, enhancement_scale):
        held = getattr(fit_cops.P, "held", None)
        s = el["Y"] @ a
        if held is not None:
            s = s + held
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
            delta = act - np.exp(el["Y"] @ a + (held if held is not None else 0.0))   # added by the template factor(s)
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
the floor do not depend on the activity parameters; {co2_text}). Temperature law held: Tmax {T0[0]},
T75 {T0[1]}, Tmin {T0[2]} K. Late window (>= {late} UTC): model = W10x5 of the surface model (+ grain at the time stamp).
Tikhonov: D^2 = sum area (s - s0)^2 / sum area (nepers; s0 = the start and prior = the emission-weighted, w = area A_orig,
LS projection of ln A_orig on the same basis). lambda = {lam:g}.
Ratio (night/shadow = H2O floor{co2_ratio}): {night} -- {ratio_decision}.
{co2_par}Optimizer start: {start}
Held exponent term: {held}

VALID in this directory:
  activity_exp_facets.npz       s (includes the held term), A = exp(s), s0, held (h; zeros without one),
                                F_min_over_qmax = {floor} A_orig, seen  (per facet, mesh order); meta: D direct, D fitter
{co2_npz}  activity_facets.dat           legacy_linear = Y coef[0] of the run; residual_linear = A - legacy_linear;
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
  fit_restart.npz "activity_linear" (= s INCLUDING the held term, not a linear activity): tools that use
  max(activity_linear, 0) (fine_tune_gpu.py, fine_tuning.py, restart_tool.py, grain_finetune.py, late_response_eval.py)
  are not valid for this fit. The restart file can continue this fit with positive_fit.py --restart-from (marker meta
  "model_extension": ratio, held term path/sha256, fit chain) or start another with --start-restart.
"""


D_CHECK_RTOL = 1e-9     # D direct (sum over facets) vs the fitter's D (Lt and the held decomposition)


def write_exp_products(out, P, r, el, seen, record):
    out = Path(out)
    a, b, T, kappa = P.unpack(r["x"])
    s = exp_state(P, el, a, b)
    A = np.exp(s)
    floor = float(os.environ.get("COPS_POS_FLOOR_RATIO", "nan"))
    night = os.environ.get("COPS_POS_NIGHT_RATIO", "?")
    held = getattr(P, "held", None)
    D = float(np.sqrt(np.sum(el["area"] * (s - P.s0) ** 2) / P.area_total))
    D_fit = float(r["D"])
    D_rel = abs(D - D_fit) / max(D, D_fit) if max(D, D_fit) > 0 else 0.0
    hinfo = record.get("held_exponent")
    meta = {"parameterization": PARAMETERIZATION_ID, "created_utc": utc_now(), "lambda": float(r["lam"]),
            "names": list(P.names), "x": [float(v) for v in r["x"]], "ascale": float(P.ascale),
            "floor_ratio": floor, "night_ratio": night, "templates": list(P.template_names),
            "run": str(getattr(P.args, "run", "")), "ratio_decision": RATIO_DECISIONS.get(current_ratio()),
            "fit_chain": record.get("fit_chain"), "grain_start": record.get("grain_start"),
            "D_direct_nepers": D, "D_fitter_nepers": D_fit, "D_direct_vs_fitter_rel": D_rel,
            "held_exponent": ({k: v for k, v in hinfo.items() if k not in ("h_B_physical", "file_meta")}
                              if hinfo else None),
            "s_definition": "s = Y a + c T" + (" + held (the held exponent term)" if held is not None else "")}
    co2 = current_co2_floor()
    co2_arrays = {}
    if co2:                     # 2026-11-08: the CO2 flux floor (keys added only when it is on)
        meta.update(co2_floor=True, co2_floor_ratio=floor,
                    co2_floor_definition=CO2_FLOOR_TEXT.format(floor=repr(current_ratio())))
        co2_arrays["F_min_CO2_over_qmax_CO2"] = floor * np.asarray(el["ACO2"], float)
    path = out / "activity_exp_facets.npz"
    with path.open("xb") as f:
        np.savez_compressed(f, s=s, A=A, s0=np.asarray(P.s0, float),
                            held=np.asarray(held, float) if held is not None else np.zeros(len(s)),
                            F_min_over_qmax=floor * np.asarray(el["AH2O"], float),
                            seen=np.asarray(seen, bool), meta_json=np.array(json.dumps(meta, sort_keys=True)),
                            **co2_arrays)
    bundle = [src["bundle"]["path"] for src in getattr(P, "template_sources", [])]
    chain = record.get("fit_chain") or {}
    if chain.get("start_restart"):
        start = "the restart %s (sha256 %s; lambda %g, final cost %.6e); grain A %s, B %s" % (
            chain["start_restart"]["path"], chain["start_restart"]["sha256"], chain["start_restart"]["lambda"],
            chain["start_restart"]["cost"], (record.get("grain_start") or {}).get("source_A"),
            (record.get("grain_start") or {}).get("source_B"))
    elif chain.get("resumed_from"):
        start = "resumed from %s (wall-time segments)%s" % (chain["resumed_from"][-1]["path"],
                                                            "; first segment started from %s" % chain["start_restart"]["path"]
                                                            if chain.get("start_restart") else "")
    else:
        start = "the projection s0 (grain A %s, B %s)" % ((record.get("grain_start") or {}).get("source_A"),
                                                           (record.get("grain_start") or {}).get("source_B"))
    if chain.get("start_set"):                         # --start-set (2026-11-08), first segment of this fit
        start += "; then --start-set %s (the start restart's values replaced: %s; %s)" % (
            ", ".join("%s=%r" % (e["name"], e["value"]) for e in chain["start_set"]),
            ", ".join("%s %r" % (e["name"], e["start_restart_value"]) for e in chain["start_set"]),
            "; ".join(sorted({e["decision"] for e in chain["start_set"]})))
    if hinfo:
        held_text = ("h from %s (key %s, sha256 %s; file array sha256 %s): s = Y a + c T + h, h held fixed (user "
                     "2026-11-07: the accumulated fine-tuning change). rms_area(h) = %.6f nepers (= D at the prior "
                     "x0). D contribution: h = B h_B + h_perp (area-weighted projection on the fitted columns), "
                     "D^2 = |B (x - x0) + h|^2_area = |Lt (x - x0 + h_B)|^2 + |h_perp|^2_area, |h_perp|_area = %.6f "
                     "nepers (constant; D^2 contribution %.6e). D direct %.9f, fitter %.9f (rel %.1e)."
                     % (hinfo["path"], hinfo["key"], hinfo["sha256"], hinfo["array_sha256"], hinfo["rms_area_nepers"],
                        hinfo["D_perp_nepers"], hinfo["D_perp_squared"], D, D_fit, D_rel))
    else:
        held_text = "none (s = Y a + c T). D direct %.9f, fitter %.9f (rel %.1e)." % (D, D_fit, D_rel)
    text = PARAMETERIZATION_TEXT.format(utc=utc_now(), pid=PARAMETERIZATION_ID, start=start, held=held_text,
                                        ratio_decision=RATIO_DECISIONS.get(current_ratio()),
                                        templates=", ".join(P.template_names) or "(none)",
                                        bundle=", ".join(bundle) or "-", night=night, floor=floor, T0=USER_T0,
                                        late=LATE_START_UTC, lam=float(r["lam"]),
                                        ngrid=len(P.late_layout.grid_rows), reason=GRID_REASON,
                                        co2_text=("CO2 flux floor F_CO2 >= %s q_max,CO2 A_CO2,run, see below" % floor
                                                  if co2 else "CO2 has no floor"),
                                        co2_par=("\n".join(textwrap.wrap("CO2 flux floor (COPS_POS_FLOOR_CO2=1): "
                                                                          + CO2_FLOOR_TEXT.format(floor=floor)
                                                                          + ". Pressures p = p_H2O + kappa p_CO2(records "
                                                                          "not at the CO2 floor) + p_CO2(records at "
                                                                          "the floor, F_min,CO2 in place of kappa Q_CO2 "
                                                                          "A_CO2,run); d p / d kappa from the records "
                                                                          "not at the floor only.", 118)) + "\n"
                                                 if co2 else ""),
                                        co2_ratio=(" = CO2 floor" if co2 else ""),
                                        co2_npz=("                                F_min_CO2_over_qmax_CO2 = %s A_CO2,run "
                                                 "(ACO2; the CO2 floor over q_max,CO2)\n" % floor if co2 else ""))
    tpath = out / "ACTIVITY_PARAMETERIZATION.txt"
    with tpath.open("x") as f:
        f.write(text)
    record["products"] = {p.name: {"sha256": sha256(p), "size_bytes": p.stat().st_size} for p in (path, tpath)}
    if held is not None:
        record["held_activity_calls"] = int(HELD["calls"])
        if HELD["calls"] == 0:
            print("positive_fit: WARNING the held term was never added to act_lin (0 calls of the patched "
                  "activity_linear): the GPU pressures did not see h", flush=True)
    record["result_summary"] = {"D_direct_nepers": D, "D_fitter": D_fit, "D_direct_vs_fitter_rel": D_rel,
                                "D_check_tolerance": D_CHECK_RTOL, "D_check_pass": D_rel <= D_CHECK_RTOL,
                                "held_exponent": bool(hinfo),
                                "total_activity_over_A_orig": float(np.sum(el["area"] * A) / np.sum(el["area"] * el["AH2O"])),
                                "A_range": [float(A.min()), float(A.max())], "kappa": float(kappa)}
    print("positive_fit: wrote %s, %s (D direct %.9f nepers, fitter %.9f, rel %.1e%s)" % (
        path, tpath, D, D_fit, D_rel, "" if D_rel <= D_CHECK_RTOL else " > %.0e: WARNING D check failed" % D_CHECK_RTOL),
        flush=True)


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
        held = getattr(problem, "held", None)
        if held is not None:
            state.arrays["activity_linear"] = np.asarray(state.arrays["activity_linear"], float) + held
        state.meta["model_extension"] = {
            "activity_linear_includes_held": held is not None,
            "parameterization": PARAMETERIZATION_ID, "policy": current_policy(),
            "night_ratio": current_ratio(), "floor_ratio": float(os.environ.get(ENV_FLOOR, "nan")),
            "co2_floor": current_co2_floor(),
            "ratio_decision": RATIO_DECISIONS.get(current_ratio()),
            "held_exponent": record.get("held_exponent_marker"),
            "fit_chain": record.get("fit_chain"),
            "activity": "A = exp(activity_linear); activity_linear = Y a + c T (exponent)",
            "warning": "not a linear activity: tools that use max(activity_linear, 0) do not apply",
            "ascale": float(problem.ascale), "positive_activity": record.get("positive_activity"),
            "start_projection_weighting": PROJECTION_WEIGHTING,
            "late_window": {k: record.get("late_window", {}).get(k) for k in ("late_start_utc", "late_observations")}}
        if record.get("start_set"):                    # --start-set (2026-11-08); also in fit_chain (resumes keep it)
            state.meta["model_extension"]["start_set"] = record["start_set"]
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
    nreg = len(P.regularization(x)) if P.lam > 0 else 0
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
