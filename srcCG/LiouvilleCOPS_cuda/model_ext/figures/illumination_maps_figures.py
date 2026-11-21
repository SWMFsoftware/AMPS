#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/illumination_maps_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/illumination_maps_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c6df23b08582bfb9fcb2b0945987be7409d8eee3b5d00d6eef6d6fdd72c92f58
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: illumination maps (day / night / cast shadow, cos SZA, the H2O source-rate factor q/q_max) at the production-rate
#   map times 10:33:10 and 10:39:26 UTC in the layouts of the full sets' production-rate figures, the production rate next
#   to q/q_max, the numbers and the cross-case identity check (2026-11-18)
# RUN: jobs/illumination_numbers.sbatch (CPU), then jobs/illumination_maps.sbatch TIME=10:33:10 / 10:39:26 (GPU), then
#   'docs' on the login node (HOW-TO-RERUN.md of figures_cluster/12M/illumination_maps_103310_103926/)
# 2026-11-18: H2O night-side fix (user decision 2026-11-18, q_min/q_max 0.2 -> 0.02 as srcCG): QRATIO = the maps' nonneg_maps.H2O_QMIN_RATIO (set_qratio; numbers json key qratio_h2o), factor classes / ticks / texts derived from it (identical for 0.2)
# ==== end of header ====
"""Illumination maps at the production-rate map times (10:33:10 and 10:39:26 UTC on 2016-09-30) in the layouts of the
full sets' production-rate figures, and the H2O production rate next to the illumination factor of its source rate.

User request (2026-11-17 CST): "check the production rate plots -- i dont see the effect of the shadow or night side
decrease of the productions -- is it there? generate illumination maps for the time intervals that corresponds to that
on the production rate maps". Lead's specification (2026-11-18 CST): three quantities, (1) the illumination class
(day / night / cast shadow), (2) cos SZA of the day facets in 0.1 classes with night and shadow in their own colours,
(3) the illumination factor q/q_max of the H2O source rate (0.2 on night and cast-shadow facets, 0.2 + 0.8 cos SZA on
day facets; optionally also on a log scale), each for the WHOLE surface and for the facets seen during the fit, in
every layout of full_set_stage2_softplus_reg/production_rate/ (templates figure, its value panels, lon/lat and 3-D
ground tracks; relief shading and the ground-track overlay as there, Rosetta impact only); one side-by-side figure per
time (softplus reg: production rate next to the illumination factor, same windows); the numbers from the data; the
cross-case identity of cos SZA and illumination.

WHAT IS SHOWN (per facet, at the map time)
  illumination  the maps' field (nonneg_maps.illumination): 0 day (cos SZA > 0, Sun ray not blocked), 1 night
                (cos SZA <= 0), 2 cast shadow (cos SZA > 0, Sun ray blocked by the nucleus)
  cos SZA       the maps job's geometry (<geometry>/<label>.cos.f8; equal to the .vtu field cos_SZA, checked bitwise
                by pub_maps_figures.load_fields)
  factor        f = q_H2O / q_max = r on night and cast-shadow facets, r + (1 - r) min(cos SZA, 1) on day facets (the
                law of nonneg_maps.source_rates = fit_cops.Model.load_kernels / prepare), r = QRATIO = the maps' own
                nonneg_maps.H2O_QMIN_RATIO (2026-11-18: 0.02 since the H2O night-side fix; the figures of 2026-11-18
                used 0.2); checked: F_H2O = q A bitwise with q = q_min (night / shadow) else q_min + (q_max - q_min)
                min(cos, 1); q_min / q_max = r; |q / q_max - f| / f and |F / (q_max A) - f| / f small (reported)
  classes       fixed (the same in every panel): cos SZA 10 classes of 0.1 (day facets); f classes of 0.1 from
                floor(r / 0.1) 0.1 to 1 (day facets; r = 0.2: 8 classes from 0.2, r = 0.02: 10 classes from 0), night
                and cast shadow (f = r exactly) the colour of the triangle below them; log10 f classes of 0.05 from
                floor(log10 r / 0.05) 0.05 to 0 (r = 0.2: 14 from -0.70, r = 0.02: 34 from -1.70), night and cast
                shadow likewise (2026-11-18: derived from r, unchanged for r = 0.2)
  surface       whole: every facet coloured; seen: the fit's seen facets coloured (also those with A = 0: illumination
                exists there), the others grey with relief shading as in the production-rate maps
  relief        as in the production-rate figures: the templates figure and its value panels with the class colours
                modulated by the relief shading (pub_maps_v2_figures.relief_modulate, --template-relief 0.7); the ground
                tracks with flat colours (--gt-relief 0)
  side by side  per time, lon/lat global + impact area: (a, b) the H2O production rate F_H2O of the reference state
                exactly as in its production-rate templates figure (classes recomputed with pub_maps_v2_figures.Scale
                and checked equal to its manifest), (c, d) f on the same facets (seen, F > 0; seen F = 0 light green,
                unseen grey), (e, f) f on the whole surface

CODE REUSE: pub_maps_v2_figures.py and pub_maps_figures.py (figures_cluster/scripts) are imported unchanged, and
through them cops_cuda.figures (AMPS, unchanged). Copied with changes (marked "copy of"): the layouts of
pub_maps_v2_figures.main (templates figure, single value panels, lon/lat and 3-D ground tracks), with the illumination
colours and fixed colour bars instead of the per-panel log classes.

Subcommands (run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python):
  numbers   CPU job (jobs/illumination_numbers.sbatch): the reference state's fields at both times
            (pub_maps_figures.load_fields: bitwise against its maps' .vtu), the factor checks, the numbers per time for
            the whole and the seen surface, the class changes between the two times, and the cross-case check (cos_SZA
            and illumination of the .vtu of every --case against the reference, the geometry folders, the seen masks)
            -> <out>/checks/illumination_numbers.{json,txt}, <out>/checks/illumination_crosscase.{json,txt}
  figures   GPU job, one time (jobs/illumination_maps.sbatch; needs the numbers output for the captions)
            -> <out>/<quantity>/{12M_templates_*, 12M_groundtrack_*_{latlon,3d}, panels/*_{a_global,b_impact_area}},
               <out>/side_by_side/, <out>/checks/illumination_<hhmmss>_{class_tables.txt,manifest.json,captions.json}
  docs      login node: <out>/README.txt (numbers, the answer), CAPTIONS.md (collect_captions.py, unchanged),
            HOW-TO-RERUN.md (the jobs' commands.sh), code/ (+ SHA256SUMS)
Products of numbers / figures are staged in <out>/<step>.incomplete_<UTC> and moved into place file by file; existing
files are never overwritten.
"""

import argparse
import datetime as dt
import json
import math
import os
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import pub_maps_v2_figures as V2  # noqa: E402  (imports pub_maps_figures and cops_cuda.figures)

M = V2.M
FIG = V2.FIG

SCHEMA = "illumination_maps_v1"
T0 = time.time()
MAP_TIMES = (("10:33:10", "t11", "2016-09-30T10:33:10"), ("10:39:26", "t12", "2016-09-30T10:39:26"))
# 2026-11-18 (user decision, H2O night-side fix): q_min / q_max of H2O is no longer the literal 0.2 here. QRATIO is the
# value of the maps' own code, nonneg_maps.H2O_QMIN_RATIO (0.02 since 2026-11-18), set by set_qratio() as soon as that
# module is loaded (numbers, figures) or read back from the numbers json (docs), and checked against the data
# (factor_checks); the factor classes and the texts derive from it (for 0.2 they are those of 2026-11-18).
QRATIO = None
CLASS_KEYS = ("day", "night", "cast_shadow")
VARIANTS = (("whole", "whole_surface"), ("seen", "seen_only"))


def say(text):
    print("illumination_maps [%6.0f s]: %s" % (time.time() - T0, text), flush=True)


def set_qratio(value):
    """QRATIO = the maps' nonneg_maps.H2O_QMIN_RATIO (a module) or a number (docs: the numbers json); 2026-11-18."""
    global QRATIO
    if not isinstance(value, (int, float)):
        if not hasattr(value, "H2O_QMIN_RATIO"):
            raise SystemExit("the maps' nonneg_maps.py has no H2O_QMIN_RATIO (written before the 2026-11-18 H2O "
                             "night-side fix): make the maps with the current code first")
        value = value.H2O_QMIN_RATIO
    QRATIO = float(value)
    if not 0.0 < QRATIO < 1.0:
        raise SystemExit("H2O q_min / q_max %r outside (0, 1)" % QRATIO)
    return QRATIO


def qr():
    """Text pieces and class layout of the factor law for QRATIO (2026-11-18)."""
    r = QRATIO
    lo_f = math.floor(r / 0.1 + 1e-9) * 0.1                 # 0.2 -> 0.2, 0.02 -> 0
    lo_l = math.floor(math.log10(r) / 0.05 + 1e-9) * 0.05   # 0.2 -> -0.70, 0.02 -> -1.70
    return {"r": "%g" % r, "day": "%g + %g cos SZA" % (r, 1.0 - r), "dec": -math.log10(r), "fac": 1.0 / r,
            "lo_f": lo_f, "n_f": int(round((1.0 - lo_f) / 0.1)), "lo_l": lo_l, "n_l": int(round(-lo_l / 0.05)),
            "logr": math.log10(r), "logr_u": ("%.2f" % math.log10(r)).replace("-", "\u2212"),
            "lo_l_u": ("%.2f" % lo_l).replace("-", "\u2212")}


def hexrgb(h):
    return np.array([int(h[k:k + 2], 16) for k in (1, 3, 5)]) / 255.0


# Okabe-Ito colours (cops_cuda.figures.PALETTE): day orange, night blue, cast shadow reddish purple
C_DAY, C_NIGHT, C_SHADOW = hexrgb("#E69F00"), hexrgb("#0072B2"), hexrgb("#CC79A7")
QSYM = "$q/q_\\mathrm{max}$"


def utc_stamp():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def central(utc_text):
    """'2026-11-18T06:04:00Z' -> the same instant in US Central (CST/CDT)."""
    try:
        from zoneinfo import ZoneInfo
        t = dt.datetime.strptime(utc_text, "%Y%m%dT%H%M%SZ").replace(tzinfo=dt.timezone.utc)
        return t.astimezone(ZoneInfo("America/Chicago")).strftime("%Y-%m-%d %H:%M %Z")
    except Exception:  # noqa: BLE001
        return utc_text


# ============================================================================================ quantities
class Quantity:
    """One illumination quantity with fixed classes; per-facet codes >= 0: class index, < 0: special category."""

    def __init__(self, key, title, cb_label, cb_short, edges, colors, specials, labels_wide, labels_narrow,
                 under=None, categorical=False, ticklabels=None):
        self.key, self.title, self.cb_label, self.cb_short = key, title, cb_label, cb_short
        self.edges = np.asarray(edges, float)
        self.colors = np.asarray(colors, float)[:, :3]
        self.n = len(self.colors)
        self.step = float(self.edges[1] - self.edges[0])
        self.specials = specials             # [(code, legend label, rgb)]
        self.under = under                   # colour of the colour bar's 'min' triangle (the factor's night / shadow)
        self.categorical = categorical
        self.labels_wide, self.labels_narrow = labels_wide, labels_narrow
        self.ticklabels = ticklabels
        ns = len(specials)
        self.palette = np.zeros((ns + self.n, 3))
        for code, _, c in specials:
            self.palette[-code - 1] = c
        self.palette[ns:] = self.colors

    def rgb(self, codes):
        ns = len(self.specials)
        return self.palette[np.where(codes >= 0, codes + ns, -codes - 1)]

    def codes(self, ill, cos, f):
        """(codes, value) per facet; value = the quantity the classes are made of."""
        day = ill == 0
        if self.key == "illumination_class":
            return ill.astype(np.int16), ill.astype(np.float64)
        if self.key == "cos_sza":
            v = cos
            codes = np.where(ill == 1, -1, np.where(ill == 2, -2, 0)).astype(np.int16)
        elif self.key == "illumination_factor":
            v = f
            codes = np.where(day, 0, -1).astype(np.int16)
        else:
            v = np.log10(f)
            codes = np.where(day, 0, -1).astype(np.int16)
        k = np.floor((v[day] - self.edges[0]) / self.step + 1e-9).astype(np.int64)
        codes[day] = np.clip(k, 0, self.n - 1)
        return codes, v

    def check(self, codes, v, ill, cos, sh, f):
        """Every day facet's value inside its class; the special categories exactly the rule's facets."""
        day = ill == 0
        out = {}
        if self.key == "illumination_class":
            out["codes_equal_illumination"] = bool(np.array_equal(codes, ill.astype(np.int16)))
        else:
            k = codes[day].astype(np.int64)
            vv = v[day]
            out["day_values_inside_their_class"] = bool(np.all((k >= 0) & (vv >= self.edges[k] - 1e-9) &
                                                               (vv <= self.edges[k + 1] + 1e-9)))
            if self.key == "cos_sza":
                out["night_code_equal_cos_le_0"] = bool(np.array_equal(codes == -1, cos <= 0.0))
                out["shadow_code_equal_cos_gt_0_and_shadow"] = bool(np.array_equal(codes == -2,
                                                                                   (cos > 0.0) & (sh != 0)))
            else:
                out["special_code_equal_night_or_shadow"] = bool(np.array_equal(codes == -1, ill != 0))
                out["factor_of_special_equal_qratio"] = bool(np.all(f[ill != 0] == QRATIO))   # was _0.2 (2026-11-18)
        out["passed"] = all(out.values())
        return out


def make_quantities():
    import matplotlib
    cos_cols = matplotlib.colormaps["YlOrBr"](np.linspace(0.92, 0.16, 10))      # grazing: dark brown, overhead: pale
    # the factor: the same 'sunlight' ramp, darkest brown = QRATIO (night / cast shadow, the colour bar's triangle) to pale
    # yellow = 1; no grey (the unseen facets are grey) and no green (seen facets with F = 0 are light green); 2026-11-18
    # test job 136360: 'cividis' had a grey middle and was replaced
    ylb = matplotlib.colormaps["YlOrBr"]
    k = qr()                                  # 2026-11-18: the factor classes from QRATIO (0.2: as before)
    fac = ylb(np.linspace(1.0, 0.12, k["n_f"] + 1))[:, :3]
    lfac = ylb(np.linspace(1.0, 0.12, k["n_l"] + 1))[:, :3]
    e_cos = np.round(np.arange(11) * 0.1, 6)
    e_f = np.round(k["lo_f"] + np.arange(k["n_f"] + 1) * 0.1, 6)
    e_lf = np.round(k["lo_l"] + np.arange(k["n_l"] + 1) * 0.05, 6)
    t_f = sorted(set([round(QRATIO, 6), 0.2, 0.4, 0.6, 0.8, 1.0]))
    s1, s2 = (0.1, 0.2) if k["lo_l"] >= -1.0 else (0.2, 0.4)
    t_l1 = list(np.round(np.arange(math.ceil(k["lo_l"] / s1 - 1e-9), 1) * s1, 6))
    t_l2 = list(np.round(np.arange(math.ceil(k["lo_l"] / s2 - 1e-9), 1) * s2, 6))
    return [
        Quantity("illumination_class", "illumination class", "illumination class", "illumination class",
                 [0, 1, 2, 3], [C_DAY, C_NIGHT, C_SHADOW], [], None, None, categorical=True,
                 ticklabels=["day", "night", "cast shadow"]),
        Quantity("cos_sza", "cos SZA (day facets)", "cos SZA of the day facets (classes of 0.1)",
                 "cos SZA (day facets)", e_cos, cos_cols,
                 [(-1, "night (cos SZA ≤ 0)", C_NIGHT), (-2, "cast shadow (cos SZA > 0, Sun ray blocked)", C_SHADOW)],
                 [0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0]),
        Quantity("illumination_factor", "illumination factor %s of the H₂O source rate" % QSYM,
                 "%s = %s (day); triangle: %s (night, cast shadow); range %s–1" % (QSYM, k["day"], k["r"], k["r"]),
                 "%s (triangle: %s, night / shadow)" % (QSYM, k["r"]), e_f, fac[1:],
                 [(-1, "night or cast shadow: %s = %s" % (QSYM, k["r"]), fac[0])],
                 t_f, t_f, under=fac[0]),
        Quantity("illumination_factor_log", "log$_{10}$ of the illumination factor %s" % QSYM,
                 "log$_{10}$ %s (classes of 0.05 decade); triangle: log$_{10}$ %s = %s (night, cast shadow)" % (
                     QSYM, k["r"], k["logr_u"]),
                 "log$_{10}$ %s (triangle: %s, night / shadow)" % (QSYM, k["logr_u"]), e_lf, lfac[1:],
                 [(-1, "night or cast shadow: log$_{10}$ %s = %s" % (QSYM, k["logr_u"]), lfac[0])],
                 t_l1, t_l2, under=lfac[0]),
    ]


def illumination_factor(cos, ill):
    return np.where(ill != 0, QRATIO, QRATIO + (1.0 - QRATIO) * np.minimum(cos, 1.0))


def factor_checks(cos, sh, ill, F, A, q_h2o, nm):
    """The factor f against the model's own q (bitwise F = q A) and the maps' illumination rule."""
    qmin, qmax = float(q_h2o[0]), float(q_h2o[1])
    night = ill != 0
    c = {"q_min": qmin, "q_max": qmax, "q_min_over_q_max": qmin / qmax}
    c["illumination_equal_rule"] = bool(np.array_equal(ill, nm.illumination(cos, sh)))
    c["night_flag_equal_source_rates_rule"] = bool(np.array_equal(night, (cos <= 0.0) | (sh != 0)))
    q = np.where(night, qmin, qmin + (qmax - qmin) * np.minimum(cos, 1.0))
    c["F_equal_q_times_A_bitwise"] = bool(np.array_equal(F, q * A))
    f = illumination_factor(cos, ill)
    c["max_rel_diff_q_over_qmax_vs_f"] = float(np.max(np.abs(q / qmax - f) / f))
    ok = A > 1e-280
    r = F[ok] / (qmax * A[ok])
    c["max_rel_diff_F_over_qmaxA_vs_f"] = float(np.max(np.abs(r - f[ok]) / f[ok]))
    c["factor_min"], c["factor_max"] = float(f.min()), float(f.max())
    c["factor_min_day"] = float(f[~night].min()) if (~night).any() else None
    c["passed"] = bool(c["illumination_equal_rule"] and c["night_flag_equal_source_rates_rule"] and
                       c["F_equal_q_times_A_bitwise"] and abs(c["q_min_over_q_max"] - QRATIO) <= 1e-15 and
                       c["max_rel_diff_q_over_qmax_vs_f"] <= 1e-14 and c["max_rel_diff_F_over_qmaxA_vs_f"] <= 1e-13)
    return c, f


# ============================================================================================ numbers
def wquantile(v, w, qs):
    o = np.argsort(v, kind="stable")
    cw = np.cumsum(w[o]) / w.sum()
    return [float(v[o][min(int(np.searchsorted(cw, q)), len(o) - 1)]) for q in qs]


def wstd(v, w):
    m = np.dot(v, w) / w.sum()
    return float(math.sqrt(max(np.dot((v - m) ** 2, w) / w.sum(), 0.0)))


def scope_numbers(m, area, ill, f, F, A, qmax):
    tot = float(area[m].sum())
    dark = ill != 0
    out = {"facets": int(m.sum()), "area_m2": tot, "classes": {}}
    for k, name in enumerate(CLASS_KEYS):
        s = m & (ill == k)
        out["classes"][name] = {"facets": int(s.sum()), "area_m2": float(area[s].sum()),
                                "area_fraction": float(area[s].sum() / tot)}
    out["night_or_shadow_area_fraction"] = float(area[m & dark].sum() / tot)
    out["mean_factor_area_weighted"] = float(np.dot(f[m], area[m]) / tot)
    d = m & ~dark
    out["mean_factor_day_area_weighted"] = float(np.dot(f[d], area[d]) / area[d].sum())
    out["mean_cos_sza_day_area_weighted"] = float(np.dot(np.minimum((f[d] - QRATIO) / (1 - QRATIO), 1.0), area[d]) /
                                                  area[d].sum())
    nz = m & (F > 0)
    if not nz.any():
        raise SystemExit("scope without F > 0 facets")
    anz = float(area[nz].sum())
    out["nonzero_F"] = {"facets": int(nz.sum()), "area_m2": anz, "area_fraction_of_scope": anz / tot,
                        "night_or_shadow_area_fraction": float(area[nz & dark].sum() / anz),
                        "class_area_fraction": {name: float(area[nz & (ill == k)].sum() / anz)
                                                for k, name in enumerate(CLASS_KEYS)},
                        "mean_factor_area_weighted": float(np.dot(f[nz], area[nz]) / anz)}
    P = float(np.dot(F[m], area[m]))
    out["production_s-1"] = P
    out["production_share"] = {name: float(np.dot(F[m & (ill == k)], area[m & (ill == k)]) / P)
                               for k, name in enumerate(CLASS_KEYS)}
    out["night_or_shadow_production_share"] = float(np.dot(F[m & dark], area[m & dark]) / P)
    out["production_weighted_factor"] = P / (qmax * float(np.dot(A[m], area[m])))
    w = area[nz]
    lgF, lgA, lgf = np.log10(F[nz]), np.log10(A[nz]), np.log10(f[nz])
    sp = {}
    for name, v in (("log10_F", lgF), ("log10_A", lgA), ("log10_f", lgf)):
        q02, q50, q98 = wquantile(v, w, (0.02, 0.5, 0.98))
        sp[name] = {"q02": q02, "q50": q50, "q98": q98, "range_02_98_decades": q98 - q02, "std_decades": wstd(v, w)}
    sp["note"] = ("area-weighted over the facets of the scope with F > 0; log10 F = log10 q_max + log10 f + log10 A "
                  "exactly, so the spread of log10 F is that of log10 A plus at most %.2f decade from f" % qr()["dec"])
    out["spread_nonzero_F"] = sp
    out["mean_A_nonzero_area_weighted_per_class"] = {
        name: (float(np.dot(A[nz & (ill == k)], area[nz & (ill == k)]) / area[nz & (ill == k)].sum())
               if (nz & (ill == k)).any() else None) for k, name in enumerate(CLASS_KEYS)}
    return out


def cmd_numbers(a):
    out = Path(a.out_dir)
    stamp = utc_stamp()
    finals = [out / "checks" / n for n in ("illumination_numbers.json", "illumination_numbers.txt",
                                           "illumination_crosscase.json", "illumination_crosscase.txt")]
    for p in finals:
        if p.exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % p)
    stage = out / ("numbers.incomplete_%s" % stamp)
    (stage / "checks").mkdir(parents=True)
    say("stage %s" % stage)
    M.TIMES = tuple((lab, utc) for _, lab, utc in MAP_TIMES)
    fields = M.load_fields(a)
    if not fields["checks"]["passed"]:
        raise SystemExit("field checks failed: %s" % fields["checks"])
    nm = M.load_module(a.maps_code, "nonneg_maps_frozen_ill")
    set_qratio(nm)                                   # 2026-11-18: q_min/q_max of the maps' code
    A, area, seen = fields["A"], fields["area"], fields["seen"]
    nf = len(A)
    res = {"schema": SCHEMA, "created_utc": stamp, "host": socket.gethostname(), "qratio_h2o": QRATIO,
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
           "script_sha256": M.sha256_file(__file__), "helper_sha256": {
               "pub_maps_v2_figures.py": M.sha256_file(V2.__file__), "pub_maps_figures.py": M.sha256_file(M.__file__),
               "figures.py": M.sha256_file(FIG.__file__), "nonneg_maps.py": M.sha256_file(a.maps_code)},
           "reference": {"key": a.ref_key, "model_name": a.model_name, "maps_dir": str(Path(a.maps_dir).resolve()),
                         "ft_dir": str(Path(a.ft_dir).resolve()), "fit_dir": str(Path(a.fit_dir).resolve()),
                         "base": str(Path(a.base).resolve()), "geometry": str(Path(a.geometry).resolve())},
           "field_checks": fields["checks"], "facets": nf, "total_area_m2": float(area.sum()),
           "seen_facets": int(seen.sum()), "seen_area_m2": float(area[seen].sum()), "times": {}}
    ills = {}
    for hhmmss, lab, utc in MAP_TIMES:
        T = fields["times"][lab]
        cos, sh = nm.geometry_time(a.geometry, lab, nf)
        ill, F = T["illumination"], T["F_H2O"]
        chk, f = factor_checks(cos, sh, ill, F, A, T["q"]["H2O"], nm)
        if not chk["passed"]:
            raise SystemExit("factor checks failed at %s: %s" % (utc, chk))
        qmax = chk["q_max"]
        r = {"utc": utc, "label": lab, "rh_AU": T["rh_AU"], "sun": T["sun"], "q_H2O": list(T["q"]["H2O"]),
             "p": T["p"], "Q_H2O_s-1": T["Q_H2O"], "factor_checks": chk,
             "whole": scope_numbers(np.ones(nf, bool), area, ill, f, F, A, qmax),
             "seen": scope_numbers(seen, area, ill, f, F, A, qmax),
             "unseen": scope_numbers(~seen, area, ill, f, F, A, qmax)}
        r["seen_production_share_of_Q_H2O"] = r["seen"]["production_s-1"] / r["whole"]["production_s-1"]
        # sample facets (the lead's examples and a few of each class among seen F > 0)
        samp = []
        rng = np.random.default_rng(20261118)
        for k in range(3):
            cand = np.nonzero(seen & (F > 0) & (ill == k))[0]
            if len(cand):
                samp.extend(int(e) for e in rng.choice(cand, min(3, len(cand)), replace=False))
        for e in (928927, 916467):
            if e < nf:
                samp.append(e)
        r["sample"] = [{"facet": e, "illumination": int(ill[e]), "cos_SZA": float(cos[e]), "shadow": int(sh[e]),
                        "seen": bool(seen[e]), "A": float(A[e]), "F_H2O": float(F[e]), "factor": float(f[e]),
                        "F_over_A": float(F[e] / A[e]) if A[e] > 0 else None, "factor_times_qmax": float(f[e] * qmax)}
                       for e in samp]
        res["times"][hhmmss] = r
        ills[hhmmss] = ill
        say("%s: whole day/night/shadow area %s; seen F>0 night+shadow area fraction %.4f, production share %.4f" % (
            utc, ["%.4e" % r["whole"]["classes"][n]["area_m2"] for n in CLASS_KEYS],
            r["seen"]["nonzero_F"]["night_or_shadow_area_fraction"], r["seen"]["night_or_shadow_production_share"]))
    # ---- class changes between the two times
    i1, i2 = ills[MAP_TIMES[0][0]], ills[MAP_TIMES[1][0]]
    ch = {}
    for scope, m in (("whole", np.ones(nf, bool)), ("seen", seen)):
        mat = [[float(area[m & (i1 == j) & (i2 == k)].sum()) for k in range(3)] for j in range(3)]
        changed = m & (i1 != i2)
        ch[scope] = {"area_m2_from_row_to_col": mat, "classes": list(CLASS_KEYS),
                     "changed_facets": int(changed.sum()), "changed_area_m2": float(area[changed].sum()),
                     "changed_area_fraction": float(area[changed].sum() / area[m].sum())}
    s1 = np.asarray(res["times"][MAP_TIMES[0][0]]["sun"])
    s2 = np.asarray(res["times"][MAP_TIMES[1][0]]["sun"])
    ang = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(s1, s2) / np.linalg.norm(s1) / np.linalg.norm(s2))))))
    res["between_times"] = {"illumination_changes": ch, "sun_direction_angle_deg": ang}
    say("between the times: Sun direction moved %.3f deg; whole-surface area changing class %.4e m^2 (%.2f %%)" % (
        ang, ch["whole"]["changed_area_m2"], 100 * ch["whole"]["changed_area_fraction"]))
    del ills, i1, i2
    # ---- cross-case identity
    cross = crosscase(a, nm, fields)
    json.dump(res, open(stage / "checks" / "illumination_numbers.json", "w"), indent=1, default=float)
    (stage / "checks" / "illumination_numbers.txt").write_text(numbers_text(res), encoding="utf-8")
    json.dump(cross, open(stage / "checks" / "illumination_crosscase.json", "w"), indent=1, default=float)
    (stage / "checks" / "illumination_crosscase.txt").write_text(crosscase_text(cross), encoding="utf-8")
    move_into_place(stage, out)
    say("done")
    return 0


def crosscase(a, nm, fields):
    """cos_SZA and illumination of every case's maps (.vtu) at both times against the reference case's; the geometry
    folders; the seen masks of the static .vtu; field data TimeValue / rh_m."""
    cases = []
    for c in a.case:
        key, _, d = c.partition("=")
        cases.append((key, Path(d)))
    keys = [k for k, _ in cases]
    if a.ref_key not in keys:
        raise SystemExit("--ref-key %s not among --case keys %s" % (a.ref_key, keys))
    cases.sort(key=lambda kd: kd[0] != a.ref_key)          # reference first
    out = {"reference": a.ref_key, "cases": {}, "times": [t for t, _, _ in MAP_TIMES], "geometry_dirs": {}}
    ref = {}
    nf = len(fields["A"])
    geoms = {}
    for key, d in cases:
        man = json.load(open(d / "maps_manifest.json"))
        c = {"maps_dir": str(d.resolve()), "model": man["model"], "geometry": man["geometry"], "base": man["base"],
             "times": {}}
        geoms.setdefault(man["geometry"], []).append(key)
        for hhmmss, lab, utc in MAP_TIMES:
            stamp = utc.replace("-", "").replace(":", "")
            vtu = d / ("%s_flux_%s.vtu" % (man["model"], stamp))
            say("cross-case %s %s: read %s" % (key, hhmmss, vtu.name))
            rd, fd, _ = nm.read_vtu(vtu)
            cos = np.array(rd["cos_SZA"], copy=True)
            ill = np.array(rd["illumination"], copy=True)
            del rd
            cg, sg = nm.geometry_time(man["geometry"], lab, nf)
            t = {"vtu": str(vtu), "TimeValue": fd.get("TimeValue", [None])[0], "rh_m": fd.get("rh_m", [None])[0],
                 "cos_SZA_equal_own_geometry": bool(np.array_equal(cos, cg)),
                 "illumination_equal_rule_of_own_geometry": bool(np.array_equal(ill, nm.illumination(cg, sg))),
                 "cos_SZA_sha256": sha256_bytes(cos.tobytes()), "illumination_sha256": sha256_bytes(ill.tobytes())}
            del cg, sg
            if key == a.ref_key:
                ref[hhmmss] = (cos, ill, t["TimeValue"], t["rh_m"])
            else:
                rc, ri, rt, rr = ref[hhmmss]
                t["cos_SZA_equal_reference"] = bool(np.array_equal(cos, rc))
                t["illumination_equal_reference"] = bool(np.array_equal(ill, ri))
                t["cos_SZA_facets_differing"] = int(np.sum(cos != rc))
                t["illumination_facets_differing"] = int(np.sum(ill != ri))
                t["TimeValue_equal_reference"] = t["TimeValue"] == rt
                t["rh_m_equal_reference"] = t["rh_m"] == rr
            c["times"][hhmmss] = t
        st = d / ("%s_static.vtu" % man["model"])
        say("cross-case %s: seen mask of %s" % (key, st.name))
        rd, _, _ = nm.read_vtu(st)
        s = np.asarray(rd["seen"]).astype(bool)
        del rd
        c["seen_facets"] = int(s.sum())
        c["seen_area_m2"] = float(fields["area"][s].sum())
        c["seen_equal_reference"] = bool(np.array_equal(s, fields["seen"]))
        c["seen_facets_differing_from_reference"] = int(np.sum(s != fields["seen"]))
        out["cases"][key] = c
    rg = fields["maps_manifest"]["geometry"]
    ref_sha = {}
    for _, lab, _ in MAP_TIMES:
        for kind in ("cos.f8", "shadow.u1"):
            ref_sha[(lab, kind)] = M.sha256_file(Path(rg) / ("%s.%s" % (lab, kind)))
    out["reference_geometry"] = rg
    for g, ks in geoms.items():
        e = {"cases": ks}
        for hhmmss, lab, _ in MAP_TIMES:
            for kind in ("cos.f8", "shadow.u1"):
                p = Path(g) / ("%s.%s" % (lab, kind))
                e["%s_%s_sha256" % (lab, kind)] = (ref_sha[(lab, kind)] if Path(g).resolve() == Path(rg).resolve()
                                                   else M.sha256_file(p))
                e["%s_%s_equal_reference_geometry" % (lab, kind)] = e["%s_%s_sha256" % (lab, kind)] == ref_sha[(lab, kind)]
        out["geometry_dirs"][g] = e
    others = [c for k, c in out["cases"].items() if k != a.ref_key]
    out["all_cos_SZA_and_illumination_identical"] = bool(all(
        t["cos_SZA_equal_reference"] and t["illumination_equal_reference"] for c in others for t in c["times"].values()))
    out["all_own_geometry_consistent"] = bool(all(
        t["cos_SZA_equal_own_geometry"] and t["illumination_equal_rule_of_own_geometry"]
        for c in out["cases"].values() for t in c["times"].values()))
    out["all_geometry_files_identical"] = bool(all(v for e in out["geometry_dirs"].values() for k, v in e.items()
                                                   if k.endswith("_equal_reference_geometry")))
    out["all_seen_identical"] = bool(all(c["seen_equal_reference"] for c in out["cases"].values()))
    return out


def sha256_bytes(b):
    import hashlib
    return hashlib.sha256(b).hexdigest()


def numbers_text(r):
    L = ["# illumination_maps_figures.py numbers (%s), %s UTC on %s, Slurm job %s" % (
        r["schema"], r["created_utc"], r["host"], r["slurm_job_id"]),
        "# reference state: %s (%s)" % (r["reference"]["key"], r["reference"]["model_name"]),
        "# maps %s; geometry %s" % (r["reference"]["maps_dir"], r["reference"]["geometry"]),
        "# field checks against the maps' .vtu passed: %s" % r["field_checks"]["passed"], ""]
    for hhmmss, t in r["times"].items():
        L.append("== %s UTC (r_h %.6f AU, q_max %.6g, q_min/q_max %.17g); factor checks passed %s "
                 "(F = q A bitwise %s; max rel |q/q_max - f| %.2e; max rel |F/(q_max A) - f| %.2e)" % (
                     t["utc"], t["rh_AU"], t["q_H2O"][1], t["factor_checks"]["q_min_over_q_max"],
                     t["factor_checks"]["passed"], t["factor_checks"]["F_equal_q_times_A_bitwise"],
                     t["factor_checks"]["max_rel_diff_q_over_qmax_vs_f"],
                     t["factor_checks"]["max_rel_diff_F_over_qmaxA_vs_f"]))
        for scope in ("whole", "seen", "unseen"):
            s = t[scope]
            L.append("  [%s surface] %d facets, %.6e m^2" % (scope, s["facets"], s["area_m2"]))
            for n in CLASS_KEYS:
                c = s["classes"][n]
                L.append("    %-12s %9d facets  %.6e m^2  %7.3f %% of the area" % (
                    n, c["facets"], c["area_m2"], 100 * c["area_fraction"]))
            L.append("    area-weighted mean q/q_max %.4f (day facets only %.4f; mean cos SZA of the day facets %.4f)" % (
                s["mean_factor_area_weighted"], s["mean_factor_day_area_weighted"], s["mean_cos_sza_day_area_weighted"]))
            nz = s["nonzero_F"]
            L.append("    F_H2O > 0 on %d facets, %.6e m^2 (%.3f %% of the area); of that area night or cast shadow "
                     "%.3f %% (day %.3f %%, night %.3f %%, shadow %.3f %%); their mean q/q_max %.4f" % (
                         nz["facets"], nz["area_m2"], 100 * nz["area_fraction_of_scope"],
                         100 * nz["night_or_shadow_area_fraction"], 100 * nz["class_area_fraction"]["day"],
                         100 * nz["class_area_fraction"]["night"], 100 * nz["class_area_fraction"]["cast_shadow"],
                         nz["mean_factor_area_weighted"]))
            L.append("    production %.6e s^-1: from day %.3f %%, night %.3f %%, cast shadow %.3f %% (night + shadow "
                     "%.3f %%); production-weighted q/q_max (production / (q_max sum A area)) %.4f" % (
                         s["production_s-1"], 100 * s["production_share"]["day"], 100 * s["production_share"]["night"],
                         100 * s["production_share"]["cast_shadow"], 100 * s["night_or_shadow_production_share"],
                         s["production_weighted_factor"]))
            sp = s["spread_nonzero_F"]
            L.append("    spread over F > 0 (area-weighted 2-98 %% range / std, decades): log10 F %.2f / %.2f, log10 A "
                     "%.2f / %.2f, log10 f %.2f / %.2f" % (
                         sp["log10_F"]["range_02_98_decades"], sp["log10_F"]["std_decades"],
                         sp["log10_A"]["range_02_98_decades"], sp["log10_A"]["std_decades"],
                         sp["log10_f"]["range_02_98_decades"], sp["log10_f"]["std_decades"]))
            L.append("    area-weighted mean A on F > 0 facets: %s" % ", ".join(
                "%s %s" % (n, "-" if v is None else "%.4e" % v)
                for n, v in s["mean_A_nonzero_area_weighted_per_class"].items()))
        L.append("  seen production / Q_H2O = %.4f" % t["seen_production_share_of_Q_H2O"])
        L.append("  sample facets: " + "; ".join("%d ill %d cos %.4f seen %s F/A %s f*q_max %.6g" % (
            e["facet"], e["illumination"], e["cos_SZA"], e["seen"],
            "-" if e["F_over_A"] is None else "%.6g" % e["F_over_A"], e["factor_times_qmax"]) for e in t["sample"]))
        L.append("")
    b = r["between_times"]
    L.append("== between the two times: Sun direction moved %.4f deg" % b["sun_direction_angle_deg"])
    for scope, c in b["illumination_changes"].items():
        L.append("  [%s] facets changing class %d, %.6e m^2 (%.3f %%); area from (row) to (col) day/night/shadow: %s" % (
            scope, c["changed_facets"], c["changed_area_m2"], 100 * c["changed_area_fraction"],
            " | ".join(" ".join("%.4e" % v for v in row) for row in c["area_m2_from_row_to_col"])))
    return "\n".join(L) + "\n"


def crosscase_text(c):
    L = ["# cross-case identity of cos_SZA and illumination (maps' .vtu), reference %s" % c["reference"],
         "all cos_SZA and illumination identical to the reference (bitwise, every facet, both times): %s" %
         c["all_cos_SZA_and_illumination_identical"],
         "every case's .vtu equal to its own geometry folder (cos bitwise, illumination = rule): %s" %
         c["all_own_geometry_consistent"],
         "all geometry folders' t11/t12 cos.f8 and shadow.u1 files identical (sha256): %s" %
         c["all_geometry_files_identical"],
         "seen masks (static .vtu) identical to the reference: %s" % c["all_seen_identical"], ""]
    for k, v in c["cases"].items():
        L.append("[%s] model %s, maps %s, geometry %s; seen %d facets (%.6e m^2), equal reference %s (%d differ)" % (
            k, v["model"], v["maps_dir"], v["geometry"], v["seen_facets"], v["seen_area_m2"], v["seen_equal_reference"],
            v["seen_facets_differing_from_reference"]))
        for t, x in v["times"].items():
            L.append("    %s: %s; cos sha256 %s..., ill sha256 %s...; own geometry %s / %s%s" % (
                t, Path(x["vtu"]).name, x["cos_SZA_sha256"][:16], x["illumination_sha256"][:16],
                x["cos_SZA_equal_own_geometry"], x["illumination_equal_rule_of_own_geometry"],
                "" if "cos_SZA_equal_reference" not in x else
                "; = reference: cos %s (%d differ), illumination %s (%d differ), TimeValue %s, rh_m %s" % (
                    x["cos_SZA_equal_reference"], x["cos_SZA_facets_differing"], x["illumination_equal_reference"],
                    x["illumination_facets_differing"], x["TimeValue_equal_reference"], x["rh_m_equal_reference"])))
    L.append("")
    for g, e in c["geometry_dirs"].items():
        L.append("geometry %s (cases %s): %s" % (g, ", ".join(e["cases"]), ", ".join(
            "%s %s" % (k.replace("_equal_reference_geometry", ""), v) for k, v in e.items()
            if k.endswith("_equal_reference_geometry"))))
    return "\n".join(L) + "\n"


def move_into_place(stage, out):
    for p in sorted(stage.rglob("*")):
        if p.is_file():
            dest = out / p.relative_to(stage)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                raise SystemExit("%s exists" % dest)
            p.rename(dest)
    for d in sorted((p for p in stage.rglob("*") if p.is_dir()), key=lambda p: -len(p.parts)):
        d.rmdir()
    stage.rmdir()


# ============================================================================================ drawing helpers
def fixed_colorbar(fig, cax, Q, label, narrow=False, orientation="horizontal", fontsize=6.8):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cm = ListedColormap(Q.colors)
    ext = "neither"
    if Q.under is not None:
        cm.set_under(Q.under)
        ext = "min"
    norm = BoundaryNorm(Q.edges, cm.N)
    if Q.categorical:
        ticks = 0.5 * (Q.edges[:-1] + Q.edges[1:])
        labels = Q.ticklabels
    else:
        ticks = Q.edges
        want = Q.labels_narrow if narrow else Q.labels_wide
        labels = [("%g" % (e + 0.0)).replace("-", "−") if any(abs(e - w) < 1e-9 for w in want) else "" for e in Q.edges]
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation=orientation, ticks=ticks,
                      spacing="uniform", extend=ext)
    (cb.ax.set_xticklabels if orientation == "horizontal" else cb.ax.set_yticklabels)(labels)
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=6.3, length=0 if Q.categorical else 2.5)
    cb.set_label(label, fontsize=fontsize)
    return cb


def surface_handles_q(Q, variant):
    from matplotlib.patches import Patch
    hs = [Patch(facecolor=c, edgecolor="none", label=lab) for _, lab, c in Q.specials]
    if variant == "seen":
        hs.append(Patch(facecolor=V2.NOT_SEEN * 0.85, edgecolor="none", label="not seen during the fit"))
    return hs


def surface_rgb(Q, codes, seen, intensity, variant, relief=0.0):
    """Per-facet displayed colour: the class colours (relief-modulated with relief > 0); seen variant: unseen grey,
    relief-shaded (as the production-rate maps)."""
    rgb = Q.rgb(codes)
    col = np.ones(len(codes), bool) if variant == "whole" else seen
    if relief > 0.0:
        rgb[col] = V2.relief_modulate(rgb[col], intensity[col], relief)
    if variant == "seen":
        rgb[~seen] = V2.NOT_SEEN[None, :] * intensity[~seen, None]
    return np.clip(rgb, 0.0, 1.0)


def same_mask_rgb(Q, codes, seen, F, intensity, relief):
    """The factor on the production-rate maps' facets: seen F > 0 class colours (relief-modulated), seen F = 0 light green
    and unseen grey, both relief-shaded (pub_maps_v2_figures.facet_rgb's surface colours)."""
    rgb = V2.NOT_SEEN[None, :] * intensity[:, None]
    z = seen & (F == 0)
    rgb[z] = V2.ZERO_SEEN[None, :] * intensity[z, None]
    c = seen & (F > 0)
    rgb[c] = Q.rgb(codes[c])
    if relief > 0.0:
        rgb[c] = V2.relief_modulate(rgb[c], intensity[c], relief)
    return np.clip(rgb, 0.0, 1.0)


def panel_table(label, Q, codes, seen, area, disp, variant, F=None):
    """Facets drawn in one panel per class / category (each facet once), with the sum check."""
    if variant == "whole":
        d, u, z = disp, disp[:0], disp[:0]
    elif variant == "seen":
        d, u, z = disp[seen[disp]], disp[~seen[disp]], disp[:0]
    else:                                         # same mask as the production-rate maps
        d = disp[seen[disp] & (F[disp] > 0)]
        z = disp[seen[disp] & (F[disp] == 0)]
        u = disp[~seen[disp]]
    rows = []
    cats = [(code, lab) for code, lab, _ in Q.specials] + [(k, None) for k in range(Q.n)]
    for code, lab in cats:
        m = codes[d] == code
        if code >= 0:
            lab = (Q.ticklabels[code] if Q.categorical else "[%g, %g%s" % (
                Q.edges[code], Q.edges[code + 1], "]" if code == Q.n - 1 else ")"))
        rows.append({"code": int(code), "label": lab.replace("$", "") if lab else lab, "facets": int(m.sum()),
                     "area_m2": float(area[d][m].sum())})
    tot = float(area[d].sum())
    out = {"panel": label, "variant": variant, "displayed_facets": int(disp.size), "coloured_facets": int(d.size),
           "coloured_area_m2": tot, "unseen_grey_facets": int(u.size), "seen_zero_F_facets": int(z.size),
           "rows": rows, "sum_check": int(sum(r["facets"] for r in rows)) == int(d.size) and
           int(d.size + u.size + z.size) == int(disp.size)}
    L = ["[%s] %s: displayed facets %d, coloured %d (%.6e m^2), unseen grey %d, seen F = 0 light green %d; counts add "
         "up: %s" % (label, variant, disp.size, d.size, tot, u.size, z.size, out["sum_check"])]
    for r in rows:
        L.append("    %-58s %10d facets %14.6e m^2 %7.3f %%" % (r["label"], r["facets"], r["area_m2"],
                                                               100 * r["area_m2"] / tot if tot > 0 else 0.0))
    if not out["sum_check"]:
        raise SystemExit("panel table check failed: %s" % label)
    return out, "\n".join(L)


def fmt_area(x):
    return "%.3e m²" % x


# ============================================================================================ figures
def cmd_figures(a):
    out = Path(a.out_dir)
    hit = [t for t in MAP_TIMES if t[0] == a.time]
    if len(hit) != 1:
        raise SystemExit("--time must be one of %s" % [t[0] for t in MAP_TIMES])
    hhmmss, lab, utc = hit[0]
    tag = hhmmss.replace(":", "")
    V2.NMIN, V2.NMAX = (int(v) for v in a.classes.split(","))
    M.MODEL_NAME = a.model_name
    only = set(a.only.split(",")) if a.only else None
    ndir = Path(a.numbers_dir) if a.numbers_dir else out / "checks"
    numbers = json.load(open(ndir / "illumination_numbers.json"))
    cross = json.load(open(ndir / "illumination_crosscase.json"))
    nt = numbers["times"][hhmmss]
    stamp = utc_stamp()
    final = out / "checks" / ("illumination_%s_manifest.json" % tag)
    if final.exists():
        raise SystemExit("%s exists (never overwritten; move it aside first)" % final)
    stage = out / ("figures_%s.incomplete_%s" % (tag, stamp))
    stage.mkdir(parents=True)
    say("stage %s, time %s" % (stage, utc))

    # ---- track (observation-only base) and fields (fit base), exactly as pub_maps_v2_figures.main
    base_tab = FIG.read_table(a.track_base)
    window = FIG.late_window(base_tab, a.bundle, None, None)
    gt = FIG.groundtrack_data(a.track_base, None, a.bundle, window, a.device)
    say("track: %d rows, impact %s UTC; window %d rows" % (len(gt["table"]["utc"]), gt["impact"]["utc"],
                                                           len(window["points"])))
    M.TIMES = ((lab, utc),)
    fields = M.load_fields(a)
    if not fields["checks"]["passed"]:
        raise SystemExit("field checks failed")
    mb = [(n[2:], x_, lo_) for n, x_, lo_ in M.multiplier_bounds(a.fit_dir)]
    lows = sorted({lo for _, _, lo in mb})
    at = [n for n, x, lo in mb if x - lo <= 1e-9 * max(1.0, abs(lo))]
    if mb and lows != [0.0]:
        M.MULT_NOTE = ("; the %d free multipliers bounded below by %s, %s at that bound" % (
            len(mb), " / ".join("%g" % v for v in lows),
            ("none" if not at else ", ".join(n.replace("slope_", "").replace("_", "-") + " deg" for n in at))))
    if not np.array_equal(fields["normal"], gt["normals"]):
        raise SystemExit("normals differ between the track mesh and the fields")
    nm = M.load_module(a.maps_code, "nonneg_maps_frozen_ill")
    set_qratio(nm)                                   # 2026-11-18: q_min/q_max of the maps' code
    if numbers.get("qratio_h2o") != QRATIO:
        raise SystemExit("numbers json q_min/q_max %r differs from the maps' code %r" % (numbers.get("qratio_h2o"), QRATIO))
    e1 = nm.read_elements(Path(a.base) / "elements.bin")
    e2 = nm.read_elements(Path(a.track_base) / "elements.bin")
    same = {k: bool(np.array_equal(e1[k], e2[k])) for k in ("centroid", "normal", "area", "ACO2")}
    del e1, e2
    if not all(same.values()):
        raise SystemExit("track base mesh differs: %s" % same)
    fields["checks"]["track_base_elements_equal_fit_base"] = same
    seen, area, A = fields["seen"], fields["area"], fields["A"]
    T = fields["times"][lab]
    F_h2o = T["F_H2O"]
    nf = len(A)
    cos, sh = nm.geometry_time(a.geometry, lab, nf)
    ill = T["illumination"]
    fchk, f = factor_checks(cos, sh, ill, F_h2o, A, T["q"]["H2O"], nm)
    if not fchk["passed"]:
        raise SystemExit("factor checks failed: %s" % fchk)
    if fchk["F_equal_q_times_A_bitwise"] != nt["factor_checks"]["F_equal_q_times_A_bitwise"]:
        raise SystemExit("numbers output disagrees with this job's checks")
    say("fields, factor checks passed; f in [%.6g, %.6g]" % (fchk["factor_min"], fchk["factor_max"]))
    meta = FIG.bundle_metadata(a.bundle)
    slope, cache_path = FIG.gravity_slope(a.bundle, a.gravity_cache, len(gt["faces"]))
    ctx = {"gt": gt, "fields": fields, "meta": meta, "cache_path": cache_path, "args": a,
           "track_base": str(Path(a.track_base).resolve()), "fit_base": str(Path(a.base).resolve()),
           "track_offset": M.gauge_offset(a.track_base), "fit_offset": M.gauge_offset(a.base)}
    plt = FIG.style()
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl, za = M.zoom_box(gt, late)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), V2.GLOBAL_W), "zoom": (zl, za, V2.ZOOM_W)}
    Qs = make_quantities()
    qcodes = {}
    qchecks = {}
    for Q in Qs:
        codes, v = Q.codes(ill, cos, f)
        qchecks[Q.key] = Q.check(codes, v, ill, cos, sh, f)
        if not qchecks[Q.key]["passed"]:
            raise SystemExit("class check failed for %s: %s" % (Q.key, qchecks[Q.key]))
        qcodes[Q.key] = codes
        del v
    say("class checks passed: %s" % {k: v["passed"] for k, v in qchecks.items()})
    QF = [Q for Q in Qs if Q.key == "illumination_factor"][0]

    def wanted(qkey, variant):
        return only is None or ("%s:%s" % (qkey, variant)) in only or qkey in only

    # ---- face images: lon/lat windows and 3-D panels (time-independent), slope / cliff images
    from matplotlib.colors import Normalize
    snorm = Normalize(0.0, 90.0)
    magma = plt.get_cmap("magma")
    slope_rgb = magma(snorm(np.clip(slope, 0, 90)))[:, :3]
    cliff_rgb = FIG.face_colors(len(F), None)
    masks = gt["masks"] or {}
    cliff_rgb[masks.get("cliff_ge_45deg", slope >= 45.0)] = FIG.CLIFF45_COLOR
    cliff_rgb[masks.get("cliff_ge_60deg", slope >= 60.0)] = FIG.CLIFF60_COLOR
    cliff_rgb = cliff_rgb * I_ll[:, None]
    faces_ll, disp_ll, imgs = {}, {}, {}
    for wname, (lon_r, lat_r, width) in wins.items():
        face, ss = V2.lonlat_faces(V, F, lon_r, lat_r, width, a.device)
        faces_ll[wname] = (face, ss)
        disp_ll[wname] = V2.displayed(face)
        imgs[("slope", wname)] = V2.shade(face, slope_rgb, ss)
        imgs[("cliff", wname)] = V2.shade(face, cliff_rgb, ss)
        say("window %s: %d displayed facets" % (wname, disp_ll[wname].size))
    del slope_rgb, cliff_rgb
    xyz = np.column_stack((table["x"], table["y"], table["z"]))
    ip = np.asarray(impact["point_m"])
    ground = gt["ground"]
    view = FIG.OrthoView(xyz[late][0] - ip)
    I3, _front = M.ortho_intensity(N, view)
    body_half = 1.04 * np.max(np.abs(np.column_stack(view.project(V - V.mean(axis=0))[:2])))
    zoom_half = 1.15 * max(np.max(np.linalg.norm(ground[late] - ip, axis=1)), 1.0)
    for name in ("late_ng_not_rg", "late_rg_not_ng"):
        if name in masks and masks[name].any():
            cen = V[F[masks[name]]].mean(axis=1)
            zoom_half = max(zoom_half, 1.08 * np.max(np.abs(np.column_stack(view.project(cen - ip)[:2]))))
    close_half = 1.3 * max(np.max(np.abs(np.column_stack(view.project(ground[late] - ip)[:2]))), 20.0)
    panels3 = ((V.mean(axis=0), body_half), (ip, zoom_half), (ip, close_half))
    pixels = a.pixels3d
    dummy = np.tile(V2.NOT_SEEN, (len(F), 1))
    r3 = []
    for k, (center, half) in enumerate(panels3):
        _, to_image, visible, face, ss = FIG.render_ortho(V, F, N, view, center, half, pixels, dummy, a.device,
                                                          return_faces=True)
        r3.append({"to_image": to_image, "visible": visible, "face": face, "ss": ss, "half": half,
                   "disp": V2.displayed(face)})
        say("3-D panel %s: %d displayed facets" % ("abc"[k], r3[-1]["disp"].size))
    del dummy

    # ---- texts
    trk = M.track_text(ctx)
    slp = M.slope_text(ctx)
    mesh = ("on the %s shape model (%s facets), longitude/latitude, outermost surface where overhangs put several "
            "surface layers at one longitude/latitude" % (Path(gt["shape"]).name, M.thousands(len(F))))
    geom_dir = str(Path(a.geometry).resolve())
    ill_def = (
        "Illumination as the model computes it (the maps' fields cos_SZA and illumination, nonneg_maps.py; the "
        "fitter's law of ../LiouvilleCOPS/tools/fit_cops.py Model.load_kernels / prepare = src/Source.cpp): cos SZA = "
        "facet normal · Sun direction (Sun = spkpos_c(SUN, et, 67P/C-G_CK, LT+S, CHURYUMOV-GERASIMENKO)); day: cos SZA "
        "> 0 and the Sun ray from the facet is not blocked; night: cos SZA ≤ 0; cast shadow: cos SZA > 0 and the ray "
        "from the facet centroid (offset 1e-6 of the longest edge along the normal) toward the Sun hits the nucleus; "
        "for all facets (geometry %s, r_h = %.6f AU)." % (geom_dir, T["rh_AU"]))
    law = ("The H2O source rate of the model is F_H2O = q_H2O A with q_H2O = q_min on night and cast-shadow facets and "
           "q_min + (q_max - q_min) cos SZA on day facets, q_max = 600 / r_h^%.7g = %.4f, q_min = %s q_max (checked "
           "on the data: F_H2O = q_H2O A bitwise on every facet); the illumination factor q_H2O/q_max is therefore %s "
           "on night and cast-shadow facets and %s on day facets (range %s to 1): at a given activity "
           "A, illumination changes F_H2O by at most a factor %g (%.2f decade)." % (
               T["p"], T["q"]["H2O"][1], qr()["r"], qr()["r"], qr()["day"], qr()["r"], qr()["fac"], qr()["dec"]))
    cr = cross
    others = [k for k in cr["cases"] if k != cr["reference"]]
    if cr["all_cos_SZA_and_illumination_identical"]:
        fitind = ("Illumination does not depend on the fit: cos SZA and the illumination class are bitwise identical on "
                  "every facet at both map times in the maps (.vtu) of all %d full-set cases (%s; geometry folders %s; "
                  "checks/illumination_crosscase.txt), so the whole-surface maps serve every case." % (
                      len(cr["cases"]), ", ".join(cr["cases"]), ", ".join(Path(g).parent.name
                                                                         for g in cr["geometry_dirs"])))
    else:
        fitind = ("Cross-case check: cos SZA / illumination differ between the cases' maps (checks/"
                  "illumination_crosscase.txt); these maps are those of %s." % cr["reference"])
    seen_note = ("the seen mask is the same in all %d cases" % len(cr["cases"]) if cr["all_seen_identical"] else
                 "the seen mask of %s; the seen masks of %s differ from it (checks/illumination_crosscase.txt)" % (
                     cr["reference"], ", ".join(k for k in others if not cr["cases"][k]["seen_equal_reference"])))
    W_, S_ = nt["whole"], nt["seen"]

    def class_numbers(s):
        return ", ".join("%s %s (%.1f %%)" % (n.replace("_", " "), fmt_area(s["classes"][n]["area_m2"]),
                                             100 * s["classes"][n]["area_fraction"]) for n in CLASS_KEYS)
    numbers_txt = (
        "At %s UTC: whole surface (%s): %s; area-weighted mean q_H2O/q_max %.3f. Seen surface (%s): %s; mean "
        "q_H2O/q_max %.3f (checks/illumination_numbers.txt)." % (
            hhmmss, fmt_area(W_["area_m2"]), class_numbers(W_), W_["mean_factor_area_weighted"],
            fmt_area(S_["area_m2"]), class_numbers(S_), S_["mean_factor_area_weighted"]))
    seen_def = (
        "Colour only on the facets seen during the fit (the fit's nonneg_facets.npz 'seen': the facet contributes to "
        "at least one fitted NG or RG observation; %s of the %s facets, %s of the area; %s), as in the production-rate "
        "maps; facets not seen during the fit are grey with relief shading. Unlike the production-rate maps, the seen "
        "facets with F_H2O = 0 (A = 0, light green there) are coloured too: illumination exists on every facet." % (
            M.thousands(seen.sum()), M.thousands(len(seen)), M.pct(np.sum(area[seen]) / np.sum(area)), seen_note))
    whole_def = "Every facet coloured (the whole surface, not only the facets seen during the fit)."
    relief_tpl = ("colours modulated by the relief shading (strength %g: darkened on slopes facing away from the light, "
                  "lightened toward white on slopes facing it, unchanged on level ground lit from 40 deg above the "
                  "horizon; hue kept; the relief light is the figure's, not the Sun); the colour bar shows the unshaded "
                  "colours" % a.template_relief)
    relief_gt = ("drawn without relief shading as the coloured facets of the production-rate ground-track maps, so "
                 "that the colour equals the colour bar's" if a.gt_relief <= 0 else
                 "colours modulated by the relief shading (strength %g)" % a.gt_relief)
    pr_steps = {}
    pr_man = Path(a.pr_dir) / "checks" / ("h2o_%s_manifest.json" % tag) if a.pr_dir else None
    if pr_man and pr_man.exists():
        pm = json.load(open(pr_man))
        pr_steps = {k: v["scale"]["step"] for k, v in pm["class_tables"].items()}

    def qtext(Q):
        if Q.key == "illumination_class":
            return ("illumination class at %s UTC on 2016-09-30: day (orange, #E69F00), night (blue, #0072B2), cast "
                    "shadow (reddish purple, #CC79A7)" % hhmmss)
        if Q.key == "cos_sza":
            return ("cosine of the solar zenith angle, cos SZA, of the day facets at %s UTC on 2016-09-30 in 10 classes "
                    "of 0.1 from 0 to 1 (matplotlib 'YlOrBr': dark brown = grazing Sun, pale yellow = Sun near the "
                    "local zenith); night facets blue (#0072B2), cast-shadow facets reddish purple (#CC79A7; their cos "
                    "SZA > 0 is not shown)" % hhmmss)
        if Q.key == "illumination_factor":
            k = qr()
            return ("illumination factor of the H2O source rate, q_H2O/q_max, at %s UTC on 2016-09-30: %s "
                    "on day facets in %d classes of 0.1 from %g to 1 (matplotlib 'YlOrBr' reversed: dark brown = low, pale "
                    "yellow = high), %s exactly on night and "
                    "cast-shadow facets (the colour of the triangle at the low end of the colour bar); range %s to 1. "
                    "It is the factor by which illumination modulates the production rate: F_H2O = q_max (q_H2O/q_max) "
                    "A" % (hhmmss, k["day"], k["n_f"], k["lo_f"], k["r"], k["r"]))
        k = qr()
        return ("log10 of the illumination factor of the H2O source rate, log10(q_H2O/q_max), at %s UTC on 2016-09-30: "
                "log10(%s) on day facets in %d classes of 0.05 decade from %.2f to 0 (matplotlib "
                "'YlOrBr' reversed: dark brown = low, pale yellow = high), log10 %s = %.3f on night and cast-shadow facets (the colour of the triangle at the low "
                "end of the colour bar); the whole range is %.2f decade%s" % (
                    hhmmss, k["day"], k["n_l"], k["lo_l"], k["r"], k["logr"], k["dec"],
                    (", i.e. %.1f classes of the production-rate maps at this time (%s decade per class in their "
                             "lon/lat panels, %s decade in their 3-D panels)" % (
                                 k["dec"] / pr_steps["latlon_global"], "%g" % pr_steps["latlon_global"],
                                 "%g" % pr_steps["3d_a"])) if pr_steps else ""))

    written, captions, ttext, tables = [], [], [], {}

    def save(fig, prefix, caption):
        written.extend(FIG.save(fig, prefix))
        plt.close(fig)
        written.append(FIG.write_caption(prefix, caption))
        captions.append((str(Path(prefix).relative_to(stage)), caption))

    def draw_map(ax, img, wname, c, ylabel, annotate):
        """copy of pub_maps_v2_figures.main draw_map (image + black track overlay of the templates figure)."""
        lon_r, lat_r, _ = wins[wname]
        ax.imshow(img, extent=(lon_r[0], lon_r[1], lat_r[0], lat_r[1]), interpolation="lanczos",
                  aspect="equal" if c == 0 else 1.0 / coslat)
        V2.latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, c, zl, za, annotate=annotate)
        V2.map_axes_setup(ax, c, lon_r, lat_r, xlabel=False, ylabel=ylabel)

    # ---- ground-track helpers (copy of pub_maps_v2_figures.main)
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    from matplotlib.cm import ScalarMappable
    alt = np.maximum(table["altitude"], 0.1)
    anorm = LogNorm(vmin=max(alt.min(), 0.1), vmax=alt.max())
    vir = plt.get_cmap("viridis")
    no_jump = step <= max(150.0, 5.0 * np.median(step))
    no_seam = (np.abs(np.diff(glon)) <= 180.0) & no_jump
    legend_extra = [Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                           label="Rosetta impact (lon %.2f°, lat %.2f°)" % (impact["lon_deg"], impact["lat_deg"])),
                    Line2D([], [], marker="s", ls="none", ms=4, mfc="white", mec="k", mew=0.6,
                           label="track start, %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                                    table["altitude"][0] / 1000.0))]

    def gt_legend(fig, anchor, hidden, shs):
        hs = list(shs)
        hs.append(Line2D([], [], color="#21918C", lw=2.6, label="ground track, late-visibility window"))
        hs.append(Line2D([], [], color="#21918C", lw=1.1, label="ground track (sub-spacecraft intercept)"))
        if hidden:
            hs.append(Line2D([], [], color="0.45", lw=0.6, ls=(0, (1, 1.5)), label="ground track, hidden"))
        hs.append(Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                         label="1-min marks from window start"))
        hs.extend(legend_extra)
        fig.legend(handles=hs, loc="lower left", bbox_to_anchor=anchor, fontsize=6.3, ncol=2, frameon=False,
                   handlelength=2.2, columnspacing=1.0)

    letters6 = "abcdef"
    man_q = {}
    for Q in Qs:
        codes = qcodes[Q.key]
        for variant, vname in VARIANTS:
            if not wanted(Q.key, variant):
                continue
            stem = "%s_%sUTC_%s" % (Q.key, tag, vname)
            sub = stage / Q.key
            vlabel = "whole surface" if variant == "whole" else "seen facets"
            sdef = whole_def if variant == "whole" else seen_def
            shs = surface_handles_q(Q, variant)
            title_v = "%s, %s UTC, %s" % (Q.title, hhmmss, vlabel)
            qt = qtext(Q)
            common = "%s %s %s %s %s" % (sdef, ill_def, law, numbers_txt, fitind)
            # per-facet colours: templates (relief) and ground tracks (flat)
            rgb_tpl = surface_rgb(Q, codes, seen, I_ll, variant, a.template_relief)
            rgb_gt = surface_rgb(Q, codes, seen, I_ll, variant, a.gt_relief)
            img_tpl, img_gt = {}, {}
            for wname in wins:
                face, ss = faces_ll[wname]
                img_tpl[wname] = V2.shade(face, rgb_tpl, ss)
                img_gt[wname] = V2.shade(face, rgb_gt, ss) if a.gt_relief != a.template_relief else img_tpl[wname]
                tab, txt = panel_table("%s %s, lon/lat %s" % (Q.key, hhmmss, wname), Q, codes, seen, area,
                                       disp_ll[wname], variant)
                tables["%s_%s_latlon_%s" % (Q.key, variant, wname)] = tab
                ttext.append(txt)
            del rgb_tpl
            if a.template_relief > 0:
                rc = V2.relief_check(Q, a.template_relief)
                tables["%s_relief_check" % Q.key] = rc

            # -- templates figure (copy of the pub_maps_v2_figures.main layout; value row = this quantity)
            row_h, hsp = 2.05, 0.30
            rows_in = 3 * row_h + 0.05
            height = (rows_in + 3 * hsp * rows_in / 4 + 1.35) / 0.995
            fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, height))
            gs = fig.add_gridspec(4, 3, width_ratios=[2.0, 1.0, 0.05], height_ratios=[row_h, 0.05, row_h, row_h],
                                  wspace=0.10, hspace=hsp, left=0.07, right=0.93, top=0.995, bottom=1.35 / height)
            letters = iter(letters6)
            axes = {}
            title_row = {"value": title_v, "slope": "effective gravitational slope", "cliff": "cliff templates"}
            for r, row in ((0, "value"), (2, "slope"), (3, "cliff")):
                for c, wname in enumerate(("global", "zoom")):
                    ax = fig.add_subplot(gs[r, c])
                    draw_map(ax, img_tpl[wname] if row == "value" else imgs[(row, wname)], wname, c, ylabel=(c == 0),
                             annotate=(row == "value"))
                    if c == 0:
                        ax.set_title(title_row[row], fontsize=8, loc="left", pad=3)
                    if r == 3:
                        ax.set_xlabel("longitude [°]")
                    ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
                    axes[(row, wname)] = ax
                cax = fig.add_subplot(gs[r, 2])
                if row == "slope":
                    cb = fig.colorbar(ScalarMappable(norm=snorm, cmap=magma), cax=cax, extend="max")
                    cb.set_label("slope [°] (> 90: overhang)", fontsize=6.5)
                    cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
                else:
                    cax.axis("off")
            for wname in ("global", "zoom"):
                fixed_colorbar(fig, V2.under_axes(fig, axes[("value", wname)], 0.33, 0.11), Q,
                               Q.cb_label if wname == "global" else Q.cb_short, narrow=(wname == "zoom"))
            fig.legend(handles=shs + V2.cliff_handles(meta) + V2.track_handles_templates(gt), loc="lower left",
                       bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3, frameon=False, handlelength=2.0,
                       columnspacing=1.2)
            cap_t = ("Illumination of the nucleus at %s UTC on 2016-09-30, the time of the production-rate map of the "
                     "same layout (full sets, production_rate/12M_templates_h2o_production_rate_%sUTC), with the "
                     "effective gravitational slope and the cliff templates, %s. Left: global; right: the impact area "
                     "(box in the left panels). (a, b) Surface: %s; %s; %s %s (c, d) Effective gravitational slope and "
                     "(e, f) cliff templates (as in the production-rate and late-visibility figures). %s %s" % (
                         hhmmss, tag, mesh, qt, vlabel, relief_tpl, common, slp, trk))
            save(fig, sub / ("12M_templates_%s" % stem), cap_t)
            # -- the two value panels as their own figures (copy of the pub_maps_v2_figures.main panel layout)
            for c, wname in enumerate(("global", "zoom")):
                figsize = (FIG.DOUBLE_COLUMN_IN, 4.6) if c == 0 else (4.4, 5.2)
                fig = plt.figure(figsize=figsize)
                ax = fig.add_axes([0.10, 0.36, 0.86, 0.60] if c == 0 else [0.15, 0.36, 0.80, 0.60])
                draw_map(ax, img_tpl[wname], wname, c, ylabel=True, annotate=True)
                ax.set_xlabel("longitude [°]")
                ax.set_title(title_v + ("" if c == 0 else ", impact area"), fontsize=8, loc="left", pad=3)
                fixed_colorbar(fig, V2.under_axes(fig, ax, 0.42, 0.12, 0.9), Q,
                               Q.cb_label if c == 0 else Q.cb_short, narrow=(c == 1))
                fig.legend(handles=shs + V2.track_handles_templates(gt), loc="lower left", bbox_to_anchor=(0.02, 0.0),
                           ncol=2, fontsize=6.0, frameon=False, handlelength=1.8, columnspacing=1.0)
                pc = ("Illumination of the nucleus at %s UTC on 2016-09-30 (time of the production-rate map), %s, %s. "
                      "Surface: %s; %s; %s %s %s" % (hhmmss, "global map" if c == 0 else
                                                     "the impact area (box in the global map)", mesh, qt, vlabel,
                                                     relief_tpl, common, trk))
                save(fig, sub / "panels" / ("12M_templates_%s_%s" % (stem, "a_global" if c == 0 else "b_impact_area")),
                     pc)
            del img_tpl

            # -- ground track, lon/lat (copy of the pub_maps_v2_figures.main layout)
            Hll = 9.6
            fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, Hll))
            ax = fig.add_axes([0.08, 0.665, 0.90, 0.315])
            ax.imshow(img_gt["global"], extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
            FIG._colored_track(ax, glon, glat, alt, anorm, vir, 2.6, mask=late & gt["outside"], seg_mask=no_seam,
                               zorder=5)
            lc = FIG._colored_track(ax, glon, glat, alt, anorm, vir, 1.0, mask=gt["outside"], seg_mask=no_seam, zorder=6)
            ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=9, mfc="#F0E442", mec="k", mew=0.6,
                    zorder=9)
            ax.plot([glon[0]], [glat[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
            ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6,
                    zorder=7)
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            ax.set_xticks(range(-180, 181, 60))
            ax.set_yticks(range(-90, 91, 30))
            ax.set_xlabel("longitude [°]")
            ax.set_ylabel("latitude [°]")
            ax.set_title(title_v, fontsize=8, loc="left", pad=3)
            ax.text(0.01, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
            fixed_colorbar(fig, V2.under_axes(fig, ax, 0.48, 0.11, 0.8), Q, Q.cb_label)
            axz = fig.add_axes([0.08, 0.17, 0.90, 0.31])
            axz.imshow(img_gt["zoom"], extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos",
                       aspect=1.0 / coslat)
            FIG._colored_track(axz, glon, glat, alt, anorm, vir, 2.6, mask=late, seg_mask=no_seam, zorder=5)
            FIG._colored_track(axz, glon, glat, alt, anorm, vir, 1.1, mask=gt["outside"], seg_mask=no_seam, zorder=6)
            FIG._minute_marks(axz, gt, late, lambda p: FIG.lonlat(p), label=False)
            axz.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=11, mfc="#F0E442", mec="k", mew=0.6,
                     zorder=9)
            axz.annotate("Rosetta impact\n%s UTC\nlon %.2f°, lat %.2f°" % (
                impact["utc"].split("T")[1], impact["lon_deg"], impact["lat_deg"]),
                xy=(impact["lon_deg"], impact["lat_deg"]), xytext=(-75, 45), textcoords="offset points",
                fontsize=6.3, ha="center", arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
                bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"}, zorder=10)
            axz.set_xlim(*zl)
            axz.set_ylim(*za)
            axz.set_xlabel("longitude [°]")
            axz.set_ylabel("latitude [°]")
            axz.text(0.01, 0.97, "(b)", transform=axz.transAxes, va="top", fontsize=8,
                     bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
            fixed_colorbar(fig, V2.under_axes(fig, axz, 0.48, 0.11, 1.0), Q, Q.cb_label)
            gt_legend(fig, (0.08, 0.0), False, shs)
            cb = fig.colorbar(lc, cax=fig.add_axes([0.70, 0.03, 0.28, 0.011]), orientation="horizontal")
            cb.set_label("spacecraft altitude [m]", fontsize=7)
            cap_ll = ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
                      "longitude/latitude (equirectangular), outermost surface where overhangs put several surface "
                      "layers at one longitude/latitude, in the layout of the production-rate figure "
                      "12M_groundtrack_h2o_production_rate_%sUTC_latlon. (a) Global map; (b) the impact area (box in a). "
                      "Surface: %s; %s; %s. %s Track: coloured by altitude (log colour bar 'spacecraft altitude'), "
                      "thick within the late-visibility window. %s" % (
                          Path(gt["shape"]).name, M.thousands(len(F)), tag, qt, vlabel, relief_gt, common, trk))
            save(fig, sub / ("12M_groundtrack_%s_latlon" % stem), cap_ll)
            del img_gt

            # -- ground track, 3-D views (copy of the pub_maps_v2_figures.main layout)
            rgb3 = surface_rgb(Q, codes, seen, I3, variant, a.gt_relief)
            H3 = 4.05
            fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, H3))
            width = 0.315
            hpan = width * FIG.DOUBLE_COLUMN_IN / H3
            axs = [fig.add_axes([0.005 + k * (width + 0.0175), 1.0 - hpan - 0.005, width, hpan]) for k in range(3)]
            lc = None
            for k, (ax, rr) in enumerate(zip(axs, r3)):
                image = V2.shade(rr["face"], rgb3, rr["ss"])
                tab, txt = panel_table("%s %s, 3-D panel (%s)" % (Q.key, hhmmss, "abc"[k]), Q, codes, seen, area,
                                       rr["disp"], variant)
                tables["%s_%s_3d_%s" % (Q.key, variant, "abc"[k])] = tab
                ttext.append(txt)
                half, to_image, visible = rr["half"], rr["to_image"], rr["visible"]
                ax.imshow(image, interpolation="lanczos", origin="upper", extent=(0, pixels, pixels, 0))
                del image
                gx, gy, _ = to_image(ground)
                tol = max(3.0, 4.0 * 2.0 * half / pixels)
                vis = visible(ground, tol) & gt["outside"]
                hidden = ~vis & gt["outside"]
                FIG._colored_track(ax, gx, gy, alt, anorm, vir, 0.6, mask=hidden, seg_mask=no_jump, zorder=4,
                                   solid_color="0.45", dotted=True)
                FIG._colored_track(ax, gx, gy, alt, anorm, vir, 2.6, mask=vis & late, seg_mask=no_jump, zorder=5)
                lc = FIG._colored_track(ax, gx, gy, alt, anorm, vir, 1.1, mask=vis, seg_mask=no_jump, zorder=6)
                if vis[0]:
                    ax.plot([gx[0]], [gy[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
                ix, iy, _ = to_image(ip[None, :])
                ax.plot([ix[0]], [iy[0]], marker="*", ms=9 if k < 2 else 11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
                if k < 2:
                    s = panels3[k + 1][1] / half * pixels / 2
                    ax.plot([ix[0] - s, ix[0] + s, ix[0] + s, ix[0] - s, ix[0] - s],
                            [iy[0] - s, iy[0] - s, iy[0] + s, iy[0] + s, iy[0] - s], color="k", lw=0.6, zorder=7)
                if k == 2:
                    FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=True)
                    ax.annotate("Rosetta impact\n%s UTC" % impact["utc"].split("T")[1], xy=(ix[0], iy[0]),
                                xytext=(-45, -55), textcoords="offset points", fontsize=6.3, ha="center",
                                arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
                                bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5,
                                      "boxstyle": "round,pad=0.25"}, zorder=10)
                else:
                    FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=False)
                ax.set_xlim(0, pixels)
                ax.set_ylim(pixels, 0)
                ax.set_xticks([])
                ax.set_yticks([])
                FIG._scale_bar(ax, half, pixels)
                ax.text(0.03, 0.97, "(%s)" % "abc"[k], transform=ax.transAxes, va="top", fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
                fixed_colorbar(fig, V2.under_axes(fig, ax, 0.08, 0.10, 0.78), Q, Q.cb_short, narrow=True)
            del rgb3
            gt_legend(fig, (0.005, 0.0), True, shs)
            cb = fig.colorbar(lc, cax=fig.add_axes([0.70, 0.07, 0.29, 0.03]), orientation="horizontal")
            cb.set_label("spacecraft altitude [m]", fontsize=7)
            cap_3d = ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
                      "orthographic views rendered on the GPU (z-buffer), in the layout of the production-rate figure "
                      "12M_groundtrack_h2o_production_rate_%sUTC_3d. (a) The whole nucleus seen down the final approach "
                      "(view direction: from the start of the late-visibility window toward the impact point); (b) the "
                      "impact region (box in a); (c) close-up of the late descent (box in b), with labelled 1-minute "
                      "marks (UTC, altitude). Scale bars in each panel. Surface: %s; %s; %s. %s Track: coloured by "
                      "altitude (log colour bar 'spacecraft altitude'), thick within the late-visibility window, dotted "
                      "grey where hidden behind the nucleus. %s" % (
                          Path(gt["shape"]).name, M.thousands(len(F)), tag, qt, vlabel, relief_gt, common, trk))
            save(fig, sub / ("12M_groundtrack_%s_3d" % stem), cap_3d)
            man_q["%s_%s" % (Q.key, variant)] = stem
            say("wrote %s (5 figures)" % stem)

    # ---- side by side: production rate | factor on the same facets | factor on the whole surface
    side = None
    if only is None or "side" in only:
        codesF = qcodes[QF.key]
        scales, img_s, sc_check = {}, {}, {}
        pm = json.load(open(pr_man)) if pr_man and pr_man.exists() else None
        for wname in ("global", "zoom"):
            face, ss = faces_ll[wname]
            disp = disp_ll[wname]
            sel = disp[seen[disp] & (F_h2o[disp] > 0)]
            sc = V2.Scale("h2o_%s/%s" % (tag, wname), F_h2o[sel], area[sel], a.low_cut_area)
            scales[wname] = sc
            tab, txt = V2.class_table("h2o_%s, lon/lat %s window (side by side)" % (tag, wname), F_h2o, seen, area,
                                      disp, sc)
            tables["side_h2o_latlon_%s" % wname] = tab
            ttext.append(txt)
            if pm is not None:
                ref = pm["class_tables"]["latlon_" + wname]
                sc_check[wname] = {
                    "edges_equal_production_rate_manifest": [float(e) for e in sc.edges] == ref["scale"]["edges_log10"],
                    "class_facets_equal": [r_["facets"] for r_ in tab["classes"]] == [r_["facets"] for r_ in
                                                                                      ref["classes"]],
                    "displayed_facets_equal": tab["displayed_facets"] == ref["displayed_facets"],
                    "manifest": str(pr_man)}
            img_s[("F", wname)] = V2.shade(face, V2.facet_rgb(F_h2o, seen, sc, I_ll, a.template_relief), ss)
            img_s[("fmask", wname)] = V2.shade(face, same_mask_rgb(QF, codesF, seen, F_h2o, I_ll, a.template_relief), ss)
            img_s[("fwhole", wname)] = V2.shade(face, surface_rgb(QF, codesF, seen, I_ll, "whole", a.template_relief),
                                                ss)
            tab, txt = panel_table("illumination_factor %s, lon/lat %s (side by side, production-rate mask)" % (
                hhmmss, wname), QF, codesF, seen, area, disp, "mask", F=F_h2o)
            tables["side_factor_mask_latlon_%s" % wname] = tab
            ttext.append(txt)
            say("side by side %s: %s; check vs production-rate manifest %s" % (wname, sc.text(), sc_check.get(wname)))
        if pm is not None and not all(all(v for k, v in c.items() if k != "manifest") for c in sc_check.values()):
            raise SystemExit("side-by-side production-rate classes differ from the production-rate figure: %s" % sc_check)
        row_h, hsp = 2.05, 0.30
        rows_in = 3 * row_h + 0.05
        height = (rows_in + 3 * hsp * rows_in / 4 + 1.35) / 0.995
        fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, height))
        gs = fig.add_gridspec(4, 3, width_ratios=[2.0, 1.0, 0.05], height_ratios=[row_h, 0.05, row_h, row_h],
                              wspace=0.10, hspace=hsp, left=0.07, right=0.93, top=0.995, bottom=1.35 / height)
        letters = iter(letters6)
        axes = {}
        titles = {"F": "H₂O production rate, %s UTC (seen facets)" % hhmmss,
                  "fmask": "illumination factor %s on the same facets (seen, F$_\\mathrm{H_2O}$ > 0)" % QSYM,
                  "fwhole": "illumination factor %s, whole surface" % QSYM}
        for r, row in ((0, "F"), (2, "fmask"), (3, "fwhole")):
            for c, wname in enumerate(("global", "zoom")):
                ax = fig.add_subplot(gs[r, c])
                draw_map(ax, img_s[(row, wname)], wname, c, ylabel=(c == 0), annotate=(row == "F"))
                if c == 0:
                    ax.set_title(titles[row], fontsize=8, loc="left", pad=3)
                if r == 3:
                    ax.set_xlabel("longitude [°]")
                ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
                axes[(row, wname)] = ax
            cax = fig.add_subplot(gs[r, 2])
            if row == "F":
                cax.axis("off")
            else:
                fixed_colorbar(fig, cax, QF, "%s (triangle: %s)" % (QSYM, qr()["r"]), narrow=True,
                               orientation="vertical", fontsize=6.5)
        cbl = "H₂O production rate [m⁻² s⁻¹] (seen facets; log classes)"
        for wname, ml in (("global", 7), ("zoom", 5)):
            V2.discrete_colorbar(fig, V2.under_axes(fig, axes[("F", wname)], 0.33, 0.11), scales[wname],
                                 cbl if wname == "global" else cbl.split(" (")[0] + " (impact area)", max_labels=ml)
        fig.legend(handles=V2.surface_handles("F_H2O") + surface_handles_q(QF, "whole") +
                   V2.track_handles_templates(gt), loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3,
                   frameon=False, handlelength=2.0, columnspacing=1.2)

        def scale_caption(sc, where):
            """copy of pub_maps_v2_figures.main scale_caption"""
            return ("%s: the seen non-zero values drawn in this panel range from %.3e to %.3e m^-2 s^-1; %s." % (
                where, sc.vmin, sc.vmax, sc.text().replace("$10^{", "10^").replace("}$", "")) + sc.below_text() +
                (" The smallest value is a denormal (subnormal) float64 number, i.e. an underflow of the model's "
                 "softplus / multiplier product below 2.2e-308, not a physical value." if sc.denormal else ""))
        S2 = nt["seen"]
        nz = S2["nonzero_F"]
        side_numbers = (
            "Numbers at %s UTC (checks/illumination_numbers.txt): of the seen facets with F_H2O > 0 (%s, %.2f %% of the "
            "seen area) %.1f %% of the area is night or cast shadow (day %.1f %%, night %.1f %%, cast shadow %.1f %%; "
            "their area-weighted mean q/q_max %.3f), and %.1f %% of the seen H2O production comes from night or "
            "cast-shadow facets; the production-weighted factor of the seen surface (production / (q_max sum A area)) "
            "is %.3f. On these facets log10 F_H2O spans %.2f decades (area-weighted 2-98 %% range; log10 A %.2f, log10 "
            "q/q_max %.2f): the activity A, not the illumination, sets the colour of the production-rate map." % (
                hhmmss, fmt_area(nz["area_m2"]), 100 * nz["area_fraction_of_scope"],
                100 * nz["night_or_shadow_area_fraction"], 100 * nz["class_area_fraction"]["day"],
                100 * nz["class_area_fraction"]["night"], 100 * nz["class_area_fraction"]["cast_shadow"],
                nz["mean_factor_area_weighted"], 100 * S2["night_or_shadow_production_share"],
                S2["production_weighted_factor"], S2["spread_nonzero_F"]["log10_F"]["range_02_98_decades"],
                S2["spread_nonzero_F"]["log10_A"]["range_02_98_decades"],
                S2["spread_nonzero_F"]["log10_f"]["range_02_98_decades"]))
        n_ = T["facets"]
        cap_s = (
            "Where night and cast shadow lower the H2O production: the production rate of %s next to the illumination "
            "factor of its source rate, at %s UTC on 2016-09-30, %s. Left: global; right: the impact area (box in the "
            "left panels). (a, b) H2O production rate F_H2O (instantaneous emitted H2O molecular flux per facet, m^-2 "
            "s^-1), exactly as in the panels (a, b) of the full set's production_rate/12M_templates_h2o_production_rate_"
            "%sUTC (the same classes, recomputed by pub_maps_v2_figures.Scale with the area cut %g and checked equal to "
            "its manifest %s): colour only on the facets seen during the fit, seen facets with F_H2O = 0 light green, "
            "unseen facets grey; %s %s (c, d) The illumination factor q_H2O/q_max on the same facets (seen, F_H2O > 0; "
            "seen facets with F_H2O = 0 light green and unseen facets grey as in a, b): %s on day facets "
            "in %d classes of 0.1 (matplotlib 'YlOrBr' reversed), %s on night and cast-shadow facets (the colour of the "
            "triangle at the low end of the colour bar). (e, f) The same factor on the whole surface. All panels: "
            "%s. Facets at this time: %s day, %s night, %s in cast shadow. %s %s %s %s %s %s" % (
                a.model_name, hhmmss, mesh, tag, a.low_cut_area, pr_man.name if pr_man else "(not checked)",
                scale_caption(scales["global"], "(a)"), scale_caption(scales["zoom"], "(b)"),
                qr()["day"], qr()["n_f"], qr()["r"], relief_tpl,          # 2026-11-18: the factor law from QRATIO
                M.thousands(n_["day"]), M.thousands(n_["night"]), M.thousands(n_["cast_shadow"]), side_numbers,
                ill_def, law, fitind, M.model_text(ctx), trk))
        save(fig, stage / "side_by_side" / ("12M_side_by_side_h2o_production_rate_vs_illumination_factor_%sUTC" % tag),
             cap_s)
        side = {"scales": {w: {"edges_log10": [float(e) for e in s_.edges], "step": s_.step, "cut_value": s_.cut_value}
                           for w, s_ in scales.items()}, "check_vs_production_rate_manifest": sc_check}
        say("wrote the side-by-side figure")

    # ---- checks, manifest, captions; move into place
    (stage / "checks").mkdir(exist_ok=True)
    with open(stage / "checks" / ("illumination_%s_class_tables.txt" % tag), "w") as fh:
        fh.write("# illumination_maps_figures.py (%s), %s UTC: facets drawn in every panel (rasterizer face index per "
                 "supersampled pixel, each facet once) per class / category; seen = the fit's nonneg_facets.npz seen\n\n"
                 % (SCHEMA, utc))
        fh.write("\n\n".join(ttext) + "\n")
    man = {"schema": SCHEMA, "created_utc": stamp, "host": socket.gethostname(),
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
           "script_sha256": M.sha256_file(__file__), "helper_sha256": {
               "pub_maps_v2_figures.py": M.sha256_file(V2.__file__), "pub_maps_figures.py": M.sha256_file(M.__file__),
               "figures.py": M.sha256_file(FIG.__file__), "nonneg_maps.py": M.sha256_file(a.maps_code)},
           "time_utc": utc, "label": lab, "field_checks": fields["checks"], "factor_checks": fchk,
           "class_checks": qchecks, "quantities": {Q.key: {"edges": [float(e) for e in Q.edges],
                                                            "colours": [V2.hexcol(c) for c in Q.colors],
                                                            "specials": [(cd, lb, V2.hexcol(c)) for cd, lb, c in
                                                                         Q.specials],
                                                            "under": None if Q.under is None else V2.hexcol(Q.under)}
                                                   for Q in Qs},
           "figures": man_q, "side_by_side": side, "tables": tables, "template_relief": a.template_relief,
           "gt_relief": a.gt_relief, "numbers_json": str((ndir / "illumination_numbers.json").resolve()),
           "crosscase_json": str((ndir / "illumination_crosscase.json").resolve()),
           "files": [str(Path(p).relative_to(stage)) for p in written], "wall_seconds": time.time() - T0}
    json.dump(man, open(stage / "checks" / ("illumination_%s_manifest.json" % tag), "w"), indent=1, default=float)
    with open(stage / "checks" / ("illumination_%s_captions.json" % tag), "w") as fh:
        json.dump([{"figure": n, "caption": c} for n, c in captions], fh, indent=1, ensure_ascii=False)
    move_into_place(stage, out)
    say("done: %d figure files" % len(written))
    return 0


# ============================================================================================ docs
def cmd_docs(a):
    out = Path(a.out_dir)
    for name in ("README.txt", "HOW-TO-RERUN.md", "CAPTIONS.md"):
        if (out / name).exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % (out / name))
    if (out / "code").exists():
        raise SystemExit("%s exists" % (out / "code"))
    nb = json.load(open(out / "checks" / "illumination_numbers.json"))
    set_qratio(nb["qratio_h2o"])                     # 2026-11-18 (numbers written before the fix have no key: refused)
    cr = json.load(open(out / "checks" / "illumination_crosscase.json"))
    mans = {t: json.load(open(out / "checks" / ("illumination_%s_manifest.json" % t.replace(":", ""))))
            for t, _, _ in MAP_TIMES}
    provs = sorted(out.glob("provenance-*"))
    # ---- code/
    code = out / "code"
    code.mkdir()
    for p in [Path(__file__), Path(V2.__file__), Path(M.__file__), Path(FIG.__file__), Path(a.maps_code),
              HERE / "collect_captions.py"] + [Path(x) for x in a.job_script]:
        shutil.copy2(p, code / p.name)
    git = next((p / "git.txt" for p in provs if (p / "git.txt").exists()), None)
    if git:
        shutil.copy2(git, code / "git.txt")
    with open(code / "SHA256SUMS", "w") as fh:
        for p in sorted(code.iterdir()):
            if p.name != "SHA256SUMS":
                fh.write("%s  %s\n" % (M.sha256_file(p), p.name))
    # ---- CAPTIONS.md (collect_captions.py, unchanged)
    subprocess.run([sys.executable, str(HERE / "collect_captions.py"), str(out), "--title",
                    "Figure captions: 12M illumination maps at 10:33:10 and 10:39:26 UTC (2016-09-30)", "--order",
                    "side_by_side,illumination_class,cos_sza,illumination_factor,illumination_factor_log"], check=True)
    # ---- README.txt
    figs = sorted(str(p.relative_to(out))[:-len("_caption.txt")] for p in out.rglob("*_caption.txt")
                  if not any(s in str(p) for s in (".incomplete_", "superseded", "provenance-", "/code/")))
    L = []
    L.append("Illumination maps of the 12M nucleus at the production-rate map times 10:33:10 and 10:39:26 UTC (2016-09-30)")
    L.append("=" * len(L[0]))
    L.append("Written %s (%s UTC) by figures_cluster/scripts/illumination_maps_figures.py docs." % (
        central(utc_stamp()), utc_stamp()))
    L.append("")
    L.append("User (2026-11-17, CST): \"check the production rate plots -- i dont see the effect of the shadow or night "
             "side decrease of the productions -- is it there? generate illumination maps for the time intervals that "
             "corresponds to that on the production rate maps\".")
    L.append("")
    L.append("ANSWER")
    L.extend(answer_lines(nb, cr, mans))
    L.append("")
    L.append("NUMBERS (computed from the data by the numbers job; checks/illumination_numbers.json / .txt)")
    for hhmmss, t in nb["times"].items():
        L.append("  %s UTC (r_h %.6f AU; q_max %.4f, q_min = %s q_max):" % (t["utc"], t["rh_AU"], t["q_H2O"][1], qr()["r"]))
        for scope, name in (("whole", "whole surface"), ("seen", "seen surface (the fit's seen facets)")):
            s = t[scope]
            L.append("    %s, %.4e m^2:" % (name, s["area_m2"]))
            L.append("      area per class: " + "; ".join("%s %.4e m^2 (%.2f %%)" % (
                n.replace("_", " "), s["classes"][n]["area_m2"], 100 * s["classes"][n]["area_fraction"])
                for n in CLASS_KEYS))
            L.append("      area-weighted mean q/q_max: %.4f (day facets alone %.4f)" % (
                s["mean_factor_area_weighted"], s["mean_factor_day_area_weighted"]))
            nz = s["nonzero_F"]
            L.append("      F_H2O > 0 on %.4e m^2 (%.2f %% of the area); of it night or cast shadow %.2f %% (night %.2f "
                     "%%, shadow %.2f %%)" % (nz["area_m2"], 100 * nz["area_fraction_of_scope"],
                                             100 * nz["night_or_shadow_area_fraction"],
                                             100 * nz["class_area_fraction"]["night"],
                                             100 * nz["class_area_fraction"]["cast_shadow"]))
            L.append("      production %.4e s^-1, from night or cast-shadow facets %.2f %% (night %.2f %%, shadow %.2f "
                     "%%); production-weighted q/q_max %.4f" % (
                         s["production_s-1"], 100 * s["night_or_shadow_production_share"],
                         100 * s["production_share"]["night"], 100 * s["production_share"]["cast_shadow"],
                         s["production_weighted_factor"]))
            sp = s["spread_nonzero_F"]
            L.append("      spread over F > 0 (area-weighted 2-98 %% range): log10 F %.2f, log10 A %.2f, log10 q/q_max "
                     "%.2f decades" % (sp["log10_F"]["range_02_98_decades"], sp["log10_A"]["range_02_98_decades"],
                                       sp["log10_f"]["range_02_98_decades"]))
        L.append("    seen production / Q_H2O = %.4f (Q_H2O %.4e s^-1)" % (t["seen_production_share_of_Q_H2O"],
                                                                         t["Q_H2O_s-1"]))
    b = nb["between_times"]
    L.append("  Between the two times the Sun direction moves by %.2f deg in the body frame; %.2f %% of the whole surface "
             "(%.2f %% of the seen surface) changes its illumination class." % (
                 b["sun_direction_angle_deg"], 100 * b["illumination_changes"]["whole"]["changed_area_fraction"],
                 100 * b["illumination_changes"]["seen"]["changed_area_fraction"]))
    L.append("")
    L.append("CROSS-CASE CHECK (checks/illumination_crosscase.txt)")
    L.append("  cos_SZA and illumination of the maps (.vtu) of %d cases (%s) at both times, against %s: bitwise "
             "identical on every facet: %s. Geometry folders (%s): t11 / t12 cos.f8 and shadow.u1 identical: %s. Seen "
             "masks identical: %s%s." % (
                 len(cr["cases"]), ", ".join(cr["cases"]), cr["reference"], cr["all_cos_SZA_and_illumination_identical"],
                 ", ".join(cr["geometry_dirs"]), cr["all_geometry_files_identical"], cr["all_seen_identical"],
                 "" if cr["all_seen_identical"] else " (" + "; ".join(
                     "%s: %d facets differ" % (k, v["seen_facets_differing_from_reference"])
                     for k, v in cr["cases"].items() if not v["seen_equal_reference"]) + ")"))
    L.append("  So the whole-surface maps serve every case%s." % (
        "; the seen-only maps too (same seen mask)" if cr["all_seen_identical"] else
        "; the seen-only maps use the seen mask of %s (the other cases' seen masks differ)" % cr["reference"]))
    L.append("")
    L.append("FOLDERS (every figure .eps + .png 300 dpi + _caption.txt; all captions in CAPTIONS.md)")
    L.append("  side_by_side/           per time: (a, b) the production rate of %s exactly as in its production-rate "
             "templates figure, (c, d) q/q_max on the same facets, (e, f) q/q_max on the whole surface" % nb[
                 "reference"]["key"])
    L.append("  illumination_class/     day / night / cast shadow")
    L.append("  cos_sza/                cos SZA of the day facets (0.1 classes), night and cast shadow own colours")
    L.append("  illumination_factor/    q/q_max = %s (night, shadow), %s (day); 0.1 classes, range %s-1" % (
        qr()["r"], qr()["day"], qr()["r"]))
    L.append("  illumination_factor_log/ log10 q/q_max, 0.05-decade classes from %.2f to 0 (to compare with the "
             "production-rate classes)" % qr()["lo_l"])
    L.append("    each: 12M_templates_<q>_<hhmmss>UTC_<whole_surface|seen_only> (value row + slope + cliffs, relief 0.7), "
             "panels/..._a_global, _b_impact_area, 12M_groundtrack_..._latlon, ..._3d (flat colours), i.e. every layout "
             "of full_set_*/production_rate/ (the slope and cliff panels c-f are not repeated: they are those of the "
             "production-rate folders)")
    L.append("  checks/                 numbers, cross-case check, class tables per panel, manifests, captions json")
    L.append("  code/                   the code that ran (SHA256SUMS; git.txt = AMPS HEAD and status at the jobs)")
    L.append("  provenance-<job>/       job script, log, git, module list, nvidia-smi -q (GPU), toolchain, pip list, "
             "commands.sh")
    L.append("")
    L.append("JOBS: %s" % "; ".join(p.name.replace("provenance-", "") for p in provs))
    if a.note:
        L.append("")
        L.append("NOTES")
        L.extend("  * " + x for x in a.note)
    L.append("")
    L.append("FIGURES (%d)" % len(figs))
    L.extend("  " + f_ for f_ in figs)
    (out / "README.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    # ---- HOW-TO-RERUN.md
    H = ["# How to rerun: 12M illumination maps at 10:33:10 and 10:39:26 UTC", "",
         "Written %s by figures_cluster/scripts/illumination_maps_figures.py docs. Code that ran: code/ (SHA256SUMS); "
         "AMPS tree at the jobs: code/git.txt. Version-controlled copies: "
         "srcCG/LiouvilleCOPS_cuda/model_ext/figures/illumination_maps_figures.py and model_ext/jobs/"
         "illumination_{numbers,maps}.sbatch (RUN_COPIES.sha256)." % central(utc_stamp()), "",
         "## A. Submit again (the tools refuse existing products: move this folder aside as *.superseded_<UTC> first)",
         "", "```bash", "cd /data001/highresolu_vt/cops_output/jobs",
         "J=$(sbatch --parsable illumination_numbers.sbatch)",
         "sbatch --dependency=afterok:$J --export=ALL,TIME=10:33:10 illumination_maps.sbatch",
         "sbatch --dependency=afterok:$J --export=ALL,TIME=10:39:26 illumination_maps.sbatch",
         "# afterwards on the login node (light): README.txt, CAPTIONS.md, HOW-TO-RERUN.md, code/",
         "%s" % " ".join(sh_quote(x) for x in a.docs_command), "```", "",
         "## B. The exact commands as run (provenance-<job>/commands.sh)", ""]
    for p in provs:
        H += ["### %s" % p.name, "", "```bash", (p / "commands.sh").read_text().strip() if (p / "commands.sh").exists()
              else "(no commands.sh)", "```", ""]
    H += ["## C. Inputs", "",
          "* reference state (%s): %s; maps %s; geometry %s" % (nb["reference"]["key"], nb["reference"]["ft_dir"],
                                                               nb["reference"]["maps_dir"], nb["reference"]["geometry"]),
          "* cross-case maps: " + "; ".join("%s %s" % (k, v["maps_dir"]) for k, v in cr["cases"].items()),
          "* environment: as in provenance-<job>/commands.sh (exports) and env.txt", ""]
    (out / "HOW-TO-RERUN.md").write_text("\n".join(H) + "\n", encoding="utf-8")
    say("docs written: README.txt, CAPTIONS.md, HOW-TO-RERUN.md, code/ (%d figures)" % len(figs))
    return 0


def sh_quote(s):
    import shlex
    return shlex.quote(str(s))


def answer_lines(nb, cr, mans):
    """The answer to the user's question, from the numbers (2026-11-18: rewritten after the numbers job 136359 showed
    that most of the active area is dark at these times)."""
    t1, t2 = MAP_TIMES[0][0], MAP_TIMES[1][0]
    a1, a2 = nb["times"][t1], nb["times"][t2]
    s1, w1, u1 = a1["seen"], a1["whole"], a1["unseen"]
    side = mans[t1].get("side_by_side") or {}
    pstep = side["scales"]["global"]["step"] if side.get("scales") else None
    L = []
    L.append("  Yes. The night-side and cast-shadow decrease is in the model and in the production-rate maps, but the "
             "maps cannot show it as a visible step:")
    k = qr()                                         # 2026-11-18: the factor law from QRATIO
    L.append("  * Every night or cast-shadow facet has F_H2O = q_min A = %s q_max A, a day facet (%s) "
             "q_max A (checked on every facet at both times: F_H2O = q_H2O A bitwise, q_min/q_max = %s). Illumination "
             "changes F_H2O by at most a factor %g (%.2f decade) at a given activity A%s." % (
                 k["r"], k["day"], k["r"], k["fac"], k["dec"],
                 (", i.e. about %.1f of the %g-decade colour classes of the lon/lat production-rate maps" % (
                     k["dec"] / pstep, pstep)) if pstep else ""))
    L.append("  * The colour of the production-rate maps is set by the activity A, not by the illumination: on the seen "
             "facets with F_H2O > 0 at %s UTC, log10 F_H2O spans %.2f decades (area-weighted 2-98 %% range), log10 A "
             "%.2f decades and log10 q/q_max only %.2f decade; and A = 0 (light green) on %.1f %% of the seen area, "
             "where illumination cannot show at all. The maps do not mark night or shadow, so a dark facet with a large "
             "A looks like a lit facet with a smaller A." % (
                 t1, s1["spread_nonzero_F"]["log10_F"]["range_02_98_decades"],
                 s1["spread_nonzero_F"]["log10_A"]["range_02_98_decades"],
                 s1["spread_nonzero_F"]["log10_f"]["range_02_98_decades"],
                 100 * (1 - s1["nonzero_F"]["area_fraction_of_scope"])))
    L.append("  * At these times most of the active surface is dark, so the decrease applies to most of the production: "
             "%.1f %% (%s UTC: %.1f %%) of the seen area with F_H2O > 0 is night or cast shadow and gives %.1f %% "
             "(%.1f %%) of the seen production. On the whole surface night or cast shadow covers %.1f %% of the area "
             "and %.1f %% of the area with F_H2O > 0, and gives %.1f %% of Q_H2O (%s UTC: %.1f %%), mostly from facets "
             "no fitted observation sees (they give %.1f %% of Q_H2O, %.1f %% of it from their night or cast-shadow "
             "facets; the data do not constrain their activity)." % (
                 100 * s1["nonzero_F"]["night_or_shadow_area_fraction"], t2,
                 100 * a2["seen"]["nonzero_F"]["night_or_shadow_area_fraction"],
                 100 * s1["night_or_shadow_production_share"], 100 * a2["seen"]["night_or_shadow_production_share"],
                 100 * w1["night_or_shadow_area_fraction"], 100 * w1["nonzero_F"]["night_or_shadow_area_fraction"],
                 100 * w1["night_or_shadow_production_share"], t2, 100 * a2["whole"]["night_or_shadow_production_share"],
                 100 * (1 - a1["seen_production_share_of_Q_H2O"]), 100 * u1["night_or_shadow_production_share"]))
    L.append("  * Hence the production-weighted q/q_max (production / (q_max sum A area)) is %.3f on the whole surface "
             "and %.3f on the seen surface at %s UTC: the model's Q_H2O is %.0f %% of what the same activity would give "
             "with the Sun overhead everywhere. The area-weighted mean q/q_max is %.3f (whole) and %.3f (seen)." % (
                 w1["production_weighted_factor"], s1["production_weighted_factor"], t1,
                 100 * w1["production_weighted_factor"], w1["mean_factor_area_weighted"],
                 s1["mean_factor_area_weighted"]))
    zt = (mans[t1].get("tables") or {}).get("side_factor_mask_latlon_zoom")
    if zt and zt["coloured_area_m2"] > 0:
        tot = zt["coloured_area_m2"]
        dark = sum(r["area_m2"] for r in zt["rows"] if r["code"] < 0)
        hi = sum(r["area_m2"] for r in zt["rows"] if r["code"] >= int(round((0.8 - qr()["lo_f"]) / 0.1)))  # f >= 0.8
        L.append("  * Where: in the side-by-side figures (side_by_side/) the seen active facets drawn in the impact-area "
                 "window at %s UTC (%.3e m^2) are %.1f %% night or cast shadow and %.1f %% at q/q_max >= 0.8 (class "
                 "table in checks/illumination_%s_class_tables.txt); in the global window the long band of seen active "
                 "facets at longitude about -60 deg is drawn in the darkest factor colours (the lowest q/q_max classes: night, "
                 "cast shadow or grazing Sun) at both times." % (
                     t1, tot, 100 * dark / tot, 100 * hi / tot, t1.replace(":", "")))
    L.append("  * The two map times are %.2f deg of Sun motion apart: %.2f %% of the surface (%.2f %% of the seen "
             "surface) changes its illumination class between them, so the two production-rate maps look alike." % (
                 nb["between_times"]["sun_direction_angle_deg"],
                 100 * nb["between_times"]["illumination_changes"]["whole"]["changed_area_fraction"],
                 100 * nb["between_times"]["illumination_changes"]["seen"]["changed_area_fraction"]))
    if cr["all_cos_SZA_and_illumination_identical"]:
        L.append("  * Illumination does not depend on the fit: cos SZA and the illumination class are bitwise identical "
                 "in the maps of all %d full-set cases, so these maps serve every case." % len(cr["cases"]))
    return L


# ============================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("numbers", "figures", "docs"):
        p = sp.add_parser(name)
        p.add_argument("--out-dir", required=True)
        if name == "docs":
            p.add_argument("--maps-code", required=True)
            p.add_argument("--job-script", action="append", default=[], help="job scripts to copy into code/")
            p.add_argument("--note", action="append", default=[], help="free-text line for README.txt NOTES (2026-11-18)")
            continue
        for x in ("--base", "--ft-dir", "--fit-dir", "--maps-dir", "--geometry", "--maps-code"):
            p.add_argument(x, required=True)
        p.add_argument("--model-name", default=M.MODEL_NAME)
        if name == "numbers":
            p.add_argument("--ref-key", required=True)
            p.add_argument("--case", action="append", required=True, help="KEY=MAPS_DIR (the reference among them)")
        else:
            for x in ("--track-base", "--bundle"):
                p.add_argument(x, required=True)
            p.add_argument("--time", required=True, help="HH:MM:SS, one of %s" % [t[0] for t in MAP_TIMES])
            p.add_argument("--pr-dir", default=None, help="the reference's full-set folder (production_rate classes "
                                                          "are checked against its checks/h2o_<hhmmss>_manifest.json)")
            p.add_argument("--numbers-dir", default=None, help="folder of illumination_numbers.json (default <out>/checks)")
            p.add_argument("--classes", default="15,20")
            p.add_argument("--low-cut-area", type=float, default=0.02)
            p.add_argument("--template-relief", type=float, default=0.7)
            p.add_argument("--gt-relief", type=float, default=0.0)
            p.add_argument("--gravity-cache", default=None)
            p.add_argument("--device", type=int, default=0)
            p.add_argument("--pixels3d", type=int, default=800)
            p.add_argument("--only", default=None, help="test: comma list of quantity[:whole|seen] and/or 'side'")
    a = ap.parse_args()
    if a.cmd == "docs":
        a.docs_command = [sys.executable] + sys.argv
    return {"numbers": cmd_numbers, "figures": cmd_figures, "docs": cmd_docs}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
