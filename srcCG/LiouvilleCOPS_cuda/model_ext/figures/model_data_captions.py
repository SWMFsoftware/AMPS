#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/model_data_captions.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/model_data_captions.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 751f2adb8a26db07220955f4033919d88fc947537c55091c837dd6b23efa767e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure captions. The _caption.txt files of a full set's model_data/ figures (model_obs_figures.py and
#   tools/pub_figures.py pressures) computed from the set's tables and the fits' result files (2026-11-13).
# RUN: login node: model_data_captions.py --model-data <set>/model_data --tag <tag> --ft-dir <c5/ft/final> --ref-dir
#   <D2> --model-text '...' --job-set4 <job> --job-pressures <job>
# ==== end of header ====
"""Captions (<figure>_caption.txt) of the model / data figures of a full figure set (model_data/: model_obs_figures.py
*_vs_COPS_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG} and tools/pub_figures.py *_pressures_{time,altitude}),
computed from the set's tables and the fits' result files (2026-11-13; the first set's captions of 2026-11-12 have the
same wording with that fit's numbers). Login node; reads small files; refuses existing caption files.

  model_data_captions.py --model-data <set>/model_data --tag softplus_smm1_5cycles --ft-dir <c5/ft/final>
      --ref-dir <clip_refit/D2/final> --model-text "the softplus activity model ... after 5 fit and fine-tuning cycles"
      --job-set4 <job> --job-pressures <job> [--bundle <bundle>]
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cops_cuda import figures as FIG  # noqa: E402
import model_obs_figures as MO       # noqa: E402


def crossing(tab):
    alt = tab["altitude"]
    o = np.argsort(-alt)
    c = MO.crossings(np.log(alt[o]), np.log(tab["pNG"][o] / tab["pRG"][o]))
    return ", ".join("%.0f m" % np.exp(a) for a, _, _ in c) or "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("--model-data", "--tag", "--ft-dir", "--ref-dir", "--model-text", "--job-set4", "--job-pressures"):
        ap.add_argument(k, required=True)
    ap.add_argument("--bundle", default="/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/bundle/"
                                        "restart_templates.npz")
    a = ap.parse_args()
    md = Path(a.model_data)
    fit = FIG.read_table(md / "tables" / ("%s_final.dat" % a.tag))
    ref = FIG.read_table(md / "tables" / "clipped_refit_reference.dat")
    ft = (Path(a.ft_dir) / "finetune_result.txt").read_text()
    g = lambda k: float(re.search(r"^%s\s+initial\s+\S+\s+final\s+(\S+)" % k, ft, re.M).group(1))   # noqa: E731
    fr = (Path(a.ref_dir) / "fit_result.txt").read_text()
    h = lambda k: float(re.search(r"^%s\s+(\S+)" % k, fr, re.M).group(1))                            # noqa: E731
    cf, nf_, rf = g("cost_data"), g("rms_ln_NG"), g("rms_ln_RG")
    cr, nr, rr = h("final_data_cost"), h("rms_ln_residual_NG"), h("rms_ln_residual_RG")
    win = FIG.late_window(fit, a.bundle)
    gaps = [(hi, lo, t0, t1) for hi, lo, t0, t1 in __import__("cliff_importance_figures").ng_gaps(fit, fit["altitude"])]
    gap = gaps[0] if gaps else None
    nng = int(np.sum((fit["refNG"] > 0) & (fit["useNG"] > 0)))
    nrg = int(np.sum((fit["refRG"] > 0) & (fit["useRG"] > 0)))
    span = "%s–%s UTC (gauge altitude %.0f–%.1f m)" % (win["start"].split("T")[1], win["stop"].split("T")[1],
                                                       win["altitude_max_m"], win["altitude_min_m"])
    cross = ("NG = RG crossings: data between %.0f m and %.0f m altitude (%s–%s UTC, a gap in the NG data during which "
             "only RG was measured; hatched), the fitted model at %s, clipped reference model at %s (open circles). " % (
                 gap[0], gap[1], gap[2].strftime("%H:%M:%S"), gap[3].strftime("%H:%M:%S"), crossing(fit),
                 crossing(ref)) if gap else "")
    base = ("Model base: the SHAP8 12M shape model (cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf, 12,364,290 facets), "
            "collisionless (Liouville) transport from the surface to the gauges, the COPS gauges 1.5 m from the spacecraft "
            "reference point toward the comet; the model at the observation times, except for the 95 observations from "
            "10:33:10 UTC, where the model is averaged over the gauges' time response (10-s running mean and 5-readout "
            "average, from a 1-s model grid) and the NG model includes the fitted grain term at the observation time. "
            "Only observation rows are drawn (the model-only 1-s grid rows are removed by "
            "figures_cluster/scripts/positive_table.py; tables/). ")
    data = ("Data: %d NG and %d RG COPS observations, all used in both fits (filled dots). Blue solid: %s (final state, "
            "its last fine tuning applied), final data cost %.4e, rms ln(model/data) NG %.4f, RG %.4f. Orange dashed: the "
            "clipped activity model refitted on the same base, for reference, final data cost %.4e, rms ln NG %.4f, RG "
            "%.4f (data costs from finetune_result.txt / fit_result.txt). Shaded: the late-visibility window %s, the "
            "interval of the late-visibility templates. Altitude: of the gauge position above the shape model. " % (
                nng, nrg, a.model_text, cf, nf_, rf, cr, nr, rr, span))
    made = "Made with figures_cluster/scripts/model_obs_figures.py (job %s)." % a.job_set4
    alt_intro = ("COPS nude-gauge (NG) and ram-gauge (RG) pressures of Rosetta's final descent on 2016-09-30 against the "
                 "altitude above the shape model (logarithmic; the descent runs from left to right), with two model fits. ")
    stem = "12M_%s" % a.tag
    caps = {
        "%s_vs_COPS_time" % stem: (
            "COPS nude-gauge (NG) and ram-gauge (RG) pressures of Rosetta's final descent on 2016-09-30 against UTC, with "
            "two model fits. Left: the whole descent, from 00:00:04 UTC (16.7 km) to the last observation at 10:39:26 "
            "UTC; right: the last 25 minutes. (a, b) NG pressure, (c, d) RG pressure [Pa], each model with its line "
            "style; (e, f) NG and RG together: NG solid, RG dashed, one colour per model, NG data black, RG data grey. "
            + cross + base + data + made),
        "%s_vs_COPS_altitude" % stem: (
            alt_intro + "(a) NG pressure, (b) RG pressure [Pa], each model with its line style; (c) NG and RG together: NG "
            "solid, RG dashed, one colour per model, NG data black, RG data grey. " + cross + base + data + made),
    }
    for g_, pnl in (("NG", "(a)"), ("RG", "(b)"), ("NG_RG", "(c)")):
        what = {"NG": "COPS NG pressure [Pa]", "RG": "COPS RG pressure [Pa]",
                "NG_RG": "COPS NG and RG pressures [Pa] (NG solid, RG dashed)"}[g_]
        caps["%s_vs_COPS_altitude_%s" % (stem, g_)] = (
            "%s of Rosetta's final descent on 2016-09-30 against the altitude above the shape model (logarithmic; the "
            "descent runs from left to right), with two model fits (panel %s of %s_vs_COPS_altitude as a separate "
            "figure). %s%s%s%s" % (what, pnl, stem, cross if g_ == "NG_RG" else "", base, data, made))
    pr = ("Data: %d NG and %d RG COPS observations, all used in both fits (black dots). Blue solid: %s (final state), "
          "final data cost %.4e. Orange dashed: the clipped activity model refitted on the same base, for reference, final "
          "data cost %.4e. Model base: the SHAP8 12M shape model (12,364,290 facets), collisionless (Liouville) transport, "
          "the COPS gauges 1.5 m from the spacecraft reference point toward the comet; for the 95 observations from "
          "10:33:10 UTC the model is averaged over the gauges' time response (10-s running mean and 5-readout average) and "
          "the NG model includes the fitted grain term. Only observation rows are drawn (tables/). The unfitted base model "
          "is not drawn (--no-base). Shaded: the late-visibility window %s. Made with tools/pub_figures.py pressures "
          "(job %s)." % (nng, nrg, a.model_text, cf, cr, span, a.job_pressures))
    caps["%s_pressures_time" % stem] = (
        "COPS NG and RG pressures of Rosetta's final descent on 2016-09-30 against UTC, with two model fits, in the format "
        "of the tracked tools/pub_figures.py pressures (unchanged). Left: the whole descent from 00:00:04 UTC; right: the "
        "last 25 minutes. (a, b) NG pressure [Pa], (c, d) RG pressure [Pa], (e, f) NG / RG ratio. " + pr)
    caps["%s_pressures_altitude" % stem] = (
        "COPS NG and RG pressures of Rosetta's final descent on 2016-09-30 against the altitude above the shape model "
        "(logarithmic), with two model fits, in the format of the tracked tools/pub_figures.py pressures (unchanged). "
        "(a) NG pressure [Pa], (b) RG pressure [Pa], (c) NG / RG ratio. " + pr)
    for name, cap in caps.items():
        if not (md / (name + ".png")).exists():
            raise SystemExit("figure %s.png missing" % name)
        p = md / (name + "_caption.txt")
        with open(p, "x") as f:
            f.write(cap + "\n")
        print("wrote %s" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
