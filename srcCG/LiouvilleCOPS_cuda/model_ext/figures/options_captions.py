#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/options_captions.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/options_captions.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: bc7c3e3db8cbbf5e7b2a87db34809af32ef1d7723f3505ee47944bc203563418
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: captions of the stage-1 option figures (12M_options_<A|B|C|D|best>_*) from summary.csv (2026-11-14; v2
#   convergence, C GPU check; 2026-11-15: --texts, continuation and P_m details, partial C GPU checks)
# RUN: called by figures_12m_nonneg_options_stage1.sh: options_captions.py --out-dir <folder> [--texts <json>]
# ==== end of header ====
"""Captions (<figure>_caption.txt) of the option-comparison figures (figures_12m_nonneg_options_stage1.sh, 2026-11-14) from
the folder's runs.json, summary.csv and best.json. Refuses existing caption files.
    options_captions.py --out-dir <folder>
2026-11-14 (v2): when summary.csv has the "convergence" column (options_table.py --diag), each fit's description adds its
convergence verdict and the D "under diagnosis" sentence is replaced by the diagnosis source; option C's text states the GPU
cross-check from each included C run's gpu_eval/gpu_check.json ("pending" while a run has none).
2026-11-15 (final): --texts <json> replaces option texts (keys A-D), COMMON, CONV (sentence on the convergence column) and
D_NOTE (the D sentence of the per-option caption) - default: the built-in texts; a fit with prev_* columns (continuation)
names the stage-1 fit it continued; a C fit with penalty_mult names its multiplier penalty; when only some C runs have a
passed GPU check, the text names them.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

OPT = {"A": "option A, sum of squares of rank 3 (A = m_b sum_k (u_k)^2, three fields)",
       "B": "option B, background plus softplus (A = m_b [A_bg A_orig + softplus], A_orig the srcCG activity map, A_bg >= the "
            "setting)",
       "C": "option C, local non-negative basis (A = m_b sum_j c_j phi_j on a 150-m node grid, c_j >= 0, smoothness penalty "
            "lambda_C = f C_ref / R_pilot; CPU fits whose COPS_model.dat comes from the design matrices, GPU cross-check "
            "pending)",
       "D": "option D, log form with a spread penalty (A = m_b exp(u), penalty f C_ref Var_w(u), C_ref = the softplus c9 data "
            "cost)"}
COMMON = ("Stage-1 fits: one fit per setting from the clipped refit D2, the six slope-bin multipliers (40-50 ... 90-180 deg) "
          "bounded below by 1, the D2 fine-tuning templates on, max_nfev 4000, no fine-tuning cycle. Reference: the bounded "
          "softplus model with slope multipliers >= 1 after 9 fit and fine-tuning cycles (c9; solid), the current best fit.")
CONV = ("Convergence verdicts: CPU diagnosis of the end points (checks/stage1_diag_20261114: scaled gradient, Gauss-Newton "
        "predicted decrease and trial steps; numbers in summary.txt); option C not diagnosed. ")
SUFFIX = {
    "_time": "against time: (a) NG, (b) RG, (c) NG and RG together (one colour per fit, NG solid, RG dashed; crossings NG = RG "
             "marked); full descent on the left, the last 25 minutes on the right",
    "_altitude": "against the gauge altitude (logarithmic, descent from left to right): (a) NG, (b) RG, (c) NG and RG together "
                 "(one colour per fit, NG solid, RG dashed; crossings NG = RG marked)",
    "_altitude_NG": "NG against the gauge altitude (logarithmic, descent from left to right)",
    "_altitude_RG": "RG against the gauge altitude (logarithmic, descent from left to right)",
    "_altitude_NG_RG": "NG and RG together against the gauge altitude (one colour per fit, NG solid, RG dashed; every crossing "
                       "NG = RG of a model marked with an open circle, the crossings listed in the legend)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--texts", default=None)
    a = ap.parse_args()
    global COMMON, CONV
    d_note = None
    if a.texts:
        tx = json.load(open(a.texts))
        OPT.update({k: v for k, v in tx.items() if k in OPT})
        COMMON, CONV = tx.get("COMMON", COMMON), tx.get("CONV", CONV)
        d_note = tx.get("D_NOTE")
    O = Path(a.out_dir)
    runs = json.load(open(O / "runs.json"))
    rows = {(r["option"], r["setting"]): r for r in csv.DictReader(open(O / "summary.csv"))}
    best = json.load(open(O / "best.json"))
    conv = any("convergence" in r and r["convergence"] for r in rows.values())
    cruns = [r for r in runs["runs"] if r["option"] == "C"]
    if cruns:
        chk = []
        for r in cruns:
            g = Path(r["root"]) / "gpu_eval" / "gpu_check.json"
            chk.append(json.load(open(g)) if g.exists() else None)
        if all(c is not None and c.get("passed") for c in chk):
            OPT["C"] = OPT["C"].replace(
                "GPU cross-check pending",
                "the GPU forward model reproduces them (one GPU evaluation per setting at the same parameters: data cost "
                "relative difference <= %.1e, no printed difference in COPS_model.dat; gpu_eval/gpu_check.json)" %
                max(c["cost_rel_diff"] for c in chk))
        elif any(c is not None and c.get("passed") for c in chk):
            ok = [(r, c) for r, c in zip(cruns, chk) if c is not None and c.get("passed")]
            OPT["C"] = OPT["C"].replace(
                "GPU cross-check pending",
                "the GPU forward model reproduces them for %s (one GPU evaluation per setting at the same parameters: data "
                "cost relative difference <= %.1e, no printed difference in COPS_model.dat; gpu_eval/gpu_check.json), not "
                "checked on the GPU for %s" % (", ".join(r["setting"] for r, _ in ok), max(c["cost_rel_diff"] for _, c in ok),
                                                ", ".join(r["setting"] for r, c in zip(cruns, chk)
                                                          if c is None or not c.get("passed"))))

    def desc(opt, setting, label=None):
        r = rows[(opt, setting)]
        f = lambda k: ("%.4f" % float(r[k])) if r.get(k) not in (None, "", "-") else "-"   # noqa: E731
        return ("%s: data cost %.4e%s, rms ln NG %s (before the gap %s, late %s), RG %s (%s, %s), %s" % (
            label or ("%s %s" % (opt, setting)), float(r["data_cost"]),
            ((", penalty %.3e" % float(r["penalty"])) + ((" incl. the multiplier penalty P_m %.3e" % float(r["penalty_mult"]))
                                                         if r.get("penalty_mult") not in (None, "") else ""))
            if r.get("penalty") not in ("", "0", "0.0") else "",
            f("rms_NG_all"), f("rms_NG_before_gap"), f("rms_NG_late"), f("rms_RG_all"), f("rms_RG_before_gap"),
            f("rms_RG_late"), (r.get("termination") or "-").replace("The maximum number of function evaluations is exceeded.",
                                                                     "stopped at max_nfev").replace(
                "`ftol` termination condition is satisfied.", "stopped on ftol").rstrip(".")) +
                ((", convergence: %s (%s)" % (r["convergence"], r["conv_note"])) if r.get("convergence") else "") +
                ((", continued from the stage-1 fit with data cost %.4e%s (%s evaluations, %s)" % (
                    float(r["prev_data_cost"]), (", penalty %.3e" % float(r["prev_penalty"]))
                    if r.get("prev_penalty") not in ("", "0", "0.0") else "", r["prev_nfev"],
                    r["prev_termination"].replace("The maximum number of function evaluations is exceeded.", "stopped at max_nfev")
                    .replace("`ftol` termination condition is satisfied.", "stopped on ftol").rstrip(".")))
                 if r.get("prev_data_cost") not in (None, "") else ""))
    data = ("Black (NG) and grey (RG) dots: the COPS observations used in the fits (934 NG, 1150 RG; NG excluded from "
            "10:24:40 UTC, the NG data gap and the late window shaded as in the other model / data figures; the 95 "
            "late-window observations compared with the gauges' 10-s / 5-readout response averages W10x5 of the 1-s grid "
            "model, the NG model including the fitted grain term). Lines: the fitted models (observation rows of each "
            "fit's COPS_model.dat).")
    written = 0
    for opt in sorted({r["option"] for r in runs["runs"]}):
        fits = [r for r in runs["runs"] if r["option"] == opt]
        for sfx, what in SUFFIX.items():
            p = O / ("12M_options_%s%s_caption.txt" % (opt, sfx))
            if not (O / ("12M_options_%s%s.png" % (opt, sfx))).exists():
                continue
            if p.exists():
                raise SystemExit("%s exists" % p)
            cap = ("Modelled COPS pressures of %s, its stage-1 settings completed so far (%s), with the bounded softplus c9 and the "
                   "COPS data, %s. %s %s Fits: %s; %s. %s%s" % (
                       OPT[opt], ", ".join(r["setting"] for r in fits), what, COMMON, data,
                       "; ".join(desc(opt, r["setting"]) for r in fits), desc("ref", "9 cycles", "softplus c9"),
                       ("Option D f 1 and f 10 stopped on ftol after 55 and 19 evaluations (under diagnosis). "
                        if opt == "D" else "") if not conv else CONV + ((d_note + " ") if (d_note and opt == "D") else ""),
                       "Table: summary.txt."))
            p.write_text(cap + "\n")
            written += 1
    for sfx, what in SUFFIX.items():
        p = O / ("12M_options_best%s_caption.txt" % sfx)
        if not (O / ("12M_options_best%s.png" % sfx)).exists():
            continue
        if p.exists():
            raise SystemExit("%s exists" % p)
        cap = ("Modelled COPS pressures of the best stage-1 setting (lowest data cost) of each non-negative activity option, "
               "the bounded softplus c9 and the clipped refit D2 (the start of the stage-1 fits), with the COPS data, %s. "
               "Options: %s. %s %s Fits: %s; %s; %s. %sTable: summary.txt." % (
                   what, "; ".join(OPT[o] for o in sorted(best)), COMMON, data,
                   "; ".join(desc(o, b["setting"], "%s %s" % (o, b["setting"])) for o, b in sorted(best.items())),
                   desc("ref", "9 cycles", "softplus c9"), desc("ref", "D2", "clipped D2"), CONV if conv else ""))
        p.write_text(cap + "\n")
        written += 1
    print("captions written: %d" % written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
