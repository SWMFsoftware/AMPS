#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-23, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/full_set_allheld_spec.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/full_set_allheld_spec.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 95786cef603f15ee057fa1f4acbf299a1c9759c4851cee48e4f47784c2bf56f7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: writes the dated specs of the all-held candidate's source-rate and full sets (2026-11-23)
# RUN: python3 full_set_allheld_spec.py (see its header)
# ==== end of header ====
"""Specs of the full figure set of the all-held candidate model (FINDINGS A23; user 2026-11-23: "a full figure set for
the all-held candidate"): option C, f 10, templates off, all six slope multipliers held at 1, new H2O law, the chain of
follow-up c (slopetest_c_q002_..._20261122T195607Z/chain, final state c5/ft/final). New dated files, the 2026-11-18 /
2026-11-20 specs are not edited:

  data/source_rate_q002_runs_20261123.json   = the 2026-11-18 runs spec + run C_f10_allheld_q002 + model text
                                               localnn_allheld (for source_rate_q002_submit.py --runs)
  data/full_sets_q002_20261123.json          = one set q002_C_f10_allheld built like full_set_q002.py build_spec's
                                               option-C entries (for full_set_q002.py --spec)

    full_set_allheld_spec.py --post-job <the chain post job>     (refuses existing files; needs the chain's maps_ft)
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import full_set_q002 as FS

S = FS.S
RUNS_OLD = S / "data/source_rate_q002_runs_20261118.json"
RUNS_NEW = S / "data/source_rate_q002_runs_20261123.json"
SPEC_NEW = S / "data/full_sets_q002_20261123.json"
KEY, SRK = "q002_C_f10_allheld", "C_f10_allheld_q002"
CHAIN = FS.NN / "localnn/slopetest_c_q002_f10_pm0.01_tploff_20261122T195607Z/chain"
ADAPTER = FS.FZ / "localnn_slopetest_20261122T033334Z/localnn_tools.py"
TEXT = ("Its H2O activity factor is the local non-negative basis (option C, templates off) with every slope multiplier "
        "held at 1, A = (sum_j c_j phi_j) exp(h): phi_j piecewise-linear hat functions of a coarse triangulation of the "
        "surface with coefficients c_j >= 0 (smoothness penalty f 10; no slope dependence and no multiplier penalty), h the "
        "accumulated per-facet fine-tuning term. The emitted H2O flux is F_H2O = q_H2O A; the code gives units only to F "
        "(m^-2 s^-1), not separately to q and A.")
NAME_SR = ("the local non-negative basis model with all slope multipliers held at 1 (option C, f 10, templates off; the "
           "candidate of FINDINGS A23: the same solution from seven starts), new H2O law (q_min 0.02 q_max), chain of 5 fit "
           "and fine-tuning cycles, final state")
TITLE = "the option-C local non-negative basis with all slope multipliers held at 1 (candidate; f 10, templates off)"
NAME = ("option C, f 10, templates off, all slope multipliers held at 1 (candidate, FINDINGS A23), chain of 5 fit and "
        "fine-tuning cycles, final state")
LABEL = "local basis C, all multipliers 1 (candidate)"
CAVEATS = ("last fit gtol (1)|fits of cycles 1-5: ftol 68, ftol 3, then gtol after 1 evaluation each; data cost 7.62008030e-3 "
           "after cycle 1, 7.62008004e-3 after cycles 2-5 (3.5e-8 relative)|all slope multipliers held at 1 (no slope "
           "dependence by construction; P_m = 0)|kappa at its lower bound")
INTRO = ("User 2026-11-23: \"a full figure set for the all-held candidate\": a FULL figure set of the candidate of FINDINGS "
         "A23 (option C, f 10, templates off, all six slope multipliers held at 1; the same solution from seven starts and "
         "through five fit and fine-tuning cycles), new H2O law q_min = 0.02 q_max, exactly like the ten new-law full sets "
         "of 2026-11-20 (../full_set_q002_*/; same tools, settings, styles and user decisions), no movies. activity/ and "
         "production_rate/ are copies (cp -p, sha256 in COPIED_FROM_source_rate_q002.txt) of ../source_rate_q002/"
         + SRK + "/; there is no old-law counterpart, so no inspection/ (as C f10 / P_m 1e-4).")

PUB_REASON = 'user 2026-11-23: "a full figure set for the all-held candidate"'
CMP = ("vs the C f10 / P_m 0.01 / templates-off set with free multipliers (../full_set_q002_C_f10_pm0p01_tploff/; the same "
       "model and law, multipliers fitted)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--post-job", required=True)
    a = ap.parse_args()
    for p in (RUNS_NEW, SPEC_NEW):
        if p.exists():
            raise SystemExit("%s exists" % p)
    c5 = CHAIN / "c5"
    ft, fit = c5 / "tools_state_ft", c5 / "tools_state_ft/fit_state"
    maps, contrib = c5 / "maps_ft", c5 / "contributions_ft"
    xc, phase = c5 / "v2_total_crosscheck_ft", c5 / "contrib_phase_ft"
    m = FS.maps_state(maps)
    if hashlib.sha256(ADAPTER.read_bytes()).hexdigest() != hashlib.sha256(
            (FS.FZ / "localnn_q002_20261118T162327Z/localnn_tools.py").read_bytes()).hexdigest():
        raise SystemExit("the adapter differs from the one of the ten sets")
    runs = json.load(open(RUNS_OLD))
    runs["doc"] += " 2026-11-23: + run %s (the all-held candidate, FINDINGS A23) and model text localnn_allheld." % SRK
    runs["runs"][SRK] = {"model": "localnn_allheld", "maps": str(maps), "old_maps": None, "model_name": NAME_SR,
                         "old_name": None}
    runs["model_texts"]["localnn_allheld"] = TEXT
    RUNS_NEW.write_text(json.dumps(runs, indent=1) + "\n")
    geom = m["geometry"]
    mname = "%s, %s" % (NAME, FS.LAW)
    pub = "publication_%s" % KEY
    env = {"KEY": KEY, "CHAIN": str(c5), "FT": str(ft), "FIT": str(fit), "MAPS": str(maps), "GEOM": geom,
           "CONTRIB": str(contrib), "XC": str(xc), "PHASE": str(phase), "MX": str(FS.MX), "CX": str(FS.CX),
           "ADAPTER": str(ADAPTER), "PUB": str(FS.F12 / pub), "PUBSTEM": "12M_%s_NG_RG_altitude" % KEY,
           "PUBTABLE": "%s_final.dat" % KEY, "MODEL_NAME": mname,
           "MODEL_LABEL": "%s, new law q_min 0.02 q_max, final state" % LABEL, "MODEL_SHORT": "%s, new law" % LABEL,
           "TITLE_PRESSURES": "SHAP8 12M: %s, new H2O law (q_min 0.02 q_max) vs COPS" % LABEL, "MODEL_TEXT": TEXT,
           "BASE": str(FS.BASE), "D2REF": str(FS.D2REF), "MCODE": str(FS.MX / "nonneg_maps.py"), "MEM_CAP": "24G",
           "SKIP": "activity,h2o", "NODE": "gpu002", "NCYC": "5"}
    fin = {"KEY": KEY, "FT": str(ft), "FIT": str(fit), "MAPS": str(maps), "CONTRIB": str(contrib), "PHASE": str(phase),
           "MODEL_NAME": mname, "BASE": str(FS.BASE), "REF_DIR": str(FS.D2REF), "CAVEATS": CAVEATS,
           "CHAIN": str(CHAIN), "NCYC": "5"}
    notes = ["Law: %s." % FS.LAW,
             "State inputs from its maps manifest: maps %s (job %s), geometry %s." % (maps, a.post_job, geom),
             "activity/ and production_rate/ and their checks/ files copied from ../source_rate_q002/%s/ "
             "(COPIED_FROM_source_rate_q002.txt); the job that made them: the copied provenance-<job>/ folder." % SRK,
             "Model-data reference (model_data/, regions): the refitted D2 with the same law (%s)." % FS.D2REF,
             "Caveats: %s." % CAVEATS.replace("|", "; "),
             "The chain's post (report, tools export / contrib / maps; degeneracy not run: all multipliers held): job %s "
             "with the fits' frozen code model_ext_frozen/localnn_slopetest_20261122T033334Z." % a.post_job,
             "'A = 0' / 'A > 0' shares (facts, captions; light green = exactly zero) do not measure switched-off surface, as "
             "in CORRECTION_20261121.md of the ten new-law sets (decision 2b): option C's bounded solver leaves switched-off "
             "activity at small positive values (A_min 2.2e8 here). Decision 1b (active = A > T = 2.114e15, 1 % of the "
             "median srcCG activity of the seen facets): 35.5 % of the seen area switched off (T x 0.1: 33.1 %, x 10: "
             "43.1 %), passes (fails above 50 %); checks/switched_off_1b.{py,txt}.",
             "The late-visibility templates of visibility/ belong to an earlier base with another gauge position (FINDINGS "
             "A20, CORRECTION_20261121.md of the ten sets)."]
    spec = {"_about": "Full figure set of the all-held candidate (user 2026-11-23); made by full_set_allheld_spec.py; "
                      "full_set_q002.py --spec this; full_set_docs.py --spec this --key %s" % KEY,
            "sets": {KEY: {
                "title": "Complete figure set (12M) of %s, new H2O law (q_min = 0.02 q_max), final state" % TITLE,
                "model_name": mname, "state": "%s (tools: %s)" % (c5 / "ft/final", ft), "chain": str(c5), "ft": str(ft),
                "fit": str(fit), "maps": str(maps), "contrib": str(contrib), "phase": str(phase), "xc": str(xc),
                "geom": geom, "srkey": SRK, "jobs_state": ["chain post (maps, contributions) job %s" % a.post_job],
                "frozen": {"maps / contributions": "model_ext_frozen/nonneg_maps_20261118T192126Z, nonneg_contrib_20261118T192126Z",
                           "adapter": str(ADAPTER)},
                "notes": notes, "comparisons": [CMP], "intro": INTRO, "pub_reason": PUB_REASON,
                "reference": "the refitted D2 with the same law (reference)", "pub": pub, "written": "2026-11-23 (CST)",
                "extra_dirs": ["inspection/        not made (no old-law counterpart of this model)"],
                "env": env, "finish_env": fin,
                "submit_command": "full_set_q002.py submit --spec %s --keys %s (bash full_set_submit.sh)" % (SPEC_NEW, KEY),
                "finish_command": "srun -p normal128cores -A general -N1 -n1 -c 8 --mem=96G -t 00:45:00 <full_set_q002.py "
                                  "finish --spec %s --key %s --job-set4 <set1/set4 job>>" % (SPEC_NEW, KEY),
                "state_commands": ["# maps / contributions: the chain's post, job %s (localnn_chain_post.sbatch, "
                                   "POST_STEPS='report tools docs')" % a.post_job]}}}
    SPEC_NEW.write_text(json.dumps(spec, indent=1) + "\n")
    print("wrote", RUNS_NEW, SPEC_NEW)
    return 0


if __name__ == "__main__":
    sys.exit(main())
