#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_anchor.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_anchor.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f7425c6c70ea4639d030897bdbb9138128a60900bc5aa6088652265442e85429
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (option C): the solution-anchored smoothness rule (user 2026-11-14 ~02:55 CST): R_pilot = x_c^T Lg x_c of
#   the 200-evaluation pilot state, C_ref = the bounded softplus c9 data cost, lambda_C = f C_ref / R_pilot (f 0.1, 1, 10)
#   -> the rule JSON given to localnn_fit.py --smooth-rule-file.
# RUN: srun 2026-11-14 09:42 UTC (job 125106) on the pilot 125071 -> runs/.../nonneg/localnn/lambda_rule_anchored_20261114.json.
# ==== end of header ====
"""Option C (local non-negative basis): the SOLUTION-ANCHORED smoothness rule (user decision via the lead, 2026-11-14
~02:55 CST):  lambda_C = f * C_ref / R_pilot,  f = 0.1, 1, 10,
  C_ref   = the data cost of the bounded softplus c9 final state (finetune_state.npz costs_final.cost_data;
            4.268806657481027e-03 = the lead's 4.2688066575e-03),
  R_pilot = x_c^T Lg x_c (x_c = c / A_s, Lg = the coarse grid's P1 stiffness; the normalization of the penalty
            lambda_C x_c^T Lg x_c) of the PILOT state: the templates-on, start-rule f = 0.1 fit (lambda_C =
            3.7558561837e-8) after 200 evaluations (the stage-1 job 125014 saved no intermediate state, so the pilot is
            the same command with --max-nfev 200, same snapshot and CPU type).
Checks: R_pilot recomputed from the state's x equals penalty / lambda of the pilot's manifest (rel < 1e-12); the pilot
ran templates on and the start rule's f = 0.1.

CLI: localnn_anchor.py --pilot-seg <pilot run root>/fit/final --grid GRID --c9-ft <c9 ft dir> --out rule.json
Output: rule.json {"rule", "description", "C_ref", "R_pilot", "pilot": {...}, "factors", "lambda_C": {"0.1": .., "1": ..,
"10": ..}, ...} -- given to localnn_fit.py --smooth-rule-file (jobs/localnn_fit.sbatch SMOOTH_RULE_FILE).
"""

import argparse
import json
import socket
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402

FACTORS = (0.1, 1.0, 10.0)
RULE = "lambda_C = f * C_ref / R_pilot (solution-anchored; user 2026-11-14 ~02:55 CST)"
C_REF_LEAD = 4.2688066575e-03


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("pilot-seg", "grid", "c9-ft", "out"):
        ap.add_argument("--" + k, required=True)
    args = ap.parse_args(argv)
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    seg = Path(args.pilot_seg).resolve()
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        x = np.asarray(z["x"], float)
        meta = json.loads(str(z["meta_json"]))
    grid = lp.load_grid(args.grid)
    if lp.sha256_file(Path(args.grid) / "grid.npz") != meta["grid"]["sha256"]:
        raise SystemExit("the pilot ran on another grid")
    n = grid["n"]
    xc = x[:n]
    R = float(xc @ (grid["Lg"] @ xc))
    lam_p = float(meta["lambda_C"])
    pen = float(meta["costs"]["penalty"])
    chk = abs(pen / lam_p / R - 1.0)
    if chk > 1e-12:
        raise SystemExit("R_pilot %r vs penalty / lambda %r" % (R, pen / lam_p))
    if meta["templates"] != "on" or meta["lambda_source"].get("factor") != 0.1 or "rule_file" in meta["lambda_source"]:
        raise SystemExit("the pilot is not the templates-on start-rule f = 0.1 fit")
    with np.load(Path(args.c9_ft) / "finetune_state.npz", allow_pickle=False) as z:
        C_ref = float(json.loads(str(z["meta_json"]))["costs_final"]["cost_data"])
    if abs(C_ref / C_REF_LEAD - 1.0) > 1.2e-11:
        raise SystemExit("C_ref %r differs from the lead's %r" % (C_ref, C_REF_LEAD))
    lam = {("%g" % f): f * C_ref / R for f in FACTORS}
    rec = {"rule": RULE, "created_utc": lp.utc_now(), "host": socket.gethostname(), "command": sys.argv,
           "description": "lambda_C = f C_ref / R_pilot, C_ref = %.12e (bounded softplus c9 data cost), R_pilot = %.12e "
                          "(x_c^T Lg x_c of the 200-evaluation pilot)" % (C_ref, R),
           "C_ref": C_ref, "C_ref_lead": C_REF_LEAD,
           "C_ref_source": {"path": str((Path(args.c9_ft) / "finetune_state.npz").resolve()),
                            "sha256": lp.sha256_file(Path(args.c9_ft) / "finetune_state.npz"),
                            "key": "meta_json costs_final.cost_data"},
           "R_pilot": R, "R_definition": "x_c^T Lg x_c, x_c = c / A_s (A_s %.17g), Lg = grid.npz Lg" % meta["A_scale"],
           "R_check_penalty_over_lambda_rel": chk,
           "pilot": {"segment": str(seg), "runroot": str(seg.parent.parent), "nfev": meta["nfev"],
                     "nfev_total": meta["nfev_total"], "stop": meta["stop"], "message": meta["message"],
                     "lambda_C": lam_p, "lambda_source": meta["lambda_source"], "templates": meta["templates"],
                     "data_cost": meta["costs"]["data_cost"], "penalty": pen, "objective": meta["costs"]["objective"],
                     "crossing_cost": meta["costs"]["crossing_cost"], "rms_ln_NG": meta["costs"]["rms_ln_NG"],
                     "rms_ln_RG": meta["costs"]["rms_ln_RG"], "slurm_job_id": meta["slurm_job_id"],
                     "localnn_state_sha256": lp.sha256_file(seg / "localnn_state.npz")},
           "grid": {"path": str(Path(args.grid).resolve()), "sha256": meta["grid"]["sha256"], "nodes": n},
           "factors": list(FACTORS), "lambda_C": lam,
           "penalty_of_the_pilot_state_with_these_lambda": {k: v * R for k, v in lam.items()},
           "code_sha256": {"localnn_anchor.py": lp.sha256_file(__file__)}}
    with open(out, "x") as f:
        json.dump(rec, f, indent=1, default=float)
    lp.say("R_pilot %.12e (pilot nfev %d, data cost %.10e); lambda_C %s" % (R, meta["nfev"], meta["costs"]["data_cost"],
                                                                          json.dumps(lam)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
