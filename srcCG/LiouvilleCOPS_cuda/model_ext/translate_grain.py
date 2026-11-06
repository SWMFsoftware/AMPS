#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/translate_grain.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/translate_grain.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0ddfb55522d31db33e90880d37b29c9619991f574bbd70d42e37c2f4980e7c89
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: carries a grain fit's restart file to another mesh/base (restart_tool translate + predictions there + the
#   grain term), so the fine tuning can run on that mesh.
# USED FOR: 12M T75 = Tmin = 133 K fit -> 20M (job 121974), -> 20M refined (122011), -> 44M (121998); 12M fits ->
#   the refined 12M mesh (122009).
# RUN: .venv/bin/python <this> --restart <fit dir>/fit_restart.npz --to-run <target run dir>
#        [--target-bundle <npz>] --out <new dir> --devices 0,...,7 --precision mixed
#   (model_ext/grain_translate_finetune.sbatch has the full command).
# ==== end of header ====
"""Carry a grain fit's restart to another mesh and give it predictions there (user request 2026-11-04: "take
fitting restart files from the last completed fit and run fine tuning stage of the 20M surface model"; answer:
the 12M fit with T75 = Tmin = 133 K fixed, translated to 20M).

    .venv/bin/python <this> --restart <fit dir>/fit_restart.npz --to-run <base B> --target-bundle <bundle on B>
        --out <new dir> [--devices 0,...,7] [--precision mixed]

1. tools/restart_tool.py translate (unchanged; no --evaluate: it uses one GPU) -> <out>/translated_no_predictions.npz
   (coefficients, T, kappa unchanged; template amplitudes rescaled by A_scale(A)/A_scale(B); a template present in
   the target bundle takes that bundle's column, recomputed on B).
2. Predictions on B with the fine tuning's own evaluator (cops_cuda.fine_tuning.SplitEvaluator, all devices, the
   same record precision as the fine tuning): surface pressures + the source fit's grain on NG (A, B and t_end of
   the source grain_term.json; the observation times are the same on every base). They are the pressures the
   grain fine tuning recomputes as its baseline, so its replay check compares like with like.
3. <out>/fit_restart.npz (new file) with those predictions and B's input.snapshot.in identity in meta.base
   (translate records only elements.bin and COPS_model.dat; the fine tuning checks all three), and <out>/grain_term.json (copy of the source's; the
   grain fine tuning reads it next to the restart). The selection flags (use_ng, use_rg, saturated) and
   meta.selection are the source fit's (translate takes the target table's, which leave the post-gap NG points
   out), after checking that point ids and observations are identical. The source files are only read.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(APP), str(APP / "tools"), str(HERE)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restart", required=True, type=Path)
    ap.add_argument("--to-run", required=True, type=Path)
    ap.add_argument("--target-bundle", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--precision", choices=("fp64", "mixed"), default="mixed")
    a = ap.parse_args()
    if a.out.exists():
        sys.exit("translate_grain: %s exists" % a.out)
    a.out.mkdir(parents=True)
    tmp = a.out / "translated_no_predictions.npz"
    cmd = [sys.executable, "-W", "ignore", str(APP / "tools/restart_tool.py"), "translate", str(a.restart),
           "--to-run", str(a.to_run), "--save", str(tmp), "--target-bundle", str(a.target_bundle), "--device", "0"]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    from cops_cuda import precision as cuda_precision
    os.environ[cuda_precision.RECORD_PRECISION_ENVIRONMENT] = "fp32" if a.precision == "mixed" else "fp64"
    from cops_cuda import fine_tuning as ftg, restart as rs
    import grain_crossing as gc
    term_src = a.restart.parent / "grain_term.json"
    term = json.loads(term_src.read_text())
    A0, B0, t_end = term["result"]["A_Pa"], term["result"]["B_Pa_per_s"], float(term["t_end_et"])
    r = rs.read(tmp)
    ev = ftg.SplitEvaluator(a.to_run, [int(d) for d in a.devices.split(",")])
    try:
        png, prg = ev.pressures(r)
        ids = ev.point_ids
    finally:
        ev.release()
    if not np.array_equal(np.asarray(ids), np.asarray(r.point_ids)):
        sys.exit("translate_grain: evaluator point order differs from the restart table")
    ets = gc.et_by_point(a.to_run)
    s = np.array([ets[int(p)] for p in ids]) - t_end
    grain = np.where(s >= 0, np.maximum(A0 - B0 * np.maximum(s, 0), 0.0), 0.0)
    new = r.copy()
    from cops_cuda import provenance
    # restart_tool translate takes the selection flags from the target base table (there the post-gap NG points are
    # not selected); the fit's selection (profile all_available: use flags and saturated NG) is the source's, on
    # the same observation points: copy it after checking that the points and observations are identical
    src = rs.read(a.restart)
    for k in ("point_ids", "ref_ng", "ref_rg"):
        if not np.array_equal(np.asarray(src.arrays[k]), np.asarray(new.arrays[k])):
            sys.exit("translate_grain: %s differs between the source restart and the target table" % k)
    selection_from_table = {k: int(np.asarray(new.arrays[k]).sum()) for k in ("use_ng", "use_rg", "saturated")}
    for k in ("use_ng", "use_rg", "saturated"):
        new.arrays[k] = np.array(src.arrays[k], copy=True)
    new.meta["selection"] = src.meta.get("selection")
    # restart_tool translate records base_elements and base_model_table only; the fine tuning also checks the
    # base input: the state is on base B, so record B's input.snapshot.in
    new.meta.setdefault("base", {})["base_input"] = provenance.file_identity(a.to_run / "input.snapshot.in")
    new.arrays["pred_ng"], new.arrays["pred_rg"] = png + grain, prg
    new.meta["history"] = list(new.meta.get("history", [])) + [{
        "action": "predictions on the target base (model_ext/translate_grain.py): fine_tuning.SplitEvaluator "
                  "surface pressures + source grain on NG", "devices": a.devices, "record_precision":
                  "fp32" if a.precision == "mixed" else "fp64", "grain": {"A_Pa": A0, "B_Pa_per_s": B0,
                                                                          "t_end_et": t_end},
        "source_grain_term": str(term_src)}]
    rs.write(a.out / "fit_restart.npz", new)
    shutil.copy2(term_src, a.out / "grain_term.json")
    use_ng = (new.use_ng > 0) & (new.ref_ng > 0)
    use_rg = (new.use_rg > 0) & (new.ref_rg > 0)
    rms = lambda m, d, u: float(np.sqrt(np.mean(np.log(m[u] / d[u]) ** 2)))  # noqa: E731
    summary = {"source_restart": str(a.restart), "to_run": str(a.to_run), "target_bundle": str(a.target_bundle),
               "T": [float(v) for v in new.T], "kappa": float(np.asarray(new.kappa).reshape(-1)[0]),
               "templates": dict(zip(new.template_names, [float(v) for v in new.template_values])),
               "translation": new.meta.get("translation"),
               "rms_ln_residual_on_target": {"NG": rms(new.pred_ng, new.ref_ng, use_ng),
                                             "RG": rms(new.pred_rg, new.ref_rg, use_rg)},
               "grain": {"A_Pa": A0, "B_Pa_per_s": B0},
               "selection": {"from_source_restart": {k: int(np.asarray(new.arrays[k]).sum())
                                                     for k in ("use_ng", "use_rg", "saturated")},
                             "target_table_had": selection_from_table}}
    (a.out / "translate_grain.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["rms_ln_residual_on_target"]), "T", summary["T"], "templates", summary["templates"],
          flush=True)


if __name__ == "__main__":
    main()
