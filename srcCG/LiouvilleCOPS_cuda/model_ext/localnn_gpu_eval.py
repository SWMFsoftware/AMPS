#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_gpu_eval.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_gpu_eval.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9cdecba69b8c81e490fc6d68e73ab2f464881f6a35b1c540ff6e060dc914b677
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (option C): GPU standard evaluation of a localnn final state through D2's path (clip_refit.py
#   --variant D2, frozen clip_20261108T234013Z): the activity as a held template (A/A_s, amplitude 1), a = 0, kappa, T held;
#   grain A, B free (D2's grain_crossing reads them from the last two optimizer entries: holding them failed, jobs 125052,
#   125123), --max-nfev 1; compare = the CPU state at the GPU's grain vs the GPU (data cost full precision, tables).
# RUN: jobs/localnn_gpu_eval.sbatch (frozen model_ext_frozen/localnn_post2_20261114T222122Z): 125183 (f10, passed), 125185, 125181.
# 2026-11-15: stage 2: prep / compare --ft-dir (a localnn fine tuning's final state A_final and its grain).
# RUN 2026-11-15: the GPU cross-check of the chains' final states (to be submitted when a GPU is free).
# ==== end of header ====
"""Option C (local non-negative basis): GPU cross-check of a localnn fit segment's FINAL per-facet activity through the
standard fitter path (lead 2026-11-14: "one short GPU cross-check per final state on gpu002, writing the standard
products").

The standard path is D2's: clip_refit.py --variant D2 (frozen model_ext_frozen/clip_20261108T234013Z) around the
unchanged grain_fit_fixed.py -> grain_fit.py -> gpu_fit.py -> fit_cops.py, GPU mixed precision (FP32 records, FP64
arithmetic), W10x5 for the late observations, grain at the time stamp, crossing weight 10. Its activity is the fitter's
max(0, Y a + A_s sum_k c_k T_k). With
    a00 .. a24 held at 0, template localnn_a = A_final / A_s held at amplitude 1, kappa and T held at the segment's
    values (fixed_parameters --fix, exact repr),
the fitter's activity is A_s (A_final / A_s) = A_final (one rounding). The GRAIN cannot be held on this path:
grain_crossing.GrainProblem.grain reads A, B from the last two entries of the optimizer vector and fixed_parameters sets
the grain start at x0[-2] (2026-11-14: jobs 125052 and 125123 failed with grain A, B held, IndexError). So grain A and B
are the two free parameters and the fit runs --max-nfev 1: scipy TRF returns its start, the fitter's own grain start
(A = the observed jump at t_end with the held values, B = 0). The fitter's products (fit_result.txt, COPS_model.dat,
late_window_model.dat, fit_restart.npz, fit_manifest.json, grain_term.*, CLIP_REFIT.txt, ...) are then the standard GPU
evaluation of the localnn surface model with that grain. compare evaluates the CPU state with the GPU's grain values
(read from its fit_restart.npz) and compares the data cost (full precision) and the tables; the surface columns do not
depend on the grain.

Subcommands
  prep     --seg SEG --base BASE --out DIR: DIR/localnn_activity.templates.npz (activity_templates.write_bundle;
           localnn_a), DIR/fix_args.txt (the --fix arguments, one per line), DIR/prep.json
  compare  --seg SEG --gpu GPUFIT --design DESIGN --base BASE --d2 D2 --point-table PT --out FILE.json: the CPU state
           with the GPU's grain vs the GPU: data cost (full precision, rel < 1e-10), pNG / pRG of COPS_model.dat (all rows,
           printed digits), late_window_model.dat (surface and grain, 10 digits); exit 1 if a check fails
"""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import localnn_problem as lp      # noqa: E402

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
TOOLS = APP.parent / "LiouvilleCOPS" / "tools"


def sha(path):
    return lp.sha256_file(path)


def cmd_prep(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    seg = Path(args.seg).resolve()
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import activity_templates as at
    import localnn_mesh as lm
    el = lm.read_elements(Path(args.base) / "elements.bin")
    with np.load(seg / "nonneg_facets.npz", allow_pickle=False) as z:
        A = np.asarray(z["A"], np.float64)
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        kappa, gA, gB = float(z["kappa"]), float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    if args.ft_dir:                          # a chain's fine-tuned state: A_final of the fine tuning of this fit
        A, gA, gB = ft_state(args.ft_dir, seg)
    A_s = float(meta["A_scale"])
    col = A / A_s
    back = A_s * col
    rt = {"max_rel_A_s_roundtrip": float(np.max(np.abs(back[A > 0] / A[A > 0] - 1.0))) if np.any(A > 0) else 0.0,
          "zero_facets_preserved": bool(np.all(back[A == 0] == 0.0))}
    bundle = out / "localnn_activity.templates.npz"
    at.write_bundle(bundle, ["localnn_a"], col[:, None], el["centroid"], el["area"],
                    {"templates": {"localnn_a": {"definition": "A_final / A_s of %s/nonneg_facets.npz (A_s %.17g); amplitude "
                                                                "1 (held) gives the localnn activity" % (seg, A_s),
                                                  "interpretation": "fixed per-facet activity for a GPU evaluation"}},
                     "source_segment": str(seg), "nonneg_facets_sha256": sha(seg / "nonneg_facets.npz")})
    fixed = ["a%02d=0" % j for j in range(25)] + ["c_localnn_a=1", "kappa=%r" % kappa, "Tmax=182.09922", "T75=136.11149",
                                                  "Tmin=133.0"]
    (out / "fix_args.txt").write_text("\n".join(fixed) + "\n")
    rep = {"segment": str(seg), "A_scale": A_s, "kappa": kappa, "grain_A_uPa": gA, "grain_B_uPa_per_100s": gB,
           "roundtrip": rt, "grain_free": "grain A, B are the free parameters (start = the fitter's own)", "bundle": {"path": str(bundle), "sha256": sha(bundle)}, "fix": fixed,
           "segment_costs": meta["costs"]}
    with open(out / "prep.json", "x") as f:
        json.dump(rep, f, indent=1, default=float)
    lp.say("prep: %s" % json.dumps(rep["roundtrip"]))
    return 0


def read_cost(path):
    for line in open(path):
        t = line.split()
        if len(t) >= 2 and t[0] == "final_data_cost":
            return float(t[1])
    raise SystemExit("no final_data_cost in %s" % path)


def ft_state(ft_dir, seg):
    """(A_final, grain A, grain B) of a localnn fine tuning of the fit segment seg."""
    ft = Path(ft_dir).resolve()
    j = json.load(open(ft / "fine_tuning.json"))
    if Path(j["fit"]["directory"]).resolve() != Path(seg).resolve():
        raise SystemExit("the fine tuning %s tuned %s, not %s" % (ft, j["fit"]["directory"], seg))
    with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
        return np.asarray(z["A_final"], np.float64), float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])


def state_rows_surface(design, seg, ft_dir=None):
    """CPU instantaneous surface pressures of every row of a fit segment (held exact) or of its fine-tuned state:
    p_fit + G (A_final - A_fit) over the late facets (localnn_ft's linear model)."""
    sN, sR, meta = lp.state_surface(design, seg)
    if ft_dir:
        import localnn_late as llate
        A_fin, _, _ = ft_state(ft_dir, seg)
        with np.load(Path(seg) / "nonneg_facets.npz", allow_pickle=False) as z:
            A_fit = np.asarray(z["A"], np.float64)
        LK = llate.load(json.load(open(Path(ft_dir) / "localnn_ft.json"))["localnn_ft"]["late_kernels"]["path"],
                        len(A_fit))
        li = LK["late_idx"]
        if np.any(A_fin[~LK["late"]] != A_fit[~LK["late"]]):
            raise SystemExit("the fine tuning changed facets outside the late kernels")
        dA = A_fin[li] - A_fit[li]
        sN, sR = sN + LK["G_NG"] @ dA, sR + LK["G_RG"] @ dA
    return sN, sR, meta


def cmd_compare(args):
    seg, gpu = Path(args.seg).resolve(), Path(args.gpu).resolve()
    with np.load(seg / "localnn_state.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        c, m, kappa = np.asarray(z["c"]), np.asarray(z["m"]), float(z["kappa"])
        gA_fit, gB_fit = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    with np.load(gpu / "fit_restart.npz", allow_pickle=False) as z:
        gm_r = json.loads(str(z["meta_json"]))
        xg = np.asarray(z["x"], float)
    names = list(gm_r["names"])
    gA, gB = float(xg[names.index("grain_A_uPa")]), float(xg[names.index("grain_B_uPa_per_100s")])
    gm = json.load(open(gpu / "fit_manifest.json"))
    opt = [o for o in gm["optimizer_results"]][0]
    design = lp.load_design(args.design)
    D = lp.Data(args.base, args.d2, args.point_table, design)
    sN, sR, _ = state_rows_surface(design, seg, args.ft_dir)
    rows = D.model_rows(sN, sR, gA, gB)
    cst = D.costs(rows["pNG"], rows["pRG"])
    cm_g = lp.read_table_cols(gpu / "COPS_model.dat", ("pNG", "pRG"))
    lw_g = lp.read_table_cols(gpu / "late_window_model.dat", ("surface_pNG_Pa", "surface_pRG_Pa", "grain_NG_Pa"))
    res = {"segment": str(seg), "ft_dir": args.ft_dir, "gpu_fit": str(gpu), "gpu_names": names, "gpu_x": [float(v) for v in xg],
           "grain_gpu": {"A_uPa": gA, "B_uPa_per_100s": gB}, "grain_fit": {"A_uPa": gA_fit, "B_uPa_per_100s": gB_fit},
           "cost_gpu_manifest_data_cost": float(opt["data_cost"]), "cost_cpu_with_gpu_grain": cst["data_cost"],
           "cost_rel_diff": abs(cst["data_cost"] / float(opt["data_cost"]) - 1.0),
           "cost_cpu_with_fit_grain_segment": meta["costs"]["data_cost"],
           "COPS_model_pNG": lp.compare(rows["pNG"], {k: v["pNG"] for k, v in cm_g.items()}, D.ipoint, "%.6e"),
           "COPS_model_pRG": lp.compare(rows["pRG"], {k: v["pRG"] for k, v in cm_g.items()}, D.ipoint, "%.6e"),
           "late_surface_pNG": lp.compare(rows["surface_NG"], {k: v["surface_pNG_Pa"] for k, v in lw_g.items()}, D.ipoint, "%.9e"),
           "late_surface_pRG": lp.compare(rows["surface_RG"], {k: v["surface_pRG_Pa"] for k, v in lw_g.items()}, D.ipoint, "%.9e"),
           "late_grain": lp.compare(rows["grain"], {k: v["grain_NG_Pa"] for k, v in lw_g.items()}, D.ipoint, "%.9e"),
           "gpu_nfev": opt.get("nfev"), "gpu_message": opt.get("message")}
    res["passed"] = bool(res["cost_rel_diff"] < 1e-10 and
                         all(res[k]["max_rel"] <= 5.0e-7 * (1 + 1e-6) for k in ("COPS_model_pNG", "COPS_model_pRG")) and
                         all(res[k]["max_rel"] <= 5.0e-10 * (1 + 1e-6) + 1e-14 for k in ("late_surface_pNG", "late_surface_pRG",
                                                                                          "late_grain")))
    with open(args.out, "x") as f:
        json.dump(res, f, indent=1, default=float)
    lp.say("GPU vs CPU: %s" % json.dumps({k: res[k] for k in res if k not in ("gpu_x", "gpu_names")}))
    return 0 if res["passed"] else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    p.add_argument("--seg", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--ft-dir", help="a localnn fine tuning of SEG: evaluate its final state (A_final, its grain)")
    c = sub.add_parser("compare")
    c.add_argument("--seg", required=True)
    c.add_argument("--gpu", required=True)
    for k in ("design", "base", "d2", "point-table"):
        c.add_argument("--" + k, required=True)
    c.add_argument("--ft-dir")
    c.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    return cmd_prep(args) if args.cmd == "prep" else cmd_compare(args)


if __name__ == "__main__":
    sys.exit(main())
