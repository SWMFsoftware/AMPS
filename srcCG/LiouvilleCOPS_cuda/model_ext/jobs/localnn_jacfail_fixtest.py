# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/jacfail/jac_fixtest.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/jacfail/jac_fixtest.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5be0cc0c17d1058f50b85325ec48f30f61b3b2b13e0bd27caa43c33b02cdca7f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the fixed jacobian_check at the failing state (7.4e-7) and two other stage-2 starts (1.6e-7, 1.7e-7).
# RUN 2026-11-15: srun 2026-11-16 16:05 UTC (jac_fixtest.out).
# ==== end of header ====
# the fixed jacobian_check at the failing state (f1 / 1e-4 c4 start) and at the stage-1 / stage-2 states of the earlier tests
import sys, json, os, numpy as np
sys.path.insert(0, "/data001/highresolu_vt/cops_output/model_ext")
import localnn_problem as lp, localnn_fit as lf, localnn_late as llate
print(lf.__file__)
R="/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5"; L=R+"/nonneg/localnn"
grid=lp.load_grid(L+"/grid_150m"); design=lp.load_design(L+"/grid_150m/design"); D=lp.Data(R+"/base", R+"/clip_refit/D2/final", R+"/late_grid_points.txt", design)
nf=len(grid["facet_node"]); LK=llate.load(L+"/late_kernels_150m.npz", nf)
f=np.load(design["meta"]["paths"]["facets_from"]); TS=f["TS"]; bins=f["bin"]
cases=[("f1_pm1e-4 c4 start (failed)", L+"/stage2_f1_pm1e-4_20261115T190017Z", 4, 8.249495752351283e-06, 1e-4),
       ("f10_pm0.01 c2 start", L+"/stage2_f10_pm0.01_20261115T190019Z", 2, 8.249495752351283e-05, 0.01),
       ("f1_pm0.01 c3 start", L+"/stage2_f1_pm0.01_20261115T190017Z", 3, 8.249495752351283e-06, 0.01)]
for lab, C, k, lam, fm in cases:
    x0=np.load(C+"/c%d/fit/start/start.npz" % k)["x0"]
    h=np.load(os.path.realpath(C+"/c%d/ft/final" % (k-1))+"/finetune_state.npz")["held_next"]
    P=lp.Problem(D, design, grid, True, 1.0, design["meta"]["A_scale"], lam, mult_penalty_f=fm, held={"h":h,"LK":LK,"bins":bins,"TS":TS})
    x=np.minimum(np.maximum(x0,P.lo),P.hi)
    chk=lf.jacobian_check(P,x)
    w=sorted(chk["per_column"].items(), key=lambda t:-t[1])[:3]
    print(lab, "passed", chk["passed"], "max_rel %.3e" % chk["max_rel_error"], "worst", [(a, "%.2e" % b) for a,b in w], flush=True)
