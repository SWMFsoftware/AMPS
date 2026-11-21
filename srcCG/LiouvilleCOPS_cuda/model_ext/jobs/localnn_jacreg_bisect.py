# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261118_jacreg/jac_bisect.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261118_jacreg/jac_bisect.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b47bfe92bcc2ce0981b36ce250a6f5e666453728885dad3b5e7df00c45239604
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the test_localnn check-5 states with the jacobian_check of one localnn snapshot (the bisect of the
#   2026-11-18 regression).
# RUN 2026-11-18: sruns 2026-11-18 08:00 UTC (bisect.out).
# ==== end of header ====
# test_localnn check 5 states (start_templates_off; multipliers 1.2..4, grain A 12, B 0.4) with the jacobian_check of
# one localnn snapshot (argv[1]); the Problem of that snapshot; the NEW-law contrib code is irrelevant (check 5 uses
# only the design route)
import sys, json, numpy as np
snap = sys.argv[1]; sys.path.insert(0, snap)
import localnn_problem as lp, localnn_fit as lf
print("snapshot", snap, lf.__file__)
R="/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5"; L=R+"/nonneg/localnn"
grid=lp.load_grid(L+"/grid_150m"); design=lp.load_design(L+"/grid_150m/design"); D=lp.Data(R+"/base", R+"/clip_refit/D2/final", R+"/late_grid_points.txt", design)
A_s=design["meta"]["A_scale"]
with np.load(L+"/grid_150m/start_templates_off/start.npz", allow_pickle=False) as z:
    xs=np.asarray(z["x0"], float); tmpl=str(z["templates"])=="on"
P=lp.Problem(D, design, grid, tmpl, 1.0, A_s, smooth_lambda=1e-6)
xs=np.minimum(np.maximum(xs, P.lo+1e-9), P.hi)
xt=xs.copy(); xt[P.im]=np.array([1.5,2.0,3.0,1.2,4.0,2.5]); xt[P.iA], xt[P.iB] = 12.0, 0.4
for lab, xx in (("start", xs), ("multipliers", xt)):
    chk=lf.jacobian_check(P, xx); w=sorted(chk["per_column"].items(), key=lambda t:-t[1])[:3]
    print(" ", lab, "max_rel %.3e" % chk["max_rel_error"], "passed", chk["passed"], "worst", [(a, "%.2e" % b) for a, b in w], flush=True)
