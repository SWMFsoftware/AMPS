# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/jacfail/jac_debug.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/jacfail/jac_debug.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: be7433e74382a94346fa5942fadacbda6f91ba819c9cbb69772f7efac799efb2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: diagnosis of the failed c4 Jacobian check of stage2_f1_pm1e-4: the worst columns are bound nodes, the one-
#   sided FD error falls as 1/h (rounding noise), the analytic columns are right.
# RUN 2026-11-16: srun 2026-11-16 16:00 UTC.
# ==== end of header ====
import sys, json, numpy as np
sys.path.insert(0, "/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2_20261115T185842Z")
import localnn_problem as lp, localnn_fit as lf, localnn_late as llate
R="/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5"; L=R+"/nonneg/localnn"; C=L+"/stage2_f1_pm1e-4_20261115T190017Z"
grid=lp.load_grid(L+"/grid_150m"); design=lp.load_design(L+"/grid_150m/design"); D=lp.Data(R+"/base", R+"/clip_refit/D2/final", R+"/late_grid_points.txt", design)
x0=np.load(C+"/c4/fit/start/start.npz")["x0"]
nf=len(grid["facet_node"]); LK=llate.load(L+"/late_kernels_150m.npz", nf)
f=np.load(design["meta"]["paths"]["facets_from"]); TS=f["TS"]; bins=f["bin"]
import os
ft=os.path.realpath(C+"/c3/ft/final"); h=np.load(ft+"/finetune_state.npz")["held_next"]
P=lp.Problem(D, design, grid, True, 1.0, design["meta"]["A_scale"], 8.249495752351283e-06, mult_penalty_f=1e-4, held={"h":h,"LK":LK,"bins":bins,"TS":TS})
x=np.minimum(np.maximum(x0,P.lo),P.hi)
chk=lf.jacobian_check(P,x); print("check max", chk["max_rel_error"])
J=P.jac(x); xc=x[P.ic]; typ=lf.node_fd_scale(xc); print("typ", typ, "held nonzero", P.held["nonzero"])
for nm in ("c2374","c0851","c1119","c1393"):
    k=P.names.index(nm); den=np.abs(J[:,k]).max(); out=[]
    f0=P.fun(x)
    for scale in (1.0, 10.0, 30.0, 100.0, 1000.0):
        hh=1e-4*max(abs(x[k]),typ)*scale
        x1,x2=x.copy(),x.copy(); x1[k]+=hh; x2[k]+=2*hh
        n2=(-3*f0+4*P.fun(x1)-P.fun(x2))/(2*hh)
        out.append("%.2e" % (np.abs(n2-J[:,k]).max()/den))
    hh=1e-4*max(abs(x[k]),typ); x1=x.copy(); x1[k]+=hh
    dp = np.abs(P.fun(x1)-f0); i=int(np.argmax(np.abs(J[:,k])))
    print(nm, "max|f| %.3e den %.3e x %.3e" % (np.abs(f0).max(), den, x[k]), "one-sided 2nd order rel err at h x 1, 10, 30, 100, 1000:", out, "| row of max |J|:", i, "rel change of r over h there %.2e" % (dp[i]/max(abs(f0[i]),1e-300)))
