# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261118_jacreg/jac_old_check_current_problem.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261118_jacreg/jac_old_check_current_problem.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 19750b2e7aa45b1b02784ce396b7e0906022225f66926ac123110a7579ba396d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the 2026-11-14 jacobian_check applied to the current Problem at the check-5 states (the analytic Jacobian is
#   unchanged: 6.8e-6 / 1.1e-6).
# RUN 2026-11-18: srun 2026-11-18 08:10 UTC (multistep_test.out).
# ==== end of header ====
# the 2026-11-14 jacobian_check (step 1e-4 max(|x|, 1) for the nodes) applied to the CURRENT localnn_problem.Problem
# (run copy) at the test_localnn check-5 states: the analytic Jacobian is unchanged / correct
import sys, types, numpy as np, ast, inspect
sys.path.insert(0, "/data001/highresolu_vt/cops_output/model_ext")
import localnn_problem as lp
src = open("/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_20261114T081912Z/localnn_fit.py").read()
tree = ast.parse(src); fn = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "jacobian_check"][0]
ns = {"np": np}; exec(compile(ast.Module(body=[fn], type_ignores=[]), "jc1114", "exec"), ns); old_check = ns["jacobian_check"]
R="/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5"; L=R+"/nonneg/localnn"
grid=lp.load_grid(L+"/grid_150m"); design=lp.load_design(L+"/grid_150m/design"); D=lp.Data(R+"/base", R+"/clip_refit/D2/final", R+"/late_grid_points.txt", design)
with np.load(L+"/grid_150m/start_templates_off/start.npz", allow_pickle=False) as z:
    xs=np.asarray(z["x0"], float); tmpl=str(z["templates"])=="on"
P=lp.Problem(D, design, grid, tmpl, 1.0, design["meta"]["A_scale"], smooth_lambda=1e-6)
xs=np.minimum(np.maximum(xs, P.lo+1e-9), P.hi)
xt=xs.copy(); xt[P.im]=np.array([1.5,2.0,3.0,1.2,4.0,2.5]); xt[P.iA], xt[P.iB] = 12.0, 0.4
for lab, xx in (("start", xs), ("multipliers", xt)):
    chk=old_check(P, xx); print("11-14 check on the current Problem:", lab, "max_rel %.3e passed %s" % (chk["max_rel_error"], chk["passed"]), flush=True)
