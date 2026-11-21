# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/postA/a_repro.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/postA/a_repro.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c06417c187e27fae3c1998c012d583ea98f6ae8fee901d3f5d377c8c7d6317c1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: diagnosis of post 136094: CPU reproduction of A vs the fit's nonneg_facets.npz A per node type and
#   NPY_DISABLE_CPU_FEATURES (node028 EPYC 9554 without it: 38 held facets differ, rel <= 3.3e-16; with it, and on
#   node018 EPYC 7763: bitwise).
# RUN 2026-11-17: sruns 136119-136122.
# ==== end of header ====
# where does the CPU reproduction of A differ from the fit's nonneg_facets.npz A (stage2r f10 / 1e-4 c5, fit on gpu002)?
import sys, json, os, platform, numpy as np
sys.path.insert(0, "/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2jc2_20261117T080730Z")
import localnn_problem as lp
seg = sys.argv[1]
z = np.load(seg + "/localnn_state.npz"); meta = json.loads(str(z["meta_json"])); c, m = z["c"], z["m"]
f = np.load(seg + "/nonneg_facets.npz"); A, TS, bins, held, u_f = f["A"], f["TS"], f["bin"], f["held"], f["u"]
grid = lp.load_grid(meta["grid"]["path"]); u = lp.facet_field(grid, c)
A_re = m[bins.astype(np.intp)] * np.maximum(u + TS, 0.0); hz = held != 0.0
A_noexp = A_re.copy(); A_re[hz] = A_re[hz] * np.exp(held[hz])
d = A_re != A
print("host", platform.node(), "cpu", open("/proc/cpuinfo").read().split("model name")[1].split("\n")[0], "NPY_DISABLE", os.environ.get("NPY_DISABLE_CPU_FEATURES", "-")[:20])
print("u bitwise", bool(np.array_equal(u, u_f)), "A mismatches", int(d.sum()), "of which held facets", int((d & hz).sum()), "held", int(hz.sum()),
      "max rel %.3e" % float(np.max(np.abs(A_re[d] / A[d] - 1))) if d.any() else 0.0)
print("fit device", (meta.get("device") or {}).get("device", "cpu"), "fit host", meta.get("host"))
