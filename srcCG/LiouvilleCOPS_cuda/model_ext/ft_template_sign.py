# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/checks/ft_template_sign_20261114/ft_template_sign.py (on Voyager: /data001/highresolu_vt/cops_output/checks/ft_template_sign_20261114/ft_template_sign.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4c054f0ad9844c702e86fa7f3ffce87a19855fbd388fc81921ba563a897ccb1a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: check. Sign of D2's ft_c1..ft_c6 templates (lead decision 1 of 2026-11-14): min / max of A_s sum T_k, facet
#   counts, where (seen, latitude / longitude). Result: >= 0 everywhere (265 facets > 0).
# USED BY: srun on normal128cores, job 125029 (checks/ft_template_sign_20261114/).
# ==== end of header ====
"""Sign check of D2's ft_c1..ft_c6 templates (held at 1 in the non-negative models; lead 2026-11-14 decision 1):
min/max of A_s sum_k T_k, facet counts > 0 / < 0 / = 0, where (slope bins, seen, latitude/longitude, the 265-facet
supports). CPU, read-only."""
import json, sys, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
import nonneg_problem as npb
from cops_cuda import reference
reference.bind()
fc = reference.load("fit_cops"); at = reference.load("activity_templates")
R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
el = fc.read_elements(str(R / "base/elements.bin"))
area = np.asarray(el["area"]); cen = np.asarray(el["centroid"])
D2 = R / "clip_refit/D2/final"
fts, src = npb.ft_templates_from(str(D2 / "fit_manifest.json"))
meta = json.loads(str(np.load(str(D2 / "fit_restart.npz"))["meta_json"]))
A_s = float(meta["enhancement_scale"])
cols = [at.select_templates(at.load_bundle(t["bundle"]["path"], el["centroid"], area), [t["name"]]) for t in fts]
T = np.column_stack(cols)
S = A_s * T.sum(axis=1)
bins = np.load(str(R / "slope_mult/bundle/slope_bins.npz"))
out = {"A_s": A_s, "templates": [t["name"] for t in fts], "source": src}
for k, t in enumerate(fts):
    c = T[:, k]
    out[t["name"]] = {"min": float(c.min()), "max": float(c.max()), "n_pos": int((c > 0).sum()), "n_neg": int((c < 0).sum())}
sup = [set(np.flatnonzero(T[:, k] > 0).tolist()) for k in range(T.shape[1])]
out["same_support_all_six"] = all(s == sup[0] for s in sup)
pos = S > 0
out["sum"] = {"min": float(S.min()), "max": float(S.max()), "n_pos": int(pos.sum()), "n_neg": int((S < 0).sum()),
              "n_zero": int((S == 0).sum()), "area_pos_m2": float(area[pos].sum()),
              "max_over_A_s": float(S.max() / A_s), "mean_pos_over_A_s": float(S[pos].mean() / A_s)}
r = np.linalg.norm(cen[pos], axis=1)
lat = np.degrees(np.arcsin(cen[pos, 2] / r)); lon = np.degrees(np.arctan2(cen[pos, 1], cen[pos, 0]))
out["where"] = {"lat_deg": [float(lat.min()), float(lat.max())], "lon_deg": [float(lon.min()), float(lon.max())],
                "centroid_mean_m": [float(v) for v in cen[pos].mean(axis=0)],
                "max_distance_from_mean_m": float(np.max(np.linalg.norm(cen[pos] - cen[pos].mean(axis=0), axis=1)))}
print("bins files", bins.files)
key = [k for k in bins.files if k in ("bin", "bin_index", "bins")]
if key:
    b = np.asarray(bins[key[0]])
    out["where"]["slope_bin_counts"] = {str(i): int(((b == i) & pos).sum()) for i in range(7)}
z = np.load(str(D2 / "nonneg_facets.npz")) if (D2 / "nonneg_facets.npz").exists() else None
sm = None
for cand in sorted((R / "nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c1/fit/final").glob("nonneg_facets.npz")):
    sm = np.load(str(cand))
if sm is not None:
    seen = np.asarray(sm["seen"]).astype(bool)
    out["where"]["seen_of_support"] = int((seen & pos).sum())
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[2], "x"), indent=1)
