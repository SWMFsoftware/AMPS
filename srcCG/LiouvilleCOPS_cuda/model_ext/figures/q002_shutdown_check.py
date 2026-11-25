# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/shutdown_check.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/slope_activity/shutdown_check/shutdown_check.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2730c52b8c935bd78d63c6890f62ead2433c0a97d4fe1ed4dc71c9f9e806a6f8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: switched-off surface of the ten fits (A < k A_u); slope means without template facets (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""2026-11-21 check of two numbers before they are reported:
1. how much surface each of the ten new-law fits switches off, without relying on A == 0 exactly (bounded solvers leave
   switched-off coefficients slightly positive): area fraction with A below k x A_u, k = 1e-6, 0.01, 0.1, with A_u the
   uniform activity that gives the fit's own seen production at 10:33:10 (A_u = Q_seen / sum_seen q a); also against one
   common A_u (median over the ten), and the smallest area holding 50 / 90 / 99 % of the seen production;
2. slope statistics of the five fits without the held-template and late fine-tuned facets.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
import q002_slope_activity as Q
met = json.load(open(Q.PL))["runs"]
nf, cen, a = Q.elements(); del cen
with np.load(Path(met[Q.A_RUN]["fit_dir"]) / "nonneg_facets.npz") as z:
    seen = z["seen"].astype(bool); bin_ = z["bin"].astype(int)
    special = (z["TS"] != 0) | (z["held"] != 0)
T = Q.map_times(); t = [k for k, v in T.items() if v["utc"].endswith("10:33:10")][0]
q, qmax, cos, sh = Q.qfac(t, T[t]["rh_m"], nf)
top = [k for k, _, _ in Q.RUNS]
w = (q * a)[seen]; atot, asn = a.sum(), a[seen].sum()
d = {"special_facets": int(special.sum()), "special_area_m2": float(a[special].sum()), "fits": {}}
Au, keep = {}, {}
for k in met:
    with np.load(Path(met[k]["ft_dir"]) / "finetune_state.npz") as z:
        A = z["A_final"].astype(float)
    As = A[seen]; a_s = a[seen]
    u = float((As * w).sum() / w.sum()); Au[k] = u; keep[k] = A.astype(np.float32)
    r = {"A_u": u, "median_A_seen_over_A_u": Q.wquant(As, a_s, [0.5])[0] / u}
    for n, f in (("1e-6", 1e-6), ("0.01", 0.01), ("0.1", 0.1)):
        r["off_seen_" + n] = float(a_s[As < f * u].sum() / asn)
        r["off_whole_" + n] = float(a[A < f * u].sum() / atot)
    o = np.argsort(-As); c = np.cumsum((As * w)[o]); c /= c[-1]; ca = np.cumsum(a_s[o]) / asn
    for p in (0.5, 0.9, 0.99):
        r["seen_area_for_%g_of_production" % p] = float(ca[min(np.searchsorted(c, p), len(ca) - 1)])
    if k in top:
        m = seen & ~special; mean = (A[m] * a[m]).sum() / a[m].sum()
        r["A_rel_no_special"] = [float((A[m & (bin_ == b)] * a[m & (bin_ == b)]).sum() / a[m & (bin_ == b)].sum() / mean)
                                 for b in range(7)]
        r["special_share_of_seen_activity"] = float((A[seen & special] * a[seen & special]).sum() / (As * a_s).sum())
    d["fits"][k] = r
c = float(np.median(list(Au.values()))); d["A_u_common"] = c
for k, A in keep.items():
    for n, f in (("1e-6", 1e-6), ("0.01", 0.01), ("0.1", 0.1)):
        d["fits"][k]["off_seen_common_" + n] = float(a[seen][A[seen] < f * c].sum() / asn)
out = Path(__file__).with_name("shutdown_check.json")
if out.exists():
    raise SystemExit("exists")
out.write_text(json.dumps(d, indent=1))
print("special facets", d["special_facets"], "area %.0f m2" % d["special_area_m2"], "common A_u %.3g" % c)
print("%-30s %9s %9s | own A_u: off seen <1e-6 <.01 <.1 | whole <1e-6 <.01 | common: <1e-6 <.01 <.1 | area for 50/90/99%%" % ("fit", "A_u", "median/Au"))
for k, r in d["fits"].items():
    print("%-30s %9.3g %9.2g | %.3f %.3f %.3f | %.3f %.3f | %.3f %.3f %.3f | %.4f %.4f %.4f" % (
        k, r["A_u"], r["median_A_seen_over_A_u"], r["off_seen_1e-6"], r["off_seen_0.01"], r["off_seen_0.1"],
        r["off_whole_1e-6"], r["off_whole_0.01"], r["off_seen_common_1e-6"], r["off_seen_common_0.01"],
        r["off_seen_common_0.1"], r["seen_area_for_0.5_of_production"], r["seen_area_for_0.9_of_production"],
        r["seen_area_for_0.99_of_production"]))
for k in top:
    r = d["fits"][k]
    print("%-30s no-special A_rel" % k, " ".join("%7.3g" % v for v in r["A_rel_no_special"]),
          " special share %.3f" % r["special_share_of_seen_activity"])
