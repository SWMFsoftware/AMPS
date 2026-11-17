#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/kappa_table.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/kappa_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ad24330fba65a1a2446d67825dfa70320474b2ff748272f4b1f1b1438f3aec05
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: kappa (CO2 scale, lower bound 0) of every stage-2 chain (first and last fit), the stage-1 fits, the softplus
#   chains and D2, with the at-bound test (2026-11-17)
# RUN: login node: kappa_table.py --spec figures/data/stage2_final_spec_20261117.json --out-dir <folder>/kappa_check
# ==== end of header ====
"""kappa (the CO2 scale on the srcCG CO2 map; lower bound 0) of the stage-2 chains and of the earlier fits (2026-11-17, lead:
"state kappa per chain (value + at bound) for all options and D2 ... say whether this was already the case before stage 2,
with evidence"). Read only, small members only (login node).
    kappa_table.py --spec <stage-2 spec.json> --out-dir <dir>
Per fit directory: A, B, D, softplus, D2 (fit_cops / nonneg fits): fit_restart.npz names / x / lo / hi (kappa) and the
"kappa <initial> <final> <sigma>" row of fit_result.txt (both must agree to 1e-12 relative); option C: the fit's
localnn_state.npz is not read - analysis/report.json kappa and bounds.kappa_at_0 (stage-2 chains: c<N>/analysis; stage 1:
<run>/analysis) and the fit_result.txt kappa row when present. "at bound": x - lo <= 1e-9 max(1, |lo|) (C: the report's
flag, x <= 1e-8 in localnn_report.py); a C fit without a report (the chains' cycle-1 fits): the fit_result.txt kappa row,
at bound when <= 1e-8 (the report's criterion). kappa = 1 is the srcCG CO2 map at face value (every fit starts there).
Writes kappa_table.txt / kappa_table.json (refuses to overwrite).
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

RUN = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
NN = RUN / "nonneg"


def from_fit(fit_dir):
    fit_dir = Path(fit_dir)
    out = {"dir": str(fit_dir)}
    t = (fit_dir / "fit_result.txt").read_text() if (fit_dir / "fit_result.txt").exists() else ""
    m = re.search(r"^kappa\s+(\S+)\s+(\S+)\s+(\S+)\s*$", t, re.M)
    if m:
        out["fit_result_kappa"] = float(m.group(2))
    rs = fit_dir / "fit_restart.npz"
    if rs.exists():
        with np.load(rs, allow_pickle=False) as z:
            names = [str(x) for x in json.loads(str(z["meta_json"]))["names"]]
            if "kappa" in names:
                i = names.index("kappa")
                x, lo, hi = float(z["x"][i]), float(z["lo"][i]), float(z["hi"][i])
                out.update({"kappa": x, "lo": lo, "hi": hi, "at_lower_bound": bool(x - lo <= 1e-9 * max(1.0, abs(lo))),
                            "source": "fit_restart.npz"})
    if "kappa" in out and "fit_result_kappa" in out:
        a, b = out["kappa"], out["fit_result_kappa"]
        out["fit_result_agrees"] = bool(abs(a - b) <= 1e-6 * max(abs(a), abs(b), 1e-300) or abs(a - b) < 1e-40)
    return out


def from_report(rep, fit_dir=None):
    d = json.load(open(rep))
    out = {"kappa": float(d["kappa"]), "lo": 0.0, "at_lower_bound": bool(d["bounds"].get("kappa_at_0")),
           "source": str(rep)}
    if fit_dir is not None:
        f = from_fit(fit_dir)
        if "fit_result_kappa" in f:
            out["fit_result_kappa"] = f["fit_result_kappa"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    O = Path(a.out_dir)
    O.mkdir(parents=True, exist_ok=True)
    for f in ("kappa_table.txt", "kappa_table.json"):
        if (O / f).exists():
            raise SystemExit("%s exists" % (O / f))
    spec = json.load(open(a.spec))
    rows = []

    def add(group, name, info):
        info.update({"group": group, "name": name})
        rows.append(info)

    # before stage 2
    add("before stage 2", "D2 (clipped refit, start of everything)", from_fit(RUN / "clip_refit/D2/final"))
    for c in (5,):
        add("before stage 2", "softplus unbounded multipliers, chain_5cycles c%d fit" % c,
            from_fit(NN / ("softplus/chain_5cycles_20261112T012553Z/c%d/fit/final" % c)))
    for c in (1, 5, 9):
        add("before stage 2", "softplus multipliers >= 1, chain_5cycles_smm1 c%d fit" % c,
            from_fit(NN / ("softplus/chain_5cycles_smm1_20261113T041215Z/c%d/fit/final" % c)))
    for rel in ("background/stage1_bg1e-3_20261114T084711Z", "background/stage1_bg1e-2_20261114T084711Z",
                "background/stage1_bg1e-1_20261114T084711Z", "sos/stage1_r3_20261114T091354Z",
                "logspread/stage1_f0.1_20261114T090251Z", "logspread/stage1_f1_20261114T090251Z",
                "logspread/stage1_f10_20261114T090251Z", "background/stage1c_bg1e-3_20261115T023450Z",
                "background/stage1c_bg1e-2_20261115T023450Z", "background/stage1c_bg1e-1_20261115T023450Z",
                "sos/stage1c_r3_20261115T023450Z", "logspread/stage1c_f0.1_20261115T023450Z",
                "logspread/stage1c_f1_20261115T023450Z", "logspread/stage1c_f10_20261115T023450Z"):
        add("before stage 2", "stage 1 " + rel.split("/")[1].rsplit("_", 1)[0] + " (" + rel.split("/")[0] + ")",
            from_fit(NN / rel / "fit/final"))
    for rel in ("localnn/stage1_anch_f0.1_20261114T094256Z", "localnn/stage1_anch_f1_20261114T222136Z",
                "localnn/stage1_anch_f10_20261114T094259Z", "localnn/stage1_anch_f10_pm0.01_cont_20261115T020243Z"):
        add("before stage 2", "stage 1 C " + rel.split("/")[1].rsplit("_", 1)[0],
            from_report(NN / rel / "analysis/report.json", NN / rel / "fit/final"))
    # stage 2: the last fit of every chain of the spec (and its first fit)
    for c in spec["chains"]:
        root = NN / c["root"]
        n = c["cycles"]
        for k in (1, n):
            fit = root / ("c%d" % k) / "fit" / "final"
            if c["kind"] == "localnn":
                rep = root / ("c%d" % k) / "analysis" / "report.json"
                if rep.exists():
                    info = from_report(rep, fit)
                elif (fit / "fit_result.txt").exists():
                    info = from_fit(fit)
                    if "fit_result_kappa" in info:
                        info.update({"kappa": info["fit_result_kappa"], "lo": 0.0,
                                     "at_lower_bound": bool(info["fit_result_kappa"] <= 1e-8)})
                    info["source"] = "fit_result.txt kappa row (no report; at bound: <= 1e-8 as the C report)"
                else:
                    continue
            else:
                if not (fit / "fit_restart.npz").exists():
                    continue
                info = from_fit(fit)
            add("stage 2 %s" % ("regularized" if c["reg"] else "unregularized"), "%s c%d fit" % (c["key"], k), info)
    L = ["# kappa (CO2 scale, lower bound 0) of the fits; at bound: x - lo <= 1e-9 (fit_restart.npz), or the option-C "
         "report's kappa_at_0; fit_result = the kappa row of fit_result.txt (final value)", "",
         "%-16s %-62s %-12s %-6s %-12s %s" % ("group", "fit", "kappa", "bound", "fit_result", "source")]
    for r in rows:
        L.append("%-16s %-62s %-12s %-6s %-12s %s" % (
            r["group"][:16], r["name"][:62], ("%.3e" % r["kappa"]) if "kappa" in r else "-",
            ("yes" if r.get("at_lower_bound") else "no") if "kappa" in r else "-",
            ("%.3e" % r["fit_result_kappa"]) if "fit_result_kappa" in r else "-", r.get("source", "-")))
    n_all = len([r for r in rows if "kappa" in r])
    n_b = len([r for r in rows if r.get("at_lower_bound")])
    mx = max((r["kappa"] for r in rows if "kappa" in r), default=float("nan"))
    L += ["", "fits with kappa: %d; at the lower bound 0: %d; largest kappa %.3e" % (n_all, n_b, mx)]
    (O / "kappa_table.txt").write_text("\n".join(L) + "\n")
    json.dump({"rows": rows, "n": n_all, "n_at_bound": n_b, "max_kappa": mx, "command": sys.argv},
              open(O / "kappa_table.json", "x"), indent=1)
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
