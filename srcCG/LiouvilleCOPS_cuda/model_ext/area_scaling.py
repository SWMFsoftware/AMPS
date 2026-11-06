#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/area_scaling.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/area_scaling.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 209444077a4a387500048cd733387802bfbdc5f54bb91ab6179bbd4adb85d9a2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: test without refitting: scales the fine-tuned activity on the facets that make the modelled pressure of
#   chosen points (main contributors) and recomputes all NG and RG pressures exactly (linear in the activity).
# USED FOR: the 12M area-scaling test (job 121943; last RG area x1.2, 4th-6th NG from the end x1.5).
# RUN: .venv/bin/python <this> --fit <fit dir> --finetune <fine-tuning dir> --run <run dir> --out <new dir>
#        [--fraction 0.9 --rg-last 1 --rg-factor 1.2 --ng-from-end 4-6 --ng-factor 1.5] (or area_scaling.sbatch).
#   No restart file is written.
# ==== end of header ====
"""Scale the fine-tuned H2O activity on the surface areas that make the modelled pressure of chosen observation
points, and compute the changed NG and RG pressures of every point (user request 2026-11-04). Nothing is refitted
and no restart file is written or changed.

    .venv/bin/python <this> --fit <fit dir> --finetune <fine-tuning dir> --run <base> --out <dir>
        [--fraction 0.9] [--rg-last 1] [--rg-factor 1.2] [--ng-from-end 4 6] [--ng-factor 1.5]
        [--devices 0,...,7] [--precision mixed]

User answers (2026-11-04): "area that contributed" = the main contributors: per target point, the smallest set of
facets whose H2O contributions make --fraction (90 %) of the modelled H2O pressure of that gauge at that point
(contribution of facet e at point p = H_p,e * A_e, A = the fine-tuned clipped activity, H = gauge pressure per unit
clipped activity at the frozen temperatures; CO2 and the grain are not area contributions and are not scaled);
no confinement (the other points change as they will); factors multiply A: RG area x1.2, NG area (union over the
NG points 4..6 from the end) x1.5, facets in both x1.2*1.5 = 1.8.

Exact: the pressures are linear in each facet's clipped activity and the factors are >= 1 on A >= 0, so
p_new = p_finetuned + dp, dp = point_changes(delta A = (f - 1) A) (cops_cuda/restart.py; the fine tuning's own
evaluator and kernels). p_finetuned = 'after' of <finetune>/pressure_comparison.txt (NG includes the tuned grain).
Products in --out: area_scaling.json (targets, sets, areas, changes), regions.npz (facet indices, factors),
tables/{finetuned,scaled}.dat (model_obs_figures.py tables), region map png.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path[:0] = [str(APP), str(APP / "tools")]


def table(path, header):
    names, rows = None, {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(header):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip() or names is None:
            continue
        v = line.split()
        rows[int(v[0])] = dict(zip(names, v))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, type=Path)
    ap.add_argument("--finetune", required=True, type=Path)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--fraction", type=float, default=0.9)
    ap.add_argument("--rg-last", type=int, default=1, help="the RG point: k-th fitted RG point from the end")
    ap.add_argument("--rg-factor", type=float, default=1.2)
    ap.add_argument("--ng-from-end", type=int, nargs=2, default=(4, 6))
    ap.add_argument("--ng-factor", type=float, default=1.5)
    ap.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    ap.add_argument("--precision", choices=("fp64", "mixed"), default="mixed")
    a = ap.parse_args()
    if a.out.exists():
        sys.exit("area_scaling: %s exists" % a.out)
    from cops_cuda import precision as cuda_precision
    os.environ[cuda_precision.RECORD_PRECISION_ENVIRONMENT] = "fp32" if a.precision == "mixed" else "fp64"
    from cops_cuda import fine_tuning as ftg, restart as rs

    r = rs.read(a.fit / "fit_restart.npz")
    st = np.load(a.finetune / "fine_tuning_state.npz")
    A = np.asarray(st["activity_final"], float)                     # fine-tuned clipped H2O activity per facet
    pc = table(a.finetune / "pressure_comparison.txt", "# iPoint")
    base = table(a.fit / "COPS_model.dat", "#iPoint")
    ids = sorted(pc)
    useNG = [i for i in ids if pc[i]["useNG"] == "1" and float(pc[i]["refNG"]) > 0]
    useRG = [i for i in ids if pc[i]["useRG"] == "1" and float(pc[i]["refRG"]) > 0]
    rg_target = [useRG[-a.rg_last]]
    k0, k1 = a.ng_from_end
    ng_targets = [useNG[-k] for k in range(k1, k0 - 1, -1)]          # chronological

    ev = ftg.SplitEvaluator(a.run, [int(d) for d in a.devices.split(",")])
    try:
        pNG_s, pRG_s = ev.pressures(r)                               # binds every part (fit state, frozen T)
        row_of = {int(p): k for k, p in enumerate(ev.point_ids)}
        area = ev.elements["area"]

        def contributions(pid, gauge):
            for part in ev.parts:
                if pid in set(int(x) for x in part.point_ids):
                    sn, sr, pn, pr = part.facet_point_shares(r, A, [pid])
                    local = int(np.nonzero(part.point_ids == pid)[0][0])
                    return (sn * pn[local]) if gauge == "NG" else (sr * abs(pr[local]))
            raise ValueError("point %d not loaded" % pid)

        def main_set(c):
            order = np.argsort(-c, kind="stable")
            cum = np.cumsum(c[order])
            total = float(cum[-1])
            n = int(np.searchsorted(cum, a.fraction * total) + 1)
            return order[:n], total, float(cum[n - 1])

        targets, sets = [], {"RG": [], "NG": []}
        for gauge, pids in (("RG", rg_target), ("NG", ng_targets)):
            for pid in pids:
                c = contributions(pid, gauge)
                idx, total, kept = main_set(c)
                sets[gauge].append(idx)
                p_after = float(pc[pid]["after" + gauge])
                targets.append({"gauge": gauge, "iPoint": pid, "UTC": pc[pid]["UTC"],
                                "altitude_m": float(base[pid]["altitude"]),
                                "seen_facets": int(np.count_nonzero(c > 0)),
                                "h2o_pressure_from_facets_Pa": total, "kept_fraction": kept / total,
                                "facets_kept": int(len(idx)), "area_kept_m2": float(area[idx].sum()),
                                "finetuned_model_Pa": p_after, "observed_Pa": float(pc[pid]["ref" + gauge])})
                print("%s point %d (%s): %d facets seen, %d facets give %.1f%% of %.4e Pa (model %.4e Pa)" % (
                    gauge, pid, pc[pid]["UTC"], targets[-1]["seen_facets"], len(idx), 100 * kept / total, total,
                    p_after), flush=True)
        rg = np.zeros(len(A), bool)
        ng = np.zeros(len(A), bool)
        for idx in sets["RG"]:
            rg[idx] = True
        for idx in sets["NG"]:
            ng[idx] = True
        f = np.ones(len(A))
        f[rg] *= a.rg_factor
        f[ng] *= a.ng_factor
        delta = (f - 1.0) * A
        dNG, dRG = ev.point_changes(r, delta)
        fac_ng, fac_rg = ev.elements["centroid"], None
    finally:
        ev.release()

    a.out.mkdir(parents=True)
    (a.out / "tables").mkdir()
    rows_ft, rows_sc, change = [], [], []
    for pid in ids:
        k = row_of[pid]
        b, t = base[pid], pc[pid]
        n0, r0 = float(t["afterNG"]), float(t["afterRG"])
        n1, r1 = n0 + dNG[k], r0 + dRG[k]
        common = [str(pid), b["UTC"], b["et"], b["altitude"]]
        rest = [t["refNG"], t["refRG"], t["useNG"], t["useRG"]]
        rows_ft.append(" ".join(common + [t["afterNG"], t["afterRG"]] + rest))
        rows_sc.append(" ".join(common + ["%.17e" % n1, "%.17e" % r1] + rest))
        change.append((pid, n1 / n0 - 1.0, r1 / r0 - 1.0))
    head = "#iPoint UTC et altitude pNG pRG refNG refRG useNG useRG"
    (a.out / "tables/finetuned.dat").write_text(
        "# fine-tuned model (after, %s/pressure_comparison.txt)\n%s\n%s\n" % (a.finetune, head, "\n".join(rows_ft)))
    (a.out / "tables/scaled.dat").write_text(
        "# fine-tuned model with the activity of the target areas scaled (area_scaling.py)\n%s\n%s\n"
        % (head, "\n".join(rows_sc)))
    ch = np.array([(c[1], c[2]) for c in change])
    nz = np.nonzero((np.abs(ch) > 1e-12).any(axis=1))[0]
    first_changed = {"iPoint": int(change[nz[0]][0]), "UTC": pc[change[nz[0]][0]]["UTC"]} if len(nz) else None
    for t in targets:
        k = [c[0] for c in change].index(t["iPoint"])
        t["relative_change_NG"], t["relative_change_RG"] = float(ch[k, 0]), float(ch[k, 1])
    summary = {"fit": str(a.fit), "finetune": str(a.finetune), "run": str(a.run), "fraction": a.fraction,
               "rg_factor": a.rg_factor, "ng_factor": a.ng_factor, "overlap_factor": a.rg_factor * a.ng_factor,
               "precision": a.precision, "targets": targets,
               "facets": {"rg_area": int(rg.sum()), "ng_area": int(ng.sum()), "both": int((rg & ng).sum())},
               "area_m2": {"rg_area": float(area[rg].sum()), "ng_area": float(area[ng].sum()),
                           "both": float(area[rg & ng].sum())},
               "first_changed_point": first_changed,
               "max_relative_change": {"NG": float(np.abs(ch[:, 0]).max()), "RG": float(np.abs(ch[:, 1]).max())},
               "points_changed_by_more_than_1e-3": int(np.count_nonzero((np.abs(ch) > 1e-3).any(axis=1))),
               "per_point_relative_change_last_30": [
                   {"iPoint": int(c[0]), "UTC": pc[c[0]]["UTC"], "NG": float(c[1]), "RG": float(c[2])}
                   for c in change[-30:]]}
    (a.out / "area_scaling.json").write_text(json.dumps(summary, indent=1))
    np.savez_compressed(a.out / "regions.npz", rg_area=np.flatnonzero(rg), ng_area=np.flatnonzero(ng),
                        factor_nonunit=np.flatnonzero(f != 1.0), factor=f[f != 1.0])
    # region map: facet centroids in body-frame longitude/latitude
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    c = fac_ng
    lon = np.degrees(np.arctan2(c[:, 1], c[:, 0]))
    lat = np.degrees(np.arcsin(c[:, 2] / np.linalg.norm(c, axis=1)))
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    for axi, zoom in zip(ax, (False, True)):
        for m, col, lab in ((rg & ~ng, "#0072B2", "RG area only (x%.2g)" % a.rg_factor),
                            (ng & ~rg, "#D55E00", "NG area only (x%.2g)" % a.ng_factor),
                            (rg & ng, "#009E73", "both (x%.2g)" % (a.rg_factor * a.ng_factor))):
            if m.any():
                axi.scatter(lon[m], lat[m], s=1, c=col, label="%s: %d facets" % (lab, m.sum()), rasterized=True)
        axi.set_xlabel("longitude [deg]")
        axi.set_ylabel("latitude [deg]")
        if zoom:
            sel = rg | ng
            axi.set_xlim(lon[sel].min() - 1, lon[sel].max() + 1)
            axi.set_ylim(lat[sel].min() - 1, lat[sel].max() + 1)
            axi.set_title("zoom")
        else:
            axi.set_xlim(-180, 180)
            axi.set_ylim(-90, 90)
            axi.set_title("whole nucleus (body frame, facet centroids)")
        axi.legend(markerscale=8, fontsize=8, loc="best")
    fig.suptitle("Scaled areas: %d%% of the modelled H2O pressure of RG point %s and NG points %s" % (
        round(100 * a.fraction), rg_target, ng_targets))
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(a.out / ("region_map." + ext), dpi=300)
    print(json.dumps({k: summary[k] for k in ("facets", "area_m2", "max_relative_change", "first_changed_point",
                                              "points_changed_by_more_than_1e-3")}), flush=True)
    for t in targets:
        print("target %s %d: NG %+.3f%%, RG %+.3f%%" % (t["gauge"], t["iPoint"], 100 * t["relative_change_NG"],
                                                        100 * t["relative_change_RG"]), flush=True)


if __name__ == "__main__":
    main()
