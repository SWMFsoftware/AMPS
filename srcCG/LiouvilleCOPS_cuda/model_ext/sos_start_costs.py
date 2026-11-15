# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/checks/sos_start_costs_20261114/sos_start_costs.py (on Voyager: /data001/highresolu_vt/cops_output/checks/sos_start_costs_20261114/sos_start_costs.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 339ed92052ec5f6d8189528c40dcd72dab25a395624813f456223c80e117fefb
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: check. Data cost of option-A starts without the GPU (lead 2026-11-14): direct facet sums with the fit's FP32
#   records (nonneg_contrib), W10x5, grain, D2's used rows. Gram 33.70346689911884 (= GPU start check 124998),
#   LM 9.519127043550094.
# USED BY: checks/sos_start_costs_20261114/run.sbatch, job 125100.
# ==== end of header ====
"""Data cost of option-A (sos, rank 3) starts WITHOUT the GPU (lead 2026-11-14): A = ||L^T y||^2 / A_s + A_s sum T
(m = 1, h = 0) from a start file, D2's kappa and grain, T held (182.09922, 136.11149, 133); the fitted forward model term
by term with the fit's FP32 records (nonneg_contrib.record_contributions, frozen), W10x5 for the late observations
(positive_problem.LateLayout), + grain for NG, floors, residuals w (ln p - ln obs) over D2's used rows (w = 1/sqrt(N)),
data cost 0.5 r.r. Validation: the Gram start must reproduce the GPU start check 3.370346689911884e+01 (job 124998).
    sos_start_costs.py <summary snapshot> <nonneg snapshot> <out.json> <start.npz> [<start.npz> ...]"""
import json, math, multiprocessing as mp, os, sys, time
from pathlib import Path
import numpy as np
SX, NX, OUTJ = sys.argv[1], sys.argv[2], sys.argv[3]
STARTS = sys.argv[4:]
for p_ in (SX, NX, "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"):
    sys.path.insert(0, p_) if p_ not in sys.path else None
sys.path.remove(SX); sys.path.insert(0, SX)          # the summary snapshot's nonneg_maps / contrib / positive_problem first
import nonneg_maps as NM, nonneg_contrib as NC, positive_problem as pp
sys.path.insert(1, NX)
import nonneg_problem as npb
R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
BASE, D2 = R / "base", R / "clip_refit/D2/final"
T0 = (182.09922, 136.11149, 133.0)
_W = {}

def work(task):
    r, ip = task
    kp, h, rec = NM.read_kernel_full(BASE / "kernels" / ("point=%04d.bin" % ip))
    face, out, _ = NC.record_contributions(rec, h, _W["A_list"], "fp32", r)
    return r, [(float(np.sum(ng)), float(np.sum(rg))) for ng, rg in out], h["et"]

def main():
    t0 = time.time()
    from cops_cuda import reference, device_store
    reference.bind()
    fc = reference.load("fit_cops"); at = reference.load("activity_templates")
    el = fc.read_elements(str(BASE / "elements.bin"))
    Y = np.ascontiguousarray(el["Y"], np.float64); area = np.asarray(el["area"])
    with np.load(D2 / "fit_restart.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"])); x = np.asarray(z["x"], float); names = list(meta["names"])
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        ref = {"NG": np.asarray(z["ref_ng"], float), "RG": np.asarray(z["ref_rg"], float)}
        use = {"NG": np.asarray(z["use_ng"]) > 0, "RG": np.asarray(z["use_rg"]) > 0}
        point_ids = np.asarray(z["point_ids"], np.int64)
    A_s = float(meta["enhancement_scale"])
    gA, gB = float(x[names.index("grain_A_uPa")]), float(x[names.index("grain_B_uPa_per_100s")])
    t_end = float(json.load(open(D2 / "grain_term.json"))["t_end_et"])
    fts, _ = npb.ft_templates_from(str(D2 / "fit_manifest.json"))
    TS = A_s * np.column_stack([at.select_templates(at.load_bundle(t["bundle"]["path"], el["centroid"], area), [t["name"]])
                                for t in fts]).sum(axis=1)
    A_list, info = [], []
    for sf in STARTS:
        with np.load(sf, allow_pickle=False) as z:
            L = np.asarray(z["L"], float); sm = json.loads(str(z["meta_json"]))
        if sm.get("d2_restart", {}).get("sha256") != NM.sha256_file(D2 / "fit_restart.npz"):
            raise SystemExit("%s: another D2 restart" % sf)
        ss = np.zeros(len(area))
        for c0 in range(0, len(area), 1 << 18):
            V = Y[c0:c0 + (1 << 18)] @ L
            ss[c0:c0 + (1 << 18)] = np.sum(V * V, axis=1)
        A = ss / A_s + TS
        A_list.append(A)
        info.append({"start_file": sf, "sha256": NM.sha256_file(sf), "refine": sm.get("refine", "none"),
                     "quality": {k: v for k, v in (sm.get("quality") or {}).items() if not isinstance(v, (list, dict))},
                     "A_min": float(A.min()), "A_sum_area": float(np.sum(A * area))})
    del Y
    fc_, header, cols, rows = NM.base_table(BASE)
    helio = fc_.header_value(header, "HeliocentricModel")
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    assert np.array_equal(ipoint, point_ids)
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, _ = NC.gauge_constants(fc_, header, TRG)
    NC._G.update({"ACO2": el["ACO2"], "T": T0, "kappa": kappa, "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio,
                  "AU": float(fc_.AU), "COS75": float(fc_.COS75), "KT": float(fc_.KT),
                  "fp32_record_fields": device_store.fp32_record_fields})
    _W["A_list"] = tuple(A_list)
    nrow = len(ipoint)
    sN = np.zeros((len(A_list), nrow)); sR = np.zeros((len(A_list), nrow)); et = np.zeros(nrow)
    tasks = sorted([(r, int(ipoint[r])) for r in range(nrow)], key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    with mp.get_context("fork").Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))) as pool:
        for r, vals, e in pool.imap_unordered(work, tasks, chunksize=1):
            for k, (n_, g_) in enumerate(vals):
                sN[k, r], sR[k, r] = n_, g_
            et[r] = e
    lay = pp.LateLayout(pp.read_point_table(str(R / "late_grid_points.txt")), pp.read_run_points(BASE), ipoint)
    W = lay.W.tocsr()
    s = et - t_end
    grain = np.where(s >= 0, np.maximum(gA * 1e-6 - gB * 1e-8 * np.maximum(s, 0.0), 0.0), 0.0)
    NNG, NRG = int(use["NG"].sum()), int(use["RG"].sum())
    for k in range(len(A_list)):
        pN, pR = sN[k].copy(), sR[k].copy()
        pN[lay.late_rows] = W @ sN[k][lay.grid_rows]
        pR[lay.late_rows] = W @ sR[k][lay.grid_rows]
        pN = np.maximum(pN + grain, 1e-300); pR = np.maximum(pR, 1e-300)
        rN = (np.log(pN[use["NG"]]) - np.log(ref["NG"][use["NG"]])) / math.sqrt(NNG)
        rR = (np.log(pR[use["RG"]]) - np.log(ref["RG"][use["RG"]])) / math.sqrt(NRG)
        info[k]["data_cost"] = 0.5 * float(rN @ rN + rR @ rR)
        info[k]["rms_ln"] = {"NG": float(np.sqrt(np.mean((rN * math.sqrt(NNG)) ** 2))),
                             "RG": float(np.sqrt(np.mean((rR * math.sqrt(NRG)) ** 2)))}
        print("%s: data cost %.15e, rms ln NG %.4f RG %.4f" % (info[k]["start_file"], info[k]["data_cost"],
                                                              info[k]["rms_ln"]["NG"], info[k]["rms_ln"]["RG"]), flush=True)
    out = {"created_utc": NM.utc_now(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
           "kappa": kappa, "grain_A_uPa": gA, "grain_B_uPa_per_100s": gB, "t_end_et": t_end, "T": T0,
           "NNG": NNG, "NRG": NRG, "starts": info, "gpu_start_check_gram_job_124998": 3.370346689911884e+01,
           "wall_seconds": time.time() - t0}
    json.dump(out, open(OUTJ, "x"), indent=1)
    return 0

if __name__ == "__main__":
    sys.exit(main())
