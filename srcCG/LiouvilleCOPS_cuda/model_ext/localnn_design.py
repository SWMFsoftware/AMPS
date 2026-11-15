#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_design.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_design.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8b360f326ae5431d526014bc3e8483df881dcd723246211bc1cb677885e8f07e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module + CLI (option C). Design matrices G^(b)[row, node] (NG, RG, 7 slope bins), template and CO2 columns
#   from the base's kernel records with the frozen nonneg_contrib arithmetic (FP32 records, FP64), and the direct facet
#   sums of D2, the softplus c9 state and two test activities; contrib_setup / direct_row used by the tests.
# RUN: jobs/localnn_build.sbatch STEP design (frozen model_ext_frozen/localnn_20261114T081912Z).
# ==== end of header ====
"""Option C (local non-negative basis, user 2026-11-14): the design matrices of the coarse hat functions, built ONCE from
the fitting base's kernel records with exactly the arithmetic of the frozen contribution code, plus the direct facet
sums that verify them. The fit is localnn_fit.py (localnn_problem.py); the grid is localnn_mesh.py.

Local, untracked model extension: nothing in the git tree is edited and no existing model_ext file is changed; the
frozen nonneg_contrib.py (model_ext_frozen/nonneg_contrib_20261112T210012Z, with its nonneg_maps.py) is imported
unchanged and its record_contributions() computes every per-record contribution used here.

LINEARITY. For fixed T (held: the original law), the modelled surface pressures are linear in the H2O activity and in
kappa (fit_cops.Model.evaluate; the frozen code per kernel record r = (point p, facet e)):
    NG_r = cNG_H bm_H wnH_r A_e u_r + kappa cNG_C bm_C wnC_r u_r
    RG_r = cRG_H,p (wpH_r A_e + bm_H wvH_r A_e u_r) + kappa cRG_C,p (phiC_r + bm_C wvC_r u_r)
(FP32 record fields Kn, KPhi, KV, cos of cops_cuda.device_store.fp32_record_fields, shadow folded, FP64 arithmetic,
F = q A with the night/shadow rule and r_h of the point). With g_r = the H2O contribution at A_e = 1 (record_contributions
with kappa = 0 and A = 1) and k_r = the CO2 contribution at kappa = 1 (A = 0), for A_e = m_b(e) (sum_j c_j phi_j(e) +
TS_e) [TS = A_s sum_k ft_ck, held at 1 (--templates on) or 0 (off)]:
    p_row = sum_b m_b ( sum_j G^(b)[row, j] c_j + Tcol^(b)[row] ) + kappa Ccol[row]
    G^(b)[row, j] = sum_{r in row, b(e_r) = b} g_r phi_j(e_r),   Tcol^(b)[row] = sum_{r, b(e_r) = b} g_r TS_e,
    Ccol[row] = sum_r k_r                                        (NG and RG separately; every base row, 1551)
The W10x5 time response of the late observations and the grain are applied in localnn_problem.py (both linear /
additive). Rounding: the design route sums g_r phi_j(e_r) per node instead of g_r A_e per record (identical up to the
order of FP64 additions).

DIRECT FACET SUMS computed in the same pass (record_contributions with the state's per-facet A): D2's final activity
max(0, activity_linear) (clip_refit/D2/final/fit_restart.npz), the bounded softplus c9 final state (finetune_state.npz
A_final), and two test activities of the basis (random c >= 0, m >= 1; without and with the templates). The tests
(test_localnn.py) and localnn_fit.py --verify compare them with the design route and with the written tables.
Checks in the pass: no negative or non-finite record contribution; every record with Kn, KPhi or KV > 0 lies on a
facet of the fit's seen mask (c9 nonneg_facets.npz seen = fit_cops.Model.seen_faces of the same base rows).

CLI: localnn_design.py --base BASE --grid GRID_DIR --slope-bins slope_bins.npz --facets-from <c9 fit>/nonneg_facets.npz
     --d2 <D2 fit dir> --c9-ft <c9 ft dir> --c9-fit <c9 fit dir> --contrib-code <frozen nonneg_contrib dir>
     --out <GRID_DIR>/design [--procs N]
Products (<out>.incomplete_<UTC> renamed at the end; never overwritten): design.npz (G_NG, G_RG [rows x 7 x nodes],
T_NG, T_RG [rows x 7], C_NG, C_RG [rows], the direct sums, the test activities' c and m, row metadata, meta_json),
design_report.json, README.txt.
"""

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
if str(APP) not in sys.path:
    sys.path.insert(0, str(APP))
SCHEMA = "localnn_design_v1"
NBINS = 7
LOG = "localnn_design:"
TEST_SEED = 20261114
_W = {}                 # worker globals (set before the pool forks)


def say(text):
    print("%s [%s] %s" % (LOG, dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(bufsize), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_module(path, name):
    path = Path(path)
    if str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))        # its own helpers (nonneg_maps) from the same frozen folder
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_grid(gdir):
    with np.load(Path(gdir) / "grid.npz", allow_pickle=False) as z:
        g = {k: np.asarray(z[k]) for k in ("node_xyz", "facet_node", "facet_w", "tri")}
        g["meta"] = json.loads(str(z["meta_json"]))
    g["n"] = len(g["node_xyz"])
    return g


def contrib_setup(base, ccode, T):
    """The frozen nonneg_contrib module with its globals set for the base (kappa 0) and the held T; returns
    (nc, NM, fit_cops, header, cols, rows, gauge-constant record, law record)."""
    nc = load_module(Path(ccode) / "nonneg_contrib.py", "nonneg_contrib_frozen")
    NM = nc.NM
    if Path(NM.__file__).resolve().parent != Path(ccode).resolve():
        raise SystemExit("nonneg_maps imported from %s, not from the frozen folder %s" % (NM.__file__, ccode))
    fc, header, cols, rows = NM.base_table(base)
    law = NM.law_guard(fc)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    if sza != "continuous":
        raise SystemExit("SZAMode %s" % sza)
    el = NM.read_elements(Path(base) / "elements.bin")
    TRG = np.array([float(r[cols["T_RG"]]) for r in rows])
    bm, cNG, cRG, gconst = nc.gauge_constants(fc, header, TRG)
    from cops_cuda import device_store
    nc._G.update({"ACO2": el["ACO2"], "T": tuple(T), "kappa": 0.0, "bm": bm, "cNG": cNG, "cRG": cRG, "helio": helio,
                  "AU": float(fc.AU), "COS75": float(fc.COS75), "KT": float(fc.KT),
                  "fp32_record_fields": device_store.fp32_record_fields})
    return {"nc": nc, "NM": NM, "fc": fc, "header": header, "cols": cols, "rows": rows, "gconst": gconst, "law": law,
            "el": el, "TRG": TRG}


def direct_row(S, base, row, ip, A, kappa):
    """Direct facet sum (NG, RG) [Pa] of one base row for a per-facet activity A and kappa (frozen arithmetic)."""
    nc, NM = S["nc"], S["NM"]
    kp, h, rec = NM.read_kernel_full(Path(base) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    old = nc._G["kappa"]
    nc._G["kappa"] = kappa
    try:
        _, ((ng, rg),), _ = nc.record_contributions(rec, h, (A,), "fp32", row)
    finally:
        nc._G["kappa"] = old
    return float(np.sum(ng)), float(np.sum(rg))


def test_activities(n, bins, fnode, fw, TS, A_s):
    """Two activities of the basis for the design-vs-direct checks: (c, m, A) without and with the templates."""
    rng = np.random.default_rng(TEST_SEED)
    out = []
    for with_t in (False, True):
        c = A_s * rng.uniform(0.0, 2.0, n)
        m = np.concatenate([[1.0], rng.uniform(1.0, 5.0, NBINS - 1)])
        u = c[fnode[:, 0]] * fw[:, 0] + c[fnode[:, 1]] * fw[:, 1] + c[fnode[:, 2]] * fw[:, 2]
        A = m[bins] * (u + TS) if with_t else m[bins] * u
        out.append((c, m, A))
    return out


def row_work(task):
    row, ip = task
    W = _W
    nc, NM = W["nc"], W["NM"]
    kp, h, rec = NM.read_kernel_full(Path(W["base"]) / "kernels" / ("point=%04d.bin" % ip))
    if kp != ip:
        raise RuntimeError("kernel file of point %d stores %d" % (ip, kp))
    nc._G["kappa"] = 0.0
    face, outs, conv = nc.record_contributions(rec, h, W["A_list"], "fp32", row)
    nc._G["kappa"] = 1.0
    _, ((kN, kR),), _ = nc.record_contributions(rec, h, (W["ZERO"],), "fp32", row)
    gN, gR = outs[0]
    nn = W["nn"]
    b = W["bins"][face]
    fn, fw = W["fnode"][face], W["fw"][face]
    res = {"row": row, "ip": ip, "records": int(len(rec)), "et": h["et"], "rh": h["rh"], "VdotB": h["VdotB"],
           "lit_threshold_adjusted": int((conv or {}).get("lit_threshold_adjusted", 0))}
    for name, g in (("NG", gN), ("RG", gR)):
        idx = np.concatenate([b.astype(np.int64) * nn + fn[:, k] for k in range(3)])
        wts = np.concatenate([g * fw[:, k] for k in range(3)])
        res["G_" + name] = np.bincount(idx, weights=wts, minlength=NBINS * nn).reshape(NBINS, nn)
        res["T_" + name] = np.bincount(b, weights=g * W["TS"][face], minlength=NBINS)
    res["C_NG"], res["C_RG"] = float(np.sum(kN)), float(np.sum(kR))
    for k, key in enumerate(W["A_keys"]):
        if k == 0:
            continue
        ng, rg = outs[k]
        res["sum_" + key] = (float(np.sum(ng)), float(np.sum(rg)))
    bad = 0
    for arr in (gN, gR, kN, kR):
        bad += int(np.sum(arr < 0) + np.sum(~np.isfinite(arr)))
    res["negative_or_nonfinite"] = bad
    nz = (rec["Kn"] > 0) | (rec["KPhi"] > 0) | (rec["KV"] > 0)
    res["records_nonzero_on_unseen"] = int(np.sum(nz & ~W["seen"][face]))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True)
    ap.add_argument("--grid", required=True)
    ap.add_argument("--slope-bins", required=True)
    ap.add_argument("--facets-from", required=True, help="nonneg_facets.npz of the c9 fit (seen, TS)")
    ap.add_argument("--d2", required=True, help="D2 fit directory (fit_restart.npz)")
    ap.add_argument("--c9-ft", required=True)
    ap.add_argument("--c9-fit", required=True)
    ap.add_argument("--contrib-code", required=True, help="frozen nonneg_contrib folder")
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    args = ap.parse_args(argv)
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    base = Path(args.base)
    ccode = Path(args.contrib_code)
    with np.load(Path(args.d2) / "fit_restart.npz", allow_pickle=False) as z:
        T = tuple(float(v) for v in np.asarray(z["T"]))
    CS = contrib_setup(base, ccode, T)
    nc, NM, fc, cols, rows, law, el = CS["nc"], CS["NM"], CS["fc"], CS["cols"], CS["rows"], CS["law"], CS["el"]
    gconst, TRG = CS["gconst"], CS["TRG"]
    nf = el["nfaces"]
    grid = load_grid(args.grid)
    nn = grid["n"]
    if len(grid["facet_node"]) != nf:
        raise SystemExit("grid has %d facets, the base %d" % (len(grid["facet_node"]), nf))
    with np.load(args.slope_bins, allow_pickle=False) as z:
        names = [str(v) for v in z["names"]]
        vals = np.asarray(z["values"])
    if vals.shape != (nf, NBINS):
        raise SystemExit("slope bins %s" % (vals.shape,))
    if not np.all(vals.sum(axis=1) == 1.0) or not np.all((vals == 0) | (vals == 1)):
        raise SystemExit("slope-bin columns are not a partition")
    bins = np.argmax(vals, axis=1).astype(np.int8)
    del vals
    with np.load(args.facets_from, allow_pickle=False) as z:
        seen = np.asarray(z["seen"], bool)
        TS = np.asarray(z["TS"], np.float64)
        fbins = np.asarray(z["bin"])
        fmeta = json.loads(str(z["meta_json"]))
    if not np.array_equal(fbins.astype(np.int8), bins):
        raise SystemExit("the slope bins differ from the facets file's")
    A_s = float(fmeta["A_scale"])
    if TS.min() < 0:
        raise SystemExit("the ft templates are negative somewhere (TS min %.3e): the wrapper is a decision" % TS.min())
    # D2 state
    d2 = Path(args.d2)
    with np.load(d2 / "fit_restart.npz", allow_pickle=False) as z:
        A_D2 = np.maximum(np.asarray(z["activity_linear"], np.float64), 0.0)
        point_ids = np.asarray(z["point_ids"], np.int64)
        kappa_D2 = float(np.asarray(z["kappa"]).reshape(-1)[0])
        tnames = json.loads(str(z["meta_json"]))["templates"]["names"]
        tcols = np.asarray(z["template_columns"])
        meta_d2 = json.loads(str(z["meta_json"]))
    ft_idx = [k for k, nme in enumerate(tnames) if nme.startswith("ft_c")]
    TS_re = A_s * np.sum(tcols[:, ft_idx], axis=1)
    del tcols
    ts_check = {"A_scale_D2": float(meta_d2["enhancement_scale"]), "A_scale_facets_file": A_s,
                "max_abs_diff_TS_recomputed_from_D2_columns": float(np.max(np.abs(TS_re - TS))),
                "max_TS": float(TS.max()), "nonzero_facets": int(np.sum(TS != 0)), "min_TS": float(TS.min())}
    del TS_re
    if ts_check["A_scale_D2"] != A_s:
        raise SystemExit("A_scale differs between D2 and the facets file")
    say("templates: " + json.dumps(ts_check))
    ipoint = np.array([int(r[cols["iPoint"]]) for r in rows])
    if not np.array_equal(ipoint, point_ids):
        raise SystemExit("D2 restart rows differ from the base rows")
    # c9 state
    with np.load(Path(args.c9_ft) / "finetune_state.npz", allow_pickle=False) as z:
        A_C9 = np.asarray(z["A_final"], np.float64)
    with np.load(Path(args.c9_fit) / "fit_restart.npz", allow_pickle=False) as z:
        T9 = tuple(float(v) for v in np.asarray(z["T"]))
        if not np.array_equal(np.asarray(z["point_ids"], np.int64), point_ids):
            raise SystemExit("c9 restart rows differ")
    if T9 != T:
        raise SystemExit("c9 temperatures %s differ from D2's %s" % (T9, T))
    tests = test_activities(nn, bins.astype(np.intp), grid["facet_node"], grid["facet_w"], TS, A_s)
    ONE = np.ones(nf)
    A_keys = ["unit", "D2", "c9", "test1", "test2"]
    _W.update({"nc": nc, "NM": NM, "base": str(base), "A_list": (ONE, A_D2, A_C9, tests[0][2], tests[1][2]),
               "A_keys": A_keys, "ZERO": np.zeros(nf), "nn": nn, "bins": bins.astype(np.intp),
               "fnode": grid["facet_node"].astype(np.int64), "fw": grid["facet_w"], "TS": TS, "seen": seen})
    tasks = [(r, int(ipoint[r])) for r in range(len(ipoint))]
    tasks.sort(key=lambda t: -int(float(rows[t[0]][cols["nVisible"]])))
    say("%d rows, %d nodes, %d processes" % (len(tasks), nn, args.procs))
    nrow = len(ipoint)
    G = {g: np.zeros((nrow, NBINS, nn)) for g in ("NG", "RG")}
    Tc = {g: np.zeros((nrow, NBINS)) for g in ("NG", "RG")}
    Cc = {g: np.zeros(nrow) for g in ("NG", "RG")}
    S = {k: np.zeros((nrow, 2)) for k in A_keys[1:]}
    meta_rows = {k: np.zeros(nrow) for k in ("et", "rh", "VdotB")}
    nrec = np.zeros(nrow, dtype=np.int64)
    bad = unseen = lit_adj = 0
    done = 0
    with mp.get_context("fork").Pool(args.procs) as pool:
        for res in pool.imap_unordered(row_work, tasks, chunksize=1):
            r = res["row"]
            for g in ("NG", "RG"):
                G[g][r] = res["G_" + g]
                Tc[g][r] = res["T_" + g]
                Cc[g][r] = res["C_" + g]
            for k in A_keys[1:]:
                S[k][r] = res["sum_" + k]
            for k in meta_rows:
                meta_rows[k][r] = res[k]
            nrec[r] = res["records"]
            bad += res["negative_or_nonfinite"]
            unseen += res["records_nonzero_on_unseen"]
            lit_adj += res["lit_threshold_adjusted"]
            done += 1
            if done % 200 == 0 or done == len(tasks):
                say("  %d / %d rows (%.0f s)" % (done, len(tasks), time.time() - t0))
    checks = {"records": int(nrec.sum()), "negative_or_nonfinite_record_contributions": bad,
              "records_with_nonzero_kernel_on_unseen_facets": unseen, "fp32_lit_threshold_adjustments": lit_adj,
              "templates": ts_check}
    say("checks: " + json.dumps(checks))
    if bad or unseen:
        raise SystemExit("negative / non-finite contributions or records on unseen facets: stop")
    # self-consistency of the two test activities (design route vs direct sum, every row)
    tc = {}
    for k, (c, m, _) in enumerate(tests):
        key = "test%d" % (k + 1)
        for gi, g in enumerate(("NG", "RG")):
            route = np.einsum("b,rbj,j->r", m, G[g], c)
            if k == 1:
                route = route + Tc[g] @ m
            d = np.abs(route / S[key][:, gi] - 1.0)
            tc["%s_%s_max_rel" % (key, g)] = float(np.max(d))
    say("design route vs direct sums (test activities, all rows): " + json.dumps(tc))
    stage.mkdir(parents=True)
    inputs = {"base_COPS_model": sha256_file(base / "COPS_model.dat"), "grid_npz": sha256_file(Path(args.grid) / "grid.npz"),
              "slope_bins": sha256_file(args.slope_bins), "facets_from": sha256_file(args.facets_from),
              "d2_fit_restart": sha256_file(d2 / "fit_restart.npz"),
              "c9_finetune_state": sha256_file(Path(args.c9_ft) / "finetune_state.npz"),
              "c9_fit_restart": sha256_file(Path(args.c9_fit) / "fit_restart.npz")}
    paths = {"base": str(base.resolve()), "grid": str(Path(args.grid).resolve()),
             "slope_bins": str(Path(args.slope_bins).resolve()), "facets_from": str(Path(args.facets_from).resolve()),
             "d2": str(d2.resolve()), "c9_ft": str(Path(args.c9_ft).resolve()), "c9_fit": str(Path(args.c9_fit).resolve())}
    code = {"localnn_design.py": sha256_file(__file__), "nonneg_contrib.py": sha256_file(ccode / "nonneg_contrib.py"),
            "nonneg_maps.py": sha256_file(ccode / "nonneg_maps.py"), "contrib_code_dir": str(ccode.resolve())}
    meta = {"schema": SCHEMA, "created_utc": utc_now(), "host": socket.gethostname(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv, "paths": paths, "inputs_sha256": inputs,
            "code": code, "law": law, "gauge_constants": gconst, "T": T, "A_scale": A_s, "bin_names": names,
            "kappa_D2": kappa_D2, "nodes": nn, "rows": nrow, "checks": checks, "design_vs_direct_test": tc,
            "precision": "records as the fit: cops_cuda.device_store.fp32_record_fields (FP32 Kn, KPhi, KV, cos; shadow "
                         "-> cos -2), FP64 arithmetic (the frozen nonneg_contrib.record_contributions)",
            "definitions": {"G_NG/G_RG[row, b, j]": "sum over the row's records r with bin(e_r) = b of g_r phi_j(e_r); g_r "
                                                    "= H2O contribution [Pa] at A_e = 1 (kappa 0)",
                            "T_NG/T_RG[row, b]": "sum_{r, bin b} g_r TS_e, TS = A_s sum_k ft_ck (c9 nonneg_facets.npz TS)",
                            "C_NG/C_RG[row]": "CO2 contribution [Pa] at kappa = 1",
                            "sum_<state>[row, (NG, RG)]": "direct facet sum of the H2O contributions of the state's "
                                                          "per-facet A at kappa 0 (add kappa * C for the total)"},
            "wall_seconds": time.time() - t0}
    np.savez(stage / "design.npz", G_NG=G["NG"], G_RG=G["RG"], T_NG=Tc["NG"], T_RG=Tc["RG"], C_NG=Cc["NG"],
             C_RG=Cc["RG"], sum_D2=S["D2"], sum_c9=S["c9"], sum_test1=S["test1"], sum_test2=S["test2"],
             test1_c=tests[0][0], test1_m=tests[0][1], test2_c=tests[1][0], test2_m=tests[1][1], iPoint=ipoint,
             records=nrec, et=meta_rows["et"], rh=meta_rows["rh"], VdotB=meta_rows["VdotB"], TRG=TRG,
             meta_json=json.dumps(meta, sort_keys=True, default=float))
    with open(stage / "design_report.json", "x") as f:
        json.dump(meta, f, indent=1, sort_keys=True, default=float)
    with open(stage / "README.txt", "x") as f:
        f.write("Design matrices of the local non-negative basis (option C), %s, %s (job %s)\n" % (
            meta["created_utc"], meta["host"], meta["slurm_job_id"]))
        f.write("Rows %d, nodes %d, bins %d, records %d; checks %s\n" % (nrow, nn, NBINS, checks["records"],
                                                                      json.dumps(checks)))
        f.write("Design route vs direct facet sums of two test activities (all rows): %s\n" % json.dumps(tc))
        f.write("Definitions: see design_report.json and the docstring of localnn_design.py.\nCOMMAND\n  %s\n"
                % " ".join(sys.argv))
    stage.rename(out)
    say("done: %s (%.0f s)" % (out, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
