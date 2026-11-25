#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_seen_gauges.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_seen_gauges.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 110b1a9481bd81b18e34b5c9d6b9f1d3db65b3eb2b0ad65b019bf8cd2cdd1886
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: per-facet visibility by NG, RG, both in the same row, lit (kernel pass) (2026-11-21)
# RUN: see OUTPUT_VISUALIZATION.md 3.21
# ==== end of header ====
"""Which facets each COPS gauge sees, per observation (2026-11-21; user: "do this analysis separately for faces seen at
least by one gauge and then for faces seen by both gauges concurrently or in different time").

One pass over the base's kernel files (kernels/point=NNNN.bin, fit_cops.read_kernel format: face, Kn, KPhi, KV, cos,
shadow). NG sees a facet in a row if Kn > 0; RG if KPhi > 0 or KV > 0 (fit_cops' "seen" = any of the three). It records per facet, for all 1551 rows (*_all) and for the fitted instants (*_obs): seen by NG, seen by RG, seen by
both in the SAME row ("concurrently"). Fitted instants of a gauge = the observation rows where that gauge has data (NG
934, RG 1150) + the 401 1-s grid rows of the late response average (both gauges); concurrent = the same facet in the same
row, where both gauges have data. "Both at different times only" = seen by both & not concurrently. *_lit = the same
with the facet lit (cos SZA > 0, no cast shadow) at that instant: a facet seen only in the dark constrains A 50 times
more weakly (q_min = 0.02 q_max). Fit-independent (geometry only).

    q002_seen_gauges.py --base <base_q002> --regions <run>/contrib_phase*/regions.npz --out <npz> [--procs 16]
"""
import argparse
import multiprocessing as mp
import struct
import time
from pathlib import Path

import numpy as np

REC = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])
G = {}


def read(fname):
    with open(fname, "rb") as f:
        version, ip = struct.unpack("<ii", f.read(8))
        if version != 1:
            raise SystemExit("%s: kernel version %d" % (fname, version))
        f.read(80)
        (nk,) = struct.unpack("<i", f.read(4))
        rec = np.fromfile(f, dtype=REC, count=nk)
    if len(rec) != nk:
        raise SystemExit("%s: truncated" % fname)
    return ip, rec


def work(chunk):
    base, nf, useN, useR = G["base"], G["nf"], G["useN"], G["useR"]
    M = {k: np.zeros(nf, bool) for k in ("ng_obs", "rg_obs", "both_obs", "ng_all", "rg_all", "both_all",
                                         "ng_lit", "rg_lit", "both_lit")}
    n = 0
    for ip in chunk:
        got, rec = read(base / "kernels" / ("point=%04d.bin" % ip))
        if got != ip:
            raise SystemExit("kernel file %d says point %d" % (ip, got))
        face = rec["face"]
        ng = rec["Kn"] > 0
        rg = (rec["KPhi"] > 0) | (rec["KV"] > 0)
        fn, fr = face[ng], face[rg]
        rowN = np.zeros(nf, bool); rowN[fn] = True
        rowR = np.zeros(nf, bool); rowR[fr] = True
        both = rowN & rowR
        M["ng_all"] |= rowN; M["rg_all"] |= rowR; M["both_all"] |= both
        lit = (rec["cos"] > 0) & (rec["shadow"] == 0)
        if useN[ip]:
            M["ng_obs"] |= rowN; M["ng_lit"][face[ng & lit]] = True
        if useR[ip]:
            M["rg_obs"] |= rowR; M["rg_lit"][face[rg & lit]] = True
        if useN[ip] and useR[ip]:
            M["both_obs"] |= both
            f = face[lit]
            M["both_lit"][f[both[f]]] = True
        n += len(rec)
    return {k: np.packbits(v) for k, v in M.items()}, n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--base", required=True)
    ap.add_argument("--regions", required=True, help="a contrib_phase regions.npz (row table: ipoint, obs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    t0 = time.time()
    base = Path(a.base)
    nf = (Path(base / "elements.bin").stat().st_size - 12) // (34 * 8)
    with np.load(a.regions, allow_pickle=False) as z:
        ipoint, obs_row, grid = z["ipoint"].astype(int), z["obs"].astype(bool), z["grid"].astype(bool)
        datN, datR = obs_row & (z["refng"] > 0), obs_row & (z["refrg"] > 0)
    useN, useR = np.zeros(ipoint.max() + 1, bool), np.zeros(ipoint.max() + 1, bool)
    useN[ipoint] = datN | grid
    useR[ipoint] = datR | grid
    if len(ipoint) != 1551 or obs_row.sum() != 1150 or datN.sum() != 934 or datR.sum() != 1150:
        raise SystemExit("unexpected row table: %d rows, %d observations" % (len(ipoint), obs_row.sum()))
    G.update(base=base, nf=nf, useN=useN, useR=useR)
    chunks = [list(c) for c in np.array_split(ipoint, a.procs * 4)]
    tot = {}
    nrec = 0
    with mp.get_context("fork").Pool(a.procs) as pool:
        for i, (res, n) in enumerate(pool.imap_unordered(work, chunks)):
            nrec += n
            for k, v in res.items():
                tot[k] = v if k not in tot else (tot[k] | v)
            print("chunk %d/%d records %d  %.0f s" % (i + 1, len(chunks), nrec, time.time() - t0), flush=True)
    np.savez(out, nfaces=nf, nrecords=nrec, nrows=len(ipoint), nobs=int(obs_row.sum()),
             definition="NG: Kn > 0; RG: KPhi > 0 or KV > 0; both: the same facet in the same row with both gauges having data", **tot)
    u = {k: int(np.unpackbits(v, count=nf).sum()) for k, v in tot.items()}
    print("done %.0f s: %d records; facets %s" % (time.time() - t0, nrec, u), flush=True)


if __name__ == "__main__":
    main()
