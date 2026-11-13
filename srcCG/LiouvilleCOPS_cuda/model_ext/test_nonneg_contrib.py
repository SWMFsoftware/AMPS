#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_nonneg_contrib.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_nonneg_contrib.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0575871f35b5a3d3afc95a8bb5b2fdf7554f192443b56e86ea2b50bfe8f14e04
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests of nonneg_contrib.py (run first by jobs/nonneg_contrib.sbatch).
# ==== end of header ====
"""Tests of nonneg_contrib.py: the W10x5 combination and the P% source-area ranking of obs_work against a brute-force
statement on synthetic data, tie order, and the gauge constants against the reference fitter's own expressions.
usage: .venv/bin/python model_ext/test_nonneg_contrib.py
"""
import sys
import types
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nonneg_contrib as C  # noqa: E402
import nonneg_maps as NM  # noqa: E402

RES = []


def check(name, ok, detail=""):
    RES.append(bool(ok))
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)


def test_obs_work():
    rng = np.random.default_rng(5)
    nf, ng = 2000, 6
    FAC = {}
    dense = []
    for j in range(ng):
        face = np.sort(rng.choice(nf, 700, replace=False)).astype(np.int32)
        c_ng = rng.exponential(1.0, 700) ** 3
        c_rg = rng.exponential(1.0, 700)
        c_rg[:5] = 0.25                                     # ties
        FAC[100 + j] = {"face": face, "c_ng": c_ng, "c_rg": c_rg, "ng": float(np.sum(c_ng)), "rg": float(np.sum(c_rg))}
        d = {"NG": np.zeros(nf), "RG": np.zeros(nf)}
        d["NG"][face], d["RG"][face] = c_ng, c_rg
        dense.append(d)
    W = sp.csr_matrix(np.array([[0.2, 0.3, 0.5, 0, 0, 0], [0, 0, 0.1, 0.4, 0.25, 0.25], [0, 0, 0, 0, 0, 1.0]]))
    area = rng.uniform(0.5, 2.0, nf)
    seen = rng.uniform(size=nf) < 0.5
    C._G.clear()
    C._G.update({"W": W, "FAC": FAC, "grid_rows": np.arange(100, 100 + ng), "nfaces": nf, "area": area, "seen": seen})
    for k in range(3):
        o = C.obs_work(k)
        for g in C.GAUGES:
            v = sum(W[k, j] * dense[j][g] for j in range(ng))
            f = np.flatnonzero(v > 0)
            order = sorted(f.tolist(), key=lambda e: (-v[e], e))      # decreasing, ties by facet index
            tot = sum(v[e] for e in order)
            ok = np.array_equal(o[g]["faces"], np.array(order)) and abs(o[g]["total_sorted"] / tot - 1) < 1e-13
            for i, p in enumerate(C.P_LEVELS):
                acc, n = 0.0, 0
                for e in order:
                    acc += v[e]
                    n += 1
                    if acc >= p / 100.0 * o[g]["total_sorted"]:
                        break
                ok &= o[g]["k"][i] == n and abs(o[g]["area"][i] - area[order[:n]].sum()) < 1e-9
                ok &= o[g]["cum_at_k"][i] >= p / 100.0 and (n == 1 or (np.cumsum(v[order])[n - 2] < p / 100.0 * tot))
            check("obs_work_%d_%s" % (k, g), ok, "k %s" % o[g]["k"])
        check("obs_work_%d_W_rowtotal" % k, abs(o["W_rows_total"]["NG"] / o["NG"]["total_face_order"] - 1) < 1e-13)


def test_constants():
    fc = NM.fit_cops_module()
    header = ["TemperatureNG=293 RGPressureModel=equilibrium_sphere SpeciesSensitivity=on (H2O 0.893, CO2 0.704)"]
    TRG = np.array([250.082, 260.5])
    bm, cNG, cRG, g = C.gauge_constants(fc, header, TRG)
    # independent statement of fit_cops.Model.prepare / crg (equilibrium_sphere): sqrt(2 pi m k T) / S
    ok = True
    for s, Ss in (("H2O", 0.893), ("CO2", 0.704)):
        m = fc.MASS[s]
        ok &= abs(cNG[s] / (fc.KB * 293.0 / Ss) - 1) < 1e-15
        ok &= np.allclose(cRG[s], np.sqrt(2.0 * np.pi * m * fc.KB * TRG) / Ss, rtol=1e-15, atol=0)
        ok &= abs(bm[s] / (np.sqrt(m / (2.0 * fc.KB)) * np.sqrt(np.pi) / 2.0) - 1) < 1e-15
    check("gauge_constants", ok, str(g))


def main():
    test_obs_work()
    test_constants()
    n = sum(1 for r in RES if not r)
    print("test_nonneg_contrib: %d checks, %d failed" % (len(RES), n))
    return 1 if n else 0


if __name__ == "__main__":
    sys.exit(main())
