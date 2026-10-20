#!/usr/bin/env python3
"""C11: compare a GPU-built fitting base with the C++-built one, file by file (either interpreter).

Usage::

    python tools/c11_compare_bases.py <C++ base dir> <GPU base dir> [--json report.json]

Reported (the evidence for the C11 validation policy, which the user approves):

* ``elements.bin``: bitwise. It is written by the unchanged C++ code in both builds.
* per kernel file (point):
  * header: bitwise geometry (SPICE on the host in both builds);
  * record set: facets only in the C++ base, facets only in the GPU base, and their
    share of the point's total Kn / KPhi / KV (a visibility decision that flipped);
  * common records:
    * ``Kn``, ``KPhi``, ``KV``: max relative difference, and the sum of |Δ| over the point
      relative to the sum of |value|;
    * ``cos`` (CosSZA, host): bitwise count;
    * shadow flag mismatches (GPU cast-shadow ray).
* ``COPS_model.dat``: max relative difference of every numeric model column (pNG, pRG,
  densities, fluxes, nVisible, altitude), observation columns bitwise.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HEADER = np.dtype([("version", "<i4"), ("point", "<i4"), ("geometry", "<f8", (10,)), ("count", "<i4")])
RECORD = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])


def read_kernel(path):
    raw = Path(path).read_bytes()
    h = np.frombuffer(raw[:HEADER.itemsize], HEADER)[0]
    rec = np.frombuffer(raw[HEADER.itemsize:], RECORD)
    if len(rec) != int(h["count"]):
        raise SystemExit("bad record count in %s" % path)
    return h, rec


def read_table(path):
    names, rows = None, []
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
        elif not line.startswith("#") and line.strip():
            rows.append(line.split())
    return names, rows


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    scale = np.maximum(np.abs(a), np.abs(b))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(scale > 0, np.abs(a - b) / scale, 0.0)
    return float(r.max()) if r.size else 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("cpu")
    p.add_argument("gpu")
    p.add_argument("--json")
    a = p.parse_args()
    cpu, gpu = Path(a.cpu), Path(a.gpu)
    out = {"cpu_base": str(cpu), "gpu_base": str(gpu)}
    out["elements_bitwise"] = (cpu / "elements.bin").read_bytes() == (gpu / "elements.bin").read_bytes()

    files = sorted(p.name for p in (cpu / "kernels").glob("point=*.bin"))
    gfiles = sorted(p.name for p in (gpu / "kernels").glob("point=*.bin"))
    if files != gfiles:
        raise SystemExit("kernel file sets differ")
    agg = {"points": len(files), "records_cpu": 0, "records_gpu": 0, "header_bitwise": 0,
           "only_cpu": 0, "only_gpu": 0, "cos_bitwise_records": 0, "common_records": 0,
           "shadow_mismatch": 0, "max_rel": {"Kn": 0.0, "KPhi": 0.0, "KV": 0.0},
           "max_point_rel_sum": {"Kn": 0.0, "KPhi": 0.0, "KV": 0.0},
           "max_flipped_share": {"Kn": 0.0, "KPhi": 0.0, "KV": 0.0},
           "points_with_flips": 0, "worst": {}}
    for name in files:
        hc, rc = read_kernel(cpu / "kernels" / name)
        hg, rg = read_kernel(gpu / "kernels" / name)
        agg["records_cpu"] += len(rc)
        agg["records_gpu"] += len(rg)
        agg["header_bitwise"] += int(hc["geometry"].tobytes() == hg["geometry"].tobytes())
        fc, fg = rc["face"], rg["face"]
        common, ic, ig = np.intersect1d(fc, fg, assume_unique=True, return_indices=True)
        only_c = np.setdiff1d(fc, fg, assume_unique=True)
        only_g = np.setdiff1d(fg, fc, assume_unique=True)
        agg["only_cpu"] += len(only_c)
        agg["only_gpu"] += len(only_g)
        agg["common_records"] += len(common)
        if len(only_c) or len(only_g):
            agg["points_with_flips"] += 1
        c, g = rc[ic], rg[ig]
        agg["cos_bitwise_records"] += int(np.count_nonzero(c["cos"] == g["cos"]))
        agg["shadow_mismatch"] += int(np.count_nonzero(c["shadow"] != g["shadow"]))
        for k in ("Kn", "KPhi", "KV"):
            r = rel(c[k], g[k])
            if r > agg["max_rel"][k]:
                agg["max_rel"][k] = r
                agg["worst"][k] = name
            total = np.sum(np.abs(rc[k])) + np.sum(np.abs(rg[k]))
            if total > 0:
                agg["max_point_rel_sum"][k] = max(agg["max_point_rel_sum"][k],
                                                  float(np.sum(np.abs(c[k] - g[k])) / (0.5 * total)))
                flipped = (np.sum(np.abs(rc[k][np.isin(fc, only_c)])) +
                           np.sum(np.abs(rg[k][np.isin(fg, only_g)])))
                agg["max_flipped_share"][k] = max(agg["max_flipped_share"][k], float(flipped / (0.5 * total)))
    out["kernels"] = agg

    names, rc = read_table(cpu / "COPS_model.dat")
    gnames, rgt = read_table(gpu / "COPS_model.dat")
    if names != gnames or len(rc) != len(rgt):
        raise SystemExit("COPS_model.dat layouts differ")
    obs = {"iPoint", "UTC", "et", "x", "y", "z", "rh_AU", "refNG", "refRG", "T_RG", "useNG", "useRG", "inside", "VdotB_RG"}
    table = {}
    for j, n in enumerate(names):
        col_c = [r[j] for r in rc]
        col_g = [r[j] for r in rgt]
        if n in obs or n == "UTC":
            table[n] = {"identical": col_c == col_g}
        else:
            table[n] = {"max_rel": rel(np.array(col_c, float), np.array(col_g, float))}
    out["model_table"] = table

    print(json.dumps({"elements_bitwise": out["elements_bitwise"], "kernels": agg,
                      "model_table": {k: v for k, v in table.items()}}, indent=1))
    if a.json:
        with open(a.json, "x") as f:
            json.dump(out, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
