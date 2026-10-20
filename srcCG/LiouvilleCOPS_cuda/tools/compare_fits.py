#!/usr/bin/env python3
"""Compare two fit result directories written by fit_cops (CPU or GPU) (either interpreter).

Usage::

    python tools/compare_fits.py <fit dir A> <fit dir B> [--json out.json]

Reported, in the terms of the provisional end-to-end criteria (VALIDATION.md):

* **termination:** optimizer message, status, nfev and first-order optimality of each fit;
* **cost:** final cost and the relative difference |ΔC| / C_A;
* **parameters:** fitted values and their difference in units of fit A's Gauss-Newton σ,
  and whether each parameter sits on a bound. Bounds are taken from fit A's manifest
  in the fitter's scaled coordinates; activity coefficients are scaled by ``ascale``;
* **predictions:** max relative difference of the fitted NG/RG pressures
  (``COPS_model.dat``, all rows);
* **identities:** whether the selection, inputs, prior and optimizer settings in the two
  manifests are the same. If they differ, the comparison is not an equivalence test.
* ``compute_backend`` of each manifest (CPU when absent).

The tool only reads. With identical inputs and evaluators (for example one vs two GPUs)
it also states whether the parameter and prediction values are bitwise identical.
"""

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

PARAM = re.compile(r"^(\S+)\s+([+-]\d\.\d+e[+-]\d+)\s+([+-]\d\.\d+e[+-]\d+)\s+(\d\.\d+e[+-]\d+)\s*$")


def read_fit(directory):
    d = Path(directory)
    manifest = json.loads((d / "fit_manifest.json").read_text())
    params, scalars = {}, {}
    for line in (d / "fit_result.txt").read_text().splitlines():
        m = PARAM.match(line)
        if m:
            params[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("final_cost", "initial_cost", "final_data_cost"):
            scalars[parts[0]] = float(parts[1])
    table = {"pNG": [], "pRG": []}
    with (d / "COPS_model.dat").open() as stream:
        names = None
        for line in stream:
            if line.startswith("#iPoint"):
                names = line[1:].split()
                continue
            if line.startswith("#") or not line.strip() or names is None:
                continue
            row = line.split()
            table["pNG"].append(float(row[names.index("pNG")]))
            table["pRG"].append(float(row[names.index("pRG")]))
    return {"dir": str(d), "manifest": manifest, "params": params, "scalars": scalars,
            "pNG": np.array(table["pNG"]), "pRG": np.array(table["pRG"]),
            "optimizer": manifest.get("optimizer_results", [{}])[0]}


def on_bound(name, value, lo, hi, ascale):
    scale = ascale if name.startswith("a") else 1.0
    lo, hi = lo * scale if lo is not None else None, hi * scale if hi is not None else None
    tol = 1e-9 * max(1.0, abs(value))
    if lo is not None and abs(value - lo) <= tol:
        return "lower"
    if hi is not None and abs(value - hi) <= tol:
        return "upper"
    return ""


def compare(A, B):
    ma, mb = A["manifest"], B["manifest"]
    same_inputs = {key: ma.get(key) == mb.get(key)
                   for key in ("selection", "inputs", "prior", "optimizer_settings")}
    pa = ma["parameters"]
    names = pa["names"]
    lo = [None if v is None or isinstance(v, str) else float(v) for v in pa["lower_bounds"]]
    hi = [None if v is None or isinstance(v, str) else float(v) for v in pa["upper_bounds"]]
    # ascale = |initial a| / |initial scaled a| for any nonzero activity coefficient
    ascale = 1.0
    for k, n in enumerate(names):
        if n.startswith("a") and pa["initial_scaled"][k]:
            ascale = A["params"][n][0] / pa["initial_scaled"][k]
            break
    rows, max_sigma = [], 0.0
    for k, n in enumerate(names):
        if n not in A["params"] or n not in B["params"]:
            continue
        _, va, sa = A["params"][n]
        _, vb, _ = B["params"][n]
        dsig = abs(va - vb) / sa if sa > 0 else (0.0 if va == vb else float("inf"))
        max_sigma = max(max_sigma, dsig)
        rows.append({"name": n, "A": va, "B": vb, "sigma_A": sa, "diff_in_sigma_A": dsig,
                     "bound_A": on_bound(n, va, lo[k], hi[k], ascale),
                     "bound_B": on_bound(n, vb, lo[k], hi[k], ascale)})
    ca, cb = A["scalars"].get("final_cost"), B["scalars"].get("final_cost")
    pred = {}
    if A["pNG"].shape == B["pNG"].shape:
        for g in ("pNG", "pRG"):
            a, b = A[g], B[g]
            ok = (a > 0) & (b > 0)
            pred[g] = {"max_rel_diff": float(np.max(np.abs(b[ok] / a[ok] - 1.0))) if ok.any() else None,
                       "bitwise": bool(a.tobytes() == b.tobytes())}
    oa, ob = A["optimizer"], B["optimizer"]
    return {
        "A": A["dir"], "B": B["dir"],
        "backend_A": ma.get("compute_backend", {}).get("kind", "CPU (fit_cops.py)"),
        "backend_B": mb.get("compute_backend", {}).get("kind", "CPU (fit_cops.py)"),
        "same_inputs": same_inputs,
        "termination": {s: {"message": o.get("message"), "status": o.get("status"),
                            "nfev": o.get("nfev"), "optimality": o.get("optimality")}
                        for s, o in (("A", oa), ("B", ob))},
        "same_termination_class": (oa.get("success") == ob.get("success") and
                                   oa.get("status") == ob.get("status")),
        "final_cost": {"A": ca, "B": cb,
                       "rel_diff": abs(cb - ca) / ca if ca and cb is not None else None},
        "active_bounds_identical": all(r["bound_A"] == r["bound_B"] for r in rows),
        "max_parameter_diff_in_sigma_A": max_sigma,
        "parameters_bitwise": all(r["A"] == r["B"] for r in rows),
        "predictions": pred, "parameters": rows,
        "only_in_A": {n: A["params"][n][1] for n in A["params"] if n not in B["params"]},
        "only_in_B": {n: B["params"][n][1] for n in B["params"] if n not in A["params"]}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--json", help="also write the full comparison here (new file)")
    args = parser.parse_args()
    result = compare(read_fit(args.a), read_fit(args.b))
    print("A: %s\n   %s\nB: %s\n   %s" % (result["A"], result["backend_A"], result["B"], result["backend_B"]))
    print("same inputs: %s" % result["same_inputs"])
    for side in ("A", "B"):
        t = result["termination"][side]
        print("termination %s: %s (status %s, nfev %s, optimality %.3g)" %
              (side, t["message"], t["status"], t["nfev"], t["optimality"] or float("nan")))
    fc = result["final_cost"]
    print("final cost A %.9e  B %.9e  rel diff %.3e" % (fc["A"], fc["B"], fc["rel_diff"] or 0.0))
    print("active bounds identical: %s; max |dx| %.3g sigma_A; parameters bitwise: %s" %
          (result["active_bounds_identical"], result["max_parameter_diff_in_sigma_A"],
           result["parameters_bitwise"]))
    print("predictions: %s" % result["predictions"])
    if result["only_in_A"] or result["only_in_B"]:
        print("parameters only in A: %s; only in B: %s" % (result["only_in_A"], result["only_in_B"]))
    print("%-18s %14s %14s %10s %8s %6s %6s" % ("param", "A", "B", "sigma_A", "d/sig", "bndA", "bndB"))
    for r in result["parameters"]:
        print("%-18s %+14.6e %+14.6e %10.3e %8.3g %6s %6s" % (r["name"], r["A"], r["B"], r["sigma_A"],
                                                             r["diff_in_sigma_A"], r["bound_A"], r["bound_B"]))
    if args.json:
        with open(args.json, "x") as stream:
            json.dump(result, stream, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
