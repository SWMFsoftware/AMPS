#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/clip_compare.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/clip_compare.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b5080e01c44bd3cf7f8502ec37ff4ed719292c9e2574152aa6025b02f25fd8c8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Compares two fits (arrays, optimizer, products) bitwise, and D1 vs D2 at the start. 2026-11-08.
# USED BY: jobs/clip_refit_bitwise.sbatch.
# ==== end of header ====
"""Comparisons for the clip refit tests (model_ext/clip_refit.py; 2026-11-08). Read only; numpy + standard library.

    clip_compare.py fits <fit dir A> <fit dir B> [--json OUT]
        bitwise comparison of two fits of the same problem (jobs/clip_refit_bitwise.sbatch: grain_fit_fixed.py vs
        clip_refit.py --variant D1 on the old base): fit_restart.npz arrays x, a, T, kappa, template_values,
        activity_linear, pred_ng, pred_rg (np.array_equal); optimizer cost, nfev, status, message; the bytes of
        COPS_model.dat, grain_term.txt, grain_components.dat, ActivityFit.txt, fit_result.txt (header lines with the
        output path excluded), pressure_residuals.dat; manifest optimizer_results cost/data_cost/D. Exit 0 iff all
        numerical items are identical.
    clip_compare.py start <D1 fit dir> <D2 fit dir> <point table> <run> [--json OUT]
        start_check_model.npz of D1 and D2 on the same base: the rows before the late window (bitwise), the late
        observations (D2 = W @ D2's grid rows, positive_problem's operator; D2 - D1 at the late rows).
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ARRAYS = ("x", "a", "T", "kappa", "template_values", "activity_linear", "pred_ng", "pred_rg", "point_ids")
FILES = ("COPS_model.dat", "grain_term.txt", "grain_components.dat", "ActivityFit.txt", "pressure_residuals.dat",
         "fit_result.txt")


def restart(path):
    with np.load(str(path), allow_pickle=False) as z:
        return {k: np.array(z[k]) for k in z.files if k != "meta_json"}, json.loads(str(z["meta_json"]))


def text_without_paths(path, dirs):
    lines = Path(path).read_text().splitlines(True)
    return "".join(l for l in lines if not any(str(d) in l for d in dirs))


def fits(a, b):
    a, b = Path(a).resolve(), Path(b).resolve()
    ra, ma = restart(a / "fit_restart.npz")
    rb, mb = restart(b / "fit_restart.npz")
    out = {"A": str(a), "B": str(b), "arrays": {}, "optimizer": {}, "files": {}, "manifest": {}}
    for k in ARRAYS:
        eq = bool(np.array_equal(ra[k], rb[k]))
        d = None
        if not eq and ra[k].dtype.kind == "f" and ra[k].shape == rb[k].shape:
            d = float(np.max(np.abs(ra[k] - rb[k]) / np.maximum(np.abs(ra[k]), 1e-300)))
        out["arrays"][k] = {"equal": eq, "max_rel_diff": d}
    for k in ("cost", "nfev", "status", "message", "optimality"):
        out["optimizer"][k] = {"A": ma["optimizer"].get(k), "B": mb["optimizer"].get(k),
                               "equal": ma["optimizer"].get(k) == mb["optimizer"].get(k)}
    for f in FILES:
        if not (a / f).is_file() or not (b / f).is_file():
            out["files"][f] = {"equal": None, "note": "missing"}
            continue
        ta, tb = text_without_paths(a / f, (a, b)), text_without_paths(b / f, (a, b))
        out["files"][f] = {"equal": ta == tb, "sha256_A": hashlib.sha256((a / f).read_bytes()).hexdigest(),
                           "sha256_B": hashlib.sha256((b / f).read_bytes()).hexdigest(),
                           "byte_identical": (a / f).read_bytes() == (b / f).read_bytes()}
    fa, fb = (json.loads((d / "fit_manifest.json").read_text())["optimizer_results"][0] for d in (a, b))
    for k in ("cost", "data_cost", "deviation_from_prior", "nfev", "message"):
        out["manifest"][k] = {"A": fa.get(k), "B": fb.get(k), "equal": fa.get(k) == fb.get(k)}
    numeric = ([v["equal"] for v in out["arrays"].values()] + [v["equal"] for v in out["optimizer"].values()] +
               [v["equal"] for v in out["files"].values() if v["equal"] is not None] +
               [v["equal"] for v in out["manifest"].values()])
    out["all_identical"] = bool(all(numeric))
    print("clip_compare fits: %s vs %s" % (a, b))
    for k, v in out["arrays"].items():
        print("  array %-16s equal %s%s" % (k, v["equal"], "" if v["equal"] else " (max rel diff %s)" % v["max_rel_diff"]))
    for k, v in out["optimizer"].items():
        print("  optimizer %-12s A %s | B %s  equal %s" % (k, v["A"], v["B"], v["equal"]))
    for k, v in out["files"].items():
        print("  file %-22s equal %s (byte-identical %s)" % (k, v["equal"], v.get("byte_identical")))
    for k, v in out["manifest"].items():
        print("  manifest %-20s A %s | B %s  equal %s" % (k, v["A"], v["B"], v["equal"]))
    print("  ALL IDENTICAL: %s" % out["all_identical"])
    return out


def start(d1, d2, table, run):
    import positive_problem as pp
    z1, z2 = np.load(Path(d1) / "start_check_model.npz"), np.load(Path(d2) / "start_check_model.npz")
    late_key = pp.time_key(pp.LATE_START_UTC)
    k1 = {int(k): j for j, k in enumerate(z1["time_key_ms"])}
    points = pp.read_run_points(run)
    ipoint = list(points)          # file (row) order of the base's COPS_model.dat
    lay = pp.LateLayout(pp.read_point_table(table), points, np.array(ipoint))
    pos2 = {int(r): j for j, r in enumerate(z2["rows"])}
    out = {"D1": str(Path(d1).resolve()), "D2": str(Path(d2).resolve())}
    early = [j for j, k in enumerate(z2["time_key_ms"]) if int(k) < late_key and int(k) in k1]
    e1 = [k1[int(z2["time_key_ms"][j])] for j in early]
    for g in ("pNG", "pRG"):
        s1, s2 = z1["surface_" + g], z2["surface_" + g]
        eq = bool(np.array_equal(s2[early], s1[e1]))
        d = float(np.max(np.abs(s2[early] - s1[e1]) / np.abs(s1[e1]))) if len(early) else 0.0
        lp = [pos2[int(r)] for r in lay.late_rows]
        gp = [pos2[int(r)] for r in lay.grid_rows]
        dw = float(np.max(np.abs(s2[lp] - lay.W @ s2[gp]) / np.abs(s2[lp])))
        l1 = [k1[int(z2["time_key_ms"][j])] for j in lp]
        dl = np.log(s2[lp]) - np.log(s1[l1])
        out[g] = {"early_rows": len(early), "early_bitwise_equal": eq, "early_max_rel_diff": d,
                  "late_rows": len(lp), "late_D2_minus_W_grid_max_rel": dw,
                  "late_ln_D2_over_D1": {"max_abs": float(np.max(np.abs(dl))), "mean": float(np.mean(dl)),
                                         "rms": float(np.sqrt(np.mean(dl * dl)))}}
        print("clip_compare start %s: %d rows before %s bitwise equal %s (max rel diff %.1e); %d late rows: "
              "|D2 - W @ grid|/D2 %.1e; ln(D2/D1) max %.3e mean %+.3e rms %.3e" % (
                  g, len(early), pp.LATE_START_UTC, eq, d, len(lp), dw, out[g]["late_ln_D2_over_D1"]["max_abs"],
                  out[g]["late_ln_D2_over_D1"]["mean"], out[g]["late_ln_D2_over_D1"]["rms"]))
    gr1, gr2 = z1["grain"], z2["grain"]
    out["grain_early_late_equal"] = bool(np.array_equal(gr2[[pos2[int(r)] for r in lay.late_rows]],
                                                        gr1[[k1[int(z2["time_key_ms"][pos2[int(r)]])]
                                                             for r in lay.late_rows]]))
    print("clip_compare start: grain at the late observations D1 == D2 (unaveraged, time stamp): %s"
          % out["grain_early_late_equal"])
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fits")
    f.add_argument("a")
    f.add_argument("b")
    f.add_argument("--json")
    s = sub.add_parser("start")
    s.add_argument("d1")
    s.add_argument("d2")
    s.add_argument("table")
    s.add_argument("run")
    s.add_argument("--json")
    a = ap.parse_args()
    out = fits(a.a, a.b) if a.cmd == "fits" else start(a.d1, a.d2, a.table, a.run)
    if a.json:
        with open(a.json, "x") as fh:
            json.dump(out, fh, indent=1)
    if a.cmd == "fits":
        return 0 if out["all_identical"] else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
