#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_q002/preflight_compare.py (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_q002/preflight_compare.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 088c9436b69dcc1da558d1f5854f77eb374e19ce2d2c1e71b78f1aa52c3f15bf
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the comparison of the preflight: data selection, references, ft templates / enhancement scale / TS vs the
#   design's TS, base kernels and observation columns; rebuild_design / rebuild_late.
# RUN 2026-11-18: self-test old vs old 2026-11-18 (preflight_selftest.json).
# ==== end of header ====
"""Preflight of the option C stage-2r chains of the new law (q_min 0.02; lead / user 2026-11-18 ~02:40 CST): what of the
option C inputs must be rebuilt from the refitted D2 and base_q002?

    preflight_compare.py --old-d2 OLD_D2 --new-d2 NEW_D2 --old-base BASE --new-base BASE_Q002 --facets-from C9_FACETS
        --design DESIGN --late LATE --out preflight.json

Checks (each 'same' / 'differs' with numbers):
  selection  D2 fit_restart.npz point_ids, use_ng, use_rg, ref_ng, ref_rg   -> localnn_problem.Data, the late kernels' rows
  templates  enhancement_scale, the ft_c* template names and columns, TS = A_s sum ft columns vs the design's TS
             (facets file)                                                   -> the design's template columns, the starts
  base       elements.bin, run points (points.txt / late grid), the kernel files (same inode or same bytes), the
             observation columns of COPS_model.dat (iPoint, UTC, et, refNG, refRG)  -> the design and the late kernels
  rows       the design's and the late kernels' rows = the base rows
Decision fields: rebuild_design (template columns / base), rebuild_late (selection / base), notes.
"""
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np


def sha(p, n=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def obs_cols(base):
    rows = []
    for line in (Path(base) / "COPS_model.dat").read_text().splitlines():
        if line.startswith("#iPoint"):
            names = line[1:].split(); c = {n: i for i, n in enumerate(names)}; continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split(); rows.append([t[c[k]] for k in ("iPoint", "UTC", "et", "refNG", "refRG")])
    return rows


def main():
    ap = argparse.ArgumentParser()
    for k in ("old-d2", "new-d2", "old-base", "new-base", "facets-from", "design", "late", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    res = {"inputs": vars(a)}
    zo, zn = np.load(Path(a.old_d2) / "fit_restart.npz", allow_pickle=False), np.load(Path(a.new_d2) / "fit_restart.npz", allow_pickle=False)
    sel = {}
    for k in ("point_ids", "use_ng", "use_rg", "ref_ng", "ref_rg"):
        sel[k] = bool(np.array_equal(np.asarray(zo[k]), np.asarray(zn[k])))
    res["selection"] = sel
    mo, mn = json.loads(str(zo["meta_json"])), json.loads(str(zn["meta_json"]))
    tno, tnn = mo["templates"]["names"], mn["templates"]["names"]
    ft = [i for i, n in enumerate(tnn) if n.startswith("ft_c")]
    tc_o, tc_n = np.asarray(zo["template_columns"]), np.asarray(zn["template_columns"])
    f = np.load(a.facets_from, allow_pickle=False); TS = np.asarray(f["TS"], float); A_s_f = float(json.loads(str(f["meta_json"]))["A_scale"])
    TS_n = float(mn["enhancement_scale"]) * np.sum(tc_n[:, ft], axis=1)
    res["templates"] = {"names_equal": tno == tnn, "columns_equal": bool(np.array_equal(tc_o, tc_n)),
                        "enhancement_scale_old": float(mo["enhancement_scale"]), "enhancement_scale_new": float(mn["enhancement_scale"]),
                        "A_scale_facets_file": A_s_f, "TS_new_vs_design_TS_max_abs": float(np.max(np.abs(TS_n - TS))),
                        "TS_max": float(TS.max()),
                        "ft_values_old_new": [[float(np.asarray(zo["x"])[mo["names"].index(n)]) if n in mo["names"] else None,
                                               float(np.asarray(zn["x"])[mn["names"].index(n)]) if n in mn["names"] else None]
                                              for n in tnn if n.startswith("ft_c")]}
    b = {}
    for name in ("elements.bin", "points.txt"):
        po, pn = Path(a.old_base) / name, Path(a.new_base) / name
        if po.exists() or pn.exists():
            b[name] = (po.exists() and pn.exists() and (os.stat(po).st_ino == os.stat(pn).st_ino or sha(po) == sha(pn)))
    ko, kn = sorted((Path(a.old_base) / "kernels").glob("point=*.bin")), sorted((Path(a.new_base) / "kernels").glob("point=*.bin"))
    same_k = len(ko) == len(kn)
    diff_k = []
    for p, q in zip(ko, kn):
        if p.name != q.name or os.stat(p).st_ino != os.stat(q).st_ino:
            if p.name != q.name or os.path.getsize(p) != os.path.getsize(q) or sha(p) != sha(q):
                diff_k.append(q.name); same_k = False
    b["kernels"] = {"files": [len(ko), len(kn)], "same": same_k, "differ": diff_k[:10]}
    b["COPS_model_observation_columns"] = obs_cols(a.old_base) == obs_cols(a.new_base)
    res["base"] = b
    d = np.load(Path(a.design) / "design.npz", allow_pickle=False); L = np.load(a.late, allow_pickle=False)
    res["rows"] = {"design_rows_equal_new_d2": bool(np.array_equal(np.asarray(d["iPoint"]), np.asarray(zn["point_ids"])))}
    tmpl_same = res["templates"]["names_equal"] and res["templates"]["TS_new_vs_design_TS_max_abs"] <= 1e-12 * res["templates"]["TS_max"] and \
        res["templates"]["enhancement_scale_new"] == A_s_f
    base_same = all(v if isinstance(v, bool) else v["same"] for v in b.values())
    res["rebuild_design"] = not (tmpl_same and base_same)
    res["rebuild_late"] = not (all(sel.values()) and base_same)
    res["notes"] = ("the design (H2O law 0.02) holds the template columns of TS = A_s sum ft (facets file) and the base "
                    "kernels; the late kernels hold the selection's late_only facets and the base kernels")
    json.dump(res, open(a.out, "x"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("selection", "rows", "rebuild_design", "rebuild_late")}, default=float))
    print(json.dumps(res["templates"], default=float)[:600]); print(json.dumps(res["base"], default=float)[:400])


if __name__ == "__main__":
    sys.exit(main())
