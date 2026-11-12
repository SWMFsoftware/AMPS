#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/test_nonneg_maps.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/test_nonneg_maps.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c942f6b0be30e0a82e3f559c95efd8adea6b793b1e69db086c05b23d27feeeec
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tests. Unit tests of nonneg_maps.py (39): the law vs fit_cops.Model and C++ EvaluateElement, mesh order vs elements.bin,
#   VTK / Tecplot writers and read-back, state reproduction (A_final, seen, slope bins).
# USED BY: nonneg_maps.sbatch (every job runs them first).
# ==== end of header ====
"""Tests of nonneg_maps.py (activity and source-rate maps of the non-negative models): the law against the reference
fitter's source text and an independent statement of Source.cpp, the illumination coding, the VTU and Tecplot writers
(read back by this module and, when /usr/bin/pvpython exists, by ParaView's own readers), the final-state activity
checks on a small synthetic state, and the map times.

usage: .venv/bin/python model_ext/test_nonneg_maps.py [--scratch DIR] [--no-pvpython]
"""
import argparse
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nonneg_maps as M  # noqa: E402

PVPYTHON = "/usr/bin/pvpython"
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)


def tetra():
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int32)
    return nodes, faces


def test_law(fc):
    rec = M.law_guard(fc)
    check("law_guard", rec["AU_m"] == 149597870700.0, rec["fit_cops"])
    rng = np.random.default_rng(1)
    n = 100000
    cos = rng.uniform(-1.0, 1.0000001, n)
    cos[:10] = [0.0, -0.0, 1.0, 1.0000001, 1e-300, -1e-300, 0.5, 0.5, -0.5, 0.99]
    sh = ((rng.uniform(size=n) < 0.3) & (cos > 0)).astype(np.uint8)
    A = rng.uniform(0, 1e20, n)
    AC = rng.uniform(0, 1e17, n)
    kappa = 0.731
    for helio in ("shap5", "shap5_1"):
        for rh in (3.8310817 * fc.AU, 3.8334296 * fc.AU, 1.3 * fc.AU):
            Fh, Fc, night, qf = M.source_rates(cos, sh, rh, helio, "continuous", A, AC, kappa, fc.AU)
            # independent statement of src/Source.cpp GetQ + EvaluateElement (continuous, cos <= 1 not clipped there)
            p = 4.2143229 if helio == "shap5" else 4.8
            rhA = rh / 149597870700.0
            out_h = np.empty(n)
            out_c = np.empty(n)
            for i in range(0, n, 997):
                c, nt = cos[i], (cos[i] <= 0.0) or (sh[i] != 0)
                qmin, qmax = 10.0 * 0.02 / rhA ** p * 600, 1.0 / rhA ** p * 600
                out_h[i] = (qmin if nt else qmin + (qmax - qmin) * min(c, 1.0)) * A[i]
                qmin, qmax = 0.1 / rhA ** 2.0 * 600, 1.0 / rhA ** 2.0 * 600
                out_c[i] = kappa * ((qmin if nt else qmin + (qmax - qmin) * min(c, 1.0)) * AC[i])
            sel = np.arange(0, n, 997)
            rel = max(np.max(np.abs(Fh[sel] / out_h[sel] - 1)), np.max(np.abs(Fc[sel] / out_c[sel] - 1)))
            check("source_rates_vs_Source.cpp_%s_%.3f" % (helio, rh / fc.AU), rel < 1e-14, "max rel %.2e" % rel)
            check("night_rule_%s_%.3f" % (helio, rh / fc.AU), np.array_equal(night, (cos <= 0) | (sh != 0)))
            # night value = qmin A exactly, day at cos=1 = qmax A
            i0 = np.nonzero(night)[0][0]
            check("night_is_qmin", Fh[i0] == qf["H2O"][0] * A[i0])
    try:
        M.source_rates(cos, sh, 3.8 * fc.AU, "shap5", "srccg_table", A, AC, 1.0, fc.AU)
        check("refuse_non_continuous", False)
    except SystemExit:
        check("refuse_non_continuous", True)


def test_illumination():
    cos = np.array([0.5, 0.5, 0.0, -0.3, 1.0, -0.0])
    sh = np.array([0, 1, 0, 0, 1, 0], dtype=np.uint8)
    ill = M.illumination(cos, sh)
    check("illumination_codes", ill.tolist() == [0, 2, 1, 1, 2, 1], str(ill.tolist()))


def pv(script):
    r = subprocess.run([PVPYTHON, "-c", script], capture_output=True, text=True, timeout=300)
    return r.returncode, r.stdout, r.stderr


def test_vtu(d, use_pv):
    nodes, faces = tetra()
    f = np.array([1.5, -2.0, 3.25e20, 0.0])
    u8 = np.array([0, 1, 2, 1], dtype=np.uint8)
    p = d / "t.vtu"
    M.write_vtu(p, nodes, faces, [("F_H2O_per_m2_s", f), ("illumination", u8)], [("TimeValue", [3600.0]),
                                                                                  ("sun_direction", [0.1, 0.2, 0.3])],
                "comment with -- dashes and <tags> & units [m^-2 s^-1]")
    rd, fd, _ = M.read_vtu(p)
    ok = (np.array_equal(rd["Points"], nodes) and np.array_equal(rd["connectivity"], faces.reshape(-1)) and
          np.array_equal(rd["offsets"], [3, 6, 9, 12]) and np.array_equal(rd["types"], [5] * 4) and
          np.array_equal(rd["F_H2O_per_m2_s"], f) and np.array_equal(rd["illumination"], u8) and
          fd["TimeValue"] == [3600.0])
    check("vtu_roundtrip", ok)
    try:
        M.write_vtu(p, nodes, faces, [], [], "x")
        check("vtu_never_overwrites", False)
    except FileExistsError:
        check("vtu_never_overwrites", True)
    M.write_pvd(d / "t.pvd", [(0.0, "t.vtu"), (3600.0, "t.vtu")])
    if use_pv:
        rc, out, err = pv(
            "from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader as R\n"
            "r=R(); r.SetFileName(%r); r.Update(); g=r.GetOutput(); cd=g.GetCellData()\n"
            "a=cd.GetArray('F_H2O_per_m2_s'); b=cd.GetArray('illumination')\n"
            "print(g.GetNumberOfPoints(), g.GetNumberOfCells(), g.GetCellType(0), [a.GetValue(i) for i in range(4)],"
            " [b.GetValue(i) for i in range(4)], g.GetFieldData().GetArray('TimeValue').GetValue(0))\n"
            "from paraview import simple\n"
            "s=simple.PVDReader(FileName=%r); print('pvd times', list(s.TimestepValues))\n" % (str(p), str(d / "t.pvd")))
        exp = "4 4 5 [1.5, -2.0, 3.25e+20, 0.0] [0, 1, 2, 1] 3600.0"
        check("vtu_pvpython", rc == 0 and exp in out and "pvd times [0.0, 3600.0]" in out,
              (out.strip() + " " + err.strip()[-300:]))
        # a file larger than the XML parser's buffers (the size-dependent failure of appended data on this cluster)
        nb = 300000
        big = np.random.default_rng(3).integers(0, 4, (nb, 3)).astype(np.int32)
        vals = np.arange(nb, dtype=np.float64) * 1.5e15
        pb = d / "big.vtu"
        M.write_vtu(pb, nodes, big, [("F", vals), ("ill", (np.arange(nb) % 3).astype(np.uint8))], [("TimeValue", [1.0])],
                    "big")
        rc, out, err = pv(
            "from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader as R\n"
            "r=R(); r.SetFileName(%r); r.Update(); g=r.GetOutput(); a=g.GetCellData().GetArray('F')\n"
            "print('big', g.GetNumberOfCells(), a.GetValue(%d), g.GetCellData().GetArray('ill').GetValue(%d))\n" % (
                str(pb), nb - 1, nb - 1))
        exp = "big %d %r %d" % (nb, float(vals[-1]), (nb - 1) % 3)
        check("vtu_pvpython_large", rc == 0 and exp in out and "ERR" not in err, out.strip() + " " + err.strip()[-300:])


def test_tecplot(d, use_pv):
    nodes, faces = tetra()
    pool = mp.get_context("fork").Pool(2)
    try:
        tv = ["X [m]", "Y [m]", "Z [m]", "A_H2O", "seen", "cos_SZA", "illumination", "F_H2O [m^-2 s^-1]"]
        tdt = ["DOUBLE"] * 4 + ["BYTE", "DOUBLE", "BYTE", "DOUBLE"]
        p = d / "t.dat"
        w = M.TecplotWriter(p, pool, "test", tv, tdt, 2, "comment line 1\nunits [m^-2 s^-1]",
                            [("model", "softplus"), ("kappa", "3e-36")])
        A = np.array([1e20, 0.0, 2.5e18, 7.0])
        seen = np.array([1, 0, 1, 0], dtype=np.uint8)
        ref = {}
        for z, t in enumerate((0.0, 3600.0), start=1):
            cos = np.array([0.1, -0.5, 0.9, 1.0]) * (z if z == 1 else -1)
            ill = M.illumination(cos, np.zeros(4, np.uint8))
            F = np.abs(cos) * 1e20 + z
            vals = {"cos_SZA": cos, "illumination": ill, "F_H2O [m^-2 s^-1]": F}
            if z == 1:
                vals.update({"A_H2O": A, "seen": seen})
                w.zone("t1", t, nodes, 4, vals, conn=faces, zone_aux=[("UTC", "2016-09-30T00:00:00")])
                for c in range(3):
                    ref[(1, tv[c])] = nodes[:, c]
                ref[(1, "A_H2O")], ref[(1, "seen")] = A, seen
                ref[(1, "connectivity")] = faces.astype(np.int64) + 1
            else:
                w.zone("t2", t, nodes, 4, vals, zone_aux=[("UTC", "2016-09-30T01:00:00")])
            for v in tv[5:]:
                ref[(z, v)] = vals[v]
        w.close()
        res = M.verify_tecplot(p, w.blocks, ref, pool)
        check("tecplot_reparse", res["passed"] and res["values"] == 4 * 3 + 4 * 2 + 12 + 4 * 3 * 2,
              json.dumps(res, sort_keys=True))
        # a deliberately wrong reference is detected
        bad = dict(ref)
        bad[(2, "F_H2O [m^-2 s^-1]")] = ref[(2, "F_H2O [m^-2 s^-1]")] * (1 + 1e-8)
        check("tecplot_reparse_detects_errors", not M.verify_tecplot(p, w.blocks, bad, pool)["passed"])
    finally:
        pool.close()
        pool.join()
    print(p.read_text()[:1500])
    if use_pv:
        rc, out, err = pv(
            "from vtkmodules.vtkIOGeometry import vtkTecplotReader as R\n"
            "r=R(); r.SetFileName(%r); r.Update(); o=r.GetOutput()\n"
            "print('blocks', o.GetNumberOfBlocks())\n"
            "for k in range(o.GetNumberOfBlocks()):\n"
            "    b=o.GetBlock(k); cd=b.GetCellData()\n"
            "    print('block', k, b.GetNumberOfPoints(), b.GetNumberOfCells(), [cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())])\n"
            "    for nm in ('A_H2O','cos_SZA','F_H2O [m^-2 s^-1]'):\n"
            "        a=cd.GetArray(nm)\n"
            "        print(' ', nm, None if a is None else [a.GetValue(i) for i in range(a.GetNumberOfTuples())])\n" % str(p))
        print(out, err[-2000:])
        z1 = ("block 0 4 4" in out and "A_H2O [1.0000000200408773e+20, 0.0, 2.499999995126612e+18, 7.0]" in out)
        check("tecplot_pvpython_zone1", rc == 0 and "blocks 2" in out and z1,
              "(informational: ParaView's vtkTecplotReader reads zone 1 (in float32); it does not implement VARSHARELIST /"
              " CONNECTIVITYSHAREZONE, so zones 2.. are NOT usable in ParaView -- use the .vtu files there)")


def synthetic_state(d, model):
    """A small final state in the layout of the chain products (fit dir + ft dir)."""
    nm = M.load_module(HERE / "nonneg_mult.py", "nonneg_mult_test_%s" % model)
    rng = np.random.default_rng(7)
    n = 400
    fit, ft = d / ("%s_fit" % model), d / ("%s_ft" % model)
    fit.mkdir()
    ft.mkdir()
    scale = 2.2e18
    u = rng.normal(0, 1e18, n)
    TS = np.where(rng.uniform(size=n) < 0.2, rng.uniform(0, 1e17, n), 0.0)
    bins = rng.integers(0, 7, n).astype(np.int8)
    m = np.array([1.0, 0.0, 16.7, 1e-40, 0.3, 0.0, 3.48])
    A_orig = np.where(rng.uniform(size=n) < 0.05, 0.0, rng.uniform(1e16, 1e18, n))
    s = 1e-5 * A_orig if model == "softplus" else None
    held = np.where(rng.uniform(size=n) < 0.1, rng.normal(0, 3, n), 0.0)
    A_fit = nm.fields(model, u, TS, s, m[bins.astype(np.intp)], np.exp(held), scale, False)[0]
    delta = np.zeros(n)
    tun = np.nonzero(A_fit > 0)[0][:20]
    delta[tun] = rng.normal(0, 0.2, len(tun))
    A_final = A_fit.copy()
    A_final[tun] = A_fit[tun] * np.exp(delta[tun])
    seen = rng.uniform(size=n) < 0.6
    seen[tun] = True
    np.savez(fit / "nonneg_facets.npz", u=u, TS=TS, A=A_fit, P=A_fit, bin=bins, seen=seen, held=held)
    params = np.zeros(n, bool)
    params[tun] = True
    np.savez(ft / "finetune_state.npz", delta=delta, held_prev=held, held_next=held + delta, A_final=A_final,
             lower_bound=np.full(n, np.nan), seen=params, late=params)
    kappa = 3.07221091e-36
    np.savez(fit / "fit_restart.npz", kappa=np.array([kappa]), template_values=np.concatenate([m, np.ones(6)]))
    json.dump({"model_extension": {"nonneg": {"config": {"model": model, "s_ratio": 1e-5 if model == "softplus" else None},
                                              "model_check": {"A_scale": scale}}},
               "parameters": {"enhancement_scale": scale}}, open(fit / "fit_manifest.json", "w"))
    (fit / "fit_result.txt").write_text("a00    +1.0e+18 -5.5e+20 1.2e+20\nkappa  +1.000000e+00 +%.6e 2.665e-02\n"
                                        "kappa 0.000000\n" % kappa)
    (ft / "finetune_result.txt").write_text("facets late_only  WINDOW_MIN 5  window x; late facets 30, with A = 0 in "
                                            "the fit (not tuned) 10, parameters %d\ncomplete yes\n" % len(tun))
    json.dump({"fit": {"directory": str(fit)},
               "modules": {"nonneg_mult.py": {"path": str(HERE / "nonneg_mult.py"),
                                              "sha256": M.sha256_file(HERE / "nonneg_mult.py")},
                           "slope_mult.py": {"path": str(HERE / "slope_mult.py"),
                                             "sha256": M.sha256_file(HERE / "slope_mult.py")}}},
              open(ft / "fine_tuning.json", "w"))
    return fit, ft, {"AH2O": A_orig, "nfaces": n}, A_final, kappa, seen, bins, (fit, ft)


def test_state(d):
    for model in ("softplus", "squared"):
        fit, ft, el, A_final, kappa, seen, bins, _ = synthetic_state(d, model)
        S = M.state_activity(ft, fit, el, 1)
        c = S["checks"]
        check("state_%s_checks_pass" % model, c["passed"], json.dumps({k: v for k, v in c.items() if k not in (
            "nonneg_mult", "slope_mult")}, sort_keys=True))
        check("state_%s_fields" % model, np.array_equal(S["A"], A_final) and S["kappa"] == kappa and
              np.array_equal(S["seen"], seen) and np.array_equal(S["bins"], bins.astype(np.uint8)))
        # the maps' seen flag is the fit's (not the fine tuning's parameter mask)
        check("state_%s_seen_is_fit_seen" % model, int(S["seen"].sum()) == int(seen.sum()) and
              c["ft_parameter_facets"] == 20 and c["ft_parameters_subset_of_fit_seen"])
        # a corrupted A_final at an untuned facet must fail the checks
        z = dict(np.load(ft / "finetune_state.npz"))
        k = int(np.nonzero((z["delta"] == 0) & (z["A_final"] > 0))[0][0])
        z["A_final"][k] *= 1.0000001
        np.savez(ft / "finetune_state.npz", **z)
        check("state_%s_detects_corruption" % model, not M.state_activity(ft, fit, el, 1)["checks"]["passed"])


def test_times():
    check("map_times", len(M.MAP_TIMES) == 13 and M.MAP_TIMES[0] == "2016-09-30T00:00:00" and
          M.MAP_TIMES[10] == "2016-09-30T10:00:00" and M.MAP_TIMES[11:] == ("2016-09-30T10:33:10", "2016-09-30T10:39:26"))
    check("seconds_after_day0", M.seconds_after_day0("2016-09-30T10:39:26") == 38366.0 and
          M.seconds_after_day0("2016-09-30T10:33:10.282") == 37990.282 and M.seconds_after_day0(M.DAY0) == 0.0)
    check("nastran_real", M.nastran_real(b"-4.569067955E+02") == -456.9067955 and M.nastran_real(b" 1.5-3") == 1.5e-3
          and M.nastran_real(b"2.0+2") == 200.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch")
    ap.add_argument("--no-pvpython", action="store_true")
    a = ap.parse_args()
    use_pv = (not a.no_pvpython) and os.path.exists(PVPYTHON)
    d = Path(a.scratch) if a.scratch else Path(tempfile.mkdtemp(prefix="test_nonneg_maps_"))
    d.mkdir(parents=True, exist_ok=False) if a.scratch else None
    fc = M.fit_cops_module()
    test_law(fc)
    test_illumination()
    test_times()
    test_vtu(d, use_pv)
    test_tecplot(d, use_pv)
    test_state(d)
    nfail = sum(1 for r in RESULTS if not r[1] and "informational" not in r[2])
    print("test_nonneg_maps: %d checks, %d failed (scratch %s)" % (len(RESULTS), nfail, d))
    if not a.scratch:
        shutil.rmtree(d, ignore_errors=True)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
