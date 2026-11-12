#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_maps.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_maps.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e4643ddfaa57d6c30e9f8fd1e4d0a1fd1cda0a3af35d0391ea5ce4829c4de065
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Activity and instantaneous source-rate maps of ONE final state of the softplus / squared models (user request
#   2026-11-11 "save the activity and source rate maps for visualization"; answers: instantaneous H2O and CO2 flux as
#   the fitter computes it (fit_cops.py: q = qmin at night or cast shadow, else qmin + (qmax - qmin) cos SZA, r_h),
#   00:00 ... 10:00, 10:33:10, 10:39:26 UTC, final state only, VTK .vtu and Tecplot ASCII). Subcommands: times, geometry
#   (builds and runs nonneg_maps_geom.cpp), validate-geometry (vs the base's kernel records), maps, verify.
# USED BY: jobs/nonneg_maps.sbatch (test 124322 on the 3-cycle states; production 124324 softplus, 124325 squared).
# ==== end of header ====
"""Activity and instantaneous source-rate maps of the non-negative H2O activity models (softplus and squared harmonics
x slope multipliers, model_ext/nonneg_*.py) of ONE final state, for visualization: VTK XML (.vtu, ParaView) and
Tecplot ASCII (FE triangle surface, cell-centred variables).

User request (2026-11-11): "save the activity and source rate maps for visualization". User answers: source rate = the
INSTANTANEOUS H2O and CO2 flux per facet exactly as the model computes it (cos SZA, the night and cast-shadow rule,
r_h) at 2016-09-30 00:00, 01:00, ..., 10:00, 10:33:10 and 10:39:26 UTC (13 times), together with the activity map;
for the FINAL state only (after the last fine tuning); both VTK (.vtu) and Tecplot ASCII.

Local, untracked model extension: nothing in the git tree is edited; the reference fitter is imported unchanged
(cops_cuda.reference.load("fit_cops")) for its constants and readers.

THE LAW (read from the code; law_guard() checks that these expressions are still the reference's, verbatim)
  ../LiouvilleCOPS/tools/fit_cops.py Model.load_kernels / prepare (= src/Source.cpp cSourceModel::GetQ /
  EvaluateElement, = cops_cuda/kernels/pressure.cu), per facet e at time t:
    rhA = r_h / AU (r_h = |Sun - nucleus| of the base's SPICE geometry [m], AU = 149597870700 m)
    p = 4.2143229 (HeliocentricModel shap5 of the base) else 4.8
    H2O: qmin = 0.2 * 600 / rhA**p, qmax = 600 / rhA**p;  CO2: qmin = 0.1 * 600 / rhA**2, qmax = 600 / rhA**2
    night = (cos SZA <= 0) | cast shadow                      (the kernel records' fields cos and shadow)
    q_s = qmin_s if night else qmin_s + (qmax_s - qmin_s) * min(cos SZA, 1)        (SZAMode continuous)
    F_H2O = q_H2O * A_H2O          (A_H2O = the model's activity, af = max(A, 0) = A for these models)
    F_CO2 = kappa * (q_CO2 * A_CO2) (A_CO2 = elements.bin ACO2 = max(0, srcCG shap5 CO2 coefficients . Y), kappa fitted)
  Units: F in m^-2 s^-1 (src/Source.h "emitted molecular flux [m^-2 s^-1]"; COPS_model.dat "units: SI (... m^-2 s^-1
  ...)"; Q_s = sum_e F_s,e area_e in s^-1). The code assigns no separate units to q and A; only their product F.
  Non-negative models (nonneg_mult.fields): A = m_b(e) q(L_e) exp(h_e); the fine tuning's final A is
  finetune_state.npz A_final (= A_fit at untuned facets, A_fit exp(delta) at the tuned ones); checked here against
  nonneg_mult.fields(..., exp(held_next)) with the fit's u, TS, bins, multipliers and s.

GEOMETRY (nonneg_maps_geom.cpp, the base builder's own method; validated by validate-geometry against the kernel
  records and the Q_H2O / Q_CO2 columns of the base at 17 observation / grid points).

Subcommands
  times              the times file of nonneg_maps_geom: 13 map times (t00..t12) + 17 validation rows (pNNNN)
  validate-geometry  nonneg_maps_geom output vs the base (mesh, BDF, Sun, cos, shadow, production)
  maps               the maps of one final state (--ft-dir c<N>/ft/final, --fit-dir c<N>/fit/final) in --out (new)
"""

import argparse
import base64
import datetime as dt
import hashlib
import inspect
import importlib.util
import json
import math
import multiprocessing as mp
import os
import re
import socket
import struct
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
for _p in (str(APP),):
    if _p not in sys.path:
        sys.path.insert(0, _p)

SCHEMA = "nonneg_maps_v1"
DAY0 = "2016-09-30T00:00:00"
MAP_TIMES = tuple(["2016-09-30T%02d:00:00" % h for h in range(11)] + ["2016-09-30T10:33:10", "2016-09-30T10:39:26"])
MAP_LABELS = tuple("t%02d" % k for k in range(len(MAP_TIMES)))
# base rows for the geometry validation: the observation points 4 s after each map hour, the exact map times
# 10:33:10 (row 1079) and 10:39:26 (row 1549), and late-window grid points
VALIDATION_ROWS = (0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 1079, 1142, 1292, 1442, 1517, 1549)
ILLUMINATION = ((0, "day", "cos SZA > 0 and the Sun ray from the facet is not blocked"),
                (1, "night", "cos SZA <= 0"),
                (2, "cast_shadow", "cos SZA > 0 and the Sun ray is blocked by the nucleus (night-side source)"))
KREC = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])
FLOAT_FMT = "%.9e"          # Tecplot ASCII: 10 significant digits (relative rounding <= 5e-10); the .vtu files are exact
PER_LINE = 8
CHUNK = 1 << 20             # values per formatting task (multiple of PER_LINE)
TEC_REL_TOL = 5.0000001e-10
LOG = "nonneg_maps:"


def say(text):
    print("%s [%s] %s" % (LOG, dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ"), text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def seconds_after_day0(utc):
    """Seconds after 2016-09-30T00:00:00 UTC of an ISO UTC string (no leap second inside this day)."""
    t = dt.datetime.strptime(utc[:19], "%Y-%m-%dT%H:%M:%S")
    frac = float("0" + utc[19:]) if len(utc) > 19 else 0.0
    return (t - dt.datetime(2016, 9, 30)).total_seconds() + frac


def fit_cops_module():
    from cops_cuda import reference
    return reference.load("fit_cops")


# ------------------------------------------------------------------------------------------------ the law
LAW_SOURCE_LINES = {
    "load_kernels": (
        'cs.append(rec["cos"]), night.append((rec["cos"] <= 0.0) | (rec["shadow"] != 0))',
        'rhA = self.rh / AU',
        'p = 4.2143229 if self.helio == "shap5" else 4.8',
        'c = np.minimum(self.cos, 1.0)',
        'if self.SZAMode == "continuous":',
        'cz = c',
        'for s, (qmin, qmax) in (("H2O", (0.2 * 600 / rhA ** p, 600 / rhA ** p)), ("CO2", (0.1 * 600 / rhA ** 2, '
        '600 / rhA ** 2))):',
        'self.Q[s] = np.where(self.night, qmin, qmin + (qmax - qmin) * cz)'),
    "prepare": (
        'ACO2 = self.ACO2[self.face]',
        'q = {"H2O": self.Q["H2O"], "CO2": self.Q["CO2"] * ACO2}',
        'pNG = pNG_H + kappa * pNG_C'),
}


def law_guard(fit_cops):
    """The source law of the reference fitter, verbatim: refuse to make maps if any expression changed."""
    missing = []
    src = {"load_kernels": inspect.getsource(fit_cops.Model.load_kernels),
           "prepare": inspect.getsource(fit_cops.Model.prepare) + inspect.getsource(fit_cops.Model.evaluate)}
    for key, lines in LAW_SOURCE_LINES.items():
        for line in lines:
            if line not in src[key]:
                missing.append("%s: %s" % (key, line))
    if float(fit_cops.AU) != 149597870700.0:
        missing.append("AU = %r" % fit_cops.AU)
    if missing:
        raise SystemExit("nonneg_maps: the reference source law changed; not making maps:\n  " + "\n  ".join(missing))
    return {"checked": "fit_cops.Model.load_kernels / prepare / evaluate expressions verbatim", "lines": LAW_SOURCE_LINES,
            "fit_cops": str(Path(fit_cops.__file__).resolve()), "fit_cops_sha256": sha256_file(fit_cops.__file__),
            "AU_m": float(fit_cops.AU)}


def q_factors(rh_m, helio, AU):
    """(qmin, qmax) of H2O and CO2 at r_h [m]: fit_cops.Model.load_kernels' expressions."""
    rhA = np.float64(rh_m) / AU
    p = 4.2143229 if helio == "shap5" else 4.8
    return {"H2O": (0.2 * 600 / rhA ** p, 600 / rhA ** p), "CO2": (0.1 * 600 / rhA ** 2, 600 / rhA ** 2),
            "p": p, "rhA": float(rhA)}


def source_rates(cos, shadow, rh_m, helio, sza_mode, A_h2o, A_co2, kappa, AU):
    """Instantaneous emitted fluxes per facet [m^-2 s^-1] exactly as the fitter forms them, and the night flag."""
    if sza_mode != "continuous":
        raise SystemExit("nonneg_maps: SZAMode %r; the GPU fits support continuous only" % sza_mode)
    qf = q_factors(rh_m, helio, AU)
    night = (cos <= 0.0) | (shadow != 0)
    cz = np.minimum(cos, 1.0)
    q = {s: np.where(night, qf[s][0], qf[s][0] + (qf[s][1] - qf[s][0]) * cz) for s in ("H2O", "CO2")}
    F_h2o = q["H2O"] * A_h2o
    F_co2 = kappa * (q["CO2"] * A_co2)
    return F_h2o, F_co2, night, qf


def illumination(cos, shadow):
    """0 day, 1 night (cos <= 0), 2 cast shadow (cos > 0 and shadowed); uint8."""
    out = np.zeros(len(cos), dtype=np.uint8)
    out[cos <= 0.0] = 1
    out[(cos > 0.0) & (shadow != 0)] = 2
    return out


# ------------------------------------------------------------------------------------------------ inputs
def read_elements(path):
    """elements.bin (version 1): centroid, normal, area, AH2O, ACO2 (copies; the 25 Y columns are dropped)."""
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        version, n, nb = struct.unpack("<iii", f.read(12))
        if version != 1 or nb != 25 or size != 12 + n * 34 * 8:
            raise SystemExit("unsupported element file %s" % path)
        d = np.fromfile(f, dtype="<f8", count=n * 34).reshape(n, 34)
    out = {"centroid": np.ascontiguousarray(d[:, 0:3]), "normal": np.ascontiguousarray(d[:, 3:6]),
           "area": np.ascontiguousarray(d[:, 6]), "AH2O": np.ascontiguousarray(d[:, 32]),
           "ACO2": np.ascontiguousarray(d[:, 33]), "nfaces": n}
    del d
    return out


def read_kernel_full(path):
    with open(path, "rb") as f:
        raw = f.read()
    version, ip = struct.unpack_from("<ii", raw, 0)
    h = struct.unpack_from("<10d", raw, 8)
    (nk,) = struct.unpack_from("<i", raw, 88)
    if version != 1 or len(raw) != 92 + nk * KREC.itemsize:
        raise SystemExit("unexpected kernel file %s" % path)
    rec = np.frombuffer(raw, dtype=KREC, offset=92, count=nk)
    return ip, {"et": h[0], "VdotB": h[1], "sc": h[2:5], "sun": h[5:8], "rh": h[8], "alt": h[9]}, rec


def base_table(base):
    fc = fit_cops_module()
    header, names, rows = fc.read_model_file(str(Path(base) / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    return fc, header, cols, rows


def read_geometry(gdir):
    gdir = Path(gdir)
    times = {}
    head = []
    with open(gdir / "times.txt") as f:
        for line in f:
            if line.startswith("#"):
                head.append(line.rstrip("\n"))
                continue
            t = line.split()
            if not t:
                continue
            times[t[0]] = {"label": t[0], "utc": t[1], "et": float(t[2]), "sun": np.array([float(v) for v in t[3:6]]),
                           "rh_m": float(t[6]), "n_cos_pos": int(t[7]), "n_shadow": int(t[8]),
                           "Q_H2O_base": float(t[9]), "Q_CO2_base": float(t[10])}
    fg = np.fromfile(gdir / "mesh_facegeom.f8", dtype="<f8")
    nf = len(fg) // 8
    fg = fg.reshape(nf, 8)
    nodes = np.fromfile(gdir / "mesh_nodes.f8", dtype="<f8").reshape(-1, 3)
    faces = np.fromfile(gdir / "mesh_faces.i4", dtype="<i4").reshape(-1, 3)
    if len(faces) != nf:
        raise SystemExit("geometry %s: %d faces vs %d face records" % (gdir, len(faces), nf))
    return {"dir": gdir, "times": times, "header": head, "facegeom": fg, "nodes": nodes, "faces": faces, "nfaces": nf}


def geometry_time(gdir, label, nf):
    cos = np.fromfile(Path(gdir) / ("%s.cos.f8" % label), dtype="<f8")
    sh = np.fromfile(Path(gdir) / ("%s.shadow.u1" % label), dtype=np.uint8)
    if len(cos) != nf or len(sh) != nf:
        raise SystemExit("geometry %s: %s has %d / %d values, expected %d" % (gdir, label, len(cos), len(sh), nf))
    return cos, sh


# ------------------------------------------------------------------------------------------------ times
def cmd_times(args):
    _, _, cols, rows = base_table(args.base)
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    lines = ["# nonneg_maps times file (%s): label UTC; tNN = map times, pNNNN = base rows (iPoint) for the geometry "
             "validation (UTC string of the base's COPS_model.dat)" % utc_now()]
    for lab, utc in zip(MAP_LABELS, MAP_TIMES):
        lines.append("%s %s" % (lab, utc))
    for r in VALIDATION_ROWS:
        row = rows[r]
        lines.append("p%04d %s" % (int(row[cols["iPoint"]]), row[cols["UTC"]]))
    with open(out, "x") as f:
        f.write("\n".join(lines) + "\n")
    say("wrote %s (%d map times, %d validation points)" % (out, len(MAP_TIMES), len(VALIDATION_ROWS)))


# ------------------------------------------------------------------------------------------------ BDF (independent)
def nastran_real(s):
    try:
        return float(s)
    except ValueError:
        t = s.replace(b" ", b"")
        for i in range(1, len(t)):
            if t[i:i + 1] in (b"-", b"+") and t[i - 1:i] not in (b"E", b"e", b"D", b"d"):
                t = t[:i] + b"E" + t[i:]
                break
        return float(t.replace(b"D", b"E").replace(b"d", b"E"))


def parse_bdf(path):
    """GRID / GRID* nodes and CTRIA3 triangles of a NASTRAN bulk data file (independent of Mesh.cpp)."""
    with open(path, "rb") as f:
        lines = f.read().split(b"\n")
    ids, xyz, tri = [], [], []
    i, n = 0, len(lines)
    while i < n:
        ln = lines[i]
        if ln.startswith(b"GRID*"):
            nxt = lines[i + 1]
            ids.append(int(ln[8:24]))
            xyz.append((nastran_real(ln[40:56]), nastran_real(ln[56:72]), nastran_real(nxt[8:24])))
            i += 2
            continue
        if ln.startswith(b"GRID"):
            ids.append(int(ln[8:16]))
            xyz.append((nastran_real(ln[24:32]), nastran_real(ln[32:40]), nastran_real(ln[40:48])))
        elif ln.startswith(b"CTRIA3"):
            t = ln.split()
            if len(t) >= 6:
                tri.append((int(t[3]), int(t[4]), int(t[5])))
            else:
                tri.append(tuple(int(ln[8 + 8 * k:16 + 8 * k]) for k in (2, 3, 4)))
        i += 1
    ids = np.array(ids, dtype=np.int64)
    xyz = np.array(xyz, dtype=np.float64)
    tri = np.array(tri, dtype=np.int64)
    order = np.argsort(ids, kind="stable")
    sid = ids[order]
    if np.any(np.diff(sid) == 0):
        raise SystemExit("BDF %s: duplicate GRID ids" % path)
    pos = np.searchsorted(sid, tri)
    if np.any(pos >= len(sid)) or np.any(sid[np.minimum(pos, len(sid) - 1)] != tri):
        raise SystemExit("BDF %s: CTRIA3 references an undefined GRID" % path)
    return xyz[order], pos


# ------------------------------------------------------------------------------------------------ validate-geometry
def cmd_validate_geometry(args):
    t0 = time.time()
    G = read_geometry(args.geometry)
    out_json = Path(args.geometry) / "geometry_validation.json"
    if out_json.exists():
        raise SystemExit("%s exists" % out_json)
    base = Path(args.base)
    fc, header, cols, rows = base_table(base)
    by_ip = {int(r[cols["iPoint"]]): r for r in rows}
    el = read_elements(base / "elements.bin")
    nf = el["nfaces"]
    rep = {"schema": SCHEMA + ".geometry_validation", "created_utc": utc_now(), "host": socket.gethostname(),
           "geometry": str(Path(args.geometry).resolve()), "base": str(base.resolve()), "checks": {}, "points": []}
    ok = True
    # 1. mesh of the driver vs elements.bin (same code: bitwise)
    fg = G["facegeom"]
    m = {"nfaces_geometry": int(G["nfaces"]), "nfaces_elements": int(nf)}
    if G["nfaces"] != nf:
        raise SystemExit("face counts differ: geometry %d, elements.bin %d" % (G["nfaces"], nf))
    for name, sl in (("centroid", slice(0, 3)), ("normal", slice(3, 6)), ("area", slice(6, 7))):
        ref = el[name].reshape(nf, -1)
        d = fg[:, sl] - ref
        m[name + "_bitwise_equal_faces"] = int(np.sum(np.all(fg[:, sl] == ref, axis=1)))
        m[name + "_max_abs_diff"] = float(np.max(np.abs(d)))
    m["passed"] = all(m[k + "_bitwise_equal_faces"] == nf for k in ("centroid", "normal", "area"))
    ok &= m["passed"]
    rep["checks"]["1_mesh_vs_elements_bin"] = m
    say("check 1 mesh vs elements.bin: centroid %d, normal %d, area %d of %d faces bitwise; passed %s" % (
        m["centroid_bitwise_equal_faces"], m["normal_bitwise_equal_faces"], m["area_bitwise_equal_faces"], nf,
        m["passed"]))
    # 2. independent BDF parse: facet order and mesh (centroid, area) vs elements.bin
    bdf = None
    for line in G["header"]:
        mm = re.search(r"shape (\S+) \(", line)
        if mm:
            bdf = mm.group(1)
    manifest = json.load(open(base / "run_manifest.json"))
    shape_dep = [d for d in manifest["dependencies"] if d["role"] == "shape_model"][0]["file"]
    b = {"bdf": bdf, "manifest_shape_model": shape_dep["path"]}
    if bdf != shape_dep["path"]:
        raise SystemExit("geometry shape %s is not the base's %s" % (bdf, shape_dep["path"]))
    if not args.skip_bdf_hash:
        b["sha256"] = sha256_file(bdf)
        b["sha256_matches_manifest"] = b["sha256"] == shape_dep["sha256"]
    t1 = time.time()
    xyz, tri = parse_bdf(bdf)
    b["parse_seconds"] = time.time() - t1
    b["grid_nodes"], b["ctria3"] = int(len(xyz)), int(len(tri))
    if len(tri) != nf:
        raise SystemExit("BDF has %d CTRIA3, elements.bin %d faces" % (len(tri), nf))
    a_, b_, c_ = xyz[tri[:, 0]], xyz[tri[:, 1]], xyz[tri[:, 2]]
    corners_drv = G["nodes"][G["faces"]]
    b["corner_coordinates_bitwise_equal_driver_faces"] = int(np.sum(np.all(
        (corners_drv[:, 0] == a_) & (corners_drv[:, 1] == b_) & (corners_drv[:, 2] == c_), axis=1)))
    cen = (a_ + b_ + c_) / 3.0
    cr = np.cross(b_ - a_, c_ - a_)
    area = 0.5 * np.sqrt(np.einsum("ij,ij->i", cr, cr))
    b["centroid_bitwise_equal_faces"] = int(np.sum(np.all(cen == el["centroid"], axis=1)))
    b["centroid_max_abs_diff_m"] = float(np.max(np.abs(cen - el["centroid"])))
    b["area_max_rel_diff"] = float(np.max(np.abs(area / el["area"] - 1.0)))
    b["normal_max_abs_diff"] = float(np.max(np.abs(cr / (2.0 * area)[:, None] - el["normal"])))
    b["passed"] = (b["corner_coordinates_bitwise_equal_driver_faces"] == nf and b["area_max_rel_diff"] < 1e-12
                   and b["centroid_max_abs_diff_m"] < 1e-9 and b.get("sha256_matches_manifest", True))
    ok &= b["passed"]
    rep["checks"]["2_bdf_independent"] = b
    del a_, b_, c_, cen, cr, area, corners_drv, xyz, tri
    say("check 2 BDF (independent parse, %.0f s): %d CTRIA3, %d corners bitwise = driver mesh; centroid %d bitwise, max "
        "|d| %.3e m; area max rel %.3e; passed %s" % (b["parse_seconds"], b["ctria3"],
                                                       b["corner_coordinates_bitwise_equal_driver_faces"],
                                                       b["centroid_bitwise_equal_faces"], b["centroid_max_abs_diff_m"],
                                                       b["area_max_rel_diff"], b["passed"]))
    # 3. per validation point: Sun, cos, shadow vs the kernel records; production vs the base table
    helio = fc.header_value(header, "HeliocentricModel")
    sza = fc.header_value(header, "SZAMode")
    tot = {"points": 0, "records": 0, "cos_equal": 0, "shadow_mismatch": 0, "night_mismatch": 0, "max_abs_dcos": 0.0}
    for lab in sorted(k for k in G["times"] if k.startswith("p")):
        T = G["times"][lab]
        ip = int(lab[1:])
        row = by_ip[ip]
        kp, h, rec = read_kernel_full(base / "kernels" / ("point=%04d.bin" % ip))
        cos, sh = geometry_time(args.geometry, lab, nf)
        face = rec["face"].astype(np.int64)
        dc = np.abs(cos[face] - rec["cos"])
        shm = (sh[face] != 0) != (rec["shadow"] != 0)
        nm_ = ((cos[face] <= 0) | (sh[face] != 0)) != ((rec["cos"] <= 0) | (rec["shadow"] != 0))
        F_h, F_c, night, qf = source_rates(cos, sh, T["rh_m"], helio, sza, el["AH2O"], el["ACO2"], 1.0, fc.AU)
        Qh_py, Qc_py = float(np.sum(F_h * el["area"])), float(np.sum(F_c * el["area"]))
        Qh_tab, Qc_tab = float(row[cols["Q_H2O"]]), float(row[cols["Q_CO2"]])
        pt = {"label": lab, "iPoint": ip, "utc": T["utc"], "records": int(len(rec)),
              "et_equal": T["et"] == h["et"], "sun_equal": bool(np.all(T["sun"] == np.array(h["sun"]))),
              "sun_max_abs_diff": float(np.max(np.abs(T["sun"] - np.array(h["sun"])))),
              "rh_equal": T["rh_m"] == h["rh"], "rh_rel_diff": float(abs(T["rh_m"] / h["rh"] - 1.0)),
              "cos_equal_records": int(np.sum(dc == 0.0)), "max_abs_dcos": float(np.max(dc)) if len(dc) else 0.0,
              "shadow_records_kernel": int(np.sum(rec["shadow"] != 0)), "shadow_records_geometry": int(np.sum(sh[face])),
              "shadow_mismatch": int(np.sum(shm)), "night_mismatch": int(np.sum(nm_)),
              "facets_cos_pos_all": int(np.sum(cos > 0)), "facets_shadow_all": int(np.sum(sh != 0)),
              "Q_H2O_table": Qh_tab, "Q_CO2_table": Qc_tab, "Q_H2O_driver": T["Q_H2O_base"],
              "Q_CO2_driver": T["Q_CO2_base"], "Q_H2O_python": Qh_py, "Q_CO2_python": Qc_py,
              "Q_H2O_rel_driver_vs_table": abs(T["Q_H2O_base"] / Qh_tab - 1.0),
              "Q_CO2_rel_driver_vs_table": abs(T["Q_CO2_base"] / Qc_tab - 1.0),
              "Q_H2O_rel_python_vs_driver": abs(Qh_py / T["Q_H2O_base"] - 1.0),
              "Q_CO2_rel_python_vs_driver": abs(Qc_py / T["Q_CO2_base"] - 1.0),
              "rh_AU_table": float(row[cols["rh_AU"]])}
        if pt["shadow_mismatch"]:
            idx = np.nonzero(shm)[0][:50]
            pt["shadow_mismatch_examples"] = [{"face": int(face[i]), "cos_kernel": float(rec["cos"][i]),
                                               "cos_geometry": float(cos[face[i]]),
                                               "shadow_kernel": int(rec["shadow"][i]),
                                               "shadow_geometry": int(sh[face[i]])} for i in idx]
        rep["points"].append(pt)
        tot["points"] += 1
        tot["records"] += pt["records"]
        tot["cos_equal"] += pt["cos_equal_records"]
        tot["shadow_mismatch"] += pt["shadow_mismatch"]
        tot["night_mismatch"] += pt["night_mismatch"]
        tot["max_abs_dcos"] = max(tot["max_abs_dcos"], pt["max_abs_dcos"])
        say("  %s %s: %d records, et/sun/rh bitwise %s/%s/%s, cos equal %d (max |d| %.2e), shadow mismatches %d "
            "(kernel %d shadowed), night mismatches %d; Q_H2O driver/table-1 %.1e, python/driver-1 %.1e; Q_CO2 %.1e, %.1e" % (
                lab, T["utc"], pt["records"], pt["et_equal"], pt["sun_equal"], pt["rh_equal"], pt["cos_equal_records"],
                pt["max_abs_dcos"], pt["shadow_mismatch"], pt["shadow_records_kernel"], pt["night_mismatch"],
                pt["Q_H2O_rel_driver_vs_table"], pt["Q_H2O_rel_python_vs_driver"], pt["Q_CO2_rel_driver_vs_table"],
                pt["Q_CO2_rel_python_vs_driver"]))
    P = rep["points"]
    tot["sun_rh_et_bitwise_all"] = all(p["et_equal"] and p["sun_equal"] and p["rh_equal"] for p in P)
    tot["max_Q_rel_driver_vs_table"] = max(max(p["Q_H2O_rel_driver_vs_table"], p["Q_CO2_rel_driver_vs_table"]) for p in P)
    tot["max_Q_rel_python_vs_driver"] = max(max(p["Q_H2O_rel_python_vs_driver"], p["Q_CO2_rel_python_vs_driver"])
                                            for p in P)
    tot["passed"] = (tot["sun_rh_et_bitwise_all"] and tot["shadow_mismatch"] == 0 and tot["night_mismatch"] == 0
                     and tot["max_abs_dcos"] <= 1e-15 and tot["max_Q_rel_driver_vs_table"] < 1e-6
                     and tot["max_Q_rel_python_vs_driver"] < 1e-12)
    tot["criteria"] = ("et, Sun direction, r_h bitwise; |dcos| <= 1e-15; 0 shadow and 0 night mismatches; Q (table, "
                       "7 printed digits) rel < 1e-6; python law vs C++ EvaluateElement rel < 1e-12")
    ok &= tot["passed"]
    rep["checks"]["3_points"] = tot
    rep["passed"] = bool(ok)
    rep["seconds"] = time.time() - t0
    with open(out_json, "x") as f:
        json.dump(rep, f, indent=1, sort_keys=True)
    lines = ["# nonneg_maps validate-geometry %s: %s" % (rep["created_utc"], "PASSED" if ok else "FAILED"),
             "# geometry %s; base %s" % (rep["geometry"], rep["base"]),
             "1 mesh vs elements.bin: " + json.dumps(m, sort_keys=True),
             "2 BDF independent: " + json.dumps(b, sort_keys=True),
             "3 points total: " + json.dumps(tot, sort_keys=True),
             "#label iPoint UTC records cos_equal max_abs_dcos shadow_kernel shadow_mismatch night_mismatch "
             "et/sun/rh_bitwise Q_H2O_rel_driver_table Q_CO2_rel_driver_table Q_rel_python_driver"]
    for p in P:
        lines.append("%s %d %s %d %d %.3e %d %d %d %s/%s/%s %.2e %.2e %.2e" % (
            p["label"], p["iPoint"], p["utc"], p["records"], p["cos_equal_records"], p["max_abs_dcos"],
            p["shadow_records_kernel"], p["shadow_mismatch"], p["night_mismatch"], p["et_equal"], p["sun_equal"],
            p["rh_equal"], p["Q_H2O_rel_driver_vs_table"], p["Q_CO2_rel_driver_vs_table"],
            max(p["Q_H2O_rel_python_vs_driver"], p["Q_CO2_rel_python_vs_driver"])))
    with open(Path(args.geometry) / "geometry_validation.txt", "x") as f:
        f.write("\n".join(lines) + "\n")
    say("validate-geometry %s (%.0f s): %d points, %d records, cos equal %d (max |d| %.2e), shadow mismatches %d, night "
        "mismatches %d, Q driver/table max rel %.2e, python/driver %.2e" % (
            "PASSED" if ok else "FAILED", rep["seconds"], tot["points"], tot["records"], tot["cos_equal"],
            tot["max_abs_dcos"], tot["shadow_mismatch"], tot["night_mismatch"], tot["max_Q_rel_driver_vs_table"],
            tot["max_Q_rel_python_vs_driver"]))
    return 0 if ok else 3


# ------------------------------------------------------------------------------------------------ writers
VTK_TYPE = {np.dtype("<f8"): "Float64", np.dtype("<f4"): "Float32", np.dtype("<i4"): "Int32",
            np.dtype("<i8"): "Int64", np.dtype("uint8"): "UInt8", np.dtype("int8"): "Int8"}


def xml_comment(text):
    return "<!--\n" + text.replace("--", "- -") + "\n-->\n"


def write_vtu(path, nodes, faces, cell_fields, field_data, comment):
    """VTK XML UnstructuredGrid, triangles (VTK type 5), every array INLINE base64 (format="binary": one base64 stream
    of the UInt64 byte count followed by the little-endian data, no compression; VTK's own layout). Inline, not
    appended: the cluster's ParaView 5.11.1 (VTK 9.2 with the RHEL expat-2.5.0-5.el9_6 reparse-deferral backport)
    cannot read appended data (raw or base64), even files written by VTK itself; inline base64 is read by every VTK
    XML reader. cell_fields: [(name, array)], field_data: [(name, list of floats)]. Returns the array descriptors."""
    nn, nc = len(nodes), len(faces)
    conn = np.ascontiguousarray(faces, dtype="<i4").reshape(-1)
    offs = np.arange(3, 3 * nc + 1, 3, dtype="<i4")
    types = np.full(nc, 5, dtype=np.uint8)
    blocks = [("Points", "Points", np.ascontiguousarray(nodes, dtype="<f8"), 3), ("Cells", "connectivity", conn, 1),
              ("Cells", "offsets", offs, 1), ("Cells", "types", types, 1)]
    for name, a in cell_fields:
        a = np.ascontiguousarray(a)
        if a.dtype not in VTK_TYPE or len(a) != nc:
            raise ValueError("cell field %s: dtype %s length %d" % (name, a.dtype, len(a)))
        blocks.append(("CellData", name, a, 1))
    desc = [{"section": sec, "name": name, "type": VTK_TYPE[a.dtype], "ncomp": ncomp, "nbytes": int(a.nbytes)}
            for sec, name, a, ncomp in blocks]

    def da(f, d, a):
        nc_ = ' NumberOfComponents="%d"' % d["ncomp"] if d["ncomp"] != 1 else ""
        f.write(('        <DataArray type="%s" Name="%s"%s format="binary">\n          ' % (
            d["type"], d["name"], nc_)).encode("ascii"))
        f.write(base64.b64encode(struct.pack("<Q", d["nbytes"]) + memoryview(a).cast("B")))
        f.write(b"\n        </DataArray>\n")
    with open(path, "xb") as f:
        x = ['<?xml version="1.0"?>\n', xml_comment(comment),
             '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n',
             '  <UnstructuredGrid>\n']
        if field_data:
            x.append('    <FieldData>\n')
            for name, vals in field_data:
                x.append('      <DataArray type="Float64" Name="%s" NumberOfTuples="%d" format="ascii">%s</DataArray>\n'
                         % (name, len(vals), " ".join("%.17g" % v for v in vals)))
            x.append('    </FieldData>\n')
        x.append('    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n      <Points>\n' % (nn, nc))
        f.write("".join(x).encode("ascii"))
        da(f, desc[0], blocks[0][2])
        f.write(b"      </Points>\n      <Cells>\n")
        for k in (1, 2, 3):
            da(f, desc[k], blocks[k][2])
        f.write(b"      </Cells>\n")
        if len(blocks) > 4:
            f.write(b"      <CellData>\n")
            for k in range(4, len(blocks)):
                da(f, desc[k], blocks[k][2])
            f.write(b"      </CellData>\n")
        f.write(b"    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n")
    return desc


def read_vtu(path):
    """The inline arrays of a write_vtu file: ({name: array}, {field data name: values}, xml head) (verification)."""
    with open(path, "rb") as f:
        d = f.read()
    rev = {v: k for k, v in VTK_TYPE.items()}
    out = {}
    for m in re.finditer(rb'<DataArray type="(\w+)" Name="([^"]+)"(?: NumberOfComponents="(\d+)")? format="binary">'
                         rb'\s*([A-Za-z0-9+/=]*)\s*</DataArray>', d):
        typ, name, ncomp = m.group(1).decode(), m.group(2).decode(), int(m.group(3) or 1)
        raw = base64.b64decode(m.group(4))
        nbytes = struct.unpack("<Q", raw[:8])[0]
        if nbytes != len(raw) - 8:
            raise ValueError("%s: array %s header %d bytes, data %d" % (path, name, nbytes, len(raw) - 8))
        a = np.frombuffer(raw, dtype=rev[typ], offset=8)
        out[name] = a.reshape(-1, ncomp) if ncomp > 1 else a
    head = d[:d.index(b"<Points>")].decode("ascii")
    fd = {m.group(1): [float(v) for v in m.group(2).split()]
          for m in re.finditer(r'<DataArray type="Float64" Name="([^"]+)" NumberOfTuples="\d+" format="ascii">([^<]*)<',
                               head)}
    return out, fd, head


def _fmt_chunk(task):
    a, kind = task
    f = FLOAT_FMT if kind == "f" else "%d"
    if kind == "c":                       # connectivity: one triangle per line
        return ("\n".join("%d %d %d" % tuple(r) for r in a.tolist()) + "\n").encode("ascii")
    vals = a.tolist()
    m = len(vals) // PER_LINE * PER_LINE
    lf = " ".join([f] * PER_LINE)
    parts = [lf % tuple(vals[k:k + PER_LINE]) for k in range(0, m, PER_LINE)]
    if m < len(vals):
        parts.append(" ".join(f % v for v in vals[m:]))
    return ("\n".join(parts) + "\n").encode("ascii")


def _parse_chunk(raw):
    return np.fromstring(raw.decode("ascii"), dtype=np.float64, sep=" ")


class TecplotWriter:
    """Tecplot ASCII FE surface: FETRIANGLE zones, BLOCK packing, X Y Z nodal, other variables cell-centred; the first
    zone writes the mesh and every variable, later zones share X..static variables (VARSHARELIST) and the connectivity
    (CONNECTIVITYSHAREZONE) of zone 1 and write the time-dependent variables only (one zone per time, STRANDID 1,
    SOLUTIONTIME = seconds after 2016-09-30T00:00:00 UTC). Byte offsets of every block are recorded (verification)."""

    def __init__(self, path, pool, title, variables, dtypes, n_static, comment, aux):
        self.path, self.pool = path, pool
        self.variables, self.dtypes, self.n_static = variables, dtypes, n_static
        self.f = open(path, "xb")
        self.blocks = []
        self.nzones = 0
        head = ["# " + line for line in comment.splitlines()]
        head.append('TITLE = "%s"' % title)
        head.append("VARIABLES = " + " ".join('"%s"' % v for v in variables))
        for k, v in aux:
            head.append('DATASETAUXDATA %s = "%s"' % (k, v.replace('"', "'")))
        self.f.write(("\n".join(head) + "\n").encode("ascii"))

    def _block(self, zone, var, a, kind):
        start = self.f.tell()
        if kind == "c":
            tasks = ((a[k:k + CHUNK // 8], kind) for k in range(0, len(a), CHUNK // 8))
        else:
            tasks = ((a[k:k + CHUNK], kind) for k in range(0, len(a), CHUNK))
        for b in self.pool.imap(_fmt_chunk, tasks):
            self.f.write(b)
        self.blocks.append({"zone": zone, "var": var, "start": start, "end": self.f.tell(), "kind": kind,
                            "n": int(a.size)})

    def zone(self, title, solution_time, nodes, nelem, values, conn=None, zone_aux=()):
        """values: {variable name: array} (zone 1: all variables; later zones: the time-dependent ones)."""
        self.nzones += 1
        nv = len(self.variables)
        first = self.nzones == 1
        h = ['ZONE T="%s", ZONETYPE=FETRIANGLE, NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK' % (title, len(nodes), nelem),
             "VARLOCATION=([4-%d]=CELLCENTERED), STRANDID=1, SOLUTIONTIME=%.17g" % (nv, solution_time),
             "DT=(%s)" % " ".join(self.dtypes)]
        if not first:
            h.append("VARSHARELIST=([1-%d]=1), CONNECTIVITYSHAREZONE=1" % (3 + self.n_static))
        for k, v in zone_aux:
            h.append('AUXDATA %s = "%s"' % (k, v))
        self.f.write(("\n".join(h) + "\n").encode("ascii"))
        z = self.nzones
        if first:
            for c in range(3):
                self._block(z, self.variables[c], np.ascontiguousarray(nodes[:, c]), "f")
        for v in self.variables[3 if first else 3 + self.n_static:]:
            a = values[v]
            self._block(z, v, a, "i" if a.dtype.kind in "ui" else "f")
        if first:
            self._block(z, "connectivity", np.ascontiguousarray(conn, dtype=np.int64) + 1, "c")

    def close(self):
        self.f.close()


def verify_tecplot(path, blocks, ref, pool):
    """Re-read every block of the Tecplot file and compare with the written arrays (floats: |x - x_ref| <=
    5e-10 |x_ref|; integers and connectivity exact)."""
    res = {"blocks": len(blocks), "values": 0, "max_rel_float": 0.0, "int_mismatch": 0, "passed": True}
    mm = np.memmap(path, dtype=np.uint8, mode="r")
    for b in blocks:
        a = ref[(b["zone"], b["var"])]
        a = a.reshape(-1)
        raw = bytes(mm[b["start"]:b["end"]])
        pieces, k0, L = [], 0, len(raw)
        while k0 < L:                       # pieces of ~64 MB cut at line ends
            c = min(k0 + (1 << 26), L)
            k1 = L if c == L else raw.rfind(b"\n", k0, c) + 1
            if k1 <= k0:
                k1 = raw.find(b"\n", c) + 1 or L
            pieces.append(raw[k0:k1])
            k0 = k1
        x = np.concatenate(list(pool.imap(_parse_chunk, pieces)))
        res["values"] += len(x)
        if len(x) != len(a):
            res["passed"] = False
            res.setdefault("count_errors", []).append([b["zone"], b["var"], len(x), len(a)])
            continue
        if b["kind"] == "f":
            ref_ = a.astype(np.float64)
            d = np.abs(x - ref_)
            bad = d > TEC_REL_TOL * np.abs(ref_)
            rel = float(np.max(np.where(ref_ != 0, d / np.where(ref_ != 0, np.abs(ref_), 1.0), d)))
            res["max_rel_float"] = max(res["max_rel_float"], rel)
            if np.any(bad):
                res["passed"] = False
                res.setdefault("float_errors", []).append([b["zone"], b["var"], int(np.sum(bad))])
        else:
            n = int(np.sum(x != a.astype(np.float64)))
            res["int_mismatch"] += n
            if n:
                res["passed"] = False
    return res


def write_pvd(path, entries):
    x = ['<?xml version="1.0"?>\n', '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n',
         '  <Collection>\n']
    for t, fname in entries:
        x.append('    <DataSet timestep="%.17g" group="" part="0" file="%s"/>\n' % (t, fname))
    x.append('  </Collection>\n</VTKFile>\n')
    with open(path, "x") as f:
        f.write("".join(x))


# ------------------------------------------------------------------------------------------------ maps
def load_module(path, name, expect_sha256=None):
    path = Path(path)
    if expect_sha256 is not None and sha256_file(path) != expect_sha256:
        raise SystemExit("%s: sha256 differs from the record %s" % (path, expect_sha256))
    if str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))         # its own helpers (slope_mult) come from the same snapshot
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_result_kappa(path):
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) >= 3 and t[0] == "kappa" and t[1].startswith(("+", "-")):
                return float(t[2])
    return None


def state_activity(ft_dir, fit_dir, el, pool_threads):
    """The final state's per-facet fields + checks: A (A_final), kappa, seen, bins, and the reproduction checks."""
    ft_dir, fit_dir = Path(ft_dir), Path(fit_dir)
    fres = (ft_dir / "finetune_result.txt").read_text()
    if not re.search(r"^complete yes$", fres, re.M):
        raise SystemExit("%s: the fine tuning is not complete" % ft_dir)
    ftj = json.load(open(ft_dir / "fine_tuning.json"))
    if Path(ftj["fit"]["directory"]).resolve() != fit_dir.resolve():
        raise SystemExit("the fine tuning %s was made from %s, not from %s" % (ft_dir, ftj["fit"]["directory"], fit_dir))
    man = json.load(open(fit_dir / "fit_manifest.json"))
    nn = man["model_extension"]["nonneg"]
    model, s_ratio = nn["config"]["model"], nn["config"]["s_ratio"]
    scale = float(man["parameters"]["enhancement_scale"])
    if scale != float(nn["model_check"]["A_scale"]):
        raise SystemExit("enhancement scale differs between parameters and model_check")
    nmrec = ftj["modules"]["nonneg_mult.py"]
    nm = load_module(nmrec["path"], "nonneg_mult_of_state", nmrec["sha256"])
    smrec = ftj["modules"].get("slope_mult.py")
    sm_file = str(Path(sys.modules["slope_mult"].__file__).resolve())
    if smrec is not None and (sm_file != str(Path(smrec["path"]).resolve()) or sha256_file(sm_file) != smrec["sha256"]):
        raise SystemExit("slope_mult imported from %s, the fine tuning used %s" % (sm_file, smrec["path"]))
    with np.load(fit_dir / "nonneg_facets.npz", allow_pickle=False) as z:
        u, TS, A_fit, bins, seen_fit = (np.asarray(z[k]) for k in ("u", "TS", "A", "bin", "seen"))
        held_fit = np.asarray(z["held"]) if "held" in z.files else np.zeros(len(u))
    # finetune_state.npz "seen" is the fine tuning's PARAMETER mask (positive_finetune.Objective(..., params, ...)),
    # not visibility; the seen flag of the maps is the fit's (nonneg_facets.npz: facet contributes to a fit point,
    # fit_cops.Model.seen_faces)
    with np.load(ft_dir / "finetune_state.npz", allow_pickle=False) as z:
        A_final, held_prev, held_next, delta, ft_params, late = (np.asarray(z[k]) for k in (
            "A_final", "held_prev", "held_next", "delta", "seen", "late"))
    with np.load(fit_dir / "fit_restart.npz", allow_pickle=False) as z:
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        tv = np.asarray(z["template_values"], np.float64)
    nb = len(nm.BIN_NAMES)
    m = tv[:nb]
    n = el["nfaces"]
    for name, a in (("u", u), ("A_final", A_final), ("bins", bins)):
        if len(a) != n:
            raise SystemExit("%s has %d facets, the base %d" % (name, len(a), n))
    s = s_ratio * el["AH2O"] if model == "softplus" else None
    mrow = m[bins.astype(np.intp)]
    chk = {"model": model, "s_ratio": s_ratio, "enhancement_scale": scale, "kappa": kappa,
           "multipliers": [float(v) for v in m], "bin_names": list(nm.BIN_NAMES),
           "bin_labels_deg": list(nm.BIN_LABELS), "nonneg_mult": nmrec,
           "slope_mult": {"path": sm_file, "sha256": sha256_file(sm_file), "record": smrec}}
    chk["held_fit_equals_held_prev"] = bool(np.array_equal(held_fit, held_prev))
    chk["held_next_equals_held_prev_plus_delta"] = bool(np.array_equal(held_next, held_prev + delta))
    gh_fit = np.exp(held_fit) if np.any(held_fit != 0.0) else None
    A_fit_re = nm.fields(model, u, TS, s, mrow, gh_fit, scale, False)[0]
    chk["A_fit_reproduced_bitwise"] = bool(np.array_equal(A_fit_re, A_fit))
    A_next = nm.fields(model, u, TS, s, mrow, np.exp(held_next), scale, False)[0]
    tuned = delta != 0.0
    chk["tuned_facets"] = int(np.sum(tuned))
    chk["A_final_equals_A_fit_at_untuned"] = bool(np.array_equal(A_final[~tuned], A_fit[~tuned]))
    chk["A_final_equals_recomputed_at_untuned"] = bool(np.array_equal(A_final[~tuned], A_next[~tuned]))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(A_final[tuned] / A_next[tuned] - 1.0)
    chk["A_final_vs_recomputed_max_rel_at_tuned"] = float(np.max(rel)) if rel.size else 0.0
    chk["A_final_min"] = float(np.min(A_final))
    chk["A_final_zero_facets"] = int(np.sum(A_final == 0.0))
    mpar = re.search(r"parameters (\d+)", fres)
    chk["ft_parameter_facets"] = int(np.sum(ft_params))
    chk["ft_parameter_facets_result_txt"] = int(mpar.group(1)) if mpar else None
    chk["ft_parameters_subset_of_fit_seen"] = bool(np.all(seen_fit[ft_params.astype(bool)]))
    chk["delta_only_at_parameters"] = bool(np.all(ft_params.astype(bool)[delta != 0.0]))
    chk["seen_fit_facets"] = int(np.sum(seen_fit))
    chk["seen_ok"] = bool(chk["ft_parameter_facets"] == chk["ft_parameter_facets_result_txt"] and
                          chk["ft_parameters_subset_of_fit_seen"] and chk["delta_only_at_parameters"])
    try:
        from numpy.core._multiarray_umath import __cpu_features__ as cpuf
        chk["numpy_cpu_features"] = {k: bool(cpuf.get(k)) for k in ("AVX2", "FMA3", "AVX512F", "AVX512_SKX")}
    except ImportError:
        chk["numpy_cpu_features"] = None
    chk["NPY_DISABLE_CPU_FEATURES"] = os.environ.get("NPY_DISABLE_CPU_FEATURES")
    kres = parse_result_kappa(fit_dir / "fit_result.txt")
    chk["kappa_fit_result_txt"] = kres
    chk["kappa_matches_fit_result"] = kres is not None and abs(kres - kappa) <= 1e-6 * abs(kappa)
    chk["passed"] = bool(chk["held_fit_equals_held_prev"] and chk["held_next_equals_held_prev_plus_delta"] and
                         chk["A_fit_reproduced_bitwise"] and chk["A_final_equals_A_fit_at_untuned"] and
                         chk["A_final_equals_recomputed_at_untuned"] and
                         chk["A_final_vs_recomputed_max_rel_at_tuned"] < 1e-12 and chk["A_final_min"] >= 0.0 and
                         chk["seen_ok"] and chk["kappa_matches_fit_result"])
    return {"A": A_final, "kappa": kappa, "seen": seen_fit.astype(bool), "bins": bins.astype(np.uint8), "late": late,
            "model": model, "checks": chk, "manifest": man, "ft_json": ftj}


def stats(a):
    a = np.asarray(a)
    return {"min": float(np.min(a)), "max": float(np.max(a)), "sum": float(np.sum(a, dtype=np.float64)),
            "dtype": str(a.dtype), "n": int(a.size)}


def human(nbytes):
    if nbytes < 1 << 20:
        return "%.1f KiB" % (nbytes / 1024.0)
    if nbytes < 1 << 30:
        return "%.1f MiB" % (nbytes / 2.0 ** 20)
    return "%.2f GiB" % (nbytes / 2.0 ** 30)


def cmd_maps(args):
    t0 = time.time()
    out = Path(args.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten; move it aside first)" % out)
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    gval = json.load(open(Path(args.geometry) / "geometry_validation.json"))
    if not gval.get("passed"):
        raise SystemExit("geometry %s did not pass validate-geometry" % args.geometry)
    base = Path(args.base)
    if Path(gval["base"]) != base.resolve():
        raise SystemExit("geometry validated against %s, not %s" % (gval["base"], base))
    fc, header, cols, rows = base_table(base)
    law = law_guard(fc)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    G = read_geometry(args.geometry)
    el = read_elements(base / "elements.bin")
    nf = el["nfaces"]
    if not (np.array_equal(G["facegeom"][:, 6], el["area"]) and np.array_equal(G["facegeom"][:, 0:3], el["centroid"])):
        raise SystemExit("geometry mesh is not the base's (area / centroid)")
    say("state: ft %s, fit %s" % (Path(args.ft_dir).resolve(), Path(args.fit_dir).resolve()))
    mrun = re.search(r"^# run: (\S+)$", (Path(args.fit_dir) / "fit_result.txt").read_text(), re.M)
    if mrun is None or Path(mrun.group(1)).resolve() != base.resolve():
        raise SystemExit("the fit %s ran on %s, not on the base %s" % (args.fit_dir, mrun and mrun.group(1), base))
    S = state_activity(args.ft_dir, args.fit_dir, el, args.procs)
    say("state checks: " + json.dumps({k: v for k, v in S["checks"].items() if k not in ("nonneg_mult",)},
                                      sort_keys=True))
    if not S["checks"]["passed"]:
        raise SystemExit("state checks failed")
    model = S["model"]
    if args.expect_model and args.expect_model != model:
        raise SystemExit("state is model %s, expected %s" % (model, args.expect_model))
    A, kappa = S["A"], S["kappa"]
    co2_act = kappa * el["ACO2"]
    area = el["area"]
    stage.mkdir(parents=True)
    nodes, faces = G["nodes"], G["faces"]
    static = [("A_H2O", A), ("CO2_activity_kappa_ACO2", co2_act), ("seen", S["seen"].astype(np.uint8)),
              ("slope_bin", S["bins"]), ("area_m2", area)]
    bins_txt = ", ".join("%d = %s deg" % (k, lab) for k, lab in enumerate(S["checks"]["bin_labels_deg"]))
    ill_txt = ", ".join("%d = %s (%s)" % t for t in ILLUMINATION)
    unit_txt = ("UNITS AND LEGEND\n"
                "F_H2O_per_m2_s, F_CO2_per_m2_s: instantaneous emitted molecular fluxes per facet [m^-2 s^-1]\n"
                "  F_H2O = q_H2O(cos SZA, night/shadow, r_h) * A_H2O, F_CO2 = kappa * q_CO2 * ACO2 (fit_cops.Model "
                "load_kernels / prepare)\n"
                "  q_s = qmin_s (night or cast shadow) else qmin_s + (qmax_s - qmin_s) min(cos SZA, 1);\n"
                "  H2O qmin = 0.2*600/rhA^p, qmax = 600/rhA^p, p = %s; CO2 qmin = 0.1*600/rhA^2, qmax = 600/rhA^2; "
                "rhA = r_h/AU\n"
                "A_H2O: the model's H2O activity factor (%s model, final state after the last fine tuning); "
                "CO2_activity_kappa_ACO2: kappa * ACO2 (kappa = %.17g).\n"
                "  The code defines F = q A in m^-2 s^-1 and assigns no separate units to q and A.\n"
                "cos_SZA: cosine of the solar zenith angle (facet normal . Sun direction)\n"
                "illumination: %s\n"
                "seen: 1 = facet contributes to at least one fit point (the fit's nonneg_facets.npz seen = "
                "fit_cops.Model.seen_faces), 0 otherwise\n"
                "slope_bin: effective-gravitational-slope bin of the multipliers: %s\n"
                "area_m2: facet area [m^2]; coordinates X Y Z [m] in the body frame 67P/C-G_CK\n"
                "TimeValue / SOLUTIONTIME: seconds after 2016-09-30T00:00:00 UTC") % (
        "%.8g" % q_factors(1.0, helio, 1.0)["p"], model, kappa, ill_txt, bins_txt)
    pool = mp.get_context("fork").Pool(args.procs)
    files, man_times = [], []
    try:
        # static VTU
        cmt = ("nonneg_maps (%s) %s model, time-independent fields. State: ft %s, fit %s.\n%s" % (
            SCHEMA, model, Path(args.ft_dir).resolve(), Path(args.fit_dir).resolve(), unit_txt))
        f_static = "%s_static.vtu" % model
        write_vtu(stage / f_static, nodes, faces, static, [("kappa", [kappa])], cmt)
        files.append(f_static)
        say("wrote %s" % f_static)
        # Tecplot
        tvars = ["X [m]", "Y [m]", "Z [m]", "A_H2O", "CO2_activity_kappa_ACO2", "seen", "slope_bin", "area [m^2]",
                 "cos_SZA", "illumination", "F_H2O [m^-2 s^-1]", "F_CO2 [m^-2 s^-1]"]
        tdt = ["DOUBLE"] * 5 + ["BYTE", "BYTE", "DOUBLE", "DOUBLE", "BYTE", "DOUBLE", "DOUBLE"]
        f_tec = "%s_maps.dat" % model
        aux = [("model", model), ("kappa", "%.17g" % kappa), ("ft_dir", str(Path(args.ft_dir).resolve())),
               ("fit_dir", str(Path(args.fit_dir).resolve())), ("schema", SCHEMA)]
        tec = TecplotWriter(stage / f_tec, pool, "%s model: activity and instantaneous source rates (13 times)" % model,
                            tvars, tdt, 5, "nonneg_maps (%s) %s model; floats %s (the .vtu files hold the exact "
                            "float64 values).\n%s\nZones: one per time; zone 1 holds the mesh and the time-independent "
                            "variables, later zones share them (VARSHARELIST, CONNECTIVITYSHAREZONE)." % (
                                SCHEMA, model, FLOAT_FMT, unit_txt), aux)
        tec_ref = {}
        pvd = []
        for lab, utc in zip(MAP_LABELS, MAP_TIMES):
            T = G["times"][lab]
            if T["utc"] != utc:
                raise SystemExit("geometry %s is %s, expected %s" % (lab, T["utc"], utc))
            cos, sh = geometry_time(args.geometry, lab, nf)
            F_h, F_c, night, qf = source_rates(cos, sh, T["rh_m"], helio, sza, A, el["ACO2"], kappa, fc.AU)
            ill = illumination(cos, sh)
            if not (np.all(np.isfinite(F_h)) and np.all(np.isfinite(F_c)) and np.all(F_h >= 0) and np.all(F_c >= 0)):
                raise SystemExit("non-finite or negative flux at %s" % utc)
            tsec = seconds_after_day0(utc)
            stamp = utc.replace("-", "").replace(":", "")
            fname = "%s_flux_%s.vtu" % (model, stamp)
            tfields = [("cos_SZA", cos), ("illumination", ill), ("F_H2O_per_m2_s", F_h), ("F_CO2_per_m2_s", F_c)]
            fd = [("TimeValue", [tsec]), ("et_s", [T["et"]]), ("rh_m", [T["rh_m"]]), ("sun_direction", list(T["sun"])),
                  ("q_H2O_min_max", list(qf["H2O"])), ("q_CO2_min_max", list(qf["CO2"])), ("kappa", [kappa])]
            write_vtu(stage / fname, nodes, faces, tfields, fd, "nonneg_maps (%s) %s model at %s UTC (TimeValue %g s "
                      "after 2016-09-30T00:00:00 UTC; r_h %.9f AU).\n%s" % (SCHEMA, model, utc, tsec, qf["rhA"],
                                                                           unit_txt))
            files.append(fname)
            pvd.append((tsec, fname))
            zone_values = {"cos_SZA": cos, "illumination": ill, "F_H2O [m^-2 s^-1]": F_h, "F_CO2 [m^-2 s^-1]": F_c}
            zaux = [("UTC", utc), ("et_s", "%.17g" % T["et"]), ("rh_m", "%.17g" % T["rh_m"]),
                    ("sun_direction", " ".join("%.17g" % v for v in T["sun"]))]
            if lab == MAP_LABELS[0]:
                zone_values.update({"A_H2O": A, "CO2_activity_kappa_ACO2": co2_act, "seen": S["seen"].astype(np.uint8),
                                    "slope_bin": S["bins"], "area [m^2]": area})
                tec.zone(utc + " UTC", tsec, nodes, nf, zone_values, conn=faces, zone_aux=zaux)
                for c, v in enumerate(tvars[:3]):
                    tec_ref[(1, v)] = nodes[:, c]
                tec_ref[(1, "connectivity")] = faces.astype(np.int64) + 1
                for v in tvars[3:8]:
                    tec_ref[(1, v)] = zone_values[v]
            else:
                tec.zone(utc + " UTC", tsec, nodes, nf, zone_values, zone_aux=zaux)
            for v in tvars[8:]:
                tec_ref[(tec.nzones, v)] = zone_values[v]
            Qh, Qc = float(np.sum(F_h * area)), float(np.sum(F_c * area))
            cnt = {name: int(np.sum(ill == k)) for k, name, _ in ILLUMINATION}
            ar = {name: float(np.sum(area[ill == k])) for k, name, _ in ILLUMINATION}
            man_times.append({"label": lab, "utc": utc, "TimeValue_s": tsec, "et_s": T["et"], "rh_m": T["rh_m"],
                              "rh_AU": qf["rhA"], "sun_direction": list(T["sun"]),
                              "q_H2O_min_max": [float(v) for v in qf["H2O"]],
                              "q_CO2_min_max": [float(v) for v in qf["CO2"]], "p": qf["p"],
                              "Q_H2O_s-1": Qh, "Q_CO2_s-1": Qc, "facets": cnt, "area_m2": ar,
                              "vtu": fname, "stats": {"cos_SZA": stats(cos), "F_H2O": stats(F_h), "F_CO2": stats(F_c)}})
            say("%s %s: Q_H2O %.6e s^-1, Q_CO2 %.6e s^-1; facets day %d night %d shadow %d; wrote %s + Tecplot zone %d" % (
                lab, utc, Qh, Qc, cnt["day"], cnt["night"], cnt["cast_shadow"], fname, tec.nzones))
        tec.close()
        files.append(f_tec)
        f_pvd = "%s_flux.pvd" % model
        write_pvd(stage / f_pvd, pvd)
        files.append(f_pvd)
        # verification: read every file back
        ver = {}
        if not args.no_verify:
            say("verify: VTU files (bitwise) and the Tecplot file (re-parse)")
            v_ok = True
            rd, fdd, _ = read_vtu(stage / f_static)
            for name, a in static:
                v_ok &= bool(np.array_equal(rd[name], a))
            v_ok &= bool(np.array_equal(rd["Points"], nodes)) and bool(np.array_equal(rd["connectivity"],
                                                                                      faces.reshape(-1)))
            for lab, utc in zip(MAP_LABELS, MAP_TIMES):
                T = G["times"][lab]
                cos, sh = geometry_time(args.geometry, lab, nf)
                F_h, F_c, _, _ = source_rates(cos, sh, T["rh_m"], helio, sza, A, el["ACO2"], kappa, fc.AU)
                rd, fdd, _ = read_vtu(stage / [t["vtu"] for t in man_times if t["label"] == lab][0])
                v_ok &= bool(np.array_equal(rd["cos_SZA"], cos) and np.array_equal(rd["F_H2O_per_m2_s"], F_h) and
                             np.array_equal(rd["F_CO2_per_m2_s"], F_c) and
                             np.array_equal(rd["illumination"], illumination(cos, sh)) and
                             np.array_equal(rd["Points"], nodes))
                v_ok &= fdd["TimeValue"][0] == seconds_after_day0(utc)
            ver["vtu_bitwise"] = bool(v_ok)
            ver["tecplot"] = verify_tecplot(stage / f_tec, tec.blocks, tec_ref, pool)
            ver["passed"] = bool(v_ok and ver["tecplot"]["passed"])
            say("verify: " + json.dumps(ver, sort_keys=True))
            if not ver["passed"]:
                raise SystemExit("verification of the written files failed")
    finally:
        pool.close()
        pool.join()
    sizes = {f: os.path.getsize(stage / f) for f in files}
    code = {p.name: sha256_file(p) for p in sorted(HERE.glob("*nonneg_maps*")) if p.is_file() and p.suffix in
            (".py", ".cpp")}
    man = {"schema": SCHEMA, "created_utc": utc_now(), "host": socket.gethostname(),
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv, "model": model,
           "state": {"ft_dir": str(Path(args.ft_dir)), "ft_dir_resolved": str(Path(args.ft_dir).resolve()),
                     "fit_dir": str(Path(args.fit_dir)), "fit_dir_resolved": str(Path(args.fit_dir).resolve()),
                     "inputs_sha256": {n: sha256_file(Path(d) / n) for d, n in (
                         (args.ft_dir, "finetune_state.npz"), (args.ft_dir, "finetune_result.txt"),
                         (args.fit_dir, "nonneg_facets.npz"), (args.fit_dir, "fit_restart.npz"),
                         (args.fit_dir, "fit_manifest.json"), (args.fit_dir, "fit_result.txt"))}},
           "base": str(base.resolve()), "geometry": str(Path(args.geometry).resolve()),
           "geometry_validation": {"passed": gval["passed"], "checks": gval["checks"]},
           "law": law, "heliocentric_model": helio, "sza_mode": sza, "state_checks": S["checks"],
           "static_stats": {name: stats(a) for name, a in static}, "times": man_times, "verification": ver,
           "files": sizes, "code_sha256": code, "code_dir": str(HERE),
           "illumination_codes": ILLUMINATION, "tecplot_float_format": FLOAT_FMT, "wall_seconds": time.time() - t0}
    with open(stage / "maps_manifest.json", "x") as f:
        json.dump(man, f, indent=1, sort_keys=True)
    write_readme(stage, man, unit_txt, sizes, bins_txt)
    stage.rename(out)
    say("maps of %s in %s (%.0f s; %s in %d files)" % (model, out, time.time() - t0, human(sum(sizes.values())),
                                                      len(sizes)))
    return 0


def write_readme(stage, man, unit_txt, sizes, bins_txt):
    c = man["state_checks"]
    L = ["Activity and instantaneous source-rate maps -- %s model, FINAL state (after the last fine tuning)" % man["model"],
         "Written %s by nonneg_maps.py (%s) on %s (Slurm job %s)." % (man["created_utc"], man["schema"], man["host"],
                                                                   man["slurm_job_id"]),
         "",
         "STATE",
         "  fine tuning: %s (-> %s)" % (man["state"]["ft_dir"], man["state"]["ft_dir_resolved"]),
         "  fit:         %s (-> %s)" % (man["state"]["fit_dir"], man["state"]["fit_dir_resolved"]),
         "  A_H2O = finetune_state.npz A_final; kappa = fit_restart.npz kappa = %.17g; multipliers %s (bins %s)" % (
             c["kappa"], " ".join("%.6g" % v for v in c["multipliers"]), " ".join(c["bin_labels_deg"])),
         "  checks: A_fit reproduced bitwise from nonneg_facets.npz %s; A_final = A_fit at the %d untuned facets %s, "
         "= m_b q exp(held_next) there %s; at the %d tuned facets max rel |A_final / (m_b q exp(held_next)) - 1| = %.2e;"
         " seen = the fit's nonneg_facets.npz seen (%d facets; the fine tuning's %d parameter facets are a subset: %s);"
         " passed %s" % (
             c["A_fit_reproduced_bitwise"], man["static_stats"]["A_H2O"]["n"] - c["tuned_facets"],
             c["A_final_equals_A_fit_at_untuned"], c["A_final_equals_recomputed_at_untuned"], c["tuned_facets"],
             c["A_final_vs_recomputed_max_rel_at_tuned"], c["seen_fit_facets"], c["ft_parameter_facets"],
             c["ft_parameters_subset_of_fit_seen"], c["passed"]),
         "  base: %s; geometry: %s (validate-geometry passed %s)" % (man["base"], man["geometry"],
                                                                     man["geometry_validation"]["passed"]),
         "",
         "THE LAW (exactly the fitter's: ../LiouvilleCOPS/tools/fit_cops.py Model.load_kernels and Model.prepare;",
         "  = src/Source.cpp cSourceModel::GetQ/EvaluateElement; = cops_cuda/kernels/pressure.cu)",
         "  night = (cos SZA <= 0) or cast shadow; SZAMode %s; HeliocentricModel %s" % (man["sza_mode"],
                                                                                         man["heliocentric_model"]),
         "  geometry: Sun = spkpos_c(SUN, et, 67P/C-G_CK, LT+S, CHURYUMOV-GERASIMENKO) (Geometry.cpp); cos SZA = normal .",
         "  Sun direction; cast shadow = ray from centroid + 1e-6 MaxEdge normal toward the Sun hits the nucleus (Source.cpp",
         "  Evaluate / gpu_liouville.cu k_shadow), only for cos SZA > 0; for ALL facets (not only those seen by the gauges).",
         "",
         unit_txt,
         "",
         "FILES (sizes)"]
    for f, s in sorted(sizes.items()):
        L.append("  %-44s %12d bytes (%s)" % (f, s, human(s)))
    L += ["  maps_manifest.json  every input (sha256), check, statistic and the code sha256",
          "",
          "  ParaView: open %s_flux.pvd (13 time steps: cos_SZA, illumination, F_H2O_per_m2_s, F_CO2_per_m2_s)" % man["model"],
          "            and %s_static.vtu (A_H2O, CO2_activity_kappa_ACO2, seen, slope_bin, area_m2); same mesh, cell" % man["model"],
          "            data (use 'Append Attributes' to combine). VTK XML with inline base64 binary arrays (format=binary,",
          "            UInt64 headers, uncompressed), exact float64; inline because the cluster's ParaView 5.11.1 cannot read",
          "            appended VTK XML data (RHEL expat-2.5.0-5.el9_6), inline base64 is read by every VTK XML reader.",
          "  Tecplot:  %s_maps.dat, ASCII, 13 FETRIANGLE zones (one per time, SOLUTIONTIME = s after 2016-09-30T00:00:00" % man["model"],
          "            UTC, STRANDID 1); zone 1 holds the mesh and the time-independent variables, zones 2-13 share them",
          "            (VARSHARELIST, CONNECTIVITYSHAREZONE). Floats %s (relative rounding <= 5e-10)." % man["tecplot_float_format"],
          "            The .dat file is for Tecplot 360. ParaView's Tecplot reader (vtkTecplotReader, ParaView 5.11) does not",
          "            implement VARSHARELIST / CONNECTIVITYSHAREZONE (zones 2-13 come out wrong there) and reads in float32:",
          "            in ParaView use the .vtu files.",
          "",
          "TIMES  (TimeValue = s after 2016-09-30T00:00:00 UTC; Q_s = sum_e F_s,e area_e; facets by illumination)",
          "#label UTC TimeValue_s rh_AU Q_H2O_s^-1 Q_CO2_s^-1 day night cast_shadow area_day_m2 area_night_m2 area_shadow_m2"]
    for t in man["times"]:
        L.append("%s %s %.0f %.9f %.6e %.6e %d %d %d %.6e %.6e %.6e" % (
            t["label"], t["utc"], t["TimeValue_s"], t["rh_AU"], t["Q_H2O_s-1"], t["Q_CO2_s-1"], t["facets"]["day"],
            t["facets"]["night"], t["facets"]["cast_shadow"], t["area_m2"]["day"], t["area_m2"]["night"],
            t["area_m2"]["cast_shadow"]))
    v = man.get("verification") or {}
    L += ["", "VERIFICATION of the written files: VTU bitwise %s; Tecplot re-parse %s (max rel %.2e)" % (
        v.get("vtu_bitwise"), (v.get("tecplot") or {}).get("passed"), (v.get("tecplot") or {}).get("max_rel_float", 0)),
          "", "COMMAND", "  " + " ".join(man["command"])]
    with open(stage / "README.txt", "x") as f:
        f.write("\n".join(L) + "\n")


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("times")
    a.add_argument("--base", required=True)
    a.add_argument("--out", required=True)
    a = sub.add_parser("validate-geometry")
    a.add_argument("--geometry", required=True)
    a.add_argument("--base", required=True)
    a.add_argument("--skip-bdf-hash", action="store_true")
    a = sub.add_parser("maps")
    a.add_argument("--ft-dir", required=True, help="c<N>/ft/final of the chain's last cycle")
    a.add_argument("--fit-dir", required=True, help="c<N>/fit/final (the fit the fine tuning started from)")
    a.add_argument("--base", required=True)
    a.add_argument("--geometry", required=True)
    a.add_argument("--out", required=True, help="new directory (e.g. <chain>/maps)")
    a.add_argument("--expect-model", choices=("softplus", "squared"))
    a.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    a.add_argument("--no-verify", action="store_true")
    args = ap.parse_args(argv)
    return {"times": cmd_times, "validate-geometry": cmd_validate_geometry, "maps": cmd_maps}[args.cmd](args) or 0


if __name__ == "__main__":
    sys.exit(main())
