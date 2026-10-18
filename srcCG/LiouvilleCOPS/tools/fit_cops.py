#!/usr/bin/env python3
"""Fit of the surface source rate and gas temperature to the COPS NG and RG observations.

The fit uses the geometric kernels of a LiouvilleCOPS run (FlightTimeRotation=none, WriteKernels=on;
files <run>/kernels/point=NNNN.bin, <run>/elements.bin, <run>/COPS_model.dat). For every element e
visible from observation point i:
   n_s  = sum_e F_s,e beta_s,e sqrt(pi)/2 Kn_e
   Phi_s = sum_e F_s,e [KPhi_e + (V.b) beta_s,e sqrt(pi)/2 KV_e],   beta=sqrt(m/(2kT))
   pNG = sum_s n_s k T_NG / S_s,   pRG = sum_s Phi_s C_RG,s(T_RG,i) / S_s
Parameters (srcCG model structure):
   H2O : F = Q_H2O(SZA, night, r_h) * max(0, sum_j a_j Y_j)       a_j: 25 spherical-harmonic coefficients (fitted)
   CO2 : F = kappa * Q_CO2(SZA, night, r_h) * A_CO2,srcCG           kappa: global CO2 scale (fitted or fixed)
   T   : analytical model  T = max(Tmin, a/cos(SZA)+b) for SZA<=75 deg, else Tmin; night/shadow: Tmin,
         a=(Tmax-T75)/(1-1/cos 75deg), b=Tmax-a                     (Tmax, T75, Tmin fitted or fixed)
         --ordered-temperature fits Tmin, T75-Tmin>=0, Tmax-T75>=0 instead
Misfit: log residuals of the used points (useNG/useRG=1 in COPS_model.dat), equal weight per gauge,
plus a Tikhonov term that ties the H2O activity map to the prior (the coefficients of the run, srcCG):
   cost = 1/2 [ (1/N_NG) sum (ln pNG_model - ln pNG_obs)^2 + (1/N_RG) sum (ln pRG_model - ln pRG_obs)^2
                + lambda D^2 ]
   D^2 = int_S (A(a) - A(a0))^2 dS / int_S A(a0)^2 dS,   A(a) = sum_j a_j Y_j (linear map, nucleus surface)
D is the rms deviation of the activity map from the prior relative to the rms prior activity. Regions not seen
from the trajectory do not affect the data term and stay at the prior.

DATA SOURCES
  Observations     refNG, refRG of <run>/COPS_model.dat: the COPS NG and RG pressures of srcCG/RosinaMeasurements.h
                   (identical to the ROSINA COPS L2 BG housekeeping pressures, PSA/PDS3 RO-C-ROSINA-2-EXT3-V2.0,
                   N2-equivalent, mbar -> Pa). Only the points with useNG/useRG = 1 are used (selection by
                   ExcludeNG, ExcludeRG, FitInterval of the base run; points without data or inside the shape model
                   are excluded).
                   --data FILE replaces these values and flags by a table (columns iPoint UTC pNG_Pa pRG_Pa useNG
                   useRG T_RG_K; data/COPS_fit_data.txt of tools/collect_cops_sources.py is the table of the base
                   run). The iPoint/UTC must be points of the run (kernels), T_RG must equal the run's, and points
                   with use=1 need p > 0 and the gauge outside the shape model; run points missing in the table
                   are not used.
  ExcludeNG        The ExcludeNG intervals of the run's input file (fit_base.in: 2016-09-30T10:24:40 .. 10:40:00,
                   from the NG filament emission range change HIGH->LOW until the impact) are NOT applied by default:
                   the NG data in them are used (with data > 0, gauge outside the shape model, inside FitInterval;
                   the base run has kernels for all points). --apply-exclude-ng blocks them (with --data: sets useNG=0
                   in these intervals; without --data the flags of the run, which have the intervals applied, are kept).
                   Without --data the default adds these points to the selection of the run.
  Saturated NG     The NG reading is saturated (held constant at 2.2633e-8 mbar, the maximum of the HIGH filament
                   emission range, filament emission current 0) in the L2 BG records 2016-09-30 10:23:46.168 ..
                   10:24:38.251, i.e. the 13 points 10:23:50 .. 10:24:38 of the srcCG table (time = end of the
                   interval). They are found in the housekeeping table (--hk-table, default data/COPS_HK_20160930.txt
                   of tools/cops_hk_table.py). By default they are USED in the fit; --exclude-ng-saturated removes them
                   (with or without --data). The choice is written to fit_result.txt and to the COPS_model.dat header.
  T_RG             T_RG column of COPS_model.dat (the COPS sensor temperature of the L2 BG housekeeping,
                   data/COPS_HK_20160930.txt, when the base run uses TemperatureRGMode = hk_table).
  Geometry         <run>/kernels/point=NNNN.bin (version 1): per visible surface element the geometric kernels Kn,
                   KPhi, KV, cos(SZA) and the shadow flag of the SPICE geometry of the base run (spacecraft position
                   and velocity, gauge pointing, Sun direction, shape model).
  Surface elements <run>/elements.bin: centroid, normal, area, the 25 spherical-harmonic basis values of the srcCG
                   activity model and the H2O/CO2 activity factors of the base run.
  Prior            the activity coefficients of the base run (srcCG tables Activity_SHAP5/Activity_SHAP5_1 in
                   src/SrcCGTables.h, or the ActivityFile) and its analytical temperature parameters.

Before fitting, the model is evaluated with the parameters of the run and compared with the pressures
computed by LiouvilleCOPS (consistency check of this implementation).

--lambda with one value: one fit, output in <out>.
--lambda with several values: L-curve; one fit per value (in parallel, --jobs) in <out>/lambda=<value>/, and the
  table <out>/lcurve.txt (data misfit vs. D, effective number of parameters, maximum-curvature lambda).

Data used in the fit (always written to <out>, before the fit; --data-only writes them and stops):
   fit_data_used.txt  the points and pressures that enter the misfit, taken from the fit problem itself (after the
                      selection of the run or --data, and --exclude-ng-saturated): iPoint UTC altitude_m pNG_used_Pa
                      pRG_used_Pa refNG_Pa refRG_Pa (-1: not used)
   fit_data_used.png  these data vs. UTC (whole descent, last 40 min) and vs. altitude; used points filled, data
                      not used open grey
   fit_profiles_used.png  only the used NG (solid) and RG (dashed) profiles vs. UTC and vs. altitude; the lines
                      break at unused points and at time gaps > 90 s

Output of each fit (default <AMPS>/cops_output/<run-name>_fit):
   fit_result.txt   parameters, uncertainties, misfit, singular values, activity diagnostics
   ActivityFit.txt  3x25 coefficients (H2O fitted, CO2=kappa*srcCG, CO=srcCG) for ActivityModel=file
   cops_fit.in      LiouvilleCOPS input that reruns the fitted model (writes <fit dir>_verify)
   COPS_model.dat   fitted aggregate pNG/pRG predictions. Source/species/production fields that this cached
                    Python evaluation does not recompute are explicitly nan; use the C++ verification table.
   fit_manifest.json versioned inputs/hashes, cache and selection identities, settings, optimizer state, and
                    independent C++-agreement/scientific-review states

Relative output paths use <AMPS>/cops_output (override: COPS_OUTPUT_ROOT).
Product inputs prefer that root, then retained application-local files.
Absolute output paths override the defaults only when outside LiouvilleCOPS.

usage: fit_cops.py <run directory> [options]
"""
import argparse
import datetime as dt
import json
import multiprocessing
import os
from pathlib import Path
import platform
import re
import struct
import sys

import matplotlib

from cops_paths import output_path, product_input
from cops_io import (FIT_SCHEMA, ProductError, canonical_sha256, file_identity,
                     selection_identity, utc_now, validate_fit_base,
                     verify_identity, write_json_atomic)

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.optimize import least_squares

KB = 1.380649E-23
AMU = 1.66053906660E-27
AU = 149597870700.0
MASS = {"H2O": 18.01528 * AMU, "CO2": 44.0095 * AMU}
SPECIES = ("H2O", "CO2")
SQRTPI = np.sqrt(np.pi)
COS75 = np.cos(np.radians(75.0))
KT = 1.0 / (1.0 - 1.0 / COS75)
NBASIS = 25
DEGREE_OF_BASIS = [0] + [1] * 3 + [2] * 5 + [3] * 7 + [4] * 9
TNAMES = ("Tmax", "T75", "Tmin")
FLOOR = 1e-300

# The active manifest is set only after the output directory has been created.
# Keeping it here lets the top-level exception handler mark an interrupted or
# failed process without converting the existing user-facing validation exits
# into tracebacks.
ACTIVE_MANIFEST = None
ACTIVE_MANIFEST_PATH = None


def save_active_manifest():
    """Atomically persist the current fit state, if a fit has started."""
    if ACTIVE_MANIFEST is not None and ACTIVE_MANIFEST_PATH is not None:
        write_json_atomic(ACTIVE_MANIFEST_PATH, ACTIVE_MANIFEST)


# ------------------------------------------------------------------------------------------- input
def read_model_file(fname):
    header, names, rows = [], None, []
    with open(fname) as f:
        for line in f:
            if line.startswith("#iPoint"):
                names = line[1:].split()
            elif line.startswith("#"):
                header.append(line[1:].strip())
            elif line.strip():
                rows.append(line.split())
    if names is None:
        sys.exit("Error: %s has no column header" % fname)
    return header, names, rows


def apply_data_file(fname, header, cols, rows):
    """replace refNG, refRG, useNG, useRG of the rows of the run by the table <fname> (tools/collect_cops_sources.py,
    data/COPS_fit_data.txt); the iPoint/UTC of the table must be those of the run, T_RG must agree (RG kernels)"""
    names, tab = None, {}
    for line in open(fname):
        if line.startswith("#iPoint"):
            names = line[1:].split()
        elif line.strip() and not line.startswith("#"):
            v = line.split()
            point = int(v[0])
            if point in tab:
                sys.exit("Error: %s has duplicate iPoint %d" % (fname, point))
            tab[point] = v
    if names is None:
        sys.exit("Error: %s has no column header (#iPoint ...)" % fname)
    t = {n: i for i, n in enumerate(names)}
    for k in ("iPoint", "UTC", "pNG_Pa", "pRG_Pa", "useNG", "useRG", "T_RG_K"):
        if k not in t:
            sys.exit("Error: %s has no column %s" % (fname, k))
    out, nmiss = [], 0
    for r in rows:
        r = list(r)
        ip = int(r[cols["iPoint"]])
        v = tab.pop(ip, None)
        if v is None:                                     # point of the run not in the table: not used
            r[cols["useNG"]] = r[cols["useRG"]] = "0"
            nmiss += 1
            out.append(r)
            continue
        if v[t["UTC"]] != r[cols["UTC"]]:
            sys.exit("Error: %s: point %d has UTC %s, the run has %s" % (fname, ip, v[t["UTC"]], r[cols["UTC"]]))
        trg, trg_run = float(v[t["T_RG_K"]]), float(r[cols["T_RG"]])
        if abs(trg - trg_run) > 1e-6 * trg_run:
            sys.exit("Error: %s: point %d has T_RG %s K, the run %s K (T_RG enters the RG kernels: rerun the base run)"
                     % (fname, ip, v[t["T_RG_K"]], r[cols["T_RG"]]))
        for g in ("NG", "RG"):
            use, p = int(v[t["use" + g]]), float(v[t["p%s_Pa" % g]])
            if use not in (0, 1):
                sys.exit("Error: %s: point %d: use%s must be 0 or 1" % (fname, ip, g))
            if use and not p > 0:
                sys.exit("Error: %s: point %d: use%s=1 but p%s=%s" % (fname, ip, g, g, v[t["p%s_Pa" % g]]))
            if use and int(r[cols["inside"]]):
                sys.exit("Error: %s: point %d: use%s=1 but the gauge is inside the shape model of the run" % (fname, ip, g))
            r[cols["ref" + g]] = v[t["p%s_Pa" % g]]
            r[cols["use" + g]] = str(use)
        out.append(r)
    if tab:
        sys.exit("Error: %s has %d points that are not in the run (e.g. %d)" % (fname, len(tab), min(tab)))
    print("fit data from %s (pressures in Pa): %d NG, %d RG points used%s" % (
        fname, sum(int(r[cols["useNG"]]) for r in out), sum(int(r[cols["useRG"]]) for r in out),
        "; %d points of the run not in the table (not used)" % nmiss if nmiss else ""))
    return out


def parse_utc(s):
    """seconds after 2016-01-01 of 'YYYY-MM-DDTh:mm:ss[.sss]' (hours not zero-padded in COPS_model.dat)"""
    d, t = s.strip().rstrip("Z").split("T")
    hh, mm, ss = t.split(":")
    return (dt.datetime.strptime(d, "%Y-%m-%d") - dt.datetime(2016, 1, 1)).total_seconds() + 3600 * int(hh) + 60 * int(mm) + float(ss)


def run_intervals(header, srcdir, run=None):
    """Read selection intervals from the archived input when available.

    Manifested runs carry ``input.snapshot.in``.  It is authoritative for the
    run configuration because the original path may later be edited or moved.
    Legacy products fall back to the filename recorded in the table header.
    """
    snapshot = os.path.join(run, "input.snapshot.in") if run else None
    if snapshot and os.path.isfile(snapshot):
        fname = snapshot
    else:
        m = re.search(r"input file: (\S+)", header[0])
        fname = m.group(1) if m else None
        if fname and not os.path.isabs(fname):
            fname = os.path.join(srcdir, fname)
    if not fname or not os.path.exists(fname):
        sys.exit("Error: the input file of the run (%s) is not found; it is needed for the ExcludeNG intervals" % fname)
    iv = {"ExcludeNG": [], "ExcludeRG": [], "FitInterval": []}
    for line in open(fname):
        mm = re.match(r"\s*(ExcludeNG|ExcludeRG|FitInterval)\s*=\s*(\S+)\s+(\S+)", line.split("#")[0])
        if mm:
            iv[mm.group(1)].append((parse_utc(mm.group(2)), parse_utc(mm.group(3)), mm.group(2), mm.group(3)))
    return fname, iv


def ng_saturated(hk_file, utc):
    """flags of the points (UTC of the run) whose NG reading is saturated: the COPS L2 BG record (table of
    tools/cops_hk_table.py) that ends at the time of the point (srcCG time = end of the BG interval, rounded to 1 s)
    has the NG reading constant at the maximum of the HIGH filament emission range (the definition of
    tools/collect_cops_sources.py; 2016-09-30 10:23:46.168-10:24:38.251, 2.2633e-8 mbar, emission current 0)"""
    names, rec = None, []
    for line in open(hk_file):
        if line.startswith("# utc_mid"):
            names = line[1:].split()
        elif line.strip() and not line.startswith("#"):
            rec.append(line.split())
    if names is None:
        sys.exit("Error: %s is not a table of tools/cops_hk_table.py" % hk_file)
    c = {n: i for i, n in enumerate(names)}
    pNG = np.array([float(r[c["pNG_mbar"]]) for r in rec])
    high = np.array([r[c["emisRangeFIL"]] == "HIGH" for r in rec])
    sat = high & (pNG > 0) & (pNG == np.max(np.where(high, pNG, -1.0)))
    stop = np.array([parse_utc(r[c["stop"]]) for r in rec])
    out = np.zeros(len(utc), dtype=bool)
    for k, u in enumerate(utc):
        t = parse_utc(u)
        j = np.argmin(np.abs(stop - t))
        out[k] = abs(stop[j] - t) < 1.0 and sat[j]
    return out


def header_value(header, key):
    for h in header:
        m = re.search(r"(?:^|\s)" + re.escape(key) + r"=(\S+)", h)
        if m:
            return m.group(1)
    sys.exit("Error: '%s' is not found in the header of COPS_model.dat" % key)


def sensitivities(header):
    for h in header:
        m = re.search(r"SpeciesSensitivity=(\S+) \(H2O ([0-9.eE+-]+), CO2 ([0-9.eE+-]+)\)", h)
        if m:
            on = m.group(1) == "on"
            return (float(m.group(2)) if on else 1.0, float(m.group(3)) if on else 1.0)
    sys.exit("Error: SpeciesSensitivity is not found in the header")


def read_activity_coefficients(srcdir, header):
    """H2O/CO2/CO coefficients used by the run: srcCG tables (SrcCGTables.h) or the ActivityFile."""
    model = header_value(header, "ActivityModel")
    if model == "file":
        fname = product_input(header_value(header, "ActivityFile"))
        txt = open(fname).read()
        txt = re.sub(r"#[^\n]*", "", txt).replace(",", " ")
        v = np.array([float(x) for x in txt.split()])
        return v[:75].reshape(3, 25)
    table = {"shap5": "Activity_SHAP5", "shap5_1": "Activity_SHAP5_1"}[model]
    txt = open(os.path.join(srcdir, "src", "SrcCGTables.h")).read()
    i = txt.index(table + "[3][25]=")
    body = txt[txt.index("{", i):txt.index("};", i)]
    v = np.array([float(x) for x in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", body)])
    return v.reshape(3, 25)


def read_elements(fname):
    """Load an already audited version-1 element cache.

    The complete cache is validated before this function is called.  The
    repeated exact-length checks below protect against replacement between the
    audit and load (and give this function safe standalone behavior in tests).
    """
    with open(fname, "rb") as f:
        header = f.read(12)
        if len(header) != 12:
            sys.exit("Error: truncated element header %s" % fname)
        version, n, nb = struct.unpack("<iii", header)
        if version != 1 or nb != NBASIS:
            sys.exit("Error: unsupported element file %s" % fname)
        raw = f.read()
        if len(raw) != n * 34 * 8:
            sys.exit("Error: element file %s has %d data bytes; expected %d" %
                     (fname, len(raw), n * 34 * 8))
        d = np.frombuffer(raw, dtype="<f8").reshape(n, 34)
    return {"centroid": d[:, 0:3], "normal": d[:, 3:6], "area": d[:, 6], "Y": np.ascontiguousarray(d[:, 7:32]),
            "AH2O": d[:, 32], "ACO2": d[:, 33]}


REC = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])


def read_kernel(fname):
    """Load one audited static kernel with exact record-length checks."""
    with open(fname, "rb") as f:
        prefix = f.read(8)
        if len(prefix) != 8:
            sys.exit("Error: truncated kernel header %s" % fname)
        version, ip = struct.unpack("<ii", prefix)
        if version != 1:
            sys.exit("Error: %s: kernel version %d; the fit needs version 1 (FlightTimeRotation=none)" % (fname, version))
        raw_header = f.read(80)
        raw_count = f.read(4)
        if len(raw_header) != 80 or len(raw_count) != 4:
            sys.exit("Error: truncated kernel header %s" % fname)
        h = struct.unpack("<10d", raw_header)
        (nk,) = struct.unpack("<i", raw_count)
        raw = f.read()
        if nk < 0 or len(raw) != nk * REC.itemsize:
            sys.exit("Error: kernel %s has %d data bytes for %d records; expected %d" %
                     (fname, len(raw), nk, max(nk, 0) * REC.itemsize))
        rec = np.frombuffer(raw, dtype=REC)
    return ip, {"et": h[0], "VdotB": h[1], "rh": h[8], "alt": h[9]}, rec


# ------------------------------------------------------------------------------------------- model
class Model:
    def __init__(self, run, header, cols, rows, elements, coef):
        self.run = run
        c = cols
        self.ipoint = np.array([int(r[c["iPoint"]]) for r in rows])
        self.useNG = np.array([int(r[c["useNG"]]) for r in rows])
        self.useRG = np.array([int(r[c["useRG"]]) for r in rows])
        self.refNG = np.array([float(r[c["refNG"]]) for r in rows])
        self.refRG = np.array([float(r[c["refRG"]]) for r in rows])
        self.TRG = np.array([float(r[c["T_RG"]]) for r in rows])
        self.pNG_run = np.array([float(r[c["pNG"]]) for r in rows])
        self.pRG_run = np.array([float(r[c["pRG"]]) for r in rows])

        if header_value(header, "TemperatureModel") != "analytical":
            sys.exit("Error: the run must use TemperatureModel=analytical")
        if header_value(header, "FlightTimeRotation") != "none":
            sys.exit("Error: the run must use FlightTimeRotation=none (version 1 kernels)")

        self.SZAMode = header_value(header, "SZAMode")
        self.helio = header_value(header, "HeliocentricModel")
        self.TNG = float(header_value(header, "TemperatureNG"))
        self.RGModel = header_value(header, "RGPressureModel")
        self.S = sensitivities(header)
        self.T0 = np.array([float(header_value(header, "Analytical" + k)) for k in TNAMES])
        self.coef = coef
        self.Y = elements["Y"]
        self.ACO2 = elements["ACO2"]
        self.nfaces = len(self.ACO2)

        # points used in the fit (and all points with kernels, for the output of the full model)
        self.sel = np.nonzero((self.useNG > 0) | (self.useRG > 0))[0]
        self.load_kernels(np.arange(len(rows)))

    def load_kernels(self, rowsel):
        pts, faces, kn, kphi, kv, cs, night, vb, rh, cnt = [], [], [], [], [], [], [], [], [], []
        for k, irow in enumerate(rowsel):
            fname = os.path.join(self.run, "kernels", "point=%04d.bin" % self.ipoint[irow])
            ip, h, rec = read_kernel(fname)
            if ip != self.ipoint[irow]:
                sys.exit("Error: %s stores iPoint %d; expected %d" %
                         (fname, ip, self.ipoint[irow]))
            if len(rec) and (np.min(rec["face"]) < 0 or np.max(rec["face"]) >= self.nfaces):
                sys.exit("Error: %s contains a face index outside [0,%d)" %
                         (fname, self.nfaces))
            n = len(rec)
            cnt.append(n)
            pts.append(np.full(n, k, dtype=np.int32))
            faces.append(rec["face"]), kn.append(rec["Kn"]), kphi.append(rec["KPhi"]), kv.append(rec["KV"])
            cs.append(rec["cos"]), night.append((rec["cos"] <= 0.0) | (rec["shadow"] != 0))
            vb.append(h["VdotB"]), rh.append(h["rh"])
        self.npt = len(rowsel)
        self.rows = rowsel
        self.pt = np.concatenate(pts)
        self.face = np.concatenate(faces)
        self.indptr = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
        self.Kn, self.KPhi, self.KV = np.concatenate(kn), np.concatenate(kphi), np.concatenate(kv)
        self.cos = np.concatenate(cs)
        self.night = np.concatenate(night)
        self.VdotB = np.array(vb)
        self.rh = np.array(rh)
        print("kernels: %d points, %d element contributions" % (self.npt, len(self.pt)))

        # production per unit activity (srcCG Comet::GetTotalProductionRateBjornNASTRAN) for each contribution
        rhA = self.rh / AU
        p = 4.2143229 if self.helio == "shap5" else 4.8
        c = np.minimum(self.cos, 1.0)
        if self.SZAMode == "continuous":
            cz = c
        else:
            ang = np.minimum(np.floor(np.degrees(np.arccos(c))), 89.0)
            cz = np.cos(np.radians(ang))
        self.Q = {}
        for s, (qmin, qmax) in (("H2O", (0.2 * 600 / rhA ** p, 600 / rhA ** p)), ("CO2", (0.1 * 600 / rhA ** 2, 600 / rhA ** 2))):
            qmin, qmax = qmin[self.pt], qmax[self.pt]
            self.Q[s] = np.where(self.night, qmin, qmin + (qmax - qmin) * cz)
        self.VB = self.VdotB[self.pt]
        self.prepare()

    def prepare(self):
        """parameter-independent factors of evaluate()"""
        c = np.where(self.night, -1.0, np.clip(self.cos, -1.0, 1.0))
        self.lit = c >= COS75                                         # analytical model: SZA <= 75 deg, lit
        self.xlit = np.where(self.lit, 1.0 / np.where(self.lit, c, 1.0) - 1.0, 0.0)
        self.bm = {s: np.sqrt(MASS[s] / (2.0 * KB)) * SQRTPI / 2.0 for s in SPECIES}
        self.cNG = {s: KB * self.TNG / S for s, S in zip(SPECIES, self.S)}
        self.cRG = {s: self.crg(s) / S for s, S in zip(SPECIES, self.S)}
        ACO2 = self.ACO2[self.face]
        q = {"H2O": self.Q["H2O"], "CO2": self.Q["CO2"] * ACO2}       # CO2: activity map of the run (scaled by kappa)
        self.wn = {s: q[s] * self.Kn for s in SPECIES}
        self.wv = {s: q[s] * self.VB * self.KV for s in SPECIES}
        self.wp = {"H2O": q["H2O"] * self.KPhi}
        self.PhiCO2 = np.bincount(self.pt, weights=q["CO2"] * self.KPhi, minlength=self.npt)

    def csr(self, data):
        """points x faces matrix with the per-contribution values 'data' (contributions are grouped by point)"""
        return sp.csr_matrix((data, self.face, self.indptr), shape=(self.npt, self.nfaces))

    def seen_faces(self):
        """elements that contribute to at least one loaded point"""
        m = np.zeros(self.nfaces, dtype=bool)
        m[self.face[(self.Kn > 0.0) | (self.KPhi > 0.0) | (self.KV > 0.0)]] = True
        return m

    def crg(self, s):
        m, kT = MASS[s], KB * self.TRG[self.rows]
        if self.RGModel == "equilibrium_sphere":
            return np.sqrt(2.0 * np.pi * m * kT)
        if self.RGModel == "calibration_srccg":
            return np.sqrt(np.pi * m * kT / 2.0)
        return kT * np.pi * np.sqrt(m / (2.0 * kT)) / 2.0

    def evaluate(self, a, T, kappa, jac=False, basis_index=None):
        """pNG, pRG for every loaded point; with jac=True also the derivatives with respect to
        a_j (j in basis_index), (Tmax, T75, Tmin) and kappa.
        beta sqrt(pi)/2 = bm_s u with u = T^-1/2; the parameter-independent factors are set in prepare()."""
        Tmax, T75, Tmin = T
        act_lin = self.Y @ a
        af = np.maximum(act_lin, 0.0)[self.face]
        Tl = Tmax + KT * (Tmax - T75) * self.xlit
        use_lit = self.lit & (Tl > Tmin)
        Te = np.where(use_lit, Tl, Tmin)
        u = 1.0 / np.sqrt(Te)
        bc = lambda w: np.bincount(self.pt, weights=w, minlength=self.npt)
        cNG, cRG, bm = self.cNG, self.cRG, self.bm

        pNG_H = cNG["H2O"] * bm["H2O"] * bc(self.wn["H2O"] * af * u)
        pRG_H = cRG["H2O"] * (bc(self.wp["H2O"] * af) + bm["H2O"] * bc(self.wv["H2O"] * af * u))
        pNG_C = cNG["CO2"] * bm["CO2"] * bc(self.wn["CO2"] * u)
        pRG_C = cRG["CO2"] * (self.PhiCO2 + bm["CO2"] * bc(self.wv["CO2"] * u))
        pNG = pNG_H + kappa * pNG_C
        pRG = pRG_H + kappa * pRG_C
        if not jac:
            return pNG, pRG

        J = {}
        # spherical-harmonic coefficients (H2O): d p / d a_j = sum_e (dp/dF_e) Q_e [act_lin_e > 0] Y_e,j
        B = (act_lin > 0.0)[:, None] * self.Y[:, basis_index]
        JaNG = cNG["H2O"] * bm["H2O"] * (self.csr(self.wn["H2O"] * u) @ B)
        JaRG = cRG["H2O"][:, None] * (self.csr(self.wp["H2O"] + bm["H2O"] * self.wv["H2O"] * u) @ B)
        for k, j in enumerate(basis_index):
            J[("a", j)] = (JaNG[:, k], JaRG[:, k])
        # temperature parameters: du/dT = -u/(2T); dT/dTmax = KT x + 1, dT/dT75 = -KT x (lit), dT/dTmin = 1 (otherwise)
        du = -0.5 * u / Te
        dl, dn_ = np.where(use_lit, du, 0.0), np.where(use_lit, 0.0, du)
        dlx = dl * self.xlit

        def sums(w):        # sum_e w (du/dT)(dT/dparam) for (Tmax, T75, Tmin)
            s_l, s_lx, s_n = bc(w * dl), bc(w * dlx), bc(w * dn_)
            return KT * s_lx + s_l, -KT * s_lx, s_n

        gN = cNG["H2O"] * bm["H2O"] * self.wn["H2O"] * af + kappa * cNG["CO2"] * bm["CO2"] * self.wn["CO2"]
        tN = sums(gN)
        tRH = sums(bm["H2O"] * self.wv["H2O"] * af)
        tRC = sums(bm["CO2"] * self.wv["CO2"])
        for k in range(3):
            J[("T", k)] = (tN[k], cRG["H2O"] * tRH[k] + kappa * cRG["CO2"] * tRC[k])
        # CO2 scale
        J[("kappa", 0)] = (pNG_C, pRG_C)
        return pNG, pRG, J


# ------------------------------------------------------------------------------------------- problem
class Problem:
    """parameter vector, residuals (data + Tikhonov) and Jacobian; the model is restricted to the fit points"""

    def __init__(self, M, el, coef, args):
        self.M, self.coef, self.args = M, coef, args
        self.fit_j = [j for j in range(NBASIS) if DEGREE_OF_BASIS[j] <= args.degree]
        self.ascale = np.max(np.abs(coef[0]))
        self.ordered = args.ordered_temperature and not args.fix_temperature

        names, x0, lo, hi = ["a%02d" % j for j in self.fit_j], [coef[0][j] / self.ascale for j in self.fit_j], [], []
        lo, hi = [-np.inf] * len(self.fit_j), [np.inf] * len(self.fit_j)
        if not args.fix_temperature:
            if self.ordered:
                Tmax, T75, Tmin = M.T0
                if not (Tmax >= T75 >= Tmin):
                    sys.exit("Error: --ordered-temperature needs Tmax >= T75 >= Tmin in the run")
                names += ["Tmin", "dT75", "dTmax"]          # Tmin, T75-Tmin, Tmax-T75
                x0 += [Tmin, T75 - Tmin, Tmax - T75]
                lo += [args.tbounds[0], 0.0, 0.0]
                hi += [args.tbounds[1], np.inf, np.inf]
            else:
                names += list(TNAMES)
                x0 += list(M.T0)
                lo += [args.tbounds[0]] * 3
                hi += [args.tbounds[1]] * 3
        if not args.fix_kappa:
            names += ["kappa"]
            x0 += [1.0]
            lo += [0.0]
            hi += [np.inf]
        self.names, self.x0, self.lo, self.hi = names, np.array(x0), np.array(lo), np.array(hi)
        self.na = len(self.fit_j)

        # Tikhonov operator: D^2 = (a-a0)_fit^T G (a-a0)_fit / N0, G = Y_fit^T diag(area) Y_fit = L L^T
        Yf = el["Y"][:, self.fit_j]
        G = Yf.T @ (el["area"][:, None] * Yf)
        self.N0 = np.sum(el["area"] * (el["Y"] @ coef[0]) ** 2)
        self.Lt = np.linalg.cholesky(G).T * self.ascale / np.sqrt(self.N0)   # r_reg = sqrt(lambda) Lt (x - x0)_a

        # data selection, equal weight per gauge (M.rows: the rows of the run whose kernels are loaded)
        sel = M.rows
        self.useNG = (M.useNG[sel] > 0) & (M.refNG[sel] > 0)
        self.useRG = (M.useRG[sel] > 0) & (M.refRG[sel] > 0)
        self.NNG, self.NRG = self.useNG.sum(), self.useRG.sum()
        self.wNG, self.wRG = 1.0 / np.sqrt(self.NNG), 1.0 / np.sqrt(self.NRG)
        self.lnNG, self.lnRG = np.log(M.refNG[sel][self.useNG]), np.log(M.refRG[sel][self.useRG])
        self.ndata = self.NNG + self.NRG
        self.lam = 0.0

    def unpack(self, x):
        a = self.coef[0].copy()
        a[self.fit_j] = x[:self.na] * self.ascale
        k = self.na
        T = self.M.T0.copy()
        if not self.args.fix_temperature:
            if self.ordered:
                Tmin, d1, d2 = x[k:k + 3]
                T = np.array([Tmin + d1 + d2, Tmin + d1, Tmin])
            else:
                T = x[k:k + 3].copy()
            k += 3
        kappa = 1.0 if self.args.fix_kappa else x[k]
        return a, T, kappa

    def regularization(self, x):
        return np.sqrt(self.lam) * (self.Lt @ (x[:self.na] - self.x0[:self.na]))

    def data_residuals(self, x):
        a, T, kappa = self.unpack(x)
        pNG, pRG = self.M.evaluate(a, T, kappa)
        return np.concatenate([self.wNG * (np.log(np.maximum(pNG[self.useNG], FLOOR)) - self.lnNG),
                               self.wRG * (np.log(np.maximum(pRG[self.useRG], FLOOR)) - self.lnRG)])

    def residuals(self, x):
        r = self.data_residuals(x)
        return np.concatenate([r, self.regularization(x)]) if self.lam > 0 else r

    def data_jacobian(self, x):
        a, T, kappa = self.unpack(x)
        pNG, pRG, J = self.M.evaluate(a, T, kappa, jac=True, basis_index=self.fit_j)
        pNG, pRG = np.maximum(pNG, FLOOR), np.maximum(pRG, FLOOR)
        col = lambda dn, dp: np.concatenate([self.wNG * dn[self.useNG] / pNG[self.useNG], self.wRG * dp[self.useRG] / pRG[self.useRG]])
        cols = [col(*J[("a", j)]) * self.ascale for j in self.fit_j]
        if not self.args.fix_temperature:
            JT = [col(*J[("T", k)]) for k in range(3)]       # d/d(Tmax, T75, Tmin)
            if self.ordered:                                  # x = (Tmin, T75-Tmin, Tmax-T75)
                JT = [JT[0] + JT[1] + JT[2], JT[0] + JT[1], JT[0]]
            cols += JT
        if not self.args.fix_kappa:
            cols.append(col(*J[("kappa", 0)]))
        return np.array(cols).T

    def reg_jacobian(self):
        R = np.zeros((self.na, len(self.x0)))
        R[:, :self.na] = np.sqrt(self.lam) * self.Lt
        return R

    def jacobian(self, x):
        Jd = self.data_jacobian(x)
        return np.vstack([Jd, self.reg_jacobian()]) if self.lam > 0 else Jd

    def deviation(self, x):
        """D: rms deviation of the linear H2O activity map from the prior, relative to the rms prior"""
        return np.linalg.norm(self.Lt @ (x[:self.na] - self.x0[:self.na]))


# ------------------------------------------------------------------------------------------- fit
P = None     # the problem (global so that forked worker processes share the kernels)


def fit_one(lam, x_start=None):
    P.lam = lam
    x_start = P.x0 if x_start is None else x_start
    res = least_squares(P.residuals, x_start, jac=P.jacobian, bounds=(P.lo, P.hi), x_scale="jac", method="trf",
                        max_nfev=P.args.max_nfev, ftol=P.args.ftol, xtol=P.args.xtol, gtol=P.args.gtol,
                        verbose=1 if P.args.jobs == 1 else 0)
    x = res.x
    rd = P.data_residuals(x)
    Jd = P.data_jacobian(x)
    R = P.reg_jacobian()
    H = Jd.T @ Jd + R.T @ R
    Hi = np.linalg.pinv(H, rcond=1e-14)
    peff = np.trace(Hi @ (Jd.T @ Jd))           # effective number of parameters resolved by the data
    v = np.sum(rd ** 2) / max(P.ndata - peff, 1.0)
    sig = np.sqrt(np.maximum(np.diag(v * Hi), 0.0))
    sv = np.linalg.svd(np.vstack([Jd, R]), compute_uv=False)
    return {"lam": lam, "x": x, "rd": rd, "cost": res.cost, "message": res.message,
            "success": bool(res.success), "nfev": res.nfev,
            "status": res.status, "optimality": res.optimality, "D": P.deviation(x), "peff": peff, "sig": sig, "sv": sv,
            "msNG": np.sum(rd[:P.NNG] ** 2), "msRG": np.sum(rd[P.NNG:] ** 2)}


def fit_worker(lam):
    r = fit_one(lam)
    print("lambda %.3e: cost %.4e  rms NG %.4f RG %.4f  D %.4f  p_eff %.2f  nfev %d  %s" % (
        lam, r["cost"], np.sqrt(r["msNG"]), np.sqrt(r["msRG"]), r["D"], r["peff"], r["nfev"], r["message"]), flush=True)
    return r


def activity_diagnostics(a, el, seen, coef):
    act = np.maximum(el["Y"] @ a, 0.0)
    act0 = np.maximum(el["Y"] @ coef[0], 0.0)
    area = el["area"]
    tot, tot0 = np.sum(area * act), np.sum(area * act0)
    return {"active_area": np.sum(area[act > 0]) / np.sum(area),
            "activity_ratio": tot / tot0,
            "activity_ratio_seen": np.sum((area * act)[seen]) / np.sum((area * act0)[seen]),
            "unseen_activity_fraction": np.sum((area * act)[~seen]) / tot if tot > 0 else np.nan,
            "unseen_activity_fraction_prior": np.sum((area * act0)[~seen]) / tot0,
            "unseen_area_fraction": np.sum(area[~seen]) / np.sum(area)}


def write_fit_data(out, P, rows, cols, msg):
    """the data that enter the misfit, taken from the fit problem itself (P.useNG/P.useRG, P.lnNG/P.lnRG over the
    loaded rows M.rows): table fit_data_used.txt and plot fit_data_used.png in <out>"""
    M = P.M
    used = {"NG": np.zeros(len(rows), dtype=bool), "RG": np.zeros(len(rows), dtype=bool)}
    pused = {"NG": np.full(len(rows), np.nan), "RG": np.full(len(rows), np.nan)}
    for g, mask, lnp in (("NG", P.useNG, P.lnNG), ("RG", P.useRG, P.lnRG)):
        k = M.rows[mask]
        used[g][k] = True
        pused[g][k] = np.exp(lnp)
    utc = [r[cols["UTC"]] for r in rows]
    alt = np.array([float(r[cols["altitude"]]) * (-1 if int(r[cols["inside"]]) else 1) for r in rows])
    ref = {g: np.array([float(r[cols["ref" + g]]) for r in rows]) for g in ("NG", "RG")}
    with open(os.path.join(out, "fit_data_used.txt"), "w") as f:
        f.write("# data used in the misfit of tools/fit_cops.py (taken from the fit problem: the points and values that enter the\n")
        f.write("# log residuals); %d NG, %d RG points; %s\n" % (P.NNG, P.NRG, msg))
        f.write("# pressures in Pa (N2-equivalent); p_used = -1: the point is not used for the gauge; altitude above the shape\n")
        f.write("# model of the run (negative: inside)\n")
        f.write("#iPoint UTC altitude_m pNG_used_Pa pRG_used_Pa refNG_Pa refRG_Pa\n")
        for k, r in enumerate(rows):
            f.write("%s %s %.3f %s %s %s %s\n" % (r[cols["iPoint"]], utc[k], alt[k],
                    "%.6e" % pused["NG"][k] if used["NG"][k] else "-1", "%.6e" % pused["RG"][k] if used["RG"][k] else "-1",
                    r[cols["refNG"]], r[cols["refRG"]]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    t = np.array([np.datetime64(dt.datetime(2016, 1, 1) + dt.timedelta(seconds=parse_utc(u))) for u in utc])
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.2), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    for g, mk, col in (("NG", "o", "#1f5fbf"), ("RG", "s", "#d62728")):
        u = used[g]
        nu = (~u) & (ref[g] > 0)
        for a, x in ((ax[0], t), (ax[1], t), (ax[2], alt)):
            a.plot(x[u], pused[g][u], mk, ms=3, color=col, label="%s used in the fit (%d)" % (g, u.sum()))
            a.plot(x[nu], ref[g][nu], mk, ms=4, mfc="none", mec="0.55", mew=0.8, label="%s data not used (%d)" % (g, nu.sum()))
    t0 = np.datetime64("2016-09-30T00:00")
    for a, (x0, x1) in zip(ax[:2], ((t0, t0 + np.timedelta64(645, "m")), (t0 + np.timedelta64(600, "m"), t0 + np.timedelta64(641, "m")))):
        a.set_xlim(x0, x1)
        a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        a.set_xlabel("UTC, 30 Sept 2016")
    ax[2].set_xscale("log")
    ax[2].set_xlim(1.0, 2.5e4)
    ax[2].set_xlabel("altitude above the shape model of the run [m]")
    for a in ax:
        a.set_yscale("log")
        a.set_ylabel("pressure [Pa]  (N2-equivalent)")
        a.grid(True, which="both", alpha=0.25)
    ax[0].legend(fontsize=8, loc="upper left", frameon=False)
    ax[0].set_title("data in the misfit of fit_cops.py: %d NG, %d RG points" % (P.NNG, P.NRG), fontsize=10)
    ax[1].set_title("last 40 minutes", fontsize=10)
    ax[2].set_title("vs. altitude; NG in ExcludeNG: %s; saturated NG: %s" % (
        "blocked" if P.args.apply_exclude_ng else "not blocked (default)",
        "excluded" if P.args.exclude_ng_saturated else "not removed (default)"), fontsize=9)
    fig.suptitle("%s   data: %s" % (M.run, os.path.abspath(P.args.data) if P.args.data else "the run (COPS_model.dat)"), fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(out, "fit_data_used.png"), dpi=110)
    plt.close(fig)

    # the used profiles only (lines break at unused points and at time gaps > 90 s)
    fig, ax = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1.2, 1]})
    sec = np.array([parse_utc(u) for u in utc])
    for g, ls, col in (("NG", "-", "#1f5fbf"), ("RG", "--", "#d62728")):
        k = np.nonzero(used[g])[0]
        cut = np.nonzero((np.diff(k) > 1) | (np.diff(sec[k]) > 90.0))[0] + 1
        for a, x in ((ax[0], t), (ax[1], alt)):
            for j, seg in enumerate(np.split(k, cut)):
                a.plot(x[seg], pused[g][seg], ls, marker=".", ms=3, lw=1.2, color=col,
                       label="%s (%d points)" % (g, len(k)) if j == 0 else None)
    ax[0].set_xlim(t0, t0 + np.timedelta64(645, "m"))
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[0].set_xlabel("UTC, 30 Sept 2016")
    ax[1].set_xscale("log")
    ax[1].set_xlim(1.0, 2.5e4)
    ax[1].set_xlabel("altitude above the shape model of the run [m]")
    for a in ax:
        a.set_yscale("log")
        a.set_ylabel("pressure [Pa]  (N2-equivalent)")
        a.grid(True, which="both", alpha=0.25)
        a.legend(fontsize=9, frameon=False)
    ax[0].set_title("COPS profiles used in the fit (only the points in the misfit)", fontsize=10)
    ax[1].set_title("vs. altitude; NG in ExcludeNG: %s; saturated NG: %s" % (
        "blocked" if P.args.apply_exclude_ng else "not blocked (default)",
        "excluded" if P.args.exclude_ng_saturated else "not removed (default)"), fontsize=10)
    fig.suptitle("%s   data: %s" % (M.run, os.path.abspath(P.args.data) if P.args.data else "the run (COPS_model.dat)"), fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(out, "fit_profiles_used.png"), dpi=110)
    plt.close(fig)
    print("fit data used: %s, %s, %s" % tuple(os.path.join(out, f) for f in ("fit_data_used.txt", "fit_data_used.png", "fit_profiles_used.png")))


def write_fit(out, r, run, header, names, rows, cols, srcdir, el, seen, r0):
    """Write one optimizer result without mislabeling inherited diagnostics.

    The cached Python model evaluates only the aggregate NG/RG pressures for a
    new parameter vector.  It does not reconstruct global production, the
    species-resolved table fields, or the below-deck diagnostic.  Those fields
    are therefore written as NaN, never copied from the parent run.  Running
    ``cops_fit.in`` through the C++ solver is the supported way to obtain them.
    """
    M, args = P.M, P.args
    a, T, kappa = P.unpack(r["x"])
    diag = activity_diagnostics(a, el, seen, P.coef)
    rmsNG, rmsRG = np.sqrt(r["msNG"]), np.sqrt(r["msRG"])

    with open(os.path.join(out, "fit_result.txt"), "w") as f:
        f.write("# LiouvilleCOPS fit of the source rate and temperature\n# run: %s\n" % run)
        f.write("# data: %d NG, %d RG points (useNG/useRG=1) from %s; misfit: log residuals, equal weight per gauge\n"
                % (P.NNG, P.NRG, os.path.abspath(args.data) if args.data else "the run (COPS_model.dat)"))
        f.write("# saturated NG points (10:23:50-10:24:38): %s\n" % ("excluded (--exclude-ng-saturated)" if args.exclude_ng_saturated
                                                                       else "not removed (default)"))
        f.write("# NG data in the ExcludeNG intervals of the run: %s\n" % ("blocked (--apply-exclude-ng)" if args.apply_exclude_ng
                                                                         else "not blocked (default)"))
        f.write("# fitted: spherical harmonics up to degree %d%s%s\n" % (args.degree, "" if args.fix_temperature else (
            ", Tmin T75-Tmin Tmax-T75 (ordered)" if P.ordered else ", Tmax T75 Tmin"), "" if args.fix_kappa else ", kappa"))
        f.write("# Tikhonov prior: H2O activity map of the run; D = rms(A-A_prior)/rms(A_prior) over the nucleus surface\n")
        f.write("lambda %.6e\n" % r["lam"])
        f.write("initial_cost %.6e\nfinal_cost %.6e\n" % (r0, r["cost"]))
        f.write("rms_ln_residual_NG %.4f  (typical factor %.3f)\nrms_ln_residual_RG %.4f  (typical factor %.3f)\n" % (rmsNG, np.exp(rmsNG), rmsRG, np.exp(rmsRG)))
        f.write("deviation_from_prior_D %.6f\n" % r["D"])
        f.write("effective_number_of_parameters %.3f (of %d)\n" % (r["peff"], len(P.x0)))
        f.write("optimizer: %s (status %d, nfev %d, first-order optimality %.3e)\n" % (r["message"], r["status"], r["nfev"], r["optimality"]))
        f.write("singular values of the weighted Jacobian incl. regularization (scaled parameters): %s\n" % " ".join("%.3e" % v for v in r["sv"]))
        f.write("condition number: %.3e\n" % (r["sv"][0] / r["sv"][-1] if r["sv"][-1] > 0 else np.inf))
        f.write("# activity diagnostics (H2O, A = max(0, sum a_j Y_j), area weighted; 'seen' = elements visible from a fit point)\n")
        for k, v in diag.items():
            f.write("%s %.6f\n" % (k, v))
        f.write("# parameter  initial  fitted  sigma (Gauss-Newton, posterior with the Tikhonov prior)\n")
        for k, n in enumerate(P.names):
            scale = P.ascale if n.startswith("a") else 1.0
            f.write("%-6s %+.6e %+.6e %.3e\n" % (n, P.x0[k] * scale, r["x"][k] * scale, r["sig"][k] * scale))
        f.write("Tmax %.6f\nT75 %.6f\nTmin %.6f\nkappa %.6f\n" % (T[0], T[1], T[2], kappa))

    with open(os.path.join(out, "ActivityFit.txt"), "w") as f:
        f.write("# activity coefficients (rows: H2O, CO2, CO), fitted by tools/fit_cops.py from %s, lambda=%.6e\n" % (run, r["lam"]))
        f.write("# H2O: fitted (degree <= %d); CO2: kappa=%.6f x run coefficients; CO: run coefficients\n" % (args.degree, kappa))
        for row in (a, kappa * P.coef[1], P.coef[2]):
            f.write(" ".join("%.15e" % v for v in row) + "\n")

    # input file that reruns the fitted model with LiouvilleCOPS
    snapshot = os.path.join(run, "input.snapshot.in")
    if os.path.isfile(snapshot):
        infile = snapshot
    else:
        infile = re.search(r"input file: (\S+)", header[0]).group(1)
        infile = infile if os.path.isabs(infile) else os.path.join(srcdir, infile)
    txt = open(infile).read()
    repl = {"ActivityModel": "file", "ActivityFile": os.path.join(out, "ActivityFit.txt"),
            "ActivityFileHeliocentricModel": M.helio, "AnalyticalTmax": "%.10g" % T[0], "AnalyticalT75": "%.10g" % T[1],
            "AnalyticalTmin": "%.10g" % T[2], "OutputDirectory": out + "_verify"}
    for k, v in repl.items():
        if re.search(r"^%s\s*=" % k, txt, flags=re.M):
            txt = re.sub(r"^%s\s*=.*$" % k, "%s = %s" % (k, v), txt, flags=re.M)
        else:
            txt += "%s = %s\n" % (k, v)
    with open(os.path.join(out, "cops_fit.in"), "w") as f:
        f.write("# generated by tools/fit_cops.py from %s\n" % infile + txt)

    # Fitted prediction table.  Geometry, observations, gauge temperature and
    # selection are inherited inputs; aggregate pressures are recomputed.  The
    # remaining source/transport diagnostics would be false if copied from the
    # parent base, so schema v2 marks them unavailable as NaN.
    pNG, pRG = MF.evaluate(a, T, kappa)
    with open(os.path.join(out, "COPS_model.dat"), "w") as f:
        f.write("# LiouvilleCOPS FIT PREDICTION (tools/fit_cops.py, lambda=%.3e) of the run %s\n" % (r["lam"], run))
        f.write("# ProductKind=fit_prediction TableSchema=liouville_cops_fit_prediction_v2 Optimizer=%s; "
                "authoritative process/cache/verification/review states are in fit_manifest.json\n" %
                ("converged" if r.get("success") else "not_converged"))
        f.write("# RecomputedColumns=pNG,pRG; UnavailableColumns=Q_H2O,Q_CO2,n_H2O,n_CO2,PhiRG_H2O,PhiRG_CO2,"
                "pNG_H2O,pNG_CO2,pRG_H2O,pRG_CO2,NGBelowDeckFraction (written as nan)\n")
        f.write("# Run ./liouville_cops cops_fit.in to obtain a complete forward table; this prediction table is not C++ verification.\n")
        for h in header[1:]:
            f.write("# " + h + "\n")
        f.write("#" + " ".join(names) + "\n")
        unavailable = ("Q_H2O", "Q_CO2", "n_H2O", "n_CO2", "PhiRG_H2O",
                       "PhiRG_CO2", "pNG_H2O", "pNG_CO2", "pRG_H2O",
                       "pRG_CO2", "NGBelowDeckFraction")
        for k, rr in enumerate(rows):
            rr = list(rr)
            for n in unavailable:
                rr[cols[n]] = "nan"
            if int(rr[cols["inside"]]) == 0:
                rr[cols["pNG"]], rr[cols["pRG"]] = "%.6e" % pNG[k], "%.6e" % pRG[k]
            else:
                rr[cols["pNG"]] = rr[cols["pRG"]] = "nan"
            f.write(" ".join(rr) + "\n")
    return diag


def max_curvature(lams, rho, eta):
    """lambda at the maximum curvature of the L-curve (log rho, log eta), parametrised by log lambda"""
    if len(lams) < 4:
        return None
    t = np.log(lams)
    x, y = np.log(rho), np.log(eta)
    dx, dy = np.gradient(x, t), np.gradient(y, t)
    ddx, ddy = np.gradient(dx, t), np.gradient(dy, t)
    kappa = (dx * ddy - dy * ddx) / np.maximum((dx ** 2 + dy ** 2) ** 1.5, 1e-300)
    k = np.argmax(kappa[1:-1]) + 1
    return lams[k], kappa


MF = None    # the model with the kernels of all points (output of the fitted model)


def optimizer_record(result):
    """JSON-safe optimizer state, separate from process/scientific status."""
    return {
        "lambda": float(result["lam"]), "success": bool(result["success"]),
        "status": int(result["status"]), "message": str(result["message"]),
        "nfev": int(result["nfev"]), "optimality": float(result["optimality"]),
        "cost": float(result["cost"]), "deviation_from_prior": float(result["D"]),
    }


def output_identities(directory, names):
    """Hash declared outputs only; undeclared files cannot acquire authority."""
    result = {}
    for name in names:
        path = Path(directory) / name
        if path.is_file():
            result[name] = file_identity(path)
    return result


def manifest_number(value):
    """Represent finite numbers as JSON numbers and open bounds explicitly."""
    value = float(value)
    if np.isposinf(value):
        return "+inf"
    if np.isneginf(value):
        return "-inf"
    if np.isnan(value):
        return "nan"
    return value


def main():
    global P, MF, ACTIVE_MANIFEST, ACTIVE_MANIFEST_PATH
    fit_started_utc = utc_now()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="LiouvilleCOPS run directory (FlightTimeRotation=none, WriteKernels=on)")
    ap.add_argument("--out", default=None, help="NEW output directory under cops_output or absolute path outside LiouvilleCOPS; existing paths are refused (default: <run-name>_fit)")
    ap.add_argument("--degree", type=int, default=4, help="max degree of the fitted spherical harmonics (0..4); higher degrees keep the run values")
    ap.add_argument("--lambda", dest="lam", type=float, nargs="+", default=[0.0],
                    help="Tikhonov weight(s) of the H2O activity prior; several values: L-curve")
    ap.add_argument("--jobs", type=int, default=None, help="parallel fits of the L-curve (default: min(#lambda, 12))")
    ap.add_argument("--fix-temperature", action="store_true", help="keep Tmax, T75, Tmin of the run")
    ap.add_argument("--ordered-temperature", action="store_true", help="fit Tmin, T75-Tmin>=0, Tmax-T75>=0 (Tmax>=T75>=Tmin)")
    ap.add_argument("--fix-kappa", action="store_true", help="keep the CO2 scale kappa=1")
    ap.add_argument("--tbounds", type=float, nargs=2, default=(30.0, 400.0), help="bounds of the temperature parameters [K]")
    ap.add_argument("--max-nfev", type=int, default=1000)
    ap.add_argument("--ftol", type=float, default=1e-10)
    ap.add_argument("--xtol", type=float, default=1e-10)
    ap.add_argument("--gtol", type=float, default=1e-10)
    ap.add_argument("--check-jacobian", action="store_true", help="compare the analytic Jacobian with finite differences and stop")
    ap.add_argument("--check-tolerance", type=float, default=1e-6, help="max relative difference allowed in the consistency check")
    ap.add_argument("--apply-exclude-ng", action="store_true",
                    help="block the NG data in the ExcludeNG intervals of the run's input file (fit_base.in: 10:24:40-10:40:00, "
                    "NG emission range HIGH->LOW until the impact); default: not blocked, these NG data are used")
    ap.add_argument("--exclude-ng-saturated", action="store_true",
                    help="do not use the saturated NG points (NG reading held at the maximum of the HIGH emission range, "
                    "2016-09-30 10:23:50-10:24:38, 13 points); default: they are used")
    ap.add_argument("--hk-table", default=None, help="COPS housekeeping table for --exclude-ng-saturated "
                    "(default: data/COPS_HK_20160930.txt of tools/cops_hk_table.py)")
    ap.add_argument("--data-only", action="store_true",
                    help="write and plot the data used in the fit (<out>/fit_data_used.txt, .png) and stop")
    ap.add_argument("--data", default=None, help="table of the fit data (pressures, useNG/useRG) instead of the values and flags "
                    "of the run, e.g. data/COPS_fit_data.txt of tools/collect_cops_sources.py")
    ap.add_argument("--allow-legacy-base", action="store_true",
                    help="explicitly accept a base without run_manifest.json after strict cache validation")
    args = ap.parse_args()
    lams = sorted(set(args.lam))
    if min(lams) < 0:
        sys.exit("Error: lambda must be >= 0")
    args.jobs = args.jobs or min(len(lams), 12)

    run = product_input(args.run.rstrip("/"))
    srcdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = str(output_path(args.out if args.out else os.path.basename(run) + "_fit"))
    if os.path.lexists(out):
        sys.exit("Error: output already exists; choose a new --out directory: %s" % out)

    # Validate every byte-level relationship before an optimizer or output
    # writer is started.  Old bases remain usable only by an explicit choice;
    # their absence of provenance is recorded rather than silently repaired.
    try:
        base_audit = validate_fit_base(run)
    except (OSError, ProductError) as exc:
        sys.exit("Error: invalid fitting base: %s" % exc)
    if base_audit["legacy_manifest"]:
        if not args.allow_legacy_base:
            sys.exit("Error: the base has no run_manifest.json. It passed no W01 provenance check; "
                     "rerun the base with the current solver or explicitly use --allow-legacy-base "
                     "to perform strict structural validation and record the limitation.")
        print("Warning: accepting an unmanifested legacy base after strict cache validation: %s" % run,
              file=sys.stderr)
    header, names, rows, cols = (base_audit[key] for key in
                                 ("headers", "names", "rows", "columns"))
    for k in ("useNG", "useRG", "T_RG"):
        if k not in cols:
            sys.exit("Error: COPS_model.dat has no column %s (rerun the model with the current code)" % k)

    coef = read_activity_coefficients(srcdir, header)
    el = read_elements(os.path.join(run, "elements.bin"))

    # the stored activity factors must equal max(0, sum a_j Y_j) of the coefficients of the run
    act = np.maximum(el["Y"] @ coef[0], 0.0)
    act_co2 = np.maximum(el["Y"] @ coef[1], 0.0)
    e = np.max(np.abs(act - el["AH2O"])) / max(np.max(np.abs(el["AH2O"])), FLOOR)
    e_co2 = np.max(np.abs(act_co2 - el["ACO2"])) / max(np.max(np.abs(el["ACO2"])), FLOOR)
    print("activity factors reproduced from the coefficients: max rel. diff H2O %.2e, CO2 %.2e" %
          (e, e_co2))
    if max(e, e_co2) > 1e-10:
        sys.exit("Error: the coefficients do not match the element file")
    if not base_audit["legacy_manifest"] and header_value(header, "ActivityModel") == "file":
        activity = [entry.get("file") for entry in base_audit["manifest"].get("dependencies", [])
                    if entry.get("role") == "activity_coefficients"]
        if len(activity) != 1:
            sys.exit("Error: manifested ActivityModel=file base has no unique activity-coefficient identity")
        try:
            verify_identity(activity[0], "base activity coefficients")
        except ProductError as exc:
            sys.exit("Error: %s" % exc)

    MF = Model(run, header, cols, rows, el, coef)

    # consistency check: the model with the parameters of the run vs. the LiouvilleCOPS pressures
    pNG0, pRG0 = MF.evaluate(coef[0], MF.T0, 1.0)
    ok = (MF.pNG_run > 0) & (MF.pRG_run > 0)
    dNG = np.max(np.abs(pNG0[ok] / MF.pNG_run[ok] - 1.0))
    dRG = np.max(np.abs(pRG0[ok] / MF.pRG_run[ok] - 1.0))
    print("consistency with LiouvilleCOPS (%d points): max rel. diff NG %.2e, RG %.2e" % (ok.sum(), dNG, dRG))
    if max(dNG, dRG) > args.check_tolerance:
        sys.exit("Error: the fit model does not reproduce the LiouvilleCOPS pressures")

    # fit data: the reference values and flags of the run, or the table --data
    if args.data:
        args.data = product_input(args.data)
        rows = apply_data_file(args.data, header, cols, rows)
        header = [header[0]] + ["fit data: %s" % args.data] + header[1:]

    # NG data in the ExcludeNG intervals of the run (fit_base.in: 10:24:40-10:40:00, NG emission range HIGH->LOW
    # until the impact): not blocked by default; --apply-exclude-ng blocks them
    infile, iv = run_intervals(header, srcdir, run)
    tt = np.array([parse_utc(r[cols["UTC"]]) for r in rows])
    inex = np.zeros(len(rows), dtype=bool)
    for a, b, _, _ in iv["ExcludeNG"]:
        inex |= (tt >= a) & (tt <= b)
    infit = np.ones(len(rows), dtype=bool)
    if iv["FitInterval"]:
        infit = np.zeros(len(rows), dtype=bool)
        for a, b, _, _ in iv["FitInterval"]:
            infit |= (tt >= a) & (tt <= b)
    exs = ", ".join("%s..%s" % (x[2], x[3]) for x in iv["ExcludeNG"]) or "none"
    rows = [list(r) for r in rows]
    n = 0
    if args.apply_exclude_ng:
        for k in np.nonzero(inex)[0]:
            n += int(rows[k][cols["useNG"]])
            rows[k][cols["useNG"]] = "0"
        msg_ex = "NG data in the ExcludeNG intervals of %s (%s) blocked (--apply-exclude-ng)%s" % (
            os.path.basename(infile), exs, "; %d points of the table removed" % n if args.data else "")
    else:
        if not args.data:              # the run flags have the intervals applied: unblock (data > 0, outside the nucleus)
            for k in np.nonzero(inex & infit)[0]:
                r = rows[k]
                if float(r[cols["refNG"]]) > 0 and not int(r[cols["inside"]]) and not int(r[cols["useNG"]]):
                    r[cols["useNG"]] = "1"
                    n += 1
        msg_ex = "NG data in the ExcludeNG intervals of %s (%s) not blocked (default; --apply-exclude-ng blocks them)%s" % (
            os.path.basename(infile), exs, "; %d points added to the selection of the run" % n if not args.data else
            "; the table is used as given")
    print(msg_ex)
    header = [header[0]] + [msg_ex] + header[1:]

    hk = product_input(args.hk_table or "data/COPS_HK_20160930.txt")
    sat = ng_saturated(hk, [r[cols["UTC"]] for r in rows])
    nsat = int(np.sum([sat[k] and int(r[cols["useNG"]]) for k, r in enumerate(rows)]))
    if args.exclude_ng_saturated:
        rows = [list(r) for r in rows]
        for k in np.nonzero(sat)[0]:
            rows[k][cols["useNG"]] = "0"
        msg = "saturated NG points excluded (--exclude-ng-saturated): %d of the %d flagged in %s were selected" % (nsat, sat.sum(), hk)
    else:
        msg = ("saturated NG points not removed (default, no --exclude-ng-saturated): %d of the %d flagged in %s are in the fit data"
               % (nsat, sat.sum(), hk))
    print(msg)
    header = [header[0]] + [msg] + header[1:]

    # model restricted to the fit points (used by the optimisation)
    M = Model(run, header, cols, rows, el, coef)
    sel = M.sel
    M.load_kernels(sel)
    P = Problem(M, el, coef, args)
    seen = M.seen_faces()
    print("fit data: %d NG points, %d RG points; %d parameters: %s" % (P.NNG, P.NRG, len(P.x0), " ".join(P.names)))

    # Create the experiment directory only after all base/cache/data checks
    # succeeded.  Existing directories are refused above, so this operation
    # cannot merge a new experiment with stale or historical artifacts.
    os.makedirs(out, exist_ok=False)
    chosen_rows = {
        "point_ids": [int(r[cols["iPoint"]]) for r in rows],
        "utc": [r[cols["UTC"]] for r in rows],
        "ref_ng": [float(r[cols["refNG"]]) for r in rows],
        "ref_rg": [float(r[cols["refRG"]]) for r in rows],
        "use_ng": [int(r[cols["useNG"]]) for r in rows],
        "use_rg": [int(r[cols["useRG"]]) for r in rows],
        "t_rg": [float(r[cols["T_RG"]]) for r in rows],
    }
    selection = selection_identity(
        chosen_rows["point_ids"], chosen_rows["utc"], chosen_rows["ref_ng"],
        chosen_rows["ref_rg"], chosen_rows["use_ng"], chosen_rows["use_rg"],
        chosen_rows["t_rg"])
    input_files = {
        "base_model_table": base_audit["model_table"],
        "base_elements": base_audit["elements"],
        "base_input": file_identity(infile),
        "housekeeping_table": file_identity(hk),
    }
    if args.data:
        input_files["replacement_fit_data"] = file_identity(args.data)
    run_manifest_path = Path(run) / "run_manifest.json"
    parent = ({"kind": "manifested_run", "manifest": file_identity(run_manifest_path)}
              if run_manifest_path.is_file() else
              {"kind": "legacy_unmanifested_run", "path": str(Path(run).resolve()),
               "accepted_explicitly": True})
    ACTIVE_MANIFEST = {
        "schema": FIT_SCHEMA,
        "schema_version": 1,
        "created_utc": fit_started_utc,
        "updated_utc": utc_now(),
        "fit_kind": "data_selection_only" if args.data_only else
                    ("single" if len(lams) == 1 else "lambda_scan"),
        "command": [sys.executable] + sys.argv,
        "cwd": str(Path.cwd().resolve()),
        "parent": parent,
        "states": {
            "process": "running",
            "optimizer": "not_run",
            "cache": "verified",
            "base_python_cpp_agreement": "verified",
            "cpp_prediction_agreement": "not_run",
            "scientific_review": "not_reviewed",
        },
        "base_cache": {
            "legacy_manifest": bool(base_audit["legacy_manifest"]),
            "kernels": base_audit["kernels"],
            "python_cpp_max_relative_difference": {"NG": float(dNG), "RG": float(dRG)},
            "consistency_tolerance": float(args.check_tolerance),
        },
        "inputs": input_files,
        "selection": selection,
        "prior": {
            "activity_coefficients_sha256": canonical_sha256(coef.tolist()),
            "activity_coefficients_shape": [3, 25],
            "temperature": {name: float(value) for name, value in zip(TNAMES, MF.T0)},
            "co2_scale": 1.0,
        },
        "objective": {
            "name": "equal_gauge_normalized_log_residuals_plus_h2o_map_tikhonov",
            "pressure_floor_Pa": FLOOR,
            "lambda_values": [float(value) for value in lams],
            "selection_policy": {
                "apply_exclude_ng": bool(args.apply_exclude_ng),
                "exclude_ng_saturated": bool(args.exclude_ng_saturated),
                "replacement_data": bool(args.data),
            },
        },
        "parameters": {
            "names": list(P.names), "degree": int(args.degree),
            "fix_temperature": bool(args.fix_temperature),
            "ordered_temperature": bool(args.ordered_temperature),
            "fix_kappa": bool(args.fix_kappa),
            "lower_bounds": [manifest_number(value) for value in P.lo],
            "upper_bounds": [manifest_number(value) for value in P.hi],
            "initial_scaled": [float(value) for value in P.x0],
        },
        "optimizer_settings": {
            "method": "scipy.optimize.least_squares/trf", "x_scale": "jac",
            "max_nfev": int(args.max_nfev), "ftol": float(args.ftol),
            "xtol": float(args.xtol), "gtol": float(args.gtol),
            "jobs": int(args.jobs),
        },
        "software": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
            "environment": {name: os.environ[name] for name in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
                            if name in os.environ},
            "fit_cops": file_identity(__file__),
            "cops_io": file_identity(Path(__file__).with_name("cops_io.py")),
            "srcCG_activity_tables": file_identity(Path(srcdir) / "src" / "SrcCGTables.h"),
        },
        "output_schema": {
            "name": "liouville_cops_fit_prediction_v2",
            "aggregate_pressures": "recomputed",
            "source_and_species_diagnostics": "unavailable_nan_until_cpp_verification",
        },
        "optimizer_results": [],
        "outputs": {},
    }
    ACTIVE_MANIFEST_PATH = str(Path(out) / "fit_manifest.json")
    save_active_manifest()
    write_fit_data(out, P, rows, cols, msg)
    if args.data_only:
        ACTIVE_MANIFEST["states"]["process"] = "complete"
        ACTIVE_MANIFEST["states"]["optimizer"] = "not_requested"
        ACTIVE_MANIFEST["updated_utc"] = utc_now()
        ACTIVE_MANIFEST["outputs"] = output_identities(
            out, ("fit_data_used.txt", "fit_data_used.png", "fit_profiles_used.png"))
        save_active_manifest()
        return
    print("elements seen from the fit points: %.1f%% of the surface area" % (100 * np.sum(el["area"][seen]) / np.sum(el["area"])))

    if args.check_jacobian:
        P.lam = max(lams)
        x = P.x0
        Ja = P.jacobian(x)
        worst = 0.0
        for k, n in enumerate(P.names):
            h = 1e-6 * max(abs(x[k]), 1.0)
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            fd = (P.residuals(xp) - P.residuals(xm)) / (2 * h)
            err = np.max(np.abs(fd - Ja[:, k])) / max(np.max(np.abs(fd)), 1e-300)
            worst = max(worst, err)
            print("  d/d%-6s analytic vs finite difference: max rel. diff %.2e" % (n, err))
        print("Jacobian check (lambda %.3e): worst %.2e" % (P.lam, worst))
        ACTIVE_MANIFEST["states"]["process"] = "complete"
        ACTIVE_MANIFEST["states"]["optimizer"] = "not_requested"
        ACTIVE_MANIFEST["jacobian_check"] = {"lambda": float(P.lam),
                                              "worst_relative_difference": float(worst)}
        ACTIVE_MANIFEST["updated_utc"] = utc_now()
        ACTIVE_MANIFEST["outputs"] = output_identities(
            out, ("fit_data_used.txt", "fit_data_used.png", "fit_profiles_used.png"))
        save_active_manifest()
        return

    rd0 = P.data_residuals(P.x0)
    r0 = 0.5 * np.sum(rd0 ** 2)
    print("initial cost %.6e  (rms ln residual NG %.3f, RG %.3f)" % (r0, np.sqrt(np.sum(rd0[:P.NNG] ** 2)), np.sqrt(np.sum(rd0[P.NNG:] ** 2))))

    if len(lams) == 1:
        r = fit_worker(lams[0])
        write_fit(out, r, run, header, names, rows, cols, srcdir, el, seen, r0)
        ACTIVE_MANIFEST["optimizer_results"] = [optimizer_record(r)]
        ACTIVE_MANIFEST["states"]["process"] = "complete"
        ACTIVE_MANIFEST["states"]["optimizer"] = "converged" if r["success"] else "not_converged"
        ACTIVE_MANIFEST["updated_utc"] = utc_now()
        ACTIVE_MANIFEST["outputs"] = output_identities(
            out, ("fit_data_used.txt", "fit_data_used.png", "fit_profiles_used.png",
                  "fit_result.txt", "ActivityFit.txt", "cops_fit.in", "COPS_model.dat"))
        save_active_manifest()
        print("results: %s" % out)
        print("verify with: ./liouville_cops %s" % os.path.join(out, "cops_fit.in"))
        return

    # L-curve: independent fits from the prior (forked processes share the kernels)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(args.jobs) as pool:
        results = pool.map(fit_worker, lams, chunksize=1)

    table = []
    result_outputs = {}
    for r in results:
        d = os.path.join(out, "lambda=%.3e" % r["lam"])
        os.makedirs(d, exist_ok=False)
        diag = write_fit(d, r, run, header, names, rows, cols, srcdir, el, seen, r0)
        result_outputs[os.path.basename(d)] = output_identities(
            d, ("fit_result.txt", "ActivityFit.txt", "cops_fit.in", "COPS_model.dat"))
        a, T, kappa = P.unpack(r["x"])
        table.append((r, diag, T, kappa))

    pos = [k for k, (r, _, _, _) in enumerate(table) if r["lam"] > 0 and r["D"] > 0]
    lam_c = None
    if len(pos) >= 4:
        lp = np.array([table[k][0]["lam"] for k in pos])
        rho = np.array([np.sqrt(table[k][0]["msNG"] + table[k][0]["msRG"]) for k in pos])
        eta = np.array([table[k][0]["D"] for k in pos])
        lam_c, curv = max_curvature(lp, rho, eta)

    with open(os.path.join(out, "lcurve.txt"), "w") as f:
        f.write("# L-curve of the Tikhonov-regularised fit; run %s\n" % run)
        f.write("# rho = sqrt(msNG + msRG) (data misfit, ln units), D = relative deviation of the H2O activity map from the prior\n")
        if lam_c is not None:
            f.write("# maximum curvature of (ln rho, ln D) at lambda = %.3e\n" % lam_c)
        f.write("#lambda rmsNG rmsRG rho D p_eff Tmax T75 Tmin kappa activity_ratio activity_ratio_seen unseen_activity_fraction active_area nfev status optimality\n")
        for r, dg, T, kappa in table:
            f.write("%.4e %.5f %.5f %.5f %.5e %.3f %.3f %.3f %.3f %.5f %.4f %.4f %.4f %.4f %d %d %.2e\n" % (
                r["lam"], np.sqrt(r["msNG"]), np.sqrt(r["msRG"]), np.sqrt(r["msNG"] + r["msRG"]), r["D"], r["peff"],
                T[0], T[1], T[2], kappa, dg["activity_ratio"], dg["activity_ratio_seen"], dg["unseen_activity_fraction"],
                dg["active_area"], r["nfev"], r["status"], r["optimality"]))
    print(open(os.path.join(out, "lcurve.txt")).read())
    ACTIVE_MANIFEST["optimizer_results"] = [optimizer_record(result) for result in results]
    if all(result["success"] for result in results):
        optimizer_state = "converged"
    elif any(result["success"] for result in results):
        optimizer_state = "mixed"
    else:
        optimizer_state = "not_converged"
    ACTIVE_MANIFEST["states"]["process"] = "complete"
    ACTIVE_MANIFEST["states"]["optimizer"] = optimizer_state
    ACTIVE_MANIFEST["updated_utc"] = utc_now()
    ACTIVE_MANIFEST["outputs"] = output_identities(
        out, ("fit_data_used.txt", "fit_data_used.png", "fit_profiles_used.png", "lcurve.txt"))
    ACTIVE_MANIFEST["outputs"]["lambda_directories"] = result_outputs
    save_active_manifest()
    print("results: %s" % out)


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        # A failure after the output directory is created must remain visibly
        # incomplete.  Preserve SystemExit's concise CLI message and keyboard
        # interrupt semantics after recording the independent process state.
        if ACTIVE_MANIFEST is not None and ACTIVE_MANIFEST.get("states", {}).get("process") == "running":
            ACTIVE_MANIFEST["states"]["process"] = "failed"
            ACTIVE_MANIFEST["updated_utc"] = utc_now()
            ACTIVE_MANIFEST["error"] = str(exc) or exc.__class__.__name__
            save_active_manifest()
        raise
