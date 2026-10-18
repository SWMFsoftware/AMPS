#!/usr/bin/env python3
"""Collect the end-of-mission (2016-09-30) COPS nude gauge (NG) and ram gauge (RG) data of all available sources
into one machine-readable JSON file, and plot them in one unit (Pa) vs. time and vs. altitude. It also writes
the data actually used in the fit of tools/fit_cops.py as a table (input of fit_cops.py --data) and plots them.

DATA SOURCES (dataset id in the JSON file: origin; local input; processing)
  PSA_L2_BG      ROSINA team, ESA PSA / PDS3 data set RO-C-ROSINA-2-EXT3-V2.0, COPS BG housekeeping files
                 (DATA/MTP35/COPS/BG), CODMAC level 2: NG and RG pressures [mbar] calibrated on board for N2 at 20 C,
                 one value per file interval (60 s, 4 s near the surface); sensor temperature, emission range and
                 autorange flags.
                 https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/ROSINA/RO-C-ROSINA-2-EXT3-V2.0/
                 Local generated table: <AMPS>/cops_output/data/COPS_HK_20160930.txt
                 (tools/cops_hk_table.py from the 1155 BG files under <COPS_DATA_ROOT>/data_COPS_L2). Time: middle of the interval.
                 Pressure 0 (no valid value) -> null. 'ng_saturated' flags the NG readings that are constant at the
                 maximum of the HIGH emission range (10:23:46-10:24:38).
  fit_data       the data used in the fit: refNG, refRG (= srcCG_table) and the flags useNG/useRG of the base run of
                 the fit (--fit-run, default output_fit_base; selection by LiouvilleCOPS: reference value > 0, gauge
                 outside the shape model, not in ExcludeNG/ExcludeRG, inside FitInterval of the run's input file).
                 The reason of each exclusion is derived from the input file of the run and checked against the
                 flags (the flags are authoritative; 'unknown' if they disagree).
                 ExcludeNG: the NG data in the ExcludeNG intervals of the run's input file (fit_base.in:
                 10:24:40-10:40:00, NG emission range HIGH->LOW until the impact) are NOT blocked by default: the
                 table has useNG = 1 for them (if the data > 0 and the gauge is outside the shape model).
                 --apply-exclude-ng keeps the flags of the run (useNG = 0, reason ExcludeNG:<begin>/<end>).
                 Saturated NG points (column ng_saturated): the points whose L2 BG record (the record that ends at the
                 time of the point) has ng_saturated, i.e. the 13 points 10:23:50 .. 10:24:38 (NG reading held at
                 2.2633e-8 mbar, filament emission current 0). By default they are INCLUDED (useNG as in the run);
                 --exclude-ng-saturated sets useNG = 0 for them in the table (reason NG_saturated).
  srcCG_table    AMPS srcCG, RosinaMeasurements.h (NudeGaugeReferenceData, RamGaugeReferenceData [mbar], 1150 points),
                 the reference data of the srcCG model. Identical to PSA_L2_BG; time = END of the BG interval.
                 Local: read from the COPS_model.dat of the SHAP5 LiouvilleCOPS run (refNG, refRG in Pa; LiouvilleCOPS
                 parses RosinaMeasurements.h with ReferenceDataUnit = mbar).
  PSA_L4         ROSINA team, ESA PSA / PDS3 data set RO-C-ROSINA-4-EXT3-V2.0, DATA/MTP35/COPS/COPS_L4_MTP35.ASC
                 (+ .LBL, LABEL/COPS_TS_TABLE.FMT), CODMAC level 4: NG_LOCAL_DENSITY [m^-3] and its uncertainty
                 (20%), derived from the NG pressure with the ionisation cross sections of the DFMS composition
                 (ROSINA Data Users Guide sect. 5.4), 2016-09-26T13:31 .. 2016-09-30T10:36:41; also DIST_COMET [km].
                 https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/ROSINA/RO-C-ROSINA-4-EXT3-V2.0/
                 Local: <COPS_DATA_ROOT>/sources/COPS_L4_MTP35.ASC; only the rows of 2016-09-30 are used. No data from 10:24:33
                 to 10:33:05; the last row (10:36:41) is ~3 min before the impact (the Users Guide, p. 38, advises
                 to ignore the large densities late in the mission close to 67P, where saturation may have occurred).
  ESA_blog_2016  ESA Rosetta blog 2016-12-15, "Rosetta's last words: science descending to a comet", COPS pressure
                 figure (credit K. Altwegg), https://blogs.esa.int/rosetta/files/2016/12/COPS_EoM.jpg.
                 Local image: <COPS_DATA_ROOT>/sources; generated CSV: <AMPS>/cops_output/data/sources.
                 Made from the data of 2016: it has NG values in the gap of the archive (10:24:42-10:33:14) and is
                 ~230 s earlier than PSA_L2_BG (fit of a time shift, 06:00-10:18, rms 0.044 in ln p). Not corrected.
  EGU2019_s02    V. Tenishev et al., "End-of-Mission ROSINA/COPS observation of the innermost coma of comet
  EGU2019_s16    67P/Churyumov-Gerasimenko", EGU General Assembly 2019 (~/Downloads/2019-EGU.pptx), slide 2 and
                 slide 16 (left, "Observed RG and NG pressure"): NG and RG vs. "Altitude" (definition not given).
                 Local images: <COPS_DATA_ROOT>/sources; generated CSVs: <AMPS>/cops_output/data/sources.
                 (tools/digitize_cops_figures.py). No time; the altitude is kept as plotted. Slides 2 and 16 are two
                 different figures (see tools/digitize_cops_figures.py): they coincide above ~4 km and differ below
                 ~3 km. Slide 16 agrees with PSA_L2_BG vs. the M004 altitude and its NG follows the ESA blog curve
                 through the archive gap. Slide 2 at 0.5-2 km has the archive values of 370-520 m larger M004
                 altitude; below ~300 m its NG is one straight segment (2.3e-8 to 8e-8 mbar). The altitude of
                 slide 2 is not stated.
  Not included (listed in the JSON under "not_included"): Heritier et al. 2017, MNRAS 469, S118 (COPS density of
  the final descent as figures only; full text not accessible here), Tzou 2017 PhD thesis (not accessible),
  PSA L2 SN/SR science files (end 10:19:26), PSA L3/L5 EXT3 (no COPS products for MTP35).

GEOMETRY of the time-tagged datasets
  Signed altitude above the SHAP5 (SHAP5_Volume.bdf) and the ESA DSK M004 (ROS_CG_M004_OSPGDLR_N_V1.BDS) surfaces
  (distance to the nearest surface point; negative inside) and the distance from the nucleus centre, linearly
  interpolated in time from the COPS_model.dat of two LiouvilleCOPS runs (SPICE: spacecraft position in
  67P/C-G_CK, kernels of cops.in). Outside the time span of the runs: null.

CONVERSIONS (ROSINA Data Users Guide RO-ROS-UG-1000, sect. 5.2)
  p [Pa] = 100 p [mbar]; all pressures are N2-equivalent gauge pressures (20 C)
  NG:  nNG_N2eq = p / (k T_NG), T_NG = 293 K, k = 1.380649e-23 J/K
       nNG_H2O  = beta_H2O nNG_N2eq, beta_H2O = 0.893 (pure H2O)
  RG:  PhiRG_N2eq = p / sqrt(2 pi m k T_RG)  (= 1/4 sqrt(8kT/(pi m)) p/(kT)), m = 18.01528 amu, T_RG = sensor temp.
       PhiRG_H2O  = beta_H2O PhiRG_N2eq
  PSA_L4 in the plots: pNG_equiv = n k 293 K (n is composition corrected; this is not the gauge reading)

OUTPUT
  data/COPS_EoM_sources.json
    title, created, generator
    conventions     units and the formulae above
    geometry_runs   the LiouvilleCOPS runs used for the geometry (directory, header lines)
    events          NG emission range switch, NG saturation, NG gap, last L4 row, last L2 record, touchdown
    datasets[]      id, kind, gauges, source {citation, url, local, ...}, notes[], columns {name: [values]}
                    (column-oriented; equal lengths within a dataset except the EGU datasets, which have separate
                    altitude grids for NG and RG; null = missing)
    not_included[]  citation, reason
  data/COPS_EoM_sources_time.png      NG (solid) and RG (dashed) pressure vs. UTC: whole descent and last 40 min
                                      (red filled markers: points used in the fit)
  data/COPS_EoM_sources_altitude.png  NG and RG pressure vs. altitude (M004 for the time-tagged data, as plotted for EGU)

  data/COPS_fit_data.txt (--fit-data)  the fit data, one row per point of the fit run:
                                      iPoint UTC altitude_SHAP5_m pNG_Pa pRG_Pa useNG useRG T_RG_K ng_saturated
                                      reason_NG reason_RG (pressures as in COPS_model.dat, Pa, 0 = no data; reason:
                                      used | no_data | inside_shape_model | NG_saturated | Exclude<NG|RG>:<begin>/<end> |
                                      outside_FitInterval | unknown).
                                      Edit useNG/useRG or the pressures and pass it to tools/fit_cops.py --data.
  data/COPS_fit_data.png              the fit data vs. UTC (whole descent, last 40 min) and vs. SHAP5 altitude; used
                                      points filled, excluded points open and coloured by reason; Exclude intervals
                                      and the saturated NG points shaded (included or excluded, see above)

  data/COPS_fit_profiles.png          the NG (solid) and RG (dashed) profiles actually used in a fit, vs. UTC (whole
                                      descent, last 40 min) and vs. M004 altitude: only the points in the misfit, read
                                      from fit_data_used.txt of tools/fit_cops.py (--fit-used, default
                                      output_fit_base_selection/fit_data_used.txt, the default selection written by
                                      fit_cops.py output_fit_base --data-only --out output_fit_base_selection;
                                      skipped if the file does not exist)

usage: collect_cops_sources.py [--shap5-run DIR] [--m004-run DIR] [--fit-run DIR] [--fit-data FILE]
                               [--apply-exclude-ng] [--exclude-ng-saturated] [--fit-used FILE] [--out FILE.json] [--plots PREFIX] [--show]
  defaults: --shap5-run output_fit_base, --m004-run output_M004 (LiouvilleCOPS runs with COPS_model.dat),
            --fit-run output_fit_base, --fit-data data/COPS_fit_data.txt
"""
import argparse
import csv
import datetime as dt
import json
import os
import re

from cops_paths import prepare_output, product_input, sources_root

import numpy as np
import matplotlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--shap5-run", default="output_fit_base")
ap.add_argument("--m004-run", default="output_M004")
ap.add_argument("--sources-dir", default=str(sources_root()),
                help="external L4/image source directory (default: <COPS_DATA_ROOT>/sources; override with COPS_SOURCES_ROOT)")
ap.add_argument("--out", default="data/COPS_EoM_sources.json", help="JSON under cops_output, or absolute path outside LiouvilleCOPS")
ap.add_argument("--plots", default="data/COPS_EoM_sources", help="PNG prefix under cops_output, or absolute path outside LiouvilleCOPS")
ap.add_argument("--show", action="store_true", help="open the plots in windows")
ap.add_argument("--fit-run", default="output_fit_base",
                help="base run of the fit (tools/fit_cops.py); its useNG/useRG flags select the fit data")
ap.add_argument("--fit-data", default="data/COPS_fit_data.txt",
                help="table of the fit data (input of tools/fit_cops.py --data)")
ap.add_argument("--fit-used", default="output_fit_base_selection/fit_data_used.txt",
                help="fit_data_used.txt of a fit (tools/fit_cops.py) for the plot of the profiles used in that fit")
ap.add_argument("--apply-exclude-ng", action="store_true",
                help="block the NG data in the ExcludeNG intervals of the fit run's input file in the fit data table "
                     "(useNG as in the run); default: not blocked (useNG=1 there)")
ap.add_argument("--exclude-ng-saturated", action="store_true",
                help="set useNG=0 in the fit data table for the saturated NG points (default: they are included)")
args = ap.parse_args()
for key in ("shap5_run", "m004_run", "fit_run", "fit_used"):
    setattr(args, key, product_input(getattr(args, key)))
if not args.show:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

KB = 1.380649e-23
AMU = 1.66053906660e-27
M_H2O = 18.01528 * AMU
T_NG = 293.0
BETA_H2O = 0.893
T0 = dt.datetime(2016, 9, 30)
SRC = os.path.abspath(os.path.expanduser(args.sources_dir))


def check_required_inputs():
    """Report every missing prerequisite before parsing data or creating outputs."""
    required = []
    run_roles = {}
    for role, directory in (("SHAP5 geometry (--shap5-run)", args.shap5_run),
                            ("M004 geometry (--m004-run)", args.m004_run),
                            ("fit selection (--fit-run)", args.fit_run)):
        model = os.path.join(directory, "COPS_model.dat")
        run_roles.setdefault(model, []).append(role)
    for model, roles in run_roles.items():
        required.append((model, ", ".join(roles)))

    required.append((product_input("data/COPS_HK_20160930.txt"), "processed COPS L2 housekeeping"))
    required.append((os.path.join(SRC, "COPS_L4_MTP35.ASC"), "external COPS L4 table"))
    for name in ("digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv",
                 "digitized_EGU2019_Tenishev_slide02_NG.csv",
                 "digitized_EGU2019_Tenishev_slide02_RG.csv",
                 "digitized_EGU2019_Tenishev_slide16_data_NG.csv",
                 "digitized_EGU2019_Tenishev_slide16_data_RG.csv"):
        required.append((product_input(os.path.join("data", "sources", name)), "digitized source curve"))

    missing = [(path, role) for path, role in required if not os.path.isfile(path)]
    if not missing:
        return
    lines = ["Error: collect_cops_sources.py cannot run because required inputs are missing:"]
    lines.extend("  - %s [%s]" % (path, role) for path, role in missing)
    lines.extend((
        "",
        "Generate the default prerequisites from the LiouvilleCOPS directory:",
        "  python3 tools/cops_hk_table.py",
        "  python3 tools/digitize_cops_figures.py",
        "  env OMP_NUM_THREADS=4 ./liouville_cops fit_base.in",
        "  env OMP_NUM_THREADS=4 ./liouville_cops m004.in",
        "  python3 tools/collect_cops_sources.py --show",
        "",
        "Alternatively, point --shap5-run, --m004-run, and --fit-run at existing compatible runs.",
    ))
    raise SystemExit("\n".join(lines))


check_required_inputs()
for key in ("out", "plots", "fit_data"):
    setattr(args, key, prepare_output(getattr(args, key)))


def parse_utc(s):
    s = s.strip().rstrip("Z")
    d, t = s.split("T")
    hh, mm, ss = t.split(":")
    sec = float(ss)
    return dt.datetime.strptime(d, "%Y-%m-%d") + dt.timedelta(hours=int(hh), minutes=int(mm), seconds=sec)


def tsec(d):
    return (d - T0).total_seconds()


def iso(d):
    return d.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]


def nn(v, fmt=None):
    """JSON value: None for nan"""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    return float("%.6g" % v) if fmt is None else float(fmt % v)


# ------------------------------------------------------------------ geometry from the LiouvilleCOPS runs
def read_run(d):
    names, rows, header = None, [], []
    for line in open(os.path.join(d, "COPS_model.dat")):
        if line.startswith("#iPoint"):
            names = line[1:].split()
        elif line.startswith("#"):
            header.append(line[1:].strip())
        elif line.strip():
            rows.append(line.split())
    c = {n: i for i, n in enumerate(names)}
    t = np.array([tsec(parse_utc(r[c["UTC"]])) for r in rows])
    alt = np.array([float(r[c["altitude"]]) for r in rows])
    inside = np.array([int(r[c["inside"]]) for r in rows])
    xyz = np.array([[float(r[c[k]]) for k in ("x", "y", "z")] for r in rows])
    shape = [h for h in header if h.startswith("ShapeModel=")][0].split()[0].split("=", 1)[1]
    out = {"t": t, "alt": np.where(inside > 0, -alt, alt), "r": np.linalg.norm(xyz, axis=1) / 1e3,
           "utc": [r[c["UTC"]] for r in rows], "refNG": np.array([float(r[c["refNG"]]) for r in rows]),
           "refRG": np.array([float(r[c["refRG"]]) for r in rows]), "shape": shape, "header": header,
           "inside": inside, "dir": d}
    for k in ("iPoint", "useNG", "useRG"):
        out[k] = np.array([int(r[c[k]]) for r in rows]) if k in c else None
    out["T_RG"] = np.array([float(r[c["T_RG"]]) for r in rows]) if "T_RG" in c else None
    out["raw"] = {k: [r[c[k]] for r in rows] for k in ("refNG", "refRG", "T_RG") if k in c}
    return out


S5 = read_run(args.shap5_run)
M4 = read_run(args.m004_run)
FR = S5 if os.path.abspath(args.fit_run) == os.path.abspath(args.shap5_run) else read_run(args.fit_run)
if FR["useNG"] is None or FR["T_RG"] is None:
    raise SystemExit("Error: %s/COPS_model.dat has no useNG/useRG/T_RG columns (rerun it with the current code)" % args.fit_run)


def geom(t):
    """signed altitude (SHAP5, M004) [m] and cometocentric distance [km] at times t [s after 2016-09-30T00:00]"""
    t = np.asarray(t, float)
    out = {}
    for key, g in (("altitude_SHAP5_m", S5), ("altitude_M004_m", M4)):
        a = np.interp(t, g["t"], g["alt"])
        a[(t < g["t"][0]) | (t > g["t"][-1])] = np.nan
        out[key] = a
    r = np.interp(t, M4["t"], M4["r"])
    r[(t < M4["t"][0]) | (t > M4["t"][-1])] = np.nan
    out["r_km"] = r
    return out


# ------------------------------------------------------------------ datasets
datasets = []

# PSA L2 BG housekeeping
hk_path = product_input("data/COPS_HK_20160930.txt")
hk = [l.split() for l in open(hk_path) if not l.startswith("#") and l.strip()]
hk_hdr = [l for l in open(hk_path) if l.startswith("# utc_mid")][0][1:].split()
hc = {n: i for i, n in enumerate(hk_hdr)}
tm = np.array([tsec(parse_utc(r[hc["utc_mid"]])) for r in hk])
pNG = np.array([float(r[hc["pNG_mbar"]]) for r in hk]) * 100.0
pRG = np.array([float(r[hc["pRG_mbar"]]) for r in hk]) * 100.0
TRG = np.array([float(r[hc["Tsensor_K"]]) for r in hk])
pNG = np.where(pNG > 0, pNG, np.nan)
pRG = np.where(pRG > 0, pRG, np.nan)
# NG saturated: reading constant at the largest NG value of the HIGH emission range
hi = np.array([r[hc["emisRangeFIL"]] == "HIGH" for r in hk])
pmax_hi = np.nanmax(np.where(hi, pNG, np.nan))
sat = hi & (pNG == pmax_hi)
g = geom(tm)
datasets.append({
    "id": "PSA_L2_BG", "kind": "archive measurement", "gauges": ["NG", "RG"],
    "source": {"citation": "ROSINA team (K. Altwegg), PSA dataset RO-C-ROSINA-2-EXT3-V2.0, COPS BG housekeeping files "
                           "(DATA/MTP35/COPS/BG), product creation 2019",
               "url": "https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/ROSINA/RO-C-ROSINA-2-EXT3-V2.0/",
               "local": "data/COPS_HK_20160930.txt (tools/cops_hk_table.py from the L2 BG files)"},
    "notes": ["NG/RG pressure = ROSINA_COPS_PRESSURE_NG/RG (N2-equivalent, 20 C); 0 in the archive (no valid value) -> null",
              "time = middle of the BG file interval; the srcCG table uses the end of the interval",
              "NG emission range HIGH -> LOW at 10:24:40; NG = 0 from 10:24:42 to 10:33:14",
              "ng_saturated: NG reading constant at the maximum of the HIGH emission range (10:23:46-10:24:38)",
              "COPS team (Users Guide p. 38): ignore the large densities late in the mission close to 67P (saturation)"],
    "columns": {
        "utc": [iso(T0 + dt.timedelta(seconds=float(s))) for s in tm],
        "utc_start": [r[hc["start"]] for r in hk], "utc_stop": [r[hc["stop"]] for r in hk],
        "altitude_SHAP5_m": [nn(v) for v in g["altitude_SHAP5_m"]], "altitude_M004_m": [nn(v) for v in g["altitude_M004_m"]],
        "r_km": [nn(v) for v in g["r_km"]],
        "pNG_Pa": [nn(v) for v in pNG], "pRG_Pa": [nn(v) for v in pRG], "T_RG_K": [nn(v) for v in TRG],
        "nNG_N2eq_m3": [nn(v) for v in pNG / (KB * T_NG)], "nNG_H2O_m3": [nn(v) for v in BETA_H2O * pNG / (KB * T_NG)],
        "PhiRG_N2eq_m2s": [nn(v) for v in pRG / np.sqrt(2 * np.pi * M_H2O * KB * TRG)],
        "PhiRG_H2O_m2s": [nn(v) for v in BETA_H2O * pRG / np.sqrt(2 * np.pi * M_H2O * KB * TRG)],
        "emission_range_NG": [r[hc["emisRangeFIL"]] for r in hk], "emission_range_RG": [r[hc["emisRangeMIC"]] for r in hk],
        "autorange_NG": [r[hc["autorangeNG"]] for r in hk], "ng_saturated": [bool(v) for v in sat]},
})

# srcCG reference table (as parsed by LiouvilleCOPS; COPS_model.dat stores it in Pa)
ts = S5["t"]
g = geom(ts)
datasets.append({
    "id": "srcCG_table", "kind": "model input table (copy of the archive values)", "gauges": ["NG", "RG"],
    "source": {"citation": "AMPS srcCG/RosinaMeasurements.h, NudeGaugeReferenceData / RamGaugeReferenceData (mbar)",
               "local": "srcCG/RosinaMeasurements.h; values via %s/COPS_model.dat" % os.path.relpath(args.shap5_run, ROOT)},
    "notes": ["identical to PSA_L2_BG; time = end of the BG file interval", "0 in the table -> null",
              "altitude_SHAP5_m < 0: the point is inside the SHAP5 surface (points 1146-1149)"],
    "columns": {"utc": S5["utc"], "altitude_SHAP5_m": [nn(v) for v in g["altitude_SHAP5_m"]],
                "altitude_M004_m": [nn(v) for v in g["altitude_M004_m"]], "r_km": [nn(v) for v in g["r_km"]],
                "pNG_Pa": [nn(v if v > 0 else np.nan) for v in S5["refNG"]],
                "pRG_Pa": [nn(v if v > 0 else np.nan) for v in S5["refRG"]]},
})

# data used in the fit: the reference values and the useNG/useRG flags of the base run of the fit
def input_intervals(run):
    """ExcludeNG, ExcludeRG, FitInterval of the input file of a run (for the exclusion reasons)"""
    m = re.search(r"input file: (\S+)", run["header"][0])
    fname = m.group(1) if m else None
    if fname and not os.path.isabs(fname):
        fname = os.path.join(ROOT, fname)
    iv = {"ExcludeNG": [], "ExcludeRG": [], "FitInterval": []}
    if fname and os.path.exists(fname):
        for line in open(fname):
            line = line.split("#")[0]
            mm = re.match(r"\s*(ExcludeNG|ExcludeRG|FitInterval)\s*=\s*(\S+)\s+(\S+)", line)
            if mm:
                iv[mm.group(1)].append((mm.group(2), mm.group(3)))
    return fname, iv


fit_input, FIV = input_intervals(FR)


def fit_reason(k, gauge):
    """why point k of the fit run is not used for the gauge ('' = used)"""
    ref = FR["ref" + gauge][k]
    t = FR["t"][k]
    inside_iv = lambda a, b: tsec(parse_utc(a)) <= t <= tsec(parse_utc(b))
    if ref <= 0:
        return "no_data"
    if FR["inside"][k]:
        return "inside_shape_model"
    for a, b in FIV["Exclude" + gauge]:
        if inside_iv(a, b):
            return "Exclude%s:%s/%s" % (gauge, a, b)
    if FIV["FitInterval"] and not any(inside_iv(a, b) for a, b in FIV["FitInterval"]):
        return "outside_FitInterval"
    return ""


reason = {gg: [fit_reason(k, gg) for k in range(len(FR["t"]))] for gg in ("NG", "RG")}
for gg in ("NG", "RG"):
    for k, rr in enumerate(reason[gg]):
        if (rr == "") != (FR["use" + gg][k] == 1):        # the flags of the run are authoritative
            reason[gg][k] = "used" if FR["use" + gg][k] else "unknown"
        elif rr == "":
            reason[gg][k] = "used"
# saturated NG points of the fit run: the L2 BG record that ends at the time of the point (srcCG time = end of the
# interval, rounded to 1 s) has ng_saturated (NG reading held at the maximum of the HIGH emission range)
hk_stop = np.array([tsec(parse_utc(r[hc["stop"]])) for r in hk])
ngsat = np.zeros(len(FR["t"]), dtype=bool)
for k, t in enumerate(FR["t"]):
    j = np.argmin(np.abs(hk_stop - t))
    ngsat[k] = abs(hk_stop[j] - t) < 1.0 and sat[j]
FR["useNG"] = FR["useNG"].copy()
# ExcludeNG intervals (fit_base.in: 10:24:40-10:40:00): not blocked by default -> the points excluded only by them are used
nunblock = 0
if not args.apply_exclude_ng:
    for k, rr in enumerate(reason["NG"]):
        if rr.startswith("ExcludeNG:"):
            FR["useNG"][k] = 1
            reason["NG"][k] = "used"
            nunblock += 1
ex_note = ("NG data in the ExcludeNG intervals of %s: %s" % (fit_input and os.path.basename(fit_input),
           "BLOCKED (--apply-exclude-ng; useNG as in the run)" if args.apply_exclude_ng else
           "NOT blocked (default; %d points with useNG=0 in the run set to useNG=1; --apply-exclude-ng blocks them)" % nunblock))
if args.exclude_ng_saturated:
    for k in np.nonzero(ngsat & (FR["useNG"] == 1))[0]:
        FR["useNG"][k] = 0
        reason["NG"][k] = "NG_saturated"
iem = np.array([float(hk[j][hc["IemisFIL_uA"]]) for j in np.nonzero(sat)[0]])
sat_note = ("NG saturated (ng_saturated = 1): %d points, %s .. %s (L2 records %s .. %s: NG reading constant at %.5e mbar, "
            "filament emission current %g..%g uA); %s" % (ngsat.sum(), FR["utc"][np.nonzero(ngsat)[0][0]], FR["utc"][np.nonzero(ngsat)[0][-1]],
            hk[np.nonzero(sat)[0][0]][hc["start"]], hk[np.nonzero(sat)[0][-1]][hc["stop"]], pmax_hi / 100.0, iem.min(), iem.max(),
            "EXCLUDED from the fit data (--exclude-ng-saturated)" if args.exclude_ng_saturated else
            "INCLUDED in the fit data (default; --exclude-ng-saturated sets useNG=0)"))
nfit = {gg: int(np.sum(FR["use" + gg])) for gg in ("NG", "RG")}
g = geom(FR["t"])
fit_note = ("selection of the fit run %s (input file %s): a point is used when the reference value is > 0, the gauge "
            "is outside the shape model, the point is not in Exclude<NG|RG> and inside FitInterval; intervals: %s"
            % (os.path.relpath(FR["dir"], ROOT), fit_input and os.path.relpath(fit_input, ROOT),
               "; ".join("%s %s %s" % (k, a, b) for k, v in FIV.items() for a, b in v) or "none"))
datasets.append({
    "id": "fit_data", "kind": "data used in the fit (tools/fit_cops.py)", "gauges": ["NG", "RG"],
    "source": {"citation": "srcCG reference table (= PSA_L2_BG) with the fit selection of LiouvilleCOPS",
               "local": "%s/COPS_model.dat; table %s" % (os.path.relpath(FR["dir"], ROOT), os.path.relpath(args.fit_data, ROOT))},
    "notes": [fit_note, ex_note, sat_note, "%d NG and %d RG points used" % (nfit["NG"], nfit["RG"]),
              "pNG_Pa/pRG_Pa: all points of the run (null = no data); useNG/useRG = 1: used in the fit",
              "the misfit is the log residual with equal weight per gauge; T_RG_K is the COPS sensor temperature "
              "used for the RG model"],
    "columns": {"iPoint": [int(v) for v in FR["iPoint"]], "utc": FR["utc"],
                "altitude_SHAP5_m": [nn(v) for v in g["altitude_SHAP5_m"]], "altitude_M004_m": [nn(v) for v in g["altitude_M004_m"]],
                "pNG_Pa": [nn(v if v > 0 else np.nan) for v in FR["refNG"]],
                "pRG_Pa": [nn(v if v > 0 else np.nan) for v in FR["refRG"]],
                "T_RG_K": [nn(v) for v in FR["T_RG"]],
                "useNG": [int(v) for v in FR["useNG"]], "useRG": [int(v) for v in FR["useRG"]],
                "ng_saturated": [int(v) for v in ngsat], "reason_NG": reason["NG"], "reason_RG": reason["RG"]},
})

with open(args.fit_data, "w") as f:
    f.write("# COPS NG and RG pressures used in the fit of the source rate and temperature (tools/fit_cops.py --data)\n")
    f.write("# written by tools/collect_cops_sources.py on %s from %s/COPS_model.dat\n" % (dt.date.today().isoformat(), FR["dir"]))
    f.write("# data: srcCG/RosinaMeasurements.h (= ROSINA COPS L2 BG pressures, PSA RO-C-ROSINA-2-EXT3-V2.0), N2-equivalent\n")
    f.write("#   gauge pressures (20 C) in Pa; time = end of the L2 BG file interval; 0 = no data\n")
    f.write("# %s\n" % fit_note)
    f.write("# %s\n" % ex_note)
    f.write("# %s\n" % sat_note)
    f.write("# used: %d NG, %d RG points of %d\n" % (nfit["NG"], nfit["RG"], len(FR["t"])))
    f.write("# use the table with: python3 tools/fit_cops.py %s --data %s\n" % (os.path.relpath(FR["dir"], ROOT), os.path.relpath(args.fit_data, ROOT)))
    f.write("#   the iPoint/UTC must be those of the run (kernels); edit useNG/useRG (0/1) or the pressures to change the fit data;\n")
    f.write("#   T_RG_K must equal the T_RG of the run (it is used in the RG kernels); altitude_SHAP5_m is for information\n")
    f.write("#   (negative: inside the shape model; such points cannot be used)\n")
    f.write("# ng_saturated: 1 = NG reading saturated (L2 BG housekeeping); fit_cops.py --exclude-ng-saturated excludes these\n")
    f.write("#   points also when useNG = 1 here\n")
    f.write("# reason_NG/reason_RG (information; not read by the fit): used | no_data | inside_shape_model | NG_saturated |\n")
    f.write("#   Exclude<NG|RG>:<UTC begin>/<UTC end> | outside_FitInterval | unknown (flag of the run not explained by its input file)\n")
    f.write("#iPoint UTC altitude_SHAP5_m pNG_Pa pRG_Pa useNG useRG T_RG_K ng_saturated reason_NG reason_RG\n")
    for k in range(len(FR["t"])):
        f.write("%d %s %.3f %s %s %d %d %s %d %s %s\n" % (FR["iPoint"][k], FR["utc"][k], FR["alt"][k], FR["raw"]["refNG"][k],
                FR["raw"]["refRG"][k], FR["useNG"][k], FR["useRG"][k], FR["raw"]["T_RG"][k], ngsat[k], reason["NG"][k], reason["RG"][k]))
print("fit data: %s (%d NG, %d RG points used; %d saturated NG points %s)" % (args.fit_data, nfit["NG"], nfit["RG"], ngsat.sum(),
      "excluded" if args.exclude_ng_saturated else "included"))

# PSA L4 NG density
l4 = [l.split(",") for l in open(os.path.join(SRC, "COPS_L4_MTP35.ASC")) if l.startswith("2016-09-30")]
t4 = np.array([tsec(parse_utc(r[0])) for r in l4])
n4 = np.array([float(r[7]) for r in l4])
s4 = np.array([float(r[8]) for r in l4])
g = geom(t4)
datasets.append({
    "id": "PSA_L4", "kind": "archive measurement (derived density)", "gauges": ["NG"],
    "source": {"citation": "ROSINA team, PSA dataset RO-C-ROSINA-4-EXT3-V2.0, DATA/MTP35/COPS/COPS_L4_MTP35.ASC "
                           "(NG_LOCAL_DENSITY, composition-corrected with DFMS abundances)",
               "url": "https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/ROSINA/RO-C-ROSINA-4-EXT3-V2.0/DATA/MTP35/COPS/",
               "local": os.path.join(SRC, "COPS_L4_MTP35.ASC") + " (+ .LBL, COPS_TS_TABLE.FMT)"},
    "notes": ["only 2016-09-30 rows; the file starts 2016-09-26", "no data 10:24:33-10:33:05; last row 10:36:41",
              "pNG_equiv_Pa = n k 293 K (for comparison with the pressures; not the gauge reading)",
              "r_km_L4 = DIST_COMET of the file"],
    "columns": {"utc": [r[0].strip() for r in l4], "altitude_SHAP5_m": [nn(v) for v in g["altitude_SHAP5_m"]],
                "altitude_M004_m": [nn(v) for v in g["altitude_M004_m"]], "r_km": [nn(v) for v in g["r_km"]],
                "r_km_L4": [float(r[2]) for r in l4],
                "nNG_m3": [nn(v) for v in n4], "nNG_sigma_m3": [nn(v) for v in s4],
                "pNG_equiv_Pa": [nn(v) for v in n4 * KB * T_NG]},
})


def read_csv(fname):
    rows = [r for r in csv.reader(l for l in open(product_input(os.path.join("data", "sources", fname))) if not l.startswith("#"))]
    return rows[0], rows[1:]


# ESA blog figure (digitized)
h, rows = read_csv("digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv")
tb = np.array([tsec(parse_utc(r[0])) for r in rows])
pb = np.array([float(r[1]) for r in rows]) * 100.0
g = geom(tb)
datasets.append({
    "id": "ESA_blog_2016", "kind": "digitized figure", "gauges": ["NG"],
    "source": {"citation": "ESA Rosetta blog, 'Rosetta's last words: science descending to a comet', 2016-12-15; "
                           "figure 'COPS measured gas pressure ... final hours' (credit K. Altwegg)",
               "url": "https://blogs.esa.int/rosetta/2016/12/15/rosettas-last-words-science-descending-to-a-comet/",
               "image": "https://blogs.esa.int/rosetta/files/2016/12/COPS_EoM.jpg",
               "local": os.path.join(SRC, "ESA_blog_2016-12-15_COPS_EoM.jpg") +
                        "; generated data/sources/digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv"},
    "notes": ["'COPS pressure' in mbar vs. UTC; the gauge is not named in the figure (NG: matches PSA_L2_BG NG)",
              "digitized by tools/digitize_cops_figures.py (dot centroids); resolution 0.002 decades, ~2 min in time",
              "values above the frame top (~2e-7 mbar) are clipped in the figure",
              "contains NG values in 10:24:42-10:33:14 where the archive (V2.0) has none (preliminary 2016 data)",
              "time offset vs. PSA_L2_BG (06:00-10:18): the blog curve is ~230 s earlier (rms 0.044 in ln p); not corrected",
              "at the end the curve is nearly vertical: time uncertainty ~ +-2 min (dot width 4 min)"],
    "columns": {"utc": [r[0] for r in rows], "altitude_SHAP5_m": [nn(v) for v in g["altitude_SHAP5_m"]],
                "altitude_M004_m": [nn(v) for v in g["altitude_M004_m"]], "r_km": [nn(v) for v in g["r_km"]],
                "pNG_Pa": [nn(v) for v in pb], "nNG_N2eq_m3": [nn(v) for v in pb / (KB * T_NG)]},
})

# EGU 2019 figures (digitized; altitude as plotted, no time)
for key, fname, slide, note in (
        ("EGU2019_s02", "EGU2019_Tenishev_slide02", 2, "slide text: NG and RG profiles 'flip' at ~360 m; a different "
         "figure from slide 16: equal above ~4 km, below ~3 km NG lower by up to a factor 1.9 and RG by 8-15%; "
         "at 0.5-2 km it has the archive values of 370-520 m larger M004 altitude; altitude not stated"),
        ("EGU2019_s16", "EGU2019_Tenishev_slide16_data", 16, "left figure 'Observed RG and NG pressure'; NG/RG cross at ~700 m; agrees with PSA_L2_BG vs. the "
         "M004 altitude (RG within 1-3%; NG above 700 m), NG follows ESA_blog_2016 through the archive gap")):
    cols = {}
    for gauge in ("NG", "RG"):
        h, rows = read_csv("digitized_%s_%s.csv" % (fname, gauge))
        cols["altitude_%s_m" % gauge] = [float(r[0]) for r in rows]
        cols["p%s_Pa" % gauge] = [nn(float(r[1]) * 100.0) for r in rows]
        if gauge == "NG":
            cols["nNG_N2eq_m3"] = [nn(float(r[1]) * 100.0 / (KB * T_NG)) for r in rows]
    datasets.append({
        "id": key, "kind": "digitized figure", "gauges": ["NG", "RG"],
        "source": {"citation": "V. Tenishev et al., 'End-of-Mission ROSINA/COPS observation of the innermost coma of comet "
                               "67P/Churyumov-Gerasimenko', EGU General Assembly 2019, slide %d" % slide,
                   "local": "~/Downloads/2019-EGU.pptx; %s; generated data/sources/digitized_%s_{NG,RG}.csv" %
                            (os.path.join(SRC, fname + ".jpg"), fname)},
        "notes": [note, "pressure axis 'Pressure' without unit; values equal the mbar values of the archive",
                  "altitude as plotted (definition not given); NG and RG have separate altitude grids",
                  "digitized dashed curves (median row per pixel column); resolution ~0.005 decades"],
        "columns": cols,
    })

not_included = [
    {"citation": "K. L. Heritier et al. 2017, MNRAS 469, S118, doi:10.1093/mnras/stx1459",
     "reason": "COPS NG density of the final descent only as figures (Figs. 4, 6, 8); full text not accessible from this "
               "machine (bot protection). Text: ~5e9 cm^-3 just before impact, +-15%, pure H2O assumed."},
    {"citation": "C.-Y. Tzou 2017, PhD thesis, Univ. Bern, 'Calibrations of ROSINA-COPS and observations at comet 67P'",
     "reason": "not accessible (bot protection); content on the final descent unknown"},
    {"citation": "PSA RO-C-ROSINA-2-EXT3-V2.0 COPS SN/SR science files", "reason": "last file ends 2016-09-30T10:19:26"},
    {"citation": "PSA RO-C-ROSINA-3-EXT3 / RO-C-ROSINA-5-EXT3", "reason": "no COPS products for MTP35"},
]

doc = {
    "title": "Rosetta ROSINA/COPS end-of-mission (2016-09-30) nude and ram gauge data from all available sources",
    "created": dt.date.today().isoformat(), "generator": "LiouvilleCOPS/tools/collect_cops_sources.py",
    "conventions": {
        "time": "UTC (spacecraft), ISO 8601", "pressure": "Pa (1 mbar = 100 Pa); N2-equivalent gauge pressure (20 C)",
        "nNG_N2eq_m3": "p / (k T_NG), T_NG = %g K, k = %.6e J/K (ROSINA Data Users Guide RO-ROS-UG-1000, sect. 5.2)" % (T_NG, KB),
        "nNG_H2O_m3": "beta_H2O * nNG_N2eq, beta_H2O = %g (pure H2O assumed)" % BETA_H2O,
        "PhiRG_N2eq_m2s": "p_RG / sqrt(2 pi m k T_RG), m = 18.01528 amu (Users Guide ram gauge flux formula), T_RG = COPS sensor temperature",
        "PhiRG_H2O_m2s": "beta_H2O * PhiRG_N2eq",
        "altitude": "signed distance to the nearest surface point of the shape model (negative: inside), interpolated in time "
                    "from LiouvilleCOPS (SPICE) runs: SHAP5 %s, M004 %s" % (S5["shape"], M4["shape"]),
        "r_km": "distance of the spacecraft (gauge) from the nucleus centre (LiouvilleCOPS M004 run)",
        "missing": "null"},
    "geometry_runs": {"SHAP5": {"run": args.shap5_run, "header": S5["header"][:3]},
                      "M004": {"run": args.m004_run, "header": M4["header"][:3]}},
    "events": {"NG_emission_range_HIGH_to_LOW": "2016-09-30T10:24:40",
               "NG_saturated": ["2016-09-30T10:23:46", "2016-09-30T10:24:38"],
               "NG_no_data": ["2016-09-30T10:24:42", "2016-09-30T10:33:14"],
               "last_L4_density": "2016-09-30T10:36:41", "last_L2_record_stop": hk[-1][hc["stop"]],
               "touchdown_ESA_DSK_M004": "2016-09-30T10:39:32.8 (LiouvilleCOPS, ESA DSK)"},
    "datasets": datasets,
    "not_included": not_included,
}
with open(args.out, "w") as f:
    json.dump(doc, f, indent=1)
print("JSON: %s (%d datasets)" % (args.out, len(datasets)))

# ------------------------------------------------------------------ plots (one unit: Pa)
STYLE = {"PSA_L2_BG": ("#1a1a1a", "archive L2 (PSA, BG housekeeping)"),
         "srcCG_table": ("#8c8c8c", "srcCG reference table"),
         "PSA_L4": ("#2a6fdb", "archive L4 NG density x k 293 K (PSA)"),
         "ESA_blog_2016": ("#e8710a", "ESA blog 2016 figure (digitized)"),
         "EGU2019_s02": ("#1e9e5a", "EGU 2019 slide 2 (digitized)"),
         "EGU2019_s16": ("#b8338a", "EGU 2019 slide 16 (digitized)"),
         "fit_data": ("#d62728", "used in the fit (small red markers; %s; see COPS_fit_data.png)" % os.path.relpath(args.fit_run, ROOT))}
D = {d["id"]: d["columns"] for d in datasets}
arr = lambda v: np.array([np.nan if x is None else x for x in v], float)
tdate = lambda c: [parse_utc(u) for u in c["utc"]]


def broken(c, *keys, gap=90.0):
    """columns of a time-tagged dataset with NaN rows inserted where consecutive samples are > gap s apart,
    so that lines are not drawn across data gaps"""
    t = np.array([tsec(parse_utc(u)) for u in c["utc"]])
    cut = np.nonzero(np.diff(t) > gap)[0] + 1
    out = [np.insert(t, cut, np.nan)]
    for k in keys:
        out.append(np.insert(arr(c[k]), cut, np.nan))
    tt = [T0 + dt.timedelta(seconds=float(x)) if np.isfinite(x) else None for x in out[0]]
    return [np.array([np.datetime64(x) if x is not None else np.datetime64("NaT") for x in tt])] + out[1:]


def gauge_legend(ax):
    from matplotlib.lines import Line2D
    h1 = ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.add_artist(h1)
    ax.legend(handles=[Line2D([], [], color="0.3", ls="-", lw=1.6, label="NG (nude gauge)"),
                       Line2D([], [], color="0.3", ls="--", lw=1.6, label="RG (ram gauge)")],
              loc="lower right", fontsize=8, frameon=False)


# plot 1: vs. time (full day and the last 40 minutes)
fig, ax = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1.25, 1]})
for a, (t0, t1) in zip(ax, ((T0, T0 + dt.timedelta(hours=10, minutes=45)),
                            (T0 + dt.timedelta(hours=10), T0 + dt.timedelta(hours=10, minutes=41)))):
    c = D["PSA_L2_BG"]
    a.plot(tdate(c), arr(c["pNG_Pa"]), "-", color=STYLE["PSA_L2_BG"][0], lw=1.4, label=STYLE["PSA_L2_BG"][1])
    a.plot(tdate(c), arr(c["pRG_Pa"]), "--", color=STYLE["PSA_L2_BG"][0], lw=1.4)
    c = D["srcCG_table"]
    a.plot(tdate(c), arr(c["pNG_Pa"]), "o", ms=2.5, mfc="none", color=STYLE["srcCG_table"][0], label=STYLE["srcCG_table"][1])
    a.plot(tdate(c), arr(c["pRG_Pa"]), "s", ms=2.5, mfc="none", color=STYLE["srcCG_table"][0])
    tt, pp = broken(D["PSA_L4"], "pNG_equiv_Pa")
    a.plot(tt, pp, "-", color=STYLE["PSA_L4"][0], lw=1.0, alpha=0.9, label=STYLE["PSA_L4"][1])
    c = D["ESA_blog_2016"]
    a.plot(tdate(c), arr(c["pNG_Pa"]), ".", ms=4, color=STYLE["ESA_blog_2016"][0], label=STYLE["ESA_blog_2016"][1])
    c = D["fit_data"]
    u = np.array(c["useNG"]) > 0
    a.plot(np.array(tdate(c))[u], arr(c["pNG_Pa"])[u], "o", ms=1.0, color=STYLE["fit_data"][0], label=STYLE["fit_data"][1])
    u = np.array(c["useRG"]) > 0
    a.plot(np.array(tdate(c))[u], arr(c["pRG_Pa"])[u], "s", ms=1.0, color=STYLE["fit_data"][0])
    a.set_xlim(t0, t1)
    a.set_yscale("log")
    a.set_xlabel("UTC, 30 Sept 2016")
    a.set_ylabel("pressure [Pa]  (N2-equivalent)")
    a.grid(True, which="both", alpha=0.25)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax[0].set_title("COPS NG and RG pressure vs. time (EGU figures have no time axis)", fontsize=10)
ax[1].set_title("last 40 minutes", fontsize=10)
z = ax[1]
for (a_, b_), lab, col in (((T0 + dt.timedelta(hours=10, minutes=23, seconds=46), T0 + dt.timedelta(hours=10, minutes=24, seconds=38)), "NG saturated", "#e8710a"),
                           ((T0 + dt.timedelta(hours=10, minutes=24, seconds=42), T0 + dt.timedelta(hours=10, minutes=33, seconds=14)), "no NG data", "0.6")):
    z.axvspan(a_, b_, color=col, alpha=0.18, lw=0)
    z.text(a_ + (b_ - a_) / 2, 2.0e-8 * 100, lab, ha="center", va="bottom", fontsize=8, rotation=90 if lab == "NG saturated" else 0)
z.axvline(T0 + dt.timedelta(hours=10, minutes=36, seconds=41), color=STYLE["PSA_L4"][0], ls=":", lw=1)
z.text(T0 + dt.timedelta(hours=10, minutes=36, seconds=41), 4e-7 * 100, " last L4", fontsize=8, color=STYLE["PSA_L4"][0])
gauge_legend(ax[0])
fig.tight_layout()
f1 = args.plots + "_time.png"
fig.savefig(f1, dpi=110)

# plot 2: vs. altitude
fig2, a = plt.subplots(figsize=(10.5, 7))
for key in ("PSA_L2_BG", "srcCG_table", "PSA_L4", "ESA_blog_2016"):
    c = D[key]
    alt = arr(c["altitude_M004_m"])
    col, lab = STYLE[key]
    if key == "PSA_L2_BG":
        a.plot(alt, arr(c["pNG_Pa"]), "-", color=col, lw=1.4, label=lab + ", M004 altitude")
        a.plot(alt, arr(c["pRG_Pa"]), "--", color=col, lw=1.4)
    elif key == "srcCG_table":
        a.plot(alt, arr(c["pNG_Pa"]), "o", ms=2.5, mfc="none", color=col, label=lab + ", M004 altitude")
        a.plot(alt, arr(c["pRG_Pa"]), "s", ms=2.5, mfc="none", color=col)
    elif key == "PSA_L4":
        tt, alt, pp = broken(c, "altitude_M004_m", "pNG_equiv_Pa")
        a.plot(alt, pp, "-", color=col, lw=1.0, label=lab + ", M004 altitude")
    else:
        a.plot(alt, arr(c["pNG_Pa"]), ".", ms=4, color=col, label=lab + ", M004 altitude")
c = D["fit_data"]
for gg, mk in (("NG", "o"), ("RG", "s")):
    u = np.array(c["use" + gg]) > 0
    a.plot(arr(c["altitude_M004_m"])[u], arr(c["p%s_Pa" % gg])[u], mk, ms=1.0, color=STYLE["fit_data"][0],
           label=STYLE["fit_data"][1] + ", M004 altitude" if gg == "NG" else None)
for key in ("EGU2019_s02", "EGU2019_s16"):
    c = D[key]
    col, lab = STYLE[key]
    a.plot(arr(c["altitude_NG_m"]), arr(c["pNG_Pa"]), "-", color=col, lw=1.6, label=lab + ", altitude as plotted")
    a.plot(arr(c["altitude_RG_m"]), arr(c["pRG_Pa"]), "--", color=col, lw=1.6)
a.set_xscale("log")
a.set_yscale("log")
a.set_xlim(1.0, 2.5e4)
a.set_xlabel("altitude [m]  (time-tagged sources: LiouvilleCOPS/SPICE, ESA DSK M004)")
a.set_ylabel("pressure [Pa]  (N2-equivalent)")
a.set_title("COPS NG and RG pressure vs. altitude, all sources", fontsize=10)
a.grid(True, which="both", alpha=0.25)
gauge_legend(a)
fig2.tight_layout()
f2 = args.plots + "_altitude.png"
fig2.savefig(f2, dpi=110)

# plot 3: the fit data (used: filled; not used: open, with the reason)
c = D["fit_data"]
tt = np.array(tdate(c))
alt5 = arr(c["altitude_SHAP5_m"])
fig3, ax3 = plt.subplots(1, 3, figsize=(19, 6.2), gridspec_kw={"width_ratios": [1.2, 1, 1]})
RCOL = {"no_data": None, "inside_shape_model": "#7f7f7f", "outside_FitInterval": "#9467bd", "unknown": "#000000",
        "NG_saturated": "#1f5fbf"}
for gg, mk, col in (("NG", "o", "#1f5fbf"), ("RG", "s", "#d62728")):
    p = arr(c["p%s_Pa" % gg])
    u = np.array(c["use" + gg]) > 0
    rs = np.array(c["reason_" + gg])
    for a, x in ((ax3[0], tt), (ax3[1], tt), (ax3[2], alt5)):
        a.plot(x[u], p[u], mk, ms=3, color=col, label="%s used (%d)" % (gg, u.sum()))
        for r in sorted(set(rs[~u])):
            m = (rs == r) & np.isfinite(p)
            if m.any():
                ec = RCOL.get(r) or ("#e8710a" if r.startswith("Exclude") else "#555555")
                a.plot(x[m], p[m], mk, ms=4, mfc="none", mec=ec, mew=0.8, label="%s not used: %s (%d)" % (gg, r, m.sum()))
for a, (t0, t1) in zip(ax3[:2], ((T0, T0 + dt.timedelta(hours=10, minutes=45)),
                                 (T0 + dt.timedelta(hours=10), T0 + dt.timedelta(hours=10, minutes=41)))):
    a.set_xlim(t0, t1)
    a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    a.set_xlabel("UTC, 30 Sept 2016 (time of the srcCG table = end of the L2 interval)")
    for k, v in FIV.items():
        for b0, b1 in v:
            if k == "ExcludeRG" or (k == "ExcludeNG" and args.apply_exclude_ng):
                a.axvspan(parse_utc(b0), parse_utc(b1), color="#e8710a", alpha=0.10, lw=0)
isat = np.nonzero(ngsat)[0]                                  # saturated NG points of the fit run (srcCG times)
sat0, sat1 = parse_utc(FR["utc"][isat[0]]) - dt.timedelta(seconds=4), parse_utc(FR["utc"][isat[-1]])
ax3[1].axvspan(sat0, sat1, color="#1f5fbf", alpha=0.25, lw=0)
ax3[1].text(sat0, 3e-7, "NG saturated (%d points; %s)\n--exclude-ng-saturated: %s " % (
    len(isat), "excluded" if args.exclude_ng_saturated else "included, default", "on" if args.exclude_ng_saturated else "off"),
    ha="right", fontsize=8, color="#1f5fbf")
ax3[2].set_xscale("log")
ax3[2].set_xlim(1.0, 2.5e4)
ax3[2].set_xlabel("altitude above SHAP5 [m] (fit run geometry)")
for a in ax3:
    a.set_yscale("log")
    a.set_ylabel("pressure [Pa]  (N2-equivalent)")
    a.grid(True, which="both", alpha=0.25)
ax3[0].legend(fontsize=7.5, loc="upper left", frameon=False)
ax3[0].set_title("COPS data used in the fit: %d NG, %d RG points (%s)" % (nfit["NG"], nfit["RG"], os.path.relpath(args.fit_run, ROOT)), fontsize=10)
ax3[1].set_title("last 40 minutes; ExcludeNG of %s %s" % (fit_input and os.path.basename(fit_input),
                 "applied (shaded)" if args.apply_exclude_ng else "not applied (default)"), fontsize=10)
ax3[2].set_title("vs. altitude", fontsize=10)
fig3.tight_layout()
f3 = os.path.splitext(args.fit_data)[0] + ".png"
fig3.savefig(f3, dpi=110)

# plot 4: the profiles actually used in a fit (fit_data_used.txt of tools/fit_cops.py: the points in the misfit)
f4 = None
if os.path.exists(args.fit_used):
    hdr, U = [], []
    for line in open(args.fit_used):
        if line.startswith("#iPoint"):
            un = line[1:].split()
        elif line.startswith("#"):
            hdr.append(line[1:].strip())
        elif line.strip():
            U.append(line.split())
    uc = {n: i for i, n in enumerate(un)}
    ut = np.array([tsec(parse_utc(r[uc["UTC"]])) for r in U])
    ualt = np.array([float(r[uc["altitude_m"]]) for r in U])
    ug = geom(ut)
    fig4, ax4 = plt.subplots(1, 3, figsize=(19, 6.2), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    for gg, ls, col in (("NG", "-", "#1f5fbf"), ("RG", "--", "#d62728")):
        p = np.array([float(r[uc["p%s_used_Pa" % gg]]) for r in U])
        k = np.nonzero(p > 0)[0]
        cut = np.nonzero((np.diff(k) > 1) | (np.diff(ut[k]) > 90.0))[0] + 1
        tt = np.array([T0 + dt.timedelta(seconds=float(s)) for s in ut])
        for a, x in ((ax4[0], tt), (ax4[1], tt), (ax4[2], ug["altitude_M004_m"])):
            for j, seg in enumerate(np.split(k, cut)):
                a.plot(x[seg], p[seg], ls, marker=".", ms=3, lw=1.2, color=col,
                       label="%s used in the fit (%d points)" % (gg, len(k)) if j == 0 else None)
    for a, (t0, t1) in zip(ax4[:2], ((T0, T0 + dt.timedelta(hours=10, minutes=45)),
                                     (T0 + dt.timedelta(hours=10), T0 + dt.timedelta(hours=10, minutes=41)))):
        a.set_xlim(t0, t1)
        a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        a.set_xlabel("UTC, 30 Sept 2016")
    ax4[2].set_xscale("log")
    ax4[2].set_xlim(1.0, 2.5e4)
    ax4[2].set_xlabel("altitude [m]  (LiouvilleCOPS/SPICE, ESA DSK M004)")
    for a in ax4:
        a.set_yscale("log")
        a.set_ylabel("pressure [Pa]  (N2-equivalent)")
        a.grid(True, which="both", alpha=0.25)
    ax4[0].legend(fontsize=9, loc="upper left", frameon=False)
    ax4[0].set_title("COPS NG and RG profiles used in the fit", fontsize=10)
    ax4[1].set_title("last 40 minutes", fontsize=10)
    sat_line = [h for h in hdr if "saturated NG points" in h]
    ex_line = [h for h in hdr if "ExcludeNG intervals" in h]
    ax4[2].set_title("vs. altitude" + ("; NG in ExcludeNG: %s" % ("blocked" if " blocked (--apply" in ex_line[0] else "not blocked")
                                       if ex_line else "") +
                     ("; saturated NG: %s" % ("excluded" if "excluded" in sat_line[0] else "not removed") if sat_line else ""),
                     fontsize=10)
    fig4.suptitle("%s" % os.path.abspath(args.fit_used), fontsize=8)
    fig4.tight_layout(rect=(0, 0, 1, 0.95))
    f4 = os.path.join(os.path.dirname(args.plots), "COPS_fit_profiles.png")
    fig4.savefig(f4, dpi=110)
else:
    print("no %s: the plot of the profiles used in the fit is skipped (run tools/fit_cops.py, e.g. with --data-only)" % args.fit_used)
print("plots: %s" % ", ".join(f for f in (f1, f2, f3, f4) if f))
if args.show:
    plt.show()
