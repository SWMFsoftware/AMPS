#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_all_fits_summary_docx.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_all_fits_summary_docx.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e610f02aabec5e384d0132424033ca0f4b44e69b47fdceaff0177673f4fdd213
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: ONE Word document summarizing all fits of the new-law round (H2O q_min 0.02 q_max): runs table, A18 summary, overview, per-run and full-set clusters, inspection, old-law appendix (2026-11-20)
# RUN: separate venv cops_output/tools/docx_venv; check / build / validate in Slurm jobs (HOW-TO-RERUN.md of figures_cluster/12M/q002_all_fits_summary/)
# CHANGED: 2026-11-21 user decisions 1b / 2b (FINDINGS A21)
# ==== end of header ====
"""ONE Word document summarizing all fits of the new-law round (H2O q_min = 0.02 q_max), 2026-11-20.

User request 2026-11-20 (relayed by the lead 14:45 CST): "prepare the full set of plots for the last round of fits with the
fixed surface source distribution; also prepare a single word doc summarizing all the fit". Supersedes, for this round,
12M_stage2_q002_comparison.docx (2026-11-19). Derived from q002_comparison_docx.py and fullsets_comparison_docx.py (same
rules: every printed number read from a file and listed in numbers_used.json; literal numbers of hand-written sentences
checked against the source files; observations from the images labelled as such and kept in a data file).
Runs (fixed order): D2 (D2_q002_warm), softplus, A, B, D, C f1/1e-4, C f10/1e-4, C f1/0.01, C f10/0.01, C f10/0.01 off.
Inputs:
  figures_cluster/12M/model_vs_observations_nonneg_options_stage2_q002/   summary.csv/.txt, multi-run and per-run figures,
                                                                          tables/, *_residuals_check.json
  figures_cluster/12M/full_set_q002_<run>/                                the full figure sets of the ten runs (discovered,
                                                                          or given by --sets-spec)
  figures_cluster/12M/source_rate_q002/<run>/inspection/                  new- vs old-law source-rate inspection
  --observations <json>                                                   observations from the images, per cluster id

    q002_all_fits_summary_docx.py check    [--sets-spec <json>]          inputs, discovered sets and plot types, gaps
    q002_all_fits_summary_docx.py build    --out <new folder> --observations <json> [--sets-spec <json>] [--jobs 8]
                                           [--min-px 800] [--quality 75]
    q002_all_fits_summary_docx.py validate --out <folder>
Interpreter: the separate venv $COPS_OUTPUT_ROOT/tools/docx_venv (python-docx, Pillow), not the project .venv.
Refuses to overwrite an existing document and to write it if a literal number of a hand-written sentence is not found.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import json
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import PIL
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = "/data001/highresolu_vt/cops_output/figures_cluster/12M"
SR = ROOT + "/source_rate_q002"
MV = ROOT + "/model_vs_observations_nonneg_options_stage2_q002"
FINDINGS = "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/model_ext/FINDINGS.md"
OUTVIS = "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/model_ext/OUTPUT_VISUALIZATION.md"
PREV_DOCS = [ROOT + "/stage2_q002_comparison/12M_stage2_q002_comparison.docx",
             ROOT + "/full_sets_stage2_comparison/12M_stage2_full_sets_comparison.docx"]
NUMTOK = re.compile(r"(?<![A-Za-z_\d.])\d+(?:[.,]\d+)*(?:e[+-]?\d+)?")
DOC_NAME = "12M_q002_all_fits_summary.docx"
DOC_DATE = "2026-11-22"
NEW_LAW = "new law q_min = 0.02 q_max"
OLD_LAW = "old law q_min = 0.2 q_max"
NO_LAW = "fit-independent, independent of the H2O law"

RUNS = [
    dict(key="D2", sr="D2_q002_warm", mv="D2", csv="D2_q002"),
    dict(key="softplus", sr="softplus_q002_reg", mv="softplus", csv="softplus_q002"),
    dict(key="A", sr="A_q002_reg", mv="A", csv="A_q002"),
    dict(key="B", sr="B_q002_reg", mv="B", csv="B_q002"),
    dict(key="D", sr="D_q002_reg", mv="D", csv="D_q002"),
    dict(key="C f1/1e-4", sr="C_f1_pm1em4_q002_reg", mv="C_f1_pm1em4", csv="C_f1_pm1em4_q002"),
    dict(key="C f10/1e-4", sr="C_f10_pm1em4_q002_reg", mv="C_f10_pm1em4", csv="C_f10_pm1em4_q002"),
    dict(key="C f1/0.01", sr="C_f1_pm0p01_q002_reg", mv="C_f1_pm0p01", csv="C_f1_pm0p01_q002"),
    dict(key="C f10/0.01", sr="C_f10_pm0p01_q002_reg", mv="C_f10_pm0p01", csv="C_f10_pm0p01_q002"),
    dict(key="C f10/0.01 off", sr="C_f10_pm0p01_tploff_q002_reg", mv="C_f10_pm0p01_tploff",
         csv="C_f10_pm0p01_tploff_q002"),
]
RK = [r["key"] for r in RUNS]
RUN = {r["key"]: r for r in RUNS}
BINS = ["[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]"]
BIN_LABEL = ["0-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-180"]
MULT = ["40_50", "50_60", "60_70", "70_80", "80_90", "90_180"]
FS_SECTIONS = ["model_data", "activity", "production_rate", "source_areas", "contribution_density", "cliffs", "regions",
               "visibility"]
SECTIONS = [("overview", "Overview: the multi-run figures (new law)"),
            ("mv", "Model vs COPS per run: new law vs old law"),
            ("model_data", "Full sets: model_data/ (new law)"),
            ("activity", "Full sets: activity/ (new law)"),
            ("production_rate", "Full sets: production_rate/ (new law)"),
            ("source_areas", "Full sets: source_areas/ (new law)"),
            ("contribution_density", "Full sets: contribution_density/ (new law)"),
            ("cliffs", "Full sets: cliffs/ (new law)"),
            ("regions", "Full sets: regions/ (new law)"),
            ("visibility", "Full sets: visibility/ (fit-independent)"),
            ("inspection", "Source-rate inspection: new-law vs old-law production rate per run")]
WHY_MISSING = {
    "inspection": ("no old-law maps exist for this chain, so no inspection figure",
                   SR + "/C_f10_pm1em4_q002_reg/README.txt"),
    "residuals": ("the comparison set has new-vs-old residual figures only for softplus, A, B, D, D2 and the best C "
                  "(one figure per C setting would repeat the model/data figures)", MV + "/README.txt"),
    "convergence": ("one fit, no cycles, so no convergence figure (a line in 12M_q002_convergence_all)",
                    MV + "/README.txt"),
    "fs_set": ("no full figure set of this run was found at the build", None),
    "fs_plot": ("this plot is not in the run's full figure set", None),
}

SETS = {}      # run -> full-set folder (or None)
STEMS = {}     # run -> model_data file stem
OBS = {}       # cluster id -> observation text
MISSING_INPUTS = []


# switched-off surface, user decision 1b of 2026-11-21 (q002_active_area.py; FINDINGS A21)
AA = "/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121/active_area.json"


MS_REP = "/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_q002_f10_pm0.01_tploff_20261122T044019Z/report_137981/slopetest.json"     # multi-start and control fits (FINDINGS A22)
MS_STA = "/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/multistart_q002_f10_pm0.01_tploff_20261122T065041Z/stationarity_137980/stationarity.json"


def aa_key(run):
    key = os.path.basename(SETS[run].rstrip("/")).replace("full_set_q002_", "")
    return "D2_q002_warm" if key == "D2" else key + "_q002_reg"


def offv(run, field="off_T"):
    return jload(AA)["runs"][aa_key(run)][field]


def off(run, field="off_T"):
    k = aa_key(run)
    return N.add(run, "share of the seen area switched off (decision 1b): %s" % field, 100.0 * offv(run, field), AA,
                 "active_area.json:runs/%s/%s" % (k, field), pct)


# ------------------------------------------------------------------------------------------------ numbers registry
class Numbers:
    """Every number printed in the document goes through add(); numbers_used.json lists them with their source."""

    def __init__(self):
        self.items, self.index, self.ctx = [], {}, "front matter"
        self.ctx_tokens, self.unverified = {}, []

    def add(self, run, quantity, value, source, locator, fmt="%s", note=None):
        text = fmt(value) if callable(fmt) else (fmt % value)
        self.ctx_tokens.setdefault(self.ctx, set()).update(NUMTOK.findall(text))
        k = (run, quantity, source, locator, text)
        if k in self.index:
            it = self.items[self.index[k]]
            if self.ctx not in it["used_in"]:
                it["used_in"].append(self.ctx)
            return text
        v = value
        if isinstance(v, float) and not math.isfinite(v):
            v = str(v)
        self.items.append(dict(id=len(self.items) + 1, run=run, quantity=quantity, value=v, text=text,
                               source=source, locator=locator, note=note, used_in=[self.ctx]))
        self.index[k] = len(self.items) - 1
        return text


N = Numbers()
_cache = {}


def tload(path):
    if path not in _cache:
        with open(path) as f:
            _cache[path] = f.read()
    return _cache[path]


def jload(path):
    if ("json", path) not in _cache:
        with open(path) as f:
            _cache[("json", path)] = json.load(f)
    return _cache[("json", path)]


def T(text, sources, fallback=None, force=False):
    """hand-written text: every numeric token must occur in one of the source files (recorded in numbers_used.json)"""
    skip = set() if force else N.ctx_tokens.get(N.ctx, set())
    for tok in NUMTOK.findall(text):
        if tok in skip:
            continue
        found = None
        for src in sources:
            if src and os.path.exists(src):
                body = tload(src)
                if tok in body or tok.replace(",", "") in body:
                    found = src
                    break
        if found:
            N.add("text", "number in hand-written text", tok, found, "substring of the source file", "%s")
        elif fallback:
            N.add("text", "number in hand-written text", tok, fallback[1], fallback[0], "%s")
        else:
            N.unverified.append(dict(ctx=N.ctx, token=tok, text=text[:160]))
    return text


def safe(fn, *a, what="", **kw):
    """an extractor that cannot find its input prints 'n/a' and is listed in missing_inputs.json"""
    try:
        v = fn(*a, **kw)
        return "n/a" if v is None else v
    except (KeyError, IndexError, OSError, TypeError, ValueError, AttributeError) as e:
        MISSING_INPUTS.append(dict(ctx=N.ctx, what=what or getattr(fn, "__name__", "?"), args=[str(x) for x in a],
                                   error="%s: %s" % (type(e).__name__, e)))
        return "n/a"


def nested(d, keys):
    for k in keys:
        d = d[k]
    return d


def pct(v):
    return "%.1f %%" % v


def pct0(v):
    return "%.0f %%" % v


def e2(v):
    return "%.2e" % v


def e3(v):
    return "%.3e" % v


def e4(v):
    return "%.4e" % v


def intc(v):
    return "{:,}".format(int(v))


def sgn_pct(v):
    return "%+.1f %%" % v


def order_text(vals, reverse=False):
    vals = {k: v for k, v in vals.items() if v is not None and v[0] is not None and v[1] != "n/a"}
    items = sorted(vals.items(), key=lambda kv: kv[1][0], reverse=reverse)
    return ", ".join("%s %s" % (k, t) for k, (v, t) in items)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ summary.csv (MV)
def srow(key):
    p = MV + "/summary.csv"
    if p not in _cache:
        with open(p) as f:
            _cache[p] = {r["key"]: r for r in csv.DictReader(f)}
    return _cache[p][key], p


def row_key(run, old=False):
    r, _ = srow(RUN[run]["csv"])
    return r["old_key"] if old else RUN[run]["csv"]


def sval(run, col, old=False):
    r, _ = srow(row_key(run, old))
    return r[col]


def sx(run, col, quantity, fmt="%s", conv=float, old=False, scale=1.0):
    key = row_key(run, old)
    r, p = srow(key)
    raw = r[col]
    if raw == "":
        return None
    v = conv(raw) if conv else raw
    if scale != 1.0:
        v = v * scale
    return N.add(run, quantity + (" (old law)" if old else ""), v, p,
                 "summary.csv row key=%s column %s%s" % (key, col, " x %g" % scale if scale != 1.0 else ""), fmt)


def resid_json_for(run):
    p = MV + "/12M_q002_vs_old_%s_residuals_check.json" % RUN[run]["mv"]
    return p if os.path.exists(p) else None


def rms(run, gauge, late, old=False):
    col = "rms_%s_%s" % (gauge, "late" if late else "all")
    s = sx(run, col, "rms ln %s (%s)" % (gauge, "late window" if late else "all fitted rows"), "%.4f", old=old)
    if s is not None:
        return s
    for p in [resid_json_for(run), MV + "/12M_q002_best_residuals_check.json"]:
        if not p:
            continue
        tab = row_key(run, old) + ".dat"
        for k, v in jload(p)["fits"].items():
            if os.path.basename(v["table"]) == tab:
                return N.add(run, "rms ln %s (%s)%s" % (gauge, "late window" if late else "all fitted rows",
                                                        " (old law)" if old else ""),
                             v[gauge]["rms_late" if late else "rms"], p,
                             "fits['%s'].%s.%s" % (k, gauge, "rms_late" if late else "rms"), "%.4f")
    return "n/a"


def table(run, old=False):
    p = MV + "/tables/" + row_key(run, old) + ".dat"
    key = ("tab", p)
    if key not in _cache:
        cols, ng, rg = None, [], []
        for line in open(p):
            if line.startswith("#iPoint"):
                cols = line[1:].split()
                continue
            if line.startswith("#") or not line.strip():
                continue
            r = line.split()
            ix_ = {n: i for i, n in enumerate(cols)}
            alt = float(r[ix_["altitude"]])
            if r[ix_["useNG"]] == "1":
                ng.append((alt, math.log(float(r[ix_["pNG"]]) / float(r[ix_["refNG"]]))))
            if r[ix_["useRG"]] == "1":
                rg.append((alt, math.log(float(r[ix_["pRG"]]) / float(r[ix_["refRG"]]))))
        _cache[key] = (ng, rg, p)
    return _cache[key]


def t_min_ng(run, old=False):
    ng, rg, p = table(run, old)
    v = min(x for a, x in ng if 10.0 <= a <= 100.0)
    return v, N.add(run, "lowest ln(model/obs) NG, gauge altitude 10-100 m" + (" (old law)" if old else ""), v, p,
                    "min over rows useNG=1, 10<=altitude<=100 m of ln(pNG/refNG)", "%+.2f")


def t_max_rg(run, old=False):
    ng, rg, p = table(run, old)
    v = max(x for a, x in rg if a < 20.0)
    return v, N.add(run, "highest ln(model/obs) RG, gauge altitude < 20 m" + (" (old law)" if old else ""), v, p,
                    "max over rows useRG=1, altitude<20 m of ln(pRG/refRG)", "%+.2f")


def t_last(run, gauge, old=False):
    ng, rg, p = table(run, old)
    rows = ng if gauge == "NG" else rg
    a, v = min(rows, key=lambda t: t[0])
    return v, N.add(run, "ln(model/obs) %s at the last observation" % gauge + (" (old law)" if old else ""), v, p,
                    "row useNG/RG=1 with the lowest altitude (%.2f m): ln(p%s/ref%s)" % (a, gauge, gauge), "%+.2f")


# ------------------------------------------------------------------------------------------------ source-rate inspection
def srp(run, rel=""):
    return os.path.join(SR, RUN[run]["sr"], rel)


def insp_path(run, t):
    p = srp(run, "inspection/checks/inspect_%s.json" % t)
    return p if os.path.exists(p) else None


def ix(run, t, keys, quantity, fmt="%s", scale=1.0):
    p = insp_path(run, t)
    if p is None:
        return None
    v = nested(jload(p), keys)
    if scale != 1.0:
        v = v * scale
    return N.add(run, quantity, v, p, "inspect_%s.json:%s%s" % (t, "/".join(keys), " x %g" % scale if scale != 1.0 else ""),
                 fmt)


def ixv(run, t, keys):
    p = insp_path(run, t)
    return None if p is None else nested(jload(p), keys)


def capx(path, pattern, run, quantity, group=1, fmt="%s", conv=str):
    m = re.search(pattern, tload(path), re.S)
    if not m:
        MISSING_INPUTS.append(dict(ctx=N.ctx, what=quantity, args=[path], error="regex no match: " + pattern))
        return "n/a"
    return N.add(run, quantity, conv(m.group(group)), path, "regex group %d: %s" % (group, pattern), fmt)


def manx(path, keys, run, quantity, fmt="%s", scale=1.0):
    v = nested(jload(path), keys)
    if scale != 1.0:
        v = v * scale
    return N.add(run, quantity, v, path, os.path.basename(path) + ":" + "/".join(str(k) for k in keys), fmt)


# ------------------------------------------------------------------------------------------------ full sets
def discover_sets(spec=None):
    """run -> full_set_q002_* folder: --sets-spec json {run key: folder} or the folder whose name holds the run's
    source-rate name, comparison-set key or summary key (exactly one match, else None)"""
    out = {}
    if spec:
        d = json.load(open(spec))
        for r in RK:
            out[r] = d.get(r)
        return out
    cands = sorted(p for p in glob.glob(ROOT + "/full_set_q002_*") if os.path.isdir(p) and not re.search(
        r"\.(superseded|incomplete|tar|provenance)", p))
    for r in RK:
        names = {RUN[r]["sr"], RUN[r]["mv"], RUN[r]["csv"], RUN[r]["sr"].replace("_q002", "").replace("_reg", ""),
                 RUN[r]["mv"] + "_reg", RUN[r]["mv"] + "_q002", RUN[r]["mv"] + "_q002_reg"}
        hits = [p for p in cands if os.path.basename(p)[len("full_set_q002_"):] in names]
        out[r] = hits[0] if len(hits) == 1 else None
    return out


def find_stem(folder):
    hits = glob.glob(os.path.join(folder, "model_data", "*_vs_COPS_altitude.png"))
    if len(hits) != 1:
        return None
    return os.path.basename(hits[0])[:-len("_vs_COPS_altitude.png")]


def fs_list(run):
    """relative PNG paths of a full set with the run's stem replaced by <STEM>"""
    folder = SETS.get(run)
    if not folder:
        return []
    stem = STEMS.get(run)
    out = []
    for sec in FS_SECTIONS:
        base = os.path.join(folder, sec)
        for dp, dn, fn in os.walk(base):
            rel = os.path.relpath(dp, folder)
            if re.search(r"(incomplete|superseded|provenance|/code|/tables)", "/" + rel):
                continue
            for f in fn:
                if f.endswith(".png"):
                    r_ = os.path.join(rel, f)
                    out.append(r_.replace(stem, "<STEM>") if stem else r_)
    return out


def fsp(run, rel):
    folder = SETS.get(run)
    if not folder:
        return None
    p = os.path.join(folder, rel.replace("<STEM>", STEMS.get(run) or "<STEM>"))
    return p if os.path.exists(p) else None


def fsc(run, rel):
    """the full set's file for a relative path (no existence filter)"""
    return os.path.join(SETS[run], rel.replace("<STEM>", STEMS.get(run) or "<STEM>"))


# plot-type table for the full sets: (regex on the rel path, title, family, params); first match wins
def fs_types():
    T_ = []

    def add(pat, title, fam, **kw):
        T_.append((re.compile("^" + pat + "$"), title, fam, kw))

    add(r"model_data/publication_NG_RG_altitude/<STEM>_NG_RG_altitude\.png",
        "Publication figure: COPS and model NG / RG against altitude", "md_pub")
    add(r"model_data/<STEM>_vs_COPS_altitude\.png", "Model vs COPS against altitude: (a) NG, (b) RG, (c) NG and RG", "md_alt")
    add(r"model_data/<STEM>_vs_COPS_altitude_NG\.png", "Model vs COPS against altitude: NG", "md_ng")
    add(r"model_data/<STEM>_vs_COPS_altitude_RG\.png", "Model vs COPS against altitude: RG", "md_rg")
    add(r"model_data/<STEM>_vs_COPS_altitude_NG_RG\.png", "Model vs COPS against altitude: NG and RG together, crossings",
        "md_ngrg")
    add(r"model_data/<STEM>_vs_COPS_time\.png", "Model vs COPS against UTC", "md_time")
    add(r"model_data/<STEM>_pressures_altitude\.png", "Pressures and NG / RG ratio against altitude", "md_press")
    add(r"model_data/<STEM>_pressures_time\.png", "Pressures and NG / RG ratio against UTC", "md_press")
    for sec, q, what0, pq in (("activity", "activity", "H2O activity factor A", "activity"),
                              ("production_rate", "h2o", "H2O production rate F_H2O", "h2o_production_rate")):
        tpre = r"12M_templates_activity" if q == "activity" else r"12M_templates_h2o_production_rate_(?P<t>\d{6})UTC"
        gpre = r"12M_groundtrack_activity" if q == "activity" else r"12M_groundtrack_h2o_production_rate_(?P<t>\d{6})UTC"
        add(sec + "/" + tpre + r"\.png", what0 + "{T}: maps with effective slope and cliff templates (six panels)",
            "map_tpl", q=q, panels=[("a", "latlon_global"), ("b", "latlon_zoom")])
        add(sec + "/" + gpre + r"_latlon\.png", what0 + "{T} with the ground track, longitude / latitude", "map_latlon",
            q=q, panels=[("a", "latlon_global"), ("b", "latlon_zoom")])
        add(sec + "/" + gpre + r"_3d\.png", what0 + "{T} with the ground track, 3-D views", "map_3d", q=q,
            panels=[("a", "3d_a"), ("b", "3d_b"), ("c", "3d_c")])
        add(sec + "/panels/" + tpre + "_a_" + pq + r"_global\.png", what0 + "{T}: panel (a) global map alone", "map_panel",
            q=q, panels=[("a", "latlon_global")])
        add(sec + "/panels/" + tpre + "_b_" + pq + r"_impact_area\.png", what0 + "{T}: panel (b) impact area alone",
            "map_panel", q=q, panels=[("b", "latlon_zoom")])
        for p, nm in (("c_slope_global", "panel (c) effective gravitational slope, global"),
                      ("d_slope_impact_area", "panel (d) effective gravitational slope, impact area"),
                      ("e_cliff_global", "panel (e) cliff templates, global"),
                      ("f_cliff_impact_area", "panel (f) cliff templates, impact area")):
            add(sec + "/panels/" + tpre + "_" + p + r"\.png", what0 + "{T}: " + nm, "geom_panel", q=q)
    add(r"source_areas/12M_source_area_union_(?P<g>NG|RG|NG_RG)\.png", "Union of the P% source areas, P = 50-90 %, {G}",
        "sa_union")
    add(r"source_areas/12M_source_area_occurrence_P(?P<P>\d\d)_(?P<g>NG|RG|NG_RG)\.png", "Source-area occurrence N_{P}, {G}",
        "sa_occ")
    add(r"contribution_density/12M_contrib_density_(?P<k>103310|103731|103909|103926|pre_gap|gap|late)_(?P<g>NG|RG|NG_RG)\.png",
        "Contribution density, {K}, {G}", "cd")
    add(r"cliffs/12M_cliff_slope_bins\.png", "Use of the effective-slope bins: area, multipliers, production", "cliff_bins")
    add(r"cliffs/12M_cliff_pressure_fraction_(?P<g>NG|RG|NG_RG)\.png", "Share of the modelled surface pressure by slope bin, {G}",
        "cliff_pf")
    add(r"regions/12M_region_share_lobe_(?P<g>NG|RG|NG_RG)\.png", "Share of the modelled surface pressure by lobe, {G}",
        "reg_lobe")
    add(r"regions/12M_region_share_distance_(?P<g>NG|RG|NG_RG)\.png",
        "Share of the modelled surface pressure by distance from the touchdown point, {G}", "reg_distance")
    add(r"visibility/12M_templates_visibility\.png", "Late-visibility templates, effective slope and cliff templates", "vis")
    add(r"visibility/12M_groundtrack_visibility_latlon\.png", "Ground track and impact point, longitude / latitude", "vis")
    add(r"visibility/12M_groundtrack_visibility_3d\.png", "Ground track and impact point, 3-D views", "vis")
    return T_


CD_NAME = {"103310": "10:33:10 UTC (start of the late window)", "103731": "10:37:31.282 UTC", "103909": "10:39:09.282 UTC",
           "103926": "10:39:26 UTC (last observation)", "pre_gap": "mean before the NG gap",
           "gap": "mean over the gap in the NG data", "late": "mean over the late window"}
FS_ORDER = ["fs_md", "map_tpl", "map_latlon", "map_3d", "map_panel", "geom_panel", "sa_union", "sa_occ", "cd", "cliff_bins",
            "cliff_pf", "reg_lobe", "reg_distance", "vis", "generic"]


def fs_classify(rel):
    for rx, title, fam, kw in fs_types():
        m = rx.match(rel)
        if m:
            gd = m.groupdict()
            t = gd.get("t")
            ttl = title.replace("{T}", " at %s:%s:%s UTC" % (t[:2], t[2:4], t[4:]) if t else "")
            ttl = ttl.replace("{G}", (gd.get("g") or "").replace("_", " and ")).replace("{P}", gd.get("P") or "")
            ttl = ttl.replace("{K}", CD_NAME.get(gd.get("k"), gd.get("k") or ""))
            d = dict(title=ttl, fam=fam, t=t, gauges=(gd.get("g") or "").split("_") if gd.get("g") else None,
                     P=int(gd["P"]) if gd.get("P") else None, key=gd.get("k"))
            d.update(kw)
            return d
    base = os.path.basename(rel)[:-4].replace("<STEM>", "run").replace("_", " ")
    return dict(title=base, fam="generic", t=None, gauges=None, P=None, key=None)


def fs_sort_key(rel):
    c = fs_classify(rel)
    sec = rel.split("/")[0]
    fam = "fs_md" if c["fam"].startswith("md") else c["fam"]
    sub = {"md_pub": 0, "md_alt": 1, "md_ng": 2, "md_rg": 3, "md_ngrg": 4, "md_time": 5, "md_press": 6}.get(c["fam"], 0)
    return (FS_SECTIONS.index(sec) if sec in FS_SECTIONS else 99, c.get("t") or "", FS_ORDER.index(fam) if fam in FS_ORDER
            else 99, sub, c.get("key") or "", c.get("P") or 0, ["NG", "RG", "NG_RG"].index("_".join(c["gauges"]))
            if c.get("gauges") and "_".join(c["gauges"]) in ("NG", "RG", "NG_RG") else 0, rel)


# ------------------------------------------------------------------------------------------------ catalogue
def catalogue():
    G = []

    def add(sec, cid, title, fam, paths, law, **kw):
        G.append(dict(section=sec, cid=cid, title=title, fam=fam, paths=paths, law=law, **kw))

    for setn, who in (("best", "best of each option (softplus, D2, A, B, best C, D)"),
                      ("C", "the five option-C chains and the softplus reference")):
        for suf, what in (("time", "against UTC"), ("altitude", "against altitude: (a) NG, (b) RG, (c) NG and RG"),
                          ("altitude_NG", "NG against altitude"), ("altitude_RG", "RG against altitude"),
                          ("altitude_NG_RG", "NG and RG together against altitude"),
                          ("residuals", "log residuals ln(model/obs) against altitude")):
            add("overview", "ov:" + setn, "New law, %s: %s" % (who, what), "ov",
                {"all": MV + "/12M_q002_%s_%s.png" % (setn, suf)}, NEW_LAW, set=setn)
    add("overview", "ov:conv", "New law: data cost per cycle of the nine stage-2 chains and the refitted D2", "ov_conv",
        {"all": MV + "/12M_q002_convergence_all.png"}, NEW_LAW)
    for suf, what, fam in (("time", "Model vs COPS against UTC", "vo"),
                           ("altitude", "Model vs COPS against altitude: (a) NG, (b) RG, (c) NG and RG", "vo"),
                           ("altitude_NG", "Model vs COPS against altitude: NG", "vo_ng"),
                           ("altitude_RG", "Model vs COPS against altitude: RG", "vo_rg"),
                           ("altitude_NG_RG", "Model vs COPS against altitude: NG and RG together", "vo_ngrg"),
                           ("residuals", "Log residuals ln(model/obs) against altitude", "vr")):
        paths = {}
        for r in RK:
            p = MV + "/12M_q002_vs_old_%s_%s.png" % (RUN[r]["mv"], suf)
            paths[r] = p if os.path.exists(p) else None
        add("mv", "mv:" + ("residuals" if fam == "vr" else "vo"), "%s, new vs old law, per run" % what, fam, paths,
            "%s (solid) vs %s (dashed)" % (NEW_LAW, OLD_LAW), suf=suf, why="residuals" if fam == "vr" else None)
    paths = {}
    for r in RK:
        p = MV + "/12M_q002_convergence_%s.png" % RUN[r]["mv"]
        paths[r] = p if os.path.exists(p) else None
    add("mv", "mv:convergence", "Data cost per cycle (fit and fine tuning), new vs old law, per chain", "cv", paths,
        "%s (solid black) vs %s (dashed grey)" % (NEW_LAW, OLD_LAW), why="convergence")
    rels = sorted(set(x for r in RK for x in fs_list(r)), key=fs_sort_key)
    for rel in rels:
        c = fs_classify(rel)
        paths = {r: fsp(r, rel) for r in RK}
        whys = {r: ("fs_set" if not SETS.get(r) else "fs_plot") for r in RK if not paths[r]}
        law = NO_LAW if c["fam"] in ("vis", "geom_panel") else NEW_LAW
        add(rel.split("/")[0], "fs:" + rel, c["title"], c["fam"], paths, law, rel=rel, whys=whys,
            **{k: v for k, v in c.items() if k not in ("title", "fam")})
    for t in ("103310", "103926"):
        rel = "inspection/12M_inspect_h2o_production_rate_new_vs_old_law_%sUTC.png" % t
        paths = {r: (srp(r, rel) if os.path.exists(srp(r, rel)) else None) for r in RK}
        add("inspection", "sr:inspection", "Inspection at %s:%s:%s UTC: new-law vs old-law F_H2O on common classes, and "
            "the new q/q_max" % (t[:2], t[2:4], t[4:]), "insp", paths,
            "%s (a, b) vs %s (c, d), common colour classes; (e, f) the new law's q/q_max" % (NEW_LAW, OLD_LAW),
            t=t, rel=rel, why="inspection")
    return G


def caption_of(png):
    return png[:-4] + "_caption.txt" if png else None


def why_of(g, r):
    w = (g.get("whys") or {}).get(r) or g.get("why")
    return WHY_MISSING[w]


# ------------------------------------------------------------------------------------------------ images and layout
PAGE = {"P": (8.5, 11.0), "L": (11.0, 8.5)}
MARGIN = 0.6
OVERHEAD = 1.55
LABEL = 0.24
PAD = 0.12


def usable(o):
    w, h = PAGE[o]
    return w - 2 * MARGIN, h - 2 * MARGIN


def layout(aspects, n):
    rmax = max(aspects) if aspects else 1.0
    best = None
    opts = [("P", 1, 1), ("L", 1, 1)] if n == 1 else [(o, c, math.ceil(n / c)) for o in ("P", "L") for c in (2, 3, 4, 5)
                                                       if c <= n]
    for o, cols, rows in opts:
        W, H = usable(o)
        hrow = (H - OVERHEAD - (0.3 if n == 1 else 0.0)) / rows - LABEL - 0.05
        w = min(W / cols - PAD, hrow / rmax)
        if n == 1:
            w = min(w, 7.0)
        score = w * w * rmax
        if best is None or score > best[4] + 1e-9:
            best = (o, cols, rows, w, score)
    return dict(orient=best[0], cols=best[1], rows=best[2], width_in=round(best[3], 3))


def prep_image(job):
    src, dst, kind, target, quality = job
    if os.path.exists(dst):
        im = Image.open(dst)
        return dst, im.size, os.path.getsize(dst)
    im = Image.open(src)
    im.load()
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[3])
        im = bg
    else:
        im = im.convert("RGB")
    w = min(im.size[0], target)
    if w < im.size[0]:
        im = im.resize((w, round(im.size[1] * w / im.size[0])), Image.LANCZOS)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + ".part"
    if kind == "png":
        im.quantize(colors=256, method=Image.Quantize.MEDIANCUT).save(tmp, "PNG", optimize=True)
    else:
        im.save(tmp, "JPEG", quality=quality, optimize=True)
    os.replace(tmp, dst)
    return dst, im.size, os.path.getsize(dst)


# ------------------------------------------------------------------------------------------------ text sources
def readmes():
    out = [MV + "/README.txt", MV + "/summary.txt"] + [srp(r, "README.txt") for r in RK]
    out += [os.path.join(SETS[r], "README.txt") for r in RK if SETS.get(r)]
    return out


def gsrc(g):
    caps = [caption_of(p) for p in g["paths"].values() if p]
    caps = [c for c in caps if c and os.path.exists(c)]
    return caps + readmes() + [FINDINGS, OUTVIS, AA]


def fig_fallback(g):
    first = [p for p in g["paths"].values() if p][0]
    return ("observation from the figure: approximate axis reading or a bin / category label shown in the figure",
            first + " and the other runs' PNG of this figure")


def front_sources():
    return readmes() + [MV + "/CAPTIONS.md", FINDINGS, OUTVIS, AA]


# ------------------------------------------------------------------------------------------------ document
def build_doc(out, groups, img, args):
    from docx import Document
    from docx.enum.section import WD_ORIENT, WD_SECTION
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt, RGBColor

    doc = Document()
    st = doc.styles
    st["Normal"].font.name = "Arial"
    st["Normal"].font.size = Pt(10)
    st["Normal"].element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    st["Normal"].paragraph_format.space_after = Pt(4)
    for name, size, before, after in (("Title", 22, 0, 6), ("Heading 1", 16, 12, 6), ("Heading 2", 13, 6, 4),
                                      ("Heading 3", 11, 6, 3)):
        s = st[name]
        s.font.name = "Arial"
        s.font.size = Pt(size)
        s.font.bold = True
        s.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)
        s.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
        s.paragraph_format.space_before = Pt(before)
        s.paragraph_format.space_after = Pt(after)
        s.paragraph_format.keep_with_next = True
    cap_style = st["Caption"]
    cap_style.font.name = "Arial"
    cap_style.font.size = Pt(9)
    cap_style.font.italic = False
    cap_style.font.color.rgb = RGBColor(0, 0, 0)
    zoom = doc.settings.element.find(qn("w:zoom"))
    if zoom is not None and zoom.get(qn("w:percent")) is None:
        zoom.set(qn("w:percent"), "100")

    def set_orient(sec, o):
        w, h = PAGE[o]
        sec.orientation = WD_ORIENT.PORTRAIT if o == "P" else WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = Inches(w), Inches(h)
        for a in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
            setattr(sec, a, Inches(MARGIN))
        sec.header_distance = Inches(0.3)
        sec.footer_distance = Inches(0.3)

    CUR = {"orient": "L"}

    def new_section(o):
        sec = doc.add_section(WD_SECTION.NEW_PAGE)
        set_orient(sec, o)
        CUR["orient"] = o
        return sec

    SRC = {"src": [], "fb": None}

    def para(text="", bold_lead=None, italic=False, size=None, style=None, after=None, force=False, fb=None, src=None):
        text = T(text, SRC["src"] if src is None else src, fb or SRC["fb"], force=force)
        p = doc.add_paragraph(style=style)
        if bold_lead:
            r = p.add_run(bold_lead + " ")
            r.bold = True
            if size:
                r.font.size = Pt(size)
        r = p.add_run(text)
        r.italic = italic
        if size:
            r.font.size = Pt(size)
        if after is not None:
            p.paragraph_format.space_after = Pt(after)
        return p

    def bullet(text, lead=None, size=None):
        text = T(text, SRC["src"], SRC["fb"])
        p = doc.add_paragraph(style="List Bullet")
        if lead:
            r = p.add_run(lead + " ")
            r.bold = True
            if size:
                r.font.size = Pt(size)
        r = p.add_run(text)
        if size:
            r.font.size = Pt(size)
        p.paragraph_format.space_after = Pt(2)
        return p

    def shade(cell, hexcol):
        tcPr = cell._tc.get_or_add_tcPr()
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear")
        shd.set(qn("w:color"), "auto")
        shd.set(qn("w:fill"), hexcol)
        tcPr.append(shd)

    def tblpr_insert(tblPr, el):
        for succ in ("w:tblLook", "w:tblCaption", "w:tblDescription"):
            s_el = tblPr.find(qn(succ))
            if s_el is not None:
                s_el.addprevious(el)
                return
        tblPr.append(el)

    def fixed(table, widths):
        tblPr = table._tbl.tblPr
        lay = OxmlElement("w:tblLayout")
        lay.set(qn("w:type"), "fixed")
        tblpr_insert(tblPr, lay)
        for i, col in enumerate(table.columns):
            col.width = Inches(widths[i])
        for row in table.rows:
            for i, cell in enumerate(row.cells):
                cell.width = Inches(widths[i])

    def cell_margins(table, twips=40):
        tblPr = table._tbl.tblPr
        mar = OxmlElement("w:tblCellMar")
        for side in ("top", "left", "bottom", "right"):
            e = OxmlElement("w:" + side)
            e.set(qn("w:w"), str(twips if side in ("left", "right") else 10))
            e.set(qn("w:type"), "dxa")
            mar.append(e)
        tblpr_insert(tblPr, mar)

    def header_repeat(row):
        trPr = row._tr.get_or_add_trPr()
        e = OxmlElement("w:tblHeader")
        e.set(qn("w:val"), "true")
        trPr.append(e)

    def data_table(header, rows, widths=None, size=7.5, first_bold=True):
        avail = usable(CUR["orient"])[0] - 0.1
        if widths is None:
            widths = [1.2] + [1.0] * (len(header) - 1)
        if sum(widths) > avail:
            widths = [w * avail / sum(widths) for w in widths]
        t = doc.add_table(rows=1 + len(rows), cols=len(header))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, h in enumerate(header):
            c = t.rows[0].cells[i]
            c.text = ""
            r = c.paragraphs[0].add_run(T(h, SRC["src"], SRC["fb"]))
            r.bold = True
            r.font.size = Pt(size)
            shade(c, "D9E2F3")
        header_repeat(t.rows[0])
        for j, row in enumerate(rows):
            for i, v in enumerate(row):
                c = t.rows[j + 1].cells[i]
                c.text = ""
                parts = v if isinstance(v, list) else [v]
                for k, part in enumerate(parts):
                    p = c.paragraphs[0] if k == 0 else c.add_paragraph()
                    p.paragraph_format.space_after = Pt(0)
                    r = p.add_run(T("-" if part is None else str(part), SRC["src"], SRC["fb"]))
                    r.font.size = Pt(size)
                    if i == 0 and first_bold:
                        r.bold = True
        fixed(t, widths)
        cell_margins(t, 50)
        doc.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    def add_footer(sec):
        p = sec.footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run("12M, new H2O law (q_min = 0.02 q_max): all fits of the last round, %s   |   page " % DOC_DATE)
        r.font.size = Pt(8)
        fld = OxmlElement("w:fldSimple")
        fld.set(qn("w:instr"), "PAGE")
        rr = OxmlElement("w:r")
        rpr = OxmlElement("w:rPr")
        sz = OxmlElement("w:sz")
        sz.set(qn("w:val"), "16")
        rpr.append(sz)
        rr.append(rpr)
        t = OxmlElement("w:t")
        t.text = "1"
        rr.append(t)
        fld.append(rr)
        p._p.append(fld)

    sec0 = doc.sections[0]
    set_orient(sec0, "L")
    add_footer(sec0)
    FR = front_sources()
    have_sets = [r for r in RK if SETS.get(r)]

    # ---------------- title, purpose, law, how to read
    N.ctx = "1. title / purpose"
    SRC["src"], SRC["fb"] = FR, ("document date / build record (generator constants, build_manifest.json)",
                                 os.path.abspath(__file__))
    doc.add_paragraph("12M, corrected H2O law: all fits of the last round", style="Title")
    para(DOC_DATE + " (US Central). SHAP8 12M shape model, COPS NG / RG data of Rosetta's final descent (2016-09-30), "
         "collisionless (Liouville) forward model.", italic=True)
    para("One document for all fits of the last fitting round, the refitted D2 (D2_q002_warm) and the nine regularized "
         "stage-2 chains refitted with the corrected H2O night-side law: the runs table and the summary of the round, "
         "the multi-run figures, the per-run model/data, residual and convergence figures (new law vs old law), one "
         "cluster per plot type of the runs' full figure sets, and the source-rate inspection figures. In every cluster "
         "the runs appear in the order D2, softplus, A, B, D, C f1/1e-4, C f10/1e-4, C f1/0.01, C f10/0.01, C f10/0.01 "
         "off, each image labelled, followed by a summary: what the plot shows (from its caption), a table of the "
         "numbers of every run (read from the files; all listed with their source in numbers_used.json), how the runs "
         "compare, and what is visible in the images (marked as observations). No new fits or figures and no physical "
         "interpretation beyond the captions, READMEs and FINDINGS.md. For this round it supersedes "
         "12M_stage2_q002_comparison.docx (2026-11-19); 12M_stage2_full_sets_comparison.docx (2026-11-17) describes the "
         "old-law states.", bold_lead="Purpose.")
    para("Every fit, map and figure before 2026-11-18 used the H2O production q_H2O = q_min on night and cast-shadow "
         "facets and q_min + (q_max - q_min) cos SZA on day facets with q_min = 0.2 q_max (not intended); srcCG has "
         "0.02. User decision 2026-11-18: use 0.02 (FINDINGS A18). All runs here were refitted on base_q002 with the "
         "NEW law q_min = 0.02 q_max; the per-run model/data, residual, convergence and inspection figures compare them "
         "with their OLD-law (q_min = 0.2 q_max) counterparts, and Appendix A lists the old-law results for context. "
         "Every figure caption below states the law of what it shows.", bold_lead="The law correction.")
    para("Run labels: C f1/1e-4 = option C, smoothness f 1, multiplier penalty P_m 1e-4 (f10: f 10; 0.01: P_m 0.01; "
         "off: templates off). A run that lacks a plot keeps its place in the grid with the reason. Plot types whose "
         "PNGs are byte-identical in every run (SHA-256) are shown once. Map colour classes are chosen per panel and per "
         "run, so equal colours do not mean equal values across runs (the inspection figures compare new and old law of "
         "one run on common classes); source-area occurrence classes, slope-bin and region shares use the same scale in "
         "every run. Full figure sets found at the build: %d of 10 (%s). Images: downscaled copies (%s); originals "
         "(300 dpi PNG and EPS) with SHA-256 in Appendix B." % (len(have_sets), ", ".join(have_sets) or "none",
                                                               args.img_note), bold_lead="How to read it.")
    para("The version of 2026-11-20 measured switched-off surface by activity exactly zero. That showed the option-C fits "
         "as active everywhere, which is wrong: their bounded solver leaves switched-off activity at tiny positive "
         "values. The numbers now follow the user's decision 1b of 2026-11-21 (FINDINGS A21; Section 2). The map images "
         "are unchanged: light green in them means exactly zero, and switched-off terrain of the C fits and the imposed "
         "background of B appear in colour. The statement that the overhang dominance of NG at touchdown is not "
         "required by the data is also corrected, and the slope statements now name the lower bound of the multipliers.",
         bold_lead="Correction (2026-11-21).")
    para("The option-C fit with templates off (C f10/0.01 off) is one local minimum of its problem, not its best solution, "
         "although it ended on ftol. From other starts the same model has four further solutions with lower objective and "
         "lower data cost, each with nearly all the seen production on another slope class. Its row, figures and maps "
         "describe that one solution. Its change against the old law compares two different solutions and is not an "
         "effect of the law, and the shift of the multipliers with the law is not established either (Section 2, "
         "'Several solutions'; FINDINGS A22). The other runs were not tested for other solutions.",
         bold_lead="Correction (2026-11-22).")

    # ---------------- runs table
    N.ctx = "2. runs table"
    SRC["fb"] = None
    doc.add_paragraph("1. Runs", style="Heading 1")
    para("Final state = the last cycle's fine tuning (c5/ft; D2: its fit, no fine tuning). Columns from "
         "model_vs_observations_nonneg_options_stage2_q002/summary.csv (new-law row and its old-law counterpart, "
         "change = new / old - 1); rms ln(model/obs) all fitted rows (late window W10x5 in brackets; D2's late values "
         "from 12M_q002_vs_old_D2_residuals_check.json); Q_H2O at 10:33:10 UTC total and of the seen facets as listed in "
         "summary.csv (its Q_source column: the run's inspection file; for C f10/1e-4 the final maps' manifest and the "
         "C report, c5 fit state). Reaching max_nfev is not convergence.", size=8.5)
    hdr = ["run", "model (summary.csv label)", "data cost new", "old law", "change", "last fit; fine tuning (nfev)",
           "rms ln NG all (late)", "rms ln RG all (late)", "multipliers 40-50 ... 90-180 deg", "kappa",
           "Q_H2O 10:33:10 total / seen [s^-1]", "flags"]
    rows = []
    for r in RK:
        rr, _ = srow(RUN[r]["csv"])
        lab = sx(r, "label", "summary.csv label", "%s", conv=None)
        term = rr["fit_termination"]
        ts = "max_nfev" if "maximum number" in term else ("ftol" if "ftol" in term else ("xtol" if "xtol" in term else term))
        sx(r, "fit_termination", "last fit termination", "%s", conv=None)
        nf = sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s)))
        seg = rr["fit_segments"]
        segt = (", %s segments" % sx(r, "fit_segments", "last fit segments", "%d", conv=lambda s: int(float(s)))
                if seg and int(float(seg)) > 1 else "")
        ft = ("; ft %s (%s)" % (sx(r, "ft_termination", "last fine tuning termination", "%s", conv=None),
                                sx(r, "ft_nfev", "last fine tuning nfev", "%d", conv=lambda s: int(float(s))))
              if rr["ft_termination"] else "; no fine tuning")
        ms = [sx(r, "mult_" + b, "multiplier " + b, "%.3g") for b in MULT]
        mtxt = ", ".join(ms) if all(m is not None for m in ms) else "none (D2: no slope multipliers)"
        flags = []
        if ts == "max_nfev":
            flags.append("last fit at max_nfev: not converged")
        if "1e-4" in r:
            flags.append("multipliers not determined (A12)")
        if sval(r, "kappa_at_bound") == "True":
            flags.append("kappa at bound")
        q = sx(r, "Q_H2O_1033", "Q_H2O at 10:33:10 [s^-1]", e2)
        qs = sx(r, "Q_H2O_seen_1033", "Q_H2O of the seen facets at 10:33:10 [s^-1]", e2)
        rows.append([r, lab, sx(r, "data_cost", "final data cost", e4), sx(r, "data_cost", "final data cost", e4, old=True),
                     sx(r, "data_cost_change", "data cost change new/old - 1", sgn_pct, scale=100.0),
                     "fit %s (%s%s)%s" % (ts, nf, segt, ft),
                     "%s (%s)" % (rms(r, "NG", False), rms(r, "NG", True)),
                     "%s (%s)" % (rms(r, "RG", False), rms(r, "RG", True)), mtxt,
                     sx(r, "kappa", "kappa (CO2 scale)", "%.2g"), "%s / %s" % (q, qs), "; ".join(flags)])
    data_table(hdr, rows, [0.75, 1.3, 0.62, 0.62, 0.5, 1.05, 0.7, 0.7, 1.05, 0.5, 0.9, 1.1], size=6.5)
    para("Parameters at bounds in the last fit (summary.csv 'at_bounds'): " + "; ".join(
        "%s: %s" % (r, sx(r, "at_bounds", "parameters at bounds (last fit)", "%s", conv=None) or "(none listed)")
        for r in RK) + ".", size=7.5)
    para("Run notes (summary.csv 'note'): " + "; ".join(
        "%s: %s" % (r, sx(r, "note", "note", "%s", conv=None)) for r in RK) + ".", size=7.5)
    if have_sets:
        para("Facts of the full figure sets (checks/facts.json of each set, full_set_facts.py; n/a: set or value not "
             "present):", size=8.5, bold_lead="")
        hdr = ["run", "seen area switched off (decision 1b)", "Q_H2O 10:33:10 [s^-1]", "from unseen facets", "seen production from >= 60 deg",
               "P50 NG union [m^2]", "small lobe late NG / RG (median)", "90-180 deg at last obs. NG / RG"]
        rows = []
        for r in RK:
            if not SETS.get(r):
                rows.append([r, "no full set"] + [""] * 6)
                continue
            rows.append([r, safe(off, r, what="switched off (decision 1b)"),
                         safe(fx, r, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"], "Q_H2O total 10:33:10", e2),
                         safe(fx, r, ["production_seen_unseen", "10:33:10", "unseen_share"], "unseen share of Q_H2O", pct,
                              100.0),
                         safe(f60, r), safe(fx, r, ["source_area_unions", "NG_P50", "area_m2"], "P50 NG union area", e2),
                         "%s / %s" % (safe(fx, r, ["regions", "small_lobe_median_late", 0], "small lobe NG late median",
                                           pct0, 100.0),
                                      safe(fx, r, ["regions", "small_lobe_median_late", 1], "small lobe RG late median",
                                           pct0, 100.0)),
                         "%s / %s" % (safe(fx, r, ["pressure_share_by_bin", "NG", "last", 6], "90-180 NG last obs.", pct0,
                                           100.0),
                                      safe(fx, r, ["pressure_share_by_bin", "RG", "last", 6], "90-180 RG last obs.", pct0,
                                           100.0))])
        data_table(hdr, rows, [0.9, 0.9, 1.0, 0.9, 1.0, 1.0, 1.1, 1.1], size=7)

    # ---------------- summary of the round
    N.ctx = "3. summary"
    new_section("L")
    doc.add_paragraph("2. Summary of the round (FINDINGS A18 and the files)", style="Heading 1")
    cost = {r: (float(sval(r, "data_cost")), sx(r, "data_cost", "final data cost", e4)) for r in RK}
    chg = {r: (float(sval(r, "data_cost_change")),
               sx(r, "data_cost_change", "data cost change new/old - 1", sgn_pct, scale=100.0)) for r in RK}
    maxn = [r for r in RK if "maximum number" in sval(r, "fit_termination")]
    better = [r for r in RK if chg[r][0] < 0]
    least = sorted([r for r in RK if chg[r][0] >= 0], key=lambda r: chg[r][0])[:3]
    para("With the intended law and the activity held fixed every old-law final state misfit the data by a factor 3.3 "
         "to 53 in data cost (A18), so every run was refitted. After the refits the final data costs are: %s. "
         "Change against the old-law counterpart: %s. Lower data cost than with the old law: %s; the smallest increases: "
         "%s. The last fits of %s stopped at max_nfev (not converged)." %
         (order_text(cost), order_text(chg), ", ".join(better) if better else "none", ", ".join(
             "%s %s" % (r, chg[r][1]) for r in least), ", ".join(maxn)), bold_lead="Fit quality.")
    q = {r: (float(sval(r, "Q_H2O_1033")), sx(r, "Q_H2O_1033", "Q_H2O at 10:33:10 [s^-1]", e2)) for r in RK}
    qs = {r: (float(sval(r, "Q_H2O_seen_1033")), sx(r, "Q_H2O_seen_1033", "Q_H2O of the seen facets at 10:33:10 [s^-1]", e2))
          for r in RK}
    ss, dk, dks = {}, {}, {}
    for r in RK:
        v = ixv(r, "103310", ["new", "seen_share_of_Q"])
        if v is not None:
            ss[r] = (v, ix(r, "103310", ["new", "seen_share_of_Q"], "seen share of Q_H2O at 10:33:10", pct0, 100.0))
            dk[r] = (ixv(r, "103310", ["new", "dark_share_of_Q"]),
                     ix(r, "103310", ["new", "dark_share_of_Q"], "night + cast-shadow share of Q_H2O at 10:33:10", pct0,
                        100.0))
            dks[r] = (ixv(r, "103310", ["new", "dark_share_of_seen_Q"]),
                      ix(r, "103310", ["new", "dark_share_of_seen_Q"],
                         "night + cast-shadow share of the seen Q_H2O at 10:33:10", pct0, 100.0))
    para("Q_H2O at 10:33:10 UTC spans %s s^-1 (%s), while the production of the facets a fitted observation sees is "
         "%s s^-1. Seen share of Q_H2O (inspection files; none for C f10/1e-4): %s. Night + cast-shadow share of "
         "Q_H2O: %s; of the seen Q_H2O: %s. Only option C with P_m 0.01 has its total mostly seen; in the other runs "
         "most of Q_H2O comes from terrain no fitted observation sees, so their totals are not set by the data (A15, "
         "A18)." % ("%s ... %s" % (min(q.values())[1], max(q.values())[1]), order_text(q), order_text(qs),
                    order_text(ss, reverse=True), order_text(dk, reverse=True), order_text(dks, reverse=True)),
         bold_lead="Source rates.")
    para("The slope bins with extra activity differ between the old-law and the new-law fit in every model with free "
         "multipliers (A18). This is not evidence of an effect of the law: one model has several local minima that differ "
         "in exactly this (A22). Last fits, new law [old law]:", bold_lead="Multiplier shift.")
    rows = []
    for r in RK:
        if sval(r, "mult_40_50") == "":
            continue
        rows.append([r] + ["%s [%s]" % (sx(r, "mult_" + b, "multiplier " + b, "%.3g"),
                                        sx(r, "mult_" + b, "multiplier " + b, "%.3g", old=True)) for b in MULT])
    data_table(["run"] + ["%s deg" % b.replace("_", "-") for b in MULT], rows, [1.2] + [1.35] * 6, size=7.5)
    if have_sets:
        big, sml, smr, p50, ovn, ovr = {}, {}, {}, {}, {}, {}
        for r in have_sets:
            try:
                big[r] = (fxv(r, ["regions", "big_lobe_NG_median_pre"]),
                          fx(r, ["regions", "big_lobe_NG_median_pre"], "big lobe NG median share before the late window",
                             pct0, 100.0))
                sml[r] = (fxv(r, ["regions", "small_lobe_median_late"])[0],
                          fx(r, ["regions", "small_lobe_median_late", 0], "small lobe NG late median", pct0, 100.0))
                smr[r] = (fxv(r, ["regions", "small_lobe_median_late"])[1],
                          fx(r, ["regions", "small_lobe_median_late", 1], "small lobe RG late median", pct0, 100.0))
                p50[r] = (fxv(r, ["source_area_unions", "NG_P50", "area_m2"]),
                          fx(r, ["source_area_unions", "NG_P50", "area_m2"], "P50 NG union area", e2))
                ovn[r] = (fxv(r, ["pressure_share_by_bin", "NG", "last"])[6],
                          fx(r, ["pressure_share_by_bin", "NG", "last", 6], "90-180 NG last obs.", pct0, 100.0))
                ovr[r] = (fxv(r, ["pressure_share_by_bin", "RG", "last"])[6],
                          fx(r, ["pressure_share_by_bin", "RG", "last", 6], "90-180 RG last obs.", pct0, 100.0))
            except (KeyError, IndexError, OSError) as e:
                MISSING_INPUTS.append(dict(ctx=N.ctx, what="facts.json summary", args=[r], error=str(e)))
        para("Before the late window the big lobe (body) gives a median %s of the modelled NG; in the late window the small "
             "lobe (head) gives NG %s and RG %s. The P50 source-area union of the 95 post-gap NG observations covers %s "
             "m^2. At the last observation (5.5 m) the 90-180 deg (overhang) bin gives %s of the modelled NG and %s of RG: "
             "the NG overhang dominance at touchdown appears in the models with activity on overhangs (multiplier or "
             "template); the run without templates has none and also misses the observed NG rise at the end, so the data "
             "do not show that it is unnecessary (FINDINGS A18, corrected by A21; checks/facts.json of each set)." % (
                 order_text(big, reverse=True), order_text(sml), order_text(smr), order_text(p50), order_text(ovn, reverse=True),
                 order_text(ovr, reverse=True)), bold_lead="Full figure sets.")
    if have_sets:
        para("A facet is active if its activity, not counting an imposed floor, exceeds T, the threshold of the project's "
             "reports: T = %s (a hundredth of the median srcCG activity of the seen facets). A fit fails if more than half "
             "of the seen area is switched off (user decision 1b, 2026-11-21; FINDINGS A21). The only imposed floor is "
             "option B's background. Only A and D pass, and they have the two highest data costs of the round." % (
                 N.add("all", "switched-off threshold T", jload(AA)["T"], AA, "active_area.json:T", "%.4g")),
             bold_lead="Switched-off surface.")
        data_table(["run", "switched off (decision 1b)", "with the floor counted", "at T/10", "at 10 T", "result"],
                   [[r, safe(off, r), safe(off, r, "off_T_floor_counted"), safe(off, r, "off_T0p1"),
                     safe(off, r, "off_T10"), "pass" if offv(r, "passes") else "fail"] for r in RK if SETS.get(r)],
                   [1.1, 1.2, 1.2, 1.0, 1.0, 0.8], size=7.5)
        para("The slope multipliers of every model that has them are bounded below by one (user decision 2026-11-12). "
             "No slope class can therefore be less active than flat terrain with the same underlying activity. A "
             "multiplier at or near that bound is not evidence that the class is as active as flat terrain, and a "
             "small production share of cliffs is not evidence that active cliffs are unnecessary: these fits did not "
             "test either.", bold_lead="Slopes.")
        msr = {r["variant"]: r for r in jload(MS_REP)["rows"]}
        mst = jload(MS_STA)
        para("Option C, f 10, P_m 0.01, templates off, was refitted from other starts (user decisions of 2026-11-21 and "
             "2026-11-22). Five solutions with all multipliers free and bounded below by one were found, and each puts "
             "nearly all the seen production on a different slope class, so the slope class is not measured. The run of "
             "this document is the one with the highest objective and data cost of the five. All five switch off most "
             "of the seen surface; with every multiplier held at one the same model does not, at a data cost below A and "
             "D. Whether cliffs are less active than flat terrain is still open: the fit with a lower bound below one "
             "did not reach its minimum.", bold_lead="Several solutions (FINDINGS A22).")
        rows = []
        for key, sk, lab in (("reference", "reference", "40-50 deg (this document's run)"),
                             ("ms_only_50_60", "ms_50_60", "50-60 deg"), ("ms_only_60_70", "ms_60_70", "60-70 deg"),
                             ("ms_only_70_80", "ms_70_80", "70-80 deg (lowest objective)"),
                             ("ms_only_80_90", "ms_80_90", "80-90 deg (lowest data cost)"),
                             ("none", "none", "all multipliers held at 1")):
            rr = msr[key]

            def num(field, fmt, scale=1.0):
                return N.add("C f10/0.01 off: " + key, "multi-start: " + field, rr[field] * scale, MS_REP,
                             "slopetest.json:rows[variant=%s]/%s" % (key, field), fmt)

            rows.append([lab, num("objective", e4), num("data_cost", e4),
                         num("top_class_share_of_seen_production", pct, 100.0) if key != "none" else "-",
                         "%s %s" % (num("off_seen_1b", pct, 100.0), "pass" if rr["passes_1b"] else "fail"),
                         num("Q_103310", e2), num("Q_seen_share", pct, 100.0), num("Q_share_over_F_E", pct, 100.0),
                         "local minimum" if mst[sk]["local_minimum"] else "nearly a minimum"])
        data_table(["slope class with the seen production", "objective", "data cost", "its share of seen production",
                    "seen area off (decision 1b)", "Q_H2O 10:33:10 [s^-1]", "seen share of Q", "Q above the energy limit",
                    "minimum test"], rows, [1.7, 0.85, 0.85, 0.9, 0.9, 0.9, 0.8, 0.9, 1.0], size=7)
    bullet("Not converged: the last fits of %s stopped at max_nfev; reaching max_nfev is not convergence." % ", ".join(maxn),
           lead="Caveats.")
    bullet("kappa at its lower bound in every run (summary.csv kappa_at_bound): the fits prefer the H2O law for the whole "
           "signal; total pressure only, so this does not measure CO2 as absent.")
    bullet("C with P_m 1e-4: the 70-80 / 80-90 deg multipliers are not determined (A12); their Q_H2O is set by unseen "
           "steep terrain times these factors, not by the data.")
    bullet("Q_H2O totals are an extrapolation of each model onto unseen terrain (A15 still applies with the new law, "
           "A18); the seen production is what the data constrain.")
    bullet("Observation from the figures: in D2's model_data figures (vs_COPS_altitude, _NG, _RG, _NG_RG, _time) the legend "
           "entry of the blue curve reads 'data cost' without the number (the captions give it); a label defect, as in the "
           "old-law D2 set.")
    bullet("C f10/1e-4: no old-law maps exist, so no inspection figure and no old-law Q_H2O; its seen Q_H2O in summary.csv "
           "is from the C report (c5 fit state). D2: one fit, no cycles, no fine tuning; its fine-tuning templates come "
           "from old-law fits (A18).")

    # ---------------- sections
    fig_no = 0
    snum = 2
    for sec, sname in SECTIONS:
        gs = [g for g in groups if g["section"] == sec]
        if not gs:
            continue
        snum += 1
        N.ctx = "section %s (intro)" % sec
        SRC["src"], SRC["fb"] = gsrc(gs[0]) + FR, None
        new_section("L")
        doc.add_paragraph("%d. %s" % (snum, sname), style="Heading 1")
        intro = SECTION_INTRO.get(sec)
        if sec in FS_SECTIONS and have_sets:
            fr = [os.path.join(SETS[r], "README.txt") for r in have_sets]
            ftxt = readme_folder_text(fr[0], sec)
            if ftxt:
                para(ftxt + " (README.txt of %s's full set; the same layout in every set)" % have_sets[0],
                     bold_lead="Folder (README.txt):")
        if intro:
            para(intro, bold_lead="About this section.")
        for k, g in enumerate(gs):
            bullet("Figure %d: %s%s" % (fig_no + 1 + k, g["title"], " (identical in all runs)" if g["identical"] else ""),
                   size=8.5)
        for gi, g in enumerate(gs):
            fig_no += 1
            g["fig"] = fig_no
            N.ctx = "Figure %d (%s)" % (fig_no, g["title"])
            SRC["src"], SRC["fb"] = gsrc(g), None
            L = g["layout"]
            new_section(L["orient"])
            doc.add_paragraph("%d.%d %s" % (snum, gi + 1, g["title"]), style="Heading 2")
            if "all" in g["paths"]:
                ftxt = os.path.relpath(g["paths"]["all"], ROOT)
            else:
                first = [p for p in g["paths"].values() if p][0]
                ftxt = os.path.relpath(first, ROOT) + " (and the other runs' files, Appendix B)"
            para("Files: " + ftxt, size=7.5, after=2, fb=("file path (Appendix B)", ROOT))
            show = g["show"]
            cols = L["cols"]
            nrow = math.ceil(len(show) / cols)
            W, H = usable(L["orient"])
            t = doc.add_table(rows=nrow, cols=cols)
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for k, who in enumerate(show):
                cell = t.rows[k // cols].cells[k % cols]
                cell.text = ""
                p0 = cell.paragraphs[0]
                p0.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p0.paragraph_format.space_after = Pt(1)
                p0.paragraph_format.keep_with_next = True
                lab = {"all": "all runs (one figure)" if g["section"] == "overview" else "all runs (identical files)"}.get(
                    who, who)
                r = p0.add_run(lab)
                r.bold = True
                r.font.size = Pt(9)
                r.font.color.rgb = RGBColor(0xA0, 0x10, 0x10)
                p1 = cell.add_paragraph()
                p1.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p1.paragraph_format.space_after = Pt(2)
                if g["paths"].get(who):
                    p1.add_run().add_picture(img[(g["id"], who)]["dst"], width=Inches(L["width_in"]))
                else:
                    why, src = why_of(g, who)
                    txt = "not produced for this run: %s%s" % (why, " (%s)" % os.path.relpath(src, ROOT) if src else "")
                    rr = p1.add_run(T(txt, [src] if src else [], ("reason of the missing plot", ROOT)))
                    rr.italic = True
                    rr.font.size = Pt(8)
            fixed(t, [W / cols] * cols)
            cell_margins(t, 40)
            cp = doc.add_paragraph(style=cap_style)
            miss = [w for w in g["show"] if w != "all" and not g["paths"].get(w)]
            who_txt = ("one multi-run figure" if g["section"] == "overview" else
                       "identical in all runs (one image shown)" if g["identical"] else
                       "all ten runs (D2, softplus, A, B, D, C f1/1e-4, C f10/1e-4, C f1/0.01, C f10/0.01, C f10/0.01 off)"
                       if not miss else "the ten runs in their fixed order, not produced for %s" % ", ".join(miss))
            rcap = cp.add_run(T("Figure %d: %s — %s; %s." % (fig_no, g["title"], who_txt, g["law"]), SRC["src"],
                                ("figure number / title", os.path.abspath(__file__))))
            rcap.font.size = Pt(9)
            rcap.bold = True
            doc.add_paragraph("Summary", style="Heading 3")
            summarize(g, para, data_table)
    # ---------------- appendix A: old-law results
    N.ctx = "Appendix A"
    SRC["src"], SRC["fb"] = FR, None
    new_section("L")
    doc.add_paragraph("Appendix A. The old-law (q_min = 0.2 q_max) results, for context", style="Heading 1")
    para("The old-law counterparts of the ten runs (summary.csv old_* rows; Q_H2O old law from the new rows' old_Q_H2O "
         "columns, i.e. the old-law blocks of the inspection files; D2's late rms from its residual-check json). These "
         "states used the unintended night-side law and are superseded (A18).", size=8.5)
    hdr = ["run", "old-law row (summary.csv label)", "data cost", "last fit; fine tuning (nfev)", "rms ln NG all (late)",
           "rms ln RG all (late)", "multipliers 40-50 ... 90-180 deg", "kappa", "Q_H2O 10:33:10 total / seen [s^-1]"]
    rows = []
    for r in RK:
        ro, _ = srow(row_key(r, True))
        term = ro["fit_termination"]
        ts = "max_nfev" if "maximum number" in term else ("ftol" if "ftol" in term else term)
        ft = ("; ft %s (%s)" % (sx(r, "ft_termination", "last fine tuning termination", "%s", conv=None, old=True),
                                sx(r, "ft_nfev", "last fine tuning nfev", "%d", conv=lambda s: int(float(s)), old=True))
              if ro["ft_termination"] else "; no fine tuning")
        ms = [sx(r, "mult_" + b, "multiplier " + b, "%.3g", old=True) for b in MULT]
        rows.append([r, sx(r, "label", "summary.csv label", "%s", conv=None, old=True),
                     sx(r, "data_cost", "final data cost", e4, old=True),
                     "fit %s (%s)%s" % (ts, sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s)), old=True),
                                        ft),
                     "%s (%s)" % (rms(r, "NG", False, True), rms(r, "NG", True, True)),
                     "%s (%s)" % (rms(r, "RG", False, True), rms(r, "RG", True, True)),
                     ", ".join(ms) if all(m is not None for m in ms) else "none",
                     sx(r, "kappa", "kappa (CO2 scale)", "%.2g", old=True) or "n/a",
                     "%s / %s" % (sx(r, "old_Q_H2O_1033", "old-law Q_H2O at 10:33:10", e2) or "n/a",
                                  sx(r, "old_Q_H2O_seen_1033", "old-law seen Q_H2O at 10:33:10", e2) or "n/a")])
    data_table(hdr, rows, [0.8, 1.7, 0.7, 1.3, 0.8, 0.8, 1.3, 0.6, 1.2], size=6.5)
    # ---------------- appendix B: sources
    N.ctx = "Appendix B"
    SRC["src"], SRC["fb"] = [], ("file path or checksum (Appendix B)", ROOT)
    new_section("L")
    doc.add_paragraph("Appendix B. Source files per figure (SHA-256)", style="Heading 1")
    para("Paths relative to %s/ (originals: .png as listed and .eps with the same name; captions <name>_caption.txt). "
         "Full SHA-256 of every source PNG. Full sets used: %s." % (ROOT, "; ".join(
             "%s: %s" % (r, os.path.relpath(SETS[r], ROOT) if SETS.get(r) else "none") for r in RK)), size=8)
    rows = []
    for g in groups:
        lines = [g["title"] + (" [identical in all runs]" if g["identical"] else "")]
        for who, p in g["paths"].items():
            if p:
                lines.append("%s: %s  sha256 %s" % (who, os.path.relpath(p, ROOT), g["sha"][who]))
            else:
                lines.append("%s: none (%s)" % (who, why_of(g, who)[0]))
        rows.append(["Figure %d" % g["fig"], lines])
    data_table(["figure", "plot and its source PNGs"], rows, [0.8, 8.8], size=6)
    # ---------------- appendix C: build
    N.ctx = "Appendix C"
    SRC["src"], SRC["fb"] = [os.path.abspath(__file__)], ("build parameter / build record", os.path.abspath(__file__))
    doc.add_paragraph("Appendix C. How this document was built", style="Heading 1")
    for txt in build_notes(args):
        bullet(txt, size=8.5)
    if N.unverified:
        with open(os.path.join(out, "numbers_unverified.json"), "w") as f:
            json.dump(N.unverified, f, indent=1, ensure_ascii=False)
        sys.exit("%d numeric tokens of hand-written text not found in a source file (numbers_unverified.json); "
                 "document not written" % len(N.unverified))
    doc.save(os.path.join(out, DOC_NAME + ".part"))
    os.replace(os.path.join(out, DOC_NAME + ".part"), os.path.join(out, DOC_NAME))


def readme_folder_text(readme, sec):
    try:
        lines = tload(readme).splitlines()
    except OSError:
        return ""
    out, on = [], False
    for ln in lines:
        if re.match(r"^%s/\s" % re.escape(sec), ln):
            on = True
            out.append(ln[len(sec) + 1:].strip())
            continue
        if on:
            if ln.startswith("    ") and ln.strip():
                out.append(ln.strip())
            else:
                break
    return " ".join(out)


SECTION_INTRO = {
    "overview": "The multi-run figures of the comparison set model_vs_observations_nonneg_options_stage2_q002/: the "
                "new-law best of each option (softplus reference, refitted D2, A, B, best C, D) and the five option-C "
                "chains with the softplus reference, each against the COPS data and as log residuals, and the data cost "
                "per cycle of all chains. One image per figure.",
    "mv": "Per run: the final state of the new-law run (solid) and of its old-law counterpart (dashed) against the COPS "
          "observations (934 NG, 1150 RG; the 95 late-window observations as the gauges' 10-s / 5-readout response "
          "average W10x5), the log residuals of both, and the data cost after every cycle's fit and fine tuning for both "
          "laws.",
    "model_data": "The model/data figures of each run's full set: the run's final state (new law) against COPS, with the "
                  "refitted D2 with the same law as the reference curve (README Notes of each set; in D2's own set both "
                  "curves are the same state).",
    "activity": "The fitted H2O activity factor A of each new-law final state on the seen facets (log classes per panel "
                "and run; light green = seen with A exactly 0, which is not the switched-off measure; grey = not seen "
                "during the fit), with the fit-independent "
                "effective gravitational slope and cliff templates. In each full set activity/ is a SHA-256-identical copy "
                "of source_rate_q002/<run>/activity/ (COPIED_FROM_source_rate_q002.txt); each run's figure is shown once.",
    "production_rate": "The H2O production rate F_H2O = q_H2O A per facet with the NEW law (q_min = 0.02 q_max on night "
                       "and cast-shadow facets, q_min + (q_max - q_min) cos SZA on day facets), same layouts as the "
                       "activity maps, relief shading under the colours. In each full set production_rate/ is a "
                       "SHA-256-identical copy of source_rate_q002/<run>/production_rate/ (COPIED_FROM_source_rate_q002.txt); "
                       "each run's figure is shown once.",
    "source_areas": "P% source areas of the 95 post-gap observations: for each observation and gauge the smallest set of "
                    "facets that supplies P % of the modelled surface pressure; union over the observations and "
                    "occurrence N_P. The close-up box is fitted per run, so the close-ups differ in extent.",
    "contribution_density": "Contribution of each facet per unit surface area to the modelled NG density and RG flux at "
                            "snapshots and as phase means, with the 50 % / 90 % source-area outlines; log colour range per "
                            "map; facets on overhangs or under other surface layers are not drawn in the lon/lat maps.",
    "cliffs": "How each state uses the effective-gravitational-slope bins of its multipliers (area, multipliers, "
              "production shares) and the share of the modelled surface pressure from each bin along the descent. The "
              "multipliers are bounded below by one, so these fits cannot show a class to be less active than flat "
              "terrain (Section 2, Slopes).",
    "regions": "Share of the modelled surface pressure from the small lobe (head, with the touchdown point) and the big "
               "lobe (body), and from distance shells around the touchdown point, along the descent.",
    "visibility": "Fit-independent figures (late-visibility templates, effective slope, cliff templates, ground track and "
                  "impact point); shown once where the PNGs are byte-identical in every set. Correction 2026-11-21: the "
                  "late-visibility templates belong to an earlier base with another gauge position and do not show which "
                  "facets each gauge saw in these fits (FINDINGS A20).",
    "inspection": "Per run and time: (a, b) the new-law F_H2O and (c, d) the old-law counterpart's F_H2O on ONE set of "
                  "colour classes (so the colours compare directly), (e, f) the new law's illumination factor q/q_max on "
                  "the whole surface. Not for C f10/1e-4 (no old-law maps). The full sets' inspection/ folders are "
                  "SHA-256-identical copies of these figures (COPIED_FROM_source_rate_q002.txt); shown once here.",
}


def build_notes(args):
    return [
        "Generator: %s (run copy; version-controlled copy model_ext/figures/q002_all_fits_summary_docx.py, SHA256 in "
        "model_ext/RUN_COPIES.sha256), derived from q002_comparison_docx.py and fullsets_comparison_docx.py; run with %s "
        "(python-docx %s, Pillow %s; separate venv, not the project .venv); observations file %s." % (
            os.path.abspath(__file__), sys.executable, args.docx_version, args.pil_version, args.observations),
        "Command: %s; Slurm job %s on %s; %s UTC." % (" ".join(sys.argv), os.environ.get("SLURM_JOB_ID", "-"),
                                                   os.uname().nodename, args.started),
        "Images: every source PNG hashed (SHA-256); a plot type whose PNGs are identical in all runs is shown once. "
        "Downscaled copies in images/: maps as JPEG (quality %d), line plots as PNG with 256 colours; width = min(native, "
        "max(%d px, 300 dpi at the printed width), 1400 px); RGBA flattened on white. Layout: per cluster the page "
        "orientation and grid (2-5 columns) that gives the largest images with all runs on one page (US Letter, margins "
        "0.6 in)." % (args.quality, args.min_px),
        "Numbers: read at build time from summary.csv, the residual-check json files, the inspection files, the full sets' "
        "facts.json, manifests, source_area_tables.txt, slope_bins json, the captions (regular expressions) and the "
        "comparison set's tables (residual extremes and last-observation values recomputed); listed in numbers_used.json "
        "with source and locator. Literal numbers of hand-written sentences are checked against the captions, READMEs, "
        "summary.txt, FINDINGS.md and OUTPUT_VISUALIZATION.md (the build refuses to write the document otherwise); numbers "
        "in 'Observation from the figures' paragraphs are recorded as figure readings; inputs an extractor could not find "
        "print 'n/a' and are listed in missing_inputs.json.",
        "Validation (validate step): the document reopened with python-docx (captions numbered and stating the law, images, "
        "placeholders, sections, tables, image and table widths against the pages, zip entries) and the OOXML schema check "
        "of the docx skill's validate.py; LibreOffice is not installed on the cluster, so no PDF; validation/ has layout "
        "previews drawn from the same geometry.",
    ]


# ------------------------------------------------------------------------------------------------ summaries
def summarize(g, para, data_table):
    F = FAMILY[g["fam"]]
    para(F["about"](g), bold_lead="What the plot shows.")
    missing = [r for r, p in g["paths"].items() if r != "all" and not p]
    if missing:
        parts = {}
        for r in missing:
            parts.setdefault(why_of(g, r), []).append(r)
        para(" ".join("Not produced for %s: %s%s." % (", ".join(rs), w, " (%s)" % os.path.relpath(s, ROOT) if s else "")
                      for (w, s), rs in parts.items()), bold_lead="Runs without this plot.",
             fb=("reason of the missing plot", ROOT))
    if g["identical"]:
        para("The PNGs of all runs are byte-identical (SHA-256 %s...): the figure does not depend on the fit or on the "
             "H2O law. One image is shown." % N.add("all", "SHA-256 of the identical PNGs", g["sha"][RK[0]],
                                                     g["paths"][RK[0]], "sha256 (equal for all runs)",
                                                     lambda v: v[:16]), bold_lead="Comparison.")
        return
    if F.get("numbers"):
        hdr, rows, w = F["numbers"](g)
        if hdr:
            data_table(hdr, rows, w)
    comp = F["compare"](g) if F.get("compare") else ""
    if comp:
        para(comp, bold_lead="How the runs compare.")
    o = OBS.get(g["cid"])
    if o:
        para(o, bold_lead="Observation from the figures:", force=True, fb=fig_fallback(g), src=[])


def first_path(g):
    return [x for x in g["paths"].values() if x][0]


def first_caption(g):
    return caption_of(first_path(g))


def cap_sentences(g, n=1):
    c = first_caption(g)
    if not c or not os.path.exists(c):
        return ""
    parts = re.split(r"(?<=[a-z0-9\)])\. (?=[A-Z(])", tload(c))
    return ". ".join(parts[:n]).rstrip(".") + "."


LAW_SENT = re.compile(r"NEW H2O law: .*?\(not intended\)")


def law_sentence(g):
    c = first_caption(g)
    m = LAW_SENT.search(tload(c)) if c and os.path.exists(c) else None
    return (m.group(0) + ".") if m else ""


# overview / per-run comparison (MV) ---------------------------------------------------------------------------------
TAB2RUN = {RUN[r]["csv"] + ".dat": r for r in RK}


def ov_fits(g):
    p = MV + "/12M_q002_%s_residuals_check.json" % g["set"]
    return [(TAB2RUN[os.path.basename(v["table"])], k, v, p) for k, v in jload(p)["fits"].items()]


def ov_about(g):
    if "residuals" in first_path(g):
        return (cap_sentences(g, 1) + " Log residuals ln(model / observation) of the runs in the table against the gauge "
                "altitude (logarithmic, descent from left to right): (a) NG, (b) RG; one point per fitted observation; for "
                "the 95 late-window observations the model is the gauges' W10x5 response average, as fitted.")
    return cap_sentences(g, 1) + " " + law_sentence(g)


def ov_numbers(g):
    hdr = ["run", "final data cost", "rms ln NG all (late)", "rms ln RG all (late)"]
    rows = []
    for r, k, v, p in ov_fits(g):
        rows.append([r, sx(r, "data_cost", "final data cost", e4)] + [
            "%s (%s)" % (N.add(r, "rms ln %s (all fitted rows)" % ga, v[ga]["rms"], p, "fits['%s'].%s.rms" % (k, ga), "%.4f"),
                         N.add(r, "rms ln %s (late window)" % ga, v[ga]["rms_late"], p, "fits['%s'].%s.rms_late" % (k, ga),
                               "%.4f")) for ga in ("NG", "RG")])
    return hdr, rows, [1.4, 1.4, 1.6, 1.6]


def ov_compare(g):
    fits = ov_fits(g)
    c = {r: (float(sval(r, "data_cost")), sx(r, "data_cost", "final data cost", e4)) for r, k, v, p in fits}
    ln = {r: (v["NG"]["rms_late"], N.add(r, "rms ln NG (late window)", v["NG"]["rms_late"], p,
                                          "fits['%s'].NG.rms_late" % k, "%.4f")) for r, k, v, p in fits}
    lr = {r: (v["RG"]["rms_late"], N.add(r, "rms ln RG (late window)", v["RG"]["rms_late"], p,
                                          "fits['%s'].RG.rms_late" % k, "%.4f")) for r, k, v, p in fits}
    return "Final data cost: %s. Late-window rms ln NG: %s; RG: %s." % (order_text(c), order_text(ln), order_text(lr))


def ovc_about(g):
    return cap_sentences(g, 2)


def ovc_numbers(g):
    hdr = ["run", "data cost after each fit (c1 ... c5)", "after each fine tuning", "final", "last fit"]
    rows = []
    for r in RK:
        rr, _ = srow(RUN[r]["csv"])
        nf = sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s)))
        if rr["history_fit_data"]:
            rows.append([r, sx(r, "history_fit_data", "data cost history, fits", "%s", conv=None),
                         sx(r, "history_ft_data", "data cost history, fine tunings", "%s", conv=None),
                         sx(r, "data_cost", "final data cost", e4),
                         "%s (%s)" % ("max_nfev" if "maximum" in rr["fit_termination"] else "ftol", nf)])
        else:
            rows.append([r, "one fit (no cycles)", "-", sx(r, "data_cost", "final data cost", e4), "max_nfev (%s)" % nf])
    return hdr, rows, [1.2, 3.0, 3.0, 1.0, 1.1]


def ovc_compare(g):
    c = {r: (float(sval(r, "data_cost")), sx(r, "data_cost", "final data cost", e4)) for r in RK}
    return "Final data cost: %s (the dotted line of the figure is the refitted D2)." % order_text(c)


def vo_about(g):
    return re.sub(r"^Modelled COPS pressures of [^:]*:", "Modelled COPS pressures of each run:", cap_sentences(g, 1)) + \
        " " + law_sentence(g)


def vo_numbers(g):
    fam = g["fam"]
    hdr = ["run", "data cost new (old; change)"]
    if fam in ("vo", "vo_ng"):
        hdr += ["rms ln NG all / late new (old)", "lowest ln NG 10-100 m new (old)", "ln NG at last obs. new (old)"]
    if fam in ("vo", "vo_rg"):
        hdr += ["rms ln RG all / late new (old)", "highest ln RG < 20 m new (old)", "ln RG at last obs. new (old)"]
    if fam == "vo_ngrg":
        hdr += ["lowest ln NG 10-100 m new (old)", "highest ln RG < 20 m new (old)"]
    rows = []
    for r in RK:
        row = [r, "%s (%s; %s)" % (sx(r, "data_cost", "final data cost", e4),
                                   sx(r, "data_cost", "final data cost", e4, old=True),
                                   sx(r, "data_cost_change", "data cost change new/old - 1", sgn_pct, scale=100.0))]
        if fam in ("vo", "vo_ng"):
            row += ["%s / %s (%s / %s)" % (rms(r, "NG", False), rms(r, "NG", True), rms(r, "NG", False, True),
                                           rms(r, "NG", True, True)),
                    "%s (%s)" % (t_min_ng(r)[1], t_min_ng(r, True)[1]),
                    "%s (%s)" % (t_last(r, "NG")[1], t_last(r, "NG", True)[1])]
        if fam in ("vo", "vo_rg"):
            row += ["%s / %s (%s / %s)" % (rms(r, "RG", False), rms(r, "RG", True), rms(r, "RG", False, True),
                                           rms(r, "RG", True, True)),
                    "%s (%s)" % (t_max_rg(r)[1], t_max_rg(r, True)[1]),
                    "%s (%s)" % (t_last(r, "RG")[1], t_last(r, "RG", True)[1])]
        if fam == "vo_ngrg":
            row += ["%s (%s)" % (t_min_ng(r)[1], t_min_ng(r, True)[1]), "%s (%s)" % (t_max_rg(r)[1], t_max_rg(r, True)[1])]
        rows.append(row)
    return hdr, rows, [1.1] + [1.6] * (len(hdr) - 1)


def vo_compare(g):
    fam = g["fam"]
    chg = {r: (float(sval(r, "data_cost_change")),
               sx(r, "data_cost_change", "data cost change new/old - 1", sgn_pct, scale=100.0)) for r in RK}
    worse = [r for r in RK if chg[r][0] > 0]
    better = [r for r in RK if chg[r][0] < 0]
    out = ["Data cost change new law vs old law: %s (worse in %d runs%s)." % (
        order_text(chg), len(worse), ", better in " + ", ".join(better) if better else "")]
    if fam in ("vo", "vo_ng", "vo_ngrg"):
        out.append("ln(model/obs) NG at the last observation, new law: %s." % order_text({r: t_last(r, "NG") for r in RK}))
    if fam in ("vo", "vo_rg", "vo_ngrg"):
        v = {r: t_last(r, "RG") for r in RK}
        vo_ = {r: t_last(r, "RG", True) for r in RK}
        lower = [r for r in RK if v[r][0] < vo_[r][0]]
        out.append("ln(model/obs) RG at the last observation, new law: %s; lower than with the old law in %s." % (
            order_text(v), "all ten runs" if len(lower) == len(RK) else (", ".join(lower) if lower else "no run")))
    return " ".join(out)


def vr_about(g):
    return (cap_sentences(g, 1) + " Log residuals ln(model / observation) of the new-law run and of its old-law "
            "counterpart against the gauge altitude (logarithmic, descent from left to right): (a) NG, (b) RG; one point "
            "per fitted observation (934 NG, 1150 RG), the late window shaded; for the 95 late-window observations the "
            "model is the gauges' W10x5 response average, as fitted. Positive: model above the data.")


def vr_numbers(g):
    hdr = ["run", "rms ln NG all (late) new", "old", "rms ln RG all (late) new", "old"]
    rows = []
    for r in RK:
        if not g["paths"][r]:
            rows.append([r, "no figure", "", "", ""])
            continue
        rows.append([r, "%s (%s)" % (rms(r, "NG", False), rms(r, "NG", True)),
                     "%s (%s)" % (rms(r, "NG", False, True), rms(r, "NG", True, True)),
                     "%s (%s)" % (rms(r, "RG", False), rms(r, "RG", True)),
                     "%s (%s)" % (rms(r, "RG", False, True), rms(r, "RG", True, True))])
    return hdr, rows, [1.3, 1.5, 1.5, 1.5, 1.5]


def vr_compare(g):
    have = [r for r in RK if g["paths"][r]]
    v = {r: t_last(r, "NG") for r in have}
    vo_ = {r: t_last(r, "NG", True) for r in have}
    rn = {r: t_last(r, "RG") for r in have}
    ro = {r: t_last(r, "RG", True) for r in have}
    lower = [r for r in have if rn[r][0] < ro[r][0]]
    return ("ln(model/obs) NG at the last observation, new law: %s; old law: %s. RG at the last observation, new law: "
            "%s; old law: %s; the new-law value is lower in %s (tables/)." %
            (order_text(v), order_text(vo_), order_text(rn), order_text(ro),
             "all %d runs shown" % len(have) if len(lower) == len(have) else (", ".join(lower) or "none")))


def cv_about(g):
    m = re.search(r"^Convergence of the stage-2 chain \S+: (.*?); solid black", tload(first_caption(g)), re.S)
    return "Convergence of each stage-2 chain: " + (m.group(1) if m else "") + ". " + (
        "Solid black: the NEW H2O law (q_min = 0.02 q_max); dashed grey: the OLD-law chain (q_min = 0.2 q_max); both "
        "against the same observations; each panel has its own data-cost axis.")


def cv_numbers(g):
    hdr = ["run", "fits c1 ... c5, new law", "fine tunings c1 ... c5, new law", "final new", "final old", "last fit new"]
    rows = []
    for r in RK:
        if not g["paths"][r]:
            rows.append([r, "no chain (one fit)", "", sx(r, "data_cost", "final data cost", e4),
                         sx(r, "data_cost", "final data cost", e4, old=True), ""])
            continue
        rr, _ = srow(RUN[r]["csv"])
        rows.append([r, sx(r, "history_fit_data", "data cost history, fits", "%s", conv=None),
                     sx(r, "history_ft_data", "data cost history, fine tunings", "%s", conv=None),
                     sx(r, "data_cost", "final data cost", e4), sx(r, "data_cost", "final data cost", e4, old=True),
                     "%s (%s)" % ("max_nfev" if "maximum" in rr["fit_termination"] else "ftol",
                                  sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s))))])
    return hdr, rows, [1.1, 2.6, 2.6, 0.9, 0.9, 1.0]


def cv_compare(g):
    drop = {}
    for r in RK:
        if not g["paths"][r]:
            continue
        rr, p = srow(RUN[r]["csv"])
        h = [float(x) for x in rr["history_fit_data"].split()]
        fin = float(rr["data_cost"])
        drop[r] = (fin / h[0] - 1.0, N.add(r, "change of the data cost from the c1 fit to the final state",
                                          100.0 * (fin / h[0] - 1.0), p,
                                          "summary.csv: data_cost / first value of history_fit_data - 1", sgn_pct))
    lower = [r for r in drop if float(sval(r, "data_cost")) < float(sval(r, "data_cost", old=True))]
    return ("Change of the data cost from the cycle-1 fit to the final state, new law: %s. Final new-law cost below the "
            "old-law one: %s." % (order_text(drop), ", ".join(lower) if lower else "none"))


# inspection (SR) ------------------------------------------------------------------------------------------------------
def tt(g):
    return "%s:%s:%s" % (g["t"][:2], g["t"][2:4], g["t"][4:])


def insp_about(g):
    t = tload(first_caption(g))
    return re.split(r"(?<=[a-z0-9\)])\. (?=[A-Z(])", t)[0] + (
        ". (a, b) F_H2O of the new-law run with the NEW law q_min = 0.02 q_max; (c, d) F_H2O of its old-law counterpart "
        "with the OLD law q_min = 0.2 q_max, both on ONE set of log10 classes per window (colours compare directly); "
        "(e, f) the new law's illumination factor q_H2O/q_max on the whole surface (0.02 on night and cast-shadow "
        "facets).")


def insp_numbers(g):
    t = g["t"]
    hdr = ["run", "Q_H2O new (old) [s^-1]", "new / old", "seen Q new (old)", "seen share new (old)",
           "night + shadow share of Q new (old)", "of the seen Q new (old)"]
    rows = []
    for r in RK:
        if not g["paths"][r]:
            rows.append([r, "no figure", "", "", "", "", ""])
            continue
        row = [r]
        for key, nm, fmt, sc in (("Q_H2O_s-1", "Q_H2O total", e3, 1.0), (None, None, None, None),
                                 ("Q_H2O_seen_s-1", "Q_H2O of the seen facets", e3, 1.0),
                                 ("seen_share_of_Q", "seen share of Q_H2O", pct, 100.0),
                                 ("dark_share_of_Q", "night + cast-shadow share of Q_H2O", pct, 100.0),
                                 ("dark_share_of_seen_Q", "night + cast-shadow share of the seen Q_H2O", pct, 100.0)):
            if key is None:
                qn_, qo = ixv(r, t, ["new", "Q_H2O_s-1"]), ixv(r, t, ["old", "Q_H2O_s-1"])
                row.append(N.add(r, "Q_H2O new / old at %s" % tt(g), qn_ / qo, insp_path(r, t),
                                 "inspect_%s.json: new/Q_H2O_s-1 / old/Q_H2O_s-1" % t, "%.3f"))
                continue
            row.append("%s (%s)" % (ix(r, t, ["new", key], nm + " at %s (new law)" % tt(g), fmt, sc),
                                    ix(r, t, ["old", key], nm + " at %s (old law)" % tt(g), fmt, sc)))
        rows.append(row)
    return hdr, rows, [1.1, 1.5, 0.7, 1.5, 1.2, 1.4, 1.3]


def insp_compare(g):
    t = g["t"]
    rat, sr = {}, {}
    for r in RK:
        if g["paths"][r]:
            qn_, qo = ixv(r, t, ["new", "Q_H2O_s-1"]), ixv(r, t, ["old", "Q_H2O_s-1"])
            rat[r] = (qn_ / qo, N.add(r, "Q_H2O new / old at %s" % tt(g), qn_ / qo, insp_path(r, t),
                                      "inspect_%s.json: new/Q_H2O_s-1 / old/Q_H2O_s-1" % t, "%.3f"))
            a, b = ixv(r, t, ["new", "Q_H2O_seen_s-1"]), ixv(r, t, ["old", "Q_H2O_seen_s-1"])
            sr[r] = (a / b, N.add(r, "seen Q_H2O new / old at %s" % tt(g), a / b, insp_path(r, t),
                                  "inspect_%s.json: new/Q_H2O_seen_s-1 / old/Q_H2O_seen_s-1" % t, "%.3f"))
    more = [r for r in rat if rat[r][0] < sr[r][0]]
    return ("Q_H2O new law / old law at %s UTC: %s; seen Q_H2O new / old: %s. The total falls more than the seen "
            "production in %d of the %d runs shown%s." % (tt(g), order_text(rat), order_text(sr), len(more), len(rat),
                                                          "" if len(more) == len(rat) else " (" + ", ".join(more) + ")"))


# full sets ------------------------------------------------------------------------------------------------------------
def fx(run, keys, quantity, fmt="%s", scale=1.0):
    p = os.path.join(SETS[run], "checks/facts.json")
    v = nested(jload(p), keys)
    if scale != 1.0:
        v = v * scale
    return N.add(run, quantity, v, p, "facts.json:" + "/".join(str(k) for k in keys), fmt)


def fxv(run, keys):
    return nested(jload(os.path.join(SETS[run], "checks/facts.json")), keys)


def f60(run):
    p = os.path.join(SETS[run], "checks/facts.json")
    d = fxv(run, ["production_share_seen_103310"])
    v = sum(d[b] for b in BINS[3:])
    return N.add(run, "share of the seen H2O production at 10:33:10 from slope >= 60 deg", 100.0 * v, p,
                 "facts.json: sum of production_share_seen_103310 [60,70) ... [90,180]", pct)


def have(g):
    return [r for r in RK if g["paths"].get(r)]


# model_data of the full sets
def fmd_about(g):
    if g["fam"] == "md_pub":
        return ("COPS NG and RG and the model NG and RG of each run against the altitude above the surface; the gap in the "
                "NG data hatched; publication format of pub_ng_rg_altitude.py (no late-window shading, user decision "
                "2026-11-12; OUTPUT_VISUALIZATION.md). New law q_min = 0.02 q_max.")
    s = cap_sentences(g, 1).replace(STEMS.get(have(g)[0]) or "<STEM>", "<STEM>")
    c = tload(first_caption(g))
    m = re.search(r"Orange dashed: (.*?)(?:, final data cost|\. )", c)
    ref = (" Reference curve (orange dashed, caption of %s): %s." % (have(g)[0], m.group(1))) if m else ""
    return s + ref + " Blue: the run's new-law final state."


def fmd_numbers(g):
    fam = g["fam"]
    hdr = ["run", "final data cost", "rms ln NG all (late)", "rms ln RG all (late)", "lowest ln NG 10-100 m",
           "highest ln RG < 20 m"]
    if fam in ("md_alt", "md_time", "md_ngrg"):
        hdr.append("model NG = RG crossing (caption)")
    rows = []
    for r in RK:
        row = [r, sx(r, "data_cost", "final data cost", e4), "%s (%s)" % (rms(r, "NG", False), rms(r, "NG", True)),
               "%s (%s)" % (rms(r, "RG", False), rms(r, "RG", True)), t_min_ng(r)[1], t_max_rg(r)[1]]
        if fam in ("md_alt", "md_time", "md_ngrg"):
            p = g["paths"].get(r)
            row.append(capx(caption_of(p), r"the fitted model at (\d+) m", r, "model NG = RG crossing altitude [m]") + " m"
                       if p and os.path.exists(caption_of(p)) else "n/a")
        rows.append(row)
    return hdr, rows, [1.2] + [1.4] * (len(hdr) - 1)


def fmd_compare(g):
    ln = {r: (float(sval(r, "rms_NG_late")) if sval(r, "rms_NG_late") else None, rms(r, "NG", True)) for r in RK}
    lr = {r: (float(sval(r, "rms_RG_late")) if sval(r, "rms_RG_late") else None, rms(r, "RG", True)) for r in RK}
    v = {r: t_min_ng(r) for r in RK}
    w = {r: t_max_rg(r) for r in RK}
    return ("Late-window rms ln NG: %s; RG: %s (D2's late values are not in summary.csv, see the table). Lowest "
            "ln(model/obs) NG between 100 and 10 m: %s. Highest ln(model/obs) RG below 20 m: %s (comparison-set tables "
            "of the same states)." % (order_text(ln), order_text(lr), order_text(v, reverse=True), order_text(w)))


# maps of the full sets
PANEL_NAME = {"latlon_global": "lon/lat global", "latlon_zoom": "lon/lat impact area", "3d_a": "3-D (a) nucleus",
              "3d_b": "3-D (b) impact region", "3d_c": "3-D (c) late descent"}


def man_path(r, g):
    return os.path.join(SETS[r], "checks/activity_manifest.json" if g["q"] == "activity" else
                        "checks/h2o_%s_manifest.json" % g["t"])


def cls_range(r, g, panel):
    p = man_path(r, g)
    ct = jload(p)["class_tables"][panel]
    lo, hi, n = ct["classes"][0]["log10_from"], ct["classes"][-1]["log10_to"], len(ct["classes"])
    s = N.add(r, "colour classes %s %s: log10 range and count" % (g["q"] + (g["t"] or ""), panel), [lo, hi, n], p,
              "class_tables/%s/classes[0].log10_from, [-1].log10_to, len" % panel,
              lambda v: "10^%g-10^%g (%d)" % tuple(v))
    return s, hi - lo


def qcap(r, g):
    p = os.path.join(SETS[r], "production_rate/12M_templates_h2o_production_rate_%sUTC_caption.txt" % g["t"])
    return capx(p, r"total H2O production Q_H2O = sum F_H2O area = (\S+) s\^-1", r,
                "Q_H2O total at %s [s^-1] (caption)" % tt(g))


def map_about(g):
    what = "H2O activity factor A" if g["q"] == "activity" else ("H2O production rate F_H2O (instantaneous emitted H2O "
                                                                "molecular flux per facet, m^-2 s^-1) at %s UTC" % tt(g))
    first = {"map_tpl": "%s of each run's new-law final state, with the effective gravitational slope and the cliff "
                        "templates, on the cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf shape model (12,364,290 facets), "
                        "longitude/latitude: left global, right the impact area" % what,
             "map_latlon": "%s of each run's new-law final state with the end-of-mission ground track of Rosetta and its "
                           "impact point, longitude/latitude: (a) global map, (b) the impact area" % what,
             "map_3d": "%s of each run's new-law final state with the ground track, orthographic views rendered on the "
                       "GPU: (a) the whole nucleus, (b) the impact region, (c) close-up of the late descent with "
                       "1-minute marks" % what,
             "map_panel": "%s of each run's new-law final state: %s" % (what, g["title"].split(": ")[-1])}[g["fam"]]
    pat = r"'seen': the facet contributes to at least one fitted NG or RG observation; ([\d,]+ of the [\d,]+ facets, [\d.]+ % of the area)"
    seen = {r: capx(caption_of(g["paths"][r]), pat, r, "seen facets during the fit (caption)") for r in have(g)
            if os.path.exists(caption_of(g["paths"][r]))}
    same = len(set(seen.values())) == 1
    law = " F_H2O = q_H2O A with the NEW law q_min = 0.02 q_max (captions)." if g["q"] != "activity" else ""
    return (first + ". Colour only on the facets seen during the fit (%s%s); light green: seen with %s exactly 0 (not the switched-off measure); grey: not "
            "seen. Classes scanned per panel and per run (captions; checks/*_class_tables.txt), so the colour scales "
            "differ between the runs.%s" % (list(seen.values())[0] if seen else "n/a",
                                             ", the same in all runs" if same else "; differs between runs",
                                             "A" if g["q"] == "activity" else "F_H2O", law))


def map_numbers(g):
    hdr = ["run"]
    if g["q"] == "activity":
        hdr += ["switched off, seen area (decision 1b)", "the same, floor counted"]
    else:
        hdr += ["Q_H2O total [s^-1] (caption)", "seen share of Q", "night + shadow share of Q / of seen Q"]
    for k, panel in g["panels"]:
        hdr.append("classes %s: range (n)" % PANEL_NAME[panel])
    if g["fam"] == "map_tpl" and g["q"] != "activity":
        hdr.append("seen production from >= 60 deg")
    rows = []
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r] + ["no figure"] + [""] * (len(hdr) - 2))
            continue
        row = [r]
        if g["q"] == "activity":
            row.append(safe(off, r, what="switched off (decision 1b)"))
            row.append(safe(off, r, "off_T_floor_counted", what="switched off, floor counted"))
        else:
            row.append(safe(qcap, r, g))
            s = ix(r, g["t"], ["new", "seen_share_of_Q"], "seen share of Q_H2O at %s" % tt(g), pct, 100.0)
            row.append(s if s else "n/a (no inspection)")
            d1 = ix(r, g["t"], ["new", "dark_share_of_Q"], "night + cast-shadow share of Q_H2O at %s" % tt(g), pct, 100.0)
            d2 = ix(r, g["t"], ["new", "dark_share_of_seen_Q"],
                    "night + cast-shadow share of the seen Q_H2O at %s" % tt(g), pct, 100.0)
            row.append("%s / %s" % (d1, d2) if d1 else "n/a")
        for k, panel in g["panels"]:
            row.append(safe(lambda: cls_range(r, g, panel)[0], what="class table %s" % panel))
        if g["fam"] == "map_tpl" and g["q"] != "activity":
            def s60():
                p = os.path.join(SETS[r], "checks/slope_bins_h2o_%s.json" % g["t"])
                v = sum(b["production_share_seen"] for b in jload(p)["bins"][3:])
                return N.add(r, "seen H2O production share from >= 60 deg at %s" % tt(g), 100.0 * v, p,
                             "sum bins[3:].production_share_seen", pct)
            row.append(safe(s60, what="slope bins"))
        rows.append(row)
    return hdr, rows, [1.1] + [1.5] * (len(hdr) - 1)


def map_compare(g):
    out = []
    hv = have(g)
    if g["q"] == "activity":
        az = {}
        for r in hv:
            try:
                az[r] = (offv(r), off(r))
            except (KeyError, OSError):
                pass
        if az:
            out.append("Share of the seen area switched off (decision 1b; the light green of the maps is exact zero only): %s."
                       % order_text(az, reverse=True))
    else:
        q = {}
        for r in hv:
            s = safe(qcap, r, g)
            if s != "n/a":
                q[r] = (float(s), s)
        ss = {}
        for r in hv:
            v = ixv(r, g["t"], ["new", "seen_share_of_Q"])
            if v is not None:
                ss[r] = (v, ix(r, g["t"], ["new", "seen_share_of_Q"], "seen share of Q_H2O at %s" % tt(g), pct, 100.0))
        out.append("Q_H2O at %s UTC: %s s^-1; seen share: %s." % (tt(g), order_text(q), order_text(ss, reverse=True)))
    k, panel = g["panels"][0]
    spans = {}
    for r in hv:
        try:
            s, d = cls_range(r, g, panel)
            spans[r] = (d, N.add(r, "decades covered by the classes, %s %s" % (g["q"], panel), d, man_path(r, g),
                                 "class_tables/%s: log10_to(last) - log10_from(first)" % panel, "%g"))
        except (KeyError, OSError):
            pass
    if spans:
        out.append("Decades spanned by the colour classes in %s: %s." % (PANEL_NAME[panel], order_text(spans)))
    return " ".join(out)


def geom_about(g):
    return cap_sentences(g, 1) + (" Effective gravitational slope: angle between the inward facet normal and the "
                                  "effective gravity; values above 90 deg are overhangs. Cliff templates: slope >= 45 "
                                  "deg (light) and >= 60 deg (dark). Geometry only: independent of the fit and of the H2O "
                                  "law.")


# source areas
def sa_val(r, gauge, P, what):
    p = os.path.join(SETS[r], "checks/source_area_tables.txt")
    pat = r"\[%s, P = %d %%\] facets with N > 0: (\d+) \(= union: (\w+)\), max N (\d+); area (\S+) m\^2" % (gauge, P)
    m = re.search(pat, tload(p))
    grp = {"facets": 1, "maxN": 3, "area": 4}[what]
    v = float(m.group(grp)) if what == "area" else int(m.group(grp))
    fmt = e2 if what == "area" else (intc if what == "facets" else "%d")
    return v, N.add(r, "source area %s P%d %s" % (gauge, P, what), v, p, "regex group %d: %s" % (grp, pat), fmt)


def sa_about(g):
    if g["fam"] == "sa_union":
        return ("U_P = the facets in the P % source area S_P of at least one of the 95 post-gap observations (10:33:10-"
                "10:39:26 UTC), P = 50 (darkest) to 90 % (lightest), each facet coloured by the smallest P whose union "
                "contains it; global map and a close-up around the impact (box fitted to the union facets). S_P(k, g): "
                "the smallest set of facets, in decreasing order of their contribution to the modelled surface pressure "
                "of observation k at gauge g, that supplies P % of it; NG grain term excluded. New law q_min = 0.02 "
                "q_max.")
    return ("N_%d = the number of the 95 post-gap observations whose %d %% source area contains the facet (classes 1, 2-4, "
            "5-9, ..., 80-94, 95, the same in every run; grey N = 0); global map and close-up around the impact (box "
            "fitted per run). New law q_min = 0.02 q_max." % (g["P"], g["P"]))


def sa_numbers(g):
    hdr, rows = ["run"], []
    for ga in g["gauges"]:
        if g["fam"] == "sa_union":
            hdr += ["%s U50 facets / area [m^2]" % ga, "%s U90 facets / area [m^2]" % ga]
        else:
            hdr += ["%s facets with N_%d > 0" % (ga, g["P"]), "%s area [m^2]" % ga, "%s max N" % ga]
    hdr.append("close-up box (caption)")
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r, "no figure"] + [""] * (len(hdr) - 2))
            continue
        row = [r]
        for ga in g["gauges"]:
            if g["fam"] == "sa_union":
                for P in (50, 90):
                    row.append(safe(lambda: "%s / %s" % (sa_val(r, ga, P, "facets")[1], sa_val(r, ga, P, "area")[1]),
                                    what="source_area_tables"))
            else:
                for wh in ("facets", "area", "maxN"):
                    row.append(safe(lambda: sa_val(r, ga, g["P"], wh)[1], what="source_area_tables"))
        row.append(capx(caption_of(g["paths"][r]), r"box in the global map: (lon \S+ deg, lat \S+ deg)", r,
                        "close-up box (lon, lat)"))
        rows.append(row)
    return hdr, rows, [1.2] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def sa_compare(g):
    out = []
    P = 50 if g["fam"] == "sa_union" else g["P"]
    for ga in g["gauges"]:
        a = {}
        for r in have(g):
            try:
                a[r] = sa_val(r, ga, P, "area")
            except (AttributeError, OSError):
                pass
        out.append("%s P%d area: %s m^2." % (ga, P, order_text(a)))
    gr = {}
    if "NG" in g["gauges"]:
        for r in have(g):
            gr[r] = capx(caption_of(g["paths"][r]), r"its share of the model NG pressure is ([\d.]+-[\d.]+, median [\d.]+)",
                         r, "NG grain-term share of the model NG pressure (range, median)")
        out.append("NG grain-term share of the model NG in these observations (range, median; excluded from the source "
                   "areas): " + "; ".join("%s %s" % (r, gr[r]) for r in have(g)) + ".")
    return " ".join(out)


# contribution density
def cd_v(r, key, ga, field, quantity, fmt, scale=1.0):
    return manx(os.path.join(SETS[r], "contribution_density/contribution_density_manifest.json"),
                ["ranges", key, ga, field], r, quantity, fmt, scale)


def cd_about(g):
    t = tload(first_caption(g))
    first = t.split(". Colour:")[0]
    return (first + ". Colour: the contribution per unit surface area as the fitted model computes it (NG: density at "
            "the nude gauge per m^2 of surface; RG: flux into the ram gauge per m^2), log scale per map; white contours: "
            "50 % (solid) and 90 % (dashed) source areas of the map; close-up box fitted per run. New law q_min = 0.02 "
            "q_max.")


def cd_numbers(g):
    hdr = ["run"]
    for ga in g["gauges"]:
        hdr += ["%s P50 facets / area [m^2]" % ga, "%s P90 facets / area [m^2]" % ga, "%s signal drawn in close-up" % ga]
    rows = []
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r, "no figure"] + [""] * (len(hdr) - 2))
            continue
        row = [r]
        for ga in g["gauges"]:
            k = g["key"]
            row.append(safe(lambda: "%s / %s" % (cd_v(r, k, ga, "P50_facets", "%s %s P50 facets" % (k, ga), intc),
                                                 cd_v(r, k, ga, "P50_area_m2", "%s %s P50 area [m^2]" % (k, ga), e2)),
                            what="cd manifest"))
            row.append(safe(lambda: "%s / %s" % (cd_v(r, k, ga, "P90_facets", "%s %s P90 facets" % (k, ga), intc),
                                                 cd_v(r, k, ga, "P90_area_m2", "%s %s P90 area [m^2]" % (k, ga), e2)),
                            what="cd manifest"))
            row.append(safe(lambda: cd_v(r, k, ga, "signal_share_drawn_zoom",
                                         "%s %s share of the signal drawn in the close-up" % (k, ga), pct0, 100.0),
                            what="cd manifest"))
        rows.append(row)
    return hdr, rows, [1.2] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def cd_compare(g):
    out = []
    for ga in g["gauges"]:
        a, d = {}, {}
        for r in have(g):
            try:
                p = os.path.join(SETS[r], "contribution_density/contribution_density_manifest.json")
                e = jload(p)["ranges"][g["key"]][ga]
                nf = cd_v(r, g["key"], ga, "P50_facets", "%s %s P50 facets" % (g["key"], ga), intc)
                a[r] = (e["P50_area_m2"], "%s (%s facet%s)" % (
                    cd_v(r, g["key"], ga, "P50_area_m2", "%s %s P50 area [m^2]" % (g["key"], ga), e2), nf,
                    "" if nf == "1" else "s"))
                d[r] = (e["signal_share_drawn_zoom"], cd_v(r, g["key"], ga, "signal_share_drawn_zoom",
                                                           "%s %s share of the signal drawn in the close-up" % (g["key"], ga),
                                                           pct0, 100.0))
            except (KeyError, OSError):
                pass
        out.append("%s: P50 area of this map %s m^2; share of the map's signal drawn in the close-up %s (the rest lies "
                   "outside the close-up, on overhangs or under other surface layers)." % (ga, order_text(a), order_text(d)))
    return " ".join(out)


# cliffs
def cliff_bins_about(g):
    s = tload(first_caption(g)).split(". (a)")[0]
    s = re.sub(r" \(the [^()]*(\([^()]*\)[^()]*)*\)", "", s, count=1)
    return (s + ": (a) share of the total and of the seen area per bin and the seen fraction of each bin, (b) the fitted "
            "multiplier m_b (log scale; 0-40 deg held at 1), (c) share of the H2O production of the seen surface (solid) "
            "and of the whole surface (hatched) per bin at 10:33:10 and 10:39:26 UTC. New law q_min = 0.02 q_max.")


def cliff_bins_numbers(g):
    hdr = ["run"] + ["%s deg: m / seen prod." % b for b in BIN_LABEL]
    rows = []
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r, "no figure"] + [""] * 6)
            continue
        p = os.path.join(SETS[r], "checks/slope_bins_h2o_103310.json")
        row = [r]
        for i in range(7):
            row.append(safe(lambda: "%s / %s" % (
                manx(p, ["bins", i, "multiplier"], r, "multiplier bin %s" % BIN_LABEL[i], "%.3g"),
                manx(p, ["bins", i, "production_share_seen"], r, "seen H2O production share 10:33:10 bin %s" % BIN_LABEL[i],
                     pct, 100.0)), what="slope_bins"))
        rows.append(row)
    return hdr, rows, [1.2] + [1.2] * 7


def cliff_bins_compare(g):
    st, al, parts = {}, {}, []
    for r in have(g):
        try:
            p = os.path.join(SETS[r], "checks/slope_bins_h2o_103310.json")
            bins = jload(p)["bins"]
        except OSError:
            continue
        v = sum(b["production_share_seen"] for b in bins[3:])
        st[r] = (v, N.add(r, "seen H2O production share 10:33:10 from >= 60 deg", 100.0 * v, p,
                          "sum bins[3:].production_share_seen", pct))
        w = sum(b["production_share_all"] for b in bins[3:])
        al[r] = (w, N.add(r, "whole-surface H2O production share 10:33:10 from >= 60 deg", 100.0 * w, p,
                          "sum bins[3:].production_share_all", pct))
        i = max(range(7), key=lambda j: bins[j]["production_share_seen"])
        parts.append("%s %s deg (%s)" % (r, BIN_LABEL[i], N.add(r, "largest bin of the seen production at 10:33:10",
                                                                100.0 * bins[i]["production_share_seen"], p,
                                                                "bins[%d].production_share_seen (max)" % i, pct0)))
    return ("Seen production from >= 60 deg at 10:33:10 (solid bars of panel c): %s; whole surface (hatched): %s. Largest "
            "bin of the seen production: %s. Panel (b) of D2 has no multipliers (all 1)." %
            (order_text(st, reverse=True), order_text(al, reverse=True), "; ".join(parts)))


def cliff_pf_about(g):
    return ("Fraction of the modelled surface pressure contributed by each effective-slope bin along the descent (stacked; "
            "one colour per bin), against the gauge altitude; late window boxed; NG gap hatched in the NG panel; NG grain "
            "term not included. New law q_min = 0.02 q_max.")


def cliff_pf_numbers(g):
    hdr = ["run"]
    for ga in g["gauges"]:
        hdr += ["%s >= 60 deg before late window (caption)" % ga, "%s >= 60 deg in late window (caption)" % ga,
                "%s mean >= 60 deg late (facts)" % ga, "%s 90-180 deg at last obs. (facts)" % ga]
    rows = []
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r, "no figure"] + [""] * (len(hdr) - 2))
            continue
        row = [r]
        relc = os.path.join(SETS[r], "cliffs/12M_cliff_pressure_fraction_NG_RG_caption.txt")
        for ga in g["gauges"]:
            pat = r"%s: slopes >= 60 deg give (\d+-\d+) %% of the surface pressure before the late window and (\d+-\d+) %% in it" % ga
            row.append(safe(capx, relc, pat, r, "%s >= 60 deg share before the late window (range)" % ga, 1) + " %")
            row.append(safe(capx, relc, pat, r, "%s >= 60 deg share in the late window (range)" % ga, 2) + " %")

            def lm():
                v = sum(fxv(r, ["pressure_share_by_bin", ga, "late_mean"])[3:])
                return N.add(r, "%s mean share >= 60 deg in the late window" % ga, 100.0 * v,
                             os.path.join(SETS[r], "checks/facts.json"),
                             "facts.json: sum pressure_share_by_bin/%s/late_mean[3:]" % ga, pct0)
            row.append(safe(lm, what="facts pressure_share_by_bin"))
            row.append(safe(fx, r, ["pressure_share_by_bin", ga, "last", 6], "90-180 deg share of %s at the last obs." % ga,
                            pct0, 100.0))
        rows.append(row)
    return hdr, rows, [1.2] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def cliff_pf_compare(g):
    out = []
    for ga in g["gauges"]:
        lm = {}
        for r in have(g):
            try:
                v = sum(fxv(r, ["pressure_share_by_bin", ga, "late_mean"])[3:])
                lm[r] = (v, N.add(r, "%s mean share >= 60 deg in the late window" % ga, 100.0 * v,
                                  os.path.join(SETS[r], "checks/facts.json"),
                                  "facts.json: sum pressure_share_by_bin/%s/late_mean[3:]" % ga, pct0))
            except (KeyError, OSError):
                pass
        out.append("%s, mean share of slopes >= 60 deg in the late window: %s." % (ga, order_text(lm, reverse=True)))
    return " ".join(out)


# regions
LOBES = ["small lobe (head, with the touchdown point)", "big lobe (body)"]
SHELLS = ["distance 0–30 m", "distance 30–100 m", "distance 100–300 m", "distance 1000–3000 m"]


def reg_stat(r, kind, ga, region, phase, idx, quantity, fmt=pct0):
    return manx(os.path.join(SETS[r], "regions/regions_figures_manifest.json"), ["stats", kind, ga, region, phase, idx], r,
                quantity, fmt, 100.0)


def reg_about(g):
    t = tload(first_caption(g))
    return (t.split(", along the descent")[0] + ", along the descent (stacked, against the gauge altitude, late window "
            "boxed; NG gap hatched in the NG panel). The region sums equal the model surface pressure; the NG grain term "
            "is not included. New law q_min = 0.02 q_max.")


def reg_numbers(g):
    hdr, rows = ["run"], []
    for ga in g["gauges"]:
        if g["fam"] == "reg_lobe":
            hdr += ["%s small lobe median before / in late window" % ga, "%s small lobe range in late window" % ga]
        else:
            hdr += ["%s median in late window: 0-30 / 30-100 / 100-300 m" % ga, "%s 1000-3000 m median before" % ga]
    for r in RK:
        if not g["paths"].get(r):
            rows.append([r, "no figure"] + [""] * (len(hdr) - 2))
            continue
        row = [r]
        for ga in g["gauges"]:
            if g["fam"] == "reg_lobe":
                L = LOBES[0]
                row.append(safe(lambda: "%s / %s" % (
                    reg_stat(r, "lobe", ga, L, "before_late_window", 2, "%s small lobe median before" % ga),
                    reg_stat(r, "lobe", ga, L, "late_window", 2, "%s small lobe median late" % ga)), what="regions"))
                row.append(safe(lambda: "%s-%s" % (
                    reg_stat(r, "lobe", ga, L, "late_window", 0, "%s small lobe min late" % ga, "%.0f"),
                    reg_stat(r, "lobe", ga, L, "late_window", 1, "%s small lobe max late" % ga)), what="regions"))
            else:
                row.append(safe(lambda: " / ".join(reg_stat(r, "distance", ga, s, "late_window", 2,
                                                            "%s %s median late" % (ga, s)) for s in SHELLS[:3]),
                                what="regions"))
                row.append(safe(lambda: reg_stat(r, "distance", ga, SHELLS[3], "before_late_window", 2,
                                                 "%s %s median before" % (ga, SHELLS[3])), what="regions"))
        rows.append(row)
    return hdr, rows, [1.2] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def reg_compare(g):
    out = []
    for ga in g["gauges"]:
        v = {}
        for r in have(g):
            try:
                st = jload(os.path.join(SETS[r], "regions/regions_figures_manifest.json"))["stats"]
                if g["fam"] == "reg_lobe":
                    v[r] = (st["lobe"][ga][LOBES[0]]["late_window"][2],
                            reg_stat(r, "lobe", ga, LOBES[0], "late_window", 2, "%s small lobe median late" % ga))
                else:
                    v[r] = (st["distance"][ga][SHELLS[1]]["late_window"][2],
                            reg_stat(r, "distance", ga, SHELLS[1], "late_window", 2, "%s %s median late" % (ga, SHELLS[1])))
            except (KeyError, OSError):
                pass
        if g["fam"] == "reg_lobe":
            out.append("%s: small-lobe median share in the late window %s." % (ga, order_text(v)))
        else:
            out.append("%s: median share from 30-100 m around the touchdown point in the late window %s." %
                       (ga, order_text(v, reverse=True)))
    return " ".join(out)


def vis_about(g):
    return cap_sentences(g, 1)


def generic_about(g):
    return cap_sentences(g, 2) or "(no caption file)"


FAMILY = {
    "ov": dict(about=ov_about, numbers=ov_numbers, compare=ov_compare),
    "ov_conv": dict(about=ovc_about, numbers=ovc_numbers, compare=ovc_compare),
    "vo": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare),
    "vo_ng": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare),
    "vo_rg": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare),
    "vo_ngrg": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare),
    "vr": dict(about=vr_about, numbers=vr_numbers, compare=vr_compare),
    "cv": dict(about=cv_about, numbers=cv_numbers, compare=cv_compare),
    "insp": dict(about=insp_about, numbers=insp_numbers, compare=insp_compare),
    "md_pub": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_alt": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_ng": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_rg": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_ngrg": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_time": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "md_press": dict(about=fmd_about, numbers=fmd_numbers, compare=fmd_compare),
    "map_tpl": dict(about=map_about, numbers=map_numbers, compare=map_compare),
    "map_latlon": dict(about=map_about, numbers=map_numbers, compare=map_compare),
    "map_3d": dict(about=map_about, numbers=map_numbers, compare=map_compare),
    "map_panel": dict(about=map_about, numbers=map_numbers, compare=map_compare),
    "geom_panel": dict(about=geom_about),
    "sa_union": dict(about=sa_about, numbers=sa_numbers, compare=sa_compare),
    "sa_occ": dict(about=sa_about, numbers=sa_numbers, compare=sa_compare),
    "cd": dict(about=cd_about, numbers=cd_numbers, compare=cd_compare),
    "cliff_bins": dict(about=cliff_bins_about, numbers=cliff_bins_numbers, compare=cliff_bins_compare),
    "cliff_pf": dict(about=cliff_pf_about, numbers=cliff_pf_numbers, compare=cliff_pf_compare),
    "reg_lobe": dict(about=reg_about, numbers=reg_numbers, compare=reg_compare),
    "reg_distance": dict(about=reg_about, numbers=reg_numbers, compare=reg_compare),
    "vis": dict(about=vis_about),
    "generic": dict(about=generic_about),
}


# ------------------------------------------------------------------------------------------------ setup / check
def amps_head():
    try:
        head = open("/data001/highresolu_vt/AMPS/.git/HEAD").read().strip()
        if head.startswith("ref: "):
            ref = head[5:]
            p = "/data001/highresolu_vt/AMPS/.git/" + ref
            if os.path.exists(p):
                return open(p).read().strip() + " (%s, read from .git; no git command run)" % ref
        return head
    except OSError as e:
        return "unknown (%s)" % e


LINE_FAMS = ("ov", "ov_conv", "vo", "vo_ng", "vo_rg", "vo_ngrg", "vr", "cv", "md_pub", "md_alt", "md_ng", "md_rg",
             "md_ngrg", "md_time", "md_press", "cliff_bins", "cliff_pf", "reg_lobe", "reg_distance")


def setup(args):
    SETS.clear()
    SETS.update(discover_sets(getattr(args, "sets_spec", None)))
    STEMS.clear()
    for r in RK:
        STEMS[r] = find_stem(SETS[r]) if SETS.get(r) else None
    groups = catalogue()
    for gi, g in enumerate(groups):
        g["id"] = gi
        g["sha"] = {who: sha256(p) for who, p in g["paths"].items() if p}
        avail = [who for who, p in g["paths"].items() if p]
        g["identical"] = (g["section"] != "overview" and len(avail) == len(RK) and len(set(g["sha"].values())) == 1)
        g["show"] = ["all"] if g["section"] == "overview" else ([RK[0]] if g["identical"] else list(RK))
    return groups


def cmd_check(args):
    groups = setup(args)
    print("full sets:")
    for r in RK:
        print("  %-15s %s  stem %s  %d PNGs" % (r, SETS.get(r), STEMS.get(r), len(fs_list(r))))
    fam_missing = [g for g in groups if g["section"] in FS_SECTIONS and g["fam"] == "generic"]
    print("clusters: %d (overview %d, mv %d, full-set %d, inspection %d); identical %d; unknown plot types %d" % (
        len(groups), sum(g["section"] == "overview" for g in groups), sum(g["section"] == "mv" for g in groups),
        sum(g["section"] in FS_SECTIONS for g in groups), sum(g["section"] == "inspection" for g in groups),
        sum(g["identical"] for g in groups), len(fam_missing)))
    for g in fam_missing:
        print("  unknown type (generic summary):", g["cid"])
    obs = json.load(open(args.observations)) if args.observations else {}
    no_obs = sorted(set(g["cid"] for g in groups if not g["identical"] and g["cid"] not in obs))
    print("clusters without an observation entry: %d" % len(no_obs))
    for c in no_obs:
        print("  ", c)
    for g in groups:
        miss = [w for w in g["show"] if w != "all" and not g["paths"].get(w)]
        if miss:
            print("  missing in %s: %s" % (g["cid"], ", ".join(miss)))
    # dry run of every extractor (numbers only, no document)
    for g in groups:
        N.ctx = g["cid"]
        F = FAMILY[g["fam"]]
        if g["identical"]:
            continue
        try:
            F["about"](g)
            if F.get("numbers"):
                F["numbers"](g)
            if F.get("compare"):
                F["compare"](g)
        except Exception as e:  # report and continue: check mode only
            print("  EXTRACTOR ERROR %s (%s): %s: %s" % (g["cid"], g["fam"], type(e).__name__, e))
    print("missing inputs (n/a): %d" % len(MISSING_INPUTS))
    for m in MISSING_INPUTS[:40]:
        print("  ", m)


def cmd_build(args):
    out = os.path.abspath(args.out)
    if os.path.exists(os.path.join(out, DOC_NAME)):
        sys.exit("refusing to overwrite %s (move the folder aside first)" % os.path.join(out, DOC_NAME))
    os.makedirs(out, exist_ok=True)
    args.started = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    import docx
    args.docx_version = getattr(docx, "__version__", "?")
    args.pil_version = PIL.__version__
    OBS.clear()
    OBS.update({k: v for k, v in json.load(open(args.observations)).items() if not k.startswith("_")})
    groups = setup(args)
    jobs, img = [], {}
    for g in groups:
        sizes = {}
        for who in g["show"]:
            if g["paths"].get(who):
                with Image.open(g["paths"][who]) as im:
                    sizes[who] = im.size
        g["native"] = sizes
        g["layout"] = layout([h / w for (w, h) in sizes.values()], len(g["show"]))
        kind = "png" if g["fam"] in LINE_FAMS else "jpg"
        for who in sizes:
            target = int(min(sizes[who][0], max(args.min_px, round(g["layout"]["width_in"] * 300)), 1400))
            base = re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.basename(g["paths"][who])[:-4])
            dst = os.path.join(out, "images", g["section"], "%03d_%s__%s.%s" % (g["id"], base, re.sub(r"[^A-Za-z0-9]", "", who),
                                                                                kind))
            jobs.append((g["paths"][who], dst, kind, target, args.quality))
            img[(g["id"], who)] = dict(src=g["paths"][who], dst=dst)
    with ProcessPoolExecutor(args.jobs) as ex:
        res = list(ex.map(prep_image, jobs))
    tot = 0
    bydst = {}
    for (src, dst, kind, target, q), (d, size, nbytes) in zip(jobs, res):
        tot += nbytes
        bydst[d] = (size, nbytes)
    for v in img.values():
        v["px"], v["bytes"] = bydst[v["dst"]]
    args.img_note = "%d images, %.1f MB, JPEG for maps and 256-colour PNG for line plots, %d-1400 px wide" % (
        len(jobs), tot / 1e6, args.min_px)
    build_doc(out, groups, img, args)
    with open(os.path.join(out, "numbers_used.json"), "w") as f:
        json.dump(dict(created_utc=args.started, generator=os.path.abspath(__file__), n=len(N.items), numbers=N.items),
                  f, indent=1, ensure_ascii=False)
    with open(os.path.join(out, "missing_inputs.json"), "w") as f:
        json.dump(MISSING_INPUTS, f, indent=1)
    man = dict(created_utc=args.started,
               finished_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S"),
               command=sys.argv, python=sys.executable, python_docx=args.docx_version, pillow=args.pil_version,
               slurm_job=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename, amps_head=amps_head(),
               generator_sha256=sha256(os.path.abspath(__file__)), observations=args.observations,
               observations_sha256=sha256(args.observations), sets=SETS, stems=STEMS,
               document=os.path.join(out, DOC_NAME), document_bytes=os.path.getsize(os.path.join(out, DOC_NAME)),
               images_bytes=tot,
               groups=[dict(fig=g["fig"], id=g["id"], cid=g["cid"], section=g["section"], title=g["title"], fam=g["fam"],
                            law=g["law"], identical=g["identical"], show=g["show"], layout=g["layout"], native=g["native"],
                            paths=g["paths"], sha256=g["sha"], missing=[w for w, p in g["paths"].items() if not p],
                            images={who: dict(dst=os.path.relpath(img[(g["id"], who)]["dst"], out),
                                              px=img[(g["id"], who)]["px"], bytes=img[(g["id"], who)]["bytes"])
                                    for who in g["show"] if g["paths"].get(who)}) for g in groups])
    with open(os.path.join(out, "build_manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    print("document %s: %.1f MB; %d figures (%d identical in all runs, %d missing-run placeholders); %d images (%.1f MB); "
          "%d numbers; %d missing inputs (n/a)" % (
              os.path.join(out, DOC_NAME), man["document_bytes"] / 1e6, len(groups), sum(g["identical"] for g in groups),
              sum(1 for g in groups for w in g["show"] if w != "all" and not g["paths"].get(w)), len(jobs), tot / 1e6,
              len(N.items), len(MISSING_INPUTS)))


def cmd_validate(args):
    import zipfile
    from docx import Document
    from docx.oxml.ns import qn
    out = os.path.abspath(args.out)
    path = os.path.join(out, DOC_NAME)
    man = json.load(open(os.path.join(out, "build_manifest.json")))
    doc = Document(path)
    heads, caps, placeholders = {}, [], 0
    for p in doc.paragraphs:
        sn = p.style.name if p.style is not None else ""
        if sn.startswith("Heading") or sn == "Title":
            heads[sn] = heads.get(sn, 0) + 1
        if sn == "Caption" and re.match(r"^Figure \d+: ", p.text):
            caps.append(p.text)
    for t in doc.tables:
        for row in t.rows:
            for c in row.cells:
                if c.text.strip().startswith(tuple(RK)) and "not produced for this run" in c.text:
                    placeholders += 1
    shapes = doc.inline_shapes
    pics = [(s.width / 914400, s.height / 914400) for s in shapes]
    problems = []
    k = exp_imgs = exp_place = 0
    for g in man["groups"]:
        n = sum(1 for w in g["show"] if g["paths"].get(w) or w == "all")
        exp_imgs += n
        exp_place += sum(1 for w in g["show"] if w != "all" and not g["paths"].get(w))
        W, H = usable(g["layout"]["orient"])
        cols = g["layout"]["cols"]
        sizes = pics[k:k + n]
        k += n
        if not sizes:
            continue
        wmax = max(w for w, h in sizes)
        hrow = max(h for w, h in sizes)
        rows = math.ceil(len(g["show"]) / cols)
        if wmax > W / cols - 0.05:
            problems.append("Figure %d: image width %.2f in > cell %.2f in" % (g["fig"], wmax, W / cols))
        if rows * (hrow + LABEL + 0.05) + OVERHEAD - 0.4 > H:
            problems.append("Figure %d: grid height > page" % g["fig"])
    body = doc.element.body
    cur, nt = [], 0
    for el in body.iterchildren():
        if el.tag == qn("w:tbl"):
            cur.append(el)
        sp = el.find(".//" + qn("w:sectPr")) if el.tag == qn("w:p") else (el if el.tag == qn("w:sectPr") else None)
        if sp is not None:
            pg, mar = sp.find(qn("w:pgSz")), sp.find(qn("w:pgMar"))
            uw = (int(pg.get(qn("w:w"))) - int(mar.get(qn("w:left"))) - int(mar.get(qn("w:right")))) / 1440.0
            for tb in cur:
                grid = tb.find(qn("w:tblGrid"))
                tw = sum(int(gc.get(qn("w:w"))) for gc in grid.findall(qn("w:gridCol"))) / 1440.0
                nt += 1
                if tw > uw + 0.02:
                    problems.append("table width %.2f in > usable %.2f in" % (tw, uw))
            cur = []
    with zipfile.ZipFile(path) as zf:
        fe = [n for n in zf.namelist() if n.endswith("/")]
        media = len([n for n in zf.namelist() if n.startswith("word/media/")])
    if fe:
        problems.append("zip folder entries: %s" % fe[:5])
    nums = json.load(open(os.path.join(out, "numbers_used.json")))
    numbered = [int(re.match(r"^Figure (\d+):", c).group(1)) for c in caps] == list(range(1, len(caps) + 1))
    law_ok = all(("q_min = 0.02 q_max" in c) or ("independent of the H2O law" in c) for c in caps)
    miss_in = json.load(open(os.path.join(out, "missing_inputs.json")))
    res = dict(document=path, bytes=os.path.getsize(path), mb=round(os.path.getsize(path) / 1e6, 2), headings=heads,
               figure_captions=len(caps), figure_groups=len(man["groups"]), captions_numbered_1_to_n=numbered,
               captions_state_the_law=law_ok, identical_groups=sum(g["identical"] for g in man["groups"]),
               inline_images=len(shapes), expected_images=exp_imgs, media_files=media, placeholders=placeholders,
               expected_placeholders=exp_place, sections=len(doc.sections), tables=len(doc.tables),
               tables_width_checked=nt, numbers_recorded=nums["n"], missing_inputs=len(miss_in), problems=problems,
               libreoffice="not available (no soffice on login or compute nodes; not on conda-forge): no PDF render")
    res["passed"] = (len(caps) == len(man["groups"]) and numbered and law_ok and len(shapes) == exp_imgs and
                     placeholders == exp_place and not problems)
    os.makedirs(os.path.join(out, "validation"), exist_ok=True)
    with open(os.path.join(out, "validation", "validation.json"), "w") as f:
        json.dump(res, f, indent=1)
    from PIL import ImageDraw, ImageFont
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf", 11)
        font2 = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf", 10)
    except OSError:
        font = font2 = ImageFont.load_default()
    dpi = 100
    sel = args.preview or sorted(set([1] + [g["fig"] for g in man["groups"]][::8] + [man["groups"][-1]["fig"]]))
    for g in [g for g in man["groups"] if g["fig"] in sel]:
        o = g["layout"]["orient"]
        pw, ph = PAGE[o]
        im = Image.new("RGB", (int(pw * dpi), int(ph * dpi)), "white")
        dr = ImageDraw.Draw(im)
        dr.rectangle([MARGIN * dpi, MARGIN * dpi, (pw - MARGIN) * dpi, (ph - MARGIN) * dpi], outline=(200, 200, 200))
        dr.text((MARGIN * dpi, MARGIN * dpi), "x.y %s" % g["title"][:110], fill=(31, 58, 95), font=font)
        W, H = usable(o)
        cols = g["layout"]["cols"]
        w_in = g["layout"]["width_in"]
        hs = [w_in * v["px"][1] / v["px"][0] for v in g["images"].values()]
        hrow = max(hs) if hs else 1.0
        y0 = MARGIN + 0.5
        for k2, who in enumerate(g["show"]):
            r, col = divmod(k2, cols)
            x = MARGIN + col * (W / cols) + (W / cols - w_in) / 2
            y = y0 + r * (hrow + LABEL + 0.05)
            dr.text((x * dpi + 2, y * dpi), "all runs" if who == "all" else who, fill=(160, 16, 16), font=font)
            if who in g["images"]:
                pic = Image.open(os.path.join(out, g["images"][who]["dst"])).convert("RGB")
                pic = pic.resize((int(w_in * dpi), int(w_in * dpi * pic.size[1] / pic.size[0])), Image.LANCZOS)
                im.paste(pic, (int(x * dpi), int((y + LABEL) * dpi)))
            else:
                dr.text((x * dpi + 2, (y + LABEL + 0.2) * dpi), "not produced for this run", fill=(90, 90, 90), font=font2)
        yc = y0 + math.ceil(len(g["show"]) / cols) * (hrow + LABEL + 0.05) + 0.05
        dr.text((MARGIN * dpi, yc * dpi), "Figure %d: %s" % (g["fig"], g["title"][:110]), fill=(0, 0, 0), font=font)
        im.save(os.path.join(out, "validation", "layout_preview_fig%03d.png" % g["fig"]))
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    c = sp.add_parser("check")
    c.add_argument("--sets-spec")
    c.add_argument("--observations")
    b = sp.add_parser("build")
    b.add_argument("--out", required=True)
    b.add_argument("--observations", required=True)
    b.add_argument("--sets-spec")
    b.add_argument("--jobs", type=int, default=8)
    b.add_argument("--min-px", type=int, default=800)
    b.add_argument("--quality", type=int, default=75)
    v = sp.add_parser("validate")
    v.add_argument("--out", required=True)
    v.add_argument("--preview", type=int, nargs="*")
    args = ap.parse_args()
    {"check": cmd_check, "build": cmd_build, "validate": cmd_validate}[args.cmd](args)


if __name__ == "__main__":
    main()
