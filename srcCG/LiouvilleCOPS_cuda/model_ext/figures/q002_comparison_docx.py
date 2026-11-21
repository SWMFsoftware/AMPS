#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-19, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_comparison_docx.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_comparison_docx.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: d6e98acfeeb6e5c2b11d37f88e7eb5e89e4b1b175bd87173c0861fd57f863a93
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Word document clustering the same plot of the ten new-law (H2O q_min 0.02 q_max) stage-2 runs, with summaries (2026-11-19)
# RUN: separate venv cops_output/tools/docx_venv (python-docx, Pillow); build + validate in a Slurm job (HOW-TO-RERUN.md of figures_cluster/12M/stage2_q002_comparison/)
# ==== end of header ====
"""Word document clustering the same plot of every run of the NEW-LAW stage-2 round (H2O q_min = 0.02 q_max), 2026-11-19.

User request 2026-11-19 (relayed by the lead): "similarly to what you have done priorly, generate word doc clastering the
same plot from each of the runs and summarizing them". Derived from fullsets_comparison_docx.py (2026-11-17; same rules:
every printed number read from a file and listed in numbers_used.json; hand-written sentences checked against the source
files; observations from the images labelled as such).
Runs (fixed order): the refitted D2 (D2_q002_warm) and the nine regularized stage-2 chains on base_q002: softplus, A, B, D,
C f1 / P_m 1e-4, C f10 / P_m 1e-4, C f1 / P_m 0.01, C f10 / P_m 0.01, C f10 / P_m 0.01 templates off.
Inputs: figures_cluster/12M/source_rate_q002/<run>/ (activity/, production_rate/, inspection/, checks/, README.txt) and
figures_cluster/12M/model_vs_observations_nonneg_options_stage2_q002/ (per-run new-vs-old-law figures, residuals,
convergence, the multi-run figures, summary.csv / .txt, tables/, *_residuals_check.json).

    q002_comparison_docx.py build    --out <new folder> [--jobs 8]
    q002_comparison_docx.py validate --out <folder>
Interpreter: the separate venv $COPS_OUTPUT_ROOT/tools/docx_venv (python-docx, Pillow), not the project .venv.
Refuses to overwrite an existing document; refuses to write it if a literal number of a hand-written sentence is not
found in a source file (numbers_unverified.json).
"""
import argparse
import csv
import datetime
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
NUMTOK = re.compile(r"(?<![A-Za-z_\d.])\d+(?:[.,]\d+)*(?:e[+-]?\d+)?")
DOC_NAME = "12M_stage2_q002_comparison.docx"
DOC_DATE = "2026-11-19"
NEW_LAW = "new law q_min = 0.02 q_max"
OLD_LAW = "old law q_min = 0.2 q_max"

# fixed order requested by the lead
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
SECTIONS = [("overview", "Overview: the multi-run figures (new law)"),
            ("model_data", "Model vs COPS per run: new law vs old law"),
            ("activity", "H2O activity factor A per run (new law)"),
            ("production_rate", "H2O production rate per run (new law)"),
            ("inspection", "Inspection: new-law vs old-law production rate per run")]
WHY_MISSING = {
    "inspection": ("no old-law maps exist for this chain, so no inspection figure",
                   SR + "/C_f10_pm1em4_q002_reg/README.txt"),
    "residuals": ("the comparison set has new-vs-old residual figures only for softplus, A, B, D, D2 and the best C "
                  "(one figure per C setting would repeat the model/data figures)", MV + "/README.txt"),
    "convergence": ("one fit, no cycles, so no convergence figure (a line in 12M_q002_convergence_all)",
                    MV + "/README.txt"),
}


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
    """hand-written text: every numeric token must occur in one of the source files (recorded in numbers_used.json);
    tokens already registered in this context are skipped unless force; fallback = (locator, source)"""
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


# summary.csv of the comparison set ------------------------------------------------------------------------------------
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
    """rms ln(model/obs): summary.csv; if empty there (D2), the run's residual-check json (recomputed from the tables)"""
    col = "rms_%s_%s" % (gauge, "late" if late else "all")
    s = sx(run, col, "rms ln %s (%s)" % (gauge, "late window" if late else "all fitted rows"), "%.4f", old=old)
    if s is not None:
        return s
    for p in [resid_json_for(run), MV + "/12M_q002_best_residuals_check.json"]:
        if not p:
            continue
        tab = (row_key(run, old)) + ".dat"
        for k, v in jload(p)["fits"].items():
            if os.path.basename(v["table"]) == tab:
                return N.add(run, "rms ln %s (%s)%s" % (gauge, "late window" if late else "all fitted rows",
                                                        " (old law)" if old else ""),
                             v[gauge]["rms_late" if late else "rms"], p,
                             "fits['%s'].%s.%s" % (k, gauge, "rms_late" if late else "rms"), "%.4f")
    return "n/a"


# the plotted observation rows (comparison set tables/) -------------------------------------------------------------------
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
            ix = {n: i for i, n in enumerate(cols)}
            alt = float(r[ix["altitude"]])
            if r[ix["useNG"]] == "1":
                ng.append((alt, math.log(float(r[ix["pNG"]]) / float(r[ix["refNG"]]))))
            if r[ix["useRG"]] == "1":
                rg.append((alt, math.log(float(r[ix["pRG"]]) / float(r[ix["refRG"]]))))
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


# source-rate sets -----------------------------------------------------------------------------------------------------
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
        return N.add(run, quantity, "n/a", path, "regex (no match): " + pattern, "%s", note="not in caption")
    return N.add(run, quantity, conv(m.group(group)), path, "regex group %d: %s" % (group, pattern), fmt)


def manx(path, keys, run, quantity, fmt="%s", scale=1.0):
    v = nested(jload(path), keys)
    if scale != 1.0:
        v = v * scale
    return N.add(run, quantity, v, path, os.path.basename(path) + ":" + "/".join(str(k) for k in keys), fmt)


def order_text(vals, reverse=False):
    items = sorted(vals.items(), key=lambda kv: kv[1][0], reverse=reverse)
    return ", ".join("%s %s" % (k, t) for k, (v, t) in items)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ catalogue
def catalogue():
    G = []

    def add(sec, title, fam, paths, law, **kw):
        G.append(dict(section=sec, title=title, fam=fam, paths=paths, law=law, **kw))

    # overview: multi-run figures, one image each
    for setn, who in (("best", "best of each option (softplus, D2, A, B, best C, D)"),
                      ("C", "the five option-C chains and the softplus reference")):
        for suf, what in (("time", "against UTC"), ("altitude", "against altitude: (a) NG, (b) RG, (c) NG and RG"),
                          ("altitude_NG", "NG against altitude"), ("altitude_RG", "RG against altitude"),
                          ("altitude_NG_RG", "NG and RG together against altitude"),
                          ("residuals", "log residuals ln(model/obs) against altitude")):
            add("overview", "New law, %s: %s" % (who, what), "ov", {"all": MV + "/12M_q002_%s_%s.png" % (setn, suf)},
                NEW_LAW, set=setn)
    add("overview", "New law: data cost per cycle of the nine stage-2 chains and the refitted D2", "ov_conv",
        {"all": MV + "/12M_q002_convergence_all.png"}, NEW_LAW)
    # per run: model / data, new vs old
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
        add("model_data", "%s, new vs old law, per run" % what, fam, paths,
            "%s (solid) vs %s (dashed)" % (NEW_LAW, OLD_LAW), suf=suf,
            why=None if fam != "vr" else "residuals")
    paths = {}
    for r in RK:
        p = MV + "/12M_q002_convergence_%s.png" % RUN[r]["mv"]
        paths[r] = p if os.path.exists(p) else None
    add("model_data", "Data cost per cycle (fit and fine tuning), new vs old law, per chain", "cv", paths,
        "%s (solid black) vs %s (dashed grey)" % (NEW_LAW, OLD_LAW), why="convergence")
    # maps
    for sec, q, tpl, times in (("activity", "activity", "12M_templates_activity", [None]),
                               ("production_rate", "h2o", "12M_templates_h2o_production_rate_%sUTC",
                                ["103310", "103926"])):
        for t in times:
            base = tpl % t if t else tpl
            gt = "12M_groundtrack_activity" if not t else "12M_groundtrack_h2o_production_rate_%sUTC" % t
            what = "H2O activity factor A" if not t else "H2O production rate F_H2O at %s:%s:%s UTC" % (t[:2], t[2:4], t[4:])
            pq = "activity" if not t else "h2o_production_rate"
            items = [(base + ".png", what + ": maps with effective slope and cliff templates (six panels)", "map_tpl",
                      [("a", "latlon_global"), ("b", "latlon_zoom")]),
                     (gt + "_latlon.png", what + " with the ground track, longitude / latitude", "map_latlon",
                      [("a", "latlon_global"), ("b", "latlon_zoom")]),
                     (gt + "_3d.png", what + " with the ground track, 3-D views", "map_3d",
                      [("a", "3d_a"), ("b", "3d_b"), ("c", "3d_c")]),
                     ("panels/" + base + "_a_%s_global.png" % pq, what + ": panel (a) global map alone", "map_panel",
                      [("a", "latlon_global")]),
                     ("panels/" + base + "_b_%s_impact_area.png" % pq, what + ": panel (b) impact area alone",
                      "map_panel", [("b", "latlon_zoom")])]
            for p, nm in (("c_slope_global", "panel (c) effective gravitational slope, global"),
                          ("d_slope_impact_area", "panel (d) effective gravitational slope, impact area"),
                          ("e_cliff_global", "panel (e) cliff templates, global"),
                          ("f_cliff_impact_area", "panel (f) cliff templates, impact area")):
                items.append(("panels/" + base + "_" + p + ".png", what + ": " + nm, "geom_panel", []))
            for rel, title, fam, panels in items:
                paths = {r: (srp(r, sec + "/" + rel) if os.path.exists(srp(r, sec + "/" + rel)) else None) for r in RK}
                add(sec, title, fam, paths, NEW_LAW if fam != "geom_panel" else
                    "geometry, independent of the H2O law (panel of the new-law figure set)", q=q, t=t, panels=panels,
                    rel=sec + "/" + rel)
    for t in ("103310", "103926"):
        rel = "inspection/12M_inspect_h2o_production_rate_new_vs_old_law_%sUTC.png" % t
        paths = {r: (srp(r, rel) if os.path.exists(srp(r, rel)) else None) for r in RK}
        add("inspection", "Inspection at %s:%s:%s UTC: new-law vs old-law F_H2O on common classes, and the new q/q_max"
            % (t[:2], t[2:4], t[4:]), "insp", paths,
            "%s (a, b) vs %s (c, d), common colour classes; (e, f) the new law's q/q_max" % (NEW_LAW, OLD_LAW),
            t=t, rel=rel, why="inspection")
    return G


def caption_of(png):
    return png[:-4] + "_caption.txt" if png else None


# ------------------------------------------------------------------------------------------------ images and layout
PAGE = {"P": (8.5, 11.0), "L": (11.0, 8.5)}
MARGIN = 0.6
OVERHEAD = 1.55    # heading, files line, caption (inches)
LABEL = 0.24       # run label above each image
PAD = 0.12         # cell padding


def usable(o):
    w, h = PAGE[o]
    return w - 2 * MARGIN, h - 2 * MARGIN


def layout(aspects, n):
    """best (orientation, cols, rows, image width in inches) for n images of the given aspect ratios (h/w), one page"""
    rmax = max(aspects)
    best = None
    if n == 1:
        opts = [("P", 1, 1), ("L", 1, 1)]
    else:
        opts = [(o, c, math.ceil(n / c)) for o in ("P", "L") for c in (2, 3, 4, 5) if c <= n]
    for o, cols, rows in opts:
        W, H = usable(o)
        wcell = W / cols
        hrow = (H - OVERHEAD - (0.3 if n == 1 else 0.0)) / rows - LABEL - 0.05
        w = min(wcell - PAD, hrow / rmax)
        if n == 1:
            w = min(w, 7.0)
        score = w * w * rmax
        if best is None or score > best[4] + 1e-9:
            best = (o, cols, rows, w, score)
    return dict(orient=best[0], cols=best[1], rows=best[2], width_in=round(best[3], 3))


def prep_image(job):
    src, dst, kind, target = job
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
        im.save(tmp, "JPEG", quality=78, optimize=True)
    os.replace(tmp, dst)
    return dst, im.size, os.path.getsize(dst)


# ------------------------------------------------------------------------------------------------ sources for text
def gsrc(g):
    caps = [caption_of(p) for p in g["paths"].values() if p]
    caps = [c for c in caps if c and os.path.exists(c)]
    readmes = [srp(r, "README.txt") for r in RK]
    return caps + [MV + "/README.txt", MV + "/summary.txt"] + readmes + [FINDINGS, OUTVIS]


def fig_fallback(g):
    first = [p for p in g["paths"].values() if p][0]
    return ("observation from the figure: approximate axis reading or a bin / category label shown in the figure",
            first + " and the other runs' PNG of this figure")


FRONT = None


def front_sources():
    return ([MV + "/README.txt", MV + "/summary.txt", MV + "/CAPTIONS.md"] + [srp(r, "README.txt") for r in RK] +
            [FINDINGS, OUTVIS])


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
        r = p.add_run("12M stage 2, new H2O law (q_min = 0.02 q_max): the ten runs plot by plot, %s   |   page " % DOC_DATE)
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

    # ---------------- title, purpose, law
    N.ctx = "1. title / purpose"
    SRC["src"], SRC["fb"] = FR, ("document date / build record (generator constants, build_manifest.json)",
                                 os.path.abspath(__file__))
    doc.add_paragraph("12M stage 2 with the corrected H2O law: the ten runs clustered plot by plot", style="Title")
    para(DOC_DATE + " (US Central). SHAP8 12M shape model, COPS NG / RG data of Rosetta's final descent (2016-09-30), "
         "collisionless (Liouville) forward model.", italic=True)
    para("This document puts the same figure of every run of the last fitting round side by side: the refitted D2 "
         "(D2_q002_warm) and the nine regularized stage-2 chains refitted with the corrected H2O night-side law. For "
         "each plot type the runs appear in one labelled grid, always in the order D2, softplus, A, B, D, C f1/1e-4, "
         "C f10/1e-4, C f1/0.01, C f10/0.01, C f10/0.01 off, followed by a summary: what the plot shows (from its "
         "caption), a table of the numbers of every run (read from the captions, summary.csv, the inspection files, the "
         "manifests and the plotted tables; all listed with their source in numbers_used.json), how the runs compare, "
         "and what is visible in the images (marked as observations). It adds no new fits or figures and no physical "
         "interpretation beyond the captions, READMEs and FINDINGS.md.", bold_lead="Purpose.")
    para("Every fit, map and figure before 2026-11-18 used the H2O production q_H2O = q_min on night and cast-shadow "
         "facets and q_min + (q_max - q_min) cos SZA on day facets with q_min = 0.2 q_max (not intended); srcCG has "
         "0.02. User decision 2026-11-18: use 0.02 (FINDINGS A18). All runs here were refitted on base_q002 with the "
         "NEW law q_min = 0.02 q_max; the per-run model/data, residual, convergence and inspection figures compare them "
         "with their OLD-law (q_min = 0.2 q_max) counterparts. Every figure caption below states the law of what it "
         "shows.", bold_lead="The law correction.")
    para("Run labels: C f1/1e-4 = option C, smoothness f 1, multiplier penalty P_m 1e-4 (f10: f 10; 0.01: P_m 0.01; "
         "off: templates off). A run that lacks a plot keeps its place in the grid with the reason. Plot types whose "
         "PNGs are byte-identical in every run (SHA256) are shown once. Colour classes of the activity and "
         "production-rate maps are chosen per panel and per run, so equal colours do not mean equal values across runs "
         "(the inspection figures compare new and old law of one run on common classes). Not yet produced for the "
         "new-law runs: the full-set categories contribution density, source areas, cliffs and regions; this document "
         "therefore has no clusters for them. Images: downscaled copies (%s); originals (300 dpi PNG and EPS) listed "
         "with SHA-256 in Appendix A." % args.img_note, bold_lead="How to read it.")

    # ---------------- cases table
    N.ctx = "2. cases table"
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

    # ---------------- summary (A18)
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
    ss = {}
    dk = {}
    for r in RK:
        v = ixv(r, "103310", ["new", "seen_share_of_Q"])
        if v is not None:
            ss[r] = (v, ix(r, "103310", ["new", "seen_share_of_Q"], "seen share of Q_H2O at 10:33:10", pct0, 100.0))
            dk[r] = (ixv(r, "103310", ["new", "dark_share_of_Q"]),
                     ix(r, "103310", ["new", "dark_share_of_Q"], "night + cast-shadow share of Q_H2O at 10:33:10", pct0,
                        100.0))
    para("Q_H2O at 10:33:10 UTC spans %s s^-1 (%s), while the production of the facets a fitted observation sees is "
         "%s s^-1. Seen share of Q_H2O (inspection files; none for C f10/1e-4): %s. Night + cast-shadow share of "
         "Q_H2O: %s. Only option C with P_m 0.01 has its total mostly seen; in the other runs most of Q_H2O comes from "
         "terrain no fitted observation sees, so their totals are not set by the data (A15, A18)." %
         ("%s ... %s" % (min(q.values())[1], max(q.values())[1]), order_text(q), order_text(qs),
          order_text(ss, reverse=True), order_text(dk, reverse=True)), bold_lead="Production.")
    para("The slope multipliers move with the law in every model with free multipliers (A18; the values are in the "
         "table of section 1); which slope bins carry extra activity is not robust to the night-side law. kappa is at "
         "its lower bound in all ten runs (summary.csv kappa_at_bound).", bold_lead="Multipliers and kappa.")
    bullet("Not converged: the last fits of %s stopped at max_nfev; reaching max_nfev is not convergence." % ", ".join(maxn),
           lead="Caveats.")
    bullet("kappa at its lower bound in every run: the fits prefer the H2O law for the whole signal; total pressure only, "
           "so this does not measure CO2 as absent.")
    bullet("C with P_m 1e-4: the 70-80 / 80-90 deg multipliers (C f1/1e-4 %s / %s, C f10/1e-4 %s / %s) are not "
           "determined (A12); their Q_H2O is set by unseen steep terrain times these factors, not by the data." % (
               sx("C f1/1e-4", "mult_70_80", "multiplier 70_80", "%.3g"),
               sx("C f1/1e-4", "mult_80_90", "multiplier 80_90", "%.3g"),
               sx("C f10/1e-4", "mult_70_80", "multiplier 70_80", "%.3g"),
               sx("C f10/1e-4", "mult_80_90", "multiplier 80_90", "%.3g")))
    bullet("Q_H2O totals are an extrapolation of each model onto unseen terrain (A15 still applies with the new law, "
           "A18); the seen production is what the data constrain.")
    bullet("C f10/1e-4: no old-law maps exist, so no inspection figure and no old-law Q_H2O; its seen Q_H2O is from the C "
           "report (c5 fit state, before the last fine tuning). D2: one fit, no cycles (no convergence figure), no fine "
           "tuning; its fine-tuning templates come from old-law fits (A18).")
    bullet("The full-set categories contribution density, source areas, cliffs and regions were not yet produced for the "
           "new-law runs.")

    # ---------------- sections
    fig_no = 0
    for si, (sec, sname) in enumerate(SECTIONS):
        N.ctx = "section %s (intro)" % sec
        gs = [g for g in groups if g["section"] == sec]
        SRC["src"], SRC["fb"] = gsrc(gs[0]) + FR, None
        new_section("L")
        doc.add_paragraph("%d. %s" % (si + 3, sname), style="Heading 1")
        para(SECTION_INTRO[sec], bold_lead="About this section.")
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
            doc.add_paragraph("%d.%d %s" % (si + 3, gi + 1, g["title"]), style="Heading 2")
            if "all" in g["paths"]:
                ftxt = os.path.relpath(g["paths"]["all"], ROOT)
            else:
                first = [p for p in g["paths"].values() if p][0]
                ftxt = os.path.relpath(first, ROOT) + " (and the other runs' files, Appendix A)"
            para("Files: " + ftxt, size=7.5, after=2)
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
                    why, src = WHY_MISSING[g["why"]]
                    rr = p1.add_run(T("not produced for this run: %s (%s)" % (why, os.path.relpath(src, ROOT)), [src]))
                    rr.italic = True
                    rr.font.size = Pt(8)
            fixed(t, [W / cols] * cols)
            cell_margins(t, 40)
            cp = doc.add_paragraph(style=cap_style)
            miss = [w for w in g["show"] if w != "all" and not g["paths"].get(w)]
            who_txt = ("one multi-run figure" if g["section"] == "overview" else
                       "identical in all runs (one image shown)" if g["identical"] else
                       "all ten runs (D2, softplus, A, B, D, C f1/1e-4, C f10/1e-4, C f1/0.01, C f10/0.01, C f10/0.01 off)"
                       if not miss else
                       "the ten runs in their fixed order, not produced for %s" % ", ".join(miss))
            rcap = cp.add_run(T("Figure %d: %s — %s; %s." % (fig_no, g["title"], who_txt, g["law"]), SRC["src"]))
            rcap.font.size = Pt(9)
            rcap.bold = True
            doc.add_paragraph("Summary", style="Heading 3")
            summarize(g, para, data_table)
    # ---------------- appendix A
    N.ctx = "Appendix A"
    SRC["src"], SRC["fb"] = [], ("file path or checksum (Appendix A)", ROOT)
    new_section("L")
    doc.add_paragraph("Appendix A. Source files per figure (SHA-256)", style="Heading 1")
    para("Paths relative to %s/ (originals: .png as listed and .eps with the same name; captions <name>_caption.txt). "
         "Full SHA-256 of every source PNG." % ROOT, size=8.5)
    rows = []
    for g in groups:
        lines = [g["title"] + (" [identical in all runs]" if g["identical"] else "")]
        for who, p in g["paths"].items():
            if p:
                lines.append("%s: %s  sha256 %s" % (who, os.path.relpath(p, ROOT), g["sha"][who]))
            else:
                lines.append("%s: none (%s)" % (who, WHY_MISSING[g["why"]][0]))
        rows.append(["Figure %d" % g["fig"], lines])
    data_table(["figure", "plot and its source PNGs"], rows, [0.8, 8.8], size=6)
    # ---------------- appendix B
    N.ctx = "Appendix B"
    SRC["src"], SRC["fb"] = [os.path.abspath(__file__)], ("build parameter / build record", os.path.abspath(__file__))
    doc.add_paragraph("Appendix B. How this document was built", style="Heading 1")
    for txt in build_notes(args):
        bullet(txt, size=8.5)
    if N.unverified:
        with open(os.path.join(out, "numbers_unverified.json"), "w") as f:
            json.dump(N.unverified, f, indent=1, ensure_ascii=False)
        sys.exit("%d numeric tokens of hand-written text not found in a source file (numbers_unverified.json); "
                 "document not written" % len(N.unverified))
    doc.save(os.path.join(out, DOC_NAME + ".part"))
    os.replace(os.path.join(out, DOC_NAME + ".part"), os.path.join(out, DOC_NAME))


SECTION_INTRO = {
    "overview": "The multi-run figures of the comparison set model_vs_observations_nonneg_options_stage2_q002/: the "
                "new-law best of each option (softplus reference, refitted D2, A, B, best C, D) and the five option-C "
                "chains with the softplus reference, each against the COPS data (time, altitude, NG, RG, NG and RG) and "
                "as log residuals, and the data cost per cycle of all chains. One image per figure. The page "
                "source_rate_q002/cluster/ (F_H2O panels (a) and (b) of every run in one tall figure) is not inserted: "
                "it repeats the production-rate panel clusters of section 6.",
    "model_data": "Per run: the final state of the new-law run (solid) and of its old-law counterpart (dashed) against "
                  "the COPS observations (934 NG, 1150 RG; the 95 late-window observations as the gauges' 10-s / "
                  "5-readout response average W10x5), the log residuals of both, and the data cost after every cycle's "
                  "fit and fine tuning for both laws.",
    "activity": "The fitted H2O activity factor A of each new-law final state on the seen facets (log classes per panel "
                "and run; light green = seen with A = 0; grey = not seen during the fit), with the fit-independent "
                "effective gravitational slope and cliff templates.",
    "production_rate": "The H2O production rate F_H2O = q_H2O A per facet at 10:33:10 and 10:39:26 UTC with the NEW law "
                       "(q_min = 0.02 q_max on night and cast-shadow facets, q_min + (q_max - q_min) cos SZA on day "
                       "facets), same layouts as the activity maps, relief shading under the colours.",
    "inspection": "Per run and time: (a, b) the new-law F_H2O and (c, d) the old-law counterpart's F_H2O on ONE set of "
                  "colour classes (so the colours compare directly), (e, f) the new law's illumination factor q/q_max "
                  "on the whole surface. Not for C f10/1e-4 (no old-law maps).",
}


def build_notes(args):
    return [
        "Generator: %s (run copy; version-controlled copy model_ext/figures/q002_comparison_docx.py, SHA256 in "
        "model_ext/RUN_COPIES.sha256), derived from fullsets_comparison_docx.py; run with %s (python-docx %s, Pillow "
        "%s; separate venv, not the project .venv)." % (os.path.abspath(__file__), sys.executable, args.docx_version,
                                                         args.pil_version),
        "Command: %s; Slurm job %s on %s; %s UTC." % (" ".join(sys.argv), os.environ.get("SLURM_JOB_ID", "-"),
                                                   os.uname().nodename, args.started),
        "Images: every source PNG hashed (SHA-256); a plot type whose PNGs are identical in all runs is shown once. "
        "Downscaled copies in images/: maps as JPEG (quality 78), line plots (model/data, residuals, convergence) as "
        "PNG with 256 colours; width = min(native, max(1000 px, 300 dpi at the printed width), 1400 px); RGBA flattened "
        "on white. Layout: per cluster the page orientation and grid (2-5 columns) that gives the largest images with "
        "all runs on one page (US Letter, margins 0.6 in).",
        "Numbers: read at build time from summary.csv, the residual-check json files, the inspection files "
        "(inspect_<hhmmss>.json), the captions (regular expressions), the map manifests (class tables), "
        "slope_bins_h2o_<hhmmss>.json and the plotted tables (tables/<run>.dat: residual extremes and the values at the "
        "last observation, recomputed); listed in numbers_used.json with source, locator and where used. Literal "
        "numbers of hand-written sentences are checked against the captions, READMEs, summary.txt, FINDINGS.md and "
        "OUTPUT_VISUALIZATION.md (the build refuses to write the document otherwise); numbers in 'Observation from the "
        "figures' paragraphs are recorded as figure readings.",
        "Validation (validate step): the document reopened with python-docx (captions, images, placeholders, sections, "
        "tables, image and table widths against the pages, zip entries), the OOXML schema check of the docx skill's "
        "validate.py; LibreOffice is not installed on the cluster (login and compute nodes), so no PDF; validation/ has "
        "layout previews drawn from the same geometry.",
    ]


# ------------------------------------------------------------------------------------------------ summaries
def summarize(g, para, data_table):
    F = FAMILY[g["fam"]]
    para(F["about"](g), bold_lead="What the plot shows.")
    missing = [r for r, p in g["paths"].items() if r != "all" and not p]
    if missing:
        why, src = WHY_MISSING[g["why"]]
        para("Not produced for %s: %s (%s)." % (", ".join(missing), why, os.path.relpath(src, ROOT)),
             bold_lead="Runs without this plot.")
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
    obs = F.get("observe")
    if obs:
        o = obs(g) if callable(obs) else obs
        if o:
            para(o, bold_lead="Observation from the figures:", force=True, fb=fig_fallback(g), src=[])


def first_caption(g):
    p = [x for x in g["paths"].values() if x][0]
    return caption_of(p)


def cap_sentences(g, n=1, stop=None):
    t = tload(first_caption(g))
    parts = re.split(r"(?<=[a-z0-9\)])\. (?=[A-Z(])", t)
    s = ". ".join(parts[:n])
    if stop:
        s = s.split(stop)[0]
    return s.rstrip(".") + "."


LAW_SENT = re.compile(r"NEW H2O law: .*?\(not intended\)")


def law_sentence(g):
    m = LAW_SENT.search(tload(first_caption(g)))
    return (m.group(0) + ".") if m else ""


# overview -------------------------------------------------------------------------------------------------------------
TAB2RUN = {RUN[r]["csv"] + ".dat": r for r in RK}


def ov_fits(g):
    p = MV + "/12M_q002_%s_residuals_check.json" % g["set"]
    out = []
    for k, v in jload(p)["fits"].items():
        out.append((TAB2RUN[os.path.basename(v["table"])], k, v, p))
    return out


def ov_about(g):
    if g["title"].endswith("against altitude") and "residuals" in first_caption(g):
        return (cap_sentences(g, 1) + " Log residuals ln(model / observation) of the runs in the table against the gauge "
                "altitude (logarithmic, descent from left to right): (a) NG, (b) RG; one point per fitted observation; for "
                "the 95 late-window observations the model is the gauges' W10x5 response average, as fitted.")
    return cap_sentences(g, 1) + " " + law_sentence(g)


def ov_numbers(g):
    hdr = ["run", "final data cost", "rms ln NG all (late)", "rms ln RG all (late)"]
    rows = []
    for r, k, v, p in ov_fits(g):
        rows.append([r, sx(r, "data_cost", "final data cost", e4),
                     "%s (%s)" % (N.add(r, "rms ln NG (all fitted rows)", v["NG"]["rms"], p, "fits['%s'].NG.rms" % k,
                                        "%.4f"),
                                  N.add(r, "rms ln NG (late window)", v["NG"]["rms_late"], p,
                                        "fits['%s'].NG.rms_late" % k, "%.4f")),
                     "%s (%s)" % (N.add(r, "rms ln RG (all fitted rows)", v["RG"]["rms"], p, "fits['%s'].RG.rms" % k,
                                        "%.4f"),
                                  N.add(r, "rms ln RG (late window)", v["RG"]["rms_late"], p,
                                        "fits['%s'].RG.rms_late" % k, "%.4f"))])
    return hdr, rows, [1.4, 1.4, 1.6, 1.6]


def ov_compare(g):
    fits = ov_fits(g)
    c = {r: (float(sval(r, "data_cost")), sx(r, "data_cost", "final data cost", e4)) for r, k, v, p in fits}
    ln = {r: (v["NG"]["rms_late"], N.add(r, "rms ln NG (late window)", v["NG"]["rms_late"], p,
                                          "fits['%s'].NG.rms_late" % k, "%.4f")) for r, k, v, p in fits}
    lr = {r: (v["RG"]["rms_late"], N.add(r, "rms ln RG (late window)", v["RG"]["rms_late"], p,
                                          "fits['%s'].RG.rms_late" % k, "%.4f")) for r, k, v, p in fits}
    return ("Final data cost: %s. Late-window rms ln NG: %s; RG: %s." % (order_text(c), order_text(ln), order_text(lr)))


def ov_observe(g):
    if g["set"] == "best" and g["fam"] == "ov":
        return ("Above the late window the curves of the six runs lie close together at this scale. In the late window "
                "D2's model NG stays below the data and below the other runs at the end, A and D put the model RG above "
                "the RG data near touchdown, and C f10/1e-4 follows the NG data most closely.")
    return ("Above the late window the six curves lie close together at this scale; the differences are in the late "
            "window, where the two P_m 1e-4 chains nearly coincide and C f10/0.01 off departs most from the RG data "
            "near touchdown.")


def ovc_about(g):
    return cap_sentences(g, 2)


def ovc_numbers(g):
    hdr = ["run", "data cost after each fit (c1 ... c5)", "after each fine tuning", "final", "last fit"]
    rows = []
    for r in RK:
        rr, _ = srow(RUN[r]["csv"])
        if rr["history_fit_data"]:
            rows.append([r, sx(r, "history_fit_data", "data cost history, fits", "%s", conv=None),
                         sx(r, "history_ft_data", "data cost history, fine tunings", "%s", conv=None),
                         sx(r, "data_cost", "final data cost", e4),
                         "%s (%s)" % ("max_nfev" if "maximum" in rr["fit_termination"] else "ftol",
                                      sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s))))])
        else:
            rows.append([r, "one fit (no cycles)", "-", sx(r, "data_cost", "final data cost", e4),
                         "max_nfev (%s)" % sx(r, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s)))])
    return hdr, rows, [1.2, 3.0, 3.0, 1.0, 1.1]


def ovc_compare(g):
    c = {r: (float(sval(r, "data_cost")), sx(r, "data_cost", "final data cost", e4)) for r in RK}
    return "Final data cost: %s (the dotted line of the figure is the refitted D2)." % order_text(c)


# per-run model / data, new vs old ------------------------------------------------------------------------------------
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
            row += ["%s (%s)" % (t_min_ng(r)[1], t_min_ng(r, True)[1]),
                    "%s (%s)" % (t_max_rg(r)[1], t_max_rg(r, True)[1])]
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
        v = {r: t_last(r, "NG") for r in RK}
        out.append("ln(model/obs) NG at the last observation, new law: %s." % order_text(v))
    if fam in ("vo", "vo_rg", "vo_ngrg"):
        v = {r: t_last(r, "RG") for r in RK}
        vo_ = {r: t_last(r, "RG", True) for r in RK}
        lower = [r for r in RK if v[r][0] < vo_[r][0]]
        out.append("ln(model/obs) RG at the last observation, new law: %s; lower than with the old law in %s." % (
            order_text(v), "all ten runs" if len(lower) == len(RK) else (", ".join(lower) if lower else "no run")))
    return " ".join(out)


def vo_observe(g):
    return ("In every run the new-law (solid) and old-law (dashed) curves lie close together above the late window at "
            "this scale; the visible differences are in the late window. The two P_m 1e-4 chains show nearly coinciding "
            "curves; in C f10/0.01 off the new-law RG falls well below the RG data in the last metres; in D2 the "
            "new-law NG stays below the data and below the old-law curve at the end.")


def vr_about(g):
    return (cap_sentences(g, 1) + " Log residuals ln(model / observation) of the new-law run and of its old-law "
            "counterpart against the gauge altitude (logarithmic, descent from left to right): (a) NG, (b) RG; one point "
            "per fitted observation (934 NG, 1150 RG), the late window shaded; for the 95 late-window observations the "
            "model is the gauges' W10x5 response average, as fitted. Positive: model above the data.")


def vr_numbers(g):
    hdr = ["run", "rms ln NG all (late) new", "old", "rms ln RG all (late) new", "old"]
    rows = []
    for r in RK:
        p = g["paths"][r]
        if not p:
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


def vr_observe(g):
    return ("D2's new-law NG residual falls to the most negative value at the end of the descent, while in the other "
            "five runs shown the NG residual rises to positive values at the last points with both laws; above the late "
            "window the new- and old-law residuals nearly coincide in every run shown.")


def cv_about(g):
    m = re.search(r"^Convergence of the stage-2 chain \S+: (.*?); solid black", tload(first_caption(g)), re.S)
    return "Convergence of each stage-2 chain: " + (m.group(1) if m else "") + ". " + (
        "Solid black: the NEW H2O law (q_min = 0.02 q_max); dashed grey: the OLD-law chain (q_min = 0.2 q_max); both "
        "against the same observations; each panel has its own data-cost axis.")


def cv_numbers(g):
    hdr = ["run", "fits c1 ... c5, new law", "fine tunings c1 ... c5, new law", "final new", "final old",
           "last fit new"]
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
        drop[r] = (fin / h[0] - 1.0, N.add(r, "change of the data cost from the c1 fit to the final state", 100.0 * (fin / h[0] - 1.0),
                                          p, "summary.csv: data_cost / first value of history_fit_data - 1", sgn_pct))
    lower = [r for r in drop if float(sval(r, "data_cost")) < float(sval(r, "data_cost", old=True))]
    return ("Change of the data cost from the cycle-1 fit to the final state, new law: %s. Final new-law cost below the "
            "old-law one: %s." % (order_text(drop), ", ".join(lower) if lower else "none"))


def cv_observe(g):
    return ("The two option-C chains with P_m 0.01 and templates on drop within the first cycle and then stay flat; "
            "C f10/0.01 off stays flat from the first fit; in A the new-law curve lies below the old-law one in every "
            "cycle; in the other chains it lies above.")


# maps ---------------------------------------------------------------------------------------------------------------
PANEL_NAME = {"latlon_global": "lon/lat global", "latlon_zoom": "lon/lat impact area", "3d_a": "3-D (a) nucleus",
              "3d_b": "3-D (b) impact region", "3d_c": "3-D (c) late descent"}


def man_path(r, g):
    return srp(r, "checks/activity_manifest.json" if g["q"] == "activity" else "checks/h2o_%s_manifest.json" % g["t"])


def cls_range(r, g, panel):
    p = man_path(r, g)
    ct = jload(p)["class_tables"][panel]
    lo, hi, n = ct["classes"][0]["log10_from"], ct["classes"][-1]["log10_to"], len(ct["classes"])
    s = N.add(r, "colour classes %s %s: log10 range and count" % (g["q"] + (g["t"] or ""), panel), [lo, hi, n], p,
              "class_tables/%s/classes[0].log10_from, [-1].log10_to, len" % panel,
              lambda v: "10^%g-10^%g (%d)" % tuple(v))
    return s, hi - lo


def tt(g):
    return "%s:%s:%s" % (g["t"][:2], g["t"][2:4], g["t"][4:])


def qcap(r, g):
    p = srp(r, "production_rate/12M_templates_h2o_production_rate_%sUTC_caption.txt" % g["t"])
    return capx(p, r"total H2O production Q_H2O = sum F_H2O area = (\S+) s\^-1", r,
                "Q_H2O total at %s [s^-1] (caption)" % tt(g))


def map_numbers(g):
    hdr = ["run"]
    if g["q"] == "activity":
        hdr += ["active area A > 0, whole / seen (inspection)"]
    else:
        hdr += ["Q_H2O total [s^-1] (caption)", "seen share of Q", "night + shadow share of Q / of seen Q"]
    for k, panel in g["panels"]:
        hdr.append("classes %s: range (n)" % PANEL_NAME[panel])
    if g["fam"] == "map_tpl" and g["q"] == "activity":
        hdr.append("multipliers 40-50 ... 90-180 deg")
    if g["fam"] == "map_tpl" and g["q"] != "activity":
        hdr.append("seen production from >= 60 deg")
    rows = []
    for r in RK:
        row = [r]
        t = "103310" if g["q"] == "activity" else g["t"]
        if g["q"] == "activity":
            a = ix(r, t, ["new", "active_area_fraction"], "active area fraction (whole) 10:33:10 inspection", pct, 100.0)
            b = ix(r, t, ["new", "active_area_fraction_seen"], "active area fraction (seen) 10:33:10 inspection", pct,
                   100.0)
            row.append("%s / %s" % (a, b) if a else "n/a (no inspection)")
        else:
            row.append(qcap(r, g))
            s = ix(r, t, ["new", "seen_share_of_Q"], "seen share of Q_H2O at %s" % tt(g), pct, 100.0)
            row.append(s if s else "n/a (no inspection)")
            d1 = ix(r, t, ["new", "dark_share_of_Q"], "night + cast-shadow share of Q_H2O at %s" % tt(g), pct, 100.0)
            d2 = ix(r, t, ["new", "dark_share_of_seen_Q"], "night + cast-shadow share of the seen Q_H2O at %s" % tt(g),
                    pct, 100.0)
            row.append("%s / %s" % (d1, d2) if d1 else "n/a")
        for k, panel in g["panels"]:
            row.append(cls_range(r, g, panel)[0])
        if g["fam"] == "map_tpl" and g["q"] == "activity":
            ms = [sx(r, "mult_" + b, "multiplier " + b, "%.3g") for b in MULT]
            row.append(", ".join(ms) if all(ms) else "none")
        if g["fam"] == "map_tpl" and g["q"] != "activity":
            p = srp(r, "checks/slope_bins_h2o_%s.json" % g["t"])
            v = sum(b["production_share_seen"] for b in jload(p)["bins"][3:])
            row.append(N.add(r, "seen H2O production share from >= 60 deg at %s" % tt(g), 100.0 * v, p,
                             "sum bins[3:].production_share_seen", pct))
        rows.append(row)
    return hdr, rows, [1.1] + [1.5] * (len(hdr) - 1)


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
    seen = {r: capx(caption_of(g["paths"][r]), pat, r, "seen facets during the fit (caption)") for r in RK
            if g["paths"][r]}
    same = len(set(seen.values())) == 1
    law = (" F_H2O = q_H2O A with the NEW law q_min = 0.02 q_max (captions)." if g["q"] != "activity" else "")
    return (first.rstrip(".") + ". Colour only on the facets seen during the fit (%s%s); light green: seen with %s = 0; "
            "grey: not seen. Classes scanned per panel and per run (captions; checks/*_class_tables.txt), so the colour "
            "scales differ between the runs.%s" % (list(seen.values())[0], ", the same in all runs" if same else
                                                     "; differs between runs", "A" if g["q"] == "activity" else "F_H2O",
                                                     law))


def map_compare(g):
    out = []
    if g["q"] == "activity":
        af = {}
        for r in RK:
            v = ixv(r, "103310", ["new", "active_area_fraction_seen"])
            if v is not None:
                af[r] = (v, ix(r, "103310", ["new", "active_area_fraction_seen"],
                               "active area fraction (seen) 10:33:10 inspection", pct, 100.0))
        out.append("Share of the seen area with A > 0 (inspection files; none for C f10/1e-4): %s." %
                   order_text(af, reverse=True))
    else:
        q = {r: (float(qcap(r, g)), qcap(r, g)) for r in RK}
        ss = {}
        for r in RK:
            v = ixv(r, g["t"], ["new", "seen_share_of_Q"])
            if v is not None:
                ss[r] = (v, ix(r, g["t"], ["new", "seen_share_of_Q"], "seen share of Q_H2O at %s" % tt(g), pct, 100.0))
        out.append("Q_H2O at %s UTC: %s s^-1; seen share: %s." % (tt(g), order_text(q), order_text(ss, reverse=True)))
        dn = {}
        for r in RK:
            p = srp(r, "production_rate/12M_templates_h2o_production_rate_%sUTC_caption.txt" % g["t"])
            dn[r] = capx(p, r"facets at this time: ([\d,]+ day, [\d,]+ night, [\d,]+ in cast shadow)", r,
                         "day / night / cast-shadow facets at %s" % tt(g))
        if len(set(dn.values())) == 1:
            out.append("The illumination classes are the same in all runs (%s; captions)." % dn[RK[0]])
        if g["fam"] == "map_tpl":
            parts = []
            for r in RK:
                p = srp(r, "checks/slope_bins_h2o_%s.json" % g["t"])
                bins = jload(p)["bins"]
                i = max(range(7), key=lambda j: bins[j]["production_share_seen"])
                parts.append("%s %s deg (%s)" % (r, BIN_LABEL[i], N.add(
                    r, "largest slope bin of the seen production at %s" % tt(g), 100.0 * bins[i]["production_share_seen"],
                    p, "bins[%d].production_share_seen (max over bins)" % i, pct0)))
            out.append("Largest slope bin of the seen H2O production (slope_bins_h2o_%s.json): %s." % (g["t"],
                                                                                                    "; ".join(parts)))
    k, panel = g["panels"][0]
    spans = {}
    for r in RK:
        s, d = cls_range(r, g, panel)
        spans[r] = (d, N.add(r, "decades covered by the classes, %s %s" % (g["q"], panel), d, man_path(r, g),
                             "class_tables/%s: log10_to(last) - log10_from(first)" % panel, "%g"))
    out.append("Decades spanned by the colour classes in %s: %s." % (PANEL_NAME[panel], order_text(spans)))
    return " ".join(out)


def map_observe(g):
    q = "A" if g["q"] == "activity" else "F_H2O"
    s = ("D2 and softplus: most of the seen surface is light green (%s = 0) with a few isolated active patches; A, B "
         "and D: smooth fields over the whole seen surface; the five option-C runs: a patchwork of polygonal cells (the "
         "hat-function basis), in the P_m 0.01 runs with large areas in the lowest classes and a few bright cells." % q)
    if g["fam"] == "map_3d":
        s += (" In the close-up of the late descent (c) softplus shows an elliptical active patch next to the impact "
              "point, D2 thin elongated active streaks.")
    if g["fam"] == "map_tpl":
        s += " Panels (c)-(f) are the same in all runs (the separate panel PNGs are byte-identical)."
    return s


def geom_about(g):
    return cap_sentences(g, 1) + (" Effective gravitational slope: angle between the inward facet normal and the "
                                  "effective gravity; values above 90 deg are overhangs. Cliff templates: slope >= 45 "
                                  "deg (light) and >= 60 deg (dark). Geometry only: independent of the fit and of the H2O "
                                  "law.")


# inspection ---------------------------------------------------------------------------------------------------------
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
    rat = {}
    for r in RK:
        if g["paths"][r]:
            qn_, qo = ixv(r, t, ["new", "Q_H2O_s-1"]), ixv(r, t, ["old", "Q_H2O_s-1"])
            rat[r] = (qn_ / qo, N.add(r, "Q_H2O new / old at %s" % tt(g), qn_ / qo, insp_path(r, t),
                                      "inspect_%s.json: new/Q_H2O_s-1 / old/Q_H2O_s-1" % t, "%.3f"))
    sr = {}
    for r in RK:
        if g["paths"][r]:
            a, b = ixv(r, t, ["new", "Q_H2O_seen_s-1"]), ixv(r, t, ["old", "Q_H2O_seen_s-1"])
            sr[r] = (a / b, N.add(r, "seen Q_H2O new / old at %s" % tt(g), a / b, insp_path(r, t),
                                  "inspect_%s.json: new/Q_H2O_seen_s-1 / old/Q_H2O_seen_s-1" % t, "%.3f"))
    more = [r for r in rat if rat[r][0] < sr[r][0]]
    return ("Q_H2O new law / old law at %s UTC: %s; seen Q_H2O new / old: %s. The total falls more than the seen "
            "production in %d of the %d runs shown%s." % (tt(g), order_text(rat), order_text(sr), len(more), len(rat),
                                                          "" if len(more) == len(rat) else " (" + ", ".join(more) + ")"))


def insp_observe(g):
    return ("Within each run the new-law (a, b) and old-law (c, d) maps show a similar spatial pattern of the seen "
            "production at this scale; panels (e, f) look the same in all runs (same time and geometry).")


FAMILY = {
    "ov": dict(about=ov_about, numbers=ov_numbers, compare=ov_compare, observe=ov_observe),
    "ov_conv": dict(about=ovc_about, numbers=ovc_numbers, compare=ovc_compare),
    "vo": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare, observe=vo_observe),
    "vo_ng": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare, observe=vo_observe),
    "vo_rg": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare, observe=vo_observe),
    "vo_ngrg": dict(about=vo_about, numbers=vo_numbers, compare=vo_compare, observe=vo_observe),
    "vr": dict(about=vr_about, numbers=vr_numbers, compare=vr_compare, observe=vr_observe),
    "cv": dict(about=cv_about, numbers=cv_numbers, compare=cv_compare, observe=cv_observe),
    "map_tpl": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_latlon": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_3d": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_panel": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "geom_panel": dict(about=geom_about),
    "insp": dict(about=insp_about, numbers=insp_numbers, compare=insp_compare, observe=insp_observe),
}


# ------------------------------------------------------------------------------------------------ build / validate
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


LINE_FAMS = ("ov", "ov_conv", "vo", "vo_ng", "vo_rg", "vo_ngrg", "vr", "cv")


def cmd_build(args):
    out = os.path.abspath(args.out)
    if os.path.exists(os.path.join(out, DOC_NAME)):
        sys.exit("refusing to overwrite %s (move the folder aside first)" % os.path.join(out, DOC_NAME))
    os.makedirs(out, exist_ok=True)
    args.started = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    import docx
    args.docx_version = getattr(docx, "__version__", "?")
    args.pil_version = PIL.__version__
    groups = catalogue()
    jobs, img = [], {}
    for gi, g in enumerate(groups):
        g["id"] = gi
        g["sha"] = {who: sha256(p) for who, p in g["paths"].items() if p}
        avail = [who for who, p in g["paths"].items() if p]
        g["identical"] = (g["section"] != "overview" and len(avail) == len(RK) and len(set(g["sha"].values())) == 1)
        if g["section"] == "overview":
            g["show"] = ["all"]
        elif g["identical"]:
            g["show"] = [RK[0]]
        else:
            g["show"] = list(RK)
        sizes = {}
        for who in g["show"]:
            if g["paths"].get(who):
                with Image.open(g["paths"][who]) as im:
                    sizes[who] = im.size
        g["native"] = sizes
        g["layout"] = layout([h / w for (w, h) in sizes.values()], len(g["show"]))
        kind = "png" if g["fam"] in LINE_FAMS else "jpg"
        for who in sizes:
            nat = sizes[who][0]
            target = int(min(nat, max(1000, round(g["layout"]["width_in"] * 300)), 1400))
            base = re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.basename(g["paths"][who])[:-4])
            dst = os.path.join(out, "images", g["section"], "%03d_%s__%s.%s" % (gi, base, re.sub(r"[^A-Za-z0-9]", "", who), kind))
            jobs.append((g["paths"][who], dst, kind, target))
            img[(gi, who)] = dict(src=g["paths"][who], dst=dst, kind=kind, target=target)
    with ProcessPoolExecutor(args.jobs) as ex:
        res = list(ex.map(prep_image, jobs))
    tot = 0
    bydst = {}
    for (src, dst, kind, target), (d, size, nbytes) in zip(jobs, res):
        tot += nbytes
        bydst[d] = (size, nbytes)
    for v in img.values():
        v["px"], v["bytes"] = bydst[v["dst"]]
    args.img_note = "%d images, %.1f MB, JPEG for maps and 256-colour PNG for line plots, 1000-1400 px wide" % (
        len(jobs), tot / 1e6)
    build_doc(out, groups, img, args)
    with open(os.path.join(out, "numbers_used.json"), "w") as f:
        json.dump(dict(created_utc=args.started, generator=os.path.abspath(__file__), n=len(N.items), numbers=N.items),
                  f, indent=1, ensure_ascii=False)
    man = dict(created_utc=args.started,
               finished_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S"),
               command=sys.argv, python=sys.executable, python_docx=args.docx_version, pillow=args.pil_version,
               slurm_job=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename, amps_head=amps_head(),
               generator_sha256=sha256(os.path.abspath(__file__)), document=os.path.join(out, DOC_NAME),
               document_bytes=os.path.getsize(os.path.join(out, DOC_NAME)), images_bytes=tot,
               groups=[dict(fig=g["fig"], id=g["id"], section=g["section"], title=g["title"], fam=g["fam"], law=g["law"],
                            identical=g["identical"], show=g["show"], layout=g["layout"], native=g["native"],
                            paths=g["paths"], sha256=g["sha"],
                            missing=[w for w, p in g["paths"].items() if not p],
                            images={who: dict(dst=os.path.relpath(img[(g["id"], who)]["dst"], out),
                                              px=img[(g["id"], who)]["px"], bytes=img[(g["id"], who)]["bytes"])
                                    for who in g["show"] if g["paths"].get(who)}) for g in groups])
    with open(os.path.join(out, "build_manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    print("document %s: %.1f MB; %d clusters / figures (%d identical in all runs, %d missing-run placeholders); %d images "
          "(%.1f MB); %d numbers" % (os.path.join(out, DOC_NAME), man["document_bytes"] / 1e6, len(groups),
                                     sum(g["identical"] for g in groups),
                                     sum(1 for g in groups for w in g["show"] if w != "all" and not g["paths"].get(w)),
                                     len(jobs), tot / 1e6, len(N.items)))


def cmd_validate(args):
    """reopen the document; count captions / images / placeholders; check image and table widths; previews"""
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
    k = 0
    exp_imgs = 0
    exp_place = 0
    for g in man["groups"]:
        n = sum(1 for w in g["show"] if g["paths"].get(w) or w == "all")
        exp_imgs += n
        exp_place += sum(1 for w in g["show"] if w != "all" and not g["paths"].get(w))
        W, H = usable(g["layout"]["orient"])
        cols = g["layout"]["cols"]
        sizes = pics[k:k + n]
        k += n
        wmax = max(w for w, h in sizes)
        hrow = max(h for w, h in sizes)
        rows = math.ceil(len(g["show"]) / cols)
        if wmax > W / cols - 0.05:
            problems.append("Figure %d: image width %.2f in > cell %.2f in" % (g["fig"], wmax, W / cols))
        if rows * (hrow + LABEL + 0.05) + OVERHEAD - 0.4 > H:
            problems.append("Figure %d: grid height > page" % g["fig"])
        g["_check"] = dict(width_in=round(wmax, 3), max_h_in=round(hrow, 3), rows=rows)
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
    res = dict(document=path, bytes=os.path.getsize(path), mb=round(os.path.getsize(path) / 1e6, 2), headings=heads,
               figure_captions=len(caps), figure_groups=len(man["groups"]), captions_numbered_1_to_n=numbered,
               captions_state_the_law=law_ok, identical_groups=sum(g["identical"] for g in man["groups"]),
               inline_images=len(shapes), expected_images=exp_imgs, media_files=media,
               placeholders=placeholders, expected_placeholders=exp_place, sections=len(doc.sections),
               tables=len(doc.tables), tables_width_checked=nt, numbers_recorded=nums["n"], problems=problems,
               libreoffice="not available (no soffice on login or compute nodes; not on conda-forge): no PDF render",
               checks=[dict(fig=g["fig"], **g["_check"]) for g in man["groups"]])
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
    for g in [g for g in man["groups"] if g["fig"] in args.preview]:
        o = g["layout"]["orient"]
        pw, ph = PAGE[o]
        im = Image.new("RGB", (int(pw * dpi), int(ph * dpi)), "white")
        dr = ImageDraw.Draw(im)
        dr.rectangle([MARGIN * dpi, MARGIN * dpi, (pw - MARGIN) * dpi, (ph - MARGIN) * dpi], outline=(200, 200, 200))
        dr.text((MARGIN * dpi, MARGIN * dpi), "x.y %s" % g["title"], fill=(31, 58, 95), font=font)
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
        im.save(os.path.join(out, "validation", "layout_preview_fig%02d.png" % g["fig"]))
    print(json.dumps({k: v for k, v in res.items() if k != "checks"}, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build")
    b.add_argument("--out", required=True)
    b.add_argument("--jobs", type=int, default=8)
    v = sp.add_parser("validate")
    v.add_argument("--out", required=True)
    v.add_argument("--preview", type=int, nargs="*", default=[1, 13, 14, 16, 19, 20, 21, 29, 36, 47, 48, 49])
    args = ap.parse_args()
    if args.cmd == "build":
        cmd_build(args)
    else:
        cmd_validate(args)


if __name__ == "__main__":
    main()
