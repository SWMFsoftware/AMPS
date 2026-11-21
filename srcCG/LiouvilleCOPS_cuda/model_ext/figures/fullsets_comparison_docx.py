#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/fullsets_comparison_docx.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/fullsets_comparison_docx.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 1f130e76e0081964023713bd89c099594fc5c75ea69358d3843bc9b615873bc9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Word document comparing the six full figure sets of the 12M stage-2 final states plot by plot (2026-11-17)
# RUN: separate venv cops_output/tools/docx_venv (python-docx, Pillow); build + validate in a Slurm job (HOW-TO-RERUN.md of figures_cluster/12M/full_sets_stage2_comparison/)
# ==== end of header ====
"""Word document comparing the six full figure sets of the 12M stage-2 final states (2026-11-17).

User request 2026-11-17: "prepare word doc that compares all the cases: group the same plot for all cases and provide
summary and that for every plot". Inputs: the six full sets under figures_cluster/12M/ (full_set_stage2_{softplus_reg,
A_reg,B_reg,D_reg,C_f1_pm1em4_reg}/, full_set_D2/), the final stage-2 comparison
(model_vs_observations_nonneg_options_stage2_final/: summary.csv, 12M_stage2_reg_best_residuals_check.json).
Every number in the document is read from a file and listed in numbers_used.json (value, formatted text, source file,
locator, where used). Images: downscaled copies of the PNGs (images/), originals referenced by path (Appendix A).

    fullsets_comparison_docx.py build    --out <new folder> [--jobs 8]
    fullsets_comparison_docx.py validate --out <folder>
Interpreter: a separate venv with python-docx and Pillow ($COPS_OUTPUT_ROOT/tools/docx_venv; not the project .venv).
Refuses to overwrite an existing document. Light on memory (< 4 GB); the image step decodes ~430 PNGs (run in a job).
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
STAGE2 = ROOT + "/model_vs_observations_nonneg_options_stage2_final"
FINDINGS = "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/model_ext/FINDINGS.md"
OUTVIS = "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/model_ext/OUTPUT_VISUALIZATION.md"
NUMTOK = re.compile(r"(?<![A-Za-z_\d.])\d+(?:[.,]\d+)*(?:e[+-]?\d+)?")
DOC_NAME = "12M_stage2_full_sets_comparison.docx"
DOC_DATE = "2026-11-17"

# case order fixed by the request: softplus, A, B, D, C, D2
CASES = [
    dict(key="softplus", label="softplus (reference)", set="full_set_stage2_softplus_reg", stem="12M_stage2_softplus_reg",
         csv="softplus_reg", resid="softplus_reg.dat", table="stage2_softplus_reg_final.dat"),
    dict(key="A", label="A sum of squares r3", set="full_set_stage2_A_reg", stem="12M_stage2_A_reg",
         csv="A_reg", resid="A_reg.dat", table="stage2_A_reg_final.dat"),
    dict(key="B", label="B background 1e-2", set="full_set_stage2_B_reg", stem="12M_stage2_B_reg",
         csv="B_reg", resid="B_reg.dat", table="stage2_B_reg_final.dat"),
    dict(key="D", label="D log spread f 0.1", set="full_set_stage2_D_reg", stem="12M_stage2_D_reg",
         csv="D_reg", resid="D_reg.dat", table="stage2_D_reg_final.dat"),
    dict(key="C", label="C local basis f 1, P_m 1e-4", set="full_set_stage2_C_f1_pm1em4_reg",
         stem="12M_stage2_C_f1_pm1em4_reg", csv="C_f1_pm1em4_reg", resid="C_f1_pm1em4_reg.dat",
         table="stage2_C_f1_pm1em4_reg_final.dat"),
    dict(key="D2", label="D2 clipped refit", set="full_set_D2", stem="12M_D2", csv="D2", resid="D2.dat",
         table="D2_final.dat"),
]
CK = [c["key"] for c in CASES]
CASE = {c["key"]: c for c in CASES}
BINS = ["[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]"]
BIN_LABEL = ["0-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-180"]
SECTIONS = ["model_data", "activity", "production_rate", "source_areas", "contribution_density", "cliffs", "regions",
            "visibility"]


# ------------------------------------------------------------------------------------------------ numbers registry
class Numbers:
    """Every number printed in the document goes through add(); numbers_used.json lists them with their source."""

    def __init__(self):
        self.items, self.index, self.ctx = [], {}, "front matter"
        self.ctx_tokens, self.unverified = {}, []

    def add(self, case, quantity, value, source, locator, fmt="%s", note=None):
        text = fmt(value) if callable(fmt) else (fmt % value)
        self.ctx_tokens.setdefault(self.ctx, set()).update(NUMTOK.findall(text))
        k = (case, quantity, source, locator, text)
        if k in self.index:
            it = self.items[self.index[k]]
            if self.ctx not in it["used_in"]:
                it["used_in"].append(self.ctx)
            return text
        v = value
        if isinstance(v, float) and not math.isfinite(v):
            v = str(v)
        self.items.append(dict(id=len(self.items) + 1, case=case, quantity=quantity, value=v, text=text,
                               source=source, locator=locator, note=note, used_in=[self.ctx]))
        self.index[k] = len(self.items) - 1
        return text


N = Numbers()
_cache = {}


def T(text, sources, fallback=None, force=False):
    """hand-written text: every numeric token must occur in one of the source files (recorded in numbers_used.json);
    tokens already registered in this context are skipped; fallback = (locator, source) for figure readings / build
    parameters; anything else is listed as unverified (the build reports it)"""
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


def jload(path):
    if path not in _cache:
        with open(path) as f:
            _cache[path] = json.load(f)
    return _cache[path]


def tload(path):
    if path not in _cache:
        with open(path) as f:
            _cache[path] = f.read()
    return _cache[path]


def setp(c, rel=""):
    return os.path.join(ROOT, CASE[c]["set"], rel)


def nested(d, keys):
    for k in keys:
        d = d[k]
    return d


def fx(c, keys, quantity, fmt="%s", scale=1.0):
    """value from <set>/checks/facts.json"""
    p = setp(c, "checks/facts.json")
    v = nested(jload(p), keys)
    if scale != 1.0:
        v = v * scale
    return N.add(c, quantity, v, p, "facts.json:" + "/".join(str(k) for k in keys), fmt)


def fxv(c, keys):
    return nested(jload(setp(c, "checks/facts.json")), keys)


def csvrow(c):
    p = STAGE2 + "/summary.csv"
    if p not in _cache:
        with open(p) as f:
            _cache[p] = {r["key"]: r for r in csv.DictReader(f)}
    return _cache[p][CASE[c]["csv"]], p


def cx(c, col, quantity, fmt="%s", conv=float):
    r, p = csvrow(c)
    v = conv(r[col]) if conv else r[col]
    return N.add(c, quantity, v, p, "summary.csv row key=%s column %s" % (CASE[c]["csv"], col), fmt)


def resid_entry(c):
    p = STAGE2 + "/12M_stage2_reg_best_residuals_check.json"
    for k, v in jload(p)["fits"].items():
        if os.path.basename(v["table"]) == CASE[c]["resid"]:
            return k, v, p
    raise KeyError(c)


def rx(c, gauge, key, quantity, fmt="%.4f"):
    k, v, p = resid_entry(c)
    return N.add(c, quantity, v[gauge][key], p, "fits['%s'].%s.%s" % (k, gauge, key), fmt)


def caption_path(c, rel):
    return setp(c, rel.replace("<STEM>", CASE[c]["stem"]).replace(".png", "_caption.txt"))


def capx(c, rel, pattern, quantity, group=1, fmt="%s", conv=str):
    p = caption_path(c, rel)
    m = re.search(pattern, tload(p), re.S)
    if not m:
        return N.add(c, quantity, "n/a", p, "regex (no match): " + pattern, "%s", note="not in caption")
    return N.add(c, quantity, conv(m.group(group)), p, "regex group %d: %s" % (group, pattern), fmt)


def manx(c, fname, keys, quantity, fmt="%s", scale=1.0, conv=None):
    p = setp(c, fname)
    v = nested(jload(p), keys)
    if conv:
        v = conv(v)
    if scale != 1.0:
        v = v * scale
    return N.add(c, quantity, v, p, os.path.basename(p) + ":" + "/".join(str(k) for k in keys), fmt)


def pct(v):
    return "%.1f %%" % v


def pct0(v):
    return "%.0f %%" % v


def g3(v):
    return "%.3g" % v


def e2(v):
    return "%.2e" % v


def e4(v):
    return "%.4e" % v


def intc(v):
    return "{:,}".format(int(v))


# residual extremes recomputed from the full sets' own model/data tables (they reproduce FINDINGS A14 where it has them)
def table_resid(c):
    p = setp(c, "model_data/tables/" + CASE[c]["table"])
    key = "tab:" + p
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


def tr_min_ng(c):
    ng, rg, p = table_resid(c)
    v = min(x for a, x in ng if 10.0 <= a <= 100.0)
    return N.add(c, "lowest ln(model/obs) NG, gauge altitude 10-100 m", v, p,
                 "min over rows useNG=1, 10<=altitude<=100 m of ln(pNG/refNG)", "%+.2f")


def tr_max_rg(c):
    ng, rg, p = table_resid(c)
    v = max(x for a, x in rg if a < 20.0)
    return N.add(c, "highest ln(model/obs) RG, gauge altitude < 20 m", v, p,
                 "max over rows useRG=1, altitude<20 m of ln(pRG/refRG)", "%+.2f")


# ------------------------------------------------------------------------------------------------ catalogue
def catalogue():
    G = []

    def add(sec, rel, title, fam, **kw):
        G.append(dict(section=sec, rel=rel, title=title, fam=fam, **kw))

    md = "model_data/"
    add("model_data", md + "publication_NG_RG_altitude/<STEM>_NG_RG_altitude.png",
        "Publication figure: COPS and model NG / RG against altitude", "md_pub")
    add("model_data", md + "<STEM>_vs_COPS_altitude.png",
        "Model vs COPS against altitude: (a) NG, (b) RG, (c) NG and RG together", "md_alt")
    add("model_data", md + "<STEM>_vs_COPS_altitude_NG.png", "Model vs COPS against altitude: NG (panel a alone)", "md_ng")
    add("model_data", md + "<STEM>_vs_COPS_altitude_RG.png", "Model vs COPS against altitude: RG (panel b alone)", "md_rg")
    add("model_data", md + "<STEM>_vs_COPS_altitude_NG_RG.png",
        "Model vs COPS against altitude: NG and RG together, NG = RG crossings (panel c alone)", "md_ngrg")
    add("model_data", md + "<STEM>_vs_COPS_time.png",
        "Model vs COPS against UTC: whole descent and the last 25 minutes", "md_time")
    add("model_data", md + "<STEM>_pressures_altitude.png",
        "Pressures and NG / RG ratio against altitude (tools/pub_figures.py format)", "md_press")
    add("model_data", md + "<STEM>_pressures_time.png",
        "Pressures and NG / RG ratio against UTC (tools/pub_figures.py format)", "md_press")
    for sec, q, tpl, times in (("activity", "activity", "12M_templates_activity", [None]),
                               ("production_rate", "h2o", "12M_templates_h2o_production_rate_%sUTC",
                                ["103310", "103926"])):
        for t in times:
            pre = sec + "/"
            base = tpl % t if t else tpl
            gt = "12M_groundtrack_activity" if not t else "12M_groundtrack_h2o_production_rate_%sUTC" % t
            what = "H2O activity factor A" if not t else "H2O production rate F_H2O at %s:%s:%s UTC" % (t[:2], t[2:4], t[4:])
            add(sec, pre + base + ".png", what + ": maps with effective slope and cliff templates (six panels)",
                "map_tpl", q=q, t=t, panels=[("a", "latlon_global"), ("b", "latlon_zoom")])
            add(sec, pre + gt + "_latlon.png", what + " with the ground track, longitude / latitude",
                "map_latlon", q=q, t=t, panels=[("a", "latlon_global"), ("b", "latlon_zoom")])
            add(sec, pre + gt + "_3d.png", what + " with the ground track, 3-D views", "map_3d", q=q, t=t,
                panels=[("a", "3d_a"), ("b", "3d_b"), ("c", "3d_c")])
            pq = "activity" if not t else "h2o_production_rate"
            add(sec, pre + "panels/" + base + "_a_%s_global.png" % pq, what + ": panel (a) global map alone",
                "map_panel", q=q, t=t, panels=[("a", "latlon_global")])
            add(sec, pre + "panels/" + base + "_b_%s_impact_area.png" % pq, what + ": panel (b) impact area alone",
                "map_panel", q=q, t=t, panels=[("b", "latlon_zoom")])
            for p, nm in (("c_slope_global", "panel (c) effective gravitational slope, global"),
                          ("d_slope_impact_area", "panel (d) effective gravitational slope, impact area"),
                          ("e_cliff_global", "panel (e) cliff templates, global"),
                          ("f_cliff_impact_area", "panel (f) cliff templates, impact area")):
                add(sec, pre + "panels/" + base + "_" + p + ".png", what + ": " + nm, "geom_panel", q=q, t=t)
    sa = "source_areas/12M_source_area_"
    for g in ("NG", "RG", "NG_RG"):
        add("source_areas", sa + "union_%s.png" % g, "Union of the P%% source areas, P = 50-90 %%, %s" % g.replace("_", " and "),
            "sa_union", gauges=g.split("_"))
    for P in (50, 60, 70, 80, 90):
        for g in ("NG", "RG", "NG_RG"):
            add("source_areas", sa + "occurrence_P%d_%s.png" % (P, g),
                "Source-area occurrence N_%d, %s" % (P, g.replace("_", " and ")), "sa_occ", P=P, gauges=g.split("_"))
    cd = "contribution_density/12M_contrib_density_"
    snap = [("103310", "10:33:10 UTC (start of the late window)"), ("103731", "10:37:31.282 UTC"),
            ("103909", "10:39:09.282 UTC"), ("103926", "10:39:26 UTC (last observation)")]
    for t, nm in snap:
        for g in ("NG", "RG", "NG_RG"):
            add("contribution_density", cd + "%s_%s.png" % (t, g),
                "Contribution density at %s, %s" % (nm, g.replace("_", " and ")), "cd", key=t, gauges=g.split("_"))
    for t, nm, gs in (("pre_gap", "mean before the NG gap", ("NG", "RG", "NG_RG")),
                      ("gap", "mean over the gap in the NG data (RG only)", ("RG",)),
                      ("late", "mean over the late window", ("NG", "RG", "NG_RG"))):
        for g in gs:
            add("contribution_density", cd + "%s_%s.png" % (t, g),
                "Contribution density, %s, %s" % (nm, g.replace("_", " and ")), "cd", key=t, gauges=g.split("_"))
    add("cliffs", "cliffs/12M_cliff_slope_bins.png", "Use of the effective-slope bins: area, multipliers, production",
        "cliff_bins")
    for g in ("NG", "RG", "NG_RG"):
        add("cliffs", "cliffs/12M_cliff_pressure_fraction_%s.png" % g,
            "Share of the modelled surface pressure by slope bin, %s" % g.replace("_", " and "), "cliff_pf",
            gauges=g.split("_"))
    for kind, nm in (("lobe", "by lobe"), ("distance", "by distance from the touchdown point")):
        for g in ("NG", "RG", "NG_RG"):
            add("regions", "regions/12M_region_share_%s_%s.png" % (kind, g),
                "Share of the modelled surface pressure %s, %s" % (nm, g.replace("_", " and ")), "reg_" + kind,
                gauges=g.split("_"))
    add("visibility", "visibility/12M_templates_visibility.png",
        "Late-visibility templates, effective slope and cliff templates (fit-independent)", "vis")
    add("visibility", "visibility/12M_groundtrack_visibility_latlon.png",
        "Ground track and impact point with the late-visibility templates, longitude / latitude", "vis")
    add("visibility", "visibility/12M_groundtrack_visibility_3d.png",
        "Ground track and impact point with the late-visibility templates, 3-D views", "vis")
    return G


def gsrc(g):
    """files whose text may hold the literal numbers of a plot's hand-written summary"""
    caps = [caption_path(c, g["rel"]) for c in CK]
    return ([x for x in caps if os.path.exists(x)] + [setp("softplus", "README.txt"), setp(g_case_readme(g), "README.txt"),
            FINDINGS, OUTVIS, STAGE2 + "/summary.txt"])


def g_case_readme(g):
    return "D2" if g["fam"].startswith("md") else "C"


def fig_fallback(g):
    return ("observation from the figure: approximate axis reading or a bin / category label shown in the figure",
            src_path("softplus", g["rel"]) + " and the other five sets' PNG of this figure")


def src_path(c, rel):
    return setp(c, rel.replace("<STEM>", CASE[c]["stem"]))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ images and layout
PAGE = {"P": (8.5, 11.0), "L": (11.0, 8.5)}
MARGIN = 0.6
OVERHEAD = 1.55    # heading, files line, caption (inches)
LABEL = 0.24       # case label above each image
PAD = 0.12         # cell padding


def usable(o):
    w, h = PAGE[o]
    return w - 2 * MARGIN, h - 2 * MARGIN


def layout(aspects, single):
    """best (orientation, cols, rows, image width in inches) for the aspect ratios (h/w) of the images"""
    rmax = max(aspects)
    best = None
    opts = [("P", 1, 1), ("L", 1, 1)] if single else [("P", 3, 2), ("P", 2, 3), ("L", 3, 2), ("L", 2, 3)]
    for o, cols, rows in opts:
        W, H = usable(o)
        wcell = W / cols
        hrow = (H - OVERHEAD - (0.3 if single else 0.0)) / rows - LABEL - 0.05
        w = min(wcell - PAD, hrow / rmax)
        if single:
            w = min(w, 7.0)
        score = w * w * rmax
        if best is None or score > best[4] + 1e-9:
            best = (o, cols, rows, w, score)
    return dict(orient=best[0], cols=best[1], rows=best[2], width_in=round(best[3], 3))


LINE_SECTIONS = ("model_data", "cliffs", "regions")


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


# ------------------------------------------------------------------------------------------------ text: per family
def order_text(vals, fmt=None, reverse=False):
    items = sorted(vals.items(), key=lambda kv: kv[1][0], reverse=reverse)
    return ", ".join("%s %s" % (k, t) for k, (v, t) in items)


def readme_title(c):
    p = setp(c, "README.txt")
    line = tload(p).splitlines()[0]
    return line.replace("Complete figure set of the 12M ", "").strip(), p


def readme_caveats(c):
    p = setp(c, "README.txt")
    out = []
    for line in tload(p).splitlines():
        s = line.strip()
        if s.startswith("* Caveat:"):
            out.append(s[2:])
        if c == "A" and "The last fit stopped at max_nfev" in s:
            out.append(s[2:])
    return out, p


def readme_folder_text(sec):
    p = setp("softplus", "README.txt")
    lines = tload(p).splitlines()
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
    return " ".join(out), p


# ------------------------------------------------------------------------------------------------ docx helpers
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
    try:
        cap_style = st["Caption"]
        cap_style.font.name = "Arial"
        cap_style.font.size = Pt(9)
        cap_style.font.italic = False
        cap_style.font.bold = False
        cap_style.font.color.rgb = RGBColor(0, 0, 0)
    except KeyError:
        cap_style = None

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

    def para(text="", bold_lead=None, italic=False, size=None, style=None, align=None, after=None):
        text = T(text, SRC["src"], SRC["fb"])
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
        if align == "center":
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
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
        """insert in schema order (CT_TblPr: ... shd, tblLayout, tblCellMar, tblLook, tblCaption, tblDescription)"""
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

    def data_table(header, rows, widths, size=8, first_bold=True):
        avail = usable(CUR["orient"])[0] - 0.1
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
                    r = p.add_run(T(str(part), SRC["src"], SRC["fb"]))
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
        r = p.add_run("12M stage-2 full figure sets: comparison of the six cases, %s   |   page " % DOC_DATE)
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

    zoom = doc.settings.element.find(qn("w:zoom"))
    if zoom is not None and zoom.get(qn("w:percent")) is None:
        zoom.set(qn("w:percent"), "100")
    sec0 = doc.sections[0]
    set_orient(sec0, "L")
    add_footer(sec0)

    # ---------------- 1. title, date, purpose
    N.ctx = "1. title / purpose"
    FRONT = [setp(c, "README.txt") for c in CK] + [FINDINGS, OUTVIS, STAGE2 + "/summary.txt",
                                                    setp("D2", "model_data/12M_D2_vs_COPS_altitude_caption.txt")]
    SRC["src"], SRC["fb"] = FRONT, ("document date / build record (generator constants, build_manifest.json)",
                                    os.path.abspath(__file__))
    doc.add_paragraph("12M stage-2 final states: the six full figure sets compared", style="Title")
    para(DOC_DATE + " (US Central). SHAP8 12M shape model, COPS NG / RG data of Rosetta's final descent (2016-09-30), "
         "collisionless (Liouville) forward model.", italic=True)
    para("This document puts the same figure of the six full figure sets of the 12M stage-2 final states side by side "
         "(the regularized stage-2 chains softplus, A, B, D and the best option-C chain, plus the clipped refit D2), so that "
         "each plot can be compared across the activity models at a glance. For every plot it gives what the plot shows "
         "(from its caption), the numbers that differ between the cases (read from the captions, the checks/ tables, "
         "the manifests and the README facts of each set; all listed with their source in numbers_used.json) and how the "
         "cases compare, including what is visible in the images. It adds no new fits or figures and no new physical "
         "interpretation: statements beyond the captions, READMEs and FINDINGS.md are marked as observations from the "
         "figures.", bold_lead="Purpose.")
    para("Case order everywhere: softplus (reference), A, B, D, C, D2. Each grid shows the six PNGs of one figure, each "
         "labelled with its case. Colour classes of the activity, production-rate and contribution-density maps are chosen "
         "per panel and per case (captions), so equal colours do not mean equal values across cases; source-area "
         "occurrence classes, slope-bin and region shares use the same scale in every case. Images are downscaled copies "
         "(%s); the originals (300 dpi PNG and EPS) are listed in Appendix A." % args.img_note,
         bold_lead="How to read it.")

    # ---------------- 2. cases table
    N.ctx = "2. cases table"
    SRC["fb"] = None
    doc.add_paragraph("1. Cases", style="Heading 1")
    para("Final state of each case = the last cycle's fine tuning (c5/ft/final; D2: its fit, no fine tuning). Data cost, "
         "terminations: summary.csv of model_vs_observations_nonneg_options_stage2_final; rms ln(model/obs): "
         "12M_stage2_reg_best_residuals_check.json there (all fitted rows; late window W10x5 in brackets); Q_H2O, seen "
         "share, kappa, multipliers: checks/facts.json of each full set (seen Q = Q_H2O x seen share). Reaching max_nfev "
         "is not convergence.", size=8.5)
    hdr = ["case", "model (README title)", "final data cost", "last fit / last fine tuning: termination (nfev)",
           "rms ln NG all (late)", "rms ln RG all (late)", "Q_H2O 10:33:10 total [s^-1]", "Q_H2O seen [s^-1]",
           "seen share", "kappa", "slope multipliers 40-50 ... 90-180 deg", "key caveats"]
    rows, cav_list = [], []
    for c in CK:
        title, tp = readme_title(c)
        N.add(c, "README title", title, tp, "first line", "%s")
        r, _ = csvrow(c)
        term = cx(c, "fit_termination", "last fit termination", "%s", conv=None)
        term_s = "max_nfev" if "maximum number" in term else ("ftol" if "ftol" in term else ("xtol" if "xtol" in term else term))
        nfev = cx(c, "fit_nfev", "last fit nfev", "%d", conv=lambda s: int(float(s)))
        seg = r.get("fit_segments") or ""
        segtxt = (", %s segments" % cx(c, "fit_segments", "last fit segments", "%d", conv=lambda s: int(float(s)))
                  if seg and int(float(seg)) > 1 else "")
        fit_t = "fit %s (%s%s)" % (term_s, nfev, segtxt)
        if r.get("ft_termination"):
            ft_t = "; fine tuning %s (%s)" % (cx(c, "ft_termination", "last fine tuning termination", "%s", conv=None),
                                             cx(c, "ft_nfev", "last fine tuning nfev", "%d", conv=lambda s: int(float(s))))
        else:
            ft_t = "; no fine tuning"
        cost = cx(c, "data_cost", "final data cost", e4)
        rng = rx(c, "NG", "rms", "rms ln NG (all fitted rows)")
        rngl = rx(c, "NG", "rms_late", "rms ln NG (late window)")
        rrg = rx(c, "RG", "rms", "rms ln RG (all fitted rows)")
        rrgl = rx(c, "RG", "rms_late", "rms ln RG (late window)")
        Q = fx(c, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"], "Q_H2O total at 10:33:10 [s^-1]", e2)
        qv = fxv(c, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"])
        sv = fxv(c, ["production_seen_unseen", "10:33:10", "seen_share"])
        qs = N.add(c, "Q_H2O of the seen facets at 10:33:10 [s^-1] (= Q_H2O x seen_share)", qv * sv,
                   setp(c, "checks/facts.json"),
                   "facts.json: production_seen_unseen/10:33:10/Q_H2O_from_vtu x .../seen_share", e2)
        ss = fx(c, ["production_seen_unseen", "10:33:10", "seen_share"], "seen share of Q_H2O at 10:33:10", pct, 100.0)
        kap = fx(c, ["kappa"], "kappa (CO2 scale)", "%.2g")
        mult = fxv(c, ["multipliers"])
        if mult:
            ms = ", ".join(fx(c, ["multipliers", i, 1], "multiplier " + mult[i][0], "%.4g") for i in range(len(mult)))
        else:
            ms = N.add(c, "slope multipliers", "none", setp(c, "README.txt"),
                       "Facts: 'No slope multipliers (all bins 1)'", "%s") + " (all bins 1)"
        cav, cp = readme_caveats(c)
        refs = []
        flags = []
        if term_s == "max_nfev":
            flags.append("last fit at max_nfev: not converged")
        for cv in cav:
            cav_list.append((c, cv, cp))
            refs.append("[%d]" % len(cav_list))
        rows.append([CASE[c]["label"], title, cost, fit_t + ft_t, "%s (%s)" % (rng, rngl), "%s (%s)" % (rrg, rrgl), Q, qs,
                     ss, kap, ms, "; ".join(flags + ["README " + " ".join(refs)] if refs else flags) or "-"])
    data_table(hdr, rows, [0.6, 1.25, 0.62, 1.05, 0.62, 0.62, 0.65, 0.62, 0.5, 0.55, 1.0, 1.25], size=7)
    para("Caveats quoted from the READMEs:", size=8.5, bold_lead="")
    for i, (c, cv, cp) in enumerate(cav_list):
        N.add(c, "README caveat [%d]" % (i + 1), cv, cp, "README.txt line", "%s")
        bullet(cv, lead="[%d] %s:" % (i + 1, c), size=8)
    para("Parameters at bounds in the last fit (summary.csv 'at_bounds'): " + "; ".join(
        "%s: %s" % (c, cx(c, "at_bounds", "parameters at bounds (last fit)", "%s", conv=None) or "(none listed)")
        for c in CK) + ".", size=7.5)

    # ---------------- 3. overall summary
    N.ctx = "3. overall summary"
    new_section("L")
    doc.add_paragraph("2. Overall summary", style="Heading 1")
    costs = {c: (float(csvrow(c)[0]["data_cost"]), cx(c, "data_cost", "final data cost", e4)) for c in CK}
    cC, cS = costs["C"][0], costs["softplus"][0]
    rel = N.add("C vs softplus", "relative data cost C / softplus - 1", 100.0 * (cC / cS - 1.0), STAGE2 + "/summary.csv",
                "computed: data_cost(C_f1_pm1em4_reg)/data_cost(softplus_reg) - 1", "%.0f %%")
    maxn = [c for c in CK if "maximum number" in csvrow(c)[0]["fit_termination"]]
    ftol = [c for c in CK if "ftol" in csvrow(c)[0]["fit_termination"]]
    for c in maxn + ftol:
        cx(c, "fit_termination", "last fit termination", "%s", conv=None)
    para("Ordered by final data cost: " + order_text(costs) + ". C has the lowest data cost (%s relative to the regularized "
         "softplus; FINDINGS A14: 47 %% below). The last fits of %s stopped at max_nfev (A in its second segment after a "
         "time limit), those of %s on ftol; every last fine tuning stopped on ftol. Reaching max_nfev is not "
         "convergence, so the ranking is that of non-converged states for %s." % (rel, ", ".join(maxn), ", ".join(ftol),
                                                                               ", ".join(maxn)),
         bold_lead="Fit quality.")
    kap = {c: (fxv(c, ["kappa"]), fx(c, ["kappa"], "kappa (CO2 scale)", "%.2g")) for c in CK}
    qs = {}
    for c in CK:
        qv = fxv(c, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"])
        sv = fxv(c, ["production_seen_unseen", "10:33:10", "seen_share"])
        qs[c] = (qv * sv, N.add(c, "Q_H2O of the seen facets at 10:33:10 [s^-1] (= Q_H2O x seen_share)", qv * sv,
                                setp(c, "checks/facts.json"),
                                "facts.json: production_seen_unseen/10:33:10/Q_H2O_from_vtu x .../seen_share", e2))
    qt = {c: (fxv(c, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"]),
              fx(c, ["production_seen_unseen", "10:33:10", "Q_H2O_from_vtu"], "Q_H2O total at 10:33:10 [s^-1]", e2))
          for c in CK}
    us = {c: (fxv(c, ["production_seen_unseen", "10:33:10", "unseen_share"]),
              fx(c, ["production_seen_unseen", "10:33:10", "unseen_share"], "unseen share of Q_H2O at 10:33:10", pct, 100.0))
          for c in CK}
    big = {c: (fxv(c, ["regions", "big_lobe_NG_median_pre"]),
               fx(c, ["regions", "big_lobe_NG_median_pre"], "big lobe NG median share before the late window (facts)",
                  pct0, 100.0)) for c in CK}
    lastng = {c: (fxv(c, ["pressure_share_by_bin", "NG", "last"])[6],
                  fx(c, ["pressure_share_by_bin", "NG", "last", 6], "90-180 deg share of NG at the last observation",
                     pct0, 100.0)) for c in CK}
    lastrg = {c: (fxv(c, ["pressure_share_by_bin", "RG", "last"])[6],
                  fx(c, ["pressure_share_by_bin", "RG", "last", 6], "90-180 deg share of RG at the last observation",
                     pct0, 100.0)) for c in CK}
    para("kappa (the CO2 scale on the srcCG CO2 map, lower bound 0) is at its lower bound in all six (%s); the gauges "
         "measure total pressure, so this says the fits prefer the H2O law for the whole signal, not that CO2 is measured "
         "absent (A14). The H2O production of the facets seen by a fitted observation is similar in all six (%s s^-1 at "
         "10:33:10), while the total Q_H2O spans %s (%s); %s of Q_H2O comes from facets no fitted observation sees. "
         "Before the late window the big lobe gives a median %s of the modelled NG. At the last observation (5.5 m) the "
         "90-180 deg (overhang) bin gives %s of the modelled NG but %s of RG. The visibility figures are identical in "
         "all six (checksums)." % (
             order_text(kap), order_text(qs),
             "%s ... %s" % (min(qt.values())[1], max(qt.values())[1]), order_text(qt), order_text(us), order_text(big),
             order_text(lastng), order_text(lastrg)),
         bold_lead="Common to all six (FINDINGS A15).")
    az = {c: (fxv(c, ["A_zero", "seen_area"]), fx(c, ["A_zero", "seen_area"], "share of the seen area with A = 0",
                                                  pct, 100.0)) for c in CK}
    st60 = {}
    for c in CK:
        d = fxv(c, ["production_share_seen_103310"])
        v = sum(d[b] for b in BINS[3:])
        st60[c] = (v, N.add(c, "share of the seen H2O production at 10:33:10 from slope >= 60 deg", 100.0 * v,
                            setp(c, "checks/facts.json"),
                            "facts.json: sum of production_share_seen_103310 [60,70) ... [90,180]", pct))
    p50 = {c: (fxv(c, ["source_area_unions", "NG_P50", "area_m2"]),
               fx(c, ["source_area_unions", "NG_P50", "area_m2"], "P50 NG source-area union area [m^2]", e2)) for c in CK}
    sm = {c: (fxv(c, ["regions", "small_lobe_median_late"])[0],
              fx(c, ["regions", "small_lobe_median_late", 0], "small lobe NG median share in the late window", pct0, 100.0))
          for c in CK}
    smr = {c: (fxv(c, ["regions", "small_lobe_median_late"])[1],
               fx(c, ["regions", "small_lobe_median_late", 1], "small lobe RG median share in the late window", pct0, 100.0))
           for c in CK}
    para("Where the activity sits differs strongly. A = 0 on %s of the seen area. The share of the seen H2O production "
         "at 10:33:10 from effective slopes >= 60 deg (14.6 %% of the seen area in every case) is %s. The P50 source-area "
         "union of the 95 post-gap NG observations covers %s m^2. In the late window the small lobe (head) gives a median "
         "%s of the modelled NG and %s of RG. A15: the data constrain the integrated signal, not where on the near "
         "surface it comes from." % (order_text(az), order_text(st60, reverse=True), order_text(p50), order_text(sm),
                                     order_text(smr)), bold_lead="What differs.")
    m78 = fx("C", ["multipliers", 3, 1], "multiplier c_slope_70_80", "%.3g")
    m89 = fx("C", ["multipliers", 4, 1], "multiplier c_slope_80_90", "%.3g")
    bullet("C's steep-bin multipliers (70-80 deg %s, 80-90 deg %s) are not determined (the level of the c-field trades "
           "against them, A12) and its Q_H2O (%s s^-1, %s from unseen facets) is driven by unseen steep terrain, not by "
           "the data." % (m78, m89, qt["C"][1], us["C"][1]), lead="Caveats.")
    bullet("Q_H2O total is not set by the COPS data (A15): it depends on how each model extends the activity to terrain "
           "no fitted observation sees; the seen production is what the data constrain.")
    bullet("Several last fits stopped at max_nfev (%s): not converged. kappa at its lower bound in every fit." % ", ".join(maxn))
    bullet("D2 is the clipped linear model without slope multipliers and without fine tuning; its model_data figures "
           "show the D2 curve on top of the identical clipped reference (README), and its slope-bin multiplier panel "
           "shows 1 in every bin.")
    bullet("Observation from the figures: in D2's model_data figures (vs_COPS_altitude, _NG, _RG, _NG_RG, _time) the "
           "legend entry of the blue curve reads 'data cost' without the number (the captions give 4.8864e-03); a label "
           "defect of the D2 figures, not a difference in the data.")
    para("FINDINGS.md A12, A14, A15 (%s) and the README Facts of each set are the sources of the statements above; this "
         "document only arranges them next to the figures." % FINDINGS, size=8.5, italic=True)

    # ---------------- 4. per section
    fig_no = 0
    for si, sec in enumerate(SECTIONS):
        N.ctx = "4.%d %s (intro)" % (si + 1, sec)
        gs = [g for g in groups if g["section"] == sec]
        SRC["src"] = ([caption_path("softplus", g["rel"]) for g in gs if os.path.exists(caption_path("softplus", g["rel"]))]
                      + FRONT)
        SRC["fb"] = None
        new_section("L")
        doc.add_paragraph("3.%d %s/" % (si + 1, sec), style="Heading 1")
        ftxt, fp = readme_folder_text(sec)
        N.add("all", "README folder description " + sec, ftxt, fp, "README.txt folder block", "%s")
        para(ftxt + " (softplus set's README; the same layout in all six sets with their own names)",
             bold_lead="Folder (README.txt):")
        para(SECTION_INTRO[sec], bold_lead="About the plots in this folder.")
        lines = []
        for g in gs:
            lines.append("Figure %d: %s%s" % (fig_no + 1 + len(lines), g["title"],
                                              " (identical in all six)" if g["identical"] else ""))
        for ln in lines:
            bullet(ln, size=8.5)
        for gi, g in enumerate(gs):
            fig_no += 1
            g["fig"] = fig_no
            N.ctx = "Figure %d (%s)" % (fig_no, g["rel"])
            SRC["src"] = gsrc(g)
            L = g["layout"]
            new_section(L["orient"])
            doc.add_paragraph("3.%d.%d %s" % (si + 1, gi + 1, g["title"]), style="Heading 2")
            pf = para("Files: <set>/%s%s" % (g["rel"], "" if g["fam"] == "md_pub" else "  (caption: same name, _caption.txt)"),
                      size=7.5, after=2)
            show = ["softplus"] if g["identical"] else CK
            cols = L["cols"] if not g["identical"] else 1
            nrow = math.ceil(len(show) / cols)
            W, H = usable(L["orient"])
            t = doc.add_table(rows=nrow, cols=cols)
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for k, c in enumerate(show):
                cell = t.rows[k // cols].cells[k % cols]
                cell.text = ""
                p0 = cell.paragraphs[0]
                p0.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p0.paragraph_format.space_after = Pt(1)
                p0.paragraph_format.keep_with_next = True
                lab = ("all six cases (identical files)" if g["identical"] else "%s" % CASE[c]["label"])
                r = p0.add_run(lab)
                r.bold = True
                r.font.size = Pt(9)
                r.font.color.rgb = RGBColor(0xA0, 0x10, 0x10)
                p1 = cell.add_paragraph()
                p1.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p1.paragraph_format.space_after = Pt(2)
                p1.add_run().add_picture(img[(g["rel"], c)]["dst"], width=Inches(L["width_in"]))
            fixed(t, [W / cols] * cols)
            cell_margins(t, 40)
            cp = doc.add_paragraph(style=cap_style) if cap_style is not None else doc.add_paragraph()
            r = cp.add_run("Figure %d: %s \u2014 %s." % (fig_no, g["title"],
                                                  "identical in all six cases (one image shown)" if g["identical"]
                                                  else "all six cases (softplus, A, B, D, C, D2)"))
            r.font.size = Pt(9)
            r.bold = True
            # summary
            doc.add_paragraph("Summary", style="Heading 3")
            summarize(g, doc, para, bullet, data_table)
    # ---------------- appendix A
    N.ctx = "Appendix A"
    SRC["src"], SRC["fb"] = [], ("file path or checksum prefix (Appendix A)", ROOT)
    new_section("L")
    doc.add_paragraph("Appendix A. Source files per figure", style="Heading 1")
    para("Paths relative to %s/ (originals: .png as listed and .eps with the same name; captions <name>_caption.txt). "
         "<set> and <STEM> per case:" % ROOT, size=8.5)
    data_table(["case", "<set>", "<STEM>"], [[CASE[c]["label"], CASE[c]["set"], CASE[c]["stem"]] for c in CK],
               [2.0, 3.6, 3.6], size=8)
    rows = []
    for g in groups:
        paths = []
        for c in CK:
            p = src_path(c, g["rel"])
            paths.append("%s: %s (sha256 %s...)" % (c, os.path.relpath(p, ROOT), g["sha"][c][:12]))
        rows.append(["Figure %d" % g["fig"], [g["title"] + (" [identical in all six]" if g["identical"] else "")] + paths])
    data_table(["figure", "plot and the six source PNGs"], rows, [0.8, 8.8], size=6.5)
    # ---------------- appendix B
    N.ctx = "Appendix B"
    SRC["src"], SRC["fb"] = [os.path.abspath(__file__)], ("build parameter / build record", os.path.abspath(__file__))
    doc.add_paragraph("Appendix B. How this document was built", style="Heading 1")
    for txt in BUILD_NOTES(args):
        bullet(txt, size=8.5)
    if N.unverified:
        with open(os.path.join(out, "numbers_unverified.json"), "w") as f:
            json.dump(N.unverified, f, indent=1, ensure_ascii=False)
        sys.exit("%d numeric tokens of hand-written text not found in a source file (numbers_unverified.json); "
                 "document not written" % len(N.unverified))
    doc.save(os.path.join(out, DOC_NAME + ".part"))
    os.replace(os.path.join(out, DOC_NAME + ".part"), os.path.join(out, DOC_NAME))


# ------------------------------------------------------------------------------------------------ section texts
SECTION_INTRO = {
    "model_data": "COPS NG and RG pressures against altitude or UTC with the model of each case (blue) and, in every "
                  "set, the clipped refit D2 as the reference (orange dashed). Model at the observation times; for the 95 "
                  "observations from 10:33:10 UTC the gauges' time response (10-s running mean, 5-readout average) and "
                  "the fitted grain term in NG, as in the fit. Late window 10:33:10-10:39:26 UTC shaded (not in the "
                  "publication figure).",
    "activity": "The fitted H2O activity factor A of each final state on the seen facets (log classes per panel; light "
                "green = seen with A = 0; grey = not seen during the fit), with the fit-independent effective "
                "gravitational slope and cliff templates. A is unitless in the code (F_H2O = q_H2O A carries the units).",
    "production_rate": "The H2O production rate F_H2O = q_H2O A per facet [m^-2 s^-1] at 10:33:10 and 10:39:26 UTC "
                       "(start and end of the late window), as the model computes it (day / night / cast-shadow rule), "
                       "with relief shading under the colours; the same layout as the activity maps.",
    "source_areas": "P% source areas of the 95 post-gap observations: for each observation and gauge the smallest set of "
                    "facets that supplies P % of the modelled surface pressure (facets in decreasing order of their "
                    "contribution). Union: the facets in S_P of at least one observation; occurrence N_P: in how many of "
                    "the 95. The close-up box is fitted per case to its own source areas, so the close-ups differ in "
                    "extent between the cases.",
    "contribution_density": "Contribution of each facet per unit surface area to the modelled NG density and RG flux, "
                            "at four snapshots and as means over three phases, with the 50 % / 90 % source-area outlines. "
                            "Log colour range per map (at most 3 decades below its 99.99th percentile); facets on "
                            "overhangs or under other surface layers are not drawn in the lon/lat maps (the share of the "
                            "signal drawn in the close-up is given per map).",
    "cliffs": "How each state uses the effective-gravitational-slope bins of its multipliers (area shares, multipliers, "
              "production shares) and which share of the modelled surface pressure comes from each bin along the descent.",
    "regions": "Share of the modelled surface pressure from the small lobe (head, with the touchdown point) and the big "
               "lobe (body), and from distance shells around the touchdown point, along the descent (NG grain term not "
               "included).",
    "visibility": "Fit-independent figures on the observation-only base: the late-visibility templates (NG-only / RG-only "
                  "facets), effective slope, cliff templates and the ground track with the Rosetta impact point. The "
                  "PNGs are byte-identical in the six sets (SHA256 checked here) and pixel-identical to the c9 set's "
                  "(CAPTIONS_COPIED.txt).",
}


def BUILD_NOTES(args):
    return [
        "Generator: %s (run copy; version-controlled copy model_ext/figures/fullsets_comparison_docx.py, SHA256 in "
        "model_ext/RUN_COPIES.sha256), run with %s (python-docx %s, Pillow %s; separate venv, not the project .venv)."
        % (os.path.abspath(__file__), sys.executable, args.docx_version, args.pil_version),
        "Command: %s; Slurm job %s on %s; %s UTC." % (" ".join(sys.argv), os.environ.get("SLURM_JOB_ID", "-"),
                                                   os.uname().nodename, args.started),
        "Images: each PNG of the six sets (85 per set) hashed (SHA256); a plot whose six PNGs are identical is shown once. "
        "Downscaled copies in images/: maps as JPEG (quality 78), line plots (model_data, cliffs, regions) as PNG with 256 "
        "colours; width = min(native, max(1000 px, 300 dpi at the printed width), 1400 px); RGBA flattened on white.",
        "Layout: per plot the page orientation and grid (3 x 2 or 2 x 3, portrait or landscape) that gives the largest "
        "images within the page (US Letter, margins 0.6 in), from the largest height/width ratio of the six PNGs.",
        "Numbers: every number in the text and tables is read from a file at build time (facts.json, summary.csv, "
        "residuals_check.json, captions by regular expression, manifests, source_area_tables.txt, slope_bins_*.json, "
        "the model/data tables) and listed in numbers_used.json with its source, locator and where it is used. Literal "
        "numbers of the hand-written sentences are checked against the plot's captions, the READMEs, FINDINGS.md, "
        "OUTPUT_VISUALIZATION.md and summary.txt (the build refuses to write the document if one is not found); "
        "numbers in the 'Observation from the figures' paragraphs are recorded as figure readings. Residual "
        "extremes are recomputed from model_data/tables/*_final.dat (ln(model/obs), rows used in the fit) and reproduce "
        "FINDINGS A14 where A14 gives them.",
        "Validation (validate step): the document reopened with python-docx, figures, captions, sections and images "
        "counted, image sizes checked against the pages; validation/validation.json. LibreOffice (soffice) is not "
        "installed on the login node or the compute nodes (checked 2026-11-17; not on conda-forge either), so no PDF was "
        "rendered; validation/ has layout previews drawn from the same geometry instead.",
    ]


# ------------------------------------------------------------------------------------------------ summaries per plot
def summarize(g, doc, para, bullet, data_table):
    fam = g["fam"]
    rel = g["rel"]
    F = FAMILY[fam]
    para(T(F["about"](g) if callable(F["about"]) else F["about"], gsrc(g)), bold_lead="What the plot shows.")
    if g["identical"]:
        para("The six PNGs are byte-identical (SHA256 %s...): the figure does not depend on the fit. One image is shown."
             % N.add("all", "SHA256 of the six identical PNGs", g["sha"]["softplus"], src_path("softplus", g["rel"]),
                     "sha256 of the file (equal for the six sets)", lambda v: v[:16]), bold_lead="Comparison.")
        extra = F.get("identical_note")
        if extra:
            para(T(extra(g) if callable(extra) else extra, gsrc(g)), bold_lead="Note.")
        return
    hdr, rows, widths = F["numbers"](g)
    if hdr:
        data_table(hdr, rows, widths, size=7.5)
    comp = F["compare"](g)
    if comp:
        para(T(comp, gsrc(g)), bold_lead="How the cases compare.")
    obs = F.get("observe")
    if obs:
        o = obs(g) if callable(obs) else obs
        if o:
            para(T(o, [], fallback=fig_fallback(g), force=True), bold_lead="Observation from the figures:", italic=False)


# model_data -------------------------------------------------------------------------------------------------------
def md_numbers(g):
    fam = g["fam"]
    cols = ["case", "final data cost"]
    rows = []
    for c in CK:
        r = [CASE[c]["label"], cx(c, "data_cost", "final data cost", e4)]
        if fam in ("md_pub", "md_alt", "md_time", "md_ng", "md_press"):
            r.append("%s (%s)" % (rx(c, "NG", "rms", "rms ln NG (all fitted rows)"),
                                  rx(c, "NG", "rms_late", "rms ln NG (late window)")))
        if fam in ("md_pub", "md_alt", "md_time", "md_rg", "md_press"):
            r.append("%s (%s)" % (rx(c, "RG", "rms", "rms ln RG (all fitted rows)"),
                                  rx(c, "RG", "rms_late", "rms ln RG (late window)")))
        if fam in ("md_pub", "md_ng", "md_press", "md_alt", "md_time"):
            r.append(tr_min_ng(c))
        if fam in ("md_pub", "md_rg", "md_press", "md_alt", "md_time"):
            r.append(tr_max_rg(c))
        if fam in ("md_alt", "md_time", "md_ngrg"):
            r.append(capx(c, g["rel"], r"the fitted model at (\d+) m", "model NG = RG crossing altitude [m]") + " m")
        rows.append(r)
    if fam in ("md_pub", "md_alt", "md_time", "md_ng", "md_press"):
        cols.append("rms ln NG all (late)")
    if fam in ("md_pub", "md_alt", "md_time", "md_rg", "md_press"):
        cols.append("rms ln RG all (late)")
    if fam in ("md_pub", "md_ng", "md_press", "md_alt", "md_time"):
        cols.append("lowest ln(model/obs) NG, 10-100 m")
    if fam in ("md_pub", "md_rg", "md_press", "md_alt", "md_time"):
        cols.append("highest ln(model/obs) RG, < 20 m")
    if fam in ("md_alt", "md_time", "md_ngrg"):
        cols.append("model NG = RG crossing (caption)")
    w = [1.6] + [1.25] * (len(cols) - 1)
    return cols, rows, w


def md_about(g):
    fam = g["fam"]
    if fam == "md_pub":
        return ("COPS NG (black) and RG (grey) and the model NG (blue) and RG (orange dashed) against the altitude above "
                "the surface; the gap in the NG data (700-297 m, 10:24:38-10:33:10 UTC) hatched; no late-window shading "
                "(publication format of pub_ng_rg_altitude.py, user decision 2026-11-12). The figure has no _caption.txt in "
                "the full sets; it is the copy of publication_stage2_<key>/ (publication_D2/ for D2), README.txt.")
    c = "softplus"
    t = tload(caption_path(c, g["rel"]))
    first = re.split(r"(?<=\))\. |(?<=[a-z])\. ", t, maxsplit=1)[0]
    first = first.replace(CASE["softplus"]["stem"], "<STEM>")
    return (first.replace("the softplus activity model with slope multipliers >= 1, regularized stage-2 chain, final "
                          "state after 5 fit and fine-tuning cycles", "the case's model") +
            ". Blue: the case's final state; orange dashed: the clipped refit D2 (reference, data cost 4.8864e-03) in "
            "every set. In D2's own figure both curves are D2 and coincide (README).")


def md_compare(g):
    fam = g["fam"]
    costs = {c: (float(csvrow(c)[0]["data_cost"]), cx(c, "data_cost", "final data cost", e4)) for c in CK}
    out = []
    if fam in ("md_pub", "md_alt", "md_time", "md_ng", "md_press"):
        lat = {c: (resid_entry(c)[1]["NG"]["rms_late"], rx(c, "NG", "rms_late", "rms ln NG (late window)")) for c in CK}
        out.append("Late-window rms ln NG from best to worst: %s." % order_text(lat))
    if fam in ("md_pub", "md_alt", "md_time", "md_rg", "md_press"):
        lat = {c: (resid_entry(c)[1]["RG"]["rms_late"], rx(c, "RG", "rms_late", "rms ln RG (late window)")) for c in CK}
        out.append("Late-window rms ln RG: %s." % order_text(lat))
    if fam in ("md_pub", "md_ng", "md_press", "md_alt", "md_time"):
        v = {}
        for c in CK:
            ng, rg, p = table_resid(c)
            v[c] = (min(x for a, x in ng if 10.0 <= a <= 100.0), tr_min_ng(c))
        if all(x[0] < 0 for x in v.values()):
            out.append("Every case is below the NG data somewhere between 100 and 10 m; the lowest ln(model/obs) there: "
                       "%s (FINDINGS A14, for the regularized chains: the capped, penalized late facets can no longer "
                       "supply the steep NG rise)." % order_text(v, reverse=True))
        else:
            out.append("Lowest ln(model/obs) NG between 100 and 10 m: %s." % order_text(v, reverse=True))
    if fam in ("md_pub", "md_rg", "md_press", "md_alt", "md_time"):
        v = {}
        for c in CK:
            ng, rg, p = table_resid(c)
            v[c] = (max(x for a, x in rg if a < 20.0), tr_max_rg(c))
        out.append("Highest ln(model/obs) RG below 20 m: %s." % order_text(v))
    if fam in ("md_alt", "md_time", "md_ngrg"):
        cr = {c: (float(capx(c, g["rel"], r"the fitted model at (\d+) m", "model NG = RG crossing altitude [m]")),
                  capx(c, g["rel"], r"the fitted model at (\d+) m", "model NG = RG crossing altitude [m]") + " m")
              for c in CK}
        dat = capx("softplus", g["rel"], r"crossings: data between (\d+ m and \d+ m)", "data NG = RG crossing interval")
        out.append("Model NG = RG crossing (open circles): %s; the data cross between %s altitude (the NG gap); the D2 "
                   "reference crosses at %s." % (order_text(cr, reverse=True), dat,
                                                 capx("softplus", g["rel"], r"clipped reference model at (\d+) m",
                                                      "D2 reference crossing [m]") + " m"))
    if fam == "md_ngrg":
        out.append("Final data cost: %s." % order_text(costs))
    return " ".join(out)


def md_observe(g):
    fam = g["fam"]
    base = ("At the scale of the figure the six model curves follow the data alike over most of the descent; the "
            "visible differences are in the last ~300 m (late window). ")
    if fam == "md_pub":
        return base + ("A and D put the model RG clearly above the RG data in the last few tens of metres; C's model NG "
                       "stays closest to the NG data in the late window; D2's model NG ends below the last NG points while "
                       "the other five rise to them at the very end.")
    if fam in ("md_alt", "md_time"):
        return base + ("D2's own figure shows one blue curve (the orange reference lies under it) and the legend entry "
                       "of the blue curve has no data-cost number (label defect; the caption gives it). In A and D the "
                       "model RG (panel b) rises above the data below a few tens of metres; just before the NG gap the "
                       "model NG of A, B and D lies below the reference curve.")
    if fam == "md_ng":
        return base + "C's NG curve is the closest to the late NG data; A and D the farthest below them between ~100 and ~10 m."
    if fam == "md_rg":
        return base + "A and D overshoot the RG data near touchdown; softplus, B, C and D2 stay closer to the RG data."
    if fam == "md_ngrg":
        return ("Crossing markers: A, D and C cross at the lower edge of the hatched gap, softplus and B inside it, D2 "
                "at the reference value (both markers coincide in D2's figure). D2's legend lacks the data-cost number.")
    if fam == "md_press":
        return base + ("In the NG / RG ratio panel A and D drop below the data ratio in the late window (a dip at about "
                       "10-20 m) and lie below it just before the gap; softplus, B and C follow the data ratio more "
                       "closely; D2 is below the data ratio in the last ~100 m.")
    return ""


# maps (activity / production) ---------------------------------------------------------------------------------------
def man_name(g):
    return "checks/activity_manifest.json" if g["q"] == "activity" else "checks/h2o_%s_manifest.json" % g["t"]


PANEL_NAME = {"latlon_global": "lon/lat global", "latlon_zoom": "lon/lat impact area", "3d_a": "3-D (a) nucleus",
              "3d_b": "3-D (b) impact region", "3d_c": "3-D (c) late descent"}


def cls_range(c, g, panel):
    p = setp(c, man_name(g))
    ct = jload(p)["class_tables"][panel]
    lo, hi, n = ct["classes"][0]["log10_from"], ct["classes"][-1]["log10_to"], len(ct["classes"])
    s = N.add(c, "colour classes %s %s: log10 range and count" % (g["q"] + (g["t"] or ""), panel),
              [lo, hi, n], p, "class_tables/%s/classes[0].log10_from, [-1].log10_to, len" % panel,
              lambda v: "10^%g-10^%g (%d)" % tuple(v))
    return s, hi - lo


def map_numbers(g):
    hdr = ["case"]
    rows = []
    if g["q"] == "activity":
        hdr += ["A = 0 on seen area", "seen facets A > 0"]
    else:
        hdr += ["Q_H2O total [s^-1] (caption)", "seen share of Q_H2O"]
    for k, panel in g["panels"]:
        hdr.append("classes %s: range (n)" % PANEL_NAME[panel])
    if g["fam"] == "map_tpl":
        hdr.append("multipliers 40-50 ... 90-180 deg")
    for c in CK:
        r = [CASE[c]["label"]]
        if g["q"] == "activity":
            r.append(fx(c, ["A_zero", "seen_area"], "share of the seen area with A = 0", pct, 100.0))
            r.append(fx(c, ["A_zero", "nonzero_seen_facets"], "seen facets with A > 0", intc))
        else:
            tt = "%s:%s:%s" % (g["t"][:2], g["t"][2:4], g["t"][4:])
            relcap = "production_rate/12M_templates_h2o_production_rate_%sUTC.png" % g["t"]
            r.append(capx(c, relcap, r"total H2O production Q_H2O = sum F_H2O area = (\S+) s\^-1",
                          "Q_H2O total at %s [s^-1]" % tt))
            r.append(fx(c, ["production_seen_unseen", tt, "seen_share"], "seen share of Q_H2O at " + tt, pct, 100.0))
        for k, panel in g["panels"]:
            r.append(cls_range(c, g, panel)[0])
        if g["fam"] == "map_tpl":
            mult = fxv(c, ["multipliers"])
            r.append(", ".join(fx(c, ["multipliers", i, 1], "multiplier " + mult[i][0], "%.3g")
                               for i in range(len(mult))) if mult else "none (all 1)")
        rows.append(r)
    w = [1.5] + [(8.3 / (len(hdr) - 1))] * (len(hdr) - 1)
    return hdr, rows, w


MODEL_SP = ("the softplus activity model with slope multipliers >= 1, regularized stage-2 chain, final state after 5 "
            "fit and fine-tuning cycles")


def map_about(g):
    c = "softplus"
    t = tload(caption_path(c, g["rel"]))
    first = t.split(". ")[0].replace(MODEL_SP + " global", "the case's model: global").replace(
        MODEL_SP + " impact", "the case's model: impact").replace(MODEL_SP, "the case's model")
    if "groundtrack" in g["rel"]:
        m = re.search(r"Surface: (.*?)\. Colour only", t)
        if m:
            first += ". Surface: " + m.group(1).replace(MODEL_SP, "the case's model")
    pat = r"'seen': the facet contributes to at least one fitted NG or RG observation; ([\d,]+ of the [\d,]+ facets, [\d.]+ % of the area)"
    seen = {c2: capx(c2, g["rel"], pat, "seen facets during the fit (caption)") for c2 in CK}
    same = len(set(seen.values())) == 1
    s = first + ". Colour only on the facets seen during the fit (%s%s); light green: seen with %s = 0; grey: not " \
                "seen. Classes scanned per panel and per case (captions; checks/*_class_tables.txt), so the colour " \
                "scales differ between the cases." % (seen["softplus"], ", the same in all six" if same else
                                                      "; differs between the cases", "A" if g["q"] == "activity" else "F_H2O")
    return s


def map_compare(g):
    out = []
    if g["q"] == "activity":
        az = {c: (fxv(c, ["A_zero", "seen_area"]), fx(c, ["A_zero", "seen_area"], "share of the seen area with A = 0",
                                                      pct, 100.0)) for c in CK}
        out.append("Share of the seen area with A = 0 (light green): %s." % order_text(az, reverse=True))
    else:
        tt = "%s:%s:%s" % (g["t"][:2], g["t"][2:4], g["t"][4:])
        relcap = "production_rate/12M_templates_h2o_production_rate_%sUTC.png" % g["t"]
        q = {c: (float(capx(c, relcap, r"total H2O production Q_H2O = sum F_H2O area = (\S+) s\^-1",
                            "Q_H2O total at %s [s^-1]" % tt)),
                 capx(c, relcap, r"total H2O production Q_H2O = sum F_H2O area = (\S+) s\^-1",
                      "Q_H2O total at %s [s^-1]" % tt)) for c in CK}
        ss = {c: (fxv(c, ["production_seen_unseen", tt, "seen_share"]),
                  fx(c, ["production_seen_unseen", tt, "seen_share"], "seen share of Q_H2O at " + tt, pct, 100.0))
              for c in CK}
        out.append("Q_H2O at %s UTC: %s s^-1; of it from seen facets: %s (README Facts / A15: the total is set by "
                   "the extension to unseen terrain, not by the data)." % (tt, order_text(q), order_text(ss)))
        dn = {c: capx(c, relcap, r"facets at this time: ([\d,]+ day, [\d,]+ night, [\d,]+ in cast shadow)",
                      "day / night / cast-shadow facets at " + tt) for c in CK}
        if len(set(dn.values())) == 1:
            out.append("The illumination classes are the same in all six (%s; captions)." % dn["softplus"])
        else:
            out.append("Illumination classes (captions): %s." % "; ".join("%s %s" % (c, dn[c]) for c in CK))
    k, panel = g["panels"][0]
    spans = {c: (cls_range(c, g, panel)[1], N.add(c, "decades covered by the classes, %s %s" % (g["q"], panel),
                                                  cls_range(c, g, panel)[1], setp(c, man_name(g)),
                                                  "class_tables/%s: log10_to(last) - log10_from(first)" % panel, "%g"))
             for c in CK}
    other = max(v[0] for c, v in spans.items() if c != "C")
    note = ("; C's classes cover far more decades than the others' (its seen non-zero values range over many orders of "
            "magnitude, captions)" if spans["C"][0] > 2 * other else "")
    out.append("Decades spanned by the colour classes in %s: %s%s." % (PANEL_NAME[panel], order_text(spans), note))
    return " ".join(out)


def map_observe(g):
    q = "A" if g["q"] == "activity" else "F_H2O"
    if g["fam"] == "map_3d":
        return ("In the close-up of the late descent (c): softplus and B show an elliptical active patch just east of "
                "the impact point under the last part of the track (B on a non-zero background); D2 shows thin elongated "
                "active streaks instead; A and D smooth fields; C large polygonal cells (its 150 m hat-function basis) "
                "with high values around the impact point. Panel (a): softplus and D2 are light green (%s = 0) almost "
                "everywhere on the seen surface." % q)
    if g["fam"] in ("map_tpl", "map_latlon", "map_panel"):
        s = ("softplus and D2: almost the whole seen surface is light green (%s = 0) with a few isolated active patches "
             "(elongated or elliptical); A and D: smooth large-scale bands over the whole seen surface; B: a nearly uniform "
             "background level with a few brighter patches and bands; C: a patchwork of polygonal cells (the 150 m "
             "hat-function basis of option C)." % q)
        if g["fam"] == "map_tpl":
            s += " Panels (c)-(f) are the same in all six (fit-independent; the separate panel PNGs are byte-identical)."
        if g["q"] != "activity":
            s += " The patterns are those of the activity maps (F_H2O = q_H2O A), with relief shading."
        return s
    return ""


def geom_about(g):
    t = tload(caption_path("softplus", g["rel"]))
    return t.split(". ")[0] + ". " + ("Effective gravitational slope: angle between the inward facet normal and the "
                                      "effective gravity (homogeneous polyhedron, 533 kg/m^3, plus the centrifugal term); "
                                      "values above 90 deg are overhangs. Cliff templates: slope >= 45 deg (light) and "
                                      ">= 60 deg (dark).")


def geom_note(g):
    return ("Fit-independent geometry; the PNG is also byte-identical to the corresponding activity/panels PNG and to the "
            "panels of the other production-rate time (same SHA256)."
            if g["section"] == "production_rate" else
            "Fit-independent geometry (effective slope from the admitted gravity cache of the genus-4 mesh, an experiment). "
            "The same image appears in the production_rate/panels of both times.")


# source areas --------------------------------------------------------------------------------------------------------
def sa_table(c, gauge, P):
    p = setp(c, "checks/source_area_tables.txt")
    t = tload(p)
    pat = r"\[%s, P = %d %%\] facets with N > 0: (\d+) \(= union: (\w+)\), max N (\d+); area (\S+) m\^2" % (gauge, P)
    m = re.search(pat, t)
    return m, p, pat


def sa_val(c, gauge, P, what):
    m, p, pat = sa_table(c, gauge, P)
    grp = {"facets": 1, "maxN": 3, "area": 4}[what]
    v = float(m.group(grp)) if what == "area" else int(m.group(grp))
    fmt = e2 if what == "area" else (intc if what == "facets" else "%d")
    return v, N.add(c, "source area %s P%d %s" % (gauge, P, what), v, p, "regex group %d: %s" % (grp, pat), fmt)


def sa_box(c, rel):
    return capx(c, rel, r"box in the global map: lon (\S+ ?- ?\S+ deg, lat \S+-\S+ deg)", "close-up box (lon, lat)")


def sa_grain(c, rel):
    return capx(c, rel, r"its share of the model NG pressure is ([\d.]+-[\d.]+, median [\d.]+)",
                "NG grain-term share of the model NG pressure in the post-gap observations (range, median)")


def sa_numbers(g):
    hdr, rows = ["case"], []
    if g["fam"] == "sa_union":
        for ga in g["gauges"]:
            hdr += ["%s U50 facets / area [m^2]" % ga, "%s U90 facets / area [m^2]" % ga]
    else:
        for ga in g["gauges"]:
            hdr += ["%s facets with N_%d > 0" % (ga, g["P"]), "%s area [m^2]" % ga, "%s max N" % ga]
    hdr.append("close-up box (caption)")
    for c in CK:
        r = [CASE[c]["label"]]
        for ga in g["gauges"]:
            if g["fam"] == "sa_union":
                for P in (50, 90):
                    r.append("%s / %s" % (sa_val(c, ga, P, "facets")[1], sa_val(c, ga, P, "area")[1]))
            else:
                r += [sa_val(c, ga, g["P"], "facets")[1], sa_val(c, ga, g["P"], "area")[1], sa_val(c, ga, g["P"], "maxN")[1]]
        r.append(sa_box(c, g["rel"]))
        rows.append(r)
    w = [1.5] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)
    return hdr, rows, w


def sa_about(g):
    if g["fam"] == "sa_union":
        return ("U_P = the facets in the P % source area S_P of at least one of the 95 post-gap observations (10:33:10-"
                "10:39:26 UTC), P = 50 (darkest) to 90 % (lightest), each facet coloured by the smallest P whose union "
                "contains it; global map and a close-up around the impact (box fitted to the union facets). S_P(k, g): "
                "the smallest set of facets, in decreasing order of their contribution to the modelled surface pressure "
                "of observation k at gauge g, that supplies P % of it; NG grain term excluded.")
    return ("N_%d = the number of the 95 post-gap observations whose %d %% source area contains the facet (classes 1, 2-4, "
            "5-9, ..., 80-94, 95, the same in every case; grey N = 0); global map and close-up around the impact (box "
            "fitted per case)." % (g["P"], g["P"]))


def sa_compare(g):
    out = []
    P = 50 if g["fam"] == "sa_union" else g["P"]
    for ga in g["gauges"]:
        a = {c: sa_val(c, ga, P, "area") for c in CK}
        out.append("%s P%d area: %s m^2." % (ga, P, order_text(a)))
        if g["fam"] == "sa_union":
            a9 = {c: sa_val(c, ga, 90, "area") for c in CK}
            out.append("%s P90 area: %s m^2." % (ga, order_text(a9)))
    lo = min(sa_val(c, g["gauges"][0], P, "area")[0] for c in CK)
    hi = max(sa_val(c, g["gauges"][0], P, "area")[0] for c in CK)
    out.append("The %s P%d areas span a factor %s between the cases (FINDINGS A15: the size of the source areas is "
               "model-dependent)." % (g["gauges"][0], P, N.add("all", "ratio max/min of %s P%d source-area union area"
                                                                 % (g["gauges"][0], P), hi / lo,
                                                                 "checks/source_area_tables.txt of the six sets",
                                                                 "computed max/min of the six areas", "%.0f")))
    gr = {c: (0, sa_grain(c, g["rel"])) for c in CK}
    if "NG" in g["gauges"]:
        out.append("NG grain-term share of the model NG in these observations (range, median; excluded from the source "
                   "areas): " + "; ".join("%s %s" % (c, gr[c][1]) for c in CK) + ".")
    return " ".join(out)


def sa_observe(g):
    return ("softplus: a compact elliptical region just east of the impact point; B: a compact ellipse at about the "
            "same place (P50) with scattered facets of the larger P over a wider area; A and D: extended areas of many "
            "degrees in longitude and latitude around the impact point, A with a further area north-east of the impact "
            "region in the global map; C: only small, thin features along steep terrain, hardly visible at this scale; "
            "D2: narrow elongated features (a north-south streak south-east of the impact point and features farther "
            "east). The close-up boxes differ between cases (table), so equal sizes on the page are not equal areas.")


# contribution density ------------------------------------------------------------------------------------------------
def cd_v(c, key, ga, field, quantity, fmt, scale=1.0):
    return manx(c, "contribution_density/contribution_density_manifest.json", ["ranges", key, ga, field], quantity, fmt,
                scale)


def cd_numbers(g):
    hdr = ["case"]
    for ga in g["gauges"]:
        hdr += ["%s P50 facets / area [m^2]" % ga, "%s P90 facets / area [m^2]" % ga, "%s signal drawn in close-up" % ga]
    rows = []
    for c in CK:
        r = [CASE[c]["label"]]
        for ga in g["gauges"]:
            k = g["key"]
            r.append("%s / %s" % (cd_v(c, k, ga, "P50_facets", "%s %s P50 facets" % (k, ga), intc),
                                  cd_v(c, k, ga, "P50_area_m2", "%s %s P50 area [m^2]" % (k, ga), e2)))
            r.append("%s / %s" % (cd_v(c, k, ga, "P90_facets", "%s %s P90 facets" % (k, ga), intc),
                                  cd_v(c, k, ga, "P90_area_m2", "%s %s P90 area [m^2]" % (k, ga), e2)))
            r.append(cd_v(c, k, ga, "signal_share_drawn_zoom", "%s %s share of the signal drawn in the close-up" % (k, ga),
                          pct0, 100.0))
        rows.append(r)
    w = [1.5] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)
    return hdr, rows, w


def cd_about(g):
    t = tload(caption_path("softplus", g["rel"]))
    first = t.split(". Colour:")[0]
    return (first + ". Colour: the contribution per unit surface area as the fitted model computes it (NG: density at "
            "the nude gauge per m^2 of surface; RG: flux into the ram gauge per m^2), log scale per map; white contours: "
            "50 % (solid) and 90 % (dashed) source areas of the map; close-up box fitted per case.")


def cd_compare(g):
    out = []
    for ga in g["gauges"]:
        a = {c: (jload(setp(c, "contribution_density/contribution_density_manifest.json"))["ranges"][g["key"]][ga]
                 ["P50_area_m2"], cd_v(c, g["key"], ga, "P50_area_m2", "%s %s P50 area [m^2]" % (g["key"], ga), e2))
             for c in CK}
        d = {c: (jload(setp(c, "contribution_density/contribution_density_manifest.json"))["ranges"][g["key"]][ga]
                 ["signal_share_drawn_zoom"], cd_v(c, g["key"], ga, "signal_share_drawn_zoom",
                                                   "%s %s share of the signal drawn in the close-up" % (g["key"], ga),
                                                   pct0, 100.0)) for c in CK}
        fa = {}
        for c in CK:
            nf = cd_v(c, g["key"], ga, "P50_facets", "%s %s P50 facets" % (g["key"], ga), intc)
            fa[c] = "%s (%s facet%s)" % (a[c][1], nf, "" if nf == "1" else "s")
        a2 = {c: (a[c][0], fa[c]) for c in CK}
        out.append("%s: P50 area of this map %s m^2; share of the map's signal drawn in the close-up %s (the rest lies "
                   "outside the close-up, on overhangs or under other surface layers)." % (ga, order_text(a2), order_text(d)))
    return " ".join(out)


def cd_observe(g):
    k = g["key"]
    if k == "103926":
        return ("At the last observation every case shows only small patches next to the impact point in the close-up; "
                "where the drawn share (table) is small, most of the map's signal comes from facets not drawn in the "
                "lon/lat projection (overhangs, other surface layers).")
    if k == "pre_gap":
        return ("Before the NG gap A and D spread their contribution over large parts of the big lobe (global maps: a broad "
                "northern band); softplus, B and D2 show a few isolated patches; C shows almost nothing at this scale.")
    return ("The pattern follows the source areas: softplus an elliptical patch east of the impact point; B a compact core "
            "with a wide low halo; A and D broad smooth areas around the impact point; C a few small features on steep "
            "terrain; D2 narrow elongated streaks.")


# cliffs --------------------------------------------------------------------------------------------------------------
def sb(c, fname, i, field, quantity, fmt=pct, scale=100.0):
    return manx(c, "checks/" + fname, ["bins", i, field], quantity, fmt, scale)


def cliff_bins_numbers(g):
    hdr = ["case"] + ["%s deg: m / seen prod." % b for b in BIN_LABEL]
    rows = []
    for c in CK:
        r = [CASE[c]["label"]]
        for i in range(7):
            m = sb(c, "slope_bins_h2o_103310.json", i, "multiplier", "multiplier bin %s" % BIN_LABEL[i], "%.3g", 1.0)
            s = sb(c, "slope_bins_h2o_103310.json", i, "production_share_seen",
                   "seen H2O production share 10:33:10 bin %s" % BIN_LABEL[i])
            r.append("%s / %s" % (m, s))
        rows.append(r)
    return hdr, rows, [1.5] + [1.18] * 7


def cliff_bins_compare(g):
    st = {}
    for c in CK:
        p = setp(c, "checks/slope_bins_h2o_103310.json")
        v = sum(b["production_share_seen"] for b in jload(p)["bins"][3:])
        st[c] = (v, N.add(c, "seen H2O production share 10:33:10 from >= 60 deg (sum of bins 60-180)", 100.0 * v, p,
                          "sum bins[3:].production_share_seen", pct))
    al = {}
    for c in CK:
        p = setp(c, "checks/slope_bins_h2o_103310.json")
        v = sum(b["production_share_all"] for b in jload(p)["bins"][3:])
        al[c] = (v, N.add(c, "whole-surface H2O production share 10:33:10 from >= 60 deg", 100.0 * v, p,
                          "sum bins[3:].production_share_all", pct))
    sa = sb("softplus", "slope_bins_h2o_103310.json", 3, "bin_share_of_seen_area", "seen-area share bin 60-70")
    return ("Panel (a) (area shares) is the same in all six (geometry). Seen production from >= 60 deg at 10:33:10 "
            "(solid bars of panel c): %s; whole surface (hatched): %s. Panel (b): D2 has no multipliers (all 1; "
            "degenerate panel), A has all six at the bound 1, C has 70-80 and 80-90 deg at ~1e5 (not determined, A12). "
            % (order_text(st, reverse=True), order_text(al, reverse=True)) + cliff_bins_dominant(g))


def cliff_bins_dominant(g):
    out = []
    for field, nm in (("production_share_seen", "seen surface (solid)"), ("production_share_all", "whole surface (hatched)")):
        parts = []
        for c in CK:
            p = setp(c, "checks/slope_bins_h2o_103310.json")
            bins = jload(p)["bins"]
            i = max(range(7), key=lambda j: bins[j][field])
            parts.append("%s %s deg (%s)" % (c, BIN_LABEL[i], N.add(c, "largest bin of the %s at 10:33:10" % nm,
                                                                    100.0 * bins[i][field], p,
                                                                    "bins[%d].%s (max over bins)" % (i, field), pct0)))
        out.append("Largest bin of the H2O production of the %s at 10:33:10 (slope_bins_h2o_103310.json): %s." %
                   (nm, "; ".join(parts)))
    return " ".join(out)


def cliff_pf_numbers(g):
    hdr = ["case"]
    for ga in g["gauges"]:
        hdr += ["%s >= 60 deg before late window (caption)" % ga, "%s >= 60 deg in late window (caption)" % ga,
                "%s mean >= 60 deg late (facts)" % ga, "%s 90-180 deg at last obs. (facts)" % ga]
    rows = []
    relc = "cliffs/12M_cliff_pressure_fraction_NG_RG.png"
    for c in CK:
        r = [CASE[c]["label"]]
        for ga in g["gauges"]:
            pat = r"%s: slopes >= 60 deg give (\d+-\d+) %% of the surface pressure before the late window and (\d+-\d+) %% in it" % ga
            r.append(capx(c, relc, pat, "%s >= 60 deg share before the late window (range)" % ga, 1) + " %")
            r.append(capx(c, relc, pat, "%s >= 60 deg share in the late window (range)" % ga, 2) + " %")
            lm = fxv(c, ["pressure_share_by_bin", ga, "late_mean"])
            r.append(N.add(c, "%s mean share >= 60 deg in the late window" % ga, 100.0 * sum(lm[3:]),
                           setp(c, "checks/facts.json"), "facts.json: sum pressure_share_by_bin/%s/late_mean[3:]" % ga,
                           pct0))
            r.append(fx(c, ["pressure_share_by_bin", ga, "last", 6], "90-180 deg share of %s at the last observation" % ga,
                        pct0, 100.0))
        rows.append(r)
    return hdr, rows, [1.5] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def cliff_pf_compare(g):
    out = []
    for ga in g["gauges"]:
        lm = {}
        for c in CK:
            v = sum(fxv(c, ["pressure_share_by_bin", ga, "late_mean"])[3:])
            lm[c] = (v, N.add(c, "%s mean share >= 60 deg in the late window" % ga, 100.0 * v, setp(c, "checks/facts.json"),
                              "facts.json: sum pressure_share_by_bin/%s/late_mean[3:]" % ga, pct0))
        out.append("%s, mean share of slopes >= 60 deg in the late window: %s." % (ga, order_text(lm, reverse=True)))
    lastng = {c: (fxv(c, ["pressure_share_by_bin", "NG", "last"])[6],
                  fx(c, ["pressure_share_by_bin", "NG", "last", 6], "90-180 deg share of NG at the last observation",
                     pct0, 100.0)) for c in CK}
    if "NG" in g["gauges"]:
        out.append("At the last observation (5.5 m) the overhang bin (90-180 deg) gives %s of the modelled NG (A15: NG "
                   "sees through the shield onto facets behind the gauge)." % order_text(lastng, reverse=True))
    if "RG" in g["gauges"]:
        lastrg = {c: (fxv(c, ["pressure_share_by_bin", "RG", "last"])[6],
                      fx(c, ["pressure_share_by_bin", "RG", "last", 6], "90-180 deg share of RG at the last observation",
                         pct0, 100.0)) for c in CK}
        out.append("RG at the last observation from 90-180 deg: %s." % order_text(lastrg))
    return " ".join(out)


def cliff_pf_observe(g):
    return ("A is dominated by 0-40 deg (grey) almost everywhere; C by 70-90 deg (light blue and green) over the whole "
            "descent; D by 40-50 deg (yellow) and 0-40 deg; softplus and B by 60-70 deg (dark blue) and 0-40 deg; D2 by "
            "0-40 deg before the late window and by 60-80 deg in it. In the NG panel the overhang bin (purple) grows in "
            "the last metres in all six.")


# regions ------------------------------------------------------------------------------------------------------------
def reg_stat(c, kind, ga, region, phase, idx, quantity, fmt=pct0):
    return manx(c, "regions/regions_figures_manifest.json", ["stats", kind, ga, region, phase, idx], quantity, fmt, 100.0)


LOBES = ["small lobe (head, with the touchdown point)", "big lobe (body)"]
SHELLS = ["distance 0–30 m", "distance 30–100 m", "distance 100–300 m", "distance 1000–3000 m"]


def reg_numbers(g):
    hdr, rows = ["case"], []
    if g["fam"] == "reg_lobe":
        for ga in g["gauges"]:
            hdr += ["%s small lobe median before / in late window" % ga, "%s small lobe range in late window" % ga]
    else:
        for ga in g["gauges"]:
            hdr += ["%s median in late window: 0-30 / 30-100 / 100-300 m" % ga, "%s 1000-3000 m median before" % ga]
    for c in CK:
        r = [CASE[c]["label"]]
        for ga in g["gauges"]:
            if g["fam"] == "reg_lobe":
                L = LOBES[0]
                r.append("%s / %s" % (reg_stat(c, "lobe", ga, L, "before_late_window", 2, "%s small lobe median before" % ga),
                                      reg_stat(c, "lobe", ga, L, "late_window", 2, "%s small lobe median late" % ga)))
                r.append("%s-%s" % (reg_stat(c, "lobe", ga, L, "late_window", 0, "%s small lobe min late" % ga, "%.0f"),
                                    reg_stat(c, "lobe", ga, L, "late_window", 1, "%s small lobe max late" % ga)))
            else:
                r.append(" / ".join(reg_stat(c, "distance", ga, s, "late_window", 2, "%s %s median late" % (ga, s))
                                    for s in SHELLS[:3]))
                r.append(reg_stat(c, "distance", ga, SHELLS[3], "before_late_window", 2,
                                  "%s %s median before" % (ga, SHELLS[3])))
        rows.append(r)
    return hdr, rows, [1.5] + [(8.2 / (len(hdr) - 1))] * (len(hdr) - 1)


def reg_about(g):
    t = tload(caption_path("softplus", g["rel"]))
    first = t.split(", along the descent")[0]
    return (first + ", along the descent (stacked, against the gauge altitude, late window boxed; NG gap hatched in the NG "
            "panel). The region sums equal the model surface pressure; the NG grain term is not included.")


def reg_compare(g):
    out = []
    for ga in g["gauges"]:
        if g["fam"] == "reg_lobe":
            v = {c: (jload(setp(c, "regions/regions_figures_manifest.json"))["stats"]["lobe"][ga][LOBES[0]]["late_window"][2],
                     reg_stat(c, "lobe", ga, LOBES[0], "late_window", 2, "%s small lobe median late" % ga)) for c in CK}
            b = {c: (jload(setp(c, "regions/regions_figures_manifest.json"))["stats"]["lobe"][ga][LOBES[1]]
                     ["before_late_window"][2],
                     reg_stat(c, "lobe", ga, LOBES[1], "before_late_window", 2, "%s big lobe median before" % ga))
                 for c in CK}
            out.append("%s: small-lobe median in the late window %s; big-lobe median before it %s." %
                       (ga, order_text(v), order_text(b, reverse=True)))
        else:
            v = {c: (jload(setp(c, "regions/regions_figures_manifest.json"))["stats"]["distance"][ga][SHELLS[1]]
                     ["late_window"][2], reg_stat(c, "distance", ga, SHELLS[1], "late_window", 2,
                                                  "%s %s median late" % (ga, SHELLS[1]))) for c in CK}
            out.append("%s: median share from 30-100 m around the touchdown point in the late window %s." %
                       (ga, order_text(v, reverse=True)))
    if g["fam"] == "reg_lobe" and "NG" in g["gauges"]:
        out.append("The README Facts give the big-lobe NG median before the late window with a different row selection "
                   "(full_set_facts.py: rows before the NG gap and observation rows), so their values can differ "
                   "slightly from the medians here (figure manifest).")
    return " ".join(out)


def reg_observe(g):
    if g["fam"] == "reg_lobe":
        return ("In all six the small-lobe share rises from a small fraction before the late window to nearly 100 % near "
                "the end. At the start of the late window it is largest in C (NG and RG) and smallest in softplus and D2 "
                "(NG).")
    return ("In all six the 1000-3000 m shell dominates above the late window and the 30-100 m and 0-30 m shells take over "
            "near the end; in C the 300-1000 m shell (green) carries a large share in the upper late window, and in D "
            "the NG panel shows a temporary return of the 300-1000 m and 100-300 m shells at about 20 m.")


# visibility -----------------------------------------------------------------------------------------------------------
def vis_about(g):
    t = tload(caption_path("softplus", g["rel"]))
    first = t.split(". ")[0]
    return first + "."


def vis_note(g):
    rel = g["rel"]
    t = tload(caption_path("softplus", rel))
    if "NG-only facets (orange," in t:
        a = capx("softplus", rel, r"NG-only facets \(orange, ([\d,]+)\)", "NG-only facets (late visibility)")
        b = capx("softplus", rel, r"RG-only facets \(blue, ([\d,]+)\)", "RG-only facets (late visibility)")
        s = "Late-visibility templates: %s NG-only and %s RG-only facets. " % (a, b)
    else:
        s = ""
    imp = capx("softplus", rel, r"to the shape model: (10:39:\S+ UTC, longitude [\d.]+°, latitude [\d.]+°)",
               "impact point (time, lon, lat)")
    return s + "Impact point (straight-line extrapolation of the last two observation rows): %s." % imp


FAMILY = {
    "md_pub": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_alt": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_ng": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_rg": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_ngrg": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_time": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "md_press": dict(about=md_about, numbers=md_numbers, compare=md_compare, observe=md_observe),
    "map_tpl": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_latlon": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_3d": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "map_panel": dict(about=map_about, numbers=map_numbers, compare=map_compare, observe=map_observe),
    "geom_panel": dict(about=geom_about, identical_note=geom_note),
    "sa_union": dict(about=sa_about, numbers=sa_numbers, compare=sa_compare, observe=sa_observe),
    "sa_occ": dict(about=sa_about, numbers=sa_numbers, compare=sa_compare, observe=sa_observe),
    "cd": dict(about=cd_about, numbers=cd_numbers, compare=cd_compare, observe=cd_observe),
    "cliff_bins": dict(about=lambda g: tload(caption_path("softplus", g["rel"])).split(". (a)")[0].replace(
        " (the softplus activity model with slope multipliers >= 1, regularized stage-2 chain, final state after 5 fit "
        "and fine-tuning cycles)", "") + ": (a) share of the total and of the seen area per bin and the seen fraction "
        "of each bin, (b) the fitted multiplier m_b (log scale; 0-40 deg held at 1), (c) share of the H2O production of "
        "the seen surface (solid) and of the whole surface (hatched) per bin at 10:33:10 and 10:39:26 UTC.",
                       numbers=cliff_bins_numbers, compare=cliff_bins_compare,
                       observe=("Panel (b) of D2 shows 1 in every bin (no multipliers); panel (a) looks the same in all "
                                "six; the hatched (whole-surface) bars of C are concentrated in the steep bins like its "
                                "solid ones, while in softplus, B, A and D2 the hatched 0-40 deg bar is the tallest.")),
    "cliff_pf": dict(about=lambda g: ("Fraction of the modelled surface pressure contributed by each effective-slope bin "
                                      "along the descent (stacked; one colour per bin: 0-40 grey, 40-50 yellow, 50-60 "
                                      "orange, 60-70 dark blue, 70-80 light blue, 80-90 green, 90-180 deg purple), "
                                      "against the gauge altitude; late window boxed; NG gap hatched in the NG panel; "
                                      "NG grain term not included."),
                     numbers=cliff_pf_numbers, compare=cliff_pf_compare, observe=cliff_pf_observe),
    "reg_lobe": dict(about=reg_about, numbers=reg_numbers, compare=reg_compare, observe=reg_observe),
    "reg_distance": dict(about=reg_about, numbers=reg_numbers, compare=reg_compare, observe=reg_observe),
    "vis": dict(about=vis_about, identical_note=vis_note),
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
            for line in open("/data001/highresolu_vt/AMPS/.git/packed-refs"):
                if line.strip().endswith(ref):
                    return line.split()[0] + " (%s, packed-refs)" % ref
        return head
    except OSError as e:
        return "unknown (%s)" % e


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
    assert len(groups) == 85, len(groups)
    # hashes, identity, aspects
    jobs, img = [], {}
    for g in groups:
        g["sha"] = {c: sha256(src_path(c, g["rel"])) for c in CK}
        g["identical"] = len(set(g["sha"].values())) == 1
        sizes = {}
        for c in (["softplus"] if g["identical"] else CK):
            with Image.open(src_path(c, g["rel"])) as im:
                sizes[c] = im.size
        g["native"] = sizes
        g["layout"] = layout([h / w for (w, h) in sizes.values()], g["identical"])
        kind = "png" if g["section"] in LINE_SECTIONS else "jpg"
        for c in sizes:
            nat = sizes[c][0]
            target = int(min(nat, max(1000, round(g["layout"]["width_in"] * 300)), 1400))
            name = re.sub(r"[^A-Za-z0-9_.-]", "_", g["rel"].replace("<STEM>", "STEM").replace("/", "__"))
            dst = os.path.join(out, "images", g["section"], "%s__%s.%s" % (name[:-4], c, kind))
            jobs.append((src_path(c, g["rel"]), dst, kind, target))
            img[(g["rel"], c)] = dict(src=src_path(c, g["rel"]), dst=dst, kind=kind, target=target)
    with ProcessPoolExecutor(args.jobs) as ex:
        res = list(ex.map(prep_image, jobs))
    tot = 0
    for (src, dst, kind, target), (d, size, nbytes) in zip(jobs, res):
        tot += nbytes
        for k, v in img.items():
            if v["dst"] == d:
                v["px"] = size
                v["bytes"] = nbytes
    args.img_note = "%d images, %.1f MB, JPEG for maps and 256-colour PNG for line plots, 1000-1400 px wide" % (
        len(jobs), tot / 1e6)
    build_doc(out, groups, img, args)
    # products next to the document
    with open(os.path.join(out, "numbers_used.json"), "w") as f:
        json.dump(dict(created_utc=args.started, generator=os.path.abspath(__file__), n=len(N.items), numbers=N.items),
                  f, indent=1, ensure_ascii=False)
    man = dict(created_utc=args.started,
               finished_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S"),
               command=sys.argv, python=sys.executable, python_docx=args.docx_version, pillow=args.pil_version,
               slurm_job=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename, amps_head=amps_head(),
               generator_sha256=sha256(os.path.abspath(__file__)), document=os.path.join(out, DOC_NAME),
               document_bytes=os.path.getsize(os.path.join(out, DOC_NAME)), images_bytes=tot,
               groups=[dict(fig=g["fig"], section=g["section"], rel=g["rel"], title=g["title"], fam=g["fam"],
                            identical=g["identical"], layout=g["layout"], native=g["native"], sha256=g["sha"],
                            images={c: dict(dst=os.path.relpath(img[(g["rel"], c)]["dst"], out),
                                            px=img[(g["rel"], c)].get("px"), bytes=img[(g["rel"], c)].get("bytes"))
                                    for c in (["softplus"] if g["identical"] else CK)}) for g in groups])
    with open(os.path.join(out, "build_manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    if N.unverified:
        with open(os.path.join(out, "numbers_unverified.json"), "w") as f:
            json.dump(N.unverified, f, indent=1, ensure_ascii=False)
        print("WARNING: %d numeric tokens of hand-written text not found in a source file (numbers_unverified.json)"
              % len(N.unverified))
    print("document %s: %.1f MB; %d plot groups (%d identical in all six); %d images (%.1f MB); %d numbers" % (
        os.path.join(out, DOC_NAME), man["document_bytes"] / 1e6, len(groups), sum(g["identical"] for g in groups),
        len(jobs), tot / 1e6, len(N.items)))


def cmd_validate(args):
    """reopen the document, count figures / captions / headings / images, check sizes against the pages, previews"""
    from docx import Document
    from docx.shared import Emu
    out = os.path.abspath(args.out)
    path = os.path.join(out, DOC_NAME)
    man = json.load(open(os.path.join(out, "build_manifest.json")))
    doc = Document(path)
    heads = {}
    caps = []
    for p in doc.paragraphs:
        sn = p.style.name if p.style is not None else ""
        if sn.startswith("Heading") or sn == "Title":
            heads[sn] = heads.get(sn, 0) + 1
        if sn == "Caption" and re.match(r"^Figure \d+: ", p.text):
            caps.append(p.text)
    shapes = doc.inline_shapes
    exp_imgs = sum(1 if g["identical"] else 6 for g in man["groups"])
    # per section: page geometry and the pictures in its body
    body = doc.element.body
    problems = []
    sect_info = []
    for s in doc.sections:
        sect_info.append(dict(orient=str(s.orientation), w_in=s.page_width / 914400, h_in=s.page_height / 914400,
                              usable_w=(s.page_width - s.left_margin - s.right_margin) / 914400,
                              usable_h=(s.page_height - s.top_margin - s.bottom_margin) / 914400))
    # pictures in order; each figure group = one section after the section intro pages
    pics = []
    for sh in shapes:
        pics.append((sh.width / 914400, sh.height / 914400))
    # walk the groups in order and compare picture sizes with the page of their layout
    k = 0
    for g in man["groups"]:
        n = 1 if g["identical"] else 6
        o = g["layout"]["orient"]
        W, H = usable(o)
        cols = 1 if g["identical"] else g["layout"]["cols"]
        rows = math.ceil(n / cols)
        sizes = pics[k:k + n]
        k += n
        wmax = max(w for w, h in sizes)
        hrow = max(h for w, h in sizes)
        grid_h = rows * (hrow + LABEL + 0.05)
        if wmax > W / cols - 0.05:
            problems.append("Figure %d: image width %.2f in > cell %.2f in" % (g["fig"], wmax, W / cols))
        if grid_h + OVERHEAD - 0.4 > H:
            problems.append("Figure %d: grid height %.2f in + overhead > page %.2f in" % (g["fig"], grid_h, H))
        g["_check"] = dict(width_in=round(wmax, 3), max_h_in=round(hrow, 3), grid_h_in=round(grid_h, 3),
                           page_usable=[W, H])
    from docx.oxml.ns import qn as _qn
    W_NS = _qn("w:tbl")
    sect_tables, cur = [], []
    for el in body.iterchildren():
        if el.tag == W_NS:
            cur.append(el)
        sp = el.find(".//" + _qn("w:sectPr")) if el.tag == _qn("w:p") else (el if el.tag == _qn("w:sectPr") else None)
        if sp is not None:
            pg = sp.find(_qn("w:pgSz"))
            mar = sp.find(_qn("w:pgMar"))
            uw = (int(pg.get(_qn("w:w"))) - int(mar.get(_qn("w:left"))) - int(mar.get(_qn("w:right")))) / 1440.0
            for tb in cur:
                grid = tb.find(_qn("w:tblGrid"))
                tw = sum(int(gc.get(_qn("w:w"))) for gc in grid.findall(_qn("w:gridCol"))) / 1440.0
                sect_tables.append(round(tw, 3))
                if tw > uw + 0.02:
                    problems.append("table width %.2f in > usable %.2f in" % (tw, uw))
            cur = []
    import zipfile
    with zipfile.ZipFile(path) as zf:
        folder_entries = [n for n in zf.namelist() if n.endswith("/")]
    if folder_entries:
        problems.append("zip folder entries: %s" % folder_entries[:5])
    nums = json.load(open(os.path.join(out, "numbers_used.json")))
    res = dict(document=path, bytes=os.path.getsize(path), mb=round(os.path.getsize(path) / 1e6, 2),
               headings=heads, figure_captions=len(caps), figure_groups=len(man["groups"]),
               identical_groups=sum(g["identical"] for g in man["groups"]), inline_images=len(shapes),
               expected_images=exp_imgs, sections=len(doc.sections), tables=len(doc.tables),
               numbers_recorded=nums["n"], problems=problems, tables_width_checked=len(sect_tables),
               captions_numbered_1_to_n=[int(re.match(r"^Figure (\d+):", c).group(1)) for c in caps] == list(
                   range(1, len(caps) + 1)),
               libreoffice="not available (no soffice on login or compute nodes; not on conda-forge): no PDF render",
               checks=[dict(fig=g["fig"], **g["_check"]) for g in man["groups"]])
    res["passed"] = (res["figure_captions"] == 85 and res["inline_images"] == exp_imgs and not problems and
                     res["captions_numbered_1_to_n"])
    os.makedirs(os.path.join(out, "validation"), exist_ok=True)
    with open(os.path.join(out, "validation", "validation.json"), "w") as f:
        json.dump(res, f, indent=1)
    # layout previews (drawn from the same geometry; NOT a Word / LibreOffice rendering)
    from PIL import ImageDraw, ImageFont
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf", 11)
        font2 = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf", 11)
    except OSError:
        font = font2 = ImageFont.load_default()
    dpi = 100
    sel = [g for g in man["groups"] if g["fig"] in args.preview]
    for g in sel:
        o = g["layout"]["orient"]
        pw, ph = PAGE[o]
        im = Image.new("RGB", (int(pw * dpi), int(ph * dpi)), "white")
        dr = ImageDraw.Draw(im)
        dr.rectangle([MARGIN * dpi, MARGIN * dpi, (pw - MARGIN) * dpi, (ph - MARGIN) * dpi], outline=(200, 200, 200))
        dr.text((MARGIN * dpi, MARGIN * dpi), "3.x.y %s" % g["title"], fill=(31, 58, 95), font=font)
        dr.text((MARGIN * dpi, MARGIN * dpi + 22), "Files: <set>/%s" % g["rel"], fill=(0, 0, 0), font=font2)
        W, H = usable(o)
        cols = 1 if g["identical"] else g["layout"]["cols"]
        y0 = MARGIN + 0.5
        cases = list(g["images"].keys())
        w_in = g["layout"]["width_in"]
        rowh = []
        for k2, c in enumerate(cases):
            r, col = divmod(k2, cols)
            src = os.path.join(out, g["images"][c]["dst"])
            pic = Image.open(src).convert("RGB")
            h_in = w_in * pic.size[1] / pic.size[0]
            rowh.append(h_in)
        hrow = max(rowh)
        for k2, c in enumerate(cases):
            r, col = divmod(k2, cols)
            x = MARGIN + col * (W / cols) + (W / cols - w_in) / 2
            y = y0 + r * (hrow + LABEL + 0.05)
            lab = "all six cases (identical files)" if g["identical"] else CASE[c]["label"]
            dr.text((x * dpi + 2, y * dpi), lab, fill=(160, 16, 16), font=font)
            pic = Image.open(os.path.join(out, g["images"][c]["dst"])).convert("RGB")
            pic = pic.resize((int(w_in * dpi), int(w_in * dpi * pic.size[1] / pic.size[0])), Image.LANCZOS)
            im.paste(pic, (int(x * dpi), int((y + LABEL) * dpi)))
        yc = y0 + math.ceil(len(cases) / cols) * (hrow + LABEL + 0.05) + 0.05
        dr.text((MARGIN * dpi, yc * dpi), "Figure %d: %s" % (g["fig"], g["title"]), fill=(0, 0, 0), font=font)
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
    v.add_argument("--preview", type=int, nargs="*", default=[1, 2, 9, 14, 36, 38, 54, 73, 76, 82, 83])
    args = ap.parse_args()
    if args.cmd == "build":
        cmd_build(args)
    else:
        cmd_validate(args)


if __name__ == "__main__":
    main()
