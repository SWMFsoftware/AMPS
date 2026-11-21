#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cliff_effect_docx.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cliff_effect_docx.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: ef24cc7c67ff5ff3a9c5a90aad8bf5d22b60a8759692decdd7fb78dd0a411d62
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Word document of the slope-at-equal-illumination analysis: the same plot of the six stage-2 cases in one
#   cluster, summaries, numbers_used.json (2026-11-17)
# RUN: separate venv cops_output/tools/docx_venv (python-docx, Pillow); build + validate in a Slurm job (HOW-TO-
#   RERUN.md of figures_cluster/12M/cliff_effect_stage2/); imports fullsets_comparison_docx.py from its directory
# ==== end of header ====
"""Word document of the slope-at-equal-illumination analysis of the six final states (2026-11-17, user: "test several
definitions: >40, >60, >80, as well as 0-40,40-50,50-60,60-70,70-80, 80-90; for the six cases individually. the task is to
find whether at the same illuminations condition the source rate depends on the slope angle; prepare the plots and a word
doc that would similarly to that one you have generated pack the same plot for all six cases in a cluster and add a
summary"). Layout and helpers from fullsets_comparison_docx.py (imported: the numbers registry N and T(), the page layout,
the image preparation). Every number printed is registered in numbers_used.json with its source file and locator; numeric
tokens of hand-written text must occur in a listed source file, otherwise the build stops (numbers_unverified.json).

    cliff_effect_docx.py build    --root <cliff_effect_stage2> [--jobs 8]
    cliff_effect_docx.py validate --root <cliff_effect_stage2>
Interpreter: $COPS_OUTPUT_ROOT/tools/docx_venv (python-docx, Pillow). Refuses to overwrite the document.
"""
import argparse
import csv
import datetime
import json
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fullsets_comparison_docx as FX  # noqa: E402

N = FX.N
DOC_NAME = "12M_stage2_cliff_effect.docx"
DOC_DATE = "2026-11-17"
CASES_JSON = os.path.join(HERE, "data", "cliff_effect_cases_20261117.json")
PLOTS = [("P1_NG", "Source rate and activity at equal illumination, NG-seen surface"),
         ("P1_RG", "Source rate and activity at equal illumination, RG-seen surface"),
         ("P1_whole", "Source rate and activity at equal illumination, whole surface"),
         ("P2_NG", "NG observation weight per m2 and its factors"),
         ("P2_RG", "RG observation weight per m2 and its factors"),
         ("P3_NG", "Matched comparison of the NG contribution per m2"),
         ("P3_RG", "Matched comparison of the RG contribution per m2"),
         ("P4_NG", "Exact decomposition of the NG signal share"),
         ("P4_RG", "Exact decomposition of the RG signal share")]
DEFS = [">40", ">60", ">80", "40-50", "50-60", "60-70", "70-80", "80-90", "90-180 (addition)"]
BINS = ["40-50", "50-60", "60-70", "70-80", "80-90", "90-180 (addition)"]
SNAP0, SNAP3 = "10:33:10 UTC (297 m)", "10:39:26 UTC (5.5 m)"
LATE = "late window (W10x5)"
FLAG = 0.10        # dagger: more than 10 % of the averaging weight in classes with an under-sampled 0-40 reference
SCRIPTS = [os.path.join(HERE, f) for f in ("cliff_effect_figures.py", "cliff_effect_data.py", "cliff_effect_docx.py")]
BUILD_FB = ("parameter of the analysis (definitions in cliff_effect_figures.py / cliff_effect_data.py)", SCRIPTS[0])


def g3(v):
    return "%.3g" % v if isinstance(v, (int, float)) and math.isfinite(v) else "-"


def jl(p):
    return FX.jload(p)


def csvrows(p):
    with open(p) as f:
        return list(csv.DictReader(f))


class Ctx:
    root = None
    cases = None


def ans_path(key):
    return os.path.join(Ctx.root, key, "answer.json")


def res(key, surface, tl, q, d):
    return jl(ans_path(key))["answer"][surface][tl][q].get(d)


def rcell(key, surface, tl, q, d, rng=True, flag=True):
    a = res(key, surface, tl, q, d)
    if a is None:
        return "-"
    src, loc = ans_path(key), "answer.%s.%s.%s.%s" % (surface, tl, q, d)
    t = N.add(key, "SZA-controlled %s ratio %s vs 0-40 (%s, %s)" % (q, d, surface, tl), a["ratio"], src, loc + ".ratio", g3)
    if rng:
        t += " (%s-%s)" % (N.add(key, "min over classes", a["min"], src, loc + ".min", g3),
                           N.add(key, "max over classes", a["max"], src, loc + ".max", g3))
    if flag and a.get("weight_fraction_undersampled_ref", 0.0) > FLAG:
        t += " †"
    return t


def vcell(key, surface, tl, q, d):
    a = res(key, surface, tl, q + "_refmin", d)
    if a is None:
        return "-"
    src, loc = ans_path(key), "answer.%s.%s.%s_refmin.%s" % (surface, tl, q, d)
    return "%s [%s]" % (N.add(key, "variant ratio (reference >= 1 %% of its area) %s %s %s" % (d, surface, tl), a["ratio"],
                               src, loc + ".ratio", g3),
                        N.add(key, "area fraction compared", a["area_fraction_compared"], src,
                              loc + ".area_fraction_compared", "%.2f"))


def build(args):
    from docx import Document
    from docx.enum.section import WD_ORIENT, WD_SECTION
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt, RGBColor
    from PIL import Image
    root = Ctx.root = os.path.abspath(args.root)
    path = os.path.join(root, DOC_NAME)
    if os.path.exists(path):
        raise SystemExit("%s exists" % path)
    spec = jl(CASES_JSON)
    CASES = Ctx.cases = [(c["key"], c["label"]) for c in spec["cases"]]
    SRC = SCRIPTS + [CASES_JSON, spec["source"]]
    # images: one cluster of six per plot type
    jobs, groups = [], []
    for pk, title in PLOTS:
        srcs = [os.path.join(root, k, pk + ".png") for k, _ in CASES]
        asp = []
        for s in srcs:
            with Image.open(s) as im:
                asp.append(im.size[1] / im.size[0])
        lay = FX.layout(asp, False)
        dsts = [os.path.join(root, "images", pk, "%s.png" % k) for k, _ in CASES]
        target = min(2400, int(round(lay["width_in"] * 450)))
        for s, d in zip(srcs, dsts):
            jobs.append((s, d, "png", target))
        groups.append(dict(key=pk, title=title, srcs=srcs, dsts=dsts, layout=lay, sha=[FX.sha256(s) for s in srcs]))
    with ProcessPoolExecutor(args.jobs) as ex:
        info = list(ex.map(FX.prep_image, jobs))
    doc = Document()
    st = doc.styles
    st["Normal"].font.name = "Arial"
    st["Normal"].font.size = Pt(10)
    st["Normal"].element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    st["Normal"].paragraph_format.space_after = Pt(4)
    for name, size in (("Title", 20), ("Heading 1", 15), ("Heading 2", 12), ("Heading 3", 10.5)):
        s = st[name]
        s.font.name = "Arial"
        s.font.size = Pt(size)
        s.font.bold = True
        s.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)
        s.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
        s.paragraph_format.keep_with_next = True
    cap = st["Caption"]
    cap.font.name = "Arial"
    cap.font.size = Pt(9)
    cap.font.italic = False
    cap.font.bold = False
    cap.font.color.rgb = RGBColor(0, 0, 0)
    zoom = doc.settings.element.find(qn("w:zoom"))          # python-docx template: w:zoom without the required percent
    if zoom is not None and zoom.get(qn("w:percent")) is None:
        zoom.set(qn("w:percent"), "100")
    cur = {"o": "P"}

    def orient(sec, o):
        w, h = FX.PAGE[o]
        sec.orientation = WD_ORIENT.PORTRAIT if o == "P" else WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = Inches(w), Inches(h)
        for a_ in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
            setattr(sec, a_, Inches(FX.MARGIN))
        cur["o"] = o

    def section(o):
        orient(doc.add_section(WD_SECTION.NEW_PAGE), o)

    orient(doc.sections[0], "P")

    def para(text, bold=None, size=None, style=None, src=None, fb=None):
        text = FX.T(text, (src or []) + SRC, fb)
        p = doc.add_paragraph(style=style)
        if bold:
            r = p.add_run(bold + " ")
            r.bold = True
            if size:
                r.font.size = Pt(size)
        r = p.add_run(text)
        if size:
            r.font.size = Pt(size)
        return p

    def bullet(text, src=None, fb=None):
        text = FX.T(text, (src or []) + SRC, fb)
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(text)
        p.paragraph_format.space_after = Pt(2)
        return p

    def table(header, rows, widths, size=7.5):
        avail = FX.usable(cur["o"])[0] - 0.1
        if sum(widths) > avail:
            widths = [w * avail / sum(widths) for w in widths]
        t = doc.add_table(rows=1 + len(rows), cols=len(header))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, h in enumerate(header):
            c = t.rows[0].cells[i]
            c.text = ""
            r = c.paragraphs[0].add_run(h)
            r.bold = True
            r.font.size = Pt(size)
        for j, row in enumerate(rows):
            for i, v in enumerate(row):
                c = t.rows[j + 1].cells[i]
                c.text = ""
                r = c.paragraphs[0].add_run(str(v))
                r.font.size = Pt(size)
        for i, col in enumerate(t.columns):
            for c in col.cells:
                c.width = Inches(widths[i])
        doc.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    # ------------------------------------------------------------------------------------------- front matter
    N.ctx = "front matter"
    doc.add_paragraph("Does the source rate depend on the slope angle at the same illumination?", style="Title")
    para("Six final states of the 12M stage-2 fits, each analysed individually: " + "; ".join(l for _, l in CASES) +
         ". " + DOC_DATE + ".", fb=BUILD_FB)
    doc.add_heading("1. Question", level=1)
    para("User, verbatim: \"test several definitions: >40, >60, >80, as well as 0-40,40-50,50-60,60-70,70-80, 80-90; for "
         "the six cases individually. the task is to find whether at the same illuminations condition the source rate "
         "depends on the slope angle; prepare the plots and a word doc that would similarly to that one you have "
         "generated pack the same plot for all six cases in a cluster and add a summary\".")
    doc.add_heading("2. Model facts", level=1)
    for t in ("LiouvilleCOPS src/Source.cpp and Liouville.h with the base settings SurfaceBC flat, SZAMode continuous, "
              "CastShadows on.",
              "H2O source rate per m2 of facet i: F_i = A_i q(mu_i), q = Qmin + (Qmax - Qmin) cos SZA on the lit side, Qmin "
              "at night (cos SZA <= 0) and in a ray-traced cast shadow.",
              "NG contribution c_i = F_i beta(T(mu_i)) sqrt(pi)/2 dS_i g_i, g_i = Kn_i/dS_i = (cos theta_i/pi) W_NG f_vis / "
              "r_i^2; RG: c_i = F_i (KPhi_i + (V.b) beta sqrt(pi)/2 KV_i).",
              "At equal SZA and shadow state q and T are equal, so at fixed illumination the slope dependence of F is the "
              "slope dependence of the fitted activity A. A contains, per slope bin, the fitted slope multiplier (D2 has "
              "none) and the fitted activity map."):
        bullet(t)
    doc.add_heading("3. Definitions and defaults", level=1)
    for t in ("Slope: the fit's effective-gravitational-slope bins (nonneg_facets.npz). Categories: 0-40 deg (reference), "
              "40-50, 50-60, 60-70, 70-80, 80-90; thresholds >40, >60, >80 (everything above, including 90-180); and, as an "
              "addition the user did not list, the fit's 90-180 bin (overhangs) separately.",
              "Illumination: SZA in 10-deg classes on the lit side plus the separate classes shadowed (cos SZA > 0, cast "
              "shadow) and night (cos SZA <= 0).",
              "Times: the contribution-density snapshots 10:33:10 (297 m), 10:37:31.282 (100 m), 10:39:09.282 (20 m), "
              "10:39:26 (5.5 m); for the contribution plots (P3, P4) also the phase means pooled over their rows: before the "
              "NG gap, the gap (RG only), the late window as the W10x5 average of the 1-s grid rows.",
              "Surface: the facets seen by the gauge at that time (NG: Kn > 0, RG: KV > 0 in the fitting base's kernels), "
              "area-weighted (main); the whole surface as the secondary panel for the source rate (P1_whole), at the map "
              "times 10:33:10 and 10:39:26.",
              "P3 matching cells: 10-deg SZA class (plus shadowed, night) x 10 % in distance (bins of a factor 1.1) x 0.1 "
              "in cos theta.",
              "SZA-controlled ratio (the answer tables): for each category, the area-weighted average over the "
              "illumination classes of the per-class ratio of the means (category / 0-40 deg), weights = the category's "
              "area per class; the range in brackets = minimum and maximum over the classes holding at least 1 % of the "
              "category's area.",
              "Under-sampled reference (flag †, an addition): a class in which the 0-40 deg reference holds less "
              "than 1 % of its area on that surface; there the per-class ratio rests on a few 0-40 deg facets and can take "
              "any size. † marks a ratio with more than 10 % of its averaging weight in such classes. The variant "
              "ratio (an addition) averages only over the classes whose reference is not under-sampled; it is a "
              "sensitivity check, not a replacement of the specified ratio.",
              "Code paths: on the gauge-seen surface F, SZA and the shadow state per facet come from the kernel records with "
              "the fit's arithmetic (nonneg_contrib.record_contributions, FP32 records); on the whole surface from the maps' "
              "flux .vtu files (the production-rate maps' code path)."):
        bullet(t, fb=BUILD_FB)
    doc.add_heading("4. Checks", level=1)
    rows = []
    for k, lab in CASES:
        src = os.path.join(root, "data", k, "cliff_data.json")
        c = jl(src)["checks"]
        d = c["rows"]
        fmax = max(d[r].get("F_records_vs_maps_max_rel", 0.0) for r in ("1079", "1549"))
        tmax = max(max(d[r]["ng_vs_contrib_phase_snapshot_rel"], d[r]["rg_vs_contrib_phase_snapshot_rel"]) for r in d)
        eq = all(d[r]["ng_equal_record_contributions"] and d[r]["rg_equal_record_contributions"] for r in d) and \
            all(v for kk, v in c.items() if kk.startswith("phase_"))
        zp = all(d[r]["F_records_vs_maps_zero_pattern_equal"] for r in ("1079", "1549"))
        mis = sum(d[r].get("illumination_class_vs_maps_mismatch_records", 0) for r in ("1079", "1549"))
        rows.append([lab, "yes" if eq else "NO",
                     N.add(k, "F records vs maps max rel", fmax, src, "checks.rows.<1079,1549>.F_records_vs_maps_max_rel",
                           "%.1e"),
                     "yes" if zp else "NO",
                     N.add(k, "illumination class mismatches", mis, src,
                           "checks.rows.<1079,1549>.illumination_class_vs_maps_mismatch_records", "%d"),
                     N.add(k, "row totals vs contribution snapshots max rel", tmax, src,
                           "checks.rows.*.{ng,rg}_vs_contrib_phase_snapshot_rel", "%.1e")])
    table(["case", "contributions equal to the fit's (all rows used)", "F records vs maps, max rel",
           "zero pattern equal", "illumination class mismatches", "row totals vs the full set's snapshots, max rel"],
          rows, [1.5, 1.3, 1.1, 0.9, 1.1, 1.4])
    para("The F differences between the two code paths are of the size of the FP32 rounding of cos SZA in the fit's "
         "records. The row totals are compared with the full sets' contribution-density snapshots (contrib_phase_data.py).",
         size=9, fb=BUILD_FB)
    # ------------------------------------------------------------------------------------------- overall answer
    doc.add_heading("5. Overall answer per case", level=1)
    allw = []
    for k, _ in CASES:
        for d in DEFS:
            a = res(k, "whole", "10:33:10", "F", d)
            if a:
                allw.append((a["ratio"], k, d))
    lo, hi = min(allw), max(allw)
    para("At the same illumination class the source rate depends on the slope in all six cases: on the whole surface at "
         "10:33:10 the SZA-controlled ratios F(category) / F(0-40 deg) range from %s (%s, %s) to %s (%s, %s). Since q and T "
         "are fixed by SZA and the shadow state, this is the slope dependence of the fitted activity: the F and A ratios "
         "agree closely (tables below). Which bins carry it, and how strongly, differs between the models." % (
             N.add(lo[1], "smallest whole-surface F ratio", lo[0], ans_path(lo[1]), "answer.whole.10:33:10.F.%s.ratio" % lo[2],
                   g3), dict(CASES)[lo[1]], lo[2],
             N.add(hi[1], "largest whole-surface F ratio", hi[0], ans_path(hi[1]), "answer.whole.10:33:10.F.%s.ratio" % hi[2],
                   g3), dict(CASES)[hi[1]], hi[2]), bold="Answer.", fb=BUILD_FB)
    para("Table entries: the SZA-controlled ratio (range over the classes); † = more than 10 % of the weight in "
         "classes with an under-sampled 0-40 deg reference. At 10:33:10 the NG and RG gauges see the same facets, so one "
         "gauge-seen column is shown.", size=9, fb=BUILD_FB)
    for k, lab in CASES:
        doc.add_heading(lab, level=2)
        same = jl(ans_path(k))["answer"]["NG_seen"][SNAP0] == jl(ans_path(k))["answer"]["RG_seen"][SNAP0]
        rows = [[d, rcell(k, "NG_seen", SNAP0, "F", d), rcell(k, "whole", "10:33:10", "F", d),
                 rcell(k, "whole", "10:39:26", "F", d), rcell(k, "whole", "10:33:10", "A", d)] for d in DEFS]
        table(["category vs 0-40", "F, gauge-seen 10:33:10" + ("" if same else " (NG)"), "F, whole surface 10:33:10",
               "F, whole surface 10:39:26", "A, whole surface 10:33:10"], rows, [1.3, 1.5, 1.5, 1.5, 1.5])
        best = max(((res(k, "whole", "10:33:10", "F", b) or {"ratio": float("nan")})["ratio"], b) for b in BINS)
        dif = max(abs(res(k, "whole", "10:33:10", "F", d)["ratio"] / res(k, "whole", "10:33:10", "A", d)["ratio"] - 1.0)
                  for d in DEFS if res(k, "whole", "10:33:10", "A", d))
        para("Whole surface, 10:33:10: >40 %s, >60 %s, >80 %s; the largest single-bin ratio is %s (%s). F and A ratios "
             "differ by at most %s %% (the SZA spread inside a 10-deg class)." % (
                 rcell(k, "whole", "10:33:10", "F", ">40", rng=False, flag=False),
                 rcell(k, "whole", "10:33:10", "F", ">60", rng=False, flag=False),
                 rcell(k, "whole", "10:33:10", "F", ">80", rng=False, flag=False),
                 rcell(k, "whole", "10:33:10", "F", best[1], rng=False, flag=False), best[1],
                 N.add(k, "max |F ratio / A ratio - 1| x 100, whole 10:33:10", 100 * dif, ans_path(k),
                       "answer.whole.10:33:10.{F,A}.*.ratio (derived)", "%.2g")), size=9.5, fb=BUILD_FB)
        para("Caveats (stage-2 final comparison summary.txt): " + "; ".join(jl(ans_path(k)).get("caveats", [])) + ".",
             size=8.5, src=[ans_path(k)])
    # ------------------------------------------------------------------------------------------- clusters
    for i, g in enumerate(groups):
        fig_no = i + 1
        N.ctx = g["key"]
        section(g["layout"]["orient"])
        doc.add_heading("%d. %s" % (5 + fig_no, g["title"]), level=1)
        cols, rws, w = g["layout"]["cols"], g["layout"]["rows"], g["layout"]["width_in"]
        t = doc.add_table(rows=rws * 2, cols=cols)
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for j, (k, lab) in enumerate(CASES):
            rr, cc = divmod(j, cols)
            c = t.rows[2 * rr].cells[cc]
            c.text = ""
            r = c.paragraphs[0].add_run("(%s) %s" % ("abcdef"[j], lab))
            r.bold = True
            r.font.size = Pt(8)
            c2 = t.rows[2 * rr + 1].cells[cc]
            c2.text = ""
            c2.paragraphs[0].add_run().add_picture(g["dsts"][j], width=Inches(w))
        cfile = os.path.join(root, CASES[0][0], g["key"] + "_caption.txt")
        capt = open(cfile).read().strip().replace(CASES[0][1], "each case").split(" Caveat: ")[0]
        p = doc.add_paragraph(style="Caption")
        p.add_run(FX.T("Figure %d: %s, all six cases (a-f). " % (fig_no, g["title"]), SRC, BUILD_FB)).bold = True
        p.add_run(FX.T(capt, [cfile] + SRC))
        section("P")
        doc.add_heading("Summary of Figure %d: %s" % (fig_no, g["title"]), level=2)
        summary(g["key"], para, table)
    # ------------------------------------------------------------------------------------------- appendix
    N.ctx = "appendix"
    section("P")
    doc.add_heading("Appendix A. Source figures (SHA-256, first 12 hex digits)", level=1)
    for i, g in enumerate(groups):
        para("Figure %d (%s): " % (i + 1, g["key"]) + "; ".join("%s %s" % (os.path.relpath(s, root), sh[:12])
                                                                 for s, sh in zip(g["srcs"], g["sha"])),
             size=7, fb=("SHA-256 of the source figure (build_manifest.json)", g["srcs"][0]))
    doc.add_heading("Appendix B. How this document was built", level=1)
    para("Data: cliff_effect_data.py per case (Slurm CPU jobs); figures and tables: cliff_effect_figures.py; this "
         "document: cliff_effect_docx.py build --root <folder> in the docx_venv. Every number printed is listed in "
         "numbers_used.json with its source file and locator. LibreOffice is not available on the cluster: the document "
         "was validated by reopening it with python-docx and with the OOXML schema validator; no PDF rendering.", size=8.5)
    if N.unverified:
        with open(os.path.join(root, "numbers_unverified.json"), "x") as f:
            json.dump(N.unverified, f, indent=1, ensure_ascii=False)
        raise SystemExit("%d numeric tokens of hand-written text not found in a source file (numbers_unverified.json); "
                         "document not written" % len(N.unverified))
    doc.save(path)
    with open(os.path.join(root, "numbers_used.json"), "x") as f:
        json.dump(N.items, f, indent=1, default=str, ensure_ascii=False)
    man = {"document": path, "built_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "groups": [dict(key=g["key"], title=g["title"], layout=g["layout"], srcs=g["srcs"], sha256=g["sha"],
                           images=g["dsts"]) for g in groups],
           "images": [[i_[0], list(i_[1]), i_[2]] for i_ in info], "numbers": len(N.items), "command": sys.argv,
           "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "cases": CASES, "cases_json": CASES_JSON}
    with open(os.path.join(root, "build_manifest.json"), "x") as f:
        json.dump(man, f, indent=1)
    print("built %s (%d figures, %d numbers)" % (path, len(groups), len(N.items)))
    return 0


def summary(pk, para, table):
    CASES = Ctx.cases
    labs = [l for _, l in CASES]
    if pk.startswith("P1"):
        if pk == "P1_whole":
            sets = [("whole", "10:33:10"), ("whole", "10:39:26")]
        else:
            g = pk.split("_")[1]
            sets = [("%s_seen" % g, SNAP0), ("%s_seen" % g, SNAP3)]
        for surf, tl in sets:
            para("SZA-controlled F ratio vs 0-40 deg, %s, %s (range over the classes; † under-sampled reference):" % (
                surf.replace("_", "-"), tl), bold="Numbers.", size=9, fb=BUILD_FB)
            table(["category"] + labs, [[d] + [rcell(k, surf, tl, "F", d) for k, _ in CASES] for d in DEFS],
                  [1.0] + [1.05] * 6, size=7)
        if pk != "P1_whole":
            surf, tl = sets[1]
            para("Variant (an addition): only classes with a 0-40 deg reference of at least 1 %% of its area; in brackets "
                 "the fraction of the category's area compared, %s, %s:" % (surf.replace("_", "-"), tl), size=9,
                 fb=BUILD_FB)
            table(["category"] + labs, [[d] + [vcell(k, surf, tl, "F", d) for k, _ in CASES] for d in DEFS],
                  [1.0] + [1.05] * 6, size=7)
        surf, tl = sets[0]
        vals = [((res(k, surf, tl, "F", ">60") or {"ratio": float("nan")})["ratio"], k) for k, _ in CASES]
        hi, lo = max(vals), min(vals)
        para("At %s the >60 deg facets emit %s (%s) to %s (%s) times the 0-40 deg rate at the same illumination class; the "
             "spread between the cases is the spread of the fitted activities. F and A ratios agree (F = A q, q fixed by the "
             "illumination class)." % (tl, rcell(lo[1], surf, tl, "F", ">60", rng=False, flag=False), dict(CASES)[lo[1]],
                                       rcell(hi[1], surf, tl, "F", ">60", rng=False, flag=False), dict(CASES)[hi[1]]),
             bold="Answer for this plot:", size=9.5, fb=BUILD_FB)
        if pk != "P1_whole":
            para("At 20 m and 5.5 m the gauge-seen surface is small and the 0-40 deg reference is under-sampled in some "
                 "classes (†); the late-snapshot ratios are then dominated by those classes and are not a robust "
                 "measure of the slope dependence; the variant and the whole-surface figure are.", size=9.5, fb=BUILD_FB)
    elif pk.startswith("P2"):
        g = pk.split("_")[1]
        rows = []
        fac = [("g", "g"), ("inv_r2", "1/r^2"), ("cos", "cos theta"), ("W", "W"), ("fvis", "visible fraction")]
        dom = {}
        for k, lab in CASES:
            src = os.path.join(Ctx.root, k, "p2_summary.json")
            z = jl(src)
            for d in (">40", ">60", ">80"):
                v = z[g][SNAP0][d]
                rows.append([lab, d] + [N.add(k, "P2 %s factor ratio %s %s" % (g, d, f), v[f], src,
                                              "%s.%s.%s.%s" % (g, SNAP0, d, f), g3) for f, _ in fac])
                if d == ">60":
                    big = max(fac[1:], key=lambda f: abs(math.log10(v[f[0]])) if v[f[0]] > 0 else -1)
                    dom[lab] = big[1]
        para("Ratios of the area-weighted means over all illumination classes, category / 0-40 deg, %s-seen surface at "
             "10:33:10 UTC (the seen facets and their geometry are the same in every case):" % g, bold="Numbers.", size=9,
             fb=BUILD_FB)
        table(["case", "category", "g", "1/r^2", "cos theta", "W", "visible fraction"], rows,
              [1.6, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0], size=7)
        para("The observation weight per m2 of the steep categories differs from that of the 0-40 deg facets mostly "
             "through %s (largest deviation from 1 among the factors for >60 deg at 10:33:10). The visible fraction is "
             "close to 1 for the seen facets: visibility acts by selecting the seen facets, not by weighting them. The "
             "geometry is a property of the shape and the trajectory; it does not enter the source rate F, only the "
             "contribution." % " / ".join(sorted(set(dom.values()))), bold="Answer for this plot:", size=9.5, fb=BUILD_FB)
    elif pk.startswith("P3"):
        g = pk.split("_")[1]
        rows = []
        for k, lab in CASES:
            src = os.path.join(Ctx.root, k, "p3.csv")
            rr = {(r["gauge"], r["definition"], r["entry"]): r for r in csvrows(src)}
            for ent in (SNAP0, LATE):
                vals = []
                for d in (">40", ">60", ">80"):
                    r = rr.get((g, d, ent))
                    if r is None or not math.isfinite(float(r["weighted_median_ratio"])):
                        vals.append("-")
                        continue
                    vals.append("%s (%s)" % (
                        N.add(k, "P3 weighted median ratio %s %s %s" % (g, d, ent), float(r["weighted_median_ratio"]),
                              src, "%s,%s,%s,weighted_median_ratio" % (g, d, ent), g3),
                        N.add(k, "P3 matched fraction %s %s %s" % (g, d, ent), float(r["contribution_fraction_matched"]),
                              src, "%s,%s,%s,contribution_fraction_matched" % (g, d, ent), "%.2f")))
                rows.append([lab, ent] + vals)
        para("Contribution-weighted median of the matched-cell ratio of the %s contribution per m2 (category / 0-40 deg; in "
             "brackets the fraction of the category's contribution in matched cells):" % g, bold="Numbers.", size=9,
             fb=BUILD_FB)
        table(["case", "time / phase", ">40", ">60", ">80"], rows, [1.7, 1.7, 1.2, 1.2, 1.2], size=7)
        para("With the illumination class, the distance and the emission angle matched, the remaining ratio is the activity "
             "ratio times the residual gauge-weight and visibility differences. Medians above 1 mean that, cell by cell, a "
             "steep facet contributes more per m2 than a flat facet at the same illumination and viewing geometry; the "
             "matched fraction says how much of the steep facets' contribution has a flat counterpart at all.",
             bold="Answer for this plot:", size=9.5, fb=BUILD_FB)
    else:
        g = pk.split("_")[1]
        rows = []
        for k, lab in CASES:
            src = os.path.join(Ctx.root, k, "p4.csv")
            rr = {(r["gauge"], r["entry"], r["category"]): r for r in csvrows(src)}
            for ent in (SNAP0, LATE):
                r = rr.get((g, ent, ">60"))
                if not r:
                    continue
                rows.append([lab, ent] + [N.add(k, "P4 %s >60 %s %s" % (g, ent, f), float(r[f]), src,
                                                "%s,%s,>60,%s" % (g, ent, f), g3)
                                          for f in ("share_over_area", "activity", "illumination", "geometry", "covariance",
                                                    "product_check")])
        para("Slope > 60 deg, %s: signal share / area share and its factors (check = activity x illumination x geometry x "
             "covariance / (share/area), 1 by construction):" % g, bold="Numbers.", size=9, fb=BUILD_FB)
        table(["case", "time / phase", "share/area", "activity", "illumination", "geometry", "covariance", "check"], rows,
              [1.5, 1.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.6], size=7)
        para("The geometry factor of a category depends only on the shape and the trajectory; the illumination factor "
             "depends on the shape, the trajectory and the temperature law, which is the same in all six fits; the activity "
             "factor and the covariance carry the differences between the models. Whether steep facets are over- or "
             "under-represented in the signal is therefore decided by the fitted activity.",
             bold="Answer for this plot:", size=9.5, fb=BUILD_FB)


def validate(args):
    from docx import Document
    root = os.path.abspath(args.root)
    doc = Document(os.path.join(root, DOC_NAME))
    caps = [p.text for p in doc.paragraphs if p.style is not None and p.style.name == "Caption"
            and re.match(r"^Figure \d+: ", p.text)]
    shapes = doc.inline_shapes
    man = json.load(open(os.path.join(root, "build_manifest.json")))
    problems = []
    k = 0
    for gi, g in enumerate(man["groups"]):
        lay = g["layout"]
        W, H = FX.usable(lay["orient"])
        for j in range(6):
            sh = shapes[k]
            w_in, h_in = sh.width / 914400, sh.height / 914400
            if w_in > W / lay["cols"] + 1e-6:
                problems.append("figure %d image %d wider than its cell (%.2f > %.2f in)" % (gi + 1, j + 1, w_in,
                                                                                            W / lay["cols"]))
            if (h_in + FX.LABEL) * lay["rows"] > H - FX.OVERHEAD + 0.3:
                problems.append("figure %d image %d too tall (%.2f in x %d rows)" % (gi + 1, j + 1, h_in, lay["rows"]))
            k += 1
    heads = {}
    for p in doc.paragraphs:
        sn = p.style.name if p.style is not None else ""
        if sn.startswith("Heading") or sn == "Title":
            heads[sn] = heads.get(sn, 0) + 1
    res_ = {"document": os.path.join(root, DOC_NAME), "captions": len(caps), "expected_captions": len(PLOTS),
            "images": len(shapes), "expected_images": 6 * len(PLOTS), "tables": len(doc.tables),
            "sections": len(doc.sections), "headings": heads, "problems": problems,
            "numbers_registered": len(json.load(open(os.path.join(root, "numbers_used.json")))),
            "size_bytes": os.path.getsize(os.path.join(root, DOC_NAME))}
    res_["passed"] = bool(res_["captions"] == len(PLOTS) and res_["images"] == 6 * len(PLOTS) and not problems)
    os.makedirs(os.path.join(root, "validation"), exist_ok=True)
    # layout previews of the clusters (drawn from the document's geometry at 60 dpi; NOT a Word rendering)
    from PIL import Image, ImageDraw
    k = 0
    for gi, g in enumerate(man["groups"]):
        lay = g["layout"]
        pw, ph = FX.PAGE[lay["orient"]]
        dpi = 60
        im = Image.new("RGB", (int(pw * dpi), int(ph * dpi)), "white")
        dr = ImageDraw.Draw(im)
        x0, y0 = FX.MARGIN * dpi, (FX.MARGIN + 0.45) * dpi
        cw = FX.usable(lay["orient"])[0] / lay["cols"] * dpi
        dr.text((x0, FX.MARGIN * dpi), "Figure %d: %s" % (gi + 1, g["title"]), fill="black")
        yrow = y0
        for j in range(6):
            sh = shapes[k + j]
            w_px, h_px = int(sh.width / 914400 * dpi), int(sh.height / 914400 * dpi)
            rr, cc = divmod(j, lay["cols"])
            if cc == 0 and j > 0:
                yrow += prev_h + FX.LABEL * dpi
            if cc == 0:
                prev_h = 0
            pic = Image.open(g["images"][j]).convert("RGB").resize((max(1, w_px), max(1, h_px)))
            im.paste(pic, (int(x0 + cc * cw), int(yrow + FX.LABEL * dpi)))
            prev_h = max(prev_h, h_px)
        k += 6
        dr.rectangle([FX.MARGIN * dpi, FX.MARGIN * dpi, (pw - FX.MARGIN) * dpi, (ph - FX.MARGIN) * dpi], outline="red")
        im.save(os.path.join(root, "validation", "layout_preview_fig%02d.png" % (gi + 1)))
    with open(os.path.join(root, "validation", "validation.json"), "w") as f:
        json.dump(res_, f, indent=1)
    print(json.dumps(res_, indent=1))
    return 0 if res_["passed"] else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--root", required=True)
    b.add_argument("--jobs", type=int, default=8)
    v = sub.add_parser("validate")
    v.add_argument("--root", required=True)
    a = ap.parse_args()
    return build(a) if a.cmd == "build" else validate(a)


if __name__ == "__main__":
    sys.exit(main())
