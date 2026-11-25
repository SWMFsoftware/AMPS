#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_top5_docx.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_top5_docx.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f625cbaac511065fc32c545180e5efb05156accbaf2b092a7660be1b756ec509
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Word document 12M_q002_top5_models.docx on the five most plausible new-law fits (summary, plots, maps,
#   section 4 comet and sources, methods) from numbers*.json, the figures, the full sets and q002_top5_content.py (2026-11-21)
# RUN: docx venv python q002_top5_docx.py --out figures_cluster/12M/q002_top5_models (jobs/q002_top5_docx.sbatch)
# CHANGED: 2026-11-21 user decisions 1b / 2b (FINDINGS A21)
# ==== end of header ====
"""Word document on the five most plausible new-law fits (2026-11-21): summary, plots, how the models differ
numerically, their assumptions, and what the fit minimizes.

User 2026-11-21: "prepare a summary word doc for the 5 models in the table; in addition to the summary plots, include a
detailed description of how those model numerically are difference from each other, what are ther assumptions, what the
fitting procedure is actuly minimizes" (the 5 = ranks 1-5 of 12M/q002_plausibility, user's choice).

Inputs (read only): <out>/numbers.json and <out>/figures/ (q002_top5_figures.py), the five full sets
12M/full_set_q002_<run>/ (activity and production-rate maps), 12M/q002_plausibility/metrics.json, and the method text
in q002_top5_content.py (next to this script; every statement there carries its source).
Writes <out>/12M_q002_top5_models.docx, <out>/images/ (the downscaled copies in the document), <out>/build_manifest.json
(inputs with SHA-256). Never overwrites.

    <docx venv>/python q002_top5_docx.py --out <dir>
"""
import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from PIL import Image
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import q002_top5_content as C  # noqa: E402

F12 = Path("/data001/highresolu_vt/cops_output/figures_cluster/12M")
DOC = "12M_q002_top5_models.docx"
FS = {"C_f10_pm0p01_q002": "C_f10_pm0p01", "C_f1_pm0p01_q002": "C_f1_pm0p01",
      "C_f10_pm0p01_tploff_q002": "C_f10_pm0p01_tploff", "A_q002": "A", "B_q002": "B"}
ORDER = list(FS)
INPUTS = {}
FIGN = [0]
TABN = [0]
PICS = [0]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def note_input(p):
    INPUTS[str(p)] = sha256(p)


# ------------------------------------------------------------------ document helpers
def set_cell_shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hexcolor)
    tcPr.append(shd)


def para(doc, text="", style=None, bold=False, italic=False, size=None, align=None, after=4, keep=False):
    p = doc.add_paragraph(style=style)
    if text:
        add_runs(p, text, bold=bold, italic=italic, size=size)
    if align:
        p.alignment = align
    p.paragraph_format.space_after = Pt(after)
    if keep:
        p.paragraph_format.keep_with_next = True
    return p


def add_runs(p, text, bold=False, italic=False, size=None):
    """Text with **bold** segments and _{sub} / ^{sup} markers."""
    import re
    for part in re.split(r"(\*\*[^*]+\*\*|_\{[^}]*\}|\^\{[^}]*\}|\^-?\d+)", text):
        if not part:
            continue
        r = None
        if part.startswith("**"):
            r = p.add_run(part[2:-2]); r.bold = True
        elif part.startswith("_{"):
            r = p.add_run(part[2:-1]); r.font.subscript = True
        elif part.startswith("^{"):
            r = p.add_run(part[2:-1]); r.font.superscript = True
        elif part.startswith("^"):
            r = p.add_run(part[1:]); r.font.superscript = True
        else:
            r = p.add_run(part)
        if bold:
            r.bold = True
        if italic:
            r.italic = True
        if size:
            r.font.size = Pt(size)
    return p


def bullet(doc, text, level=0):
    p = doc.add_paragraph(style="List Bullet" if level == 0 else "List Bullet 2")
    add_runs(p, text)
    p.paragraph_format.space_after = Pt(2)
    return p


def equation(doc, text, label=None):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_runs(p, text)
    for r in p.runs:
        r.font.name = "Cambria Math"
        r._element.rPr.rFonts.set(qn("w:eastAsia"), "Cambria Math")
        r.font.size = Pt(10.5)
    if label:
        r = p.add_run("    (%s)" % label)
        r.font.size = Pt(10)
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(4)
    return p


def table(doc, header, rows, widths_cm, font=8.5, caption=None):
    if caption:
        TABN[0] += 1
        para(doc, "Table %d. %s" % (TABN[0], caption), italic=True, size=9, keep=True, after=2)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        add_runs(c.paragraphs[0], h, bold=True, size=font)
        set_cell_shade(c, "E7ECF2")
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            add_runs(cells[i].paragraphs[0], str(v), size=font)
    for row in t.rows:
        for i, w in enumerate(widths_cm):
            row.cells[i].width = Cm(w)
    # repeat the header row on every page
    trPr = t.rows[0]._tr.get_or_add_trPr()
    h = OxmlElement("w:tblHeader")
    h.set(qn("w:val"), "true")
    trPr.append(h)
    para(doc, "", after=4)
    return t


def image(doc, src, out_dir, width_cm, caption, max_px=1700, jpeg=False):
    note_input(src)
    im = Image.open(src)
    im.load()
    if im.mode in ("RGBA", "P", "LA"):
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im.convert("RGBA"), mask=im.convert("RGBA").split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    scale = min(1.0, max_px / max(im.size))
    if scale < 1:
        im = im.resize((round(im.width * scale), round(im.height * scale)), Image.LANCZOS)
    src_p = Path(src)
    prefix = next((q.name.replace("full_set_q002_", "") + "_" for q in src_p.parents if q.name.startswith("full_set_q002_")), "")
    dst = out_dir / (prefix + src_p.stem + (".jpg" if jpeg else ".png"))
    if dst.exists():
        raise SystemExit("%s exists" % dst)
    if jpeg:
        im.save(dst, "JPEG", quality=82, optimize=True)
    else:
        im.save(dst, "PNG", optimize=True)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(str(dst), width=Cm(width_cm))
    PICS[0] += 1
    FIGN[0] += 1
    c = para(doc, "Figure %d. %s" % (FIGN[0], caption), italic=True, size=8.5, after=8)
    return FIGN[0]


def heading(doc, text, level):
    h = doc.add_heading(text, level=level)
    h.paragraph_format.keep_with_next = True
    return h


def page_break(doc):
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


def fmt(v, d=3):
    return ("%%.%dg" % d) % v


def sci(v, d=2):
    m, e = ("%%.%de" % d % v).split("e")
    return "%s×10^{%d}" % (m, int(e))


# ------------------------------------------------------------------ build
def build(out):
    num = json.load(open(out / "numbers.json"))
    note_input(out / "numbers.json")
    R = num["runs"]
    FIG = out / "figures"
    IMG = out / "images"
    IMG.mkdir()

    doc = Document()
    # python-docx's template has <w:zoom w:val="bestFit"/> without the w:percent the OOXML schema requires
    zoom = doc.settings.element.find(qn("w:zoom"))
    if zoom is not None and zoom.get(qn("w:percent")) is None:
        zoom.set(qn("w:percent"), "100")
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(10.5)
    for lvl, size in ((1, 15), (2, 12.5), (3, 11)):
        s = doc.styles["Heading %d" % lvl]
        s.font.name = "Calibri"
        s.font.size = Pt(size)
        s.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    for side in ("left_margin", "right_margin"):
        setattr(sec, side, Cm(2.0))
    sec.top_margin = sec.bottom_margin = Cm(2.0)
    W = 17.0                                                     # text width, cm

    # title
    para(doc, C.TITLE, bold=True, size=17, after=6)
    para(doc, C.SUBTITLE, italic=True, size=10.5, after=10)
    for line in C.FRONT:
        para(doc, line, size=10, after=3)

    heading(doc, "1. Summary", 1)
    NS = json.load(open(out / "numbers_sources.json"))
    note_input(out / "numbers_sources.json")
    summ = C.summary(R)
    summ = summ[:-1] + C.summary_sources(NS) + summ[-1:]          # before the caveats bullet
    for blk in summ:
        (bullet(doc, blk[1]) if blk[0] == "b" else para(doc, blk[1]))

    heading(doc, "2. The five fits at a glance", 1)
    hdr, rows, widths, cap = C.glance_table(R, ORDER)
    table(doc, hdr, rows, widths, font=7.8, caption=cap)
    hdr, rows, widths, cap = C.settings_table(R, ORDER)
    table(doc, hdr, rows, widths, font=7.8, caption=cap)

    heading(doc, "3. Summary plots", 1)
    para(doc, C.PLOTS_INTRO)
    for stem, width, cap in C.summary_figures(R):
        image(doc, FIG / (stem + ".png"), IMG, width, cap)
    heading(doc, "3.6 Activity and production-rate maps of each fit", 2)
    para(doc, C.MAPS_INTRO)
    for k in ORDER:
        fs = F12 / ("full_set_q002_%s" % FS[k])
        heading(doc, "%s" % R[k]["label"], 3)
        t = doc.add_table(rows=1, cols=2)
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for j, (rel, what) in enumerate((("activity/12M_groundtrack_activity_latlon.png", "activity"),
                                         ("production_rate/12M_groundtrack_h2o_production_rate_103310UTC_latlon.png",
                                          "production"))):
            src = fs / rel
            note_input(src)
            im = Image.open(src).convert("RGB")
            s = min(1.0, 1300 / max(im.size))
            im = im.resize((round(im.width * s), round(im.height * s)), Image.LANCZOS)
            dst = IMG / ("%s_%s.jpg" % (FS[k], what))
            im.save(dst, "JPEG", quality=82, optimize=True)
            cell = t.rows[0].cells[j]
            cell.width = Cm(W / 2)
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
            cell.paragraphs[0].add_run().add_picture(str(dst), width=Cm(W / 2 - 0.4))
            PICS[0] += 1
        FIGN[0] += 1
        para(doc, "Figure %d. %s" % (FIGN[0], C.map_caption(R[k], fs)), italic=True, size=8.5, after=8)

    def render(blocks):
        for b in blocks:
            kind = b[0]
            if kind == "p":
                para(doc, b[1])
            elif kind == "b":
                bullet(doc, b[1])
            elif kind == "b2":
                bullet(doc, b[1], level=1)
            elif kind == "eq":
                equation(doc, b[1], b[2] if len(b) > 2 else None)
            elif kind in ("h2", "h3"):
                heading(doc, b[1], 2)
            elif kind == "table":
                table(doc, b[1], b[2], b[3], font=b[4] if len(b) > 4 else 8, caption=b[5] if len(b) > 5 else None)
            elif kind == "fig":
                src = Path(b[1]) if Path(b[1]).is_absolute() else out / b[1]
                cap = b[3]
                if cap is None:
                    cap = src.with_name(src.stem + "_caption.txt").read_text().strip().replace("five most plausible fits", "five selected fits")
                    note_input(src.with_name(src.stem + "_caption.txt"))
                image(doc, src, IMG, b[2], cap, max_px=1800, jpeg=True)
            elif kind == "grid":
                items, wcm, cap = b[1], b[2], b[3]
                t = doc.add_table(rows=(len(items) + 1) // 2, cols=2)
                t.alignment = WD_TABLE_ALIGNMENT.CENTER
                for n, (src, label) in enumerate(items):
                    note_input(src)
                    im = Image.open(src).convert("RGB")
                    sc = min(1.0, 1400 / max(im.size))
                    im = im.resize((round(im.width * sc), round(im.height * sc)), Image.LANCZOS)
                    pre = next((q.name.replace("full_set_q002_", "") for q in Path(src).parents
                                if q.name.startswith("full_set_q002_")), "x")
                    dst = IMG / ("%s_%s.jpg" % (pre, Path(src).stem))
                    if dst.exists():
                        raise SystemExit("%s exists" % dst)
                    im.save(dst, "JPEG", quality=82, optimize=True)
                    cell = t.rows[n // 2].cells[n % 2]
                    cell.width = Cm(W / 2)
                    cp = cell.paragraphs[0]
                    cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    add_runs(cp, "(%s) %s" % ("abcdefgh"[n], label), bold=True, size=8)
                    pp = cell.add_paragraph()
                    pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    pp.add_run().add_picture(str(dst), width=Cm(wcm))
                    PICS[0] += 1
                FIGN[0] += 1
                para(doc, "Figure %d. %s" % (FIGN[0], cap), italic=True, size=8.5, after=8)
            else:
                raise SystemExit("unknown block %r" % kind)

    page_break(doc)
    heading(doc, "4. The comet and its sources as seen by the five fits", 1)
    para(doc, C.SOURCES_SECTION_INTRO)
    render(C.sources(R, NS, F12))
    sa = out / "slope_activity"
    for q in ("pass2/slope_activity.json", "matched/slope_matched.json", "seen_gauges.npz",
              "matched/slope_statistics.csv"):
        note_input(sa / q)
    aa = F12 / "q002_decisions_20261121/active_area.json"
    note_input(aa)
    render(C.slope_section(json.load(open(sa / "pass2/slope_activity.json")),
                           json.load(open(sa / "matched/slope_matched.json")), out, json.load(open(aa))))

    page_break(doc)
    ms_rep = Path('/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_q002_f10_pm0.01_tploff_20261122T044019Z/report_137981/slopetest.json')
    ms_sta = Path('/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/multistart_q002_f10_pm0.01_tploff_20261122T065041Z/stationarity_137980/stationarity.json')
    note_input(ms_rep)
    note_input(ms_sta)
    last = C.multistart_section(json.load(open(ms_rep)), json.load(open(ms_sta)))
    for sec_title, blocks in C.methods(R, ORDER, num["ten_fits"]) + [last]:
        heading(doc, sec_title[0], sec_title[1])
        render(blocks)

    heading(doc, "Appendix: sources", 1)
    para(doc, C.SOURCES_INTRO)
    for p_, h in sorted(INPUTS.items()):
        para(doc, "%s  (SHA-256 %s)" % (p_, h[:16] + "…"), size=7.5, after=1)
    for s in C.EXTRA_SOURCES:
        para(doc, s, size=7.5, after=1)

    dst = out / DOC
    if dst.exists():
        raise SystemExit("%s exists" % dst)
    doc.save(dst)
    man = {"created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
           "generator": str(Path(__file__).resolve()), "generator_sha256": sha256(Path(__file__)),
           "content": str(Path(C.__file__).resolve()), "content_sha256": sha256(C.__file__),
           "figures": FIGN[0], "tables": TABN[0], "inputs": INPUTS}
    (out / "build_manifest.json").write_text(json.dumps(man, indent=1))
    print("wrote", dst, "figures", FIGN[0], "tables", TABN[0], "size %.1f MB" % (dst.stat().st_size / 1e6))


def validate(out):
    d = Document(out / DOC)
    caps = [p.text for p in d.paragraphs if p.text.startswith("Figure ")]
    pics = len(d.inline_shapes)
    print("reopened: %d paragraphs, %d tables, %d inline images, %d figure captions" % (
        len(d.paragraphs), len(d.tables), pics, len(caps)))
    if len(caps) != FIGN[0] or pics != PICS[0]:
        raise SystemExit("validation: %d captions (expected %d), %d images (expected %d)" % (
            len(caps), FIGN[0], pics, PICS[0]))
    return {"paragraphs": len(d.paragraphs), "tables": len(d.tables), "images": pics, "figure_captions": len(caps)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    build(out)
    v = validate(out)
    (out / "validation.json").write_text(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
