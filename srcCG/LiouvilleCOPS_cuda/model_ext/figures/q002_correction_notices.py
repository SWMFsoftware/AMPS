#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_correction_notices.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_correction_notices.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 658663544338c500c5ccb494f79c474ed90fa9556789412d680e9fec02d6118a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: correction notices on older documents and figure sets (decision 2b)
# RUN: jobs/q002_correction_notices.sbatch
# ==== end of header ====
"""Decision 2b of 2026-11-21 (FINDINGS A21): correction notices on the older documents and figure sets.

* the three older Word documents: a notice paragraph after the title. The original is kept as *.superseded_<UTC>
  (byte-identical, the file its build_manifest.json describes); text and images are otherwise unchanged (checked).
* a CORRECTION_20261121.md next to those documents, in the ten full sets of the round, in source_rate_q002, in the
  old-law full sets and in the five-fit report's slope_activity folder. No manifest, facts.json, figure or table is
  edited.

    q002_correction_notices.py --root <figures_cluster/12M> --stamp <UTC>
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import docx
from docx.shared import RGBColor

NOTE = "CORRECTION_20261121.md"
HEAD = "Correction notice (added 2026-11-21; user decision 2b, FINDINGS A21). "
ZERO = ("'A = 0' and 'A > 0' shares in %s do not measure switched-off surface: option C's bounded solver leaves "
        "switched-off activity at tiny positive values, and option B keeps an imposed background. ")
BOUND = ("The slope multipliers of these models are bounded below by 1, so nothing here shows that steep terrain is as "
         "active as flat terrain, or that active cliffs are unnecessary: that was not tested. ")
MASKS = ("The late-visibility templates (NG-only / RG-only facets, 5,924 / 29,210) belong to an earlier base with another "
         "gauge position and do not show what each gauge saw in these fits (FINDINGS A20). ")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def texts(d):
    return [p.text for p in d.paragraphs], [c.text for t in d.tables for r in t.rows for c in r.cells], len(d.inline_shapes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--stamp", required=True)
    a = ap.parse_args()
    R = Path(a.root)
    rows = list(csv.DictReader(open(R / "q002_decisions_20261121/active_area.csv")))
    table = ["| run | switched off (decision 1b) | with the floor counted | result |", "|---|---|---|---|"]
    for r in rows:
        table.append("| %s | %.1f %% | %.1f %% | %s |" % (r["run"], 100 * float(r["off_T"]),
                     100 * float(r["off_T_floor_counted"]), "pass" if r["passes"] == "True" else "fail"))
    new_law = ("Share of the seen area switched off in the new-law fits (active = activity without an imposed floor above "
               "T = 2.114e15, 1 % of the median srcCG activity of the seen facets; more than 50 % off fails):")
    cur = ("The current documents are q002_all_fits_summary/12M_q002_all_fits_summary.docx and "
           "q002_top5_models/12M_q002_top5_models.docx (figures_cluster/12M). ")
    docs = {
        "stage2_q002_comparison/12M_stage2_q002_comparison.docx":
            (ZERO % "this document" + "With the measure now decided the option-C fits switch off 66-94 % of the seen "
             "surface and B 96 %; only A and D pass. " + BOUND + "This document is superseded. " + cur, True),
        "full_sets_stage2_comparison/12M_stage2_full_sets_comparison.docx":
            (ZERO % "this document" + "These are old-law fits (q_min = 0.2 q_max); their switched-off shares were not "
             "recomputed. " + BOUND + MASKS + "For the current law: " + cur, False),
        "cliff_effect_stage2/12M_stage2_cliff_effect.docx":
            (ZERO % "this document" + "These are old-law fits (q_min = 0.2 q_max); their switched-off shares were not "
             "recomputed. " + BOUND + "Within about 150 m the models' only way to vary activity with slope is the fitted "
             "multiplier, so the slope dependence shown here is not independent evidence (FINDINGS A20). " + cur, False),
    }
    record = {}
    for rel in docs:                               # nothing is written unless every target is as expected
        f = R / rel
        bak = f.with_name(f.name + ".superseded_" + a.stamp)
        if bak.exists() or not f.exists() or (f.parent / NOTE).exists():
            raise SystemExit("unexpected state of %s" % f.parent)
    for rel, (text, with_table) in docs.items():
        f = R / rel
        bak = f.with_name(f.name + ".superseded_" + a.stamp)
        tmp = f.with_name(f.stem + ".notice_tmp.docx")
        d = docx.Document(str(f))
        before = texts(d)
        ps = d.paragraphs
        title = next(i for i, p in enumerate(ps) if p.text.strip())
        p = ps[title + 1].insert_paragraph_before("", style="Normal")
        for t, bold in ((HEAD, True), (text.strip(), False)):
            run = p.add_run(t)
            run.bold = bold
            run.font.color.rgb = RGBColor(0xA0, 0x1C, 0x1C)
        d.save(str(tmp))
        after = texts(docx.Document(str(tmp)))
        i = title + 1
        if after[0][:i] + after[0][i + 1:] != before[0] or after[1:] != before[1:]:
            raise SystemExit("%s: more than the notice changed" % f)
        h0 = sha(f)
        f.rename(bak)
        tmp.rename(f)
        rec = dict(original=bak.name, original_sha256=h0, with_notice_sha256=sha(f), images=before[2])
        record[rel] = rec
        L = ["# Correction, 2026-11-21 (user decision 2b; FINDINGS A21)", "", text.strip(), ""]
        if with_table:
            L += [new_law, ""] + table + [""]
        L += ["The document now carries this notice after its title; nothing else in it changed. The file that "
              "build_manifest.json describes is `%s` (SHA-256 `%s`); with the notice: `%s`." % (
                  bak.name, h0, rec["with_notice_sha256"]), ""]
        (f.parent / NOTE).write_text("\n".join(L))
        print("notice in", rel, "images", before[2])

    q002 = ["# Correction, 2026-11-21 (user decision 2b; FINDINGS A21)", "",
            (ZERO % "this folder (facts, inspection files, captions, README)").strip(),
            "Light green in the maps means exactly zero only: switched-off terrain of the C fits and B's imposed "
            "background are drawn in colour. The maps, tables and facts here are unchanged (regenerating them was "
            "option 2c, not chosen).", "", new_law, ""] + table + ["", BOUND.strip(), "", MASKS.strip(), "",
            "The statement that the NG overhang share at touchdown is 'not required by the data' goes too far: the fit "
            "without templates has none and also misses the NG rise of the last ~20 m.", "",
            "Numbers: `figures_cluster/12M/q002_decisions_20261121/active_area.json`.", ""]
    sets = sorted(p for p in R.glob("full_set_q002_*") if p.is_dir() and ".superseded" not in p.name
                  and ".incomplete" not in p.name) + [R / "source_rate_q002"]
    old = sorted(p for p in R.glob("full_set_stage2*") if p.is_dir() and ".superseded" not in p.name
                 and ".incomplete" not in p.name)
    old.append(R / "full_set_D2")                     # the old-law document's sixth set
    oldtext = ["# Correction, 2026-11-21 (user decision 2b; FINDINGS A21)", "",
               (ZERO % "this folder").strip(),
               "These are old-law fits (q_min = 0.2 q_max); their switched-off shares were not recomputed.", "",
               BOUND.strip(), "", MASKS.strip(), ""]
    slope = ["# Correction, 2026-11-21 (FINDINGS A21)", "",
             "The switched-off numbers in this folder (shutdown_check/, the off_1e-6 / off_1pct entries of pass2/ and "
             "matched/, the columns of slope_statistics.csv) use A < k A_u, which decision 1b replaced. Current shares, "
             "also per slope class: ../../q002_decisions_20261121/active_area.json.", "",
             "Matched ratios of 1.00-1.07 at 60-80 deg are the multipliers' imposed lower bound 1 or within 7 % of it, "
             "not a measurement that such cliffs are as active as flat terrain. The other statistics are unaffected.", ""]
    for folder, L in [(p, q002) for p in sets] + [(p, oldtext) for p in old] + [(R / "q002_top5_models/slope_activity", slope)]:
        n = folder / NOTE
        if n.exists():
            raise SystemExit("%s exists" % n)
        n.write_text("\n".join(L))
        record[str(n.relative_to(R))] = dict(sha256=sha(n))
        print("note", n.relative_to(R))
    out = R / "q002_decisions_20261121/correction_notices.json"
    with open(out, "x") as f:
        json.dump(record, f, indent=1)


if __name__ == "__main__":
    main()
