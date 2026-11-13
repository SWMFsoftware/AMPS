#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/collect_captions.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/collect_captions.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e64ce7819f658e1266590340a6694889b23ef10d0da5db560f3cf3d584b5298a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure helper: all <figure>_caption.txt of a folder into one CAPTIONS.md.
# RUN: login node: collect_captions.py <folder> [--order ...] [--skip ...] [--title ...]
# ==== end of header ====
"""Collect the per-figure captions (<figure>_caption.txt) of a figure folder into one Markdown file (user request
2026-11-12: "add captions to all plots in a separate file").

    collect_captions.py <figure folder> [--out CAPTIONS.md] [--title "..."]

Every <name>_caption.txt below the folder (skipping superseded*, *.incomplete_*, provenance-*, code*) becomes a section
with the figure files (<name>.eps, <name>.png) and the caption text verbatim. Sections are ordered by folder (in the
order of --order, then alphabetically) and file name. Refuses to overwrite an existing output file.
"""
import argparse
import sys
from pathlib import Path

SKIP = ("superseded", ".incomplete_", "provenance-", "code")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--out", default="CAPTIONS.md")
    ap.add_argument("--title", default="Figure captions")
    ap.add_argument("--order", default="activity,production_rate,source_areas,cliffs")
    ap.add_argument("--skip", default="", help="comma list of further path parts to skip (e.g. a preview folder)")
    a = ap.parse_args()
    root = Path(a.folder)
    out = root / a.out
    if out.exists():
        sys.exit("%s exists" % out)
    order = a.order.split(",")
    caps = []
    for p in root.rglob("*_caption.txt"):
        rel = p.relative_to(root)
        if any(s in part for part in rel.parts for s in SKIP + tuple(x for x in a.skip.split(",") if x)):
            continue
        top = rel.parts[0] if len(rel.parts) > 1 else ""
        caps.append(((order.index(top) if top in order else len(order), top, str(rel)), rel, p))
    caps.sort()
    lines = ["# %s" % a.title, "", "Folder: `%s`. %d figures; each also has its own `<figure>_caption.txt` next to it." % (
        root, len(caps)), ""]
    section = None
    for (_, top, _), rel, p in caps:
        if top != section:
            lines += ["## %s" % (top or "(top level)"), ""]
            section = top
        stem = str(rel)[:-len("_caption.txt")]
        lines += ["### `%s`" % stem, "", "Files: `%s.eps`, `%s.png`" % (stem, stem), "",
                  p.read_text(encoding="utf-8").strip(), ""]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("%s: %d captions" % (out, len(caps)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
