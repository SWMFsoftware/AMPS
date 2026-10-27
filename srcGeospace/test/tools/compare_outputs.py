#!/usr/bin/env python3
"""Compare the numerical outputs of two test run roots (made by run_srcearth_tests.sh).

Usage:
    compare_outputs.py REF_ROOT NEW_ROOT [--out report.csv] [--subdir test_output]

For every file below REF_ROOT/test_output that also exists below NEW_ROOT/test_output
the script reports one of

    IDENTICAL   byte-for-byte equal
    IDENTICAL_AFTER_PATHS  equal after replacing run-root paths and the application
                directory name (srcEarth/srcGeospace) by placeholders
    NONDETERMINISTIC  known to differ between two runs of the unmodified reference
    NUMERIC     same token structure; reports the maximum absolute and relative
                difference over all numeric tokens (text files only)
    DIFFERENT   different token structure or binary content
    MISSING     present only in the reference root

Plots (png/eps/pdf), logs and JSON/text files that only embed absolute paths or
timings are skipped by default because they always differ between run roots.
The comparison never modifies either root.
"""

import argparse
import csv
import math
import os
import re
import sys

SKIP_SUFFIXES = (".png", ".eps", ".pdf", ".log", ".json", ".gz", ".npz", ".pkl")
# Files that differ between two runs of the SAME unmodified reference executable
# (measured: AMPS writes the Mode3D AMR mesh cache "mesh->msh"/"*.mesh->bin" with
# run-dependent content).  They are reported as NONDETERMINISTIC, not as differences.
NONDETERMINISTIC_SUBSTRINGS = ("mesh->",)
TOKEN = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?|[^\s,;]+")


def _float(token):
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


ROOT_SUBST = []   # (absolute run-root path, placeholder) pairs, set in main()


# Any absolute path prefix in front of a run root's test_output/ or application link,
# so that roots that were moved after the run still compare equal.
ROOT_PREFIX = re.compile(r"(?:/[^\s,;\"'<>|]+)?/(test_output|srcEarth|srcGeospace)/")
# The application directory was renamed from srcEarth to srcGeospace; outputs of the
# two trees name it differently in commands and manifests.
APP_NAME = re.compile(r"\bsrc(?:Earth|Geospace)\b")


def normalize(text):
    # Run-root paths embedded in outputs (manifests, summaries, generated decks) differ
    # between roots by construction; replace them before comparing.
    for root, tag in ROOT_SUBST:
        text = text.replace(root, tag)
    text = ROOT_PREFIX.sub(r"<RUN_ROOT>/\1/", text)
    return APP_NAME.sub("<APP>", text)


def compare_text(a_path, b_path):
    with open(a_path, "r", errors="replace") as fa, open(b_path, "r", errors="replace") as fb:
        a_text = normalize(fa.read())
        b_text = normalize(fb.read())
    if a_text == b_text:
        return "IDENTICAL_AFTER_PATHS", 0.0, 0.0, "differs only in run-root paths / application name"
    a_tokens = TOKEN.findall(a_text)
    b_tokens = TOKEN.findall(b_text)
    if len(a_tokens) != len(b_tokens):
        return "DIFFERENT", float("nan"), float("nan"), f"token count {len(a_tokens)} vs {len(b_tokens)}"
    max_abs = 0.0
    max_rel = 0.0
    n_numeric_diff = 0
    for ta, tb in zip(a_tokens, b_tokens):
        if ta == tb:
            continue
        fa_, fb_ = _float(ta), _float(tb)
        if fa_ is None or fb_ is None:
            return "DIFFERENT", float("nan"), float("nan"), f"token '{ta}' vs '{tb}'"
        if math.isnan(fa_) and math.isnan(fb_):
            continue
        n_numeric_diff += 1
        d = abs(fa_ - fb_)
        scale = max(abs(fa_), abs(fb_))
        max_abs = max(max_abs, d)
        if scale > 0.0:
            max_rel = max(max_rel, d / scale)
    return "NUMERIC", max_abs, max_rel, f"{n_numeric_diff} numeric tokens differ"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref_root")
    ap.add_argument("new_root")
    ap.add_argument("--subdir", default="test_output")
    ap.add_argument("--out", default=None, help="optional CSV report")
    ap.add_argument("--all", action="store_true", help="do not skip plots/logs/json")
    args = ap.parse_args()

    rows = compare_roots(args.ref_root, args.new_root, args.subdir, args.all)
    counts = {}
    for r in rows:
        counts[r[1]] = counts.get(r[1], 0) + 1
    for r in rows:
        if r[1] != "IDENTICAL":
            print(f"{r[1]:9s} {r[0]}  max_abs={r[2]} max_rel={r[3]} {r[4]}")
    print("summary:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["file", "status", "max_abs", "max_rel", "note"])
            w.writerows(rows)
    return 0


def compare_roots(ref_root, new_root, subdir="test_output", include_all=False):
    """Rows (relative file, status, max_abs, max_rel, note), sorted by file."""
    ref = os.path.join(ref_root, subdir)
    new = os.path.join(new_root, subdir)
    ROOT_SUBST.clear()
    for root in (ref_root, new_root):
        ROOT_SUBST.append((os.path.abspath(root), "<RUN_ROOT>"))
        ROOT_SUBST.append((os.path.realpath(root), "<RUN_ROOT>"))
    # Longest first: one root may be a prefix of the other (double, double_check).
    ROOT_SUBST.sort(key=lambda item: len(item[0]), reverse=True)
    rows = []
    for dirpath, _, files in os.walk(ref):
        for name in sorted(files):
            if not include_all and name.lower().endswith(SKIP_SUFFIXES):
                continue
            a = os.path.join(dirpath, name)
            rel = os.path.relpath(a, ref)
            b = os.path.join(new, rel)
            if not os.path.exists(b):
                rows.append((rel, "MISSING", "", "", ""))
                continue
            with open(a, "rb") as fa, open(b, "rb") as fb:
                da, db = fa.read(), fb.read()
            if da == db:
                rows.append((rel, "IDENTICAL", 0.0, 0.0, ""))
                continue
            if any(t in name for t in NONDETERMINISTIC_SUBSTRINGS):
                rows.append((rel, "NONDETERMINISTIC", "", "", "differs between reference runs too"))
                continue
            if b"\0" in da[:4096] or b"\0" in db[:4096]:
                rows.append((rel, "DIFFERENT", "", "", "binary"))
                continue
            status, mabs, mrel, note = compare_text(a, b)
            rows.append((rel, status, mabs, mrel, note))
    rows.sort()
    return rows


if __name__ == "__main__":
    sys.exit(main())
