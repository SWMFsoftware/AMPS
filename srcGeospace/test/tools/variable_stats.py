#!/usr/bin/env python3
"""Per-variable difference statistics of Tecplot outputs between two run roots.

Usage:
    variable_stats.py REF_ROOT NEW_ROOT [--subdir test_output] [--pattern REGEX] [--out stats.csv]

compare_outputs.py reports one maximum per file, which for the float build is dominated
by single penumbral locations and by diagnostic columns.  This script reads every
Tecplot ASCII file (*.dat with a VARIABLES line and POINT zones) present in both roots
and reports, for every (test directory, file name, variable), the number of values, the
fraction that are identical as printed, and the median, 90th/99th percentile and maximum
of the relative difference |a-b|/max(|a|,|b|) over the values that differ.  Values are
compared as printed, so differences below the printed precision are not visible.
Neither root is modified.
"""

import argparse
import csv
import os
import re
import sys

VARS = re.compile(r'VARIABLES\s*=\s*(.*)', re.IGNORECASE)


def read_tecplot(path):
    names = None
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            m = VARS.match(s)
            if m:
                names = [v.strip().strip('"') for v in re.split(r'[,\s]+', m.group(1)) if v.strip()]
                continue
            if s[0].isalpha():
                continue
            try:
                vals = [float(t) for t in s.replace(",", " ").split()]
            except ValueError:
                continue
            rows.append(vals)
    return names, rows


def percentile(sorted_vals, q):
    if not sorted_vals:
        return 0.0
    k = min(len(sorted_vals) - 1, max(0, int(round(q * (len(sorted_vals) - 1)))))
    return sorted_vals[k]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref_root")
    ap.add_argument("new_root")
    ap.add_argument("--subdir", default="test_output")
    ap.add_argument("--pattern", default=None, help="only files whose relative path matches")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rows = collect(args.ref_root, args.new_root, args.subdir, args.pattern)
    print("| test | file | variable | values | identical | median rel | p90 rel | p99 rel | max rel |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print("| %s | %s | %s | %d | %.3f | %.2e | %.2e | %.2e | %.2e |" % r)
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["test", "file", "variable", "values", "identical_fraction",
                        "median_rel_of_differing", "p90_rel", "p99_rel", "max_rel"])
            w.writerows(rows)
    return 0


def collect(ref_root, new_root, subdir="test_output", pattern=None):
    """Rows (test, file, variable, values, identical fraction, median, p90, p99, max)
    for every variable that differs in at least one value."""
    ref = os.path.join(ref_root, subdir)
    new = os.path.join(new_root, subdir)
    pat = re.compile(pattern) if pattern else None

    # (test, file name, variable) -> [n, n_identical, list of rel diffs]
    acc = {}
    for dirpath, _, files in os.walk(ref):
        for name in files:
            if not name.endswith(".dat"):
                continue
            a = os.path.join(dirpath, name)
            rel = os.path.relpath(a, ref)
            if pat and not pat.search(rel):
                continue
            b = os.path.join(new, rel)
            if not os.path.exists(b):
                continue
            na, ra = read_tecplot(a)
            nb, rb = read_tecplot(b)
            if not na or na != nb or len(ra) != len(rb):
                continue
            test = rel.split(os.sep)[0]
            for va, vb in zip(ra, rb):
                if len(va) != len(na) or len(vb) != len(na):
                    continue
                for j, var in enumerate(na):
                    e = acc.setdefault((test, name, var), [0, 0, []])
                    e[0] += 1
                    x, y = va[j], vb[j]
                    if x == y:
                        e[1] += 1
                        continue
                    scale = max(abs(x), abs(y))
                    e[2].append(abs(x - y) / scale if scale > 0 else 0.0)

    rows = []
    for (test, name, var), (n, nid, d) in sorted(acc.items()):
        if nid == n:
            continue   # coordinates and values that never differ
        d.sort()
        rows.append((test, name, var, n, nid / n, percentile(d, 0.5), percentile(d, 0.9),
                     percentile(d, 0.99), d[-1]))
    return rows


if __name__ == "__main__":
    sys.exit(main())
