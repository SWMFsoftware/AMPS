#!/usr/bin/env python3
"""Join test_runner reports of several run roots into one comparison table.

Usage:
    summarize.py --root cpu=<run-root-1> --root gpu=<run-root-2> [--root ...] [--out summary.md]

Each root must contain test_runner_report.json written by srcEarth/test/test_runner.py
(or by the srcEarth runner for older roots).
Tests are matched by their command line.  The table lists, per test, the actual status
(P/F, '!' marks a mismatch with the list expectation) and the elapsed time of every
root, plus the time ratio of every root relative to the first one.
"""

import argparse
import json
import os
import sys


def load(root):
    path = os.path.join(root, "test_runner_report.json")
    with open(path) as f:
        rows = json.load(f)
    return {r["command"]: r for r in rows}


def short(cmd, width=70):
    cmd = cmd.replace("python3 ", "").replace("./", "")
    cmd = cmd.replace("srcEarth/test/", "").replace("srcGeospace/test/", "")
    return cmd if len(cmd) <= width else cmd[: width - 3] + "..."


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", action="append", required=True, help="label=path")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    roots = []
    for item in args.root:
        label, path = item.split("=", 1)
        roots.append((label, load(path)))

    commands = []
    for _, rep in roots:
        for c in rep:
            if c not in commands:
                commands.append(c)

    head = ["#", "test"]
    for label, _ in roots:
        head += [f"{label} status", f"{label} s"]
    for label, _ in roots[1:]:
        head.append(f"{label}/{roots[0][0]}")
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    totals = [0.0 for _ in roots]
    for i, cmd in enumerate(commands, 1):
        row = [str(i), "`" + short(cmd) + "`"]
        times = []
        for k, (_, rep) in enumerate(roots):
            r = rep.get(cmd)
            if r is None:
                row += ["-", "-"]
                times.append(None)
                continue
            mark = "" if r.get("matched_reference", True) else "!"
            row += [f"{r['actual']}{mark} (exp {r['expected']})", f"{r['elapsed_s']:.1f}"]
            times.append(r["elapsed_s"])
            totals[k] += r["elapsed_s"]
        for t in times[1:]:
            row.append(f"{t / times[0]:.2f}" if t is not None and times[0] else "-")
        lines.append("| " + " | ".join(row) + " |")
    tot = ["", "**total**"]
    for t in totals:
        tot += ["", f"{t:.0f}"]
    for t in totals[1:]:
        tot.append(f"{t / totals[0]:.2f}" if totals[0] else "-")
    lines.append("| " + " | ".join(tot) + " |")
    text = "\n".join(lines) + "\n"
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
