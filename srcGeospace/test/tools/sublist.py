#!/usr/bin/env python3
"""Write a subset of a test_runner list file.

Usage:
    sublist.py LIST [--only REGEX]... [--exclude REGEX]... [--np N] [--nt M]
               [--np-nt-exclude REGEX] [-o OUT]

A command entry (a line starting with P, F or a {P,F,...} loop marker) is kept when it
matches at least one --only pattern (or no --only is given) and no --exclude pattern.
The entry's attached lines go with it: the '! runner: ...' directives and the
'for $x=...' loop header that precede it, and the 'last pass:' line that follows it.
Comments and blank lines are always kept, so line structure stays readable.  The input
list is never modified; the result goes to OUT (default: standard output), which should
lie outside the source tree, e.g. in the run root.

--np/--nt replace the values of the -np/--np and -nt/--nt tokens that a kept command
already has (commands without them are unchanged), except in commands matching
--np-nt-exclude (default: the C13 and C18 runners, which choose their own rank
counts).  This is how a GPU pass runs srcEarth's tests with, e.g., one rank and 32
threads per GPU (test/tools/run_srcearth_tests.sh --np/--nt).

Examples:
    # the GPU validation scope: every test except SA01-SA03, one rank with 32 threads
    sublist.py srcEarth/test/list --exclude 'test/SA0[123]/' --np 1 --nt 32 -o <run-root>/list_gpu
"""

import argparse
import re
import sys

COMMAND = re.compile(r"^(?:[PF]\s|\{[PF](?:,[PF])*\}\s)")
ATTACHED_BEFORE = re.compile(r"^(?:!\s*runner:|for\s+\$)")
LAST_PASS = re.compile(r"^!?last pass:")


def rewrite_np_nt(command, np_value, nt_value, exclude):
    """Replace the values of existing -np/--np and -nt/--nt tokens; a command without
    them, or matching `exclude`, is returned unchanged."""
    if (np_value is None and nt_value is None) or (exclude is not None and exclude.search(command)):
        return command
    out = command
    for names, value in ((("-np", "--np"), np_value), (("-nt", "--nt"), nt_value)):
        if value is None:
            continue
        for name in names:
            out = re.sub(r"(?<!\S)(" + re.escape(name) + r")(\s+|=)(\d+)(?!\S)",
                         lambda m: f"{m.group(1)}{m.group(2)}{value}", out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("list")
    ap.add_argument("--only", action="append", default=[], help="keep matching commands")
    ap.add_argument("--exclude", action="append", default=[], help="drop matching commands")
    ap.add_argument("--np", type=int, default=None, help="new value of existing -np/--np tokens")
    ap.add_argument("--nt", type=int, default=None, help="new value of existing -nt/--nt tokens")
    ap.add_argument("--np-nt-exclude", default=r"test/C1[38]/",
                    help="commands matching this are never rewritten (default: %(default)s)")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()
    only = [re.compile(p) for p in args.only]
    exclude = [re.compile(p) for p in args.exclude]
    np_nt_exclude = re.compile(args.np_nt_exclude) if args.np_nt_exclude else None
    rewritten = 0

    out = []
    pending = []          # attached lines (and interleaved comments) before a command
    dropping_tail = False
    kept = dropped = 0
    with open(args.list) as f:
        for line in f:
            text = line.rstrip("\n")
            if dropping_tail and LAST_PASS.match(text):
                continue
            dropping_tail = False
            if ATTACHED_BEFORE.match(text):
                pending.append((True, line))
                continue
            if COMMAND.match(text):
                keep = (not only or any(p.search(text) for p in only)) and \
                       not any(p.search(text) for p in exclude)
                for attached, buffered in pending:
                    if keep or not attached:
                        out.append(buffered)
                pending = []
                if keep:
                    new = rewrite_np_nt(text, args.np, args.nt, np_nt_exclude)
                    if new != text:
                        rewritten += 1
                    out.append(new + ("\n" if line.endswith("\n") else ""))
                    kept += 1
                else:
                    dropping_tail = True
                    dropped += 1
                continue
            if pending:
                pending.append((False, line))
            else:
                out.append(line)
    out.extend(buffered for _, buffered in pending)
    data = "".join(out)
    if args.out:
        with open(args.out, "w") as f:
            f.write(data)
    else:
        sys.stdout.write(data)
    print(f"sublist: kept {kept} command entries, dropped {dropped}"
          + (f", -np/-nt rewritten in {rewritten}" if args.np is not None or args.nt is not None else ""),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
