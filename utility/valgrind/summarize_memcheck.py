#!/usr/bin/env python3
"""Group valgrind memcheck errors from one or more log files.

Usage: summarize_memcheck.py [--all-frames] LOG...

Each error is reduced to its kind (e.g. "Conditional jump or move depends on
uninitialised value(s)"), the first stack frame that lies in AMPS source code
(a frame with a file:line outside /usr and the MPI/system libraries) and, when
--track-origins=yes was used, the first AMPS frame of the origin.  Errors with
the same key are counted once per log and listed with the logs that had them.
The exit status is 1 if any error was found and 0 otherwise.
"""
import re
import sys
from collections import OrderedDict

ERR_START = re.compile(r"^==\d+== (\S.*)$")
FRAME = re.compile(r"^==\d+==\s+(?:at|by) 0x[0-9A-Fa-f]+: (.*?)(?: \((.*)\))?$")
ORIGIN = re.compile(r"^==\d+==\s+Uninitialised value was created by (.*)$")
ADDR = re.compile(r"^==\d+==\s+Address 0x[0-9A-Fa-f]+ is (.*)$")

KINDS = (
    "Conditional jump or move depends on uninitialised value",
    "Use of uninitialised value of size",
    "Syscall param",
    "Invalid read of size",
    "Invalid write of size",
    "Invalid free",
    "Mismatched free",
    "Source and destination overlap",
    "Argument 'size' of function",
    "Uninitialised byte(s) found during client check request",
)

SYSTEM = ("/usr/", "vg_replace", "memcheck_preload.c", "libmpi", "libopen-", "libc.so", "libpthread", "ld-2.", "libgfortran")


def is_amps(frame_loc):
    if frame_loc is None:
        return False
    if any(s in frame_loc for s in SYSTEM):
        return False
    return ":" in frame_loc and not frame_loc.startswith("in ")


def parse(path):
    errors = []
    cur = None
    section = "stack"
    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            m = ERR_START.match(line)
            if m and any(m.group(1).startswith(k) for k in KINDS):
                cur = {"kind": m.group(1), "stack": [], "origin": [], "addr": ""}
                errors.append(cur)
                section = "stack"
                continue
            if cur is None:
                continue
            m = ORIGIN.match(line)
            if m:
                cur["origin_what"] = m.group(1)
                section = "origin"
                continue
            m = ADDR.match(line)
            if m:
                cur["addr"] = m.group(1)
                section = "addr"
                continue
            m = FRAME.match(line)
            if m:
                if section == "stack":
                    cur["stack"].append((m.group(1), m.group(2)))
                elif section == "origin":
                    cur["origin"].append((m.group(1), m.group(2)))
                continue
            if re.match(r"^==\d+==\s*$", line):
                # a blank line ends the stack or origin section of the error
                if section in ("origin", "addr") or (section == "stack" and cur["stack"]):
                    if section == "stack":
                        section = "after-stack"
                    else:
                        cur = None
    return errors


def first_amps(frames):
    for fn, loc in frames:
        if is_amps(loc):
            return "%s (%s)" % (fn, loc)
    if frames:
        fn, loc = frames[0]
        return "%s (%s)" % (fn, loc or "?")
    return "?"


def main(argv):
    all_frames = "--all-frames" in argv
    logs = [a for a in argv if a != "--all-frames"]
    groups = OrderedDict()
    for log in logs:
        for e in parse(log):
            kind = re.sub(r"\(s\)", "", e["kind"])
            kind = re.sub(r"0x[0-9A-Fa-f]+", "0x", kind)
            key = (kind, first_amps(e["stack"]), first_amps(e["origin"]) if e["origin"] else e.get("origin_what", e["addr"]))
            g = groups.setdefault(key, {"count": 0, "logs": OrderedDict(), "example": e})
            g["count"] += 1
            g["logs"][log] = g["logs"].get(log, 0) + 1
    for n, (key, g) in enumerate(groups.items(), 1):
        kind, where, origin = key
        print("[%d] %s" % (n, kind))
        print("    at:     %s" % where)
        print("    origin: %s" % origin)
        print("    count:  %d in %s" % (g["count"], ", ".join("%s(%d)" % (l, c) for l, c in g["logs"].items())))
        if all_frames:
            for fn, loc in g["example"]["stack"][:12]:
                print("      at %s (%s)" % (fn, loc))
            if g["example"]["origin"]:
                print("      origin:")
                for fn, loc in g["example"]["origin"][:8]:
                    print("      by %s (%s)" % (fn, loc))
    print("%d distinct error(s)" % len(groups))
    return 1 if groups else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
