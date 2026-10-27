#!/usr/bin/env python3
"""inventory.py -- which srcGeospace source files and headers are used, and by what.

Usage:
    srcGeospace/scripts/inventory.py --stage <AMPS build stage of srcGeospace>
                                     [--out doc/INVENTORY.md] [--list-unused] [--remove]

A file with extension .cpp .cc .c .h .hh .hpp .dfn is classified as

  amps-compiled  compiled by the AMPS build: <stage>/build/main/<path>.o exists
  amps-included  a dependency of an amps-compiled file, as reported by the compiler
                 (the stage's own compile command with -MM, run in <stage>/build/main)
  own-build      compiled by another build in the tree (test/*/run_test.sh,
                 tools/*/Makefile) or #included from such a file
                 (quoted #include closure; system and external headers are ignored)
  referenced     not built, but named in a test, script, list, data or document file
                 (a non-source text file mentions its file name)
  unused         none of the above

--remove deletes the unused files (the inventory written with --out lists them).  The
nvc++ GPU build compiles gpu/EngineDevice.cpp, which the CPU AMPS build compiles too, so
it is covered by amps-compiled/amps-included.
"""

import argparse
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict

SOURCE_EXT = (".cpp", ".cc", ".c")
HEADER_EXT = (".h", ".hh", ".hpp", ".dfn")
INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"]+)"', re.M)
SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def rel(path):
    return os.path.relpath(os.path.realpath(path), SRC)


def universe():
    out = []
    for root, _, files in os.walk(SRC):
        for name in files:
            if name.endswith(SOURCE_EXT + HEADER_EXT):
                out.append(rel(os.path.join(root, name)))
    return sorted(out)


def amps_compiled(stage):
    main = os.path.join(stage, "build", "main")
    out = set()
    for root, _, files in os.walk(main):
        for name in files:
            if not name.endswith(".o"):
                continue
            stem = os.path.relpath(os.path.join(root, name[:-2]), main)
            for ext in SOURCE_EXT:
                if os.path.isfile(os.path.join(SRC, stem + ext)):
                    out.add(stem + ext)
    return out


def compile_template(stage):
    """The stage's compile command for a main/ source, with -c/-o removed."""
    with open(os.path.join(stage, "build.log"), errors="replace") as f:
        for line in f:
            if " -c " in line and "gridless/" in line and line.split()[0].endswith(("mpicxx", "g++", "c++")):
                toks = shlex.split(line)
                outtoks, skip = [], 0
                for t in toks:
                    if skip:
                        skip -= 1
                        continue
                    if t == "-o":
                        skip = 1
                        continue
                    if t == "-c" or t.endswith(SOURCE_EXT):
                        continue
                    outtoks.append(t)
                return outtoks
    raise RuntimeError("no main/ compile command found in build.log")


def amps_included(stage, compiled):
    main = os.path.join(stage, "build", "main")
    cmd = compile_template(stage)
    deps = set()
    for src in sorted(compiled):
        r = subprocess.run(cmd + ["-MM", "-MG", src], cwd=main, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"-MM failed for {src}: {r.stderr[:400]}")
        text = r.stdout.replace("\\\n", " ")
        for tok in text.split(":", 1)[1].split():
            path = os.path.normpath(os.path.join(main, tok)) if not os.path.isabs(tok) else tok
            path = os.path.realpath(path)
            if path.startswith(os.path.realpath(main) + os.sep):
                r_ = os.path.relpath(path, os.path.realpath(main))
                if os.path.isfile(os.path.join(SRC, r_)):
                    deps.add(r_)
    return deps - compiled


def resolve(includer, name, search):
    for base in [os.path.dirname(includer)] + search:
        p = os.path.normpath(os.path.join(SRC, base, name))
        if os.path.isfile(p):
            return rel(p)
    return None


def include_closure(roots, search):
    seen, todo = set(), list(roots)
    while todo:
        f = todo.pop()
        if f in seen:
            continue
        seen.add(f)
        try:
            text = open(os.path.join(SRC, f), errors="replace").read()
        except OSError:
            continue
        for name in INCLUDE.findall(text):
            r = resolve(f, name, search)
            if r and r not in seen:
                todo.append(r)
    return seen


def own_build_roots(files):
    """Sources named by a build file (CMakeLists.txt, run_test.sh, Makefile)."""
    roots = defaultdict(set)
    builds = []
    for root, _, names in os.walk(SRC):
        for n in names:
            if n in ("CMakeLists.txt", "Makefile") or (n.endswith(".sh") and root.startswith(os.path.join(SRC, "test"))):
                builds.append(os.path.join(root, n))
    by_name = defaultdict(list)
    for f in files:
        by_name[os.path.basename(f)].append(f)
    for b in builds:
        text = open(b, errors="replace").read()
        bdir = rel(os.path.dirname(b))
        for tok in re.findall(r"[\w./$\-{}\"]+\.(?:cpp|cc|c)\b", text):
            name = os.path.basename(tok.strip('"'))
            cands = by_name.get(name, [])
            # Prefer a candidate below the build file's directory, then any.
            local = [c for c in cands if c.startswith(bdir + "/") or bdir == "."]
            for c in (local or cands):
                roots[rel(b)].add(c)
    return roots


def name_references(files, exclude=()):
    """File name -> set of non-source text files that mention it (the inventory's own
    output file is excluded, so a listed file does not count as referenced)."""
    exclude = {os.path.realpath(e) for e in exclude}
    names = {os.path.basename(f) for f in files}
    refs = defaultdict(set)
    pattern = re.compile(r"[\w.\-]+\.(?:cpp|cc|c|h|hh|hpp|dfn)\b")
    for root, _, ns in os.walk(SRC):
        for n in ns:
            if n.endswith(SOURCE_EXT + HEADER_EXT):
                continue
            p = os.path.join(root, n)
            if os.path.realpath(p) in exclude:
                continue
            try:
                data = open(p, "rb").read()
            except OSError:
                continue
            if b"\0" in data[:8192]:
                continue
            for m in set(pattern.findall(data.decode("utf-8", "replace"))):
                if m in names:
                    refs[m].add(rel(p))
    return refs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--list-unused", action="store_true")
    ap.add_argument("--remove", action="store_true")
    args = ap.parse_args()

    files = universe()
    compiled = amps_compiled(args.stage)
    included = amps_included(args.stage, compiled)
    roots = own_build_roots(files)
    own = {}
    search = {"test": [".", "gridless", "util"], "tools": ["."]}
    for build, rs in roots.items():
        top = build.split("/")[0]
        dirs = search.get(top, ["."])
        if top == "test":
            # Unit-test scripts put their own directory (and stubs/) first on the
            # include path (-I"$script_dir"), so stub headers there shadow AMPS headers.
            bdir = os.path.dirname(build)
            dirs = [bdir, os.path.join(bdir, "stubs")] + dirs + ["3d"]
        for f in include_closure(rs, dirs):
            own.setdefault(f, build)
    refs = name_references(files, [args.out] if args.out else [])

    rows = []
    for f in files:
        if f in compiled:
            status, why = "amps-compiled", "AMPS build"
        elif f in included:
            status, why = "amps-included", "compiler dependency of the AMPS build"
        elif f in own:
            status, why = "own-build", own[f]
        elif os.path.basename(f) in refs:
            r = sorted(refs[os.path.basename(f)])
            status, why = "referenced", ", ".join(r[:3]) + (" ..." if len(r) > 3 else "")
        else:
            status, why = "unused", ""
        rows.append((f, status, why))

    counts = defaultdict(int)
    for _, s, _ in rows:
        counts[s] += 1
    unused = [f for f, s, _ in rows if s == "unused"]
    print("inventory:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)
    if args.list_unused:
        print("\n".join(unused))
    if args.out:
        with open(args.out, "w") as f:
            f.write(render(rows, counts, args))
            if args.remove and unused:
                f.write("## Removed by `inventory.py --remove`\n\n"
                        "These files were not compiled, included or referenced by any build, "
                        "test, script or document (they remain in srcEarth):\n\n")
                f.write("".join(f"* `{u}`\n" for u in unused))
    if args.remove:
        for u in unused:
            os.remove(os.path.join(SRC, u))
        print(f"removed {len(unused)} unused files", file=sys.stderr)
    return 0


def render(rows, counts, args):
    by_dir = defaultdict(list)
    for f, s, w in rows:
        by_dir[os.path.dirname(f) or "."].append((f, s, w))
    out = ["# srcGeospace source inventory", "",
           "Generated by `scripts/inventory.py --stage <AMPS build stage>`; see the script header",
           "for the classification rules.  Status counts: " +
           ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) + ".", ""]
    for d in sorted(by_dir):
        out += [f"## {d}", "", "| file | status | used by |", "|---|---|---|"]
        for f, s, w in sorted(by_dir[d]):
            out.append(f"| `{os.path.basename(f)}` | {s} | {w} |")
        out.append("")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.exit(main())
