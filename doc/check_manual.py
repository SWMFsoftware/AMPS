#!/usr/bin/env python3
"""Check the actual AMPS.tex include graph and its local assets/references.

Only files reached from the main document are build dependencies. Historical
research notes are deliberately not treated as current manual chapters. This
checker requires no AMPS checkout outside doc/ and performs no network access.
It complements, rather than replaces, a real LaTeX build and visual review.
"""

import collections
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent


def inspect():
    seen, visiting, errors = {}, set(), []

    def visit(relative):
        p = (ROOT / relative).resolve()
        try:
            p.relative_to(ROOT)
        except ValueError:
            errors.append("input outside manual directory: " + str(relative))
            return
        if p in visiting:
            errors.append("cyclic input: " + str(relative))
            return
        if p in seen:
            return
        if not p.is_file():
            errors.append("missing input: " + str(relative))
            return
        text = p.read_text(encoding="utf-8")
        # Listing text is an example, not TeX to recursively execute/resolve.
        text = re.sub(r"\\begin\{lstlisting\}.*?\\end\{lstlisting\}", "",
                      text, flags=re.S)
        text = re.sub(r"(?m)(?<!\\)%.*$", "", text)
        seen[p] = text
        visiting.add(p)
        for name in re.findall(r"\\input\{([^}]+)\}", text):
            q = Path(name)
            if not q.suffix:
                q = q.with_suffix(".tex")
            visit(q)
        for name in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", text):
            q = (ROOT / name).resolve()
            try:
                q.relative_to(ROOT)
            except ValueError:
                errors.append("graphic outside manual directory: " + name)
                continue
            if not q.is_file():
                errors.append("missing graphic: " + name)
        visiting.remove(p)

    visit(Path("AMPS.tex"))
    combined = "\n".join(seen.values())
    labels = re.findall(r"\\label\{([^}]+)\}", combined)
    for name, count in collections.Counter(labels).items():
        if count != 1:
            errors.append("duplicate label: " + name)
    references = re.findall(r"\\(?:ref|pageref|eqref)\{([^}]+)\}", combined)
    errors.extend("unknown reference: " + x for x in references if x not in labels)
    citations = re.findall(r"\\cite(?:[a-z]*)?(?:\[[^]]*\])?\{([^}]+)\}", combined)
    bibliography = set(re.findall(r"\\bibitem(?:\[[^]]*\])?\{([^}]+)\}", combined))
    for group in citations:
        errors.extend("unknown citation: " + x.strip() for x in group.split(",")
                      if x.strip() not in bibliography)
    required = ["Makefile", "README.md", "build_manual.py", "check_manual.py",
                "AMPS_clean_sweep_white.jpg"]
    errors.extend("missing build asset: " + x for x in required if not (ROOT / x).is_file())
    return seen, errors


if __name__ == "__main__":
    files, failures = inspect()
    if failures:
        print("\n".join("ERROR: " + x for x in failures), file=sys.stderr)
        sys.exit(1)
    print("PASS manual structure: {} reachable TeX files; local assets, "
          "unique labels and resolved citations/references".format(len(files)))
