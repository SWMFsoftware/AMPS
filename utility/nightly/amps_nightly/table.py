"""Reader for MakefileTest/Table: test names, their Family= and the result files each block mentions.

The block grammar follows utility/TestScripts/table_parser.cpp: a block runs from a line "#>"
to a line "<#"; keyword lines (Name=, Family=, ...) come before the first section header
(Compile=, Rundir=, Run=, Check=), after which every line belongs to that section.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set

from .registry import Registry

DEFAULT_TABLE = Path(__file__).resolve().parents[3] / "MakefileTest" / "Table"
_SECTIONS = ("Compile=", "Rundir=", "Run=", "Check=")
_DIFF_FILE = re.compile(r"\btest_([A-Za-z0-9][A-Za-z0-9_.+=-]*?)(?:\$\([^)]*\))?\.diff\b")


@dataclass
class TableTest:
    name: str                 # as used for make targets: test_<name>_compile, ...
    family: str               # "" when the block has no Family= line
    selector: str             # the @@(...) machine selector of Name=, "" for all machines
    line: int                 # line of "#>" (1-based)
    family_lines: int = 0     # number of Family= lines in the block
    result_files: Set[str] = field(default_factory=set)  # other test_<x>.diff names written by the block


def read_table(path: Optional[Path] = None) -> List[TableTest]:
    path = Path(path) if path else DEFAULT_TABLE
    tests: List[TableTest] = []
    current: Optional[TableTest] = None
    section = False
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        line = raw.rstrip()
        if line == "#>":
            current, section = TableTest(name="", family="", selector="", line=number), False
            continue
        if current is None:
            continue
        if line == "<#":
            if current.name:
                current.result_files.discard(current.name)
                tests.append(current)
            current = None
            continue
        if line in _SECTIONS:
            section = True
        if section:
            current.result_files.update(_DIFF_FILE.findall(line))
            continue
        if line.startswith("Name="):
            value = line[len("Name="):]
            name, _, selector = value.partition("@@")
            current.name = name.strip().replace("test/", "")
            current.selector = selector.strip()
        elif line.startswith("Family="):
            current.family = line[len("Family="):].strip()
            current.family_lines += 1
    return tests


def family_by_name(tests: List[TableTest]) -> Dict[str, str]:
    return {t.name: t.family for t in tests if t.family}


def check_table(tests: List[TableTest], registry: Registry) -> List[str]:
    """Problems that would put a test in the wrong (or no) family on the report."""
    problems: List[str] = []
    families: Dict[str, str] = {}
    for test in tests:
        where = "Table line {} ({})".format(test.line, test.name)
        if test.family_lines == 0:
            problems.append(where + ": no Family= line")
            continue
        if test.family_lines > 1:
            problems.append(where + ": more than one Family= line")
        family = registry.family(test.family)
        if family is None:
            problems.append(where + ": Family={} is not defined in families.json".format(test.family))
        elif family["group"] != "legacy":
            problems.append(where + ": Family={} belongs to group {}, not legacy".format(test.family, family["group"]))
        previous = families.setdefault(test.name, test.family)
        if previous != test.family:
            problems.append(where + ": Family={} differs from another block of the same test ({})".format(test.family, previous))
    return problems
