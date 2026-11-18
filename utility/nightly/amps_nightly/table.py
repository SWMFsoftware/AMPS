"""Reader for MakefileTest/Table: test names, their Family= and the result files each block mentions.

The block grammar follows utility/TestScripts/table_parser.cpp: a block runs from a line "#>"
to a line "<#"; keyword lines (Name=, Family=, ...) come before the first section header
(Compile=, Rundir=, Run=, Check=), after which every line belongs to that section.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .registry import Registry

DEFAULT_TABLE = Path(__file__).resolve().parents[3] / "MakefileTest" / "Table"
_SECTIONS = ("Compile=", "Rundir=", "Run=", "Check=")
_DIFF_FILE = re.compile(r"\btest_([A-Za-z0-9][A-Za-z0-9_.+=-]*?)(?:\$\([^)]*\))?\.diff\b")
_RUN_DIR = re.compile(r"\$[({]TEST([A-Za-z0-9][A-Za-z0-9_.+=-]*?)DIR[)}]")    # ${TEST<name>DIR}: run_test_<name>
# the SWMF trees a block builds or runs in, in the workspace: SWMF/..., cd SWMF-IH/SWMF; (not SWMF_data, not SWMF.exe)
_SWMF_TREE = re.compile(r"(?<![\w./-])SWMF(?:-IH)?(?=/|;|\s|$)")
# a make call with a fixed number of jobs (make -j8, $(MAKE) -j16 amps); -j1 (serial on purpose) is allowed
_FIXED_JOBS = re.compile(r"(?:\bmake\b|\$[({]MAKE[)}]).*?\s-j\s*([0-9]+)\b")
COMPILE_JOBS = "-j$(TESTCOMPILEJOBS)"


@dataclass
class TableTest:
    name: str                 # as used for make targets: test_<name>_compile, ...
    family: str               # "" when the block has no Family= line
    selector: str             # the @@(...) machine selector of Name=, "" for all machines
    line: int                 # line of "#>" (1-based)
    family_lines: int = 0     # number of Family= lines in the block
    result_files: Set[str] = field(default_factory=set)  # other test_<x>.diff names written by the block
    uses: Set[str] = field(default_factory=set)          # other tests whose run directory (${TEST<x>DIR}) it uses
    type: str = "amps"        # Type= of the block: amps, or swmf (an SWMF-coupled test)
    swmf_tree: bool = False   # a section uses the workspace's SWMF or SWMF-IH tree
    fixed_jobs: List[Tuple[int, str]] = field(default_factory=list)   # (line, text): make -j<N>, N > 1

    @property
    def workspace_only(self) -> bool:
        """It builds in the workspace's SWMF trees (or is an SWMF test): it cannot compile in a copy of the workspace
        (makefile_runner compile lanes)."""
        return self.type == "swmf" or self.swmf_tree


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
                current.uses.discard(current.name)
                tests.append(current)
            current = None
            continue
        if line in _SECTIONS:
            section = True
        if section:
            current.result_files.update(_DIFF_FILE.findall(line))
            current.uses.update(_RUN_DIR.findall(line))
            current.swmf_tree = current.swmf_tree or bool(_SWMF_TREE.search(line))
            jobs = _FIXED_JOBS.search(line)
            if jobs and int(jobs.group(1)) > 1:
                current.fixed_jobs.append((number, line.lstrip(">").strip()))
            continue
        if line.startswith("Name="):
            value = line[len("Name="):]
            name, _, selector = value.partition("@@")
            current.name = name.strip().replace("test/", "")
            current.selector = selector.strip()
        elif line.startswith("Family="):
            current.family = line[len("Family="):].strip()
            current.family_lines += 1
        elif line.startswith("Type="):
            current.type = line[len("Type="):].strip()
    return tests


# Variants of a test: a Table block named <test>--<variant> runs <test> in another way (its own build, run and check).
# Make targets cannot hold parentheses (make reads x(y) as an archive member), so the Table name carries the variant
# as a suffix; results, logs and pages show it as <test>(<variant>). A variant can need a machine capability: it is
# SKIPPED on machines whose file does not declare it.
VARIANTS = {"valgrind": "valgrind"}        # variant -> machine capability it requires ("" for none)


def _variant(name: str) -> str:
    base, sep, variant = name.rpartition("--")
    return variant if sep and base and variant in VARIANTS else ""


def display_name(name: str) -> str:
    """'Mars-ions--valgrind' -> 'Mars-ions(valgrind)'; other names are returned unchanged."""
    variant = _variant(name)
    return "{}({})".format(name[:-len(variant) - 2], variant) if variant else name


def table_name(display: str) -> str:
    """The inverse of display_name: 'Mars-ions(valgrind)' -> 'Mars-ions--valgrind'."""
    match = re.match(r"^(.+)\(([^()]+)\)$", display)
    return "{}--{}".format(match.group(1), match.group(2)) if match and match.group(2) in VARIANTS else display


def required_capabilities(name: str) -> List[str]:
    """The machine capabilities a Table test needs (Table name or display name)."""
    capability = VARIANTS.get(_variant(table_name(name)), "")
    return [capability] if capability else []


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
    for test in tests:
        for number, text in test.fixed_jobs:
            problems.append("Table line {} ({}): '{}' fixes the number of make jobs: use {} (the nightly sets it from "
                            "compile_jobs of the machine file; 8 otherwise)".format(number, test.name, text, COMPILE_JOBS))
    return problems
