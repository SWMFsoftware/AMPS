#!/usr/bin/env python3
"""Aggregate GCC JSON coverage for SWCME production sources only.

Header-only production code appears in several translation units.  Summing
those reports would inflate both covered and total counts, so this checker
merges executable lines by (source,line) and branch outcomes by
(source,line,branch-index), retaining the largest counter seen in any unit.
Generated fixtures and validation sources are intentionally outside the
production-source selector.
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def production_source(path: Path, root: Path) -> bool:
    """Return true only for top-level SWCME implementation/header files."""
    try:
        relative=path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return len(relative.parts)==1 and (
        relative.name=="swcme3d.cpp" or
        (relative.name.startswith("swcme") and relative.suffix==".hpp"))


def run_gcov(gcno_files: Iterable[Path], work: Path) -> List[Path]:
    reports: List[Path]=[]
    for index,gcno in enumerate(gcno_files):
        # GCC can emit a zero-length notes placeholder for a translation unit
        # whose callable validation code is removed from the coverage image by
        # the linker. Such a file contains no counters and `gcov` rejects it as
        # malformed. It cannot contribute production coverage, so ignore the
        # empty placeholder while still treating every nonempty gcov error as
        # fatal; this is not an exclusion of executable physics code.
        if gcno.stat().st_size==0:
            print(f"COV01 note: skipped empty gcno placeholder {gcno.name}")
            continue
        unit=work/f"unit-{index:04d}"
        unit.mkdir()
        proc=subprocess.run(
            ["gcov","--json-format","--branch-probabilities",
             "--branch-counts",str(gcno.resolve())],
            cwd=unit,text=True,stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,check=False)
        if proc.returncode!=0:
            raise RuntimeError(f"gcov failed for {gcno}:\n{proc.stdout}")
        generated=list(unit.glob("*.gcov.json.gz"))
        if len(generated)!=1:
            raise RuntimeError(
                f"gcov produced {len(generated)} JSON reports for {gcno}")
        reports.extend(generated)
    return reports


def aggregate(reports: Iterable[Path], root: Path) -> Dict[str,object]:
    line_counts: Dict[Tuple[str,int],int]={}
    branch_counts: Dict[Tuple[str,int,int],int]={}
    files=set()
    for report in reports:
        with gzip.open(report,"rt") as stream:
            data=json.load(stream)
        for entry in data.get("files",[]):
            source=Path(entry["file"])
            if not source.is_absolute():
                # GCC records paths relative to the directory in which the
                # instrumented compiler was invoked (`test/`), whereas this
                # script receives the repository root.  Resolve against that
                # build directory first; retain a root-relative fallback for
                # toolchains that emit paths in repository coordinates.
                build_relative=(root/"test"/source).resolve()
                root_relative=(root/source).resolve()
                source=build_relative if build_relative.exists() else root_relative
            if not production_source(source,root):
                continue
            relative=source.resolve().relative_to(root.resolve()).as_posix()
            files.add(relative)
            for line in entry.get("lines",[]):
                number=int(line["line_number"])
                key=(relative,number)
                line_counts[key]=max(line_counts.get(key,0),int(line.get("count",0)))
                # GCC's JSON order is stable within one source line.  Using the
                # branch index plus max-count merges duplicate header reports
                # without pretending that repeated translation units are new
                # branch outcomes.
                for branch_index,branch in enumerate(line.get("branches",[])):
                    bkey=(relative,number,branch_index)
                    branch_counts[bkey]=max(
                        branch_counts.get(bkey,0),int(branch.get("count",0)))

    covered_lines=sum(count>0 for count in line_counts.values())
    covered_branches=sum(count>0 for count in branch_counts.values())
    total_lines=len(line_counts)
    total_branches=len(branch_counts)
    return {
        "schema_version":1,
        "production_files":sorted(files),
        "line":{"covered":covered_lines,"total":total_lines,
                "percent":100.0*covered_lines/total_lines if total_lines else 0.0},
        "branch":{"covered":covered_branches,"total":total_branches,
                  "percent":100.0*covered_branches/total_branches
                  if total_branches else 0.0},
    }


def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--gcno-dir",type=Path,required=True)
    parser.add_argument("--root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--min-line",type=float,default=72.0)
    parser.add_argument("--min-branch",type=float,default=40.0)
    args=parser.parse_args()

    gcno=sorted(args.gcno_dir.glob("*.gcno"))
    if not gcno:
        raise SystemExit("COV01: no .gcno files found")
    with tempfile.TemporaryDirectory(prefix="swcme-cov01-gcov-") as tmp:
        result=aggregate(run_gcov(gcno,Path(tmp)),args.root)
    result["thresholds"]={"line_percent":args.min_line,
                          "branch_percent":args.min_branch}
    line=float(result["line"]["percent"])
    branch=float(result["branch"]["percent"])
    result["status"]="PASS" if (
        line>=args.min_line and branch>=args.min_branch) else "FAIL"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(f"COV01 production line coverage: {line:.2f}% "
          f"({result['line']['covered']}/{result['line']['total']})")
    print(f"COV01 production branch coverage: {branch:.2f}% "
          f"({result['branch']['covered']}/{result['branch']['total']})")
    print(f"COV01 report: {args.output}")
    return 0 if result["status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())
