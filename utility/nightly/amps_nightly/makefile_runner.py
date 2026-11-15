"""The built-in runner for the Makefile.test tests (MakefileTest/Table).

Per test, the four legacy stages run as separate make targets, each with the machine's time limit for
that stage, and their output goes to one log per test (logs/<runner>/tests/<test>.log):
  phase A, one test after another in the workspace (compiling reconfigures the tree):
      make test_<X>_compile, make test_<X>_rundir
  phase B, up to `lanes` tests at once, each in its own run_test_<X> directory:
      make test_<X>_run, make test_<X>_check
This is the split of the legacy scheduler, with time limits and per-test logs added. Status, stage and
reason come from the .diff result files exactly as for imported legacy nights; .diff contents are not kept.
"""

from __future__ import annotations

import concurrent.futures
import shutil
import subprocess
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import legacy, procs
from .machines import Machine
from .plan import RUN, Item, makefile_manifest
from .registry import UNASSIGNED
from .table import DEFAULT_TABLE, TableTest, family_by_name, read_table

TOOL_PARSER = DEFAULT_TABLE.parents[1] / "utility" / "TestScripts" / "table_parser.cpp"


def manifest(ws: Path, machine: Machine, env: Dict[str, str], log) -> List[Tuple[str, str, str]]:
    """(name, family, exception code) of the tests selected in this workspace.

    Uses the workspace's own `make test_manifest`. A commit made before that target existed gets the
    list from this tool's table_parser run on the workspace's own Table (families from this tool's Table).
    """
    try:
        return makefile_manifest(ws)
    except (RuntimeError, OSError, subprocess.SubprocessError) as error:
        log("make test_manifest unavailable in the tested commit ({}); using this tool's table_parser".format(error))
    work = ws / ".nightly-manifest"
    work.mkdir(exist_ok=True)
    subprocess.run(["g++", "-std=c++14", str(TOOL_PARSER), "-o", "table_parser"], cwd=str(work), check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.copy(str(ws / "MakefileTest" / "Table"), str(work / "Table"))
    subprocess.run(["./table_parser", "2"], cwd=str(work), check=True, stdout=subprocess.DEVNULL,
                   env=dict(env, AMPS_TABLE_HOST=machine.table_host))
    done = subprocess.run(["make", "-s", "-f", "Makefile.test", "-f", str(ws / "Makefile.conf"), "test_manifest"],
                          cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                          timeout=300)
    if done.returncode != 0:
        raise RuntimeError("test_manifest from this tool's table_parser failed: " + done.stderr.strip()[-300:])
    families = family_by_name(read_table())
    tests = []
    for line in done.stdout.splitlines():
        parts = line.split("|")
        if len(parts) == 5 and parts[0]:
            tests.append((parts[0], parts[1] or families.get(parts[0], UNASSIGNED), parts[3]))
    return tests


def _rundir_done(ws: Path, name: str) -> bool:
    """The test's own result file says its run directory is ready (exit codes of legacy targets are not reliable)."""
    for path in ws.glob("test_{}*.diff".format(name)):
        if legacy.split_result_file(path.name)[0] != name:
            continue
        try:
            if "rundir... done" in path.read_text(errors="replace")[-2000:]:
                return True
        except OSError:
            pass
    return False


def run(ws: Path, machine: Machine, items: Sequence[Item], env: Dict[str, str], logs: Path, tag: str,
        deadline: float, log: Callable[[str], None]) -> List[Dict]:
    """Run the planned Makefile.test tests in ws; return result rows (with log_src = the per-test log)."""
    make_vars = ["{}={}".format(k, v) for k, v in sorted((machine.make_vars or {}).items())]
    tests_dir = logs / "tests"
    tests_dir.mkdir(parents=True, exist_ok=True)
    limits = machine.timeouts_s
    timing: Dict[str, float] = {}
    over: Dict[str, str] = {}          # test -> "stage exceeded its time limit of N s"
    not_started: List[str] = []

    def stage(name: str, which: str) -> bool:
        target = "test_{}_{}".format(name, which)
        limit = min(limits[which], max(1.0, deadline - time.time()))
        outcome = procs.run(["make", target] + make_vars, ws, env, tests_dir / (name + ".log"), limit,
                            "{}/{}".format(tag, name), header=target)
        timing[name] = timing.get(name, 0.0) + outcome.elapsed_s
        if outcome.timed_out:
            over[name] = "{} stage exceeded {} ({:.0f} s)".format(
                which, "the night's time limit" if limit < limits[which] else "its time limit", limit)
        return outcome.exit_code == 0 and not outcome.timed_out

    # the shared library only: the legacy target test_compile also compiles every test
    procs.run(procs.shell("mkdir -p share/lib && cd share/Library/src && make LIB"), ws, env, logs / "setup.log",
              limits["compile"], tag + "/setup", header="shared library (share/Library/src: make LIB)")
    planned = [i.test for i in items if i.action == RUN]
    log("Makefile.test: {} tests to run, {} lanes".format(len(planned), machine.lanes))
    ready = []
    for name in planned:                                   # phase A: one at a time
        if time.time() >= deadline:
            not_started.append(name)
            continue
        stage(name, "compile")
        if name not in over:
            stage(name, "rundir")
        if _rundir_done(ws, name):
            ready.append(name)
    with concurrent.futures.ThreadPoolExecutor(max_workers=machine.lanes) as pool:   # phase B: lanes
        def run_and_check(name: str) -> None:
            if time.time() >= deadline:
                not_started.append(name)
                return
            stage(name, "run")
            stage(name, "check")
        list(pool.map(run_and_check, ready))

    results = {r["test"]: r for r in legacy.collect_results(legacy.read_diff_files(ws), read_table())}
    rows = []
    for item in items:
        if item.action != RUN:
            continue
        name = item.test
        result = results.pop(name, None)
        if name in not_started:
            row = {"status": "SKIPPED", "stage": "", "reason": "not started: the night's time limit was reached"}
        elif result is None:
            row = {"status": "FAILED", "stage": "compile", "reason": "no result file test_{}*.diff was written".format(name)}
        else:
            row = {"status": result["status"], "stage": result.get("stage", ""), "reason": result.get("reason", "")}
            if name in over:
                row.update(status="FAILED", stage=over[name].split()[0], reason=over[name])
        row.update(test=name, family=item.family, duration_s=round(timing.get(name, 0.0), 1),
                   log_src=str(tests_dir / (name + ".log")))
        rows.append(row)
    if results:
        log("Makefile.test: ignored result files of tests that were not run: " + ", ".join(sorted(results)))
    return rows
