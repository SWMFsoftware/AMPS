"""The built-in runner for the Makefile.test tests (MakefileTest/Table).

Per test, the four legacy stages run as separate make targets, each with the machine's time limit for
that stage, and their output goes to one log per test (logs/<runner>/tests/<test>.log):
  phase A, `compile_lanes` tests at once, each lane one test after another in its own tree (compiling
  reconfigures the tree; the files Config.pl rewrites are reset to their post-install state before each test,
  see CONFIG_FILES). Lane 1 is the workspace; the others are copies of it (compile_lane), and a test compiled
  in a copy has its run directory and result files moved into the workspace after its rundir stage. The SWMF
  tests compile in the workspace, where their SWMF trees are; a test that uses the run directory of another
  test compiles in the same tree, after it (compile_groups). Each compile runs make -j<compile_jobs>
  (TESTCOMPILEJOBS):
      make test_<X>_compile, make test_<X>_rundir
  phase B, up to `lanes` tests at once in the workspace, each in its own run_test_<X> directory:
      make test_<X>_run, make test_<X>_check
  (a test whose Table block uses the run directory of another test, ${TEST<Y>DIR}, runs after the run of Y)
This is the split of the legacy scheduler, with time limits and per-test logs added. Status, stage and
reason come from the .diff result files exactly as for imported legacy nights; .diff contents are not kept.
"""

from __future__ import annotations

import collections
import concurrent.futures
import os
import re
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import legacy, procs
from .machines import COMPILE_JOBS_VAR, Config, Machine
from .plan import RUN, Item, makefile_manifest
from .registry import UNASSIGNED
from .table import DEFAULT_TABLE, TableTest, display_name, family_by_name, read_table

TOOL_PARSER = DEFAULT_TABLE.parents[1] / "utility" / "TestScripts" / "table_parser.cpp"


def workspace_table(ws: Path) -> List[TableTest]:
    """The tests of the tested workspace's own MakefileTest/Table (this tool's Table if it has none).

    Result files are attributed to tests by name, so the names must come from the commit under test: a test
    that exists only there would otherwise be folded into the test whose name is its prefix
    (test_unit_test_cutcell_tetra.diff into unit_test) and reported as having written no result file.
    """
    path = ws / "MakefileTest" / "Table"
    return read_table(path if path.is_file() else None)


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
    # Makefile.conf first: the compiler-specific rules of Makefile.test (ifeq ($(COMPILE.c),gcc) ... SKIP or Ref=)
    # are evaluated while it is read, so COMPILE.c must already be defined
    done = subprocess.run(["make", "-s", "-f", str(ws / "Makefile.conf"), "-f", "Makefile.test", "test_manifest"],
                          cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                          timeout=300)
    if done.returncode != 0:
        raise RuntimeError("test_manifest from this tool's table_parser failed: " + done.stderr.strip()[-300:])
    families = family_by_name(workspace_table(ws))
    tests = []
    for line in done.stdout.splitlines():
        parts = line.split("|")
        if len(parts) == 5 and parts[0]:
            tests.append((parts[0], parts[1] or families.get(parts[0], UNASSIGNED), parts[3]))
    return tests


def test_log(tests_dir: Path, name: str) -> Path:
    """The log of one test, named as the test is shown: tests/Mars-ions(valgrind).log for Mars-ions--valgrind."""
    return tests_dir / (display_name(name) + ".log")


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


_ABORT = re.compile(r"^\S+: exit: line=|MPI_ABORT was invoked|exited on signal|Segmentation fault|"
                    r"terminate called after throwing|Floating point exception")


def _step_lines(log: Path, header: str) -> List[str]:
    """The output lines of the last `header` step (e.g. "test_X_run") in a per-test log."""
    try:
        lines = log.read_text(errors="replace").splitlines()
    except OSError:
        return []
    start = max((n for n, line in enumerate(lines) if line.startswith("### ") and line.rstrip().endswith(header)),
                default=-1)
    step: List[str] = []
    for line in lines[start + 1:] if start >= 0 else []:
        if line.startswith("### ") and not line.startswith(("### cwd", "### $")):
            break
        step.append(line)
    return step


def run_error(log: Path, header: str) -> str:
    """How the program of a run step aborted (AMPS exit(), MPI_Abort, a signal), or "".

    AMPS' exit() used to abort with the code 0, so mpirun and make reported success; the check that
    followed compared whatever an unfinished run had left (Periodic-BC: "extra numbers in ...")."""
    found = next((l for l in _step_lines(log, header) if _ABORT.search(l)), "")
    return " ".join(found.split())[:300]


# Compiler and linker diagnostics. "Error: " anywhere in a line also matched the source lines gcc quotes under a
# warning (`  130 | sprintf(msg,"Error: The neutral binary file ...`) and warnings that quote a message
# ("... AMPS:Error: ..."): a test that builds no ./amps (the SWMF tests) was then taken as not built and skipped.
_SOURCE_EXCERPT = re.compile(r"^\s*\d*\s*\|")                         # gcc's quoted source line or its caret line
_COMPILER_ERROR = re.compile(r"^\s*\S+?:\d+(?::\d+)?: (?:fatal )?error: "   # file:line[:col]: error: (gcc, clang)
                             r"|^\s*(?:[\w.+/-]+: )?fatal error: "          # g++: fatal error: ...
                             r"|^\s*Error: "                                # gfortran (after its file:line: line)
                             r"|\*\*\* No rule to make target")
_LINKER_ERROR = re.compile(r"undefined reference to|collect2: error|\bld: cannot find|ld returned \d+ exit status")


def build_error(log: Path, header: str) -> str:
    """The first compiler or linker error of the last `header` step in a per-test log, or ""."""
    step = [l for l in _step_lines(log, header) if not _SOURCE_EXCERPT.match(l)]
    for pattern in (_COMPILER_ERROR, _LINKER_ERROR):
        found = next((l for l in step if pattern.search(l)), "")
        if found:
            return " ".join(found.split())[:300]
    return ""


# The files Config.pl rewrites. A Config.pl option given by one test's compile rule stays in them for every test
# configured after it in the same workspace: -link-option, -f-link-option and -cpp-link-option replace lines of the
# tracked Makefile (which the SWMF tests copy into SWMF/PT/AMPS), -compiler-option appends to Makefile.local, and
# -spice-path, -model-data-path, -cplr-data-path, -openmp, -avx, ... go to Makefile.local and .amps.conf, where
# -application does not reset them. On 2026-11-15 the gtest options of the unit tests had replaced the machine's
# -f-link-option=-lmpi_cxx for every later test, and the Earth tests were built with const-field-bc's model data
# path baked in. Every enabled AMPS-only test configures its own application, so each test starts from the state
# the workspace had after the install.
CONFIG_FILES = ("Makefile", "Makefile.local", "Makefile.conf", ".amps.conf", ".general.conf")


def config_snapshot(ws: Path) -> Dict[str, Optional[bytes]]:
    """The contents of the Config.pl files of ws (None: the file does not exist)."""
    return {name: (ws / name).read_bytes() if (ws / name).is_file() else None for name in CONFIG_FILES}


def config_restore(ws: Path, saved: Dict[str, Optional[bytes]]) -> List[str]:
    """Put the Config.pl files of ws back to `saved`; return the names of the files that had changed.

    Unchanged files are not rewritten, so their modification times (make dependencies) stay."""
    changed = []
    for name, data in saved.items():
        path = ws / name
        if (path.read_bytes() if path.is_file() else None) == data:
            continue
        if data is None:
            path.unlink()
        else:
            path.write_bytes(data)
        changed.append(name)
    return changed


# Compile lanes. On 2026-11-17 the 89 compiles of passion ran one after another for 185 of the runner's 255 minutes
# (make -j8 on 36 cores), and the 4 run lanes then took 69 minutes. A compile changes the tree it runs in (Config.pl,
# srcTemp, ./amps), so tests compile at once only in separate trees: copies of the configured workspace.
# Not copied: what the tests make in the workspace (run directories, result files, the build tree, the executable,
# the SWMF trees, which only the SWMF tests use, and they compile in the workspace itself).
_NOT_COPIED = re.compile(r"^(?:run_test_.*|test_.*\.diff|SWMF|SWMF-IH|srcTemp|amps|\.nightly-manifest)$")
_REWRITE_MAX_BYTES = 1 << 20


def compile_lane(ws: Path, lane: Path) -> List[str]:
    """Make lane a copy of the configured workspace ws, in which the tests compile as in ws; return the files in
    which the workspace's path was replaced by the lane's (relative to lane).

    Symbolic links stay links (SWMF_data, output). The git metadata is shared with ws (no build step uses git).
    Config.pl writes the absolute path of the tree into Makefile.def (PTDIR, MYDIR, used by the Table rules as
    $(PTDIR)/data/input/...): every text file of the copy that holds it gets the lane's path instead."""
    if lane.exists():
        shutil.rmtree(str(lane))
    lane.parent.mkdir(parents=True, exist_ok=True)
    top = str(ws)
    shutil.copytree(top, str(lane), symlinks=True,
                    ignore=lambda where, names: [n for n in names if where == top and _NOT_COPIED.match(n)])
    old, new = top.encode(), str(lane).encode()
    rewritten = []
    for root, dirs, files in os.walk(str(lane)):
        dirs[:] = [d for d in dirs if d != ".git"]
        for name in files:
            path = Path(root, name)
            try:
                if path.is_symlink() or path.stat().st_size > _REWRITE_MAX_BYTES:
                    continue
                data = path.read_bytes()
            except OSError:
                continue
            if old in data and b"\0" not in data:
                path.write_bytes(data.replace(old, new))
                rewritten.append(str(path.relative_to(lane)))
    return sorted(rewritten)


def hand_over(lane: Path, ws: Path, name: str) -> None:
    """Move what the compile and rundir stages of test `name` made in a compile lane into the workspace, where its
    run and check stages run: its run directory (TEST<name>DIR = run_test_<name>) and its result files."""
    source = lane / ("run_test_" + name)
    if source.exists() or source.is_symlink():
        target = ws / source.name
        if target.is_symlink() or target.is_file():
            target.unlink()
        elif target.exists():
            shutil.rmtree(str(target))
        shutil.move(str(source), str(target))
    for path in lane.glob("test_{}*.diff".format(name)):
        if legacy.split_result_file(path.name)[0] == name:
            shutil.move(str(path), str(ws / path.name))


def lane_dirs(ws: Path, count: int) -> List[Path]:
    """The trees of compile lanes 2..count: <runner directory>/lane<k>/AMPS, next to the workspace <runner
    directory>/AMPS, so the orchestrator removes them with it."""
    return [ws.parent / "lane{}".format(k) / ws.name for k in range(2, count + 1)]


def make_lanes(ws: Path, machine: Machine, shareable: int, log: Callable[[str], None]) -> List[Path]:
    """The copies of ws for compile lanes 2..compile_lanes: no more than the compile units that can compile in a
    copy (`shareable`), and fewer when a copy cannot be made (too little space, an error): the tests then compile in
    the lanes there are."""
    lanes: List[Path] = []
    for number, lane in enumerate(lane_dirs(ws, min(machine.compile_lanes, shareable + 1)), 2):
        free = shutil.disk_usage(str(ws)).free / 2 ** 30
        if free < machine.min_free_gb:
            log("Makefile.test: compile lanes {}-{} not made: only {:.1f} GB free (min_free_gb = {:g})".format(
                number, machine.compile_lanes, free, machine.min_free_gb))
            break
        started = time.time()
        try:
            rewritten = compile_lane(ws, lane)
        except OSError as error:                 # shutil.Error included
            log("Makefile.test: compile lane {} not made ({}): the tests compile in {} lane{}".format(
                number, " ".join(str(error).split())[:300], len(lanes) + 1, "" if not lanes else "s"))
            shutil.rmtree(str(lane), ignore_errors=True)
            break
        log("Makefile.test: compile lane {}: {}, a copy of the workspace ({:.0f} s; the workspace's path replaced "
            "in {})".format(number, lane, time.time() - started, ", ".join(rewritten) or "no file"))
        lanes.append(lane)
    return lanes


def run_order(ready: Sequence[str], table: Sequence[TableTest]) -> Tuple[List[str], Dict[str, List[str]]]:
    """The order of phase B and, per test, the tests whose run stage it waits for.

    A Table block that uses the run directory of another test (${TEST<other>DIR}) needs that test's run:
    Moon-restart copies the restart files the run of Moon writes into run_test_Moon (2026-11-18 its run started
    12 s after Moon's in another lane and found none). Such a test waits for the other test's run stage. In the
    order every test comes after the tests it waits for, and the lanes take the tests in this order, so the
    tests a lane waits for have started in other lanes: waiting cannot block them. A dependency cycle is not
    waited for.
    """
    uses = {t.name: t.uses for t in table}
    planned = set(ready)
    order: List[str] = []
    waits: Dict[str, List[str]] = {}
    state: Dict[str, str] = {}

    def visit(name: str) -> None:
        state[name] = "visiting"
        waits[name] = []
        for other in sorted(uses.get(name, ())):
            if other not in planned or state.get(other) == "visiting":
                continue
            if other not in state:
                visit(other)
            waits[name].append(other)
        state[name] = "done"
        order.append(name)

    for name in ready:
        if name not in state:
            visit(name)
    return order, waits


def compile_groups(planned: Sequence[str], table: Sequence[TableTest]) -> List[List[str]]:
    """The planned tests in compile units, each compiled in one tree, one test after another: a test that uses
    another test's run directory (${TEST<x>DIR}) is in the unit of that test, after it.

    Such a pair was written for one tree compiling in Table order. Moon-restart is configured with
    input/moon-restart.input, which includes input/moon.input, which includes moon.input.spicekernels from the
    root of the tree: only Config.pl -application=Moon (the compile of Moon) puts it there. In a compile lane of
    its own Moon-restart did not build ("Cannot find file moon.input.spicekernels", 2026-11-18)."""
    uses = {t.name: t.uses for t in table}
    first = {name: number for number, name in enumerate(planned)}
    parent = {name: name for name in planned}

    def root(name: str) -> str:
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for name in planned:
        for other in uses.get(name, ()):
            if other in parent:
                parent[root(name)] = root(other)
    units: Dict[str, List[str]] = collections.OrderedDict()
    for name in planned:
        units.setdefault(root(name), []).append(name)
    groups = [run_order(members, table)[0] for members in units.values()]
    return sorted(groups, key=lambda group: min(first[n] for n in group))


# OpenMP threads per MPI rank of the Makefile.test runs (unless the environment sets OMP_NUM_THREADS). Without it an
# OpenMP build starts one thread per core in every rank: on 2026-11-15 fast-wave-openmp ran 8 ranks x 38 threads on
# 36 cores next to the other lanes, the load reached 135, and magnetic-reconnection (31 min on a quiet machine) hit the
# 2-hour run limit. 2 threads still exercise the threaded code paths.
OMP_THREADS = 2


def run(ws: Path, machine: Machine, items: Sequence[Item], env: Dict[str, str], logs: Path, tag: str,
        deadline: float, log: Callable[[str], None], config: Optional[Config] = None) -> List[Dict]:
    """Run the planned Makefile.test tests in ws; return result rows (with log_src = the per-test log)."""
    env = dict(env)
    env.setdefault("OMP_NUM_THREADS", str(OMP_THREADS))
    jobs = machine.compile_jobs_for(config)
    make_vars = ["{}={}".format(k, v) for k, v in sorted((machine.make_vars or {}).items())] + \
        ["{}={}".format(COMPILE_JOBS_VAR, jobs)]
    tests_dir = logs / "tests"
    tests_dir.mkdir(parents=True, exist_ok=True)
    limits = machine.timeouts_s
    timing: Dict[str, float] = {}
    over: Dict[str, str] = {}          # test -> "stage exceeded its time limit of N s"
    not_started: List[str] = []
    no_space: Dict[str, str] = {}      # test -> why it was not started (too little free space)
    aborted: Dict[str, str] = {}       # test -> how its run aborted (even when make reported success)
    not_built: Dict[str, str] = {}     # test -> the build error (the Table's "-@($(MAKE) amps)" ignores it)

    def stage(name: str, which: str, tree: Path = ws) -> bool:
        target = "test_{}_{}".format(name, which)
        limit = min(limits[which], max(1.0, deadline - time.time()))
        outcome = procs.run(["make", target] + make_vars, tree, env, test_log(tests_dir, name), limit,
                            "{}/{}".format(tag, name), header=target)
        timing[name] = timing.get(name, 0.0) + outcome.elapsed_s
        if outcome.timed_out:
            over[name] = "{} stage exceeded {} ({:.0f} s)".format(
                which, "the night's time limit" if limit < limits[which] else "its time limit", limit)
        return outcome.exit_code == 0 and not outcome.timed_out

    # the shared library only: the legacy target test_compile also compiles every test
    procs.run(procs.shell("mkdir -p share/lib && cd share/Library/src && make LIB"), ws, env, logs / "setup.log",
              limits["compile"], tag + "/setup", header="shared library (share/Library/src: make LIB)")
    planned = [i.make_name for i in items if i.action == RUN]       # make names (test_<name>_compile, ...)
    table = workspace_table(ws)
    workspace_only = {t.name for t in table if t.workspace_only}
    groups = compile_groups(planned, table)
    only_ws = (lambda group: any(n in workspace_only for n in group))
    trees = [ws] + make_lanes(ws, machine, sum(1 for g in groups if not only_ws(g)), log)
    log("Makefile.test: {} tests to run; compile: {} lane{}, make -j{} each; run: {} lanes".format(
        len(planned), len(trees), "" if len(trees) == 1 else "s", jobs, machine.lanes))
    installed = {tree: config_snapshot(tree) for tree in trees}

    def reset_config(name: str, tree: Path) -> None:
        changed = config_restore(tree, installed[tree])
        if changed:
            with open(str(test_log(tests_dir, name)), "a") as stream:
                stream.write("\n### configuration reset to the post-install state: {} (changed by the tests "
                             "configured before this one)\n".format(", ".join(changed)))

    # one tree: the Table order, as before the compile lanes
    queue = collections.deque(groups if len(trees) > 1 else [[n] for n in planned])
    queue_lock = threading.Lock()

    def next_group(tree: Path) -> Optional[List[str]]:
        """The workspace takes the units that only it can compile first; a copy takes the next other unit."""
        with queue_lock:
            fits = only_ws if tree == ws else (lambda g: not only_ws(g))
            group = next((g for g in queue if fits(g)), None)
            if group is None and tree == ws and queue:
                group = queue[0]
            if group is not None:
                queue.remove(group)
            return group

    built = set()

    def compile_and_rundir(name: str, tree: Path) -> None:
        if time.time() >= deadline:
            not_started.append(name)
            return
        free = shutil.disk_usage(str(ws)).free / 2 ** 30
        if free < machine.min_free_gb:
            no_space[name] = "not started: only {:.1f} GB free (min_free_gb = {:g})".format(free, machine.min_free_gb)
            return
        amps = tree / "amps"
        if amps.exists() or amps.is_symlink():
            amps.unlink()                  # an executable left by an earlier test must not stand in for this one
        reset_config(name, tree)           # nor the Config.pl options an earlier test set
        stage(name, "compile", tree)
        if name not in over and not amps.exists():
            error = build_error(test_log(tests_dir, name), "test_{}_compile".format(name))
            if error:
                not_built[name] = "the build failed: " + error
        if name not in over and name not in not_built:
            stage(name, "rundir", tree)
        if tree != ws:
            hand_over(tree, ws, name)
        if name not in not_built and _rundir_done(ws, name):
            built.add(name)

    def compile_lane_worker(tree: Path) -> None:                # phase A: one test after another in each tree
        while True:
            group = next_group(tree)
            if group is None:
                return
            for name in group:
                try:
                    compile_and_rundir(name, tree)
                except OSError as error:
                    not_built[name] = "the nightly could not hand over its build in {}: {}".format(tree, error)

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(trees)) as pool:
        list(pool.map(compile_lane_worker, trees))
    ready = [name for name in planned if name in built]
    # the run stages that build in the workspace (the SWMF tests copy its tracked Makefile) see the installed state,
    # not the last compiled test's
    config_restore(ws, installed[ws])
    # The check stages run one at a time (the runs of the other lanes go on meanwhile). Every check runs in the
    # workspace, and share/Scripts/DiffNum.pl keeps the text of the two files it compares in the fixed names
    # _diffnum_file1_ and _diffnum_file2_ of the current directory: two checks at once overwrite each other's files
    # and report differences that are not there (2026-11-15: CG-dust-test-tetra-cut-mesh, "34017c34017,34472",
    # while Periodic-BC-NonUniform's check ran).
    check_lock = threading.Lock()
    order, waits = run_order(ready, table)
    for name in order:
        if waits[name]:
            log("Makefile.test: {} runs after {} (it uses their run directories)".format(name, ", ".join(waits[name])))
    run_done = {name: threading.Event() for name in order}
    with concurrent.futures.ThreadPoolExecutor(max_workers=machine.lanes) as pool:   # phase B: lanes
        def run_and_check(name: str) -> None:
            try:
                for other in waits[name]:
                    if not run_done[other].is_set():
                        with open(str(test_log(tests_dir, name)), "a") as stream:
                            stream.write("\n### {} waiting for test_{}_run: this test uses run_test_{}\n".format(
                                time.strftime("%Y-%m-%d %H:%M:%S"), other, other))
                        run_done[other].wait(max(0.0, deadline - time.time()))
                if time.time() >= deadline:
                    not_started.append(name)
                    return
                stage(name, "run")
            finally:
                run_done[name].set()          # also when it did not start: nothing waits forever
            error = run_error(test_log(tests_dir, name), "test_{}_run".format(name))
            if error:
                aborted[name] = "the run aborted: " + error
            with check_lock:
                stage(name, "check")
        list(pool.map(run_and_check, order))

    results = {r["test"]: r for r in legacy.collect_results(legacy.read_diff_files(ws), workspace_table(ws))}
    rows = []
    for item in items:
        if item.action != RUN:
            continue
        name = item.make_name
        result = results.pop(item.test, None)          # collect_results names the tests as they are shown
        if name in not_started:
            row = {"status": "SKIPPED", "stage": "", "reason": "not started: the night's time limit was reached"}
        elif name in no_space:
            row = {"status": "SKIPPED", "stage": "", "reason": no_space[name]}
        elif name in not_built:
            row = {"status": "FAILED", "stage": "compile", "reason": not_built[name]}
        elif name in aborted and name not in over:
            row = {"status": "FAILED", "stage": "run", "reason": aborted[name]}
        elif result is None:
            row = {"status": "FAILED", "stage": "compile", "reason": "no result file test_{}*.diff was written".format(name)}
        else:
            row = {"status": result["status"], "stage": result.get("stage", ""), "reason": result.get("reason", "")}
            if name in over:
                row.update(status="FAILED", stage=over[name].split()[0], reason=over[name])
        row.update(test=item.test, family=item.family, duration_s=round(timing.get(name, 0.0), 1),
                   log_src=str(test_log(tests_dir, name)))
        rows.append(row)
    if results:
        log("Makefile.test: ignored result files of tests that were not run: " + ", ".join(sorted(results)))
    return rows
