"""`amps-nightly run`: the nightly test procedure on one machine (ORCHESTRATOR.md).

For each enabled configuration of the machine, the runners of runners.json run one after another
(phase 1; max_parallel > 1 is reserved for the concurrent scheduler of phase 2). Each runner gets its own
workspace (workspace.py), its commands run under a time limit with every process tracked (procs.py),
its report is read into result rows (readers.py, makefile_runner.py), and its logs are copied into the
night's run directory. `amps-nightly pack` (transfer.py) then sends each run directory to the server:

    {workdir}/<COMPILER>/<runner>/      one folder per compiler, as the old ~/Tmp_AMPS_test/GNU/AMPS:
        AMPS/                           the runner's installed copy (worktree + share, util, BATL, SWMF_data link;
                                        Config.pl calls of the machine file), replaced the next night
        out/  tmp/                      the runner's output directory and TMPDIR
    {workdir}/nights/<YYYYMMDD-HHMM>/
        orchestrator.log
        runs/<column>/       run.json results.jsonl runners.json logs/<runner>/... DONE   (packed and sent)

Every step's output in a log is cut to its first 256 KB and last 768 KB (trim_log): some tests print
hundreds of thousands of lines, and the start and the end of a step are what explain a failure.
"""

from __future__ import annotations

import collections
import datetime as dt
import fcntl
import fnmatch
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

from . import makefile_runner, procs, transfer, workspace
from .machines import Config, Machine
from .plan import RUN, SKIP, Item, SourcePlan, apply_rules, machine_vars, plan
from .readers import READERS, ReadContext, ResultsUnavailable, one_line
from .records import make_run_id, unique_rows, utc_text, write_run_dir
from .registry import UNASSIGNED, Registry
from .runners import Routine, Source

BUILTIN_PLACEHOLDERS = {"python", "ws", "out", "tmp", "jobs", "mpiexec", "table_host"}
MAX_LOG_BYTES = 50 * 2 ** 20                     # larger logs keep their first 10 MB and last 40 MB
TRIM_HEAD, TRIM_TAIL = 256 * 2 ** 10, 768 * 2 ** 10   # bytes kept of each step of a log (trim_log)
_PLACEHOLDER = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")
_KEEP_ENV = ("HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ", "SSH_AUTH_SOCK")
WORKSPACE_MARK = ".amps-nightly"     # in <workdir>/<COMPILER>/<runner>/: only such folders are ever deleted


class NightError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- placeholders
def _texts(source: Source) -> List[str]:
    texts: List[str] = []
    for command in list(source.prepare) + ([source.command] if source.command else []):
        texts.extend([command] if isinstance(command, str) else command)
    texts.extend(source.env.values())
    texts.extend(source.logs)
    texts.append(source.results.get("file", ""))
    return texts


def placeholders(source: Source) -> set:
    return {name for text in _texts(source) for name in _PLACEHOLDER.findall(text)}


def fill(text: str, values: Dict[str, str], quote: bool = False) -> str:
    def replace(match: "re.Match[str]") -> str:
        name = match.group(1)
        if name not in values:
            raise NightError("unknown placeholder {{{}}}".format(name))
        return shlex.quote(values[name]) if quote else values[name]
    return _PLACEHOLDER.sub(replace, text)


def command_argv(command: Any, values: Dict[str, str]) -> List[str]:
    if isinstance(command, str):
        return procs.shell(fill(command, values, quote=True))
    return [fill(part, values) for part in command]


# ----------------------------------------------------------------------------- logs
def copy_log(source: Path, target: Path) -> None:
    """Copy a log into the archive; a log over MAX_LOG_BYTES keeps its head and tail."""
    target.parent.mkdir(parents=True, exist_ok=True)
    size = source.stat().st_size
    if size <= MAX_LOG_BYTES:
        shutil.copyfile(str(source), str(target))
        return
    head, tail = 10 * 2 ** 20, 40 * 2 ** 20
    with open(str(source), "rb") as reader, open(str(target), "wb") as writer:
        writer.write(reader.read(head))
        writer.write("\n\n[... {} bytes omitted by the nightly archive ...]\n\n".format(size - head - tail).encode())
        reader.seek(size - tail)
        shutil.copyfileobj(reader, writer)


def trim_log(path: Path, head: int = TRIM_HEAD, tail: int = TRIM_TAIL) -> int:
    """Cut each step of a log to its first `head` and last `tail` bytes; return the bytes removed.

    A step starts at a "### " header line (procs.run writes one per command); a log without headers is
    one step. Lines are read in pieces of at most 64 KB, so even a runaway log is never read whole.
    """
    path = Path(path)
    if path.stat().st_size <= head + tail:
        return 0
    part = path.with_name("." + path.name + ".trim")
    removed = 0
    with open(str(path), "rb") as src, open(str(part), "wb") as out:
        kept, lines, held, skipped, skipped_lines = 0, collections.deque(), 0, 0, 0
        for line in iter(lambda: src.readline(65536), b""):
            if line.startswith(b"### "):
                if skipped:
                    out.write("[... amps-nightly left out {} lines ({:.1f} MB) of this step ...]\n".format(
                        skipped_lines, skipped / 2 ** 20).encode())
                    removed += skipped
                out.writelines(lines)
                out.write(line)
                kept, lines, held, skipped, skipped_lines = 0, collections.deque(), 0, 0, 0
            elif kept < head:
                out.write(line)
                kept += len(line)
            else:
                lines.append(line)
                held += len(line)
                while held > tail:
                    old = lines.popleft()
                    held -= len(old)
                    skipped += len(old)
                    skipped_lines += 1
        if skipped:
            out.write("[... amps-nightly left out {} lines ({:.1f} MB) of this step ...]\n".format(
                skipped_lines, skipped / 2 ** 20).encode())
            removed += skipped
        out.writelines(lines)
    if removed:
        os.replace(str(part), str(path))
    else:
        part.unlink()
    return removed


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+=-]+", "_", name)[:120] or "log"


# ----------------------------------------------------------------------------- the night
class Night:
    def __init__(self, machine: Machine, routine: Routine, registry: Registry, configs: Sequence[Config],
                 only: Sequence[str] = (), echo: Callable[[str], None] = print,
                 now: Optional[dt.datetime] = None, tests: Sequence[str] = (), keep_local: bool = False,
                 commit: str = "", note: str = ""):
        self.machine, self.routine, self.registry = machine, routine, registry
        self.configs, self.only, self.echo = list(configs), set(only), echo
        self.tests = list(tests)          # --test globs: limit the Makefile.test runner (trials and debugging)
        self.keep_local = keep_local      # --no-upload: `amps-nightly pack` never sends these results
        self.pin = commit                 # --commit: test this commit instead of the branch head
        self.note = note                  # --note: added to every run's notes (shown on the night's page)
        self.started = (now or dt.datetime.now(dt.timezone.utc)).astimezone(dt.timezone.utc)
        self.local_start = self.started.astimezone()
        self.workdir = Path(machine.workdir).expanduser()
        self.root = self.workdir / "nights" / self.local_start.strftime("%Y%m%d-%H%M")
        self.deadline = time.time() + machine.night_limit_h * 3600
        self.mirrors = workspace.Mirrors(self.workdir / "mirrors", routine.checkout)
        self.log_file = self.root / "orchestrator.log"

    def runner_dir(self, config: Config, runner: str) -> Path:
        """<workdir>/<COMPILER>/<runner>: the runner's installed copy (AMPS/), its out/ and tmp/."""
        return self.workdir / config.folder / runner

    def free_gb(self) -> float:
        return shutil.disk_usage(str(self.workdir)).free / 2 ** 30

    def log(self, message: str) -> None:
        line = "{} {}".format(dt.datetime.now().strftime("%H:%M:%S"), message)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(str(self.log_file), "a") as stream:
            stream.write(line + "\n")
        self.echo(line)

    # ------------------------------------------------------------------ environment
    def base_env(self, config: Config, tmp: Path) -> Dict[str, str]:
        env = {k: os.environ[k] for k in _KEEP_ENV if k in os.environ}
        env["PATH"] = ":".join(self.machine.search_path)
        env["TMPDIR"] = str(tmp)
        return env

    def values(self, ws: Path, out: Path, tmp: Path, source: Source) -> Dict[str, str]:
        values = dict(machine_vars(self.machine))
        if self.machine.core_profile:
            values["core_profile"] = str(ws / self.machine.core_profile)
        values.update(python=self.machine.python, ws=str(ws), out=str(out), tmp=str(tmp),
                      jobs=str(min(source.jobs or self.machine.cores, self.machine.cores)),
                      mpiexec=self.machine.mpiexec, table_host=self.machine.table_host)
        return values

    # ------------------------------------------------------------------ main
    def run(self) -> List[Path]:
        self.workdir.mkdir(parents=True, exist_ok=True)
        lock = open(str(self.workdir / "lock"), "a")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            lock.close()
            raise NightError("another amps-nightly run is still running on this machine ({})".format(self.workdir / "lock"))
        running = self.workdir / transfer.RUNNING_FILE
        try:
            running.write_text(self.root.name + "\n")
            self._clean_old_nights()
            self.log("night {} on {}: configurations {}; runners {}".format(
                self.root.name, self.machine.machine, ", ".join(c.config for c in self.configs),
                ", ".join(s.id for s in self.routine.sources if s.enabled and (not self.only or s.id in self.only))))
            checkout_error = ""
            try:
                linked = {c.path: machine_vars(self.machine).get(c.link_var, "")
                          for c in self.routine.checkout.companions if c.link_var}
                self.mirrors.update(self.log, self.pin, {p: d for p, d in linked.items() if d})
            except (workspace.WorkspaceError, subprocess.SubprocessError, OSError) as error:
                checkout_error = str(error)
                self.log("ERROR: checkout failed: " + checkout_error)
            done = [self.run_config(config, checkout_error) for config in self.configs]
            self.log("night finished; 'amps-nightly pack' sends the results to the server")
            return done
        finally:
            if running.exists():
                running.unlink()
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()

    def _clean_old_nights(self) -> None:
        # the previous night's installed copies (the old routine started with rm -rf GNU)
        for config in self.machine.configs:
            folder = self.workdir / config.folder
            for old in sorted(folder.iterdir()) if folder.is_dir() else []:
                if (old / WORKSPACE_MARK).is_file():
                    shutil.rmtree(str(old), ignore_errors=True)
        nights = self.workdir / "nights"
        for old in sorted(nights.iterdir()) if nights.is_dir() else []:
            if not old.is_dir() or old == self.root:
                continue
            if any(r.is_dir() and not (r / transfer.PACKED_FILE).exists() for r in (old / "runs").glob("*")):
                # results 'amps-nightly pack' has not sent yet: keep them (and the log), drop the workspaces
                for child in old.iterdir():
                    if child.name not in ("runs", "orchestrator.log"):
                        shutil.rmtree(str(child), ignore_errors=True) if child.is_dir() else child.unlink()
            elif self.machine.keep_workspaces != "always":
                shutil.rmtree(str(old), ignore_errors=True)
        for mirror in (self.workdir / "mirrors").glob("*.git"):
            subprocess.run(["git", "worktree", "prune"], cwd=str(mirror), stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)

    def run_config(self, config: Config, checkout_error: str) -> Path:
        run_id = make_run_id(self.machine.machine, config.config, self.started)
        rundir = self.root / "runs" / config.column
        rundir.mkdir(parents=True, exist_ok=True)
        rows: List[Dict[str, Any]] = []
        runners: List[Dict[str, Any]] = []
        self.log("--- {} ({})".format(config.column, run_id))
        (rundir / transfer.PARTIAL_RUN_FILE).write_text(json.dumps(
            self._run_record(config, run_id, self.started, False), indent=2) + "\n")
        if checkout_error:
            rows.append({"family": UNASSIGNED, "test": "checkout", "status": "FAILED", "stage": "setup",
                         "reason": checkout_error, "source": "orchestrator"})
        else:
            plans = plan(self.machine, config, self.routine, self.registry, self.local_start.date())
            for sp in plans:
                if self.only and sp.source.id not in self.only:
                    continue
                info: Dict[str, Any] = {"runner": sp.source.id, "title": sp.source.title}
                if not sp.active:
                    info.update(state="not run", note=sp.note)
                    self.log("{}: not run ({})".format(sp.source.id, sp.note))
                elif sp.items and all(i.action == SKIP for i in sp.items):
                    rows.extend(self._skipped(sp, sp.items[0].reason))
                    info.update(state="skipped", note=sp.items[0].reason)
                    self.log("{}: skipped ({})".format(sp.source.id, sp.items[0].reason))
                elif time.time() >= self.deadline:
                    reason = "not started: the night's time limit of {:g} h was reached".format(self.machine.night_limit_h)
                    rows.extend(self._skipped(sp, reason))
                    info.update(state="skipped", note=reason)
                    self.log("{}: {}".format(sp.source.id, reason))
                elif self.free_gb() < self.machine.min_free_gb:
                    reason = "not started: only {:.1f} GB free in {} (min_free_gb = {:g})".format(
                        self.free_gb(), self.workdir, self.machine.min_free_gb)
                    rows.extend(self._skipped(sp, reason))
                    info.update(state="skipped", note=reason)
                    self.log("{}: {}".format(sp.source.id, reason))
                else:
                    new_rows = self.run_runner(sp, config, run_id, rundir, info)
                    rows.extend(new_rows)
                runners.append(info)
                self._write_progress(rundir, rows, runners)
        run = self._run_record(config, run_id, dt.datetime.now(dt.timezone.utc), not checkout_error)
        rows = unique_rows(rows)
        (rundir / "logs").mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(self.log_file), str(rundir / "logs" / "orchestrator.log"))
        (rundir / "runners.json").write_text(json.dumps(runners, indent=2) + "\n")
        write_run_dir(rundir, run, rows, done=True)
        for name in (transfer.PROGRESS_FILE, transfer.PARTIAL_RUN_FILE):
            if (rundir / name).exists():
                (rundir / name).unlink()
        if self.keep_local:
            (rundir / transfer.PACKED_FILE).write_text("not sent: kept on this machine (amps-nightly run --no-upload)\n")
        self.log("{}: {} results ({} failed) in {}{}".format(
            config.column, len(rows), sum(1 for r in rows if r["status"] == "FAILED"), rundir,
            " (kept local: --no-upload)" if self.keep_local else ""))
        return rundir

    def _run_record(self, config: Config, run_id: str, finished: dt.datetime, complete: bool) -> Dict[str, Any]:
        return {"schema_version": 1, "run_id": run_id, "machine": self.machine.machine, "config": config.config,
                "compiler": config.compiler, "toolchain": config.toolchain, "branch": self.routine.checkout.branch,
                "commit": self.mirrors.commits.get("", ""), "commit_inferred": False,
                "started": utc_text(self.started), "finished": utc_text(finished), "complete": complete,
                "producer": "amps-nightly run", "night": self.local_start.strftime("%Y-%m-%d"),
                "notes": "; ".join(["{} {}".format(p, c[:12]) for p, c in self.mirrors.commits.items() if p]
                                   + ([self.note] if self.note else []))}

    # ------------------------------------------------------------------ one runner
    def run_runner(self, sp: SourcePlan, config: Config, run_id: str, rundir: Path, info: Dict[str, Any]) -> List[Dict]:
        source = sp.source
        home = self.runner_dir(config, source.id)
        ws, out, tmp = home / "AMPS", home / "out", home / "tmp"
        if home.exists():
            workspace.remove(ws, self.mirrors)          # also forgets its git worktrees
            shutil.rmtree(str(home), ignore_errors=True)
        home.mkdir(parents=True, exist_ok=True)
        (home / WORKSPACE_MARK).write_text("installed copy of runner {} for night {}; replaced by the next night\n".format(
            source.id, self.root.name))
        logs = rundir / "logs" / source.id
        runner_log = logs / "runner.log"
        for path in (out, tmp, logs):
            path.mkdir(parents=True, exist_ok=True)
        tag = "{}/{}".format(run_id, source.id)
        env = self.base_env(config, tmp)
        values = self.values(ws, out, tmp, source)
        env.update({k: fill(v, values) for k, v in source.env.items()})
        end = min(time.time() + source.time_limit_min * 60, self.deadline)
        began = time.time()
        info.update(state="ran", started=utc_text(dt.datetime.now(dt.timezone.utc)), workspace=str(ws),
                    time_limit_min=source.time_limit_min)
        self.log("{}: start (limit {} min, workspace {})".format(source.id, source.time_limit_min, ws))
        rows: List[Dict[str, Any]] = []
        outcome = None
        stopped = ""                                    # why the runner did not get to report normally
        try:
            workspace.create(ws, self.mirrors, self.machine, values, self.log)
            if source.workspace == "configured":
                ok, failed = workspace.configure(ws, self.machine, config, env, runner_log, tag + "/configure",
                                                 max(60.0, end - time.time()))
                if not ok:
                    stopped = "configuring the workspace failed: " + failed
        except (workspace.WorkspaceError, OSError, subprocess.SubprocessError) as error:
            stopped = "workspace setup failed: {}".format(error)
        items = list(sp.items)
        if stopped:
            rows = self._failed(sp, "setup", stopped, runner_log)
        elif source.kind == "makefile":
            try:
                found = makefile_runner.manifest(ws, self.machine, env, self.log)
                items = [Item(family or UNASSIGNED, name, SKIP if code == "SKIP" else RUN,
                              "ExeptionCode=SKIP in MakefileTest/Table for this compiler" if code == "SKIP" else "")
                         for name, family, code in found]
                apply_rules(items, source, self.machine, config, self.routine, self.local_start.date())
                if self.tests:
                    items = [i for i in items if any(fnmatch.fnmatchcase(i.test, g) for g in self.tests)]
                rows = makefile_runner.run(ws, self.machine, items, env, logs, tag, end, self.log)
                rows += [{"family": i.family, "test": i.test, "status": "SKIPPED", "reason": i.reason, "log_src": ""}
                         for i in items if i.action == SKIP]
            except (RuntimeError, OSError, subprocess.SubprocessError) as error:
                rows = self._failed(sp, "setup", "cannot list the Makefile.test tests: {}".format(error), runner_log)
        else:
            for command in source.prepare:
                outcome = procs.run(workspace.with_setup(command_argv(command, values), config.setup), ws, env,
                                    runner_log, max(1.0, end - time.time()), tag, header="prepare")
                if outcome.exit_code != 0 or outcome.timed_out:
                    stopped = "build step failed ({}): {}".format(
                        "time limit" if outcome.timed_out else "exit {}".format(outcome.exit_code),
                        command if isinstance(command, str) else " ".join(command))
                    break
            if stopped:
                rows = self._failed(sp, "build", stopped, runner_log)
            else:
                outcome = procs.run(workspace.with_setup(command_argv(source.command, values), config.setup), ws, env,
                                    runner_log, max(1.0, end - time.time()), tag, header="command")
                rows = self._read(sp, values, outcome, runner_log)
        procs.sweep(tag)
        rows = self._finish_rows(sp, rows, config, runner_log)
        self._archive(source, values, rows, rundir, logs)
        failed = sum(1 for r in rows if r["status"] == "FAILED")
        info.update(finished=utc_text(dt.datetime.now(dt.timezone.utc)), elapsed_s=round(time.time() - began, 1),
                    exit_code=outcome.exit_code if outcome else None, timed_out=bool(outcome and outcome.timed_out),
                    cpu_s=round(outcome.cpu_s, 1) if outcome else None,
                    peak_rss_mb=round(outcome.peak_rss_mb, 1) if outcome else None, problem=stopped,
                    counts={s: sum(1 for r in rows if r["status"] == s) for s in ("PASS", "FAILED", "SKIPPED")},
                    log="logs/{}/runner.log".format(source.id))
        keep = self.machine.keep_workspaces == "always" or (self.machine.keep_workspaces == "failures" and failed)
        if not keep:
            workspace.remove(ws, self.mirrors)
            shutil.rmtree(str(home), ignore_errors=True)
        info["workspace_kept"] = bool(keep)
        self.log("{}: {} PASS, {} FAILED, {} SKIPPED in {:.0f} min{}".format(
            source.id, info["counts"]["PASS"], failed, info["counts"]["SKIPPED"], info["elapsed_s"] / 60,
            "; " + stopped if stopped else ""))
        return rows

    def _read(self, sp: SourcePlan, values: Dict[str, str], outcome: procs.Outcome, runner_log: Path) -> List[Dict]:
        source = sp.source
        ctx = ReadContext(runner=source.id, file=fill(source.results.get("file", ""), values),
                          out=Path(values["out"]), exit_code=outcome.exit_code)
        try:
            rows = READERS[source.results["reader"]](ctx)
        except ResultsUnavailable as error:
            last = _last_error(runner_log)
            reason = "{}: {}".format("time limit of {} min exceeded".format(source.time_limit_min) if outcome.timed_out
                                     else "exited with {}".format(outcome.exit_code), last or "no report ({})".format(error))
            return self._failed(sp, "run", reason, runner_log, listed_only=False)
        if outcome.timed_out:
            rows.append({"test": source.id + ": time limit", "status": "FAILED", "stage": "run",
                         "reason": "time limit of {} min exceeded; the results above are partial".format(source.time_limit_min),
                         "log_src": str(runner_log)})
        elif outcome.exit_code not in (0, None) and not any(r["status"] == "FAILED" for r in rows):
            rows.append({"test": source.id + ": exit status", "status": "FAILED", "stage": "run",
                         "reason": "exited with {} although no test failed: {}".format(
                             outcome.exit_code, _last_error(runner_log) or "see the runner log"),
                         "log_src": str(runner_log)})
        reported = {r["test"] for r in rows}
        for item in sp.items:
            if item.listed and item.action == RUN and item.test not in reported:
                rows.append({"test": item.test, "family": item.family, "status": "FAILED", "stage": "run",
                             "reason": "no result (the runner stopped early)", "log_src": str(runner_log)})
        return rows

    # ------------------------------------------------------------------ rows
    def _skipped(self, sp: SourcePlan, reason: str) -> List[Dict]:
        items = [i for i in sp.items if i.listed] or [Item(sp.source.family or UNASSIGNED, sp.source.id, SKIP)]
        return [{"family": i.family, "test": i.test, "status": "SKIPPED", "reason": i.reason or reason,
                 "source": sp.source.id} for i in items]

    def _failed(self, sp: SourcePlan, stage: str, reason: str, log: Path, listed_only: bool = True) -> List[Dict]:
        items = [i for i in sp.items if i.listed and i.action == RUN]
        if not items:
            items = [Item(sp.source.family or UNASSIGNED, sp.source.id, RUN)]
        rows = [{"family": i.family, "test": i.test, "status": "FAILED", "stage": stage, "reason": reason,
                 "log_src": str(log)} for i in items]
        return rows + [{"family": i.family, "test": i.test, "status": "SKIPPED", "reason": i.reason, "log_src": ""}
                       for i in sp.items if i.listed and i.action == SKIP]

    def _finish_rows(self, sp: SourcePlan, rows: List[Dict], config: Config, runner_log: Path) -> List[Dict]:
        """Families, switched-off tests and the record fields of every row of one runner."""
        source = sp.source
        today = self.local_start.date()
        for r in rows:
            r["family"] = r.get("family") or source.family or self.registry.family_for_suite(r["test"])
            r["source"] = source.id
            entry = self.routine.disabled_entry(source.id, r["test"], self.machine.machine, config.config, today)
            if entry and r["status"] != "SKIPPED":
                r.update(status="SKIPPED", stage="", reason=entry.describe() + " (its result was ignored)")
            r.setdefault("stage", "")
            r.setdefault("reason", "")
        return rows

    def _archive(self, source: Source, values: Dict[str, str], rows: List[Dict], rundir: Path, logs: Path) -> None:
        """Copy the runner's logs into the run directory and point every row at its log."""
        out = Path(values["out"])
        copied: Dict[str, str] = {}
        for pattern in source.logs:
            for match in sorted(Path("/").glob(fill(pattern, values).lstrip("/"))):
                if match.is_file():
                    try:
                        relative = match.relative_to(out)
                    except ValueError:
                        relative = Path(match.name)
                    target = logs / "files" / relative
                    copy_log(match, target)
                    copied[str(match)] = str(target.relative_to(rundir))
        for r in rows:
            src = r.pop("log_src", "") or ""
            path = Path(src) if src else None
            if path is not None and path.is_file():
                try:
                    r["log"] = str(path.relative_to(rundir))
                except ValueError:
                    if str(path) not in copied:
                        target = logs / "tests" / (_safe(path.stem) + path.suffix)
                        n = 2
                        while target.exists():
                            target = logs / "tests" / "{}.{}{}".format(_safe(path.stem), n, path.suffix)
                            n += 1
                        copy_log(path, target)
                        copied[str(path)] = str(target.relative_to(rundir))
                    r["log"] = copied[str(path)]
            else:
                r["log"] = str((logs / "runner.log").relative_to(rundir)) if (logs / "runner.log").exists() else ""
        removed = sum(trim_log(p) for p in sorted(logs.rglob("*")) if p.is_file())
        if removed:
            self.log("{}: left out {:.0f} MB from the middle of long steps in its logs".format(source.id, removed / 2 ** 20))

    def _write_progress(self, rundir: Path, rows: List[Dict], runners: List[Dict]) -> None:
        """Rows and runner records so far, readable if the orchestrator itself dies (no DONE yet)."""
        with open(str(rundir / transfer.PROGRESS_FILE), "w") as stream:
            for r in rows:
                stream.write(json.dumps({k: v for k, v in r.items() if k != "log_src"}, sort_keys=True) + "\n")
        (rundir / "runners.json").write_text(json.dumps(runners, indent=2) + "\n")


def _last_error(log: Path) -> str:
    try:
        lines = log.read_text(errors="replace").splitlines()[-300:]
    except OSError:
        return ""
    for line in reversed(lines):
        if ("ERROR" in line or "Error" in line or "Traceback" in line) and not line.startswith("###"):
            return one_line(line, 200)
    return ""


def describe(machine: Machine, routine: Routine, registry: Registry, configs: Sequence[Config],
             only: Sequence[str] = ()) -> List[str]:
    """The commands a night would run, with placeholders filled (amps-nightly run --dry-run)."""
    lines = []
    night = Night(machine, routine, registry, configs, only, echo=lambda _: None)
    for config in configs:
        lines.append("{} (column {}): each runner is installed in {}/<runner>/AMPS".format(
            config.config, config.column, night.workdir / config.folder))
        for sp in plan(machine, config, routine, registry, night.local_start.date()):
            if only and sp.source.id not in only:
                continue
            if not sp.active:
                lines.append("  {:<10} not run: {}".format(sp.source.id, sp.note))
                continue
            if sp.items and all(i.action == SKIP for i in sp.items):
                lines.append("  {:<10} skipped: {}".format(sp.source.id, sp.items[0].reason))
                continue
            home = night.runner_dir(config, sp.source.id)
            values = night.values(home / "AMPS", home / "out", home / "tmp", sp.source)
            lines.append("  {:<10} {} {} copy, limit {} min".format(sp.source.id, home / "AMPS", sp.source.workspace,
                                                                  sp.source.time_limit_min))
            if sp.source.workspace == "configured":
                for call in config.config_pl:
                    lines.append("      install: ./Config.pl " + call)
                for line in config.after_install:
                    lines.append("      install: " + line)
            if sp.source.kind == "makefile":
                lines.append("      make test_<X>_compile/_rundir (one at a time), _run/_check ({} lanes) "
                             "for {} tests".format(machine.lanes, sum(1 for i in sp.items if i.action == RUN)))
                continue
            for command in sp.source.prepare:
                lines.append("      prepare: " + " ".join(shlex.quote(a) for a in command_argv(command, values)))
            lines.append("      command: " + " ".join(shlex.quote(a) for a in command_argv(sp.source.command, values)))
            lines.append("      results: {} {}".format(sp.source.results["reader"], fill(sp.source.results["file"], values)))
    return lines
