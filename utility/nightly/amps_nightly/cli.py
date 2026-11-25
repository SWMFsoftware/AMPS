"""amps-nightly command line.

    amps-nightly install [--root DIR]            put amps_nightly_test.sh in ~/bin and its line in crontab
    amps-nightly uninstall                       take them out again
    amps-nightly check                           Table, families.json, runners.json and machines/ agree
    amps-nightly machines                        machine files, their columns and last uploads
    amps-nightly plan [--machine M] [--list]     what the runners run (and skip, and why) on a machine
    amps-nightly run [--dry-run] [--only R]      script 1: the nightly test procedure on this machine (runners.json)
    amps-nightly pack [--no-send]                script 2: tar each finished run with its logs, send it to the server
    amps-nightly check-table                     every Table block has one registered Family=
    amps-nightly import-legacy-archive ...       server archive (AMPS_TEST_RESULTS) -> run directories
    amps-nightly import-legacy-run ...           legacy work tree (*.diff) -> one run directory
    amps-nightly ingest [TARBALL ...]            incoming/ (or the given tarballs) -> results.db + www/archive/
    amps-nightly report                          results.db -> www/index.html + www/family/
    amps-nightly update [TARBALL ...]            script 3: ingest, then report, then prune
    amps-nightly prune                           gzip archived logs after 14 nights, delete them after 400 days
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import shutil
import socket
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from . import install, legacy, report, store as store_mod, transfer
from .checks import check_all
from .machines import declared_columns, find_machine, load_machines
from .orchestrator import Night, NightError, describe
from .plan import RUN, SKIP, duplicates, plan
from .records import RecordError, parse_time
from .registry import Registry, RegistryError
from .runners import Routine, RoutineError
from .table import check_table, read_table

FALLBACK_ROOT = Path.home() / "AMPSTEST" / "nightly"


def server_root(machines_dir: Optional[str] = None) -> Path:
    """$AMPS_NIGHTLY_ROOT, else the "server_root" of this host's machine file, else ~/AMPSTEST/nightly."""
    if os.environ.get("AMPS_NIGHTLY_ROOT"):
        return Path(os.environ["AMPS_NIGHTLY_ROOT"]).expanduser()
    machines, _ = load_machines(machines_dir)
    machine = find_machine(machines, hostname=socket.gethostname()) or find_machine(machines, hostname=socket.getfqdn())
    if machine is not None and machine.server_root:
        return Path(machine.server_root).expanduser()
    return FALLBACK_ROOT


def _paths(args: argparse.Namespace) -> Tuple[Path, Path]:
    root = Path(args.root).expanduser() if args.root else server_root(getattr(args, "machines", None))
    www = Path(args.www).expanduser() if getattr(args, "www", None) else root / "www"
    return root, www


def _folder_map(items: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    mapping = {}
    for item in items:
        old, sep, new = item.partition("=")
        machine, slash, config = new.partition("/")
        if not sep or not slash or not old or not machine or not config:
            raise SystemExit("--map expects LEGACY_FOLDER=MACHINE/CONFIG, got {!r}".format(item))
        mapping[old] = (machine, config)
    return mapping


def _date(text: Optional[str]) -> Optional[dt.date]:
    return dt.datetime.strptime(text, "%Y-%m-%d").date() if text else None


# ----------------------------------------------------------------------------- commands
def cmd_check_table(args: argparse.Namespace) -> int:
    registry = Registry.load(args.families)
    tests = read_table(args.table)
    problems = check_table(tests, registry)
    for problem in problems:
        print(problem)
    print("{} Table blocks, {} problem{}".format(len(tests), len(problems), "" if len(problems) == 1 else "s"))
    return 1 if problems else 0


def cmd_check(args: argparse.Namespace) -> int:
    errors, warnings, summary = check_all(args.table, args.families, args.runners, args.machines)
    for line in summary:
        print(line)
    for line in warnings:
        print("warning:", line)
    for line in errors:
        print("ERROR:", line)
    print("OK" if not errors else "{} error{}".format(len(errors), "" if len(errors) == 1 else "s"))
    return 1 if errors else 0


def cmd_machines(args: argparse.Namespace) -> int:
    machines, problems = load_machines(args.machines)
    if args.var:
        # for scripts: only the value on stdout (amps_nightly_test.sh reads swmf_data this way)
        machine = (find_machine(machines, args.machine) if args.machine else
                   find_machine(machines, hostname=socket.gethostname()) or
                   find_machine(machines, hostname=socket.getfqdn()))
        if machine is None:
            print("no machine file for {}".format(args.machine or socket.gethostname()), file=sys.stderr)
            return 1
        value = machine.vars.get(args.var)
        if value is None:
            print("{} defines no var {!r}".format(machine.path, args.var), file=sys.stderr)
            return 1
        print(value)
        return 0
    root, _ = _paths(args)
    last: Dict[str, str] = {}
    if (root / "results.db").exists():
        with store_mod.Store(root / "results.db") as st:
            last = {c["column_id"]: c["last_finished"] for c in st.columns()}
    for m in machines:
        print("{}{}  ({})".format(m.machine, "" if m.enabled else "  [switched off]", m.path))
        if m.description:
            print("    " + m.description)
        print("    hosts {} · Table host {} · starts {} · {} lane{} · uploads to {}".format(
            ", ".join(m.hostnames), m.table_host, m.start, m.lanes, "" if m.lanes == 1 else "s", m.upload))
        print("    Makefile.test compiles: {} at once, make -j{}{}".format(
            m.compile_lanes, m.compile_jobs_for(), "".join(
                " ({} -j{})".format(c.config, c.compile_jobs) for c in m.configs if c.compile_jobs)))
        print("    capabilities: {} · src/test profile: {}".format(", ".join(m.capabilities) or "none",
                                                                  m.core_profile or "none"))
        for c in m.configs:
            print("    column {:<24} {}{}  last upload: {}{}".format(
                c.column, c.compiler, " · " + c.toolchain if c.toolchain else "",
                last.get(c.column, "never"), "" if c.enabled else "  [switched off]"))
    for problem in problems:
        print("ERROR:", problem)
    if not machines and not problems:
        print("no machine files in", args.machines or "utility/nightly/machines")
    return 1 if problems else 0


def _machine(args: argparse.Namespace):
    machines, problems = load_machines(args.machines)
    for problem in problems:
        print("ERROR:", problem)
    if args.machine:
        machine = find_machine(machines, args.machine)
    else:
        machine = (find_machine(machines, hostname=socket.gethostname())
                   or find_machine(machines, hostname=socket.getfqdn()))
    if machine is None:
        print("no machine file for {}; known: {}".format(args.machine or socket.gethostname(),
                                                        ", ".join(m.machine for m in machines) or "none"))
    return machine


def cmd_run(args: argparse.Namespace) -> int:
    machine = _machine(args)
    if machine is None:
        return 1
    if args.config:
        configs = [machine.config(args.config)]
        if configs[0] is None:
            print("{} has no configuration {!r}".format(machine.machine, args.config))
            return 1
    else:
        configs = machine.active_configs()
        if not configs:
            print("{} is switched off (\"enabled\": false) or has no enabled configuration; nothing to run".format(
                machine.machine))
            return 0
    routine, registry = Routine.load(args.runners), Registry.load(args.families)
    unknown = [r for r in args.only if routine.source(r) is None]
    if unknown:
        print("unknown runner(s): {}; known: {}".format(", ".join(unknown), ", ".join(s.id for s in routine.sources)))
        return 1
    if args.dry_run:
        print("\n".join(describe(machine, routine, registry, configs, args.only)))
        return 0
    try:
        runs = Night(machine, routine, registry, configs, args.only, tests=args.test, keep_local=args.no_upload,
                     commit=args.commit, note=args.note).run()
    except NightError as error:
        print("error:", error)
        return 3
    for path in runs:
        print("run directory:", path)
    if not args.no_upload:
        print("next: 'amps-nightly pack' sends the results to {}".format(machine.upload or "the server"))
    return 0


def cmd_pack(args: argparse.Namespace) -> int:
    machine = _machine(args)
    if machine is None:
        return 1
    workdir = Path(machine.workdir).expanduser()
    outbox = workdir / "outbox"
    status = 0
    finished, abandoned = transfer.runs_to_pack(workdir)
    busy = transfer.running_night(workdir)
    if busy:
        print("a run is in progress (night {}); its results are packed after it finishes".format(busy))
    if abandoned:
        try:
            routine = Routine.load(args.runners)
        except RoutineError:
            routine = None

        def family_of(runner: str) -> str:
            source = routine.source(runner) if routine else None
            return source.family if source else ""

        for rundir in abandoned:
            try:
                transfer.finalize(rundir, family_of)
                print("completed {}: its orchestrator stopped before the end; the run is marked incomplete".format(rundir))
                finished.append(rundir)
            except (transfer.TransferError, RecordError, OSError) as error:
                print("ERROR:", error)
                status = 1
    for rundir in finished:
        try:
            tarball = transfer.pack(rundir, outbox)
            print("packed {} -> {} ({:.1f} MB)".format(rundir, tarball, tarball.stat().st_size / 2 ** 20))
        except (transfer.TransferError, OSError) as error:
            print("ERROR: packing {} failed: {}".format(rundir, error))
            status = 1
    waiting = transfer.unsent(outbox)
    if args.no_send:
        for tarball in waiting:
            print("not sent (--no-send):", tarball)
    else:
        for tarball in waiting:
            try:
                transfer.send(tarball, machine.upload)
                print("sent {} to {}".format(tarball.name, machine.upload))
            except transfer.TransferError as error:
                print("ERROR: {}; it stays in {} and is sent next time".format(error, outbox))
                status = 1
    transfer.clean_outbox(outbox)
    if not finished and not waiting:
        print("nothing to pack or send: every run under {} has been sent".format(workdir / "nights"))
    return status


def _confirm(yes: bool) -> bool:
    if yes:
        return True
    if not sys.stdin.isatty():
        print("nothing changed: run again with --yes to make these changes")
        return False
    return input("Make these changes? [y/N] ").strip().lower() in ("y", "yes")


def cmd_install(args: argparse.Namespace) -> int:
    machine = _machine(args)
    if machine is None:
        return 1
    bin_dir = Path(args.bin).expanduser()
    old = install.read_conf(bin_dir / install.CONF_NAME)
    server = bool(machine.upload) and not transfer.is_remote(machine.upload)
    settings = dict(old)
    settings["AMPSTEST"] = str(Path(args.root or old.get("AMPSTEST") or machine.server_root or
                                    Path(machine.workdir).expanduser().parent).expanduser())
    settings["TESTDIR"] = str(Path(machine.workdir).expanduser())
    settings["BRANCH"] = args.branch or old.get("BRANCH") or "master"
    settings["UPDATE_SITE"] = args.update_site or old.get("UPDATE_SITE") or ("yes" if server else "no")
    settings["PYTHON"] = args.python or old.get("PYTHON") or "python3"
    start = args.time or machine.start or "20:00"
    try:
        crontab = install.read_crontab()
        plan = install.plan_install(bin_dir, settings, start, crontab)
    except (install.InstallError, ValueError, OSError) as error:
        print("error:", error)
        return 1
    print("amps-nightly install for {} (machine file {}):".format(machine.machine, machine.path))
    for line in plan.lines:
        print("  " + line)
    data = machine.vars.get("swmf_data")
    if data and not Path(data).expanduser().exists():
        print("  note: SWMF_data ({}) does not exist yet; the first night clones it (about 4 GB)".format(data))
    print("  note: each night downloads AMPS {} from GitHub; utility/nightly must be on that branch".format(
        settings["BRANCH"]))
    if not plan.changes:
        print("nothing to change")
        return 0
    if args.dry_run or not _confirm(args.yes):
        return 0 if args.dry_run else 1
    try:
        install.apply(plan, crontab)
    except (install.InstallError, OSError) as error:
        print("error:", error)
        return 1
    print("installed: the nightly starts at {} every day; its log is {}/logs/nightly-<date>.log".format(
        start, settings["AMPSTEST"]))
    return 0


def cmd_uninstall(args: argparse.Namespace) -> int:
    bin_dir = Path(args.bin).expanduser()
    old = install.read_conf(bin_dir / install.CONF_NAME)
    backups = Path(old.get("AMPSTEST") or Path.home()).expanduser() / "backup"
    try:
        crontab = install.read_crontab()
    except install.InstallError as error:
        print("error:", error)
        return 1
    plan = install.plan_uninstall(bin_dir, crontab, backups)
    if not plan.changes:
        print("nothing to remove: no {} in the crontab or in {}".format(install.SCRIPT_NAME, bin_dir))
        return 0
    print("amps-nightly uninstall (the nightly's data in {} stays):".format(old.get("AMPSTEST", "AMPSTEST")))
    for line in plan.lines:
        print("  " + line)
    if args.dry_run or not _confirm(args.yes):
        return 0 if args.dry_run else 1
    try:
        install.apply(plan, crontab)
    except (install.InstallError, OSError) as error:
        print("error:", error)
        return 1
    print("removed")
    return 0


def cmd_relocate(args: argparse.Namespace) -> int:
    root, www = _paths(args)
    with store_mod.Store(root / "results.db") as st:
        moved = store_mod.relocate(www, st)
    print("moved {} run directories into {}".format(moved, www / "archive"))
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    machine = _machine(args)
    if machine is None:
        return 1
    configs = [machine.config(args.config)] if args.config else machine.configs
    if None in configs:
        print("{} has no configuration {!r}".format(machine.machine, args.config))
        return 1
    routine, registry = Routine.load(args.runners), Registry.load(args.families)
    amps_dir = Path(args.amps_dir).expanduser() if args.amps_dir else None
    status = 0
    if not machine.enabled:
        print("{} is switched off (\"enabled\": false): nothing runs there; the plan below applies once it is "
              "switched on".format(machine.machine))
    for config in configs:
        plans = plan(machine, config, routine, registry, amps_dir=amps_dir, table=args.table)
        print("\n{} · {} (column {}){}".format(machine.machine, config.compiler, config.column,
                                               "" if config.enabled else "  [configuration switched off]"))
        for sp in plans:
            if not sp.active:
                print("  {:<14} not run   {}".format(sp.source.id, sp.note))
                continue
            skipped = sp.count(SKIP)
            print("  {:<14} {:>4} run{}   {}".format(sp.source.id, sp.count(RUN),
                                                   ", {} skipped".format(skipped) if skipped else "", sp.note))
            if args.list:
                for item in sp.items:
                    print("      {:<4} {:<20} {}{}".format(item.action.upper(), item.family, item.test,
                                                         "   -- " + item.reason if item.reason else ""))
        runs = sum(sp.count(RUN) for sp in plans)
        skips = sum(sp.count(SKIP) for sp in plans)
        print("  total: {} test{} run, {} reported as skipped".format(runs, "" if runs == 1 else "s", skips))
        for problem in duplicates(plans):
            print("  ERROR:", problem)
            status = 1
    return status


def cmd_import_archive(args: argparse.Namespace) -> int:
    root, _ = _paths(args)
    out = Path(args.out).expanduser() if args.out else root / "incoming" / "legacy-archive"
    out.mkdir(parents=True, exist_ok=True)
    known: List[str] = []
    db = root / "results.db"
    if db.exists():
        with store_mod.Store(db) as st:
            known = st.run_ids()
    commits = legacy.CommitIndex(Path(args.history).expanduser(), args.branch) if args.history else None
    result = legacy.import_archive(Path(args.archive).expanduser(), out, read_table(args.table), _folder_map(args.map),
                                   commits, _date(args.since), _date(args.until), known)
    print("wrote {} run directories to {}".format(len(result["written"]), out))
    print("skipped {} nights that repeat the previous night's files, {} already imported, {} empty folders".format(
        len(result["duplicate"]), len(result["existing"]), len(result["empty"])))
    if not commits:
        print("note: no --history repository given, so commit IDs of archived nights are unknown")
    return 0


def cmd_import_run(args: argparse.Namespace) -> int:
    root, _ = _paths(args)
    out = Path(args.out).expanduser() if args.out else root / "incoming" / "legacy-run"
    amps = Path(args.amps_dir).expanduser()
    table = Path(args.table) if args.table else amps / "MakefileTest" / "Table"
    path = legacy.import_workdir(amps, out, read_table(table), args.machine, args.config)
    print("wrote", path)
    return 0


def cmd_ingest(args: argparse.Namespace) -> int:
    root, www = _paths(args)
    root.mkdir(parents=True, exist_ok=True)
    incoming = root / "incoming"
    for name in getattr(args, "tarballs", []) or []:
        source = Path(name).expanduser()
        if not source.is_file() or not source.name.endswith(transfer.SUFFIX):
            print("ERROR: {} is not a {} file".format(source, transfer.SUFFIX))
            return 1
        incoming.mkdir(parents=True, exist_ok=True)
        if source.resolve().parent != incoming.resolve():
            part = incoming / ("." + source.name + ".part")
            shutil.copyfile(str(source), str(part))
            os.replace(str(part), str(incoming / source.name))
    with store_mod.Store(root / "results.db") as st:
        accepted, rejected = store_mod.ingest(root, www, st)
    print("ingested {} run{}".format(len(accepted), "" if len(accepted) == 1 else "s"))
    for directory, reason in rejected:
        print("REJECTED {}: {}".format(directory, reason))
    return 1 if rejected else 0


def cmd_report(args: argparse.Namespace) -> int:
    root, www = _paths(args)
    db = root / "results.db"
    if not db.exists():
        print("no results store at {}; run 'amps-nightly ingest' first".format(db))
        return 1
    now = parse_time(args.now) if args.now else None
    machines, problems = load_machines(args.machines)
    for problem in problems:
        print("warning: ignoring", problem)
    with store_mod.Store(db) as st:
        summary = report.build(st, Registry.load(args.families), www, now=now,
                               window_days=args.window_days, stale_hours=args.stale_hours,
                               declared=declared_columns(machines) if machines else None)
    print("wrote {} ({} columns, {} families, {} family pages)".format(
        www / "index.html", len(summary["columns"]), len(summary["families"]), summary["family_pages"]))
    print("wrote {}; updated the previous/next links of {} night pages".format(
        www / "archive" / "index.html", summary["night_pages_updated"]))
    if summary["waiting"]:
        print("no results yet for:", ", ".join(summary["waiting"]))
    if summary["hidden"]:
        print("not shown (no enabled machine file):", ", ".join(summary["hidden"]))
    return 0


def cmd_prune(args: argparse.Namespace) -> int:
    _, www = _paths(args)
    compressed, deleted = store_mod.prune(www, args.compress_after_days, args.delete_after_days)
    print("archive: compressed the logs of {} nights, deleted the logs of {} nights".format(compressed, deleted))
    return 0


def cmd_update(args: argparse.Namespace) -> int:
    status = cmd_ingest(args)
    status = cmd_report(args) or status
    cmd_prune(args)
    return status


# ----------------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="amps-nightly", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")
    sub.required = True

    def server(p: argparse.ArgumentParser, www: bool = True) -> None:
        p.add_argument("--root", help="server data directory: results.db, incoming/, www/ (default: "
                                      "$AMPS_NIGHTLY_ROOT, else this host's server_root, else ~/AMPSTEST/nightly)")
        if www:
            p.add_argument("--www", help="web directory for the pages (default: ROOT/www)")

    def config_files(p: argparse.ArgumentParser, table: bool = True) -> None:
        if table:
            p.add_argument("--table", help="Table path (default: MakefileTest/Table of this checkout)")
        p.add_argument("--families", help="families.json path (default: utility/nightly/families.json)")
        p.add_argument("--runners", help="runners.json path (default: utility/nightly/runners.json)")
        p.add_argument("--machines", help="machine files directory (default: utility/nightly/machines)")

    p = sub.add_parser("install", help="install ~/bin/amps_nightly_test.sh, its settings and its crontab line "
                                        "(shows the changes and asks first)")
    p.add_argument("--root", help="AMPSTEST: the nightly's logs, SWMF_data and (on the server) results "
                                  "(default: as installed, else the machine file's server_root, else the "
                                  "parent of its workdir)")
    p.add_argument("--branch", help="AMPS branch to download and test (default: as installed, else master)")
    p.add_argument("--update-site", choices=("yes", "no"),
                   help="also update the results and pages here (default: yes if the machine uploads to a local directory)")
    p.add_argument("--python", help="Python 3.8+ that runs amps-nightly (default: as installed, else python3)")
    p.add_argument("--time", help="start time HH:MM (default: the machine file's start)")
    p.add_argument("--machine", help="machine name (default: the machine file listing this host)")
    p.add_argument("--machines", help="machine files directory (default: utility/nightly/machines)")
    p.add_argument("--bin", default="~/bin", help="where the script goes (default: ~/bin)")
    p.add_argument("--dry-run", action="store_true", help="only show the changes")
    p.add_argument("--yes", action="store_true", help="make the changes without asking")
    p.set_defaults(func=cmd_install)

    p = sub.add_parser("uninstall", help="remove the crontab line and ~/bin/amps_nightly_test.sh (data stays)")
    p.add_argument("--bin", default="~/bin", help="where the script is (default: ~/bin)")
    p.add_argument("--dry-run", action="store_true", help="only show the changes")
    p.add_argument("--yes", action="store_true", help="make the changes without asking")
    p.set_defaults(func=cmd_uninstall)

    p = sub.add_parser("check", help="validate the Table, families.json, runners.json and machines/ together")
    config_files(p)
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("machines", help="list the machine files, their columns and last uploads")
    p.add_argument("--machines", help="machine files directory (default: utility/nightly/machines)")
    p.add_argument("--var", metavar="NAME", help="print only this var of this host's machine file (for scripts)")
    p.add_argument("--machine", help="with --var: this machine instead of the one listing this host")
    server(p, www=False)
    p.set_defaults(func=cmd_machines)

    p = sub.add_parser("plan", help="show what the nightly routine runs and skips on a machine")
    p.add_argument("--machine", help="machine name (default: the machine file listing this host)")
    p.add_argument("--config", help="only this configuration")
    p.add_argument("--amps-dir", help="configured AMPS tree for the exact Makefile.test list "
                                      "(default: this checkout, if it was configured on that machine)")
    p.add_argument("--list", action="store_true", help="list every test")
    config_files(p)
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("run", help="run the nightly test procedure of runners.json on this machine")
    p.add_argument("--machine", help="machine name (default: the machine file listing this host)")
    p.add_argument("--config", help="only this configuration")
    p.add_argument("--only", action="append", default=[], metavar="RUNNER", help="only this runner (repeatable)")
    p.add_argument("--test", action="append", default=[], metavar="GLOB",
                   help="only Makefile.test tests matching this glob (repeatable; for trials)")
    p.add_argument("--no-upload", action="store_true",
                   help="keep the results on this machine: 'amps-nightly pack' will not send them (for trials)")
    p.add_argument("--commit", default="", metavar="SHA",
                   help="test this commit of the branch instead of its head (amps_nightly_test.sh passes the "
                        "commit it downloaded, so the tool and the tested code are the same)")
    p.add_argument("--note", default="", help="text added to the notes of every run (shown on the night's page)")
    p.add_argument("--dry-run", action="store_true", help="print the commands each runner would execute")
    config_files(p, table=False)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("pack", help="tar every finished run with its logs into WORKDIR/outbox and send it to "
                                     "the machine file's upload destination")
    p.add_argument("--machine", help="machine name (default: the machine file listing this host)")
    p.add_argument("--no-send", action="store_true", help="only write the tarballs")
    config_files(p, table=False)
    p.set_defaults(func=cmd_pack)

    p = sub.add_parser("relocate-archive", help="move runs stored under www/runs/ into www/archive/YYYY/MM/DD/<column>/")
    server(p)
    p.set_defaults(func=cmd_relocate)

    p = sub.add_parser("check-table", help="verify Family= lines in MakefileTest/Table")
    p.add_argument("--table", help="Table path (default: MakefileTest/Table of this checkout)")
    p.add_argument("--families", help="families.json path (default: utility/nightly/families.json)")
    p.set_defaults(func=cmd_check_table)

    p = sub.add_parser("import-legacy-archive", help="convert the legacy server archive into run directories")
    p.add_argument("--archive", required=True, help="AMPS_TEST_RESULTS directory (YYYY/MM/DD/<folder>/test_amps.res)")
    p.add_argument("--map", action="append", default=[], metavar="FOLDER=MACHINE/CONFIG",
                   help="legacy folder name to column, e.g. amps-gpu_gnu=passion/gnu (repeatable)")
    p.add_argument("--history", help="git repository with full AMPS history, used to infer each night's commit")
    p.add_argument("--branch", default="master", help="branch the legacy nightly checked out (default: master)")
    p.add_argument("--table", help="Table that provides Family= (default: this checkout's MakefileTest/Table)")
    p.add_argument("--since", help="first archive date to import, YYYY-MM-DD")
    p.add_argument("--until", help="last archive date to import, YYYY-MM-DD")
    p.add_argument("--out", help="output directory (default: ROOT/incoming/legacy-archive)")
    server(p, www=False)
    p.set_defaults(func=cmd_import_archive)

    p = sub.add_parser("import-legacy-run", help="convert the *.diff files of a legacy work tree into a run directory")
    p.add_argument("--amps-dir", required=True, help="legacy work tree, e.g. ~/Tmp_AMPS_test/GNU/AMPS")
    p.add_argument("--machine", required=True)
    p.add_argument("--config", required=True, help="compiler configuration, e.g. gnu")
    p.add_argument("--table", help="Table path (default: AMPS_DIR/MakefileTest/Table)")
    p.add_argument("--out", help="output directory (default: ROOT/incoming/legacy-run)")
    server(p, www=False)
    p.set_defaults(func=cmd_import_run)

    p = sub.add_parser("ingest", help="move complete uploads from ROOT/incoming into the results store")
    server(p)
    p.add_argument("tarballs", nargs="*", metavar="TARBALL", help="run tarballs to add (copied into ROOT/incoming first)")
    p.set_defaults(func=cmd_ingest)

    for name, func, text in (("report", cmd_report, "regenerate the report pages"),
                             ("update", cmd_update, "ingest, then regenerate the report pages")):
        p = sub.add_parser(name, help=text)
        server(p)
        p.add_argument("--families", help="families.json path (default: utility/nightly/families.json)")
        p.add_argument("--machines", help="machine files directory; its enabled configurations are the columns "
                                          "(default: utility/nightly/machines)")
        p.add_argument("--window-days", type=int, default=60,
                       help="without machine files: hide columns without a run in this many days")
        p.add_argument("--stale-hours", type=float, default=36, help="flag columns whose newest run is older")
        p.add_argument("--now", help=argparse.SUPPRESS)
        if name == "update":
            p.add_argument("tarballs", nargs="*", metavar="TARBALL",
                           help="run tarballs to add (default: what 'amps-nightly pack' sent to ROOT/incoming)")
        p.set_defaults(func=func, compress_after_days=14, delete_after_days=400)

    p = sub.add_parser("prune", help="gzip archived logs after 14 nights, delete them after 400 days")
    server(p)
    p.add_argument("--compress-after-days", type=int, default=14)
    p.add_argument("--delete-after-days", type=int, default=400)
    p.set_defaults(func=cmd_prune)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (RegistryError, RoutineError) as error:
        print("error:", error, file=sys.stderr)
        return 2
