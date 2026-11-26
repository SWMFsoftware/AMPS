"""`amps-nightly install` / `uninstall`: the nightly cron entry of a test machine
(the successors of `Config.pl -set-test` / `-rm-test` and utility/TestScripts/SetNightlyTest.pl).

install
    copies utility/nightly/bin/amps_nightly_test.sh to ~/bin (unchanged, so it can be compared with the
    repository's version each night), writes its settings to ~/bin/amps_nightly_test.conf and adds one
    crontab line per routine the machine runs, below a marker comment: the nightly's at the machine file's
    "start" every day, and one for each routine enabled in its "routines" at that routine's "start"
    (e.g. "Sat 10:00": weekly), which passes --routine <name> to the script. Run again, it updates the copy and
    keeps the settings.
uninstall
    removes the crontab lines and the two files; the nightly's data (AMPSTEST) stays.

Both show every change first and need a yes (or --yes); the crontab is saved before it is changed.
"""

from __future__ import annotations

import datetime as dt
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .machines import parse_start

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "amps_nightly_test.sh"
SCRIPT_NAME = "amps_nightly_test.sh"
CONF_NAME = "amps_nightly_test.conf"
MARKER = "# amps-nightly: AMPS nightly testing (added by 'amps-nightly install'; remove with 'amps-nightly uninstall')"
SETTINGS = ("AMPSTEST", "TESTDIR", "BRANCH", "UPDATE_SITE", "PYTHON", "REPOSITORY", "SWMF_DATA_REPOSITORY",
            "SWMF_DATA_BRANCH", "KEEP_LOG_DAYS", "LOCK_WAIT_H")
CONF_HEADER = """# Settings of ~/bin/amps_nightly_test.sh (written by 'amps-nightly install'; see the script's header).
# AMPSTEST: the nightly's logs and SWMF_data (and on the server: results, archive, web pages);
# TESTDIR: where the tests run (the machine file's workdir); BRANCH: the AMPS branch downloaded and tested;
# UPDATE_SITE: yes on the results server; PYTHON: the Python 3.8+ that runs amps-nightly.
"""


class InstallError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- settings file
def read_conf(path: Path) -> Dict[str, str]:
    """KEY=VALUE lines of a settings file (shell syntax, as the script sources it)."""
    values: Dict[str, str] = {}
    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        return values
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, raw = line.partition("=")
        try:
            words = shlex.split(raw, comments=True)
        except ValueError:
            continue
        if key.strip() in SETTINGS:
            values[key.strip()] = words[0] if words else ""
    return values


def conf_text(settings: Dict[str, str]) -> str:
    return CONF_HEADER + "".join("{}={}\n".format(k, shlex.quote(settings[k])) for k in SETTINGS if k in settings)


# ----------------------------------------------------------------------------- crontab
def read_crontab() -> str:
    done = subprocess.run(["crontab", "-l"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if done.returncode != 0:
        if "no crontab" in done.stderr:
            return ""
        raise InstallError("crontab -l failed: " + done.stderr.strip())
    return done.stdout


def write_crontab(text: str) -> None:
    done = subprocess.run(["crontab", "-"], input=text, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    if done.returncode != 0:
        raise InstallError("crontab - failed: " + done.stderr.strip())


def cron_line(start: str, script: Path, cron_log: Path, routine: str = "nightly") -> str:
    """The crontab line that starts the script at `start`: "HH:MM" every day, "<Day> HH:MM" on that day of each week;
    a routine other than the nightly is passed to the script (--routine <name>)."""
    day, hour, minute = parse_start(start)
    return "{:02d} {:02d} * * {}  {}{} >> {} 2>&1".format(
        minute, hour, "*" if day is None else day, script, "" if routine == "nightly" else " --routine " + routine,
        cron_log)


def without_entry(text: str) -> Tuple[str, List[str]]:
    """The crontab without the marker line and the active lines that start the script."""
    kept, removed = [], []
    for line in text.splitlines():
        ours = line.strip() == MARKER or (SCRIPT_NAME in line and not line.lstrip().startswith("#"))
        (removed if ours else kept).append(line)
    return "".join(l + "\n" for l in kept), removed


# ----------------------------------------------------------------------------- plans
@dataclass
class Plan:
    lines: List[str] = field(default_factory=list)          # what will be done, for the person to read
    files: Dict[Path, Optional[str]] = field(default_factory=dict)    # path -> new text (None = delete)
    copy: Optional[Tuple[Path, Path]] = None                # script source -> target
    crontab: Optional[str] = None                           # new crontab text
    backup: Optional[Path] = None                           # where the old crontab is saved
    changes: bool = False


def plan_install(bin_dir: Path, settings: Dict[str, str], start: str, crontab: str,
                 source: Path = SCRIPT, now: Optional[dt.datetime] = None,
                 routines: Sequence[Tuple[str, str]] = ()) -> Plan:
    """start: the nightly's start time; routines: (routine, start) of the other routines the machine runs."""
    now = now or dt.datetime.now()
    bin_dir = Path(bin_dir).expanduser()
    target, conf = bin_dir / SCRIPT_NAME, bin_dir / CONF_NAME
    root = Path(settings["AMPSTEST"]).expanduser()
    plan = Plan()
    if not source.is_file():
        raise InstallError("{} not found".format(source))
    if target.is_file() and target.read_bytes() == source.read_bytes():
        plan.lines.append("keep     {} (same as {})".format(target, source))
    else:
        plan.lines.append("{}  {}  (copy of {})".format(
            "replace " if target.exists() else "write   ", target, source))
        plan.copy = (source, target)
    text = conf_text(settings)
    old_conf = conf.read_text() if conf.is_file() else None
    if old_conf == text:
        plan.lines.append("keep     {}".format(conf))
    else:
        plan.lines.append("{}  {}:".format("replace " if old_conf is not None else "write   ", conf))
        plan.lines += ["             " + l for l in text.splitlines() if not l.startswith("#")]
        plan.files[conf] = text
    lines = [cron_line(start, target, root / "logs" / "cron.log")] + \
        [cron_line(when, target, root / "logs" / "cron.log", name) for name, when in routines]
    rest, removed = without_entry(crontab)
    new = rest + MARKER + "\n" + "".join(line + "\n" for line in lines)
    if removed == [MARKER] + lines:
        for line in lines:
            plan.lines.append("keep     the crontab entry: " + line)
    else:
        for old in removed:
            plan.lines.append("crontab  remove: " + old)
        plan.lines.append("crontab  add:    " + MARKER)
        for line in lines:
            plan.lines.append("                 " + line)
        plan.backup = root / "backup" / "crontab-{}.txt".format(now.strftime("%Y%m%d-%H%M%S"))
        plan.lines.append("         (the current crontab is saved to {} first)".format(plan.backup))
        plan.crontab = new
    plan.changes = bool(plan.copy or plan.files or plan.crontab is not None)
    return plan


def plan_uninstall(bin_dir: Path, crontab: str, backup_dir: Path, now: Optional[dt.datetime] = None) -> Plan:
    now = now or dt.datetime.now()
    bin_dir = Path(bin_dir).expanduser()
    plan = Plan()
    rest, removed = without_entry(crontab)
    if removed:
        for old in removed:
            plan.lines.append("crontab  remove: " + old)
        plan.backup = Path(backup_dir) / "crontab-{}.txt".format(now.strftime("%Y%m%d-%H%M%S"))
        plan.lines.append("         (the current crontab is saved to {} first)".format(plan.backup))
        plan.crontab = rest
    for path in (bin_dir / SCRIPT_NAME, bin_dir / CONF_NAME):
        if path.exists():
            plan.lines.append("delete   {}".format(path))
            plan.files[path] = None
    plan.changes = bool(plan.files or plan.crontab is not None)
    return plan


def apply(plan: Plan, old_crontab: str, writer: Callable[[str], None] = write_crontab) -> None:
    if plan.copy:
        source, target = plan.copy
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name("." + target.name + ".new")
        shutil.copyfile(str(source), str(part))
        part.chmod(0o755)
        os.replace(str(part), str(target))           # a running copy keeps its old file
    for path, text in plan.files.items():
        if text is None:
            path.unlink()
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    if plan.crontab is not None:
        if plan.backup:
            plan.backup.parent.mkdir(parents=True, exist_ok=True)
            plan.backup.write_text(old_crontab)
        writer(plan.crontab)
