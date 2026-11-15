"""Convert results of the legacy nightly scripts (utility/TestScripts) into run directories.

Legacy semantics: every test writes test_<name>[<ref>].diff. Each make stage appends a marker
("test_<name>_compile...", "test_<name>_compile... done", rundir, run, check) and the check
stage overwrites the file with DiffNum.pl output, so an EMPTY file means the test passed.
Tests may also write extra result files (test_<name>--DUST.diff, test_<name>1.diff, ...); those
are folded into the owning test, which fails if any of its files is non-empty.

Two sources are supported:
  * the server archive  AMPS_TEST_RESULTS/YYYY/MM/DD/<folder>/test_amps.res  (ls -ltr + head -100)
  * a legacy work tree  Tmp_AMPS_test/<COMPILER>/AMPS/*.diff  right after a nightly run
"""

from __future__ import annotations

import bisect
import collections
import datetime as dt
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .records import make_run_id, utc_text, write_run_dir
from .registry import UNASSIGNED
from .table import TableTest

_MARK = re.compile(r"^test_(?P<name>.+)_(?P<stage>compile|rundir|run|check)\.\.\.(?P<done> done)?$")
_NEXT = {"compile": "rundir", "rundir": "run", "run": "check", "check": ""}
_INCOMPLETE = {
    "compile": "compilation failed",
    "rundir": "run directory setup failed",
    "run": "run did not finish (crash, abort or time limit)",
    "check": "result check did not finish",
}
_NEVER = {
    "rundir": "run directory was never prepared",
    "run": "never executed: the run stage did not start",
    "check": "result check never ran",
}
_MISSING = re.compile(r"(?:ERROR in DiffNum\.pl:\s*)?(\S+?) does not exist|could not open\s+(\S+)", re.I)
_CHECKOUT = re.compile(r"AMPS was checked out on (.+)")
_EXTRA_SUFFIX = re.compile(r"^(--|_|[0-9]|-[0-9]|-[^-]*=)")
_TZ_HOURS = {"UTC": 0, "GMT": 0, "EST": -5, "EDT": -4, "CST": -6, "CDT": -5,
             "MST": -7, "MDT": -6, "PST": -8, "PDT": -7}


# ----------------------------------------------------------------------------- one .diff file
def split_result_file(filename: str) -> Tuple[str, str]:
    """'test_Saturn-MOP[5].diff' -> ('Saturn-MOP', '5')."""
    base = filename.rsplit("/", 1)[-1]
    if base.startswith("test_"):
        base = base[len("test_"):]
    if base.endswith(".diff"):
        base = base[: -len(".diff")]
    base = re.sub(r"\{[^}]*\}$", "", base)       # an unexpanded {TEST<name>-REF} left in some old file names
    match = re.match(r"^(.*)\[([^\]]*)\]$", base)
    return (match.group(1), match.group(2)) if match else (base, "")


def _short(text: str, limit: int = 160) -> str:
    text = " ".join(text.split())
    return text if len(text) <= limit else text[: limit - 1] + "…"


def evaluate_diff(size: int, text: str) -> Tuple[str, str, str]:
    """(status, stage, reason) of one legacy result file of the given byte size."""
    if size == 0:
        return "PASS", "check", ""
    marks, content = [], []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        match = _MARK.match(line)
        (marks if match else content).append(match or line)
    if not content:
        if not marks:
            return "FAILED", "check", "result file is not empty but holds only white space"
        last = marks[-1]
        stage = last.group("stage")
        if last.group("done"):
            following = _NEXT[stage]
            if not following:
                return "FAILED", "check", "result check finished without a verdict"
            return "FAILED", following, _NEVER[following]
        return "FAILED", stage, _INCOMPLETE[stage]
    last = marks[-1] if marks else None
    if last is not None and not last.group("done") and last.group("stage") != "check":
        # the stage printed something (banners, messages) but never reported completion
        return "FAILED", last.group("stage"), _INCOMPLETE[last.group("stage")]
    stage = "check"
    if last is not None:
        stage = last.group("stage")
        if last.group("done") and _NEXT[stage]:
            stage = _NEXT[stage]
    informative = [line for line in content if re.search(r"[A-Za-z0-9]", line)]
    if not informative:
        return "FAILED", stage, "result file holds only separator lines"
    for line in informative:
        missing = _MISSING.search(line)
        if missing:
            target = (missing.group(1) or missing.group(2) or "").rsplit("/", 1)[-1]
            return "FAILED", stage, "output not produced: " + _short(target, 120)
    if informative[0].startswith("DiffNum.pl"):
        reference = informative[0].split()[-1].rsplit("/", 1)[-1]
        differing = sum(1 for line in content if line.startswith("<"))
        truncated = text.count("\n") >= 99          # archive and importer keep only the first 100 lines
        if differing:
            what = "{}{} line{} differ".format("at least " if truncated else "", differing, "" if differing == 1 else "s")
        else:
            what = "values differ"
        return "FAILED", stage, "{} from reference {}".format(what, _short(reference, 100))
    return "FAILED", stage, _short(informative[0])


# ----------------------------------------------------------------------------- many files -> results
def collect_results(files: "collections.OrderedDict[str, Tuple[int, str]]",
                    table: Sequence[TableTest]) -> List[Dict[str, object]]:
    """Fold legacy result files (file name -> (size, text)) into one result per test."""
    family = {t.name: t.family for t in table if t.family}
    owner_of_extra: Dict[str, str] = {}
    for test in table:
        for extra in test.result_files:
            owner_of_extra.setdefault(extra, test.name)
    known = sorted(family, key=len, reverse=True)

    groups: "collections.OrderedDict[str, Dict[str, list]]" = collections.OrderedDict()
    for filename, (size, text) in files.items():
        name, _ = split_result_file(filename)
        if name in family:
            owner, primary = name, True
        elif name in owner_of_extra:
            owner, primary = owner_of_extra[name], False
        else:
            # test_<owner>--DUST, test_<owner>1, test_<owner>_x, test_<owner>-2, test_<owner>-test=0-wave;
            # a plain "-word" suffix is not accepted: it usually names a different (renamed or removed) test
            owner = next((k for k in known if name.startswith(k) and _EXTRA_SUFFIX.match(name[len(k):])), name)
            primary = owner == name
        group = groups.setdefault(owner, {"primary": [], "extra": []})
        group["primary" if primary else "extra"].append((filename, size, text))

    results = []
    for owner, group in groups.items():
        if group["primary"]:
            filename, size, text = group["primary"][0]
            status, stage, reason = evaluate_diff(size, text)
        else:
            status, stage, reason = "PASS", "check", ""
        extras = [(f, s, t) for f, s, t in group["extra"] if s]
        if extras and status == "PASS":
            f, s, t = extras[0]
            _, stage, why = evaluate_diff(s, t)
            status = "FAILED"
            reason = "{} of {} extra result files differ ({}: {})".format(
                len(extras), len(group["extra"]), split_result_file(f)[0], why)
        # The .diff contents are only read to decide status, stage and a one-line reason;
        # they are deliberately not stored.
        results.append({
            "family": family.get(owner, UNASSIGNED),
            "test": owner,
            "status": status,
            "stage": stage if status != "PASS" else "",
            "reason": reason,
            "source": "makefile",
        })
    return results


def parse_res(text: str) -> "collections.OrderedDict[str, Tuple[int, str]]":
    """Parse test_amps.res: 'ls -ltr *diff' lines, a line of '=', then 'head -100 *diff'."""
    parts = re.split(r"^=+[ \t]*$", text, maxsplit=1, flags=re.M)
    listing_text, heads = parts[0], (parts[1] if len(parts) > 1 else "")
    if heads.startswith("\n"):          # rest of the separator line
        heads = heads[1:]
    sizes: "collections.OrderedDict[str, int]" = collections.OrderedDict()
    for line in listing_text.splitlines():
        parts = line.split()
        if len(parts) >= 9 and parts[-1].endswith(".diff") and parts[4].isdigit():
            sizes[parts[-1]] = int(parts[4])
    bodies: Dict[str, List[str]] = {}
    current: Optional[str] = None
    for line in heads.splitlines():
        header = re.match(r"^==> (.+) <==$", line)
        if header:
            current = header.group(1)
            bodies[current] = []
        elif current is not None:
            bodies[current].append(line)
    if not bodies and len(sizes) == 1:
        bodies[next(iter(sizes))] = heads.splitlines()
    return collections.OrderedDict((name, (size, "\n".join(bodies.get(name, []))))
                                   for name, size in sizes.items())


# ----------------------------------------------------------------------------- times and commits
def parse_legacy_date(text: str) -> Optional[dt.datetime]:
    """'Thu 22 Oct 2026 10:30:01 PM CDT' or 'Thu Oct 22 22:30:01 CDT 2026' -> aware datetime."""
    words = text.split()
    zone = next((w for w in words if w in _TZ_HOURS), None)
    if zone is None:
        return None
    rest = " ".join(w for w in words if w != zone)
    for fmt in ("%a %d %b %Y %I:%M:%S %p", "%a %b %d %H:%M:%S %Y", "%a %d %b %Y %H:%M:%S"):
        try:
            naive = dt.datetime.strptime(rest, fmt)
        except ValueError:
            continue
        return naive.replace(tzinfo=dt.timezone(dt.timedelta(hours=_TZ_HOURS[zone])))
    return None


def legacy_checkout_time(log: Path) -> Optional[dt.datetime]:
    try:
        with open(str(log), "rb") as stream:
            head = stream.read(4096).decode("utf-8", "replace")
    except OSError:
        return None
    match = _CHECKOUT.search(head)
    return parse_legacy_date(match.group(1).strip()) if match else None


class CommitIndex:
    """Commit that was the tip of a branch at a given time (first-parent history, committer dates)."""

    def __init__(self, git_dir: Path, ref: str = "master"):
        out = subprocess.run(["git", "-C", str(git_dir), "log", "--first-parent", "--format=%ct %H", ref],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True).stdout
        pairs = sorted((int(t), h) for t, h in (line.split() for line in out.splitlines() if line.strip()))
        self._times = [t for t, _ in pairs]
        self._shas = [h for _, h in pairs]

    def at(self, when: dt.datetime) -> str:
        index = bisect.bisect_right(self._times, int(when.timestamp())) - 1
        return self._shas[index] if index >= 0 else ""


def _git_head(path: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""
    return out if re.match(r"^[0-9a-f]{40}$", out) else ""


def _mtime(path: Path) -> dt.datetime:
    return dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc)


# ----------------------------------------------------------------------------- importers
def _run_record(machine: str, config: str, finished: dt.datetime, started: Optional[dt.datetime],
                commit: str, inferred: bool, producer: str, toolchain: str, notes: str) -> Dict[str, object]:
    return {
        "schema_version": 1,
        "run_id": make_run_id(machine, config, finished),
        "machine": machine, "config": config, "compiler": config.upper(),
        "toolchain": toolchain, "branch": "master",
        "commit": commit, "commit_inferred": inferred,
        "started": utc_text(started) if started else "",
        "finished": utc_text(finished), "complete": True,
        "producer": producer, "notes": notes,
    }


def import_archive(archive: Path, out_dir: Path, table: Sequence[TableTest],
                   folders: Dict[str, Tuple[str, str]], commits: Optional[CommitIndex] = None,
                   since: Optional[dt.date] = None, until: Optional[dt.date] = None,
                   skip_run_ids: Iterable[str] = ()) -> Dict[str, List[str]]:
    """Convert AMPS_TEST_RESULTS/YYYY/MM/DD/<folder>/test_amps.res into run directories under out_dir.

    folders maps legacy folder names to (machine, config); other folders are split at the last '_'.
    Nights whose checkout time equals the previous imported night of the same column are copies
    of stale files (the legacy scp kept sending the last results) and are skipped.
    """
    skip = set(skip_run_ids)
    report: Dict[str, List[str]] = {"written": [], "duplicate": [], "existing": [], "empty": []}
    last_checkout: Dict[Tuple[str, str], dt.datetime] = {}
    days = sorted(p for p in Path(archive).glob("[0-9][0-9][0-9][0-9]/[0-9][0-9]/[0-9][0-9]") if p.is_dir())
    for day in days:
        date = dt.date(int(day.parent.parent.name), int(day.parent.name), int(day.name))
        if (since and date < since) or (until and date > until):
            continue
        for folder in sorted(p for p in day.iterdir() if p.is_dir()):
            res = folder / "test_amps.res"
            label = "{}/{}".format(day.relative_to(archive), folder.name)
            if not res.is_file() or res.stat().st_size == 0:
                if folder.exists() and any(folder.iterdir()):
                    report["empty"].append(label)
                continue
            machine, config = folders.get(folder.name) or _split_folder(folder.name)
            checkout = legacy_checkout_time(folder / "test_amps.log")
            key = (machine, config)
            if checkout is not None and last_checkout.get(key) == checkout:
                report["duplicate"].append(label)
                continue
            if checkout is not None:
                last_checkout[key] = checkout
            finished = _mtime(res)
            if abs((finished.date() - date).days) > 2:      # file touched later: fall back to the folder date
                finished = dt.datetime(date.year, date.month, date.day, 14, 0, tzinfo=dt.timezone.utc)
            commit = commits.at(checkout) if (commits and checkout) else ""
            run = _run_record(machine, config, finished, checkout, commit, bool(commit), "import-legacy-archive",
                              "legacy nightly scripts (utility/TestScripts)",
                              "imported from {}".format(label))
            if run["run_id"] in skip or (Path(out_dir) / str(run["run_id"])).exists():
                report["existing"].append(label)
                continue
            results = collect_results(parse_res(res.read_text(errors="replace")), table)
            write_run_dir(Path(out_dir) / str(run["run_id"]), run, results)
            report["written"].append(str(run["run_id"]))
    return report


def read_diff_files(directory: Path) -> "collections.OrderedDict[str, Tuple[int, str]]":
    """The *.diff result files of a legacy tree: name -> (size, first 100 lines), as the archive kept them."""
    files: "collections.OrderedDict[str, Tuple[int, str]]" = collections.OrderedDict()
    for path in sorted(Path(directory).glob("*.diff")):
        with open(str(path), errors="replace") as stream:
            files[path.name] = (path.stat().st_size, "".join(line for _, line in zip(range(100), stream)))
    return files


def import_workdir(amps_dir: Path, out_dir: Path, table: Sequence[TableTest], machine: str, config: str,
                   toolchain: str = "legacy nightly scripts (utility/TestScripts)") -> Path:
    """Convert the *.diff files left in a legacy work tree (e.g. ~/Tmp_AMPS_test/GNU/AMPS) into a run directory."""
    amps_dir = Path(amps_dir)
    files = read_diff_files(amps_dir)
    newest = max((_mtime(amps_dir / name) for name in files), default=None)
    if not files:
        raise ValueError("no *.diff files in {}".format(amps_dir))
    commit = _git_head(amps_dir)
    run = _run_record(machine, config, newest, legacy_checkout_time(amps_dir / "test_amps.log"), commit, False,
                      "import-legacy-run", toolchain, "imported from {}".format(amps_dir))
    return write_run_dir(Path(out_dir) / str(run["run_id"]), run, collect_results(files, table))


def _split_folder(name: str) -> Tuple[str, str]:
    machine, _, config = name.rpartition("_")
    return (machine.lower(), config.lower()) if machine and config else (name.lower(), "default")
