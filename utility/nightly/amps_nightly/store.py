"""Results store (SQLite) and ingestion of uploaded run directories on the server.

Server layout (root: the server_root of the server's machine file, /data/vtenishe/AMPSTEST on passion):

    <root>/incoming/<run-id>.tar.gz           uploads from `amps-nightly pack` (hidden names are still arriving)
    <root>/incoming/[<anything>/]<run-id>/   run directories; only directories containing DONE are taken
    <root>/rejected/<name>-<time>/           uploads that failed validation, with REJECTED.txt
    <root>/results.db                        every run and every test result, kept forever
    <www>/archive/YYYY/MM/DD/<column>/        accepted runs: results, logs and the night summary index.html
                                              (a run of a routine other than the nightly: <column>-<routine>/)
    <www>/index.html, <www>/family/...        report pages (see report.py)

Every run belongs to a routine: "nightly" (runners.json), or another (runners.<routine>.json, e.g. "long").
A results.db made before routines existed gets the column, and its runs are all nightly.
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .records import DONE_FILE, NIGHTLY, RecordError, parse_time, read_run_dir
from .transfer import SUFFIX, TransferError, unpack

_SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id          TEXT PRIMARY KEY,
    column_id       TEXT NOT NULL,
    machine         TEXT NOT NULL,
    config          TEXT NOT NULL,
    compiler        TEXT NOT NULL DEFAULT '',
    toolchain       TEXT NOT NULL DEFAULT '',
    branch          TEXT NOT NULL DEFAULT '',
    commit_sha      TEXT NOT NULL DEFAULT '',
    commit_inferred INTEGER NOT NULL DEFAULT 0,
    started         TEXT NOT NULL DEFAULT '',
    finished        TEXT NOT NULL,
    complete        INTEGER NOT NULL,
    producer        TEXT NOT NULL DEFAULT '',
    notes           TEXT NOT NULL DEFAULT '',
    path            TEXT NOT NULL DEFAULT '',
    ingested        TEXT NOT NULL,
    routine         TEXT NOT NULL DEFAULT 'nightly'
);
CREATE TABLE IF NOT EXISTS results (
    run_id       TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    family       TEXT NOT NULL,
    test         TEXT NOT NULL,
    status       TEXT NOT NULL CHECK (status IN ('PASS', 'FAILED', 'SKIPPED')),
    stage        TEXT NOT NULL DEFAULT '',
    reason       TEXT NOT NULL DEFAULT '',
    source       TEXT NOT NULL DEFAULT '',
    duration_s   REAL,
    cases_total  INTEGER,
    cases_failed INTEGER,
    log          TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (run_id, family, test)
);
CREATE INDEX IF NOT EXISTS runs_by_column ON runs (column_id, finished);
CREATE INDEX IF NOT EXISTS results_by_status ON results (status, family, test);
"""


class Store:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.db = sqlite3.connect(str(self.path))
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys = ON")
        self.db.executescript(_SCHEMA)
        if "routine" not in {row[1] for row in self.db.execute("PRAGMA table_info(runs)")}:
            # a store from before the routines: every run in it is a nightly run
            with self.db:
                self.db.execute("ALTER TABLE runs ADD COLUMN routine TEXT NOT NULL DEFAULT 'nightly'")
        self.db.execute("CREATE INDEX IF NOT EXISTS runs_by_routine ON runs (column_id, routine, finished)")
        self.db.execute("PRAGMA user_version = 2")

    def close(self) -> None:
        self.db.close()

    def __enter__(self) -> "Store":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    # ------------------------------------------------------------------ writes
    def has_run(self, run_id: str) -> bool:
        return self.db.execute("SELECT 1 FROM runs WHERE run_id = ?", (run_id,)).fetchone() is not None

    def run_ids(self) -> List[str]:
        return [r[0] for r in self.db.execute("SELECT run_id FROM runs")]

    def add_run(self, run: Dict[str, Any], results: List[Dict[str, Any]], path: str = "",
                now: Optional[dt.datetime] = None) -> None:
        ingested = (now or dt.datetime.now(dt.timezone.utc)).astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        with self.db:
            self.db.execute(
                "INSERT INTO runs (run_id, column_id, machine, config, compiler, toolchain, branch, commit_sha,"
                " commit_inferred, started, finished, complete, producer, notes, path, ingested, routine)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (run["run_id"], run["column"], run["machine"], run["config"], run["compiler"], run["toolchain"],
                 run["branch"], run["commit"], int(run["commit_inferred"]), run["started"], run["finished"],
                 int(run["complete"]), run["producer"], run["notes"], path, ingested, run.get("routine", NIGHTLY)))
            self.db.executemany(
                "INSERT INTO results (run_id, family, test, status, stage, reason, source, duration_s,"
                " cases_total, cases_failed, log) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [(run["run_id"], r["family"], r["test"], r["status"], r["stage"], r["reason"], r["source"],
                  r["duration_s"], r["cases_total"], r["cases_failed"], r["log"]) for r in results])

    # ------------------------------------------------------------------ reads
    def columns(self, routine: str = NIGHTLY) -> List[Dict[str, Any]]:
        """One row per column with a run of this routine, with its newest run time."""
        rows = self.db.execute(
            "SELECT column_id, machine, config, MAX(finished) AS last_finished, COUNT(*) AS runs"
            " FROM runs WHERE routine = ? GROUP BY column_id ORDER BY machine, config", (routine,)).fetchall()
        return [dict(r) for r in rows]

    def routines(self) -> List[str]:
        """The routines other than the nightly that have runs, by name."""
        return [r[0] for r in self.db.execute(
            "SELECT DISTINCT routine FROM runs WHERE routine != ? ORDER BY routine", (NIGHTLY,))]

    def latest_runs(self, column: str, count: int = 2, routine: str = NIGHTLY) -> List[Dict[str, Any]]:
        rows = self.db.execute("SELECT * FROM runs WHERE column_id = ? AND routine = ? ORDER BY finished DESC LIMIT ?",
                               (column, routine, count)).fetchall()
        return [dict(r) for r in rows]

    def latest_results(self, column: str, routine: str, since: str = "") -> Dict[Tuple[str, str], Dict[str, Any]]:
        """(family, test) -> its newest result in this column among the runs of `routine` that finished at `since`
        (UTC text) or later, with that run under "run": run_id, finished, commit_sha, commit_inferred, path.

        A periodic routine need not run every test each time (a run of a few tests sent by hand), so its report
        shows each test's own newest result."""
        rows = self.db.execute(
            "SELECT r.*, u.finished, u.commit_sha, u.commit_inferred, u.path FROM results r"
            " JOIN runs u ON u.run_id = r.run_id WHERE u.column_id = ? AND u.routine = ? AND u.finished >= ?"
            " ORDER BY u.finished", (column, routine, since)).fetchall()
        latest: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for row in rows:
            row = dict(row)
            row["run"] = {k: row.pop(k) for k in ("finished", "commit_sha", "commit_inferred", "path")}
            row["run"]["run_id"] = row["run_id"]
            latest[(row["family"], row["test"])] = row
        return latest

    def results(self, run_id: str) -> List[Dict[str, Any]]:
        rows = self.db.execute("SELECT * FROM results WHERE run_id = ? ORDER BY family, test", (run_id,)).fetchall()
        return [dict(r) for r in rows]

    def all_runs(self) -> List[Dict[str, Any]]:
        """Every stored run, grouped by column (machine, config) and oldest first within a column."""
        rows = self.db.execute("SELECT * FROM runs ORDER BY machine, config, column_id, finished, run_id").fetchall()
        return [dict(r) for r in rows]

    def status_counts(self) -> Dict[str, Dict[str, int]]:
        """run_id -> {status: number of result rows}, as on each night's summary page."""
        counts: Dict[str, Dict[str, int]] = {}
        for r in self.db.execute("SELECT run_id, status, COUNT(*) AS n FROM results GROUP BY run_id, status"):
            counts.setdefault(r["run_id"], {})[r["status"]] = r["n"]
        return counts

    def last_passes(self, column: str, routine: str = NIGHTLY) -> Dict[Tuple[str, str], Dict[str, Any]]:
        """(family, test) -> newest PASS in this column and routine: run_id, finished, commit_sha, commit_inferred."""
        rows = self.db.execute(
            "SELECT r.family, r.test, u.run_id, u.finished, u.commit_sha, u.commit_inferred"
            " FROM results r JOIN runs u ON u.run_id = r.run_id"
            " WHERE u.column_id = ? AND u.routine = ? AND r.status = 'PASS' ORDER BY u.finished",
            (column, routine)).fetchall()
        return {(r["family"], r["test"]): dict(r) for r in rows}


# ----------------------------------------------------------------------------- archive
def night_of(run: Dict[str, Any]) -> str:
    """Local date on which the run's night started: run.json "night", else the start (or end) in local time."""
    if run.get("night"):
        return run["night"]
    return parse_time(run["started"] or run["finished"]).astimezone().strftime("%Y-%m-%d")


def archive_path(www: Path, run: Dict[str, Any]) -> Path:
    """archive/YYYY/MM/DD/<column>, like the legacy AMPS_TEST_RESULTS/YYYY/MM/DD/<folder>; a second run of
    the same column on the same night gets <column>.2, and so on. A run of another routine than the nightly
    goes to <column>-<routine> (archive/2026/11/28/passion-gnu-long)."""
    base = Path("archive") / night_of(run).replace("-", "/")
    routine = run.get("routine") or NIGHTLY
    folder = run["column"] if routine == NIGHTLY else "{}-{}".format(run["column"], routine)
    relative, number = base / folder, 2
    while (Path(www) / relative).exists():
        relative, number = base / "{}.{}".format(folder, number), number + 1
    return relative


def write_index(directory: Path) -> None:
    """The night's summary page in its archive folder (imported lazily: report imports this module)."""
    from .report import write_run_index
    try:
        write_run_index(Path(directory))
    except (OSError, ValueError, RecordError):
        pass


def relocate(www: Path, store: "Store") -> int:
    """One-time move of runs stored under www/runs/YYYY/MM/DD/<run-id> into the archive layout."""
    moved = 0
    for row in store.db.execute("SELECT run_id, path, column_id, started, finished FROM runs ORDER BY finished").fetchall():
        if not row["path"].startswith("runs/"):
            continue
        source = Path(www) / row["path"]
        if not source.is_dir():
            continue
        run = {"column": row["column_id"], "started": row["started"], "finished": row["finished"]}
        try:
            night = json.loads((source / "run.json").read_text()).get("night")
            if night:
                run["night"] = night
        except (OSError, ValueError):
            pass
        relative = archive_path(www, run)
        (Path(www) / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(Path(www) / relative))
        with store.db:
            store.db.execute("UPDATE runs SET path = ? WHERE run_id = ?", (str(relative), row["run_id"]))
        write_index(Path(www) / relative)
        moved += 1
    old = Path(www) / "runs"
    if old.is_dir() and not any(p.is_file() for p in old.rglob("*")):
        shutil.rmtree(str(old))
    return moved


def prune(www: Path, compress_after_days: int = 14, delete_after_days: int = 400,
          today: Optional[dt.date] = None) -> Tuple[int, int]:
    """Gzip the logs of nights older than compress_after_days (browsers only download .gz, so recent logs stay
    plain), delete the logs of nights older than delete_after_days; results and summaries are kept.
    Returns (nights compressed, nights whose logs were deleted)."""
    import gzip
    today = today or dt.date.today()
    compressed = deleted = 0
    for night in sorted(Path(www).glob("archive/*/*/*/*")):
        logs = night / "logs"
        if not logs.is_dir():
            continue
        try:
            date = dt.date(*(int(x) for x in night.parts[-4:-1]))
        except ValueError:
            continue
        age = (today - date).days
        if age > delete_after_days:
            shutil.rmtree(str(logs))
            deleted += 1
        elif age > compress_after_days:
            changed = False
            for path in logs.rglob("*"):
                if path.is_file() and path.suffix != ".gz":
                    with open(str(path), "rb") as source, gzip.open(str(path) + ".gz", "wb") as target:
                        shutil.copyfileobj(source, target)
                    path.unlink()
                    changed = True
            compressed += changed
        else:
            continue
        write_index(night)
    return compressed, deleted


# ----------------------------------------------------------------------------- ingestion
def _candidates(incoming: Path) -> List[Path]:
    found = []
    if not incoming.is_dir():
        return found
    for first in sorted(incoming.iterdir()):
        if not first.is_dir() or first.name.startswith("."):
            continue
        if (first / DONE_FILE).is_file():
            found.append(first)
            continue
        found.extend(sorted(p for p in first.iterdir() if p.is_dir() and (p / DONE_FILE).is_file()))
    return found


def _tarballs(incoming: Path) -> List[Path]:
    if not incoming.is_dir():
        return []
    return sorted(p for p in incoming.glob("*" + SUFFIX) if p.is_file() and not p.name.startswith("."))


def _ingest_one(directory: Path, www: Path, store: "Store", now: dt.datetime) -> str:
    """Validate one run directory, move it into the archive and add it to the store; return its run ID."""
    run, results = read_run_dir(directory)
    if store.has_run(run["run_id"]):
        raise RecordError("run {} was already ingested".format(run["run_id"]))
    relative = archive_path(www, run)
    target = www / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(directory), str(target))
    try:
        store.add_run(run, results, str(relative), now)
    except Exception:
        shutil.move(str(target), str(directory))
        raise
    write_index(target)
    return run["run_id"]


def _reject(item: Path, root: Path, reason: str, now: dt.datetime) -> str:
    bin_dir = root / "rejected" / "{}-{}".format(item.name, now.strftime("%Y%m%dT%H%M%SZ"))
    bin_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        if item.is_dir():
            shutil.move(str(item), str(bin_dir))
        else:
            bin_dir.mkdir()
            shutil.move(str(item), str(bin_dir / item.name))
        (bin_dir / "REJECTED.txt").write_text(reason + "\n")
    except OSError as move_error:
        reason += " (could not move to rejected/: {})".format(move_error)
    return reason


def ingest(root: Path, www: Path, store: "Store", now: Optional[dt.datetime] = None) -> Tuple[List[str], List[Tuple[str, str]]]:
    """Move every complete upload (tarball or run directory) into the store and the archive;
    return (accepted run IDs, [(upload, reason)])."""
    root, www = Path(root), Path(www)
    now = now or dt.datetime.now(dt.timezone.utc)
    incoming = root / "incoming"
    accepted: List[str] = []
    rejected: List[Tuple[str, str]] = []
    for tarball in _tarballs(incoming):
        staging = incoming / ".unpack" / tarball.name[:-len(SUFFIX)]
        try:
            accepted.append(_ingest_one(unpack(tarball, staging), www, store, now))
            tarball.unlink()
        except (TransferError, RecordError, OSError, sqlite3.Error) as error:
            rejected.append((str(tarball), _reject(tarball, root, str(error), now)))
        finally:
            shutil.rmtree(str(staging), ignore_errors=True)
    shutil.rmtree(str(incoming / ".unpack"), ignore_errors=True)
    for directory in _candidates(incoming):
        try:
            accepted.append(_ingest_one(directory, www, store, now))
        except (RecordError, OSError, sqlite3.Error) as error:
            rejected.append((str(directory), _reject(directory, root, str(error), now)))
    return accepted, rejected
