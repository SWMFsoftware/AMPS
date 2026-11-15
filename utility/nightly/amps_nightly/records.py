"""Run-directory record format shared by every producer (runner, importers) and the server.

A run directory holds the results of ONE column (machine + compiler configuration) for one night:

    <run-id>/run.json        run metadata (machine, config, commit, times, completeness)
    <run-id>/results.jsonl   one JSON object per test
    <run-id>/logs/...        optional per-test logs referenced by results
    <run-id>/DONE            written last; the server ignores directories without it
"""

from __future__ import annotations

import datetime as dt
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

SCHEMA_VERSION = 1
STATUSES = ("PASS", "FAILED", "SKIPPED")
RUN_FILE = "run.json"
RESULTS_FILE = "results.jsonl"
DONE_FILE = "DONE"

_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,63}$")
_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,127}$")
_COMMIT = re.compile(r"^[0-9a-f]{7,40}$")
_STAGE = re.compile(r"^[a-z][a-z-]{0,19}$")
_CONTROL = re.compile(r"[\x00-\x08\x0b-\x1f\x7f]")


class RecordError(ValueError):
    """A run directory or record does not satisfy the format."""


# ----------------------------------------------------------------------------- times
def parse_time(value: Any) -> dt.datetime:
    """ISO-8601 time with an explicit offset (or Z)."""
    if not isinstance(value, str) or not value.strip():
        raise RecordError("time must be a non-empty ISO-8601 string")
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = dt.datetime.fromisoformat(text)
    except ValueError:
        raise RecordError("not an ISO-8601 time: {!r}".format(value))
    if parsed.tzinfo is None:
        raise RecordError("time needs a UTC offset: {!r}".format(value))
    return parsed


def utc_text(when: dt.datetime) -> str:
    return when.astimezone(dt.timezone.utc).replace(microsecond=0).strftime("%Y-%m-%dT%H:%M:%SZ")


def column_id(machine: str, config: str) -> str:
    return "{}-{}".format(machine, config)


def make_run_id(machine: str, config: str, finished: dt.datetime) -> str:
    return "{}-{}-{}".format(finished.astimezone(dt.timezone.utc).strftime("%Y%m%dT%H%MZ"), machine, config)


# ----------------------------------------------------------------------------- fields
def _text(raw: Dict[str, Any], key: str, limit: int, default: str = "", multiline: bool = False) -> str:
    value = raw.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, str):
        raise RecordError("{} must be a string".format(key))
    if not multiline and ("\n" in value or "\r" in value):
        raise RecordError("{} must be a single line".format(key))
    value = _CONTROL.sub("", value)
    if len(value) > limit:
        value = value[: limit - 1] + "…"
    return value


def _name(raw: Dict[str, Any], key: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not _NAME.match(value):
        raise RecordError("{} must match {}: got {!r}".format(key, _NAME.pattern, value))
    return value


def _count(raw: Dict[str, Any], key: str) -> Optional[int]:
    value = raw.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise RecordError("{} must be a non-negative integer".format(key))
    return value


def normalize_run(raw: Any) -> Dict[str, Any]:
    if not isinstance(raw, dict):
        raise RecordError("run.json must hold an object")
    if raw.get("schema_version") != SCHEMA_VERSION:
        raise RecordError("run.json schema_version must be {}".format(SCHEMA_VERSION))
    run_id = raw.get("run_id")
    if not isinstance(run_id, str) or not _RUN_ID.match(run_id):
        raise RecordError("run_id must match {}".format(_RUN_ID.pattern))
    if not isinstance(raw.get("complete"), bool):
        raise RecordError("complete must be true or false")
    commit = raw.get("commit") or ""
    if commit and (not isinstance(commit, str) or not _COMMIT.match(commit)):
        raise RecordError("commit must be a 7-40 character lowercase hex SHA")
    run = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "machine": _name(raw, "machine"),
        "config": _name(raw, "config"),
        "compiler": _text(raw, "compiler", 40),
        "toolchain": _text(raw, "toolchain", 200),
        "branch": _text(raw, "branch", 100),
        "commit": commit,
        "commit_inferred": bool(raw.get("commit_inferred", False)),
        "started": utc_text(parse_time(raw["started"])) if raw.get("started") else "",
        "finished": utc_text(parse_time(raw.get("finished"))),
        "complete": raw["complete"],
        "producer": _text(raw, "producer", 60),
        "notes": _text(raw, "notes", 500),
    }
    night = raw.get("night") or ""
    if night:
        try:
            dt.datetime.strptime(str(night), "%Y-%m-%d")
        except ValueError:
            raise RecordError("night must be a date YYYY-MM-DD (local date on which the night started)")
        run["night"] = str(night)
    run["column"] = column_id(run["machine"], run["config"])
    return run


def normalize_result(raw: Any, where: str = "") -> Dict[str, Any]:
    try:
        if not isinstance(raw, dict):
            raise RecordError("result must be an object")
        test = raw.get("test")
        if not isinstance(test, str) or not test.strip() or len(test) > 200 or _CONTROL.search(test) or "\n" in test:
            raise RecordError("test must be a single-line name of at most 200 characters")
        status = raw.get("status")
        if status not in STATUSES:
            raise RecordError("status must be one of {}".format(", ".join(STATUSES)))
        stage = raw.get("stage") or ""
        if stage and (not isinstance(stage, str) or not _STAGE.match(stage)):
            raise RecordError("stage must be a short lowercase word")
        duration = raw.get("duration_s")
        if duration is not None:
            if isinstance(duration, bool) or not isinstance(duration, (int, float)) or not math.isfinite(duration) or duration < 0:
                raise RecordError("duration_s must be a non-negative number")
            duration = float(duration)
        log = _text(raw, "log", 300)
        if log and (log.startswith("/") or "\\" in log or ".." in log.split("/")):
            raise RecordError("log must be a relative path inside the run directory")
        cases_total, cases_failed = _count(raw, "cases_total"), _count(raw, "cases_failed")
        if cases_total is not None and cases_failed is not None and cases_failed > cases_total:
            raise RecordError("cases_failed exceeds cases_total")
        return {
            "family": _name(raw, "family"),
            "test": test.strip(),
            "status": status,
            "stage": stage,
            "reason": _text(raw, "reason", 500),
            "source": _text(raw, "source", 40),
            "duration_s": duration,
            "cases_total": cases_total,
            "cases_failed": cases_failed,
            "log": log,
        }
    except RecordError as error:
        raise RecordError("{}{}".format(where + ": " if where else "", error))


# ----------------------------------------------------------------------------- run directories
def _write_atomic(path: Path, text: str) -> None:
    tmp = path.with_name("." + path.name + ".tmp")
    tmp.write_text(text)
    os.replace(str(tmp), str(path))


def write_run_dir(path: Path, run: Dict[str, Any], results: Iterable[Dict[str, Any]], done: bool = True) -> Path:
    """Validate and write a run directory; DONE is written last so readers never see a partial run."""
    path = Path(path)
    run = normalize_run(run)
    rows = [normalize_result(r, "result {}".format(i + 1)) for i, r in enumerate(results)]
    _check_unique(rows)
    path.mkdir(parents=True, exist_ok=True)
    stored = {k: v for k, v in run.items() if k != "column"}
    _write_atomic(path / RUN_FILE, json.dumps(stored, indent=2, sort_keys=True) + "\n")
    _write_atomic(path / RESULTS_FILE, "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows))
    if done:
        _write_atomic(path / DONE_FILE, run["run_id"] + "\n")
    return path


def read_run_dir(path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    path = Path(path)
    try:
        run = normalize_run(json.loads((path / RUN_FILE).read_text()))
    except (OSError, ValueError) as error:
        raise RecordError("{}: {}".format(RUN_FILE, error))
    rows: List[Dict[str, Any]] = []
    try:
        lines = (path / RESULTS_FILE).read_text().splitlines()
    except OSError as error:
        raise RecordError("{}: {}".format(RESULTS_FILE, error))
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        try:
            raw = json.loads(line)
        except ValueError as error:
            raise RecordError("{} line {}: {}".format(RESULTS_FILE, number, error))
        rows.append(normalize_result(raw, "{} line {}".format(RESULTS_FILE, number)))
    _check_unique(rows)
    return run, rows


def _check_unique(rows: List[Dict[str, Any]]) -> None:
    seen = set()
    for row in rows:
        key = (row["family"], row["test"])
        if key in seen:
            raise RecordError("duplicate result for {} / {}".format(*key))
        seen.add(key)


def unique_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Tidy the rows of one run so that it holds one row per (family, test): long texts are shortened,
    odd stages become "run", and a later duplicate gets its runner's name added to the test name."""
    def one_line(text: Any, limit: int) -> str:
        words = " ".join(str(text or "").split())
        return words if len(words) <= limit else words[:limit - 3] + "..."

    seen = set()
    for r in rows:
        r["test"] = one_line(r["test"], 180) or "?"
        r["reason"] = one_line(r.get("reason", ""), 480)
        if r.get("stage") and not re.match(r"^[a-z][a-z-]{0,19}$", r["stage"]):
            r["stage"] = "run"
        r.pop("log_src", None)
        key = (r["family"], r["test"])
        if key in seen:
            r["test"] = "{} [{}]".format(r["test"], r.get("source", "?"))
        seen.add((r["family"], r["test"]))
    return rows
