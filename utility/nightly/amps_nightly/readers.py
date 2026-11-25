"""Readers: turn the report a test runner already writes into result rows.

A reader gets a ReadContext and returns rows {test, status, stage, reason, duration_s, cases_total,
cases_failed, log_src[, family]}; log_src is the path of that test's own log file (or ""), which the
orchestrator copies into the archive. A reader raises ResultsUnavailable when the report is missing or
cannot be read. Adding a runner with a new report format means adding one function to READERS.
"""

from __future__ import annotations

import glob
import json
import xml.etree.ElementTree as ElementTree
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional


class ResultsUnavailable(Exception):
    """The runner's report is missing or unreadable."""


@dataclass
class ReadContext:
    runner: str
    file: str               # results.file with placeholders filled (may be a glob for junit)
    out: Path               # the runner's output directory
    exit_code: Optional[int]


Row = Dict[str, Any]
_STATUS = {"PASS": "PASS", "PASSED": "PASS", "OK": "PASS", "P": "PASS",
           "FAIL": "FAILED", "FAILED": "FAILED", "F": "FAILED", "ERROR": "FAILED",
           "SKIP": "SKIPPED", "SKIPPED": "SKIPPED"}


def one_line(text: Any, limit: int = 300) -> str:
    words = " ".join(str(text or "").split())
    return words if len(words) <= limit else words[:limit - 3] + "..."


def status_of(raw: Any) -> str:
    return _STATUS.get(str(raw or "").strip().upper(), "FAILED")


def stage_of(message: str) -> str:
    """Failing stage guessed from a runner's message."""
    text = message.lower()
    if any(w in text for w in ("did not build", "build fail", "compil", "config.pl", "make amps", "link")):
        return "build"
    if any(w in text for w in ("timed out", "time limit", "timeout", "crash", "signal", "did not finish")):
        return "run"
    return "check"


def row(test: str, status: str, reason: str = "", stage: str = "", duration: Any = None, cases_total: Any = None,
        cases_failed: Any = None, log_src: Any = "", family: str = "") -> Row:
    result = {"test": test, "status": status, "reason": one_line(reason) if status != "PASS" or reason else "",
              "stage": stage if status == "FAILED" else "",
              "duration_s": float(duration) if isinstance(duration, (int, float)) and duration >= 0 else None,
              "cases_total": cases_total if isinstance(cases_total, int) else None,
              "cases_failed": cases_failed if isinstance(cases_failed, int) else None,
              "log_src": str(log_src or "")}
    if status == "FAILED" and not result["stage"]:
        result["stage"] = stage_of(result["reason"])
    if family:
        result["family"] = family
    return result


def _json(path: Path) -> Any:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ResultsUnavailable("cannot read {}: {}".format(path, error))


def _last_line(path: Any, words: tuple = ("FAIL", "ERROR", "Error", "error:")) -> str:
    if not path:
        return ""
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError:
        return ""
    for line in reversed(lines[-400:]):
        if any(w in line for w in words) and not line.startswith("###"):
            return one_line(line)
    return ""


# ----------------------------------------------------------------------------- nested checks
# src/test suites whose one registered case runs a whole nested suite. Its receipt lists every check,
# and each check becomes a row of its own family. The suite's own row keeps the suite's verdict; these
# rows only make the individual checks visible, with their own history.
#   receipt               path of the nested receipt, relative to the suite's artifact_dir
#   schema                the receipt's "schema" value; anything else is reported, never guessed at
#   family                family of the check rows (registered in families.json)
#   out_of_scope_reason   reason shown for each ID the receipt lists in outside_selected_scope
NESTED_CHECKS: Dict[str, Dict[str, str]] = {
    "CORE-CHEMICAL-DATABASE-001": {
        "receipt": "shared-protocol/result.json",
        "schema": "amps-chemical-database-software-tests-v1",
        "family": "core-chemical-database",
        "out_of_scope_reason": "not run by the nightly: --protocol-only leaves out the checks that need a pinned "
                               "database release (--database-lock) or the original source files; see "
                               "test/chemical_database/README.md",
    },
}
_CHILD_STATUSES = ("PASS", "FAIL", "ERROR", "SKIPPED")


def nested_check_rows(suite: Dict[str, Any], spec: Dict[str, str], base: Path) -> List[Row]:
    """Rows of the checks listed in a suite's nested receipt (see NESTED_CHECKS).

    base is the report's output directory, for a relative artifact_dir. No receipt is expected when the
    suite was skipped or the harness reported ERROR, since its own row already says why. Otherwise a
    missing, unreadable or self-contradictory receipt gives one FAILED row "<suite> checks", so the page
    never shows fewer checks without saying so.
    """
    suite_id = str(suite.get("id"))
    family = spec["family"]
    raw = str(suite.get("status", "")).upper()

    def unavailable(why: str) -> List[Row]:
        return [row(suite_id + " checks", "FAILED", "nested check results unavailable: " + why,
                    stage="check", family=family)]

    artifact = Path(str(suite.get("artifact_dir") or ""))
    receipt = (artifact if artifact.is_absolute() else base / artifact) / spec["receipt"]
    if not suite.get("artifact_dir") or not receipt.is_file():
        return [] if raw in ("SKIP", "SKIPPED", "ERROR") else unavailable("no receipt {}".format(receipt))
    try:
        data = json.loads(receipt.read_text())
    except (OSError, ValueError) as error:
        return unavailable("cannot read {}: {}".format(receipt, error))
    if not isinstance(data, dict) or data.get("schema") != spec["schema"]:
        return unavailable("{} is not a {} receipt".format(receipt, spec["schema"]))
    records, counts = data.get("tests"), data.get("counts")
    excluded = data.get("outside_selected_scope", [])
    if (not isinstance(records, list) or not all(isinstance(r, dict) for r in records)
            or not isinstance(counts, dict) or not isinstance(excluded, list)
            or not all(isinstance(e, str) and e for e in excluded)):
        return unavailable("{} lacks its check list, counts or out-of-scope list".format(receipt))
    tally = {s: sum(1 for r in records if r.get("status") == s) for s in _CHILD_STATUSES}
    if (any(r.get("status") not in _CHILD_STATUSES for r in records) or tally != counts
            or data.get("executed") != len(records)):
        return unavailable("{}: check records contradict its counts".format(receipt))
    names = [str(r.get("test") or "") for r in records] + excluded
    if not all(names) or len(set(names)) != len(names):
        return unavailable("{}: a check name is empty or repeated".format(receipt))
    rows = []
    for record in records:
        status = str(record["status"])
        reason = record.get("reason") or ""
        if status == "ERROR":
            reason = "error: " + str(reason or "check reported ERROR")
        rows.append(row(str(record["test"]), status_of(status), reason, log_src=record.get("log") or "",
                        family=family))
    rows += [row(name, "SKIPPED", spec["out_of_scope_reason"], family=family) for name in excluded]
    return rows


# ----------------------------------------------------------------------------- readers
def src_test_json(ctx: ReadContext) -> List[Row]:
    """src/test/run_tests.py --json: one row per suite, with its case counts.

    Suites in NESTED_CHECKS also give one row per nested check (nested_check_rows).
    """
    data = _json(Path(ctx.file))
    base = Path(data.get("output_directory") or Path(ctx.file).parent)
    rows = []
    for suite in data.get("results", []):
        cases = suite.get("cases") or []
        failed = sum(1 for c in cases if str(c.get("status", "")).upper() in ("FAIL", "ERROR"))
        raw = str(suite.get("status", "")).upper()
        message = one_line(suite.get("message"))
        reason = ("harness: " + message) if raw == "ERROR" else message
        log = suite.get("log") or ""
        rows.append(row(str(suite.get("id")), status_of(raw), reason, duration=suite.get("seconds"),
                        cases_total=len(cases) or None, cases_failed=failed if cases else None,
                        log_src=(log if not log or Path(log).is_absolute() else base / log)))
        spec = NESTED_CHECKS.get(str(suite.get("id")))
        if spec:
            rows += nested_check_rows(suite, spec, base)
    for item in data.get("not_executed_suites", []) or []:
        sid = item.get("id") if isinstance(item, dict) else str(item)
        rows.append(row(sid, "FAILED", "not executed by run_tests.py", stage="run"))
    return rows


def sep_orchestrator_json(ctx: ReadContext) -> List[Row]:
    """tools/sep_test_orchestrator.py run_all_test.json: one row per record (build checks included)."""
    data = _json(Path(ctx.file))
    records = data.get("records", [])
    ids = [r.get("test_id") for r in records]
    rows = []
    for record in records:
        test = str(record.get("test_id"))
        if ids.count(test) > 1:
            test = "{}/{}".format(record.get("step"), test)
        stage = "build" if record.get("kind") == "check" and record.get("step") == "build" else ""
        message = record.get("message") or ""
        if str(record.get("status", "")).upper() == "ERROR":
            message = "error: " + str(message or "runner reported ERROR")
        rows.append(row(test, status_of(record.get("status")), message, stage=stage, log_src=record.get("log", "")))
    return rows


_INTERPRETERS = {"python", "sh", "bash", "env", "perl", "time"}      # also python3, python3.8


def earth_report(ctx: ReadContext) -> List[Row]:
    """srcEarth/test/test_runner.py <prefix>.json: one row per list entry."""
    data = _json(Path(ctx.file))
    if not isinstance(data, list):
        raise ResultsUnavailable("{} is not a list of entries".format(ctx.file))

    def name(entry: Dict[str, Any]) -> str:
        # the script, after an interpreter, its options and VAR=value settings ("python3 srcEarth/test/SA01/..."
        # was named "python3 (line 233)")
        words = str(entry.get("command") or "").split()
        i = 0
        while i < len(words) - 1 and (Path(words[i]).name.rstrip("0123456789.") in _INTERPRETERS or "=" in words[i]
                                      or (i and words[i].startswith("-"))):
            i += 1
        first = words[i] if words else "entry"
        parts = Path(first).parts
        return parts[2] if len(parts) > 3 and parts[:2] == ("srcEarth", "test") else Path(first).stem

    names = [name(e) for e in data]
    rows = []
    for entry, test in zip(data, names):
        if names.count(test) > 1:
            test = "{} (line {})".format(test, entry.get("line_no"))
        if entry.get("not_started"):          # the runner was stopped (its time limit) before this entry started
            rows.append(row(test, "SKIPPED", "not started: the runner was stopped before this entry"))
            continue
        passed = entry.get("actual") == "P"
        reason = ""
        if not passed:
            if entry.get("interrupted"):
                reason = "stopped after {:.0f} s: the runner was stopped".format(entry.get("elapsed_s") or 0)
            elif entry.get("timed_out"):
                reason = "timed out after {:.0f} s".format(entry.get("elapsed_s") or 0)
            else:
                reason = _last_line(entry.get("log_file")) or "exit code {}".format(entry.get("exit_code"))
            if entry.get("expected") == "F":
                reason = "expected to fail (srcEarth/test/list): " + reason
        elif entry.get("expected") == "F":
            reason = "passes although srcEarth/test/list expects it to fail"
        rows.append(row(test, "PASS" if passed else "FAILED", reason,
                        stage="run" if entry.get("timed_out") or entry.get("interrupted") else "",
                        duration=entry.get("elapsed_s"),
                        log_src=entry.get("log_file", "")))
    return rows


def moon_summary(ctx: ReadContext) -> List[Row]:
    """srcMoon/test/run_tests.py summary.json (+ <out>/<id>/result.json): one row per U-test."""
    path = Path(ctx.file)
    data = _json(path)
    rows = []
    for item in data.get("results", []):
        record = item if isinstance(item, dict) else {"id": str(item)}
        tid = str(record.get("id"))
        detail = path.parent / tid / "result.json"
        if detail.is_file():
            try:
                record = dict(json.loads(detail.read_text()), **{k: v for k, v in record.items() if v})
            except ValueError:
                pass
        raw = str(record.get("status", "")).upper()
        reason = record.get("reason") or record.get("message") or ""
        if raw == "ERROR":
            reason = "error: " + str(reason or "test reported ERROR")
        log = next((p for p in (path.parent / tid / "stdout.log", detail) if p.is_file()), "")
        rows.append(row(tid, status_of(raw), reason, log_src=log))
    return rows


def enceladus_receipt(ctx: ReadContext) -> List[Row]:
    """srcEnceladus/test/run_tests.py test_receipt.json: one row per result (PASS or FAIL)."""
    data = _json(Path(ctx.file))
    rows = []
    for item in data.get("results", []):
        detail = [l for l in str(item.get("detail") or "").splitlines() if l.strip()]
        first = next((l for l in detail if l.startswith(("FAIL:", "ERROR:"))), "")
        failing = len([l for l in detail if l.startswith(("FAIL:", "ERROR:"))])
        reason = "{}{}".format(first, " (and {} more)".format(failing - 1) if failing > 1 else "") if first else \
            (detail[-1] if detail else "")
        rows.append(row(str(item.get("test")), status_of(item.get("status")), reason))
    return rows


def junit(ctx: ReadContext) -> List[Row]:
    files = sorted(glob.glob(ctx.file, recursive=True))
    if not files:
        raise ResultsUnavailable("no JUnit file matches {}".format(ctx.file))
    rows = []
    for path in files:
        try:
            root = ElementTree.parse(path).getroot()
        except (OSError, ElementTree.ParseError) as error:
            raise ResultsUnavailable("cannot read {}: {}".format(path, error))
        for case in root.iter("testcase"):
            test = ".".join(x for x in (case.get("classname"), case.get("name")) if x)
            bad = case.find("failure") if case.find("failure") is not None else case.find("error")
            if bad is not None:
                status, reason = "FAILED", bad.get("message") or (bad.text or "").strip()
            elif case.find("skipped") is not None:
                status, reason = "SKIPPED", case.find("skipped").get("message") or ""
            else:
                status, reason = "PASS", ""
            try:
                seconds = float(case.get("time", ""))
            except ValueError:
                seconds = None
            rows.append(row(test, status, reason, duration=seconds))
    return rows


def exit_code(ctx: ReadContext) -> List[Row]:
    ok = ctx.exit_code == 0
    return [row(ctx.runner, "PASS" if ok else "FAILED", "" if ok else "exit code {}".format(ctx.exit_code), stage="run")]


def amps_results(ctx: ReadContext) -> List[Row]:
    """A runner that writes the orchestrator's own result rows, one JSON object per line."""
    try:
        lines = Path(ctx.file).read_text().splitlines()
    except OSError as error:
        raise ResultsUnavailable("cannot read {}: {}".format(ctx.file, error))
    rows = []
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        try:
            item = json.loads(line)
        except ValueError as error:
            raise ResultsUnavailable("{} line {}: {}".format(ctx.file, number, error))
        rows.append(row(str(item.get("test")), status_of(item.get("status")), item.get("reason", ""),
                        stage=item.get("stage", ""), duration=item.get("duration_s"),
                        cases_total=item.get("cases_total"), cases_failed=item.get("cases_failed"),
                        log_src=item.get("log", ""), family=item.get("family", "")))
    return rows


READERS: Dict[str, Callable[[ReadContext], List[Row]]] = {
    "src-test-json": src_test_json,
    "sep-orchestrator-json": sep_orchestrator_json,
    "earth-report": earth_report,
    "moon-summary": moon_summary,
    "enceladus-receipt": enceladus_receipt,
    "junit": junit,
    "exit-code": exit_code,
    "amps-results": amps_results,
}
