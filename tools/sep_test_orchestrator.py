#!/usr/bin/env python3
"""Build srcSEP or srcSEP3D natively, run all of its test procedures, and
report how many tests ran, passed, failed, and were skipped, with the list of
failed tests (name, runner, message, and test-log location).

This is the shared implementation behind ``srcSEP/test/run_all_test.sh`` and
``srcSEP3D/test/run_all_test.sh`` (each wrapper calls it with ``--app``).

Design
------
* A *step* is one command (or the native build) from the application's step
  table.  Steps are either ``test`` steps, which run many tests and write a
  report, or ``check`` steps (the build, initialization previews), which are a
  single pass/fail outcome taken from the exit status.
* After a test step finishes, its *result reader* converts the report that
  step's runner already writes (JSON where available, otherwise the runner's
  stable text format) into uniform :class:`TestRecord` objects with a status
  in PASS/FAIL/SKIP/ERROR, the runner's name, the message, and a test-log
  location.  Runners are not re-implemented; the reader only interprets them.
* Every failed test gets a log location, in this order of preference:
  1. the runner's own per-test log or diagnostic file;
  2. the test's section copied from a shared runner log, or its captured
     output, into ``OUTPUT_DIR/failed-tests/<step>/<test>.log`` (each copied
     file starts with a header naming its source and line range);
  3. the step log and the line holding the test's failure line (used where a
     runner genuinely cannot separate tests, e.g. one MPI process running
     several native tests).
* Fail-safe accounting: a test step whose report is missing or unreadable
  contributes one ERROR record ``<step>:results-unavailable``; a step that
  exits nonzero although its report contains no failing test contributes one
  ERROR record ``<step>:exit-status``.  A failure therefore cannot look green.
* Counting: "test executions" is the sum over all steps.  "Unique tests"
  de-duplicates records that carry the same (namespace, test ID); only steps
  that genuinely run the same tests share a namespace (for srcSEP3D, the
  stage-1 registry tests that both run-tests-all and pd-mover execute).  The
  worst status wins (FAIL > ERROR > SKIP > PASS).

Outputs (OUTPUT_DIR defaults to test_output/run_all_test/<app>-<UTC time>)
-------
* ``NN-<step>.log``: complete terminal output of each step;
* ``run_all_test.log`` (or ``--log``): everything printed to the terminal;
* ``run_all_test.json``: all steps, counts, and every test record;
* ``failed_tests.txt``: the failed-test list in plain text;
* ``failed-tests/<step>/<test>.log``: copied per-test logs (see above);
* optional JUnit XML (``--junit FILE``).

Exit status: 0 when every executed check passed and no test FAILed or
ERRORed; 1 otherwise; 2 for usage or preflight errors.

Requires Python 3.8 or newer and only the standard library.
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import json
import os
import re
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, TextIO, Tuple
from xml.sax.saxutils import escape as _xml_escape

AMPS_ROOT = Path(__file__).resolve().parents[1]

STATUSES = ("PASS", "FAIL", "SKIP", "ERROR")
# Worst status wins when the same test is seen more than once.
SEVERITY = {"PASS": 0, "SKIP": 1, "ERROR": 2, "FAIL": 3}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class TestRecord:
    """One test (or check) outcome in the uniform vocabulary.

    test_id    name as the runner reports it (scoped runners use scope/ID)
    status     PASS, FAIL, SKIP, or ERROR (normalized by normalize_status)
    message    runner message or failure reason, single line
    step       step that ran it; runner = command/script that produced it
    log, line  test-log location; line is 1-based, 0 = whole file
    evidence   other per-test files (JSON/XML reports, artifact directories)
    namespace  de-duplication scope for "unique tests" (see module docstring)
    kind       "test" or "check"
    """
    test_id: str
    status: str
    message: str
    step: str
    runner: str
    log: str = ""
    line: int = 0
    evidence: List[str] = field(default_factory=list)
    namespace: str = ""
    kind: str = "test"


class ResultsUnavailable(Exception):
    """Raised by a reader when its step's report is missing or unreadable."""


@dataclass
class StepContext:
    """What a result reader needs to know about one finished step."""
    step: "Step"
    output_dir: Path          # OUTPUT_DIR of the whole run
    step_log: Path            # complete terminal output of this step
    exit_code: int
    started: float           # time.time() when the step started
    failed_dir: Path          # OUTPUT_DIR/failed-tests/<step>

    def record(self, test_id: str, status: str, message: str = "", *,
               log: str = "", line: int = 0, evidence: Sequence[str] = (),
               runner: str = "") -> TestRecord:
        """Build a TestRecord with this step's defaults filled in."""
        return TestRecord(
            test_id=test_id, status=normalize_status(status),
            message=one_line(message), step=self.step.name,
            runner=runner or self.step.runner, log=log or str(self.step_log),
            line=line, evidence=[str(item) for item in evidence],
            namespace=self.step.namespace or self.step.name,
            kind=self.step.kind)


Reader = Callable[[StepContext], List[TestRecord]]


@dataclass
class Step:
    """One entry of an application's step table.

    commands  argv lists run in order from the AMPS root; the first nonzero
              exit stops the step.  ``pre`` (optional) runs first with the
              step log and returns an exit code (used for the build preflight).
    reader    converts the step's report into TestRecords (test steps only).
    runner    label shown for this step's tests (script or command name).
    long      skipped by --quick.
    make_dirs directories created before the step (native --test-json parents).
    """
    name: str
    description: str
    commands: List[List[str]]
    kind: str = "test"
    reader: Optional[Reader] = None
    runner: str = ""
    namespace: str = ""
    long: bool = False
    env: Dict[str, str] = field(default_factory=dict)
    pre: Optional[Callable[[TextIO], int]] = None
    make_dirs: List[Path] = field(default_factory=list)


@dataclass
class StepRun:
    step: Step
    status: str               # PASS/FAIL for checks and test steps, SKIP, DRY
    exit_code: Optional[int]
    seconds: float
    log: str
    records: List[TestRecord] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def normalize_status(raw: object) -> str:
    """Map the runners' spellings onto PASS/FAIL/SKIP/ERROR.

    passed/pass/ok -> PASS; failed/fail/failure -> FAIL; skipped/skip ->
    SKIP; errors/error/MISSING (pd-mover) -> ERROR.  Anything unknown is an
    ERROR, so an unexpected runner state is never counted as a pass."""
    text = str(raw).strip().upper()
    if text in ("PASS", "PASSED", "OK"):
        return "PASS"
    if text in ("FAIL", "FAILED", "FAILURE"):
        return "FAIL"
    if text in ("SKIP", "SKIPPED"):
        return "SKIP"
    return "ERROR"


def one_line(text: object, limit: int = 300) -> str:
    """Collapse whitespace and bound the length of a message for reports."""
    value = " ".join(str(text or "").split())
    return value if len(value) <= limit else value[:limit - 3] + "..."


def utc_stamp() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def read_json(path: Path) -> dict:
    """Load a report; any problem becomes ResultsUnavailable for the step."""
    if not path.is_file():
        raise ResultsUnavailable(f"report not found: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError) as error:
        raise ResultsUnavailable(f"cannot read report {path}: {error}")


def find_line(path: Path, pattern: str) -> int:
    """1-based number of the first line matching regex pattern (0 if none)."""
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            regex = re.compile(pattern)
            for number, line in enumerate(stream, start=1):
                if regex.search(line):
                    return number
    except OSError:
        pass
    return 0


def safe_name(test_id: str) -> str:
    """File-system-safe name for a copied per-test log."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", test_id) or "test"


def write_test_log(ctx: StepContext, test_id: str, body: str, *,
                   source: str, first: int = 0, last: int = 0) -> str:
    """Write OUTPUT_DIR/failed-tests/<step>/<test>.log and return its path.

    The header names the test, step, runner, and the source the body was
    copied from (with its line range when it is a section of a shared log),
    so the extracted file can always be traced back."""
    ctx.failed_dir.mkdir(parents=True, exist_ok=True)
    target = ctx.failed_dir / (safe_name(test_id) + ".log")
    span = f" lines {first}-{last}" if first else ""
    header = (f"# test:   {test_id}\n# step:   {ctx.step.name}\n"
              f"# runner: {ctx.step.runner}\n# source: {source}{span}\n\n")
    target.write_text(header + body, encoding="utf-8")
    return str(target)


def worst(records: Iterable[TestRecord]) -> Dict[Tuple[str, str], TestRecord]:
    """De-duplicate by (namespace, test_id), keeping the worst status."""
    merged: Dict[Tuple[str, str], TestRecord] = {}
    for record in records:
        key = (record.namespace, record.test_id)
        if key not in merged or SEVERITY[record.status] > SEVERITY[merged[key].status]:
            merged[key] = record
    return merged


def count(records: Iterable[TestRecord]) -> Dict[str, int]:
    totals = {status: 0 for status in STATUSES}
    for record in records:
        totals[record.status] += 1
    totals["TOTAL"] = sum(totals[status] for status in STATUSES)
    return totals


# ---------------------------------------------------------------------------
# Result readers (one per runner report format)
# ---------------------------------------------------------------------------

def read_check(ctx: StepContext) -> List[TestRecord]:
    """Check steps (build, init previews): one outcome from the exit code."""
    status = "PASS" if ctx.exit_code == 0 else "FAIL"
    message = "" if ctx.exit_code == 0 else f"exit code {ctx.exit_code}"
    return [ctx.record(ctx.step.name, status, message)]


def reader_srcsep_run_tests(report_dir: Path) -> Reader:
    """srcSEP/test/run_tests.py: report_dir/srcsep-tests.json.

    results[] = {id, status, message, artifacts}.  All tests share one
    report_dir/test-run.log in which each test starts with a timestamped
    ``] RUN ... --test ID`` (native registry) or ``--case ID`` (validation
    case) line; a failed test's sections are copied to its own log.  Its
    report_dir/individual/<ID>/<ID>.json is listed as evidence."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(report_dir / "srcsep-tests.json")
        shared = report_dir / "test-run.log"
        lines: List[str] = []
        if shared.is_file():
            lines = shared.read_text(encoding="utf-8", errors="replace").splitlines(True)
        records = []
        for item in payload.get("results", []):
            test_id = str(item.get("id", "?"))
            status = normalize_status(item.get("status"))
            evidence = [str(path) for path in
                        (report_dir / "individual" / test_id / f"{test_id}.json",)
                        if path.is_file()]
            log, line = str(shared), 0
            if status in ("FAIL", "ERROR"):
                log, line = extract_marked_sections(ctx, test_id, shared, lines)
            records.append(ctx.record(test_id, status, item.get("message", ""),
                                      log=log, line=line, evidence=evidence))
        return records
    return read


def extract_marked_sections(ctx: StepContext, test_id: str, shared: Path,
                            lines: List[str]) -> Tuple[str, int]:
    """Copy every ``] RUN ... --test|--case ID`` section of a shared log.

    A section runs from its RUN line to the line before the next RUN line.
    Returns the copied file (or the shared log and 0 when no marker exists)."""
    start = re.compile(r"\] RUN .*--(?:test|case) " + re.escape(test_id) + r"(?:\s|$)")
    any_run = re.compile(r"^\[[^\]]+\] RUN ")
    sections: List[Tuple[int, int]] = []
    index = 0
    while index < len(lines):
        if start.search(lines[index]):
            end = index + 1
            while end < len(lines) and not any_run.search(lines[end]):
                end += 1
            sections.append((index, end))
            index = end
        else:
            index += 1
    if not sections:
        return str(shared), 0
    body = "".join("".join(lines[a:b]) + "\n" for a, b in sections)
    path = write_test_log(ctx, test_id, body, source=str(shared),
                          first=sections[0][0] + 1, last=sections[-1][1])
    return path, 0


def reader_srcsep3d_run_tests(report_dir: Path) -> Reader:
    """srcSEP3D/test/run_tests.py: report_dir/srcsep3d-tests.json.

    results[] = {test_id, status, message, log}; ``log`` is the runner's own
    per-test log (report_dir/logs/<ID>.log).  Reports from runner versions
    without it fall back to the step log line ``[ID] FAIL|ERROR``."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(report_dir / "srcsep3d-tests.json")
        records = []
        for item in payload.get("results", []):
            test_id = str(item.get("test_id", "?"))
            log, line = str(item.get("log") or ""), 0
            if not log:
                log = str(ctx.step_log)
                line = find_line(ctx.step_log, r"^\[" + re.escape(test_id) + r"\] (FAIL|ERROR)")
            evidence = [str(path) for path in (report_dir / f"{test_id}.json",)
                        if path.is_file()]
            records.append(ctx.record(test_id, item.get("status"),
                                      item.get("message", ""), log=log,
                                      line=line, evidence=evidence))
        return records
    return read


def reader_pd_mover(report: Path, script: str,
                    namespaces: Optional[Dict[str, str]] = None) -> Reader:
    """run_parallel_diffusion_mover_tests.py --json REPORT.

    suites[].results[] = {test_id, status (incl. MISSING), suite, detail,
    log, line}; ``log`` is the saved suite output and ``line`` the result
    line.  The runner label names the suite; ``namespaces`` maps suites that
    re-run tests counted elsewhere (e.g. stage-1) onto a shared namespace."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(report)
        records = []
        for suite in payload.get("suites", []):
            for item in suite.get("results", []):
                suite_name = str(item.get("suite") or suite.get("name") or "")
                record = ctx.record(
                    str(item.get("test_id", "?")), item.get("status"),
                    item.get("detail", ""), log=str(item.get("log") or ""),
                    line=int(item.get("line") or 0),
                    runner=f"{script} [{suite_name}]")
                record.namespace = (namespaces or {}).get(
                    suite_name, f"{ctx.step.name}/{suite_name}")
                records.append(record)
        return records
    return read


def reader_native_json(report: Path) -> Reader:
    """Native ``amps --test ... --test-json REPORT``: results[] = {id, status,
    message, artifacts}.  All selected tests share one MPI process, so their
    output is the step log; the test's own ``[ID]`` line is located there and
    its artifacts are listed as evidence."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(report)
        records = []
        for item in payload.get("results", []):
            test_id = str(item.get("id", "?"))
            line = find_line(ctx.step_log, r"\[" + re.escape(test_id) + r"\]")
            evidence = [str(path) for path in item.get("artifacts", []) or []]
            evidence.insert(0, str(report))
            records.append(ctx.record(test_id, item.get("status"),
                                      item.get("message", ""), line=line,
                                      evidence=evidence))
        return records
    return read


def reader_coupled(output_dir: Path) -> Reader:
    """run_coupled_sep_corona.py --output-dir DIR: DIR/summary.json.

    results[] = {id, scope, status, message, output}; the runner already
    writes run_directory/failures/<scope>/<id>.log for each failed case, which
    is used directly; otherwise the captured output is copied."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(output_dir / "summary.json")
        run_dir = Path(str(payload.get("run_directory") or output_dir))
        records = []
        for item in payload.get("results", []):
            scope, raw_id = str(item.get("scope", "")), str(item.get("id", "?"))
            test_id = f"{scope}/{raw_id}" if scope else raw_id
            status = normalize_status(item.get("status"))
            log = ""
            if status in ("FAIL", "ERROR"):
                diagnostic = run_dir / "failures" / scope / f"{raw_id}.log"
                if diagnostic.is_file():
                    log = str(diagnostic)
                elif item.get("output"):
                    log = write_test_log(ctx, test_id, str(item["output"]),
                                         source=f"{output_dir / 'summary.json'} (captured output)")
            records.append(ctx.record(test_id, status, item.get("message", ""),
                                      log=log, evidence=item.get("artifacts") or []))
        return records
    return read


def reader_reduced_front(output_dir: Path) -> Reader:
    """run_reduced_shock_front.py --output-dir DIR: DIR/summary.json.
    results[] = {scope, test_id, status, message, log, report}; ``log`` is the
    runner's phase log for that test and ``report`` its report, if any."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(output_dir / "summary.json")
        records = []
        for item in payload.get("results", []):
            scope, raw_id = str(item.get("scope", "")), str(item.get("test_id", "?"))
            records.append(ctx.record(
                f"{scope}/{raw_id}" if scope else raw_id, item.get("status"),
                item.get("message", ""), log=str(item.get("log") or ""),
                evidence=[item["report"]] if item.get("report") else []))
        return records
    return read


def reader_positive_shock(output_dir: Path) -> Reader:
    """validate_positive_shock_example.py --output-root DIR: DIR/summary.json.
    checks[] = {run, name, status, detail, log}; ``log`` is that run's
    execution.log.  Each check is one test, named run/check."""
    def read(ctx: StepContext) -> List[TestRecord]:
        payload = read_json(output_dir / "summary.json")
        return [ctx.record(f"{item.get('run', '')}/{item.get('name', '?')}",
                           item.get("status"), item.get("detail", ""),
                           log=str(item.get("log") or ""))
                for item in payload.get("checks", [])]
    return read


def reader_coronal_cme(results_json: Path) -> Reader:
    """src/models/sep_coronal_cme ``make test``: build/test-results/results.json.

    tests[] = {id, status, message, output}.  The report lives at a fixed
    path, so a report older than this step is stale (e.g. the build failed
    before the tests ran) and is treated as unavailable.  A failed test's
    captured output is copied to its own log."""
    def read(ctx: StepContext) -> List[TestRecord]:
        if results_json.is_file() and results_json.stat().st_mtime < ctx.started - 1:
            raise ResultsUnavailable(f"report is older than this step: {results_json}")
        payload = read_json(results_json)
        records = []
        for item in payload.get("tests", []):
            test_id = str(item.get("id", "?"))
            status = normalize_status(item.get("status"))
            log = ""
            if status in ("FAIL", "ERROR") and item.get("output"):
                log = write_test_log(ctx, test_id, str(item["output"]),
                                     source=f"{results_json} (captured output)")
            records.append(ctx.record(test_id, status, item.get("message", ""),
                                      log=log, evidence=[str(results_json)]))
        return records
    return read


_UNITTEST_NAME = re.compile(r"^(test\w*) \(([^)]+)\)")
_UNITTEST_VERDICT = re.compile(
    r"\.\.\. (ok|FAIL|ERROR|skipped|expected failure|unexpected success)\b")
_UNITTEST_BLOCK = re.compile(r"^(FAIL|ERROR): (test\w*) \(([^)]+)\)")


def read_unittest(ctx: StepContext) -> List[TestRecord]:
    """Python unittest with verbosity 2 (step log text).

    Each test prints ``test_x (module.Class)`` and, possibly after a
    docstring line, ``... ok|FAIL|ERROR|skipped``.  Failure details follow in
    blocks that start with ``FAIL: test_x (...)`` after a ``=====`` line;
    each block is copied to that test's log."""
    lines = ctx.step_log.read_text(encoding="utf-8", errors="replace").splitlines()
    records: Dict[str, TestRecord] = {}
    pending: Optional[Tuple[str, int]] = None
    for number, line in enumerate(lines, start=1):
        name = _UNITTEST_NAME.match(line)
        if name:
            pending = (f"{name.group(2).split('.')[-1]}.{name.group(1)}"
                       if "." in name.group(2) else name.group(1), number)
        verdict = _UNITTEST_VERDICT.search(line)
        if verdict and pending:
            word = verdict.group(1)
            status = {"ok": "PASS", "expected failure": "PASS", "FAIL": "FAIL",
                      "ERROR": "ERROR", "skipped": "SKIP",
                      "unexpected success": "FAIL"}[word]
            records[pending[0]] = ctx.record(pending[0], status, line.strip(),
                                             line=pending[1])
            pending = None
    # Copy each failure block and attach it to its test.
    for number, line in enumerate(lines):
        block = _UNITTEST_BLOCK.match(line)
        if not block:
            continue
        end = number + 1
        while end < len(lines) and not lines[end].startswith(("=" * 20, "-" * 20)):
            end += 1
        key = (f"{block.group(3).split('.')[-1]}.{block.group(2)}"
               if "." in block.group(3) else block.group(2))
        path = write_test_log(ctx, key, "\n".join(lines[number:end]) + "\n",
                              source=str(ctx.step_log), first=number + 1, last=end)
        record = records.get(key) or ctx.record(key, block.group(1))
        record.log, record.line = path, 0
        records[key] = record
    if not records:
        raise ResultsUnavailable("no unittest result lines in the step log")
    return list(records.values())


_BRACKET = re.compile(r"^\[([A-Za-z0-9_.:/-]+)\] (PASS|FAIL|SKIP|ERROR)\b(.*)$")
_NAMED = re.compile(r"^([A-Za-z0-9_.-]+) (PASS|FAIL|SKIP|ERROR): (.*)$")


def read_bracket_text(ctx: StepContext) -> List[TestRecord]:
    """Text suites (sep_corona_swcme ``make test``): ``[ID] STATUS ...`` lines
    and ``NAME STATUS: ...`` verdict lines (``BG3D-4 FAIL: ...``).  These
    binaries share the step log; each record points to its own line there."""
    records: Dict[str, TestRecord] = {}
    with ctx.step_log.open(encoding="utf-8", errors="replace") as stream:
        for number, line in enumerate(stream, start=1):
            line = line.rstrip("\n")
            match = _BRACKET.match(line) or _NAMED.match(line)
            if not match:
                continue
            record = ctx.record(match.group(1), match.group(2),
                                match.group(3).strip(), line=number)
            previous = records.get(record.test_id)
            if previous is None or SEVERITY[record.status] > SEVERITY[previous.status]:
                records[record.test_id] = record
    if not records:
        raise ResultsUnavailable("no [ID] STATUS result lines in the step log")
    return list(records.values())


_MAKE_IGNORED = re.compile(r"^make(?:\[\d+\])?: \[.*\] Error (\d+) \(ignored\)")


def reader_make_recipe(binary_pattern: str) -> Reader:
    """``make -i ... test`` whose recipe runs several test binaries in turn.

    A make recipe stops at its first failing command and ``-k`` does not
    change that, so the step runs make with ``-i`` (ignore recipe errors):
    every binary then runs.  make echoes each recipe command before running
    it; lines matching ``binary_pattern`` (e.g. ``build/test_bg3d1``) split
    the step log into per-binary segments.  Within a segment the
    ``[ID] STATUS`` / ``NAME STATUS:`` lines become records (see
    read_bracket_text).  Because -i hides failures from the exit status, a
    segment that ends in make's ``Error N (ignored)`` without any FAIL/ERROR
    record becomes one ERROR record named after the binary; a segment with
    no result lines and no error counts as one PASS for that binary.  Every
    binary is therefore counted at least once and no failure is lost."""
    echo = re.compile(binary_pattern)

    def read(ctx: StepContext) -> List[TestRecord]:
        records: Dict[str, TestRecord] = {}
        segments: List[Tuple[str, int, List[TestRecord], int]] = []
        current: Optional[List] = None  # [binary, line, records, error_code]
        with ctx.step_log.open(encoding="utf-8", errors="replace") as stream:
            for number, raw in enumerate(stream, start=1):
                line = raw.rstrip("\n")
                started = echo.match(line)
                if started:
                    if current:
                        segments.append(tuple(current))
                    current = [started.group(1) if started.groups() else line,
                               number, [], 0]
                    continue
                ignored = _MAKE_IGNORED.match(line)
                if ignored and current:
                    current[3] = int(ignored.group(1))
                    continue
                match = _BRACKET.match(line) or _NAMED.match(line)
                if match:
                    record = ctx.record(match.group(1), match.group(2),
                                        match.group(3).strip(), line=number)
                    if current:
                        current[2].append(record)
                    else:
                        # Result lines before the first binary (e.g. an
                        # architecture check run as a prerequisite).
                        segments.append(("", number, [record], 0))
        if current:
            segments.append(tuple(current))
        for binary, line, found, error in segments:
            if binary and not found:
                found = [ctx.record(binary, "ERROR" if error else "PASS",
                                    f"exited with status {error}" if error
                                    else "exited 0 without per-test result lines",
                                    line=line)]
            elif binary and error and not any(r.status in ("FAIL", "ERROR") for r in found):
                found = found + [ctx.record(binary, "ERROR",
                                            f"exited with status {error}", line=line)]
            for record in found:
                previous = records.get(record.test_id)
                if previous is None or SEVERITY[record.status] > SEVERITY[previous.status]:
                    records[record.test_id] = record
        if not records:
            raise ResultsUnavailable("no test binaries or result lines in the step log")
        return list(records.values())
    return read


# ---------------------------------------------------------------------------
# Native build (AGENTS.md preflight + the user-specified sequence)
# ---------------------------------------------------------------------------

def running_processes_in_checkout() -> List[str]:
    """make/amps/mpiexec processes of this user whose cwd is inside AMPS_ROOT.

    A native rebuild removes build/ and ./amps, so it must not start while
    another build or run uses this checkout."""
    found = []
    uid = os.getuid()
    for entry in os.listdir("/proc"):
        if not entry.isdigit() or int(entry) == os.getpid():
            continue
        try:
            if os.stat(f"/proc/{entry}").st_uid != uid:
                continue
            name = Path(f"/proc/{entry}/comm").read_text().strip()
            if name not in ("make", "amps", "mpiexec", "mpirun"):
                continue
            cwd = os.readlink(f"/proc/{entry}/cwd")
        except OSError:
            continue
        if cwd == str(AMPS_ROOT) or cwd.startswith(str(AMPS_ROOT) + os.sep):
            found.append(f"pid {entry} {name} in {cwd}")
    return found


def build_preflight(log: TextIO) -> int:
    """Confirm the root, refuse while the checkout is in use, rm -rf build."""
    log.write(f"preflight: AMPS root {AMPS_ROOT}\n")
    busy = running_processes_in_checkout()
    if busy:
        log.write("preflight: refusing to remove build/ while running:\n  "
                  + "\n  ".join(busy) + "\n")
        return 2
    completed = subprocess.run(["rm", "-rf", "--", "build"], cwd=str(AMPS_ROOT),
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
    log.write(f"$ rm -rf -- build\n{completed.stdout}[exit {completed.returncode}]\n")
    return completed.returncode


def build_step(application: str) -> Step:
    return Step(
        name="build", kind="check", runner="native build",
        description=(f"clean native build: preflight, rm -rf -- build, "
                     f"./Config.pl -application={application}, make -j"),
        commands=[["./Config.pl", f"-application={application}"], ["make", "-j"]],
        pre=build_preflight)


# ---------------------------------------------------------------------------
# Application step tables
# ---------------------------------------------------------------------------

SELFTEST = "tools/test_sep_test_orchestrator.py"


def srcsep_steps(out: Path, args: argparse.Namespace) -> List[Step]:
    """srcSEP: steps from the srcSEP section of srcSEP3D/LOG.md.  New srcSEP
    test procedures must be added here (with a reader) when developed."""
    run_tests = out / "run-tests"
    pd_json = out / "pd-mover.json"
    return [
        build_step("test/sep_parker_spiral__field_line"),
        Step("run-tests",
             "srcSEP/test/run_tests.py --all against the new executable: every "
             "test in the linked registry plus the CV/IV/XM/OV/EV validation "
             "cases, each in an isolated process (does not rebuild AMPS)",
             [["python3", "srcSEP/test/run_tests.py", "--amps", "./amps", "--all",
               "--output-dir", str(run_tests)]],
             reader=reader_srcsep_run_tests(run_tests),
             runner="srcSEP/test/run_tests.py", namespace="srcsep-registry"),
        Step("pd-mover",
             "parallel-diffusion / mover runner: library PD13-*, PDB01-PDB05, "
             "CLI01-CLI06, COEF01-COEF06",
             [["python3", "srcSEP/test/run_parallel_diffusion_mover_tests.py",
               "--json", str(pd_json)]],
             reader=reader_pd_mover(pd_json,
                                    "srcSEP/test/run_parallel_diffusion_mover_tests.py",
                                    {"library": "parallel-diffusion-library"}),
             runner="srcSEP/test/run_parallel_diffusion_mover_tests.py"),
        Step("orchestrator-selftest",
             "self-test of this orchestrator's result readers and log extraction",
             [["python3", SELFTEST]], reader=read_unittest, runner=SELFTEST),
    ]


def srcsep3d_steps(out: Path, args: argparse.Namespace) -> List[Step]:
    """srcSEP3D: steps from srcSEP3D/LOG.md.  New srcSEP3D test procedures
    must be added here (with a reader) when developed."""
    ranks = str(args.ranks)
    large = str(args.ranks_large)
    stage1 = "srcsep3d-registry"   # stage-1 tests are run by several steps
    pd_json = out / "pd-mover.json"
    coronal = AMPS_ROOT / "src/models/sep_coronal_cme/build/test-results/results.json"

    def native(name: str, description: str, nranks: str, test_args: List[str],
               long: bool = False) -> Step:
        report = out / name / "native.json"
        return Step(name, description,
                    [["mpiexec", "-n", nranks, "./amps", *test_args,
                      "--expect-mpi-ranks", nranks, "--test-json", str(report),
                      "--artifact-directory", str(out / name / "artifacts")]],
                    reader=reader_native_json(report), runner="./amps (native)",
                    long=long, make_dirs=[out / name])

    def init(name: str, deck: str) -> Step:
        return Step(name, f"initialization-only run of {deck} ({large} ranks); "
                    "writes the mesh/Parker-line preview products",
                    [["mpiexec", "-n", large, "./amps", "--input", deck,
                      "--initialization-only", "--initialization-output-dir",
                      str(out / name)]], kind="check", runner="./amps (native)")

    return [
        build_step("sep3d"),
        init("init-full", "srcSEP3D/examples/sep3d_analytic_parker.in"),
        init("init-tube", "srcSEP3D/examples/sep3d_analytic_parker_active_tube.in"),
        native("sccm3d", f"native SCCM3D01-SCCM3D07, zero steps ({ranks} ranks)", ranks,
               [*sum((["--test", f"SCCM3D0{i}"] for i in range(1, 8)), []),
                "--test-input", "srcSEP3D/examples/sep3d_analytic_parker.in",
                "--test-steps", "0"]),
        native("sep-corona-suite",
               f"native --test-suite sep-corona on the active-tube deck ({ranks} ranks)",
               ranks, ["--test-suite", "sep-corona", "--test-input",
                       "srcSEP3D/examples/sep3d_analytic_parker_active_tube.in",
                       "--test-steps", "0"]),
        Step("coupled-runner",
             "run_coupled_sep_corona.py: shared SEP/corona gates plus every "
             "live native sep-corona descriptor",
             [["python3", "srcSEP3D/test/run_coupled_sep_corona.py", "--amps", "./amps",
               "--ranks", ranks, "--test-input",
               "srcSEP3D/examples/sep3d_analytic_parker_active_tube.in",
               "--test-steps", "0", "--output-dir", str(out / "coupled-runner")]],
             reader=reader_coupled(out / "coupled-runner"),
             runner="srcSEP3D/test/run_coupled_sep_corona.py"),
        native("swcme-coupled",
               f"native sep-corona suite on the SWCME SSE 20 Rs-1 AU deck, "
               f"20 steps ({large} ranks)", large,
               ["--test-suite", "sep-corona", "--test-input",
                "srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in",
                "--test-steps", "20"], long=True),
        Step("run-tests-all",
             "srcSEP3D/test/run_tests.py --all: standalone, R0-R2, C01-C05, "
             "M/B/T/P/A/O/V and production-build tests (rebuilds only test/stage1)",
             [["python3", "srcSEP3D/test/run_tests.py", "--all", "--amps-source", ".",
               "--make-config", "Makefile.conf", "--rebuild",
               "--output-dir", str(out / "run-tests-all")]],
             env={"MAKEFLAGS": "-j16"},
             reader=reader_srcsep3d_run_tests(out / "run-tests-all"),
             runner="srcSEP3D/test/run_tests.py", namespace=stage1),
        Step("pd-mover",
             "parallel-diffusion / mover runner --all: library PD13-*, CFG3D17, "
             "R3D01, BLDL3D05, full stage-1 suite",
             [["python3", "srcSEP3D/test/run_parallel_diffusion_mover_tests.py",
               "--all", "--json", str(pd_json)]],
             reader=reader_pd_mover(
                 pd_json, "srcSEP3D/test/run_parallel_diffusion_mover_tests.py",
                 {"library": "parallel-diffusion-library", "stage1": stage1,
                  "stage1-all": stage1, "bldl3d05": stage1}),
             runner="srcSEP3D/test/run_parallel_diffusion_mover_tests.py"),
        Step("reduced-front",
             "run_reduced_shock_front.py: reduced shock-front shared, portable "
             "and native tests (uses this build; no --rebuild-native)",
             [["python3", "srcSEP3D/test/run_reduced_shock_front.py", "--amps", "./amps",
               "--output-dir", str(out / "reduced-front")]],
             reader=reader_reduced_front(out / "reduced-front"),
             runner="srcSEP3D/test/run_reduced_shock_front.py", long=True),
        Step("positive-shock",
             "validate_positive_shock_example.py (~13 min): 1-/4-rank positive "
             "1-AU reduced-front qualification",
             [["python3", "srcSEP3D/test/validate_positive_shock_example.py",
               "--amps", "./amps", "--output-root", str(out / "positive-shock")]],
             reader=reader_positive_shock(out / "positive-shock"),
             runner="srcSEP3D/test/validate_positive_shock_example.py", long=True),
        Step("viewer-tests", "test_shock_front_viewer.py unit tests (no AMPS run)",
             [["python3", "srcSEP3D/test/test_shock_front_viewer.py"]],
             reader=read_unittest, runner="srcSEP3D/test/test_shock_front_viewer.py"),
        Step("corona-swcme",
             "make -i -C src/models/sep_corona_swcme test: shared Corona-SWCME "
             "model tests incl. BG3D-4 (UNQUALIFIED; its FAIL is expected "
             "evidence); -i runs every test binary of the recipe even after "
             "one fails, and the reader attributes failures per binary",
             [["make", "-i", "-C", "src/models/sep_corona_swcme", "-j16", "test"]],
             reader=reader_make_recipe(r"^build/(test_\w+)$"),
             runner="src/models/sep_corona_swcme make test",
             # test_bg3d4_piston_refinement alone runs longer than 10 minutes.
             long=True),
        Step("coronal-cme",
             "make -C src/models/sep_coronal_cme test: coronal-CME regression",
             [["make", "-C", "src/models/sep_coronal_cme", "-j16", "test"]],
             reader=reader_coronal_cme(coronal),
             runner="src/models/sep_coronal_cme make test"),
        Step("orchestrator-selftest",
             "self-test of this orchestrator's result readers and log extraction",
             [["python3", SELFTEST]], reader=read_unittest, runner=SELFTEST),
    ]


APPS = {
    "srcsep": ("srcSEP", srcsep_steps),
    "srcsep3d": ("srcSEP3D", srcsep3d_steps),
}


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

class Console:
    """Terminal output that is also appended to the run log."""

    def __init__(self, run_log: Optional[Path]):
        self.stream = None
        if run_log is not None:
            run_log.parent.mkdir(parents=True, exist_ok=True)
            self.stream = run_log.open("a", encoding="utf-8")

    def write(self, text: str) -> None:
        sys.stdout.write(text)
        sys.stdout.flush()
        if self.stream:
            self.stream.write(text)
            self.stream.flush()

    def line(self, text: str = "") -> None:
        self.write(text + "\n")


def run_command(argv: List[str], log: TextIO, console: Console, stream: bool,
                env: Dict[str, str]) -> int:
    """Run one command from the AMPS root; its stdout+stderr go to the step
    log and, with --stream, live to the console (and so the run log)."""
    log.write(f"$ {shlex.join(argv)}\n")
    log.flush()
    environment = dict(os.environ, **env)
    try:
        process = subprocess.Popen(argv, cwd=str(AMPS_ROOT), env=environment,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   universal_newlines=True, errors="replace", bufsize=1)
    except OSError as error:
        log.write(f"cannot start command: {error}\n")
        return 127
    assert process.stdout is not None
    # "with process" closes the pipe and reaps the child on every path.
    with process:
        try:
            for text in process.stdout:
                log.write(text)
                if stream:
                    console.write(text)
            return process.wait()
        except KeyboardInterrupt:
            process.terminate()
            raise


def execute(step: Step, index: int, out: Path, console: Console,
            args: argparse.Namespace) -> StepRun:
    """Run one step, then read its results (see module docstring)."""
    log_path = out / f"{index:02d}-{step.name}.log"
    console.line(f"==> [{step.name}] " + " && ".join(shlex.join(c) for c in step.commands))
    for directory in step.make_dirs:
        directory.mkdir(parents=True, exist_ok=True)
    started_wall, started = time.monotonic(), time.time()
    code = 0
    with log_path.open("w", encoding="utf-8") as log:
        if step.pre is not None:
            code = step.pre(log)
            if args.stream:
                console.write(log_path.read_text(encoding="utf-8", errors="replace"))
        for argv in step.commands:
            if code != 0:
                break
            code = run_command(argv, log, console, args.stream, step.env)
            log.write(f"[exit {code}]\n")
    seconds = time.monotonic() - started_wall
    ctx = StepContext(step, out, log_path, code, started,
                      out / "failed-tests" / step.name)
    if step.kind == "check":
        records = read_check(ctx)
    else:
        try:
            records = step.reader(ctx) if step.reader else read_check(ctx)
        except ResultsUnavailable as error:
            records = [ctx.record(f"{step.name}:results-unavailable", "ERROR",
                                  f"{error}; exit code {code}")]
        if code != 0 and not any(r.status in ("FAIL", "ERROR") for r in records):
            records.append(ctx.record(
                f"{step.name}:exit-status", "ERROR",
                f"step exited {code} but its report shows no failing test"))
    failed = any(r.status in ("FAIL", "ERROR") for r in records)
    run = StepRun(step, "FAIL" if failed else "PASS", code, seconds,
                  str(log_path), records)
    totals = count(records)
    if step.kind == "check":
        console.line(f"    {run.status} ({seconds:.0f}s)  log: {log_path}")
    else:
        console.line(f"    {run.status} ({seconds:.0f}s) tests={totals['TOTAL']} "
                     f"pass={totals['PASS']} fail={totals['FAIL']} "
                     f"skip={totals['SKIP']} error={totals['ERROR']}  log: {log_path}")
    return run


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def failure_lines(record: TestRecord) -> List[str]:
    """Plain-text entry for one failed test (used on screen and in files)."""
    location = record.log + (f":{record.line}" if record.line else "")
    lines = [f"  {record.status:<5}  {record.test_id}",
             f"         runner:  {record.runner} [step {record.step}]",
             f"         message: {record.message or '(none)'}",
             f"         log:     {location or '(none)'}"]
    for item in record.evidence[:3]:
        lines.append(f"         report:  {item}")
    return lines


def summarize(app_label: str, runs: List[StepRun], console: Console, out: Path,
              args: argparse.Namespace, started_utc: str) -> int:
    """Print and write the summary; return the process exit status."""
    tests = [r for run in runs for r in run.records if r.kind == "test"]
    checks = [r for run in runs for r in run.records if r.kind == "check"]
    executions = count(tests)
    unique = count(worst(tests).values())
    check_totals = count(checks)

    console.line()
    console.line(f"=== {app_label} run_all_test summary ===")
    console.line(f"  {'step':<22}{'result':<7}{'tests':>6}{'pass':>6}{'fail':>6}"
                 f"{'skip':>6}{'error':>7}{'time':>8}")
    for run in runs:
        if run.status in ("SKIP", "DRY"):
            reason = "--quick" if run.status == "SKIP" else "dry run"
            console.line(f"  {run.step.name:<22}{run.status:<7}  ({reason})")
            continue
        if run.step.kind == "check":
            console.line(f"  {run.step.name:<22}{run.status:<7}{'(check)':>6}"
                         f"{'':>31}{run.seconds:>7.0f}s")
            continue
        totals = count(run.records)
        console.line(f"  {run.step.name:<22}{run.status:<7}{totals['TOTAL']:>6}"
                     f"{totals['PASS']:>6}{totals['FAIL']:>6}{totals['SKIP']:>6}"
                     f"{totals['ERROR']:>7}{run.seconds:>7.0f}s")
    console.line()
    console.line(f"  test executions: total={executions['TOTAL']} "
                 f"pass={executions['PASS']} fail={executions['FAIL']} "
                 f"skip={executions['SKIP']} error={executions['ERROR']}")
    console.line(f"  unique tests:    total={unique['TOTAL']} pass={unique['PASS']} "
                 f"fail={unique['FAIL']} skip={unique['SKIP']} error={unique['ERROR']}")
    console.line(f"  checks:          total={check_totals['TOTAL']} "
                 f"pass={check_totals['PASS']} fail={check_totals['FAIL']}")

    failed = [r for r in checks + tests if r.status in ("FAIL", "ERROR")]
    text = [f"Failed tests and checks ({len(failed)}):"]
    for record in failed:
        text.extend(failure_lines(record))
    if not failed:
        text.append("  none")
    console.line()
    for line in text:
        console.line(line)
    if any(r.step == "build" and r.status != "PASS" for r in checks):
        console.line("NOTE: the native build failed; native steps used a stale or missing ./amps.")

    if not args.dry_run:
        (out / "failed_tests.txt").write_text("\n".join(text) + "\n", encoding="utf-8")
        report = {
            "schema": "sep-run-all-test-v1", "application": app_label,
            "amps_root": str(AMPS_ROOT), "started_utc": started_utc,
            "completed_utc": utc_stamp(), "output_dir": str(out),
            "counts": {"test_executions": executions, "unique_tests": unique,
                       "checks": check_totals},
            "steps": [{"name": run.step.name, "kind": run.step.kind,
                       "status": run.status, "exit_code": run.exit_code,
                       "seconds": round(run.seconds, 3), "log": run.log,
                       "commands": run.step.commands,
                       "counts": count(run.records)} for run in runs],
            "records": [asdict(r) for run in runs for r in run.records],
        }
        (out / "run_all_test.json").write_text(json.dumps(report, indent=2) + "\n",
                                               encoding="utf-8")
        console.line(f"\nreports: {out / 'run_all_test.json'}  {out / 'failed_tests.txt'}")
        if args.junit:
            write_junit(args.junit, app_label, runs)
            console.line(f"junit:   {args.junit}")
    return 1 if failed else 0


def write_junit(path: Path, app_label: str, runs: List[StepRun]) -> None:
    """Minimal JUnit XML: one testsuite per step, one testcase per record."""
    parts = [f'<?xml version="1.0" encoding="UTF-8"?>\n<testsuites name="{_xml_escape(app_label)}">']
    for run in runs:
        totals = count(run.records)
        parts.append(f'  <testsuite name="{_xml_escape(run.step.name)}" tests="{totals["TOTAL"]}" '
                     f'failures="{totals["FAIL"]}" errors="{totals["ERROR"]}" '
                     f'skipped="{totals["SKIP"]}" time="{run.seconds:.3f}">')
        for record in run.records:
            parts.append(f'    <testcase classname="{_xml_escape(record.runner)}" '
                         f'name="{_xml_escape(record.test_id)}">')
            detail = _xml_escape(f"{record.message} (log: {record.log})")
            if record.status == "FAIL":
                parts.append(f'      <failure message="{detail}"/>')
            elif record.status == "ERROR":
                parts.append(f'      <error message="{detail}"/>')
            elif record.status == "SKIP":
                parts.append("      <skipped/>")
            parts.append("    </testcase>")
        parts.append("  </testsuite>")
    parts.append("</testsuites>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def make_parser(app: Optional[str]) -> argparse.ArgumentParser:
    epilog = ""
    if app:
        label, factory = APPS[app]
        placeholder = argparse.Namespace(ranks="$SEP3D_TEST_RANKS",
                                         ranks_large="$SEP3D_TEST_RANKS_LARGE")
        steps = factory(Path("OUTPUT_DIR"), placeholder)
        epilog = f"{label} steps (in order; [long] = skipped by --quick):\n" + "".join(
            f"  {s.name:<22}{'[long] ' if s.long else ''}{'[check] ' if s.kind == 'check' else ''}"
            f"{s.description}\n" for s in steps)
        epilog += ("\nNew tests must be added to this step table "
                   "(tools/sep_test_orchestrator.py) with a result reader.\n")
    parser = argparse.ArgumentParser(
        prog=f"{APPS[app][0]}/test/run_all_test.sh" if app else None,
        description=__doc__.split("\n\n")[0], epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--app", choices=sorted(APPS), required=app is None,
                        default=app, help=argparse.SUPPRESS if app else "application")
    parser.add_argument("--skip-build", action="store_true",
                        help="reuse the existing ./amps and build/ (must match the app)")
    parser.add_argument("--quick", action="store_true", help="skip [long] steps")
    parser.add_argument("--only", action="append", default=[], metavar="STEP",
                        help="run only this step (repeatable; see --list)")
    parser.add_argument("--list", action="store_true",
                        help="list the steps and their commands, then exit")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the commands without running anything")
    parser.add_argument("--output-dir", type=Path,
                        help="logs and products (default: "
                             "test_output/run_all_test/<app>-<UTC time>)")
    parser.add_argument("--stream", action="store_true",
                        help="also show each step's complete output live")
    parser.add_argument("--log", type=Path, metavar="FILE",
                        help="run log of all terminal output "
                             "(default: OUTPUT_DIR/run_all_test.log)")
    parser.add_argument("--junit", type=Path, metavar="FILE",
                        help="also write JUnit XML")
    parser.add_argument("--ranks", type=int,
                        default=int(os.environ.get("SEP3D_TEST_RANKS", "4")),
                        help="MPI ranks for 4-rank cases (env SEP3D_TEST_RANKS)")
    parser.add_argument("--ranks-large", type=int,
                        default=int(os.environ.get("SEP3D_TEST_RANKS_LARGE", "10")),
                        help="MPI ranks for large cases (env SEP3D_TEST_RANKS_LARGE)")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Pre-scan --app so --help can describe that application's steps.
    app = None
    for index, item in enumerate(argv):
        if item.startswith("--app="):
            app = item.split("=", 1)[1]
        elif item == "--app" and index + 1 < len(argv):
            app = argv[index + 1]
    args = make_parser(app if app in APPS else None).parse_args(argv)
    label, factory = APPS[args.app]
    started_utc = utc_stamp()
    out = (args.output_dir or AMPS_ROOT / "test_output" / "run_all_test"
           / f"{args.app}-{started_utc}").expanduser().resolve()
    steps = factory(out, args)
    known = {step.name for step in steps}
    unknown = [name for name in args.only if name not in known]
    if unknown:
        print(f"unknown step(s): {', '.join(unknown)}; known: {', '.join(sorted(known))}",
              file=sys.stderr)
        return 2
    if args.list:
        for step in steps:
            flags = ("[long] " if step.long else "") + ("[check] " if step.kind == "check" else "")
            print(f"{step.name:<22}{flags}" + " && ".join(shlex.join(c) for c in step.commands))
        return 0

    selected = [s for s in steps if not args.only or s.name in args.only]
    if args.skip_build:
        selected = [s for s in selected if s.name != "build"]
    console = Console(None if args.dry_run else
                      (args.log or out / "run_all_test.log").expanduser().resolve())
    if not args.dry_run:
        out.mkdir(parents=True, exist_ok=True)
    console.line(f"{label} run_all_test: AMPS root {AMPS_ROOT}")
    console.line(f"output: {out}")
    runs: List[StepRun] = []
    index = 0
    for step in selected:
        if step.long and args.quick:
            runs.append(StepRun(step, "SKIP", None, 0.0, ""))
            continue
        if args.dry_run:
            console.line(f"==> [{step.name}] "
                         + " && ".join(shlex.join(c) for c in step.commands))
            runs.append(StepRun(step, "DRY", None, 0.0, ""))
            continue
        index += 1
        runs.append(execute(step, index, out, console, args))
    return summarize(label, runs, console, out, args, started_utc)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        sys.exit(130)
