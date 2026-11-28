#!/usr/bin/env python3
"""Run standalone plume-model unittest modules and record one outcome per test.

``run_tests.py`` starts this script once per group of ``standalone/tests``
modules (in a fresh interpreter, working directory ``srcEnceladus``, exactly as
``python3 -m unittest standalone.tests.<module>`` would run them) and turns the
JSON it writes into one receipt row per test.

Usage (from ``srcEnceladus``)::

    python3 test/standalone_test_receipt.py --output outcomes.json test_physics [test_fitting ...]

Output: ``{"collected": N, "tests": [{"name", "id", "module", "status",
"detail", "duration_s"}, ...]}``, in loader order (module order as given,
then classes and methods sorted by name).

* ``name`` is ``S_<module without test_>.<Class>.<method>``. Removing the
  ``S_`` and putting back ``test_`` gives the unittest selector, e.g.
  ``S_physics.PhysicsTest.test_x`` is
  ``python3 -m unittest standalone.tests.test_physics.PhysicsTest.test_x``.
* ``status`` is ``PASS``, ``FAIL`` or ``SKIP``. A unittest error counts as
  ``FAIL``; an expected failure counts as ``PASS``, an unexpected success as
  ``FAIL``; any failing subtest fails its test.
* ``detail`` starts with one summary line (``PASS (…s)``, ``SKIP: <reason>``,
  ``FAIL: <exception>`` or ``ERROR: <exception>``), which the AMPS nightly
  shows as the reason, followed for failures by the traceback and the
  test's captured output.

Nothing is hidden. Every collected test gets a row. A test that never ran
because a ``setUpClass``/``setUpModule`` fixture was skipped or failed
inherits that fixture's outcome. A fixture problem that explains no
unrecorded test (for example a failing ``tearDownClass``) gets its own
``FAIL`` row (``S_<module>.<Class>.tearDownClass``). A module that cannot be
imported, or contains no tests, is one ``FAIL`` row ``S_<module>``. The exit
status is 1 when any row is ``FAIL``, else 0.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = "standalone.tests."
FAILED_IMPORT_PREFIX = "unittest.loader._FailedTest."
# _ErrorHolder description of a class/module fixture event: "setUpClass (standalone.tests.test_x.SomeTest)"
FIXTURE_EVENT = re.compile(r"^(?P<method>\w+) \((?P<parent>[\w.]+)\)$")


def row_name(test_id: str) -> str:
    """standalone.tests.test_physics.PhysicsTest.test_x -> S_physics.PhysicsTest.test_x."""

    if not test_id.startswith(PACKAGE + "test_"):
        raise ValueError("not a standalone/tests test id: " + test_id)
    return "S_" + test_id[len(PACKAGE + "test_"):]


def module_of(test_id: str) -> str:
    return test_id[len(PACKAGE):].split(".", 1)[0]


def iter_tests(suite: unittest.TestSuite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from iter_tests(item)
        else:
            yield item


def first_line(text: str) -> str:
    return next((line.strip() for line in text.splitlines() if line.strip()), "")


class RecordingResult(unittest.TestResult):
    """TestResult that keeps one PASS/FAIL/SKIP outcome per test id."""

    def __init__(self) -> None:
        super().__init__()
        # Capture each test's stdout/stderr; TestResult appends it to the
        # formatted traceback of a failing test, so passing tests stay quiet.
        self.buffer = True
        self.outcomes: dict[str, dict] = {}
        self.fixture_events: list[dict] = []
        self._started: dict[str, float] = {}

    # -- bookkeeping -------------------------------------------------------
    def startTest(self, test) -> None:
        super().startTest(test)
        self._started[test.id()] = time.perf_counter()

    def _duration(self, test_id: str):
        started = self._started.get(test_id)
        return None if started is None else round(time.perf_counter() - started, 6)

    def _set(self, test, status: str, detail: str) -> None:
        test_id = test.id()
        event = FIXTURE_EVENT.match(test_id) if isinstance(test, unittest.suite._ErrorHolder) else None
        if event:
            self.fixture_events.append({"method": event.group("method"), "parent": event.group("parent"),
                                        "status": status, "detail": detail})
            return
        previous = self.outcomes.get(test_id)
        if previous is not None and previous["status"] == "FAIL":
            # A test with several failing subtests (or a failure plus a
            # cleanup error) keeps every message.
            previous["detail"] += "\n" + detail
            return
        self.outcomes[test_id] = {"status": status, "detail": detail, "duration_s": self._duration(test_id)}

    def _failure(self, kind: str, test, err) -> str:
        # Summary from the exception itself (type and first message line):
        # a multi-line assertion diff must not become the one-line reason.
        exc_type, exc_value = err[0], err[1]
        message = first_line(str(exc_value))
        summary = exc_type.__name__ + (": " + message if message else "")
        return "{}: {}\n{}".format(kind, summary, self._exc_info_to_string(err, test))

    # -- unittest outcome hooks -------------------------------------------
    def addSuccess(self, test) -> None:
        super().addSuccess(test)
        self._set(test, "PASS", "PASS ({:.3f} s)".format(self._duration(test.id()) or 0.0))

    def addFailure(self, test, err) -> None:
        super().addFailure(test, err)
        self._set(test, "FAIL", self._failure("FAIL", test, err))

    def addError(self, test, err) -> None:
        super().addError(test, err)
        self._set(test, "FAIL", self._failure("ERROR", test, err))

    def addSkip(self, test, reason) -> None:
        super().addSkip(test, reason)
        self._set(test, "SKIP", "SKIP: {}".format(reason))

    def addExpectedFailure(self, test, err) -> None:
        super().addExpectedFailure(test, err)
        self._set(test, "PASS", "PASS (expected failure)\n" + self._exc_info_to_string(err, test))

    def addUnexpectedSuccess(self, test) -> None:
        super().addUnexpectedSuccess(test)
        self._set(test, "FAIL", "FAIL: unexpected success of a test marked expectedFailure")

    def addSubTest(self, test, subtest, err) -> None:
        super().addSubTest(test, subtest, err)
        if err is not None:
            kind = "FAIL" if issubclass(err[0], test.failureException) else "ERROR"
            detail = self._failure(kind, subtest, err)
            # Record on the parent test: subtests are not separate rows.
            self._set(test, "FAIL", "{} [subtest {}]\n{}".format(
                detail.splitlines()[0], subtest._subDescription(), detail))


def collect(modules: list[str]):
    """Load the modules; return (tests, failed-import rows)."""

    loader = unittest.TestLoader()
    tests, failed = [], []
    for module in modules:
        loaded = list(iter_tests(loader.loadTestsFromName(PACKAGE + module)))
        broken = [test for test in loaded if test.id().startswith(FAILED_IMPORT_PREFIX)]
        if broken or loader.errors:
            # The loader keeps the import traceback in loader.errors and in
            # the placeholder test's _exception.
            message = "\n".join(loader.errors) or "\n".join(
                repr(getattr(test, "_exception", "unknown import error")) for test in broken)
            failed.append({"name": "S_" + module[len("test_"):], "id": PACKAGE + module, "module": module,
                           "status": "FAIL", "detail": "ERROR: cannot import {}\n{}".format(module, message),
                           "duration_s": None})
            loader.errors.clear()
            continue
        if not loaded:
            failed.append({"name": "S_" + module[len("test_"):], "id": PACKAGE + module, "module": module,
                           "status": "FAIL", "detail": "FAIL: no tests were collected from " + module,
                           "duration_s": None})
            continue
        tests.extend(loaded)
    return tests, failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("--output", type=Path, required=True, help="JSON file to write")
    parser.add_argument("--root", type=Path, default=ROOT,
                        help="directory holding the standalone/tests package (default: srcEnceladus; for self-tests)")
    parser.add_argument("modules", nargs="+", help="standalone/tests module names (test_<name>)")
    args = parser.parse_args(argv)
    for module in args.modules:
        if not re.fullmatch(r"test_\w+", module):
            parser.error("not a standalone/tests module name: " + module)

    # Match "python3 -m unittest" started in srcEnceladus: the package root
    # (not this script's directory) is the first import location.
    sys.path[0] = str(args.root.resolve())
    tests, rows = collect(args.modules)
    result = RecordingResult()
    unittest.TestSuite(tests).run(result)

    for test in tests:
        test_id = test.id()
        outcome = result.outcomes.get(test_id)
        if outcome is None:
            # Not run: inherit the class/module fixture event that covers it.
            event = next((event for event in result.fixture_events
                          if test_id.startswith(event["parent"] + ".")), None)
            if event is not None:
                event["explains"] = True
                outcome = {"status": event["status"], "duration_s": None,
                           "detail": "{} [{} of {}]\n{}".format(first_line(event["detail"]), event["method"],
                                                               event["parent"], event["detail"])}
            else:
                outcome = {"status": "FAIL", "duration_s": None,
                           "detail": "FAIL: the test was collected but never run"}
        rows.append({"name": row_name(test_id), "id": test_id, "module": module_of(test_id), **outcome})
    for event in result.fixture_events:
        if not event.get("explains"):
            name = "S_" + event["parent"][len(PACKAGE + "test_"):] + "." + event["method"] \
                if event["parent"].startswith(PACKAGE + "test_") else "S_fixture." + event["method"]
            rows.append({"name": name, "id": event["parent"] + "." + event["method"],
                         "module": module_of(event["parent"]) if event["parent"].startswith(PACKAGE) else "",
                         "status": "FAIL" if event["status"] == "FAIL" else event["status"],
                         "detail": event["detail"], "duration_s": None})

    names = [row["name"] for row in rows]
    if len(names) != len(set(names)):
        raise SystemExit("duplicate standalone row names: {}".format(sorted({n for n in names if names.count(n) > 1})))
    args.output.write_text(json.dumps({"collected": len(tests), "modules": args.modules, "tests": rows},
                                      indent=2) + "\n", encoding="utf-8")
    return 1 if any(row["status"] == "FAIL" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
