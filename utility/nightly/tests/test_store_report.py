import datetime as dt
import json
import tempfile
import unittest
from pathlib import Path

from amps_nightly import report
from amps_nightly.records import RecordError, normalize_result, read_run_dir, write_run_dir
from amps_nightly.registry import Registry
from amps_nightly.store import Store, ingest

UTC = dt.timezone.utc
SHA = ["1111111aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "2222222bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", "3333333ccccccccccccccccccccccccccccccccc"]


def run(machine, config, finished, commit, inferred=False):
    return {"schema_version": 1, "run_id": "{}-{}-{}".format(finished.replace(":", ""), machine, config),
            "machine": machine, "config": config, "compiler": config.upper(), "toolchain": "gcc 13",
            "branch": "master", "commit": commit, "commit_inferred": inferred, "finished": finished, "complete": True}


def res(family, test, status, **extra):
    row = {"family": family, "test": test, "status": status}
    row.update(extra)
    return row


class Records(unittest.TestCase):
    def test_round_trip_and_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_run_dir(Path(tmp, "r"), run("m", "gnu", "2026-11-13T10:00:00Z", SHA[0]),
                                 [res("fam", "t1", "PASS"), res("fam", "t2", "FAILED", reason="x", stage="run")])
            meta, rows = read_run_dir(path)
            self.assertTrue((path / "DONE").is_file())
        self.assertEqual(meta["column"], "m-gnu")
        self.assertEqual([r["status"] for r in rows], ["PASS", "FAILED"])
        with self.assertRaises(RecordError):
            normalize_result(res("fam", "t", "BROKEN"))
        with self.assertRaises(RecordError):
            normalize_result(res("fam", "t", "FAILED", log="../../etc/passwd"))


class StoreAndReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.www = self.root / "www"
        incoming = self.root / "incoming"
        # passion-gnu: night 1 all pass, night 2 t2 fails (new failure) and s1 skipped
        write_run_dir(incoming / "a", run("passion", "gnu", "2026-11-12T11:00:00Z", SHA[0], inferred=True),
                      [res("pic-field-solver", "t1", "PASS"), res("pic-field-solver", "t2", "PASS"),
                       res("core-pic", "CORE-PIC-NATIVE-001", "PASS", cases_total=95)])
        write_run_dir(incoming / "b", run("passion", "gnu", "2026-11-13T11:00:00Z", SHA[1]),
                      [res("pic-field-solver", "t1", "PASS"),
                       res("pic-field-solver", "t2", "FAILED", stage="check", reason="3 lines differ",
                           excerpt="==> test_t2.diff <==\n< 1"),
                       res("core-pic", "CORE-PIC-NATIVE-001", "FAILED", cases_total=95, cases_failed=1),
                       res("pic-field-solver", "s1", "SKIPPED", reason="no AVX-512")])
        # pleiades-intel: one old run (stale)
        write_run_dir(incoming / "c", run("pleiades", "intel", "2026-11-01T11:00:00Z", SHA[2]),
                      [res("pic-field-solver", "t1", "PASS")])
        # malformed upload
        bad = incoming / "bad"
        bad.mkdir()
        (bad / "run.json").write_text("{}")
        (bad / "results.jsonl").write_text("")
        (bad / "DONE").write_text("x")
        # unfinished upload (no DONE) must be left alone
        write_run_dir(incoming / "partial", run("passion", "gnu", "2026-11-14T11:00:00Z", SHA[2]), [], done=False)

    def tearDown(self):
        self.tmp.cleanup()

    def test_ingest_and_report(self):
        now = dt.datetime(2026, 11, 13, 12, 0, tzinfo=UTC)
        with Store(self.root / "results.db") as st:
            accepted, rejected = ingest(self.root, self.www, st, now)
            self.assertEqual(len(accepted), 3)
            self.assertEqual(len(rejected), 1)
            self.assertTrue((self.root / "incoming" / "partial").is_dir())
            self.assertTrue(any((self.root / "rejected").glob("bad-*/REJECTED.txt")))
            self.assertEqual(sorted(str(p.parent.relative_to(self.www)) for p in self.www.glob("archive/*/*/*/*/DONE")),
                             ["archive/2026/11/01/pleiades-intel", "archive/2026/11/12/passion-gnu",
                              "archive/2026/11/13/passion-gnu"])
            self.assertTrue((self.www / "archive/2026/11/13/passion-gnu/index.html").is_file())   # night summary
            # same run uploaded again is rejected, not duplicated
            write_run_dir(self.root / "incoming" / "again", run("passion", "gnu", "2026-11-13T11:00:00Z", SHA[1]), [])
            accepted, rejected = ingest(self.root, self.www, st, now)
            self.assertEqual((accepted, len(rejected)), ([], 1))

            summary = report.build(st, Registry.load(), self.www, now=now, tz=UTC)
        self.assertEqual(summary["columns"], ["passion-gnu", "pleiades-intel"])
        self.assertEqual(summary["families"], ["pic-field-solver", "core-pic"])
        index = (self.www / "index.html").read_text()
        self.assertIn("NEW FAILURE", index)
        self.assertIn("No upload for 12 days", index)            # pleiades column is stale
        self.assertIn(report.tally_html(94, 0, 1), index)        # suite cases are counted individually
        self.assertIn(report.tally_html(1, 1, 1), index)
        self.assertIn('<span class="n n-fail">1</span><span class="sl" aria-hidden="true">/</span>'
                      '<span class="n n-skip zero">0</span><span class="sl" aria-hidden="true">/</span>'
                      '<span class="n n-pass">94</span>', index)   # #fail/#skip/#pass, a 0 faded
        self.assertIn("1 failed, 0 skipped, 94 passed", index)   # in words for the tooltip and screen readers
        self.assertIn('href="family/pic-field-solver/passion-gnu.html"', index)

        page = (self.www / "family" / "pic-field-solver" / "passion-gnu.html").read_text()
        self.assertLess(page.index(">t2<"), page.index(">s1<"))   # failures first, then skips, then passes
        self.assertLess(page.index(">s1<"), page.index(">t1<"))
        self.assertIn("1 day before this run", page)              # t2 last passed the night before
        self.assertIn(SHA[0], page)                                # ... at the first night's commit
        self.assertIn("Inferred from the checkout time", page)    # first night's commit was inferred
        self.assertIn("3 lines differ", page)
        self.assertNotIn("test_t2.diff", page)                     # .diff contents are never stored or shown
        with Store(self.root / "results.db") as st:
            self.assertNotIn("excerpt", st.results(st.latest_runs("passion-gnu", 1)[0]["run_id"])[0])
        self.assertIn("never", page)                               # s1 never passed
        self.assertTrue((self.www / "family" / "core-pic" / "pleiades-intel.html").is_file())
        self.assertTrue((self.www / "assets" / "report.css").is_file())

    def test_prune_compresses_then_deletes_logs(self):
        from amps_nightly.store import prune
        for day, name in (("2026/11/10", "a"), ("2026/10/01", "b"), ("2025/01/01", "c")):
            night = self.www / "archive" / day / "passion-gnu"
            write_run_dir(night, run("passion", "gnu", day.replace("/", "-") + "T11:00:00Z", SHA[0]),
                          [res("pic-field-solver", "t1", "FAILED", log="logs/core/runner.log")])
            (night / "logs" / "core").mkdir(parents=True)
            (night / "logs" / "core" / "runner.log").write_text("log " + name)
        self.assertEqual(prune(self.www, today=dt.date(2026, 11, 14)), (1, 1))
        base = self.www / "archive"
        self.assertTrue((base / "2026/11/10/passion-gnu/logs/core/runner.log").is_file())
        self.assertTrue((base / "2026/10/01/passion-gnu/logs/core/runner.log.gz").is_file())
        self.assertIn('href="logs/core/runner.log.gz"', (base / "2026/10/01/passion-gnu/index.html").read_text())
        self.assertFalse((base / "2025/01/01/passion-gnu/logs").exists())
        self.assertTrue((base / "2025/01/01/passion-gnu/results.jsonl").is_file())

    def test_machine_files_decide_the_columns(self):
        now = dt.datetime(2026, 11, 13, 12, 0, tzinfo=UTC)
        declared = [{"column": "newbox-gnu", "machine": "newbox", "config": "gnu", "compiler": "GNU", "toolchain": "GCC 13"},
                    {"column": "passion-gnu", "machine": "passion", "config": "gnu", "compiler": "GNU", "toolchain": ""}]
        with Store(self.root / "results.db") as st:
            ingest(self.root, self.www, st, now)
            summary = report.build(st, Registry.load(), self.www, now=now, tz=UTC, declared=declared)
            self.assertEqual(len(st.columns()), 2)                  # the hidden column's history is kept
        self.assertEqual((summary["columns"], summary["waiting"], summary["hidden"]),
                         (["passion-gnu"], ["newbox-gnu"], ["pleiades-intel"]))
        index = (self.www / "index.html").read_text()
        self.assertLess(index.index("newbox"), index.index("passion"))   # machine-file order
        self.assertIn("— no results yet", index)
        self.assertIn("newbox · GNU: no results yet", index)
        self.assertIn("no enabled machine file declares them: pleiades-intel", index)
        self.assertFalse((self.www / "family" / "pic-field-solver" / "newbox-gnu.html").exists())
        self.assertFalse((self.www / "family" / "pic-field-solver" / "pleiades-intel.html").exists())
        page = (self.www / "family" / "pic-field-solver" / "passion-gnu.html").read_text()
        self.assertNotIn("pleiades", page)


if __name__ == "__main__":
    unittest.main()
