import datetime as dt
import io
import json
import sqlite3
import tarfile
import tempfile
import unittest
from pathlib import Path

from amps_nightly import report
from amps_nightly.records import (RecordError, make_run_id, normalize_result, normalize_run, read_run_dir,
                                  write_run_dir)
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

    def test_archive_index_and_night_navigation(self):
        from amps_nightly.store import write_index
        now = dt.datetime(2026, 11, 13, 12, 0, tzinfo=UTC)
        base = self.www / "archive"
        first, second = base / "2026/11/12/passion-gnu", base / "2026/11/13/passion-gnu"
        with Store(self.root / "results.db") as st:
            ingest(self.root, self.www, st, now)
            self.assertEqual(report.build(st, Registry.load(), self.www, now=now, tz=UTC)["night_pages_updated"], 3)
            page = (first / "index.html").read_text()
            self.assertIn('rel="next" href="../../13/passion-gnu/index.html"', page)
            self.assertIn("Next night · 2026-11-13", page)
            self.assertIn("No previous night", page)
            page = (second / "index.html").read_text()
            self.assertIn('rel="prev" href="../../12/passion-gnu/index.html"', page)
            self.assertIn("No next night", page)
            self.assertIn('href="../../../../index.html">All nights', page)
            other = (base / "2026/11/01/pleiades-intel/index.html").read_text()     # columns are separate
            self.assertIn("No previous night", other)
            self.assertIn("No next night", other)
            # nothing new: no night page is rewritten
            self.assertEqual(report.build(st, Registry.load(), self.www, now=now, tz=UTC)["night_pages_updated"], 0)
            # a new night: it links back and the night before gains its next link; older nights stay as they are
            before = (first / "index.html").read_text()
            write_run_dir(self.root / "incoming" / "d", run("passion", "gnu", "2026-11-14T11:30:00Z", SHA[2]),
                          [res("pic-field-solver", "t1", "PASS")])
            ingest(self.root, self.www, st, now)
            self.assertEqual(report.build(st, Registry.load(), self.www, now=now, tz=UTC)["night_pages_updated"], 2)
            self.assertIn('rel="next" href="../../14/passion-gnu/index.html"', (second / "index.html").read_text())
            self.assertIn('rel="prev" href="../../13/passion-gnu/index.html"',
                          (base / "2026/11/14/passion-gnu/index.html").read_text())
            self.assertEqual((first / "index.html").read_text(), before)
        # any other writer of a night page (ingest, prune, relocate) keeps its links
        write_index(second)
        self.assertIn('rel="next" href="../../14/passion-gnu/index.html"', (second / "index.html").read_text())

        archive = (base / "index.html").read_text()
        self.assertIn("4 stored runs in 2 columns", archive)
        self.assertLess(archive.index("2026-11-14"), archive.index("2026-11-13"))     # newest first
        self.assertIn('href="2026/11/13/passion-gnu/index.html"', archive)
        self.assertIn('href="2026/11/01/pleiades-intel/index.html"', archive)          # every stored column
        self.assertLess(archive.index('id="passion-gnu"'), archive.index('id="pleiades-intel"'))   # latest run first
        # night of 2026-11-13: t2 and CORE-PIC-NATIVE-001 failed, s1 skipped, t1 passed (result rows)
        row = archive[archive.index('href="2026/11/13/passion-gnu/index.html"'):]
        self.assertIn('<td class="num">2</td><td class="num">1</td><td class="num">1</td>', row[:row.index("</tr>")])
        self.assertIn('href="archive/index.html">All nights', (self.www / "index.html").read_text())
        self.assertIn('href="../../archive/index.html">All nights',
                      (self.www / "family" / "pic-field-solver" / "passion-gnu.html").read_text())

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


# the runs table before routines existed (results.db of November 2026)
OLD_RUNS = """CREATE TABLE runs (run_id TEXT PRIMARY KEY, column_id TEXT NOT NULL, machine TEXT NOT NULL,
    config TEXT NOT NULL, compiler TEXT NOT NULL DEFAULT '', toolchain TEXT NOT NULL DEFAULT '',
    branch TEXT NOT NULL DEFAULT '', commit_sha TEXT NOT NULL DEFAULT '', commit_inferred INTEGER NOT NULL DEFAULT 0,
    started TEXT NOT NULL DEFAULT '', finished TEXT NOT NULL, complete INTEGER NOT NULL,
    producer TEXT NOT NULL DEFAULT '', notes TEXT NOT NULL DEFAULT '', path TEXT NOT NULL DEFAULT '',
    ingested TEXT NOT NULL);
INSERT INTO runs VALUES ('20261113T1100Z-passion-gnu', 'passion-gnu', 'passion', 'gnu', 'GNU', '', 'master', '', 0,
    '', '2026-11-13T11:00:00Z', 1, 'amps-nightly run', '', 'archive/2026/11/13/passion-gnu', '2026-11-13T12:00:00Z');
"""


def long_run(finished, commit, night):
    raw = run("passion", "gnu", finished, commit)
    raw.update(routine="long", night=night, run_id=raw["run_id"] + "-long")
    return raw


class RoutineRecords(unittest.TestCase):
    def test_run_id_and_record(self):
        when = dt.datetime(2026, 11, 28, 16, 5, tzinfo=UTC)
        self.assertEqual(make_run_id("passion", "gnu", when), "20261128T1605Z-passion-gnu")
        self.assertEqual(make_run_id("passion", "gnu", when, "long"), "20261128T1605Z-passion-gnu-long")
        old = run("passion", "gnu", "2026-11-13T11:00:00Z", SHA[0])          # a record from before routines
        self.assertEqual(normalize_run(old)["routine"], "nightly")
        self.assertEqual(normalize_run(dict(old, routine="long"))["routine"], "long")
        for bad in ("Long", "long run", 7):
            with self.assertRaises(RecordError):
                normalize_run(dict(old, routine=bad))

    def test_old_store_and_old_tarball(self):
        """A results.db and a tarball made before routines: every run in them is nightly."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            db = sqlite3.connect(str(root / "results.db"))
            db.executescript(OLD_RUNS)
            db.commit()
            db.close()
            # run.json as write_run_dir wrote it before routines: every field but "routine"
            raw = {k: v for k, v in normalize_run(run("passion", "gnu", "2026-11-14T11:00:00Z", SHA[1])).items()
                   if k not in ("routine", "column")}
            (root / "incoming").mkdir()
            with tarfile.open(str(root / "incoming" / (raw["run_id"] + ".tar.gz")), "w:gz") as tar:
                for name, text in (("run.json", json.dumps(raw)),
                                   ("results.jsonl", json.dumps(res("pic-field-solver", "t1", "PASS")) + "\n"),
                                   ("DONE", raw["run_id"] + "\n")):
                    info = tarfile.TarInfo("{}/{}".format(raw["run_id"], name))
                    info.size = len(text.encode())
                    tar.addfile(info, io.BytesIO(text.encode()))
            with Store(root / "results.db") as st:
                accepted, rejected = ingest(root, root / "www", st)
                self.assertEqual((accepted, rejected), ([raw["run_id"]], []))
                self.assertEqual([r["routine"] for r in st.latest_runs("passion-gnu", 5)], ["nightly", "nightly"])
                self.assertEqual(st.routines(), [])
                self.assertEqual(st.columns()[0]["runs"], 2)
            self.assertTrue((root / "www/archive/2026/11/14/passion-gnu/DONE").is_file())


class RoutineReport(unittest.TestCase):
    """A long routine next to the nightly: separate pages, each test's newest long result, its own staleness."""

    DECLARED = [{"column": "passion-gnu", "machine": "passion", "config": "gnu", "compiler": "GNU", "toolchain": ""}]

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.www = self.root / "www"
        incoming = self.root / "incoming"
        write_run_dir(incoming / "n1", run("passion", "gnu", "2026-11-27T11:00:00Z", SHA[0]),
                      [res("pic-field-solver", "t1", "PASS"), res("pic-field-solver", "t2", "PASS")])
        write_run_dir(incoming / "n2", run("passion", "gnu", "2026-11-28T11:00:00Z", SHA[1]),
                      [res("pic-field-solver", "t1", "PASS"), res("pic-field-solver", "t2", "FAILED", reason="drift")])
        write_run_dir(incoming / "l1", long_run("2026-11-21T22:00:00Z", SHA[0], "2026-11-21"),
                      [res("sep-long", "L1", "PASS"), res("sep-long", "L2", "FAILED", reason="differs"),
                       res("sep-long", "L3", "PASS")])
        # the newest long run: newer than every nightly run, and only L1 was in it (a run of one test sent by hand)
        write_run_dir(incoming / "l2", long_run("2026-11-29T20:00:00Z", SHA[2], "2026-11-29"),
                      [res("sep-long", "L1", "FAILED", reason="flux differs", log="logs/legacy-long/tests/L1.log")])
        (incoming / "l2" / "logs" / "legacy-long" / "tests").mkdir(parents=True)
        (incoming / "l2" / "logs" / "legacy-long" / "tests" / "L1.log").write_text("log of L1\n")
        with Store(self.root / "results.db") as st:
            accepted, rejected = ingest(self.root, self.www, st, dt.datetime(2026, 11, 29, 21, 0, tzinfo=UTC))
        self.assertEqual((len(accepted), rejected), (4, []))

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, now, stale_after_h=192, machines=True):
        with Store(self.root / "results.db") as st:
            return report.build(st, Registry.load(), self.www, now=now, tz=UTC,
                                declared=self.DECLARED if machines else None,
                                routines={"long": [dict(self.DECLARED[0], stale_after_h=stale_after_h)]}
                                if machines else None)

    def test_nightly_and_long_results_are_separate(self):
        summary = self.build(dt.datetime(2026, 11, 29, 21, 0, tzinfo=UTC))
        self.assertEqual((summary["columns"], summary["families"]), (["passion-gnu"], ["pic-field-solver"]))
        self.assertEqual(summary["routines"], {"long": {"columns": ["passion-gnu"], "waiting": [],
                                                        "families": ["sep-long"]}})
        self.assertTrue((self.www / "archive/2026/11/29/passion-gnu-long/DONE").is_file())
        self.assertIn("Routine long, started on 2026-11-29",
                      (self.www / "archive/2026/11/29/passion-gnu-long/index.html").read_text())
        index = (self.www / "index.html").read_text()
        nightly, long_ = index.split('<section class="rsec"')
        # the nightly table: the newest nightly run, although a long run finished after it
        self.assertIn("at " + report.commit_html(SHA[1], False), nightly)
        self.assertNotIn("sep-long", nightly)
        self.assertIn("NEW FAILURE", nightly)                                    # t2, nightly 1 -> nightly 2
        page = (self.www / "family/pic-field-solver/passion-gnu.html").read_text()
        self.assertIn("drift", page)
        self.assertNotIn(".long.html", page)
        # the long section: L1 from the newest long run, L2 and L3 from the one before
        self.assertIn("<h2>Long validation runs (periodic)</h2>", long_)
        self.assertIn('href="family/sep-long/passion-gnu.long.html"', long_)
        self.assertIn(report.tally_html(1, 0, 2), long_)
        self.assertNotIn("No upload", long_)
        page = (self.www / "family/sep-long/passion-gnu.long.html").read_text()
        self.assertIn("Last run and its commit", page)
        rows = {name: page[page.index(">{}<".format(name)):] for name in ("L1", "L2", "L3")}
        self.assertIn("29 Nov 2026", rows["L1"][:rows["L1"].index("</div></div>")])
        self.assertIn(SHA[2][:7], rows["L1"].split('class="trow')[0])
        self.assertIn("8 days before this run", rows["L1"].split('class="trow')[0])      # passed in the first long run
        self.assertIn('href="../../archive/2026/11/29/passion-gnu-long/logs/legacy-long/tests/L1.log"', page)
        self.assertIn("21 Nov 2026", rows["L3"].split('class="trow')[0])
        self.assertIn(SHA[0][:7], rows["L3"].split('class="trow')[0])
        self.assertIn("this run", rows["L3"].split('class="trow')[0])
        self.assertFalse((self.www / "family/sep-long/passion-gnu.html").exists())
        # previous/next links stay within a routine; the archive index has one table per column and routine
        nav = json.loads((self.www / "archive/2026/11/29/passion-gnu-long/nav.json").read_text())
        self.assertEqual((nav["previous"]["href"], nav["previous"]["night"], nav["next"]),
                         ("../../21/passion-gnu-long/index.html", "2026-11-21", None))
        nav = json.loads((self.www / "archive/2026/11/28/passion-gnu/nav.json").read_text())
        self.assertEqual((nav["previous"]["night"], nav["next"]), ("2026-11-27", None))
        archive = (self.www / "archive/index.html").read_text()
        self.assertIn('id="passion-gnu-long"', archive)
        self.assertIn("GNU · gcc 13 · routine long", archive)

    def test_staleness_of_each_routine(self):
        # 4 days after the last nightly run, 3 days after the last long run: only the nightly column is flagged
        self.build(dt.datetime(2026, 12, 2, 21, 0, tzinfo=UTC))
        nightly, long_ = (self.www / "index.html").read_text().split('<section class="rsec"')
        self.assertIn("No upload for 4 days", nightly)
        self.assertNotIn("No upload", long_)
        self.assertIn("older than 192 h is flagged", long_)
        # 9 days after the last long run: flagged after its stale_after_h (192 h), or a shorter one
        self.build(dt.datetime(2026, 12, 8, 21, 0, tzinfo=UTC))
        self.assertIn("No upload for 9 days", (self.www / "index.html").read_text().split('<section class="rsec"')[1])
        self.build(dt.datetime(2026, 12, 2, 21, 0, tzinfo=UTC), stale_after_h=48)
        self.assertIn("No upload for 3 days", (self.www / "index.html").read_text().split('<section class="rsec"')[1])

    def test_without_machine_files(self):
        summary = self.build(dt.datetime(2026, 11, 29, 21, 0, tzinfo=UTC), machines=False)
        self.assertEqual(summary["routines"]["long"]["columns"], ["passion-gnu"])
        self.assertTrue((self.www / "family/sep-long/passion-gnu.long.html").is_file())

    def test_declared_routine_without_results(self):
        with Store(self.root / "results.db") as st:
            summary = report.build(st, Registry.load(), self.www, now=dt.datetime(2026, 11, 29, 21, 0, tzinfo=UTC),
                                   tz=UTC, declared=self.DECLARED,
                                   routines={"weekly": [dict(self.DECLARED[0], stale_after_h=192)]})
        self.assertEqual(summary["routines"]["weekly"], {"columns": [], "waiting": ["passion-gnu"], "families": []})
        self.assertIn("No results yet from passion · GNU.", (self.www / "index.html").read_text())


if __name__ == "__main__":
    unittest.main()
