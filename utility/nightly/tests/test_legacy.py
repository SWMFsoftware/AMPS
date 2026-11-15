import collections
import datetime as dt
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from amps_nightly import legacy
from amps_nightly.records import read_run_dir
from amps_nightly.table import TableTest

CDT = dt.timezone(dt.timedelta(hours=-5))


def table(*specs):
    """specs: (name, family, extra result-file names)"""
    return [TableTest(name=n, family=f, selector="", line=i + 1, family_lines=1, result_files=set(x))
            for i, (n, f, x) in enumerate(specs)]


class EvaluateDiff(unittest.TestCase):
    def test_empty_file_passes(self):
        self.assertEqual(legacy.evaluate_diff(0, ""), ("PASS", "check", ""))

    def test_compile_started_but_not_done(self):
        status, stage, reason = legacy.evaluate_diff(30, "test_X_compile...\n")
        self.assertEqual((status, stage), ("FAILED", "compile"))
        self.assertIn("compilation failed", reason)

    def test_rundir_done_means_run_never_started(self):
        text = "test_X_compile...\ntest_X_compile... done\ntest_X_rundir...\ntest_X_rundir... done\n"
        status, stage, reason = legacy.evaluate_diff(len(text), text)
        self.assertEqual((status, stage), ("FAILED", "run"))
        self.assertIn("never executed", reason)

    def test_run_started_not_done(self):
        text = "test_X_rundir... done\ntest_X_run...\n"
        self.assertEqual(legacy.evaluate_diff(len(text), text)[:2], ("FAILED", "run"))

    def test_numerical_differences(self):
        text = ("DiffNum.pl -a=1e-30 -r=1e-12 run_test_CG/PT/plots/test_CG.dat output/test_CG/test_CG.ref_np8\n"
                "3n3\n< AMPS: 0 5.1e+02\n< AMPS: 1 5.2e+02\n")
        status, stage, reason = legacy.evaluate_diff(len(text), text)
        self.assertEqual((status, stage), ("FAILED", "check"))
        self.assertEqual(reason, "2 lines differ from reference test_CG.ref_np8")

    def test_missing_output(self):
        text = "ERROR in DiffNum.pl: run_test_A/PT/plots/test_A.dat does not exist\n"
        self.assertEqual(legacy.evaluate_diff(len(text), text),
                         ("FAILED", "check", "output not produced: test_A.dat"))

    def test_banner_printed_by_unfinished_run(self):
        text = ("test_E_rundir... done\ntest_E_run...\n========================================\n"
                "Running 3D forward tracking\n========================================\n")
        self.assertEqual(legacy.evaluate_diff(len(text), text),
                         ("FAILED", "run", "run did not finish (crash, abort or time limit)"))

    def test_truncated_listing_is_a_lower_bound(self):
        text = "DiffNum.pl a b/ref_np8\n" + "".join("< {}\n".format(i) for i in range(99))
        self.assertTrue(legacy.evaluate_diff(len(text), text)[2].startswith("at least 99 lines differ"))

    def test_whitespace_only(self):
        self.assertEqual(legacy.evaluate_diff(3, "\n \n")[0], "FAILED")


class ResultFiles(unittest.TestCase):
    def test_split_result_file(self):
        self.assertEqual(legacy.split_result_file("test_Saturn-MOP[5].diff"), ("Saturn-MOP", "5"))
        self.assertEqual(legacy.split_result_file("test_CG.diff"), ("CG", ""))

    def test_parse_res_with_heads(self):
        res = ("-rw-rw-r-- 1 u g 0 Oct 22 22:39 test_A.diff\n"
               "-rw-rw-r-- 1 u g 41 Oct 22 22:42 test_B[2].diff\n"
               "===============================================\n"
               "==> test_A.diff <==\n\n==> test_B[2].diff <==\ntest_B_compile...\n")
        files = legacy.parse_res(res)
        self.assertEqual(list(files), ["test_A.diff", "test_B[2].diff"])
        self.assertEqual(files["test_A.diff"][0], 0)
        self.assertEqual(files["test_B[2].diff"], (41, "test_B_compile..."))

    def test_parse_res_single_file_has_no_head_header(self):
        res = "-rw-r--r-- 1 u g 18 Oct 22 22:39 test_A.diff\n=====\ntest_A_compile...\n"
        self.assertEqual(legacy.parse_res(res)["test_A.diff"], (18, "test_A_compile..."))

    def test_collect_folds_extra_files_into_their_test(self):
        tests = table(("CG-dust-test", "comet-dust", ()), ("CG", "comet-dust", ()),
                      ("const-field-bc", "pic-field-solver", ("field-bc-constE_dirichlet",)),
                      ("Earth", "earth", ()))
        files = collections.OrderedDict([
            ("test_CG.diff", (0, "")),
            ("test_CG-dust-test.diff", (0, "")),
            ("test_CG-dust-test--DUST.diff", (90, "DiffNum.pl x y/ref_np8\n2n2\n< 1\n")),
            ("test_const-field-bc.diff", (0, "")),
            ("test_field-bc-constE_dirichlet.diff", (0, "")),
            ("test_Earth-renamed-test.diff", (20, "test_Earth-renamed-test_run...\n")),
            ("test_const-field-bc-test=0-wave.diff", (0, "")),
        ])
        got = {r["test"]: r for r in legacy.collect_results(files, tests)}
        self.assertEqual(got["CG"]["status"], "PASS")
        self.assertEqual(got["CG-dust-test"]["status"], "FAILED")
        self.assertIn("1 of 1 extra result files differ", got["CG-dust-test"]["reason"])
        # .diff contents decide the status but are not stored
        self.assertNotIn("excerpt", got["CG-dust-test"])
        self.assertNotIn("< 1", repr(got))
        self.assertEqual(got["const-field-bc"]["status"], "PASS")
        self.assertNotIn("field-bc-constE_dirichlet", got)
        self.assertNotIn("const-field-bc-test=0-wave", got)      # "-name=..." suffix is an extra file
        # a "-word" suffix is a different test, not an extra file of "Earth"
        self.assertEqual(got["Earth-renamed-test"]["family"], "unassigned")
        self.assertNotIn("Earth", got)


class Times(unittest.TestCase):
    def test_both_date_formats(self):
        want = dt.datetime(2026, 10, 22, 22, 30, 1, tzinfo=CDT)
        self.assertEqual(legacy.parse_legacy_date("Thu 22 Oct 2026 10:30:01 PM CDT"), want)
        self.assertEqual(legacy.parse_legacy_date("Thu Oct 22 22:30:01 CDT 2026"), want)
        self.assertIsNone(legacy.parse_legacy_date("yesterday"))

    def test_commit_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t",
                       GIT_COMMITTER_EMAIL="t@t")
            run = lambda *a, **kw: subprocess.run(["git", "-C", tmp] + list(a), check=True, stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE, env=dict(env, **kw), universal_newlines=True)
            run("init", "-q")
            run("symbolic-ref", "HEAD", "refs/heads/master")     # older git has no 'init -b'
            for day in ("2026-10-20T12:00:00-05:00", "2026-10-22T12:00:00-05:00"):
                run("commit", "-q", "--allow-empty", "-m", day, GIT_COMMITTER_DATE=day, GIT_AUTHOR_DATE=day)
            shas = run("log", "--format=%H").stdout.split()
            index = legacy.CommitIndex(Path(tmp))
            self.assertEqual(index.at(dt.datetime(2026, 10, 21, 22, 30, tzinfo=CDT)), shas[1])
            self.assertEqual(index.at(dt.datetime(2026, 10, 22, 22, 30, tzinfo=CDT)), shas[0])
            self.assertEqual(index.at(dt.datetime(2026, 10, 1, tzinfo=CDT)), "")


class ImportArchive(unittest.TestCase):
    def test_archive_import_skips_repeated_nights(self):
        res = ("-rw-rw-r-- 1 u g 0 Oct 21 22:39 test_A.diff\n-rw-rw-r-- 1 u g 18 Oct 21 22:40 test_B.diff\n"
               "=====\n==> test_A.diff <==\n\n==> test_B.diff <==\ntest_B_compile...\n")
        with tempfile.TemporaryDirectory() as tmp:
            archive, out = Path(tmp, "arch"), Path(tmp, "out")
            for day, checkout in (("2026/10/22", "Wed 21 Oct 2026 10:30:01 PM CDT"),
                                  ("2026/10/23", "Wed 21 Oct 2026 10:30:01 PM CDT"),   # stale copy
                                  ("2026/10/24", "Fri 23 Oct 2026 10:30:01 PM CDT")):
                folder = archive / day / "amps-gpu_gnu"
                folder.mkdir(parents=True)
                (folder / "test_amps.res").write_text(res)
                (folder / "test_amps.log").write_text("AMPS was checked out on {}\nmore\n".format(checkout))
            got = legacy.import_archive(archive, out, table(("A", "fam-a", ()), ("B", "fam-b", ())),
                                        {"amps-gpu_gnu": ("passion", "gnu")})
            self.assertEqual(len(got["written"]), 2)
            self.assertEqual(len(got["duplicate"]), 1)
            run, results = read_run_dir(out / got["written"][0])
            self.assertEqual((run["machine"], run["config"], run["column"]), ("passion", "gnu", "passion-gnu"))
            self.assertEqual(run["started"], "2026-10-22T03:30:01Z")
            self.assertEqual({r["test"]: r["status"] for r in results}, {"A": "PASS", "B": "FAILED"})


if __name__ == "__main__":
    unittest.main()
