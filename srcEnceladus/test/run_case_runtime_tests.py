#!/usr/bin/env python3
"""Actual C++ case/catalog resolver exercised with independent synthetic data."""
import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from case_fixtures import make_fixture
from source_case_fixtures import fixture_rates


class CaseRuntimeTests(unittest.TestCase):
    binary=None
    def invoke(self,deck,root,case,ok=True):
        p=subprocess.run([str(self.binary),"case_resolve",str(deck),str(root),case],text=True,capture_output=True)
        self.assertEqual(p.returncode==0,ok,p.stdout+p.stderr)
        return p.stdout+p.stderr

    def test_all_28_measurement_basis_cases_preserve_nonuniform_rate_totals(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,pub=make_fixture(Path(t))
            for folder in sorted((pub/"cases").iterdir()):
                with self.subTest(case=folder.name):
                    text=self.invoke(deck,root,folder.name)
                    rate=float(re.search(r"Q_s-1=([\deE+.-]+)",text).group(1))
                    self.assertTrue(math.isclose(rate,math.fsum(fixture_rates(folder.name.split("--")[1]).values()),rel_tol=1e-12))
            listing=subprocess.run([str(self.binary),"case_list",str(deck),str(root)],text=True,capture_output=True)
            self.assertEqual(listing.returncode,0,listing.stderr)
            self.assertEqual(listing.stdout.count("READY "),28)

    def test_cli_aliases_and_rejections(self):
        case="E3--tenishev2014_hybrid_8"
        for alias in ("-i","-input","--input","--input="):
            argv=[alias+"example.in"] if alias.endswith("=") else [alias,"example.in"]
            p=subprocess.run([str(self.binary),"case_cli",*argv,"--case="+case],capture_output=True)
            self.assertEqual(p.returncode,0,p.stderr)
        for argv in (("-case",),("-case","E3"),("-case","E3---basis"),("-case",case,"--case",case),
                     ("-i","a","-i","b"),("--list-cases","-case",case),("--fit-strategy","independent")):
            with self.subTest(argv=argv):
                p=subprocess.run([str(self.binary),"case_cli",*argv],capture_output=True)
                self.assertNotEqual(p.returncode,0)

    def test_missing_selector_unknown_case_and_missing_data_do_not_fallback(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,_=make_fixture(Path(t),all_cases=False)
            self.assertIn("selector",self.invoke(deck,root,"",False))
            self.assertIn("unknown case",self.invoke(deck,root,"E5--tenishev2014_hybrid_8",False))
            self.assertIn("model-data",self.invoke(deck,root/"missing","E3--tenishev2014_hybrid_8",False))

    def test_modified_source_is_not_accepted_under_old_hash(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,pub=make_fixture(Path(t),all_cases=False)
            p=next((pub/"cases").glob("*/sources.in"));p.write_text(p.read_text()+"! changed\n")
            self.assertIn("SHA-256 mismatch",self.invoke(deck,root,"E3--tenishev2014_hybrid_8",False))

    def test_symlink_escape_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,pub=make_fixture(Path(t),all_cases=False)
            p=next((pub/"cases").glob("*/sources.in"));outside=Path(t)/"outside.in"
            shutil.copy2(p,outside);p.unlink();p.symlink_to(outside)
            self.assertIn("escapes",self.invoke(deck,root,"E3--tenishev2014_hybrid_8",False))

    def test_identical_bundles_at_different_roots_have_same_fingerprint(self):
        with tempfile.TemporaryDirectory() as t:
            directory=Path(t);a=directory/"machine a";a.mkdir()
            deck,root,_,pub=make_fixture(a,all_cases=False)
            moved=directory/"machine b";moved.mkdir();shutil.copytree(root,moved/"data")
            shutil.copy2(deck,moved/"cases.in")
            x=self.invoke(deck,root,"E3--tenishev2014_hybrid_8")
            y=self.invoke(moved/"cases.in",moved/"data","E3--tenishev2014_hybrid_8")
            self.assertEqual(re.search(r"fingerprint=(\w+)",x).group(1),re.search(r"fingerprint=(\w+)",y).group(1))

    def test_shape_totals_required_and_absolute_totals_forbidden(self):
        for mode in ("shape","absolute"):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as t:
                deck,root,_,_=make_fixture(Path(t),all_cases=False,mode=mode)
                self.invoke(deck,root,"E3--tenishev2014_hybrid_8")
                s=deck.read_text()
                if mode=="shape":s=s.replace("total_source_rate.H2O_s-1=1e28\n","")
                else:s=s.replace("#end case","total_source_rate.H2O_s-1=1e28\n#end case")
                deck.write_text(s)
                self.invoke(deck,root,"E3--tenishev2014_hybrid_8",False)

    def test_case_overrides_domain_and_inactive_controls_are_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,_=make_fixture(Path(t),all_cases=False)
            original=deck.read_text()
            deck.write_text(original.replace("#end case","mesh.domain_max_m=2000000,1260500,1260500\n#end case"))
            self.assertIn("domain_max_m=2000000",self.invoke(deck,root,"E3--tenishev2014_hybrid_8"))
            for text in (original.replace("#end case","courant.CO2=0.5\n#end case"),
                         original.replace("#end case","fit_strategy=independent\n#end case"),
                         original.replace("-1260500,-1260500,-1260500","-100,-100,-100")):
                deck.write_text(text);self.invoke(deck,root,"E3--tenishev2014_hybrid_8",False)

    def test_original_raw_files_are_not_runtime_dependencies(self):
        with tempfile.TemporaryDirectory() as t:
            deck,root,_,publication=make_fixture(Path(t),all_cases=False)
            shutil.rmtree(Path(t)/"raw");shutil.rmtree(Path(t)/"observations")
            self.invoke(deck,root,"E3--tenishev2014_hybrid_8")
            # Explicit model-only/context requests must remain selected even
            # though 'none' is not a transported species. The catalog is the
            # authority; changing its metadata does not change asset hashes.
            catalog=publication/"enceladus_cases.in"
            text=catalog.read_text().replace("comparison_species=H2O","comparison_species=none")
            # A model-only request intentionally has no density/uncertainty
            # column mapping; keep the unchanged native validation contract.
            text=re.sub(r"^\s*(value_column|uncertainty_column)=.*\n","",text,flags=re.M)
            catalog.write_text(text)
            self.assertIn("observations=1",self.invoke(deck,root,"E3--tenishev2014_hybrid_8"))


def main():
    p=argparse.ArgumentParser();p.add_argument("--executable",type=Path,required=True)
    a=p.parse_args();CaseRuntimeTests.binary=a.executable.resolve()
    return 0 if unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CaseRuntimeTests)).wasSuccessful() else 1

if __name__=="__main__":raise SystemExit(main())
