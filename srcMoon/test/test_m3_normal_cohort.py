#!/usr/bin/env python3
"""Administrative native-receipt controls, not observations or PIC evidence."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from run_m3_normal_cohort import _checks,python_uniform,run_startup

class ReceiptControls(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)
    def write(self,rows,summary='M3 fixture checks=2 status=PASS'):
        with (self.path/'semantic.checks.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['check','status']);writer.writerows(rows)
        (self.path/'native.log').write_text(summary+'\n')
    def test_complete_exact_named_contract_is_evaluated(self):
        self.write([('a','PASS'),('b','PASS')]);self.assertEqual(len(_checks(self.path,{'a','b'},r'M3 fixture checks=(\d+) status=(PASS|FAIL|ERROR)',{'returncode':0})),2)
    def test_empty_or_duplicate_rows_cannot_make_a_green_receipt(self):
        for rows in ([],[('a','PASS'),('a','PASS')]):
            self.write(rows)
            with self.assertRaises(ValueError):_checks(self.path,{'a','b'},r'M3 fixture checks=(\d+) status=(PASS|FAIL|ERROR)',{'returncode':0})
    def test_signal_and_arbitrary_exit_never_qualify_even_with_matching_text(self):
        self.write([('a','FAIL'),('b','PASS')],'M3 fixture checks=2 status=FAIL')
        for code in (-11,139,3,77):
            with self.assertRaises(ValueError):_checks(self.path,{'a','b'},r'M3 fixture checks=(\d+) status=(PASS|FAIL|ERROR)',{'returncode':code})
    def test_evaluated_failure_and_bad_exit_are_distinct(self):
        self.write([('a','FAIL'),('b','PASS')],'M3 fixture checks=2 status=FAIL')
        self.assertEqual(_checks(self.path,{'a','b'},r'M3 fixture checks=(\d+) status=(PASS|FAIL|ERROR)',{'returncode':1})[0]['status'],'FAIL')
        with self.assertRaises(ValueError):_checks(self.path,{'a','b'},r'M3 fixture checks=(\d+) status=(PASS|FAIL|ERROR)',{'returncode':0})
    def test_missing_normal_build_is_retained_as_skipped(self):
        result=run_startup(self.path/'missing',30,None)
        self.assertEqual(result['status'],'SKIPPED');self.assertEqual(json.loads((self.path/'missing/result.json').read_text())['status'],'SKIPPED')
    def test_malformed_normal_proof_is_retained_as_error(self):
        result=run_startup(self.path/'bad',30,[])
        self.assertEqual(result['status'],'ERROR');self.assertTrue((self.path/'bad/infrastructure_exception.log').is_file())
    def test_independent_uniform_is_open_and_full_child_key_distinct(self):
        a,b=python_uniform(73,0),python_uniform(73,1)
        self.assertTrue(0<a<1 and 0<b<1);self.assertNotEqual(a,b)

if __name__=='__main__':unittest.main()
