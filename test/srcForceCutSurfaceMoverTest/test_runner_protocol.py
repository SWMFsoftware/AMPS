#!/usr/bin/env python3
"""Reporting-only fixtures: fabricated rows are NEVER physical/observational evidence."""
import tempfile,unittest
from pathlib import Path
from run import score_native_report

class NativeProtocol(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='force-tracer-PROTOCOL-ONLY-');self.addCleanup(self.temp.cleanup)
        self.report=Path(self.temp.name)/'PROTOCOL_ONLY.report'
    def fixture(self,count=124,failures=0,mode='fixed'):
        rows=[('FAIL' if i<failures else 'PASS')+' order='+str(1+i%2)+' PROTOCOL ONLY row '+str(i)+' value=0 expected=0' for i in range(count)]
        rows.append(f'SUMMARY checks={count} passed={count-failures} failed={failures}')
        self.report.write_text('\n'.join(rows)+'\n')
    def score(self,code=0,stdout='',mode='fixed'):
        return score_native_report(self.report,stdout,code,mode)[0]
    def test_complete_fixed_protocol(self):self.fixture();self.assertEqual(self.score(),'PASS')
    def test_complete_adaptive_protocol(self):self.fixture(134);self.assertEqual(self.score(mode='adaptive'),'PASS')
    def test_real_failed_criterion_preserves_fail(self):self.fixture(failures=1);self.assertEqual(self.score(1),'FAIL')
    def test_runtime_abort_is_error_even_if_exit_one(self):self.fixture(failures=1);self.assertEqual(self.score(1,'AMPS: exit: deliberate protocol fixture'),'ERROR')
    def test_missing_report_is_error(self):self.assertEqual(self.score(),'ERROR')
    def test_truncated_report_is_error(self):self.fixture();self.report.write_text('\n'.join(self.report.read_text().splitlines()[:-1]));self.assertEqual(self.score(),'ERROR')
    def test_removed_checks_are_error(self):self.fixture(102);self.assertEqual(self.score(),'ERROR')
    def test_contradictory_exit_is_error(self):self.fixture();self.assertEqual(self.score(1),'ERROR')
    def test_corrupted_nan_metric_is_error(self):self.fixture();self.report.write_text(self.report.read_text().replace('value=0','value=nan',1));self.assertEqual(self.score(),'ERROR')
    def test_count_contradiction_is_error(self):self.fixture();self.report.write_text(self.report.read_text().replace('failed=0','failed=1'));self.assertEqual(self.score(),'ERROR')
if __name__=='__main__':unittest.main()
