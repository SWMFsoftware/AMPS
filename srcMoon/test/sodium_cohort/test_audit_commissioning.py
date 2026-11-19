#!/usr/bin/env python3
"""Analytical helpers/negative provenance only; no fake native PASS receipts."""
import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_commissioning as a


class AnalyticalClock(unittest.TestCase):
    """Declared arithmetic fixtures, not physical rates or observations."""
    def fixture(self, when='1'):
        with localcontext() as context:
            context.prec = 80
            rate = -Decimal.from_float(a.numerical.source_uniform(731993, 'ANALYTICAL-NOT-NATIVE', 0)).ln() / Decimal(when)
        return rate, [{'original_parent_id':'0', 'status':'REACTED',
                       'relative_event_time_s':when, 'channel_id':'ANALYTICAL',
                       'initial_represented_number':'1'}]
    def evaluate(self, rate, rows, horizon='2'):
        return a.first_hit_comparisons(rows, rate, Decimal(horizon), 731993, 'ANALYTICAL-NOT-NATIVE')[0]
    def test_exact_analytical_exponential(self):
        rate, rows=self.fixture();r=self.evaluate(rate,rows)
        self.assertEqual(r['reference_first_event_s'],'1')
        self.assertEqual(r['normalized_timing_error'],'0')
        self.assertEqual(r['single_level_timing_status'],'PASS')
    def test_exact_existing_budget_boundary(self):
        rate,rows=self.fixture();rows[0]['relative_event_time_s']='1.004'
        self.assertEqual(self.evaluate(rate,rows)['single_level_timing_status'],'PASS')
    def test_no_threshold_relaxation(self):
        rate,rows=self.fixture();rows[0]['relative_event_time_s']='1.0040000000000000000001'
        self.assertEqual(self.evaluate(rate,rows)['single_level_timing_status'],'FAIL')
    def test_endpoint_first_hit_is_reaction(self):
        rate,rows=self.fixture('1');self.assertEqual(self.evaluate(rate,rows,'1')['expected_status'],'REACTED')
    def test_future_first_hit_is_survival(self):
        rate,rows=self.fixture('3');rows[0].update(status='SURVIVED',relative_event_time_s='')
        self.assertEqual(self.evaluate(rate,rows)['single_level_timing_status'],'PASS')
    def test_missing_expected_reaction_is_evaluated_fail(self):
        rate,rows=self.fixture();rows[0].update(status='SURVIVED',relative_event_time_s='')
        self.assertEqual(self.evaluate(rate,rows)['membership_status'],'FAIL')
    def test_zero_hazard_has_no_invented_lifetime(self):
        _,rows=self.fixture();rows[0].update(status='SURVIVED',relative_event_time_s='')
        result=self.evaluate(Decimal(0),rows)
        self.assertEqual(result['reference_first_event_s'],'')
        self.assertEqual(result['single_level_timing_status'],'PASS')
    def test_invalid_rate_and_horizon_refused(self):
        for rate,horizon in [(Decimal(-1),Decimal(1)),(Decimal('NaN'),Decimal(1)),(Decimal(1),Decimal(0))]:
            with self.subTest(values=(rate,horizon)),self.assertRaises(ValueError):
                a.first_hit_comparisons([],rate,horizon,731993,'ANALYTICAL')


class ParentRecords(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='commissioning-parser-only-');self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)/'parents.csv'
        self.rows=[dict(zip(a.q.PARENT_COLUMNS,('0','SURVIVED','','','1')))]
        self.receipt={'original_parents':1,'physical_horizon_s':'1','initial_represented_number':'1'}
    def write(self):
        with self.path.open('w',newline='')as stream:
            w=csv.DictWriter(stream,fieldnames=a.q.PARENT_COLUMNS);w.writeheader();w.writerows(self.rows)
    def read(self):self.write();return a.read_parents(self.path,self.receipt,'ANALYTICAL')
    def test_small_survivor_cohort_needs_no_statistical_plan(self):
        self.assertEqual(len(self.read()),1)
    def test_complete_range_required(self):
        self.rows[0]['original_parent_id']='1'
        with self.assertRaises(ValueError):self.read()
    def test_duplicate_ids_refused(self):
        self.receipt['original_parents']=2;self.rows*=2
        with self.assertRaises(ValueError):self.read()
    def test_censoring_refused(self):
        self.rows[0]['status']='ESCAPED'
        with self.assertRaises(ValueError):self.read()
    def test_invented_survivor_event_refused(self):
        self.rows[0]['relative_event_time_s']='0.5'
        with self.assertRaises(ValueError):self.read()
    def test_wrong_source_channel_refused(self):
        self.rows[0].update(status='REACTED',relative_event_time_s='.5',channel_id='UNKNOWN')
        with self.assertRaises(ValueError):self.read()


class NegativeBinding(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='commissioning-negative-only-');self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.binding=self.root/'binding.json'
    def run_audit(self):return a.run(self.binding,self.root/'audit')
    def test_missing_binding_is_skip(self):
        r=self.run_audit();self.assertEqual(r['status'],'SKIPPED')
        self.assertEqual(r['refinement_status'],'NOT_EVALUATED')
        self.assertEqual(r['statistical_qualification_status'],'NOT_EVALUATED')
    def test_malformed_binding_is_error(self):
        self.binding.write_text('{}');self.assertEqual(self.run_audit()['status'],'ERROR')
    def test_duplicate_binding_key_is_error(self):
        self.binding.write_text('{"schema":"x","schema":"y"}')
        self.assertEqual(self.run_audit()['status'],'ERROR')
    def test_corrupt_profile_is_error_before_any_owner(self):
        profile=self.root/'profile.json';profile.write_text('ADMINISTRATIVE NOT A PROFILE')
        budget=self.root/'budget.json';budget.write_text('ADMINISTRATIVE NOT AN APPROVAL')
        self.binding.write_text(json.dumps({'schema':a.BINDING_SCHEMA,
            'profile':{'path':str(profile),'sha256':'0'*64},
            'budget':{'path':str(budget),'sha256':a.ledger.digest(budget)},'owners':[]}))
        self.assertEqual(self.run_audit()['status'],'ERROR')
    def test_no_actual_owner_is_skip_not_fake_pass(self):
        profile=self.root/'profile.json';profile.write_text('ADMINISTRATIVE')
        budget=self.root/'budget.json';budget.write_text('ADMINISTRATIVE')
        self.binding.write_text(json.dumps({'schema':a.BINDING_SCHEMA,
            'profile':{'path':str(profile),'sha256':a.ledger.digest(profile)},
            'budget':{'path':str(budget),'sha256':a.ledger.digest(budget)},'owners':[]}))
        self.assertEqual(self.run_audit()['status'],'SKIPPED')
    def test_output_nonoverwrite(self):
        self.run_audit();p=self.root/'audit/result.json';before=p.read_bytes()
        with self.assertRaises(FileExistsError):self.run_audit()
        self.assertEqual(p.read_bytes(),before)


if __name__=='__main__':unittest.main(verbosity=2)
