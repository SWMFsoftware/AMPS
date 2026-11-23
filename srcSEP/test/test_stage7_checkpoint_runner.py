#!/usr/bin/env python3
"""Isolated evidence-corruption controls; never native physics qualification."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import run_stage7_checkpoint_tests as gate


class Stage7EvidenceControls(unittest.TestCase):
    def plan(self,p,job):
        # A complete synthetic report catalog is required even in reader
        # controls. Missing receipts for unrelated roles remain explicit ERROR;
        # only the chosen case's defect/status is asserted by these unit tests.
        job.update(ranks=2)
        key=lambda j:(j['id'],j['role'],j['threads'],j['ranks'])
        others=[dict(id=i,role=r,threads=t,ranks=n,cwd=str(p/f'missing-{i}-{r}'),argv=['unit-only'])
                for i,r,t,n in sorted(gate.expected_job_keys(True)) if (i,r,t,n)!=key(job)]
        return dict(hashes={},jobs=[job]+others,fixtures_only=True,acceptance={})
    def ledger(self,**changes):
        # Exact manufactured open-system arithmetic: 4 particles, 10 J birth,
        # +2 J prescribed work, 8 J retained and 4 J returned event energy.
        v=dict(epoch_s=12300,born=4,retained=3,returned=1,domain_exit=0,geometry_retired=0,
            born_kinetic_j=10,retained_birth_j=6,returned_birth_j=4,domain_birth_j=0,retired_birth_j=0,
            retained_kinetic_j=8,returned_event_j=4,domain_event_j=0,retired_event_j=0,
            background_work_j=2,wave_work_j=0)
        v.update(changes)
        return 'REDUCED_SOURCE_LEDGER '+' '.join(f'{k}={x}' for k,x in v.items())+'\n'

    def test_wrong_signed_work_birth_or_physical_endpoint_fails(self):
        self.assertEqual(gate.audit_budget(self.ledger(),[12300]),0)
        for changes in [dict(background_work_j=-2),dict(retained_birth_j=5),dict(epoch_s=12301),dict(retained=2)]:
            with self.assertRaises(gate.GateFailure):gate.audit_budget(self.ledger(**changes),[12300])

    def test_empty_short_driver_cannot_qualify(self):
        v={name:0 for name in ['born','retained','returned','born_kinetic_j','retained_birth_j','returned_birth_j',
                              'retained_kinetic_j','returned_event_j','background_work_j']}
        with self.assertRaises(gate.GateFailure):gate.audit_budget(self.ledger(**v),[12300])

    def test_missing_receipts_and_contradictory_reports_are_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);job=dict(id='CP05C01',role='routine',threads=0,cwd=str(p),argv=['unit-only'])
            plan=self.plan(p,job)
            with contextlib.redirect_stdout(io.StringIO()):gate.collect(plan,p)
            self.assertEqual(json.loads((p/'summary.json').read_text())['results'][0]['status'],'ERROR')
            (p/'process-result.json').write_text(json.dumps(dict(argv=job['argv'],exit_code=0)))
            (p/'results.json').write_text(json.dumps(dict(results=[dict(id=job['id'],status='FAIL')])))
            with contextlib.redirect_stdout(io.StringIO()):gate.collect(plan,p)
            self.assertEqual(json.loads((p/'summary.json').read_text())['results'][0]['status'],'ERROR')

    def test_missing_artifact_is_error_and_wrong_bytes_are_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);job=dict(id='CP05F01',role='resume',threads=0,cwd=str(p),argv=['unit-only'],prefix=str(p/'physics'))
            plan=self.plan(p,job)
            (p/'process-result.json').write_text(json.dumps(dict(argv=job['argv'],exit_code=0)))
            (p/'results.json').write_text(json.dumps(dict(results=[dict(id=job['id'],status='PASS')])))
            with contextlib.redirect_stdout(io.StringIO()):gate.collect(plan,p)
            self.assertEqual(json.loads((p/'summary.json').read_text())['results'][0]['status'],'ERROR')
            (p/'physics.expected').write_bytes(b'unchanged reference process')
            (p/'physics.resumed').write_bytes(b'wrong resumed state')
            with contextlib.redirect_stdout(io.StringIO()):gate.collect(plan,p)
            self.assertEqual(json.loads((p/'summary.json').read_text())['results'][0]['status'],'FAIL')

    def test_source_identity_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);source=p/'physics.cpp';source.write_text('original')
            plan=dict(hashes={str(source):gate.digest(source)},jobs=[],acceptance={})
            source.write_text('changed')
            with self.assertRaisesRegex(ValueError,'identity changed'):gate.collect(plan,p)

    def test_prepared_only_or_missing_roles_cannot_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with self.assertRaisesRegex(ValueError,'catalog'):
                gate.collect(dict(hashes={},jobs=[],fixtures_only=False,acceptance={}),p)


if __name__=='__main__':unittest.main()
