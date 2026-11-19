#!/usr/bin/env python3
"""Infrastructure refusal/artifact checks; these are NOT lunar science tests.

Temporary predecessor JSONs exercise status protocol only. They are never
submitted as qualification evidence, input observations or a frozen reference.
The critical check is that even all-green fixture predecessors cannot turn a
missing production kinetic model into a formal PASS or an invented native run.
"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import integration_m3 as owner
import run_all_test as aggregate
import test_runner as global_owner


class FormalOwnerControls(unittest.TestCase):
    def test_literal_global_registration(self):
        self.assertEqual(aggregate.validate_registry(),[])
        self.assertEqual(global_owner.audit_architecture(),[])
        self.assertEqual({x[0] for x in aggregate.LINKED_M3_TESTS},{'I21','I22','I23','I24'})

    def test_missing_capabilities_stay_skipped_for_each_owner(self):
        for identity in ('I21','I22','I23','I24'):
            with self.subTest(identity=identity),tempfile.TemporaryDirectory() as raw:
                output=Path(raw)
                for number in range(13,21):
                    folder=output/f'I{number:02d}';folder.mkdir()
                    (folder/'result.json').write_text(json.dumps({'status':'PASS',
                        'scope':'synthetic orchestration fixture, NOT scientific evidence'}))
                value=owner.run(identity,owner.ROOT/'amps',
                    owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',output)
                self.assertEqual(value['status'],'SKIPPED')
                self.assertFalse(value['criteria_evaluated'])
                folder=output/identity
                self.assertFalse((folder/'receipt.txt').exists())
                self.assertFalse((folder/'comparison.png').exists())
                manifest=json.loads((folder/'run_manifest.json').read_text())
                self.assertEqual(manifest['execution'],'not-run')
                self.assertEqual(manifest['mpi_ranks'],8)
                self.assertTrue(manifest['input_snapshot_available'])
                normalized=aggregate.load_runner_result(identity,
                    owner.TEST/'integration-tests'/identity/'test.py',
                    folder/'result.json',folder/'stdout.log',77)
                self.assertEqual(normalized['status'],'SKIPPED')

    def test_missing_input_is_absent_not_fabricated(self):
        with tempfile.TemporaryDirectory() as raw:
            output=Path(raw)
            value=owner.run('I21',owner.ROOT/'amps',output/'missing.in',output)
            self.assertEqual(value['status'],'SKIPPED')
            self.assertFalse((output/'I21/application.in').exists())
            normalized=aggregate.load_runner_result('I21',owner.TEST/'integration-tests/I21/test.py',
                output/'I21/result.json',output/'I21/stdout.log',77)
            self.assertEqual(normalized['status'],'SKIPPED')

    def test_approved_protocol_is_recorded_without_stale_permission_blockers(self):
        with tempfile.TemporaryDirectory() as raw:
            output=Path(raw)
            value=owner.run('I22',owner.ROOT/'amps',
                owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',output)
            self.assertEqual(value['status'],'SKIPPED')
            self.assertEqual(value['approved_numerical_bias_limit'],.002)
            self.assertFalse(any(reason in value['blocking_reasons'] for reason in
                ('numerical error budget requires review before execution',
                 'full statistical family/stream/rare-branch/rate-estimator design requires review')))
            self.assertTrue((output/'I22/approved_qualification_addendum.json').is_file())
            self.assertFalse((output/'I22/receipt.txt').exists())

    def test_invalid_predecessor_is_error(self):
        with tempfile.TemporaryDirectory() as raw:
            output=Path(raw);(output/'I13').mkdir()
            (output/'I13/result.json').write_text('{"status":"invented-success"}')
            value=owner.run('I21',owner.ROOT/'amps',
                owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',output)
            self.assertEqual(value['status'],'ERROR')

    def test_capability_boolean_cannot_manufacture_execution(self):
        with tempfile.TemporaryDirectory() as raw:
            output=Path(raw);caps=output/'unsupported-activation.json'
            value=json.loads(owner.CAPABILITIES.read_text())
            value['capabilities']['I21']['missing_formal_capabilities']=[]
            caps.write_text(json.dumps(value))
            with patch.object(owner,'CAPABILITIES',caps):
                result=owner.run('I21',owner.ROOT/'amps',
                    owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',output)
            self.assertEqual(result['status'],'ERROR')
            self.assertFalse(result['criteria_evaluated'])


if __name__=='__main__':unittest.main()
