#!/usr/bin/env python3
"""Infrastructure refusal/artifact checks; these are NOT lunar science tests.

Temporary predecessor JSONs exercise status protocol only. They are never
submitted as qualification evidence, input observations or a frozen reference.
The critical check is that even all-green fixture predecessors cannot turn a
missing production kinetic model into a formal PASS or an invented native run.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import integration_m3 as owner
import run_all_test as aggregate
import test_runner as global_owner


class FormalOwnerControls(unittest.TestCase):
    def test_i24_expired_checkout_and_changed_evidence_are_distinct(self):
        # Exercise the build-receipt admission code only. These administrative
        # bytes and patched production identity cannot launch a physical case.
        with patch.object(sys,'path',[str(owner.TEST/'integration-tests/I24'),*sys.path]):
            import run_pic_ground_capture as pic
        with tempfile.TemporaryDirectory() as work:
            root=Path(work)
            def record(path,content):
                path.parent.mkdir(exist_ok=True,parents=True)
                path.write_text(content)
                return dict(path=str(path),sha256=owner.sha(path))
            library=record(root/'scratch/build/libAMPS.a','fixture-library')
            saved=record(root/'generated/pic.h','fixture-ABI')
            generated=record(root/'scratch/build/pic.h','fixture-ABI')
            path=root/'normal_build_evidence.json'
            recipe=record(root/'recipe.json','administrative fixture')
            recipe.update(packet_layout='conservative-products',reference_problem='homogeneous-zero-acceleration',
                          mapping=dict(record(root/'mapping.json','fixture'),states=pic.MAPPING),
                          profile=record(root/'profile.json','fixture'))
            payload=dict(schema='srcMoon-actual-normal-restore-producer-v1',status='PASS',
                source_qualification_recipe=recipe, production_source_hashes={},
                build_log=record(root/'build.log','EXIT: 0\nEXIT: 0\n'),
                executable=record(root/'amps','fixture'),library=library,
                retained_generated_snapshots={saved['path']:saved['sha256']})
            path.write_text(json.dumps(payload))
            with patch.object(pic,'production_source_identity',return_value={}), \
                    patch.object(pic,'native_sources',side_effect=lambda value:value):
                self.assertEqual(pic.checked_build(path),root/'scratch/build')
                Path(generated['path']).write_text('altered header')
                with self.assertRaisesRegex(ValueError,'changed configured build member'):
                    pic.checked_build(path)
                Path(library['path']).unlink()
                with self.assertRaisesRegex(ValueError,'changed configured build member'):
                    pic.checked_build(path)
                Path(generated['path']).write_text('fixture-ABI')
                with self.assertRaisesRegex(pic.MissingPrerequisite,'disposable configured build member'):
                    pic.checked_build(path)
                Path(saved['path']).write_text('altered retained snapshot')
                with self.assertRaisesRegex(ValueError,'changed retained generated ABI'):
                    pic.checked_build(path)

    def test_i24_pic_registration_requires_explicit_sodium_lock(self):
        # Validate the actual aggregate CLI before it configures or runs a case.
        command = [sys.executable, '-B', str(owner.TEST/'run_all_test.py'), '--skip-build',
                   '--only', 'I24', '--i24-pic-producer', 'absent-producer.json']
        attempt = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(attempt.returncode, 2)
        self.assertIn('--i24-pic-producer requires --i24-sodium-lock', attempt.stderr)

    def test_i24_report_keeps_stronger_child_status(self):
        # These temporary fixtures exercise bookkeeping only. No MPI is
        # launched and no fixture survives this test as scientific evidence.
        for statuses, expected in ((('PASS','PASS'),'SKIPPED'),(('PASS','FAIL'),'FAIL'),
                                   (('FAIL','ERROR'),'ERROR'),(('ERROR','FAIL'),'ERROR')):
            with self.subTest(statuses=statuses), tempfile.TemporaryDirectory() as work:
                root = Path(work)
                iterator = iter(statuses)
                def issue_control(command, cwd, log, *, output):
                    child = Path(command[command.index('--output-dir')+1])
                    child.mkdir()
                    status = next(iterator)
                    (child/'result.json').write_text(json.dumps({'status':status,
                        'reason':'temporary orchestration fixture'}))
                    index = {}
                    if status == 'FAIL':
                        if child.name == 'sodium-ground-capture':
                            (child/'native_receipt.json').write_text('{"fixture":true}')
                        else:
                            csv = child/'fixture-events-rank-0.csv'
                            csv.write_text('temporary orchestration fixture')
                            index = {'result_files':{csv.name:owner.sha(csv)}}
                    (child/'run_manifest.json').write_text(json.dumps(index))
                    return {'reason':None,'returncode':{'PASS':0,'FAIL':1,'ERROR':2}[status]}
                with patch.object(owner.CampaignBudget,'run',side_effect=issue_control):
                    result = owner.run('I24',owner.ROOT/'amps',owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',
                                       root,i24_lock=root/'fixture-lock',i24_pic_producer=root/'fixture-build')
                self.assertEqual(result['formal_status'],'SKIPPED / NOT QUALIFIED')
                self.assertEqual(result['status'],expected)
                self.assertEqual(tuple(row['status'] for row in result['supporting_checks']),statuses)
                # Run the final summary's real schema/count/severity validator.
                aggregate.validate_supporting_results(result,result['status'])
                self.assertEqual((root/'I24/receipt.txt').is_file(),'FAIL' in statuses)

    def test_i24_aggregate_retains_both_independent_children(self):
        # Exercise both real child CLIs via the public aggregate. One broken
        # input at a time tests that either child can return ERROR without
        # losing the other child's SKIP and its separately hashed evidence.
        for damaged in ('none','lock','producer'):
            with self.subTest(damaged=damaged), tempfile.TemporaryDirectory() as work:
                root = Path(work)
                lock, producer = root/'database-lock.json', root/'producer.json'
                if damaged == 'lock':
                    lock.write_text('{}')
                elif damaged == 'producer':
                    producer.write_text('{}')
                command = [sys.executable, '-B', str(owner.TEST/'run_all_test.py'), '--only','I24',
                           '--skip-build','--i24-sodium-lock',str(lock),'--i24-pic-producer',str(producer),
                           '--i24-pic-refinement',
                           '--test-runtime-seconds','60','--output-dir',str(root/'out')]
                attempt = subprocess.run(command, capture_output=True, text=True,timeout=90)
                folder = root/'out/I24'
                result = json.loads((folder/'result.json').read_text())
                manifest = json.loads((folder/'run_manifest.json').read_text())
                expected = ['ERROR' if damaged=='lock' else 'SKIPPED',
                            'ERROR' if damaged=='producer' else 'SKIPPED']
                self.assertEqual([r['status'] for r in result['supporting_checks']], expected)
                self.assertEqual(result['status'],'SKIPPED' if damaged=='none' else 'ERROR')
                self.assertEqual(result['formal_status'],'SKIPPED / NOT QUALIFIED')
                self.assertEqual(attempt.returncode,0 if damaged=='none' else 2,attempt.stdout+attempt.stderr)
                for key, status in zip(('i24_sodium_provider','i24_pic_commissioning'),expected):
                    entry = manifest[key]
                    self.assertEqual(json.loads(Path(entry['result']).read_text())['status'],status)
                    self.assertEqual(entry['result_sha256'],owner.sha(entry['result']))
                    self.assertEqual(entry['manifest_sha256'],owner.sha(entry['manifest']))
                    self.assertTrue(Path(entry['log']).is_file())
                    self.assertEqual('--refinement' in entry['supervision']['command'],
                                     key=='i24_pic_commissioning')
                # Neither invalid configuration evaluated a physical criterion.
                self.assertFalse((folder/'receipt.txt').exists())

    def test_i24_absent_and_corrupt_requested_database_have_different_statuses(self):
        # No science data are synthesized. An empty file is deliberately an
        # invalid release lock, while the absent file exercises a dependency
        # SKIP on the same public owner and aggregator.
        for exists,expected,code in ((False,'SKIPPED',77),(True,'ERROR',2)):
            with self.subTest(existing_corrupt_lock=exists),tempfile.TemporaryDirectory() as work:
                folder=Path(work)
                lock=folder/'database-lock.json'
                if exists:lock.touch()
                result=owner.run('I24',owner.ROOT/'amps',owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',
                    folder,i24_lock=lock)
                self.assertEqual(result['supporting_checks'][0]['status'],expected)
                self.assertEqual(result['status'],expected)
                self.assertFalse((folder/'I24/receipt.txt').exists())
                aggregate_result=aggregate.load_runner_result('I24',owner.TEST/'integration-tests/I24/test.py',
                    folder/'I24/result.json',folder/'I24/stdout.log',code)
                self.assertEqual(aggregate_result['status'],expected)

    def test_i23_pic_registration_requires_explicit_helium_profile(self):
        # Validate the actual aggregate CLI before it configures or runs a case.
        command = [sys.executable, '-B', str(owner.TEST/'run_all_test.py'), '--skip-build',
                   '--only', 'I23', '--i23-pic-producer', 'absent-producer.json']
        attempt = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(attempt.returncode, 2)
        self.assertIn('--i23-pic-producer requires --i23-helium-profile', attempt.stderr)
        command = [sys.executable, '-B', str(owner.TEST/'run_all_test.py'), '--skip-build',
                   '--only', 'I24', '--i23-helium-profile', 'absent-profile.json']
        attempt = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(attempt.returncode, 2)
        self.assertIn('--i23-helium-profile requires I23 in regression mode', attempt.stderr)

    def test_i23_report_keeps_stronger_child_status(self):
        # Bookkeeping fixtures only: no MPI, no physical evidence survives.
        for statuses, expected in ((('PASS','PASS'),'SKIPPED'),(('PASS','FAIL'),'FAIL'),
                                   (('FAIL','ERROR'),'ERROR'),(('SKIPPED','PASS'),'SKIPPED')):
            with self.subTest(statuses=statuses), tempfile.TemporaryDirectory() as work:
                root = Path(work)
                iterator = iter(statuses)
                def issue_control(command, cwd, log, *, output):
                    child = Path(command[command.index('--output-dir')+1])
                    child.mkdir()
                    status = next(iterator)
                    (child/'result.json').write_text(json.dumps({'status':status,
                        'reason':'temporary orchestration fixture'}))
                    index = {}
                    if status == 'FAIL':
                        if child.name == 'helium-electron-bath':
                            (child/'native_receipt.json').write_text('{"fixture":true}')
                        else:
                            csv = child/'fixture-events-rank-0.csv'
                            csv.write_text('temporary orchestration fixture')
                            index = {'result_files':{csv.name:owner.sha(csv)}}
                    (child/'run_manifest.json').write_text(json.dumps(index))
                    return {'reason':None,'returncode':{'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]}
                with patch.object(owner.CampaignBudget,'run',side_effect=issue_control):
                    result = owner.run('I23',owner.ROOT/'amps',owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',
                                       root,i23_profile=root/'fixture-profile',i23_pic_producer=root/'fixture-build')
                self.assertEqual(result['formal_status'],'SKIPPED / NOT QUALIFIED')
                self.assertEqual(result['status'],expected)
                self.assertEqual(tuple(row['status'] for row in result['supporting_checks']),statuses)
                aggregate.validate_supporting_results(result,result['status'])
                self.assertEqual((root/'I23/receipt.txt').is_file(),'FAIL' in statuses)

    def test_i23_absent_and_corrupt_profile_have_different_statuses(self):
        # The real child CLI authenticates the profile before compiling. An
        # absent profile is a dependency SKIP; wrong bytes are an ERROR.
        for exists,expected,code in ((False,'SKIPPED',77),(True,'ERROR',2)):
            with self.subTest(existing_corrupt_profile=exists),tempfile.TemporaryDirectory() as work:
                folder=Path(work)
                profile=folder/'profile.json'
                if exists:profile.write_text('{}')
                result=owner.run('I23',owner.ROOT/'amps',owner.ROOT/'srcMoon/examples/m3_reaction_rates.in',
                    folder,i23_profile=profile,wall_seconds=120)
                self.assertEqual(result['supporting_checks'][0]['status'],expected)
                self.assertEqual(result['status'],expected)
                self.assertEqual(result['formal_status'],'SKIPPED / NOT QUALIFIED')
                self.assertFalse((folder/'I23/receipt.txt').exists())
                aggregate_result=aggregate.load_runner_result('I23',owner.TEST/'integration-tests/I23/test.py',
                    folder/'I23/result.json',folder/'I23/stdout.log',code)
                self.assertEqual(aggregate_result['status'],expected)

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
