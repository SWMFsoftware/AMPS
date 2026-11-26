#!/usr/bin/env python3
"""Administrative family controls and independent endpoint comparison tests."""
from pathlib import Path
import sys
import tempfile
import unittest

import m3_qualification as q
import run_m3_campaign as campaign
from campaign_budget import CampaignBudget

sys.path.insert(0, str(Path(__file__).resolve().parent/'sodium_cohort'))
from audit import compare_parent_layouts


class FamilyControls(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.plan = self.root/'campaign.json'
        self.output = self.root/'output'

    def family(self, *cases):
        campaign.write_json(self.plan, {'schema': campaign.SCHEMA,
            'approved_addendum_sha256': q.APPROVED_SHA256, 'cases': list(cases)})
        return q.sha(self.plan)

    def input_path(self, name, text):
        path = self.root/name
        path.write_text(text)
        return path

    def test_frozen_order_and_relative_input_resolution(self):
        a, b = [self.input_path(name, name) for name in ('one.json', 'two.json')]
        digest = self.family(
            {'name': 'first', 'kind': 'i21-cohort', 'inputs': {'plan': {'path': a.name, 'sha256': q.sha(a)}}},
            {'name': 'second', 'kind': 'i21-cohort', 'inputs': {'plan': {'path': b.name, 'sha256': q.sha(b)}}},
        )
        cases = campaign.read_family(self.plan, digest)
        self.assertEqual([case['name'] for case in cases], ['first','second'])
        self.assertEqual(cases[1]['inputs']['plan']['path'], str(b))

    def test_i23_electron_bath_cases_bind_profile_and_producer(self):
        # Administrative identities only: the commands are formed, not run.
        profile = self.input_path('profile.json', 'administrative profile identity')
        producer = self.input_path('producer.json', 'administrative producer identity')
        digest = self.family(
            {'name': 'he-ensemble', 'kind': 'i23-helium-ensemble',
             'inputs': {'profile': campaign.descriptor(profile)}},
            {'name': 'he-pic', 'kind': 'i23-helium-pic',
             'inputs': {'profile': campaign.descriptor(profile), 'producer': campaign.descriptor(producer)}},
        )
        ensemble, pic = campaign.read_family(self.plan, digest)
        self.assertEqual((ensemble['test_id'], pic['test_id']), ('I23', 'I23'))
        command = campaign.case_command(ensemble, self.root/'a', self.root, 60)
        self.assertEqual(command[0], str(campaign.RUNNERS['i23-helium-ensemble']))
        self.assertEqual(command[command.index('--profile')+1], str(profile))
        command = campaign.case_command(pic, self.root/'b', self.root, 60)
        self.assertIn('--refinement', command)
        self.assertEqual(command[command.index('--producer-manifest')+1], str(producer))
        # A PIC case without its producer is refused before any launch.
        digest = self.family({'name': 'he-pic', 'kind': 'i23-helium-pic',
                              'inputs': {'profile': campaign.descriptor(profile)}})
        with self.assertRaisesRegex(ValueError, 'exactly these pinned input roles'):
            campaign.read_family(self.plan, digest)

    def test_i24_total_capture_cases_bind_own_case_and_producer(self):
        # Administrative identities only: the commands are formed, not run.
        case = campaign.descriptor(campaign.TOTAL_CAPTURE_CASE)
        producer = self.input_path('producer.json', 'administrative producer identity')
        digest = self.family(
            {'name': 'na-total-ensemble', 'kind': 'i24-total-capture-ensemble', 'inputs': {'case': case}},
            {'name': 'na-total-pic', 'kind': 'i24-total-capture-pic',
             'inputs': {'case': case, 'producer': campaign.descriptor(producer)}},
        )
        ensemble, pic = campaign.read_family(self.plan, digest)
        self.assertEqual((ensemble['test_id'], pic['test_id']), ('I24', 'I24'))
        command = campaign.case_command(ensemble, self.root/'a', self.root, 60)
        self.assertEqual(command[0], str(campaign.RUNNERS['i24-total-capture-ensemble']))
        self.assertNotIn('--lock', command)
        command = campaign.case_command(pic, self.root/'b', self.root, 60)
        self.assertIn('--refinement', command)
        self.assertEqual(command[command.index('--producer-manifest')+1], str(producer))
        # A different case file cannot be substituted for the runner's own.
        foreign = self.input_path('case.json', 'administrative foreign case identity')
        digest = self.family({'name': 'na-total-ensemble', 'kind': 'i24-total-capture-ensemble',
                              'inputs': {'case': campaign.descriptor(foreign)}})
        other, = campaign.read_family(self.plan, digest)
        with self.assertRaisesRegex(ValueError, 'its own pinned case file'):
            campaign.case_command(other, self.root/'c', self.root, 60)

    def test_duplicate_physical_matrix_rejected(self):
        binding = {'plan': campaign.descriptor(self.input_path('cohort.json','administrative identity'))}
        digest = self.family(*({'name': name, 'kind': 'i21-cohort', 'inputs': binding}
                               for name in ('first','second')))
        with self.assertRaisesRegex(ValueError, 'duplicate physical matrix'):
            campaign.read_family(self.plan, digest)

    def test_bad_outer_pin_records_error_before_launch(self):
        self.family()
        result = campaign.run(self.plan, '0'*64, self.output)
        self.assertEqual(result['status'], 'ERROR')
        self.assertEqual(result['formal'], [])
        self.assertIn('changed', result['reason'])

    def test_missing_preflight_input_is_skipped_without_physical_execution(self):
        source = self.root/'absent.json'
        digest = self.family({'name': 'dark', 'kind': 'i21-cohort', 'inputs':
                             {'plan': {'path': str(source), 'sha256': '0'*64}}})
        result = campaign.run(self.plan, digest, self.output)
        self.assertEqual(result['status'], 'SKIPPED')
        self.assertEqual(result['formal_counts']['SKIPPED'], 4)
        self.assertEqual(result['independent_check_counts']['PASS'], 0)
        self.assertFalse((self.output/'cases').exists())
        self.assertIn('unavailable', Path(result['cases'][0]['log']).read_text())
        self.assertEqual(q.strict_object(self.output/'progress.json'), result)

    def test_changed_input_error_survives_missing_neighbour(self):
        # One altered pin must be reported even if its case cannot otherwise
        # run because a different dependency is unavailable.
        wrong = self.input_path('inventory.json', '{}')
        digest = self.family({'name': 'capture', 'kind': 'i24-sodium-ensemble', 'inputs':
                             {'lock': {'path': 'missing', 'sha256': '0'*64},
                              'mass_inventory': {'path': str(wrong), 'sha256': '0'*64}}})
        result = campaign.run(self.plan, digest, self.output)
        self.assertEqual(result['status'], 'ERROR')
        self.assertIn('changed', result['cases'][0]['reason'])
        self.assertFalse((self.output/'cases').exists())

    def test_result_exit_disagreement_is_error(self):
        # Administrative subprocess fixture for the supervisor, with no
        # executable model, particle files or independent science artifacts.
        child = self.root/'child'
        child.mkdir()
        campaign.write_json(child/'result.json', {'status':'PASS'})
        campaign.write_json(child/'run_manifest.json', {})
        runner = self.input_path('exit.py', 'raise SystemExit(1)\n')
        with self.assertRaisesRegex(ValueError, 'result and exit disagree'):
            campaign.supervise({}, [str(runner)], child, self.root/'process.log',
                               self.root/'resources', CampaignBudget(8,1), self.root)

    def test_malformed_independent_child_receipts_are_rejected(self):
        for result in ({'checks': []}, {'checks': {'first': {'budget': 'PASS'}}},
                       {'checks': {'first': [True]}}):
            with self.assertRaisesRegex(ValueError, 'invalid PIC layout or check status'):
                campaign.checks_from_case('i24-sodium-pic', self.root, result)
        child = self.root/'independent-audit'
        child.mkdir()
        campaign.write_json(child/'result.json', {'status': 'FAIL', 'checks': []})
        with self.assertRaisesRegex(ValueError, 'owner and independent audit disagree'):
            campaign.checks_from_case('i21-cohort', self.root, {
                'execution': 'actual-native-cohort-and-independent-audit', 'status': 'PASS', 'checks': []})

    def test_artifact_limit_is_shared_across_cases(self):
        self.input_path('previous-case.log', 'a'*4096)
        budget = CampaignBudget(8, 1, artifact_bytes=1024)
        receipt = budget.run([sys.executable, '-B', '-c', 'pass'], self.root,
            self.root/'process.log', output=self.root/'resources', artifact_root=self.root)
        self.assertEqual(receipt['status'], 'ERROR')
        self.assertIn('artifact byte limit', receipt['reason'])

    def test_independent_mpi_parent_comparison(self):
        # Small deterministic administrative histories test order independence
        # and distinguish an omitted parent from a changed physical channel.
        hit = {'original_parent_id':'2', 'status':'REACTED', 'channel_id':'declared',
               'relative_event_time_s':'10'}
        survivor = {'original_parent_id':'4', 'status':'SURVIVED', 'channel_id':'',
                    'relative_event_time_s':''}
        result = compare_parent_layouts([hit,survivor], [survivor,dict(hit,relative_event_time_s='10.5')])
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['maximum_serialized_time_difference_s'], .5)
        result = compare_parent_layouts([hit,survivor], [dict(hit,channel_id='other')])
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['changed_fate_or_channel_parent_ids'], ['2'])
        self.assertEqual(result['unpaired_parent_ids'], ['4'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
