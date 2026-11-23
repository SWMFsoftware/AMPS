#!/usr/bin/env python3
"""Software negative controls for immutable full-application launch evidence."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_standalone_application as runner
import run_backend_admission_tests as backend


class ApplicationPlan(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.input = self.root/'source'
        self.input.write_text('software identity fixture\n')
        runner.write_json(self.root/'build-sha256.json', {str(self.input): runner.digest(self.input)})
        self.plan = {'schema': 'srcsep-full-application-plan-v1', 'jobs': runner.jobs_for(self.root)}
        self.save()

    def save(self):
        runner.write_json(self.root/'launch-plan.json', self.plan)

    def test_exact_complete_catalog(self):
        self.assertEqual([j['mode'] for j in runner.check_prepared(self.root)['jobs']], list(runner.audit.MODES))

    def test_missing_mover_rejects(self):
        self.plan['jobs'].pop(); self.save()
        with self.assertRaisesRegex(ValueError, 'catalog'):
            runner.check_prepared(self.root)

    def test_short_run_rejects(self):
        self.plan['jobs'][0]['command'] = self.plan['jobs'][0]['command'].replace('400', '4')
        self.save()
        with self.assertRaisesRegex(ValueError, 'catalog'):
            runner.check_prepared(self.root)

    def test_changed_source_rejects(self):
        self.input.write_text('changed software identity\n')
        with self.assertRaisesRegex(ValueError, 'changed'):
            runner.check_prepared(self.root)

    def test_foreign_receipt_rejects(self):
        job = self.plan['jobs'][0]; directory = Path(job['cwd']); directory.mkdir()
        runner.write_json(directory/'process-result.json', {'exit_code': 0, 'argv': ['/bin/true'], 'cwd': job['cwd']})
        with self.assertRaisesRegex(ValueError, 'another application invocation'):
            runner.collect(self.root)

    def test_preflight_errors_preserve_all_movers(self):
        runner.report_error(self.root, 'missing configured native build')
        rows = json.loads((self.root/'verified-results.json').read_text())['results']
        self.assertEqual([row['test_id'] for row in rows], ['SAPP01', 'SAPP02', 'SAPP03'])
        self.assertTrue(all(row['status'] == 'ERROR' for row in rows))

    def test_paths_are_shell_quoted(self):
        import shlex
        root = self.root/"spaces ' and $(literal)"
        job = runner.jobs_for(root)[0]
        self.assertEqual(shlex.split(job['command']), job['argv'])
        child = shlex.split(job['argv'][-1])
        self.assertEqual(child[3], str(root/'amps'))

    def test_preparation_is_not_execution(self):
        fresh = self.root/'fresh'
        output = io.StringIO()
        # Expected child reports are fixture data, not this unittest stream's
        # verdicts. Capture and assert them so the routine reader sees one
        # unambiguous actual result for each negative control.
        with contextlib.redirect_stdout(output), \
                patch.object(runner, 'execute', side_effect=AssertionError('must not launch')), \
                patch.object(runner, 'prepare'), patch.object(runner, 'check_prepared'), \
                patch('sys.argv', ['runner', '--build-root', str(self.root),
                    '--output-dir', str(fresh), '--prepare-only']):
            self.assertEqual(runner.main(), 0)
        self.assertIn('PREPARED, not executed:', output.getvalue())
        self.assertFalse((fresh/'verified-results.json').exists())

    def test_existing_application_evidence_is_immutable(self):
        marker = self.root/'verified-results.json'
        marker.write_text('retained prior qualification\n')
        output = io.StringIO()
        with contextlib.redirect_stderr(output), \
                patch.object(runner, 'prepare', side_effect=AssertionError('must not prepare over evidence')), \
                patch('sys.argv', ['runner', '--build-root', str(self.root), '--output-dir', str(self.root)]):
            self.assertEqual(runner.main(), 1)
        self.assertIn('choose a fresh application evidence directory:', output.getvalue())
        self.assertEqual(marker.read_text(), 'retained prior qualification\n')

    def test_existing_backend_evidence_is_immutable(self):
        marker = self.root/'summary.json'
        marker.write_text('retained prior backend evidence\n')
        output = io.StringIO()
        with contextlib.redirect_stderr(output), \
                patch.object(backend.subprocess, 'run', side_effect=AssertionError('must not launch over evidence')), \
                patch('sys.argv', ['runner', '--native', '--amps', '/bin/true', '--output-dir', str(self.root)]):
            self.assertEqual(backend.main(), 1)
        self.assertIn('choose a fresh backend evidence directory:', output.getvalue())
        self.assertEqual(marker.read_text(), 'retained prior backend evidence\n')

    def test_malformed_backend_plan_reports_each_required_error(self):
        # /bin/true is only a readable executable identity in this software
        # negative control. The malformed plan must reject BEFORE any MPI
        # launch, and this check makes no native transport qualification claim.
        (self.root/'prepared.json').write_text('{}\n')
        output = io.StringIO()
        with contextlib.redirect_stdout(output), \
                patch.object(backend.subprocess, 'run', side_effect=AssertionError('must not launch malformed plan')), \
                patch('sys.argv', ['runner', '--native', '--collect', '--amps', '/bin/true', '--output-dir', str(self.root)]):
            self.assertEqual(backend.main(), 1)
        rows = json.loads((self.root/'summary.json').read_text())['results']
        self.assertEqual([r['test_id'] for r in rows], ['SBACK06', 'SBACK07', 'SBACK08'])
        self.assertTrue(all(r['status']=='ERROR' for r in rows))
        self.assertEqual(output.getvalue().count("] ERROR 'hashes'"), 3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
