#!/usr/bin/env python3
"""Negative controls for the stage-6 audit, using a completed real campaign.

Only isolated copies of synthetic-verification output are deliberately damaged.
These are checker tests, not simulation reruns or physical reference data. The
actual process receipts, immutable executable and campaign output stay intact.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET

import verify_standalone_application as audit


class CorruptedCampaign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.jobs = json.loads((cls.root/'launch-plan.json').read_text())['jobs']
        cls.job = next(job for job in cls.jobs if job['mode'] == 'fte-dmumu')
        cls.hashes = json.loads((cls.root/'build-sha256.json').read_text())
        # First establish that the unmodified actual case passes. Otherwise a
        # negative control could succeed because of an unrelated missing input.
        audit.verify_mode(cls.root, 'fte-dmumu', cls.job, cls.hashes)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='stage6-audit-negative-')
        self.root_copy = Path(self.temporary.name)
        self.case = self.root_copy/'fte-dmumu'
        for name in ['process-result.json', 'application.rank0.log',
                     'PT/plots/amps.FieldLines.out=0.dat', 'PT/plots/amps.FieldLines.out=0.dat.sampling']:
            destination = self.case/name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(self.root/'fte-dmumu'/name, destination)

    def tearDown(self):
        self.temporary.cleanup()

    def check(self, message, hashes=None):
        with self.assertRaisesRegex(audit.GateFailure, message):
            audit.verify_mode(self.root_copy, 'fte-dmumu', self.job,
                              self.hashes if hashes is None else hashes)

    def change_log(self, transform):
        path = self.case/'application.rank0.log'
        path.write_text(transform(path.read_text()))

    def test_wrong_physical_endpoint(self):
        self.change_log(lambda text: re.sub(r'^(REDUCED_SOURCE_LEDGER )epoch_s=1.20000000000000000e\+05',
            r'\g<1>epoch_s=1.19999000000000000e+05', text, count=1, flags=re.M))
        self.check('physical endpoint changed')

    def test_shortened_sampling_window(self):
        self.change_log(lambda text: text.replace('current sample length:400', 'current sample length:399', 1))
        self.check('sampling cadence changed')

    def test_missing_completed_step(self):
        self.change_log(lambda text: re.sub(r'^REDUCED_SOURCE_LEDGER[^\n]*\n', '', text, count=1, flags=re.M))
        self.check('incomplete physical campaign')

    def test_birth_energy_in_empty_population(self):
        # The first sampled population is empty. A one-joule retained birth
        # measure must fail with zero allowed residual, not pass via a floor.
        self.change_log(lambda text: re.sub(r'retained_birth_j=[^\s]+', 'retained_birth_j=1', text, count=1))
        self.check('immutable birth energy does not close')

    def test_segment_sampling_dt(self):
        path = self.case/'PT/plots/amps.FieldLines.out=0.dat.sampling'
        lines = path.read_text().splitlines()
        index = next(i for i, line in enumerate(lines) if line and not line.startswith('#'))
        fields = lines[index].split(); fields[6] = '299'; lines[index] = ' '.join(fields)
        path.write_text('\n'.join(lines)+'\n')
        self.check('physical sampling dt changed')

    def test_inactive_shock_title(self):
        path = self.case/'PT/plots/amps.FieldLines.out=0.dat'
        path.write_text(path.read_text().replace('shock model=reduced-SSE/Parker;', 'shock model=inactive-legacy;'))
        self.check('inactive shock authority')

    def test_nonfinite_active_field(self):
        path = self.case/'PT/plots/amps.FieldLines.out=0.dat'
        lines = path.read_text().splitlines()
        names = re.findall(r'"([^"]*)"', next(line for line in lines if line.startswith('VARIABLES=')))
        index = next(i for i, line in enumerate(lines) if line.strip() and line.lstrip()[0] in '-+.0123456789')
        fields = lines[index].split(); fields[names.index('Bz [nT]')] = 'nan'; lines[index] = ' '.join(fields)
        path.write_text('\n'.join(lines)+'\n')
        self.check('nonfinite active background primitive')

    def test_changed_executable_identity(self):
        changed = dict(self.hashes)
        changed[str(Path(self.job['cwd']).parent/'amps')] = '0'*64
        self.check('executable or physical input changed', changed)

    def test_changed_source_identity(self):
        # Old successful output cannot validate a subsequently edited checkout.
        # Corrupt only the COPY of the expected digest, never the source file.
        changed = dict(self.hashes)
        source = next(name for name in changed if name.endswith('/srcSEP/parker_mover.cpp'))
        changed[source] = '0'*64
        self.check('recorded source or build input changed', changed)

    def test_missing_process_receipt(self):
        (self.case/'process-result.json').unlink()
        with self.assertRaises(FileNotFoundError):
            audit.verify_mode(self.root_copy, 'fte-dmumu', self.job, self.hashes)

    def test_failed_process_receipt(self):
        # Deliberate corruption of the copied receipt is a negative fixture;
        # this is never saved as the original process's actual return code.
        path = self.case/'process-result.json'
        receipt = json.loads(path.read_text()); receipt['exit_code'] = 1
        path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(RuntimeError, 'application exited with code 1'):
            audit.verify_mode(self.root_copy, 'fte-dmumu', self.job, self.hashes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-dir', type=Path, required=True)
    args = parser.parse_args()
    CorruptedCampaign.root = args.campaign_dir.resolve()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CorruptedCampaign)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    failures = dict((test.id().split('.')[-1], message) for test, message in result.failures+result.errors)
    names = unittest.defaultTestLoader.getTestCaseNames(CorruptedCampaign)
    records = [{'test_id': name, 'status': 'FAIL' if name in failures else 'PASS',
                'message': failures.get(name, 'intentional corrupted output rejected')}
               for name in names] if result.testsRun == len(names) else [
                   {'test_id': 'negative-control-prerequisite', 'status': 'ERROR',
                    'message': str(result.errors)}]
    root = CorruptedCampaign.root
    (root/'negative-control-results.json').write_text(json.dumps({'scope': 'audit controls, not physical simulations',
                                                                'results': records}, indent=2)+'\n')
    xml = ET.Element('testsuite', name='stage6-audit-negative-controls', tests=str(len(records)))
    for row in records:
        case = ET.SubElement(xml, 'testcase', name=row['test_id'])
        if row['status'] != 'PASS':
            ET.SubElement(case, 'failure', message=row['message'])
    ET.ElementTree(xml).write(root/'negative-controls.xml', encoding='unicode')
    return int(not result.wasSuccessful())


if __name__ == '__main__':
    raise SystemExit(main())
