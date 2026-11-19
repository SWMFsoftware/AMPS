"""Real-source lexical controls; no manufactured observation or model truth.

The original author ZIP is the independent byte oracle. Missing retained source
is SKIPPED. These tests do not assign a retention law, error confidence level,
temperature, sample identifier or physical-validation PASS.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'utility/chemical_database'))
import native_surface_sources as native

DEFAULT = ROOT/'test_output/srcMoon/chemical-input-acquisition-20261119/surface/acquisition-20261119T062501Z/raw/As_published_in_Comm_Earth_Env.zip'
SOURCE = Path(os.environ.get('AMPS_SURFACE_SOURCE_ZIP', str(DEFAULT)))


class NativeSurfaceSource(unittest.TestCase):
    requires_authoritative_source = True
    @classmethod
    def setUpClass(cls):
        if not SOURCE.is_file():
            raise unittest.SkipTest('retained authoritative surface ZIP missing; set AMPS_SURFACE_SOURCE_ZIP')

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='native-surface-source-controls-')
        self.work = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_producer_identity_and_every_native_member_byte(self):
        self.assertEqual(SOURCE.stat().st_size, 22673)
        self.assertEqual(hashlib.md5(SOURCE.read_bytes()).hexdigest(), native.SOURCE_MD5)
        records = native.reproduce(SOURCE, self.work/'result')
        with zipfile.ZipFile(SOURCE) as archive:
            self.assertEqual({row['output'] for row in records}, set(archive.namelist()))
            for row in records:
                original = archive.read(row['native_archive_member'])
                self.assertEqual((self.work/'result'/row['output']).read_bytes(), original)
                self.assertEqual(row['output_sha256'], hashlib.sha256(original).hexdigest())

    def test_measured_and_modelled_rows_are_kept_distinct(self):
        records = native.reproduce(SOURCE, self.work/'result')
        measured = [row for row in records if row['source_role']=='measured-sputter-yield']
        models = [row for row in records if row['source_role']=='model-output-sputter-yield']
        self.assertEqual(len(measured), 4)
        self.assertEqual(sum(row['row_count'] for row in measured), 17)
        self.assertEqual(len(models), 14)
        self.assertEqual(sum(row['row_count'] for row in models), 2033)

    def test_no_unmeasured_surface_physics_is_invented(self):
        for row in native.reproduce(SOURCE, self.work/'result'):
            for key in ('uncertainty_confidence_level','temperature_K',
                        'incoming_ion_retention_fraction','outgoing_charge_distribution'):
                self.assertIsNone(row[key])
            self.assertIn('unresolved', row['material_id_status'])
            self.assertEqual(row['status'], 'SOURCE_REPRODUCTION_ONLY_NOT_PHYSICS_QUALIFIED')

    def test_source_bytes_and_mtime_remain_unchanged(self):
        before = (SOURCE.read_bytes(), SOURCE.stat().st_mtime_ns)
        native.reproduce(SOURCE, self.work/'result')
        self.assertEqual(before, (SOURCE.read_bytes(), SOURCE.stat().st_mtime_ns))

    def test_same_size_corruption_refuses_before_output_creation(self):
        damaged = self.work/'changed.zip'
        payload = bytearray(SOURCE.read_bytes()); payload[-1] ^= 1; damaged.write_bytes(payload)
        with self.assertRaisesRegex(ValueError, 'digest changed'):
            native.reproduce(damaged, self.work/'result')
        self.assertFalse((self.work/'result').exists())

    def test_changed_size_refuses_before_output_creation(self):
        damaged = self.work/'truncated.zip'; damaged.write_bytes(SOURCE.read_bytes()[:-1])
        with self.assertRaisesRegex(ValueError, 'byte length changed'):
            native.reproduce(damaged, self.work/'result')
        self.assertFalse((self.work/'result').exists())

    def test_existing_destination_and_its_contents_are_preserved(self):
        destination = self.work/'result'; destination.mkdir(); sentinel = destination/'user'; sentinel.write_text('existing output')
        with self.assertRaises(FileExistsError):
            native.reproduce(SOURCE, destination)
        self.assertEqual(sentinel.read_text(), 'existing output')

    def test_replay_is_deterministic(self):
        first = native.reproduce(SOURCE, self.work/'first')
        second = native.reproduce(SOURCE, self.work/'second')
        self.assertEqual(first, second)
        self.assertEqual((self.work/'first/reproduction.json').read_bytes(),
                         (self.work/'second/reproduction.json').read_bytes())

    def test_cli_missing_source_is_skip_without_outputs(self):
        result = subprocess.run([sys.executable, '-B', str(ROOT/'utility/chemical_database/native_surface_sources.py'),
                                 '--source-zip', str(self.work/'missing.zip'),
                                 '--output-dir', str(self.work/'result')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 77)
        self.assertEqual(json.loads(result.stdout)['status'], 'SKIPPED')
        self.assertFalse((self.work/'result').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
