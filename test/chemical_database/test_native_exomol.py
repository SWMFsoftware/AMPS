"""Independent native-token checks; measurement data are never manufactured.

Small parser controls are explicitly administrative. Optional real-source tests
consume only the exact packet selected by AMPS_EXOMOL_PACKET; no newest-packet
search, source interpolation or physical-provider activation occurs.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("tested_native_exomol", ROOT/'utility/chemical_database/native_exomol.py')
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


class ExoMolParser(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_signed_administrative_tokens_are_not_clipped(self):
        source = self.root/'administrative.photo'
        source.write_text('100.000 -1.25E-9\n100.020 0.000E+00\n100.021 2.5E-9\n')
        result = native.project(source, self.root/'projection.csv', 'photo')
        with (self.root/'projection.csv').open() as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows[1:], [line.split() for line in source.read_text().splitlines()])
        self.assertEqual(result['source_checks']['negative_sample_count'], 1)
        self.assertFalse(result['source_checks']['nonnegative_probability_table'])
        self.assertFalse(result['physical_activation'])

    def test_nonincreasing_native_grid_refused_without_sorting(self):
        source = self.root/'administrative.photo'
        source.write_text('100 1\n99 2\n')
        with self.assertRaisesRegex(ValueError, 'nonincreasing'):
            native.project(source, self.root/'projection.csv', 'photo')
        self.assertFalse((self.root/'projection.csv').exists())

    def test_nonfinite_source_is_not_numeric_input(self):
        source = self.root/'administrative.photo'
        source.write_text('100 NaN\n')
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            native.project(source, self.root/'projection.csv', 'photo')

    def test_existing_output_is_not_overwritten(self):
        source = self.root/'administrative.photo';source.write_text('100 1\n')
        output = self.root/'existing.csv';output.write_text('ORIGINAL')
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            native.project(source, output, 'photo')
        self.assertEqual(output.read_text(), 'ORIGINAL')

    def test_json_max_energy_does_not_impose_a_state_cutoff(self):
        # Administrative schema control: the metadata label is insufficient
        # to identify a physical upper bound. These are not measured levels.
        source = self.root/'administrative.states';source.write_text('1 3.0 1 0\n')
        definition = self.root/'administrative.def.json'
        definition.write_text(json.dumps({'dataset': {'version': 'administrative', 'states': {
            'states_file_fields': [{'name': name} for name in ('ID', 'E', 'gtot', 'J')],
            'number_of_states': 1, 'max_energy': 2.0}}}))
        result = native.project(source, self.root/'states.csv', 'states', definition)
        self.assertEqual(result['source_checks']['maximum_energy_native_cm-1'], '3.0')
        self.assertEqual(result['source_checks']['definition_max_energy_native_cm-1'], '2.0')
        self.assertNotIn('maximum_energy_matches_definition', result['source_checks'])
        self.assertIn('not an inferred maximum', result['source_checks']['definition_max_energy_semantics'])


class ExoMolRealSource(unittest.TestCase):
    requires_authoritative_source = True
    @classmethod
    def setUpClass(cls):
        selected = os.environ.get('AMPS_EXOMOL_PACKET')
        if not selected:
            raise unittest.SkipTest('explicit authoritative ExoMol packet not selected; set AMPS_EXOMOL_PACKET')
        cls.packet = Path(selected).resolve()
        manifest = json.loads((cls.packet/'acquisition_manifest.json').read_text())
        cls.records = {row['path']: row for row in manifest['records']}
        for row in cls.records.values():
            if row['status'] != 'PASS' or not row['complete'] or row['http_status'] != 200:
                raise ValueError('authoritative source packet includes failed/incomplete requests')
            if hashlib.sha256((cls.packet/row['path']).read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('authoritative source packet changed')

    def test_real_water_samples_and_gaps_survive_exactly(self):
        with tempfile.TemporaryDirectory() as directory:
            for name in ('h2o-dtu-423.15K.photo', 'h2o-dtu-573.15K.photo'):
                with self.subTest(source=name):
                    result = native.project(self.packet/name, Path(directory)/(name+'.csv'), 'photo')
                    self.assertEqual(result['row_count'], 7933)
                    self.assertTrue(result['native_token_identity'])
                    self.assertFalse(result['physical_activation'])

    def test_real_oh_states_match_native_definition(self):
        with tempfile.TemporaryDirectory() as directory:
            result = native.project(self.packet/'oh-mythos.states.bz2', Path(directory)/'states.csv', 'states',
                                    self.packet/'oh-mythos.def.json')
            self.assertEqual(result['source_checks']['declared_states'], result['row_count'])
            self.assertEqual(result['row_count'], 10606)
            self.assertEqual(len(result['columns']), 16)
            self.assertIn('not molecular formation energy', ' '.join(result['limitations']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
