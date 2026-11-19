"""Production-backed Apollo lexical checks, never synthetic physical evidence.

Pure format fixtures contain named strings/empty cells, not observations.
The separately marked authoritative class consumes only the fixed, real
Jones2020 deposit and skips when that dependency is absent. All outputs live
in owned temporary directories; source bytes and mtimes must stay unchanged.
"""
import csv
import hashlib
import importlib.util
import io
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('native_apollo_tpd', ROOT / 'utility/chemical_database/native_apollo_tpd.py')
N = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(N)
DEFAULT_SOURCE = ROOT / 'test_output/srcMoon/m3-physics-literature-20261119T072531Z/surface/jones2020-native'


def rows(payload):
    return list(csv.DictReader(io.StringIO(payload.decode('utf-8'), newline='')))


def independent_original(payload):
    """Simple independent byte identity oracle, not an adsorption evaluator."""
    data = rows(payload)
    result = []
    # The native output is ordered. One linear traversal avoids quadratic
    # rescanning of real large files in the global development test routine.
    for _, group in itertools.groupby(data, key=lambda row: row['source_line_1based']):
        parts = list(group)
        result.append('\t'.join(r['native_token'] for r in parts))
        result.append({'LF': '\n', 'CRLF': '\r\n', 'CR': '\r', 'NONE': ''}[parts[0]['native_line_ending']])
    return ''.join(result).encode('utf-8')


class AdministrativeLexicalControls(unittest.TestCase):
    """These strings test the public lexical boundary, not laboratory physics."""
    def test_independent_roundtrip_all_endings_and_empty_tokens(self):
        raw = b'header\tmissing header\r\n\t-2.500E-03\t\r\n\nlast\ttoken\rfinal\t'
        projected = N.lexical_csv(raw)
        self.assertEqual(independent_original(projected), raw)
        self.assertEqual(N.restore_native(projected), raw)
        self.assertEqual([r['native_token'] for r in rows(projected) if r['source_line_1based'] == '2'], ['', '-2.500E-03', ''])

    def test_csv_special_characters_and_unicode_are_not_altered(self):
        raw = 'original,"label"\tΩ\vnot-a-line\n'.encode('utf-8')
        projected = N.lexical_csv(raw)
        self.assertEqual(independent_original(projected), raw)
        self.assertEqual(N.restore_native(projected), raw)
        self.assertEqual(len(rows(projected)), 2)

    def test_empty_source_has_no_invented_line(self):
        self.assertEqual(rows(N.lexical_csv(b'')), [])
        self.assertEqual(N.restore_native(N.lexical_csv(b'')), b'')

    def test_no_missing_header_field_is_fabricated(self):
        projected = rows(N.lexical_csv(b'one\ttwo\n1\t2\t3\n'))
        self.assertEqual([r['native_token'] for r in projected if r['line_role'] == 'native-header'], ['one', 'two'])
        self.assertEqual([r['source_field_1based'] for r in projected if r['source_line_1based'] == '2'], ['1', '2', '3'])

    def test_malformed_quote_and_nonbytes_are_explicit_errors(self):
        with self.assertRaises(N.NativeApolloError):
            N.restore_native(','.join(N.CSV_COLUMNS).encode()+b'\n1,1,native-header,"unterminated,NONE')
        for method in (N.restore_native, N.lexical_csv):
            with self.subTest(method=method.__name__), self.assertRaises(N.NativeApolloError):
                method('text-instead-of-required-bytes')

    def test_restore_rejects_changed_schema(self):
        with self.assertRaises(N.NativeApolloError):
            N.restore_native(N.lexical_csv(b'header\n').replace(b'native_token', b'invented_value', 1))

    def test_restore_rejects_gaps_duplicate_fields_and_ending_conflicts(self):
        original = N.lexical_csv(b'a\tb\nc\td\n')
        header = list(N.CSV_COLUMNS)
        controls = [('source_field_1based', '1'), ('native_line_ending', 'CR'), ('source_line_1based', '3'), ('line_role', 'native-data')]
        for key, value in controls:
            altered = rows(original)
            altered[1][key] = value
            stream = io.StringIO(newline='')
            writer = csv.DictWriter(stream, fieldnames=header, lineterminator='\n')
            writer.writeheader(); writer.writerows(altered)
            with self.subTest(key=key), self.assertRaises(N.NativeApolloError):
                N.restore_native(stream.getvalue().encode())

    def test_missing_source_cli_is_explicit_skip_without_output(self):
        with tempfile.TemporaryDirectory() as name:
            output = Path(name) / 'new'
            result = subprocess.run([sys.executable, '-B', str(Path(N.__file__)), '--source-dir', str(Path(name) / 'absent'), '--output-dir', str(output)], capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 77, result.stderr)
            self.assertEqual(json.loads(result.stdout)['status'], 'SKIPPED')
            self.assertFalse(output.exists())


class AuthoritativeApolloSource(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        cls.source = Path(os.environ.get('AMPS_APOLLO_TPD_SOURCE_DIR', str(DEFAULT_SOURCE)))
        if not cls.source.is_dir() or any(not (cls.source/name).is_file() for name in N.SOURCE_FILES):
            raise unittest.SkipTest('fixed authoritative Jones2020 source collection unavailable')
        cls.before = {name: ((cls.source/name).read_bytes(), (cls.source/name).stat().st_mtime_ns) for name in N.SOURCE_FILES}
        cls.temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temporary.name) / 'reproduction'
        cls.records = N.reproduce(cls.source, cls.output)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'temporary'):
            cls.temporary.cleanup()

    def clone_source(self, target):
        shutil.copytree(str(self.source), str(target))
        return target

    def test_all_fourteen_native_objects_match_producer_identities_and_copies(self):
        self.assertEqual(len(N.SOURCE_FILES), 14)
        for name, expected in N.SOURCE_FILES.items():
            raw = (self.source/name).read_bytes()
            with self.subTest(name=name):
                self.assertEqual(len(raw), expected['bytes'])
                self.assertEqual(hashlib.md5(raw).hexdigest(), expected['md5'])
                self.assertEqual(hashlib.sha256(raw).hexdigest(), expected['sha256'])
                self.assertEqual((self.output/'native'/name).read_bytes(), raw)

    def test_all_thirteen_tables_reconstruct_independently_to_original_bytes(self):
        self.assertEqual(len(self.records), 13)
        for record in self.records:
            with self.subTest(name=record['source_filename']):
                projected = (self.output/record['output']).read_bytes()
                self.assertEqual(independent_original(projected), (self.source/record['source_filename']).read_bytes())
                self.assertEqual(len(rows(projected)), record['row_count'])
                self.assertEqual(hashlib.sha256(projected).hexdigest(), record['output_sha256'])

    def test_actual_figure3_missing_tab_is_not_repaired(self):
        r = next(r for r in self.records if r['source_filename'] == 'figure3.txt')
        self.assertEqual(r['native_header_tokens'], ['Coverage', 'Energy 2SD'])
        data = rows((self.output/r['output']).read_bytes())
        self.assertEqual(sum(row['source_line_1based'] == '1' for row in data), 2)
        self.assertEqual(sum(row['source_line_1based'] == '2' for row in data), 3)

    def test_blank_signed_and_bound_tokens_are_preserved(self):
        r = next(r for r in self.records if r['source_filename'] == 'figure6.txt')
        data = rows((self.output/r['output']).read_bytes())
        self.assertTrue(any(row['native_token'] == '' for row in data))
        self.assertTrue(any(row['native_token'] == '-1.71E-04' for row in data))
        # Independent temperature axes survive as two separate native header
        # positions; no row-wise pairing/interpolation is assigned by this API.
        header = [row['native_token'] for row in data if row['source_line_1based'] == '1']
        self.assertEqual(header.count('Temperature (K)'), 2)

    def test_measured_inverted_and_model_roles_remain_distinct(self):
        records = {r['source_filename']: r for r in self.records}
        self.assertEqual(records['figure8top.txt']['source_role'], 'lunar-cycle-model-prediction-not-observation')
        self.assertEqual(records['figure3.txt']['source_role'], 'model-inverted-coverage-dependent-activation-energy')
        self.assertEqual(records['figure2.txt']['source_role'], 'measured-TPD-at-distinct-water-exposures')
        self.assertEqual(records['figure4.txt']['source_role'], 'mixed-measured-TPD-and-fitted-component-simulations')
        self.assertTrue(all(r['physical_activation'] is False for r in self.records))
        self.assertTrue(all(r['native_units']['energy_coverage_distribution'] is None for r in self.records))

    def test_source_bytes_and_mtimes_unchanged(self):
        for name, (raw, mtime) in self.before.items():
            with self.subTest(name=name):
                self.assertEqual((self.source/name).read_bytes(), raw)
                self.assertEqual((self.source/name).stat().st_mtime_ns, mtime)

    def test_corruption_and_truncation_refused_before_any_output(self):
        for kind in ('same-size', 'truncated'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as name:
                source = self.clone_source(Path(name)/'source')
                path = source/'figure3.txt'
                raw = path.read_bytes()
                path.write_bytes((b'X'+raw[1:]) if kind == 'same-size' else raw[:-1])
                output = Path(name)/'new'
                with self.assertRaises(N.NativeApolloError):
                    N.reproduce(source, output)
                self.assertFalse(output.exists())

    def test_missing_member_is_unavailable_before_output(self):
        with tempfile.TemporaryDirectory() as name:
            source = self.clone_source(Path(name)/'source')
            (source/'figure3.txt').unlink()
            output = Path(name)/'new'
            with self.assertRaises(FileNotFoundError):
                N.reproduce(source, output)
            self.assertFalse(output.exists())

    def test_extra_member_and_symlink_refused_before_output(self):
        for kind in ('extra', 'symlink'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as name:
                source = self.clone_source(Path(name)/'source')
                if kind == 'extra':
                    (source/'unreviewed.txt').write_text('administrative fixture only')
                else:
                    (source/'figure3.txt').unlink()
                    (source/'figure3.txt').symlink_to(self.source/'figure3.txt')
                output = Path(name)/'new'
                with self.assertRaises(N.NativeApolloError):
                    N.reproduce(source, output)
                self.assertFalse(output.exists())

    def test_existing_output_and_sentinel_are_preserved(self):
        with tempfile.TemporaryDirectory() as name:
            output = Path(name)/'exists'; output.mkdir()
            (output/'sentinel').write_bytes(b'preserve')
            with self.assertRaises(FileExistsError):
                N.reproduce(self.source, output)
            self.assertEqual((output/'sentinel').read_bytes(), b'preserve')

    def test_reproduction_is_deterministic(self):
        with tempfile.TemporaryDirectory() as name:
            target = Path(name)/'replay'
            records = N.reproduce(self.source, target)
            self.assertEqual(records, self.records)
            for record in records:
                self.assertEqual((target/record['output']).read_bytes(), (self.output/record['output']).read_bytes())


if __name__ == '__main__':
    unittest.main(verbosity=2)
