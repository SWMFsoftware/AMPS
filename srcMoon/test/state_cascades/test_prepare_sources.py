#!/usr/bin/env python3
"""Bounded source-integrity controls using an actual acquired NIST packet.

Negative controls change only copied manifests or malformed metadata fixtures.
Neither the original raw source nor any canonical validation data is edited.
The controls evaluate software/source bookkeeping, never atomic lifetimes or
agreement with observations.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal
import json
from pathlib import Path
import sys
import time
import unittest

try:
    from . import prepare_sources as source
except ImportError:
    import prepare_sources as source


def suite(packet, output):
    class SourceChecks(unittest.TestCase):
        @classmethod
        def setUpClass(cls):
            # Parse the actual complete packet once, preserving the separate
            # offline replay and its native hashes as artifacts of this run.
            cls.processing = output / 'processing'
            cls.processing.mkdir()
            cls.qa = source.audit(packet, cls.processing)
            cls.receipt = source.verify_source_packet(packet)
            cls.levels = source.native_csv(packet / 'raw/oi_levels_complete.csv', source.LEVEL_FIELDS)
            cls.lines = source.native_csv(packet / 'raw/oi_lines_complete.csv', source.LINE_FIELDS)

        def copied_manifest(self, name, change):
            # The raw-directory symlink is read-only by procedure; the only
            # mutated object is a new manifest in this control's own folder.
            path = output / name
            path.mkdir()
            (path / 'raw').symlink_to(packet / 'raw', target_is_directory=True)
            value = json.loads((packet / 'acquisition.json').read_text())
            change(value)
            (path / 'acquisition.json').write_text(json.dumps(value, indent=2) + '\n')
            return path

        def test_complete_native_html_csv_counts(self):
            self.assertEqual(len(self.lines), source.metadata(packet / 'raw/oi_lines_complete.html', 'Lines')['row_count'])
            self.assertEqual(len(self.levels), source.metadata(packet / 'raw/oi_levels_complete.html', 'Levels')['row_count'])

        def test_every_acquired_file_identity_verified(self):
            expected = source.fixed_requests() + source.crosscheck_requests(self.levels)
            self.assertEqual({r['file']: r['url'] for r in self.receipt['requests']}, dict(expected))
            for row in self.receipt['requests']:
                path = packet / 'raw' / row['file']
                self.assertEqual(source.digest(path), row['sha256'])
                self.assertEqual(source.digest(path, 'md5'), row['md5'])

        def test_full_native_values_and_missing_A_preserved(self):
            with (self.processing / 'processed/oi_transitions.csv').open(newline='') as stream:
                processed = list(csv.DictReader(stream))
            self.assertEqual(processed, self.lines)
            self.assertGreater(sum(not r['Aki(s^-1)'] for r in processed), 0)
            self.assertEqual(sum(not r['Aki(s^-1)'] for r in processed), self.qa['complete_export_missing_A_rows'])

        def test_multiple_multipoles_for_same_lower_level_retained(self):
            upper = source.selected_levels(self.levels)['O_I_1D_2']
            exits = source.select_exits(self.lines, upper)
            for angular in ('2', '1'):
                pair = [r for r in exits if r['term_i'] == '3P' and r['J_i'] == angular]
                self.assertEqual({r['Type'] for r in pair}, {'M1', 'E2'})

        def test_selected_source_uncertainties_remain_missing(self):
            for upper in source.selected_levels(self.levels).values():
                self.assertEqual(upper['Uncertainty (cm-1)'], '')
            self.assertIn('No explicit information on uncertainties',
                          source.metadata(packet / 'raw/oi_lines_complete.html', 'Lines')['native_prose'])

        def test_accuracy_bounds_are_source_codes_not_sigma(self):
            codes = source.accuracy_codes(packet / 'raw/lines_help.html')
            self.assertEqual(Decimal(codes['B+']['fractional_bound']), Decimal('.07'))
            self.assertEqual(Decimal(codes['C+']['fractional_bound']), Decimal('.18'))
            self.assertEqual(codes['B+']['relation'], '≤')

        def test_source_integrity_does_not_qualify_lifetime_or_branching(self):
            self.assertEqual(self.qa['status'], 'PASS')
            self.assertEqual(self.qa['radiative_cascade_qualification'], 'SKIPPED / NOT QUALIFIED')
            self.assertFalse(self.qa['physical_activation'])
            for state in self.qa['selected_states']:
                self.assertFalse(state['sum_is_complete_decay_rate'])
                self.assertIsNone(state['total_decay_rate_s_1'])
                self.assertIsNone(state['lifetime_s'])
                self.assertIsNone(state['normalized_branching_probabilities'])

        def test_absent_pair_is_not_a_zero_rate(self):
            with (self.processing / 'processed/selected_lower_level_coverage.csv').open(newline='') as stream:
                absent = [r for r in csv.DictReader(stream) if r['listed_component_count'] == '0']
            self.assertTrue(absent)
            self.assertTrue(all(r['source_status'] == 'NO_STORED_TRANSITION_NOT_A_ZERO_RATE' for r in absent))

        def test_fine_structure_descendants_remain_explicit(self):
            graph = json.loads((self.processing / 'processed/stored_cascade_graph.json').read_text())
            nodes = [value['source_level'] for value in graph['nodes'].values()]
            self.assertEqual({(r['Term'], r['J']) for r in nodes},
                             {('1D', '2'), ('1S', '0'), ('3P', '2'), ('3P', '1'), ('3P', '0')})
            fine = [r for r in graph['components'] if r['term_k'] == '3P']
            self.assertEqual(len(fine), 4)
            self.assertTrue(all(r['Aki(s^-1)'] for r in fine))
            self.assertFalse(graph['complete_physical_cascade'])
            self.assertIsNone(graph['physical_lifetime_s'])
            self.assertTrue(all(value['physical_lifetime_s'] is None for value in graph['nodes'].values()))

        def test_wrong_expected_sha_is_refused(self):
            path = self.copied_manifest('wrong-sha-control',
                lambda value: value['requests'][0].update(sha256='0' * 64))
            with self.assertRaisesRegex(ValueError, 'hash/size mismatch'):
                source.verify_source_packet(path)

        def test_wrong_expected_md5_is_refused(self):
            path = self.copied_manifest('wrong-md5-control',
                lambda value: value['requests'][0].update(md5='0' * 32))
            with self.assertRaisesRegex(ValueError, 'hash/size mismatch'):
                source.verify_source_packet(path)

        def test_filtered_query_is_refused(self):
            path = self.copied_manifest('filtered-query-control',
                lambda value: value['requests'][0].update(url=value['requests'][0]['url'] + '?min_str=1'))
            with self.assertRaisesRegex(ValueError, 'request set differs'):
                source.verify_source_packet(path)

        def test_missing_request_is_refused(self):
            path = self.copied_manifest('missing-request-control',
                lambda value: value['requests'].pop())
            with self.assertRaisesRegex(ValueError, 'request set differs'):
                source.verify_source_packet(path)

        def test_duplicate_manifest_key_is_refused(self):
            path = output / 'duplicate-json-control.json'
            path.write_text('{"status":"COMPLETE","status":"ERROR"}\n')
            with self.assertRaisesRegex(ValueError, 'duplicate JSON field'):
                source.strict_json(path)

        def test_unknown_native_schema_is_refused(self):
            path = output / 'malformed-header-control.csv'
            path.write_text('NOT_A_NIST_EXPORT\n')
            with self.assertRaisesRegex(ValueError, 'schema changed'):
                source.native_csv(path, source.LINE_FIELDS)

    return unittest.defaultTestLoader.loadTestsFromTestCase(SourceChecks)


def main(argv=None, default_packet=None, entrypoint=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-packet', type=Path, required=default_packet is None,
                        default=default_packet)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = dict(schema='srcMoon-NIST-OI-source-controls-v1', status='ERROR',
        scope='source/parser/provenance controls only',
        formal_cascade_status='SKIPPED / NOT QUALIFIED',
        runner=str(entrypoint or Path(__file__).resolve()),
        log=str(output / 'checks.log'))
    packet = args.source_packet.resolve()
    before = {}
    script_hashes = {str(Path(path).resolve()): source.digest(path)
                     for path in (source.__file__, __file__, entrypoint) if path}
    try:
        if not packet.is_dir():
            result.update(status='SKIPPED', reason='actual acquired source packet unavailable',
                          counts=dict(PASS=0, FAIL=0, ERROR=0, SKIPPED=1))
            (output / 'checks.log').write_text('SKIPPED: actual acquired source packet unavailable; no controls executed.\n')
        else:
            before = {str(p.relative_to(packet)): source.digest(p)
                      for p in (packet / 'raw').iterdir() if p.is_file()}
            checks = suite(packet, output)
            ids = [case.id() for case in checks]
            with (output / 'checks.log').open('w') as stream:
                actual = unittest.TextTestRunner(stream=stream, verbosity=2).run(checks)
            outcomes = dict.fromkeys(ids, 'PASS')
            for status, records in [('SKIPPED', actual.skipped), ('FAIL', actual.failures), ('ERROR', actual.errors)]:
                for case, _detail in records:
                    if case.id() in outcomes:
                        outcomes[case.id()] = status
                    else:
                        # A class-setup failure prevents all declared methods
                        # from executing. Retain that one actual error and
                        # explicitly classify the unexecuted methods SKIPPED.
                        outcomes = {identity: 'SKIPPED' for identity in ids}
                        outcomes[case.id()] = status
            result.update(status='ERROR' if actual.errors else 'FAIL' if actual.failures else 'SKIPPED' if actual.skipped else 'PASS',
                counts={status: sum(s == status for s in outcomes.values()) for status in ('PASS', 'FAIL', 'ERROR', 'SKIPPED')},
                outcomes=outcomes,
                errors=[dict(test=str(case), detail=detail) for case, detail in actual.errors],
                failures=[dict(test=str(case), detail=detail) for case, detail in actual.failures])
            after = {str(p.relative_to(packet)): source.digest(p)
                     for p in (packet / 'raw').iterdir() if p.is_file()}
            if before != after:
                raise RuntimeError('original source bytes changed during controls')
        if {path: source.digest(path) for path in script_hashes} != script_hashes:
            raise RuntimeError('source/control scripts changed during execution')
    except (OSError, ValueError, KeyError, RuntimeError, csv.Error) as error:
        result.update(status='ERROR', reason=str(error))
        if not result.get('counts'):
            result['counts'] = dict(PASS=0, FAIL=0, ERROR=0, SKIPPED=0)
        result['counts']['ERROR'] += 1
        with (output / 'checks.log').open('a') as stream:
            stream.write('ERROR: ' + str(error) + '\n')
    # Missing/corrupt packages still need a manifest and an actual result, so
    # a global runner never has to infer scientific status from a traceback.
    acquisition = dict(path=str(packet / 'acquisition.json'))
    try:
        acquisition['sha256'] = source.digest(packet / 'acquisition.json')
    except OSError as error:
        acquisition['unavailable'] = str(error)
    (output / 'run_manifest.json').write_text(json.dumps(dict(
        source_packet=str(packet), original_raw_hashes=before,
        source_hashes=script_hashes, acquisition=acquisition,
        python=sys.version, canonical_data_modified=False), indent=2) + '\n')
    result['elapsed_seconds'] = time.monotonic() - started
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return dict(PASS=0, FAIL=1, ERROR=2, SKIPPED=77)[result['status']]


if __name__ == '__main__':
    raise SystemExit(main())
