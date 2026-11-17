#!/usr/bin/env python3
"""Native-operator runner CONTROL tests, not replacement native physics.

No MPI executable is mocked or launched here. Temporary CSV records below
test reporting/schema/error classes only. The separate supporting runner must
still execute the actual production I15/I20 code in the combined campaign.
No fixture is a lunar observation, numerical solution or validation artifact.
"""
from __future__ import annotations

from contextlib import ExitStack
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import integration_common as integration
import run_m2_operator_checks as checks
from m2_illumination import setting


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def illumination_protocol_rows(root):
    queries = []
    rows = []
    for index in range(8):
        query = {'query_id': 'protocol-query-'+str(index), 'site_id': 'protocol-site-'+str(index),
                 'latitude_deg': 0, 'longitude_east_deg': index, 'time_system': 'TDB', 'epoch': 0}
        queries.append(query)
        rows.append(dict(query, input_epoch='0', input_time_system='TDB', ephemeris_time_tdb_s='0',
                         terrain_visible_area_fraction='.5', earth_visible_area_fraction='1',
                         jointly_visible_area_fraction='.5', projected_incidence='.25',
                         sampled_any_terrain_visible='1', sampled_any_jointly_visible='1', ray_count='64'))
    source_id = 'protocol-source-clock'
    source_utc = '2009-01-24T00:00:00'
    source_query = dict(queries[-1], query_id=source_id, site_id='protocol-source-site',
                        time_system='UTC', epoch=source_utc)
    queries.append(source_query)
    rows.append(dict(rows[-1], query_id=source_id, site_id='protocol-source-site',
                     input_epoch=source_utc, input_time_system='UTC', ephemeris_time_tdb_s='123'))
    write_csv(root/'receipt.txt.I15.rank0000.csv', rows)
    # 123 is a pure report-protocol number, NEVER a claimed UTC conversion.
    # The REAL native runner obtains the conversion from a separate CSPICE
    # query. These temporary fixtures only test identity/error classification.
    (root/'receipt.txt').write_text('test_id=I15\nmpi_ranks=8\nruntime_psd_geometry_hook=1\nsource_begin_epoch_tdb_s=123\n')
    return {'queries': queries, 'ray_count': 64, 'source_epoch_query_id': source_id,
            'source_start_utc': source_utc}, rows


def thermal_protocol_rows(root):
    queries = [];rows = []
    for index in range(8):
        depth = 0 if index%2 == 0 else .02
        query = {'query_id': 'protocol-query-'+str(index), 'site_id': 'protocol-site-'+str(index),
                 'logical_surface_element': index, 'depth_m': depth,
                 'time_support_id': 'protocol-only-not-a-lunar-period'}
        queries.append(query)
        rows.append(dict(query, sample_weight='1', frame='MOON_ME_DE421', model_physical_temperature_K='100',
                         ephemeris_time_tdb_s='123',
                         model_bolometric_brightness_temperature_K=str(100*.95**.25) if depth == 0 else '',
                         absorbed_solar_flux_W_m2='0', model_latitude_deg='0', model_longitude_east_deg='0'))
        # Preserve the actual native spelling; the values exercise report
        # protocol only, never pretend to be output of a thermal solution.
        rows[-1]['absorbed_solar_flux_W_m-2'] = rows[-1].pop('absorbed_solar_flux_W_m2')
        # Report-protocol fixture values only: actual native law evaluation
        # remains a separate compiled run. Surface/depth kind must survive
        # serialization even when all checkpoint writes are disabled.
        rows[-1].update(surface_retention_timescale_s='2550240' if depth == 0 else '',
                        na_surface_residence_time_s=str(math.exp(1.85*1.602176634e-19/(1.380649e-23*100))/1e13) if depth == 0 else '',
                        na_surface_thermal_retention='1' if depth == 0 else '',
                        ar_surface_residence_time_s=str(6.62607015e-34/(1.380649e-23*100)*math.exp(28000/(8.31446261815324*100))) if depth == 0 else '',
                        ar_surface_thermal_retention='0' if depth == 0 else '')
    write_csv(root/'receipt.txt.I20_thermal_model.csv', rows)
    (root/'receipt.txt').write_text('test_id=I20\nmpi_ranks=8\nthermal_model_rows=8\nthermal_spinup_completed_cycles=1\n'
        'surface_kinetics_export_schema=native-species-residence-retention-v1\n'
        'thermal_initialization_mode=periodic-spinup\nthermal_initialization_et_s=123\n'
        'thermal_spinup_end_tdb_s=123\ninitial_thermal_surface_elements=8\n'
        'initial_thermal_layer_count=2\ninitial_thermal_state_file='+str(root/'receipt.txt.I20_initial_thermal_state.csv')+'\n')
    (root/'application.in').write_text('thermal_infrared_emissivity = 0.95\n'
        'cold_trap_retention_timescale_s = 2550240\nsodium_activation_energy_ev = 1.85\n'
        'sodium_attempt_frequency_hz = 10000000000000\nargon_desorption_energy_j_mol = 28000\n'
        'thermal_layer_count = 2\n')
    # Temporary model-state protocol records, never a claimed heat solution
    # or a lunar observation. The separate MPI fixture writes actual columns.
    initial_rows=[{'logical_surface_element':element,'depth_index':index,
                   'depth_m':depth,'kinetic_temperature_K':100}
                  for element in range(8) for index,depth in ((-1,0),(0,.01),(1,.03))]
    write_csv(root/'receipt.txt.I20_initial_thermal_state.csv',initial_rows)
    return {'queries': queries}, rows


class ManifestPolicy(unittest.TestCase):
    def fixture(self, root, stack):
        # Actual hash operations run on tiny protocol bytes. No kernel, DEM
        # or observed value is fabricated and used as physical evidence.
        executable = root/'not-an-executable';executable.write_text('PROTOCOL ONLY')
        kernel = root/'fixture.tls';kernel.write_text('HASH FIXTURE, NOT A SPICE KERNEL')
        image = root/'fixture.img';image.write_text('HASH FIXTURE, NOT A DEM')
        label = root/'fixture.lbl';label.write_text('HASH FIXTURE, NOT A PDS LABEL')
        application = root/'application.in'
        application.write_text('surface_geometry = lola\nlola_image_file = '+str(image)+
                               '\nlola_label_file = '+str(label)+'\n')
        stack.enter_context(patch.object(integration, 'SPICE_ROOT', root))
        stack.enter_context(patch.object(integration, 'SPICE_FILES', ['fixture.tls']))
        stack.enter_context(patch.object(integration, 'git_value', return_value='PROTOCOL-IDENTITY-ONLY'))
        stack.enter_context(patch.object(integration.subprocess, 'run', return_value=SimpleNamespace(stdout='fixture compiler\n')))
        authority = stack.enter_context(patch.object(integration, 'dataset_package_preflight',
            return_value=(False, 'PROTOCOL D02 PARTIAL', {})))
        return executable, application, kernel, image, label, authority

    def test_formal_default_still_refuses_partial_d02(self):
        for test in ('I15', 'I20'):
            with self.subTest(test=test), tempfile.TemporaryDirectory(prefix='m2-operator-metadata-fixture-') as temp, ExitStack() as stack:
                root = Path(temp)
                executable, application, _, _, _, authority = self.fixture(root, stack)
                with self.assertRaises(integration.MissingPrerequisite):
                    integration.manifest(test, executable, [], 1, application, 8)
                authority.assert_called_once_with('D02')

    def test_support_omits_unconsumed_observations_keeps_consumed_hashes(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-metadata-fixture-') as temp, ExitStack() as stack:
            root = Path(temp)
            executable, application, kernel, image, label, authority = self.fixture(root, stack)
            result = integration.manifest('I15', executable, [], 1, application, 8,
                                          require_scientific_data=False)
            authority.assert_not_called()
            self.assertEqual(result['test_id'], 'I15')
            self.assertEqual(result['external_data']['fixture.tls'], integration.sha256(kernel))
            for path in (image, label):
                self.assertEqual(result['external_data'][str(path)], integration.sha256(path))
            self.assertEqual(result['inputs'][str(application)], integration.sha256(application))


class FixturePreparation(unittest.TestCase):
    def prepare(self, root):
        def manifest(test, amps, command, seed, application, ranks, **kwargs):
            self.assertFalse(kwargs['require_scientific_data'])
            return {'test_id': test, 'inputs': {}, 'scope': 'PREPARATION PROTOCOL ONLY'}
        with patch.object(checks, 'manifest', side_effect=manifest):
            return checks.prepare_jobs(root/'not-launched', checks.REPO/'srcMoon/examples/m2_development.in', root, 8)

    def test_query_inputs_strip_i19_only_controls_and_keep_formal_skipped(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-preparation-fixture-') as temp:
            jobs, contracts = self.prepare(Path(temp))
            for test, job in jobs.items():
                application = job['directory']/'application.in'
                self.assertNotIn('m2_verification_case =', application.read_text())
                self.assertEqual(setting(application, 'illumination_disc_radial_points'), '4')
                self.assertEqual(setting(application, 'illumination_disc_azimuth_points'), '16')
                identity = json.loads((job['directory']/'run_manifest.json').read_text())
                self.assertEqual(identity['formal_test_status'], 'SKIPPED')
                self.assertEqual(identity['execution_scope'], 'supporting-native-operator-wiring')
                self.assertEqual(len(contracts[test]['queries']), 9 if test == 'I15' else 8)
                if test == 'I15':
                    self.assertEqual(sum(row['time_system'] == 'TDB' and row['epoch'] == 0
                                         for row in contracts[test]['queries']), 8)
                    self.assertEqual(contracts[test]['source_start_utc'],
                                     setting(checks.REPO/'srcMoon/examples/m2_development.in', 'simulation_start_utc'))

    def test_millisecond_replay_preserves_thermal_convergence_limits(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-preparation-fixture-') as temp:
            jobs, contracts = self.prepare(Path(temp))
            application = jobs['I20']['directory']/'application.in'
            self.assertEqual(setting(application, 'thermal_time_step_s'), '0.001')
            self.assertEqual(float(setting(application, 'thermal_spinup_temperature_tolerance_k')), .01)
            self.assertEqual(float(setting(application, 'thermal_spinup_energy_tolerance_j_m2')), 10000)
            start = checks._utc(setting(application, 'thermal_spinup_cycle_start_utc'))
            end = checks._utc(setting(application, 'thermal_spinup_cycle_end_utc'))
            self.assertEqual((end-start).total_seconds(), .001)
            self.assertEqual(setting(application,'simulation_start_utc'),setting(application,'thermal_spinup_cycle_end_utc'))
            self.assertEqual(setting(application,'thermal_initialization_mode'),'periodic-spinup')
            self.assertTrue(contracts['I20']['no_cycle_limit_or_threshold_relaxation'])


class ReportingSemantics(unittest.TestCase):
    def test_unapplied_initial_state_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,_=thermal_protocol_rows(root)
            path=root/'receipt.txt'
            path.write_text(path.read_text().replace('thermal_initialization_mode=periodic-spinup',
                'thermal_initialization_mode=declared-uniform'))
            with self.assertRaises(checks.FixtureFailure):checks.evaluate_thermal(root,contract)

    def test_recorded_state_must_match_actual_t0_query(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,_=thermal_protocol_rows(root)
            path=root/'receipt.txt.I20_initial_thermal_state.csv'
            rows=checks._rows(path);rows[0]['kinetic_temperature_K']='101'
            write_csv(path,rows)
            with self.assertRaises(checks.FixtureFailure):checks.evaluate_thermal(root,contract)

    def test_missing_full_grid_initial_state_is_error(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,_=thermal_protocol_rows(root)
            path=root/'receipt.txt.I20_initial_thermal_state.csv'
            write_csv(path,checks._rows(path)[:-1])
            with self.assertRaises(ValueError):checks.evaluate_thermal(root,contract)

    def test_surface_kinetics_independent_of_checkpoint_writes(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,_=thermal_protocol_rows(root)
            self.assertEqual(checks.evaluate_thermal(root,contract)['native_thermal_rows'],8)

    def test_resolved_wrong_residence_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,rows=thermal_protocol_rows(root)
            rows[0]['na_surface_residence_time_s']='1'
            write_csv(root/'receipt.txt.I20_thermal_model.csv',rows)
            with self.assertRaises(checks.FixtureFailure):checks.evaluate_thermal(root,contract)

    def test_retention_classification_miss_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,rows=thermal_protocol_rows(root)
            rows[0]['na_surface_thermal_retention']='0'
            write_csv(root/'receipt.txt.I20_thermal_model.csv',rows)
            with self.assertRaises(checks.FixtureFailure):checks.evaluate_thermal(root,contract)

    def test_subsurface_retention_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root=Path(temp);contract,rows=thermal_protocol_rows(root)
            rows[1]['ar_surface_thermal_retention']='0'
            write_csv(root/'receipt.txt.I20_thermal_model.csv',rows)
            with self.assertRaises(checks.FixtureFailure):checks.evaluate_thermal(root,contract)

    def test_installed_provider_and_native_source_clock_are_evaluated(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, _ = illumination_protocol_rows(root)
            result = checks.evaluate_illumination(root, contract)
            self.assertEqual(result['runtime_psd_geometry_hook'], 1)
            self.assertEqual(result['source_begin_epoch_tdb_s'], 123)
            self.assertEqual(result['expected_source_begin_epoch_tdb_s'], 123)

    def test_uninstalled_provider_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, _ = illumination_protocol_rows(root)
            (root/'receipt.txt').write_text('runtime_psd_geometry_hook=0\nsource_begin_epoch_tdb_s=123\n')
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_illumination(root, contract)

    def test_wrong_provider_start_clock_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, _ = illumination_protocol_rows(root)
            (root/'receipt.txt').write_text('runtime_psd_geometry_hook=1\nsource_begin_epoch_tdb_s=0\n')
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_illumination(root, contract)

    def test_nonfinite_provider_clock_is_report_error(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, _ = illumination_protocol_rows(root)
            (root/'receipt.txt').write_text('runtime_psd_geometry_hook=1\nsource_begin_epoch_tdb_s=NaN\n')
            with self.assertRaises(ValueError):
                checks.evaluate_illumination(root, contract)

    def test_native_epoch_echo_mismatch_is_identity_error(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = illumination_protocol_rows(root)
            rows[0]['input_time_system'] = 'UTC'
            write_csv(root/'receipt.txt.I15.rank0000.csv', rows)
            with self.assertRaises(ValueError):
                checks.evaluate_illumination(root, contract)

    def test_resolved_probability_bound_miss_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = illumination_protocol_rows(root)
            rows[0]['jointly_visible_area_fraction'] = '.75'
            write_csv(root/'receipt.txt.I15.rank0000.csv', rows)
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_illumination(root, contract)

    def test_resolved_occurrence_identity_miss_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = illumination_protocol_rows(root)
            rows[0]['sampled_any_jointly_visible'] = '0'
            write_csv(root/'receipt.txt.I15.rank0000.csv', rows)
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_illumination(root, contract)

    def test_nonbinary_output_is_report_error(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = illumination_protocol_rows(root)
            rows[0]['sampled_any_jointly_visible'] = '2'
            write_csv(root/'receipt.txt.I15.rank0000.csv', rows)
            with self.assertRaises(ValueError):
                checks.evaluate_illumination(root, contract)

    def test_subsurface_brightness_misclassification_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = thermal_protocol_rows(root)
            rows[1]['model_bolometric_brightness_temperature_K'] = '100'
            write_csv(root/'receipt.txt.I20_thermal_model.csv', rows)
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_thermal(root, contract)

    def test_resolved_grey_operator_miss_is_failure(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp);contract, rows = thermal_protocol_rows(root)
            rows[0]['model_bolometric_brightness_temperature_K'] = '1'
            write_csv(root/'receipt.txt.I20_thermal_model.csv', rows)
            with self.assertRaises(checks.FixtureFailure):
                checks.evaluate_thermal(root, contract)

    def test_missing_native_dependency_is_skipped_not_false_pass(self):
        with tempfile.TemporaryDirectory(prefix='m2-operator-report-fixture-') as temp:
            root = Path(temp)
            result = checks.run_checks(root/'absent-amps', root/'absent-input', root/'out')
            self.assertEqual(result['status'], 'SKIPPED')
            self.assertEqual(result['counts']['PASS'], 0)
            self.assertEqual(result['formal_i15_status'], 'SKIPPED')
            self.assertEqual(result['formal_i20_status'], 'SKIPPED')


if __name__ == '__main__':
    unittest.main(verbosity=2)
