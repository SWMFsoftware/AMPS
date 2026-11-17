#!/usr/bin/env python3
"""I20 infrastructure/row-arithmetic tests; NEVER observational evidence.

All inputs are tiny temporary syntax/arithmetic fixtures, not synthetic lunar
observations. These tests validate metadata isolation, stable preselection,
explicit units/depth/time support and honest failure semantics. No fixture is
written to the canonical data tree or rendered as a scientific comparison.
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import m2_polar as polar


def candidate(index, *, element=7, edge=0):
    return {'site_id': 'metadata-fixture-'+str(index), 'native_product_id': 'FIXTURE-NOT-DATA',
            'native_line': index+1, 'native_sample': 2, 'latitude_deg': '80',
            'longitude_east_deg': str(index), 'hemisphere': 'N',
            'permanent_shadow': '0', 'edge_ambiguous': str(edge),
            'logical_surface_element': str(element), 'd03_quality_flag': 'GOOD',
            'collocation_distance_m': '2'}


def arithmetic_rows():
    sites = {'arithmetic-fixture': {'logical_surface_element': '7'}}
    model, reference, criteria = [], [], {}
    values = {'any_solar_disc_visibility_occurrence_fraction': .25,
              'permanent_shadow_classification': 0,
              'annual_maximum_surface_temperature_K': 100,
              'annual_average_temperature_at_2cm_K': 80}
    for quantity, (units, depth) in polar.QUANTITIES.items():
        identity = {'site_id': 'arithmetic-fixture', 'quantity': quantity,
                    'units': units, 'depth_m': str(depth), 'frame': polar.FRAME,
                    'time_support_id': 'explicit-arithmetic-support',
                    'temperature_kind': 'physical_temperature_K' if units == 'K' else ''}
        reference.append(dict(identity, reference_value=str(values[quantity]), uncertainty_value='',
                              quality_flag='GOOD', source_product_id='ARITHMETIC-FIXTURE-NOT-OBSERVATION',
                              source_record='0'))
        model.append(dict(identity, model_value=str(values[quantity]), logical_surface_element='7'))
        criteria[quantity] = {'maximum_absolute_error': 0}
    return sites, reference, model, criteria


def thermal_fixture():
    sites = {'arithmetic-fixture': {'logical_surface_element': '7'}}
    samples = []
    for quantity, depth in (('surface', '0'), ('depth', '.02')):
        for index, value in enumerate((10, 30)):
            samples.append({'query_id': quantity+str(index), 'site_id': 'arithmetic-fixture',
                            'logical_surface_element': '7', 'ephemeris_time_tdb_s': str(index),
                            'depth_m': depth, 'sample_weight': str(index+1),
                            'time_support_id': 'arithmetic-support', 'frame': polar.FRAME,
                            'model_physical_temperature_K': str(value),
                            'model_bolometric_brightness_temperature_K': str(value-1) if depth == '0' else ''})
    supports = []
    for name, quantity, depth, aggregation in (
            ('surface', 'annual_maximum_surface_temperature_K', 0, 'maximum'),
            ('depth', 'annual_average_temperature_at_2cm_K', .02, 'weighted-mean')):
        supports.append({'quantity': quantity, 'depth_m': depth, 'aggregation': aggregation,
                         'coverage_qualified': True, 'sampling_evidence': {'fixture': 'ARITHMETIC ONLY'},
                         'temperature_kind': 'physical_temperature_K', 'start_tdb_s': 0,
                         'stop_tdb_s': 1, 'endpoint_convention': 'closed',
                         'time_support_id': 'arithmetic-support',
                         'queries_by_site': {'arithmetic-fixture': [name+'0', name+'1']},
                         'sample_weights': {name+'0': 1, name+'1': 2}})
    return sites, samples, supports


def illumination_fixture():
    return [{'site_id': 'arithmetic-fixture', 'frame': polar.FRAME,
             'time_support_id': 'arithmetic-support', 'logical_surface_element': '7',
             'epoch_count': '4', 'any_disc_occurrence_fraction': '.25',
             'terrain_any_disc_occurrence_fraction': '.25',
             'mean_visible_area_fraction': '.1', 'geometric_psr': '0'}]


def retention_fixture():
    """Native-column ARITHMETIC fixture, not a physical residence model."""
    sites, samples, supports = thermal_fixture()
    for row in samples:
        row.update({key: '' for key in polar.SURFACE_KINETICS_FIELDS})
        if row['depth_m'] == '0':
            row.update(surface_retention_timescale_s='2', na_surface_residence_time_s='1',
                       na_surface_thermal_retention='0', ar_surface_residence_time_s='inf',
                       ar_surface_thermal_retention='1')
    return sites, samples, supports, illumination_fixture()


def cohort_fixture(*, shared=False):
    """Tiny identity/accounting fixture; values are NEVER observations."""
    sites = {'first-metadata-site': {'logical_surface_element': '7'}}
    if shared:
        sites['second-metadata-site'] = {'logical_surface_element': '7'}
    rows, hazards = [], []
    for cohort, site in enumerate(sites.values()):
        for species in ('Na', 'Ar'):
            for time in (0, 1):
                rows.append({'species': species, 'cohort': str(cohort), 'element': site['logical_surface_element'],
                             'elapsed_s': str(time), 'initial_atoms': '1',
                             'first_survival': '1' if time == 0 else '.5',
                             'restricted_first_residence_s': '0' if time == 0 else '.75',
                             'surface_fraction': '1' if time == 0 else '.5',
                             'gas_fraction': '0' if time == 0 else '.5',
                             'escape_fraction': '0', 'numerical_loss_fraction': '0'})
            hazards.append({'species': species, 'cohort': str(cohort), 'element': site['logical_surface_element'],
                            'begin_elapsed_s': '0', 'dt_s': '1', 'temperature_K': '30',
                            'total_release_hazard_s-1': '1', 'thermal_release_hazard_s-1': '1',
                            'photon_release_hazard_s-1': '0'})
    return sites, rows, hazards


def kinetic_runtime_fixture(root):
    """Actual files for runtime-IDENTITY testing, never a native/physics run."""
    from m2_protocol import cases
    protocol = polar._json(polar.HERE/'M2_I19_I20_PROTOCOL.json')
    case = next(row for row in cases('cluster') if row['case'] == 'cohort' and
                row['dt_divisor'] == row['weight_divisor'] == 4 and row['seed_index'] == 0)
    sites, rows, hazards = cohort_fixture()
    horizon = protocol['physical_duration_s']
    for row in rows:
        if row['elapsed_s'] == '1':
            row['elapsed_s'] = str(horizon)
    for row in hazards:
        row['dt_s'] = str(horizon)
    def csv_file(name, values):
        path = root/name
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)
        return path
    application = root/'application.in'
    from integration_common import SPICE_FILES,parse_receipt
    spice=root/'SPICE-METADATA-FIXTURE';physical_data={}
    for name in SPICE_FILES:
        path=spice/'Kernels'/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text('PROTOCOL METADATA ONLY: NOT A VALID SPICE KERNEL\n');physical_data[name]=polar.sha256(path)
    application.write_text('#begin application: moon\n'
        'spice_path = '+str(spice)+'\nsurface_geometry = sphere\n'
        'surface_temperature_mode = thermal-model\nthermal_initial_temperature_k = 30\n'
        'thermal_initialization_mode = declared-uniform\n'
        'm2_verification_record_initial_thermal_state = true\n'
        'simulation_start_utc = 2000-01-01T12:00:00Z\n'
        'cold_trap_retention_timescale_s = '+str(protocol['source_and_inventory_controls']['retention_reporting_timescale_s'])+'\n'
        'm2_verification_profile = cluster\nm2_verification_case = cohort\n'
        'm2_verification_dt_divisor = 4\nm2_verification_weight_divisor = 4\n'
        'm2_verification_surface_elements = 7\nm2_verification_cohort_packets = 32\n'
        'm2_verification_source_integral_tolerance = 1e-8\n#end application\n')
    receipt = root/'receipt.txt'
    initial_state=root/'receipt.txt.I19_initial_thermal_state.csv'
    initial_state.write_text('logical_surface_element,depth_index,depth_m,kinetic_temperature_K\n'+''.join(
        f'{element},{layer},{depth},30\n' for element in range(8) for layer,depth in ((-1,0),(0,.1),(1,.3))))
    receipt.write_text('test_id=I19\nmpi_ranks=8\nverification_case=cohort\n'
                       'seed='+str(case['seed'])+'\nduration_s='+str(horizon)+'\ninitial_epoch_tdb_s=123\n'
                       'initialized_na_weight_atoms=0.03125\ninitialized_ar_weight_atoms=0.03125\n'
                       'thermal_initialization_mode=declared-uniform\nthermal_initialization_et_s=123\n'
                       'initial_thermal_state_csv='+str(initial_state)+'\ninitial_thermal_state_rows=24\n'
                       'initial_thermal_state_element_count=8\ninitial_thermal_state_layer_count=2\n')
    with receipt.open('a') as stream:
        stream.write('thermal_initialization_spice_kernel_count='+str(len(SPICE_FILES))+'\n'+''.join(
            'thermal_initialization_spice_kernel_'+str(index)+'_file='+str(spice/'Kernels'/name)+'\n'
            for index,name in enumerate(SPICE_FILES)))
    manifest = root/'run_manifest.json'
    from m2_initial_state import physics_identity_sha256
    identity = {'test_id': 'I19', 'execution_scope': 'cluster-numerical-design', 'case': case,
                'started_utc':'2026-11-02T00:00:00Z','physics_identity_sha256':physics_identity_sha256(application),
                'mpi_ranks': 8, 'openmp_threads': 1, 'random_seed': case['seed'],
                'source_revision': 'PROTOCOL-FIXTURE-NOT-A-NATIVE-BUILD',
                'external_data':physical_data,
                'executable': {'sha256': 'a'*64},
                'inputs': {str(application): polar.sha256(application)},
                'build_configuration_hashes': {'FIXTURE-NOT-A-BUILD': 'b'*64},
                'protocol': {'sha256': polar.sha256(polar.HERE/'M2_I19_I20_PROTOCOL.json')}}
    manifest.write_text(json.dumps(identity))
    budgets = [{'species': species, 'elapsed_s': str(elapsed), 'species_elapsed_s': str(elapsed),
                'epoch_tdb_s': str(123+elapsed), 'actual_epoch_tdb_s': str(123+elapsed)}
               for species in ('Na', 'Ar') for elapsed in (0, horizon)]
    paths = dict(application=application, receipt=receipt, manifest=manifest,
                 cohorts=csv_file('cohorts.csv', rows), hazards=csv_file('hazards.csv', hazards),
                 budget=csv_file('budget.csv', budgets),initial_state=initial_state)
    descriptors = {key: {'path': str(path), 'sha256': polar.sha256(path)} for key, path in paths.items()}
    record = dict(seed_index=0, **descriptors)
    from m2_initial_state_evidence import read_state,physical_external_data
    state_qualification={'approved_utc':'2026-11-01T00:00:00Z','state_sha256':polar.sha256(initial_state),
        'shape':read_state(initial_state,8,2),'initialization_mode':'declared-uniform','ephemeris_time_tdb_s':123,
        'physics_identity_sha256':physics_identity_sha256(application),'executable_sha256':'a'*64,
        'qualification_sha256':'f'*64,'frozen_surface_elements':[7],
        'physical_external_data':physical_external_data(identity,application,parse_receipt(receipt)),
        'numerical_controls':{'m2_verification_cohort_packets':32,'m2_verification_source_integral_tolerance':1e-8}}
    qualified = {'native_evidence': [descriptors['receipt']], 'native_source_evidence':list(descriptors.values()),
        'predecessor_and_build_evidence':{'initial_state_qualification':state_qualification},
        'applied_initial_state':{'native_case_count':288,'qualification_sha256':'f'*64,
            'cases':{case['key']:{'path':str(initial_state),'sha256':polar.sha256(initial_state)}}}}
    runtime = {key: identity[key] for key in ('mpi_ranks', 'openmp_threads', 'executable', 'build_configuration_hashes')}
    return record, sites, protocol, runtime, application, qualified


def metadata_graph(root):
    """Create a tiny PROTOCOL graph with no readable observational science.

    The fake prerequisite labels exist only under TemporaryDirectory and the
    package authority is explicitly patched by each test. This graph tests
    hash/chronology/selection handling; it cannot qualify canonical D02/D03.
    """
    def write(name, value):
        path = root/name;path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) if isinstance(value, dict) else value)
        return path
    def descriptor(path):
        return {'path': str(path), 'sha256': polar.sha256(path)}
    review = {'approved': True, 'reviewer': 'INFRASTRUCTURE FIXTURE ONLY',
              'approved_utc': '2026-11-01T00:00:00Z', 'before_model_evaluation': True}
    collocation = write('collocation.json', {'status': 'APPROVED', 'selection_frozen': True,
                                            'review': review, 'maximum_distance_m': 3})
    candidates = root/'candidates.csv'
    row = candidate(0)
    with candidates.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row));writer.writeheader();writer.writerow(row)
    selection = polar.select_sites(polar._csv(candidates), maximum_distance_m=3)
    registry = write('sites.json', dict(selection, schema='srcMoon-I20-native-site-registry-v1',
                                       status='READY', selection_frozen=True, review=review,
                                       candidates=descriptor(candidates), collocation_contract=descriptor(collocation)))
    evidence = write('proof.txt', 'PROTOCOL ONLY; no physical or observational evidence')
    operator = write('operator.json', {'status': 'PASS', 'execution': 'executed',
                                      'site_registry_sha256': polar.sha256(registry),
                                      'collocation_contract_sha256': polar.sha256(collocation),
                                      'scorer_sha256': polar.sha256(Path(polar.__file__)),
                                      'evidence_files': [descriptor(evidence)]})
    # These values are deliberately not CSV. A metadata preflight must never
    # attempt to parse/score reference or predicted science values.
    not_science = write('not-science.csv', 'DO NOT OPEN AS OBSERVATIONS')
    queries = write('queries.csv', 'GEOMETRY PROTOCOL PLACEHOLDER')
    _, _, supports = thermal_fixture()
    for support in supports:
        support['sampling_evidence'] = descriptor(evidence)
    contract = write('polar.json', {'schema': polar.SCHEMA, 'status': 'APPROVED',
                                    'selection_frozen': True, 'review': review, 'frame': polar.FRAME,
                                    'coordinate_convention': 'planetocentric-east-positive',
                                    'site_registry': descriptor(registry), 'candidates': descriptor(candidates),
                                    'collocation_contract': descriptor(collocation),
                                    'reference_rows': descriptor(not_science), 'thermal_queries': descriptor(queries),
                                    'operator_qualification': descriptor(operator),
                                    'initial_state_qualification':descriptor(evidence),
                                    'application_template_sha256': 'a'*64, 'thermal_supports': supports,
                                    'criteria': {q: {'maximum_absolute_error': 0} for q in polar.QUANTITIES}})
    numerical = write('numeric.json', {'status': 'PASS', 'scope': 'fixed full numerical design; not observational agreement'})
    native = write('native.json', {'schema': 'srcMoon-I20-production-export-v1', 'status': 'PASS',
                                  'execution': 'executed', 'started_utc': '2026-11-02T00:00:00Z',
                                  'comparison_contract_sha256': polar.sha256(contract),
                                  'application_template_sha256': 'a'*64,
                                  'source_revision': 'PROTOCOL-NOT-A-REAL-BUILD',
                                  'input': descriptor(evidence), 'executable': descriptor(evidence),
                                  'runtime_manifest': descriptor(evidence), 'runtime_receipt': descriptor(evidence),
                                  'initial_state':descriptor(evidence),
                                  'model_rows': descriptor(not_science), 'thermal_samples': descriptor(not_science),
                                  'illumination_model': descriptor(not_science),
                                  'numerical_qualification': descriptor(numerical),
                                  'native_evidence': [descriptor(evidence)], 'kinetic_runs': ['protocol-only']})
    for test in polar.PREDECESSORS:
        write(test+'/result.json', {'status': 'PASS', 'scope': 'PROTOCOL-FIXTURE-ONLY'})
        write(test+'/run_manifest.json', {'mpi_ranks': 1, 'executable': descriptor(evidence)})
        write(test+'/receipt.txt', 'test_id='+test+'\nmpi_ranks=1\n')
    return contract, native


class SelectionProtocol(unittest.TestCase):
    def test_order_independent_approved_hash_selection(self):
        values = [candidate(index) for index in range(8)]
        first = polar.select_sites(values, maximum_distance_m=3)
        second = polar.select_sites(list(reversed(values)), maximum_distance_m=3)
        self.assertEqual(first, second)
        self.assertEqual(len(first['sites']), 3)
        self.assertEqual([row['selection_sha256'] for row in first['sites']],
                         sorted(row['selection_sha256'] for row in first['sites']))

    def test_shared_model_elements_preserved(self):
        result = polar.select_sites([candidate(i) for i in range(3)], maximum_distance_m=3)
        self.assertTrue(result['shared_model_support'])
        self.assertEqual(len(result['sites']), 3)
        self.assertEqual({row['logical_surface_element'] for row in result['sites']}, {'7'})

    def test_real_missing_strata_preserved(self):
        result = polar.select_sites([candidate(0)], maximum_distance_m=3)
        self.assertEqual(sum(row['selected'] for row in result['strata']), 1)
        self.assertGreater(sum(row['missing'] for row in result['strata']), 0)

    def test_edges_distance_and_quality_are_predeclared(self):
        values = [candidate(0), candidate(1, edge=1), candidate(2), candidate(3)]
        values[2]['collocation_distance_m'] = '4'
        values[3]['d03_quality_flag'] = 'BAD'
        result = polar.select_sites(values, maximum_distance_m=3)
        self.assertEqual([row['site_id'] for row in result['sites']], ['metadata-fixture-0'])

    def test_selection_refuses_residual_inputs(self):
        row = dict(candidate(0), model_value='arbitrary')
        with self.assertRaises(ValueError):
            polar.select_sites([row], maximum_distance_m=3)

    def test_selection_refuses_duplicate_native_pixel(self):
        with self.assertRaises(ValueError):
            polar.select_sites([candidate(0), candidate(0)], maximum_distance_m=3)


class ComparisonArithmetic(unittest.TestCase):
    def test_exact_rows_pass_with_missing_uncertainty_still_missing(self):
        sites, reference, model, criteria = arithmetic_rows()
        result = polar.compare_rows(reference, model, sites, criteria)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['evaluated_rows'], 4)
        self.assertTrue(all(row['uncertainty_value'] == '' for row in result['comparisons']))
        self.assertIn('arithmetic only', result['scope'])

    def test_evaluated_limit_miss_fails(self):
        sites, reference, model, criteria = arithmetic_rows()
        model[-1]['model_value'] = '81'
        result = polar.compare_rows(reference, model, sites, criteria)
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(len(result['failures']), 1)
        self.assertEqual(criteria[model[-1]['quantity']]['maximum_absolute_error'], 0)

    def test_missing_criterion_is_unavailable_not_auto_acceptance(self):
        sites, reference, model, criteria = arithmetic_rows()
        criteria.pop(model[-1]['quantity'])
        with self.assertRaises(polar.MissingEvidence):
            polar.compare_rows(reference, model, sites, criteria)

    def test_missing_or_extra_row_cannot_shrink_family(self):
        sites, reference, model, criteria = arithmetic_rows()
        with self.assertRaises(ValueError):
            polar.compare_rows(reference, model[:-1], sites, criteria)
        with self.assertRaises(ValueError):
            polar.compare_rows(reference, model+[copy.deepcopy(model[-1])], sites, criteria)

    def test_depth_quantity_time_frame_and_kind_are_distinct(self):
        for column, value in (('depth_m', '0'), ('time_support_id', 'different'),
                              ('frame', 'J2000'), ('temperature_kind', 'bolometric_brightness_temperature_K')):
            sites, reference, model, criteria = arithmetic_rows()
            model[-1][column] = value
            with self.subTest(column=column), self.assertRaises(ValueError):
                polar.compare_rows(reference, model, sites, criteria)

    def test_flagged_reference_and_wrong_model_support_rejected(self):
        sites, reference, model, criteria = arithmetic_rows()
        reference[0]['quality_flag'] = 'BAD'
        with self.assertRaises(ValueError):
            polar.compare_rows(reference, model, sites, criteria)
        reference[0]['quality_flag'] = 'GOOD';model[0]['logical_surface_element'] = '8'
        with self.assertRaises(ValueError):
            polar.compare_rows(reference, model, sites, criteria)


class ModelTemporalOperator(unittest.TestCase):
    def test_explicit_max_and_weighted_mean_of_model_samples(self):
        sites, samples, supports = thermal_fixture()
        result = polar.aggregate_thermal_samples(samples, supports, sites)
        self.assertEqual(result[0]['model_value'], 30)
        self.assertEqual(result[1]['model_value'], 70/3)
        self.assertEqual(result[1]['depth_m'], .02)

    def test_instantaneous_sample_does_not_imply_annual_coverage(self):
        sites, samples, supports = thermal_fixture()
        supports[0]['coverage_qualified'] = False
        with self.assertRaises(polar.MissingEvidence):
            polar.aggregate_thermal_samples(samples, supports, sites)

    def test_model_weights_and_complete_membership_are_fixed(self):
        sites, samples, supports = thermal_fixture()
        samples[0]['sample_weight'] = '2'
        with self.assertRaises(ValueError):
            polar.aggregate_thermal_samples(samples, supports, sites)
        sites, samples, supports = thermal_fixture()
        with self.assertRaises(ValueError):
            polar.aggregate_thermal_samples(samples[1:], supports, sites)

    def test_subsurface_temperature_not_brightness(self):
        sites, samples, supports = thermal_fixture()
        supports[1]['temperature_kind'] = 'bolometric_brightness_temperature_K'
        with self.assertRaises(ValueError):
            polar.aggregate_thermal_samples(samples, supports, sites)

    def test_explicit_endpoint_convention_is_enforced(self):
        sites, samples, supports = thermal_fixture()
        supports[0]['endpoint_convention'] = 'left-closed-right-open'
        with self.assertRaises(ValueError):
            polar.aggregate_thermal_samples(samples, supports, sites)

    def test_native_illumination_occurrence_not_visible_area(self):
        sites, samples, supports = thermal_fixture()
        illumination = illumination_fixture()
        rows = polar.assemble_model_rows(samples, illumination, supports, sites)
        values = {row['quantity']: row['model_value'] for row in rows}
        self.assertEqual(values['any_solar_disc_visibility_occurrence_fraction'], .25)
        illumination[0]['any_disc_occurrence_fraction'] = '.1'
        with self.assertRaises(ValueError):
            polar.assemble_model_rows(samples, illumination, supports, sites)

    def test_earth_occultation_does_not_create_geometric_psr(self):
        sites, samples, supports = thermal_fixture()
        illumination = illumination_fixture()
        illumination[0]['any_disc_occurrence_fraction'] = '0'
        rows = polar.assemble_model_rows(samples, illumination, supports, sites)
        values = {row['quantity']: row['model_value'] for row in rows}
        self.assertEqual(values['any_solar_disc_visibility_occurrence_fraction'], 0)
        self.assertEqual(values['permanent_shadow_classification'], 0)

    def test_terrain_mask_cannot_be_inferred_from_joint_visibility(self):
        sites, samples, supports = thermal_fixture()
        illumination = illumination_fixture();del illumination[0]['terrain_any_disc_occurrence_fraction']
        with self.assertRaises(polar.MissingEvidence):
            polar.assemble_model_rows(samples, illumination, supports, sites)


class NativeRetentionDiagnostics(unittest.TestCase):
    def check(self, fixture=None, **kwargs):
        sites, samples, supports, sky = fixture or retention_fixture()
        return polar.surface_retention_samples(samples, supports, sites, sky,
                                               retention_timescale_s=2, **kwargs)

    def test_native_zero_release_limit_and_dual_classification_reported(self):
        fixture = retention_fixture();fixture[3][0]['any_disc_occurrence_fraction'] = '0'
        result = self.check(fixture, native_export_schema=polar.KINETICS_EXPORT_SCHEMA)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(len(result['rows']), 4)
        argon = [row for row in result['rows'] if row['species'] == 'Ar']
        self.assertTrue(all(row['surface_residence_time_s'] == 'inf' and row['thermal_retention'] == 1 for row in argon))
        self.assertTrue(all(row['terrain_only_geometric_psr'] == 0 and
                            row['benchmark_any_disc_occurrence_fraction'] == 0 and
                            not row['permanent_retention_established'] for row in result['rows']))
        # JSON remains standard, with the infinite limit represented as text.
        json.dumps(result, allow_nan=False)

    def test_evaluated_retention_invariant_miss_is_fail(self):
        fixture = retention_fixture();fixture[1][0]['na_surface_thermal_retention'] = '1'
        result = self.check(fixture)
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(len(result['failures']), 1)

    def test_native_zero_residence_is_instantaneous_release_limit(self):
        fixture = retention_fixture();fixture[1][0]['na_surface_residence_time_s'] = '0'
        self.assertEqual(self.check(fixture)['status'], 'PASS')

    def test_exact_threshold_uses_native_greater_or_equal_definition(self):
        fixture = retention_fixture()
        fixture[1][0].update(na_surface_residence_time_s='2', na_surface_thermal_retention='1')
        self.assertEqual(self.check(fixture)['status'], 'PASS')

    def test_missing_old_capability_skips_but_broken_declared_schema_errors(self):
        fixture = retention_fixture();del fixture[1][0]['na_surface_residence_time_s']
        with self.assertRaises(polar.MissingEvidence):
            self.check(fixture)
        with self.assertRaises(ValueError):
            self.check(fixture, native_export_schema=polar.KINETICS_EXPORT_SCHEMA)

    def test_malformed_residence_and_binary_state_are_errors(self):
        for field, value in (('na_surface_residence_time_s', '-1'),
                             ('na_surface_residence_time_s', 'NaN'),
                             ('na_surface_residence_time_s', '-inf'),
                             ('na_surface_thermal_retention', '2')):
            fixture = retention_fixture();fixture[1][0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.check(fixture)

    def test_subsurface_state_cannot_claim_surface_kinetics(self):
        fixture = retention_fixture();fixture[1][2]['na_surface_residence_time_s'] = '1'
        with self.assertRaises(ValueError):
            self.check(fixture)

    def test_site_time_depth_identity_and_explicit_timescale_enforced(self):
        for field, value in (('site_id', 'different'), ('ephemeris_time_tdb_s', '2'),
                             ('surface_retention_timescale_s', '3'), ('depth_m', '.02')):
            fixture = retention_fixture();fixture[1][0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.check(fixture)

    def test_missing_surface_query_cannot_shrink_report(self):
        fixture = retention_fixture();fixture[1].pop(0)
        with self.assertRaises(ValueError):
            self.check(fixture)

    def test_scorer_preserves_fail_error_skip_before_reference_access(self):
        # This patches ONLY metadata authority to exercise the score control
        # flow. Every model row is labelled arithmetic fixture, the would-be
        # reference is deliberately unreadable text, and no observations or
        # native computation are represented by this temporary test graph.
        for variant, expected in (('label-invariant', 'FAIL'), ('invalid-value', 'ERROR'),
                                  ('older-exporter', 'SKIPPED')):
            with self.subTest(variant=variant), tempfile.TemporaryDirectory(prefix='i20-score-control-fixture-') as temp:
                root = Path(temp);sites, samples, supports, sky = retention_fixture()
                protocol = polar._json(polar.HERE/'M2_I19_I20_PROTOCOL.json')
                timescale = protocol['source_and_inventory_controls']['retention_reporting_timescale_s']
                for row in samples:
                    if row['depth_m'] == '0':
                        row['surface_retention_timescale_s'] = str(timescale)
                if variant == 'label-invariant':
                    samples[0]['na_surface_thermal_retention'] = '1'
                elif variant == 'invalid-value':
                    samples[0]['na_surface_residence_time_s'] = 'NaN'
                def csv_file(name, values):
                    path = root/name
                    with path.open('w', newline='') as stream:
                        writer = csv.DictWriter(stream, fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)
                    return str(path)
                application = root/'application.in'
                application.write_text('cold_trap_retention_timescale_s = '+str(timescale)+'\n')
                runtime = {'test_id': 'I20', 'mpi_ranks': 8, 'source_revision': 'ARITHMETIC-FIXTURE',
                           'started_utc': '2026-11-02T00:00:00Z', 'executable': {'sha256': 'a'*64},
                           'inputs': {str(application): polar.sha256(application)}}
                manifest = root/'runtime.json';manifest.write_text(json.dumps(runtime))
                receipt = root/'receipt.txt'
                receipt.write_text('test_id=I20\nmpi_ranks=8\nthermal_model_rows=4\n'+
                    ('' if variant == 'older-exporter' else 'surface_kinetics_export_schema='+polar.KINETICS_EXPORT_SCHEMA+'\n'))
                reference = root/'reference-not-readable.csv';reference.write_text('FORBIDDEN: NOT OBSERVATIONS')
                gate = {'id': 'I20', 'status': 'PASS', 'scope': 'preflight-only', 'sites': sites,
                        'evidence':{'initial_state_qualification':{'approved_utc':'2026-11-01T00:00:00Z'}},
                        'contract': {'thermal_supports': supports},
                        'native_manifest': {'executable': runtime['executable'],
                                            'source_revision': runtime['source_revision'],
                                            'started_utc': runtime['started_utc']},
                        'files': {'model_rows': csv_file('model.csv', polar.assemble_model_rows(samples, sky, supports, sites)),
                                  'thermal_samples': csv_file('samples.csv', samples),
                                  'illumination_model': csv_file('sky.csv', sky),
                                  'runtime_manifest': str(manifest), 'runtime_receipt': str(receipt),
                                  'input': str(application), 'reference_rows': str(reference),'initial_state':str(reference)}}
                with patch.object(polar, 'preflight', return_value=gate), patch.object(polar, 'compare_rows') as compare, \
                     patch('m2_initial_state_evidence.verify_surface_initial_state',return_value={'scope':'CONTROL-FIXTURE-NOT-PHYSICS'}):
                    result = polar.score(None, root, root/'out', render=False)
                self.assertEqual(result['status'], expected, result)
                self.assertFalse(result['observational_scoring'])
                compare.assert_not_called()
                self.assertFalse((root/'out/comparison.csv').exists())
                self.assertFalse((root/'out/comparison.png').exists())
                self.assertEqual((root/'out/surface_retention.csv').exists(), expected == 'FAIL')


class CohortSiteBinding(unittest.TestCase):
    def check(self, fixture=None, **kwargs):
        sites, rows, hazards = fixture or cohort_fixture()
        controls = dict(declared_elements=[r['logical_surface_element'] for r in sites.values()],
                        horizon_s=1, initial_epoch_s=10, seed_index=0)
        controls.update(kwargs)
        return polar.bind_cohort_sites(rows, hazards, sites, **controls)

    def test_shared_surface_elements_keep_distinct_native_site_labels(self):
        diagnostics = self.check(cohort_fixture(shared=True))
        self.assertEqual(len(diagnostics), 8)
        self.assertEqual({row['site_id'] for row in diagnostics}, {'first-metadata-site', 'second-metadata-site'})
        self.assertEqual({row['logical_surface_element'] for row in diagnostics}, {7})
        self.assertEqual({row['ephemeris_time_tdb_s'] for row in diagnostics}, {10, 11})
        self.assertTrue(all(row['origin_depth_m'] == 0 and row['frame'] == polar.FRAME for row in diagnostics))

    def test_declared_input_order_cannot_replace_frozen_native_sites(self):
        with self.assertRaises(ValueError):
            self.check(declared_elements=[8])

    def test_native_species_site_and_hazard_identity_mismatches_are_errors(self):
        for which, field, value in ((1, 'species', 'He'), (1, 'element', '8'), (1, 'cohort', '1'),
                                     (2, 'species', 'He'), (2, 'element', '8'), (2, 'cohort', '1')):
            fixture = cohort_fixture();fixture[which][0][field] = value
            with self.subTest(which=which, field=field), self.assertRaises(ValueError):
                self.check(fixture)

    def test_complete_species_and_initial_final_times_required(self):
        for remove in (0, 1, 2):
            fixture = cohort_fixture();fixture[1].pop(remove)
            with self.subTest(remove=remove), self.assertRaises(ValueError):
                self.check(fixture)

    def test_duplicate_site_checkpoints_and_hazard_intervals_are_errors(self):
        for which in (1, 2):
            fixture = cohort_fixture();fixture[which].append(dict(fixture[which][0]))
            with self.subTest(which=which), self.assertRaises(ValueError):
                self.check(fixture)

    def test_hazard_interval_cannot_extend_past_actual_run(self):
        fixture = cohort_fixture();fixture[2][0]['dt_s'] = '2'
        with self.assertRaises(ValueError):
            self.check(fixture)

    def test_initial_stock_identity_and_pre_sweep_weight_normalization(self):
        fixture = cohort_fixture();fixture[1][1]['initial_atoms'] = '2'
        with self.assertRaises(ValueError):
            self.check(fixture)
        with self.assertRaises(ValueError):
            self.check(expected_initial_atoms={'Na': 2, 'Ar': 1})
        self.assertEqual(len(self.check(expected_initial_atoms={'Na': 1, 'Ar': 1})), 4)


class NativeKineticRuntimeBinding(unittest.TestCase):
    def validate(self, root, fixture):
        record, sites, protocol, runtime, application, qualified = fixture
        return polar.validate_kinetic_run(record, root, sites, protocol, runtime, application, qualified)

    def test_qualified_site_seed_clock_and_build_files_are_retained(self):
        with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
            root = Path(temp);fixture = kinetic_runtime_fixture(root)
            paths, diagnostics, identity = self.validate(root, fixture)
            self.assertEqual(len(diagnostics), 4)
            self.assertEqual(identity['initial_epoch_tdb_s'], 123)
            self.assertEqual(identity['mpi_ranks'], 8)
            self.assertEqual(set(paths), {'cohorts', 'hazards', 'application', 'receipt', 'manifest', 'budget', 'initial_state'})

    def test_unqualified_copy_of_native_bytes_cannot_replace_i19_evidence(self):
        with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
            root = Path(temp);fixture = kinetic_runtime_fixture(root)
            copied = root/'unqualified-copy.csv';copied.write_bytes(Path(fixture[0]['cohorts']['path']).read_bytes())
            fixture[0]['cohorts'] = {'path': str(copied), 'sha256': polar.sha256(copied)}
            with self.assertRaises(ValueError):
                self.validate(root, fixture)

    def test_missing_native_manifest_is_missing_capability(self):
        with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
            root = Path(temp);fixture = kinetic_runtime_fixture(root);fixture[0].pop('manifest')
            with self.assertRaises(polar.MissingEvidence):
                self.validate(root, fixture)

    def test_changed_manifest_bytes_are_error_before_native_values(self):
        with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
            root = Path(temp);fixture = kinetic_runtime_fixture(root)
            Path(fixture[0]['manifest']['path']).write_text('{}')
            with self.assertRaises(ValueError):
                self.validate(root, fixture)

    def test_wrong_layout_or_build_cannot_borrow_qualified_kinetics(self):
        for key, value in (('mpi_ranks', 4), ('openmp_threads', 2),
                           ('executable', {'sha256': 'c'*64}),
                           ('build_configuration_hashes', {'OTHER-BUILD': 'd'*64})):
            with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
                root = Path(temp);fixture = kinetic_runtime_fixture(root);fixture[3][key] = value
                with self.subTest(key=key), self.assertRaises(ValueError):
                    self.validate(root, fixture)

    def test_missing_qualified_native_evidence_is_skipped(self):
        with tempfile.TemporaryDirectory(prefix='i20-runtime-identity-fixture-') as temp:
            root = Path(temp);fixture = kinetic_runtime_fixture(root);fixture[5].clear()
            with self.assertRaises(polar.MissingEvidence):
                self.validate(root, fixture)


class FailureSemantics(unittest.TestCase):
    def test_missing_contract_skips_without_science_outputs(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            root = Path(temp)
            result = polar.score(None, root, root/'out', render=False)
            self.assertEqual(result['status'], 'SKIPPED')
            self.assertFalse(result['observational_scoring'])
            self.assertFalse((root/'out/comparison.csv').exists())
            self.assertFalse((root/'out/comparison.png').exists())
            self.assertFalse((root/'out/comparison.eps').exists())

    def test_json_duplicate_and_nonfinite_evidence_is_error(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            path = Path(temp)/'contract.json'
            for text in ('{"schema":"x","schema":"y"}', '{"value":NaN}'):
                path.write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    polar._json(path)

    def test_hash_mismatch_is_not_missing_evidence(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            path = Path(temp)/'metadata.txt';path.write_text('protocol fixture only')
            with self.assertRaises(ValueError):
                polar._descriptor({'path': str(path), 'sha256': '0'*64}, Path(temp), 'fixture')

    def test_unreviewed_registry_cannot_be_implicitly_frozen(self):
        with self.assertRaises(polar.MissingEvidence):
            polar._review({'status': 'READY', 'selection_frozen': True}, 'fixture')

    def test_existing_scientific_result_is_never_overwritten(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            path = Path(temp)/'result.json'
            polar._write_new_json(path, {'fixture': 1})
            with self.assertRaises(ValueError):
                polar._write_new_json(path, {'fixture': 2})

    def test_duplicate_optional_input_is_error_not_missing(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            root = Path(temp);template = root/'receipt.in'
            template.write_text('m2_formal_contract_file = first.json\n'
                                'm2_formal_contract_file = second.json\n')
            with self.assertRaises(ValueError):
                polar.prepare_input(template, root/'application.in', root)
            self.assertFalse((root/'application.in').exists())

    def test_malformed_contract_is_error_without_science_outputs(self):
        with tempfile.TemporaryDirectory(prefix='i20-parser-fixture-') as temp:
            root = Path(temp);contract = root/'contract.json'
            contract.write_text('{invalid-json')
            result = polar.score(contract, root, root/'out', render=False)
            self.assertEqual(result['status'], 'ERROR')
            self.assertFalse((root/'out/comparison.csv').exists())

    def test_formal_staging_preserves_mesh_and_source_physics(self):
        # Metadata/launcher fixture only. The patched gate stops all real
        # data processing; the assertion protects the meaningful production
        # boundary between bounded I18 fixtures and formal I20 input staging.
        with tempfile.TemporaryDirectory(prefix='i20-input-fixture-') as temp:
            root = Path(temp);contract = root/'contract.json';contract.write_text('{}')
            query = root/'queries.csv';query.write_text('GEOMETRY FIXTURE ONLY\n')
            template = root/'source.in'
            template.write_text('#begin application: moon\n'
                                'm2_formal_contract_file = contract.json\n'
                                'surface_temperature_mode = thermal-model\n'
                                'mesh_profile = full\n'
                                'mesh_near_surface_cell_size_m = 12345\n'
                                'argon_source_enabled = false\n'
                                'surface_cea_file = original.cea\n'
                                'surface_tecplot_file = original.dat\n'
                                '#end application\n')
            from m2_initial_state import physics_identity_sha256
            gate = {'status': 'PASS', 'scope': 'preflight-only', 'contract_path': str(contract),
                    'files': {'thermal_queries': str(query)},
                    'contract': {'application_template_sha256': polar.sha256(template)},
                    'evidence': {'initial_state_qualification':{'physics_identity_sha256':physics_identity_sha256(template)}}}
            destination = root/'output/application.in'
            with patch.object(polar, 'preflight', return_value=gate):
                polar.prepare_input(template, destination, root)
            from m2_illumination import setting
            for key, value in (('mesh_profile', 'full'), ('mesh_near_surface_cell_size_m', '12345'),
                               ('argon_source_enabled', 'false')):
                self.assertEqual(setting(destination, key), value)


class MetadataOrchestration(unittest.TestCase):
    def gate(self, root, contract, native):
        # This patched authority is confined to testing orchestration. It
        # never changes canonical metadata or produces a formal result file.
        with patch('integration_common.dataset_package_preflight', return_value=(True, 'PROTOCOL ONLY', {})):
            with patch.object(polar,'_qualified_initial_state',return_value={'scope':'METADATA AUTHORITY FIXTURE ONLY'}):
                return polar.preflight(contract, root, model_manifest_path=native)

    def test_preflight_does_not_read_reference_or_model_values(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            result = self.gate(root, contract, native)
            self.assertEqual(result['status'], 'PASS', result)
            self.assertFalse(result['observational_scoring'])
            self.assertEqual(result['scope'], 'preflight-only')

    def test_posthoc_freeze_is_error_before_residual_access(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            value = json.loads(native.read_text());value['started_utc'] = '2026-10-31T00:00:00Z'
            native.write_text(json.dumps(value))
            self.assertEqual(self.gate(root, contract, native)['status'], 'ERROR')

    def test_supporting_predecessor_cannot_be_promoted(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            (root/'I15/result.json').write_text('{"status":"PASS","scope":"supporting-development"}')
            result = self.gate(root, contract, native)
            self.assertEqual(result['status'], 'SKIPPED')
            self.assertIn('formal linked', result['reason'])

    def test_wrong_template_binding_is_error(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            value = json.loads(native.read_text());value['application_template_sha256'] = 'b'*64
            native.write_text(json.dumps(value))
            self.assertEqual(self.gate(root, contract, native)['status'], 'ERROR')

    def test_missing_native_export_skips_after_reference_preflight(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            native.unlink()
            self.assertEqual(self.gate(root, contract, native)['status'], 'SKIPPED')

    def test_metadata_only_launch_route_needs_no_future_output_hash(self):
        with tempfile.TemporaryDirectory(prefix='i20-metadata-fixture-') as temp:
            root = Path(temp);contract, native = metadata_graph(root)
            native.unlink()
            with patch('integration_common.dataset_package_preflight', return_value=(True, 'PROTOCOL ONLY', {})):
                with patch.object(polar,'_qualified_initial_state',return_value={'scope':'METADATA AUTHORITY FIXTURE ONLY'}):
                    result = polar.preflight(contract, root, require_model=False)
            self.assertEqual(result['status'], 'PASS', result)
            self.assertNotIn('native_manifest', result)
            self.assertFalse(result['observational_scoring'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
