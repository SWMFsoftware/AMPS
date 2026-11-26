#!/usr/bin/env python3
"""Controls for the I22 atomic photoionization pilot (approved 2026-11-26).

Covers the verifier-side atomic target table, the source-bound input preparer
and crosswalk binding, the frozen case generator, the domain and speed limits
that set the approved cases, and the I22 plan/owner wiring. Real D04 and
profile files are used only when present (otherwise the test is skipped);
no native run is involved and nothing here is physical evidence.
"""
import hashlib
import importlib.util
import json
import math
from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

TEST = Path(__file__).resolve().parent
ROOT = TEST.parents[1]
I22 = TEST / 'integration-tests/I22'
sys.path.insert(0, str(TEST))
sys.path.insert(0, str(TEST / 'sodium_cohort'))
import atomic_targets  # noqa: E402
import m3_qualification as q  # noqa: E402
import verify_numerical as numerical  # noqa: E402
import verify_ledger as ledger  # noqa: E402
import prepare_physical_case as preparation  # noqa: E402

D04 = Path('/data/vtenishe/moon_validation_data/photochemistry')
PROFILE = ROOT / 'test_output/srcMoon/m3-followup-MfThQfny/atomic-family/profile.json'
APPROVAL = I22 / 'reference/atomic_pilot_approval_20261126.json'
ELECTRON_MASS_KG = Decimal('9.1093837139E-31')   # CODATA 2022, as in the profile constants
ELECTRON_VOLT_J = Decimal('1.602176634E-19')     # exact SI
AMU_KG = Decimal('1.66053906892E-27')            # CODATA 2022 atomic mass constant
ISOTOPE_U = {'He': Decimal('4.00260325413'), 'Ne': Decimal('19.9924401753'), 'Ar': Decimal('39.9623831237')}
# Moon build domain (srcMoon/main_lib.cpp: xMaxDomain=400, yMaxDomain=30,
# BlockCells=30,5,5, R=1737.4 km): x in [-436.67 R, 30 R], |y|,|z| <= 38.89 R.
R_MOON = Decimal('1737400')
DOMAIN_X_MAX = 30 * R_MOON
DOMAIN_YZ = Decimal(400) * (1 + Decimal(5) / 30) * Decimal(5) / 30 / 2 * R_MOON


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def d04_rows():
    import csv
    path = D04 / 'processed/phidrates/photo_rates_1au.csv'
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


class AtomicTargetTable(unittest.TestCase):
    def test_table_and_scopes(self):
        self.assertEqual(atomic_targets.nucleus('He'), '4He')
        self.assertEqual(atomic_targets.products('Ar'), 'Ar+ + e')
        self.assertEqual(atomic_targets.element_of_nucleus('20Ne'), 'Ne')
        for bad in ('K', 'H2O', ''):
            with self.assertRaises(ValueError):
                atomic_targets.require(bad)
        with self.assertRaises(ValueError):
            atomic_targets.element_of_nucleus('39K')
        atomic_targets.check_test('I21', 'Na')
        atomic_targets.check_test('I22', 'He')
        for test_id, element in (('I21', 'He'), ('I22', 'Na'), ('I23', 'He')):
            with self.assertRaises(ValueError):
                atomic_targets.check_test(test_id, element)

    def test_receipt_element(self):
        # A receipt without source_parent predates the pilot and is Na by construction.
        self.assertEqual(atomic_targets.receipt_element({'source_channels': [{'id': 'PHIDRATES:Na:00:na_e'}]}), 'Na')
        self.assertEqual(atomic_targets.receipt_element(
            {'source_parent': 'He', 'source_channels': [{'id': 'PHIDRATES:He:00:he_e'}]}), 'He')
        for receipt in ({'source_parent': 'He', 'source_channels': [{'id': 'PHIDRATES:Ne:00:ne_e'}]},
                        {'source_parent': 'Xe'},
                        {'source_parent': 'Ar', 'source_channels': []}):
            with self.assertRaises(ValueError):
                atomic_targets.receipt_element(receipt)


@unittest.skipUnless((D04 / 'qa_report.json').is_file(), 'canonical D04 package unavailable')
class D04Bindings(unittest.TestCase):
    def reference(self):
        return {'solar_endpoint': 'quiet', 'total_acceleration_macro': 'TotalParticleAcceleration_default(a,s,p,x,v,n);',
                'point_sun_visible': 1, 'sun_position_lso_m': ['1', '0', '0'], 'source_position_lso_m': ['0', '0', '0'],
                'source_reference_distance_m': '1', 'source_reference_distance_identity': 'ANALYTICAL-UNIT-DISTANCE',
                'geometry_basis_epoch_et_s': '0'}

    def test_numerical_source_rate_selects_the_parent_family(self):
        table = {r['species']: r for r in d04_rows() if r['solar_activity'] == '0.0'}
        for element in ('Na', 'He', 'Ne', 'Ar'):
            rate, channel = numerical.source_rate(self.reference(), {}, element)
            self.assertEqual(channel, table[element]['reaction_id'])
            self.assertEqual(rate, Decimal.from_float(float(table[element]['rate_1au_s-1'])))
        with self.assertRaises(ValueError):
            numerical.source_rate(self.reference(), {}, 'H2O')

    def test_ledger_photon_source_carries_its_element(self):
        base = D04 / 'processed/phidrates'
        pins = {k: {'path': str(base / n), 'sha256': hashlib.sha256((base / n).read_bytes()).hexdigest()}
                for k, n in (('rates', 'photo_rates_1au.csv'), ('thresholds', 'thresholds.csv'))}
        from fractions import Fraction as F
        constants = {'elementary charge': (F('1.602176634e-19'), 'C')}
        for element, channel in (('He', 'PHIDRATES:He:00:he_e'), ('Ar', 'PHIDRATES:Ar:00:ar_e')):
            value = ledger.load_photon_source({'photon_source': {'reaction_id': channel, 'solar_activity': '1', **pins}},
                                              base, {}, constants, 'source-bound-experimental')
            self.assertEqual(value['element'], element)
        with self.assertRaises(ValueError):
            ledger.load_photon_source({'photon_source': {'reaction_id': 'PHIDRATES:H2O:00:h2o_e', 'solar_activity': '0', **pins}},
                                      base, {}, constants, 'source-bound-experimental')


class CaseDesign(unittest.TestCase):
    def test_generator_reproduces_frozen_cases(self):
        generator = load('i22_cases', I22 / 'prepare_atomic_pilot_cases.py')
        files = generator.render()
        frozen = I22 / 'reference/atomic_pilot'
        self.assertEqual(sorted(files), sorted(p.name for p in frozen.iterdir()))
        for name, text in files.items():
            self.assertEqual((frozen / name).read_text(), text, name)
        self.assertEqual(generator.sha(generator.PROTOCOL), json.loads(files['He-quiet.json'])['controls']['transport_policy_evidence_sha256'])
        # Approved counts with the later He-quiet revision (16,384 -> 8,192) applied.
        counts = generator.parent_counts()
        self.assertEqual(counts, {'He-quiet': 8192, 'He-active': 8192, 'Ne-quiet': 2048, 'Ne-active': 512,
                                  'Ar-quiet': 256, 'Ar-active': 128})
        for name, count in counts.items():
            self.assertEqual(json.loads(files[name + '.json'])['controls']['original_parents'], str(count))
        # Revision b: production builds the photon history over each mover
        # substep, so the substep cap must equal the coarsest common step
        # (horizon/8) or step refinement leaves the exposure grid unchanged
        # (the campaign-2 Ne/Ar plateau). He keeps the I21 value 8192 = 65536/8.
        caps = {}
        for name in counts:
            element, endpoint = name.split('-')
            step = json.loads(files[name + '.json'])['controls']['physical_step_s']
            for order in ('first', 'second'):
                lines = [l for l in files[f'{element}_{order}_{endpoint}.in'].splitlines()
                         if l.split('=', 1)[0].strip() == 'tracing_maximum_substep_s']
                self.assertEqual(lines, ['tracing_maximum_substep_s = ' + step], name)
            caps[name] = step
        self.assertEqual(caps, {'He-quiet': '8192', 'He-active': '8192', 'Ne-quiet': '32768', 'Ne-active': '32768',
                                'Ar-quiet': '131072', 'Ar-active': '131072'})

    @unittest.skipUnless((D04 / 'qa_report.json').is_file(), 'canonical D04 package unavailable')
    def test_ions_stay_in_the_domain_and_electrons_below_the_cap(self):
        """Independent worst case for the approved horizons and the 3.4e6 m/s cap."""
        approval = json.loads(APPROVAL.read_text())['cases']
        rows = {(r['species'], r['solar_activity']): r for r in d04_rows()}
        for element in ('He', 'Ne', 'Ar'):
            for endpoint, activity in (('quiet', '0.0'), ('active', '1.0')):
                name = f'{element}-{endpoint}'
                # Energy above threshold <= mean excess + |physical - printed| < mean + 1e-3 eV.
                energy = (Decimal(rows[(element, activity)]['excess_energy_eV']) + Decimal('0.001')) * ELECTRON_VOLT_J
                momentum = (2 * ELECTRON_MASS_KG * energy).sqrt()
                ion_speed = momentum / (ISOTOPE_U[element] * AMU_KG)
                electron_speed = momentum / ELECTRON_MASS_KG
                horizon = Decimal(approval['horizons_s'][name])
                reach_x = Decimal(10**7) + (ion_speed + 1) * horizon
                self.assertLess(reach_x, DOMAIN_X_MAX - Decimal('1.1e7'), name)  # approval: at least 1.1e7 m
                self.assertLess(ion_speed * horizon, DOMAIN_YZ, name)
                # Moon heliocentric speed (~3e4 m/s) included as a conservative margin.
                self.assertLess(electron_speed + Decimal('3.1e4'), Decimal('3.4e6'), name)
                relativistic = Decimal('0.75') * (Decimal('3.4e6') / Decimal(299792458)) ** 2
                self.assertLess(relativistic, Decimal('1e-4'))


@unittest.skipUnless(PROFILE.is_file() and (D04 / 'qa_report.json').is_file(), 'atomic-family profile or D04 absent')
class InputPreparation(unittest.TestCase):
    def setUp(self):
        self.tool = I22 / 'prepare_atomic_pilot_inputs.py'
        self.digest = hashlib.sha256(PROFILE.read_bytes()).hexdigest()

    def run_tool(self, output, digest=None):
        return subprocess.run([sys.executable, '-B', str(self.tool), '--profile', str(PROFILE),
                               '--profile-sha256', digest or self.digest, '--output-dir', str(output)],
                              capture_output=True, text=True, timeout=300)

    def test_preparation_and_crosswalk_binding(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'inputs'
            done = self.run_tool(output)
            self.assertEqual(done.returncode, 0, done.stderr + done.stdout)
            profile = json.loads((output / 'profile.json').read_text())
            self.assertEqual(sorted(profile['atomic_photoionization_crosswalks']), ['Ar', 'He', 'Ne'])
            self.assertEqual(profile['derived_from']['sha256'], self.digest)
            self.assertIs(profile['physical_activation'], False)
            for element, isotope in (('He', '4He'), ('Ne', '20Ne'), ('Ar', '40Ar')):
                row, _ = preparation.authenticate_crosswalk(profile, output / 'profile.json',
                                                            output / f'{element}_source_energy_crosswalk.csv')
                self.assertEqual((row['parent_state_id'], row['ion_state_id']), (isotope + ':ground', isotope + '+:ground'))
                self.assertTrue(row['physical_threshold_source_identity'].startswith('GROUND-ION-' + isotope + '@sha256:'))
            # The Na atomic_extension still binds; an unlisted path is refused.
            na = Path(profile['atomic_extension']['crosswalk_csv']['path'])
            self.assertEqual(preparation.authenticate_crosswalk(profile, output / 'profile.json', na)[0]['ion_state_id'],
                             '23Na+:ground')
            stray = Path(temporary) / 'stray.csv'
            stray.write_bytes((output / 'He_source_energy_crosswalk.csv').read_bytes())
            with self.assertRaises(ValueError):
                preparation.authenticate_crosswalk(profile, output / 'profile.json', stray)
            # A listed crosswalk under the wrong parent is refused.
            swapped = json.loads(json.dumps(profile))
            listed = swapped['atomic_photoionization_crosswalks']
            listed['He'], listed['Ne'] = listed['Ne'], listed['He']
            with self.assertRaises(ValueError):
                preparation.authenticate_crosswalk(swapped, output / 'profile.json', output / 'He_source_energy_crosswalk.csv')

    def test_refusals_preserve_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            wrong = self.run_tool(Path(temporary) / 'wrong', '0' * 64)
            self.assertEqual(wrong.returncode, 2)
            self.assertEqual(json.loads((Path(temporary) / 'wrong/result.json').read_text())['status'], 'ERROR')
            existing = Path(temporary) / 'existing'
            existing.mkdir()
            self.assertEqual(self.run_tool(existing).returncode, 2)
            self.assertEqual(list(existing.iterdir()), [])


class PlanWiring(unittest.TestCase):
    def plan(self, test_id):
        return {'schema': 'srcMoon-M3-cohort-plan-v1', 'approved_addendum_sha256': q.APPROVED_SHA256,
                'test_id': test_id, 'cohort_id': 'atomic-I22-pilot-fixture', 'original_parents': 128,
                'initial_represented_number': 2.625, 'horizon_s': 1048576.0, 'reference_total_rate_s': 3.5e-7,
                'source_channel_rates_s': [{'id': 'PHIDRATES:Ar:00:ar_e', 'rate_s': 3.5e-7}], 'branch_claims': 6,
                'seed': 731993, 'primary_estimator': 'integrated-exposure', 'test_kind': 'trajectory-integrated-hazard'}

    def test_trajectory_plans_belong_to_I21_and_I22(self):
        message = 'belong to I21 and the I22 atomic pilot'
        for test_id in ('I21', 'I22'):
            try:
                q.validate_plan(self.plan(test_id))
            except ValueError as error:
                self.assertNotIn(message, str(error))
        for test_id in ('I23', 'I24'):
            with self.assertRaisesRegex(ValueError, message):
                q.validate_plan(self.plan(test_id))

    def test_producer_and_owner_options(self):
        import m3_qualification_executor as executor
        self.assertEqual(executor.PHYSICAL_PRODUCERS['I22'], 'normal-source-bound-atomic-photoionization-cohort-v1')
        self.assertEqual(executor.PHYSICAL_PRODUCERS['I21'], 'normal-source-bound-sodium-cohort-v1')
        with tempfile.TemporaryDirectory() as temporary:
            for test_id, option in (('I21', '--i22-cohort-plan'), ('I22', '--i21-cohort-plan')):
                done = subprocess.run([sys.executable, '-B', str(TEST / 'integration_m3.py'), test_id, '--amps',
                                       str(Path(temporary) / 'amps'), '--output-dir', str(Path(temporary) / test_id),
                                       option, str(Path(temporary) / 'plan.json'),
                                       option + '-sha256', '0' * 64], capture_output=True, text=True, timeout=120)
                self.assertEqual(done.returncode, 2, done.stderr)
                self.assertIn('must match the owner test ID', done.stderr)

    def test_build_mappings(self):
        build = load('m3_build_i22', TEST / 'build_m3_source_qualification.py')
        self.assertEqual(list(build.NEON_PHOTOIONIZATION_MAPPING), ['Na', 'Ne', 'Ne_plus'])
        self.assertEqual(list(build.ARGON_PHOTOIONIZATION_MAPPING), ['Na', 'Ar', 'AR_PLUS'])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'mapping.json'
            for good in (build.HELIUM_IMPACT_MAPPING, build.NEON_PHOTOIONIZATION_MAPPING, build.ARGON_PHOTOIONIZATION_MAPPING):
                path.write_text(json.dumps(good))
                self.assertEqual(build.validate_mapping(path), good)
            path.write_text(json.dumps({'Na': '23Na:ground', 'Ar': '40Ar:ground', 'Ne_plus': '20Ne+:ground'}))
            with self.assertRaises(ValueError):
                build.validate_mapping(path)


if __name__ == '__main__':
    unittest.main()
