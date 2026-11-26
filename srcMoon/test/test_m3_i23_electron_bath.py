#!/usr/bin/env python3
"""I23 electron-bath admission and protocol controls (infrastructure only).

These checks bind the 2026-11-25 approval record, the committed admission
generator and the frozen evidence files that the production composer
authenticates. They exercise refusal semantics and reproducibility; they are
not lunar science tests and never create observations or reference values.
Source-dependent checks SKIP when the authenticated profile or D05 package is
absent on this host.
"""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

TEST = Path(__file__).resolve().parent
I23 = TEST / 'integration-tests/I23'
sys.path.insert(0, str(I23))
import prepare_helium_admission as tool  # noqa: E402

CASE = json.loads((I23 / 'reference/helium_case.json').read_text())
APPROVAL = json.loads((I23 / 'reference/helium_admission_approval_20261125.json').read_text())
DATA_ROOT = Path('/data/vtenishe/moon_validation_data')


def sources_available():
    profile = Path(CASE['profile_path'])
    return profile.is_file() and tool.sha(profile) == CASE['profile_sha256'] and \
        (DATA_ROOT / 'photochemistry/raw/electron_impact/voronov/cfit.dat').is_file()


class ApprovalRecord(unittest.TestCase):
    def test_approval_scope_is_explicit_and_non_qualifying(self):
        self.assertEqual(APPROVAL['approval_status'], 'APPROVED')
        self.assertEqual(APPROVAL['approval_date'], '2026-11-25')
        self.assertEqual(APPROVAL['bath']['electron_densities_m3'], ['0', '1e6', '1e7'])
        self.assertEqual(APPROVAL['bath']['thermal_energy_eV'], '10')
        self.assertEqual(APPROVAL['bath']['frame'], 'LSO')
        self.assertEqual(APPROVAL['products']['formal_I23_status'], 'SKIPPED / NOT QUALIFIED')
        self.assertIn('no He+ velocity', APPROVAL['admission']['observable_scope'])
        self.assertIn('No observational agreement', APPROVAL['status_limits'])

    def test_case_pins_frozen_admission_files(self):
        for key in ('crosswalk', 'model_error'):
            record = CASE['admission_files'][key]
            self.assertEqual(tool.sha(I23 / record['path']), record['sha256'])

    def test_pic_input_template_has_no_hidden_helium_controls(self):
        # The runner appends the explicit He/He+ keys from the case JSON.
        text = (I23 / 'reference/pic_helium_impact.in').read_text()
        for key in ('helium_ion_characteristic_speed_m_s', 'helium_ion_to_sodium_reference_weight_ratio',
                    'helium_to_sodium_reference_weight_ratio', 'electron_density_m3'):
            self.assertNotIn(key + ' =', text)
        self.assertIn('simulation_environment_mode = frozen', text)
        self.assertIn('tracing_clock_mode = common-time', text)


@unittest.skipUnless(sources_available(), 'authenticated helium profile or D05 package unavailable on this host')
class AdmissionGenerator(unittest.TestCase):
    def test_generator_reproduces_frozen_evidence_byte_for_byte(self):
        with tempfile.TemporaryDirectory() as work:
            result = tool.run(Path(CASE['profile_path']), DATA_ROOT, Path(work) / 'out')
            self.assertEqual(result['status'], 'PASS', result)
            for key in ('crosswalk', 'model_error'):
                frozen = I23 / CASE['admission_files'][key]['path']
                self.assertEqual((Path(work) / 'out' / frozen.name).read_bytes(), frozen.read_bytes())

    def test_central_recoil_tolerance_is_refused_without_relaxation(self):
        # The composer evaluates the bound at the upper slope. A tolerance at
        # the central value must be refused by the generator, not widened.
        changed = copy.deepcopy(APPROVAL)
        changed['admission']['recoil']['maximum_heavy_velocity_error_m_s'] = '1073'
        changed['admission']['recoil']['maximum_heavy_kinetic_error_eV'] = '0.0239'
        with self.assertRaisesRegex(ValueError, 'below the composer upper-end bound'):
            tool.derive(CASE, changed, Path(CASE['profile_path']), DATA_ROOT)

    def test_derived_bounds_match_the_approved_values(self):
        crosswalk, model = tool.derive(CASE, APPROVAL, Path(CASE['profile_path']), DATA_ROOT)
        self.assertAlmostEqual(float(crosswalk['source_threshold_difference_eV']), 0.0126, delta=5e-5)
        self.assertLess(float(crosswalk['evaluation_rounding_allowance_J']),
                        1e-15 * float(crosswalk['rate_fit']['threshold_J']))
        self.assertAlmostEqual(float(model['conditional_mean']['admitted_maximum_absolute_log_slope_error']), 0.120,
                               delta=5e-4)
        self.assertLessEqual(float(model['recoil']['upper_end']['rms_speed_bound_m_s']), 1095)
        self.assertLessEqual(float(model['recoil']['upper_end']['kinetic_bound_stationary_eV']), 0.0249)
        self.assertGreater(float(model['recoil']['upper_end']['rms_speed_bound_m_s']),
                           float(model['recoil']['central']['rms_speed_bound_m_s']))
        self.assertEqual(model['mechanical']['heavy_particle_speed_max_m_s'], '3e6')

    def test_unapproved_record_is_skipped(self):
        changed = copy.deepcopy(APPROVAL)
        changed['approval_status'] = 'PROPOSED'
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / 'approval.json'
            path.write_text(json.dumps(changed))
            original = tool.APPROVAL
            try:
                tool.APPROVAL = path
                result = tool.run(Path(CASE['profile_path']), DATA_ROOT, Path(work) / 'out')
            finally:
                tool.APPROVAL = original
            self.assertEqual(result['status'], 'SKIPPED')


if __name__ == '__main__':
    unittest.main()
