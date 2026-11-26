#!/usr/bin/env python3
"""Controls for the approved I24 total-capture inputs and independent oracle.

No MPI process or particle is launched. The parser negative control edits a
private copy of the retained NIST response; nothing here is observational
evidence. The generator and oracle are those used by run_total_capture.py.
"""
from decimal import Decimal as D
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

TEST = Path(__file__).resolve().parent
ROOT = TEST.parents[1]
I24 = TEST / 'integration-tests/I24'
sys.path.insert(0, str(I24))
sys.path.insert(0, str(TEST / 'sodium_cohort'))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


generator = load('i24_total_generator', I24 / 'prepare_total_capture.py')
runner = load('i24_total_runner', I24 / 'run_total_capture.py')


def sources_present(case):
    return all((ROOT / entry['path']).exists() for entry in case['sources'].values())


class TotalCaptureInputs(unittest.TestCase):
    def setUp(self):
        self.case = json.loads(generator.CASE.read_text())
        if not sources_present(self.case):
            self.skipTest('retained I24 total-capture sources are unavailable on this machine')

    def test_regeneration_matches_frozen_reviewed_bytes(self):
        outputs = generator.generate(self.case)
        for name, data in outputs.items():
            self.assertEqual(data, (I24 / 'reference' / name).read_bytes())
            self.assertEqual(generator.sha256(data).hexdigest(), self.case['derived_files'][name])

    def test_changed_source_identity_is_refused(self):
        case = json.loads(json.dumps(self.case))
        case['sources']['nist_asd_h_i']['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'changed source identity'):
            generator.generate(case)

    def test_unapproved_record_skips(self):
        with tempfile.TemporaryDirectory() as work:
            record = Path(work) / 'approval.json'
            record.write_text(json.dumps({'status': 'PENDING'}))
            with patch.object(generator, 'APPROVAL', record):
                with self.assertRaises(generator.NotApproved):
                    generator.generate(self.case)

    def test_statistical_weight_must_equal_2J_plus_1(self):
        raw = (ROOT / self.case['sources']['nist_asd_h_i']['path']).read_text()
        damaged = raw.replace('"=""3/2""",4,', '"=""3/2""",2,', 1)
        self.assertNotEqual(raw, damaged)
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / 'levels.csv'
            path.write_text(damaged)
            with self.assertRaisesRegex(ValueError, 'statistical weight'):
                generator.hydrogen_rows(path)


class TotalCaptureOracle(unittest.TestCase):
    def setUp(self):
        self.case = json.loads(generator.CASE.read_text())
        inventory = ROOT / self.case['sources']['mass_inventory']['path']
        if not inventory.is_file():
            self.skipTest('retained NIST mass inventory is unavailable on this machine')
        self.oracle = runner.Oracle(self.case, json.loads(inventory.read_text()),
                                    I24 / 'reference/shingal_table4.csv', I24 / 'reference/hydrogen_levels.csv')

    @runner.precise  # expectations at the oracle's 50-digit context
    def test_shell_fractions_use_printed_rows_and_normalize(self):
        # Shingal 1986 Table 4, 4.0 keV/u: 0.22 | 17.73 + 37.27 | 2.25 + 5.47 + 3.20 A^2.
        total = D('0.22') + D('55.00') + D('10.92')
        self.assertEqual(self.oracle.fractions(D(4000)), [D('0.22') / total, D('55.00') / total, D('10.92') / total])
        between = self.oracle.fractions(D(2200))
        self.assertAlmostEqual(float(sum(between)), 1.0, places=30)
        self.assertTrue(all(a < b < c or a > b > c for a, b, c in zip(self.oracle.fractions(D(1000)), between,
                                                                       self.oracle.fractions(D(4000)))))

    def test_morgan_law_returns_knots_and_interpolates_log_log(self):
        self.assertEqual(self.oracle.sigma(D(3000)), D('75') * D('1e-20'))
        value = self.oracle.sigma(D(4000))
        self.assertTrue(D('66e-20') < value < D('75e-20'))
        with self.assertRaises(ValueError):
            self.oracle.sigma(D(900))

    def test_shell_centroids_lie_inside_their_nist_levels(self):
        e = self.oracle.e
        self.assertEqual(self.oracle.excitation[0], 0)
        self.assertTrue(D('10.19880615024') * e < self.oracle.excitation[1] < D('10.19885151459') * e)
        self.assertTrue(D('12.0874936591') * e < self.oracle.excitation[2] < D('12.0875115582') * e)

    def test_keyed_uniform_reproduces_the_source_clock_draw(self):
        from verify_numerical import source_uniform
        for parent in range(16):
            self.assertEqual(runner.keyed_uniform(731993, 'control', parent, 1, 0),
                             source_uniform(731993, 'control', parent))

    @runner.precise  # expectations at the oracle's 50-digit context
    def test_zero_angle_products_conserve_momentum_and_energy(self):
        g = self.oracle.speed(3000)
        for shell in (1, 2, 3):
            ref = self.oracle.products(g, shell)
            momentum = [self.oracle.ion * a + self.oracle.h * b for a, b in zip(ref['ion'], ref['hydrogen'])]
            incident = [self.oracle.proton * g * c for c in (D('0.6'), D('0.8'), D(0))]
            self.assertTrue(all(abs(m - i) <= abs(i) * D('1e-40') for m, i in zip(momentum, incident)))
            after = (self.oracle.ion * sum(x * x for x in ref['ion']) + self.oracle.h * sum(x * x for x in ref['hydrogen'])) / 2
            self.assertLess(abs(after - self.oracle.proton * g * g / 2 - ref['q']), self.oracle.proton * g * g * D('1e-40'))


if __name__ == '__main__':
    unittest.main()
