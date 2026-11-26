#!/usr/bin/env python3
"""Controls for the I21 decision-1 oracles, case generator and gravity recipe.

The oracle checks use small synthetic plan/receipt dictionaries to test the
independent reference code against brute-force evaluations; no native run is
involved and nothing here is physical evidence. Approved 2026-11-26.
"""
import importlib.util
from decimal import Decimal
from pathlib import Path
import sys
import unittest

TEST = Path(__file__).resolve().parent
ROOT = TEST.parents[1]
sys.path.insert(0, str(TEST))
sys.path.insert(0, str(TEST / 'sodium_cohort'))
import verify_numerical as numerical  # noqa: E402
import m3_qualification as q  # noqa: E402

try:
    import mpmath
except ImportError:  # pragma: no cover - the maintained interpreter has mpmath
    mpmath = None


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AU = 149597870700
R_MOON = 1737400.0
SPHERE = [{'centre_lso_m': ['0', '0', '0'], 'radius_m': repr(R_MOON)}]


def shadow_fixture():
    velocity = ['0', '122.0703125', '0']
    plan = {'horizon_s': 65536.0, 'trajectory_reference': {
        'oracle': 'moon-shadow-crossing-zero-force-straight-line-v1', 'source_velocity_lso_ms': velocity}}
    receipt = {'case_kind': 'trajectory-integrated-hazard', 'photon_source_geometry': 'spice-sun',
               'source_position_lso_m': ['-4000000', '-4000000', '0'], 'sun_position_lso_m': [str(AU), '0', '0'],
               'trajectory_reference': {'zero_force': True, 'lunar_point_mass_gravity': False, 'terrain_faces': 0,
                                        'source_velocity_lso_ms': velocity, 'opaque_spheres': SPHERE}}
    return plan, receipt


def kepler_fixture(position=('0', '10000000', '0'), velocity=('0', '0', '600')):
    velocity = list(velocity)
    plan = {'horizon_s': 16384.0, 'trajectory_reference': {
        'oracle': 'lunar-point-mass-kepler-v1', 'source_velocity_lso_ms': velocity}}
    receipt = {'case_kind': 'trajectory-lunar-gravity', 'photon_source_geometry': 'spice-sun',
               'source_position_lso_m': list(position), 'sun_position_lso_m': [str(AU), '0', '0'],
               'trajectory_reference': {'zero_force': False, 'lunar_point_mass_gravity': True,
                                        'lunar_gm_m3_s2': '4902800000000', 'terrain_faces': 0,
                                        'source_velocity_lso_ms': velocity, 'opaque_spheres': SPHERE}}
    return plan, receipt


@unittest.skipIf(mpmath is None, 'mpmath unavailable; run with /usr/bin/python3')
class ShadowOracle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        plan, receipt = shadow_fixture()
        cls.oracle = numerical.ShadowCrossingExposure(plan, receipt, Decimal('7.26e-6'))
        cls.mp = mpmath.MPContext()
        cls.mp.dps = 30

    def brute_depth(self, t):
        o = self.oracle; mp = self.mp
        def lit_rate(s):
            point = [x + v * s for x, v in zip(o.start, o.velocity)]
            if numerical._dark(o.mp, point, o.sun, o.centre, o.radius):
                return 0
            return o.strength / o.mp.fsum((a - b) ** 2 for a, b in zip(o.sun, point))
        nodes = [0, o.entry, o.exit, t] if t > o.exit else [0, min(t, o.entry)]
        return mp.quad(lambda s: lit_rate(mp.mpf(s)), sorted(set(n for n in nodes if n <= t)))

    def test_entry_and_exit_bracket_the_geometric_shadow(self):
        # A point Sun at 1 AU on +x: the shadow is slightly wider than the
        # sphere, so entry precedes (4e6-R)/v and exit follows (4e6+R)/v.
        v = 8e6 / 65536
        self.assertLess(float(self.oracle.entry), (4e6 - R_MOON) / v)
        self.assertGreater(float(self.oracle.exit), (4e6 + R_MOON) / v)
        self.assertLess(float(self.oracle.exit - self.oracle.entry) / (2 * R_MOON / v) - 1, 1e-4)

    def test_closed_form_exposure_matches_direct_integration(self):
        for t in (self.oracle.entry / 2, self.oracle.horizon):
            closed = self.oracle.first_depth if t > self.oracle.exit else self.oracle.strength * (
                self.oracle.F(t) - self.oracle.F(0))
            if t > self.oracle.exit:
                closed = self.oracle.end_depth
            self.assertLess(abs(closed / self.brute_depth(t) - 1), 1e-12)

    def test_first_hit_inverts_depth_and_never_lands_in_the_dark(self):
        for depth in ('0.01', '0.1', str(self.oracle.first_depth + 0.01), '0.26'):
            hit = self.oracle.first_hit(Decimal(depth))
            self.assertIsNotNone(hit)
            t = self.oracle.mp.mpf(str(hit[0]))
            self.assertFalse(self.oracle.in_dark(hit[0]))
            if t <= self.oracle.entry:
                value = self.oracle.strength * (self.oracle.F(t) - self.oracle.F(0))
            else:
                value = self.oracle.first_depth + self.oracle.strength * (self.oracle.F(t) - self.oracle.F(self.oracle.exit))
            self.assertLess(abs(value - self.oracle.mp.mpf(depth)), self.oracle.mp.mpf(10) ** -35)
        self.assertIsNone(self.oracle.first_hit(Decimal('0.5')))


@unittest.skipIf(mpmath is None, 'mpmath unavailable; run with /usr/bin/python3')
class KeplerOracle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        plan, receipt = kepler_fixture()
        cls.oracle = numerical.KeplerExposure(plan, receipt, Decimal('7.26e-6'))

    def test_two_body_invariants(self):
        o = self.oracle; mp = o.mp
        def invariants(t):
            x = o.position(mp.mpf(t)); h = mp.mpf(1) / 2 ** 20
            xp = o.position(mp.mpf(t) + h); xm = o.position(mp.mpf(t) - h)
            v = [(a - b) / (2 * h) for a, b in zip(xp, xm)]
            r = mp.sqrt(mp.fsum(c * c for c in x))
            energy = mp.fsum(c * c for c in v) / 2 - o.gm / r
            momentum = [x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0]]
            return energy, momentum
        e0, h0 = invariants(1)
        for t in (4096, 9000, 16384):
            e, h = invariants(t)
            self.assertLess(abs(e / e0 - 1), 1e-20)
            self.assertLess(max(abs(a - b) for a, b in zip(h, h0)) / abs(h0[0]), 1e-20)

    def test_kepler_position_matches_an_independent_ode_solution(self):
        # Classical RK4 at a 0.5 s step (8,192 steps to 4,096 s): truncation
        # error ~(h*omega)^4*r is far below the 1e-3 m check; independent of
        # the oracle's universal-variable Kepler solution.
        gm = 4902800000000.0
        def rhs(y):
            r3 = (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 1.5
            return [y[3], y[4], y[5], -gm * y[0] / r3, -gm * y[1] / r3, -gm * y[2] / r3]
        y = [0.0, 1.0e7, 0.0, 0.0, 0.0, 600.0]; h = 0.5
        for _ in range(8192):
            k1 = rhs(y); k2 = rhs([a + h / 2 * b for a, b in zip(y, k1)])
            k3 = rhs([a + h / 2 * b for a, b in zip(y, k2)]); k4 = rhs([a + h * b for a, b in zip(y, k3)])
            y = [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]
        actual = self.oracle.position(self.oracle.mp.mpf(4096))
        self.assertLess(max(abs(float(a) - b) for a, b in zip(actual, y[:3])), 1e-3)

    def test_first_hit_inverts_the_integrated_exposure(self):
        depth = Decimal('0.05')
        t = self.oracle.mp.mpf(str(self.oracle.first_hit(depth)[0]))
        self.assertLess(abs(self.oracle.depth(t) - self.oracle.mp.mpf('0.05')), self.oracle.mp.mpf(10) ** -24)
        self.assertIsNone(self.oracle.first_hit(Decimal('0.2')))

    def test_shadowed_orbit_is_refused(self):
        plan, receipt = kepler_fixture(position=('-10000000', '0', '0'), velocity=('0', '600', '0'))
        with self.assertRaises(q.MissingPrerequisite):
            numerical.KeplerExposure(plan, receipt, Decimal('7.26e-6'))


class CaseGeneratorAndRecipe(unittest.TestCase):
    def test_case_files_reproduce_byte_for_byte(self):
        generator = load('i21_decision1_cases', TEST / 'integration-tests/I21/prepare_decision1_cases.py')
        files = generator.render()
        for name, text in files.items():
            self.assertEqual((generator.FROZEN / name).read_text(), text, name)
        self.assertEqual(sorted(p.name for p in generator.FROZEN.iterdir()), sorted(files))

    def test_gravity_reference_recipe_selects_the_production_term(self):
        recipe = load('compile_recipe_prepare', ROOT / 'utility/chemical_database/state_profiles/compile_recipe/prepare.py')
        template = (ROOT / 'input/moon.input').read_text()
        text = recipe.reference_template(template, 'lunar-point-mass-gravity')
        self.assertIn('Moon::LunarPointMassReferenceAcceleration(a,s,p,x,v,n)', text)
        self.assertIn('define _MOON_COHORT_LUNAR_GRAVITY_REFERENCE_ 1', text)
        self.assertNotIn('_MOON_COHORT_ZERO_FORCE_REFERENCE_', text)
        with self.assertRaises(ValueError):
            recipe.reference_template(text, 'homogeneous-zero-acceleration')
        zero = recipe.reference_template(template, 'homogeneous-zero-acceleration')
        self.assertIn('PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)', zero)


if __name__ == '__main__':
    unittest.main()
