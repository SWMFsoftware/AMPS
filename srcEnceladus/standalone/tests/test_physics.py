import math
import copy
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fitting import load_config, sources_from_config  # noqa: E402
from plume_model import (KB, POINT_SOURCE_LAW, SURFACE_FLUX_LAW, Plume,
                         angular_grid, line_of_sight_design_matrix, line_source,
                         outward_ray_design_matrix, radial_moment,
                         spherical_source)  # noqa: E402


MASS_H2O = 2.991507624607566e-26


class PhysicsTests(unittest.TestCase):
    def test_radial_moments_match_direct_quadrature(self):
        for power in (2, 3):
            for speed_ratio in (0.0, 0.8, 3.0, 8.0):
                for cosine in (-0.9, -0.1, 0.4, 1.0):
                    # Complete the square before numerical quadrature.  Direct
                    # integration can be O(1e-31), far below a conventional
                    # absolute tolerance, and then is not an accurate oracle.
                    shift = speed_ratio*cosine
                    scale = math.exp(-speed_ratio**2*(1-cosine**2))
                    expected = scale*quad(
                        lambda t: t**power*math.exp(-(t-shift)**2),
                        0, np.inf, epsabs=1e-13, epsrel=1e-11)[0]
                    actual = float(radial_moment(power, cosine, speed_ratio))
                    self.assertLess(abs(actual/expected-1.0), 1e-7)

    def test_flux_integrates_to_unit_production_for_both_laws(self):
        position, normal, axis = spherical_source(252100, -80, 40, 25, 120)
        mu, wmu, phi, wphi = angular_grid(192)
        north = np.array([0.0, 0.0, 1.0])
        east = np.cross(north, normal)
        east /= np.linalg.norm(east)
        tangent_north = np.cross(normal, east)
        rays = (mu[:, None, None]*normal
                + np.sqrt(1-mu[:, None, None]**2)
                * (np.cos(phi)[None, :, None]*tangent_north
                   + np.sin(phi)[None, :, None]*east))
        distance = 1.0e6
        points = position+rays.reshape(-1, 3)*distance
        weights = (wmu[:, None]*wphi[None, :]).ravel()
        for law in (POINT_SOURCE_LAW, SURFACE_FLUX_LAW):
            plume = Plume(tuple(position), tuple(normal), tuple(axis), 190, 700,
                          MASS_H2O, law)
            integral = np.sum(plume.radial_flux_per_rate(points)*distance**2*weights)
            self.assertAlmostEqual(integral, 1.0, places=8)

    def test_inverse_square_and_horizon(self):
        position, normal, axis = spherical_source(252100, -90, 0, 0, 0)
        plume = Plume(tuple(position), tuple(normal), tuple(axis), 180, 600,
                      MASS_H2O, POINT_SOURCE_LAW)
        n1 = plume.density_per_rate(position+normal*1e5)
        n2 = plume.density_per_rate(position+normal*2e5)
        self.assertAlmostEqual(float(n1/n2), 4.0, places=12)
        self.assertEqual(float(plume.density_per_rate(position-normal*1e5)), 0.0)

    def test_surface_flux_has_projected_area_factor(self):
        position = np.array([252100.0, 0.0, 0.0])
        normal = axis = np.array([1.0, 0.0, 0.0])
        tangent_point = position+np.array([0.0, 1.0e5, 0.0])
        point = Plume(tuple(position), tuple(normal), tuple(axis), 180, 600,
                      MASS_H2O, POINT_SOURCE_LAW)
        surface = Plume(tuple(position), tuple(normal), tuple(axis), 180, 600,
                        MASS_H2O, SURFACE_FLUX_LAW)
        self.assertGreater(float(point.density_per_rate(tangent_point)), 0)
        self.assertEqual(float(surface.density_per_rate(tangent_point)), 0)

    def test_signed_latitude_and_west_longitude(self):
        position, normal, _ = spherical_source(1, -45, 90, 0, 0)
        expected = np.array([0, -math.sqrt(0.5), -math.sqrt(0.5)])
        np.testing.assert_allclose(position, expected, atol=2e-15)
        np.testing.assert_allclose(normal, expected, atol=2e-15)

    def test_line_source_is_length_normalized(self):
        points = [[252100, 0, 0], [252100, 1000, 0], [252100, 4000, 0]]
        stripe = line_source(points, 180, 350, MASS_H2O)
        np.testing.assert_allclose(stripe.weights, [0.125, 0.5, 0.375])
        self.assertAlmostEqual(sum(stripe.weights), 1.0)

    def test_published_configuration_has_all_sources(self):
        config = load_config(ROOT/"config"/"tenishev2014_stage0.json")
        sources, labels = sources_from_config(config)
        self.assertEqual(len(sources), 15)
        self.assertEqual(sum(len(s.components) for s in sources[8:]), 102)
        self.assertEqual(config["source_law"], POINT_SOURCE_LAW)
        self.assertEqual(len(labels), 17)

    def test_source_law_is_selectable_in_configuration(self):
        config = load_config(ROOT/"config"/"tenishev2014_stage0.json")
        alternative = copy.deepcopy(config)
        alternative["source_law"] = SURFACE_FLUX_LAW
        sources, _ = sources_from_config(alternative)
        self.assertTrue(all(source.source_law == SURFACE_FLUX_LAW
                            for source in sources[:8]))
        self.assertTrue(all(component.source_law == SURFACE_FLUX_LAW
                            for stripe in sources[8:]
                            for component in stripe.components))

    def test_infinite_chord_power_law_matches_analytic_integral(self):
        reference = 252100.0
        impact = np.array([300000.0, 500000.0, 1000000.0])
        tangent = np.column_stack((impact, np.zeros((len(impact), 2))))
        direction = np.tile([0.0, 1.0, 0.0], (len(impact), 1))
        matrix = line_of_sight_design_matrix(
            [], tangent, direction, reference, order=48)
        expected = np.pi*reference**2/impact
        np.testing.assert_allclose(matrix[:, 0], expected, rtol=2e-14)
        reverse = line_of_sight_design_matrix(
            [], tangent, -direction, reference, order=48)
        np.testing.assert_allclose(reverse, matrix, rtol=0, atol=1e-8)
        observer_coordinate = -np.array([2.0e5, 4.0e5, 8.0e5])
        finite = line_of_sight_design_matrix(
            [], tangent, direction, reference, order=48,
            observer_coordinate_from_tangent_m=observer_coordinate)
        finite_expected = (reference**2/impact
                           *(np.pi/2-np.arctan(observer_coordinate/impact)))
        np.testing.assert_allclose(finite[:, 0], finite_expected, rtol=2e-14)

        outward_origin = np.column_stack((impact, np.zeros((len(impact), 2))))
        outward_direction = np.tile([1.0, 0.0, 0.0], (len(impact), 1))
        outward = outward_ray_design_matrix(
            [], outward_origin, outward_direction, reference, order=48)
        np.testing.assert_allclose(
            outward[:, 0], reference**2/impact, rtol=2e-14)


if __name__ == "__main__":
    unittest.main()
