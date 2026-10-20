"""Small deterministic regressions for the metre-unit HEAL corruption repair."""
import sys
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from gravity_topology import require_outward_closed_body, ORIENTATION_POLICY
from gravity_cache import make_definition, write_cache, load_cache


def icosphere(level=3):
    """Build the same outward triangulation at every physical radius."""
    phi = (1 + np.sqrt(5.)) / 2
    v = np.array([(0, a, b * phi) for a in (-1, 1) for b in (-1, 1)] +
                 [(a, b * phi, 0) for a in (-1, 1) for b in (-1, 1)] +
                 [(b * phi, 0, a) for a in (-1, 1) for b in (-1, 1)], dtype=float)
    v /= np.linalg.norm(v, axis=1)[:, None]
    f = ConvexHull(v).simplices.copy()
    inward = np.einsum("ij,ij->i", v[f[:, 0]],
                      np.cross(v[f[:, 1]], v[f[:, 2]])) < 0
    f[inward] = f[inward][:, [1, 0, 2]]
    for _ in range(level):
        verts, split, refined = list(v), {}, []
        def midpoint(i, j):
            key = tuple(sorted((i, j)))
            if key not in split:
                p = v[i] + v[j]
                split[key] = len(verts)
                verts.append(p / np.linalg.norm(p))
            return split[key]
        for a, b, c in f:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        v, f = np.asarray(verts), np.asarray(refined)
    return v, f


class TopologyTests(unittest.TestCase):
    def setUp(self):
        self.v, self.f = icosphere(0)

    def test_closed_outward_scale_and_translation(self):
        small = require_outward_closed_body(self.v, self.f)
        large = require_outward_closed_body(self.v * 1000 + 1e6, self.f)
        self.assertEqual(large["euler_characteristic"], 2)
        self.assertTrue(large["vertex_links_manifold"])
        self.assertAlmostEqual(large["signed_volume_m3"] / small["signed_volume_m3"] / 1e9, 1.)

    def test_refuse_invalid_inputs(self):
        flipped = self.f.copy()
        flipped[0] = flipped[0, ::-1]
        invalid = self.f.copy()
        invalid[0, 0] = len(self.v)
        degenerate = self.f.copy()
        degenerate[0, 0] = degenerate[0, 1]
        nan = self.v.copy()
        nan[0, 0] = np.nan
        for vertices, faces in ((self.v, flipped), (self.v, self.f[:-1]),
                                (self.v, self.f[:, ::-1]), (self.v, invalid),
                                (self.v, degenerate), (nan, self.f),
                                (self.v, self.f.astype(float)),
                                (np.r_[self.v, [[0, 0, 0]]], self.f)):
            with self.subTest(faces=faces.shape):
                with self.assertRaises(ValueError):
                    require_outward_closed_body(vertices, faces)

    def test_refuse_disconnected_body(self):
        with self.assertRaisesRegex(ValueError, "connected body"):
            require_outward_closed_body(np.r_[self.v, self.v + 4],
                                        np.r_[self.f, self.f + len(self.v)])

    def test_old_heal_cache_cannot_match_new_definition(self):
        centroid = np.array([[1., 0., 0.]])
        normal, area, omega = centroid.copy(), np.ones(1), np.zeros(3)
        current = make_definition(centroid, normal, area,
                                  {"sha256": "fixture", "size_bytes": 1}, 12, 20, 533, .1)
        self.assertEqual(current["orientation_policy"], ORIENTATION_POLICY)
        old = dict(current)
        old.pop("orientation_policy")
        old["algorithm"] = "homogeneous_polyhedron_self_gravity_at_outward_offset_v1"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.npz"
            write_cache(path, centroid, normal, area, old, -centroid, omega, {})
            with self.assertRaisesRegex(ValueError, "physical inputs"):
                load_cache(path, centroid, normal, area, current, omega)


try:
    import polyhedral_gravity as pg
except ImportError:
    pg = None


@unittest.skipIf(pg is None, "optional polyhedral-gravity package is not on PYTHONPATH")
class PackageRegressionTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("COPS_GRAVITY_REFERENCE"),
                         "set COPS_GRAVITY_REFERENCE to the C00 gravity fixture")
    def test_c00_surface_targets(self):
        from build_physical_activity_templates import read_nastran_surface
        import json
        fixture = Path(os.environ["COPS_GRAVITY_REFERENCE"])
        manifest = json.loads((fixture.parent / "gravity_fixtures.json").read_text())
        specification = manifest["exact_788k"]
        v, f = read_nastran_surface(specification["shape"]["path"])
        require_outward_closed_body(v, f)
        with np.load(fixture, allow_pickle=False) as archive:
            # Keep the default developer check small; the release acceptance
            # can request all 832 recorded targets explicitly.
            count = int(os.environ.get("COPS_GRAVITY_REFERENCE_TARGETS", "8"))
            selected = np.unique(np.linspace(0, len(archive["points_offset"]) - 1,
                                             count, dtype=int))
            points = archive["points_offset"][selected]
            expected = archive["acceleration_offset"][selected]
        poly = pg.Polyhedron(polyhedral_source=(v, f),
                             density=specification["density_kg_m3"],
                             normal_orientation=pg.NormalOrientation.OUTWARDS,
                             integrity_check=pg.PolyhedronIntegrity.DISABLE,
                             metric_unit=pg.MetricUnit.METER)
        result = pg.GravityEvaluable(polyhedron=poly)(computation_points=points.tolist(),
                                                      parallel=False)
        actual = np.asarray([row[1] for row in result])
        # Parallel package reductions need not be bitwise identical. Compare
        # vector-relative error, avoiding fragile per-component zero divisions.
        relative = np.linalg.norm(actual - expected, axis=1) / np.linalg.norm(expected, axis=1)
        self.assertLessEqual(float(np.max(relative)), 1e-12)

    def test_documented_heal_defect_reproduced_on_pinned_package(self):
        if pg.__version__ != "3.3":
            self.skipTest("defect reproduction specifically targets pinned release 3.3")
        v, f = icosphere(3)
        healed = pg.Polyhedron(polyhedral_source=(1000 * v, f), density=533,
                               normal_orientation=pg.NormalOrientation.OUTWARDS,
                               integrity_check=pg.PolyhedronIntegrity.HEAL,
                               metric_unit=pg.MetricUnit.METER)
        # This is deliberately the only HEAL invocation in this test suite.
        # The exact count depends on triangulation and floating-point rounding.
        self.assertGreater(np.count_nonzero(np.any(np.asarray(healed.faces) != f, axis=1)), 0)

    def test_both_builder_paths_preserve_audited_geometry(self):
        from build_physical_activity_templates import effective_slope, validate_gravity_resolution
        v, f = icosphere(0)
        v *= 1000
        centroid = np.mean(v[f], axis=1)
        normal = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        slope, gravity, _, metadata = effective_slope(
            centroid, normal, v, f, 533, np.zeros(3), .1, .25, 100, False)
        self.assertTrue(np.all(np.sum(gravity * normal, axis=1) < 0))
        self.assertEqual(metadata["orientation_policy"], ORIENTATION_POLICY)
        comparison = validate_gravity_resolution(
            centroid, normal, np.zeros(3), .1, 533, slope, v, f, 1e-6, False)
        self.assertEqual(comparison["topology"], metadata["topology"])
        self.assertLessEqual(comparison["slope_difference_deg_max"], 1e-6)

    def test_metre_kilometre_scale_similarity_and_direction(self):
        # A scale-independent physical identity: at homologous positions and
        # equal density, acceleration scales linearly with body radius. This
        # checks the repaired path without depending on HEAL's unstable count.
        v, f = icosphere(3)
        points = np.array([[1.2, .1, .2], [-.2, 1.3, -.1], [.1, .2, -1.4]])
        def acceleration(scale, density=533):
            require_outward_closed_body(v * scale, f)
            poly = pg.Polyhedron(polyhedral_source=(v * scale, f), density=density,
                                 normal_orientation=pg.NormalOrientation.OUTWARDS,
                                 integrity_check=pg.PolyhedronIntegrity.DISABLE,
                                 metric_unit=pg.MetricUnit.METER)
            result = pg.GravityEvaluable(polyhedron=poly)(
                computation_points=(points * scale).tolist(), parallel=False)
            return np.array([row[1] for row in result])
        unit, large = acceleration(1), acceleration(1000)
        np.testing.assert_allclose(large, 1000 * unit, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(acceleration(1000, 1066), 2 * large, rtol=1e-14)
        self.assertTrue(np.all(np.sum(large * points, axis=1) < 0))


if __name__ == "__main__":
    unittest.main()
