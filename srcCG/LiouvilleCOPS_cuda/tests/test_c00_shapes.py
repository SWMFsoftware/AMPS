"""The synthetic C00 gravity shapes are closed, consistently and outwardly oriented."""

import math
from pathlib import Path
import sys

import numpy as np
import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT / "tools"))
import c00_capture_reference as c00  # noqa: E402


@pytest.mark.parametrize("name,shape,volume", [
    ("cube", lambda: c00.cube(1000.0), 1e9),
    ("tetrahedron", lambda: c00.tetrahedron(1000.0), 1000.0 ** 3 / (6.0 * math.sqrt(2.0))),
])
def test_polyhedra_are_closed_outward_with_exact_volume(name, shape, volume):
    v, f = shape()
    topology = c00.mesh_topology(v, f)
    assert topology["consistently_oriented_closed_manifold"]
    assert topology["euler_characteristic"] == 2
    assert topology["signed_volume_m3"] == pytest.approx(volume, rel=1e-12)


@pytest.mark.parametrize("level", range(4))
def test_icosphere_refinement_converges_to_sphere_volume(level):
    v, f = c00.icosphere(1000.0, level)
    topology = c00.mesh_topology(v, f)
    assert topology["faces"] == 20 * 4 ** level
    assert topology["consistently_oriented_closed_manifold"]
    assert np.allclose(np.linalg.norm(v, axis=1), 1000.0, rtol=1e-14)
    sphere = 4.0 / 3.0 * math.pi * 1000.0 ** 3
    assert 0.0 < topology["signed_volume_m3"] < sphere


def test_probe_points_include_exact_singular_locations():
    v, f = c00.cube(1000.0)
    labels, points = c00.probe_points(v, f)
    assert len(labels) == len(points)
    vertex = points[labels.index("on_vertex_0")]
    assert any(np.array_equal(vertex, row) for row in v)


def test_rotation_matrix_is_orthonormal():
    r = c00.rotation_matrix(0.3, -0.7, 1.1)
    assert np.allclose(r @ r.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(r) == pytest.approx(1.0, abs=1e-15)
