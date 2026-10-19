"""C03 exact FP64 polyhedral gravity versus the pinned package (GPU).

The reference is ``polyhedral_gravity`` 3.3 from the locked CUDA environment,
evaluated with ``DISABLE`` on audited outward meshes, exactly as the C00
reference. Budgets (VALIDATION.md, C03 amendment):

* acceleration and potential: ``|d| <= 16 eps M`` where the package result is
  finite. M is the term magnitude accumulated by the kernel, so the bound holds
  under cancellation (far field, symmetric centres);
* points where the package returns NaN must be NaN on the GPU too (same branch
  structure); they are never "fixed".
"""

import math
from pathlib import Path
import sys

import numpy as np
import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT / "tools"))
import c00_capture_reference as c00  # noqa: E402

from cops_cuda import gravity, gravity_body, memory  # noqa: E402
from cops_cuda.errors import CapacityError, IdentityError  # noqa: E402

pg = pytest.importorskip("polyhedral_gravity")
EPS = np.finfo(float).eps
BUDGET = 16.0 * EPS
DENSITY = 533.0
SHAPES = {
    "cube": lambda: c00.cube(1000.0),
    "tetrahedron": lambda: c00.tetrahedron(1000.0),
    "icosphere_l0": lambda: c00.icosphere(1000.0, 0),
    "icosphere_l2": lambda: c00.icosphere(1000.0, 2),
    "icosphere_l3": lambda: c00.icosphere(1000.0, 3),
}


def package(v, f, points, density=DENSITY):
    poly = pg.Polyhedron(polyhedral_source=(v, f), density=density,
                         normal_orientation=pg.NormalOrientation.OUTWARDS,
                         integrity_check=pg.PolyhedronIntegrity.DISABLE,
                         metric_unit=pg.MetricUnit.METER)
    evaluable = pg.GravityEvaluable(polyhedron=poly)
    result = evaluable(computation_points=np.asarray(points).tolist(), parallel=False)
    return (np.array([r[0] for r in result]), np.array([r[1] for r in result]), evaluable)


@pytest.fixture(scope="module")
def device():
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("no CUDA device")
    return cp


def evaluator(cp, body, **kwargs):
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    pool.set_limit(size=0)
    return gravity.GravityEvaluator(body, **kwargs)


@pytest.mark.parametrize("name", sorted(SHAPES))
def test_precomputed_face_terms_are_bitwise_the_package_state(name):
    v, f = SHAPES[name]()
    _, _, evaluable = package(v, f, [[1e5, 0.0, 0.0]])
    _, segments, normals, segment_normals = evaluable.__getstate__()
    G, N, n = gravity.face_terms(v, f)
    assert np.asarray(segments).tobytes() == G.tobytes()
    assert np.asarray(normals).tobytes() == N.tobytes()
    assert np.asarray(segment_normals).tobytes() == n.tobytes()


@pytest.mark.parametrize("variant", ["fold", "split"])
@pytest.mark.parametrize("name", sorted(SHAPES))
def test_gravity_matches_package_within_term_magnitude_budget(device, name, variant):
    v, f = SHAPES[name]()
    labels, points = c00.probe_points(v, f)
    p_ref, g_ref, _ = package(v, f, points)
    body = gravity_body.GravityBody(v, f, DENSITY)
    with evaluator(device, body, split_chunks=5) as ev:
        p_gpu, g_gpu, M = ev.evaluate(points, variant=variant, with_magnitude=True)
    finite = np.isfinite(g_ref).all(axis=1) & np.isfinite(p_ref)
    assert np.array_equal(finite, np.isfinite(g_gpu).all(axis=1) & np.isfinite(p_gpu)), \
        "NaN pattern differs from the package"
    err = np.linalg.norm(g_gpu - g_ref, axis=1)[finite] / M[finite]
    assert err.max() <= BUDGET, (labels[int(np.argmax(err))], err.max() / EPS)
    # Potential: the same face terms times sigma_p h_p, bounded by M times the
    # largest h_p, which is at most the distance to the farthest vertex.
    reach = np.max(np.linalg.norm(v[None, :, :] - points[:, None, :], axis=2), axis=1)
    perr = np.abs(p_gpu - p_ref)[finite] / (M * reach)[finite]
    assert perr.max() <= BUDGET


def test_fold_is_invariant_to_chunking_and_repeatable(device):
    v, f = SHAPES["icosphere_l3"]()
    _, points = c00.probe_points(v, f)
    body = gravity_body.GravityBody(v, f, DENSITY)
    with evaluator(device, body) as ev:
        p1, g1 = ev.evaluate(points)
        p2, g2 = ev.evaluate(points)
    with evaluator(device, body, target_chunk=7, max_pairs_per_launch=7 * 37) as ev:
        p3, g3 = ev.evaluate(points)
    for a, b in ((g1, g2), (g1, g3), (p1, p2), (p1, p3)):
        assert a.tobytes() == b.tobytes()


def test_split_is_repeatable_and_close_to_fold(device):
    v, f = SHAPES["icosphere_l3"]()
    _, points = c00.probe_points(v, f)
    body = gravity_body.GravityBody(v, f, DENSITY)
    with evaluator(device, body, split_chunks=9, max_pairs_per_launch=500) as ev:
        _, a, M = ev.evaluate(points, variant="split", with_magnitude=True)
        _, b = ev.evaluate(points, variant="split")
        _, c = ev.evaluate(points, variant="fold")
    assert a.tobytes() == b.tobytes()
    finite = np.isfinite(c).all(axis=1)          # exact on-vertex probes may be NaN
    assert np.array_equal(finite, np.isfinite(a).all(axis=1))
    assert (np.linalg.norm(a - c, axis=1)[finite] / M[finite]).max() <= BUDGET


def test_physical_invariances(device):
    v, f = SHAPES["icosphere_l3"]()
    labels, points = c00.probe_points(v, f)
    shift = np.array([1234.5, -987.25, 456.125])
    rotation = c00.rotation_matrix(0.3, -0.7, 1.1)
    body = gravity_body.GravityBody(v, f, DENSITY)
    with evaluator(device, body) as ev:
        _, g, M = ev.evaluate(points, with_magnitude=True)
    with evaluator(device, gravity_body.GravityBody(v + shift, f, DENSITY)) as ev:
        _, g_t = ev.evaluate(points + shift)
    with evaluator(device, gravity_body.GravityBody(v @ rotation.T, f, DENSITY)) as ev:
        _, g_r = ev.evaluate(points @ rotation.T)
    with evaluator(device, gravity_body.GravityBody(v, f, 2.0 * DENSITY)) as ev:
        _, g_2 = ev.evaluate(points)
    finite = np.isfinite(g).all(axis=1) & np.isfinite(g_t).all(axis=1) & np.isfinite(g_r).all(axis=1)
    scale = np.max(np.linalg.norm(g[finite], axis=1))
    # Translation/rotation change rounding only (coordinates are re-rounded).
    assert np.max(np.abs(g_t - g)[finite]) / scale < 1e-9
    assert np.max(np.abs(g_r - g @ rotation.T)[finite]) / scale < 1e-9
    # Density enters only through the final prefix: exact up to 1 rounding.
    assert np.max(np.abs(g_2 - 2.0 * g)[finite] / np.abs(2.0 * g)[finite].clip(1e-300)) <= 2 * EPS
    # Far field tends to a point mass of the polyhedron's volume.
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
    gm = gravity.GRAVITATIONAL_CONSTANT * DENSITY * volume
    for i, label in enumerate(labels):
        if label.startswith("far_50"):
            r2 = float(np.dot(points[i], points[i]))
            assert abs(np.linalg.norm(g[i]) * r2 / gm - 1.0) < 2e-3
            # attraction points toward the body
            assert np.dot(g[i], points[i]) < 0.0


def test_inadmissible_bodies_are_refused():
    v, f = c00.cube(1000.0)
    flipped = f.copy()
    flipped[3] = flipped[3][[1, 0, 2]]
    with pytest.raises(IdentityError, match="consistently oriented"):
        gravity_body.GravityBody(v, flipped, DENSITY)
    with pytest.raises(IdentityError, match="inward"):
        gravity_body.GravityBody(v, f[:, [1, 0, 2]], DENSITY)
    with pytest.raises(IdentityError, match="closed"):
        gravity_body.GravityBody(v, f[:-1], DENSITY)
    with pytest.raises(IdentityError, match="outside the vertex table"):
        gravity_body.GravityBody(v, np.where(f == 7, 8, f), DENSITY)


def test_capacity_refusal_before_upload(device):
    cp = device
    v, f = SHAPES["icosphere_l3"]()
    body = gravity_body.GravityBody(v, f, DENSITY)
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    pool.set_limit(size=0)
    total = cp.cuda.runtime.memGetInfo()[1]
    budget = memory.DeviceBudget(0, capacity=(int(0.2 * total) + (1 << 20), total))
    try:
        with pytest.raises(CapacityError):
            gravity.GravityEvaluator(body, budget=budget, target_chunk=1 << 20)
        assert pool.total_bytes() == 0
    finally:
        budget.close()


def test_interior_of_refined_sphere_converges(device):
    errors = []
    for level in (1, 2, 3):
        v, f = c00.icosphere(1000.0, level)
        point = np.array([[100.0, 200.0, -150.0]])
        with evaluator(device, gravity_body.GravityBody(v, f, DENSITY)) as ev:
            _, g = ev.evaluate(point)
        a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
        volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
        radius = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
        del radius   # the inscribed polyhedron is not a sphere: convergence only
        analytic = 4.0 / 3.0 * math.pi * gravity.GRAVITATIONAL_CONSTANT * DENSITY * np.linalg.norm(point)
        errors.append(abs(np.linalg.norm(g) / analytic - 1.0))
    assert errors[0] > errors[1] > errors[2]
