"""Mixed-precision gravity (experiment, 2026-10-21): exact FP64 near faces, FP32 far quadrature.

The far-field terms replace the exact plane sum by a 3-point quadrature of the face integral
of 1/r; against the exact FP64 field the relative acceleration error must fall with the near
radius and stay far below what the slope thresholds (45, 60 degrees) could notice.
"""

import numpy as np
import pytest

import c00_capture_reference as c00
from test_derive_surface import write_bdf

cp = pytest.importorskip("cupy")
pytestmark = pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")


def body_and_points(tmp_path):
    from cops_cuda import gravity_body
    v, f = c00.icosphere(1000.0, 4)
    path = tmp_path / "body.bdf"
    write_bdf(path, v, f, extra_cards=False)
    body = gravity_body.GravityBody.from_bdf(path, 533.0)
    V, F = body.vertices, body.faces
    c = V[F].mean(axis=1)
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    n /= np.linalg.norm(n, axis=1)[:, None]
    rng = np.random.default_rng(1)
    idx = rng.choice(len(F), 600, replace=False)
    far = rng.normal(size=(200, 3))
    far *= (rng.uniform(1100.0, 5000.0, 200) / np.linalg.norm(far, axis=1))[:, None]
    return body, np.vstack([c[idx] + 0.1 * n[idx], far])


def test_mixed_matches_fp64(tmp_path, clean_pool):
    from cops_cuda import gravity
    body, pts = body_and_points(tmp_path)
    with gravity.GravityEvaluator(body, 0) as ev:
        p0, g0 = ev.evaluate(pts)
    errors = []
    for k in (2.0, 6.0):
        with gravity.GravityEvaluator(body, 0, precision="mixed", near_factor=k) as ev:
            p1, g1 = ev.evaluate(pts)
            near = ev.near_counts
            ident = ev.identity()
        rel = np.linalg.norm(g1 - g0, axis=1) / np.linalg.norm(g0, axis=1)
        cosang = np.sum(g1 * g0, axis=1) / np.linalg.norm(g1, axis=1) / np.linalg.norm(g0, axis=1)
        assert np.degrees(np.arccos(np.clip(cosang, -1, 1))).max() < 1e-3
        assert np.max(np.abs(p1 - p0) / np.abs(p0)) < 1e-5
        assert near[:600].min() >= 1                          # every surface target has exact near faces
        assert ident["precision"]["policy"] == gravity.MIXED_POLICY
        errors.append(rel.max())
    assert errors[0] < 1e-4 and errors[1] < 2e-6 and errors[1] < errors[0]


def test_mixed_chunking_is_bitwise(tmp_path, clean_pool):
    """The FP32 tile sums are fixed by the face order, so launch splits do not change results."""
    from cops_cuda import gravity
    body, pts = body_and_points(tmp_path)
    with gravity.GravityEvaluator(body, 0, precision="mixed") as ev:
        _, a = ev.evaluate(pts)
    with gravity.GravityEvaluator(body, 0, precision="mixed", target_chunk=97,
                                  max_pairs_per_launch=97 * 100) as ev:      # 100 -> 64 faces per launch
        _, b = ev.evaluate(pts)
    assert a.tobytes() == b.tobytes()
