"""Target-split gravity (c04_gravity_equivalence.py gpu --devices): bitwise the one-device field.

Each target's field is a fixed-order sum over all body faces, independent of the other
targets and of the target chunking, so contiguous blocks evaluated separately (one per
device, in spawned processes) concatenate to exactly the one-device result.
"""

import argparse

import numpy as np
import pytest

import c00_capture_reference as c00
from test_derive_surface import write_bdf

cp = pytest.importorskip("cupy")
pytestmark = pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")


def body_file(tmp_path):
    v, f = c00.icosphere(1000.0, 3)
    path = tmp_path / "body.bdf"
    write_bdf(path, v, f, extra_cards=False)
    return path


def targets(n=5000, seed=3):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    return d / np.linalg.norm(d, axis=1)[:, None] * rng.uniform(1000.5, 5000.0, n)[:, None]


def test_blocks_and_chunking_do_not_change_the_field(tmp_path, clean_pool):
    from cops_cuda import gravity, gravity_body
    body = gravity_body.GravityBody.from_bdf(body_file(tmp_path), 533.0)
    pts = targets()
    with gravity.GravityEvaluator(body, device_id=0) as ev:
        p0, g0, m0 = ev.evaluate(pts, variant="fold", with_magnitude=True)
    with gravity.GravityEvaluator(body, device_id=0, target_chunk=1 << 10) as ev:
        parts = [ev.evaluate(pts[a:b], variant="fold", with_magnitude=True)
                 for a, b in ((0, 1777), (1777, 5000))]
    assert np.concatenate([p[0] for p in parts]).tobytes() == p0.tobytes()
    assert np.concatenate([p[1] for p in parts]).tobytes() == g0.tobytes()
    assert np.concatenate([p[2] for p in parts]).tobytes() == m0.tobytes()


def test_tool_split_equals_one_device(tmp_path, clean_pool):
    import c04_gravity_equivalence as tool
    from cops_cuda import gravity, gravity_body
    shape = body_file(tmp_path)
    pts = targets(3000)
    devices = [0, 1] if cp.cuda.runtime.getDeviceCount() > 1 else [0, 0]
    args = argparse.Namespace(gravity_shape=str(shape), density=533.0)
    p, g, m, partition = tool.evaluate_split(args, pts, devices, tmp_path)
    body = gravity_body.GravityBody.from_bdf(shape, 533.0)
    with gravity.GravityEvaluator(body, device_id=0) as ev:
        p0, g0, m0 = ev.evaluate(pts, variant="fold", with_magnitude=True)
    assert p.tobytes() == p0.tobytes() and g.tobytes() == g0.tobytes() and m.tobytes() == m0.tobytes()
    assert partition["bounds"] == [0, 1500, 3000]
    assert not list(tmp_path.glob("part_*.npz"))
