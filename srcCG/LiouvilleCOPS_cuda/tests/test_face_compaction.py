"""Face compaction (SCALING_PLAN.md 3.5): per-face device data only for faces in the records.

* the device faces are exactly the union of the loaded records' faces, in ascending order;
* pressures, every Jacobian column and ``seen_faces`` are BITWISE those without compaction
  (FP64 sequential and tree, FP32 records), including a template basis and a basis rebind;
* the restart kernels (point changes, facet shares) give bitwise the uncompacted results;
* the round-trip digest check still proves the uploaded bytes (faces mapped back);
* the per-face device memory shrinks with the face count.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda import formats, host_model
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator, masks  # noqa: E402
from cops_cuda.device_store import DeviceDataset  # noqa: E402

cpu_model = contract.cpu_model


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def subset_model(env, rows):
    fit_cops = env["fit_cops"]
    run = str(env["synthetic"]["run"])
    header, names, table = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {n: i for i, n in enumerate(names)}
    M = host_model.header_model(run, header, cols, table, env["elements"], env["coef"])
    M.load_kernels(np.asarray(rows, dtype=int))
    return M


def record_faces(run, point_ids):
    faces = []
    for p in point_ids:
        path = run / "kernels" / ("point=%04d.bin" % p)
        raw = np.fromfile(path, dtype=formats.KERNEL_RECORD_DTYPE,
                          offset=formats.KERNEL_HEADER_DTYPE.itemsize)
        faces.append(raw["face"])
    return np.unique(np.concatenate(faces))


def test_dataset_faces_and_roundtrip(cpu_model, clean_pool):
    run = cpu_model["synthetic"]["run"]
    points = list(cpu_model["synthetic"]["point_ids"][2:4])      # an empty point and 118 records
    with DeviceDataset.load(run, point_ids=points, device_id=0, mesh_fields=("ACO2", "area"),
                            compact_faces=True) as ds:
        want = record_faces(run, points)
        assert np.array_equal(ds.local_faces, want) and ds.nfaces_device == len(want)
        assert len(want) < ds.inputs.nfaces                    # something was dropped
        assert ds.arrays["mesh.ACO2"].size == len(want)
        assert same(ds.arrays["mesh.area"].get(), ds.inputs.elements["area"][want])
        ds.verify_roundtrip()                                    # face bytes mapped back
        assert ds.audit()["face_out_of_range"] == 0
        assert ds.identity()["face_compaction"]["device_faces"] == len(want)
        assert "compact.map" not in ds.arrays


@pytest.mark.parametrize("mode,precision", [("sequential", "fp64"), ("tree", "fp64"), ("tree", "fp32")])
def test_results_are_bitwise_uncompacted(cpu_model, clean_pool, mode, precision):
    if precision == "fp32" and mode != "tree":
        pytest.skip("fp32 records: fused tree mode only")
    M = subset_model(cpu_model, [2, 3])
    run = cpu_model["synthetic"]["run"]
    nf = len(cpu_model["elements"]["area"])
    rng = np.random.default_rng(3)
    basis = rng.normal(size=(nf, 2))
    request = bk.ColumnRequest(range(25), range(2), temperature=True, kappa=True)
    states = [bk.ParameterState(cpu_model["coef"][0] * s, M.T0, k, b)
              for s, k, b in ((1.0, 1.0, [0.2, -0.1]), (0.7, 0.4, [1.0, 0.5]))]
    out = {}
    for compact in (False, True):
        with evaluator.CudaBackend(M, run, None, 1.0, mode=mode, records_per_block=17,
                                   record_precision=precision, compact_faces=compact) as be:
            assert (be.dataset.local_faces is not None) == compact
            res = [be.evaluate(bk.ParameterState(s.a, s.T, s.kappa)) for s in states]
            be.rebind_basis(basis, 1.7)
            res += [be.evaluate(s, request) for s in states]
            out[compact] = (res, masks.seen_positive(be.dataset),
                            be.plan.total_bytes, be.dataset.nfaces_device)
    (a, seen_a, bytes_a, nf_a), (b, seen_b, bytes_b, nf_b) = out[False], out[True]
    assert same(seen_a, seen_b) and nf_b < nf_a and bytes_b < bytes_a
    for x, y in zip(a, b):
        assert same(x.pNG, y.pNG) and same(x.pRG, y.pRG)
        for key in x.columns:
            assert same(x.columns[key][0], y.columns[key][0]), key
            assert same(x.columns[key][1], y.columns[key][1]), key


def test_restart_kernels_are_bitwise_uncompacted(cpu_model, clean_pool, monkeypatch, tmp_path):
    import test_restart
    from cops_cuda import restart as rs
    from cops_cuda import evaluator as ev
    _, _, r = test_restart.make_restart(cpu_model, tmp_path, test_restart.template_config())
    run = cpu_model["synthetic"]["run"]
    frozen = np.maximum(r.activity_linear, 0.0)
    delta = np.where(np.arange(len(frozen)) % 4 == 1, 0.3 * frozen + 1e-3, 0.0)
    out = {}
    for flag in ("0", "1"):
        monkeypatch.setenv(ev.COMPACT_ENVIRONMENT, flag)
        se = rs.StateEvaluator(run, row_indices=np.arange(2, 4))
        try:
            p = se.pressures(r)
            shares = se.facet_point_shares(r, delta, se.point_ids)
            change = se.point_changes(r, delta)
            out[flag] = (p, shares, change, se._backend.dataset.nfaces_device)
        finally:
            se.release()
    (p0, s0, c0, n0), (p1, s1, c1, n1) = out["0"], out["1"]
    assert n1 < n0
    assert all(same(u, v) for u, v in zip(p0, p1))
    assert all(same(u, v) for u, v in zip(s0, s1))
    assert all(same(u, v) for u, v in zip(c0, c1))
