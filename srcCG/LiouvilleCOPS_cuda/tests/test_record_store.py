"""FP32 record store (cops_cuda/record_store.py, tools/c15_record_store.py; SCALING_PLAN P1).

* the store reproduces the loader's FP32 conversion: device arrays loaded from the store are
  BITWISE those converted from the kernel files (with and without face compaction), so every
  evaluation is bitwise the same (one GPU and two GPUs);
* the trust chain refuses a tampered store file, a kernel file changed after the build, a store
  of another base, an incomplete store, and FP64 use.
"""

import json
import os
import shutil

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda import record_store
from cops_cuda.errors import CacheFormatError, IdentityError
import synthetic_run
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator  # noqa: E402
from cops_cuda.device_store import DeviceDataset  # noqa: E402

cpu_model = contract.cpu_model
RECORD_ARRAYS = ("kernel.face", "kernel.kn", "kernel.kphi", "kernel.kv", "kernel.cos")


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def test_build_manifest(tmp_path):
    info = synthetic_run.make_run(tmp_path / "run")
    m = record_store.build(info["run"], tmp_path / "store", threads=3)
    recorded = json.loads((info["run"] / "run_manifest.json").read_text())["outputs"]["kernels"]
    assert m["kernel_set_sha256"] == recorded["set_sha256"] and m["kernel_set_verified_against_run_manifest"]
    assert m["records"] == recorded["records"] and len(m["points"]) == len(info["point_ids"])
    for e in m["points"]:
        assert (tmp_path / "store" / e["store_file"]).stat().st_size == 20 * e["records"]
    with pytest.raises(FileExistsError):
        record_store.build(info["run"], tmp_path / "store")


@pytest.mark.parametrize("compact", [False, True])
def test_store_load_is_bitwise_the_kernel_load(tmp_path, clean_pool, compact):
    info = synthetic_run.make_run(tmp_path / "run")
    record_store.build(info["run"], tmp_path / "store", threads=2)
    points = info["point_ids"][1:5]
    kw = dict(point_ids=points, device_id=0, mesh_fields=("ACO2",), record_precision="fp32",
              compact_faces=compact, chunk_records=37)
    with DeviceDataset.load(info["run"], **kw) as a:
        da, fa, ka = a.device_digests(), a.local_faces, (a.kernel_set_sha256, a.kernel_entries)
    with DeviceDataset.load(info["run"], record_store=tmp_path / "store", loader_threads=3, **kw) as b:
        db, fb, kb = b.device_digests(), b.local_faces, (b.kernel_set_sha256, b.kernel_entries)
        b.verify_roundtrip()
        assert b.record_source == "record_store" and b.identity()["record_store"]["files_sha256_verified"] == 4
    assert all(da[k] == db[k] for k in RECORD_ARRAYS) and da == db
    assert (fa is None and fb is None) or np.array_equal(fa, fb)
    assert ka == kb


def test_store_evaluations_are_bitwise(cpu_model, tmp_path, clean_pool):
    run = cpu_model["synthetic"]["run"]
    store = tmp_path / "store"
    record_store.build(run, store)
    P = contract.make_problem(cpu_model, contract.CONFIGS[-1])
    request = bk.ColumnRequest.from_problem(P)
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 6))
    out = []
    for source in (None, store):
        with evaluator.from_problem(P, run, mode="tree", record_precision="fp32",
                                    record_store=source) as be:
            out.append((be.evaluate(state), be.evaluate(state, request)))
    (p0, j0), (p1, j1) = out
    assert same(p0.pNG, p1.pNG) and same(p0.pRG, p1.pRG)
    for key in request.keys():
        assert same(j0.columns[key][0], j1.columns[key][0]) and same(j0.columns[key][1], j1.columns[key][1])


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs two CUDA devices")
def test_two_gpus_read_the_store_from_the_environment(cpu_model, tmp_path, clean_pool, monkeypatch):
    from cops_cuda.multi_gpu import MultiGpuBackend
    from test_gpu_model import make_gpu_model
    run = cpu_model["synthetic"]["run"]
    store = tmp_path / "store"
    record_store.build(run, store)
    monkeypatch.setenv("COPS_CUDA_RECORD_PRECISION", "fp32")
    G = make_gpu_model(cpu_model)
    P = contract.make_problem(dict(cpu_model, model=G), contract.CONFIGS[0])
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 2))
    res = []
    for env in ("", str(store)):
        monkeypatch.setenv("COPS_CUDA_RECORD_STORE", env)
        with MultiGpuBackend(G._spec, G, None, 1.0, devices=(0, 1)) as two:
            res.append(two.evaluate(state))
    assert same(res[0].pNG, res[1].pNG) and same(res[0].pRG, res[1].pRG)


def test_refusals(tmp_path, clean_pool):
    info = synthetic_run.make_run(tmp_path / "run")
    store = tmp_path / "store"
    record_store.build(info["run"], store)
    points = info["point_ids"][:2]
    load = lambda **kw: DeviceDataset.load(info["run"], point_ids=points, device_id=0, mesh_fields="none",
                                           record_store=store, **kw)
    with pytest.raises(ValueError):
        load(record_precision="fp64")
    # a tampered store file
    name = record_store.store_file_name(points[0])
    data = bytearray((store / name).read_bytes())
    data[-1] ^= 1
    (store / name).write_bytes(bytes(data))
    with pytest.raises(CacheFormatError, match="SHA-256"):
        load(record_precision="fp32")
    shutil.rmtree(store)
    record_store.build(info["run"], store)
    # a kernel file touched after the build
    path = synthetic_run.kernel_path(info["run"], points[1])
    st = path.stat()
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 1000))
    with pytest.raises(IdentityError, match="changed after"):
        load(record_precision="fp32")
    # another base, and an incomplete store
    other = synthetic_run.make_run(tmp_path / "other", seed=11)
    with pytest.raises(IdentityError):
        DeviceDataset.load(other["run"], point_ids=other["point_ids"][:1], device_id=0, mesh_fields="none",
                           record_precision="fp32", record_store=store)
    (store / record_store.MANIFEST).unlink()
    with pytest.raises(CacheFormatError, match="incomplete"):
        load(record_precision="fp32")
