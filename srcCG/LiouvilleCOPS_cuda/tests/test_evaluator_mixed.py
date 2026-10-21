"""Mixed precision evaluator (experiment, 2026-10-21): FP32 kernel records, FP64 arithmetic.

* Pressures and every Jacobian column agree with the FP64 tree mode within 1e-6 of the
  summed magnitudes (each FP32 input carries a relative error <= 2^-24; terms combine at most
  four inputs), and the evaluation is deterministic.
* The seen-facet set and the visibility masks are identical (the loader keeps zero-ness).
* The loader refuses a record whose value FP32 cannot hold (underflow to zero).
* The device plan per record is 21 bytes of records and no stored factors.
"""

import shutil

import numpy as np
import pytest

from cops_cuda import backend as bk
import synthetic_run
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator, masks  # noqa: E402
from c06_jacobian_budget import jacobian_magnitudes  # noqa: E402
from test_evaluator import magnitudes, fresh_pool  # noqa: E402

cpu_model = contract.cpu_model
TOL = 1e-6


@pytest.mark.parametrize("config", contract.CONFIGS, ids=[c["name"] for c in contract.CONFIGS])
def test_mixed_matches_fp64(cpu_model, config):
    fresh_pool()
    P = contract.make_problem(cpu_model, config)
    run = cpu_model["synthetic"]["run"]
    fit_cops = cpu_model["fit_cops"]
    request = bk.ColumnRequest.from_problem(P)
    for seed in (4, 9):
        x = contract.perturbed(P, seed)
        state = bk.ParameterState.from_problem(P, x)
        mNG, mRG = magnitudes(P, x)
        with evaluator.from_problem(P, run, mode="tree", records_per_block=17) as ref:
            want, want_j = ref.evaluate(state), ref.evaluate(state, request)
            seen_ref = masks.seen_positive(ref.dataset)
            ref_bytes = ref.plan.total_bytes + ref.dataset.plan.total_bytes
        fresh_pool()
        with evaluator.from_problem(P, run, mode="tree", records_per_block=17,
                                    record_precision="fp32") as gpu:
            got, got_j = gpu.evaluate(state), gpu.evaluate(state, request)
            assert gpu.evaluate(state).pNG.tobytes() == got.pNG.tobytes()     # deterministic
            assert np.array_equal(masks.seen_positive(gpu.dataset), seen_ref)
            ident = gpu.identity()
            mixed_bytes = gpu.plan.total_bytes + gpu.dataset.plan.total_bytes
        fresh_pool()
        assert np.all(np.abs(got.pNG - want.pNG) <= TOL * mNG)
        assert np.all(np.abs(got.pRG - want.pRG) <= TOL * mRG)
        for key, mags in jacobian_magnitudes(P, P.M, x, request, fit_cops).items():
            for g, w, m in zip(got_j.columns[key], want_j.columns[key], mags):
                assert np.all(np.abs(g - w) <= TOL * m), key
        assert ident["precision"]["record_storage"] == "fp32"
        assert 0 < ident["precision"]["record_conversion"]["max_relative_rounding"] <= 2.0 ** -24
        assert mixed_bytes < ref_bytes


def test_fp32_loader_refuses_underflow(tmp_path, clean_pool):
    from cops_cuda.device_store import DeviceDataset
    from cops_cuda.errors import CacheFormatError
    info = synthetic_run.make_run(tmp_path / "run")
    run = info["run"]
    point = info["point_ids"][0]
    geometry, records = info["kernels"][point]
    records = records.copy()
    records["Kn"][0] = 1e-50                          # nonzero in FP64, 0 in FP32
    synthetic_run.write_kernel(synthetic_run.kernel_path(run, point), 1, point, geometry, records)
    synthetic_run.write_manifest(run, info["point_ids"])
    with DeviceDataset.load(run, point_ids=[point], device_id=0, mesh_fields="none") as ds:
        assert ds.record_precision == "fp64"
    with pytest.raises(CacheFormatError, match="FP32"):
        DeviceDataset.load(run, point_ids=[point], device_id=0, mesh_fields="none",
                           record_precision="fp32")
    shutil.rmtree(run)
