"""MPI transport of the multi-GPU backend (cops_cuda/transport.py; SCALING_PLAN P3, invariant I2).

Real ``mpirun`` jobs on the local GPUs (tests/mpi_driver.py; several ranks may share a GPU):

* host activity, 2 ranks: pressures, every Jacobian column (also after a basis rebind) and
  ``seen_faces`` are BITWISE one GPU with host activity;
* device activity, 1, 2 and 3 ranks (rank 0 in process, devices 0,1,0), FP64 and FP32 records,
  streamed: BITWISE one GPU with device activity;
* B29 (no template basis at bind), 2 ranks, both activity modes: BITWISE one GPU;
* a request failing on the workers is a BackendError on rank 0 that names the rank; closing
  collects the stale replies and the same ranks serve the next bind correctly.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("mpi4py")
cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
MPIRUN = shutil.which("mpirun")
if MPIRUN is None:
    pytest.skip("no mpirun", allow_module_level=True)

from cops_cuda import evaluator  # noqa: E402
import mpi_driver  # noqa: E402

DRIVER = Path(__file__).resolve().parent / "mpi_driver.py"


def launch(n, case, run, out, activity, env=None, timeout=600):
    environ = dict(os.environ, **(env or {}))
    cmd = [MPIRUN, "-n", str(n), "--oversubscribe", "--bind-to", "none", sys.executable, "-W", "ignore",
           str(DRIVER), case, str(run), str(out), activity]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=environ,
                       timeout=timeout)
    assert r.returncode == 0, (r.stdout[-3000:], r.stderr[-3000:])
    return np.load(out)


def reference(run, activity, precision="fp64"):
    """One resident GPU, in this process."""
    env, G, P = mpi_driver.environment(run)
    with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale, records_per_block=17,
                               activity=activity, record_precision=precision) as be:
        return dict(mpi_driver.sequence(be, G, P))


def assert_bitwise(want, got):
    names = [k for k in got.files if k != "info"]
    assert sorted(names) == sorted(want)
    for name in names:
        assert np.asarray(want[name]).tobytes() == got[name].tobytes(), name


def test_two_ranks_host_activity_are_bitwise_one_gpu(synthetic, clean_pool, tmp_path):
    run = synthetic["run"]
    want = reference(run, "host")
    got = launch(2, "evaluate", run, tmp_path / "out.npz", "host")
    assert_bitwise(want, got)
    info = json.loads(str(got["info"]))
    assert info["transport"] == "mpi" and info["ranks"] == 2 and info["activity"] == "host"
    assert info["worker_activity"] == ["host", "host"]


@pytest.mark.parametrize("ranks", [1, 2, 3])
def test_ranks_device_activity_are_bitwise_one_gpu(synthetic, clean_pool, tmp_path, ranks):
    run = synthetic["run"]
    want = reference(run, "device")
    got = launch(ranks, "evaluate", run, tmp_path / "out.npz", "device")
    assert_bitwise(want, got)
    info = json.loads(str(got["info"]))
    assert info["ranks"] == ranks and info["worker_activity"] == ["device"] * ranks
    ndev = cp.cuda.runtime.getDeviceCount()
    assert info["devices"] == [r % ndev for r in range(ranks)]
    assert [p["rows"][0] for p in info["partition"]][0] == 0


def test_streamed_fp32_ranks_are_bitwise_one_gpu(synthetic, clean_pool, tmp_path, monkeypatch):
    run = synthetic["run"]
    want = reference(run, "device", "fp32")
    env = {"COPS_CUDA_RECORD_PRECISION": "fp32", evaluator.RESIDENT_ENVIRONMENT: "60",
           evaluator.STREAM_CHUNK_ENVIRONMENT: "40"}
    got = launch(2, "evaluate", run, tmp_path / "out.npz", "device", env)
    assert_bitwise(want, got)
    info = json.loads(str(got["info"]))
    assert all(n > 0 for n in info["streamed_records"])


def test_worker_failure_is_a_backend_error_and_ranks_serve_again(synthetic, clean_pool, tmp_path):
    run = synthetic["run"]
    got = launch(2, "failure", run, tmp_path / "out.npz", "device")
    message = str(got["message"])
    assert "GPU worker" in message and "failed during rebind_basis" in message and "Traceback" in message
    assert bool(got["broken"])
    env, G, P = mpi_driver.environment(run)
    from cops_cuda import backend as bk
    import test_backend_contract as contract
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 4))
    with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale, records_per_block=17,
                               activity="device") as be:
        want = be.evaluate(state)
    assert want.pNG.tobytes() == got["pNG"].tobytes() and want.pRG.tobytes() == got["pRG"].tobytes()


@pytest.mark.parametrize("activity", ["host", "device"])
def test_no_basis_bind_is_bitwise_one_gpu(synthetic, clean_pool, tmp_path, activity):
    from cops_cuda import backend as bk
    import test_backend_contract as contract
    run = synthetic["run"]
    env, G, P = mpi_driver.environment(run, 0)
    request = bk.ColumnRequest(P.fit_j, (), temperature=True, kappa=True)
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 5))
    with evaluator.CudaBackend(G, run, None, 1.0, records_per_block=17, activity=activity) as be:
        r = be.evaluate(state, request)
    want = {"pNG": r.pNG, "pRG": r.pRG}
    for key in sorted(r.columns):
        want["%s.%d.NG" % key], want["%s.%d.RG" % key] = r.columns[key]
    assert_bitwise(want, launch(2, "no_basis", run, tmp_path / "out.npz", activity))
