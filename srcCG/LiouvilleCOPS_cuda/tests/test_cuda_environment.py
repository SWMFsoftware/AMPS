"""C01 device checks that must hold on every visible GPU (CUDA environment only)."""

from pathlib import Path
import sys

import pytest

cp = pytest.importorskip("cupy")

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
sys.path.insert(0, str(APP_ROOT / "tools"))
import cuda_doctor  # noqa: E402
from cops_cuda import environment  # noqa: E402

DEVICES = list(range(cp.cuda.runtime.getDeviceCount()))
pytestmark = pytest.mark.skipif(not DEVICES, reason="no CUDA device")


def test_exactly_one_cupy_distribution():
    assert len(environment.cupy_distributions()) == 1


def test_devices_have_stable_uuids():
    inventory = environment.device_inventory()
    uuids = [d["uuid"] for d in inventory]
    assert all(u and u.startswith("GPU-") for u in uuids)
    assert len(set(uuids)) == len(uuids)


@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("check", ["check_fp64_kernel", "check_async", "check_reduction",
                                   "check_memory"])
def test_device_check(device, check):
    with cp.cuda.Device(device):
        ok, detail = getattr(cuda_doctor, check)(cp, device)
    assert ok, detail


def test_worker_thread_environment_limits_host_pools():
    env = environment.worker_thread_environment(1)
    assert env["OMP_NUM_THREADS"] == "1" and env["OPENBLAS_NUM_THREADS"] == "1"
