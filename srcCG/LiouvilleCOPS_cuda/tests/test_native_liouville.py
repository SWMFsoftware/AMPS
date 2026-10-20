"""C11: liouville_cops_cuda passes the reference self-test and its GPU-vs-CPU self-checks (GPU).

The binary is built by ``make -C native``; the test is skipped when it is absent.
``Mode=selftest`` of the unchanged LiouvilleCOPS driver exercises the GPU
``EvaluateKernels`` (analytic sphere densities and fluxes, spacecraft velocity) and
the GPU cast shadows (convex sphere: no shadow), plus the start-up self-checks
of the GPU BVH traversal and NG weight against the CPU classes.
"""

import os
from pathlib import Path
import subprocess

import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
BINARY = APP_ROOT / "native" / "liouville_cops_cuda"
REFERENCE = APP_ROOT.parent / "LiouvilleCOPS"


@pytest.mark.skipif(not BINARY.is_file(), reason="native/liouville_cops_cuda not built")
def test_selftest_passes_on_the_gpu():
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("no CUDA device")
    env = dict(os.environ, COPS_CUDA_DEVICE="0")
    result = subprocess.run([str(BINARY), "selftest.in"], cwd=str(REFERENCE), env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
    out = result.stdout
    assert result.returncode == 0, out[-2000:]
    assert "SelfTest: ALL PASSED (0 failures)" in out
    assert "liouville_cops_cuda: GPU kernels on device 0" in out
    assert "0 face mismatches" in out                      # GPU FirstHit == cBVH::FirstHit
    assert "self-check NGWeight GPU vs CPU" in out
    assert "FAIL" not in out
