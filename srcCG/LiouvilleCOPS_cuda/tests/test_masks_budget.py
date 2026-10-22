"""Mask scratch inside the device budget (SHAP8 20M chain, 2026-10-21).

The per-facet mask array of ``masks.visibility_masks`` (uint32 per facet: 80 MB at 20M facets)
was allocated outside every plan and exceeded the 64 MB slack of the dataset's budget.
``memory.scratch_reservation`` reserves such arrays for the call; up to half the slack they stay
in the slack as before. With a tiny slack the synthetic run reproduces the 20M situation.
"""

import numpy as np
import pytest

from cops_cuda.errors import CapacityError
import synthetic_run

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import masks  # noqa: E402
from cops_cuda.memory import DeviceBudget, MemoryPlan, scratch_reservation  # noqa: E402


def test_scratch_beyond_the_slack_is_reserved(clean_pool):
    budget = DeviceBudget(0, slack_bytes=2048)            # pool blocks are 512 B
    try:
        with pytest.raises(CapacityError):                  # the old behaviour: outside the plan
            with budget.allocation_guard(MemoryPlan("unplanned")):
                cp.zeros(4096, dtype=cp.uint32)
        with scratch_reservation(budget, "test", {"bits": (4096, np.uint32)}):
            assert len(budget.reservations) == 1
            bits = cp.zeros(4096, dtype=cp.uint32)
            del bits
        assert not budget.reservations
        with scratch_reservation(budget, "test", {"small": (64, np.uint8)}):
            assert not budget.reservations                  # within half the slack: no reservation
    finally:
        budget.close()


def test_visibility_masks_with_a_small_slack(tmp_path, clean_pool):
    info = synthetic_run.make_run(tmp_path / "run", nfaces=4000)
    points = info["point_ids"]
    want, _ = masks.visibility_masks(info["run"], points)
    budget = DeviceBudget(0, slack_bytes=4096)              # nfaces * 4 B = 16000 B > 4096
    try:
        got, meta = masks.visibility_masks(info["run"], points, budget=budget)
        assert not budget.reservations
    finally:
        budget.close()
    assert meta["face_counts"]["ng_seen"] > 0
    for name in want:
        assert np.array_equal(want[name], got[name]), name
