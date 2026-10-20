"""Shared pytest fixtures: import paths, synthetic runs, CuPy pool hygiene."""

from pathlib import Path
import sys

import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
for path in (APP_ROOT, APP_ROOT / "tests", APP_ROOT / "tools"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import synthetic_run  # noqa: E402


@pytest.fixture
def synthetic(tmp_path):
    """A fresh valid synthetic base run (see tests/synthetic_run.py)."""
    return synthetic_run.make_run(tmp_path / "run")


@pytest.fixture
def clean_pool():
    """Empty, unlimited CuPy pools before and after a GPU test.

    ``DeviceBudget`` requires an empty pool and sets a pool limit. A failing
    test must not leave either behind for the next test.
    """
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("no CUDA device")

    def reset():
        pool = cp.get_default_memory_pool()
        pool.free_all_blocks()
        pool.set_limit(size=0)
        cp.get_default_pinned_memory_pool().free_all_blocks()

    reset()
    yield cp
    reset()
