"""Device-memory planning, headroom policy and tiling (ROADMAP C02 item 3).

Plan arithmetic is pure Python and runs anywhere. Measurement and enforcement
(``DeviceBudget``) import CuPy only when used, so they need the CUDA
environment.

Why an explicit plan
--------------------
CPU resident memory is not a GPU memory estimate. The device also holds:

* the CUDA context of each process, about 0.3 GB measured on the A6000;
* loaded compiled modules;
* the CuPy pool, which keeps freed blocks cached and rounds each block up to
  512 bytes;
* the scratch and partial-sum buffers of the reductions.

Every device array a store allocates is therefore declared first as a
``MemoryItem`` with its exact byte count. The sum is checked against measured
free memory before anything is allocated. If the dataset does not fit, the
refusal names every array and its size instead of failing halfway through a
multi-gigabyte upload.

Headroom policy
---------------
``DEFAULT_HEADROOM_FRACTION`` (0.20) of the device's *total* memory is never
planned (a run may set another fraction in [0.05, 0.5] through
``$COPS_CUDA_HEADROOM_FRACTION``; ``headroom_fraction_default``). This is the initial ROADMAP policy. It leaves room for the
driver, for fragmentation, for an unexpected library workspace, and for a
second process to initialise a context. A plan must fit within
``free - headroom``, both at the moment of the first measurement and at the
moment of each reservation, because another process may have allocated in
between.

Enforcement
-----------
``DeviceBudget`` sets the CuPy default pool limit of its device to the sum of
the reserved, rounded array sizes plus a declared slack for small incidental
temporaries. An allocation outside the plan then raises CuPy's
``OutOfMemoryError`` at the allocation site instead of quietly consuming the
headroom. ``CapacityError`` (a ``MemoryError``) wraps it with the plan.

Tiling
------
A Jacobian over the full contribution set must never form a
contribution-by-parameter array (245M x 30 x 8 B is about 59 GB).
``plan_tiles`` splits observation points into contiguous groups bounded in
record count and the requested parameter columns into groups bounded in
width. The partial-sum scratch of the deterministic segmented reduction then
fits a declared budget. Tiles split only at point boundaries. A single point
larger than the record bound becomes its own tile, and the C06 reduction
splits it over several blocks internally, combining their partial sums in a
fixed order.
"""

import math
import os
from collections import namedtuple

from .errors import CapacityError


DEFAULT_HEADROOM_FRACTION = 0.20
HEADROOM_ENVIRONMENT = "COPS_CUDA_HEADROOM_FRACTION"
HEADROOM_RANGE = (0.05, 0.5)


def headroom_fraction_default():
    """The process's headroom fraction: ``$COPS_CUDA_HEADROOM_FRACTION`` or the default 0.20.

    A run setting (``gpu_fit.py --headroom-fraction``), introduced 2026-10-20 for the SHAP5 3M
    mesh, whose two-GPU evaluator plan exceeded the 20% policy by 0.8 GiB while 8.7 GiB of
    the device stayed actually free. The environment carries it into spawned multi-GPU
    workers. Values outside ``HEADROOM_RANGE`` are refused.
    """
    text = os.environ.get(HEADROOM_ENVIRONMENT)
    if text is None or text == "":
        return DEFAULT_HEADROOM_FRACTION
    value = float(text)
    if not HEADROOM_RANGE[0] <= value <= HEADROOM_RANGE[1]:
        raise ValueError("%s=%s outside %s" % (HEADROOM_ENVIRONMENT, text, HEADROOM_RANGE))
    return value

POOL_GRANULE = 512                 # CuPy rounds every pool block up to 512 B
DEFAULT_SLACK_BYTES = 64 << 20     # incidental temporaries (scalars, small reductions)
DEFAULT_MODULE_RESERVE = 32 << 20  # compiled CUDA modules of one process (declared)

CATEGORIES = ("static", "basis", "partial_sums", "scratch", "slack", "modules",
              "context")
# Items in these categories are allocated by CuPy's pool and counted in its
# limit.  "modules" and "context" are owned by the driver and only reduce the
# free memory the plan is checked against.
POOLED_CATEGORIES = frozenset(("static", "basis", "partial_sums", "scratch", "slack"))


def pool_bytes(nbytes):
    """Device bytes the CuPy pool actually reserves for one array."""
    nbytes = int(nbytes)
    if nbytes < 0:
        raise ValueError("negative allocation size")
    return -(-nbytes // POOL_GRANULE) * POOL_GRANULE


def human(nbytes):
    """Short binary-prefixed size for messages."""
    value = float(nbytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(value) < 1024.0 or unit == "TiB":
            return ("%.0f %s" % (value, unit)) if unit == "B" else ("%.2f %s" % (value, unit))
        value /= 1024.0


MemoryItem = namedtuple("MemoryItem", "name category nbytes detail")


class MemoryPlan:
    """Declared device allocations of one store or one evaluation phase."""

    def __init__(self, label):
        self.label = label
        self.items = []

    def add(self, name, category, nbytes, detail=""):
        if category not in CATEGORIES:
            raise ValueError("unknown memory category %r" % category)
        nbytes = int(nbytes)
        if nbytes < 0:
            raise ValueError("negative size for %s" % name)
        if any(item.name == name for item in self.items):
            raise ValueError("memory item %r declared twice in %s" % (name, self.label))
        self.items.append(MemoryItem(name, category, nbytes, detail))
        return self

    def add_array(self, name, category, shape, dtype, detail=""):
        """Declare one array by shape and dtype (avoids hand-computed sizes)."""
        import numpy as np
        count = 1
        for extent in (shape if isinstance(shape, tuple) else (shape,)):
            if int(extent) < 0:
                raise ValueError("negative extent for %s" % name)
            count *= int(extent)
        return self.add(name, category, count * np.dtype(dtype).itemsize,
                        detail or "%s %s" % (tuple(shape) if isinstance(shape, tuple)
                                             else (shape,), np.dtype(dtype).name))

    def extend(self, other, prefix=""):
        for item in other.items:
            self.add(prefix + item.name, item.category, item.nbytes, item.detail)
        return self

    @property
    def pooled_bytes(self):
        """Bytes the CuPy pool will hold, including its 512-byte rounding."""
        return sum(pool_bytes(item.nbytes) for item in self.items
                   if item.category in POOLED_CATEGORIES)

    @property
    def driver_bytes(self):
        return sum(item.nbytes for item in self.items
                   if item.category not in POOLED_CATEGORIES)

    @property
    def total_bytes(self):
        return self.pooled_bytes + self.driver_bytes

    def table(self):
        """Human-readable plan for error messages and logs."""
        lines = ["memory plan %s:" % self.label]
        for item in self.items:
            lines.append("  %-28s %-12s %12s  %s" % (item.name, item.category,
                                                    human(item.nbytes), item.detail))
        lines.append("  %-28s %-12s %12s" % ("total (pool-rounded)", "", human(self.total_bytes)))
        return "\n".join(lines)

    def as_dict(self):
        return {"label": self.label, "pooled_bytes": self.pooled_bytes,
                "driver_bytes": self.driver_bytes, "total_bytes": self.total_bytes,
                "items": [{"name": item.name, "category": item.category,
                           "bytes": item.nbytes, "pool_bytes": pool_bytes(item.nbytes)
                           if item.category in POOLED_CATEGORIES else item.nbytes,
                           "detail": item.detail} for item in self.items]}


def check_fits(plan, free_bytes, total_bytes, headroom_fraction, already_reserved=0):
    """Raise ``CapacityError`` unless ``plan`` fits the headroom policy.

    ``already_reserved`` is the planned memory of earlier reservations that
    ``free_bytes`` does not yet reflect (reserved but not allocated).
    """
    if not 0.0 <= headroom_fraction < 1.0:
        raise ValueError("headroom fraction must be in [0, 1)")
    headroom = int(math.ceil(headroom_fraction * total_bytes))
    usable = int(free_bytes) - headroom - int(already_reserved)
    if plan.total_bytes > usable:
        raise CapacityError(
            "%s\ndoes not fit: needs %s; device free %s, headroom %.0f%% of total %s = %s, "
            "already reserved %s, usable %s" %
            (plan.table(), human(plan.total_bytes), human(free_bytes),
             100.0 * headroom_fraction, human(total_bytes), human(headroom),
             human(already_reserved), human(max(usable, 0))))
    return usable - plan.total_bytes


class Reservation:
    """A plan admitted by a ``DeviceBudget``; returned to it on release."""

    def __init__(self, budget, plan):
        self.budget, self.plan, self.active = budget, plan, True

    def release(self):
        if self.active:
            self.active = False
            self.budget._release(self)


class DeviceBudget:
    """Admission control and pool-limit enforcement for one device.

    The budget assumes that it owns the CuPy default pool of its device, so
    the pool must be empty when the budget is created. It measures free and
    total memory once with ``cudaMemGetInfo`` at creation, after the
    context exists, so this process's context is already excluded. Each
    reservation measures again.

    ``capacity=(free, total)`` replaces the measurement in tests. Pool
    enforcement still applies, so a test can exercise refusals without
    filling a 48 GB device.
    """

    def __init__(self, device_id=0, headroom_fraction=None,
                 slack_bytes=DEFAULT_SLACK_BYTES, capacity=None):
        import cupy
        self.cupy = cupy
        self.device_id = int(device_id)
        self.headroom_fraction = float(headroom_fraction if headroom_fraction is not None
                                       else headroom_fraction_default())
        self.slack_bytes = pool_bytes(slack_bytes)
        self.injected_capacity = capacity
        self.reservations = []
        with cupy.cuda.Device(self.device_id):
            self.pool = cupy.get_default_memory_pool()
            if self.pool.total_bytes() != 0:
                raise RuntimeError("device %d CuPy pool already holds %d bytes; a DeviceBudget "
                                   "must own the pool from the start" %
                                   (self.device_id, self.pool.total_bytes()))
            self.initial_free, self.total = self._measure()
        self._apply_limit()

    def _measure(self):
        if self.injected_capacity is not None:
            free, total = self.injected_capacity
            return int(free), int(total)
        free, total = self.cupy.cuda.runtime.memGetInfo()
        return int(free), int(total)

    @property
    def reserved_pooled(self):
        return sum(r.plan.pooled_bytes for r in self.reservations)

    @property
    def reserved_total(self):
        return sum(r.plan.total_bytes for r in self.reservations)

    def _apply_limit(self):
        # CuPy treats a limit of 0 as "unlimited", so the slack is always part
        # of the limit; an empty budget therefore allows only small temporaries.
        with self.cupy.cuda.Device(self.device_id):
            self.pool.set_limit(size=self.reserved_pooled + self.slack_bytes)

    def reserve(self, plan):
        """Admit ``plan`` or raise ``CapacityError`` before any allocation."""
        with self.cupy.cuda.Device(self.device_id):
            # Ledger check against the first measurement: the sum of all plans
            # must fit what was free when the budget took the device.
            check_fits(plan, self.initial_free - self.slack_bytes, self.total,
                       self.headroom_fraction, already_reserved=self.reserved_total)
            # Live check: reservations not yet allocated are not reflected in
            # the current free memory, so add them back explicitly.
            free, _ = self._measure()
            outstanding = max(0, self.reserved_pooled - self.pool.total_bytes())
            check_fits(plan, free - max(0, self.slack_bytes - self.pool.total_bytes()),
                       self.total, self.headroom_fraction, already_reserved=outstanding)
        reservation = Reservation(self, plan)
        self.reservations.append(reservation)
        self._apply_limit()
        return reservation

    def _release(self, reservation):
        self.reservations.remove(reservation)
        with self.cupy.cuda.Device(self.device_id):
            # Return cached blocks to the driver so that the measured device
            # memory, not only the pool's bookkeeping, goes back down.
            self.pool.free_all_blocks()
        self._apply_limit()

    def allocation_guard(self, plan):
        """Context manager converting CuPy OOM into ``CapacityError`` with the plan."""
        return _AllocationGuard(self, plan)

    def state(self):
        """Measured and planned memory for manifests and leak checks."""
        with self.cupy.cuda.Device(self.device_id):
            free, total = self.cupy.cuda.runtime.memGetInfo()
            return {"device_id": self.device_id, "total_bytes": int(total),
                    "free_bytes": int(free), "initial_free_bytes": self.initial_free,
                    "injected_capacity": self.injected_capacity is not None,
                    "headroom_fraction": self.headroom_fraction,
                    "slack_bytes": self.slack_bytes,
                    "reserved_total_bytes": self.reserved_total,
                    "pool_limit_bytes": int(self.pool.get_limit()),
                    "pool_used_bytes": int(self.pool.used_bytes()),
                    "pool_total_bytes": int(self.pool.total_bytes())}

    def close(self):
        """Drop every reservation and remove the pool limit."""
        for reservation in list(self.reservations):
            reservation.release()
        with self.cupy.cuda.Device(self.device_id):
            self.pool.free_all_blocks()
            self.pool.set_limit(size=0)


class _AllocationGuard:
    def __init__(self, budget, plan):
        self.budget, self.plan = budget, plan

    def __enter__(self):
        return self

    def __exit__(self, kind, value, traceback):
        if kind is not None and issubclass(kind, self.budget.cupy.cuda.memory.OutOfMemoryError):
            raise CapacityError("device %d allocation outside or beyond the plan: %s\n%s" %
                                (self.budget.device_id, value, self.plan.table())) from value
        return False


# ----------------------------------------------------------------------- tiling
TilePlan = namedtuple("TilePlan", "point_tiles column_groups records_per_block "
                                  "partial_sum_bytes oversized_points")


def point_tiles(counts, max_records):
    """Contiguous point ranges ``[p0, p1)`` with at most ``max_records`` records.

    A point above the bound forms its own tile and is listed separately, so
    the caller knows that its reduction is split across blocks.
    """
    if max_records <= 0:
        raise ValueError("max_records must be positive")
    tiles, oversized = [], []
    start, filled = 0, 0
    for index, count in enumerate(counts):
        count = int(count)
        if count > max_records:
            if index > start:
                tiles.append((start, index))
            tiles.append((index, index + 1))
            oversized.append(index)
            start, filled = index + 1, 0
            continue
        if filled + count > max_records and index > start:
            tiles.append((start, index))
            start, filled = index, 0
        filled += count
    if start < len(counts):
        tiles.append((start, len(counts)))
    return tiles, oversized


def column_groups(ncols, max_columns):
    """Contiguous column ranges ``[c0, c1)`` of at most ``max_columns`` each."""
    if ncols < 0 or max_columns <= 0:
        raise ValueError("invalid column grouping")
    return [(start, min(start + max_columns, ncols)) for start in range(0, ncols, max_columns)]


def partial_sum_bytes(counts, records_per_block, ncols, gauges=2):
    """FP64 bytes of per-block partial sums for one tile and one column group.

    Each point is reduced by ``ceil(count / records_per_block)`` blocks (at
    least one). Each block writes ``gauges * ncols`` partial sums, which a
    second fixed-order pass combines. That two-pass shape makes the
    reduction independent of scheduling.
    """
    blocks = sum(max(1, -(-int(count) // records_per_block)) for count in counts)
    return blocks * gauges * max(ncols, 1) * 8


def plan_tiles(counts, ncols, scratch_budget_bytes, records_per_block=1 << 16,
               max_tile_records=1 << 26, gauges=2):
    """Choose observation tiles and column groups within a scratch budget.

    Starts from all columns in one group and the given tile record bound. It
    halves the column-group width until the largest tile's partial sums fit,
    and halves the tile bound when even single columns do not fit. Only
    partial sums are budgeted. The C06 kernels are fused and keep no
    per-record-per-column intermediates, so no other scratch scales with
    records times columns.
    """
    counts = [int(c) for c in counts]
    if scratch_budget_bytes <= 0:
        raise ValueError("scratch budget must be positive")
    tile_bound = max_tile_records
    while True:
        tiles, oversized = point_tiles(counts, tile_bound)
        width = max(ncols, 1)
        while True:
            need = max(partial_sum_bytes(counts[p0:p1], records_per_block, width, gauges)
                       for p0, p1 in tiles) if tiles else 0
            if need <= scratch_budget_bytes or width == 1:
                break
            width = max(1, width // 2)
        if need <= scratch_budget_bytes:
            return TilePlan(tiles, column_groups(ncols, width) if ncols else [],
                            records_per_block, need, oversized)
        if tile_bound <= records_per_block:
            raise CapacityError("scratch budget %s cannot hold the partial sums of one "
                                "block-sized tile and one column" % human(scratch_budget_bytes))
        tile_bound = max(records_per_block, tile_bound // 2)
