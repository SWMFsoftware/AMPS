"""C02 gate on synthetic runs: layout, identities, release cycles, refusals (GPU).

Gate wording (ROADMAP C02): "load/release cycles preserve data and bounded
memory; malformed caches, wrong mesh order, and insufficient capacity are
detected without silent changes".  Each test below covers one clause; every
refusal test also checks that nothing stays allocated afterwards.
"""

import json

import numpy as np
import pytest

from cops_cuda import device_store, memory, reference
from cops_cuda.errors import CacheFormatError, CapacityError, IdentityError
import synthetic_run


def host_expectation(synthetic, point_ids):
    """Concatenated records/offsets in the requested order (independent path)."""
    records = [synthetic["kernels"][p][1] for p in point_ids]
    offsets = np.concatenate([[0], np.cumsum([len(r) for r in records])]).astype(np.int64)
    joined = np.concatenate(records) if records else np.zeros(0, dtype=synthetic_run.RECORD)
    geometry = np.array([synthetic["kernels"][p][0] for p in point_ids])
    return joined, offsets, geometry


def assert_layout(ds, synthetic, point_ids):
    joined, offsets, geometry = host_expectation(synthetic, point_ids)
    a = ds.arrays
    assert np.array_equal(a["kernel.face"].get(), joined["face"])
    for device, field in (("kn", "Kn"), ("kphi", "KPhi"), ("kv", "KV"), ("cos", "cos")):
        got = a["kernel." + device].get()
        # bitwise, including signed zeros
        assert got.tobytes() == np.ascontiguousarray(joined[field]).tobytes(), device
    assert np.array_equal(a["kernel.shadow"].get(), joined["shadow"].astype(np.uint8))
    assert np.array_equal(a["kernel.offsets"].get(), offsets)
    assert np.array_equal(a["kernel.geometry"].get(), geometry.T)
    assert a["kernel.offsets"].dtype == np.int64 and a["kernel.face"].dtype == np.int32


def pool_is_empty(cp):
    pool = cp.get_default_memory_pool()
    return pool.used_bytes() == 0 and pool.total_bytes() == 0


def test_full_load_matches_files_and_manifest(clean_pool, synthetic):
    cp = clean_pool
    with device_store.DeviceDataset.load(synthetic["run"], mesh_fields="all") as ds:
        assert_layout(ds, synthetic, synthetic["point_ids"])
        el = synthetic["elements"]
        assert np.array_equal(ds.arrays["mesh.Y"].get(), el[:, 7:32].T)
        assert np.array_equal(ds.arrays["mesh.centroid"].get(), el[:, 0:3].T)
        assert np.array_equal(ds.arrays["mesh.ACO2"].get(), el[:, 33])
        ident = ds.identity()
        assert ident["kernels"]["verified_against_run_manifest"]
        assert ident["device_audit"] == dict.fromkeys(device_store.AUDIT_COUNTERS, 0)
        ds.verify_roundtrip()
        # The pool holds exactly the planned arrays (rounded) plus the audit's
        # small counter array, which was already released.
        assert cp.get_default_memory_pool().used_bytes() <= ds.plan.pooled_bytes
    assert pool_is_empty(cp)


def test_subset_in_caller_order_with_tiny_chunks(clean_pool, synthetic):
    """Records crossing chunk and file boundaries keep their exact order."""
    order = [synthetic["point_ids"][i] for i in (4, 0, 2, 5)]
    with device_store.DeviceDataset.load(synthetic["run"], point_ids=order,
                                         chunk_records=7) as ds:
        assert_layout(ds, synthetic, order)
        assert not ds.identity()["kernels"]["verified_against_run_manifest"]
        ds.verify_roundtrip()
    assert pool_is_empty(clean_pool)


def test_load_release_cycles_are_bitwise_and_bounded(clean_pool, synthetic):
    cp = clean_pool
    digests, free_after = [], []
    for _ in range(3):
        with device_store.DeviceDataset.load(synthetic["run"], chunk_records=50) as ds:
            digests.append(ds.verify_roundtrip())
        assert pool_is_empty(cp)
        assert cp.get_default_memory_pool().get_limit() == 0   # private budget closed
        free_after.append(cp.cuda.runtime.memGetInfo()[0])
    assert digests[0] == digests[1] == digests[2]
    # Device free memory returns to the same level (no growth across cycles).
    assert max(free_after) - min(free_after) <= 2 << 20


def test_expected_file_digests(clean_pool, synthetic):
    run, points = synthetic["run"], synthetic["point_ids"]
    good = {"point=%04d.bin" % p: synthetic_run.sha256(synthetic_run.kernel_path(run, p))
            for p in points}
    with device_store.DeviceDataset.load(run, expected_kernel_digests=good) as ds:
        assert ds.identity()["kernels"]["files_verified_against_expected"] == len(points)
    bad = dict(good)
    bad["point=%04d.bin" % points[1]] = "0" * 64
    with pytest.raises(IdentityError, match="differs from the expected"):
        device_store.DeviceDataset.load(run, expected_kernel_digests=bad)
    assert pool_is_empty(clean_pool)


# ------------------------------------------------------------ malformed caches
@pytest.mark.parametrize("damage, error, text", [
    ("truncate", CacheFormatError, "size mismatch"),
    ("missing", CacheFormatError, "missing"),
    ("geometry", IdentityError, "differs between table and kernel"),
    ("duplicate_across_chunks", CacheFormatError, "duplicate face"),
    ("nan_late", CacheFormatError, "non-finite"),
])
def test_malformed_cache_refused_without_leak(clean_pool, synthetic, damage, error, text):
    run, points = synthetic["run"], synthetic["point_ids"]
    path = synthetic_run.kernel_path(run, points[3])
    raw = bytearray(path.read_bytes())
    if damage == "truncate":
        path.write_bytes(bytes(raw[:-40]))
    elif damage == "missing":
        path.unlink()
    elif damage == "geometry":
        raw[8:16] = np.float64(1.0).tobytes()          # et of the header
        path.write_bytes(bytes(raw))
    elif damage == "duplicate_across_chunks":
        # Copy face of record 0 into record 5; with chunk_records=3 the two
        # records are validated in different chunks.
        raw[92 + 5 * 40:92 + 5 * 40 + 4] = raw[92:96]
        path.write_bytes(bytes(raw))
    elif damage == "nan_late":
        last = len(raw) - 40
        raw[last + 12:last + 20] = np.float64(np.inf).tobytes()   # KPhi of the last record
        path.write_bytes(bytes(raw))
    with pytest.raises(error, match=text):
        device_store.DeviceDataset.load(run, chunk_records=3)
    assert pool_is_empty(clean_pool)
    assert clean_pool.get_default_memory_pool().get_limit() == 0


def test_manifest_set_digest_mismatch(clean_pool, synthetic):
    run = synthetic["run"]
    manifest = json.loads((run / "run_manifest.json").read_text())
    manifest["outputs"]["kernels"]["set_sha256"] = "f" * 64
    (run / "run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IdentityError, match="kernel set sha256"):
        device_store.DeviceDataset.load(run)
    assert pool_is_empty(clean_pool)


# ------------------------------------------------------------ wrong mesh order
def test_reordered_mesh_is_refused(clean_pool, synthetic):
    run = synthetic["run"]
    original_identity = device_store.RunInputs(run).mesh_identity
    table = synthetic["elements"].copy()
    table[[0, 1]] = table[[1, 0]]                     # swap two facets
    (run / "elements.bin").unlink()
    synthetic_run.write_elements(run / "elements.bin", table)
    # Manifest still describes the original mesh -> refused.
    with pytest.raises(IdentityError, match="does not match run_manifest"):
        device_store.DeviceDataset.load(run)
    # Even a manifest rewritten for the new file is refused against the
    # identity a template bundle / gravity cache carries.
    synthetic_run.write_manifest(run, synthetic["point_ids"])
    with pytest.raises(IdentityError, match="mesh identity"):
        device_store.DeviceDataset.load(run, expected_mesh_identity=original_identity)
    assert pool_is_empty(clean_pool)


def test_basis_for_another_mesh_is_refused(clean_pool, synthetic):
    run, nfaces = synthetic["run"], synthetic["nfaces"]
    inputs = device_store.RunInputs(run)
    matrix = np.ones((nfaces, 3))
    with pytest.raises(IdentityError, match="built for mesh"):
        device_store.DeviceDataset.load(run, enhancement_basis={
            "matrix": matrix, "mesh_identity": "0" * 64})
    with pytest.raises(TypeError):
        device_store.DeviceDataset.load(run, enhancement_basis={
            "matrix": matrix.astype(np.float32), "mesh_identity": inputs.mesh_identity})
    with device_store.DeviceDataset.load(run, enhancement_basis={
            "matrix": matrix, "mesh_identity": inputs.mesh_identity, "scale": 2.0}) as ds:
        assert np.array_equal(ds.arrays["basis.enhancement_basis"].get(), matrix.T)
        assert ds.identity()["enhancement_basis"]["columns"] == 3
    assert pool_is_empty(clean_pool)


def test_mesh_identity_uses_reference_definition(synthetic):
    templates = reference.load("activity_templates")
    el = synthetic["elements"]
    assert device_store.RunInputs(synthetic["run"]).mesh_identity == \
        templates.mesh_identity(el[:, 0:3], el[:, 6])


# --------------------------------------------------------- insufficient capacity
def test_capacity_refused_before_any_allocation(clean_pool, synthetic):
    cp = clean_pool
    free, total = cp.cuda.runtime.memGetInfo()
    # A device "with" 1 MiB free: the synthetic dataset needs less, but the
    # 20% headroom of the real total plus slack does not leave 1 MiB.
    budget = memory.DeviceBudget(0, capacity=(1 << 20, total))
    try:
        with pytest.raises(CapacityError, match="does not fit"):
            device_store.DeviceDataset.load(synthetic["run"], budget=budget)
        assert cp.get_default_memory_pool().total_bytes() == 0
        assert not budget.reservations
    finally:
        budget.close()


def test_allocation_beyond_plan_is_refused(clean_pool):
    cp = clean_pool
    budget = memory.DeviceBudget(0, slack_bytes=1 << 20)
    try:
        plan = memory.MemoryPlan("small").add_array("x", "static", (1024,), np.float64)
        reservation = budget.reserve(plan)
        with budget.allocation_guard(plan):
            x = cp.empty(1024, dtype=np.float64)          # planned
        with pytest.raises(CapacityError, match="beyond the plan"):
            with budget.allocation_guard(plan):
                cp.empty(4 << 20, dtype=np.uint8)          # 4 MiB, not planned
        del x
        reservation.release()
        assert cp.get_default_memory_pool().total_bytes() == 0
    finally:
        budget.close()


# --------------------------------------------------------------- device audit
@pytest.mark.parametrize("field, index, value, counter", [
    ("kernel.face", 0, -1, "face_out_of_range"),
    ("kernel.kv", 3, np.nan, "nonfinite_value"),
    ("kernel.shadow", 2, 7, "invalid_shadow"),
    ("kernel.offsets", 1, 10 ** 6, "decreasing_offsets"),
])
def test_device_audit_detects_corruption(clean_pool, synthetic, field, index, value, counter):
    with device_store.DeviceDataset.load(synthetic["run"]) as ds:
        ds.arrays[field][index] = value
        with pytest.raises(CacheFormatError, match="device audit failed"):
            ds.audit()
        assert ds.audit_result[counter] >= 1
    assert pool_is_empty(clean_pool)
