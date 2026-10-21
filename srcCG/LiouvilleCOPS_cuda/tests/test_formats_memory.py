"""C02 host-side contracts: binary formats and memory/tiling arithmetic (no GPU)."""

import numpy as np
import pytest

from cops_cuda import formats, memory, precision, reference
from cops_cuda.errors import CacheFormatError, CapacityError, IdentityError
import synthetic_run


# ---------------------------------------------------------------- formats
def test_layouts_equal_reference_definitions():
    cops_io = reference.load("cops_io")
    fit_cops = reference.load("fit_cops")
    assert formats.KERNEL_HEADER_DTYPE.itemsize == cops_io.KERNEL1_HEADER.size
    assert formats.KERNEL_RECORD_DTYPE.itemsize == cops_io.KERNEL1_RECORD.size
    assert formats.KERNEL_RECORD_DTYPE == fit_cops.REC
    assert formats.ELEMENT_HEADER_DTYPE.itemsize == cops_io.ELEMENT_HEADER.size


def test_elements_and_kernels_parse_like_reference(synthetic):
    fit_cops = reference.load("fit_cops")
    run = synthetic["run"]
    fields, identity = formats.load_elements(run / "elements.bin")
    ref = fit_cops.read_elements(str(run / "elements.bin"))
    for name in ("centroid", "normal", "area", "Y", "AH2O", "ACO2"):
        assert np.array_equal(fields[name], ref[name]), name
    assert identity["sha256"] == synthetic_run.sha256(run / "elements.bin")
    for point in synthetic["point_ids"]:
        path = synthetic_run.kernel_path(run, point)
        geometry, records, digest = formats.read_kernel(path, point, synthetic["nfaces"])
        ip, header, ref_records = fit_cops.read_kernel(str(path))
        assert ip == point and np.array_equal(records, ref_records)
        assert geometry[1] == header["VdotB"] and geometry[8] == header["rh"]
        assert digest == synthetic_run.sha256(path)


def _damage(path, offset, data):
    raw = bytearray(path.read_bytes())
    raw[offset:offset + len(data)] = data
    path.write_bytes(bytes(raw))


@pytest.mark.parametrize("damage, error, text", [
    ("truncate", CacheFormatError, "size mismatch"),
    ("append", CacheFormatError, "size mismatch"),
    ("version", CacheFormatError, "kernel version"),
    ("point", IdentityError, "stores iPoint"),
    ("face_range", CacheFormatError, "out-of-range face"),
    ("duplicate", CacheFormatError, "duplicate face"),
    ("nan", CacheFormatError, "non-finite"),
    ("shadow", CacheFormatError, "invalid shadow"),
])
def test_each_kernel_defect_is_named(synthetic, damage, error, text):
    point, nfaces = synthetic["point_ids"][0], synthetic["nfaces"]
    path = synthetic_run.kernel_path(synthetic["run"], point)
    raw = path.read_bytes()
    if damage == "truncate":
        path.write_bytes(raw[:-1])
    elif damage == "append":
        path.write_bytes(raw + b"\0")
    elif damage == "version":
        _damage(path, 0, np.int32(2).tobytes())
    elif damage == "point":
        _damage(path, 4, np.int32(point + 1).tobytes())
    elif damage == "face_range":
        _damage(path, 92, np.int32(nfaces).tobytes())
    elif damage == "duplicate":
        face0 = raw[92:96]
        _damage(path, 92 + 40, face0)
    elif damage == "nan":
        _damage(path, 92 + 4, np.float64(np.nan).tobytes())
    elif damage == "shadow":
        _damage(path, 92 + 36, np.int32(2).tobytes())
    with pytest.raises(error, match=text):
        formats.read_kernel(path, point, nfaces)


def test_element_defects(synthetic):
    path = synthetic["run"] / "elements.bin"
    raw = path.read_bytes()
    path.write_bytes(raw + b"\0")
    with pytest.raises(CacheFormatError, match="size mismatch"):
        formats.load_elements(path)
    area = bytearray(raw)
    area[12 + 6 * 8:12 + 7 * 8] = np.float64(0.0).tobytes()
    path.write_bytes(bytes(area))
    with pytest.raises(CacheFormatError, match="non-positive area"):
        formats.load_elements(path)


def test_set_digest_matches_manifest_definition(synthetic):
    import json
    run = synthetic["run"]
    entries = [("point=%04d.bin" % p, synthetic_run.sha256(synthetic_run.kernel_path(run, p)))
               for p in synthetic["point_ids"]]
    manifest = json.loads((run / "run_manifest.json").read_text())
    assert formats.kernel_set_digest(entries) == manifest["outputs"]["kernels"]["set_sha256"]


# -------------------------------------------------------------- precision
def test_precision_policy_refuses_fp32_and_fast_math():
    with pytest.raises(TypeError):
        precision.require_device_dtype("kn", np.float32)
    with pytest.raises(ValueError):
        precision.require_device_dtype("unknown_field", np.float64)
    with pytest.raises(ValueError):
        precision.check_nvrtc_options(("--fmad=false", "--use_fast_math"))
    with pytest.raises(ValueError):
        precision.check_nvrtc_options(("--std=c++14",))
    assert precision.check_nvrtc_options(precision.REFERENCE_NVRTC_OPTIONS)


# ----------------------------------------------------------------- memory
def test_plan_rounding_and_categories():
    plan = memory.MemoryPlan("t")
    plan.add("a", "static", 1)
    plan.add_array("b", "static", (3, 5), np.float64)        # 120 B -> 512 B
    plan.add("ctx", "context", 1000)
    assert plan.pooled_bytes == 1024 and plan.driver_bytes == 1000
    assert plan.total_bytes == 2024
    with pytest.raises(ValueError):
        plan.add("a", "static", 8)
    with pytest.raises(ValueError):
        plan.add("c", "bogus", 8)


def test_check_fits_applies_headroom():
    plan = memory.MemoryPlan("t").add("x", "static", 700 * 512)
    total = 1000 * 512
    memory.check_fits(plan, free_bytes=total, total_bytes=total, headroom_fraction=0.2)
    with pytest.raises(CapacityError, match="does not fit"):
        memory.check_fits(plan, free_bytes=total, total_bytes=total, headroom_fraction=0.31)
    with pytest.raises(CapacityError):
        memory.check_fits(plan, free_bytes=total, total_bytes=total, headroom_fraction=0.2,
                          already_reserved=200 * 512)


def test_point_tiles_cover_in_order_and_respect_bound():
    counts = [5, 0, 7, 30, 2, 2, 9, 1]
    tiles, oversized = memory.point_tiles(counts, 10)
    flat = [p for p0, p1 in tiles for p in range(p0, p1)]
    assert flat == list(range(len(counts)))
    assert oversized == [3]
    for p0, p1 in tiles:
        assert sum(counts[p0:p1]) <= 10 or (p1 - p0 == 1 and p0 in oversized)
    assert memory.column_groups(7, 3) == [(0, 3), (3, 6), (6, 7)]


def test_plan_tiles_shrinks_columns_then_tiles():
    counts = [100000] * 50
    plan = memory.plan_tiles(counts, ncols=32, scratch_budget_bytes=1 << 20,
                             records_per_block=1 << 12)
    assert plan.partial_sum_bytes <= 1 << 20
    assert [c for g in plan.column_groups for c in range(*g)] == list(range(32))
    assert [p for t in plan.point_tiles for p in range(*t)] == list(range(50))
    with pytest.raises(CapacityError):
        memory.plan_tiles(counts, ncols=1, scratch_budget_bytes=8, records_per_block=1 << 12)


def test_headroom_run_setting(monkeypatch):
    """gpu_fit.py --headroom-fraction travels through the environment into every budget."""
    monkeypatch.delenv(memory.HEADROOM_ENVIRONMENT, raising=False)
    assert memory.headroom_fraction_default() == memory.DEFAULT_HEADROOM_FRACTION == 0.20
    monkeypatch.setenv(memory.HEADROOM_ENVIRONMENT, "0.15")
    assert memory.headroom_fraction_default() == 0.15
    for bad in ("0.01", "0.7"):
        monkeypatch.setenv(memory.HEADROOM_ENVIRONMENT, bad)
        with pytest.raises(ValueError):
            memory.headroom_fraction_default()
