"""Vectorized fitting-base validation equals the reference (2026-10-21, SHAP8 12M).

The audit returned for a valid base is identical, and every corruption the
reference rejects is rejected with the identical first message, including
files with several faults in different records.
"""

import numpy as np
import pytest

import synthetic_run
from cops_cuda import base_validation, reference

reference.bind()
io = reference.load("cops_io")


def outcome(fn, *args):
    try:
        return "ok", fn(*args)
    except io.ProductError as exc:
        return "error", str(exc)


def both(run):
    return (outcome(io.validate_fit_base, run),
            outcome(base_validation.validate_fit_base, io, run))


def test_valid_base_identical(synthetic):
    ref, fast = both(synthetic["run"])
    assert ref[0] == "ok" and ref == fast
    assert io.validate_kernel.__module__ == io.__name__          # workers restored


KERNEL_FAULTS = {
    "range": [(7, "face", 10 ** 6)],
    "negative": [(0, "face", -1)],
    "duplicate": [(9, "face", None)],                 # None: copy record 2's face
    "nan": [(4, "KV", np.nan)],
    "inf": [(4, "cos", -np.inf)],
    "shadow": [(6, "shadow", 2)],
    "nan_before_duplicate": [(9, "face", None), (5, "Kn", np.nan)],
    "duplicate_before_shadow": [(3, "face", None), (8, "shadow", -1)],
    "same_record_range_and_nan": [(5, "face", 99999), (5, "KPhi", np.nan)],
    "same_record_duplicate_and_shadow": [(5, "face", None), (5, "shadow", 7)],
}


@pytest.mark.parametrize("name", sorted(KERNEL_FAULTS))
def test_kernel_faults_identical(synthetic, name):
    run, point = synthetic["run"], synthetic["point_ids"][0]
    geometry, records = synthetic["kernels"][point]
    records = records.copy()
    for r, field, value in KERNEL_FAULTS[name]:
        records[field][r] = records["face"][2] if value is None else value
    synthetic_run.write_kernel(synthetic_run.kernel_path(run, point), 1, point, geometry, records)
    synthetic_run.write_manifest(run, synthetic["point_ids"])
    ref, fast = both(run)
    assert ref[0] == "error" and ref == fast


@pytest.mark.parametrize("name", ["header_nan", "version", "point", "count", "truncated"])
def test_kernel_header_faults_identical(synthetic, name):
    run, point = synthetic["run"], synthetic["point_ids"][1]
    geometry, records = synthetic["kernels"][point]
    geometry, version, count = geometry.copy(), 1, None
    if name == "header_nan":
        geometry[3] = np.nan
    elif name == "version":
        version = 2
    elif name == "count":
        count = synthetic["nfaces"] + 1
    elif name == "truncated":
        records = records[:-1]
        count = len(records) + 1
    synthetic_run.write_kernel(synthetic_run.kernel_path(run, point), version,
                               point + 1 if name == "point" else point, geometry, records, count)
    synthetic_run.write_manifest(run, synthetic["point_ids"])
    ref, fast = both(run)
    assert ref[0] == "error" and ref == fast


@pytest.mark.parametrize("name", ["nan", "area_zero", "area_negative", "nan_after_area"])
def test_element_faults_identical(synthetic, name):
    run = synthetic["run"]
    table = synthetic["elements"].copy()
    if name == "nan":
        table[17, 20] = np.nan
    elif name == "area_zero":
        table[30, 6] = 0.0
    elif name == "area_negative":
        table[31, 6] = -1.0
    else:
        table[40, 6] = 0.0
        table[50, 2] = np.inf
    synthetic_run.write_elements(run / "elements.bin", table)
    ref, fast = both(run)
    assert ref[0] == "error" and ref == fast


def test_element_chunks(synthetic, monkeypatch):
    """Faults found across read chunks report the global face index."""
    monkeypatch.setattr(base_validation, "ELEMENT_CHUNK", 7)
    test_valid_base_identical(synthetic)
    table = synthetic["elements"].copy()
    table[100, 6] = -2.0
    synthetic_run.write_elements(synthetic["run"] / "elements.bin", table)
    ref, fast = both(synthetic["run"])
    assert ref[0] == "error" and ref == fast


def test_install_replaces_only_fit_cops_entry():
    fit_cops = reference.load("fit_cops")
    original = fit_cops.validate_fit_base
    try:
        info = base_validation.install(fit_cops)
        assert info["policy"] == base_validation.POLICY
        assert io.validate_fit_base is original                    # reference module untouched
    finally:
        fit_cops.validate_fit_base = original


def test_serial_and_prefetch_identical(synthetic):
    run = synthetic["run"]
    assert (base_validation.validate_fit_base(io, run, threads=1) ==
            base_validation.validate_fit_base(io, run, threads=8) == io.validate_fit_base(run))
