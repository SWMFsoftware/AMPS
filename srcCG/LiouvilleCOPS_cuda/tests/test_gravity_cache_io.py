"""C04 GPU-cache writer: CPU-loader compatibility without CPU code changes (CPU only)."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT / "tools"))
import c00_capture_reference as c00  # noqa: E402
from test_derive_surface import write_bdf  # noqa: E402

from cops_cuda import gravity_body, gravity_cache_io, reference  # noqa: E402


def make_body(tmp_path, level=2, name="body.bdf"):
    v, f = c00.icosphere(1000.0, level)
    path = tmp_path / name
    write_bdf(path, v, f, extra_cards=False)
    return gravity_body.GravityBody.from_bdf(path, 533.0)


def test_gpu_cache_passes_the_cpu_reuse_checks(synthetic, tmp_path):
    gravity_cache = reference.load("gravity_cache")
    topology = reference.load("gravity_topology")
    fit_cops = reference.load("fit_cops")
    run = synthetic["run"]
    body = make_body(tmp_path)
    el = fit_cops.read_elements(str(run / "elements.bin"))
    g = np.random.default_rng(1).normal(0.0, 1e-4, el["centroid"].shape)
    omega = np.array([1e-8, -2e-8, 1.4e-4])
    out = tmp_path / "gpu_cache.npz"
    gravity_cache_io.write_gpu_cache(out, run, body, g, omega, 0.1, (0.005, 0.003, 240), 0.25,
                                     {"backend": "test"}, [], 1.0)
    # The builder's reuse path: load_cache against its own definition ...
    definition = gravity_cache.make_definition(el["centroid"], el["normal"], el["area"], body.source,
                                               len(body.vertices), len(body.faces), 533.0, 0.1)
    cached = gravity_cache.load_cache(out, el["centroid"], el["normal"], el["area"], definition, omega)
    assert cached["self_gravity"].tobytes() == g.tobytes()
    calc = cached["metadata"]["calculation"]
    # ... then exactly the builder's metadata checks.
    V, F = reference.load("build_physical_activity_templates").read_nastran_surface(body.source["path"])
    assert calc["orientation_policy"] == topology.ORIENTATION_POLICY
    assert calc["topology"] == json.loads(json.dumps(topology.require_outward_closed_body(V, F)))
    assert float(calc["offset_halving_slope_difference_deg_max"]) <= 0.25
    assert calc["integrity_check"] == "DISABLE after topology audit"
    # The selector's provenance gate reads the cache definition.
    assert definition["algorithm"] == gravity_cache.ALGORITHM
    assert definition["orientation_policy"] == topology.ORIENTATION_POLICY
    # No CPU package identity is fabricated.
    assert not {"package_version", "package_binary", "package_commit"} & set(calc)
    assert calc["backend"].startswith("cuda_polyhedral_gravity_v1")
    # A cache is exclusive.
    with pytest.raises(ValueError):
        gravity_cache_io.write_gpu_cache(out, run, body, g, omega, 0.1, (0.005, 0.003, 240), 0.25,
                                         {}, [], 1.0)


def test_cache_for_another_body_is_refused(synthetic, tmp_path):
    gravity_cache = reference.load("gravity_cache")
    fit_cops = reference.load("fit_cops")
    run = synthetic["run"]
    body = make_body(tmp_path)
    other = make_body(tmp_path, level=1, name="other.bdf")
    el = fit_cops.read_elements(str(run / "elements.bin"))
    g = -1e-4 * el["normal"]                # nonzero (the CPU writer refuses zero gravity)
    out = tmp_path / "gpu_cache.npz"
    gravity_cache_io.write_gpu_cache(out, run, body, g, np.zeros(3), 0.1, (0.0, 0.0, 240), 0.25,
                                     {}, [], 1.0)
    wrong = gravity_cache.make_definition(el["centroid"], el["normal"], el["area"], other.source,
                                          len(other.vertices), len(other.faces), 533.0, 0.1)
    with pytest.raises(ValueError):
        gravity_cache.load_cache(out, el["centroid"], el["normal"], el["area"], wrong, np.zeros(3))


def test_offset_audit_above_tolerance_is_refused(synthetic, tmp_path):
    fit_cops = reference.load("fit_cops")
    run = synthetic["run"]
    body = make_body(tmp_path)
    g = np.zeros(fit_cops.read_elements(str(run / "elements.bin"))["centroid"].shape)
    with pytest.raises(ValueError, match="offset-halving"):
        gravity_cache_io.write_gpu_cache(tmp_path / "x.npz", run, body, g, np.zeros(3), 0.1,
                                         (0.3, 0.2, 240), 0.25, {}, [], 1.0)
