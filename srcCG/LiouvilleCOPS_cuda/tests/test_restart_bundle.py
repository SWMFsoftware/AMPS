"""C05 restart bundle (cops_cuda/restart_bundle.py) on a synthetic run (CPU only).

The bundle must be accepted by the *unchanged* fitter's loader and selector
(including the gravity provenance gate for the cliff columns), and its
columns must be exactly the inputs: the visibility masks and ``slope >= th``
from the admitted cache, with band facets taken from the resolver.
"""

import numpy as np
import pytest

from cops_cuda import gravity_cache_io, reference, restart_bundle
from test_gravity_cache_io import make_body


def vis_inputs(nf, rng):
    ng, rg = rng.random(nf) < 0.3, rng.random(nf) < 0.4
    vis = {"late_ng_not_rg": ng & ~rg, "late_rg_not_ng": rg & ~ng, "ng_seen": ng, "rg_seen": rg}
    meta = {"interval_utc": ["2016-09-30T10:33:10", "2016-09-30T10:39:26"],
            "selected_point_count": 2, "selected_point_ids": [1, 2],
            "definitions": {k: "test " + k for k in vis},
            "face_counts": {k: int(v.sum()) for k, v in vis.items()},
            "kernel_set_sha256": "0" * 64, "backend": "test"}
    return vis, meta


def admitted_cache(synthetic, tmp_path, rng):
    gravity_cache = reference.load("gravity_cache")
    fit_cops = reference.load("fit_cops")
    run = synthetic["run"]
    body = make_body(tmp_path)
    el = fit_cops.read_elements(str(run / "elements.bin"))
    # Inward gravity with tangential parts: slopes spread over 0..90 deg.
    g = -1e-4 * el["normal"] + rng.normal(0.0, 1.2e-4, el["normal"].shape)
    omega = np.array([1e-8, -2e-8, 1.4e-4])
    path = tmp_path / "gpu_cache.npz"
    gravity_cache_io.write_gpu_cache(path, run, body, g, omega, 0.1, (0.005, 0.003, 240), 0.25,
                                     {"backend": "test"}, [], 1.0)
    definition = gravity_cache.make_definition(el["centroid"], el["normal"], el["area"], body.source,
                                               len(body.vertices), len(body.faces), 533.0, 0.1)
    cached = gravity_cache.load_cache(path, el["centroid"], el["normal"], el["area"], definition, omega)
    return run, body, path, definition, cached, omega, el


def test_bundle_is_accepted_by_the_unchanged_fitter(synthetic, tmp_path):
    templates = reference.load("activity_templates")
    rng = np.random.default_rng(5)
    run, body, path, definition, cached, omega, el = admitted_cache(synthetic, tmp_path, rng)
    vis, meta = vis_inputs(len(el["area"]), rng)
    out = tmp_path / "restart" / "bundle.npz"
    calls = []
    payload = restart_bundle.build_bundle(
        out, run, body, cached, definition, path, (omega, ["e0", "e1", "e2"], 0.0, ["spice"]),
        vis, meta, (45.0, 60.0), lambda faces: calls.append(faces) or cached["slope"][faces],
        write_map=False)
    names = ["late_ng_not_rg", "late_rg_not_ng", "cliff_ge_45deg", "cliff_ge_60deg"]
    assert payload["names"] == names
    bundle = templates.load_bundle(out, el["centroid"], el["area"])
    values = templates.select_templates(bundle, names)        # includes the provenance gate
    assert np.array_equal(values[:, 0] > 0, vis["late_ng_not_rg"])
    assert np.array_equal(values[:, 1] > 0, vis["late_rg_not_ng"])
    slope = cached["slope"]
    assert np.array_equal(values[:, 2] > 0, slope >= 45.0)
    assert np.array_equal(values[:, 3] > 0, slope >= 60.0)
    assert 0 < np.count_nonzero(values[:, 3]) < np.count_nonzero(values[:, 2]) < len(slope)
    gravity = bundle["metadata"]["gravity"]
    assert gravity["integrity_check"] == "DISABLE after topology audit"
    assert gravity["backend"].startswith("cuda_polyhedral_gravity_v1")
    assert gravity["cache"]["definition"] == definition
    assert gravity["density_kg_m3"] == 533.0 and gravity["evaluation_offset_m"] == 0.1
    assert (out.with_suffix(".json")).is_file()
    # Exclusive output.
    with pytest.raises(ValueError):
        restart_bundle.build_bundle(out, run, body, cached, definition, path,
                                    (omega, [], 0.0, []), vis, meta, write_map=False)


def test_cache_for_another_topology_is_refused(synthetic, tmp_path):
    rng = np.random.default_rng(6)
    run, body, path, definition, cached, omega, el = admitted_cache(synthetic, tmp_path, rng)
    vis, meta = vis_inputs(len(el["area"]), rng)
    other = make_body(tmp_path, level=1, name="other.bdf")
    with pytest.raises(ValueError, match="topology"):
        restart_bundle.build_bundle(tmp_path / "x.npz", run, other, cached, definition, path,
                                    (omega, [], 0.0, []), vis, meta, write_map=False)


def test_band_facets_take_the_resolver_classification():
    slope = np.array([10.0, 45.0 - 5e-7, 45.0 + 5e-7, 60.0, 70.0])
    # The resolver disagrees with the GPU on both side-of-45 facets and on 60.
    masks, resolved = restart_bundle.cliff_columns(
        slope, (45.0, 60.0), lambda faces: np.array([46.0, 44.0, 59.0])[:len(faces)])
    assert sorted(resolved) == [1, 2, 3]
    assert masks[45.0].tolist() == [False, True, False, True, True]
    assert masks[60.0].tolist() == [False, False, False, False, True]
    with pytest.raises(ValueError, match="no resolver"):
        restart_bundle.cliff_columns(slope, (45.0,))
    assert restart_bundle.cliff_name(52.5) == "cliff_ge_52p5deg"
