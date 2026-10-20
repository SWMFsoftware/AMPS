"""Publication figures (cops_cuda/figures.py, tools/pub_figures.py).

* the late-visibility default equals c05_write_bundle's ``--late-start/--late-stop``;
* geometry on a cube: nearest surface points, the impact extrapolation, and the GPU ray
  intercept (``kernels/render.cu`` ``ray_first_hit``);
* the GPU z-buffer keeps the nearest face per pixel;
* pressure figures are written as .eps and .png, and existing files are refused.
"""

import datetime as dt
import re

import numpy as np
import pytest

from cops_cuda import figures

APP_ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def cube():
    v = np.array([[x, y, z] for x in (-1.0, 1.0) for y in (-1.0, 1.0) for z in (-1.0, 1.0)])
    idx = {tuple(p): k for k, p in enumerate(v)}
    faces = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            others = [a for a in range(3) if a != axis]
            corners = []
            for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                p = [0.0, 0.0, 0.0]
                p[axis], p[others[0]], p[others[1]] = sign, a, b
                corners.append(idx[tuple(p)])
            faces += [[corners[0], corners[1], corners[2]], [corners[0], corners[2], corners[3]]]
    return v, np.array(faces, np.int64)


def test_default_window_matches_the_bundle_writer():
    text = (APP_ROOT / "tools" / "c05_write_bundle.py").read_text(encoding="utf-8")
    start = re.search(r'"--late-start", default="([^"]+)"', text).group(1)
    stop = re.search(r'"--late-stop", default="([^"]+)"', text).group(1)
    assert figures.DEFAULT_LATE_WINDOW == (start, stop)


def table(n=40):
    t0 = dt.datetime(2016, 9, 30, 10, 30, 0)
    time = np.array([t0 + dt.timedelta(seconds=20 * k) for k in range(n)])
    alt = np.linspace(500.0, 0.5, n)
    ref = np.exp(-alt / 300.0) * 1e-5
    return {"utc": [t.strftime("%Y-%m-%dT%H:%M:%S") for t in time], "time": time, "altitude": alt,
            "iPoint": np.arange(n, dtype=float), "inside": np.zeros(n), "refNG": ref, "refRG": 0.5 * ref,
            "useNG": (np.arange(n) % 5 != 0).astype(float), "useRG": np.ones(n), "pNG": 1.1 * ref,
            "pRG": 0.45 * ref, "et": np.arange(n) * 20.0,
            "x": np.zeros(n), "y": np.zeros(n), "z": 1.0 + alt}


def test_late_window_rows_and_altitudes():
    tb = table()
    w = figures.late_window(tb, start="2016-09-30T10:33:10", stop="2016-09-30T10:39:26")
    sel = (tb["time"] >= w["t0"]) & (tb["time"] <= w["t1"])
    assert w["points"] == [int(p) for p in tb["iPoint"][sel]]
    assert w["altitude_min_m"] == tb["altitude"][sel].min()
    assert w["altitude_max_m"] == tb["altitude"][sel].max()
    with pytest.raises(ValueError):
        figures.late_window(tb, start="2016-09-30T11:00:00", stop="2016-09-30T11:10:00")


def test_nearest_points_and_impact_on_a_cube():
    v, f = cube()
    xyz = np.array([[0.2, 0.3, 5.0], [3.0, 3.0, 0.0]])
    pts, dist, _ = figures.ground_track(v, f, xyz, k=12)
    np.testing.assert_allclose(pts[0], [0.2, 0.3, 1.0], atol=1e-12)
    np.testing.assert_allclose(dist[0], 4.0, atol=1e-12)
    np.testing.assert_allclose(pts[1], [1.0, 1.0, 0.0], atol=1e-12)        # an edge
    tb = table()
    tb["x"][:] = 0.25
    tb["y"][:] = 0.0
    tb["z"][-2:] = [1.9, 1.5]                                                # 0.5 m above the top face
    impact = figures.impact_point(v, f, tb, radius=5.0)
    np.testing.assert_allclose(impact["point_m"], [0.25, 0.0, 1.0], atol=1e-12)
    assert impact["seconds_after_last_row"] == pytest.approx(0.5 / (0.4 / 20.0))


cp = pytest.importorskip("cupy")
gpu = pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")


@gpu
def test_intercept_toward_the_origin(clean_pool):
    v, f = cube()
    xyz = np.array([[10.0, 0.0, 0.0], [0.0, -4.0, 4.0], [0.0, 0.0, 0.5]])  # last one is inside
    pts, face = figures.intercept_track(v, f, xyz)
    np.testing.assert_allclose(pts[0], [1.0, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(pts[1], [0.0, -1.0, 1.0], atol=1e-12)
    assert face[0] >= 0 and face[1] >= 0
    # from inside, the first hit along -z is the bottom face
    np.testing.assert_allclose(pts[2], [0.0, 0.0, -1.0], atol=1e-12)


@gpu
def test_zbuffer_keeps_the_nearest_face(clean_pool):
    # two squares (two triangles each) over the same pixels, at depths 5 and 2
    u = np.array([0.0, 10.0, 10.0, 0.0] * 2)
    v = np.array([0.0, 0.0, 10.0, 10.0] * 2)
    w = np.array([5.0] * 4 + [2.0] * 4)
    faces = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]], np.int64)
    face, depth = figures.rasterize(u, v, w, faces, 12, 12)
    covered = face >= 0
    assert covered[:10, :10].all() and not covered[10:, :].any() and not covered[:, 10:].any()
    assert set(np.unique(face[covered])) <= {2, 3}
    np.testing.assert_allclose(depth[covered], 2.0)
    face2, _ = figures.rasterize(u, v, w, faces, 12, 12, draw=[1, 1, 0, 0])
    assert set(np.unique(face2[face2 >= 0])) <= {0, 1}


def test_pressure_figures_are_written_once(tmp_path):
    tb = table()
    w = figures.late_window(tb)
    paths = figures.pressure_figures(tb, [("fit", tb)], tb, w, tmp_path / "p")
    names = sorted(p.split("/")[-1] for p in paths)
    assert names == ["p_altitude.eps", "p_altitude.png", "p_time.eps", "p_time.png"]
    for p in paths:
        data = open(p, "rb").read(8)
        assert data.startswith(b"%!PS") if p.endswith(".eps") else data.startswith(b"\x89PNG")
    with pytest.raises(FileExistsError):
        figures.pressure_figures(tb, [("fit", tb)], tb, w, tmp_path / "p")
