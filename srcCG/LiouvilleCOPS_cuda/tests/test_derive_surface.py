"""tools/derive_outward_surface.py: exact global winding reversal (CPU only)."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT / "tools"))
import c00_capture_reference as c00  # noqa: E402

from cops_cuda import reference  # noqa: E402


def write_bdf(path, v, f, extra_cards=True):
    """Large-field GRID* cards (fixed columns) and free-field CTRIA3, like Gmsh."""
    with path.open("w") as out:
        out.write("$ Created by test\n")
        # GRID* fields are 16 characters wide: print 9 significant digits.
        for i, (x, y, z) in enumerate(v, start=1):
            out.write("GRID*   %-16d%-16d%-16.9e%-16.9e*N%d\r\n" % (i, 0, x, y, i))
            out.write("*N%-6d%-16.9e\r\n" % (i, z))
        for k, (a, b, c) in enumerate(f, start=1):
            out.write("CTRIA3  %d       1       %d       %d       %d       \r\n" % (k, a + 1, b + 1, c + 1))
        if extra_cards:
            out.write("CTETRA  99      1       1       2       3       5\r\n")


def run(src, dst):
    return subprocess.run([sys.executable, str(APP_ROOT / "tools" / "derive_outward_surface.py"),
                           str(src), str(dst)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)


def test_inward_cube_is_reversed_exactly(tmp_path):
    builder = reference.load("build_physical_activity_templates")
    v, f = c00.cube(1000.0)
    src = tmp_path / "inward.bdf"
    write_bdf(src, v + np.array([0.125, -3.5, 7.25]), f[:, [1, 0, 2]])   # inward winding
    V0, F0 = builder.read_nastran_surface(src)
    dst = tmp_path / "derived" / "outward.bdf"
    result = run(src, dst)
    assert result.returncode == 0, result.stderr
    V1, F1 = builder.read_nastran_surface(dst)
    assert V1.tobytes() == V0.tobytes()
    assert np.array_equal(F1, F0[:, [0, 2, 1]])
    record = json.loads(Path(str(dst) + ".provenance.json").read_text())
    assert all(record["verification"].values())
    assert record["topology_output"]["signed_volume_m3"] > 0.0
    # Exclusive: a second run refuses to replace the output.
    assert run(src, dst).returncode != 0


def test_outward_or_open_sources_are_refused(tmp_path):
    v, f = c00.cube(1000.0)
    outward = tmp_path / "outward.bdf"
    write_bdf(outward, v, f)
    assert run(outward, tmp_path / "x.bdf").returncode != 0          # already outward
    open_mesh = tmp_path / "open.bdf"
    write_bdf(open_mesh, v, f[1:, [1, 0, 2]])
    assert run(open_mesh, tmp_path / "y.bdf").returncode != 0        # not closed
