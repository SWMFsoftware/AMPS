"""Provenance and reference-binding contracts (CPU only; runs in both interpreters)."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import provenance as pv  # noqa: E402
from cops_cuda import reference  # noqa: E402


def test_tree_identity_tracks_content_and_ignores_caches(tmp_path):
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "__pycache__").mkdir()
    (tmp_path / "__pycache__" / "a.cpython-311.pyc").write_bytes(b"cache")
    first = pv.tree_identity(tmp_path)
    assert first["file_count"] == 1
    (tmp_path / "__pycache__" / "b.pyc").write_bytes(b"more cache")
    assert pv.tree_identity(tmp_path)["tree_sha256"] == first["tree_sha256"]
    (tmp_path / "a.py").write_text("x = 2\n")
    assert pv.tree_identity(tmp_path)["tree_sha256"] != first["tree_sha256"]


def test_exclusive_writes_refuse_existing_products(tmp_path):
    target = tmp_path / "product.json"
    identity = pv.write_json_exclusive(target, {"value": 1})
    assert identity["sha256"] == pv.sha256_file(target)
    with pytest.raises(FileExistsError):
        pv.write_json_exclusive(target, {"value": 2})
    assert json.loads(target.read_text()) == {"value": 1}
    pv.make_new_directory(tmp_path / "run")
    with pytest.raises(FileExistsError):
        pv.make_new_directory(tmp_path / "run")


def test_array_hash_includes_dtype_and_shape():
    np = pytest.importorskip("numpy")
    a = np.zeros(4, dtype="<f8")
    assert pv.array_sha256(a) != pv.array_sha256(a.reshape(2, 2))
    assert pv.array_sha256(a) != pv.array_sha256(a.astype(">f8"))


def test_reference_binding_loads_modules_from_the_named_tree():
    root = reference.reference_root()
    module = reference.load("cops_selection")
    assert root / "tools" in Path(module.__file__).resolve().parents


def test_reference_tree_must_contain_the_reference_tools(tmp_path):
    with pytest.raises(FileNotFoundError):
        reference.reference_root(tmp_path)


def test_rebinding_a_different_reference_is_refused(tmp_path):
    fake = tmp_path / "fake"
    for name in reference.REQUIRED_FILES:
        (fake / name).parent.mkdir(parents=True, exist_ok=True)
        (fake / name).write_text("")
    code = ("import sys; sys.path.insert(0, %r)\n"
            "from cops_cuda import reference\n"
            "reference.bind()\n"
            "try:\n    reference.bind(%r)\nexcept RuntimeError:\n    print('refused')\n"
            % (str(APP_ROOT), str(fake)))
    output = subprocess.run([sys.executable, "-c", code], text=True,
                            stdout=subprocess.PIPE, check=True).stdout
    assert output.strip() == "refused"
