"""Content identities, exclusive product writes, and source-tree snapshots.

Every product written by this application is identified by content (SHA-256
and size), never by its filename.  Products are created exclusively: an
existing path is an error, so a rerun can never silently merge with or replace
an earlier reference.  The helpers are Python 3.8 compatible because the C00
capture runs under the CPU reference interpreter.
"""

import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys


CHUNK = 8 * 1024 * 1024

# Directories that hold generated or environment state, not source.  They are
# excluded from tree identities so that a cache refresh or a virtualenv does
# not change the recorded source revision.
TREE_EXCLUDED_DIRS = frozenset(("__pycache__", ".git", ".venv", ".pytest_cache",
                                "build", "RCS"))
TREE_EXCLUDED_SUFFIXES = (".pyc", ".swp", ".swo", "~")


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def utc_stamp():
    """Compact UTC timestamp used in unique run-directory names."""
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path, role=None):
    """Absolute path, size and SHA-256 of one file (the path is informative only)."""
    path = Path(path).resolve()
    identity = {"path": str(path), "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path)}
    if role is not None:
        identity["role"] = role
    return identity


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_sha256(value):
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def array_sha256(array, chunk_rows=1 << 20):
    """Hash an array's dtype, shape and C-order bytes (dtype includes byte order).

    A non-contiguous array (a memory-mapped column view) is hashed in row chunks: the same
    digest without a full contiguous copy."""
    import numpy as np
    array = np.asanyarray(array)
    digest = hashlib.sha256()
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(repr(tuple(array.shape)).encode("ascii"))
    if array.flags.c_contiguous or array.ndim == 0:
        digest.update(np.ascontiguousarray(array).tobytes(order="C"))
    else:
        for start in range(0, array.shape[0], chunk_rows):
            digest.update(np.ascontiguousarray(array[start:start + chunk_rows]).tobytes(order="C"))
    return digest.hexdigest()


def write_json_exclusive(path, value):
    """Write JSON to a new file; refuse to replace anything that exists."""
    path = Path(path)
    text = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with path.open("x") as stream:
        stream.write(text)
    return file_identity(path)


def make_new_directory(path):
    """Create a directory that must not exist yet (parents may exist)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.mkdir(exist_ok=False)
    return path.resolve()


def git_state(path):
    """Revision and dirty-file list of the git work tree containing ``path``."""
    path = Path(path).resolve()

    def run(*args):
        completed = subprocess.run(("git",) + args, cwd=str(path), text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return completed.stdout.strip() if completed.returncode == 0 else None

    top = run("rev-parse", "--show-toplevel")
    if top is None:
        return {"available": False}
    status = run("status", "--porcelain", "--untracked-files=all", "--", ".")
    return {"available": True, "toplevel": top, "head": run("rev-parse", "HEAD"),
            "describe": run("describe", "--always", "--dirty"),
            "status_porcelain": status.splitlines() if status else []}


def tree_identity(root, exclude_dirs=()):
    """Per-file identities of every source file below ``root``.

    The aggregate digest covers relative paths and file digests in sorted
    order, so it changes when any tracked or untracked source file changes,
    including files a version-control status would not report.
    """
    root = Path(root).resolve()
    skip = TREE_EXCLUDED_DIRS | frozenset(exclude_dirs)
    files = []
    for directory, subdirs, names in os.walk(str(root)):
        subdirs[:] = sorted(d for d in subdirs if d not in skip)
        for name in sorted(names):
            if name.endswith(TREE_EXCLUDED_SUFFIXES):
                continue
            path = Path(directory) / name
            if path.is_symlink() or not path.is_file():
                continue
            files.append({"path": str(path.relative_to(root)),
                          "size_bytes": path.stat().st_size,
                          "sha256": sha256_file(path)})
    digest = hashlib.sha256()
    for entry in files:
        digest.update(entry["path"].encode("utf-8") + b"\0" +
                      entry["sha256"].encode("ascii") + b"\n")
    return {"root": str(root), "file_count": len(files),
            "tree_sha256": digest.hexdigest(), "files": files}


def process_context():
    """Interpreter, host and thread-environment facts for a manifest."""
    threads = {name: os.environ.get(name) for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "COPS_OUTPUT_ROOT", "CUDA_VISIBLE_DEVICES",
        "CUDA_PATH")}
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except AttributeError:
        affinity = None
    return {"utc": utc_now(), "host": socket.gethostname(),
            "python": sys.version, "executable": sys.executable,
            "argv": list(sys.argv), "cwd": os.getcwd(),
            "cpu_affinity": affinity, "environment": threads}
