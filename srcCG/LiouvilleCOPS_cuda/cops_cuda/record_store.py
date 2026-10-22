"""FP32 device-ready record store (SCALING_PLAN.md 3.1; derived product, user decision Q4).

A base's transport-kernel records, converted once to the FP32 device layout of the mixed
precision experiment (``device_store.fp32_record_fields``: face int32; Kn, KPhi, KV, cos
float32; shadow folded into cos) and written as one file per observation point:

    <store>/point=NNNN.f32      face[n] kn[n] kphi[n] kv[n] cos[n]   (little endian, 20 n bytes)
    <store>/store_manifest.json

Why: a fit then reads 20 instead of 40 bytes per record and skips the per-record Python
conversion and validation of the loader (12M: 715 s per bind), and any rank can read exactly
its rows.

Trust chain:
* ``build`` validates every kernel file as the loader does (header against the table, record
  checks, SHA-256 of the whole file) and records, per point, the kernel file's size, mtime and
  SHA-256 and the store file's SHA-256. The kernel-set digest must equal ``run_manifest.json``.
* ``open_store`` refuses a store whose base ``run_manifest.json``, mesh identity or face count
  differ, and ``check_rows`` refuses a point whose kernel file changed (size or mtime) since the
  store was built. A loader verifies the SHA-256 of every store file it reads.
* The CPU-format kernel files remain the reference; the store is never read by the CPU tools
  and can be deleted and rebuilt at any time.
"""

import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import formats, provenance
from .errors import CacheFormatError, IdentityError

SCHEMA = "liouville_cops_cuda.record_store"
SCHEMA_VERSION = 1
FIELDS = (("face", "<i4"), ("kn", "<f4"), ("kphi", "<f4"), ("kv", "<f4"), ("cos", "<f4"))
BYTES_PER_RECORD = 20
MANIFEST = "store_manifest.json"


def store_file_name(point):
    return "point=%04d.f32" % int(point)


def split_fields(buffer, count):
    """Field views of one point's store bytes (no copies)."""
    out, offset = {}, 0
    for name, dtype in FIELDS:
        out[name] = np.frombuffer(buffer, dtype=dtype, count=count, offset=offset)
        offset += 4 * count
    return out


def _convert_point(inputs, scan, index, cos75, out):
    from .device_store import fp32_record_fields
    path = scan.paths[index]
    point = scan.point_ids[index]
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise IdentityError("%s changed while it was read" % path)
    head = formats.KERNEL_HEADER_DTYPE.itemsize
    if data[:head] != scan.headers[index]:
        raise IdentityError("%s changed between scan and conversion" % path)
    count = int(scan.counts[index])
    if len(data) != head + count * formats.KERNEL_RECORD_DTYPE.itemsize:
        raise CacheFormatError("%s size differs from its header" % path)
    records = np.frombuffer(data, dtype=formats.KERNEL_RECORD_DTYPE, count=count, offset=head)
    face_seen = np.zeros(inputs.nfaces, dtype=bool)
    formats.validate_records(records, inputs.nfaces, path, 0, face_seen)
    conversion = {"max_relative_rounding": 0.0}
    fields = fp32_record_fields(records, cos75, conversion)
    target = out / store_file_name(point)
    digest = hashlib.sha256()
    with open(target, "xb") as handle:
        for name, dtype in FIELDS:
            chunk = np.ascontiguousarray(fields[name], dtype=dtype).tobytes()
            digest.update(chunk)
            handle.write(chunk)
    return {"row": int(scan.rows[index]), "point": int(point), "kernel_file": path.name,
            "kernel_size": int(after.st_size), "kernel_mtime_ns": int(after.st_mtime_ns),
            "kernel_sha256": hashlib.sha256(data).hexdigest(), "records": count,
            "store_file": target.name, "store_sha256": digest.hexdigest(),
            "conversion": conversion}


def build(run, out, threads=8, progress=None):
    """Convert every kernel file of ``run`` into a new store directory ``out``."""
    from . import reference
    from .device_store import KernelScan, RunInputs
    run, out = Path(run).resolve(), Path(out).resolve()
    inputs = RunInputs(run)
    rows = inputs.resolve_rows(None, None)
    scan = KernelScan(inputs, rows)
    cos75 = float(reference.load("fit_cops").COS75)
    out.mkdir(parents=True, exist_ok=False)
    entries = [None] * len(rows)
    with ThreadPoolExecutor(max(1, int(threads)), thread_name_prefix="cops-store") as pool:
        futures = [pool.submit(_convert_point, inputs, scan, i, cos75, out) for i in range(len(rows))]
        for i, future in enumerate(futures):
            entries[i] = future.result()
            if progress:
                progress(i + 1, len(rows))
    kernel_set = formats.kernel_set_digest([(e["kernel_file"], e["kernel_sha256"]) for e in entries])
    recorded = (inputs.manifest or {}).get("outputs", {}).get("kernels", {})
    if recorded.get("set_sha256") and recorded["set_sha256"] != kernel_set:
        raise IdentityError("kernel set %s differs from run_manifest.json %s" %
                            (kernel_set, recorded["set_sha256"]))
    conversion = {"max_relative_rounding": max(e["conversion"]["max_relative_rounding"] for e in entries),
                  "lit_threshold_adjusted": sum(e["conversion"].get("lit_threshold_adjusted", 0)
                                                for e in entries)}
    from . import precision
    manifest = {
        "schema": SCHEMA, "schema_version": SCHEMA_VERSION, "created_utc": provenance.utc_now(),
        "base_run": str(run),
        "run_manifest": provenance.file_identity(run / "run_manifest.json"),
        "mesh_identity": inputs.mesh_identity, "nfaces": int(inputs.nfaces),
        "kernel_set_sha256": kernel_set,
        "kernel_set_verified_against_run_manifest": bool(recorded.get("set_sha256")),
        "record_precision": "fp32", "policy": precision.MIXED_POLICY_NAME, "cos75": cos75,
        "layout": {"fields": [n for n, _ in FIELDS], "dtypes": [d for _, d in FIELDS],
                   "order": "field-major within each point file", "bytes_per_record": BYTES_PER_RECORD},
        "records": int(sum(e["records"] for e in entries)), "points": entries,
        "conversion": conversion,
        "tool": provenance.file_identity(Path(__file__), "module"),
    }
    with open(out / MANIFEST, "x") as handle:          # written last: its presence marks completion
        json.dump(manifest, handle, indent=1)
    return manifest


class Store:
    """An opened store, checked against the base it is used with."""

    def __init__(self, path, inputs):
        self.path = Path(path).resolve()
        manifest_path = self.path / MANIFEST
        if not manifest_path.is_file():
            raise CacheFormatError("no %s in %s (incomplete or not a record store)" % (MANIFEST, self.path))
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("schema") != SCHEMA or manifest.get("schema_version") != SCHEMA_VERSION:
            raise CacheFormatError("%s is not a version-%d record store" % (self.path, SCHEMA_VERSION))
        current = provenance.file_identity(Path(inputs.run) / "run_manifest.json")
        if manifest["run_manifest"].get("sha256") != current.get("sha256"):
            raise IdentityError("record store %s was built for another base (run_manifest.json differs)"
                                % self.path)
        if manifest["mesh_identity"] != inputs.mesh_identity or manifest["nfaces"] != inputs.nfaces:
            raise IdentityError("record store %s was built for another mesh" % self.path)
        self.manifest = manifest
        self.by_point = {e["point"]: e for e in manifest["points"]}

    def check_rows(self, scan):
        """Entries of the scanned rows; refuse a point whose kernel file changed since the build."""
        out = []
        for index, point in enumerate(scan.point_ids):
            entry = self.by_point.get(int(point))
            if entry is None:
                raise IdentityError("record store has no point %d" % point)
            if entry["records"] != int(scan.counts[index]):
                raise IdentityError("record store point %d has %d records; the base has %d" %
                                    (point, entry["records"], int(scan.counts[index])))
            st = scan.paths[index].stat()
            if (st.st_size, st.st_mtime_ns) != (entry["kernel_size"], entry["kernel_mtime_ns"]):
                raise IdentityError("%s changed after the record store was built (rebuild the store)"
                                    % scan.paths[index])
            out.append(entry)
        return out

    def read_point(self, entry):
        """(fields dict, raw bytes) of one point, verified against the recorded SHA-256."""
        path = self.path / entry["store_file"]
        data = path.read_bytes()
        if len(data) != BYTES_PER_RECORD * entry["records"]:
            raise CacheFormatError("%s has %d bytes; expected %d" %
                                   (path, len(data), BYTES_PER_RECORD * entry["records"]))
        if hashlib.sha256(data).hexdigest() != entry["store_sha256"]:
            raise CacheFormatError("%s differs from its recorded SHA-256" % path)
        return split_fields(data, entry["records"])


def environment_store():
    """``$COPS_CUDA_RECORD_STORE`` (inherited by the multi-GPU workers), or None."""
    value = os.environ.get("COPS_CUDA_RECORD_STORE", "")
    return value or None
