#!/usr/bin/env python3
"""Strict I/O, hashing, and manifest helpers for LiouvilleCOPS products.

This module contains no scientific model and deliberately makes no choices
about observations, calibration, parameter bounds, or physical validity.  Its
job is narrower: identify bytes, reject internally inconsistent caches/tables,
and record exactly which artifacts were used.  Large files are hashed in
bounded chunks so provenance does not require loading them into memory.

Manifest states are intentionally independent.  A completed process is not
automatically a converged optimizer, a verified cache, a C++-reproduced fit,
or a scientifically reviewed result.
"""

import datetime as _datetime
import hashlib
import json
import math
import os
from pathlib import Path
import struct
from collections import Counter


RUN_SCHEMA = "liouville_cops.run_manifest"
FIT_SCHEMA = "liouville_cops.fit_manifest"
SCHEMA_VERSION = 1
MODEL_COLUMNS = (
    "iPoint", "UTC", "et", "x", "y", "z", "altitude", "inside",
    "rh_AU", "VdotB_RG", "nVisible", "Q_H2O", "Q_CO2", "n_H2O",
    "n_CO2", "PhiRG_H2O", "PhiRG_CO2", "pNG_H2O", "pNG_CO2", "pNG",
    "pRG_H2O", "pRG_CO2", "pRG", "refNG", "refRG",
    "NGBelowDeckFraction", "T_RG", "useNG", "useRG",
)
ELEMENT_HEADER = struct.Struct("<iii")
KERNEL1_HEADER = struct.Struct("<ii10di")
KERNEL1_RECORD = struct.Struct("<i4di")


class ProductError(ValueError):
    """A product is incomplete, corrupt, ambiguous, or mutually inconsistent."""


def utc_now():
    """UTC timestamp with timezone, suitable for an audit record."""
    return _datetime.datetime.now(_datetime.timezone.utc).isoformat()


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    """Return the SHA-256 of *path* using bounded memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while True:
            block = stream.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value):
    """Hash a JSON-compatible value after deterministic serialization."""
    data = json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def file_identity(path, role=None):
    """Return an absolute path, byte size, and streaming SHA-256 identity."""
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise ProductError("required file is missing: %s" % path)
    result = {"path": str(path), "size_bytes": path.stat().st_size,
              "sha256": sha256_file(path)}
    if role is not None:
        result["role"] = role
    return result


def write_json_atomic(path, value):
    """Replace a JSON document atomically after fully writing its sibling temp file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.%d" % os.getpid())
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def load_manifest(path, schema=None):
    """Read and minimally validate a versioned manifest."""
    path = Path(path)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ProductError("cannot read manifest %s: %s" % (path, exc)) from None
    if not isinstance(value, dict):
        raise ProductError("manifest is not a JSON object: %s" % path)
    if value.get("schema_version") != SCHEMA_VERSION:
        raise ProductError("unsupported manifest schema version in %s: %r" %
                           (path, value.get("schema_version")))
    if schema is not None and value.get("schema") != schema:
        raise ProductError("unexpected manifest schema in %s: %r" %
                           (path, value.get("schema")))
    return value


def verify_identity(identity, label=None):
    """Verify that a file still matches an identity stored in a manifest."""
    label = label or identity.get("role") or identity.get("path", "file")
    actual = file_identity(identity["path"])
    if actual["size_bytes"] != identity.get("size_bytes"):
        raise ProductError("%s size changed since the manifest was written" % label)
    if actual["sha256"] != identity.get("sha256"):
        raise ProductError("%s SHA-256 changed since the manifest was written" % label)


def read_model_table(path, required_columns=()):
    """Read a COPS model table and reject duplicate/malformed observations.

    Numeric physics columns are not interpreted here.  The checks are format
    invariants: one declared schema, exact row width, unique integer point IDs,
    and required columns.  Equal UTC values are allowed because simultaneous
    records may be legitimate.  Returning strings preserves the
    serialized values used by existing readers.
    """
    path = Path(path)
    headers, names, rows = [], None, []
    try:
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if line.startswith("#iPoint"):
                    if names is not None:
                        raise ProductError("%s has more than one column header" % path)
                    names = line[1:].split()
                elif line.startswith("#"):
                    headers.append(line[1:].strip())
                elif line.strip():
                    if names is None:
                        raise ProductError("%s:%d has data before the column header" %
                                           (path, line_number))
                    row = line.split()
                    if len(row) != len(names):
                        raise ProductError("%s:%d has %d fields; expected %d" %
                                           (path, line_number, len(row), len(names)))
                    rows.append(row)
    except OSError as exc:
        raise ProductError("cannot read %s: %s" % (path, exc)) from None
    if names is None:
        raise ProductError("%s has no #iPoint column header" % path)
    if len(set(names)) != len(names):
        raise ProductError("%s has duplicate column names" % path)
    missing = [name for name in required_columns if name not in names]
    if missing:
        raise ProductError("%s is missing columns: %s" % (path, ", ".join(missing)))
    columns = {name: index for index, name in enumerate(names)}
    point_ids = []
    for row_number, row in enumerate(rows, 1):
        try:
            point_ids.append(int(row[columns["iPoint"]]))
        except (KeyError, ValueError):
            raise ProductError("%s data row %d has an invalid iPoint" %
                               (path, row_number)) from None
    duplicate_ids = sorted(value for value, count in Counter(point_ids).items() if count > 1)
    if duplicate_ids:
        raise ProductError("%s has duplicate iPoint values (first: %d)" %
                           (path, duplicate_ids[0]))
    # Equal UTC values are intentionally allowed: simultaneous records could
    # be legitimate.  iPoint is the authoritative unique key.
    if not rows:
        raise ProductError("%s has no data rows" % path)
    return headers, names, rows, columns


def _finite(values, label):
    if not all(math.isfinite(value) for value in values):
        raise ProductError("non-finite value in %s" % label)


def validate_elements(path):
    """Validate the complete version-1 element cache and return its identity."""
    path = Path(path)
    size = path.stat().st_size if path.is_file() else -1
    with path.open("rb") as stream:
        header = stream.read(ELEMENT_HEADER.size)
        if len(header) != ELEMENT_HEADER.size:
            raise ProductError("truncated elements header: %s" % path)
        version, nfaces, nbasis = ELEMENT_HEADER.unpack(header)
        if version != 1 or nfaces <= 0 or nbasis != 25:
            raise ProductError("unsupported elements header in %s: version=%d faces=%d basis=%d" %
                               (path, version, nfaces, nbasis))
        expected = ELEMENT_HEADER.size + nfaces * 34 * 8
        if size != expected:
            raise ProductError("elements size mismatch in %s: %d bytes; expected %d" %
                               (path, size, expected))
        digest = hashlib.sha256(header)
        for face in range(nfaces):
            raw = stream.read(34 * 8)
            digest.update(raw)
            values = struct.unpack("<34d", raw)
            _finite(values, "%s face %d" % (path, face))
            if values[6] <= 0.0:
                raise ProductError("non-positive area in %s face %d" % (path, face))
    return {"path": str(path.resolve()), "size_bytes": size,
            "sha256": digest.hexdigest(), "version": version,
            "nfaces": nfaces, "nbasis": nbasis}


def validate_kernel(path, expected_point, nfaces):
    """Validate one complete static kernel, including record face indices."""
    path = Path(path)
    size = path.stat().st_size if path.is_file() else -1
    with path.open("rb") as stream:
        raw_header = stream.read(KERNEL1_HEADER.size)
        if len(raw_header) != KERNEL1_HEADER.size:
            raise ProductError("truncated kernel header: %s" % path)
        unpacked = KERNEL1_HEADER.unpack(raw_header)
        version, point = unpacked[:2]
        values, count = unpacked[2:12], unpacked[12]
        if version != 1:
            raise ProductError("%s has kernel version %d; static fitting needs version 1" %
                               (path, version))
        if point != expected_point:
            raise ProductError("%s stores iPoint %d; expected %d" %
                               (path, point, expected_point))
        if count < 0 or count > nfaces:
            raise ProductError("%s has invalid record count %d for %d faces" %
                               (path, count, nfaces))
        expected_size = KERNEL1_HEADER.size + count * KERNEL1_RECORD.size
        if size != expected_size:
            raise ProductError("kernel size mismatch in %s: %d bytes; expected %d" %
                               (path, size, expected_size))
        _finite(values, "%s header" % path)
        digest = hashlib.sha256(raw_header)
        seen = set()
        for record_number in range(count):
            raw = stream.read(KERNEL1_RECORD.size)
            digest.update(raw)
            face, kn, kphi, kv, cosine, shadow = KERNEL1_RECORD.unpack(raw)
            if not 0 <= face < nfaces:
                raise ProductError("%s record %d has out-of-range face %d" %
                                   (path, record_number, face))
            if face in seen:
                raise ProductError("%s contains duplicate face %d" % (path, face))
            seen.add(face)
            _finite((kn, kphi, kv, cosine), "%s record %d" % (path, record_number))
            if shadow not in (0, 1):
                raise ProductError("%s record %d has invalid shadow flag %d" %
                                   (path, record_number, shadow))
    return {"point": point, "records": count, "size_bytes": size,
            "sha256": digest.hexdigest(), "header": values}


def validate_fit_base(run):
    """Validate table/elements/all kernels and return a compact cache audit.

    The combined kernel digest includes each filename and file digest in point
    order.  This identifies the set without placing thousands of hashes in the
    manifest.  Exact size checks reject truncation and appended bytes.
    """
    run = Path(run).resolve()
    model_path = run / "COPS_model.dat"
    headers, names, rows, columns = read_model_table(
        model_path, required_columns=("iPoint", "UTC", "et", "x", "y", "z",
                                     "altitude", "rh_AU", "VdotB_RG", "pNG",
                                     "pRG", "refNG", "refRG", "T_RG",
                                     "useNG", "useRG"))
    point_ids = [int(row[columns["iPoint"]]) for row in rows]
    elements = validate_elements(run / "elements.bin")
    kernel_dir = run / "kernels"
    if not kernel_dir.is_dir():
        raise ProductError("kernel directory is missing: %s" % kernel_dir)
    expected_names = {"point=%04d.bin" % point for point in point_ids}
    actual_names = {path.name for path in kernel_dir.glob("point=*.bin")}
    missing, extra = sorted(expected_names - actual_names), sorted(actual_names - expected_names)
    if missing or extra:
        detail = []
        if missing:
            detail.append("missing %d (first %s)" % (len(missing), missing[0]))
        if extra:
            detail.append("unexpected %d (first %s)" % (len(extra), extra[0]))
        raise ProductError("kernel/table point mismatch: " + "; ".join(detail))
    kernel_set = hashlib.sha256()
    records, total_bytes = 0, 0
    for row, point in zip(rows, point_ids):
        path = kernel_dir / ("point=%04d.bin" % point)
        info = validate_kernel(path, point, elements["nfaces"])
        records += info["records"]
        total_bytes += info["size_bytes"]
        kernel_set.update(path.name.encode("ascii"))
        kernel_set.update(info["sha256"].encode("ascii"))
        # These are serialized copies of the same geometry.  The table has
        # fewer printed digits, so compare at a tolerance set by its format.
        header_values = info["header"]
        comparisons = (("et", header_values[0], 2e-6),
                       ("VdotB_RG", header_values[1], 2e-6),
                       ("x", header_values[2], 2e-6),
                       ("y", header_values[3], 2e-6),
                       ("z", header_values[4], 2e-6),
                       ("rh_AU", header_values[8] / 149597870700.0, 2e-8),
                       ("altitude", header_values[9], 2e-6))
        for name, cached, relative_tolerance in comparisons:
            table_value = float(row[columns[name]])
            tolerance = relative_tolerance * max(1.0, abs(cached))
            if abs(table_value - cached) > tolerance:
                raise ProductError("point %d %s differs between table and kernel (%r vs %r)" %
                                   (point, name, table_value, cached))
    manifest_path = run / "run_manifest.json"
    manifest = None
    legacy = not manifest_path.is_file()
    if not legacy:
        manifest = load_manifest(manifest_path, RUN_SCHEMA)
        if manifest.get("states", {}).get("process") != "complete":
            raise ProductError("base run manifest does not mark the process complete")
        recorded = manifest.get("outputs", {})
        for key, actual in (("model_table", file_identity(model_path)),
                            ("elements", elements)):
            old = recorded.get(key)
            if old and (old.get("size_bytes") != actual["size_bytes"] or
                        old.get("sha256") != actual["sha256"]):
                raise ProductError("%s does not match run_manifest.json" % key)
        old_kernel = recorded.get("kernels", {})
        if old_kernel and old_kernel.get("set_sha256") != kernel_set.hexdigest():
            raise ProductError("kernel set does not match run_manifest.json")
    return {
        "run": str(run), "legacy_manifest": legacy, "manifest": manifest,
        "model_table": file_identity(model_path), "elements": elements,
        "kernels": {"version": 1, "count": len(point_ids),
                    "records": records, "size_bytes": total_bytes,
                    "set_sha256": kernel_set.hexdigest()},
        "point_ids": point_ids, "headers": headers, "names": names,
        "rows": rows, "columns": columns,
    }


def selection_identity(point_ids, utc, ref_ng, ref_rg, use_ng, use_rg, t_rg):
    """Hash the exact selected observations without imposing a selection policy."""
    rows = [[int(ip), str(time), float(ng), float(rg), int(ung), int(urg), float(temp)]
            for ip, time, ng, rg, ung, urg, temp in
            zip(point_ids, utc, ref_ng, ref_rg, use_ng, use_rg, t_rg)]
    return {"sha256": canonical_sha256(rows), "rows": len(rows),
            "used_ng": sum(row[4] for row in rows),
            "used_rg": sum(row[5] for row in rows)}


def manifest_summary(path):
    """Return display-only product/status information for plot legends."""
    path = Path(path)
    run_manifest = path.parent / "run_manifest.json"
    fit_manifest = path.parent / "fit_manifest.json"
    # Lambda-scan products live one level below their shared fit manifest.
    if not fit_manifest.is_file() and path.parent.name.startswith("lambda="):
        candidate = path.parent.parent / "fit_manifest.json"
        if candidate.is_file():
            fit_manifest = candidate
    if fit_manifest.is_file():
        value = load_manifest(fit_manifest, FIT_SCHEMA)
        return "fit", value.get("states", {})
    if run_manifest.is_file():
        value = load_manifest(run_manifest, RUN_SCHEMA)
        return value.get("run_kind", "forward"), value.get("states", {})
    return "legacy-unmanifested", {}
