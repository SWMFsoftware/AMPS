#!/usr/bin/env python3
"""Versioned per-facet activity-template bundles for LiouvilleCOPS fits.

The bundle is a NumPy ``.npz`` file, but its contents are deliberately strict:
no pickled objects, one canonical metadata JSON string, a fixed field list,
and an exact SHA-256 identity of the centroid/area table from ``elements.bin``.
This prevents a localized template made for one shape resolution or facet
ordering from being applied to another mesh merely because the face counts
happen to agree.
"""

import hashlib
import json
from pathlib import Path
import re

import numpy as np


FORMAT = "liouville_cops.activity_templates"
VERSION = 1
REQUIRED_KEYS = frozenset(("format", "format_version", "names", "values",
                           "mesh_identity", "metadata_json"))


def mesh_identity(centroid, area):
    """Hash exact little-endian float64 centroid/area records in face order."""
    centroid = np.asarray(centroid, dtype="<f8")
    area = np.asarray(area, dtype="<f8")
    if centroid.ndim != 2 or centroid.shape[1] != 3:
        raise ValueError("centroid must have shape (nfaces,3)")
    if area.shape != (len(centroid),):
        raise ValueError("area must have shape (nfaces,)")
    table = np.empty((len(area), 4), dtype="<f8")
    table[:, :3] = centroid
    table[:, 3] = area
    return hashlib.sha256(table.tobytes(order="C")).hexdigest()


def _scalar(array, name):
    value = np.asarray(array)
    if value.shape != ():
        raise ValueError("template bundle field %s must be scalar" % name)
    return value.item()


def validate_name(name):
    """Return a stable parameter-safe name or reject ambiguous spelling."""
    if not isinstance(name, str) or not re.fullmatch(r"[a-z][a-z0-9_]{0,63}", name):
        raise ValueError("template names must match [a-z][a-z0-9_]{0,63}: %r" % name)
    return name


def write_bundle(path, names, values, centroid, area, metadata):
    """Write a new compressed bundle; existing paths are never replaced."""
    path = Path(path)
    if path.exists():
        raise ValueError("template bundle already exists: %s" % path)
    names = [validate_name(str(name)) for name in names]
    if len(set(names)) != len(names) or not names:
        raise ValueError("template names must be nonempty and unique")
    values = np.asarray(values, dtype="<f8")
    if values.shape != (len(area), len(names)):
        raise ValueError("template matrix shape %s does not match (%d,%d)" %
                         (values.shape, len(area), len(names)))
    if not np.all(np.isfinite(values)):
        raise ValueError("template matrix contains nonfinite values")
    if np.any(values < 0.0):
        raise ValueError("physical enhancement templates must be nonnegative")
    payload = dict(metadata)
    payload.update({"format": FORMAT, "format_version": VERSION,
                    "names": names, "face_count": int(len(area)),
                    "mesh_identity": mesh_identity(centroid, area)})
    metadata_json = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    path.parent.mkdir(parents=True, exist_ok=True)
    # Passing an open file prevents numpy from silently appending another
    # suffix and makes the caller's path the exact manifested artifact.
    with path.open("xb") as stream:
        np.savez_compressed(stream, format=np.array(FORMAT),
                            format_version=np.array(VERSION, dtype="<i8"),
                            names=np.asarray(names, dtype="U64"), values=values,
                            mesh_identity=np.array(payload["mesh_identity"]),
                            metadata_json=np.array(metadata_json))
    return payload


def load_bundle(path, centroid, area):
    """Load and fully validate a bundle against the active element cache."""
    path = Path(path)
    with np.load(str(path), allow_pickle=False) as archive:
        keys = frozenset(archive.files)
        if keys != REQUIRED_KEYS:
            raise ValueError("template bundle fields are %s; expected %s" %
                             (sorted(keys), sorted(REQUIRED_KEYS)))
        file_format = str(_scalar(archive["format"], "format"))
        version = int(_scalar(archive["format_version"], "format_version"))
        if file_format != FORMAT or version != VERSION:
            raise ValueError("unsupported template bundle %s version %s" %
                             (file_format, version))
        names = [validate_name(str(value)) for value in archive["names"].tolist()]
        if not names or len(set(names)) != len(names):
            raise ValueError("template names must be nonempty and unique")
        values = np.asarray(archive["values"], dtype=float)
        if values.shape != (len(area), len(names)):
            raise ValueError("template matrix has shape %s; active mesh needs (%d,%d)" %
                             (values.shape, len(area), len(names)))
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("template matrix must be finite and nonnegative")
        expected_identity = mesh_identity(centroid, area)
        saved_identity = str(_scalar(archive["mesh_identity"], "mesh_identity"))
        if saved_identity != expected_identity:
            raise ValueError("template mesh identity %s does not match active mesh %s" %
                             (saved_identity, expected_identity))
        metadata_text = str(_scalar(archive["metadata_json"], "metadata_json"))
    try:
        metadata = json.loads(metadata_text)
    except json.JSONDecodeError as exc:
        raise ValueError("invalid template metadata JSON: %s" % exc) from None
    if (metadata.get("format") != FORMAT or metadata.get("format_version") != VERSION or
            metadata.get("names") != names or metadata.get("face_count") != len(area) or
            metadata.get("mesh_identity") != expected_identity):
        raise ValueError("template metadata disagrees with the binary bundle fields")
    return {"path": str(path.resolve()), "names": names, "values": values,
            "mesh_identity": expected_identity, "metadata": metadata}


def select_templates(bundle, requested):
    """Return requested columns in caller order, rejecting duplicates."""
    requested = [validate_name(name) for name in requested]
    if len(set(requested)) != len(requested):
        raise ValueError("the same activity template was requested more than once")
    lookup = {name: index for index, name in enumerate(bundle["names"])}
    missing = [name for name in requested if name not in lookup]
    if missing:
        raise ValueError("template(s) not in %s: %s" %
                         (bundle["path"], ", ".join(missing)))
    # Selection, rather than loading, is the boundary for this check: the
    # geometry-only NG column remains scientifically valid in old bundles.
    # Their cliff columns do not. A correct geometry hash alone cannot detect
    # HEAL's in-memory face reversals, so require the repaired numerical policy
    # even when an old bundle has been renamed or moved to another machine.
    definitions = bundle.get("metadata", {}).get("templates", {})
    gravity_requested = any(
        name.startswith("cliff_") or
        "effective_gravitational_slope" in str(definitions.get(name, {}).get("definition", ""))
        for name in requested)
    if gravity_requested:
        from gravity_cache import ALGORITHM
        from gravity_topology import ORIENTATION_POLICY
        gravity = bundle.get("metadata", {}).get("gravity", {})
        definition = gravity.get("cache", {}).get("definition", {})
        if (definition.get("algorithm") != ALGORITHM or
                definition.get("orientation_policy") != ORIENTATION_POLICY or
                gravity.get("integrity_check") != "DISABLE after topology audit"):
            raise ValueError(
                "cliff template lacks repaired gravity orientation provenance; "
                "HEAL-era cliff products are invalid and must be regenerated "
                "with DISABLE after an independent topology audit")
    return np.ascontiguousarray(bundle["values"][:, [lookup[name] for name in requested]])
