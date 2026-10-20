#!/usr/bin/env python3
"""Strict reusable per-facet gravity caches for physical activity templates.

The expensive quantity is homogeneous-polyhedron *self gravity* evaluated at
the active model facet centroids displaced along their normals.  It depends on
the modeled facet geometry, gravity polyhedron, density, and evaluation offset,
but not on the observation epoch.  Centrifugal/effective gravity and slope are
also saved as an auditable reference realization; they are recomputed cheaply
from cached self gravity if a later run supplies a different SPICE angular
velocity.

Every array is bound to an exact centroid/normal/area ordering and protected by
both the NPZ container CRC and an explicit SHA-256.  A cache therefore cannot
silently migrate to a different shape resolution or facet ordering.
"""

import hashlib
import json
from pathlib import Path

import numpy as np


FORMAT = "liouville_cops.physical_gravity_cache"
VERSION = 1
from gravity_topology import ORIENTATION_POLICY

# v1 caches may contain silently HEAL-flipped faces. Changing both fields
# prevents a corrected builder from accepting those values on hash equality.
ALGORITHM = "homogeneous_polyhedron_self_gravity_at_outward_offset_v2_disable_after_topology_audit"
REQUIRED_KEYS = frozenset((
    "format", "format_version", "geometry_identity", "definition_json",
    "self_gravity_m_s2", "reference_omega_body_rad_s",
    "reference_centrifugal_m_s2", "reference_effective_gravity_m_s2",
    "reference_slope_deg", "metadata_json",
))


def _little_endian(array):
    """Return a contiguous little-endian float64 representation for hashing."""
    return np.ascontiguousarray(np.asarray(array, dtype="<f8"))


def _array_sha256(array):
    return hashlib.sha256(_little_endian(array).tobytes(order="C")).hexdigest()


def geometry_identity(centroid, normal, area):
    """Identify the exact geometry and face ordering required by a cache."""
    centroid = _little_endian(centroid)
    normal = _little_endian(normal)
    area = _little_endian(area)
    if centroid.ndim != 2 or centroid.shape[1] != 3:
        raise ValueError("centroid must have shape (nfaces,3)")
    if normal.shape != centroid.shape or area.shape != (len(centroid),):
        raise ValueError("normal/area do not match the centroid face table")
    table = np.empty((len(area), 7), dtype="<f8")
    table[:, 0:3] = centroid
    table[:, 3:6] = normal
    table[:, 6] = area
    return hashlib.sha256(table.tobytes(order="C")).hexdigest()


def make_definition(centroid, normal, area, gravity_shape_identity,
                    gravity_vertices, gravity_triangles, density, offset):
    """Build the path-independent set of inputs that determines self gravity."""
    shape = {
        "size_bytes": int(gravity_shape_identity["size_bytes"]),
        "sha256": str(gravity_shape_identity["sha256"]),
    }
    return {
        "algorithm": ALGORITHM,
        "orientation_policy": ORIENTATION_POLICY,
        "geometry_identity": geometry_identity(centroid, normal, area),
        "face_count": int(len(area)),
        "gravity_shape": shape,
        "gravity_shape_vertices": int(gravity_vertices),
        "gravity_shape_triangles": int(gravity_triangles),
        "density_kg_m3": float(density),
        "evaluation_offset_m": float(offset),
    }


def derive_effective_fields(centroid, normal, self_gravity, omega):
    """Derive centrifugal/effective acceleration and gravitational slope."""
    centroid = np.asarray(centroid, dtype=float)
    normal = np.asarray(normal, dtype=float)
    self_gravity = np.asarray(self_gravity, dtype=float)
    omega = np.asarray(omega, dtype=float)
    if (centroid.ndim != 2 or centroid.shape[1] != 3 or
            normal.shape != centroid.shape or self_gravity.shape != centroid.shape or
            omega.shape != (3,)):
        raise ValueError("gravity field arrays have inconsistent shapes")
    centrifugal = -np.cross(omega, np.cross(omega, centroid))
    effective = self_gravity + centrifugal
    magnitude = np.linalg.norm(effective, axis=1)
    if (not np.all(np.isfinite(effective)) or not np.all(np.isfinite(normal)) or
            np.any(magnitude <= 0.0)):
        raise ValueError("nonfinite/zero effective gravity at a modeled facet")
    cosine = np.sum(-normal * effective, axis=1) / magnitude
    slope = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return centrifugal, effective, slope


def _scalar(archive, name):
    value = np.asarray(archive[name])
    if value.shape != ():
        raise ValueError("gravity cache field %s must be scalar" % name)
    return value.item()


def write_cache(path, centroid, normal, area, definition, self_gravity, omega,
                calculation_metadata):
    """Write a new cache and its reference epoch-dependent derived fields."""
    path = Path(path)
    if path.exists():
        raise ValueError("gravity cache already exists: %s" % path)
    expected_geometry = geometry_identity(centroid, normal, area)
    if definition.get("geometry_identity") != expected_geometry:
        raise ValueError("gravity definition does not match the supplied geometry")
    self_gravity = _little_endian(self_gravity)
    if self_gravity.shape != np.asarray(centroid).shape or not np.all(np.isfinite(self_gravity)):
        raise ValueError("self gravity must be a finite (nfaces,3) array")
    omega = _little_endian(omega)
    centrifugal, effective, slope = derive_effective_fields(
        centroid, normal, self_gravity, omega)
    centrifugal = _little_endian(centrifugal)
    effective = _little_endian(effective)
    slope = _little_endian(slope)
    arrays = {
        "self_gravity_m_s2": self_gravity,
        "reference_omega_body_rad_s": omega,
        "reference_centrifugal_m_s2": centrifugal,
        "reference_effective_gravity_m_s2": effective,
        "reference_slope_deg": slope,
    }
    metadata = {
        "format": FORMAT,
        "format_version": VERSION,
        "definition": definition,
        "array_sha256": {name: _array_sha256(value)
                          for name, value in arrays.items()},
        "calculation": calculation_metadata,
    }
    definition_json = json.dumps(definition, sort_keys=True, separators=(",", ":"))
    metadata_json = json.dumps(metadata, sort_keys=True, separators=(",", ":"))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        np.savez_compressed(
            stream, format=np.array(FORMAT),
            format_version=np.array(VERSION, dtype="<i8"),
            geometry_identity=np.array(expected_geometry),
            definition_json=np.array(definition_json),
            metadata_json=np.array(metadata_json), **arrays)
    return {"metadata": metadata, "self_gravity": self_gravity,
            "centrifugal": centrifugal, "effective": effective, "slope": slope,
            "reference_omega_reused": True}


def load_cache(path, centroid, normal, area, expected_definition, omega):
    """Load an exact self-gravity cache and derive fields for ``omega``."""
    path = Path(path)
    with np.load(str(path), allow_pickle=False) as archive:
        keys = frozenset(archive.files)
        if keys != REQUIRED_KEYS:
            raise ValueError("gravity cache fields are %s; expected %s" %
                             (sorted(keys), sorted(REQUIRED_KEYS)))
        if (str(_scalar(archive, "format")) != FORMAT or
                int(_scalar(archive, "format_version")) != VERSION):
            raise ValueError("unsupported gravity cache format/version")
        saved_geometry = str(_scalar(archive, "geometry_identity"))
        definition_text = str(_scalar(archive, "definition_json"))
        metadata_text = str(_scalar(archive, "metadata_json"))
        arrays = {name: _little_endian(archive[name]) for name in REQUIRED_KEYS
                  if name.endswith("_m_s2") or name.endswith("_rad_s") or
                  name == "reference_slope_deg"}
    expected_geometry = geometry_identity(centroid, normal, area)
    if saved_geometry != expected_geometry:
        raise ValueError("gravity cache geometry identity does not match the active mesh")
    try:
        definition = json.loads(definition_text)
        metadata = json.loads(metadata_text)
    except json.JSONDecodeError as exc:
        raise ValueError("invalid gravity cache JSON: %s" % exc) from None
    if definition != expected_definition:
        raise ValueError("gravity cache physical inputs do not match the requested inputs")
    if (metadata.get("format") != FORMAT or metadata.get("format_version") != VERSION or
            metadata.get("definition") != definition):
        raise ValueError("gravity cache metadata disagrees with its binary identity fields")

    nfaces = len(area)
    shapes = {
        "self_gravity_m_s2": (nfaces, 3),
        "reference_omega_body_rad_s": (3,),
        "reference_centrifugal_m_s2": (nfaces, 3),
        "reference_effective_gravity_m_s2": (nfaces, 3),
        "reference_slope_deg": (nfaces,),
    }
    for name, shape in shapes.items():
        value = arrays[name]
        if value.shape != shape or not np.all(np.isfinite(value)):
            raise ValueError("gravity cache %s is not a finite array of shape %s" %
                             (name, shape))
        if metadata.get("array_sha256", {}).get(name) != _array_sha256(value):
            raise ValueError("gravity cache checksum mismatch for %s" % name)

    reference_omega = arrays["reference_omega_body_rad_s"]
    reference_centrifugal, reference_effective, reference_slope = \
        derive_effective_fields(centroid, normal, arrays["self_gravity_m_s2"],
                                reference_omega)
    for name, saved, derived in (
            ("reference_centrifugal_m_s2", arrays["reference_centrifugal_m_s2"],
             reference_centrifugal),
            ("reference_effective_gravity_m_s2",
             arrays["reference_effective_gravity_m_s2"], reference_effective),
            ("reference_slope_deg", arrays["reference_slope_deg"], reference_slope)):
        if not np.array_equal(saved, _little_endian(derived)):
            raise ValueError("gravity cache derived reference field is inconsistent: %s" % name)

    omega = _little_endian(omega)
    same_omega = np.array_equal(omega, reference_omega)
    if same_omega:
        centrifugal = arrays["reference_centrifugal_m_s2"]
        effective = arrays["reference_effective_gravity_m_s2"]
        slope = arrays["reference_slope_deg"]
    else:
        centrifugal, effective, slope = derive_effective_fields(
            centroid, normal, arrays["self_gravity_m_s2"], omega)
    return {"metadata": metadata,
            "self_gravity": arrays["self_gravity_m_s2"],
            "centrifugal": centrifugal, "effective": effective, "slope": slope,
            "reference_omega_reused": same_omega}
