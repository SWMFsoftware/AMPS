"""GPU-produced gravity caches that the unchanged CPU workflow reuses (ROADMAP C04).

The user's constraint (2026-10-19): nothing in ``../LiouvilleCOPS`` changes its
workflow. A GPU self-gravity field therefore becomes a CPU-readable product
only through the CPU's own code:

* ``gravity_cache.make_definition``: the path-independent definition (algorithm
  ``…_v2_disable_after_topology_audit``, orientation policy, exact facet
  geometry identity, shape hash, density, offset). It is backend-independent,
  because both backends evaluate the same mathematical definition.
* ``gravity_cache.write_cache``: the exclusive, hashed NPZ writer.

The *calculation* metadata says truthfully that the field was computed by the
CUDA translation. It also carries the values that the CPU builder's reuse
path and the template selector check:

* ``orientation_policy``: ``gravity_topology.ORIENTATION_POLICY``;
* ``topology``: exactly ``gravity_topology.require_outward_closed_body(V, F)``
  for the evaluated body, as kept in ``GravityBody.reference_audit``;
* ``integrity_check``: ``"DISABLE after topology audit"``. The GPU never alters
  the input winding, the meaning of DISABLE, and admits the body only after
  that audit;
* the offset-halving audit of the builder's own 512-facet sample, computed by
  the GPU.

No CPU package identity is copied: ``package_version``, ``package_binary`` and
``package_commit`` are absent, and ``backend`` and ``evidence`` name the CUDA
code, the device and the equivalence report that admitted the field.
"""

import numpy as np

from . import provenance, reference

INTEGRITY_CHECK = "DISABLE after topology audit"
BACKEND = "cuda_polyhedral_gravity_v1 (LiouvilleCOPS_cuda; GPLv3 translation of polyhedral-gravity 3.3)"


def offset_sample(n):
    """The builder's deterministic offset-halving sample (``effective_slope``)."""
    return np.unique(np.linspace(0, n - 1, min(512, n), dtype=int))


def write_gpu_cache(path, run, body, self_gravity, omega, offset, offset_halving,
                    offset_tolerance, evaluator_identity, evidence, seconds):
    """Write a CPU-format gravity cache for a GPU field; returns the writer's record.

    ``offset_halving`` is ``(max_deg, p95_deg, sample_count)`` from the GPU
    evaluation of the builder's sample at ``offset`` and ``offset/2``.
    """
    fit_cops = reference.load("fit_cops")
    gravity_cache = reference.load("gravity_cache")
    topology = reference.load("gravity_topology")
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal, area = elements["centroid"], elements["normal"], elements["area"]
    if self_gravity.shape != centroid.shape:
        raise ValueError("self-gravity shape %s does not match the facet table %s" %
                         (self_gravity.shape, centroid.shape))
    if body.source is None:
        raise ValueError("the gravity body must come from a file (shape identity)")
    definition = gravity_cache.make_definition(
        centroid, normal, area, body.source, len(body.vertices), len(body.faces),
        body.density, offset)
    max_deg, p95_deg, count = offset_halving
    if max_deg > offset_tolerance:
        raise ValueError("offset-halving slope change %.6g deg exceeds the tolerance %.6g deg" %
                         (max_deg, offset_tolerance))
    calculation = {
        "integrity_check": INTEGRITY_CHECK,
        "orientation_policy": topology.ORIENTATION_POLICY,
        "topology": body.reference_audit,
        "backend": BACKEND,
        "formulation": "Tsoulis homogeneous polyhedron; FP64 CUDA translation of ESA polyhedral-gravity 3.3",
        "evaluator": evaluator_identity,
        "evidence": evidence,
        "offset_halving_sample_count": int(count),
        "offset_halving_slope_difference_deg_max": float(max_deg),
        "offset_halving_slope_difference_deg_p95": float(p95_deg),
        "offset_tolerance_deg": float(offset_tolerance),
        "evaluation_seconds": float(seconds),
        "created_utc": provenance.utc_now(),
    }
    return gravity_cache.write_cache(path, centroid, normal, area, definition,
                                     np.ascontiguousarray(self_gravity, dtype=np.float64),
                                     np.asarray(omega, dtype=np.float64), calculation)
