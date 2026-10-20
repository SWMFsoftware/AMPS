"""Geometry-only visibility masks on the GPU (ROADMAP C05; CUDA environment).

``late_masks(run, start, stop)`` returns, for the observation points of a
closed UTC interval, the per-facet union over points of four same-point
predicates (``kernels/masks.cu``):

* ``ng_seen``, ``rg_seen``: diagnostics. Their counts are what the CPU builder
  records as ``ng_seen_face_count`` / ``rg_seen_face_count``.
* ``late_ng_not_rg``: Kn > 0 and KPhi == 0 and KV == 0 at the same point. This is
  the reference template of
  ``build_physical_activity_templates.late_gauge_template``.
* ``late_rg_not_ng``: (KPhi != 0 or KV != 0) and not Kn > 0 at the same point.
  These are the facets only RG sees (user decision 2026-10-19, same interval).

Point selection is the builder's: table rows with ``start <= UTC <= stop``
(``cops_selection.parse_utc``) and ``inside == 0``, in table order. The masks
are computed from the original FP64 kernels with integer OR, so they are
bit-identical to the CPU predicates; ``tools/c05_mask_check.py`` proves it on
the campaign data.

Response-weighted maps stay CPU diagnostics, written by the reference builder.
They never define a mask (ROADMAP §3.2).
"""

import numpy as np

from . import cuda_sources, reference
from .device_store import DeviceDataset

MASK_BITS = (("ng_seen", 1), ("rg_seen", 2), ("late_ng_not_rg", 4), ("late_rg_not_ng", 8))


def late_points(inputs, start, stop):
    """Point ids selected exactly like ``late_gauge_template`` (table order)."""
    selection = reference.load("cops_selection")
    t0, t1 = selection.parse_utc(start), selection.parse_utc(stop)
    cols = inputs.columns
    if "inside" not in cols:
        raise ValueError("table lacks the 'inside' column required for mask selection")
    points = [int(row[cols["iPoint"]]) for row in inputs.rows
              if t0 <= selection.parse_utc(row[cols["UTC"]]) <= t1 and not int(row[cols["inside"]])]
    if not points:
        raise ValueError("late interval selects no outside-nucleus observation rows")
    return points


def late_masks(run, start, stop, device_id=0, budget=None, expected_kernel_digests=None):
    """GPU masks over the interval's points; returns (masks dict, metadata)."""
    from .device_store import RunInputs
    inputs = RunInputs(run)
    points = late_points(inputs, start, stop)
    masks, meta = visibility_masks(run, points, device_id, budget, expected_kernel_digests)
    meta["interval_utc"] = [start, stop]
    return masks, meta


def visibility_masks(run, points, device_id=0, budget=None, expected_kernel_digests=None):
    """GPU masks (MASK_BITS) united over the given observation points; returns (masks, metadata)."""
    import cupy
    points = [int(p) for p in points]
    if not points:
        raise ValueError("no observation points selected")
    with DeviceDataset.load(run, point_ids=points, device_id=device_id, mesh_fields="none",
                            budget=budget, expected_kernel_digests=expected_kernel_digests) as ds:
        nf = ds.inputs.nfaces
        with cupy.cuda.Device(device_id):
            bits = cupy.zeros(nf, dtype=cupy.uint32)
            a = ds.arrays
            if ds.nrecords:
                cuda_sources.function("masks", "visibility_union", device_id)(
                    (max(1, min(65535, -(-ds.nrecords // 256))),), (256,),
                    (a["kernel.face"], a["kernel.kn"], a["kernel.kphi"], a["kernel.kv"],
                     np.int64(ds.nrecords), bits))
            host = bits.get()
            del bits
        identity = ds.identity()
    masks = {name: (host & bit) != 0 for name, bit in MASK_BITS}
    meta = {"selected_point_count": len(points),
            "selected_point_ids": points,
            "definitions": {
                "late_ng_not_rg": "union_over_points_of_Kn_gt_0_and_same_point_KPhi_eq_0_and_KV_eq_0",
                "late_rg_not_ng": "union_over_points_of_KPhi_ne_0_or_KV_ne_0_and_same_point_not_Kn_gt_0",
                "ng_seen": "union_over_points_of_Kn_gt_0",
                "rg_seen": "union_over_points_of_KPhi_ne_0_or_KV_ne_0"},
            "face_counts": {name: int(np.count_nonzero(m)) for name, m in masks.items()},
            "kernel_set_sha256": identity["kernels"]["loaded_set_sha256"],
            "mesh_identity": identity["mesh_identity"], "backend": "cuda masks.cu (integer atomicOr)",
            "cuda_source": cuda_sources.source_identity("masks")}
    return masks, meta


def seen_positive(dataset):
    """``fit_cops.Model.seen_faces`` on a resident dataset (host bool array per facet).

    Predicate (reference): facet ``face[r]`` is seen if ``Kn > 0 | KPhi > 0 | KV > 0``
    for any record r of the loaded points (``kernels/masks.cu: seen_positive``).
    """
    import cupy
    ds = dataset
    with cupy.cuda.Device(ds.device_id):
        seen = cupy.zeros(ds.inputs.nfaces, dtype=cupy.uint8)
        a = ds.arrays
        if ds.nrecords:
            cuda_sources.function("masks", "seen_positive", ds.device_id)(
                (max(1, min(65535, -(-ds.nrecords // 256))),), (256,),
                (a["kernel.face"], a["kernel.kn"], a["kernel.kphi"], a["kernel.kv"],
                 np.int64(ds.nrecords), seen))
        host = seen.get().astype(bool)
        del seen
    return host
