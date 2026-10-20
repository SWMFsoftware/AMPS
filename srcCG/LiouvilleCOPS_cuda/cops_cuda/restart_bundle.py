"""Restart template bundle: NG-only, RG-only and cliff masks (ROADMAP C05 item 4).

User decision of 2026-10-19: a restart reuses persisted products, namely the
admitted gravity cache and a bundle of the NG/RG/cliff masks. ``build_bundle``
writes that bundle with the CPU's own writer (``activity_templates.write_bundle``)
and checks that the unchanged fitter's ``load_bundle``/``select_templates``
accept every column. ``tools/c05_write_bundle.py`` is the command-line glue
(GPU masks, SPICE omega, cache loading). This module holds no I/O beyond the
bundle itself, so it is testable on synthetic runs.

Columns
-------
* ``late_ng_not_rg`` and ``late_rg_not_ng``: the GPU visibility masks
  (``masks.late_masks``), proven bit-identical to the CPU predicates;
* ``cliff_ge_<th>deg``: ``slope >= th`` from the admitted cache's slopes.

Threshold band (VALIDATION.md): facets with ``|slope - th| <= BAND`` are
re-resolved by ``resolve_band(faces) -> slopes``. The tool passes the CPU
reference package for this. Those facets take the resolver's classification,
and the metadata lists every re-resolved facet.

The gravity section of the metadata mirrors the builder's
(``build_physical_activity_templates``): the cache definition (``algorithm``,
``orientation_policy``) and the cache's calculation metadata
(``integrity_check``, ``topology``, ``backend``, evidence). The selector's
provenance gate therefore applies unchanged.
"""

import datetime as dt
import json
from pathlib import Path

import numpy as np

from . import provenance, reference

BAND = 1e-6
VISIBILITY = ("late_ng_not_rg", "late_rg_not_ng")


def cliff_name(threshold):
    """The builder's column name: ``cliff_ge_45deg``, ``cliff_ge_52p5deg``."""
    return ("cliff_ge_%gdeg" % threshold).replace(".", "p")


def cliff_columns(slope, thresholds, resolve_band=None):
    """``{threshold: mask}`` plus ``{face: resolved slope}`` for the band facets."""
    slope = np.asarray(slope, dtype=np.float64)
    band = sorted(set(int(i) for th in thresholds
                      for i in np.flatnonzero(np.abs(slope - th) <= BAND)))
    resolved = {}
    if band:
        if resolve_band is None:
            raise ValueError("%d facets lie within %.1g deg of a threshold and no "
                             "resolver was given" % (len(band), BAND))
        values = np.asarray(resolve_band(band), dtype=np.float64)
        if values.shape != (len(band),) or not np.isfinite(values).all():
            raise ValueError("band resolver returned %s for %d facets" % (values.shape, len(band)))
        resolved = {face: float(v) for face, v in zip(band, values)}
    masks = {}
    for th in thresholds:
        mask = slope >= th
        for face, value in resolved.items():
            if abs(slope[face] - th) <= BAND:
                mask[face] = value >= th
        masks[th] = mask
    return masks, resolved


def build_bundle(out, run, body, cached, definition, gravity_cache_path, omega_info,
                 vis, vis_meta, thresholds=(45.0, 60.0), resolve_band=None,
                 writer=None, write_map=True):
    """Write ``out`` (.npz) and its ``.json`` (and ``.png``); returns the payload.

    ``cached`` is ``gravity_cache.load_cache``'s result for ``definition`` (built
    with ``gravity_cache.make_definition``, which fixes density and offset);
    ``omega_info`` is ``(omega, epochs, variation, command)`` as the builder
    obtains it; ``vis``/``vis_meta`` come from ``masks.late_masks``.
    """
    builder = reference.load("build_physical_activity_templates")
    templates = reference.load("activity_templates")
    topology = reference.load("gravity_topology")
    fit_cops = reference.load("fit_cops")
    out = Path(out).resolve()
    run = Path(run).resolve()
    json_path, map_path = out.with_suffix(".json"), out.with_suffix(".png")
    for path in (out, json_path, map_path):
        if path.exists():
            raise ValueError("refusing to replace %s" % path)
    calc = cached["metadata"].get("calculation", {})
    # The builder's own reuse conditions (CPU builder, gravity cache reuse).
    if (calc.get("orientation_policy") != topology.ORIENTATION_POLICY or
            calc.get("topology") != json.loads(json.dumps(body.reference_audit))):
        raise ValueError("gravity cache lacks matching independently audited topology")
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, area = elements["centroid"], elements["area"]
    slope = np.asarray(cached["slope"], dtype=np.float64)
    if slope.shape != (len(area),):
        raise ValueError("cache slopes do not match the mesh")
    omega, epochs, variation, command = omega_info
    # Density and offset are those the cache was defined (and loaded) for.
    density, offset = definition["density_kg_m3"], definition["evaluation_offset_m"]

    names = list(VISIBILITY)
    columns = [np.asarray(vis[name], dtype=bool).astype(float) for name in VISIBILITY]
    template_meta = {}
    for name in VISIBILITY:
        template_meta[name] = {
            "definition": vis_meta["definitions"][name], "interval_utc": vis_meta["interval_utc"],
            "selected_point_count": vis_meta["selected_point_count"],
            "selected_point_ids": vis_meta["selected_point_ids"],
            "face_count": vis_meta["face_counts"][name],
            "area_fraction": builder.area_fraction(np.asarray(vis[name], dtype=bool), area),
            "computed_by": vis_meta["backend"], "kernel_set_sha256": vis_meta["kernel_set_sha256"],
            "interpretation": "geometry-defined candidate region; fitted amplitude is not proof "
                              "of a physical source"}
    template_meta["late_ng_not_rg"]["ng_seen_face_count"] = vis_meta["face_counts"]["ng_seen"]
    template_meta["late_ng_not_rg"]["rg_seen_face_count"] = vis_meta["face_counts"]["rg_seen"]

    masks, resolved = cliff_columns(slope, thresholds, resolve_band)
    for th in thresholds:
        label = cliff_name(th)
        names.append(label)
        columns.append(masks[th].astype(float))
        in_band = [int(i) for i in np.flatnonzero(np.abs(slope - th) <= BAND)]
        template_meta[label] = {
            "definition": "effective_gravitational_slope_deg >= %.17g" % th, "threshold_deg": th,
            "face_count": int(np.count_nonzero(masks[th])),
            "area_fraction": builder.area_fraction(masks[th], area),
            "threshold_margin_deg": BAND,
            "faces_within_threshold_margin_zero_based": in_band,
            "threshold_band_resolved_by_cpu_reference": {str(f): resolved[f] for f in in_band},
            "interpretation": "cliff/steep-terrain hypothesis; amplitude fitted nonnegative; "
                              "no production law asserted"}
    gravity = cached["self_gravity"]
    norm = np.linalg.norm(gravity, axis=1)
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_run": str(run),
        "base_run_manifest": provenance.file_identity(run / "run_manifest.json"),
        "elements": provenance.file_identity(run / "elements.bin"), "templates": template_meta,
        "gravity": {
            "shape": body.source, "shape_vertices": int(len(body.vertices)),
            "shape_triangles": int(len(body.faces)), "density_kg_m3": density,
            "homogeneous_density_assumption": True, "evaluation_offset_m": offset,
            "acceleration_definition": "self_gravity - omega_cross_omega_cross_r",
            "omega_body_frame_rad_s": [float(v) for v in omega], "omega_sample_epochs_utc": epochs,
            "omega_max_sample_deviation_rad_s": variation,
            "omega_command": [str(c) for c in command],
            "slope_definition": "angle_between_inward_facet_normal_and_effective_gravity",
            "slope_deg_min_median_max": [float(slope.min()), float(np.median(slope)),
                                         float(slope.max())],
            "gravity_m_s2_min_max": [float(norm.min()), float(norm.max())],
            "cache": {"status": "reused", "file": provenance.file_identity(gravity_cache_path),
                      "definition": definition,
                      "reference_omega_reused": bool(cached["reference_omega_reused"])},
            **calc},
        "restart_product": {"decision": "user 2026-10-19: restart reuses gravity cache + "
                                        "NG/RG/cliff masks",
                            "writer": writer or provenance.file_identity(Path(__file__), "module")},
        "scientific_scope": {"species": "H2O",
                             "templates_are_predictors_not_established_causal_laws": True,
                             "late_ng_transition_requires_profile_sensitivity_test": True},
    }
    payload = templates.write_bundle(out, names, np.column_stack(columns), centroid, area, metadata)
    with json_path.open("x") as stream:
        stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if write_map:
        builder.write_map(map_path, centroid, area,
                          [("late NG-visible / RG-not-visible mask", columns[0], "viridis", (0.0, 1.0)),
                           ("late RG-visible / NG-not-visible mask", columns[1], "viridis", (0.0, 1.0)),
                           ("effective gravitational slope [deg]", slope, "magma", (0.0, 90.0))] +
                          [(n, v, "viridis", (0.0, 1.0)) for n, v in zip(names[2:], columns[2:])])
    # The unchanged fitter's loader and selector must accept every column.
    bundle = templates.load_bundle(out, centroid, area)
    templates.select_templates(bundle, names)
    return payload


def build_visibility_bundle(out, run, vis, vis_meta, writer=None):
    """Bundle with only the two visibility columns (no gravity needed); returns the payload.

    Used to start NG-only/RG-only fits before a gravity cache exists. The full restart
    bundle (with the cliff columns) is written later by ``build_bundle``. Both bundles
    carry identical visibility columns.
    """
    builder = reference.load("build_physical_activity_templates")
    templates = reference.load("activity_templates")
    fit_cops = reference.load("fit_cops")
    out, run = Path(out).resolve(), Path(run).resolve()
    json_path = out.with_suffix(".json")
    for path in (out, json_path):
        if path.exists():
            raise ValueError("refusing to replace %s" % path)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, area = elements["centroid"], elements["area"]
    names = list(VISIBILITY)
    columns = [np.asarray(vis[n], dtype=bool).astype(float) for n in VISIBILITY]
    template_meta = {}
    for name in VISIBILITY:
        template_meta[name] = {
            "definition": vis_meta["definitions"][name], "interval_utc": vis_meta["interval_utc"],
            "selected_point_count": vis_meta["selected_point_count"],
            "selected_point_ids": vis_meta["selected_point_ids"],
            "face_count": vis_meta["face_counts"][name],
            "area_fraction": builder.area_fraction(np.asarray(vis[name], dtype=bool), area),
            "computed_by": vis_meta["backend"], "kernel_set_sha256": vis_meta["kernel_set_sha256"],
            "interpretation": "geometry-defined candidate region; fitted amplitude is not proof "
                              "of a physical source"}
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_run": str(run), "base_run_manifest": provenance.file_identity(run / "run_manifest.json"),
        "elements": provenance.file_identity(run / "elements.bin"), "templates": template_meta,
        "restart_product": {"decision": "visibility-only bundle (no gravity); cliff columns follow in the "
                                        "full restart bundle",
                            "writer": writer or provenance.file_identity(Path(__file__), "module")},
    }
    payload = templates.write_bundle(out, names, np.column_stack(columns), centroid, area, metadata)
    with json_path.open("x") as stream:
        stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    bundle = templates.load_bundle(out, centroid, area)
    templates.select_templates(bundle, names)
    return payload

