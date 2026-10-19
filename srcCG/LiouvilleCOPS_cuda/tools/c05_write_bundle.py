#!/usr/bin/env python3
"""Write the restart template bundle: NG-only, RG-only and cliff masks (C05; CUDA env).

User decision of 2026-10-19: a restart reuses persisted products, the admitted
gravity cache and a bundle of the NG/RG/cliff masks. This tool writes that
bundle with the CPU's own writer (``activity_templates.write_bundle``) and
maps (``build_physical_activity_templates.write_map``). The unchanged CPU
fitter can therefore select every column.

Columns:
* ``late_ng_not_rg`` and ``late_rg_not_ng``: GPU visibility masks
  (``cops_cuda.masks``), proven bit-identical to the CPU predicates
  (``tools/c05_mask_check.py``);
* ``cliff_ge_45deg`` and ``cliff_ge_60deg``: ``slope >= threshold`` from the admitted
  gravity cache, loaded with ``gravity_cache.load_cache`` against the definition
  rebuilt exactly as the builder builds it, and the builder's SPICE ω.

Threshold band (VALIDATION.md): a facet with ``|slope - threshold| <= 1e-6°``
is recomputed with the CPU reference package (DISABLE after the topology
audit) and takes the CPU classification. The re-resolved facets are listed in
the metadata.

The gravity section of the metadata mirrors the builder's, so the selector's
provenance gate applies unchanged: the cache definition (algorithm,
orientation policy) and the cache's calculation metadata (``integrity_check``,
``orientation_policy``, ``topology``, ``backend``, evidence).
"""

import argparse
import datetime as dt
import json
from pathlib import Path
import sys

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
BAND = 1e-6


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run", required=True)
    parser.add_argument("--gravity-cache", required=True, help="admitted cache (CPU format)")
    parser.add_argument("--gravity-shape", required=True)
    parser.add_argument("--out", required=True, help="new bundle .npz (must not exist)")
    parser.add_argument("--density", type=float, default=533.0)
    parser.add_argument("--offset", type=float, default=0.1)
    parser.add_argument("--late-start", default="2016-09-30T10:33:10")
    parser.add_argument("--late-stop", default="2016-09-30T10:39:26")
    parser.add_argument("--cliff-threshold", type=float, action="append", default=None)
    parser.add_argument("--device", type=int, default=0)
    args = parser.parse_args()
    thresholds = args.cliff_threshold or [45.0, 60.0]

    from cops_cuda import gravity_body, masks, provenance, reference
    builder = reference.load("build_physical_activity_templates")
    gravity_cache = reference.load("gravity_cache")
    templates = reference.load("activity_templates")
    fit_cops = reference.load("fit_cops")
    run = Path(args.run).resolve()
    out = Path(args.out).resolve()
    json_path, map_path = out.with_suffix(".json"), out.with_suffix(".png")
    for path in (out, json_path, map_path):
        if path.exists():
            raise SystemExit("refusing to replace %s" % path)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal, area = elements["centroid"], elements["normal"], elements["area"]

    # Visibility masks (GPU, integer OR).
    vis, vis_meta = masks.late_masks(run, args.late_start, args.late_stop, device_id=args.device)

    # Gravity: the builder's omega and definition, the admitted cache.
    body = gravity_body.GravityBody.from_bdf(args.gravity_shape, args.density)
    epochs = [args.late_start, builder.utc_midpoint(args.late_start, args.late_stop), args.late_stop]
    omega, variation, command = builder.spice_omega(
        (reference.bind() / "spice_angular_velocity").resolve(), run / "input.snapshot.in", epochs)
    definition = gravity_cache.make_definition(centroid, normal, area, body.source,
                                               len(body.vertices), len(body.faces),
                                               args.density, args.offset)
    cached = gravity_cache.load_cache(args.gravity_cache, centroid, normal, area, definition, omega)
    calc = cached["metadata"].get("calculation", {})
    if (calc.get("orientation_policy") != reference.load("gravity_topology").ORIENTATION_POLICY or
            calc.get("topology") != json.loads(json.dumps(body.reference_audit))):
        raise SystemExit("gravity cache lacks matching independently audited topology")
    slope = np.array(cached["slope"], dtype=np.float64)

    # Threshold band: re-resolve with the CPU reference package.
    band_faces = sorted(set(int(i) for th in thresholds
                            for i in np.flatnonzero(np.abs(slope - th) <= BAND)))
    resolved = {}
    if band_faces:
        import polyhedral_gravity as pg
        poly = pg.Polyhedron(polyhedral_source=(body.vertices, body.faces), density=args.density,
                             normal_orientation=pg.NormalOrientation.OUTWARDS,
                             integrity_check=pg.PolyhedronIntegrity.DISABLE,
                             metric_unit=pg.MetricUnit.METER)
        points = centroid[band_faces] + args.offset * normal[band_faces]
        result = pg.GravityEvaluable(polyhedron=poly)(computation_points=points.tolist(),
                                                      parallel=True)
        g_cpu = np.asarray([r[1] for r in result])
        s_cpu = gravity_cache.derive_effective_fields(centroid[band_faces], normal[band_faces],
                                                      g_cpu, omega)[2]
        for k, face in enumerate(band_faces):
            resolved[face] = float(s_cpu[k])
    names = ["late_ng_not_rg", "late_rg_not_ng"]
    columns = [vis["late_ng_not_rg"].astype(float), vis["late_rg_not_ng"].astype(float)]
    template_meta = {}
    for name, key in (("late_ng_not_rg", "late_ng_not_rg"), ("late_rg_not_ng", "late_rg_not_ng")):
        template_meta[name] = {
            "definition": vis_meta["definitions"][key], "interval_utc": vis_meta["interval_utc"],
            "selected_point_count": vis_meta["selected_point_count"],
            "selected_point_ids": vis_meta["selected_point_ids"],
            "face_count": vis_meta["face_counts"][key],
            "area_fraction": builder.area_fraction(vis[key], area),
            "computed_by": vis_meta["backend"], "kernel_set_sha256": vis_meta["kernel_set_sha256"],
            "interpretation": "geometry-defined candidate region; fitted amplitude is not proof of a physical source"}
    template_meta["late_ng_not_rg"]["ng_seen_face_count"] = vis_meta["face_counts"]["ng_seen"]
    template_meta["late_ng_not_rg"]["rg_seen_face_count"] = vis_meta["face_counts"]["rg_seen"]
    for th in thresholds:
        label = ("cliff_ge_%gdeg" % th).replace(".", "p")
        mask = slope >= th
        for face, s_cpu in resolved.items():
            if abs(slope[face] - th) <= BAND:
                mask[face] = s_cpu >= th
        names.append(label)
        columns.append(mask.astype(float))
        template_meta[label] = {
            "definition": "effective_gravitational_slope_deg >= %.17g" % th, "threshold_deg": th,
            "face_count": int(np.count_nonzero(mask)), "area_fraction": builder.area_fraction(mask, area),
            "threshold_margin_deg": BAND,
            "faces_within_threshold_margin_zero_based": [int(i) for i in np.flatnonzero(np.abs(slope - th) <= BAND)],
            "threshold_band_resolved_by_cpu_reference": {str(f): resolved[f] for f in resolved
                                                         if abs(slope[f] - th) <= BAND},
            "interpretation": "cliff/steep-terrain hypothesis; amplitude fitted nonnegative; no production law asserted"}
    gravity = cached["self_gravity"]
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_run": str(run), "base_run_manifest": provenance.file_identity(run / "run_manifest.json"),
        "elements": provenance.file_identity(run / "elements.bin"), "templates": template_meta,
        "gravity": {
            "shape": body.source, "shape_vertices": int(len(body.vertices)),
            "shape_triangles": int(len(body.faces)), "density_kg_m3": args.density,
            "homogeneous_density_assumption": True, "evaluation_offset_m": args.offset,
            "acceleration_definition": "self_gravity - omega_cross_omega_cross_r",
            "omega_body_frame_rad_s": [float(v) for v in omega], "omega_sample_epochs_utc": epochs,
            "omega_max_sample_deviation_rad_s": variation, "omega_command": [str(c) for c in command],
            "slope_definition": "angle_between_inward_facet_normal_and_effective_gravity",
            "slope_deg_min_median_max": [float(slope.min()), float(np.median(slope)), float(slope.max())],
            "gravity_m_s2_min_max": [float(np.min(np.linalg.norm(gravity, axis=1))),
                                     float(np.max(np.linalg.norm(gravity, axis=1)))],
            "cache": {"status": "reused", "file": provenance.file_identity(args.gravity_cache),
                      "definition": definition,
                      "reference_omega_reused": bool(cached["reference_omega_reused"])},
            **calc},
        "restart_product": {"decision": "user 2026-10-19: restart reuses gravity cache + NG/RG/cliff masks",
                            "writer": provenance.file_identity(Path(__file__), "tool")},
        "scientific_scope": {"species": "H2O",
                             "templates_are_predictors_not_established_causal_laws": True,
                             "late_ng_transition_requires_profile_sensitivity_test": True},
    }
    payload = templates.write_bundle(out, names, np.column_stack(columns), centroid, area, metadata)
    with json_path.open("x") as stream:
        stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    builder.write_map(map_path, centroid, area,
                      [("late NG-visible / RG-not-visible mask", columns[0], "viridis", (0.0, 1.0)),
                       ("late RG-visible / NG-not-visible mask", columns[1], "viridis", (0.0, 1.0)),
                       ("effective gravitational slope [deg]", slope, "magma", (0.0, 90.0))] +
                      [(n, v, "viridis", (0.0, 1.0)) for n, v in zip(names[2:], columns[2:])])
    # The unchanged fitter's selector must accept every column.
    bundle = templates.load_bundle(out, centroid, area)
    templates.select_templates(bundle, names)
    print("wrote %s, %s, %s; columns %s; band facets re-resolved %d" %
          (out, json_path, map_path, names, len(resolved)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
