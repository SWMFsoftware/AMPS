#!/usr/bin/env python3
"""Write the restart template bundle: NG-only, RG-only and cliff masks (C05; CUDA env).

User decision of 2026-10-19: a restart reuses persisted products, the admitted
gravity cache and a bundle of the NG/RG/cliff masks. This tool is the glue: it
computes the GPU masks, obtains the builder's SPICE omega and loads the cache;
``cops_cuda.restart_bundle.build_bundle`` then writes the bundle with the CPU's
own writer (``activity_templates.write_bundle``) and maps
(``build_physical_activity_templates.write_map``). The unchanged CPU fitter can
therefore select every column.

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
from pathlib import Path
import sys

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))


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

    from cops_cuda import gravity_body, masks, provenance, reference, restart_bundle
    builder = reference.load("build_physical_activity_templates")
    gravity_cache = reference.load("gravity_cache")
    fit_cops = reference.load("fit_cops")
    run = Path(args.run).resolve()
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

    def resolve_band(faces):
        """Threshold band: the CPU reference package, DISABLE after the topology audit."""
        import polyhedral_gravity as pg
        poly = pg.Polyhedron(polyhedral_source=(body.vertices, body.faces), density=args.density,
                             normal_orientation=pg.NormalOrientation.OUTWARDS,
                             integrity_check=pg.PolyhedronIntegrity.DISABLE,
                             metric_unit=pg.MetricUnit.METER)
        points = centroid[faces] + args.offset * normal[faces]
        result = pg.GravityEvaluable(polyhedron=poly)(computation_points=points.tolist(),
                                                      parallel=True)
        g_cpu = np.asarray([r[1] for r in result])
        return gravity_cache.derive_effective_fields(centroid[faces], normal[faces], g_cpu, omega)[2]

    payload = restart_bundle.build_bundle(
        args.out, run, body, cached, definition, args.gravity_cache,
        (omega, epochs, variation, command), vis, vis_meta, thresholds, resolve_band,
        writer=provenance.file_identity(Path(__file__), "tool"))
    out = Path(args.out).resolve()
    resolved = sum(len(m.get("threshold_band_resolved_by_cpu_reference", {}))
                   for m in payload["templates"].values())
    print("wrote %s, %s, %s; columns %s; band facets re-resolved %d" %
          (out, out.with_suffix(".json"), out.with_suffix(".png"), payload["names"], resolved))
    return 0


if __name__ == "__main__":
    sys.exit(main())
