#!/usr/bin/env python3
"""Continuous slope template: sin of the effective gravitational slope (user decision 2026-10-22).

    .venv/bin/python tools/c16_slope_template.py --run <base> --bundle <restart template bundle> --out <new.npz>

Writes an activity-template bundle in the reference format (``../LiouvilleCOPS/tools/activity_templates.py``
``write_bundle``; loaded by ``fit_cops.py --activity-template BUNDLE NAME``) with the column

    slope_sin = sin(min(theta, 90 deg))

where theta is the effective gravitational slope (self-gravity + centrifugal) of the gravity cache named by
``--bundle`` (or ``--gravity-cache``): the downslope share of the effective gravity, which drives mass
wasting. Overhangs (theta > 90 deg, where dust cannot stay either) count as fully steep. The template
replaces the binary ``cliff_ge_60deg`` mask with a continuous one; cos(solar zenith angle) is not applied
here, because the fitted activity already multiplies the production law of each observation's SZA.
Columns named by ``--copy`` (default ``late_ng_not_rg``) are copied from ``--bundle`` unchanged, so one
bundle serves the slope-only and the late-NG-only + slope fits.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))

from cops_cuda import figures as fg  # noqa: E402
from cops_cuda import formats, provenance, reference  # noqa: E402

DEFINITION = "sin(min(effective_gravitational_slope_deg, 90) * pi / 180)"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", required=True, help="the fitting base (elements.bin, run_manifest.json)")
    p.add_argument("--bundle", required=True, help="restart template bundle of this base (gravity cache, late templates)")
    p.add_argument("--gravity-cache", help="admitted gravity cache (default: the one named by --bundle)")
    p.add_argument("--copy", action="append", default=None, help="bundle columns to copy (default late_ng_not_rg)")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    reference.bind()
    templates = reference.load("activity_templates")
    run = Path(args.run)
    el, _ = formats.load_elements(run / "elements.bin", lazy=True)
    centroid, area = np.asarray(el["centroid"]), np.asarray(el["area"])
    slope, cache = fg.gravity_slope(args.bundle, args.gravity_cache, nfaces=len(area))
    if slope is None:
        raise SystemExit("no gravity cache found for %s" % args.bundle)
    source = templates.load_bundle(args.bundle, centroid, area)
    copy = args.copy or ["late_ng_not_rg"]
    columns = [templates.select_templates(source, copy)] if copy else []
    values_slope = np.sin(np.radians(np.minimum(slope, 90.0)))
    values = np.column_stack(columns + [values_slope]) if columns else values_slope[:, None]
    names = list(copy) + ["slope_sin"]
    meta_src = fg.bundle_metadata(args.bundle)
    gravity = dict(meta_src.get("gravity") or {})
    recorded = ((gravity.get("cache") or {}).get("file") or {}).get("sha256")
    used = provenance.sha256_file(cache)
    if recorded != used:
        raise SystemExit("gravity cache %s (sha256 %s) is not the one %s records (%s)" % (cache, used, args.bundle, recorded))
    w = area / area.sum()
    metadata = {
        "created_utc": provenance.utc_now(),
        "base_run": str(run.resolve()),
        "base_run_manifest": provenance.file_identity(run / "run_manifest.json", "run_manifest"),
        "elements": provenance.file_identity(run / "elements.bin", "elements"),
        # the source bundle's gravity provenance, verbatim: the reference loader checks its repaired
        # orientation policy for any template defined from the effective gravitational slope
        "gravity": gravity,
        "source_bundle": provenance.file_identity(args.bundle, "template_bundle"),
        "templates": dict({n: dict(meta_src.get("templates", {}).get(n, {}), copied_from=str(args.bundle))
                           for n in copy},
                          slope_sin={"definition": DEFINITION,
                                     "policy": "user 2026-10-22: continuous slope template sin(slope) replacing the "
                                               "binary cliff_ge_60deg mask; cos(SZA) enters through the production "
                                               "law, not the template",
                                     "overhangs": "theta > 90 deg counts as 90 deg (fully steep)",
                                     "area_weighted_mean": float(w @ values_slope),
                                     "area_fraction_above_0p5": float(w @ (values_slope > 0.5)),
                                     "overhang_facets": int(np.sum(slope > 90.0))}),
        "writer": provenance.file_identity(Path(__file__), "tool"),
    }
    payload = templates.write_bundle(args.out, names, values, centroid, area, metadata)
    Path(str(args.out)[:-4] + ".json").write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote %s: %s; slope_sin area-weighted mean %.3f, overhangs %d facets" % (
        args.out, names, metadata["templates"]["slope_sin"]["area_weighted_mean"],
        metadata["templates"]["slope_sin"]["overhang_facets"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
