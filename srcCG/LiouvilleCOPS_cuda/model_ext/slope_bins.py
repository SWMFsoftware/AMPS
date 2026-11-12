#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/slope_bins.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/slope_bins.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e96e05cf5d45e08a78355efc482ea6dac59905ccccd3944842a088a8332883a4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Builds the template bundle of the 7 effective-gravitational-slope bins ([0,40) ... [80,90), [90,180] deg)
#   whose amplitudes are the multipliers of the slope-multiplier model (user decisions 2026-11-10: "define multiple
#   cliff templates (0-40, ..., 80-90) ... activity multiplier"; answer: 7th bin 90-180), with checks (gravity-cache
#   sha256, mesh identity, slope >= 60 == cliff_ge_60deg, the lead's bin counts).
# USED BY: jobs/slope_mult_cpu.sbatch STEP=bins (CPU job 123624 ->
#   runs/shap8_012m_lategrid_dz1.5/slope_mult/bundle/slope_bins.npz, also used by the bump and nonneg models).
# ==== end of header ====
"""Template bundle of the slope bins for the slope-multiplier model (user decisions 2026-11-10; model_ext/slope_mult.py).

    .venv/bin/python -W ignore slope_bins.py --run <base> --bundle <restart template bundle> [--gravity-cache NPZ] \\
        --out <new .npz> [--expect-counts lead]

Writes an activity-template bundle in the reference format (activity_templates.write_bundle; loaded by fit_cops.py
--activity-template BUNDLE NAME) with the 7 bin indicators of the effective gravitational slope theta
(reference_slope_deg of the admitted gravity cache named by --bundle, or --gravity-cache):
    slope_00_40 = 1[theta < 40], slope_40_50 = 1[40 <= theta < 50], ..., slope_80_90 = 1[80 <= theta < 90],
    slope_90_180 = 1[90 <= theta <= 180]
(slope_mult.bin_index). The fitter fits one nonnegative amplitude per column; slope_mult.py's GPU patch reads these
amplitudes as the MULTIPLIERS m_b of A = max(1e-5 A_orig, m_b L) (they are not additive templates).

Checks (printed and stored in the metadata "slope_multiplier_bins"; the script stops on a failure):
* the gravity cache's SHA-256 is the one the source bundle records (as tools/c16_slope_template.py);
* the cache has one slope per facet of the base; slopes finite and in [0, 180];
* mesh identity: the base's elements (centroid, area) = the source bundle's mesh identity;
* (theta >= 60) == the source bundle's cliff_ge_60deg column > 0, facet by facet (the lead's identity check);
* the bin counts equal the lead's numbers (--expect-counts lead; default on);
* facets within 1e-6 deg of each inner bin edge (reported; the cliff template's threshold margin);
* the base's A_orig (AH2O): its SHA-256 (slope_mult checks it in every GPU process) and the area with A_orig = 0.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(APP)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import slope_mult as sm  # noqa: E402

EDGE_MARGIN_DEG = 1e-6


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", required=True)
    p.add_argument("--bundle", required=True)
    p.add_argument("--gravity-cache", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--expect-counts", default="lead", choices=("lead", "none"))
    args = p.parse_args(argv)
    from cops_cuda import figures as fg
    from cops_cuda import formats, provenance, reference
    reference.bind()
    templates = reference.load("activity_templates")
    run = Path(args.run).resolve()
    el, el_identity = formats.load_elements(run / "elements.bin", lazy=True)
    centroid, area = np.asarray(el["centroid"]), np.asarray(el["area"])
    ah2o = np.ascontiguousarray(np.asarray(el["AH2O"], dtype=np.float64))
    n = len(area)
    mesh = templates.mesh_identity(centroid, area)
    slope, cache = fg.gravity_slope(args.bundle, args.gravity_cache, nfaces=n)
    if slope is None:
        raise SystemExit("slope_bins: no gravity cache for %s" % args.bundle)
    meta_src = fg.bundle_metadata(args.bundle)
    gravity = dict(meta_src.get("gravity") or {})
    recorded = ((gravity.get("cache") or {}).get("file") or {}).get("sha256")
    used = provenance.sha256_file(cache)
    if recorded != used:
        raise SystemExit("slope_bins: gravity cache %s (sha256 %s) is not the one %s records (%s)"
                         % (cache, used, args.bundle, recorded))
    source = templates.load_bundle(args.bundle, centroid, area)       # checks the mesh identity
    if source["mesh_identity"] != mesh:
        raise SystemExit("slope_bins: mesh identity differs")
    cliff = templates.select_templates(source, ["cliff_ge_60deg"])[:, 0] > 0.0
    bins = sm.bin_index(slope)
    cliff_equal = bool(np.array_equal(slope >= 60.0, cliff))
    cliff_diff = int(np.sum((slope >= 60.0) != cliff))
    counts = np.bincount(bins, minlength=len(sm.BIN_NAMES))
    lead_equal = [int(v) for v in counts] == list(sm.LEAD_COUNTS)
    w = area / area.sum()
    edges = sm.BIN_EDGES_DEG[1:-1]
    margin = {"%g" % e: int(np.sum(np.abs(slope - e) <= EDGE_MARGIN_DEG)) for e in edges}
    zero = ah2o == 0.0
    per_bin = []
    for k, name in enumerate(sm.BIN_NAMES):
        sel = bins == k
        per_bin.append({"name": name, "interval_deg": sm.BIN_LABELS[k], "facets": int(sel.sum()),
                        "area_fraction": float(w[sel].sum()), "area_m2": float(area[sel].sum()),
                        "a_orig_zero_facets": int(np.sum(sel & zero)),
                        "a_orig_zero_area_fraction_of_bin": float(area[sel & zero].sum() / max(area[sel].sum(), 1e-300))})
    values = sm.indicator_basis(bins)
    sm_meta = {"names": list(sm.BIN_NAMES), "edges_deg": list(sm.BIN_EDGES_DEG), "labels": list(sm.BIN_LABELS),
               "convention": "[lo, hi) for every bin, the last one [90, 180] closed (slope_mult.bin_index)",
               "nfaces": int(n), "counts": [int(v) for v in counts], "lead_counts": list(sm.LEAD_COUNTS),
               "counts_equal_lead": bool(lead_equal), "per_bin": per_bin,
               "bins_sha256": sm.sha256_bytes(bins), "bins_dtype": "int8",
               "ah2o_sha256": sm.sha256_bytes(ah2o), "ah2o_source": "elements.bin column 32 (AH2O) of %s" % run,
               "a_orig_zero_facets": int(zero.sum()), "a_orig_zero_area_fraction": float(w[zero].sum()),
               "mesh_identity": mesh, "slope_source": {"path": str(cache), "sha256": used, "key": "reference_slope_deg",
                                                      "array_sha256": sm.sha256_bytes(np.asarray(slope, np.float64))},
               "cliff_ge_60deg_check": {"equal": cliff_equal, "facets_differing": cliff_diff,
                                        "cliff_facets": int(cliff.sum()), "source_bundle": str(Path(args.bundle).resolve())},
               "edge_margin_deg": EDGE_MARGIN_DEG, "facets_within_margin_of_edge": margin,
               "slope_min_max_deg": [float(slope.min()), float(slope.max())],
               "policy": sm.POLICY}
    print("slope_bins: %d facets; counts %s (lead %s: %s); (slope >= 60) == cliff_ge_60deg: %s (%d differ, cliff %d); "
          "within %g deg of an edge %s; A_orig = 0 on %d facets (%.4f of the area); slope %.6g .. %.6g deg"
          % (n, [int(v) for v in counts], list(sm.LEAD_COUNTS), "EQUAL" if lead_equal else "DIFFERENT", cliff_equal,
             cliff_diff, int(cliff.sum()), EDGE_MARGIN_DEG, margin, int(zero.sum()), float(w[zero].sum()),
             float(slope.min()), float(slope.max())), flush=True)
    for b in per_bin:
        print("slope_bins:   %-13s %-9s facets %8d  area fraction %.6f  A_orig = 0: %d facets" % (
            b["name"], b["interval_deg"], b["facets"], b["area_fraction"], b["a_orig_zero_facets"]), flush=True)
    if not cliff_equal:
        raise SystemExit("slope_bins: (slope >= 60) differs from cliff_ge_60deg on %d facets: not the bundle's slope"
                         % cliff_diff)
    if args.expect_counts == "lead" and not lead_equal:
        raise SystemExit("slope_bins: bin counts %s differ from the lead's %s" % (list(counts), list(sm.LEAD_COUNTS)))
    definitions = {}
    for k, name in enumerate(sm.BIN_NAMES):
        lo, hi = sm.BIN_EDGES_DEG[k], sm.BIN_EDGES_DEG[k + 1]
        rel = "<=" if k == len(sm.BIN_NAMES) - 1 else "<"
        definitions[name] = {
            "definition": "%g <= effective_gravitational_slope_deg %s %g" % (lo, rel, hi),
            "interpretation": "slope-bin indicator; its fitted amplitude is the MULTIPLIER m of A = max(1e-5 A_orig, "
                              "m L) on this bin (slope_mult.py), not an additive template",
            "face_count": int(counts[k]), "area_fraction": float(w[bins == k].sum()),
            "policy": "user 2026-11-10 (slope-multiplier model)"}
    metadata = {
        "created_utc": provenance.utc_now(),
        "base_run": str(run),
        "base_run_manifest": provenance.file_identity(run / "run_manifest.json", "run_manifest")
        if (run / "run_manifest.json").is_file() else None,
        "elements": provenance.file_identity(run / "elements.bin", "elements"),
        "gravity": gravity,
        "source_bundle": provenance.file_identity(args.bundle, "template_bundle"),
        "templates": definitions,
        "slope_multiplier_bins": sm_meta,
        "writer": provenance.file_identity(Path(__file__), "tool"),
        "writer_module": provenance.file_identity(Path(sm.__file__), "tool"),
    }
    payload = templates.write_bundle(args.out, list(sm.BIN_NAMES), values, centroid, area, metadata)
    Path(str(args.out)[:-4] + ".json").write_text(json.dumps(payload, indent=2) + "\n")
    print("slope_bins: wrote %s (sha256 %s)" % (args.out, hashlib.sha256(Path(args.out).read_bytes()).hexdigest()),
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
