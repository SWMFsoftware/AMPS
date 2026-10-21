#!/usr/bin/env python3
"""Full-field CPU-vs-GPU gravity equivalence (C04 acceptance; CUDA env).

Decision of 2026-10-19 (ROADMAP C04 item 2). The gravity model is tested at full-field
scale by evaluating the same homogeneous body at the same targets with both
implementations:

* **CPU:** the repaired reference builder
  (``build_physical_activity_templates.py``, DISABLE after the topology audit),
  which writes ``<out>.gravity.npz``;
* **GPU:** ``cops_cuda.gravity`` (``fold`` variant, the package's summation order).

The body is the derived outward surface of ``SHAP5_Volume.bdf`` (100,000
faces; see ``tools/derive_outward_surface.py``). The targets are all 787,754
facets of the 788k base run at ``centroid + 0.1 m * normal``.

Subcommands:

``gpu``
    Evaluate the GPU field and save ``gpu_field.npz`` and ``gpu_manifest.json`` in a new
    ``$COPS_OUTPUT_ROOT/cuda_c04/equivalence_<UTC>/``. The field contains the
    self-gravity, term magnitude M, potential, and the 512-sample half-offset field of
    the offset-halving audit. ``--devices 0,1`` splits the targets into contiguous
    blocks, one per GPU, in spawned processes. Each target is a fixed-order sum over all
    body faces, so the field is bitwise the one-device field
    (``tests/test_gravity_split.py``). The manifest records the partition.
``compare``
    Load the CPU cache through the reference ``gravity_cache.load_cache``, against the
    definition rebuilt exactly as the builder builds it, so a cache for another body,
    offset, density or facet order is refused. Then compare:

    * per facet: ``|dg| <= 16 eps M`` and ``<= 1e-12`` relative;
    * slopes: within 1e-9°, using the cache's reference omega;
    * 45°/60° cliff masks: identical after CPU re-resolution of the 1e-6° band;
    * the CPU's offset-halving audit value;
    * optionally the template bundle's cliff columns.

    With ``--gpu-788k-field`` (the C03 ``full_field.npz``) it also reports the corrected
    100k-vs-788k slope comparison that replaces the retracted HEAL-era "rejection".

Example::

    .venv/bin/python tools/c04_gravity_equivalence.py gpu \\
        --run $COPS_OUTPUT_ROOT/output_fit_shap5_788k_base \\
        --gravity-shape .../shape_models/derived/SHAP5_Volume_surface_outward.bdf
    .venv/bin/python tools/c04_gravity_equivalence.py compare --dir <equivalence dir> \\
        --cpu-cache $COPS_OUTPUT_ROOT/<bundle>.gravity.npz --cpu-bundle $COPS_OUTPUT_ROOT/<bundle>.npz \\
        --gpu-788k-field $COPS_OUTPUT_ROOT/cuda_c03/gravity_20261019T171225Z/full_field.npz
"""

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import gravity, gravity_body, manifest, provenance, reference  # noqa: E402

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"
EPS = float(np.finfo(float).eps)
M_BUDGET, REL_BUDGET, SLOPE_BUDGET, BAND = 16.0 * EPS, 1e-12, 1e-9, 1e-6
THRESHOLDS = (45.0, 60.0)


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def offset_sample(n):
    """The builder's deterministic offset-halving sample (effective_slope)."""
    return np.unique(np.linspace(0, n - 1, min(512, n), dtype=int))


def _split_worker(shape, density, device, points, part_path, queue, precision="fp64", near_factor=6.0,
                  allow_genus=0):
    """One device's contiguous block of targets (spawned process; own CUDA context)."""
    try:
        body = gravity_body.GravityBody.from_bdf(shape, density, allow_genus)
        last = [0.0]

        def progress(done, total, seconds):
            if seconds - last[0] >= 300.0 or done == total:
                last[0] = seconds
                log("  device %d: %d/%d (%.1f%%) %.0f s" % (device, done, total, 100.0 * done / total, seconds))
        with gravity.GravityEvaluator(body, device_id=device, precision=precision, near_factor=near_factor) as ev:
            potential, g, M = ev.evaluate(points, variant="fold", progress=progress, with_magnitude=True)
            near = ev.near_counts if ev.near_counts is not None else np.zeros(0, dtype=np.int64)
        np.savez(part_path, potential=potential, g=g, M=M, near=near)
        queue.put((device, None))
    except BaseException as exc:                          # reported to the parent
        queue.put((device, "%s: %s" % (exc.__class__.__name__, exc)))


def evaluate_split(args, points, devices, out):
    """Targets split into contiguous blocks, one per device, evaluated in parallel.

    Each target's field is a fixed-order sum over all body faces that does not depend on the
    other targets, so the concatenated result is bitwise the single-device result
    (tests/test_gravity_split.py). Returns potential, g, M and the partition record.
    """
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    bounds = np.linspace(0, len(points), len(devices) + 1).astype(int)
    queue = ctx.Queue()
    procs, parts = [], []
    for k, d in enumerate(devices):
        part = out / ("part_%d.npz" % k)
        parts.append(part)
        p = ctx.Process(target=_split_worker, args=(args.gravity_shape, args.density, d,
                                                    points[bounds[k]:bounds[k + 1]], str(part), queue,
                                                    getattr(args, "precision", "fp64"),
                                                    getattr(args, "near_factor", 6.0),
                                                    getattr(args, "allow_genus", 0)))
        p.start()
        procs.append(p)
    errors = [queue.get() for _ in procs]
    for p in procs:
        p.join()
    failed = [e for e in errors if e[1] is not None]
    if failed:
        raise RuntimeError("gravity split failed: %s" % failed)
    blocks = [np.load(str(part)) for part in parts]
    potential = np.concatenate([b["potential"] for b in blocks])
    g = np.concatenate([b["g"] for b in blocks])
    M = np.concatenate([b["M"] for b in blocks])
    near = np.concatenate([b["near"] for b in blocks])
    for part in parts:
        part.unlink()
    return potential, g, M, {"devices": devices, "bounds": bounds.tolist(),
                             "near_faces_per_target": ({"mean": float(near.mean()), "max": int(near.max()),
                                                        "min": int(near.min())} if near.size else None)}


def stage_gpu(args):
    fit_cops = reference.load("fit_cops")
    run = Path(args.run).resolve()
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal = elements["centroid"], elements["normal"]
    body = gravity_body.GravityBody.from_bdf(args.gravity_shape, args.density, args.allow_genus)
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c04" /
                                        ("equivalence_" + provenance.utc_stamp()))
    log("directory %s; body %d faces; %d targets" % (out, body.nfaces, len(centroid)))
    # The builder's exact target expression (effective_slope).
    points = centroid + args.offset * normal
    sample = offset_sample(len(centroid))
    half = centroid[sample] + 0.5 * args.offset * normal[sample]
    last = [0.0]

    def progress(done, total, seconds):
        if seconds - last[0] >= 60.0 or done == total:
            last[0] = seconds
            log("  %d/%d (%.1f%%) %.0f s" % (done, total, 100.0 * done / total, seconds))
    t0 = time.perf_counter()
    devices = [int(d) for d in args.devices.split(",")] if args.devices else [args.device]
    partition = None
    options = {"precision": args.precision, "near_factor": args.near_factor}
    if len(devices) == 1:
        with gravity.GravityEvaluator(body, device_id=devices[0], **options) as ev:
            potential, g, M = ev.evaluate(points, variant="fold", progress=progress,
                                          with_magnitude=True)
            identity = ev.identity()             # near-face statistics of the full target set
            _, g_half, M_half = ev.evaluate(half, variant="fold", with_magnitude=True)
    else:
        potential, g, M, partition = evaluate_split(args, points, devices, out)
        with gravity.GravityEvaluator(body, device_id=devices[0], **options) as ev:
            _, g_half, M_half = ev.evaluate(half, variant="fold", with_magnitude=True)
            identity = ev.identity()
    seconds = time.perf_counter() - t0
    np.savez_compressed(out / "gpu_field.npz", self_gravity_m_s2=g, term_magnitude_m_s2=M,
                        potential_m2_s2=potential, offset_sample=sample,
                        half_offset_self_gravity_m_s2=g_half,
                        half_offset_term_magnitude_m_s2=M_half)
    record = {"schema": "liouville_cops_cuda.c04_gpu_field", "schema_version": 1,
              "created_utc": provenance.utc_now(), "run": str(run),
              "elements": provenance.file_identity(run / "elements.bin", "elements"),
              "gravity_shape": provenance.file_identity(args.gravity_shape, "gravity_shape"),
              "density_kg_m3": args.density, "offset_m": args.offset,
              "targets": int(len(points)), "seconds": seconds,
              "pairs_per_second": float(len(points)) * body.nfaces / seconds,
              "evaluator": identity,
              "field": provenance.file_identity(out / "gpu_field.npz", "gpu_field"),
              "devices": devices, "target_partition": partition,
              "manifest": manifest.backend_manifest("cuda_polyhedral_gravity_v1", devices[0])}
    provenance.write_json_exclusive(out / "gpu_manifest.json", record)
    log("GPU field done in %.0f s (%.3g pairs/s): %s" % (seconds, record["pairs_per_second"], out))
    return 0


def stage_compare(args):
    fit_cops = reference.load("fit_cops")
    gravity_cache = reference.load("gravity_cache")
    directory = Path(args.dir).resolve()
    gpu_record = json.loads((directory / "gpu_manifest.json").read_text())
    run = Path(gpu_record["run"])
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal, area = elements["centroid"], elements["normal"], elements["area"]
    gpu = np.load(str(directory / "gpu_field.npz"))
    shape_path = gpu_record["gravity_shape"]["path"]
    builder = reference.load("build_physical_activity_templates")
    V, F = builder.read_nastran_surface(shape_path)

    # Rebuild the definition exactly as the builder does; load_cache refuses a
    # cache made for another body, offset, density or facet order.
    definition = gravity_cache.make_definition(
        centroid, normal, area, provenance.file_identity(shape_path), len(V), len(F),
        gpu_record["density_kg_m3"], gpu_record["offset_m"])
    with np.load(str(args.cpu_cache), allow_pickle=False) as raw:
        omega = np.asarray(raw["reference_omega_body_rad_s"], dtype=float)
    cpu = gravity_cache.load_cache(args.cpu_cache, centroid, normal, area, definition, omega)
    g_cpu, slope_cpu = cpu["self_gravity"], cpu["slope"]
    calc = cpu["metadata"].get("calculation", {})
    g, M = gpu["self_gravity_m_s2"], gpu["term_magnitude_m_s2"]
    log("loaded CPU cache %s (%s)" % (args.cpu_cache, calc.get("orientation_policy")))

    diff = np.linalg.norm(g - g_cpu, axis=1)
    rel = diff / np.linalg.norm(g_cpu, axis=1)
    per_m = diff / M
    slope = gravity_cache.derive_effective_fields(centroid, normal, g, omega)[2]
    dslope = np.abs(slope - slope_cpu)
    masks = {}
    for th in THRESHOLDS:
        band = np.abs(slope_cpu - th) <= BAND
        gpu_mask = slope >= th
        # Policy: facets inside the band take the CPU classification.
        resolved = np.where(band, slope_cpu >= th, gpu_mask)
        masks["%g" % th] = {
            "band_facets": int(band.sum()),
            "disagreements_before_band_resolution": int(np.count_nonzero(gpu_mask != (slope_cpu >= th))),
            "identical_after_band_resolution": bool(np.array_equal(resolved, slope_cpu >= th)),
            "area_fraction": float(np.sum(area[slope_cpu >= th]) / np.sum(area)),
            "face_count": int(np.count_nonzero(slope_cpu >= th))}
    # Offset-halving audit through the GPU, compared with the CPU's recorded value.
    sample = gpu["offset_sample"]
    s_half = gravity_cache.derive_effective_fields(
        centroid[sample], normal[sample], gpu["half_offset_self_gravity_m_s2"], omega)[2]
    halving_gpu = float(np.max(np.abs(slope[sample] - s_half)))
    halving_cpu = calc.get("offset_halving_slope_difference_deg_max")
    checks = {
        "all_finite": bool(np.isfinite(g).all() and np.isfinite(g_cpu).all()),
        "acceleration": {"rel_max": float(rel.max()), "rel_median": float(np.median(rel)),
                         "per_M_max_eps": float(per_m.max() / EPS),
                         "cancellation_M_over_g_max": float(np.max(M / np.linalg.norm(g_cpu, axis=1))),
                         "bitwise_facets": int(np.count_nonzero(np.all(g == g_cpu, axis=1)))},
        "slope": {"abs_max_deg": float(dslope.max()), "abs_median_deg": float(np.median(dslope))},
        "cliff_masks": masks,
        "offset_halving": {"gpu_max_deg": halving_gpu, "cpu_recorded_max_deg": halving_cpu}}
    ok = (checks["all_finite"] and rel.max() <= REL_BUDGET and per_m.max() <= M_BUDGET and
          dslope.max() <= SLOPE_BUDGET and
          all(m["identical_after_band_resolution"] for m in masks.values()) and
          halving_cpu is not None and abs(halving_gpu - float(halving_cpu)) <= SLOPE_BUDGET)

    if args.cpu_bundle:
        templates = reference.load("activity_templates")
        bundle = templates.load_bundle(args.cpu_bundle, centroid, area)
        cols = {}
        for th in THRESHOLDS:
            name = "cliff_ge_%gdeg" % th
            if name in bundle["names"]:
                column = templates.select_templates(bundle, [name])[:, 0] > 0.0
                cols[name] = bool(np.array_equal(column, slope_cpu >= th))
        checks["cpu_bundle_cliff_columns_match_cache_slopes"] = cols
        ok = ok and all(cols.values())

    if args.gpu_788k_field:
        # By-product: corrected resolution comparison (100k body vs exact 788k body).
        with np.load(str(args.gpu_788k_field), allow_pickle=False) as full:
            slope_788 = np.asarray(full["slope_deg"])
        d = np.abs(slope_cpu - slope_788)
        checks["resolution_100k_vs_788k"] = {
            "note": "same targets; gravity body 100k (CPU) vs exact 788k (GPU, C03); physical, not numerical",
            "slope_difference_deg": {"median": float(np.median(d)),
                                     "p95": float(np.percentile(d, 95.0)), "max": float(d.max())},
            "historical_gate_p95_deg": 5.0,
            "cliff_area_fraction_100k": {k: v["area_fraction"] for k, v in masks.items()},
            "cliff_area_fraction_788k": {"%g" % th: float(np.sum(area[slope_788 >= th]) / np.sum(area))
                                         for th in THRESHOLDS},
            "classification_disagreement_area_fraction": {
                "%g" % th: float(np.sum(area[(slope_cpu >= th) != (slope_788 >= th)]) / np.sum(area))
                for th in THRESHOLDS}}

    report = {"schema": "liouville_cops_cuda.c04_equivalence", "schema_version": 1,
              "created_utc": provenance.utc_now(), "gpu_record": gpu_record["field"],
              "cpu_cache": provenance.file_identity(args.cpu_cache, "cpu_gravity_cache"),
              "cpu_calculation": {k: calc.get(k) for k in ("orientation_policy", "package_version",
                                                           "evaluation_seconds", "parallelization")},
              "budgets": {"per_M": "16 eps M", "rel": REL_BUDGET, "slope_deg": SLOPE_BUDGET,
                          "band_deg": BAND},
              "checks": checks, "pass": bool(ok)}
    provenance.write_json_exclusive(directory / ("compare_%s.json" % provenance.utc_stamp()), report)
    log("compare %s: rel max %.3g, %.3f eps*M, slope max %.3g deg, masks %s" %
        ("PASS" if ok else "FAIL", rel.max(), per_m.max() / EPS, dslope.max(),
         {k: v["identical_after_band_resolution"] for k, v in masks.items()}))
    return 0 if ok else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="stage", required=True)
    g = sub.add_parser("gpu")
    g.add_argument("--run", required=True)
    g.add_argument("--gravity-shape", required=True)
    g.add_argument("--density", type=float, default=533.0)
    g.add_argument("--offset", type=float, default=0.1)
    g.add_argument("--device", type=int, default=0)
    g.add_argument("--devices", default=None,
                   help="comma-separated devices: targets split into contiguous blocks (bitwise the one-device field)")
    g.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    g.add_argument("--precision", choices=("fp64", "mixed"), default="fp64",
                   help="fp64 (default): the exact reference; mixed (experiment, 2026-10-21): exact FP64 "
                        "terms of faces within --near-factor circumradii, FP32 quadrature for all others")
    g.add_argument("--near-factor", type=float, default=6.0, help="mixed precision near radius in circumradii")
    g.add_argument("--allow-genus", type=int, default=0,
                   help="admit a closed oriented body of genus up to this value (experiment; default 0)")
    c = sub.add_parser("compare")
    c.add_argument("--dir", required=True)
    c.add_argument("--cpu-cache", required=True)
    c.add_argument("--cpu-bundle")
    c.add_argument("--gpu-788k-field")
    args = parser.parse_args()
    return stage_gpu(args) if args.stage == "gpu" else stage_compare(args)


if __name__ == "__main__":
    sys.exit(main())
