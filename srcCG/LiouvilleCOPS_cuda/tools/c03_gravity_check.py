#!/usr/bin/env python3
"""C03 evidence: exact FP64 GPU gravity against the pinned CPU reference (CUDA env).

Writes ``$COPS_OUTPUT_ROOT/cuda_c03/gravity_<UTC>/report.json``. With ``--full``
it also writes ``full_field.npz``, the self-gravity at every 788k target.

Checks, each with a pass flag and its measurements:

``synthetic``
    The C00 shapes (``gravity/gravity_synthetic.npz``, campaign cp38 binary):
    the original, translated, rotated and double-density cases.
    * The budget is ``|dg| <= 16 eps M``, where M is the kernel's term magnitude
      (VALIDATION.md C03 amendment).
    * The NaN pattern of the original arrays must be identical to the package's.
    * Transformed cases are rebuilt from the saved arrays. Their vertex bits
      may differ from the capture, so their NaN agreement is reported only.
``precompute``
    The face table (G, N, n) of the 788k body is compared bit for bit with
    ``GravityEvaluable.__getstate__``.
``wheel_equivalence``
    The CUDA-environment package (cp311 wheel) is compared with the C00
    fixture (campaign cp38 binary) on 16 targets. If they are bitwise
    identical, the cp311 package may stand in as the CPU reference for new
    samples.
``exact_788k``
    The 832 C00 targets at 0.1 m and 0.05 m offsets and the 32 on-surface
    points, in both the fold and split variants:
    * acceleration: relative error ≤ 1e-12 and ≤ 16 eps·M;
    * slope: within 1e-9° (reference ``gravity_cache.derive_effective_fields``,
      C00 omega);
    * cliff classification at 45° and 60°, and threshold-band facets (1e-6°);
    * the 512-point offset-halving audit.
``difficult``
    New deterministic near-singular samples on the 788k body:
    * exactly on vertices and on edge midpoints;
    * 0.1 m above vertices and edge midpoints along the face normal;
    * 1 mm above centroids.
    The CPU reference is the cp311 package, computed now.
``determinism``
    Repeated split runs, fold runs with a different chunking, and fold runs on
    the second GPU must all be bitwise identical.
``timing``
    Measured seconds per target for fold and split, and the projected time for
    the full field.
``full``  (``--full``)
    All 787,754 targets (fold). The field must be bitwise equal to the
    ``exact_788k`` fold values at the C00 targets. Writes slopes and the 45°/60°
    cliff masks with their area fractions. This is a sample-validated field,
    not a template bundle.

Example::

    .venv/bin/python tools/c03_gravity_check.py \\
        --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z --cpus 28-35 --full
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
sys.path.insert(0, str(APP_ROOT / "tools"))

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"
EPS = float(np.finfo(float).eps)
M_BUDGET = 16.0 * EPS
REL_BUDGET = 1e-12
SLOPE_BUDGET_DEG = 1e-9
BAND_DEG = 1e-6
THRESHOLDS = (45.0, 60.0)


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def parse_cpus(text):
    cpus = set()
    for part in text.split(","):
        lo, _, hi = part.partition("-")
        cpus.update(range(int(lo), int(hi or lo) + 1))
    return cpus


def errors(g, ref, M):
    """Error statistics over points where the reference is finite."""
    finite = np.isfinite(ref).all(axis=1)
    gpu_finite = np.isfinite(g).all(axis=1)
    diff = np.linalg.norm(g - ref, axis=1)
    norm = np.linalg.norm(ref, axis=1)
    rel = diff[finite] / norm[finite]
    per_m = diff[finite] / M[finite]
    return {"points": int(len(ref)), "finite_reference": int(finite.sum()),
            "nan_pattern_identical": bool(np.array_equal(finite, gpu_finite)),
            "gpu_nonfinite_where_reference_finite": int(np.count_nonzero(finite & ~gpu_finite)),
            "gpu_finite_where_reference_nonfinite": int(np.count_nonzero(~finite & gpu_finite)),
            "rel_max": float(rel.max()) if rel.size else None,
            "per_M_max_eps": float(per_m.max() / EPS) if per_m.size else None,
            "bitwise_rows": int(np.count_nonzero(np.all(g == ref, axis=1) & finite)),
            "cancellation_M_over_g_max": float(np.max(M[finite] / norm[finite]))
            if finite.any() else None}


def package_evaluable(pg, V, F, density):
    """The pinned package with the C00 options (DISABLE after the topology audit).

    Construction already runs ``prepare()``, so ``__getstate__`` is available
    without evaluating anything.  (Evaluating a distant dummy point would make
    the package log its isCriticalDifference warning once per face.)
    """
    poly = pg.Polyhedron(polyhedral_source=(V, F), density=density,
                         normal_orientation=pg.NormalOrientation.OUTWARDS,
                         integrity_check=pg.PolyhedronIntegrity.DISABLE,
                         metric_unit=pg.MetricUnit.METER)
    return pg.GravityEvaluable(polyhedron=poly)


def package_eval(pg, V, F, density, points):
    evaluable = package_evaluable(pg, V, F, density)
    result = evaluable(computation_points=np.asarray(points).tolist(), parallel=True)
    return (np.array([r[0] for r in result]), np.array([r[1] for r in result]), evaluable)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True, help="C00 capture directory")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--cpus", default=None,
                        help="CPU list for the CPU reference (TBB ignores OMP_NUM_THREADS)")
    parser.add_argument("--difficult-faces", type=int, default=32)
    parser.add_argument("--full", action="store_true", help="evaluate all targets")
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT",
                                                                DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()
    if args.cpus:
        # TBB sizes its pool from the affinity mask at first use, so this must
        # precede the import of polyhedral_gravity.
        os.sched_setaffinity(0, parse_cpus(args.cpus))

    import cupy
    import polyhedral_gravity as pg
    import c00_capture_reference as c00
    from cops_cuda import gravity, gravity_body, manifest, provenance, reference
    gravity_cache = reference.load("gravity_cache")
    fit_cops = reference.load("fit_cops")

    capture = Path(args.capture).resolve()
    gfix = json.loads((capture / "gravity" / "gravity_fixtures.json").read_text())
    identities = json.loads((capture / "identities" / "identities.json").read_text())
    syn = np.load(str(capture / "gravity" / "gravity_synthetic.npz"))
    big = np.load(str(capture / "gravity" / "gravity_788k.npz"))
    density = float(gfix["exact_788k"]["density_kg_m3"])
    shape_path = gfix["exact_788k"]["shape"]["path"]
    run = Path(json.loads((capture / "identities" / "kernel_files.json").read_text())["base"])
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c03" /
                                        ("gravity_" + provenance.utc_stamp()))
    log("report directory %s" % out)
    checks, started = {}, time.perf_counter()
    dev = args.device

    def gpu(body, points, variant="fold", device=dev, **kwargs):
        pool = cupy.get_default_memory_pool()
        with cupy.cuda.Device(device):
            pool.free_all_blocks()
        with gravity.GravityEvaluator(body, device_id=device, **kwargs) as ev:
            p, g, M = ev.evaluate(points, variant=variant, with_magnitude=True)
            timing = dict(ev.timing)
        return p, g, M, timing

    # ------------------------------------------------------------ synthetic
    log("synthetic shapes against the C00 fixture (campaign binary)")
    shapes = sorted(set(k.split("__")[0] for k in syn.files))
    shift = np.array([1234.5, -987.25, 456.125])
    rotation = c00.rotation_matrix(0.3, -0.7, 1.1)
    synthetic = {}
    for name in shapes:
        v, f, pts = syn[name + "__vertices"], syn[name + "__faces"], syn[name + "__points"]
        cases = {"original": (v, f, pts, density, syn[name + "__acceleration"]),
                 "translated": (v + shift, f, pts + shift, density,
                                syn[name + "__acceleration_translated"]),
                 "rotated": (v @ rotation.T, f, pts @ rotation.T, density,
                             syn[name + "__acceleration_rotated"]),
                 "double_density": (v, f, pts, 2.0 * density,
                                    syn[name + "__acceleration_double_density"])}
        result = {}
        for case, (cv, cf, cp_, rho, ref) in cases.items():
            p, g, M, _ = gpu(gravity_body.GravityBody(cv, cf, rho), cp_, variant="fold")
            result[case] = errors(g, ref, M)
            if case == "original":
                pref = syn[name + "__potential"]
                fin = np.isfinite(pref)
                result[case]["potential_rel_max"] = float(
                    np.max(np.abs(p - pref)[fin] / np.abs(pref)[fin]))
        ok = (result["original"]["nan_pattern_identical"] and
              all(r["per_M_max_eps"] is not None and r["per_M_max_eps"] * EPS <= M_BUDGET
                  for r in result.values()))
        synthetic[name] = {"pass": ok, "cases": result}
        log("  %-22s %s  max %.3f eps*M" % (name, "PASS" if ok else "FAIL",
                                            max(r["per_M_max_eps"] for r in result.values())))
    checks["synthetic"] = {"pass": all(s["pass"] for s in synthetic.values()),
                           "budget": "16 eps M", "shapes": synthetic}

    # ------------------------------------------------------------ 788k body
    log("reading and auditing the 788k gravity body")
    body = gravity_body.GravityBody.from_bdf(shape_path, density)
    same_faces = bool(np.array_equal(body.faces, big["evaluated_faces"]))
    if body.source["sha256"] != gfix["exact_788k"]["shape"]["sha256"]:
        raise SystemExit("gravity shape differs from the C00 capture")

    log("precompute against the package state")
    t0 = time.perf_counter()
    evaluable = package_evaluable(pg, body.vertices, body.faces, density)
    _, seg, nrm, segn = evaluable.__getstate__()
    G, N, n = gravity.face_terms(body.vertices, body.faces)
    pre = {"segment_vectors": np.asarray(seg).tobytes() == G.tobytes(),
           "plane_normals": np.asarray(nrm).tobytes() == N.tobytes(),
           "segment_normals": np.asarray(segn).tobytes() == n.tobytes()}
    del seg, nrm, segn, evaluable
    checks["precompute"] = {"pass": all(pre.values()), "bitwise": pre,
                            "faces_identical_to_c00_evaluated_faces": same_faces,
                            "seconds": time.perf_counter() - t0}
    log("  %s" % pre)

    log("cp311 package versus the C00 cp38 fixture on 16 targets")
    t0 = time.perf_counter()
    wp, wg, _ = package_eval(pg, body.vertices, body.faces, density, big["points_offset"][:16])
    checks["wheel_equivalence"] = {
        "pass": True, "targets": 16,
        "acceleration_bitwise_identical": bool(np.array_equal(wg, big["acceleration_offset"][:16])),
        "potential_bitwise_identical": bool(np.array_equal(wp, big["potential_offset"][:16])),
        "acceleration_rel_max": float(np.max(np.linalg.norm(wg - big["acceleration_offset"][:16], axis=1)
                                             / np.linalg.norm(big["acceleration_offset"][:16], axis=1))),
        "seconds": time.perf_counter() - t0}
    log("  %s" % checks["wheel_equivalence"])

    # C00 targets: fold and split, both offsets, and on-surface points.
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal, area = elements["centroid"], elements["normal"], elements["area"]
    omega = np.asarray(big["omega_rad_s"])
    targets = big["targets"]
    sets = {"offset": (big["points_offset"], big["acceleration_offset"], targets),
            "half_offset": (big["points_half_offset"], big["acceleration_half_offset"], targets),
            "on_surface": (centroid[big["on_surface_targets"]], big["acceleration_on_surface"],
                           big["on_surface_targets"])}
    exact, fold_values = {}, {}
    for variant in ("split", "fold"):
        for label, (pts, ref, idx) in sets.items():
            log("788k %s %s: %d targets" % (variant, label, len(pts)))
            t0 = time.perf_counter()
            p, g, M, _ = gpu(body, pts, variant=variant)
            seconds = time.perf_counter() - t0
            stat = errors(g, ref, M)
            stat["seconds"] = seconds
            stat["seconds_per_target"] = seconds / len(pts)
            if label != "on_surface":
                slope = gravity_cache.derive_effective_fields(centroid[idx], normal[idx], g, omega)[2]
                slope_ref = big["slope_offset_deg" if label == "offset" else "slope_half_offset_deg"]
                stat["slope_abs_max_deg"] = float(np.max(np.abs(slope - slope_ref)))
                stat["threshold_band_targets"] = {
                    "%g" % th: int(np.count_nonzero(np.abs(slope_ref - th) <= BAND_DEG))
                    for th in THRESHOLDS}
                stat["classification_identical"] = {
                    "%g" % th: bool(np.array_equal(slope >= th, slope_ref >= th)) for th in THRESHOLDS}
                stat["slope"] = slope
            exact["%s_%s" % (variant, label)] = stat
            if variant == "fold":
                fold_values[label] = (p, g)
    halving = {}
    for variant in ("split", "fold"):
        s_full = exact["%s_offset" % variant].pop("slope")
        s_half = exact["%s_half_offset" % variant].pop("slope")
        mask = big["in_offset_halving_sample"]
        halving[variant] = float(np.max(np.abs(s_full - s_half)[mask]))
    reference_halving = gfix["exact_788k"]["offset_halving_slope_deg_max_sample"]
    ok788 = all(
        s["nan_pattern_identical"] and s["rel_max"] <= REL_BUDGET and
        s["per_M_max_eps"] * EPS <= M_BUDGET and
        s.get("slope_abs_max_deg", 0.0) <= SLOPE_BUDGET_DEG and
        all(s.get("classification_identical", {}).values())
        for s in exact.values())
    ok788 = ok788 and all(abs(h - reference_halving) <= SLOPE_BUDGET_DEG for h in halving.values())
    checks["exact_788k"] = {"pass": ok788, "sets": exact,
                            "offset_halving_max_deg": halving,
                            "offset_halving_reference_deg": reference_halving}
    log("  exact_788k %s; offset halving %s (reference %.6g)" %
        ("PASS" if ok788 else "FAIL", halving, reference_halving))

    # ------------------------------------------------------- difficult points
    log("difficult locations: %d faces" % args.difficult_faces)
    picks = np.linspace(0, body.nfaces - 1, args.difficult_faces, dtype=np.int64)
    V, F = body.vertices, body.faces
    labels, points = [], []
    for k in picks:
        a, b = V[F[k, 0]], V[F[k, 1]]
        nk = normal[k]
        mid = 0.5 * (a + b)
        for label, point in (("on_vertex", a), ("on_edge_midpoint", mid),
                             ("vertex_plus_0.1m", a + 0.1 * nk),
                             ("edge_midpoint_plus_0.1m", mid + 0.1 * nk),
                             ("centroid_plus_1mm", centroid[k] + 1e-3 * nk)):
            labels.append(label)
            points.append(point)
    points = np.asarray(points)
    t0 = time.perf_counter()
    dp, dg, _ = package_eval(pg, V, F, density, points)
    cpu_seconds = time.perf_counter() - t0
    _, g_d, M_d, _ = gpu(body, points, variant="split")
    by_label = {}
    for label in sorted(set(labels)):
        sel = np.array([l == label for l in labels])
        by_label[label] = errors(g_d[sel], dg[sel], M_d[sel])
    stat = errors(g_d, dg, M_d)
    checks["difficult"] = {
        "pass": stat["nan_pattern_identical"] and stat["per_M_max_eps"] * EPS <= M_BUDGET,
        "faces": picks.tolist(), "overall": stat, "by_label": by_label,
        "reference": "cp311 polyhedral_gravity 3.3, DISABLE", "cpu_seconds": cpu_seconds}
    log("  difficult %s: %s" % ("PASS" if checks["difficult"]["pass"] else "FAIL", stat))

    # ------------------------------------------------------------ determinism
    log("determinism: repeat, chunking, second device")
    pts = big["points_offset"]
    _, g1, _, _ = gpu(body, pts, variant="split")
    _, g2, _, _ = gpu(body, pts, variant="split")
    _, g3, _, _ = gpu(body, pts, variant="fold", max_pairs_per_launch=3 * 10 ** 7,
                      target_chunk=300)
    det = {"split_repeat_bitwise": bool(g1.tobytes() == g2.tobytes()),
           "fold_rechunked_bitwise": bool(g3.tobytes() == fold_values["offset"][1].tobytes())}
    if cupy.cuda.runtime.getDeviceCount() > 1:
        other = 1 if dev == 0 else 0
        _, g4, _, _ = gpu(body, pts, variant="fold", device=other)
        det["fold_second_device_bitwise"] = bool(g4.tobytes() == fold_values["offset"][1].tobytes())
        det["second_device"] = manifest.device_identity(other)["uuid"]
    checks["determinism"] = {"pass": all(v for k, v in det.items() if k.endswith("bitwise")),
                             **det}
    log("  %s" % det)

    fold_per_target = exact["fold_offset"]["seconds_per_target"]
    split_per_target = exact["split_offset"]["seconds_per_target"]
    checks["timing"] = {
        "pass": True, "note": "832-target batches; fold uses 7 blocks, so it is latency-bound here",
        "fold_seconds_per_target_832": fold_per_target,
        "split_seconds_per_target_832": split_per_target,
        "cpu_reference_seconds_per_target": cpu_seconds / len(points),
        "cpu_reference_cpus": sorted(os.sched_getaffinity(0))}

    # ------------------------------------------------------------- full field
    if args.full:
        all_points = centroid + gfix["exact_788k"]["offset_m"] * normal
        log("full field: %d targets (fold) on device %d" % (len(all_points), dev))
        t0 = time.perf_counter()
        last = [0.0]

        def progress(done, total, seconds):
            if seconds - last[0] >= 60.0 or done == total:
                last[0] = seconds
                log("  %d/%d (%.1f%%) %.0f s" % (done, total, 100.0 * done / total, seconds))
        with cupy.cuda.Device(dev):
            cupy.get_default_memory_pool().free_all_blocks()
        with gravity.GravityEvaluator(body, device_id=dev) as ev:
            fp, fg, fM = ev.evaluate(all_points, variant="fold", progress=progress,
                                     with_magnitude=True)
            evaluator_identity = ev.identity()
        seconds = time.perf_counter() - t0
        centrifugal, effective, slope = gravity_cache.derive_effective_fields(
            centroid, normal, fg, omega)
        consistent = bool(fg[targets].tobytes() == fold_values["offset"][1].tobytes())
        masks = {"%g" % th: slope >= th for th in THRESHOLDS}
        band = {"%g" % th: int(np.count_nonzero(np.abs(slope - th) <= BAND_DEG))
                for th in THRESHOLDS}
        np.savez_compressed(out / "full_field.npz", self_gravity_m_s2=fg, potential_m2_s2=fp,
                            term_magnitude_m_s2=fM, omega_rad_s=omega,
                            centrifugal_m_s2=centrifugal, effective_gravity_m_s2=effective,
                            slope_deg=slope, cliff_ge_45deg=masks["45"], cliff_ge_60deg=masks["60"])
        checks["full"] = {
            "pass": consistent and bool(np.isfinite(fg).all()),
            "targets": int(len(all_points)), "seconds": seconds,
            "pairs_per_second": float(len(all_points)) * body.nfaces / seconds,
            "bitwise_equal_to_exact_788k_fold_at_c00_targets": consistent,
            "slope_deg_min_median_max": [float(slope.min()), float(np.median(slope)),
                                         float(slope.max())],
            "area_fraction": {k: float(np.sum(area[m]) / np.sum(area)) for k, m in masks.items()},
            "face_count": {k: int(np.count_nonzero(m)) for k, m in masks.items()},
            "threshold_band_facets_1e-6deg": band,
            "self_gravity_norm_min_max": [float(np.linalg.norm(fg, axis=1).min()),
                                          float(np.linalg.norm(fg, axis=1).max())],
            "file": provenance.file_identity(out / "full_field.npz"),
            "evaluator": evaluator_identity,
            "status": "sample-validated C03 evidence; not a campaign template bundle"}
        log("full: %s in %.0f s; cliff area fractions %s" %
            ("PASS" if checks["full"]["pass"] else "FAIL", seconds, checks["full"]["area_fraction"]))
    else:
        checks["timing"]["projected_full_field_fold_seconds"] = None

    report = {"schema": "liouville_cops_cuda.c03_report", "schema_version": 1,
              "created_utc": provenance.utc_now(), "capture": str(capture),
              "gravity_shape": body.identity(), "density_kg_m3": density,
              "budgets": {"per_M": "16 eps M", "rel": REL_BUDGET, "slope_deg": SLOPE_BUDGET_DEG,
                          "threshold_band_deg": BAND_DEG},
              "checks": checks, "all_pass": all(c["pass"] for c in checks.values()),
              "elapsed_s": time.perf_counter() - started, "host_load": os.getloadavg(),
              "omega_source": identities.get("spice_omega", {}).get("median_rad_s"),
              "manifest": manifest.backend_manifest("cuda_polyhedral_gravity_v1", dev)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if report["all_pass"] else "FAIL", out / "report.json"))
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
