#!/usr/bin/env python3
"""Write a GPU self-gravity field as a CPU-format gravity cache (C04; CUDA env).

Input is a GPU field directory, either:
* ``cuda_c04/equivalence_<UTC>`` (``gpu_field.npz``: 788k targets, any body, plus
  the half-offset sample); or
* a C03 report directory (``full_field.npz`` for the 788k body). The half-offset
  values of the builder's sample are then re-evaluated on the GPU.

ω is obtained exactly as the builder obtains it: the reference
``spice_omega`` with the builder's three epochs of the declared late interval
and the run's ``input.snapshot.in``. The cache is written with
``cops_cuda.gravity_cache_io.write_gpu_cache`` (CPU writer, truthful GPU
metadata). Optional ``--evidence`` links the equivalence report that admitted
the field.

The unchanged CPU builder then reuses the cache with
``--gravity-cache <path>``; see RUN.md.
"""

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--field-dir", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--gravity-shape", required=True)
    parser.add_argument("--out", required=True, help="new cache .npz (must not exist)")
    parser.add_argument("--density", type=float, default=533.0)
    parser.add_argument("--offset", type=float, default=0.1)
    parser.add_argument("--offset-tolerance", type=float, default=0.25)
    parser.add_argument("--late-start", default="2016-09-30T10:33:10")
    parser.add_argument("--late-stop", default="2016-09-30T10:39:26")
    parser.add_argument("--evidence", action="append", default=[],
                        help="report(s) that admitted this field (identified by content)")
    parser.add_argument("--device", type=int, default=0)
    args = parser.parse_args()

    from cops_cuda import gravity, gravity_body, gravity_cache_io, provenance, reference
    builder = reference.load("build_physical_activity_templates")
    gravity_cache = reference.load("gravity_cache")
    fit_cops = reference.load("fit_cops")
    run = Path(args.run).resolve()
    field_dir = Path(args.field_dir).resolve()
    body = gravity_body.GravityBody.from_bdf(args.gravity_shape, args.density)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal = elements["centroid"], elements["normal"]
    sample = gravity_cache_io.offset_sample(len(centroid))

    if (field_dir / "gpu_field.npz").is_file():
        record = json.loads((field_dir / "gpu_manifest.json").read_text())
        if record["gravity_shape"]["sha256"] != body.source["sha256"]:
            raise SystemExit("GPU field was computed for another gravity shape")
        field = np.load(str(field_dir / "gpu_field.npz"))
        g = field["self_gravity_m_s2"]
        if not np.array_equal(field["offset_sample"], sample):
            raise SystemExit("GPU field's offset sample differs from the builder's")
        g_half = field["half_offset_self_gravity_m_s2"]
        evaluator_identity, seconds = record["evaluator"], record["seconds"]
    else:
        report = json.loads((field_dir / "report.json").read_text())
        if report["gravity_shape"]["source"]["sha256"] != body.source["sha256"]:
            raise SystemExit("C03 field was computed for another gravity shape")
        with np.load(str(field_dir / "full_field.npz"), allow_pickle=False) as full:
            g = np.asarray(full["self_gravity_m_s2"])
        half = centroid[sample] + 0.5 * args.offset * normal[sample]
        with gravity.GravityEvaluator(body, device_id=args.device) as ev:
            _, g_half = ev.evaluate(half, variant="fold")
        evaluator_identity = report["checks"]["full"]["evaluator"]
        seconds = report["checks"]["full"]["seconds"]

    # The builder's omega: reference spice_omega at its three epochs.
    epochs = [args.late_start, builder.utc_midpoint(args.late_start, args.late_stop), args.late_stop]
    omega, variation, command = builder.spice_omega(
        (reference.bind() / "spice_angular_velocity").resolve(), run / "input.snapshot.in", epochs)
    slope = gravity_cache.derive_effective_fields(centroid[sample], normal[sample], g[sample], omega)[2]
    slope_half = gravity_cache.derive_effective_fields(centroid[sample], normal[sample], g_half, omega)[2]
    diff = np.abs(slope - slope_half)
    evidence = [provenance.file_identity(path, "admission_evidence") for path in args.evidence]
    t0 = time.perf_counter()
    written = gravity_cache_io.write_gpu_cache(
        Path(args.out), run, body, g, omega, args.offset,
        (float(diff.max()), float(np.percentile(diff, 95.0)), len(sample)),
        args.offset_tolerance, evaluator_identity, evidence, seconds)
    sidecar = Path(args.out).with_suffix(".provenance.json")
    provenance.write_json_exclusive(sidecar, {
        "schema": "liouville_cops_cuda.c04_gpu_cache", "schema_version": 1,
        "created_utc": provenance.utc_now(), "cache": provenance.file_identity(args.out),
        "field_dir": str(field_dir), "gravity_shape": body.source, "run": str(run),
        "omega_rad_s": [float(v) for v in omega], "omega_sample_deviation_rad_s": float(variation),
        "omega_command": [str(c) for c in command],
        "offset_halving_deg_max": float(diff.max()), "evidence": evidence,
        "write_seconds": time.perf_counter() - t0,
        "array_sha256": written["metadata"]["array_sha256"]})
    print("wrote %s (offset halving max %.6g deg); provenance %s" % (args.out, diff.max(), sidecar))
    return 0


if __name__ == "__main__":
    sys.exit(main())
