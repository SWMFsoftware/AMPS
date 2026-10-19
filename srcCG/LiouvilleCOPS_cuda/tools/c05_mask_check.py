#!/usr/bin/env python3
"""C05 evidence: GPU visibility masks against the CPU predicates (CUDA env).

For the campaign's late interval (default 2016-09-30T10:33:10–10:39:26) it
computes the four GPU masks (``cops_cuda.masks``) and compares them bit for bit
with:

* ``late_ng_not_rg`` from the reference builder's own function
  ``late_gauge_template``, called unchanged. Its counts of NG-seen and RG-seen
  facets are compared as well;
* ``late_rg_not_ng`` from the independent C00b CPU replica
  (``joint/templates.npz``) and from a NumPy replica computed here over the
  same points;
* the per-point predicates' union, recomputed in NumPy.

It also runs a control interval with both gauges active (the C00 negative
case). Writes ``$COPS_OUTPUT_ROOT/cuda_c05/masks_<UTC>/report.json`` and
``masks.npz``.
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
DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"
LATE = ("2016-09-30T10:33:10", "2016-09-30T10:39:26")
CONTROL = ("2016-09-30T09:00:00", "2016-09-30T09:10:00")


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def numpy_replica(fit_cops, run, points, nfaces):
    out = {k: np.zeros(nfaces, dtype=bool) for k in ("ng_seen", "rg_seen", "late_ng_not_rg", "late_rg_not_ng")}
    for p in points:
        _, _, rec = fit_cops.read_kernel(str(Path(run) / "kernels" / ("point=%04d.bin" % p)))
        ng = rec["Kn"] > 0.0
        rg = (rec["KPhi"] != 0.0) | (rec["KV"] != 0.0)
        out["ng_seen"][rec["face"][ng]] = True
        out["rg_seen"][rec["face"][rg]] = True
        out["late_ng_not_rg"][rec["face"][ng & ~rg]] = True
        out["late_rg_not_ng"][rec["face"][rg & ~ng]] = True
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run", default=os.path.join(os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT),
                                                      "output_fit_shap5_788k_base"))
    parser.add_argument("--capture", help="C00 capture with the joint stage (RG-only replica)")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    from cops_cuda import manifest, masks, provenance, reference
    fit_cops = reference.load("fit_cops")
    cops_io = reference.load("cops_io")
    builder = reference.load("build_physical_activity_templates")
    run = Path(args.run).resolve()
    audit_header, names, rows, _ = cops_io.read_model_table(run / "COPS_model.dat")
    nfaces = len(fit_cops.read_elements(str(run / "elements.bin"))["area"])
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c05" /
                                        ("masks_" + provenance.utc_stamp()))
    checks, saved = {}, {}
    for label, (start, stop) in (("late", LATE), ("control", CONTROL)):
        t0 = time.perf_counter()
        gpu, meta = masks.late_masks(run, start, stop, device_id=args.device)
        gpu_seconds = time.perf_counter() - t0
        replica = numpy_replica(fit_cops, run, meta["selected_point_ids"], nfaces)
        # The reference builder's own function, unchanged.
        t0 = time.perf_counter()
        ref_mask, _, _, ref_meta = builder.late_gauge_template(
            str(run), audit_header, names, rows, nfaces, start, stop)
        cpu_seconds = time.perf_counter() - t0
        check = {
            "points": meta["selected_point_count"],
            "points_equal_reference": ref_meta["selected_point_ids"] == meta["selected_point_ids"],
            "late_ng_not_rg_equals_reference_builder": bool(np.array_equal(gpu["late_ng_not_rg"], ref_mask > 0.0)),
            "ng_seen_count_equals_reference": meta["face_counts"]["ng_seen"] == ref_meta["ng_seen_face_count"],
            "rg_seen_count_equals_reference": meta["face_counts"]["rg_seen"] == ref_meta["rg_seen_face_count"],
            "all_masks_equal_numpy_replica": all(np.array_equal(gpu[k], replica[k]) for k in gpu),
            "face_counts": meta["face_counts"], "gpu_seconds": gpu_seconds,
            "cpu_reference_seconds": cpu_seconds}
        if label == "late" and args.capture:
            c00b = np.load(str(Path(args.capture) / "joint" / "templates.npz"))
            check["late_rg_not_ng_equals_c00b_replica"] = bool(
                np.array_equal(gpu["late_rg_not_ng"], c00b["late_rg_not_ng"]))
            check["late_ng_not_rg_equals_c00b"] = bool(np.array_equal(gpu["late_ng_not_rg"], c00b["late_ng_not_rg"]))
        check["pass"] = all(v for k, v in check.items() if isinstance(v, bool))
        checks[label] = check
        for k, m in gpu.items():
            saved["%s__%s" % (label, k)] = m
        saved["%s__point_ids" % label] = np.asarray(meta["selected_point_ids"], dtype=np.int32)
        log("%s: %s" % (label, check))
    with (out / "masks.npz").open("xb") as stream:
        np.savez_compressed(stream, **saved)
    report = {"schema": "liouville_cops_cuda.c05_masks", "schema_version": 1,
              "created_utc": provenance.utc_now(), "run": str(run), "checks": checks,
              "all_pass": all(c["pass"] for c in checks.values()),
              "masks_file": provenance.file_identity(out / "masks.npz"),
              "manifest": manifest.backend_manifest("cuda_masks_v1", args.device)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; %s" % ("PASS" if report["all_pass"] else "FAIL", out))
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
