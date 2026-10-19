#!/usr/bin/env python3
"""C02 evidence on the campaign data: load, verify, cycle, refuse (CUDA env).

Writes one exclusive report directory
``$COPS_OUTPUT_ROOT/cuda_c02/load_<UTC>/report.json``. Each check produces a
pass/fail entry with its measurements:

``subset``
    The 47 C00 fit points. Each file is compared with the C00 per-file
    digests, the device arrays are compared bit for bit with the CPU
    reference reader ``fit_cops.read_kernel`` (an independent code path),
    and a device round trip is run.
``cycles``
    Repeated load/release of the subset. Device digests must be identical
    in every cycle, the CuPy pool must be empty after each release, and
    device free memory must return to its starting level.
``full``  (``--full``)
    All 1150 points in table order. The loaded set digest must equal the
    ``set_sha256`` in ``run_manifest.json``. Then the device audit, a round
    trip, and measured versus planned memory, including the per-process
    usage reported by nvidia-smi (CUDA context plus pool).
``capacity``
    The full plan against a capacity just below it. ``CapacityError`` must
    be raised before any allocation.
``second_device``
    The subset on the other GPU, so that the device identity is recorded per
    device.

Nothing here modifies the run directory. Reads are sequential and
single-threaded, so the tool can run beside the CPU campaign; contention is
recorded in the report.

Example::

    .venv/bin/python tools/c02_device_load.py \\
        --capture $COPS_OUTPUT_ROOT/cuda_reference/c00_20261019T064536Z --full
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import device_store, manifest, memory, provenance, reference  # noqa: E402
from cops_cuda.errors import CapacityError  # noqa: E402

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"


def log(message):
    print("[%s] %s" % (provenance.utc_now()[11:19], message), flush=True)


def process_device_memory():
    """Device memory nvidia-smi attributes to this process, in bytes, per GPU UUID."""
    try:
        out = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,gpu_uuid,used_memory",
                              "--format=csv,noheader,nounits"], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    result = {}
    for line in out.splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) == 3 and parts[0] == str(os.getpid()):
            result[parts[1]] = int(parts[2]) << 20
    return result


def cpu_reference_compare(ds, run, point_ids):
    """Compare the device arrays with fit_cops.read_kernel output, bitwise."""
    fit_cops = reference.load("fit_cops")
    a, start, mismatches = ds.arrays, 0, []
    host = {name: a["kernel." + name].get() for name, _ in device_store.RECORD_FIELDS}
    offsets = a["kernel.offsets"].get()
    for index, point in enumerate(point_ids):
        ip, header, rec = fit_cops.read_kernel(str(Path(run) / "kernels" /
                                                   ("point=%04d.bin" % point)))
        end = start + len(rec)
        if ip != point or offsets[index] != start or offsets[index + 1] != end:
            mismatches.append((point, "offsets"))
        for name, field in device_store.RECORD_FIELDS:
            expected = rec[field].astype(np.uint8) if name == "shadow" else rec[field]
            if host[name][start:end].tobytes() != np.ascontiguousarray(expected).tobytes():
                mismatches.append((point, name))
        if a["kernel.geometry"][1, index].item() != header["VdotB"]:
            mismatches.append((point, "VdotB"))
        start = end
    return {"points": len(point_ids), "records": int(start), "mismatches": mismatches,
            "pass": not mismatches and start == ds.nrecords}


def free_memory(device_id):
    import cupy
    with cupy.cuda.Device(device_id):
        return int(cupy.cuda.runtime.memGetInfo()[0])


def pool_state(device_id):
    import cupy
    with cupy.cuda.Device(device_id):
        pool = cupy.get_default_memory_pool()
        return {"used": int(pool.used_bytes()), "total": int(pool.total_bytes())}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", default=None,
                        help="base run (default: the run recorded in the C00 capture)")
    parser.add_argument("--capture", required=True, help="C00 capture directory")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--cycles", type=int, default=3)
    parser.add_argument("--full", action="store_true", help="also load all points")
    parser.add_argument("--chunk-records", type=int, default=device_store.DEFAULT_CHUNK_RECORDS)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT",
                                                                DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    import cupy
    capture = Path(args.capture).resolve()
    kernel_files = json.loads((capture / "identities" / "kernel_files.json").read_text())
    fixtures = json.loads((capture / "fit" / "fit_fixtures.json").read_text())
    run = Path(args.base or kernel_files["base"]).resolve()
    expected = {entry["file"]: entry["sha256"] for entry in kernel_files["kernels"]}
    subset = [int(p) for p in fixtures["subset_point_ids"]]
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c02" /
                                        ("load_" + provenance.utc_stamp()))
    log("report directory %s" % out)
    checks, started = {}, time.perf_counter()
    device = args.device

    # ------------------------------------------------------------- subset
    log("subset: %d C00 points on device %d" % (len(subset), device))
    with device_store.DeviceDataset.load(run, point_ids=subset, device_id=device,
                                         expected_kernel_digests=expected,
                                         chunk_records=args.chunk_records) as ds:
        compare = cpu_reference_compare(ds, run, subset)
        ds.verify_roundtrip()
        identity = ds.identity()
        pooled = pool_state(device)
    checks["subset"] = {"pass": compare["pass"] and identity["kernels"][
        "files_verified_against_expected"] == len(subset),
        "cpu_reference_compare": compare, "identity": identity,
        "pool_while_loaded": pooled, "pool_after_release": pool_state(device)}
    log("subset: %s, %d records, cpu-reference mismatches %d" %
        ("PASS" if checks["subset"]["pass"] else "FAIL", compare["records"],
         len(compare["mismatches"])))

    # ------------------------------------------------------------- cycles
    baseline = free_memory(device)
    digests, frees, pools = [], [], []
    for cycle in range(args.cycles):
        with device_store.DeviceDataset.load(run, point_ids=subset, device_id=device,
                                             chunk_records=args.chunk_records) as ds:
            digests.append(ds.verify_roundtrip())
        frees.append(free_memory(device))
        pools.append(pool_state(device))
    drift = max(abs(f - baseline) for f in frees) if frees else 0
    checks["cycles"] = {
        "pass": all(d == digests[0] for d in digests) and
        all(p == {"used": 0, "total": 0} for p in pools) and drift <= (2 << 20),
        "cycles": args.cycles, "identical_digests": all(d == digests[0] for d in digests),
        "free_bytes_baseline": baseline, "free_bytes_after_each": frees,
        "max_free_drift_bytes": drift, "pool_after_each": pools}
    log("cycles: %s (%d cycles, free-memory drift %d bytes)" %
        ("PASS" if checks["cycles"]["pass"] else "FAIL", args.cycles, drift))

    # --------------------------------------------------------------- full
    full_plan = None
    if args.full:
        log("full: all points in table order on device %d" % device)
        before = free_memory(device)
        with device_store.DeviceDataset.load(run, device_id=device,
                                             chunk_records=args.chunk_records) as ds:
            full_plan = ds.plan
            loaded_free = free_memory(device)
            nvsmi = process_device_memory()
            log("full: loaded %s in %.1f s; verifying round trip" %
                (ds.summary(), ds.timing["load_total_s"]))
            ds.verify_roundtrip()
            ident = ds.identity()
            backend = manifest.backend_manifest("c02_device_dataset", device, ds)
            pooled = pool_state(device)
        after = free_memory(device)
        records = ident["kernels"]["records"]
        checks["full"] = {
            "pass": ident["kernels"]["verified_against_run_manifest"] and
            ident["device_audit"] == dict.fromkeys(device_store.AUDIT_COUNTERS, 0) and
            pool_state(device) == {"used": 0, "total": 0},
            "identity": ident,
            "memory": {"planned_total_bytes": full_plan.total_bytes,
                       "pool_bytes_while_loaded": pooled,
                       "device_free_before": before, "device_free_loaded": loaded_free,
                       "device_free_after_release": after,
                       "measured_device_use_bytes": before - loaded_free,
                       "nvidia_smi_process_bytes": nvsmi},
            "throughput": {"records_per_s": records / ident["timing_s"]["kernel_stream_s"],
                           "disk_bytes_per_s": records * 40 /
                           ident["timing_s"]["kernel_stream_s"]}}
        # The manifest of the loaded full dataset, as a C06 product would carry it.
        provenance.write_json_exclusive(out / "backend_manifest.json", backend)
        log("full: %s; planned %s, measured %s, stream %.1f s, round trip %.1f s" %
            ("PASS" if checks["full"]["pass"] else "FAIL", memory.human(full_plan.total_bytes),
             memory.human(before - loaded_free), ident["timing_s"]["kernel_stream_s"],
             ident["timing_s"].get("roundtrip_verify_s", float("nan"))))

    # ----------------------------------------------------------- capacity
    inputs = device_store.RunInputs(run)
    scan = device_store.KernelScan(inputs, inputs.resolve_rows())
    plan = device_store.dataset_plan(inputs.nfaces, len(scan.rows), scan.nrecords,
                                     device_store.MESH_FIELD_SETS["evaluator"])
    total = cupy.cuda.Device(device).mem_info[1]
    # Free memory such that free - 20% headroom - slack is one byte short.
    fake_free = plan.total_bytes + int(np.ceil(0.2 * total)) + memory.DEFAULT_SLACK_BYTES - 1
    budget = memory.DeviceBudget(device, capacity=(fake_free, total))
    t0 = time.perf_counter()
    try:
        device_store.DeviceDataset.load(run, device_id=device, budget=budget)
        refused, message = False, None
    except CapacityError as exc:
        refused, message = True, str(exc)
    refusal_s = time.perf_counter() - t0
    pool_after = pool_state(device)
    budget.close()
    checks["capacity"] = {"pass": refused and pool_after == {"used": 0, "total": 0},
                          "plan_bytes": plan.total_bytes, "injected_free_bytes": fake_free,
                          "refusal_seconds": refusal_s, "message": message,
                          "pool_after": pool_after}
    log("capacity: %s (refused in %.1f s, nothing allocated)" %
        ("PASS" if checks["capacity"]["pass"] else "FAIL", refusal_s))

    # ------------------------------------------------------ second device
    if cupy.cuda.runtime.getDeviceCount() > 1:
        other = 1 - device if device in (0, 1) else 0
        with device_store.DeviceDataset.load(run, point_ids=subset, device_id=other,
                                             expected_kernel_digests=expected) as ds:
            digests_other = ds.verify_roundtrip()
            device_info = manifest.device_identity(other)
        checks["second_device"] = {
            "pass": digests_other == digests[0] if digests else True,
            "device": device_info, "identical_to_first_device": digests_other == digests[0]
            if digests else None}
        log("second device %s: %s" % (device_info["uuid"],
                                      "PASS" if checks["second_device"]["pass"] else "FAIL"))

    report = {"schema": "liouville_cops_cuda.c02_report", "schema_version": 1,
              "created_utc": provenance.utc_now(), "run": str(run),
              "capture": str(capture), "checks": checks,
              "all_pass": all(c["pass"] for c in checks.values()),
              "elapsed_s": time.perf_counter() - started,
              "host_load": os.getloadavg(),
              "manifest": manifest.backend_manifest("c02_device_dataset", device)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if report["all_pass"] else "FAIL", out / "report.json"))
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
