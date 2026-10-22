#!/usr/bin/env python3
"""Build an FP32 record store of a fitting base (SCALING_PLAN.md 3.1, P1).

    .venv/bin/python tools/c15_record_store.py --run BASE --out STORE [--threads 8]

The store holds the base's kernel records in the FP32 device layout of the mixed-precision
experiment, one file per point, with a manifest that ties every file to its kernel file (size,
mtime, SHA-256) and the base's run_manifest.json (``cops_cuda/record_store.py``). Use it with
``gpu_fit.py --precision mixed --record-store STORE``. The CPU-format kernel files are not
changed and remain the reference; the store can be deleted and rebuilt.
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True, help="fitting base directory")
    ap.add_argument("--out", required=True, help="new store directory (must not exist)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    from cops_cuda import record_store, reference
    reference.bind()
    t0 = time.perf_counter()
    last = [0.0]

    def progress(done, total):
        now = time.perf_counter()
        if now - last[0] > 30 or done == total:
            last[0] = now
            print("[%6.0f s] %d/%d points" % (now - t0, done, total), flush=True)

    manifest = record_store.build(args.run, args.out, threads=args.threads, progress=progress)
    size = manifest["records"] * record_store.BYTES_PER_RECORD
    print("wrote %s: %d points, %d records, %.1f GB, kernel set %s (run_manifest verified: %s), "
          "max relative rounding %.3g, lit adjustments %d, %.0f s" %
          (args.out, len(manifest["points"]), manifest["records"], size / 1e9,
           manifest["kernel_set_sha256"][:16], manifest["kernel_set_verified_against_run_manifest"],
           manifest["conversion"]["max_relative_rounding"],
           manifest["conversion"]["lit_threshold_adjusted"], time.perf_counter() - t0), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
