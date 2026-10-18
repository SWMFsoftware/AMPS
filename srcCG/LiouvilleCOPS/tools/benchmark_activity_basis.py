#!/usr/bin/env python3
"""Benchmark the Stage-A2 Python harmonic table on a real elements.bin.

The tool reads only the version-1 element header and centroid columns.  It does
not modify the cache and does not write a high-order basis into elements.bin.
The exact dense-table allocation is reported separately from process maximum
RSS, whose delta is allocator/platform dependent.
"""

import argparse
import json
from pathlib import Path
import resource
import struct
import time

import numpy as np

from activity_basis import BASIS_ID, coefficient_count, evaluate_cartesian
from cops_paths import prepare_output, product_input


def read_centroids(path):
    """Memory-map centroid columns after strict version/size validation."""
    path = Path(path)
    with path.open("rb") as stream:
        header = stream.read(12)
    if len(header) != 12:
        raise ValueError("truncated element header: %s" % path)
    version, faces, legacy_basis = struct.unpack("<iii", header)
    if version != 1 or faces <= 0 or legacy_basis != 25:
        raise ValueError("unsupported element header in %s: version=%d faces=%d basis=%d" %
                         (path, version, faces, legacy_basis))
    expected = 12 + faces * 34 * 8
    if path.stat().st_size != expected:
        raise ValueError("element size mismatch in %s: %d; expected %d" %
                         (path, path.stat().st_size, expected))
    records = np.memmap(path, mode="r", dtype="<f8", offset=12,
                        shape=(faces, 34))
    return records[:, :3], {"version": version, "faces": faces,
                            "legacy_basis": legacy_basis,
                            "size_bytes": expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("elements", nargs="?",
                        default="output_fit_full_default_base/elements.bin",
                        help="version-1 elements.bin (default: %(default)s)")
    parser.add_argument("--min-degree", type=int, default=5)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--repeat", type=int, default=3,
                        help="number of timed constructions (default: %(default)s)")
    parser.add_argument("--json", help="optional JSON output under COPS_OUTPUT_ROOT")
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")

    path = product_input(args.elements)
    centroids, header = read_centroids(path)
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    elapsed_samples = []
    checksum = maximum = table_bytes = None
    for _ in range(args.repeat):
        started = time.perf_counter()
        table = evaluate_cartesian(centroids, args.min_degree, args.max_degree)
        # The reductions force every output page to be materialized before the
        # sample ends and provide repeatability diagnostics without saving it.
        checksum = float(np.sum(table, dtype=np.float64))
        maximum = float(np.max(np.abs(table)))
        elapsed_samples.append(time.perf_counter() - started)
        table_bytes = int(table.nbytes)
        del table
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {
        "basis_id": BASIS_ID,
        "elements": str(Path(path).resolve()),
        "element_header": header,
        "min_degree": args.min_degree,
        "max_degree": args.max_degree,
        "modes": coefficient_count(args.min_degree, args.max_degree),
        "table_shape": [header["faces"], coefficient_count(args.min_degree,
                                                             args.max_degree)],
        "table_bytes": table_bytes,
        "table_mib": table_bytes / (1024.0 ** 2),
        "repeat": args.repeat,
        "elapsed_samples_seconds": elapsed_samples,
        "elapsed_median_seconds": float(np.median(elapsed_samples)),
        "elapsed_minimum_seconds": min(elapsed_samples),
        "max_rss_delta_kib": max(0, int(rss_after - rss_before)),
        "checksum": checksum,
        "maximum_absolute_value": maximum,
        "cache_modified": False,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.json:
        destination = Path(prepare_output(args.json))
        destination.write_text(text + "\n", encoding="utf-8")
        print("result: %s" % destination)


if __name__ == "__main__":
    main()
