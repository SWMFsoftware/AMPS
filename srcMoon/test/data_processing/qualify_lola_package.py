#!/usr/bin/env python3
"""Reproducibly qualify the immutable D01/D02 LOLA raw files.

This script never edits ``raw/`` and never derives model values.  It verifies
the exact products, detached-label byte contracts, and SHA-256 digests, then
writes package-level provenance/QA metadata beside the immutable directory.
D01 can be READY independently; D02 remains PARTIAL until a separately
reviewed script creates and verifies the independent comparison table required
by I15/I20.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re


CANONICAL_ROOT = Path("/data/vtenishe/moon_validation_data")

# Frozen after downloading from the exact PDS/PGDA paths in the roadmap.  A
# mismatch is data corruption or a product-version change and must be reviewed;
# this script never updates expected hashes from the observed files.
EXPECTED = {
    "LDEM_4.IMG": "972d66d7f330abca62abc21aae270661516c00dd1a881ec102b94df9563656b5",
    "LDEM_4.LBL": "b5a197b4d417e3b5e0e84f8fe15f9a1d2203b8eb5d7213236c0ffbc86e940cc6",
    "illumination/AVGVISIB_65N_240M_201608.IMG": "4d96400a8609cfb87004010cc3affbefcb2064c4434750945a93bddc6c64254c",
    "illumination/AVGVISIB_65N_240M_201608.LBL": "ac839d42ea50ec34314d0c68bc3553d74a2719c88acca1ac7cf3b572a1742d71",
    "illumination/AVGVISIB_65S_240M_201608.IMG": "4b983fb1340c7b625a36d253e881874271c330083d2fdcc7277dc2f8ac24a438",
    "illumination/AVGVISIB_65S_240M_201608.LBL": "10b0d6d881713d6bd4d61a30cc64da6d53934e2c9a8cc702b36a5cd758751b03",
    "illumination/LPSR_65N_240M_201608.IMG": "b693690d9e7d1e31ef22cc90d1fedc06b78c7bf27d0f310e6dce69fa353fe128",
    "illumination/LPSR_65N_240M_201608.LBL": "0aaa186b0ff5003a32c7e6db18bd742a11c8ca22d18a8fdb82b47bd719b609a0",
    "illumination/LPSR_65S_240M_201608.IMG": "81eac93985118aefbe9ff5cc605018ae51708ebdcd9b80c895bdc14bfc0bb1ff",
    "illumination/LPSR_65S_240M_201608.LBL": "e546d988e903c64a77da73d28ea8058caf10953366f6d133d137911c1d9c5a3d",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def label_integer(text: str, key: str) -> int:
    matches = re.findall(rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)", text)
    if len(matches) != 1:
        raise ValueError(f"expected one integer {key} in detached label")
    return int(matches[0])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=CANONICAL_ROOT)
    args = parser.parse_args()
    package = args.data_root.resolve() / "lola"
    raw = package / "raw"
    records: list[dict[str, object]] = []
    failures: list[str] = []
    failures_by_dataset: dict[str, list[str]] = {"D01": [], "D02": []}
    for relative, expected in EXPECTED.items():
        dataset = "D02" if relative.startswith("illumination/") else "D01"
        path = raw / relative
        if not path.is_file():
            message = f"missing immutable raw file: {path}"
            failures.append(message)
            failures_by_dataset[dataset].append(message)
            continue
        actual = sha256(path)
        records.append({"path": f"raw/{relative}", "bytes": path.stat().st_size,
                        "sha256": actual, "expected_sha256": expected,
                        "hash_matches": actual == expected})
        if actual != expected:
            message = f"SHA-256 mismatch: {path}"
            failures.append(message)
            failures_by_dataset[dataset].append(message)

    # Validate every available detached binary image from its own native label
    # instead of assuming all products share LDEM dimensions.
    for label in sorted(raw.rglob("*.LBL")):
        relative_label = label.relative_to(raw).as_posix()
        dataset = "D02" if relative_label.startswith("illumination/") else "D01"
        image = label.with_suffix(".IMG")
        if not image.is_file():
            message = f"label has no matching IMG: {label}"
            failures.append(message)
            failures_by_dataset[dataset].append(message)
            continue
        text = label.read_text(encoding="ascii")
        lines = label_integer(text, "LINES")
        samples = label_integer(text, "LINE_SAMPLES")
        bits = label_integer(text, "SAMPLE_BITS")
        expected_bytes = lines * samples * bits // 8
        if bits != 16 or image.stat().st_size != expected_bytes:
            message = (
                f"binary size/sample contract mismatch: {image} "
                f"({image.stat().st_size} != {expected_bytes})")
            failures.append(message)
            failures_by_dataset[dataset].append(message)

    d01_names = {"LDEM_4.IMG", "LDEM_4.LBL"}
    d01_ok = not failures_by_dataset["D01"] and all(
        any(row["path"] == f"raw/{name}" and row["hash_matches"] for row in records)
        for name in d01_names
    )
    # D02 raw acquisition is complete for both AVGVISIB and LPSR, but
    # independent processed sampling/projection/edge QA is intentionally not
    # declared READY here.
    d02_raw_ok = all(any(row["path"] == f"raw/{name}" and row["hash_matches"]
                        for row in records)
                    for name in EXPECTED if name.startswith("illumination/"))

    provenance = {
        "schema": "moon-validation-lola-provenance-v1",
        "generated_by": "AMPS/srcMoon/test/data_processing/qualify_lola_package.py",
        "datasets": {
            "D01": {
                "title": "LOLA LDEM_4 global topography",
                "source": "NASA PDS Geosciences LOLA GDR",
                "source_urls": [
                    "https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.IMG",
                    "https://imbrium.mit.edu/DATA/LOLA_GDR/CYLINDRICAL/IMG/LDEM_4.LBL",
                ],
                "product_id": "LDEM_4",
                "frame": "MEAN EARTH/POLAR AXIS OF DE421",
                "status": "READY" if d01_ok else "ERROR",
            },
            "D02": {
                "title": "LOLA polar average visibility reference",
                "source": "NASA PDS/PGDA LOLA polar illumination",
                "source_directory": "https://imbrium.mit.edu/EXTRAS/ILLUMINATION/IMG/",
                "documentation_url": "https://pgda.gsfc.nasa.gov/products/69",
                "products": ["AVGVISIB_65N_240M_201608",
                             "AVGVISIB_65S_240M_201608"],
                "status": "PARTIAL" if d02_raw_ok else "ERROR",
                "limitation": "independent processed comparison table and edge QA not yet frozen",
            },
        },
        "raw_files": records,
    }
    qa = {
        "schema": "moon-validation-lola-qa-v1",
        "overall_status": "ERROR" if failures else "PARTIAL",
        "dataset_status": {"D01": "READY" if d01_ok else "ERROR",
                           "D02": "PARTIAL" if d02_raw_ok else "ERROR"},
        "checks": {"frozen_hashes": not any("SHA-256" in item for item in failures),
                   "label_image_byte_contracts": not any("binary size" in item for item in failures)},
        "failures": failures,
    }
    readme = """# LOLA validation package (D01/D02)

`raw/` is immutable. Run this deterministic qualification step from the AMPS
root:

```bash
python3 -B srcMoon/test/data_processing/qualify_lola_package.py
```

D01 is READY only when its frozen hashes and detached-label binary contract
pass. U08 performs the independent native-pixel and production geometry
checks; I14 remains predecessor-gated until I13 passes. D02 remains PARTIAL
after raw acquisition until an independent, reviewed spatial/temporal
comparison product and edge-selection QA are generated for I15/I20. Directory
presence alone is not readiness. The native labels control units, projection,
scale, offset, and body-fixed frame.
"""

    package.mkdir(parents=True, exist_ok=True)
    (package / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (package / "qa_report.json").write_text(
        json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (package / "raw_inventory.json").write_text(
        json.dumps({"schema": "moon-validation-raw-inventory-v1",
                    "files": records}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    (package / "qa_report.txt").write_text(
        "overall_status: " + qa["overall_status"] + "\n"
        + "D01: " + qa["dataset_status"]["D01"] + "\n"
        + "D02: " + qa["dataset_status"]["D02"] + "\n"
        + ("failures: none\n" if not failures else "failures:\n- " + "\n- ".join(failures) + "\n"),
        encoding="utf-8")
    (package / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps(qa, indent=2, sort_keys=True))
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
