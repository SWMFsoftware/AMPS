#!/usr/bin/env python3
"""Download the frozen CDAW height-time files used by VP04."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import urllib.request


# Each URL is a digital leading-edge measurement file linked by the CDAW CME
# catalog.  Byte counts and SHA-256 digests freeze the exact manual traces used
# in this campaign; a catalog correction must be reviewed as a data revision.
FILES = (
    ("20130411.072406.w360h.v0861.p085g.yht", "2013_04", 2100, "fa0643e479fdc83d3f92a0640c69a0d0f9c74dc27b0b9133e36bf4de663e43f1"),
    ("20130522.132550.w360h.v1466.p287g.yht", "2013_05", 1694, "8e042fc21cfa656363cfa8127a4de9b44bb953fa4facace9d54df2614aa3bc36"),
    ("20140106.080005.w360h.v1402.p274g.yht", "2014_01", 1578, "e4e0cf247549e1213f1907b5b31be350b99c8bdfbe99e52b3a30dd0a09b8e187"),
    ("20140107.182405.w360h.v1830.p231g.yht", "2014_01", 1520, "a1c4c5db1f2aadab4f5f0a634f8553d7f5adfcd0cf3b19c2d26dfc27ef474a69"),
    ("20170906.122405.w360h.v1571.p201g.yht", "2017_09", 1578, "8b5905fa8acb2f8206da7885b43d8070b6a1491ae08c1b7c1b7bf4f851527445"),
    ("20170910.160005.w360h.v3163.p262g.yht", "2017_09", 1230, "058f28a870129fea6c9adc310ccb0ad37368006ca9d32fa8c76225f9f8e334bd"),
)


def url(name: str, month: str) -> str:
    return f"https://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/{month}/yht/{name}"


def verify(path: Path, expected_size: int, expected_sha256: str) -> None:
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if len(data) != expected_size or digest != expected_sha256:
        raise RuntimeError(f"content mismatch for {path}: bytes={len(data)}, sha256={digest}")


def download(source_url: str, destination: Path) -> None:
    """Use a same-directory temporary file so final rename is atomic."""

    temporary = destination.with_name(destination.name + ".part")
    request = urllib.request.Request(source_url, headers={"User-Agent": "SWCME-VP04/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
            output.flush(); os.fsync(output.fileno())
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "data/raw")
    args = parser.parse_args(); args.data_dir.mkdir(parents=True, exist_ok=True)
    receipt = []
    for name, month, size, digest in FILES:
        path = args.data_dir / name
        if not path.exists():
            print(f"Downloading {name} ...")
            download(url(name, month), path)
        verify(path, size, digest)
        receipt.append({"name": name, "bytes": size, "sha256": digest, "source_url": url(name, month)})
        print(f"  PASS {digest}  {path}")
    (args.data_dir / "download_receipt.json").write_text(json.dumps({
        "schema_version": 1, "validation_id": "VP04",
        "verified_utc": datetime.now(timezone.utc).isoformat(), "files": receipt,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP04 data error: {error}", file=sys.stderr)
        raise SystemExit(2)
