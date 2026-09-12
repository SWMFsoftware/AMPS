#!/usr/bin/env python3
"""Acquire and verify the immutable Helios archive used by VP03.

The VP03 package owns this entry point so it remains runnable independently.
It may reuse a byte-identical VP01 or VP02 cache through a hard link, but it
never trusts that cache until all frozen size, MD5, and SHA-256 checks pass.
"""

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

FILE = {
    "name": "corefit.gz",
    "url": "https://zenodo.org/api/records/1009506/files/corefit.gz/content",
    "bytes": 131_088_147,
    "md5": "3706615fd835c9cbdd4e5b59a8431e30",
    "sha256": "5659d9408fe3f0c75507ab416d58cd51b91ccb243fbad9dd488721b7dd437a08",
}


def md5_hasher():
    """Support vendor Python builds that predate ``usedforsecurity``."""

    try:
        return hashlib.md5(usedforsecurity=False)
    except TypeError:
        return hashlib.md5()


def digests(path: Path) -> tuple[str, str, int]:
    md5, sha256, size = md5_hasher(), hashlib.sha256(), 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
            size += len(block)
    return md5.hexdigest(), sha256.hexdigest(), size


def verify(path: Path) -> None:
    """Fail closed when an input differs from the published archive."""

    expected = (str(FILE["md5"]), str(FILE["sha256"]), int(FILE["bytes"]))
    observed = digests(path)
    if observed != expected:
        raise RuntimeError(f"checksum mismatch for {path}: observed={observed}")


def acquire(destination: Path) -> None:
    """Install verified bytes atomically, optionally from a sibling cache."""

    validation_dir = Path(__file__).resolve().parents[1]
    caches = [validation_dir / case / "data/raw/corefit.gz" for case in ("vp01", "vp02")]
    temporary = destination.with_name(destination.name + ".part")
    try:
        source = next((candidate for candidate in caches if candidate.is_file()), None)
        if source is not None:
            verify(source)
            try:
                os.link(source, temporary)
            except OSError:
                shutil.copyfile(source, temporary)
        else:
            request = urllib.request.Request(str(FILE["url"]), headers={"User-Agent": "SWCME-VP03/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Download and verify VP03 data")
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "data/raw")
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)
    destination = args.data_dir / str(FILE["name"])
    if not destination.exists():
        acquire(destination)
    verify(destination)
    md5, sha256, size = digests(destination)
    (args.data_dir / "download_receipt.json").write_text(json.dumps({
        "schema_version": 1, "validation_id": "VP03",
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "file": {"name": destination.name, "bytes": size, "md5": md5, "sha256": sha256},
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS {sha256}  {destination}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP03 data error: {error}", file=sys.stderr)
        raise SystemExit(2)
