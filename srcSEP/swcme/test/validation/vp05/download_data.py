#!/usr/bin/env python3
"""Download and verify the frozen HELIO4CAST LineupCAT table for VP05."""

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
    "name": "HELIO4CAST_multipoint_v30.csv",
    "url": "https://helioforecast.space/static/sync/lineups/HELIO4CAST_multipoint_v30.csv",
    "bytes": 50_510,
    "sha256": "d55bdfce703606d68c90681d47486807b1d51582dff0f757e33660186c9183b9",
}


def verify(path: Path) -> None:
    """Pin the campaign to exact v3.0 bytes rather than a mutable URL label."""

    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if len(data) != FILE["bytes"] or digest != FILE["sha256"]:
        raise RuntimeError(f"content mismatch for {path}: bytes={len(data)}, sha256={digest}")


def download(destination: Path) -> None:
    temporary = destination.with_name(destination.name + ".part")
    request = urllib.request.Request(str(FILE["url"]), headers={"User-Agent": "SWCME-VP05/1.0"})
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
    destination = args.data_dir / str(FILE["name"])
    if not destination.exists():
        download(destination)
    verify(destination)
    (args.data_dir / "download_receipt.json").write_text(json.dumps({
        "schema_version": 1, "validation_id": "VP05",
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "file": {"name": destination.name, "bytes": FILE["bytes"], "sha256": FILE["sha256"], "source_url": FILE["url"]},
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS {FILE['sha256']}  {destination}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP05 data error: {error}", file=sys.stderr)
        raise SystemExit(2)
