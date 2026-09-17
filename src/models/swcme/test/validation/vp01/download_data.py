#!/usr/bin/env python3
"""Download and checksum the immutable observational inputs for VP01.

The large third-party archives deliberately remain outside source control.  A
fresh checkout can recreate the exact input bytes from Zenodo, and this script
refuses to use either a truncated download or a silently replaced remote file.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import urllib.request


# Values below are copied from the Zenodo record and independently augmented
# with SHA-256 after the first verified download.  MD5 checks the repository's
# published checksum; SHA-256 gives the campaign a stronger content identity.
FILES = (
    {
        "name": "corefit.gz",
        "url": "https://zenodo.org/api/records/1009506/files/corefit.gz/content",
        "bytes": 131_088_147,
        "md5": "3706615fd835c9cbdd4e5b59a8431e30",
        "sha256": "5659d9408fe3f0c75507ab416d58cd51b91ccb243fbad9dd488721b7dd437a08",
    },
    {
        "name": "source_code.zip",
        "url": "https://zenodo.org/api/records/1009506/files/source_code.zip/content",
        "bytes": 890_360,
        "md5": "83e28b40c4fea3cdc2cff00a5c527a98",
        "sha256": "a666bda15fad59a65f4c4b9a09a617a89802b3ccd63828a78d22b32e6633a841",
    },
)
RECORD_URL = "https://zenodo.org/api/records/1009506"


def md5_hasher():
    """Create an MD5 object on both modern and legacy Python/OpenSSL builds.

    Python builds linked to newer OpenSSL providers accept
    ``usedforsecurity=False``.  That flag correctly describes this use: MD5 is
    retained only to compare with the checksum published by Zenodo, while the
    independently frozen SHA-256 digest is the security-strength content
    identity.  Older vendor Python builds expose ``openssl_md5`` without that
    keyword and raise ``TypeError`` before hashing any bytes.  Falling back only
    for that signature mismatch preserves support for those systems without
    hiding a real provider or policy error from either call.
    """

    try:
        return hashlib.md5(usedforsecurity=False)
    except TypeError:
        return hashlib.md5()


def digests(path: Path) -> tuple[str, str, int]:
    """Return MD5, SHA-256, and byte count without loading an archive in RAM."""

    md5 = md5_hasher()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
            size += len(block)
    return md5.hexdigest(), sha256.hexdigest(), size


def verify(path: Path, specification: dict[str, object]) -> None:
    """Reject data unless every independent content check is exact."""

    observed_md5, observed_sha256, observed_size = digests(path)
    errors = []
    if observed_size != specification["bytes"]:
        errors.append(f"bytes={observed_size}, expected {specification['bytes']}")
    if observed_md5 != specification["md5"]:
        errors.append(f"md5={observed_md5}, expected {specification['md5']}")
    if observed_sha256 != specification["sha256"]:
        errors.append(
            f"sha256={observed_sha256}, expected {specification['sha256']}"
        )
    if errors:
        raise RuntimeError(f"{path}: " + "; ".join(errors))


def download(url: str, destination: Path) -> None:
    """Download atomically so an interruption never looks like valid input."""

    temporary = destination.with_name(destination.name + ".part")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SWCME-VP01-validation/1.0 (+https://zenodo.org)"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open(
            "wb"
        ) as output:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                output.write(block)
            output.flush()
            os.fsync(output.fileno())
        temporary.replace(destination)
    except BaseException:
        # The temporary pathname is validation-owned and cannot contain a prior
        # accepted input.  Removing it keeps retry behavior deterministic.
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download and verify the Helios inputs used by SWCME VP01"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "raw",
        help="destination directory (default: validation/vp01/data/raw)",
    )
    parser.add_argument(
        "--corefit-only",
        action="store_true",
        help="skip the small author-source archive; corefit.gz is always fetched",
    )
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)

    selected = FILES[:1] if args.corefit_only else FILES
    verified_files = []
    for specification in selected:
        destination = args.data_dir / str(specification["name"])
        if destination.exists():
            print(f"Verifying existing {destination.name} ...")
        else:
            print(f"Downloading {destination.name} from Zenodo ...")
            download(str(specification["url"]), destination)
        verify(destination, specification)
        print(f"  PASS {specification['sha256']}  {destination}")
        observed_md5, observed_sha256, observed_size = digests(destination)
        verified_files.append(
            {
                "name": destination.name,
                "bytes": observed_size,
                "md5": observed_md5,
                "sha256": observed_sha256,
                "source_url": specification["url"],
            }
        )

    # The live record is provenance, not a numerical input, so it is archived
    # for audit but intentionally excluded from the immutable checksum gate.
    record_path = args.data_dir / "zenodo_record_1009506.json"
    download(RECORD_URL, record_path)
    json.loads(record_path.read_text(encoding="utf-8"))
    print(f"  saved current Zenodo metadata: {record_path}")
    receipt = {
        "schema_version": 1,
        "validation_id": "VP01",
        "record_doi": "10.5281/zenodo.1009506",
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "files": verified_files,
    }
    (args.data_dir / "download_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"VP01 data error: {error}", file=sys.stderr)
        raise SystemExit(2)
