#!/usr/bin/env python3
"""Require byte-identical legacy model output before and after Stage 1."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def files(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): digest(path)
        for path in sorted(root.glob("*.dat"))
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    baseline = files(args.baseline)
    candidate = files(args.candidate)
    if not baseline:
        raise AssertionError("legacy baseline contains no .dat files")
    if baseline != candidate:
        missing = sorted(baseline.keys() - candidate.keys())
        extra = sorted(candidate.keys() - baseline.keys())
        changed = sorted(
            name for name in baseline.keys() & candidate.keys()
            if baseline[name] != candidate[name]
        )
        raise AssertionError(
            f"legacy output differs: missing={missing} extra={extra} changed={changed}"
        )
    print(f"PASS Stage-1 legacy equivalence files={len(baseline)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
