#!/usr/bin/env python3
"""Generate the dedicated I08 Na/Na+ no-SPICE AMPS input.

The I08 fixture must follow the production Moon configuration byte-for-byte
except for the two declared controls needed to exercise the ion mover:
`SpeciesList=Na,NA_PLUS` and `SPICE=off`.  Generating it from the canonical
input avoids maintaining a second hand-edited physical model.  Exact-match and
single-replacement checks fail closed when the production input changes, so a
future edit cannot silently alter the fixture contract.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"expected exactly one {old!r}, found {count}")
    return text.replace(old, new, 1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    text = arguments.source.read_text(encoding="utf-8")
    text = replace_once(text, "SpeciesList=Na  !, NA_Plus !lkdsjflksdj",
                        "SpeciesList=Na,NA_PLUS  ! dedicated I08 ion-mover fixture")
    text = replace_once(text, "  SPICE=on !on, off",
                        "  SPICE=off ! I08 isolates uniform-field Lorentz motion")
    text = replace_once(text, "  #include moon.input.spicekernels",
                        "  ! No SPICE kernel include: I08 is explicitly no-SPICE.")
    arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
