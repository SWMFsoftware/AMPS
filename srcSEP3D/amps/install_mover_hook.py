#!/usr/bin/env python3
"""Verify that AMPS' runtime particle-mover dispatch is enabled for srcSEP3D.

Generic PIC keeps legacy compile-time mover selection enabled by default
(``_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ = _PIC_MODE_ON_`` in
src/pic/picGlobal.dfn).  srcSEP3D opts into the function-pointer path from its
AMPS deck: input/sep3d.input contains

    define _PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ _PIC_MODE_OFF_

and ampsConfig.pl appends ``#undef``/``#define`` for it to the generated
build/pic/picGlobal.dfn.  ``amps_init()`` then registers the one AMPS adapter
callback (Parker or focused transport, as selected by the srcSEP3D input
parser) through ``PIC::Mover::SetUserDefinedParticleMover``.

This script no longer edits any file.  It checks the *effective* (last)
definition in the generated header and fails if the deck did not select
pointer mode.  It is called by ``make -C srcSEP3D audit-production-symbols``.
The historical file name is kept because build and test scripts reference it;
``--check`` is accepted for compatibility and is the only mode.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


def check(path: Path) -> None:
    """Verify the effective (last) dispatch-mode definition."""
    text = path.read_text(encoding="utf-8")
    definitions = re.findall(
        r"^\s*#define\s+_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_\s+[^\n]+$",
        text,
        flags=re.MULTILINE,
    )
    if not definitions:
        raise RuntimeError("particle mover dispatch mode has no definition")
    if not definitions[-1].rstrip().endswith("_PIC_MODE_OFF_"):
        raise RuntimeError(
            "the effective particle mover dispatch mode is not pointer mode: "
            f"{definitions[-1]}; add 'define "
            "_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ _PIC_MODE_OFF_' to the "
            "#General block of the AMPS input deck and reconfigure"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pic_global", type=Path,
                        help="configured build/pic/picGlobal.dfn")
    parser.add_argument(
        "--check", action="store_true",
        help="accepted for compatibility; checking is the only mode")
    args = parser.parse_args()
    try:
        check(args.pic_global)
    except (OSError, RuntimeError) as exc:
        print(f"install_mover_hook.py: {exc}", file=sys.stderr)
        return 2
    print(f"[srcSEP3D] verified runtime mover dispatch in {args.pic_global}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
