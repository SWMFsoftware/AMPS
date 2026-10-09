#!/usr/bin/env python3
"""Enable AMPS' runtime particle-mover dispatch in generated picGlobal.dfn.

Generic PIC keeps legacy compile-time mover selection enabled by default.
srcSEP3D deliberately opts into the function-pointer path, then registers its
one AMPS adapter callback during ``amps_init()``.  The callback reads the
already validated immutable run configuration and dispatches to the selected
Parker or focused-transport implementation.

This idempotent hook operates only on the generated definition header.  AMPS'
configuration scripts may append input-derived overrides after the template's
include guard, so the srcSEP3D selection is appended last and its check always
examines the effective (last) definition.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys


FINAL_BEGIN = "// BEGIN srcSEP3D final mover selection (R01)"
FINAL_END = "// END srcSEP3D final mover selection (R01)"
FINAL_SELECTION = f"""{FINAL_BEGIN}
// ampsConfig.pl appends input-selected macro overrides after the template's
// include guard. Put the application authority after those generated lines.
// The legacy mover macro remains intact for compatibility and diagnostics;
// generic PIC's dispatcher bypasses it only for this opted-in application.
#undef _PIC_PARTICLE_MOVER_LEGACY_SETTINGS_
#define _PIC_PARTICLE_MOVER_LEGACY_SETTINGS_ _PIC_MODE_OFF_
{FINAL_END}
"""


def _without_final_selection(text: str) -> str:
    pattern = re.compile(
        rf"\n?{re.escape(FINAL_BEGIN)}.*?{re.escape(FINAL_END)}\n?",
        re.DOTALL,
    )
    return pattern.sub("\n", text)


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
    if not definitions[-1].endswith("_PIC_MODE_OFF_"):
        raise RuntimeError(
            "the effective particle mover dispatch mode is not pointer mode: "
            f"{definitions[-1]}"
        )


def install(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    text = _without_final_selection(text)
    if "_PIC_PARTICLE_MOVER_LEGACY_SETTINGS_" not in text:
        raise RuntimeError(
            "generated PIC definitions do not provide the runtime mover "
            "dispatch switch; regenerate from the current source tree"
        )
    # This must be the final definition in the complete generated file.
    text = text.rstrip() + "\n\n" + FINAL_SELECTION
    path.write_text(text, encoding="utf-8")
    check(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pic_global", type=Path,
                        help="configured build/pic/picGlobal.dfn")
    parser.add_argument(
        "--check", action="store_true",
        help="verify that the effective final definition enables pointer mode")
    args = parser.parse_args()
    try:
        if args.check:
            check(args.pic_global)
        else:
            install(args.pic_global)
    except (OSError, RuntimeError) as exc:
        print(f"install_mover_hook.py: {exc}", file=sys.stderr)
        return 2
    action = "verified" if args.check else "enabled"
    print(f"[srcSEP3D] {action} runtime mover dispatch in {args.pic_global}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
