#!/usr/bin/env python3
"""Install srcSEP3D's single production mover in generated picGlobal.dfn.

AMPS compiles ``pic_mover.cpp`` before it compiles the application archive, so
the selected callback must be visible in the generated PIC definition header.
This idempotent application hook performs that configuration step explicitly;
it does not edit the canonical source template unless the caller deliberately
passes that path.  Normal use targets ``AMPS/build/pic/picGlobal.dfn`` after
AMPS configuration and before compilation.

The forward declarations are required because picGlobal.dfn is included near
the beginning of pic.h.  Merely replacing the macro (the old prototype's
behavior) leaves ``SEP3D::AMPS::Movers`` undeclared while pic_mover.cpp is
compiled.  The declaration below exactly matches the implementation in
amps_particle_adapter.cpp and therefore also gives the compiler a signature
check at every macro call site.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


DEFAULT = (
    "#define _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,LocalTimeStep,node) "
    "PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder"
    "(ptr,LocalTimeStep,node);"
)
SELECTED = (
    "#define _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,LocalTimeStep,node) "
    "SEP3D::AMPS::Movers::MoveParticle(ptr,LocalTimeStep,node);"
)
DECLARATIONS = """// BEGIN srcSEP3D generated mover declaration (R01)
template <class T> class cTreeNodeAMR;
namespace PIC { namespace Mesh { class cDataBlockAMR; } }
namespace SEP3D { namespace AMPS { namespace Movers {
int MoveParticle(long int ptr, double dtTotal,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
} } }
// END srcSEP3D generated mover declaration (R01)
"""


def install(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    if SELECTED not in text:
        if DEFAULT not in text:
            raise RuntimeError(
                "expected default AMPS mover macro was not found; refusing an "
                "ambiguous configuration edit"
            )
        text = text.replace(DEFAULT, SELECTED, 1)
    if "BEGIN srcSEP3D generated mover declaration" not in text:
        guard = "#define _PIC_GLOBAL_DEFINITIONS_H_\n"
        if guard not in text:
            raise RuntimeError("picGlobal.dfn include guard was not found")
        text = text.replace(guard, guard + "\n" + DECLARATIONS, 1)
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pic_global", type=Path,
                        help="configured build/pic/picGlobal.dfn")
    args = parser.parse_args()
    try:
        install(args.pic_global)
    except (OSError, RuntimeError) as exc:
        print(f"install_mover_hook.py: {exc}", file=sys.stderr)
        return 2
    print(f"[srcSEP3D] selected production mover in {args.pic_global}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
