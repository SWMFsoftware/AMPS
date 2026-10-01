#!/usr/bin/env python3
"""Validate the required source members of an AMPS SEP overlay.

The release archive is an overlay, not a build-tree snapshot.  This small gate
keeps the AMPS-core splitter fix in the same auditable source set as srcSEP3D
and the shared coronal-CME model, while rejecting common generated products.
It intentionally performs no copying and never modifies the inspected tree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable, Tuple


# BLDL3D09 reads this literal tuple with ast.literal_eval.  A required member
# differs from an allowlist entry: its absence must fail packaging because an
# old destination copy would otherwise remain active after overlay extraction.
TOP_LEVEL_REQUIRED: Tuple[str, ...] = (
    "srcSEP3D/core/domain_geometry.h",
    "srcSEP3D/core/domain_geometry.cpp",
    "srcSEP3D/test/individual-test/domain_geometry_probe.cpp",
    "srcSEP3D/test/individual-test/output_boundary_probe.cpp",
    "src/pic/pic.h",
    "src/pic/ecsim/domain_bc.h",
    "tools/test_pic_header_guards.py",
    "srcSEP3D/examples/sep3d_analytic_parker_corner_sphere.in",
    "src/pic/pic_particle_spliting.cpp",
    "src/models/sep_coronal_cme/README.md",
    "src/models/sep_coronal_cme/makefile",
    "srcSEP3D/main.cpp",
    "srcSEP3D/main_lib.cpp",
    "srcSEP3D/makefile",
    "srcSEP3D/validation/coronal_cme_application_test.cpp",
    "srcSEP3D/validation/coronal_cme_application_test.h",
    "src/models/sep_common/sep_coherent_transport.h",
    "src/models/sep_common/sep_coherent_transport.cpp",
    "src/models/sep_coronal_cme/include/sep_coronal_cme/research_extensions.h",
    "src/models/sep_coronal_cme/src/research_extensions.cpp",
    "src/models/sep_coronal_cme/test/test_stage13.py",
    "src/models/sep_coronal_cme/test/test_stage14.py",
    "src/models/sep_coronal_cme/tools/release/qualification.py",
    "src/models/sep_coronal_cme/tools/qualify_release.py",
    "src/models/sep_coronal_cme/tools/research_stage14.py",
    "src/models/sep_coronal_cme/docs/STAGE13_RELEASE_QUALIFICATION.md",
    "src/models/sep_coronal_cme/docs/STAGE14_RESEARCH_EXTENSIONS.md",
)

FORBIDDEN_DIRECTORY_NAMES: Tuple[str, ...] = (
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
)
FORBIDDEN_SUFFIXES: Tuple[str, ...] = (
    ".o",
    ".a",
    ".pyc",
    ".pyo",
)


def generated_members(root: Path) -> Iterable[str]:
    """Yield repository-relative generated members in deterministic order."""
    found = []
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if any(part in FORBIDDEN_DIRECTORY_NAMES for part in relative.parts):
            found.append(relative.as_posix())
        elif path.is_file() and path.suffix in FORBIDDEN_SUFFIXES:
            found.append(relative.as_posix())
        elif path.is_dir() and path.name == "build" and (
                "sep_coronal_cme" in relative.parts or
                "srcSEP3D" in relative.parts):
            found.append(relative.as_posix() + "/")
    return sorted(set(found))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="check required and generated members in a SEP overlay")
    parser.add_argument("root", nargs="?", type=Path, default=Path.cwd(),
                        help="AMPS overlay root (default: current directory)")
    parser.add_argument(
        "--allow-generated", action="store_true",
        help="check required sources only; useful for a developer build tree")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        print("ERROR: overlay root is not a directory: {}".format(root),
              file=sys.stderr)
        return 2

    missing = [member for member in TOP_LEVEL_REQUIRED
               if not (root / member).is_file()]
    generated = [] if args.allow_generated else list(generated_members(root))
    if missing:
        print("ERROR: missing required release members:", file=sys.stderr)
        for member in missing:
            print("  {}".format(member), file=sys.stderr)
    if generated:
        print("ERROR: generated members must not enter the source archive:",
              file=sys.stderr)
        for member in generated:
            print("  {}".format(member), file=sys.stderr)
    if missing or generated:
        return 1
    print("SEP overlay hygiene PASS: {} required members are present".format(
        len(TOP_LEVEL_REQUIRED)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
