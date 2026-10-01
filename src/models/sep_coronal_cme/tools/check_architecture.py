#!/usr/bin/env python3
"""Enforce the dependency boundary of the stand-alone SCCM physics library."""

from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(r"(?i)(\bmpi\b|mpi\.h|tecplot|\bpic(?:::|/|\.h)|amps|srcsep|swcme)")
# AMPS builds may leave sep_*.o next to the neutral module's sources. Only
# C/C++ source/header files belong to the text audit; archive dependencies
# still receive the separate nm check below. Never decode object/archive bytes
# as text or mask invalid UTF-8 in a genuine source with errors='ignore'.
SOURCE_SUFFIXES = {".h", ".hpp", ".hh", ".hxx", ".c", ".cc", ".cpp",
                   ".cxx", ".inc", ".ipp", ".tpp", ".ixx", ".cppm"}


def source_files(folder: Path, pattern: str = "*"):
    """Select regular source files, including future nested neutral sources."""
    return (path for path in sorted(folder.rglob(pattern))
            if path.is_file() and path.suffix.lower() in SOURCE_SUFFIXES)


def source_dependency_error(root: Path):
    """Return a source-boundary violation, or None when the text audit closes.

    Keeping this filesystem audit separate from linking permits regressions
    with binary build artifacts and forbidden sources in isolated temporary
    trees. The production gate still performs all symbol and public ABI checks.
    """
    common = root.parent / "sep_common"
    scans = ((root/"include", "*", False), (root/"src", "*", False),
             (common, "sep_*", True))
    for folder, pattern, neutral in scans:
        for path in source_files(folder, pattern):
            label = str(path.relative_to(common if neutral else root))
            if neutral:
                label = "sep_common/"+label
            try:
                text = path.read_text(encoding="utf-8")
            except (OSError, UnicodeError) as error:
                return f"cannot read UTF-8 source {label}: {error}"
            if neutral:
                if "sep_coronal_cme" in text:
                    return f"neutral dependency points upward: {label}"
            else:
                for number, line in enumerate(text.splitlines(), 1):
                    if line.lstrip().startswith("#include") and FORBIDDEN.search(line):
                        return f"forbidden include {label}:{number}: {line}"
    return None

def fail(message: str) -> int:
    print(f"[ARCHSCCM01] FAIL {message}", file=sys.stderr)
    return 1

def main() -> int:
    # Catch compile-time dependencies at their source with a file/line result.
    violation = source_dependency_error(ROOT)
    if violation:
        return fail(violation)
    common = ROOT.parent / "sep_common"
    archive = ROOT / "build" / "libsep_coronal_cme.a"
    if not archive.exists():
        return fail("library archive is absent; run make lib first")
    undefined = subprocess.run(["nm", "-u", str(archive)], check=True,
        text=True, stdout=subprocess.PIPE).stdout
    for line in undefined.splitlines():
        if FORBIDDEN.search(line):
            return fail(f"forbidden unresolved symbol: {line.strip()}")
    # A clean external translation unit detects application-only types leaked
    # through the public ABI.
    consumer = '''#include "sep_coronal_cme/configuration_parser.h"
#include "sep_coronal_cme/discontinuity_transport.h"
#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_coronal_cme/field_line_reduction.h"
#include "sep_coronal_cme/mhd_jump_solver.h"
#include "sep_coronal_cme/model_configuration.h"
#include "sep_coronal_cme/particle_source.h"
#include "sep_coronal_cme/runtime_integration.h"
#include "sep_coronal_cme/shock_provider.h"
#include "sep_coronal_cme/source_surface_coupling.h"
#include "sep_coronal_cme/turbulence_transport.h"
#include "sep_field_line_bundle_io.h"
#include "sep_field_line_exchange.h"
int main() {
  SEP::CoronalCME::ModelConfiguration c;
  SEP::CoronalCME::MhdPrimitiveState plasma;
  SEP::CoronalCME::DirectionalWaveState wave;
  SEP::CoronalCME::SurfacePatch patch;
  SEP::CoronalCME::LongitudeMap map;
  SEP::FieldLine::FieldLineSet lines;
  return c.schemaVersion + static_cast<int>(plasma.massDensityKgM3 +
      wave.totalJPerM3 + patch.areaM2 + map.forwardJacobian +
      lines.lines.size()) - 1;
}
'''
    with tempfile.TemporaryDirectory(prefix="sccm-architecture-") as directory:
        source = Path(directory) / "consumer.cpp"
        source.write_text(consumer, encoding="utf-8")
        command = ["g++", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
            "-Werror", "-I", str(ROOT / "include"), "-I", str(common),
            str(source), str(archive), "-o", str(Path(directory) / "consumer")]
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        if result.returncode:
            return fail("public-header consumer did not compile:\n" + result.stdout)
    print("[ARCHSCCM01] PASS dependency-light C++17 archive and public headers")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
