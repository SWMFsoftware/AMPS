#!/usr/bin/env python3
"""Enforce the dependency boundary of the stand-alone SCCM physics library."""

from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(r"(?i)(\bmpi\b|mpi\.h|tecplot|\bpic(?:::|/|\.h)|amps|srcsep|swcme)")

def fail(message: str) -> int:
    print(f"[ARCHSCCM01] FAIL {message}", file=sys.stderr)
    return 1

def main() -> int:
    # Catch compile-time dependencies at their source with a file/line result.
    for folder in (ROOT / "include", ROOT / "src"):
        for path in sorted(folder.rglob("*")):
            if path.suffix not in {".h", ".hpp", ".cpp", ".inc"}:
                continue
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if line.lstrip().startswith("#include") and FORBIDDEN.search(line):
                    return fail(f"forbidden include {path.relative_to(ROOT)}:{number}: {line}")
    common = ROOT.parent / "sep_common"
    for path in sorted(common.glob("sep_status*")):
        if "sep_coronal_cme" in path.read_text(encoding="utf-8"):
            return fail(f"neutral dependency points upward: {path.name}")
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
#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_coronal_cme/mhd_jump_solver.h"
#include "sep_coronal_cme/model_configuration.h"
#include "sep_coronal_cme/source_surface_coupling.h"
#include "sep_coronal_cme/turbulence_transport.h"
int main() {
  SEP::CoronalCME::ModelConfiguration c;
  SEP::CoronalCME::MhdPrimitiveState plasma;
  SEP::CoronalCME::DirectionalWaveState wave;
  SEP::CoronalCME::SurfacePatch patch;
  SEP::CoronalCME::LongitudeMap map;
  return c.schemaVersion + static_cast<int>(plasma.massDensityKgM3 +
      wave.totalJPerM3 + patch.areaM2 + map.forwardJacobian) - 1;
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
