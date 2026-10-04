#!/usr/bin/env python3
"""Dependency/public-header audit for the shared Corona--SWCME model."""

from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_INCLUDE = re.compile(
    r"^\s*#\s*include\s*[<\"](?:mpi(?:\.h)?|pic(?:\.h)?|SEP3D(?:\.h)?|amps)",
    re.IGNORECASE,
)


def fail(message):
    print("[ARCHCSWC01] FAIL " + message)
    return 1


def main():
    for folder in (ROOT / "include", ROOT / "src"):
        for source in sorted(folder.rglob("*")):
            if not source.is_file() or source.suffix not in {".h", ".hpp", ".cpp"}:
                continue
            text = source.read_text(encoding="utf-8")
            for number, line in enumerate(text.splitlines(), 1):
                if FORBIDDEN_INCLUDE.search(line):
                    return fail(f"forbidden application/MPI include {source}:{number}")

    with tempfile.TemporaryDirectory(prefix="corona-swcme-public-") as temporary:
        source = Path(temporary) / "consumer.cpp"
        output = Path(temporary) / "consumer.o"
        source.write_text(
            '#include "sep_corona_swcme/cme_event.h"\n'
            '#include "sep_corona_swcme/ambient_state.h"\n'
            '#include "sep_corona_swcme/surface_shock.h"\n'
            '#include "sep_corona_swcme/sheath_model.h"\n'
            "int main() {\n"
            "  SEP::CoronaSwcme::AmbientPrimitive primitive;\n"
            "  SEP::CoronaSwcme::EventSupport support;\n"
            "  return primitive.magneticSector + (support.endS > 0);\n"
            "}\n",
            encoding="utf-8",
        )
        command = [
            "g++", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-I" + str(ROOT / "include"),
            "-I" + str(ROOT.parent / "sep_coronal_cme" / "include"),
            "-I" + str(ROOT.parent / "sep_common"),
            "-c", str(source), "-o", str(output),
        ]
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            return fail("public-header consumer did not compile:\n" + result.stdout)

    archive = ROOT / "build" / "libsep_corona_swcme.a"
    result = subprocess.run(["nm", "-u", str(archive)], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            check=False)
    if result.returncode:
        return fail("cannot inspect shared archive:\n" + result.stdout)
    forbidden = [line for line in result.stdout.splitlines()
                 if "MPI_" in line or "PIC::" in line or "_ZN3PIC" in line]
    if forbidden:
        return fail("archive has application/MPI symbols: " + "; ".join(forbidden))
    print("[ARCHCSWC01] PASS dependency-light public headers and archive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
