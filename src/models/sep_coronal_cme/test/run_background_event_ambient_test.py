#!/usr/bin/env python3
"""Build/run the canonical ambient binding without adding a core dependency.

Only this integration executable links the existing composite ambient factory.
The coronal archive and its public headers remain standalone, particle-free,
and free of composite/PIC/MPI references. Physical source and regional closures
are not enabled or qualified by this background-only interoperability test.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    built = subprocess.run(["make", "-j8", "build/background_event_ambient_tests"], cwd=ROOT)
    if built.returncode:
        raise SystemExit(built.returncode)
    raise SystemExit(subprocess.run([str(ROOT / "build/background_event_ambient_tests")], cwd=ROOT).returncode)
