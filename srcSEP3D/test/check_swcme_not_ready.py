#!/usr/bin/env python3
"""srcSEP3D: native check that every SWCME example deck stops as "not ready".

User decision 2026-10-12: the canonical SWCME model is not ready to be used in
srcSEP3D.  Any '#begin physics: swcme' block must stop the run, before mesh
allocation, with "the SWCME model is not ready yet to be used" at the file
and line of that block.  The six SWCME example decks are kept unchanged as
rejected fixtures; this runner launches the built ./amps on each of them
(initialization-only, one rank) and checks

  SWCMENR.<deck>  non-zero exit status, the typed diagnosis, and the
                  '<deck>:<line>' of the deck's first '#begin physics: swcme'
                  (found independently from the deck text).

It replaces the run_all steps that used to run these decks (init-full,
init-tube, sccm3d, sep-corona-suite, coupled-runner, swcme-coupled) until the
SWCME model is released.  Output lines '[ID] STATUS message' are read by
tools/sep_test_orchestrator.py (read_bracket_text); the exit status is 1 when
any check fails.

Usage (from the AMPS root, after the native srcSEP3D build):
  python3 srcSEP3D/test/check_swcme_not_ready.py --amps ./amps \\
      --output-dir test_output/swcme-not-ready
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DECKS = [
    "srcSEP3D/examples/sep3d_analytic_parker.in",
    "srcSEP3D/examples/sep3d_analytic_parker_active_tube.in",
    "srcSEP3D/examples/sep3d_analytic_parker_corner_sphere.in",
    "srcSEP3D/examples/sep3d_swcme_20rs_1au.in",
    "srcSEP3D/examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in",
    "srcSEP3D/examples/sep3d_swcme_sse_mesh_background_20rs_1au.in",
]
DIAGNOSIS = "the SWCME model is not ready yet to be used"


def swcme_line(deck: Path) -> int:
    for number, line in enumerate(deck.read_text().splitlines(), start=1):
        if line.startswith("#begin physics: swcme"):
            return number
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--amps", default=str(ROOT / "amps"))
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output-dir", default=str(ROOT / "test_output/swcme-not-ready"))
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    failures = 0
    for relative in DECKS:
        deck = ROOT / relative
        test_id = "SWCMENR." + deck.stem
        line = swcme_line(deck)
        log = out / (deck.stem + ".log")
        command = [args.mpiexec, "-n", "1", str(Path(args.amps).resolve()),
                   "--input", str(deck), "--initialization-only",
                   "--initialization-output-dir", str(out / deck.stem)]
        try:
            with log.open("w") as stream:
                code = subprocess.run(command, cwd=out, stdout=stream,
                                      stderr=subprocess.STDOUT,
                                      timeout=args.timeout).returncode
        except subprocess.TimeoutExpired:
            print(f"[{test_id}] ERROR timeout; log {log}", flush=True)
            failures += 1
            continue
        text = log.read_text(errors="replace")
        where = f"{deck.name}:{line}:"
        if line and code != 0 and DIAGNOSIS in text and where in text:
            print(f"[{test_id}] PASS stopped at {where} with the not-ready "
                  f"diagnosis (exit {code})", flush=True)
        else:
            print(f"[{test_id}] FAIL exit {code}; expected '{DIAGNOSIS}' at "
                  f"{where}; log {log}", flush=True)
            failures += 1
    print(f"SWCME not-ready summary: {len(DECKS) - failures} passed, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
