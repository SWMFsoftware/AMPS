#!/usr/bin/env python3
"""Test 1 of srcHartley2/CHEMISTRY_VARIANTS.md (section 10): rates of the
SA1982_VH OH-loss variant.

An independent Python evaluation of the legacy formula is compared with the
C++ module that AMPS uses (srcHartley2/CometChemistryVariants.cpp), at the six
I3 case velocities of input/i3.in, at grid points, mid-points and both clamps.

The Python side shares no code with AMPS.  It takes OHTAUS from the legacy
Fortran DATA statement recorded verbatim in CHEMISTRY_VARIANTS.md (section 1),
not from the CSV file, so the test also checks that the CSV holds that table.
It reads tau_FUV, S and the ion fraction from the legacy code's documented
values (4.24e5 s, 0.833e5 s, 0.05) and checks that input/i3.conf and the CSV
use the same values.

Usage: test_chemistry_variants.py PATH_TO_test_chemistry_variants_BINARY
"""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "srcHartley2" / "i3" / "reference" / "oh_lifetime_sa1982.csv"
NOTE = ROOT / "srcHartley2" / "CHEMISTRY_VARIANTS.md"
I3_CONF = ROOT / "input" / "i3.conf"
I3_IN = ROOT / "input" / "i3.in"
BINARY: Path | None = None

TAU_FUV_S = 4.24e5   # legacy: taum0 = 4.24d5
SCALE_S = 0.833e5    # legacy: tauoh = ohtaus(ivel)*0.833d5
ION_FRACTION = 0.05  # decision 9 of CHEMISTRY_VARIANTS.md


def legacy_ohtaus() -> list[float]:
    """OHTAUS from the verbatim DATA statement in the design note."""
    text = NOTE.read_text(encoding="utf-8")
    block = text[text.index("DATA OHTAUS/"):text.index("DATA OHLON4/")]
    values = [float(x) for x in re.findall(r"\d+\.\d+", block.split("/", 1)[1])]
    assert len(values) == 61, len(values)
    return values


def reference_rates(velocity: float) -> tuple[float, int, float, float]:
    """Independent evaluation: grid -60 + 2(i-1) km/s, linear, clamped."""
    table = legacy_ohtaus()
    grid = [-60.0 + 2.0 * i for i in range(61)]
    if velocity <= grid[0]:
        ohtaus, clamped = table[0], int(velocity < grid[0])
    elif velocity >= grid[-1]:
        ohtaus, clamped = table[-1], int(velocity > grid[-1])
    else:
        i = int((velocity - grid[0]) // 2.0)
        i = min(i, 59)
        t = (velocity - grid[i]) / (grid[i + 1] - grid[i])
        ohtaus, clamped = table[i] * (1.0 - t) + table[i + 1] * t, 0
    k_d = 1.0 / TAU_FUV_S + 1.0 / (SCALE_S * ohtaus)
    return ohtaus, clamped, k_d, ION_FRACTION * k_d


def case_velocities() -> dict[int, float]:
    text = I3_IN.read_text(encoding="utf-8")
    result = {}
    for match in re.finditer(r"#begin case:\s*(\d+)(.*?)#end case", text, re.S):
        found = re.search(r"brightness\.heliocentric_radial_velocity_km_s\s*=\s*(\S+)",
                          match.group(2))
        result[int(match.group(1))] = float(found.group(1))
    return result


def run_binary(*arguments: str) -> str:
    return subprocess.run([str(BINARY), *arguments], check=True, text=True,
                          capture_output=True).stdout


class ChemistryVariantRates(unittest.TestCase):
    def test_csv_holds_the_legacy_table(self) -> None:
        rows = [line.split(",") for line in TABLE.read_text().splitlines()
                if line[:1] in "-0123456789"]
        self.assertEqual([float(v) for v, _ in rows], [-60.0 + 2 * i for i in range(61)])
        self.assertEqual([float(o) for _, o in rows], legacy_ohtaus())
        self.assertIn(f"tau_fuv_s,{TAU_FUV_S:.2e}".replace("e+05", "e5"), TABLE.read_text())

    def test_build_input_uses_the_legacy_constants(self) -> None:
        text = I3_CONF.read_text(encoding="utf-8")
        record = re.search(r"ChemistryRateTable=OH_TAU_SA1982 \| MODEL=SA1982_OH_LIFETIME"
                           r" \| SCALE_S=(\S+) \| ION_FRACTION=(\S+)", text)
        self.assertIsNotNone(record)
        self.assertEqual(float(record.group(1)), SCALE_S)
        self.assertEqual(float(record.group(2)), ION_FRACTION)

    def test_self_test(self) -> None:
        output = run_binary("--self-test", str(TABLE))
        self.assertIn("PASS test_chemistry_variants failures=0", output, output)

    def test_rates_match_independent_evaluation(self) -> None:
        cases = case_velocities()
        self.assertEqual(sorted(cases), [1, 2, 3, 4, 5, 6])
        velocities = list(cases.values()) + [-60.0, -59.0, -40.0, -39.0, 0.0, 1.0,
                                             58.5, 60.0, -61.0, 75.0, -100.0]
        output = run_binary("--rates", str(TABLE), repr(SCALE_S), repr(ION_FRACTION),
                            *[repr(v) for v in velocities]).splitlines()
        self.assertEqual(output[-1], "sha256 " + hashlib.sha256(TABLE.read_bytes()).hexdigest())
        for line, velocity in zip(output[:-1], velocities):
            v, ohtaus, clamped, k_d, k_ion = line.split()
            expected = reference_rates(velocity)
            with self.subTest(v_h=velocity):
                self.assertEqual(float(v), velocity)
                self.assertAlmostEqual(float(ohtaus), expected[0], delta=1e-12 * expected[0])
                self.assertEqual(int(clamped), expected[1])
                self.assertLessEqual(abs(float(k_d) / expected[2] - 1.0), 1e-12)
                self.assertLessEqual(abs(float(k_ion) / expected[3] - 1.0), 1e-12)

    def test_case_values_of_the_design_note(self) -> None:
        # Section 3 of CHEMISTRY_VARIANTS.md lists 1/k_d per case; they must be
        # what the independent evaluation gives (rounded to 1 s).
        expected_lifetimes = {1: 180002, 2: 167381, 3: 152533, 4: 121100, 5: 149074, 6: 170396}
        for case, velocity in case_velocities().items():
            with self.subTest(case=case):
                self.assertEqual(round(1.0 / reference_rates(velocity)[2]),
                                 expected_lifetimes[case])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    BINARY = Path(sys.argv.pop(1)).resolve()
    unittest.main()
