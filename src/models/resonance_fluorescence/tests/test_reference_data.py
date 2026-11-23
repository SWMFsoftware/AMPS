#!/usr/bin/env python3
"""Verify provenance, normalized grids, and generated core reference data."""

from __future__ import annotations

import builtins
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest


MODULE = Path(__file__).resolve().parents[1]
REFERENCE = MODULE / "reference" / "g_factors"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Two valid float summation algorithms.  The built-in sum() of Python <= 3.11
# adds left to right; Python >= 3.12 uses compensated summation.  Integer
# sums are exact in both and are returned unchanged, so only float rounding
# differs.
def _left_to_right_sum(values, start=0):
    total = start
    for value in values:
        total = total + value
    return total


def _compensated_sum(values, start=0):
    values = list(values)
    if isinstance(start, int) and all(isinstance(value, int) for value in values):
        return builtins.sum(values, start)
    total = float(start)
    compensation = 0.0
    for value in values:
        value = float(value)
        step = total + value
        if abs(total) >= abs(value):
            compensation += (total - step) + value
        else:
            compensation += (value - step) + total
        total = step
    return total + compensation


def _fresh_generator(name: str, summation=None):
    """Load an independent copy of the generator; optionally bind ``sum``."""
    spec = importlib.util.spec_from_file_location(
        name, MODULE / "tools/generate_reference_tables.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if summation is not None:
        module.sum = summation
    return module


class ReferenceDataTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads((REFERENCE / "manifest.json").read_text())

    def test_manifest_hashes_match_every_authoritative_file(self) -> None:
        models = self.manifest["models"]
        h = models["H_LYA_IKL2020"]
        oh = models["OH_AX00_SA1988_LOWELL_UNQUENCHED"]
        self.assertEqual(sha256(REFERENCE / h["source_file"]), h["source_sha256"])
        self.assertEqual(sha256(REFERENCE / oh["raw_file"]), oh["raw_sha256"])
        self.assertEqual(
            sha256(REFERENCE / oh["normalized_file"]), oh["normalized_sha256"]
        )
        for name in (
            "O_I_1304_KVB2022_DIRECT_J2",
            "O_I_1304_SUMER2001_DIRECT_J2",
            "O_I_1304_SUMER2001_LYB_BOWEN_J2",
        ):
            model = models[name]
            self.assertEqual(
                sha256(REFERENCE / model["source_file"]), model["source_sha256"], name
            )
        for path, digest in models["O_I_1304_KVB2022_DIRECT_J2"]["extracted_files"].items():
            self.assertEqual(sha256(REFERENCE / path), digest, path)
        nist = self.manifest["atomic_data"]["NIST_ASD_O_I"]
        for kind in ("ultraviolet", "infrared"):
            self.assertEqual(
                sha256(REFERENCE / nist[f"{kind}_file"]), nist[f"{kind}_sha256"], kind
            )
        generated = self.manifest["generated"]
        self.assertEqual(
            sha256((REFERENCE / generated["file"]).resolve()), generated["sha256"]
        )

    def test_lowell_table_is_complete_and_not_extrapolated(self) -> None:
        oh = self.manifest["models"]["OH_AX00_SA1988_LOWELL_UNQUENCHED"]
        with (REFERENCE / oh["normalized_file"]).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 1201)
        velocity = [float(row["heliocentric_radial_velocity_km_s"]) for row in rows]
        self.assertEqual((velocity[0], velocity[-1]), (-60.0, 60.0))
        self.assertTrue(
            all(abs(right - left - 0.1) < 1.0e-10 for left, right in zip(velocity, velocity[1:]))
        )
        self.assertTrue(all(row["api_extrapolated"] == "false" for row in rows))
        self.assertTrue(
            all(float(row["oh00_energy_g_1au_erg_s-1_molecule-1"]) > 0 for row in rows)
        )
        with (REFERENCE / oh["raw_file"]).open() as stream:
            raw = [json.loads(line) for line in stream]
        self.assertEqual(len(raw), 1201)
        self.assertTrue(all(row["response"]["extrapolated"] is False for row in raw))

    def test_oxygen_tables_are_reviewed_and_consistent(self) -> None:
        sys.path.insert(0, str(MODULE / "tools"))
        try:
            import generate_reference_tables as generator
        finally:
            sys.path.pop(0)
        oxygen = generator.oxygen_tables()
        # 0.352-au Killen et al. values imply 8.41e9 photons cm^-2 s^-1 for
        # the triplet at 1 AU; SOLSTICE gives 8.12e9 on 2009-03-01.
        self.assertAlmostEqual(oxygen["reference_flux_m2"] / 8.12e13, 1.0, delta=0.05)
        self.assertAlmostEqual(
            oxygen["reference_flux_m2"] * 1.0e-4,
            self.manifest["models"]["O_I_1304_KVB2022_DIRECT_J2"][
                "reference_triplet_flux_1au_photon_cm-2_s-1"
            ],
            delta=1.0,
        )
        self.assertEqual(len(oxygen["kvb_rate"]), 961)
        self.assertEqual(len(oxygen["sumer_direct"]), 2001)
        self.assertEqual(len(oxygen["bowen_rate"]), 2001)
        self.assertTrue(all(value > 0.0 for value in oxygen["kvb_rate"]))
        self.assertTrue(all(value > 0.0 for value in oxygen["sumer_direct"]))
        self.assertTrue(all(value > 0.0 for value in oxygen["bowen_rate"]))
        # 3d 3D -> 3p 3P (1128.7 nm) branching from NIST A values.
        for _, _, _, branching in oxygen["bowen"]:
            self.assertAlmostEqual(branching, 0.288, delta=0.002)
        # Line-core agreement of the two independent solar profiles: with the
        # same triplet flux the 1-AU rates agree within 30% near rest.
        reference = oxygen["reference_flux_m2"]
        for velocity in (-5.8, 0.4, 7.1):
            kvb = oxygen["kvb_rate"][round((velocity + 48.0) / 0.1)]
            sumer = oxygen["sumer_direct"][round((velocity + 100.0) / 0.1)] * reference
            self.assertLess(abs(sumer / kvb - 1.0), 0.30, velocity)
        # Mean triplet photon wavelength between 130.2 and 130.6 nm.
        wavelength = generator.PLANCK_J_S * generator.LIGHT_SPEED_M_S / oxygen["mean_energy"]
        self.assertTrue(130.2e-9 < wavelength < 130.6e-9)

    def test_generated_include_is_current(self) -> None:
        subprocess.run(
            [sys.executable, str(MODULE / "tools/generate_reference_tables.py"), "--check"],
            check=True,
        )

    def test_generated_include_is_independent_of_float_summation(self) -> None:
        # --check compares text exactly, so the generated bits must not depend
        # on the interpreter.  A float sum() inside the generator would make
        # the two renders below differ (1-3 ulp in the O I tables, found on
        # 2026-11-23 between Python 3.8 and 3.14), on any interpreter.
        reference = _fresh_generator("_gen_reference").render()
        for name, summation in (
            ("left-to-right", _left_to_right_sum),
            ("compensated", _compensated_sum),
        ):
            rendered = _fresh_generator("_gen_" + name.replace("-", "_"), summation).render()
            self.assertTrue(
                rendered == reference,
                f"generated include changes with {name} float summation; "
                "use math.fsum for every float sum in generate_reference_tables.py",
            )


if __name__ == "__main__":
    unittest.main()
