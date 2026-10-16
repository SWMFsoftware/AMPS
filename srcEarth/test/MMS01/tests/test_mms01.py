"""Data-independent MMS01 policy regressions.

The external 180 MB CDF is exercised by the registered QUICK workflow.  These small
tests protect profile selection and source/gate constants even on build hosts where
the observational bundle is intentionally absent.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import mms01_common as common  # noqa: E402
from prepare_MMS01 import _angular_weights  # noqa: E402


class MMS01PolicyTests(unittest.TestCase):
    def test_quick_is_predeclared_first_middle_last(self) -> None:
        self.assertEqual(common.PROFILES["QUICK"], (0, 19, 38))
        self.assertEqual(common.PROFILES["COMPLETE"], tuple(range(39)))
        self.assertTrue(set(common.PROFILES["QUICK"]) <= set(common.PROFILES["COMPLETE"]))

    def test_science_dimensions_close_exactly(self) -> None:
        self.assertEqual(common.EXPECTED_SWEEPS,
                         common.EXPECTED_DISTRIBUTIONS * common.SWEEPS_PER_DISTRIBUTION)
        self.assertEqual((common.ANODES, common.ENERGIES), (16, 63))

    def test_documented_constants_are_frozen(self) -> None:
        self.assertEqual(common.HPCA_HPLUS_MASS_EV_S2_CM2, 1.0453453e-12)
        self.assertEqual(len(common.ANODE_CORRECTION), common.ANODES)
        self.assertEqual(common.MAX_RELATIVE_ERROR, 1.0e-6)

    def test_documented_hpca_angular_response_is_frozen(self) -> None:
        # CMAD Table 5-6 contains duplicated elevation centroids because the two
        # top-hat anode groups observe azimuths separated by 180 degrees.  Retain
        # that geometry; sorting or replacing it by equally spaced bins changes the
        # calibrated full-sky integral even if all response factors remain positive.
        elevation = np.asarray([
            123.75, 101.25, 78.75, 56.25, 33.75, 11.25, 11.25, 33.75,
            56.25, 78.75, 101.25, 123.75, 146.25, 168.75, 168.75, 146.25,
        ])
        weights = _angular_weights(elevation)
        self.assertAlmostEqual(float(weights[0]), 0.04902707507499252, places=15)
        self.assertAlmostEqual(float(weights[5]), 0.0026922795684574045, places=15)
        # Sixteen spacecraft azimuth samples make one full-sky distribution.
        self.assertAlmostEqual(float(weights.sum() * 16.0), 9.05578279133411, places=14)

    def test_source_identities_are_sha256(self) -> None:
        for byte_count, digest in common.SOURCE_FILES.values():
            self.assertGreater(byte_count, 0)
            self.assertEqual(len(digest), 64)
            int(digest, 16)

    def test_cdf_reader_bootstrap_is_version_and_hash_locked(self) -> None:
        text = common.CDF_READER_REQUIREMENTS.read_text(encoding="ascii")
        self.assertIn(f"cdflib=={common.CDFLIB_VERSION}", text)
        self.assertIn(
            "sha256:49af97acc328c586ac5b7c27fd8e67bccf24af82c5bd8c37d8cfe048a0c1752a",
            text,
        )


if __name__ == "__main__":
    unittest.main()
