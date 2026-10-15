import csv
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CSPICE = Path("/home/vtenishe/SPICE/cspice")
KERNELS = Path("/home/vtenishe/SPICE/Kernels")


@unittest.skipUnless((CSPICE / "lib/cspice.a").is_file(),
                     "local CSPICE installation is unavailable")
class SpiceGroundTrackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="spice-ground-track-test-")
        cls.binary = Path(cls.temporary.name) / "spice_ground_tracks"
        subprocess.run([
            "g++", "-std=c++11", "-O2", "-Wall", "-Wextra", "-pedantic",
            "-I" + str(CSPICE / "include"),
            str(ROOT / "tools/spice_ground_tracks.cpp"),
            str(CSPICE / "lib/cspice.a"),
            str(CSPICE / "lib/csupport.a"), "-lm", "-o", str(cls.binary)
        ], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def calculate(self, flyby, utc, pck, spk):
        output = Path(self.temporary.name) / (flyby.lower() + ".csv")
        subprocess.run([
            str(self.binary), str(KERNELS / "NAIF/naif0012.tls"),
            str(KERNELS / "CASSINI" / pck),
            str(KERNELS / "CASSINI" / spk), flyby, utc,
            "-5", "5", "0.1", str(output)
        ], check=True)
        with open(output, newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 101)
        return rows

    def assert_geometry(self, rows):
        for row in rows:
            spacecraft = np.array([float(row[name]) for name in
                                   ("spacecraft_x_m", "spacecraft_y_m",
                                    "spacecraft_z_m")])
            subpoint = np.array([float(row[name]) for name in
                                ("subpoint_x_m", "subpoint_y_m",
                                 "subpoint_z_m")])
            self.assertAlmostEqual(np.linalg.norm(spacecraft - subpoint),
                                   float(row["altitude_m"]), delta=0.02)
            self.assertLess(np.linalg.norm(subpoint), 260000.0)
            self.assertGreaterEqual(float(row["west_lon_planetocentric_deg"]), 0)
            self.assertLess(float(row["west_lon_planetocentric_deg"]), 360)

    def test_e3_geometry_and_closest_sample(self):
        rows = self.calculate("E3", "2008-03-12T19:06:12",
                              "cpck10Mar2008.tpc",
                              "080228AP_SCPSE_08050_08076.bsp")
        self.assert_geometry(rows)
        closest = min(rows, key=lambda row: float(row["altitude_m"]))
        self.assertAlmostEqual(float(closest["time_from_reference_s"]), -0.5,
                               delta=0.11)
        self.assertAlmostEqual(float(closest["altitude_m"]), 52431.2, delta=10.0)
        at_reference = rows[50]
        self.assertAlmostEqual(float(at_reference["west_lon_planetocentric_deg"]),
                               135.090528456, places=6)
        self.assertAlmostEqual(float(at_reference["lat_planetocentric_deg"]),
                               -20.481302187, places=6)

    def test_e5_geometry_and_closest_sample(self):
        rows = self.calculate("E5", "2008-10-09T19:06:43",
                              "cpck16Sep2008.tpc",
                              "081010AP_SCPSE_08282_08315.bsp")
        self.assert_geometry(rows)
        closest = min(rows, key=lambda row: float(row["altitude_m"]))
        self.assertAlmostEqual(float(closest["time_from_reference_s"]), -3.3,
                               delta=0.11)
        self.assertAlmostEqual(float(closest["altitude_m"]), 24959.7, delta=10.0)
        at_reference = rows[50]
        self.assertAlmostEqual(float(at_reference["west_lon_planetocentric_deg"]),
                               97.420177269, places=6)
        self.assertAlmostEqual(float(at_reference["lat_planetocentric_deg"]),
                               -39.402571381, places=6)


if __name__ == "__main__":
    unittest.main()
