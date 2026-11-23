"""Independent geometry and tracing controls for the observation-side screen."""

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from reproduce import (
    calibrated_positions, compile_trace, coordinates, digest, ray_to_shell,
    RSUN_M, RSS_RSUN, SHELL_RSUN, trace,
)


class MappingTests(unittest.TestCase):
    def test_corrected_wcs_paired_provenance_and_native_order(self):
        # Fabricated bytes verify the package contract in isolation. A
        # calibration product cannot be mistaken for a different raw raster
        # or silently reorder its time-dependent exposure columns.
        with tempfile.TemporaryDirectory() as scratch:
            folder = Path(scratch)
            raster_file = folder / "test.fits"
            raster_file.write_bytes(b"unit-test raw raster marker")
            radiance = folder / "test.sav"
            radiance.write_bytes(b"unit-test SolarSoft marker")
            table = folder / "test.wcs.csv"
            table.write_text("x_index_zero,y_index_zero,x_arcsec,y_arcsec\n"
                             "1,0,-980,23\n0,0,-981,22\n")
            record = {"source_sha256": digest(raster_file), "wcs_sha256": digest(table),
                      "radiance_sha256": digest(radiance)}
            manifest = {"rasters": {raster_file.name: record}}
            x, y = calibrated_positions(folder, raster_file, (2, 1), manifest)
            np.testing.assert_array_equal(x[:, 0], [-981, -980])
            np.testing.assert_array_equal(y[:, 0], [22, 23])
            record["source_sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "calibration provenance mismatch"):
                calibrated_positions(folder, raster_file, (2, 1), manifest)

    def test_limb_ray_geometry(self):
        # A ray 980 arcsec from disk centre at 1 AU crosses two sides of
        # the 1.05-Rsun shell. Due east must decrease Carrington longitude.
        impact, (near, far) = ray_to_shell(-980.0, 0.0, 1.0, 0.0, 290.0)
        self.assertGreater(impact, 1.0)
        for point in (near, far):
            self.assertAlmostEqual(np.linalg.norm(point) / RSUN_M, SHELL_RSUN, places=13)
            lat, lon = coordinates(*point)
            self.assertAlmostEqual(lat, 0.0)
            self.assertLess(lon, 290.0)
        self.assertLess(np.linalg.norm(near - far) / RSUN_M, 1.0)
        self.assertTrue(near[0] * math.cos(math.radians(290))
                        + near[1] * math.sin(math.radians(290))
                        > far[0] * math.cos(math.radians(290))
                        + far[1] * math.sin(math.radians(290)))
        self.assertEqual(ray_to_shell(0, 0, 1, 0, 290)[1], ())
        self.assertEqual(ray_to_shell(-1200, 0, 1, 0, 290)[1], ())

    def test_dipole_open_closed_and_source_surface_polarity(self):
        # A pure axial dipole provides an open positive and an open negative
        # line as well as a closed line. This prevents an all-closed result
        # triggered by a sign, endpoint, or source-surface integration bug.
        mu = -1.0 + (np.arange(180) + 0.5) / 90.0
        field = np.tile(mu[:, np.newaxis], (1, 360))
        latitudes = (65, -65, 0)
        longitude = math.radians(173)
        seeds = []
        for latitude in latitudes:
            lat = math.radians(latitude)
            seeds.append(RSUN_M * SHELL_RSUN * np.array([
                math.cos(lat) * math.cos(longitude), math.cos(lat) * math.sin(longitude),
                math.sin(lat),
            ]))
        with tempfile.TemporaryDirectory() as scratch:
            executable = Path(scratch) / "trace_pfss"
            compile_trace(executable)
            mean, lines = trace(executable, field, seeds, 1, 0.01)
        self.assertAlmostEqual(mean, 0.0, places=17)
        self.assertEqual([line[0] for line in lines], ["open", "open", "closed"])
        # The analytic PFSS dipole poloidal flux is proportional to
        # (r^2 + 2 Rss^3/r) cos^2(latitude), constant along each line.
        ratio = (SHELL_RSUN ** 2 + 2.0 * RSS_RSUN ** 3 / SHELL_RSUN) / (3.0 * RSS_RSUN ** 2)
        for latitude, (_, x, y, z, radial_field) in zip(latitudes[:2], lines[:2]):
            lat, lon = coordinates(x, y, z)
            target = math.copysign(math.degrees(math.acos(
                math.sqrt(ratio) * math.cos(math.radians(latitude)))), latitude)
            self.assertAlmostEqual(np.linalg.norm([x, y, z]) / RSUN_M, RSS_RSUN, places=12)
            self.assertAlmostEqual(lon, 173.0, places=10)
            self.assertAlmostEqual(lat, target, places=5)
            self.assertGreater(lat * radial_field, 0.0)


if __name__ == "__main__":
    unittest.main()
