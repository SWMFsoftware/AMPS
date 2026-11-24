"""Check the release-member lineage and the limits of the sodium physics evidence."""

from decimal import Decimal
import hashlib
import json
import os
import unittest

import test_protocol


class SodiumGroundCaptureBenchmark(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_SODIUM_SOURCE_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("sodium source release was not explicitly selected")
        cls.db = test_protocol.CONTEXT.python_database()
        files = cls.db.describe("sodium-ground-capture-benchmark-20261124")["files"]
        cls.entries = {entry["id"]: entry for entry in files}

    def test_native_citation_and_knots_match_earlier_independent_decode(self):
        path = self.db.source_path("iaea-d119990:release-2024.1")
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         "418a70144dce1719c38ed4d72d5994ecf370c9bb8f5778722c689475df618537")
        lines = path.read_text().splitlines()
        header = json.loads(lines[0])
        self.assertEqual(header["refs"], {"B176": {"doi": "10.1103/physreva.70.050702"}})
        self.assertNotIn("unc_perc", header)
        # The legacy CSV was independently decoded from the HTML/Bokeh
        # arrays. Compare decimal values so different formatting is harmless.
        release_rows = tuple(tuple(map(Decimal, line.split())) for line in lines[2:])
        self.assertEqual(len(release_rows), 12)
        old_rows = self.db.table("collisiondb_119990_embedded_arrays.csv").rows
        self.assertEqual(release_rows, tuple(tuple(map(Decimal, row)) for row in old_rows))
        # Users often look up the pre-existing table directly. Its current
        # catalogue record must expose the figure crosswalk with the old knots.
        table = self.db.describe("collisiondb_119990_embedded_arrays.csv")
        self.assertEqual(table["source_figure_review_dataset"], "sodium-ground-capture-benchmark-20261124")
        self.assertIn("does not resolve the outgoing H(1s) partial", table["notes"])

    def test_source_scopes_and_permissions_remain_bound_to_their_own_sources(self):
        chart = self.entries["nifs-chart:38917-knoop"]["quantity_metadata"]
        paper = self.entries["shingal1986-na:ground-partial-benchmark"]["quantity_metadata"]
        self.assertIn("Figure 2", chart["locator"])
        self.assertIn("sums over hydrogen final states", chart["scope"])
        self.assertEqual(chart["same_archive_knots"], 12)
        self.assertIn("0.22 Å²", paper["assessment"])
        self.assertIn("2.2E-17 cm²", paper["assessment"])
        self.assertIn("No 4 keV/u H(1s) uncertainty reported", paper["uncertainty"])
        followup = self.entries["shingal1987-na:ground-partial-benchmark"]["quantity_metadata"]
        self.assertIn("4.00 keV projectile laboratory", followup["assessment"])
        self.assertIn("4.073E-17 cm²", followup["assessment"])
        self.assertIn("overlap in authorship", followup["uncertainty"])
        # The follow-up tabulates E_lab in keV, whereas A34 tabulates keV/u.
        # A client must retain the original energy basis and null uncertainties.
        for entry, energy, unit, sigma in ((paper, "4.0", "keV/u", "0.22"),
                                           (followup, "4.00", "keV", "0.4073")):
            point = entry["benchmark_point"]
            self.assertEqual((point["projectile_energy"], point["energy_unit"]), (energy, unit))
            self.assertEqual(point["energy_frame"], "target-rest laboratory")
            self.assertEqual((point["cross_section"], point["area_unit"]), (sigma, "Å²"))
            self.assertEqual(Decimal(point["cross_section_cm2"]), Decimal(sigma) * Decimal("1e-16"))
            self.assertIsNone(point["uncertainty"])
        for identifier in ("nifs-chart:38917-knoop", "shingal1986-na:ground-partial-benchmark",
                           "shingal1987-na:ground-partial-benchmark"):
            entry = self.entries[identifier]
            self.assertEqual(entry["policy"], "metadata-only")
            self.assertFalse(entry["physical_activation"])
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as refused:
                self.db.source_path(identifier)
            self.assertEqual(refused.exception.code, "PermissionRequired")
