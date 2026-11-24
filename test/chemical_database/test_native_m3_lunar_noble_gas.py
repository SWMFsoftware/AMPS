"""Verify the three published neutral rates, native metadata and public bindings."""
from decimal import Decimal
import hashlib
import os
from types import SimpleNamespace
import unittest

import test_protocol


class LunarNobleGas(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if not test_protocol.CONTEXT or os.environ.get("AMPS_CHECK_M3_LUNAR_NOBLE_GAS") != "1":
            raise unittest.SkipTest("lunar noble-gas release was not selected")
        cls.db = test_protocol.CONTEXT.python_database()
        cls.table = "peschel2026:quiet-sun-neutral-loss"
        cls.entries = {row["id"]: row for row in cls.db.describe("peschel2026-lunar-noble-gas-loss")["files"]}

    def test_printed_rates_in_both_interfaces(self):
        # Read independently from the typeset page. The powers of ten below
        # use ordinary scientific notation; the CSV also retains the printed
        # mantissa of 0.525 for helium for a round trip to the paper.
        rows = self.db.table(self.table).rows
        self.assertEqual([(r[0], Decimal(r[3])) for r in rows],
                         [("He", Decimal("5.25e-8")), ("Ne", Decimal("5.9e-7")),
                          ("40Ar", Decimal("3.06e-7"))])
        self.assertEqual([r[1] for r in rows], ["0.525", "5.9", "3.06"])
        identity = dict(self.db.identity)
        fixture = SimpleNamespace(root=self.db.root, fields=identity, export_path=identity["export_path"])
        native = test_protocol.CONTEXT.native(fixture, "table", self.table, "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in rows])

    def test_publisher_pdf_and_scientific_scope(self):
        pdf = self.db.source_path("peschel2026:article")
        self.assertEqual(hashlib.sha256(pdf.read_bytes()).hexdigest(),
                         "2aa3fc35bf93027e5879c904250c840e0a3db2ee37fd4330f5bdd5feffd978ae")
        source = self.entries["peschel2026:article"]
        self.assertEqual(source["policy"], "aas-peschel-2026-cc-by-4.0")
        metadata = source["quantity_metadata"]
        self.assertIn("quiet-Sun", metadata["spectrum_label"])
        self.assertEqual(metadata["spatial_domain"], "illuminated exosphere")
        self.assertIn("1.12e7", metadata["source_flux_note"])
        self.assertIn("1.12e-7", metadata["source_flux_note"])
        self.assertIn("neon", metadata["limitation"])
        self.assertIsNone(metadata["measurement_uncertainty"])
        self.assertIsNone(metadata["daughter_state_probabilities"])
        self.assertTrue(all(row["physical_activation"] is False for row in self.entries.values()))
        self.assertEqual(tuple(self.entries[self.table + ":csv"]["derived_from"]), ("peschel2026:article",))
