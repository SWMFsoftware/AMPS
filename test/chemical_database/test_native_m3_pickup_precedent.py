"""Check licensed model-frequency provenance and the independently printed medians."""

from decimal import Decimal
import hashlib
import os
from types import SimpleNamespace
import unittest

import test_protocol


class LunarPickupPrecedent(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_PICKUP_PRECEDENT") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("optional Poppe 2022 release checks were not selected")
        cls.db = test_protocol.CONTEXT.python_database()
        files = cls.db.describe("poppe2022-lunar-pickup-rates")["files"]
        cls.files = {item["id"]: item for item in files}

    def test_printed_table_1_medians_and_isotope_rows(self):
        # Independently read from journal Table 1 (p. 4). Its three significant
        # digits differ from the fuller precision of the deposited ASCII file.
        # Round only the values being compared, never the exported CSV.
        expected = (
            ("He", "1.12e-7", "4.48e-9", "5.00e-11"),
            ("Ne", "3.86e-7", "5.84e-10", "1.50e-8"),
            ("Na", "7.59e-6", "8.40e-7", "5.00e-7"),
            ("Ar36", "5.58e-7", "2.37e-7", "1.00e-7"),
            ("Ar", "5.58e-7", "2.37e-7", "1.00e-7"),
        )
        table = self.db.table("poppe2022:atomic-1au-rates")
        self.assertEqual(len(table.rows), len(expected))
        for row, (name, *printed) in zip(table.rows, expected):
            with self.subTest(name=name):
                self.assertEqual(row[1], name)
                for archived, reference in zip(row[3::3], printed):
                    self.assertEqual(Decimal(format(Decimal(archived), ".3g")), Decimal(reference))

    def test_projection_reproduces_original_source_tokens_in_both_readers(self):
        path = self.db.source_path("poppe2022:ionization-rates")
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                         "a208a94f318bbe2c63c59e6e28cfd15bf78901b7f7560480c8298bb235b7826e")
        native_lines = path.read_text(encoding="ascii").splitlines()
        table = self.db.table("poppe2022:atomic-1au-rates")
        self.assertEqual(tuple(int(row[0]) for row in table.rows), (9, 15, 16, 23, 26))
        for line_number, *fields in table.rows:
            self.assertEqual(fields, native_lines[int(line_number) - 1].split())
        identity = dict(self.db.identity)
        fixture = SimpleNamespace(root=self.db.root, fields=identity, export_path=identity["export_path"])
        response = test_protocol.CONTEXT.native(fixture, "table", "poppe2022:atomic-1au-rates", "table-access")
        self.assertTrue(response["query_ok"], response)
        self.assertEqual(response["rows"], [list(row) for row in table.rows])

    def test_rights_model_status_and_uncertainty(self):
        source = self.files["poppe2022:ionization-rates"]
        self.assertEqual(source["policy"], "zenodo-cc-by-4.0")
        self.assertEqual(source["quantity_metadata"]["data_doi"], "10.5281/zenodo.6994648")
        self.assertIsNone(source["quantity_metadata"]["measurement_uncertainty"])
        self.assertIsNone(source["quantity_metadata"]["daughter_state_probabilities"])
        self.assertEqual(self.files["poppe2022:atomic-1au-rates:csv"]["derived_from"],
                         ("poppe2022:ionization-rates",))
        self.assertTrue(all(row["physical_activation"] is False for row in self.files.values()))
        self.assertEqual(self.files["poppe2022:article"]["policy"], "metadata-only")
        with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as denied:
            self.db.source_path("poppe2022:article")
        self.assertEqual(denied.exception.code, "PermissionRequired")
