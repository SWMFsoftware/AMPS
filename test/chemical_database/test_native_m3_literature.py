"""Compare installed literature tables to independently read primary values."""

import os
from types import SimpleNamespace
import unittest

import test_protocol


class LiteratureRelease(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_LITERATURE_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit r7 literature-catalogue test not selected")
        cls.database = test_protocol.CONTEXT.python_database()

    def test_water_onsets_and_signal_normalization(self):
        # Cross-check the exact water subsection in the original printed
        # source and the independent PMC6624710 native text. The OH+ onset
        # has no uncertainty in angstrom in that sentence.
        onsets = self.database.table("dibeler1966-water:observed-onsets")
        self.assertEqual(onsets.rows, (("H2O+", "984.5", "1", "12.593"),
                                      ("OH+", "687", "", "18.05")))
        signals = self.database.table("dibeler1966-water:relative-ion-signals")
        self.assertEqual(signals.rows, (("H2O+", "584", "21.23", "1.0"),
                                       ("OH+", "584", "21.23", "0.4"),
                                       ("O+", "584", "21.23", "0.05")))

    def test_cui_branch_uncertainties_and_cpp_read(self):
        # These native numbers were verified against MDPI HTML Table 1,
        # independently of the recipe's PDF text layer. Match source roles,
        # the explicit missing observation, and significant trailing zeros.
        identifier = "cui2025-he-capture:table1"
        table = self.database.table(identifier)
        self.assertEqual(table.rows, (
            ("0.5", "100", "", "99.27", "", "", "0.73"),
            ("1.25", "99.13", "0.66", "98.62", "0.87", "0.45", "1.38"),
            ("2.5", "97.75", "0.51", "97.71", "2.25", "0.33", "2.29"),
            ("3.75", "97.32", "0.68", "96.87", "2.68", "0.43", "3.13"),
            ("5", "96.89", "0.42", "96.19", "3.11", "0.27", "3.81"),
        ))
        identity = dict(self.database.identity)
        fixture = SimpleNamespace(root=self.database.root, fields=identity, export_path=identity["export_path"])
        result = test_protocol.CONTEXT.native(fixture, "table", identifier, "table-access")
        self.assertTrue(result["open_ok"], result)
        self.assertTrue(result["query_ok"], result)
        self.assertEqual(result["rows"], [list(row) for row in table.rows])

    def test_kim_neutral_fits_and_held_payload(self):
        fit = self.database.table("kim1992-neutral-ionization:table1")
        self.assertEqual(fit.rows, (
            ("H", "13.61", "0.7576", "-5.521", "5.867", "2.948"),
            ("He", "24.58", "0.7326", "-5.117", "3.295", "2.468"),
        ))
        with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as refused:
            self.database.source_path("m3-literature-20261124:signalhelium-14173434")
        self.assertEqual(refused.exception.code, "PermissionRequired")
