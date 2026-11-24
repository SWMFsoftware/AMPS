"""Compare the r10 shared release to the printed primary table cells."""

import os
from types import SimpleNamespace
import unittest

import test_protocol


class LaboratoryBenchmarks(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_LAB_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit r10 laboratory source test not selected")
        cls.database = test_protocol.CONTEXT.python_database()

    def test_wang_condition_table_retains_units_uncertainty_and_native_gaps(self):
        # Independent transcription of the printed Table 2, page 10.
        # The missing humid entries must stay blank through both APIs.
        identifier = "wang2023-quartz-oh:table2"
        table = self.database.table(identifier)
        self.assertEqual(table.rows, (
            ("1.0", "dry helium", "11", "3", "10", "3", "0.011", "0.004"),
            ("1.0", "humid helium", "", "", "93", "26", "", ""),
            ("3.5", "dry helium", "7", "4", "7", "4", "0.010", "0.008"),
            ("3.5", "humid helium", "79", "45", "55", "30", "0.014", "0.008"),
            ("5.2", "dry helium", "17", "6", "15", "4", "0.011", "0.004"),
            ("5.2", "humid helium", "130", "27", "197", "28", "0.006", "0.001"),
        ))
        self.assertIn("diffusive_flux_in_1e14_cm_minus2_s_minus1_native", table.columns)
        self.assertIn("thermal_flux_in_1e16_cm_minus2_s_minus1_native", table.columns)
        identity = dict(self.database.identity)
        fixture = SimpleNamespace(root=self.database.root, fields=identity, export_path=identity["export_path"])
        native = test_protocol.CONTEXT.native(fixture, "table", identifier, "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in table.rows])

    def test_wang_resolution_rows_remain_a_paired_diagnostic(self):
        table = self.database.table("wang2023-quartz-oh:table3")
        self.assertEqual(table.rows, (
            ("10", "5", "3", "7", "4", "7", "4", "0.010", "0.008"),
            ("50", "6", "2", "18", "9", "8", "3", "0.023", "0.014"),
        ))
        files = self.database.describe("m3-lab-benchmarks-20261124")["files"]
        row, = (row for row in files if row["id"] == "wang2023-quartz-oh:table3:csv")
        self.assertIn("one experiment", row["quantity_metadata"]["conditions"])

    def test_sutton_kay_experiment_and_sutton_report_remain_distinct(self):
        preprint = self.database.table("sutton-kay1972-he:table3")
        self.assertEqual(preprint.rows, (
            ("200", "4p", "4s", "16.00", "+10.0", "-14.0"),
            ("200", "4d", "4s", "1.10", "+00.4", "-00.6"),
            ("200", "4f", "4s", "0.74", "+00.2", "-00.2"),
        ))
        # The previously published single-author report gives a separate
        # relative-error estimate from the same experimental programme.
        report = self.database.table("sutton1972-he:excitation-ratios")
        self.assertEqual(report.rows[0][3:5], ("19", "30"))

    def test_retrieved_metadata_cannot_grant_numerical_payload_access(self):
        files = self.database.describe("m3-lab-benchmarks-20261124")["files"]
        held = {row["id"] for row in files if row["policy"] == "metadata-only"}
        self.assertEqual(held, {"water1991-pepico:ntrs", "bellm2007-he-ion-states:crossref",
                                "dogan1998-he-ion-line:crossref", "durr2006-he-recoil:crossref"})
        for identifier in held:
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as refused:
                self.database.source_path(identifier)
            self.assertEqual(refused.exception.code, "PermissionRequired")
