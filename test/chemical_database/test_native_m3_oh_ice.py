"""Check the independently read rows in the IUPAC HO + ice data sheet."""

import os
from types import SimpleNamespace
import unittest

import test_protocol


class EvaluatedIceUptake(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_OH_ICE_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("real IUPAC OH/ice source release was not selected")
        cls.database = test_protocol.CONTEXT.python_database()

    def test_native_conditions_and_two_error_conventions(self):
        # Read independently from the experimental and preferred-value
        # columns on page 9075. The broader temperature range belongs to
        # the evaluators and no initial error is supplied for fresh ice.
        table = self.database.table("iupac2010-oh-ice:uptake")
        self.assertEqual(table.columns, ("evidence", "surface_condition", "temperature_min_K",
                                        "temperature_max_K", "uptake_gamma", "absolute_gamma_plus_minus",
                                        "delta_log_gamma_plus_minus"))
        self.assertEqual(table.rows, (
            ("experimental", "fresh", "205", "230", "0.1", "", ""),
            ("experimental", "conditioned", "205", "230", "0.03", "0.02", ""),
            ("evaluated", "steady_state_clean_ice", "200", "230", "0.03", "", "0.5"),
        ))
        identity = dict(self.database.identity)
        fixture = SimpleNamespace(root=self.database.root, fields=identity, export_path=identity["export_path"])
        native = test_protocol.CONTEXT.native(fixture, "table", "iupac2010-oh-ice:uptake", "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in table.rows])

    def test_evaluation_and_primary_experiment_are_attributed_once(self):
        source, derived = sorted(self.database.describe("iupac2010-oh-ice")["files"],
                                 key=lambda record: record["role"], reverse=True)
        self.assertEqual(source["id"], "iupac2010-oh-ice:evaluation")
        self.assertEqual(list(derived["derived_from"]), [source["id"]])
        for record in (source, derived):
            meta = record["quantity_metadata"]
            self.assertEqual(meta["doi"], "10.5194/acp-10-9059-2010")
            self.assertEqual(meta["primary_experiment_doi"], "10.1021/jp952142z")
            self.assertIn("1 Torr", meta["method"])
            self.assertIn("Irreversible gas-phase loss", meta["limitation"])
            self.assertFalse(record["physical_activation"])
