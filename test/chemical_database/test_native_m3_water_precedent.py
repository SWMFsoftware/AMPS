"""Check that published water model inputs retain their scope in the shared API."""

import os
import unittest

import test_protocol


class WaterModelPrecedents(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_WATER_PRECEDENT_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("water model precedent packet was not selected")
        database = test_protocol.CONTEXT.python_database()
        records = database.describe("water-photofragment-model-precedents-20261124")["files"]
        cls.records = {item["id"]: item for item in records}

    def test_pinned_lunar_branch_table_and_electron_output_scope(self):
        water = self.records["water-model:H2O.jl"]["quantity_metadata"]
        self.assertEqual(water["commit"], "62e670f9d0344ced1c1414a86b40fbba43b373ef")
        self.assertEqual(water["photon_wavelength_bin_angstrom"], ("0.0", "662.0"))
        # Check the association with each ion. Swapping O+ and H+ leaves
        # the sum unchanged but would materially alter fragment production.
        self.assertEqual([(b["products"], b["threshold_eV"], b["fraction"]) for b in water["branches"]],
                         [("H2O(+) + e-", "12.60", "0.63"), ("H + OH(+) + e-", "18.11", "0.31"),
                          ("H2 + O(+) + e-", "18.65", "0.043"), ("H(+) + OH + e-", "18.72", "0.017")])
        self.assertIsNone(water["uncertainty"])
        dispatch = self.records["water-model:photodestruction_logic.jl"]["quantity_metadata"]
        self.assertEqual(dispatch["returned_heavy_products"], 2)
        self.assertFalse(dispatch["returned_electron"])
        self.assertFalse(dispatch["single_photon_electron_heavy_budget_demonstrated"])

    def test_comet_comparison_preserves_unreported_ion_speeds(self):
        comet = self.records["water-model:shinnaka2017.pdf"]["quantity_metadata"]
        neutral = comet["neutral_channel"]
        self.assertEqual((neutral["H_km_s"], neutral["OH_km_s"], neutral["fraction"]),
                         ("17.5", "1.05", "0.670"))
        self.assertEqual(comet["ionization"]["fraction"], "0.059")
        self.assertIsNone(comet["ionization"]["product_speeds_km_s"])
        self.assertFalse(comet["ionization"]["state_resolved"])
        self.assertEqual(len(self.records), 6)
        for record in self.records.values():
            self.assertEqual(record["policy"], "metadata-only")
            self.assertFalse(record["physical_activation"])
