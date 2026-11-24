"""Verify that published lunar OH assumptions retain their species, event and scope."""

import os
import unittest

import test_protocol


class LunarSurfaceModelPrecedents(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if not test_protocol.CONTEXT or os.environ.get("AMPS_CHECK_LUNAR_SURFACE_PRECEDENT_RELEASE") != "1":
            raise unittest.SkipTest("lunar OH surface-precedent release was not selected")
        cls.db = test_protocol.CONTEXT.python_database()
        cls.dataset = "lunar-oh-surface-model-precedents"
        cls.entries = {entry["id"]: entry for entry in cls.db.describe(cls.dataset)["files"]}

    def test_jones_retention_is_conditional_on_return(self):
        record = self.entries[self.dataset + ":jones-k7"]["quantity_metadata"]
        # The authors divide the gas-phase photofragments into escaping and
        # returning groups, then apply k7 only to impacts by the second group.
        self.assertEqual(record["modeled_radical_escape_fraction"], 0.5)
        self.assertEqual(record["impact_sticking_probability"], 1.0)
        self.assertFalse(record["activation_barrier_assumed"])
        self.assertEqual(record["sticking_temperature_dependence"], "constant")
        self.assertIsNone(record["energy_accommodation_coefficient"])
        self.assertIn("chemisorption", record["event"])

    def test_thermal_launch_and_residence_preserve_the_two_models(self):
        smolka = self.entries[self.dataset + ":smolka-launch"]["quantity_metadata"]
        self.assertEqual(smolka["energy_accommodation_coefficient"], 1.0)
        self.assertEqual(smolka["launch_speed_weight"], "v^3 exp(-m v^2/(2 k_B T))")
        self.assertIsNone(smolka["impact_sticking_probability"])
        archived = self.entries[self.dataset + ":laferriere-hops"]["quantity_metadata"]
        self.assertEqual(archived["software_doi"], "10.6084/m9.figshare.28548857.v1")
        self.assertEqual(set(archived["species"]), {"OH", "H2O"})
        self.assertEqual((archived["desorption_barrier_eV"], archived["attempt_frequency_s_inv"],
                          archived["coverage_fraction"]), (0.6, 10**13, 0.01))
        self.assertEqual(archived["residence_time_s"], "1 / (f A exp(-E_d/(k_B T)))")
        self.assertEqual(archived["launch_speed_weight"], "v^2 exp(-m v^2/(2 k_B T))")
        self.assertIsNone(archived["energy_accommodation_coefficient"])

    def test_release_preserves_evidence_and_access_classes(self):
        self.assertEqual(len(self.entries), 3)
        for key, record in self.entries.items():
            with self.subTest(source=key):
                self.assertEqual(record["quantity_metadata"]["evidence_class"], "published model assumption")
                self.assertFalse(record["quantity_metadata"]["measured_lunar_oh_scattering"])
                self.assertFalse(record["physical_activation"])
                self.assertEqual(record["policy"], "metadata-only")
                with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as error:
                    self.db.source_path(key)
                self.assertEqual(error.exception.code, "PermissionRequired")
