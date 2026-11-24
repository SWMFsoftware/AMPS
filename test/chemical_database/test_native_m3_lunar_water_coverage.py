"""Check which Schörghofer H2O model parameters and mechanisms survive publication."""
import os
import unittest

import test_protocol


class LunarWaterCoverageRelease(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if not test_protocol.CONTEXT or os.environ.get("AMPS_CHECK_LUNAR_WATER_COVERAGE_RELEASE") != "1":
            raise unittest.SkipTest("Schörghofer lunar-water source checks were not selected")
        db = test_protocol.CONTEXT.python_database()
        dataset = "schorghofer-lunar-water-models-v123"
        cls.records = {entry["logical_name"]: entry for entry in db.describe(dataset)["files"]}
        cls.database = db

    def test_activated_water_sites_and_coverage_branches(self):
        # Compare with the actual tagged module constants and the final S=S1
        # assignment. The other multilayer prescription S2 is evaluated but
        # does not supply the return value in v1.2.3.
        record = self.records["desorptionrates.f90"]
        self.assertEqual(record["sha256"], "0c15a23cc629a096913085e7df12f5039bb9e90003696a0efbb0f54b29e33e9f")
        m = record["quantity_metadata"]
        self.assertEqual((m["monolayer_molecules_m2"], m["Ec_eV"], m["W_eV"], m["Eice_eV"]),
                         (1e19, 0.65, 0.22, 0.529))
        self.assertEqual(m["attempt_frequency_prefactor_s_inv"], 2.2e17)
        self.assertEqual(m["submonolayer_domain"], "0 < v <= 1")
        self.assertEqual(m["submonolayer_formula"], "nu theta exp(-Ec b) v^(W b) / (1+W b)")
        self.assertEqual(m["multilayer_gamma"], "1/v^2")
        self.assertEqual(m["zero_coverage_flux"], 0)
        self.assertIsNone(m["incident_sticking_probability"])
        self.assertIsNone(m["kinetic_energy_accommodation"])

    def test_illuminated_water_loss_and_release_scope(self):
        # The model follows an illuminated thermal hop and removes the neutral.
        # The inactive source comment listing an OH frequency cannot supply a
        # downstream radical law. The archive should make this scope visible.
        loss = self.records["montecarlo.f90"]["quantity_metadata"]
        self.assertEqual(loss["loss_frequency_s_inv_1au"], 1.26e-5)
        self.assertEqual(loss["illumination"], "in-flight loss enabled for Q > 0")
        self.assertIn("half-normal", loss["launch_law"])
        self.assertEqual(loss["loss_probability"], "for integrated hazard h <= 0.2 use h, otherwise use 1-exp(-h)")
        self.assertIsNone(loss["daughter_velocities"])
        self.assertIsNone(loss["daughter_states"])
        self.assertEqual(set(self.records), {"desorptionrates.f90", "montecarlo.f90"})
        for record in self.records.values():
            m = record["quantity_metadata"]
            self.assertEqual(m["species"], "H2O")
            self.assertFalse(m["qualified_for_incident_OH"])
            self.assertFalse(record["physical_activation"])
            self.assertEqual(m["software_doi"], "10.5281/zenodo.15161597")
            self.assertEqual(record["policy"], "metadata-only")
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as error:
                self.database.source_path(record["id"])
            self.assertEqual(error.exception.code, "PermissionRequired")
