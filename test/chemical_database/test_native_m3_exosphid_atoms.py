"""Check the scope and native values of the EXOSPHID model's three atomic records."""

import os
import unittest

import test_protocol


class ExosphidAtomicPhotoionization(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_EXOSPHID_ATOMS_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("EXOSPHID atomic source checks require a selected real release")
        cls.db = test_protocol.CONTEXT.python_database()
        cls.entries = {r["id"]: r for r in cls.db.describe("exosphid-noble-gas-photoionization")["files"]}

    def test_solar_frequencies_and_single_channel(self):
        # Independent readback of the producer's qr, ar, tsh, wv, and rp
        # literals. Require the original species and precision to survive.
        rows = {"He": ("5.57e-8", "16.80e-8", "24.59", "504.27"),
                "Ne": ("1.72e-7", "5.74e-7", "21.57", "574.94"),
                "Ar": ("3.11e-7", "7.15e-7", "17.56", "706.00")}
        self.assertEqual(set(self.entries), {"exosphid-atoms:" + atom for atom in rows})
        for atom, expected in rows.items():
            with self.subTest(atom=atom):
                facts = self.entries["exosphid-atoms:" + atom]["quantity_metadata"]
                self.assertEqual(tuple(facts[k] for k in ("quiet_rate_s_inv", "active_rate_s_inv",
                                                          "threshold_eV", "wavelength_cutoff_angstrom")), expected)
                self.assertEqual(facts["wavelength_range_angstrom"], ("0", expected[3]))
                self.assertEqual(facts["species_tokens"], (atom, atom + "(+)", "e(-)"))
                self.assertEqual((facts["reaction_type"], facts["reaction_probability"]), ("SPI", "1.00"))

    def test_argon_ground_threshold_and_admission(self):
        # Prevent a future importer from silently replacing the author's
        # threshold with NIST's value or interpreting an unspecified charge
        # state as a distribution over electronic levels.
        ar = self.entries["exosphid-atoms:Ar"]["quantity_metadata"]
        self.assertEqual(ar["threshold_eV"], "17.56")
        self.assertEqual(ar["ground_threshold_comparison"]["ground_ionization_energy_eV"], "15.7596119")
        self.assertIn("physics.nist.gov/", ar["ground_threshold_comparison"]["source_uri"])
        for key, entry in self.entries.items():
            self.assertIsNone(entry["quantity_metadata"]["residual_ion_state_probabilities"])
            self.assertIsNone(entry["quantity_metadata"]["measurement_uncertainty"])
            self.assertFalse(entry["physical_activation"])
            self.assertEqual(entry["policy"], "metadata-only")
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as denied:
                self.db.source_path(key)
            self.assertEqual(denied.exception.code, "PermissionRequired")
