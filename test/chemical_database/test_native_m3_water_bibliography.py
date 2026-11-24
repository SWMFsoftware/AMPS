"""Protect the source-specific water ion facts read from ASMS 1975 page 114."""

import os
import unittest

import test_protocol


class WaterPhotofragmentEvidence(unittest.TestCase):
    requires_authoritative_source = True

    def test_fractions_match_printed_page_and_original_remains_held(self):
        if os.environ.get("AMPS_CHECK_WATER_YIELDS_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("real 21.2 eV water photofragment evidence was not selected")
        db = test_protocol.CONTEXT.python_database()
        files = {file["id"]: file for file in db.describe("kronebusch-water-photofragmentation")["files"]}
        entry = files["kronebusch-water:asms1975"]
        evidence = entry["quantity_metadata"]
        # The scan has three zeros after the decimal point in the O+ fraction.
        # Confirmed on the printed page independently of Poppler's OCR text.
        self.assertEqual(dict(evidence["assessed_fraction_by_detected_ion"]),
                         {"OH+": "0.29", "H+": "0.05", "O+": "0.0009"})
        self.assertEqual((evidence["parent"], evidence["wavelength_angstrom"], evidence["photon_energy_eV"]),
                         ("isolated H2O", "584", "21.2"))
        self.assertEqual(evidence["denominator"], "total ionization")
        self.assertIsNone(evidence["fraction_uncertainty"])
        self.assertEqual(evidence["collection_discriminates_above_ion_kinetic_energy_eV"], "1")
        self.assertEqual(entry["sha256"], "40d1c49380c06395bf22ce0a228d8ec1efe21fb110955db4f9ce5a1687352678")
        self.assertIn("kronebusch-water:crossref", files)
        for record in files.values():
            self.assertFalse(record["physical_activation"])
            self.assertEqual(record["policy"], "metadata-only")
            self.assertNotIn("table", record)
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as held:
                db.source_path(record["id"])
            self.assertEqual(held.exception.code, "PermissionRequired")
