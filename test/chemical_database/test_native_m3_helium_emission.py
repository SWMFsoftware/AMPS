"""Check the two source-specific helium emission anchors in a selected release."""

import os
import unittest

import test_protocol


class HeliumIonEmission(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_HELIUM_EMISSION") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit helium ion emission source checks not selected")
        cls.database = test_protocol.CONTEXT.python_database()
        dataset = cls.database.describe("he-electron-ion-emission-20261124")
        cls.sources = {row["id"]: row for row in dataset["files"]}

    def test_reported_measurements_match_independent_published_abstracts(self):
        # Read independently from the publisher abstracts and checked against
        # the printed results. Each paper names its own wavelength unit and
        # error. No common error or significance level is inferred.
        for author, wavelength, unit, value, error in (
            ("bloemen", "30.4", "nm", "5.6", "1.4"),
            ("forand", "304", "angstrom", "6.2", "1.7"),
        ):
            with self.subTest(author=author):
                q = self.sources["he-ion-emission:" + author]["quantity_metadata"]
                self.assertEqual((q["incident_energy_eV"], q["wavelength"], q["wavelength_unit"],
                                  q["cross_section"], q["reported_absolute_error"], q["cross_section_unit"]),
                                 ("200", wavelength, unit, value, error, "1e-19 cm2"))
                self.assertIsNone(q["error_confidence_level"])
                self.assertIsNone(q["direct_population_fraction"])

    def test_forand_curve_scope_and_theoretical_role(self):
        q = self.sources["he-ion-emission:forand"]["quantity_metadata"]
        # The Python reader exposes source arrays as immutable tuples.
        self.assertEqual(q["published_curve"], dict(incident_energy_range_eV=("80", "1000"),
                                                    number_of_points=20, normalization_energy_eV="200"))
        self.assertEqual((q["cascade_correction"], q["polarization_correction"]),
                         ("unapplied", "unapplied"))
        theory = self.sources["he-ion-emission:zatsarinny"]["quantity_metadata"]
        self.assertIn("2p and 2s", theory["observable"])
        self.assertNotIn("cross_section", theory)

    def test_original_article_access_is_held_for_all_three_records(self):
        self.assertEqual(set(self.sources), {"he-ion-emission:" + name for name in
                                             ("bloemen", "forand", "zatsarinny")})
        for identifier, row in self.sources.items():
            with self.subTest(identifier=identifier):
                self.assertEqual(row["policy"], "metadata-only")
                self.assertFalse(row["physical_activation"])
                with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as denied:
                    self.database.source_path(identifier)
                self.assertEqual(denied.exception.code, "PermissionRequired")
