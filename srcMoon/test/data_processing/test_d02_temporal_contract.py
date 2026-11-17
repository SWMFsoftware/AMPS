#!/usr/bin/env python3
"""D02 metadata/error protocol checks; fixtures are never observations.

Tiny temporary files exercise hashing, label parsing and fail-closed status.
No model, DN science value, reference visibility, or observational uncertainty
is manufactured or scored. The real audit always uses committed native hashes.
These checks are automatically invoked by the existing M2 package-audit CLI.
"""
import hashlib
from pathlib import Path
import tempfile
import unittest

import audit_d02_temporal_contract as audit
from prepare_d02_reference import illumination_operator_metadata


class NativeContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="amps-d02-protocol-")
        self.package = Path(self.temporary.name) / "lola"
        self.raw = self.package / "raw/illumination"
        self.raw.mkdir(parents=True)
        self.name = "AVGVISIB_65N_240M_201608"
        self.label = self.raw / (self.name + ".LBL")
        self.image = self.raw / (self.name + ".IMG")
        self.label.write_text('''/* METADATA PROTOCOL FIXTURE, NEVER SCIENCE EVIDENCE */
PRODUCT_ID = "AVGVISIB_65N_240M_201608"
START_TIME = 2009-07-13T17:33:17.246
STOP_TIME = 2013-07-18T14:12:13
DESCRIPTION = "by any fraction of the solar disc"
LINES = 2
LINE_SAMPLES = 2
SAMPLE_BITS = 16
SAMPLE_TYPE = LSB_INTEGER
SCALING_FACTOR = 0.00004
OFFSET = 0.0
UNIT = NONE
MAP_SCALE = 240 <m/pix>
MAP_PROJECTION_TYPE = "POLAR STEREOGRAPHIC"
POSITIVE_LONGITUDE_DIRECTION = "EAST"
KEYWORD_LATITUDE_TYPE = "PLANETOCENTRIC"
COORDINATE_SYSTEM_NAME = "MEAN EARTH/POLAR AXIS OF DE421"
''')
        # Only its byte length is examined; these bytes supply no DN values.
        self.image.write_bytes(bytes(8))

    def tearDown(self):
        self.temporary.cleanup()

    def hashes(self):
        return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (self.label, self.image)}

    def inspect(self):
        return audit.inspect(self.package, self.hashes())

    def replace(self, before, after):
        text = self.label.read_text()
        self.assertIn(before, text)
        self.label.write_text(text.replace(before, after))

    def test_native_time_fields_cannot_supply_averaging_epochs(self):
        result = self.inspect()
        self.assertEqual(result["status"], "SKIPPED")
        self.assertEqual(result["native_integrity_status"], "PASS")
        self.assertTrue(all(v is None for v in result["temporal_contract"].values()))
        product = result["products"][0]
        self.assertEqual(product["native_start_time"], "2009-07-13T17:33:17.246")
        self.assertFalse(product["label_times_used_as_simulation_bounds"])
        self.assertEqual(product["quantity"], "fraction_of_sampled_times_any_solar_disc_visible")

    def test_changed_bytes_are_error(self):
        frozen = self.hashes()
        self.image.write_bytes(bytes(10))
        result = audit.inspect(self.package, frozen)
        self.assertEqual(result["status"], "ERROR")
        self.assertIn("SHA-256", result["reason"])

    def test_missing_native_file_is_skipped(self):
        frozen = self.hashes()
        self.image.unlink()
        result = audit.inspect(self.package, frozen)
        self.assertEqual(result["status"], "SKIPPED")
        self.assertEqual(result["missing_files"], [str(self.image)])

    def test_inconsistent_byte_dimensions_are_error(self):
        self.replace("LINES = 2", "LINES = 3")
        self.assertEqual(self.inspect()["status"], "ERROR")

    def test_unknown_frame_is_error(self):
        self.replace("MEAN EARTH/POLAR AXIS OF DE421", "unspecified frame")
        self.assertEqual(self.inspect()["status"], "ERROR")

    def test_other_observable_is_not_visibility_occurrence(self):
        self.replace("by any fraction of the solar disc", "mean visible solar disc area")
        self.assertEqual(self.inspect()["status"], "ERROR")

    def test_duplicate_native_scalar_is_error(self):
        with self.label.open("a") as stream:
            stream.write("START_TIME = 2000-01-01T00:00:00\n")
        self.assertEqual(self.inspect()["status"], "ERROR")

    def test_absent_units_are_error(self):
        self.replace("240 <m/pix>", "240")
        self.assertEqual(self.inspect()["status"], "ERROR")

    def test_inspection_preserves_all_input_bytes(self):
        before = self.hashes()
        self.inspect()
        self.assertEqual(self.hashes(), before)

    def test_spatial_qa_cannot_invent_exact_product_sampling(self):
        """A spatial transform must not freeze general publication defaults."""
        operator = illumination_operator_metadata()
        self.assertEqual(operator["status"], "UNRESOLVED")
        for key in ("reference_cycle_years", "reference_time_step_hours",
                    "simulation_start_epoch", "simulation_stop_epoch",
                    "horizon_azimuth_spacing_deg"):
            self.assertIsNone(operator[key], key)
        self.assertFalse(operator["documentation_context"][
            "exact_201608_temporal_contract_established"])


if __name__ == "__main__":
    unittest.main(buffer=True)
