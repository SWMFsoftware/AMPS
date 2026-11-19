"""Native-source preservation tests; analytical TSV fixtures are not observations."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "utility/chemical_database/new_surface_apollo15221.py"
SPEC = importlib.util.spec_from_file_location("surface_apollo15221", PATH)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class Protocol(unittest.TestCase):
    """Exercise corruption, missingness and exact token/column preservation."""

    def test_signed_measurement_and_model_columns_stay_separate(self):
        raw = ("\t".join(M.SCHEMAS["water-tpd"])+"\n300\t-1E-15\t1.2E-8\t+1.20E-8\t-2.000E-5\n").encode()
        self.assertEqual(M.decode(raw, "water-tpd")[1],
                         ["300", "-1E-15", "1.2E-8", "+1.20E-8", "-2.000E-5"])

    def test_repeated_coverage_is_not_deduplicated(self):
        rows = M.decode(b"Fractional Coverage\tEnergy (J/mol)\n0.99999\t110400\n0.99999\t110500\n", "coverage-energy")
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[1][0], rows[2][0])

    def test_explicit_blank_cell_survives(self):
        rows = M.decode(b"Energy (J/mol)\tProbability\n110300\t\n", "site-distribution")
        self.assertEqual(rows[1], ["110300", ""])

    def test_changed_units_are_not_silently_converted(self):
        with self.assertRaisesRegex(ValueError, "header changed"):
            M.decode(b"Energy (kJ/mol)\tProbability\n110.3\t1E-14\n", "site-distribution")

    def test_ragged_source_is_error(self):
        with self.assertRaisesRegex(ValueError, "ragged"):
            M.decode(b"Energy (J/mol)\tProbability\n110300\n", "site-distribution")

    def test_blank_row_is_not_silently_omitted(self):
        with self.assertRaisesRegex(ValueError, "ragged"):
            M.decode(b"Energy (J/mol)\tProbability\n\n110300\t1E-14\n", "site-distribution")

    def test_missing_data_is_error(self):
        with self.assertRaisesRegex(ValueError, "no data"):
            M.decode(b"Energy (J/mol)\tProbability\n", "site-distribution")

    def test_hash_refusal_precedes_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp)/"native.txt", Path(tmp)/"view.csv"
            source.write_bytes(b"Energy (J/mol)\tProbability\n110300\t1E-14\n")
            with self.assertRaisesRegex(ValueError, "SHA256"):
                M.project(source, output, kind="site-distribution", expected_sha256="0"*64)
            self.assertFalse(output.exists())

    def test_exact_round_trip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp)/"native.txt", Path(tmp)/"view.csv"
            raw=b"Energy (J/mol)\tProbability\r\n110300\t6.25814E-14\r\n"
            source.write_bytes(raw); digest=hashlib.sha256(raw).hexdigest()
            result=M.project(source, output, kind="site-distribution", expected_sha256=digest)
            self.assertEqual(result["row_count"], 1)
            self.assertEqual(source.read_bytes(), raw)
            with output.open(newline="") as stream:
                self.assertEqual(list(csv.reader(stream)), M.decode(raw, "site-distribution"))
            with self.assertRaises(FileExistsError):
                M.project(source, output, kind="site-distribution", expected_sha256=digest)


if __name__ == "__main__":
    unittest.main()
