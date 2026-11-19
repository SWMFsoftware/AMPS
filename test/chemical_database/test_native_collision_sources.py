"""Administrative parser fixtures only: no synthetic cross-section observations."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("native_collision_sources", ROOT / "utility/chemical_database/native_collision_sources.py")
N = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(N)


def page(rows='<tr><td>pcf(1)</td><td>+2.500e+01</td></tr><tr><td>A2</td><td>-1.0000e-12</td></tr>'):
    # Clearly administrative metadata fixture, not a physical reaction or data.
    return '<h2>DataSet D1</h2><h3>\\[\\mathrm{fixture} \\rightarrow \\mathrm{fixture}\\]</h3><table><tr><td>Columns</td><td>native typed units</td></tr></table><table class="fit-coeff-table">' + rows + '</table>'


class NativeCollisionMetadata(unittest.TestCase):
    def test_exact_decimal_and_order(self):
        value = N.parse_iaea_fit_page(page(), "D1")
        self.assertEqual(value['parameters'], [('pcf(1)', '+2.500e+01'), ('A2', '-1.0000e-12')])
        self.assertIn(b'+2.500e+01', N.parameters_csv(value))
        self.assertFalse(value['physical_activation'])

    def test_duplicate_labels_refused(self):
        with self.assertRaises(N.NativeCollisionSourceError):
            N.parse_iaea_fit_page(page('<tr><td>A</td><td>1</td></tr><tr><td>A</td><td>2</td></tr>'))

    def test_nonfinite_refused(self):
        for token in ('NaN', 'Infinity', '-Infinity'):
            with self.subTest(token=token), self.assertRaises(N.NativeCollisionSourceError):
                N.parse_iaea_fit_page(page('<tr><td>A</td><td>' + token + '</td></tr>'))

    def test_invalid_decimal_refused(self):
        with self.assertRaises(N.NativeCollisionSourceError):
            N.parse_iaea_fit_page(page('<tr><td>A</td><td>not-a-number</td></tr>'))

    def test_missing_or_extra_cells_refused(self):
        for cells in ('<td>A</td>', '<td>A</td><td>1</td><td>extra</td>'):
            with self.subTest(cells=cells), self.assertRaises(N.NativeCollisionSourceError):
                N.parse_iaea_fit_page(page('<tr>' + cells + '</tr>'))

    def test_wrong_dataset_identity_refused(self):
        with self.assertRaises(N.NativeCollisionSourceError):
            N.parse_iaea_fit_page(page(), 'D2')

    def test_error_document_refused(self):
        with self.assertRaises(N.NativeCollisionSourceError):
            N.parse_iaea_fit_page('<title>HTTP500 OSError</title>')

    def test_missing_and_duplicate_table_refused(self):
        for data in (page().replace('fit-coeff-table', 'other-table'), page() + '<table class="fit-coeff-table"></table>'):
            with self.subTest(data=data), self.assertRaises(N.NativeCollisionSourceError):
                N.parse_iaea_fit_page(data)

    def test_missing_reaction_refused(self):
        with self.assertRaises(N.NativeCollisionSourceError):
            N.parse_iaea_fit_page(page().replace('\\rightarrow', 'no-arrow'))

    def test_missing_formula_is_explicit_unqualified(self):
        result = N.parse_iaea_fit_page(page())
        self.assertIn('NATIVE_FIT_FORMULA_NOT_RENDERED', result['missing_metadata'])
        self.assertIn('NATIVE_FIT_DOMAIN_NOT_SUPPLIED', result['missing_metadata'])
        self.assertEqual(result['scientific_status'], 'SOURCE_INPUT_NOT_PHYSICAL_QUALIFICATION')
        self.assertEqual(result['native_metadata']['columns'], 'native typed units')

    def test_sha_refusal_precedes_output(self):
        with tempfile.TemporaryDirectory() as name:
            source, output = Path(name) / 'source.html', Path(name) / 'out.csv'
            source.write_text(page())
            with self.assertRaises(N.NativeCollisionSourceError):
                N.project(source, '0' * 64, output, 'D1')
            self.assertFalse(output.exists())

    def test_malformed_sha_refused(self):
        with tempfile.TemporaryDirectory() as name:
            source, output = Path(name) / "source.html", Path(name) / "out.csv"
            source.write_text(page())
            for digest in (None, 123, "bad"):
                with self.subTest(digest=digest), self.assertRaises(N.NativeCollisionSourceError):
                    N.project(source, digest, output, "D1")
            self.assertFalse(output.exists())

    def test_reproduction_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as name:
            source, output = Path(name) / 'source.html', Path(name) / 'out.csv'
            source.write_text(page())
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            result = N.project(source, digest, output, 'D1')
            self.assertEqual(result['row_count'], 2)
            self.assertEqual(output.read_bytes(), N.parameters_csv(N.parse_iaea_fit_page(page())))
            before = output.read_bytes()
            with self.assertRaises(FileExistsError):
                N.project(source, digest, output, 'D1')
            self.assertEqual(output.read_bytes(), before)


if __name__ == '__main__':
    unittest.main(verbosity=2)
