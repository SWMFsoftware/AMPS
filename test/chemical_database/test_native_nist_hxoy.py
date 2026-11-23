"""Parser controls and explicitly selected authoritative NIST article checks.

Administrative HTML fixtures below are never observations or thermochemical
reference values.  The actual-source class is SKIPPED unless the caller supplies
AMPS_NIST_HXOY_SOURCE and AMPS_NIST_HXOY_SHA256; it verifies source/view identity,
not physical accuracy, uncertainty selection or model qualification.
"""
import csv
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("native_nist_hxoy", ROOT / "utility/chemical_database/native_nist_hxoy.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture():
    return ('PMC7339710 10.6028/jres.121.005' + ''.join(
        '<section id="' + name + '"><table><tr><th rowspan="2">E<sup>−1</sup></th>'
        '<td colspan="2">A &amp; B<a href="#fn">[a]</a></td></tr>'
        '<tr><td>  literal\n space </td></tr></table></section>' for name in module.TABLES))


class ParserControls(unittest.TestCase):
    def test_exact_cells_spans_and_footnotes(self):
        rows = module.cells(fixture(), "T3")
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0][3:7], ("th", "2", "", "E−1"))
        self.assertEqual(rows[1][5:], ("2", "A & B[a]", 'A &amp; B<a href="#fn">[a]</a>'))
        self.assertEqual(rows[2][6], "  literal\n space ")

    def test_missing_duplicate_and_malformed_refused(self):
        for value in (fixture().replace('id="T3"', 'id="missing"'), fixture() + fixture(),
                      fixture().replace('</td>', '', 1), fixture().replace('rowspan="2"', 'rowspan="0"'),
                      fixture().replace('</table>', '', 1), fixture().replace('rowspan="2"', 'rowspan=2')):
            with self.subTest(value=value[:60]), self.assertRaises(ValueError):
                module.cells(value, "T3")

    def test_transactional_hash_and_wrong_article_guards(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); source = base / 'source.html'; source.write_text(fixture())
            output = base / 'output'
            with self.assertRaises(ValueError): module.prepare(source, '0' * 64, output)
            self.assertFalse(output.exists())
            source.write_text(fixture().replace('PMC7339710', 'unrelated'))
            with self.assertRaises(ValueError): module.prepare(source, module.sha256(source), output)
            self.assertFalse(output.exists())

    def test_csv_roundtrip_and_nonoverwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); source = base / 'source.html'; source.write_text(fixture())
            out = base / 'view'; module.prepare(source, module.sha256(source), out)
            with (out / 'T3_cells.csv').open(newline='') as stream:
                rows = list(csv.reader(stream))
            self.assertEqual(tuple(rows[0]), module.COLUMNS)
            self.assertEqual(tuple(map(tuple, rows[1:])), module.cells(fixture(), 'T3'))
            before = (out / 'T3_cells.csv').read_bytes()
            with self.assertRaises(ValueError): module.prepare(source, module.sha256(source), out)
            self.assertEqual((out / 'T3_cells.csv').read_bytes(), before)


@unittest.skipUnless(os.environ.get('AMPS_NIST_HXOY_SOURCE') and os.environ.get('AMPS_NIST_HXOY_SHA256'),
                     'explicit authoritative source/hash unavailable')
class AuthoritativeSource(unittest.TestCase):
    requires_authoritative_source = True

    def test_native_cells_exactly_retained(self):
        source = Path(os.environ['AMPS_NIST_HXOY_SOURCE'])
        expected = os.environ['AMPS_NIST_HXOY_SHA256']
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), expected)
        document = source.read_text()
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'views'
            receipt = module.prepare(source, expected, out)
            for item in receipt['tables']:
                with (out / item['path']).open(newline='') as stream:
                    rows = list(csv.DictReader(stream))
                # Independently inspect the original table region and inner cell
                # markup.  No energy interpretation is copied from the adapter.
                section = document.split('id="' + item['table_id'] + '"', 1)[1]
                body = section.split('<table', 1)[1].split('>', 1)[1].split('</table>', 1)[0]
                original = re.findall(r'<(?:td|th)\b[^>]*>(.*?)</(?:td|th)>', body, re.S)
                self.assertEqual([r['native_html'] for r in rows], original)
                self.assertEqual(module.sha256(out / item['path']), item['sha256'])
                self.assertGreater(len(rows), 20)


if __name__ == '__main__':
    unittest.main()
