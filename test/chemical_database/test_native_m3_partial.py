"""Check published r8 values against the original workbooks and reports."""

import os
from types import SimpleNamespace
import unittest
from xml.etree import ElementTree
from zipfile import ZipFile

import test_protocol


class PartialAtomicRelease(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_PARTIAL_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit r8 atomic-source test not selected")
        cls.database = test_protocol.CONTEXT.python_database()

    def test_workbook_all_projected_cells_match_native_xml(self):
        # Read native cells independently of the publisher's XLSX decoder.
        # Compare the complete workbook hash as well as each projected value.
        workbook = self.database.source_path("kalha-yeh-lindau-v2:workbook")
        self.assertEqual(test_protocol.digest(workbook),
                         "95770ea450655b6a26598dfe6b34a64787a558fdb933c1c13db2329da86570b5")
        table = self.database.table("kalha-yeh-lindau-v2:he-ne-na-ar")
        self.assertEqual(table.row_count, 192)
        with ZipFile(workbook) as archive:
            sheets = ElementTree.fromstring(archive.read("xl/workbook.xml")).findall(".//{*}sheet")
            relationships = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
            worksheet_targets = {entry.attrib["Id"]: entry.attrib["Target"] for entry in relationships}
            relationship_id = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            names = ("2_He", "10_Ne", "11_Na", "18_Ar")
            cell_values = {}
            for sheet in sheets:
                name = sheet.attrib["name"]
                if name in names:
                    target = worksheet_targets[sheet.attrib[relationship_id]]
                    xml = ElementTree.fromstring(archive.read("xl/" + target))
                    cell_values[name] = {cell.attrib["r"]: cell.findtext("{*}v", "")
                                         for cell in xml.findall(".//{*}c")}
            shared_strings = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            words = ["".join(entry.itertext()) for entry in shared_strings.findall("{*}si")]
        columns = {"2_He": "B", "10_Ne": "BCD", "11_Na": "BCDE", "18_Ar": "BCDE"}
        self.assertEqual({(row[4], row[5], row[6]) for row in table.rows},
                         {(sheet, f"A{index}", f"{letter}{index}")
                          for sheet, letters in columns.items()
                          for index in range(2, 18) for letter in letters})
        for symbol, orbital, energy, sigma, sheet, energy_cell, sigma_cell in table.rows:
            values = cell_values[sheet]
            self.assertEqual(symbol, sheet.split("_", 1)[1])
            self.assertEqual(orbital, words[int(values[sigma_cell[0] + "1"])])
            self.assertEqual((energy, sigma), (values[energy_cell], values[sigma_cell]))
        self.assertIn(("He", "1s", "21.2", "", "2_He", "A4", "B4"), table.rows)

    def test_sutton_ratios_uncertainty_and_cpp_reader(self):
        # Values independently read from the printed page and Table 2.
        identifier = "sutton1972-he:excitation-ratios"
        table = self.database.table(identifier)
        self.assertEqual(table.rows, (
            ("200", "4p", "4s", "19", "30", "17", "19", "21"),
            ("200", "4d", "4s", "0.85", "20", "0.81", "0.85", "0.88"),
            ("200", "4f", "4s", "0.78", "10", "0.73", "0.78", "0.79"),
        ))
        identity = dict(self.database.identity)
        fixture = SimpleNamespace(root=self.database.root, fields=identity, export_path=identity["export_path"])
        result = test_protocol.CONTEXT.native(fixture, "table", identifier, "table-access")
        self.assertTrue(result["open_ok"], result)
        self.assertTrue(result["query_ok"], result)
        self.assertEqual(result["rows"], [list(row) for row in table.rows])

    def test_bhatia_table_only_author_hybrid_values(self):
        table = self.database.table("bhatia2022-he:hybrid-table2")
        self.assertEqual(table.rows, (
            ("0.1", "7.3300"), ("0.2", "7.1544"), ("0.3", "6.8716"),
            ("0.4", "6.4951"), ("0.5", "6.0461"), ("0.6", "5.5925"),
            ("0.7", "5.0120"), ("0.8", "4.4740"), ("0.9", "3.9649"),
            ("1.0", "3.4654"), ("1.1", "3.0206"), ("1.3", "2.2561"),
            ("1.4", "1.9821"), ("1.5", "1.6817"), ("1.6", "1.6329"),
        ))

    def test_bhatia_excited_ratios_preserve_printed_exponents(self):
        table = self.database.table("bhatia2022-he:excited-table8")
        # The official HTML confirms the two decade jumps. Preserve the
        # source values and leave their physical interpretation open.
        self.assertEqual(table.rows, (
            ("0.1", "8.0146", "−3", "1.2794", "−2"),
            ("0.2", "8.6667", "−3", "1.3318", "−2"),
            ("0.3", "9.6384", "−3", "1.4112", "−1"),
            ("0.4", "1.0774", "−2", "1.5065", "−1"),
            ("0.5", "1.1935", "−2", "1.6097", "−2"),
            ("0.6", "1.3075", "−2", "1.7206", "−2"),
            ("0.7", "1.4252", "−2", "1.8497", "−2"),
            ("0.8", "1.5551", "−2", "2.0062", "−2"),
        ))

    def test_embargoed_archive_and_water_lead_are_metadata_only(self):
        archive = "wang2026-he:sciencedb-metadata"
        dataset = self.database.describe("m3-atomic-partial-20261124")
        row, = (row for row in dataset["files"] if row["id"] == archive)
        self.assertEqual(row["quantity_metadata"]["archive_conditions_of_access"], "EMBARGO")
        self.assertFalse(row["physical_activation"])
        for identifier in (archive, "suarez2018-water:publication-metadata"):
            with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as refused:
                self.database.source_path(identifier)
            self.assertEqual(refused.exception.code, "PermissionRequired")
