"""Stage-A5 reporting tests independent of mission kernels and SPICE."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))

from cops_fit_reporting import (REPORT_SCHEMA, identity_label, report_identity,
                                write_activity_products,
                                write_pressure_products,
                                write_report_descriptor)
from cops_io import file_identity
from plot_cops_model import Run


class FitReportingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.payload = {
            "profile": "all_available", "selection_sha256": "a" * 64,
            "enhancement_enabled": True, "enhancement_degree": 5,
            "legacy_degree": 4, "lambda": 1.0e-3,
        }
        self.identity = report_identity(self.payload)
        self.label = identity_label(self.payload)

    def tearDown(self):
        self.temporary.cleanup()

    def test_tables_figures_and_descriptor_share_one_identity(self):
        centroid = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                             [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
        area = np.arange(1.0, 7.0)
        baseline = np.linspace(0.6, 1.1, 6)
        residual = np.array([-0.8, -0.1, 0.0, 0.1, 0.2, 0.3])
        total = baseline + residual
        activity_files, summary = write_activity_products(
            self.root, self.identity, self.label, centroid, area,
            np.array([1, 1, 1, 0, 0, 0], dtype=bool), baseline, residual,
            total, np.maximum(total, 0.0),
            [(5, np.array([0.2, -0.1]), residual, 0.05)])

        pressure_files = write_pressure_products(
            self.root, self.identity, self.label, [10, 11, 12, 13],
            ["2016-09-30T10:22:00", "2016-09-30T10:23:55",
             "2016-09-30T10:24:20", "2016-09-30T10:25:00"],
            np.array([900.0, 500.0, 150.0, 80.0]),
            {"NG": [1.1, 2.2, 3.3, 4.4], "RG": [1.0, 2.0, 3.0, 4.0]},
            {"NG": [1.0, 2.0, 3.0, 4.0], "RG": [1.0, 2.0, 3.0, 4.0]},
            {"NG": [1, 0, 1, 1], "RG": [1, 1, 1, 1]},
            [0, 0, 1, 0],
            {"NG": ["selected", "transition", "selected", "selected"],
             "RG": ["selected"] * 4},
            [("2016-09-30T10:23:40", "2016-09-30T10:24:40")])
        products = activity_files + pressure_files
        products.append(write_report_descriptor(self.root, self.identity,
                                                self.payload, products))

        self.assertEqual(set(products), {
            "activity_facets.dat", "activity_degree_power.dat",
            "activity_template_power.dat",
            "activity_region_diagnostics.dat", "activity_maps.png",
            "pressure_residuals.dat", "pressure_comparison_time.png",
            "pressure_comparison_altitude.png", "pressure_residuals.png",
            "fit_report.json",
        })
        for name in ("activity_facets.dat", "activity_degree_power.dat",
                     "activity_template_power.dat",
                     "activity_region_diagnostics.dat", "pressure_residuals.dat"):
            header = (self.root / name).read_text().splitlines()[0]
            self.assertIn("ReportSchema=" + REPORT_SCHEMA, header)
            self.assertIn("ReportIdentity=" + self.identity, header)
        descriptor = json.loads((self.root / "fit_report.json").read_text())
        self.assertEqual(descriptor["report_identity"], self.identity)
        self.assertEqual(descriptor["identity_payload"], self.payload)
        self.assertEqual(summary["degree_count"], 1)

        rows = [line.split() for line in
                (self.root / "pressure_residuals.dat").read_text().splitlines()
                if line and not line.startswith("#")]
        ng = [row for row in rows if row[3] == "NG"]
        self.assertEqual(ng[1][4], "0")
        self.assertEqual(ng[1][8], "nan")
        self.assertEqual(ng[2][5], "1")
        self.assertTrue(np.isfinite(float(ng[2][8])))

    def test_plot_reader_verifies_fit_manifest_and_report_bytes(self):
        selection = "b" * 64
        report = self.root / "fit_report.json"
        report.write_text(json.dumps({
            "schema": REPORT_SCHEMA, "schema_version": 1,
            "report_identity": self.identity,
            "identity_payload": self.payload, "products": [],
        }) + "\n")
        table = self.root / "COPS_model.dat"
        table.write_text(
            "# ProductKind=fit_prediction TableSchema=liouville_cops_fit_prediction_v2\n"
            "# FitProfile=all_available SelectionSHA256=%s ReportSchema=%s ReportIdentity=%s\n"
            "# ActivityBasisID=liouville_cops.real_orthonormal_cs_v1 "
            "EnhancementDegree=5 Lambda=0.001 Regularization=test\n"
            "#iPoint UTC altitude inside pNG pRG refNG refRG useNG useRG\n"
            "0 2016-09-30T10:20:00 1000 0 1e-8 2e-8 1e-8 2e-8 1 1\n" %
            (selection, REPORT_SCHEMA, self.identity))
        manifest = {
            "schema": "liouville_cops.fit_manifest", "schema_version": 1,
            "states": {"process": "complete", "optimizer": "converged"},
            "selection": {"profile": "all_available", "sha256": selection},
            "activity_enhancement": {
                "enabled": True,
                "basis_id": "liouville_cops.real_orthonormal_cs_v1"},
            "optimizer_results": [{"result_directory": ".",
                                   "report_identity": self.identity}],
            "outputs": {"COPS_model.dat": file_identity(table),
                        "fit_report.json": file_identity(report)},
        }
        (self.root / "fit_manifest.json").write_text(json.dumps(manifest) + "\n")
        loaded = Run(str(table))
        self.assertEqual(loaded.report_identity, self.identity)
        self.assertIn("all_available", loaded.display_name())

        report.write_text(report.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "fit_report.json does not match"):
            Run(str(table))


if __name__ == "__main__":
    unittest.main()
