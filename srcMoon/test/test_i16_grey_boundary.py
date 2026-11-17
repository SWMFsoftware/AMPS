#!/usr/bin/env python3
"""I16 receipt/arithmetic fixtures only; no lunar observations are read.

These tests prove that the scorer uses explicitly labeled brightness and
rejects inconsistent native columns. They do not qualify D03 or I16 physics.
"""
import math
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys

from integration_common import i16_model_brightness_k,prepare_i16_input,sha256
sys.path.insert(0,str(Path(__file__).resolve().parent / "common"))
import moon_testlib


class GreyReceipt(unittest.TestCase):
    def row(self, kinetic=320.0, emissivity=81.0 / 256.0):
        return {
            "production_kinetic_temperature_K": str(kinetic),
            "production_bolometric_brightness_temperature_K": str(
                emissivity ** 0.25 * kinetic),
            "bolometric_emissivity": str(emissivity),
        }

    def test_independent_exact_three_quarter_root(self):
        # Fourth root(81/256)=3/4: the comparable observation is 240 K,
        # while the callback's 320 K is a different physical quantity.
        self.assertEqual(i16_model_brightness_k(self.row(),81.0/256.0),240.0)

    def test_unit_emissivity(self):
        self.assertEqual(i16_model_brightness_k(self.row(173.0,1.0),1.0),173.0)

    def test_direct_kinetic_column_cannot_pass(self):
        row=self.row()
        row["production_bolometric_brightness_temperature_K"]="320"
        with self.assertRaises(ValueError):
            i16_model_brightness_k(row,81.0/256.0)

    def test_stale_ambiguous_receipt_is_rejected(self):
        with self.assertRaises(KeyError):
            i16_model_brightness_k({"production_temperature_k":"320"},0.81)

    def test_row_emissivity_must_match_native_receipt(self):
        with self.assertRaises(ValueError):
            i16_model_brightness_k(self.row(),0.81)

    def test_invalid_temperatures(self):
        for field in ("production_kinetic_temperature_K",
                      "production_bolometric_brightness_temperature_K"):
            for value in (0.0,-1.0,math.nan,math.inf):
                with self.subTest(field=field,value=value):
                    row=self.row();row[field]=str(value)
                    with self.assertRaises(ValueError):
                        i16_model_brightness_k(row,81.0/256.0)

    def test_invalid_emissivity(self):
        for value in (0.0,-0.1,1.01,math.nan,math.inf):
            with self.subTest(value=value):
                row=self.row();row["bolometric_emissivity"]=str(value)
                with self.assertRaises(ValueError):
                    i16_model_brightness_k(row,value)


class StagingReceipt(unittest.TestCase):
    def stage(self, explicit):
        # These bytes test only the staging protocol's hash binding. They are
        # not scientific records, not canonical D03, and are never scored or
        # supplied to an initialized model. READY is a temporary parser fixture.
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root=Path(temporary.name)
        data=root/"_m2_validation_data/D03-GCP";data.mkdir(parents=True)
        table=data/"temperature_table.csv";table.write_text("PROTOCOL ONLY\n")
        benchmark=data/"benchmark_points.csv";benchmark.write_text("NO OBSERVATIONS\n")
        (data/"d03_gcp_qa.json").write_text(json.dumps({
            "status":"READY","processed":{
                "temperature_table":{"sha256":sha256(table)},
                "benchmark_points":{"sha256":sha256(benchmark)}}}))
        text=("#begin application: moon\n"
              "surface_temperature_mode = thermal-model\n"
              "diviner_temperature_table_file = none\n"
              "thermal_infrared_emissivity = 0.81\n")
        if explicit:
            text+="diviner_map_kinetic_interpretation = isothermal-grey\n"
        text+="#end application\n"
        template=root/"template.in";template.write_text(text)
        destination=root/"application.in"
        prepare_i16_input(template,destination,root)
        return destination,json.loads((root/"diviner_boundary_interpretation.json").read_text())

    def test_legacy_template_is_explicitly_staged_with_policy(self):
        destination,receipt=self.stage(False)
        self.assertIn("diviner_map_kinetic_interpretation = isothermal-grey",destination.read_text())
        self.assertIn("Approved grey-facet boundary",destination.read_text())
        self.assertEqual(receipt["staged_input_sha256"],sha256(destination))
        self.assertEqual(receipt["comparison_quantity"],"bolometric_brightness_temperature_K")

    def test_existing_explicit_declaration_is_not_duplicated(self):
        destination,receipt=self.stage(True)
        self.assertEqual(destination.read_text().count("diviner_map_kinetic_interpretation ="),1)
        self.assertEqual(receipt["interpretation"],"isothermal-grey")


class FailedProbeProvenance(unittest.TestCase):
    def test_furnished_bytes_hashed_even_when_metric_child_fails(self):
        # This executable-shaped child is a control-plane fixture, not an
        # AMPS executable and not a source of any physical reference. Only the
        # compile adapter is substituted; the real runner executes the child,
        # reads its inventory, hashes those bytes and retains its FAIL status.
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            kernel=root/"protocol kernel bytes.dat"
            kernel.write_bytes(b"PROTOCOL FIXTURE; NOT A SPICE KERNEL\n")
            child=root/"failing-protocol-child"
            child.write_text("#!"+sys.executable+"\n"
                "print('DEPENDENCY_SPICE_COUNT 1')\n"
                "print("+repr("DEPENDENCY_SPICE_FILE "+str(kernel))+" )\n"
                "print('FAIL CONTROL FIXTURE; not a physical metric')\n"
                "raise SystemExit(1)\n")
            child.chmod(0o700)
            with patch.object(moon_testlib,"compile_probe",return_value=(child,{"status":"PASS","scope":"protocol-only-no-build"})):
                result=moon_testlib.run_probe_case(
                    {"probe_modes":["prescribed-temperature"]},root,root/"output")
            self.assertEqual(result["status"],"FAIL")
            artifact=root/"output/U10/kernel_hashes.json"
            record=json.loads(artifact.read_text())
            self.assertEqual(record["files"],[{"path":str(kernel),"sha256":sha256(kernel)}])
            self.assertEqual(result["probes"][0]["returncode"],1)


if __name__=="__main__":
    unittest.main()
