#!/usr/bin/env python3
"""Independent toy-byte/selection controls; never native observations.

The 5x5 signed-little-endian images test decoding, explicit units/projection,
padding and gap handling. Frozen production hashes are replaced ONLY through
the Python fixture hook; the public CLI cannot select these artificial bytes.
All files remain temporary and no model or observational score is produced.
"""
import csv
import json
from pathlib import Path
import struct
import tempfile
import unittest

import prepare_i15_spatial_support as processor


class NativeRegionalPreparation(unittest.TestCase):
    def fixture(self, root):
        data=root/"data";package=data/"lola";raw=package/"raw/illumination";raw.mkdir(parents=True)
        (package/"README.md").write_text("CONTROL FIXTURE - NOT AUTHORITATIVE DATA")
        (package/"provenance.json").write_text('{"control_fixture":true}')
        (package/"qa_report.json").write_text('{"dataset_status":{"D02":"PARTIAL"}}')
        psr="LPSR_65N_240M_201608";visibility="AVGVISIB_65N_240M_201608"
        expected={}
        for product,multiplier in ((psr,1),(visibility,.00004)):
            label=f'''PRODUCT_ID = {product}
LINES = 5
LINE_SAMPLES = 5
SAMPLE_BITS = 16
SAMPLE_TYPE = LSB_INTEGER
MAP_PROJECTION_TYPE = "POLAR STEREOGRAPHIC"
COORDINATE_SYSTEM_NAME = "MEAN EARTH/POLAR AXIS OF DE421"
POSITIVE_LONGITUDE_DIRECTION = EAST
KEYWORD_LATITUDE_TYPE = PLANETOCENTRIC
UNIT = NONE
MAP_SCALE = 240 <m/pix>
A_AXIS_RADIUS = 1737.4 <km>
B_AXIS_RADIUS = 1737.4 <km>
C_AXIS_RADIUS = 1737.4 <km>
CENTER_LATITUDE = 90 <deg>
CENTER_LONGITUDE = 0 <deg>
LINE_PROJECTION_OFFSET = 2
SAMPLE_PROJECTION_OFFSET = 2
SCALING_FACTOR = {multiplier}
OFFSET = 0
MINIMUM_LATITUDE = 65 <deg>
MAXIMUM_LATITUDE = 90 <deg>
'''
            (raw/(product+".LBL")).write_text(label)
            values=[12500]*25 if product==visibility else [1]+[0]*24
            (raw/(product+".IMG")).write_bytes(struct.pack("<25h",*values))
            for suffix in (".LBL",".IMG"):expected[product+suffix]=processor.digest(raw/(product+suffix))
        provenance=root/"time-support-control.txt";provenance.write_text("PROTOCOL FIXTURE; NO REAL TIME CONTRACT")
        value=dict(schema=processor.SCHEMA,status="APPROVED",frozen_before_model_evaluation=True,
            review=dict(approved=True,reviewer="PROTOCOL FIXTURE",approved_utc="2026-01-01T00:00:00Z"),
            area_measure="native-polar-stereographic-map-plane",selection_model_independent=True,
            selection_basis="fixed control rectangle",time_support_id="CONTROL-ONLY",
            time_support_provenance={"path":str(provenance),"sha256":processor.digest(provenance)},
            domains=[dict(native_product_id=psr,visibility_product_id=visibility,source_id="CONTROL-PAIR",
                line_min=1,line_max=5,sample_min=1,sample_max=5,
                scoring_line_min=2,scoring_line_max=4,scoring_sample_min=2,scoring_sample_max=4)])
        contract=root/"selection.json";contract.write_text(json.dumps(value))
        return data,contract,expected,value,raw

    def test_signed_native_scale_projection_and_boundary_are_reproduced(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
            result=processor.prepare(contract,data,root/"output",expected_hashes=expected)
            self.assertEqual(result["status"],"PASS")
            self.assertEqual(result["formal_i15_status"],"SKIPPED")
            with (root/"output/site_reference.csv").open() as stream:rows=list(csv.DictReader(stream))
            self.assertEqual(len(rows),25)
            self.assertTrue(all(float(r["reference_average_occurrence"])==.5 for r in rows))
            self.assertTrue(all(float(r["pixel_area_m2"])==57600 for r in rows))
            centre=next(r for r in rows if r["native_line"]=="3" and r["native_sample"]=="3")
            self.assertEqual(float(centre["latitude_deg"]),90)
            near=next(r for r in rows if r["native_line"]=="2" and r["native_sample"]=="2")
            self.assertEqual(near["psr_boundary_ambiguous"],"1")
            self.assertEqual(float(near["psr_resolution_allowance_m2"]),57600)
            self.assertEqual(centre["logical_surface_element"],"")
            provenance=json.loads((root/"output/provenance.json").read_text())
            self.assertTrue(provenance["control_fixture_not_science"])
            self.assertFalse((root/"output/comparison.png").exists())

    def test_unreviewed_or_missing_time_metadata_skips_before_native_decode(self):
        for mode in ("review","timing"):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
                if mode=="review":value["review"]["approved"]=False
                else:value.pop("time_support_provenance")
                contract.write_text(json.dumps(value))
                result=processor.prepare(contract,data,root/"output",expected_hashes=expected)
                self.assertEqual(result["status"],"SKIPPED")
                self.assertFalse((root/"output/site_reference.csv").exists())

    def test_unknown_area_convention_cannot_become_true_surface_area(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
            value["area_measure"]="true-spherical-area-by-squaring-map-scale"
            contract.write_text(json.dumps(value))
            self.assertEqual(processor.prepare(contract,data,root/"output",expected_hashes=expected)["status"],"SKIPPED")

    def test_changed_native_bytes_are_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
            (raw/(value["domains"][0]["native_product_id"]+".IMG")).write_bytes(b"changed")
            result=processor.prepare(contract,data,root/"output",expected_hashes=expected)
            self.assertEqual(result["status"],"ERROR");self.assertIn("hash mismatch",result["reason"])

    def test_native_uninterpreted_value_is_preserved_as_gap_not_filled(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
            path=raw/(value["domains"][0]["visibility_product_id"]+".IMG")
            path.write_bytes(struct.pack("<25h",-5,*([12500]*24)));expected[path.name]=processor.digest(path)
            result=processor.prepare(contract,data,root/"output",expected_hashes=expected)
            self.assertEqual(result["status"],"SKIPPED")
            with (root/"output/site_reference.csv").open() as stream:row=next(csv.DictReader(stream))
            self.assertEqual(row["native_visibility_dn"],"-5")
            self.assertEqual(row["reference_average_occurrence"],"")
            self.assertEqual(row["quality"],"NATIVE_VALUE_UNINTERPRETED")

    def test_regions_need_padding_and_cannot_overlap(self):
        for mode in ("padding","overlap"):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
                if mode=="padding":value["domains"][0]["scoring_line_min"]=1
                else:value["domains"].append(value["domains"][0].copy())
                contract.write_text(json.dumps(value))
                result=processor.prepare(contract,data,root/"output",expected_hashes=expected)
                self.assertEqual(result["status"],"ERROR")

    def test_canonical_output_and_existing_evidence_are_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,contract,expected,value,raw=self.fixture(root)
            with self.assertRaisesRegex(ValueError,"outside"):
                processor.prepare(contract,data,data/"processed",expected_hashes=expected)
            out=root/"output";out.mkdir();(out/"preserve.txt").write_text("unchanged")
            with self.assertRaisesRegex(ValueError,"overwrite"):
                processor.prepare(contract,data,out,expected_hashes=expected)
            self.assertEqual((out/"preserve.txt").read_text(),"unchanged")


if __name__=="__main__":unittest.main(verbosity=2)
