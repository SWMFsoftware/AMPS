#!/usr/bin/env python3
"""I15 metadata/control/arithmetic tests, never observational evidence.

Every artificial file is confined to TemporaryDirectory. The fake executable
is explicitly non-executable text, no MPI or physics is invoked, and plots
are disabled. The assertions test refusal of contradictory provenance and
exact arithmetic joins; a fixture's result is not an I15 qualification.
"""
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import m2_illumination as illumination


class IlluminationContract(unittest.TestCase):
    def fixture(self, root, *, earth_eclipsed_site=None):
        """Return a fully connected CONTROL-PLANE fixture, not lunar data."""
        record=lambda p:{"path":str(p),"sha256":illumination.digest(p)}
        def document(name, value):
            p=root/name;p.write_text(json.dumps(value));return p
        def table(name, data):
            p=root/name
            with p.open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator="\n")
                writer.writeheader();writer.writerows(data)
            return p
        doc=root/"NOT-AUTHORITATIVE.txt";doc.write_text("PROTOCOL FIXTURE; NO OBSERVATIONS OR PHYSICS\n")
        exe=root/"NOT-AN-EXECUTABLE";exe.write_text("NEVER EXECUTE THIS CONTROL FIXTURE\n")
        for name in ("terrain.img","terrain.lbl"):(root/name).write_text("PROTOCOL ONLY")
        path=root/"reference.json";template=root/"source.in"
        template.write_text("\n".join(["#begin application: moon",
            "m2_formal_contract_file = "+str(path),"surface_geometry = lola",
            "surface_mesh_resolution_m = 1","lola_product_id = FIXTURE",
            "lola_image_file = "+str(root/"terrain.img"),"lola_label_file = "+str(root/"terrain.lbl"),
            "spice_path = "+str(root),"illumination_ray_origin_offset_m = 0.01",
            "illumination_terrain_solar_model = uniform-disc","illumination_aberration_correction = NONE",
            "illumination_disc_radial_points = 2","illumination_disc_azimuth_points = 8",
            "thermal_solar_disc_radius_m = 10","thermal_opaque_earth_radius_m = 1",
            "surface_cea_file = unused.cea","surface_tecplot_file = unused.dat",
            "surface_inventory_checkpoint_file = none","#end application: moon",""]))
        queries=[];sites=[];models=[]
        for line in range(1,4):
            for sample in range(1,4):
                site=f"CONTROL-{line}-{sample}"
                eclipsed=site==earth_eclipsed_site
                sites.append(dict(site_id=site,reference_average_occurrence="0" if eclipsed else "0.5",reference_psr="0",
                    psr_boundary_ambiguous="0",pixel_area_m2="1",psr_resolution_allowance_m2="0",
                    logical_surface_element="0",time_support_id="CONTROL-ONLY",native_product_id="FIXTURE",
                    native_line=str(line),native_sample=str(sample),source_id="CONTROL-SOURCE",
                    latitude_deg=str(line),longitude_east_deg=str(sample)))
                for epoch in (0,1):
                    key=f"{site}-t{epoch}"
                    queries.append(dict(zip(illumination.GEOMETRY_FIELDS,
                        (key,site,str(line),str(sample),str(epoch),"TDB","0","radial",
                         "true" if earth_eclipsed_site else "false"))))
                    models.append(dict(query_id=key,site_id=site,ephemeris_time_tdb_s=str(epoch),
                        input_epoch=str(epoch),input_time_system="TDB",latitude_deg=str(line),longitude_east_deg=str(sample),
                        terrain_visible_area_fraction="0.2" if (epoch or eclipsed) else "0",
                        earth_visible_area_fraction="0" if eclipsed else "1",
                        jointly_visible_area_fraction="0.2" if (epoch and not eclipsed) else "0",
                        projected_incidence="0.1" if (epoch and not eclipsed) else "0",
                        sampled_any_terrain_visible="1" if (epoch or eclipsed) else "0",
                        sampled_any_jointly_visible="0" if eclipsed else str(epoch),ray_count="16"))
        geometry=table("geometry.csv",queries)
        epochs=table("epoch-control.csv",[dict(query_id=r["query_id"],reference_visible=r["sampled_any_jointly_visible"],
                                                source_id="CONTROL-SOURCE") for r in models])
        site_table=table("site-control.csv",sites)
        timing_site="CONTROL-1-1" if earth_eclipsed_site else "CONTROL-2-2"
        timing=table("timing-control.csv",[dict(event_id="CONTROL-EVENT",before_query_id=timing_site+"-t0",
            after_query_id=timing_site+"-t1",reference_epoch_tdb_s="0.5",event_kind="sunrise",
            maximum_timing_error_s="1",source_id="CONTROL-SOURCE")])
        support=document("spatial-control.json",dict(status="READY",frozen_before_model_evaluation=True,
            scope="resolved-regional-native-cell-coverage",area_measure="native-polar-stereographic-map-plane",
            native_pixel_width_m=1,native_pixel_height_m=1,timing_tolerances_s={"CONTROL-EVENT":1},
            domains=[dict(native_product_id="FIXTURE",line_min=1,line_max=3,sample_min=1,sample_max=3,
                scoring_line_min=2,scoring_line_max=2,scoring_sample_min=2,scoring_sample_max=2)]))
        geometry_hash=hashlib.sha256(illumination.geometry_bytes(queries)).hexdigest()
        assessment=document("assessment-control.json",dict(status="PASS",any_disc_sliver_error_qualified=True,
            temporal_event_brackets_qualified=True,geometry_queries_sha256=geometry_hash,
            method="closed-form-bound",evidence=[record(doc)]))
        runs=[]
        for radial,azimuth in ((1,4),(2,8)):
            runtime=root/f"run-{radial}.in"
            runtime.write_text(template.read_text().replace("radial_points = 2",f"radial_points = {radial}").replace(
                "azimuth_points = 8",f"azimuth_points = {azimuth}"))
            receipt=root/f"receipt-{radial}.txt"
            receipt.write_text(f"test_id=I15\nmpi_ranks=8\nquery_count={len(queries)}\ndisc_radial_points={radial}\ndisc_azimuth_points={azimuth}\n")
            manifest=document(f"manifest-{radial}.json",dict(test_id="I15",mpi_ranks=8,openmp_threads=1,
                executable={"sha256":illumination.digest(exe)},inputs={str(runtime):illumination.digest(runtime),str(geometry):geometry_hash}))
            outputs=table(f"output-{radial}.csv",[dict(row,ray_count=str(radial*azimuth)) for row in models])
            runs.append(dict(receipt=record(receipt),manifest=record(manifest),model_outputs=[record(outputs)]))
        proof=document("resolution-control.json",dict(status="PASS",angular_and_temporal_resolution_qualified=True,
            source_input_sha256=illumination.digest(template),executable_sha256=illumination.digest(exe),
            production_source_signature_sha256=hashlib.sha256(b"{}").hexdigest(),geometry_queries_sha256=geometry_hash,
            disc_radial_points=2,disc_azimuth_points=8,native_runs=runs,independent_resolution_assessment=record(assessment)))
        operator=dict(frame="MOON_ME_DE421",time_system="TDB",sample_count=2,
            epoch_list_endpoint_convention="two EXPLICIT ARITHMETIC fixture epochs",cadence_description="fixture only",
            observer_height_m=0,height_direction="radial",aberration_correction="NONE",
            earth_occultation="opaque-sphere" if earth_eclipsed_site else "excluded",
            solar_radius_m=10,solar_radiance_model="uniform",time_support_id="CONTROL-ONLY")
        if earth_eclipsed_site:operator["opaque_earth_radius_m"]=1
        contract=dict(schema="srcMoon-I15-reference-contract-v1",status="READY",frozen_before_model_evaluation=True,
            scope="TEMPORARY CONTROL FIXTURE; NOT SCIENCE",selection_role="hold-out",independent_source_ids=["CONTROL-SOURCE"],
            review=dict(approved=True,reviewer="PROTOCOL FIXTURE",date_utc="2026-01-01T00:00:00Z"),
            operator=operator,operator_source_locators={k:dict(document_index=0,locator="ARITHMETIC FIXTURE ONLY") for k in operator},
            authoritative_documents=[record(doc)],geometry_queries=record(geometry),epoch_reference=record(epochs),
            site_reference=record(site_table),timing_reference=record(timing),resolution_qualification=record(proof),spatial_support=record(support))
        path.write_text(json.dumps(contract))
        return dict(root=root,path=path,template=template,exe=exe,contract=contract,proof=proof,geometry=geometry,
                    queries=queries,models=models,site_table=site_table,sites=sites,support=support,table=table,record=record)

    def stage(self, fixture):
        out=fixture["root"]/"stage output"
        with patch("build_m2_qualification.production_source_identity",return_value={}):
            illumination.stage_input(fixture["template"],out/"application.in",executable=fixture["exe"])
        fixture["table"]("stage output/receipt.txt.I15.rank0000.csv",fixture["models"])
        return out

    def score(self, fixture, out):
        receipt=dict(query_count=len(fixture["models"]),body_fixed_frame="MOON_ME_DE421",
                     aberration_correction="NONE",disc_radial_points=2,disc_azimuth_points=8)
        with patch("build_m2_qualification.production_source_identity",return_value={}):
            return illumination.score(out,receipt,executable=fixture["exe"],render=False)

    def test_missing_exact_metadata_is_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));f["contract"]["operator"]["cadence_description"]=None
            f["path"].write_text(json.dumps(f["contract"]))
            with self.assertRaises(illumination.Unavailable):illumination.load_contract(f["path"])

    def test_corrupt_authoritative_file_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));(f["root"]/"NOT-AUTHORITATIVE.txt").write_text("changed control bytes")
            with self.assertRaisesRegex(ValueError,"hash mismatch"):illumination.load_contract(f["path"])

    def test_same_site_cannot_change_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));f["queries"][1]["latitude_deg"]="2"
            f["table"]("geometry.csv",f["queries"]);f["contract"]["geometry_queries"]=f["record"](f["geometry"])
            f["path"].write_text(json.dumps(f["contract"]))
            with self.assertRaisesRegex(ValueError,"changing native geometry"):illumination.load_contract(f["path"])

    def test_equal_tdb_spellings_cannot_double_weight_an_epoch(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));f["queries"][1]["epoch"]="0.000"
            f["table"]("geometry.csv",f["queries"]);f["contract"]["geometry_queries"]=f["record"](f["geometry"])
            f["path"].write_text(json.dumps(f["contract"]))
            with self.assertRaisesRegex(ValueError,"duplicated"):illumination.load_contract(f["path"])

    def test_utc_identity_preserves_leap_spelling_without_posix_conversion(self):
        a=illumination.epoch_value(dict(time_system="UTC",epoch="2016-12-31T23:59:60.0Z"))
        b=illumination.epoch_value(dict(time_system="UTC",epoch="2016-12-31T23:59:60.000Z"))
        self.assertEqual(a,b)
        with self.assertRaises(ValueError):illumination.epoch_value(dict(time_system="UTC",epoch="2016-02-31T00:00:00Z"))

    def test_duplicate_headers_ragged_csv_and_duplicate_json_are_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"broken"
            for text in ("id,id\na,b\n","id,value\na\n","id,value\na,b,c\n"):
                p.write_text(text)
                with self.assertRaises(ValueError):illumination.rows(p,("id",))
            p.write_text('{"status":"PASS","status":"FAIL"}')
            with self.assertRaises(ValueError):illumination.json_document(p)

    def test_resolution_proof_must_bind_actual_input_executable_and_grid(self):
        for key in ("source_input_sha256","executable_sha256","geometry_queries_sha256","disc_radial_points"):
            with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
                f=self.fixture(Path(tmp));proof=json.loads(f["proof"].read_text());proof[key]=99 if key=="disc_radial_points" else "0"*64
                f["proof"].write_text(json.dumps(proof));f["contract"]["resolution_qualification"]=f["record"](f["proof"])
                f["path"].write_text(json.dumps(f["contract"]))
                with patch("build_m2_qualification.production_source_identity",return_value={}):
                    with self.assertRaisesRegex(ValueError,"resolution proof differs"):
                        illumination.stage_input(f["template"],Path(tmp)/"out/application.in",executable=f["exe"])

    def test_standalone_pass_boolean_is_not_resolution_qualification(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));f["proof"].write_text('{"status":"PASS","angular_and_temporal_resolution_qualified":true}')
            f["contract"]["resolution_qualification"]=f["record"](f["proof"]);f["path"].write_text(json.dumps(f["contract"]))
            with self.assertRaises(illumination.Unavailable):illumination.load_contract(f["path"])

    def test_rehashed_native_proof_still_needs_consistent_actual_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));proof=json.loads(f["proof"].read_text())
            bad=[dict(row,ray_count="999") for row in f["models"]]
            output=f["table"]("output-1.csv",bad)
            proof["native_runs"][0]["model_outputs"]=[f["record"](output)]
            f["proof"].write_text(json.dumps(proof));f["contract"]["resolution_qualification"]=f["record"](f["proof"])
            f["path"].write_text(json.dumps(f["contract"]))
            with patch("build_m2_qualification.production_source_identity",return_value={}):
                with self.assertRaisesRegex(ValueError,"ray count differs"):
                    illumination.stage_input(f["template"],Path(tmp)/"out/application.in",executable=f["exe"])

    def test_epoch_miss_is_fail_without_relaxing_99_percent_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp))
            values=illumination.rows(f["root"]/"epoch-control.csv",("query_id","reference_visible","source_id"))
            values[0]["reference_visible"]="1"
            file=f["table"]("epoch-control.csv",values)
            f["contract"]["epoch_reference"]=f["record"](file);f["path"].write_text(json.dumps(f["contract"]))
            out=self.stage(f);status,metrics,_=self.score(f,out)
            self.assertEqual(status,"FAIL")
            self.assertLess(metrics["epoch_visibility_agreement"],.99)

    def test_timing_failure_remains_separate_from_correct_epoch_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp))
            events=illumination.rows(f["root"]/"timing-control.csv",("event_id",))
            events[0]["maximum_timing_error_s"]="0.1"
            timing=f["table"]("timing-control.csv",events)
            support=json.loads(f["support"].read_text());support["timing_tolerances_s"]["CONTROL-EVENT"]=.1
            f["support"].write_text(json.dumps(support));f["contract"]["spatial_support"]=f["record"](f["support"])
            f["contract"]["timing_reference"]=f["record"](timing);f["path"].write_text(json.dumps(f["contract"]))
            out=self.stage(f);status,metrics,_=self.score(f,out)
            self.assertEqual(status,"FAIL");self.assertEqual(metrics["epoch_visibility_agreement"],1)
            self.assertFalse(metrics["timing_pass"])

    def test_arithmetic_occurrence_is_not_mean_area_and_quoted_paths_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));out=self.stage(f);status,metrics,_=self.score(f,out)
            self.assertEqual(status,"PASS") # arithmetic/protocol only, no native model was run
            data=illumination.rows(out/"illumination_model.csv",("any_disc_occurrence_fraction",))
            self.assertEqual(float(data[0]["any_disc_occurrence_fraction"]),.5)
            self.assertEqual(float(data[0]["mean_visible_area_fraction"]),.1)
            self.assertEqual(Path(illumination.setting(out/"application.in","m2_illumination_queries_file")),out/"illumination_geometry_queries.csv")
            self.assertFalse((out/"comparison.png").exists())
            from common.write_comparison_artifacts import read_comparison
            self.assertTrue(read_comparison(out/"comparison.csv"))

    def test_native_coordinates_and_epochs_must_match_exact_query(self):
        for key in ("latitude_deg","input_epoch","ephemeris_time_tdb_s"):
            with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
                f=self.fixture(Path(tmp));out=self.stage(f);f["models"][0][key]="99"
                f["table"]("stage output/receipt.txt.I15.rank0000.csv",f["models"])
                with self.assertRaises(ValueError):self.score(f,out)

    def test_native_area_boolean_mask_and_grid_inconsistency_are_error(self):
        for key,value in (("sampled_any_jointly_visible","1"),("sampled_any_terrain_visible","1"),
                          ("jointly_visible_area_fraction","0.5"),("ray_count","15")):
            with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
                f=self.fixture(Path(tmp));out=self.stage(f);f["models"][0][key]=value
                f["table"]("stage output/receipt.txt.I15.rank0000.csv",f["models"])
                with self.assertRaises(ValueError):self.score(f,out)

    def test_earth_eclipsed_history_does_not_create_geometric_psr(self):
        # This site is a predeclared Boolean/arithmetic control: terrain sees
        # the entire sampled solar history, while Earth's mask blocks it at
        # every listed epoch. No physical simulation/observation is generated.
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp),earth_eclipsed_site="CONTROL-2-2")
            out=self.stage(f);status,metrics,_=self.score(f,out)
            self.assertEqual(status,"PASS") # Control-plane arithmetic only.
            self.assertEqual(metrics["epoch_visibility_agreement"],1.)
            self.assertEqual(metrics["psr_interior_mismatch_count"],0)
            data=illumination.rows(out/"illumination_model.csv",("site_id",))
            site=next(row for row in data if row["site_id"]=="CONTROL-2-2")
            self.assertEqual(float(site["any_disc_occurrence_fraction"]),0.)
            self.assertEqual(float(site["joint_any_disc_occurrence_fraction"]),0.)
            self.assertEqual(float(site["terrain_any_disc_occurrence_fraction"]),1.)
            self.assertEqual(float(site["mean_visible_area_fraction"]),0.)
            self.assertAlmostEqual(float(site["mean_terrain_visible_area_fraction"]),.2)
            self.assertEqual(int(site["geometric_psr"]),0)
            self.assertEqual(site["earth_occultation_policy"],"opaque-sphere")
            self.assertEqual(site["geometric_psr_scope"],"terrain-dark-over-declared-sampled-epochs")
            self.assertFalse((out/"comparison.png").exists())

    def test_missing_terrain_flag_is_not_inferred_from_joint_or_area(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));out=self.stage(f)
            missing=[{key:value for key,value in row.items() if key!="sampled_any_terrain_visible"}
                     for row in f["models"]]
            f["table"]("stage output/receipt.txt.I15.rank0000.csv",missing)
            with self.assertRaisesRegex(ValueError,"missing illumination CSV fields"):
                self.score(f,out)

    def test_missing_temporal_support_cannot_publish_a_permanent_shadow_result(self):
        for field in ("epoch_list_endpoint_convention","sample_count","time_support_id"):
            with self.subTest(field=field),tempfile.TemporaryDirectory() as tmp:
                f=self.fixture(Path(tmp));f["contract"]["operator"][field]=None
                f["path"].write_text(json.dumps(f["contract"]))
                with self.assertRaises(illumination.Unavailable):self.stage(f)
                self.assertFalse((Path(tmp)/"stage output/illumination_model.csv").exists())

    def test_sparse_cells_cannot_qualify_psr_area(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp))
            with self.assertRaisesRegex(ValueError,"sparse/incomplete|complete support"):
                illumination.spatial_support(f["support"],f["sites"][:-1])

    def test_arbitrary_large_psr_allowance_cannot_hide_interior_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            f=self.fixture(Path(tmp));f["sites"][4]["psr_resolution_allowance_m2"]="100000000"
            with self.assertRaisesRegex(ValueError,"native-boundary"):
                illumination.spatial_support(f["support"],f["sites"])


if __name__=="__main__":unittest.main(verbosity=2)
