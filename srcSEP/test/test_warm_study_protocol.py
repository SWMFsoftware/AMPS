#!/usr/bin/env python3
"""Bounded software-protocol tests; manufactured records are not physics references."""
import copy
from decimal import Decimal
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.dont_write_bytecode=True
PHYSICS=Path(__file__).resolve().parents[1]/"validation/physics"
sys.path.insert(0,str(PHYSICS))
import run_warm_parallel_study as protocol
import plot_warm_parallel_study as plotter
import warm_parallel_maxwellian as warm
import mpmath as mp


class StudyProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storage=tempfile.TemporaryDirectory(prefix="srcsep-warm-protocol-tests-")
        cls.directory=Path(cls.storage.name)
        # Execute the real API for a deliberately different explicit N. A
        # hard-coded native128-grid copy cannot pass this capture/provenance gate.
        cls.extracted=protocol.capture_grid(cls.directory/"grid","c++",["4","1e-8","1e-2"])
        cls.input_bytes=(PHYSICS/"controlled_ep_maxwellian_baseline.json").read_bytes()
        cls.case=json.loads(cls.input_bytes)

    @classmethod
    def tearDownClass(cls):cls.storage.cleanup()

    def make_protocol(self,name="case",points=None):
        directory=self.directory/name;directory.mkdir(exist_ok=True)
        protocol.atomic_bytes(directory/"raw/input.json",self.input_bytes,True)
        result={"input":{"archive":"raw/input.json","sha256":protocol.digest(self.input_bytes)},
                "solver":protocol.archive(PHYSICS/"warm_parallel_maxwellian.py",directory),
                "base_digits":40,"verification_digits":70,"sampling":"explicit-kdi",
                "requested_helicities":["Right"],"requested_k_di":points or [".01",".02",".03"]}
        return directory,result

    def make_study(self,directory,prepared):
        # This fixture tests evidence schemas, missing-mode gaps and provenance
        # rejection. These numeric fields are intentionally manufactured and
        # never retained as a scientific reference or solver validation result.
        rows=[]
        for index,point in enumerate(prepared["requested_k_di"]):
            if index==1:
                rows.append({"status":"UNRESOLVED","spatial_helicity":"Right","k_di":point,
                             "reason":"manufactured missing-mode interval for protocol testing"})
            else:
                rows.append({"status":"RESOLVED","spatial_helicity":"Right","k_di":point,"requested_k_di":point,
                             "serialized_significant_digits":75,
                             "omega_real_radps":"1","omega_imaginary_raw_radps":"0","frequency_numerical_error_radps":"1e-30",
                             "phase_speed_relative_mhd_error":".005" if index==0 else ".02",
                             "group_speed_relative_mhd_error":".005" if index==0 else ".02",
                             "damping_over_real_frequency_raw":"0","damping_resolved":False,
                             "weak_damping_hermitian_energy_over_magnetic_energy":"2"})
        value={"schema":"srcsep-controlled-warm-parallel-study-v1","fixture_scope":"manufactured software protocol, not physics evidence",
               "input":self.case,"input_sha256":prepared["input"]["sha256"],
               "implementation_sha256":prepared["solver"]["sha256"],"implementation_archive":prepared["solver"]["archive"],
               "base_digits":40,"verification_digits":70,"requested_k_di":prepared["requested_k_di"],"requested_helicities":["Right"],
               "python_version":protocol.platform.python_version(),"mpmath_version":mp.__version__,
               "algorithm":"manufactured software-protocol fixture","root_isolation":"not physics evidence",
               "physical_production_qualified":False,"domain_limits_selected":False,"rows":rows}
        path=directory/"manufactured-study.json";protocol.atomic_json(path,value);return path,value

    def test_capture_records_actual_api_grid_and_archived_sources(self):
        grid=self.extracted["grid"]
        self.assertEqual(grid["bins"],4);self.assertEqual(len(grid["centers_hex"]),4);self.assertEqual(len(grid["edges_hex"]),5)
        for relative,record in self.extracted["source_files"].items():
            self.assertEqual(protocol.digest((self.directory/"grid"/record["archive"]).read_bytes()),record["sha256"],relative)
        self.assertTrue(self.extracted["compiler_version"].strip())

    def test_native_coordinate_retains_exact_binary_input_and_raw_spelling(self):
        grid=self.extracted["grid"]
        mapped=protocol.coordinates(grid,warm.Plasma(self.case),"native-grid",None,85)
        self.assertEqual(len(mapped["requested_k_di"]),6)
        with mp.workdps(120):
            numerator,denominator=float.fromhex(grid["centers_hex"][0]).as_integer_ratio()
            expected=mp.mpf(numerator)/denominator*warm.Plasma(self.case).scales()["ion_inertial_length_m"]
            self.assertLess(abs(mp.mpf(mapped["requested_k_di"][1])-expected),mp.mpf("1e-84")*expected)
        raw=[".00100","0.002000","3e-3"]
        explicit=protocol.coordinates(grid,warm.Plasma(self.case),"explicit-kdi",raw,85)
        self.assertEqual(explicit["requested_k_di"],raw)
        self.assertNotIn("requested_k_di",explicit["native_grid_map"][0])

    def test_changed_grid_point_or_missing_explicit_values_rejected(self):
        grid=self.extracted["grid"]
        points=protocol.coordinates(grid,warm.Plasma(self.case),"native-grid",None,85)["requested_k_di"]
        points[0]=".25"
        with self.assertRaises(ValueError):protocol.coordinates(grid,warm.Plasma(self.case),"native-grid",points,85)
        for bad in (None,[1,2],[],["nan"],["1","1"]):
            with self.assertRaises(ValueError):protocol.coordinates(grid,warm.Plasma(self.case),"explicit-kdi",bad,85)

    def test_import_preserves_gap_and_never_infers_contour_across_unknown(self):
        directory,prepared=self.make_protocol("good");path,_=self.make_study(directory,prepared)
        combined=protocol.combine_studies([path],prepared,directory)
        self.assertEqual(combined["coverage"]["Right"],{"RESOLVED":2,"UNRESOLVED":1})
        self.assertEqual(combined["reference_percentage_sampling_brackets"],[])
        self.assertEqual(combined["rows"][1]["status"],"UNRESOLVED")
        self.assertTrue((directory/combined["implementation_archive"]).is_file())

    def test_import_rejects_wrong_source_missing_rows_and_overlap(self):
        directory,prepared=self.make_protocol("bad");path,valid=self.make_study(directory,prepared)
        changed=copy.deepcopy(valid);changed["implementation_sha256"]="0"*64
        with self.assertRaises(ValueError):protocol.validate_study(changed,prepared)
        changed=copy.deepcopy(valid);changed["rows"].pop()
        with self.assertRaises(ValueError):protocol.validate_study(changed,prepared)
        with self.assertRaises(ValueError):protocol.combine_studies([path,path],prepared,directory)
        self.assertFalse((directory/"study.json").exists())

    def test_reported_coordinate_cannot_override_raw_request(self):
        directory,prepared=self.make_protocol("coordinate-tamper");path,value=self.make_study(directory,prepared)
        value["rows"][0]["k_di"]="123"
        with self.assertRaisesRegex(ValueError,"coordinate"):
            protocol.validate_study(value,prepared)
        protocol.atomic_json(path,value)
        with self.assertRaisesRegex(ValueError,"coordinate"):plotter.load_study(path)

    def test_mixed_runtime_parts_cannot_inherit_first_parts_version(self):
        directory,prepared=self.make_protocol("versions");prepared["requested_helicities"]=["Right","Left"]
        first,value=self.make_study(directory,prepared);second=copy.deepcopy(value)
        second["requested_helicities"]=["Left"];second["mpmath_version"]="manufactured different runtime"
        for row in second["rows"]:row["spatial_helicity"]="Left"
        protocol.atomic_json(directory/"left.json",second)
        with self.assertRaisesRegex(ValueError,"runtime"):
            protocol.combine_studies([first,directory/"left.json"],prepared,directory)

    def test_plot_labels_are_input_owned_and_unresolved_coverage_survives(self):
        directory,prepared=self.make_protocol("plot");path,value=self.make_study(directory,prepared)
        case=copy.deepcopy(self.case);case["magnetic_field_t"]="8e-9";case["proton_temperature_k"]="240000"
        caption=plotter.equilibrium_caption(case)
        self.assertIn("B=8",caption);self.assertIn("Tp=2.4000e+5",caption)
        metadata=plotter.render(path,directory/"figure")
        self.assertEqual(metadata["coverage"]["Right"]["UNRESOLVED"],1)
        self.assertEqual(metadata["physical_production_qualified"],False)
        for extension in ("png","pdf","plot.json"):self.assertTrue((directory/("figure."+extension)).is_file())
        missing=copy.deepcopy(value);missing["rows"].pop();protocol.atomic_json(directory/"missing.json",missing)
        with self.assertRaises(ValueError):plotter.load_study(directory/"missing.json")

    def test_archived_plotter_replays_without_repository_imports(self):
        directory,prepared=self.make_protocol("archive-replay");path,_=self.make_study(directory,prepared)
        archived=protocol.archive(PHYSICS/"plot_warm_parallel_study.py",directory)
        completed=subprocess.run([sys.executable,str(directory/archived["archive"]),"--study",str(path),
                                  "--output-prefix",str(directory/"replay")],cwd=str(self.directory),capture_output=True,text=True)
        self.assertEqual(completed.returncode,0,completed.stderr)
        self.assertEqual(json.loads((directory/"replay.plot.json").read_text())["plotter_sha256"],archived["sha256"])

    def test_positive_resolved_damping_cannot_underflow_silently(self):
        directory,prepared=self.make_protocol("underflow");path,value=self.make_study(directory,prepared)
        value["rows"][0]["damping_resolved"]=True;value["rows"][0]["damping_over_real_frequency_raw"]="1e-400"
        protocol.atomic_json(path,value)
        with self.assertRaisesRegex(ValueError,"float range"):plotter.render(path,directory/"underflow")
        self.assertFalse((directory/"underflow.png").exists())


if __name__=="__main__":unittest.main(verbosity=2)
