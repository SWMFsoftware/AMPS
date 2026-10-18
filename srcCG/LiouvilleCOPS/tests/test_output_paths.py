"""Output-root regression tests; real script/C++ checks are opt-in.

env COPS_RUN_INTEGRATION=1 python3 -m unittest discover -s tests -v
"""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))
import cops_paths


class PathTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cops_paths_")
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.app = self.base / "AMPS/srcCG/LiouvilleCOPS"
        self.app.mkdir(parents=True)
        self.root = self.base / "AMPS/cops_output"
        patch = mock.patch.object(cops_paths, "APP_ROOT", self.app)
        patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch.dict(os.environ)
        patch.start()
        self.addCleanup(patch.stop)
        for name in ("COPS_OUTPUT_ROOT", "COPS_DATA_ROOT", "COPS_L2_ROOT", "COPS_SOURCES_ROOT"):
            os.environ.pop(name, None)

    def test_default_root_and_no_creation_on_resolution(self):
        self.assertEqual(cops_paths.output_root(), self.root)
        self.assertEqual(cops_paths.output_path("plots/a.png"), self.root / "plots/a.png")
        self.assertFalse(self.root.exists())

    def test_override_and_nested_output_parent(self):
        override = self.base / "custom"
        os.environ["COPS_OUTPUT_ROOT"] = str(override)
        self.assertEqual(cops_paths.output_root(), override)
        out = cops_paths.prepare_output("data/sources/a.csv")
        self.assertEqual(Path(out), override / "data/sources/a.csv")
        self.assertTrue(Path(out).parent.is_dir())

    def test_external_data_defaults_and_overrides(self):
        expected = Path("/data/vtenishe/rosetta_data_end_of_mission_cops")
        self.assertEqual(cops_paths.data_root(), expected)
        self.assertEqual(cops_paths.cops_l2_root(), expected / "data_COPS_L2")
        self.assertEqual(cops_paths.sources_root(), expected / "sources")
        custom = self.base / "mission_data"
        os.environ["COPS_DATA_ROOT"] = str(custom)
        self.assertEqual(cops_paths.cops_l2_root(), custom / "data_COPS_L2")
        self.assertEqual(cops_paths.sources_root(), custom / "sources")

    def test_absolute_output_is_honored(self):
        target = self.base / "explicit/a.png"
        self.assertEqual(cops_paths.prepare_output(target), str(target))
        self.assertFalse(self.root.exists())

    def test_outputs_inside_application_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "must be outside"):
            cops_paths.output_path(self.app / "new-output")
        with self.assertRaisesRegex(ValueError, "must be outside"):
            cops_paths.output_path("../srcCG/LiouvilleCOPS/new-output")
        self.assertFalse((self.app / "new-output").exists())

    def test_product_prefers_new_over_legacy(self):
        legacy = self.app / "run/COPS_model.dat"
        legacy.parent.mkdir()
        legacy.write_text("old")
        self.assertEqual(cops_paths.product_input("run/COPS_model.dat"), str(legacy))
        new = self.root / "run/COPS_model.dat"
        new.parent.mkdir(parents=True)
        new.write_text("new")
        self.assertEqual(cops_paths.product_input("run/COPS_model.dat"), str(new))
        self.assertEqual(cops_paths.product_input(legacy), str(legacy))

    def test_explicit_dot_input_bypasses_product_search(self):
        self.assertEqual(cops_paths.product_input("./fit_base.in"), str(Path.cwd() / "fit_base.in"))

    def test_missing_input_points_to_new_root(self):
        self.assertEqual(cops_paths.product_input("missing-run"), str(self.root / "missing-run"))

    def test_help_does_not_create_output_directory(self):
        env = dict(os.environ, COPS_OUTPUT_ROOT=str(self.root))
        for script in ("cops_hk_table.py", "digitize_cops_figures.py", "collect_cops_sources.py",
                       "plot_cops_model.py", "plot_cops_model_altitude.py", "plot_lcurve.py",
                       "rebuild_cops_boundary.py", "fit_cops.py"):
            with self.subTest(script=script):
                result = subprocess.run([sys.executable, str(APP / "tools" / script), "--help"],
                                        env=env, text=True, capture_output=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertFalse(self.root.exists())

    def test_collector_preflights_missing_generated_inputs(self):
        env = dict(os.environ, COPS_OUTPUT_ROOT=str(self.root))
        result = subprocess.run(
            [sys.executable, str(APP / "tools/collect_cops_sources.py")],
            env=env, text=True, capture_output=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("required inputs are missing", result.stderr)
        self.assertIn("output_fit_base/COPS_model.dat", result.stderr)
        self.assertIn("output_M004/COPS_model.dat", result.stderr)
        self.assertIn("python3 tools/cops_hk_table.py", result.stderr)
        self.assertIn("./liouville_cops m004.in", result.stderr)
        self.assertFalse(self.root.exists())

    def test_hk_output_and_empty_input_guard(self):
        raw = self.base / "raw"
        raw.mkdir()
        (raw / "bg_test.tab").write_text(
            'START_TIME = "2016-09-30T00:00:00.000"\n'
            'STOP_TIME = "2016-09-30T00:00:04.000"\n'
            '"ROSINA_COPS_SENSOR_TEMP","","20","C"\n')
        env = dict(os.environ, COPS_OUTPUT_ROOT=str(self.root))
        command = [sys.executable, str(APP / "tools/cops_hk_table.py"), str(raw)]
        result = subprocess.run(command, env=env, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.root / "data/COPS_HK_20160930.txt"
        original = target.read_bytes()
        (raw / "bg_test.tab").unlink()
        result = subprocess.run(command, env=env, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(target.read_bytes(), original)


@unittest.skipUnless(os.environ.get("COPS_RUN_INTEGRATION") == "1", "opt-in script/solver integration")
class ProductIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cops_products_")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.out = self.work / "products"
        self.env = dict(os.environ, COPS_OUTPUT_ROOT=str(self.out), OMP_NUM_THREADS="2",
                        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        raw = self.work / "raw"
        raw.mkdir()
        record = (
            'START_TIME = "2016-09-30T00:00:00.000"\n'
            'STOP_TIME = "2016-09-30T00:00:04.000"\n'
            '"ROSINA_COPS_PRESSURE_NG","","1e-9","MBAR"\n'
            '"ROSINA_COPS_PRESSURE_RG","","2e-9","MBAR"\n'
            '"ROSINA_COPS_SENSOR_TEMP","","20","C"\n'
            '"ROSINA_COPS_ELECTR_TEMP","","20","C"\n'
            '"ROSINA_COPS_EMISSION_I_FILAMENT","","1","UA"\n'
            '"ROSINA_COPS_EMISSION_I_MICROTIP","","1","UA"\n'
            '"ROSINA_COPS_FIL_EMISSION_R_FLAG","HIGH","",""\n'
            '"ROSINA_COPS_MIC_EMISSION_R_FLAG","HIGH","",""\n'
            '"ROSINA_COPS_AUTORANGE_NG_FLAG","ON","",""\n'
            '"ROSINA_COPS_AUTORANGE_RG_FLAG","ON","",""\n')
        (raw / "bg_fixture_0.tab").write_text(record)
        (raw / "bg_fixture_1.tab").write_text(
            record.replace('START_TIME = "2016-09-30T00:00:00.000"',
                           'START_TIME = "2016-09-30T00:00:04.000"')
                  .replace('STOP_TIME = "2016-09-30T00:00:04.000"',
                           'STOP_TIME = "2016-09-30T00:01:04.000"'))
        self.run_command([sys.executable, str(APP / "tools/cops_hk_table.py"), str(raw)])

    def run_command(self, command, cwd=APP):
        result = subprocess.run(command, env=self.env, cwd=cwd, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def make_base(self, name="output_small"):
        source = self.work / (name + ".in")
        source.write_text((APP / "fit_base.in").read_text() +
                          "\nFirstPoint = 0\nLastPoint = 1\nOutputDirectory = " + name + "\n"
                          "NGAcceptance = full\nNGHeadOnWeight = off\nThreads = 2\n")
        self.run_command([str(APP / "liouville_cops"), str(source)])
        base = self.out / name
        manifest = json.loads((base / "run_manifest.json").read_text())
        self.assertEqual(manifest["schema"], "liouville_cops.run_manifest")
        self.assertEqual(manifest["states"]["process"], "complete")
        self.assertEqual(manifest["states"]["cache"], "verified")
        self.assertTrue((base / "input.snapshot.in").is_file())
        return base

    def test_plot_digitize_collect_into_output_root(self):
        source_files = [p for p in cops_paths.sources_root().rglob("*") if p.is_file()]
        old = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
        base = self.make_base()
        self.run_command([sys.executable, str(APP / "tools/digitize_cops_figures.py"), "--check"], cwd=self.work)
        self.assertEqual(len(list((self.out / "data/sources").glob("*.csv"))), 5)
        self.run_command([sys.executable, str(APP / "tools/collect_cops_sources.py"),
                          "--shap5-run", "output_small", "--m004-run", "output_small",
                          "--fit-run", "output_small"], cwd=self.work)
        for name in ("COPS_EoM_sources.json", "COPS_EoM_sources_time.png", "COPS_EoM_sources_altitude.png",
                     "COPS_fit_data.txt", "COPS_fit_data.png"):
            self.assertTrue((self.out / "data" / name).is_file(), name)
        scan = self.out / "scan"
        scan.mkdir()
        (scan / "lcurve.txt").write_text(
            "# maximum curvature point: lambda = 1e-2\n"
            "#lambda rmsNG rmsRG rho D Tmax T75 Tmin kappa activity_ratio activity_ratio_seen unseen_activity_fraction status\n"
            "1e-3 0.1 0.2 0.3 0.4 180 140 130 1 1 1 0.2 1\n"
            "1e-2 0.2 0.3 0.4 0.5 181 141 131 1.1 1.1 1.1 0.3 0\n")
        for script, source, filename in (
                ("plot_cops_model.py", "output_small", "model.png"),
                ("plot_cops_model_altitude.py", "output_small", "altitude.png"),
                ("plot_lcurve.py", "scan", "lcurve.png")):
            self.run_command([sys.executable, str(APP / "tools" / script), source,
                              "--save", "plots/" + filename], cwd=self.work)
            self.assertTrue((self.out / "plots" / filename).is_file())
        self.assertEqual(old, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files})

    def test_cpp_fit_and_generated_verification_paths(self):
        base = self.make_base()
        self.assertTrue((base / "COPS_model.dat").is_file())
        self.run_command([sys.executable, str(APP / "tools/fit_cops.py"), "output_small",
                          "--data-only", "--out", "selection_small"], cwd=self.work)
        self.assertTrue((self.out / "selection_small/fit_data_used.txt").is_file())
        selection_manifest = json.loads((self.out / "selection_small/fit_manifest.json").read_text())
        self.assertEqual(selection_manifest["states"]["process"], "complete")
        self.assertEqual(selection_manifest["states"]["optimizer"], "not_requested")
        self.assertEqual(selection_manifest["states"]["cache"], "verified")
        self.run_command([sys.executable, str(APP / "tools/fit_cops.py"), "output_small",
                          "--max-nfev", "1", "--lambda", "0.001"], cwd=self.work)
        fitted = self.out / "output_small_fit"
        self.assertTrue((fitted / "fit_result.txt").is_file())
        fit_manifest = json.loads((fitted / "fit_manifest.json").read_text())
        self.assertEqual(fit_manifest["states"]["process"], "complete")
        self.assertEqual(fit_manifest["states"]["cache"], "verified")
        self.assertEqual(fit_manifest["states"]["cpp_prediction_agreement"], "not_run")
        lines = (fitted / "COPS_model.dat").read_text().splitlines()
        names = next(line[1:].split() for line in lines if line.startswith("#iPoint "))
        row = next(line.split() for line in lines if line and not line.startswith("#"))
        values = dict(zip(names, row))
        for name in ("Q_H2O", "Q_CO2", "n_H2O", "n_CO2", "PhiRG_H2O",
                     "PhiRG_CO2", "pNG_H2O", "pNG_CO2", "pRG_H2O", "pRG_CO2",
                     "NGBelowDeckFraction"):
            self.assertEqual(values[name], "nan", name)
        self.run_command([str(APP / "liouville_cops"), "output_small_fit/cops_fit.in"])
        verified = self.out / "output_small_fit_verify"
        self.assertTrue((verified / "COPS_model.dat").is_file())
        self.run_command([sys.executable, str(APP / "tools/validate_cops_run.py"),
                          "--fit", "output_small_fit", "--verification",
                          "output_small_fit_verify"], cwd=self.work)
        fit_manifest = json.loads((fitted / "fit_manifest.json").read_text())
        self.assertEqual(fit_manifest["states"]["cpp_prediction_agreement"], "verified")

    def test_geometry_only_matches_full_geometry(self):
        full = self.make_base("output_geometry_full")
        source = self.work / "geometry_only.in"
        source.write_text((APP / "fit_base.in").read_text() +
                          "\nFirstPoint = 0\nLastPoint = 1\nOutputDirectory = output_geometry_only\n"
                          "WriteKernels = off\nGeometryOnly = on\nThreads = 2\n")
        self.run_command([str(APP / "liouville_cops"), str(source)])
        geometry = self.out / "output_geometry_only"
        self.assertFalse((geometry / "kernels").exists())
        self.assertFalse((geometry / "elements.bin").exists())

        def table(path):
            lines = path.read_text().splitlines()
            columns = next(line[1:].split() for line in lines if line.startswith("#iPoint "))
            rows = [dict(zip(columns, line.split())) for line in lines if line and not line.startswith("#")]
            return lines, rows

        full_lines, full_rows = table(full / "COPS_model.dat")
        geometry_lines, geometry_rows = table(geometry / "COPS_model.dat")
        self.assertTrue(any("GeometryOnly=off" in line for line in full_lines))
        self.assertTrue(any("GeometryOnly=on" in line for line in geometry_lines))
        self.assertEqual(len(geometry_rows), 2)
        for actual, expected in zip(geometry_rows, full_rows):
            for key in ("iPoint", "UTC", "x", "y", "z", "altitude", "inside", "rh_AU"):
                self.assertEqual(actual[key], expected[key], key)
            self.assertEqual(float(actual["pNG"]), 0.0)
            self.assertEqual(float(actual["pRG"]), 0.0)

    def test_cpp_refuses_to_overwrite_manifested_run(self):
        base = self.make_base("output_no_overwrite")
        protected = {name: hashlib.sha256((base / name).read_bytes()).hexdigest()
                     for name in ("COPS_model.dat", "elements.bin", "run_manifest.json")}
        result = subprocess.run([str(APP / "liouville_cops"),
                                 str(base / "input.snapshot.in")],
                                env=self.env, cwd=APP, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("output directory already exists", result.stdout)
        self.assertEqual(protected, {
            name: hashlib.sha256((base / name).read_bytes()).hexdigest()
            for name in protected})


if __name__ == "__main__":
    unittest.main()
