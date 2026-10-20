"""Standard-library tests; no CSPICE, mission data, or scientific runs required.

Run: python3 -m unittest discover -s tests -p 'test_rebuild_cops_boundary.py' -v
"""

import contextlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "rebuild_cops_boundary.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("rebuild_cops_boundary", SCRIPT)
rebuild = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rebuild)

INPUT = """# Test input
Mode = observations
ShapeModel = mesh.bdf
ObservationFile = observations.h
KernelDirectory = spice
Kernel = first.tls
Kernel = second.tf
TemperatureRGMode = hk_table
TemperatureRGTable = hk.txt
TemperatureRGTableColumn = 6
TemperatureModel = analytical
AnalyticalTmax = 182.09922
AnalyticalT75 = 136.11149
AnalyticalTmin = 133
SurfaceTemperatureCorrectionFactor = 1
FlightTimeRotation = none
SurfaceBC = rough
SurfaceBC = flat # last scalar wins in C++
WriteKernels = off
OutputDirectory = old
Threads = 8
FirstPoint = 0
LastPoint = 1
PointStep = 1
ExcludeNG = start1 stop1
ExcludeNG = start2 stop2
"""

FAKE_SOLVER = r'''
import pathlib, struct, sys
values = {}
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    active = line.split("#", 1)[0].strip()
    if active:
        key, value = active.split("=", 1)
        values[key.strip()] = value.strip()
out = pathlib.Path(values["OutputDirectory"])
(out / "kernels").mkdir()
table = "# FlightTimeRotation=none\n# SurfaceBC=" + values["SurfaceBC"] + "\n#iPoint UTC\n"
for ip in range(int(values["FirstPoint"]), min(int(values["LastPoint"]), 1) + 1, int(values["PointStep"])):
    table += str(ip) + " fake\n"
    header = struct.pack("<ii10di", 1, ip, *([0.] * 10), 1)
    record = struct.pack("<i4di", 0, 1., 1., 1., 1., 0)
    (out / "kernels" / ("point=%04d.bin" % ip)).write_bytes(header + record)
(out / "COPS_model.dat").write_text(table)
(out / "elements.bin").write_bytes(struct.pack("<iii", 1, 1, 25) + bytes(272))
print("Observations 'fake': 2 points, first .. last")
print("solver stderr is logged", file=sys.stderr)
'''


class RebuildTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cops_boundary_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        patch = mock.patch.object(rebuild, "ROOT", self.root)
        patch.start()
        self.addCleanup(patch.stop)
        (self.root / "fit_base.in").write_text(INPUT)
        (self.root / "spice").mkdir()
        for name in ("mesh.bdf", "observations.h", "hk.txt", "spice/first.tls", "spice/second.tf"):
            (self.root / name).write_text("fixture")
        self.executable = self.root / "liouville_cops"
        self.executable.write_text("#!" + sys.executable + "\n" + FAKE_SOLVER)
        self.executable.chmod(0o755)
        self.out = self.root / "new_rough"

    def run_main(self, extra=()):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = rebuild.main(["--surface-bc", "rough", "--out", str(self.out)] + list(extra))
        return status, stdout.getvalue(), stderr.getvalue()

    def test_preserves_temperature_intervals_and_normalizes_paths(self):
        text, values = rebuild.prepare_input(INPUT, "rough", self.out, 2)
        self.assertEqual(values["SurfaceBC"], ["rough"])
        self.assertEqual(values["WriteKernels"], ["on"])
        self.assertEqual(values["Threads"], ["2"])
        self.assertEqual(values["AnalyticalTmax"], ["182.09922"])
        self.assertEqual(values["ExcludeNG"], ["start1 stop1", "start2 stop2"])
        self.assertEqual(values["Kernel"], [str(self.root / "spice/first.tls"), str(self.root / "spice/second.tf")])
        self.assertEqual(values["ObservationFile"], [str(self.root / "observations.h")])
        self.assertIn("# Test input", text)

    def test_preserves_threads_without_override(self):
        _, values = rebuild.prepare_input(INPUT, "flat", self.out, None)
        self.assertEqual(values["Threads"], ["8"])

    def test_shape_model_and_format_are_explicit_overrides(self):
        """A resolution change must reach the archived solver input exactly."""
        selected = self.root / "shape_models" / "shap5_199k.bdf"
        text, values = rebuild.prepare_input(
            INPUT, "flat", self.out, 4, selected, "nastran")
        self.assertEqual(values["ShapeModel"], [str(selected)])
        self.assertEqual(values["ShapeModelFormat"], ["nastran"])
        self.assertEqual(text.count("ShapeModel = " + str(selected)), 1)
        self.assertIn("# Replaced by rebuild tool: ShapeModel = mesh.bdf", text)
        without_shape = "\n".join(
            line for line in INPUT.splitlines()
            if not line.startswith("ShapeModel =")) + "\n"
        _, replacement = rebuild.prepare_input(
            without_shape, "flat", self.out, 4, selected, "nastran")
        self.assertEqual(replacement["ShapeModel"], [str(selected)])

    def test_shape_format_validation_is_not_inferred_from_filename(self):
        """Do not guess BDF/DSK semantics from an extension or model name."""
        with self.assertRaisesRegex(ValueError, "nastran or dsk"):
            rebuild.prepare_input(INPUT, "flat", self.out, None,
                                  self.root / "mesh.any", "automatic")

    def test_rejects_incompatible_models_and_invalid_ranges(self):
        for key, value in (("Mode", "selftest"), ("FlightTimeRotation", "exact"),
                           ("TemperatureModel", "constant"),
                           ("SurfaceTemperatureCorrectionFactor", "2"),
                           ("PointStep", "0"), ("FirstPoint", "-1"),
                           ("LastPoint", "-1"), ("Threads", "0")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                rebuild.prepare_input(INPUT + key + " = " + value + "\n", "rough", self.out, None)

    def test_rejects_nonliteral_paths(self):
        for path in ("~/shape", "$DATA/shape", '"shape"', "path with space"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                rebuild.prepare_input(INPUT + "ShapeModel = " + path + "\n", "rough", self.out, None)

    def test_dry_run_creates_nothing_and_does_not_require_executable(self):
        self.executable.unlink()
        status, output, _ = self.run_main(["--dry-run"])
        self.assertEqual(status, 0)
        self.assertIn("SurfaceBC = rough", output)
        self.assertFalse(self.out.exists())

    def test_existing_directory_and_dangling_symlink_are_refused(self):
        self.out.mkdir()
        marker = self.out / "keep"
        marker.write_text("unchanged")
        self.assertEqual(self.run_main()[0], 1)
        self.assertEqual(marker.read_text(), "unchanged")
        link = self.root / "dangling"
        link.symlink_to(self.root / "absent")
        self.out = link
        self.assertEqual(self.run_main()[0], 1)
        self.assertTrue(link.is_symlink())

    def test_missing_external_file_is_rejected_before_creation(self):
        (self.root / "hk.txt").unlink()
        self.assertEqual(self.run_main()[0], 1)
        self.assertFalse(self.out.exists())

    def test_missing_shape_override_is_rejected_before_creation(self):
        status, _, error = self.run_main([
            "--shape-model", "absent.bdf", "--shape-model-format", "nastran"])
        self.assertEqual(status, 1)
        self.assertIn("Required data file is missing", error)
        self.assertFalse(self.out.exists())

    def test_rejects_nonpositive_thread_option(self):
        self.assertEqual(self.run_main(["--threads", "0"])[0], 1)
        self.assertFalse(self.out.exists())

    def test_end_to_end_with_fixture_solver(self):
        original = (self.root / "fit_base.in").read_bytes()
        status, _, error = self.run_main(["--threads", "2"])
        self.assertEqual(status, 0, error)
        metadata = json.loads((self.out / "rebuild.json").read_text())
        self.assertEqual(metadata["status"], "complete")
        self.assertEqual(metadata["verified_points"], 2)
        self.assertEqual(metadata["cwd"], str(self.root))
        self.assertEqual((self.root / "fit_base.in").read_bytes(), original)
        self.assertEqual((self.out / "source.in").read_bytes(), original)
        self.assertIn("solver stderr is logged", (self.out / "rebuild.log").read_text())

    def test_end_to_end_records_selected_shape(self):
        selected = self.root / "selected.bdf"
        selected.write_text("fixture")
        status, _, error = self.run_main([
            "--shape-model", str(selected), "--shape-model-format", "nastran"])
        self.assertEqual(status, 0, error)
        metadata = json.loads((self.out / "rebuild.json").read_text())
        self.assertEqual(metadata["shape_model"], str(selected))
        self.assertEqual(metadata["shape_model_format"], "nastran")
        self.assertTrue(metadata["shape_model_overridden"])
        generated = rebuild.parse_input((self.out / "boundary.in").read_text())
        self.assertEqual(generated["ShapeModel"], [str(selected)])

    def test_solver_failure_retains_failed_manifest(self):
        self.executable.write_text("#!" + sys.executable + "\nimport sys\nsys.exit(7)\n")
        status, _, _ = self.run_main()
        self.assertEqual(status, 1)
        metadata = json.loads((self.out / "rebuild.json").read_text())
        self.assertEqual(metadata["status"], "failed")
        self.assertEqual(metadata["solver_returncode"], 7)
        self.assertTrue((self.out / "boundary.in").exists())

    def test_failed_launch_retains_failed_manifest(self):
        self.executable.write_text("not an executable format")
        self.assertEqual(self.run_main()[0], 1)
        metadata = json.loads((self.out / "rebuild.json").read_text())
        self.assertEqual(metadata["status"], "failed")

    def test_validation_detects_truncated_or_wrong_version_kernel(self):
        self.assertEqual(self.run_main()[0], 0)
        values = rebuild.parse_input((self.out / "boundary.in").read_text())
        kernel = self.out / "kernels/point=0000.bin"
        data = kernel.read_bytes()
        for bad in (data[:20], data[:-1], struct.pack("<i", 2) + data[4:]):
            kernel.write_bytes(bad)
            with self.assertRaises(ValueError):
                rebuild.verify_outputs(self.out, values, "rough")

    def test_validation_detects_incomplete_table(self):
        self.assertEqual(self.run_main()[0], 0)
        values = rebuild.parse_input((self.out / "boundary.in").read_text())
        table = self.out / "COPS_model.dat"
        table.write_text(table.read_text().replace("1 fake\n", ""))
        with self.assertRaises(ValueError):
            rebuild.verify_outputs(self.out, values, "rough")


@unittest.skipUnless(os.environ.get("COPS_RUN_INTEGRATION") == "1",
                     "set COPS_RUN_INTEGRATION=1 to use the local solver and mission data")
class RealSolverTests(unittest.TestCase):
    def test_boundary_and_surface_temperature_dependencies(self):
        """Three one-point runs, isolated in /tmp; never touch retained outputs."""
        root = SCRIPT.parents[1]
        template = (root / "fit_base.in").read_text()
        template = template.replace("TemperatureRGMode = hk_table", "TemperatureRGMode = constant")
        template = "\n".join(
            line for line in template.splitlines()
            if not line.strip().startswith(("TemperatureRGTable =", "TemperatureRGTableColumn ="))) + "\n"
        # Turn off velocity correction to isolate the T-independent RG flux.
        template += ("\nFirstPoint = 0\nLastPoint = 0\nPointStep = 1\n"
                     "SpacecraftVelocity = off\nNGAcceptance = full\nNGHeadOnWeight = off\n"
                     "AnalyticalTmax = 182.09922\nAnalyticalT75 = 136.11149\nAnalyticalTmin = 133\n")
        with tempfile.TemporaryDirectory(prefix="cops_boundary_physics_") as temp:
            work = Path(temp)
            base = work / "base.in"
            base.write_text(template)
            warm = work / "warm.in"
            warm.write_text(template + "\nAnalyticalTmax = 364.19844\nAnalyticalT75 = 272.22298\nAnalyticalTmin = 266\n")
            outputs = {}
            for name, bc, infile in (("flat", "flat", base), ("rough", "rough", base),
                                      ("warm", "flat", warm)):
                out = work / name
                result = subprocess.run(
                    [sys.executable, str(SCRIPT), "--input", str(infile), "--surface-bc", bc,
                     "--out", str(out), "--threads", "2"],
                    cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    timeout=180)
                self.assertEqual(result.returncode, 0, result.stdout)
                outputs[name] = out
            kernel = "kernels/point=0000.bin"
            self.assertNotEqual((outputs["flat"] / kernel).read_bytes(),
                                (outputs["rough"] / kernel).read_bytes())
            self.assertEqual((outputs["flat"] / kernel).read_bytes(),
                             (outputs["warm"] / kernel).read_bytes())
            elements = (outputs["flat"] / "elements.bin").read_bytes()
            for name in ("rough", "warm"):
                self.assertEqual(elements, (outputs[name] / "elements.bin").read_bytes())
            rows = {}
            for name, out in outputs.items():
                lines = (out / "COPS_model.dat").read_text().splitlines()
                columns = next(line[1:].split() for line in lines if line.startswith("#iPoint "))
                data = next(line.split() for line in lines if line and not line.startswith("#"))
                rows[name] = dict(zip(columns, data))
            for key in ("pNG", "n_H2O", "n_CO2"):
                ratio = float(rows["warm"][key]) / float(rows["flat"][key])
                self.assertTrue(math.isclose(ratio, 1 / math.sqrt(2), rel_tol=1e-8), (key, ratio))
            for key in ("pRG", "PhiRG_H2O", "PhiRG_CO2", "Q_H2O", "Q_CO2", "T_RG"):
                self.assertTrue(math.isclose(float(rows["warm"][key]), float(rows["flat"][key]), rel_tol=1e-8), key)
            for key in ("Q_H2O", "Q_CO2"):
                self.assertEqual(rows["rough"][key], rows["flat"][key])


if __name__ == "__main__":
    unittest.main()
