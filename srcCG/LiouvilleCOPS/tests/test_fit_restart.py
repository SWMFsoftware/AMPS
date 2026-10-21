"""CPU restart round trips and optional interoperability with the untouched CUDA reader.

No GPU, SPICE kernels, or mission run is required. Synthetic pressure values
come from the production CPU evaluator. The sibling is imported read-only
only for the interoperability test; CPU runtime has no sibling dependency.
"""
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
import subprocess
import unittest
from unittest.mock import patch

import numpy as np
from test_joint_template_fit import joint_problem
from test_a4_enhancement_fit import synthetic_problem
import fit_cops
import fit_restart as rs
from cops_io import file_identity
import run_fit_profile


def fixture(root, degree=None):
    p, el, _ = joint_problem() if degree is None else synthetic_problem(degree)
    p.args.no_save_restart = False
    x = p.x0.copy()
    if p.na:
        x[0] += .1
    if p.nb:
        x[p.bslice] = np.linspace(.01, .05, p.nb)
    a, b, t, k = p.unpack(x)
    ng, rg = p.M.evaluate(a, t, k, b=b, enhancement_basis=p.enhancement_basis,
                          enhancement_scale=p.enhancement_scale)
    root.mkdir(parents=True, exist_ok=True)
    lines = ["#iPoint UTC altitude refNG refRG useNG useRG pNG pRG"]
    for i in range(p.M.npt):
        lines.append("%d 2016-09-30T10:00:%02d 100 %.17g %.17g 1 1 %.17g %.17g" %
                     (i, i, p.M.refNG[i], p.M.refRG[i], ng[i], rg[i]))
    (root / "COPS_model.dat").write_text("\n".join(lines)+"\n")
    manifest = {"inputs": {}, "selection": {}, "restart": {"save": True},
                "states": {"process": "complete"}}
    (root / "fit_manifest.json").write_text(json.dumps(manifest))
    result = dict(x=x, lam=0., success=False, status=0, nfev=2,
                  optimality=1., cost=2., message="synthetic evaluation limit")
    state = rs.from_fit(p, result, root, manifest, el)
    return p, el, result, state, manifest


class RestartTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cops_restart_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_roundtrip_and_exact_start_for_all_modes(self):
        for degree in (None, 4, 5, 8):
            p, el, result, state, _ = fixture(self.root/str(degree), degree)
            path = self.root/(str(degree)+".npz")
            prior = p.x0.copy()
            rs.write(path, state)
            read = rs.read(path)
            for key in rs.ARRAY_KEYS:
                np.testing.assert_array_equal(read.arrays[key], state.arrays[key])
            start, report = rs.start_vector(p, read)
            self.assertEqual(start.tobytes(), result["x"].tobytes())
            np.testing.assert_array_equal(prior, p.x0)
            self.assertEqual(report["from_restart"], p.names)
            rs.check_mesh(read, el)
            with self.assertRaises(FileExistsError):
                rs.write(path, state)

    def test_corruption_mesh_and_nonfinite_rejected(self):
        _, el, _, state, _ = fixture(self.root)
        path = self.root/"state.npz"
        rs.write(path, state)
        with np.load(path, allow_pickle=False) as data:
            payload = {k: np.array(data[k]) for k in data.files}
        payload["a"][0] += 1
        np.savez_compressed(self.root/"bad.npz", **payload)
        with self.assertRaisesRegex(ValueError, "corrupt"):
            rs.read(self.root/"bad.npz")
        bad = deepcopy(el)
        bad["centroid"][0, 0] += 1
        with self.assertRaisesRegex(ValueError, "mesh differs"):
            rs.check_mesh(state, bad)
        state.a[0] = np.nan
        rs.write(self.root/"nan.npz", state)
        with self.assertRaisesRegex(ValueError, "physical array"):
            rs.read(self.root/"nan.npz")

    def test_ordered_temperature_missing_names_bounds_and_scale(self):
        p, el, _, state, _ = fixture(self.root)
        args = deepcopy(p.args)
        args.ordered_temperature = True
        ordered = fit_cops.Problem(p.M, el, p.coef, args)
        x, _ = rs.start_vector(ordered, state)
        np.testing.assert_array_equal(ordered.unpack(x)[2], state.T)
        # Missing template names keep their new problem start; renamed fields
        # do not accidentally acquire an amplitude by column position.
        extra = deepcopy(p)
        extra.names[extra.na] = "c_new_template"
        x, report = rs.start_vector(extra, state)
        self.assertEqual(x[extra.na], extra.x0[extra.na])
        self.assertIn("c_new_template", report["kept_problem_start"])
        p.hi[p.names.index("Tmax")] = 180
        x, report = rs.start_vector(p, state)
        self.assertEqual(x[p.names.index("Tmax")], 180)
        self.assertEqual(report["clipped"][0]["name"], "Tmax")
        p.enhancement_scale *= 2
        x, report = rs.start_vector(p, state)
        np.testing.assert_array_equal(x[p.bslice], state.template_values / 2)
        self.assertEqual(report["enhancement_scale_ratio"], .5)

    def test_template_sources_are_verified_and_explicit_arguments_win(self):
        _, _, _, state, _ = fixture(self.root)
        bundle = self.root/"templates.npz"
        bundle.write_bytes(b"synthetic bundle identity only")
        state.meta["templates"]["sources"] = [
            dict(name=n, bundle=file_identity(bundle)) for n in state.template_names]
        args, names = rs.restore_template_arguments(state, [])
        self.assertEqual(names, state.template_names)
        self.assertEqual(args[0], [str(bundle), state.template_names[0]])
        bundle.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "bundle changed"):
            rs.restore_template_arguments(state, [])
        self.assertEqual(rs.restore_template_arguments(state, [["new", "other"]]),
                         ([["new", "other"]], False))

    def test_save_each_lambda_history_optout_and_external_owner(self):
        p, el, result, previous, manifest = fixture(self.root/"base")
        p.restart_state = previous
        manifest["restart"]["loaded"] = {"path": "earlier.npz", "sha256": "fixture"}
        live = self.root/"fit_manifest.json"
        with patch.multiple(fit_cops, P=p, ACTIVE_MANIFEST=manifest,
                            ACTIVE_MANIFEST_PATH=str(live)):
            for lam in (0., .1):
                dest = self.root/("lambda=%.3e" % lam)
                fixture(dest)
                result["lam"] = lam
                products = fit_cops.save_fit_restart(dest, result, el)
                self.assertIn("fit_restart.npz", products)
                state = rs.read(dest/"fit_restart.npz")
                self.assertEqual(state.meta["settings"]["lambda"], lam)
                self.assertEqual(len(state.meta["history"]), 2)
                self.assertEqual(state.meta["fit_manifest"],
                                 file_identity(dest/"fit_restart_manifest.json"))
            p.args.no_save_restart = True
            self.assertEqual(fit_cops.save_fit_restart(self.root/"off", result, el), {})
            self.assertFalse((self.root/"off").exists())
            p.args.no_save_restart = False
            manifest["compute_backend"] = {"restart": {"save": True}}
            self.assertEqual(fit_cops.save_fit_restart(self.root/"cuda", result, el), {})
            self.assertFalse((self.root/"cuda").exists())

    def test_fit_one_passes_restart_without_changing_prior(self):
        p, _, _, _, _ = fixture(self.root)
        p.args.max_nfev, p.args.jobs = 1, 1
        p.args.ftol = p.args.xtol = p.args.gtol = 1e-10
        p.restart_start = p.x0.copy()
        p.restart_start[0] += .2
        prior = p.x0.copy()
        class StopAtOptimizer(Exception):
            pass
        def inspect(fun, start, **kwargs):
            np.testing.assert_array_equal(start, p.restart_start)
            raise StopAtOptimizer
        with patch.multiple(fit_cops, P=p, least_squares=inspect):
            with self.assertRaises(StopAtOptimizer):
                fit_cops.fit_one(.1)
        np.testing.assert_array_equal(p.x0, prior)

    def test_profile_runner_forwards_options(self):
        args = run_fit_profile.build_parser().parse_args(
            ["--profile", "all_available", "--restart-from", "/tmp/old.npz", "--no-save-restart"])
        cmd = run_fit_profile.fit_command(args, "all_available", "new")
        self.assertIn("--no-save-restart", cmd)
        self.assertEqual(cmd[cmd.index("--restart-from")+1], "/tmp/old.npz")

    def test_cuda_reader_writer_interoperability_read_only(self):
        sibling = Path(fit_cops.__file__).resolve().parents[2]/"LiouvilleCOPS_cuda"
        if not (sibling/"cops_cuda/restart.py").is_file():
            self.skipTest("optional read-only CUDA sibling absent")
        sys.path.insert(0, str(sibling))
        self.addCleanup(lambda: sys.path.remove(str(sibling)))
        from cops_cuda import restart as cuda
        for degree in (None, 4, 5, 8):
            p, el, result, state, manifest = fixture(self.root/str(degree), degree)
            cpu_path, gpu_path = self.root/(str(degree)+"cpu.npz"), self.root/(str(degree)+"gpu.npz")
            rs.write(cpu_path, state)
            gpu_read = cuda.read(cpu_path)
            np.testing.assert_array_equal(cuda.start_vector(p, gpu_read)[0], result["x"])
            gpu_state = cuda.from_fit(p, result, self.root/str(degree), manifest, el)
            cuda.write(gpu_path, gpu_state)
            cpu_read = rs.read(gpu_path)
            np.testing.assert_array_equal(rs.start_vector(p, cpu_read)[0], result["x"])
            for key in rs.ARRAY_KEYS:
                np.testing.assert_array_equal(cpu_read.arrays[key], state.arrays[key])

    def test_cli_fit_resume_scan_and_optout(self):
        """Run real entry points on a tiny binary base, not mocked fit writers.

        Reuse the sibling's read-only synthetic binary-format generator when
        available. All generated inputs/outputs live in this test's /tmp root.
        """
        app = Path(fit_cops.__file__).resolve().parents[1]
        sibling_tests = app.parent/"LiouvilleCOPS_cuda/tests"
        if not (sibling_tests/"synthetic_run.py").is_file():
            self.skipTest("optional CUDA synthetic binary fixture absent")
        import importlib.util
        spec = importlib.util.spec_from_file_location("restart_synthetic_run", sibling_tests/"synthetic_run.py")
        synthetic = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(synthetic)
        base = synthetic.make_run(self.root/"run", nfaces=80, empty_point=-1)
        run = base["run"]
        header, names, rows = fit_cops.read_model_file(run/"COPS_model.dat")
        cols = {n: i for i, n in enumerate(names)}
        coef = fit_cops.read_activity_coefficients(str(app), header)
        table = base["elements"]
        table[:, 32] = np.maximum(table[:, 7:32]@coef[0], 0)
        table[:, 33] = np.maximum(table[:, 7:32]@coef[1], 0)
        synthetic.write_elements(run/"elements.bin", table)
        el = fit_cops.read_elements(run/"elements.bin")
        model = fit_cops.Model(str(run), header, cols, rows, el, coef)
        ng, rg = model.evaluate(coef[0], model.T0, 1.)
        for i, row in enumerate(rows):
            for gauge, pressure in (("NG", ng[i]), ("RG", rg[i])):
                row[cols["p"+gauge]] = "%.17g" % pressure
                row[cols["ref"+gauge]] = "%.17g" % (pressure*1.03)
        (run/"COPS_model.dat").write_text(
            "".join("# "+h+"\n" for h in header)+"#"+" ".join(names)+"\n"+
            "".join(" ".join(row)+"\n" for row in rows))
        (run/"input.snapshot.in").write_text("ActivityModel = shap5\n")
        synthetic.write_manifest(run, base["point_ids"])
        manifest_path = run/"run_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["effective_settings"] = {"BodyFrame": ["67P/C-G_CK"]}
        manifest_path.write_text(json.dumps(manifest))
        hk = self.root/"hk.txt"
        hk.write_text("# utc_mid pNG_mbar emisRangeFIL stop\n"
                      "2016-09-30T00:00:00 1e-9 LOW 2016-09-30T00:00:00\n")
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
                   OMP_NUM_THREADS="1", MPLCONFIGDIR=str(self.root/"mpl"))
        common = [sys.executable, str(app/"tools/fit_cops.py"), str(run),
                  "--degree", "0", "--fix-temperature", "--fix-kappa",
                  "--max-nfev", "2", "--hk-table", str(hk)]
        def launch(name, extra):
            dest = self.root/name
            proc = subprocess.run(common+["--out", str(dest)]+extra, env=env,
                                  text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stdout)
            return dest
        first = launch("first", [])
        saved = rs.read(first/"fit_restart.npz")
        scan = launch("scan", ["--restart-from", str(first/"fit_restart.npz"),
                               "--lambda", "0", ".1", "--jobs", "2"])
        scan_manifest = json.loads((scan/"fit_manifest.json").read_text())
        np.testing.assert_array_equal(scan_manifest["restart"]["x_start_scaled"], saved.x)
        for name in ("lambda=0.000e+00", "lambda=1.000e-01"):
            state = rs.read(scan/name/"fit_restart.npz")
            self.assertEqual(len(state.meta["history"]), 2)
            self.assertIn("fit_restart.npz", scan_manifest["outputs"]["lambda_directories"][name])
        off = launch("off", ["--restart-from", str(first/"fit_restart.npz"), "--no-save-restart"])
        self.assertFalse((off/"fit_restart.npz").exists())
        diagnostic = launch("data", ["--data-only"])
        self.assertFalse((diagnostic/"fit_restart.npz").exists())


if __name__ == "__main__":
    unittest.main()
