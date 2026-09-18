#!/usr/bin/env python3
"""Dependency-light contract tests for the D03 native release gate."""

from __future__ import annotations

import argparse
import copy
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "test" / "run_native_integration.py"
EXAMPLE = ROOT / "test" / "native_integration_manifest.example.json"


def _load_runner():
    specification = importlib.util.spec_from_file_location(
        "srcsep_d03_native_runner", RUNNER)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load D03 native runner")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class NativeIntegrationContractTests(unittest.TestCase):
    def setUp(self):
        self.runner = _load_runner()
        self.manifest = json.loads(EXAMPLE.read_text(encoding="utf-8"))

    def test_example_has_complete_mover_mpi_refresh_restart_coverage(self):
        self.assertEqual(self.runner.validate_manifest(self.manifest), [])

    def test_missing_second_mpi_decomposition_is_rejected(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["cases"] = [
            case for case in manifest["cases"]
            if case["id"] != "fte-mfp-mpi4"]
        errors = self.runner.validate_manifest(manifest)
        self.assertTrue(any("fte-mfp needs rank 1" in error for error in errors),
                        errors)

    def test_restart_needs_shared_exact_artifact_group(self):
        manifest = copy.deepcopy(self.manifest)
        resume = next(case for case in manifest["cases"]
                      if case["kind"] == "resume")
        resume["artifacts"][0]["equivalence_group"] = "different-final"
        errors = self.runner.validate_manifest(manifest)
        self.assertTrue(any("shared artifact equivalence_group" in error
                            for error in errors), errors)

    def test_every_case_must_prove_particle_dispatch(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["cases"][0]["expect"].pop("minimum_particle_dispatches")
        errors = self.runner.validate_manifest(manifest)
        self.assertTrue(any("minimum_particle_dispatches >= 1" in error
                            for error in errors), errors)

    def test_malformed_refresh_threshold_is_reported_not_raised(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["cases"][0]["expect"]["minimum_source_state_id"] = "two"
        errors = self.runner.validate_manifest(manifest)
        self.assertTrue(any("minimum_source_state_id >= 2" in error
                            for error in errors), errors)

    def test_example_is_explicitly_non_runnable(self):
        with tempfile.TemporaryDirectory(prefix="srcsep-d03-template-") as tmp:
            temporary = Path(tmp)
            amps = temporary / "AMPS"
            (amps / "srcSEP").mkdir(parents=True)
            (amps / "srcSEP" / "makefile").write_text(
                "# template gate must return before make\n", encoding="utf-8")
            (amps / "Makefile.conf").write_text(
                "# presence is sufficient for this pre-build test\n", encoding="utf-8")
            manifest = temporary / "manifest.json"
            manifest.write_text(json.dumps(self.manifest), encoding="utf-8")
            arguments = argparse.Namespace(
                amps_source=str(amps), make_config=str(amps / "Makefile.conf"),
                manifest=str(manifest), rebuild=False, no_build=False,
                timeout=None, release=False,
                output_dir=str(temporary / "evidence"))
            with contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.runner.run_gate(arguments), 0)
            evidence = json.loads((temporary / "evidence" /
                "native_integration_evidence.json").read_text(encoding="utf-8"))
            self.assertEqual(evidence["status"], "SKIP")
            self.assertIn("template_only", evidence["message"])
            self.assertEqual(evidence["build"], {})

    def test_runtime_source_wires_post_return_counter_to_mpi_summary(self):
        """Guard the source boundary while the configured build is external.

        This does not claim native execution; the release campaign remains the
        authoritative proof. It catches accidental removal of the post-return
        increment or rank reduction early enough to keep the resulting D03
        failure diagnostic actionable.
        """
        runtime = (ROOT / "production_mover_runtime.cpp").read_text(
            encoding="utf-8")
        main = (ROOT / "main.cpp").read_text(encoding="utf-8")
        implementation_call = runtime.index(
            "g_selected_implementation(ptr,dtTotal,startNode)")
        counter_increment = runtime.index("g_completed_dispatches.fetch_add")
        self.assertLess(implementation_call, counter_increment)
        self.assertIn("SEP::Mover::CompletedDispatchCount()", main)
        self.assertIn("MPI_Reduce(&localCompletedDispatches", main)
        self.assertIn('" completed_particle_dispatches="', main)

    def test_mover_aliases_are_reserved_for_runner_selection(self):
        for arguments in (["--particle-mover", "parker"],
                          ["--particle-mover=parker"],
                          ["--mover", "fte-dmumu"],
                          ["--mover=fte-dmumu"],
                          ["--sep-mover", "fte-mfp"],
                          ["--sep-mover=fte-mfp"]):
            self.assertTrue(self.runner._contains_mover_option(arguments))
        self.assertFalse(self.runner._contains_mover_option(
            ["--total-iterations", "20"]))

    def test_orchestrator_records_complete_synthetic_campaign(self):
        """Exercise evidence plumbing without presenting it as native physics.

        The tiny shell executable emits the metadata grammar and deterministic
        files expected by the Python layer. It does not call a transport kernel;
        the configured D03 release campaign remains the only native claim.
        """
        with tempfile.TemporaryDirectory(prefix="srcsep-d03-orchestrator-") as tmp:
            temporary = Path(tmp)
            amps = temporary / "AMPS"
            (amps / "srcSEP").mkdir(parents=True)
            (amps / "srcSEP" / "makefile").write_text(
                "# --no-build orchestrator fixture\n", encoding="utf-8")
            make_config = amps / "Makefile.conf"
            make_config.write_text("# fixture\n", encoding="utf-8")
            executable = amps / "amps"
            executable.write_text(
                "#!/bin/sh\n"
                "mover=unknown\n"
                "while test $# -gt 0; do\n"
                "  if test \"$1\" = --particle-mover; then shift; mover=$1; fi\n"
                "  shift\n"
                "done\n"
                "test \"$D03_CASE_DIR\" = \"$PWD\" || exit 9\n"
                "case $mover in\n"
                "  parker) runfp=1111;;\n"
                "  fte-dmumu) runfp=2222;;\n"
                "  fte-mfp) runfp=3333;;\n"
                "  *) exit 8;;\n"
                "esac\n"
                "printf 'state:%s\\n' \"$mover\" > final_state.bin\n"
                "printf 'checkpoint\\n' > checkpoint.dat\n"
                "printf 'RunConfiguration fingerprint=%s\\n' \"$runfp\"\n"
                "echo 'SWCME configuration fingerprint=cafebabe'\n"
                "echo 'final_source_state_id=2 completed_particle_dispatches=1 mpi_consensus=pass'\n",
                encoding="utf-8")
            executable.chmod(0o755)
            launcher = temporary / "fake-mpi"
            launcher.write_text(
                "#!/bin/sh\nshift\nexec \"$@\"\n", encoding="utf-8")
            launcher.chmod(0o755)

            manifest = copy.deepcopy(self.manifest)
            manifest.pop("template_only")
            manifest["executable"] = str(executable)
            manifest["common_args"] = []
            manifest["mpi_launcher"] = [str(launcher), "{ranks}"]
            manifest["environment"] = {"D03_CASE_DIR": "{case_dir}"}
            for case in manifest["cases"]:
                case["args"] = []
            manifest_path = temporary / "site-manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            output = temporary / "evidence"
            arguments = argparse.Namespace(
                amps_source=str(amps), make_config=str(make_config),
                manifest=str(manifest_path), rebuild=False, no_build=True,
                timeout=10.0, release=False, output_dir=str(output))
            with contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.runner.run_gate(arguments), 0)
            evidence = json.loads((output /
                "native_integration_evidence.json").read_text(encoding="utf-8"))
            self.assertEqual(evidence["status"], "PASS")
            self.assertEqual(len(evidence["cases"]), len(manifest["cases"]))
            self.assertTrue(all(case["status"] == "PASS"
                                for case in evidence["cases"]))
            self.assertTrue(all(group["status"] == "PASS"
                                for group in evidence["artifact_groups"].values()))

    def test_missing_native_tree_skips_development_but_fails_release(self):
        with tempfile.TemporaryDirectory(prefix="srcsep-d03-contract-") as tmp:
            temporary = Path(tmp)
            missing = temporary / "missing-amps"
            base = dict(
                amps_source=str(missing),
                make_config=str(missing / "Makefile.conf"),
                manifest=None, rebuild=False, no_build=False, timeout=None)
            development = argparse.Namespace(
                **base, release=False, output_dir=str(temporary / "development"))
            release = argparse.Namespace(
                **base, release=True, output_dir=str(temporary / "release"))
            with contextlib.redirect_stdout(io.StringIO()), \
                 contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(self.runner.run_gate(development), 0)
                self.assertEqual(self.runner.run_gate(release), 2)
            dev_evidence = json.loads((temporary / "development" /
                "native_integration_evidence.json").read_text(encoding="utf-8"))
            release_evidence = json.loads((temporary / "release" /
                "native_integration_evidence.json").read_text(encoding="utf-8"))
            self.assertEqual(dev_evidence["status"], "SKIP")
            self.assertEqual(release_evidence["status"], "INCOMPLETE")


if __name__ == "__main__":
    unittest.main(verbosity=2)
