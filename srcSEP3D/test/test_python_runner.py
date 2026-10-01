#!/usr/bin/env python3
"""Dependency-free unit tests for the srcSEP3D Python runner contract.

These tests exercise selection and reporting in memory.  They do not compile
the C++ suite or invoke AMPS, which keeps RUN3D01 fast enough for --routine.
Subprocess behavior and real C++ reports remain covered by HARN/BLD tests.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


RUNNER_PATH = Path(__file__).with_name("run_tests.py")
SPEC = importlib.util.spec_from_file_location("srcsep3d_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class SelectionTests(unittest.TestCase):
    def parse(self, *tokens: str):
        return runner._parser().parse_args(list(tokens))

    def test_individual_ids_are_case_insensitive_and_deduplicated(self):
        selected = runner._select(self.parse(
            "--test", "lay01", "--test", "LAY01", "--test", "BLDL3D02"))
        self.assertEqual([item.test_id for item in selected],
                         ["BLDL3D02", "LAY01"])

    def test_groups_expand_in_stable_id_order(self):
        selected = runner._select(self.parse("--group", "harn"))
        self.assertEqual([item.test_id for item in selected],
                         ["HARN01", "HARN02", "HARN03", "HARN04"])

    def test_routine_excludes_outer_shell_probes(self):
        selected = runner._select(self.parse("--routine"))
        ids = {item.test_id for item in selected}
        self.assertIn("RUN3D01", ids)
        self.assertNotIn("HARN02-EXITCODE", ids)
        self.assertNotIn("HARN03-EXITCODE", ids)

    def test_all_contains_every_public_definition(self):
        selected = runner._select(self.parse("--all"))
        self.assertEqual({item.test_id for item in selected}, set(runner.BY_ID))

    def test_conflicting_or_unknown_selection_is_usage_error(self):
        with self.assertRaises(runner.RunnerError):
            runner._select(self.parse("--routine", "--all"))
        with self.assertRaises(runner.RunnerError):
            runner._select(self.parse("--test", "NO_SUCH_TEST"))


class ReportTests(unittest.TestCase):
    def test_json_and_junit_have_consistent_counts(self):
        results = [
            runner.Result("PASS01", "UNIT", "PASS", "ok", 0.1, []),
            runner.Result("SKIP01", "UNIT", "SKIP", "not configured", 0.0, []),
        ]
        with tempfile.TemporaryDirectory(prefix="srcsep3d-runner-") as tmp:
            destination = Path(tmp)
            runner._write_reports(results, destination)
            payload = json.loads(
                (destination / "srcsep3d-tests.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["totals"],
                             {"pass": 1, "fail": 0, "skip": 1, "error": 0})
            self.assertEqual(payload["exit_code"], 0)
            xml = (destination / "srcsep3d-tests.xml").read_text(encoding="utf-8")
            self.assertIn('tests="2"', xml)
            self.assertIn('skipped="1"', xml)


class ProductionBuildRoutingTests(unittest.TestCase):
    def test_configured_build_delegates_to_enclosing_amps_root(self):
        """BLDL3D01 must not compile AMPS-dependent sources in srcSEP3D."""
        with tempfile.TemporaryDirectory(prefix="srcsep3d-production-") as tmp:
            amps_root = Path(tmp)
            config = amps_root / "Makefile.conf"
            (amps_root / "Makefile").write_text(
                "# enclosing AMPS build fixture\n", encoding="utf-8")
            config.write_text("# configuration fixture\n", encoding="utf-8")
            args = runner._parser().parse_args([
                "--test", "BLDL3D01",
                "--amps-source", str(amps_root),
                "--make-config", str(config),
            ])
            definition = runner.BY_ID["BLDL3D01"]

            # The subprocess is mocked because this unit test verifies routing,
            # while BLDL3D01 itself performs the real enclosing build.
            with mock.patch.object(
                    runner, "_run_command", return_value=(0, "", 0.25)) as call:
                result = runner._check_production_build(definition, args)

            command = call.call_args.args[0]
            self.assertEqual(result.status, "PASS")
            self.assertIn("strict-production", command)
            self.assertIn(f"AMPS_ROOT={amps_root.resolve()}", command)
            self.assertIn(f"AMPS_CONFIG={config.resolve()}", command)


class ApplicationIsolationTests(unittest.TestCase):
    def test_runner_never_constructs_a_srcsep_sibling_path(self):
        """srcSEP3D tests may use shared models, but never the srcSEP app."""
        source = RUNNER_PATH.read_text(encoding="utf-8")
        self.assertNotIn('ROOT.parent / "srcSEP"', source)
        self.assertNotIn("ROOT.parent / 'srcSEP'", source)


class PicHeaderGuardTests(unittest.TestCase):
    def make_fixture(self, directory: Path):
        """Only source headers are installed: no tools directory or helper."""
        pic = directory / "src/pic"
        (pic / "ecsim").mkdir(parents=True)
        (pic / "pic.h").write_text(
            "#ifndef _PIC_\n#define _PIC_\n"
            "//include headers for individual physical models\n"
            '#include "ecsim/domain_bc.h"\n'
            '#include "gyro/gyro_mover.h"\n'
            "#endif // _PIC_\n", encoding="utf-8")
        (pic / "ecsim/domain_bc.h").write_text(
            "#ifndef AMPS_PIC_ECSIM_DOMAIN_BC_H_INCLUDED\n"
            "#define AMPS_PIC_ECSIM_DOMAIN_BC_H_INCLUDED\n"
            "#pragma once\n"
            "class cDomainBC {};\n"
            "#endif // AMPS_PIC_ECSIM_DOMAIN_BC_H_INCLUDED\n", encoding="utf-8")
        args = runner._parser().parse_args([
            "--test", "BLDL3D12", "--amps-source", str(directory)])
        return pic, args

    def test_inline_check_uses_selected_source_without_external_tool(self):
        with tempfile.TemporaryDirectory(prefix="srcsep3d-guards-") as tmp:
            directory = Path(tmp)
            pic, args = self.make_fixture(directory)
            outcomes = [(0, "", 0.01), (1, "default argument", 0.01),
                        (0, "", 0.01), (0, "", 0.01), (1, "redefinition", 0.01)]
            # Only compiler execution is mocked. Source discovery, guard
            # analysis, fixture creation and both negative controls are real.
            with mock.patch.object(runner, "_run_command", side_effect=outcomes) as call:
                result = runner._run_source(runner.BY_ID["BLDL3D12"], args,
                                            directory / "reports")
            self.assertEqual(result.status, "PASS", result.message)
            self.assertFalse((directory / "tools").exists())
            self.assertEqual(call.call_count, 5)
            for item in call.call_args_list:
                self.assertEqual(item.args[0][0], args.cxx)
                self.assertIn("-fsyntax-only", item.args[0])
            work = directory / "reports/pic-header-guards-probe"
            self.assertEqual((work / "unguarded-domain.log").read_text(), "redefinition")
            self.assertIn("class cDomainBC", (pic / "ecsim/domain_bc.h").read_text())

    def test_early_pic_guard_is_rejected_before_compilation(self):
        with tempfile.TemporaryDirectory(prefix="srcsep3d-guards-") as tmp:
            directory = Path(tmp)
            pic, args = self.make_fixture(directory)
            text = (pic / "pic.h").read_text()
            text = text.replace("//include headers", "#endif\n//include headers")
            (pic / "pic.h").write_text(text)
            with mock.patch.object(runner, "_run_command") as call:
                result = runner._run_source(runner.BY_ID["BLDL3D12"], args,
                                            directory / "reports")
            self.assertEqual(result.status, "FAIL")
            self.assertIn("tail lies outside _PIC_", result.message)
            call.assert_not_called()

    def test_negative_control_cannot_silently_pass(self):
        with tempfile.TemporaryDirectory(prefix="srcsep3d-guards-") as tmp:
            directory = Path(tmp)
            _, args = self.make_fixture(directory)
            with mock.patch.object(runner, "_run_command", return_value=(0, "", 0.01)):
                result = runner._run_source(runner.BY_ID["BLDL3D12"], args,
                                            directory / "reports")
            self.assertEqual(result.status, "FAIL")
            self.assertIn("unguarded-model-tail compile contract failed", result.message)


if __name__ == "__main__":
    unittest.main(verbosity=2)
