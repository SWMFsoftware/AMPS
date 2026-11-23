#!/usr/bin/env python3
"""Focused tests for the dependency-light Python test orchestrator.

The fixture deliberately contains both supported analytical evidence forms:
an ordinary registry metric and a CSV solution series.  The analytical values
are literal independent fixture data, not values computed by a srcSEP kernel,
so this test can detect plotting/report regressions without validating a
production numerical routine against itself.
"""

from __future__ import annotations

import importlib.util
import io
import json
import hashlib
import os
from pathlib import Path
import subprocess
import re
import sys
import tempfile
import textwrap
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "test" / "run_tests.py"
CASE_RUNNER = ROOT / "validation" / "run_case.py"


def _load_runner_module():
    """Load the runner by path without requiring test/ to be a package."""
    specification = importlib.util.spec_from_file_location(
        "srcsep_python_test_runner", RUNNER)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import runner from {RUNNER}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _load_cross_model_runner():
    """Load the case scorer while preserving its sibling-module imports."""
    case_directory = ROOT / "validation" / "cases"
    sys.path.insert(0, str(case_directory))
    try:
        specification = importlib.util.spec_from_file_location(
            "srcsep_cross_model_case_runner",
            case_directory / "cross_model_case_runner.py")
        if specification is None or specification.loader is None:
            raise RuntimeError("cannot import cross-model case runner")
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


class PythonTestRunnerTests(unittest.TestCase):
    def test_new_source_family_aliases_use_existing_strict_make_runners(self):
        """Catalog aliases select audited runners; they do not rebuild AMPS.

        Assert the real Make recipe and existing source runner, so a misspelled
        alias cannot pass a table-only test while remaining unreachable.
        """
        runner = _load_runner_module()
        makefile = (ROOT / "makefile").read_text()
        for suite, target, script in (
                ("analytic-radial", "test-analytic-radial-background-unit", "run_analytic_radial_background_tests.sh"),
                ("analytic-shock-epoch", "test-analytic-shock-epoch-unit", "run_analytic_shock_epoch_tests.sh"),
                ("finite-amplitude-study", "test-finite-amplitude-study-unit", "run_finite_amplitude_study_tests.sh"),
                ("material-tube-kinematics", "test-material-tube-kinematics-unit", "run_material_tube_kinematics_tests.sh"),
                ("standard-batsrus-waves", "test-standard-batsrus-wave-unit", "run_standard_batsrus_wave_tests.sh"),
                ("flux-tube-state", "test-flux-tube-state-unit", "run_flux_tube_state_tests.sh"),
                ("polarized-wkb", "test-polarized-wkb-unit", "run_polarized_wkb_tests.sh"),
                ("maxwell-juttner-flux", "test-maxwell-juttner-flux-unit", "run_maxwell_juttner_flux_tests.sh")):
            with self.subTest(suite=suite):
                self.assertEqual(runner.SOURCE_SUITES[suite], target)
                # These two tracked POSIX scripts have mode 0644. Their
                # strict numerical runner must be invoked through Make's
                # shell so a source archive preserves the same behavior as
                # the checkout. Keep the exact target/script assertion and
                # failure propagation; only account for the required launcher.
                launcher = "$(SHELL) " if script in (
                    "run_standard_batsrus_wave_tests.sh",
                    "run_maxwell_juttner_flux_tests.sh") else ""
                self.assertIn(target + ":\n\t" + launcher + "./test/" + script, makefile)
                self.assertTrue((ROOT / "test" / script).is_file())

    def test_fct_source_guard_discards_compile_and_execution_mutations(self):
        """Actual shell gate rejects changed source, even after a child PASS.

        A recording fake compiler works only inside a disposable mirrored
        fixture. This tests the source/evidence validity boundary without
        recompiling physics or mutating any real checkout file.
        """
        source_script=ROOT/"test/run_frozen_collision_transaction_tests.sh"
        script=source_script.read_text()
        originals=re.findall(r'"((?:util/|test/|\.\./src/models/sep_common/)[^"\n]+)"',script)
        protected=ROOT/"util/sep_frozen_collision_transaction.h"
        original_hash=hashlib.sha256(protected.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory(prefix="srcsep-source-guard-") as temporary:
            temporary=Path(temporary)
            app=temporary/"tree/srcSEP"
            for relative in set(originals):
                target=(app/relative).resolve()
                self.assertIn(temporary,target.parents)
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_bytes((ROOT/relative).read_bytes())
            runner=app/"test/run_frozen_collision_transaction_tests.sh"
            runner.chmod(0o755)
            mutated=app/"util/sep_frozen_collision_transaction.h"
            fixture_original=mutated.read_bytes()
            marker=temporary/"executed"
            compiler=temporary/"fake-cxx"
            compiler.write_text("#!"+sys.executable+"\n"+textwrap.dedent('''\
                import os
                from pathlib import Path
                import sys
                target=Path(os.environ["SRCSEP_GUARD_TARGET"])
                output=Path(sys.argv[sys.argv.index("-o")+1])
                if os.environ["SRCSEP_GUARD_PHASE"]=="compile":
                    target.write_text(target.read_text()+"\\n// simulated concurrent edit\\n")
                program=("#!"+sys.executable+"\\nimport os\\nfrom pathlib import Path\\n"
                    "Path(os.environ['SRCSEP_GUARD_MARKER']).write_text('executed')\\n"
                    "target=Path(os.environ['SRCSEP_GUARD_TARGET'])\\n"
                    "target.write_text(target.read_text()+'\\\\n// simulated execution edit\\\\n')\\n"
                    "print('[FCT01] PASS fixture child verdict; source guard must discard it')\\n")
                output.write_text(program)
                output.chmod(0o755)
                '''),encoding="utf-8")
            compiler.chmod(0o755)
            for phase,diagnostic in (("compile","compilation"),("execute","execution")):
                with self.subTest(phase=phase):
                    mutated.write_bytes(fixture_original);marker.unlink(missing_ok=True)
                    environment=dict(os.environ,CXX=str(compiler),PYTHONDONTWRITEBYTECODE="1",
                        SRCSEP_GUARD_TARGET=str(mutated),SRCSEP_GUARD_MARKER=str(marker),SRCSEP_GUARD_PHASE=phase)
                    child=subprocess.run([str(runner)],cwd=str(temporary),env=environment,text=True,
                        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False)
                    self.assertEqual(child.returncode,2,child.stdout)
                    self.assertIn("source changed during "+diagnostic,child.stdout)
                    self.assertEqual(marker.exists(),phase=="execute")
                    if phase=="execute":self.assertIn("[FCT01] PASS",child.stdout)
            self.assertEqual(hashlib.sha256(protected.read_bytes()).hexdigest(),original_hash)

    def test_native_cwd_isolation_covers_discovery_and_relative_input_assets(self):
        """A real child deliberately writes PIC names; srcSEP stays unchanged.

        This is an infrastructure fault fixture, not model validation. Writing
        before --list-tests returns makes a discovery regression observable as
        well as an execution regression. Relative CSV paths intentionally use
        the historic native srcSEP base, not the manifest's parent directory.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-native-cwd-") as temporary:
            root = Path(temporary)
            source, caller, output = root / "source", root / "caller", root / "evidence"
            (source / "PT").mkdir(parents=True)
            caller.mkdir()
            protected = source / "PT" / "thread=0.log"
            protected.write_text("source sentinel\n", encoding="utf-8")
            (source / "fixtures").mkdir()
            asset = source / "particles.csv"
            asset.write_text("id,mu\n8,0.25\n", encoding="utf-8")
            manifest = source / "fixtures" / "case.txt"
            manifest_bytes = ("srcsep-cv01-native-args-v1\n--initial-particles\n"
                              "particles.csv\n--dt-s\n0.125\n")
            manifest.write_text(manifest_bytes, encoding="utf-8")
            trace = root / "calls.jsonl"
            executable = root / "amps"
            executable.write_text("#!" + sys.executable + "\n" + textwrap.dedent('''\
                import json
                import os
                from pathlib import Path
                import sys
                args=sys.argv[1:]
                cwd=Path.cwd()
                (cwd/"PT").mkdir(exist_ok=True)
                (cwd/"PT/thread=0.log").write_text("native scratch")
                (cwd/"mesh.cache").write_text("native cache")
                record={"cwd":str(cwd),"arguments":args}
                if "--test-input" in args:
                    lines=Path(args[args.index("--test-input")+1]).read_text().splitlines()
                    leaf=Path(lines[lines.index("--initial-particles")+1])
                    record["asset_path"]=str(leaf)
                    record["asset_bytes"]=leaf.read_text()
                with Path(os.environ["SRCSEP_CWD_TRACE"]).open("a") as stream:
                    stream.write(json.dumps(record)+"\\n")
                if "--list-tests" in args:
                    print("PCLI01 | cli | deterministic | infrastructure only")
                    raise SystemExit(0)
                report={"schema":"srcsep-component-tests-v1","exit_code":0,
                    "totals":{"passed":1,"failed":0,"skipped":0,"errors":0},
                    "results":[{"id":"PCLI01","status":"PASS","message":"cwd fixture",
                        "elapsed_seconds":0,"seed":None,"configuration":[],"metrics":[],
                        "artifacts":[str(cwd/"mesh.cache")]}]}
                Path(args[args.index("--test-json")+1]).write_text(json.dumps(report))
                Path(args[args.index("--test-junit")+1]).write_text("<testsuite/>")
                '''), encoding="utf-8")
            executable.chmod(0o755)
            previous = Path.cwd()
            try:
                # Caller-owned --amps and --output-dir retain caller-relative
                # semantics while forwarded native paths retain source base.
                os.chdir(caller)
                with mock.patch.object(runner, "ROOT", source), \
                        mock.patch.object(runner, "_validation_case_ids", return_value=set()), \
                        mock.patch.dict(os.environ, {"SRCSEP_CWD_TRACE":str(trace),
                                                     "PYTHONDONTWRITEBYTECODE":"1"}), \
                        mock.patch("sys.stdout", new_callable=io.StringIO):
                    for selection, label in ((["--test","PCLI01"],"selected"),(["--all"],"all")):
                        self.assertEqual(runner.main([
                            "--amps","../amps",*selection,"--output-dir","../evidence/"+label,
                            "--plot","none","--","--test-input","fixtures/case.txt"]),0)
            finally:
                os.chdir(previous)
            records = [json.loads(line) for line in trace.read_text().splitlines()]
            self.assertEqual(len(records),3)  # selected + all discovery + all case
            directories = [Path(record["cwd"]) for record in records]
            self.assertEqual(len(set(directories)),3)
            for record, directory in zip(records,directories):
                self.assertIn(output,directory.parents)
                self.assertEqual(record["asset_bytes"],asset.read_text())
                self.assertIn(directory,Path(record["asset_path"]).parents)
                provenance=json.loads((directory/"native-execution.json").read_text())
                self.assertEqual(provenance["native_relative_path_base"],str(source))
                self.assertTrue(any(item["original"]==str(asset) for item in provenance["input_artifacts"]))
                self.assertTrue((directory/"PT/thread=0.log").is_file())
            self.assertEqual(protected.read_text(),"source sentinel\n")
            self.assertFalse((source/"mesh.cache").exists())
            self.assertEqual(manifest.read_text(),manifest_bytes)
            self.assertEqual(asset.read_text(),"id,mu\n8,0.25\n")

    def test_native_unknown_nested_schema_fails_without_launch_or_path_guess(self):
        """An undeclared nested path base is a hard infrastructure error."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-native-schema-") as temporary:
            root=Path(temporary)
            source=root/"source";source.mkdir()
            (source/"case.txt").write_text("new-unknown-schema\nrelative/data.csv\n")
            args=SimpleNamespace(model_args=["--test-input=case.txt"])
            with mock.patch.object(runner,"ROOT",source), \
                    mock.patch.object(runner.subprocess,"Popen") as launch:
                with self.assertRaisesRegex(runner.RunnerError,"no declared nested-path isolation contract"):
                    runner._native_scope(args,root/"evidence","unknown")
                launch.assert_not_called()
            self.assertEqual((source/"case.txt").read_text(),"new-unknown-schema\nrelative/data.csv\n")

    def test_linked_validation_native_launch_and_discovery_own_their_scratch(self):
        """Exercise the shared CV/IV/XM launcher through real subprocesses.

        The fixture publishes literal minimal CSV and PASS records only. It
        proves executable selection and file ownership, with no numerical
        acceptance claim for the named case.
        """
        def load(name, path):
            specification=importlib.util.spec_from_file_location(name,path)
            module=importlib.util.module_from_spec(specification)
            specification.loader.exec_module(module)
            return module
        linked=load("srcsep_linked_scope_fixture",ROOT/"validation/cases/linked_case_common.py")
        cases=load("srcsep_case_scope_fixture",CASE_RUNNER)
        with tempfile.TemporaryDirectory(prefix="srcsep-linked-cwd-") as temporary:
            root=Path(temporary);source=root/"source";source.mkdir()
            (source/"history.csv").write_text("time_s,rate\n0,1\n")
            output=root/"evidence";output.mkdir()
            executable=root/"amps"
            executable.write_text("#!"+sys.executable+"\n"+textwrap.dedent('''\
                import json
                from pathlib import Path
                import sys
                args=sys.argv[1:]
                Path("PT").mkdir(exist_ok=True)
                Path("PT/thread=0.log").write_text("owned scratch")
                if "--list-tests" in args:
                    print("CV02 | controlled-analytical | deterministic")
                    raise SystemExit(0)
                cid=args[args.index("--test")+1]
                lines=Path(args[args.index("--test-input")+1]).read_text().splitlines()
                assert Path(lines[lines.index("--source-history-csv")+1]).read_text()=="time_s,rate\\n0,1\\n"
                model=Path(args[args.index("--test-output-dir")+1])/(cid+"_model.csv")
                model.write_text("fixture,value\\n1,2\\n")
                Path(args[args.index("--test-json")+1]).write_text(json.dumps(
                    {"results":[{"id":cid,"status":"PASS"}]}))
                Path(args[args.index("--test-junit")+1]).write_text("<testsuite/>")
                '''),encoding="utf-8")
            executable.chmod(0o755)
            discovery=output/"discovery";discovery.mkdir()
            cases._verify_linked_registry(executable,["CV02"],10,discovery)
            with mock.patch("sys.stdout",new_callable=io.StringIO):
                result=linked.run_linked_model(case_id="CV02",arguments=["--source-history-csv","history.csv"],
                    source_root=source,output_dir=output,executable=executable,timeout=10)
            self.assertFalse((source/"PT").exists())
            self.assertTrue((discovery/"PT/thread=0.log").exists())
            self.assertIn("WORKING_DIRECTORY:",result["log"].read_text())
            records=list((output/"native/native-work").glob("*/native-execution.json"))
            self.assertEqual(len(records),1)
            provenance=json.loads(records[0].read_text())
            self.assertEqual(provenance["native_relative_path_base"],str(source))
            self.assertTrue((records[0].parent/"PT/thread=0.log").exists())
            self.assertTrue(any(item["original"]==str(source/"history.csv")
                                for item in provenance["input_artifacts"]))

    def test_help_contains_annotated_selection_and_plot_examples(self):
        """Keep the command-line quick-start available to archive users."""
        completed = subprocess.run(
            [sys.executable, str(RUNNER), "--help"], cwd=str(ROOT), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        self.assertEqual(completed.returncode, 0, completed.stdout)
        # These anchors cover the distinct native, source-only, report-only,
        # MPI, plotting, and forwarded-model-argument workflows.  Exact spacing
        # is intentionally not asserted because argparse may wrap descriptions
        # differently when a Python distribution changes terminal width.
        for expected in (
                "Discover the tests registered",
                "--test PARK07",
                "--amps ../amps --validation-case CV01",
                "--validation-case XM01 --validation-case XM02",
                "validation/cases/XM02/publication_input.json",
                "validation/cases/XM03/publication_input.json",
                "validation/cases/XM03/input/earth_shock_thermal_source.csv",
                "validation/cases/XM03/reference/liu_figure12_earth_observations.csv",
                "Do not add --case-input",
                "No external production CSV is used by XM02",
                "XM03 needs no external model CSV",
                "--group parker --group fte-dmumu",
                "--routine",
                "--all",
                "D01, D02, and D03PRE",
                "--suite controlled-analytical",
                "--suite d03-native-integration",
                "--mpi-np 4",
                "--amps-source",
                "--make-config",
                "--native-manifest",
                "--rebuild",
                "--release",
                "--from-json results.json",
                "--plot none",
                "--mover parker --turbulence-source prescribed",
                "Outputs:"):
            self.assertIn(expected, completed.stdout)

    def test_registry_list_parser_ignores_headers_and_deduplicates(self):
        runner = _load_runner_module()
        # This reproduces the real C++ registry preamble and column heading.
        # The historical parser accepted ``ID`` and later asked AMPS to run an
        # unknown test called "id", terminating a monolithic --all command.
        listing = """Registered srcSEP standalone component tests
ID | group | class | initialization | build modes | description
PARK01 | parker | deterministic
FTED08 | fte-dmumu | extended
PARK01 | parker | duplicate
"""
        self.assertEqual(runner._parse_list_output(listing),
                         ["FTED08", "PARK01"])

    def test_native_all_cli_uses_only_requested_configured_build(self):
        """Regress both native commands that reached 3-D-only Make targets.

        Launch the installed runner as a separate process, rather than mocking
        its native helper. A tiny protocol fixture advertises one deterministic
        test and writes literal JSON; it verifies orchestration, not transport
        physics. A recording make on PATH rejects implicit builds and accepts
        only the enclosing clean/strict-production route under --rebuild.
        Substituting the 3-D runner must fail for both test/stage1 and
        clean-standalone; no real configured AMPS build occurs in this fixture.
        """
        with tempfile.TemporaryDirectory(prefix="srcsep-native-all-cli-") as tmp:
            temporary = Path(tmp)
            executable = temporary / "amps"
            calls = temporary / "native-calls.jsonl"
            make_called = temporary / "make-called"
            amps_source = temporary / "AMPS"
            application_makefile = amps_source / "srcSEP" / "makefile"
            application_makefile.parent.mkdir(parents=True)
            application_makefile.write_text("# protocol fixture\n", encoding="utf-8")
            make_config = amps_source / "Makefile.conf"
            make_config.write_text("# protocol fixture\n", encoding="utf-8")
            # sys.executable is the interpreter actually running this test;
            # it avoids assuming a particular Python alias on the host.
            executable.write_text("#!" + sys.executable + "\n" + textwrap.dedent('''\
                import json
                import os
                from pathlib import Path
                import sys
                args = sys.argv[1:]
                with Path(os.environ["SRCSEP_RUNNER_FIXTURE_CALLS"]).open("a") as log:
                    log.write(json.dumps(args) + "\\n")
                if args == ["--list-tests"]:
                    print("Registered srcSEP standalone component tests")
                    print("ID | group | class | description")
                    print("PCLI01 | cli | deterministic | protocol fixture")
                    raise SystemExit(0)
                if "--test" not in args or args[args.index("--test") + 1] != "PCLI01":
                    raise SystemExit("unexpected native fixture command")
                report = {
                    "schema": "srcsep-component-tests-v1", "exit_code": 0,
                    "totals": {"passed": 1, "failed": 0, "skipped": 0, "errors": 0},
                    "results": [{"id": "PCLI01", "status": "PASS",
                        "message": "runner protocol fixture; no physics claim",
                        "elapsed_seconds": 0.0, "seed": None,
                        "configuration": [], "metrics": [], "artifacts": []}],
                }
                Path(args[args.index("--test-json") + 1]).write_text(json.dumps(report))
                Path(args[args.index("--test-junit") + 1]).write_text(
                    '<testsuite tests="1" failures="0" errors="0"/>')
                '''), encoding="utf-8")
            executable.chmod(0o755)
            make = temporary / "make"
            make.write_text("#!" + sys.executable + "\n" + textwrap.dedent('''\
                import json
                import os
                from pathlib import Path
                import sys
                args = sys.argv[1:]
                with Path(os.environ["SRCSEP_RUNNER_FIXTURE_MAKE_CALLED"]).open("a") as log:
                    log.write(json.dumps(args) + "\\n")
                for target in ("test/stage1", "clean-standalone"):
                    if target in args:
                        print("make: *** No rule to make target '" + target + "'. Stop.")
                        raise SystemExit(2)
                if os.environ["SRCSEP_RUNNER_FIXTURE_ALLOW_REBUILD"] == "1":
                    root = Path(os.environ["SRCSEP_RUNNER_FIXTURE_AMPS_ROOT"])
                    if args == ["-C", str(root), "clean"]:
                        raise SystemExit(0)
                    if args == ["--no-print-directory", "-f", str(root / "srcSEP/makefile"),
                                "strict-production", "AMPS_ROOT=" + str(root),
                                "AMPS_CONFIG=" + str(root / "Makefile.conf")]:
                        raise SystemExit(0)
                print("unexpected implicit or misrouted make command: " + str(args))
                raise SystemExit(2)
                '''), encoding="utf-8")
            make.chmod(0o755)
            environment = dict(os.environ)
            environment.update({
                "PATH": str(temporary) + os.pathsep + environment.get("PATH", ""),
                "PYTHONDONTWRITEBYTECODE": "1",
                "SRCSEP_RUNNER_FIXTURE_CALLS": str(calls),
                "SRCSEP_RUNNER_FIXTURE_MAKE_CALLED": str(make_called),
                "SRCSEP_RUNNER_FIXTURE_AMPS_ROOT": str(amps_source),
            })
            for rebuild in (False, True):
                with self.subTest(rebuild=rebuild):
                    calls.unlink(missing_ok=True)
                    make_called.unlink(missing_ok=True)
                    environment["SRCSEP_RUNNER_FIXTURE_ALLOW_REBUILD"] = str(int(rebuild))
                    output = temporary / "test_output" / ("rebuilt" if rebuild else "all")
                    command = [sys.executable, str(RUNNER), "--amps", str(executable),
                               "--all", "--output-dir", str(output)]
                    if rebuild:
                        command.extend(["--rebuild", "--amps-source", str(amps_source),
                                        "--make-config", str(make_config)])
                    completed = subprocess.run(command,
                        cwd=str(ROOT), env=environment, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        timeout=30, check=False)
                    self.assertEqual(completed.returncode, 0, completed.stdout)
                    self.assertEqual(make_called.exists(), rebuild, completed.stdout)
                    native_calls = [json.loads(line) for line in calls.read_text().splitlines()]
                    self.assertEqual(len(native_calls), 2)
                    self.assertEqual(native_calls[0], ["--list-tests"])
                    self.assertEqual(native_calls[1][:2], ["--test", "PCLI01"])
                    report = json.loads((output / "srcsep-tests.json").read_text())
                    self.assertEqual(report["totals"], {
                        "passed": 1, "failed": 0, "skipped": 0, "errors": 0})
                    manifest = json.loads((output / "run_manifest.json").read_text())
                    self.assertEqual(len(manifest["command"]), 3 if rebuild else 1)
                    self.assertEqual(manifest["command"][-1][0], str(executable))
                    if rebuild:
                        make_calls = [json.loads(line) for line in make_called.read_text().splitlines()]
                        self.assertEqual(make_calls, [item[1:] for item in manifest["command"][:2]])

    def test_all_rebuild_uses_enclosing_configured_build_commands(self):
        """Make newly linked registry IDs visible before --all discovery.

        The clean/build commands are intentionally tested without invoking an
        actual AMPS tree.  Their exact routing is the important contract: both
        commands operate on the enclosing source root, and strict-production
        receives absolute AMPS_ROOT/AMPS_CONFIG values.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-all-rebuild-") as tmp:
            amps_source = Path(tmp) / "AMPS"
            (amps_source / "srcSEP").mkdir(parents=True)
            application_makefile = amps_source / "srcSEP" / "makefile"
            application_makefile.write_text("strict-production:\n\t@true\n",
                                             encoding="utf-8")
            make_config = amps_source / "Makefile.conf"
            make_config.write_text("# configured fixture\n", encoding="utf-8")
            arguments = SimpleNamespace(
                amps_source=amps_source, make_config=make_config,
                timeout=42.0)
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                commands = runner._rebuild_linked_application(
                    arguments, Path(tmp) / "rebuild.log")

        self.assertEqual(commands[0], [
            "make", "-C", str(amps_source.resolve()), "clean"])
        self.assertEqual(commands[1][:5], [
            "make", "--no-print-directory", "-f",
            str(application_makefile.resolve()), "strict-production"])
        self.assertIn(f"AMPS_ROOT={amps_source.resolve()}", commands[1])
        self.assertIn(f"AMPS_CONFIG={make_config.resolve()}", commands[1])
        self.assertEqual(run.call_count, 2)

    def test_all_accepts_rebuild_and_runs_build_before_discovery(self):
        """Preserve the documented --all --rebuild one-liner."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-all-rebuild-main-") as tmp:
            temporary = Path(tmp)
            executable = temporary / "amps"
            executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            executable.chmod(0o755)
            output = temporary / "output"
            events: list[str] = []

            def fake_rebuild(args, log_path):
                del args, log_path
                events.append("rebuild")
                return [["make", "clean"], ["make", "strict-production"]]

            def fake_all(args, output_dir, log_path):
                del args, log_path
                events.append("discover-and-run")
                report_path = output_dir / "srcsep-tests.json"
                report_path.write_text(json.dumps({
                    "schema": "srcsep-component-tests-v1", "exit_code": 0,
                    "totals": {"passed": 1, "failed": 0,
                               "skipped": 0, "errors": 0},
                    "results": [{"id": "D01", "status": "PASS",
                                 "message": "fixture", "metrics": [],
                                 "configuration": [], "artifacts": []}],
                }), encoding="utf-8")
                return 0, report_path, [[str(executable), "--test", "D01"]]

            # Keep multiple-context syntax compatible with the Python 3.8
            # interpreter supported by this source tree.
            with mock.patch.object(runner, "_rebuild_linked_application",
                                    side_effect=fake_rebuild), \
                  mock.patch.object(runner, "_run_all_tests",
                                    side_effect=fake_all), \
                  mock.patch("sys.stdout", io.StringIO()):
                exit_code = runner.main([
                    "--amps", str(executable), "--all", "--rebuild",
                    "--amps-source", str(temporary), "--make-config",
                    str(temporary / "Makefile.conf"), "--plot", "none",
                    "--output-dir", str(output)])

            self.assertEqual(exit_code, 0)
            self.assertEqual(events, ["rebuild", "discover-and-run"])
            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["command"][0], ["make", "clean"])
            self.assertEqual(manifest["command"][-1][-2:], ["--test", "D01"])

    def test_isolated_all_continues_and_lists_failed_and_error_tests(self):
        """Retain GOOD01 after BAD02 errors and FAIL03 reports failure.

        A fake process layer keeps this unit test independent of AMPS while
        exercising the real per-ID paths, synthetic ERROR record, aggregate
        JSON/JUnit generation, printed results, and actionable final lists.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-isolated-all-") as tmp:
            temporary = Path(tmp)
            executable = temporary / "amps"
            executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            executable.chmod(0o755)
            output = temporary / "output"
            output.mkdir()
            arguments = SimpleNamespace(
                amps=str(executable), model_args=[], mpi_np=None,
                mpiexec="mpiexec", timeout=None)
            called: list[list[str]] = []

            def fake_run(command, cwd, log_path, environment=None,
                         timeout_s=None):
                del cwd, log_path, environment, timeout_s
                command = list(command)
                called.append(command)
                identifier = command[command.index("--test") + 1]
                if identifier == "BAD02":
                    # Model a fatal signal/abort: no JSON exists for this ID.
                    return 134
                result_status = "FAIL" if identifier == "FAIL03" else "PASS"
                result_message = ("analytical tolerance exceeded"
                                  if result_status == "FAIL" else
                                  "fixture passed")
                report_path = Path(command[command.index("--test-json") + 1])
                report_path.write_text(json.dumps({
                    "schema": "srcsep-component-tests-v1",
                    "exit_code": 1 if result_status == "FAIL" else 0,
                    "totals": {"passed": int(result_status == "PASS"),
                               "failed": int(result_status == "FAIL"),
                               "skipped": 0, "errors": 0},
                    "results": [{
                        "id": identifier, "status": result_status,
                        "message": result_message, "elapsed_seconds": 0.01,
                        "seed": 7, "configuration": [], "metrics": [],
                        "artifacts": [],
                    }],
                }), encoding="utf-8")
                return 1 if result_status == "FAIL" else 0

            stdout = io.StringIO()
            with mock.patch.object(runner, "_discover_all_ids",
                                    return_value=["BAD02", "FAIL03", "GOOD01"]), \
                  mock.patch.object(runner, "_validation_case_ids",
                                    return_value=set()), \
                  mock.patch.object(runner, "_run_streaming",
                                    side_effect=fake_run), \
                  mock.patch("sys.stdout", stdout):
                exit_code, report_path, commands = runner._run_all_tests(
                    arguments, output, output / "test-run.log")
                # main() calls the final summary only after plot/manifest output;
                # invoke the same helper here because this focused unit calls
                # the lower-level isolated campaign function directly.
                runner._print_summary(json.loads(
                    report_path.read_text(encoding="utf-8")))

            self.assertEqual(exit_code, 2)
            self.assertEqual(len(called), 3)
            self.assertEqual(commands, called)
            self.assertEqual([command[command.index("--test") + 1]
                              for command in called],
                             ["BAD02", "FAIL03", "GOOD01"])
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["totals"], {
                "passed": 1, "failed": 1, "skipped": 0, "errors": 1})
            self.assertTrue((output / "srcsep-tests.xml").is_file())
            rendered = stdout.getvalue()
            self.assertIn("RESULT BAD02: ERROR", rendered)
            self.assertIn("RESULT FAIL03: FAIL", rendered)
            self.assertIn("RESULT GOOD01: PASS", rendered)
            self.assertIn(
                "Overall test summary: TOTAL=3 PASS=1 FAIL=1 SKIP=0 ERROR=1",
                rendered)
            self.assertIn("Failed tests (1):\n  - FAIL03: analytical tolerance exceeded",
                          rendered)
            self.assertIn("Error tests (1):\n  - BAD02: process exited 134",
                          rendered)

    def test_source_suites_record_every_outer_process_outcome(self):
        """Keep Make/launcher failures visible in console, JSON, and JUnit.

        A source suite can fail before a C++ registry exists.  These six
        outcomes exercise the B04 wrapper contract without invoking a compiler:
        success, ordinary compile failure, timeout, missing launcher, fatal
        signal, and an explicit optional-prerequisite skip.
        """
        runner = _load_runner_module()
        suites = ["fixture-pass", "fixture-compile", "fixture-timeout",
                  "fixture-missing", "fixture-signal", "fixture-skip"]
        targets = {name: "target-" + name.split("fixture-", 1)[1]
                   for name in suites}
        arguments = SimpleNamespace(
            suites=suites, keep_going=True, timeout=0.01)

        def fake_run(command, cwd, log_path, environment=None,
                     timeout_s=None, captured_output=None):
            del cwd, log_path, environment, timeout_s
            target = list(command)[-1]
            if captured_output is not None and target != "target-missing":
                captured_output.append(f"{target}: diagnostic fixture\n")
                if target == "target-skip":
                    captured_output.append("SRCSEP_SUITE_RESULT=SKIP\n")
            if target == "target-timeout":
                raise runner.RunnerError("test command exceeded 0.01 seconds")
            if target == "target-missing":
                raise runner.RunnerError(
                    "cannot start test command make: executable not found")
            return {
                "target-pass": 0,
                "target-compile": 2,
                "target-signal": -11,
                "target-skip": 0,
            }[target]

        with tempfile.TemporaryDirectory(prefix="srcsep-source-results-") as tmp:
            output = Path(tmp)
            with mock.patch.dict(runner.SOURCE_SUITES, targets, clear=False), \
                  mock.patch.object(runner, "_run_streaming",
                                    side_effect=fake_run):
                exit_code, report_path, commands = runner._run_source_suites(
                    arguments, output, output / "test-run.log")

            self.assertEqual(exit_code, 2)
            self.assertEqual(len(commands), 6)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["totals"], {
                "passed": 1, "failed": 1, "skipped": 1, "errors": 3})
            by_id = {item["id"]: item for item in report["results"]}
            expected_status = {
                "SUITE-FIXTURE-PASS": "PASS",
                "SUITE-FIXTURE-COMPILE": "FAIL",
                "SUITE-FIXTURE-TIMEOUT": "ERROR",
                "SUITE-FIXTURE-MISSING": "ERROR",
                "SUITE-FIXTURE-SIGNAL": "ERROR",
                "SUITE-FIXTURE-SKIP": "SKIP",
            }
            self.assertEqual({key: value["status"] for key, value in by_id.items()},
                             expected_status)
            for identifier, record in by_id.items():
                self.assertIn("command", record, identifier)
                self.assertIn("diagnostic_excerpt", record, identifier)
                self.assertIn("started_utc", record, identifier)
                self.assertIn("ended_utc", record, identifier)
                self.assertGreaterEqual(record["elapsed_seconds"], 0.0)
            self.assertIn("executable not found",
                          by_id["SUITE-FIXTURE-MISSING"]["diagnostic_excerpt"])
            self.assertIn("SIGSEGV", by_id["SUITE-FIXTURE-SIGNAL"]["message"])

            junit = (output / "srcsep-tests.xml").read_text(encoding="utf-8")
            self.assertIn('tests="6" failures="1" errors="3" skipped="1"',
                          junit)
            for identifier in expected_status:
                self.assertIn(f'name="{identifier}"', junit)
            self.assertIn("command: make target-compile", junit)
            self.assertIn("target-compile: diagnostic fixture", junit)

            stdout = io.StringIO()
            with mock.patch("sys.stdout", stdout):
                for record in report["results"]:
                    runner._print_result(record)
                runner._print_summary(report)
            rendered = stdout.getvalue()
            for identifier, status in expected_status.items():
                self.assertIn(f"RESULT {identifier}: {status}", rendered)
            self.assertIn(
                "Overall test summary: TOTAL=6 PASS=1 FAIL=1 SKIP=1 ERROR=3",
                rendered)

    def test_real_compile_error_is_a_structured_suite_failure(self):
        """Compile invalid C++ and retain the compiler diagnostic in reports."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-real-compile-fail-") as tmp:
            temporary = Path(tmp)
            (temporary / "broken.cpp").write_text(
                "int main( { return 0; }\n", encoding="utf-8")
            (temporary / "makefile").write_text(
                "target-compile:\n"
                "\t$(CXX) -std=c++17 -fsyntax-only broken.cpp\n",
                encoding="utf-8")
            output = temporary / "output"
            output.mkdir()
            arguments = SimpleNamespace(
                suites=["fixture-real-compile"], keep_going=True,
                timeout=30.0)
            with mock.patch.dict(
                    runner.SOURCE_SUITES,
                    {"fixture-real-compile": "target-compile"}, clear=False), \
                  mock.patch.object(runner, "ROOT", temporary):
                exit_code, report_path, _ = runner._run_source_suites(
                    arguments, output, output / "test-run.log")
            report = json.loads(report_path.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 1)
        self.assertEqual(report["totals"], {
            "passed": 0, "failed": 1, "skipped": 0, "errors": 0})
        record = report["results"][0]
        self.assertEqual(record["id"], "SUITE-FIXTURE-REAL-COMPILE")
        self.assertEqual(record["status"], "FAIL")
        self.assertIn("error:", record["diagnostic_excerpt"].lower())
        self.assertEqual(record["command"], "make target-compile")

    def test_source_suite_fail_fast_records_later_suites_as_skipped(self):
        """A fail-fast run must account for choices it deliberately blocks."""
        runner = _load_runner_module()
        arguments = SimpleNamespace(
            suites=["fixture-first", "fixture-second", "fixture-third"],
            keep_going=False, timeout=None)
        mapping = {
            "fixture-first": "target-first",
            "fixture-second": "target-second",
            "fixture-third": "target-third",
        }
        called = []

        def fail_first(command, cwd, log_path, environment=None,
                       timeout_s=None, captured_output=None):
            del cwd, log_path, environment, timeout_s
            called.append(list(command))
            if captured_output is not None:
                captured_output.append("compiler failed\n")
            return 2

        with tempfile.TemporaryDirectory(prefix="srcsep-source-blocked-") as tmp:
            output = Path(tmp)
            with mock.patch.dict(runner.SOURCE_SUITES, mapping, clear=False), \
                  mock.patch.object(runner, "_run_streaming",
                                    side_effect=fail_first):
                exit_code, report_path, _ = runner._run_source_suites(
                    arguments, output, output / "test-run.log")
            report = json.loads(report_path.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 1)
        self.assertEqual(len(called), 1)
        self.assertEqual(report["totals"], {
            "passed": 0, "failed": 1, "skipped": 2, "errors": 0})
        blocked = [item for item in report["results"]
                   if item["status"] == "SKIP"]
        self.assertTrue(all("SUITE-FIXTURE-FIRST" in item["message"]
                            for item in blocked))

    def test_xm02_publication_plot_label_names_source_and_exact_figure(self):
        """XM02 retains its detached-figure publication attribution."""
        runner = _load_cross_model_runner()
        case = json.loads((ROOT / "validation" / "cases" / "XM02" /
                           "input.json").read_text(encoding="utf-8"))
        label = runner._publication_plot_label(case)
        self.assertIn("Reference:", label)
        self.assertIn(case["reference"]["plot_citation"], label)
        self.assertIn("Figure 7", label)

    def test_xm03_plot_is_individual_legended_and_source_free(self):
        """Each XM03 time is a standalone publication panel without citation text."""
        runner = _load_cross_model_runner()
        case = json.loads((ROOT / "validation" / "cases" / "XM03" /
                           "input.json").read_text(encoding="utf-8"))
        model = []
        reference = []
        instrument_points = (
            ("ACE/EPAM", 0.3, 0.2, 0.5, 30.0),
            ("GOES-13/EPEAD", 3.0, 2.0, 5.0, 3.0),
            ("SOHO/ERNE", 30.0, 20.0, 50.0, 0.3),
        )
        for elapsed, time_factor in ((4.0, 1.0), (12.0, 2.0), (36.0, 0.5)):
            for energy, intensity in ((0.1, 100.0), (1.0, 10.0),
                                      (10.0, 1.0), (100.0, 0.1)):
                model.append({
                    "elapsed_hours": str(elapsed),
                    "energy_mev": str(energy),
                    "relative_differential_intensity": str(time_factor * intensity),
                })
            for instrument, energy, low, high, intensity in instrument_points:
                reference.append({
                    "elapsed_hours": str(elapsed),
                    "instrument": instrument,
                    "effective_energy_mev": str(energy),
                    "energy_low_mev": str(low),
                    "energy_high_mev": str(high),
                    "differential_intensity_pfu_per_mev": str(
                        time_factor * intensity),
                })

        captured = {}

        def capture(figure, output, stem, formats, **kwargs):
            axis = figure.axes[0]
            legend = axis.get_legend()
            captured[stem] = {
                "title": axis.get_title(),
                "legend": [text.get_text() for text in legend.get_texts()],
                "figure_text": " ".join(text.get_text() for text in figure.texts),
            }
            return [Path(output) / f"{stem}.{extension}" for extension in formats]

        with mock.patch.object(runner, "_save_figure", side_effect=capture):
            paths = runner._xm03_plot(
                case, Path("unused"), model, reference, 1.0, ("png", "eps"))

        self.assertEqual(len(paths), 6)
        self.assertEqual(set(captured), {
            "XM03_earth_observation_comparison_04h",
            "XM03_earth_observation_comparison_12h",
            "XM03_earth_observation_comparison_36h",
        })
        for record in captured.values():
            self.assertIn("after CME launch", record["title"])
            self.assertEqual(set(record["legend"]), {
                "srcSEP model",
                "ACE/EPAM observations",
                "GOES-13/EPEAD observations",
                "SOHO/ERNE observations",
            })
            visible_text = (record["title"] + " " + record["figure_text"]).lower()
            for forbidden in ("liu", "doi", "figure 12", "reference:", "source"):
                self.assertNotIn(forbidden, visible_text)

        with mock.patch.object(runner, "_save_figure", side_effect=capture) as save:
            outputs = runner._xm03_plot(case, Path("unused"), model, reference,
                                        1.0, ("png", "eps"), case_id="OV01")
        self.assertEqual([path.name for path in outputs], [
            "OV01_earth_observation_comparison.png",
            "OV01_earth_observation_comparison.eps",
        ])
        self.assertEqual(len(save.call_args.args[0].axes), 3)

    def test_xm03_uses_one_global_amplitude_for_all_observations(self):
        """Prevent accidental per-time or per-instrument normalization.

        The synthetic linked spectrum is exactly one fifth of every reference
        value. A correct global log-space nuisance fit must recover five while
        preserving zero residual and full coverage across all three times.
        """
        runner = _load_cross_model_runner()
        case = json.loads((ROOT / "validation" / "cases" / "XM03" /
                           "input.json").read_text(encoding="utf-8"))
        model = []
        reference = []
        for elapsed, time_factor in ((4.0, 1.0), (12.0, 10.0), (36.0, 100.0)):
            for energy, spectral_factor in ((1.0, 2.0), (10.0, 0.2)):
                observed = time_factor * spectral_factor
                model.append({
                    "elapsed_hours": str(elapsed),
                    "energy_mev": str(energy),
                    "relative_differential_intensity": str(observed / 5.0),
                    "effective_sample_count": "1000",
                })
                reference.append({
                    "elapsed_hours": str(elapsed),
                    "instrument": "fixture",
                    "energy_low_mev": str(0.9 * energy),
                    "energy_high_mev": str(1.1 * energy),
                    "effective_energy_mev": str(energy),
                    "differential_intensity_pfu_per_mev": str(observed),
                })
        metrics, scale, comparison = runner._xm03_score(case, model, reference)
        self.assertAlmostEqual(scale, 5.0, places=12)
        by_name = {item["name"]: item["value"] for item in metrics}
        self.assertEqual(by_name["observation_point_coverage"], 1.0)
        self.assertAlmostEqual(by_name["global_log10_intensity_rmse"], 0.0)
        self.assertAlmostEqual(by_name["log10_intensity_correlation"], 1.0)
        for row in comparison:
            self.assertAlmostEqual(
                float(row["model_scaled_pfu_per_mev"]),
                float(row["differential_intensity_pfu_per_mev"]), places=12)

    def test_xm02_and_xm03_duration_keys_match_their_registered_inputs(self):
        """Exercise argument construction before a costly linked execution.

        XM02 retains the historical ``duration_hours`` field. XM03 names its
        06:00-UTC origin explicitly because its spectral clock begins at the
        later CME launch. Constructing both native vectors here detects a key
        rename applied to the wrong case before AMPS is started.
        """
        runner = _load_cross_model_runner()
        xm02_path = ROOT / "validation" / "cases" / "XM02" / "input.json"
        xm03_path = ROOT / "validation" / "cases" / "XM03" / "input.json"
        xm02 = json.loads(xm02_path.read_text(encoding="utf-8"))
        xm03 = json.loads(xm03_path.read_text(encoding="utf-8"))
        xm02_arguments = runner._xm02_arguments(xm02)
        xm03_arguments = runner._xm03_arguments(xm03, xm03_path)
        self.assertEqual(
            xm02_arguments[xm02_arguments.index("--duration-s") + 1],
            "158400.0")
        self.assertEqual(
            xm03_arguments[xm03_arguments.index("--duration-s") + 1],
            "158400.0")

    def test_ov01_alias_uses_the_canonical_xm03_source(self):
        runner = _load_cross_model_runner()
        cases = ROOT / "validation" / "cases"
        ov01_path = cases / "OV01" / "input.json"
        case = json.loads(ov01_path.read_text(encoding="utf-8"))
        args = runner._xm03_arguments(case, ov01_path)
        self.assertEqual(Path(args[args.index("--source-history-csv") + 1]),
                         cases / "XM03" / "input" / "earth_shock_thermal_source.csv")
        case["physics"]["source_history_case"] = "OV01"
        with self.assertRaisesRegex(ValueError, "registered OV01"):
            runner._xm03_arguments(case, ov01_path)

    def test_existing_report_produces_png_eps_and_manifests(self):
        with tempfile.TemporaryDirectory(prefix="srcsep-python-runner-") as tmp:
            temporary = Path(tmp)
            solution_csv = temporary / "FTED08_solution.csv"
            # Conventional column names exercise the pointwise overlay path.
            # Values are fixed fixture evidence for y=x^2 at three coordinates.
            solution_csv.write_text(
                "x,numerical,analytical\n"
                "0,0,0\n"
                "1,1.02,1\n"
                "2,3.96,4\n",
                encoding="utf-8")
            report = {
                "schema": "srcsep-component-tests-v1",
                "exit_code": 0,
                "totals": {"passed": 2, "failed": 0,
                           "skipped": 0, "errors": 0},
                "results": [
                    {
                        "id": "PARK01",
                        "status": "PASS",
                        "message": "controlled analytical moment",
                        "elapsed_seconds": 0.01,
                        "seed": 1234,
                        "configuration": ["fixture=true"],
                        "metrics": [{
                            "name": "relative_error",
                            "value": 0.01,
                            "tolerance": 0.05,
                            "comparison": "relative_error <= tolerance",
                            "units": "1",
                        }],
                        "artifacts": [],
                    },
                    {
                        "id": "FTED08",
                        "status": "PASS",
                        "message": "controlled solution series",
                        "elapsed_seconds": 0.01,
                        "seed": 5678,
                        "configuration": ["fixture=true"],
                        "metrics": [],
                        "artifacts": [solution_csv.name],
                    },
                ],
            }
            report_path = temporary / "fixture-results.json"
            report_path.write_text(json.dumps(report), encoding="utf-8")
            output = temporary / "output"

            completed = subprocess.run(
                [sys.executable, str(RUNNER), "--from-json", str(report_path),
                 "--output-dir", str(output), "--formats", "png,eps"],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 0, completed.stdout)
            # The actionable lists intentionally close the transcript, after
            # plot generation and the output-directory announcement. A clean
            # report still names both categories and marks them empty.
            self.assertGreater(completed.stdout.rfind("Failed tests (0):"),
                               completed.stdout.rfind("Results:"))
            self.assertGreater(completed.stdout.rfind("Error tests (0):"),
                               completed.stdout.rfind("Failed tests (0):"))
            self.assertTrue(completed.stdout.rstrip().endswith(
                "Error tests (0):\n  none"), completed.stdout)

            manifest = json.loads(
                (output / "analytical_plot_manifest.json").read_text(
                    encoding="utf-8"))
            by_id = {entry["test_id"]: entry for entry in manifest["plots"]}
            self.assertEqual(by_id["PARK01"]["kind"], "metric-acceptance")
            self.assertEqual(by_id["FTED08"]["kind"], "solution-series")
            for test_id in ("PARK01", "FTED08"):
                for suffix in ("png", "eps"):
                    figure = output / "plots" / f"{test_id}_comparison.{suffix}"
                    self.assertTrue(figure.is_file(), figure)
                    self.assertGreater(figure.stat().st_size, 0)

            run_manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(run_manifest["exit_code"], 0)
            self.assertTrue(run_manifest["report_sha256"])

    def test_validation_mode_forwards_the_selected_linked_executable(self):
        """Prevent CV01 from regressing to a Python-compiled model driver."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=["CV01"],
                case_input=None,
                timeout=42.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[:4], [
                sys.executable,
                str(ROOT / "validation" / "run_case.py"),
                "--amps",
                "/opt/amps/bin/srcsep-amps",
            ])
            self.assertIn("CV01", command)
            self.assertIn("--timeout", command)
            run.assert_called_once()

    def test_cv02_cv05_share_the_linked_application_selection_path(self):
        """Require every new controlled case to use the requested AMPS binary."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=["CV02", "CV03", "CV04", "CV05"],
                case_input=None,
                timeout=120.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in ("CV02", "CV03", "CV04", "CV05"):
                self.assertIn(case_id, command)
            self.assertNotIn("controlled_transport_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_cv06_cv12_share_the_linked_application_selection_path(self):
        """Keep advanced physics cases attached to the selected application.

        The case-local source harness is useful for compile testing, but it is
        not scientific evidence. This assertion prevents future convenience
        changes from bypassing validation/run_case.py or replacing ``--amps``
        with advanced_validation_models.cpp.
        """
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-validation-command-") as tmp:
            output = Path(tmp)
            selected = [f"CV{index:02d}" for index in range(6, 13)]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps",
                validation_all=False,
                validation_cases=selected,
                case_input=None,
                timeout=600.0,
            )
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            self.assertNotIn("advanced_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_iv01_iv06_use_the_selected_linked_application(self):
        """Integrated cases must never be relabelled standalone evidence."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-iv-command-") as tmp:
            output = Path(tmp)
            selected = [f"IV{index:02d}" for index in range(1, 7)]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps", validation_all=False,
                validation_cases=selected, case_input=None, timeout=600.0)
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            self.assertNotIn("integrated_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_xm01_xm03_use_the_selected_linked_application(self):
        """Cross-model orchestration must retain the production binary gate."""
        runner = _load_runner_module()
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-command-") as tmp:
            output = Path(tmp)
            selected = ["XM01", "XM02", "XM03"]
            arguments = SimpleNamespace(
                amps="/opt/amps/bin/srcsep-amps", validation_all=False,
                validation_cases=selected, case_input=None, timeout=900.0)
            with mock.patch.object(runner, "_run_streaming", return_value=0) as run:
                _, _, command = runner._run_validation_cases(
                    arguments, output, output / "runner.log")
            self.assertEqual(command[2:4], ["--amps", "/opt/amps/bin/srcsep-amps"])
            for case_id in selected:
                self.assertIn(case_id, command)
            # XM02/XM03 resolve their only input from the registry.  The
            # top-level runner must not manufacture an --input override while
            # assembling a multi-case command.
            self.assertNotIn("--input", command)
            self.assertNotIn("cross_model_validation_models.cpp", " ".join(command))
            run.assert_called_once()

    def test_xm02_xm03_reject_command_line_input_overrides(self):
        """Keep each publication comparison tied to its registered input.

        This check occurs before executable validation, so a deliberately
        missing binary can be used to prove that input-policy diagnostics are
        deterministic and do not depend on the local AMPS installation.
        """
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-fixed-input-") as tmp:
            temporary = Path(tmp)
            alternate = temporary / "alternate.json"
            alternate.write_text("{}\n", encoding="utf-8")
            for case_id in ("XM02", "XM03"):
                completed = subprocess.run(
                    [sys.executable, str(RUNNER), "--amps",
                     str(temporary / "missing-amps"), "--validation-case", case_id,
                     "--case-input", str(alternate), "--output-dir",
                     str(temporary / case_id)],
                    cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, check=False)
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn(
                    f"{case_id} uses its single registered publication-derived input",
                    completed.stdout)

    def test_direct_case_runner_enforces_fixed_xm_input(self):
        """Apply the same no-override rule below the convenience front end."""
        with tempfile.TemporaryDirectory(prefix="srcsep-xm-direct-input-") as tmp:
            temporary = Path(tmp)
            alternate = temporary / "alternate.json"
            alternate.write_text("{}\n", encoding="utf-8")
            completed = subprocess.run(
                [sys.executable, str(CASE_RUNNER), "--amps",
                 str(temporary / "missing-amps"), "--case", "XM02", "--input",
                 str(alternate), "--output-dir", str(temporary / "output")],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn(
                "XM02 uses its single registered publication-derived input",
                completed.stdout)

    def test_validation_case_rejects_a_missing_linked_executable(self):
        """Never replace unavailable application evidence with a local build."""
        with tempfile.TemporaryDirectory(prefix="srcsep-missing-amps-") as tmp:
            missing = Path(tmp) / "amps-does-not-exist"
            completed = subprocess.run(
                [sys.executable, str(CASE_RUNNER), "--amps", str(missing),
                 "--case", "CV01", "--output-dir", str(Path(tmp) / "out")],
                cwd=str(ROOT), text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=False)
            self.assertEqual(completed.returncode, 2, completed.stdout)
            self.assertIn("linked srcSEP/AMPS executable is missing", completed.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
