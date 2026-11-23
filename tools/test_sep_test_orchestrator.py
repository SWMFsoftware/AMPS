#!/usr/bin/env python3
"""Self-test for tools/sep_test_orchestrator.py (run_all_test.sh backend).

Each result reader is exercised on a small synthetic fixture written in the
exact format its runner produces (JSON report fields, shared-log markers,
unittest verbose text, ``[ID] STATUS`` lines).  Expected counts, statuses,
messages and extracted log contents are written into the fixtures and
asserted literally; nothing is computed by the code under test.  The
step-level guards (missing report, nonzero exit without a failing record),
the status vocabulary, de-duplication into "unique tests", and the summary's
exit status and report files are covered as well.

Run:  python3 tools/test_sep_test_orchestrator.py      (verbosity 2)
Needs no AMPS build; it is registered as the ``orchestrator-selftest`` step.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("sep_test_orchestrator",
                                               _HERE / "sep_test_orchestrator.py")
orch = importlib.util.module_from_spec(_SPEC)
# dataclasses resolves string annotations through sys.modules, so the module
# must be registered before it executes.
sys.modules[_SPEC.name] = orch
_SPEC.loader.exec_module(orch)


class Fixture(unittest.TestCase):
    """Temporary output directory plus a context factory for one step."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def ctx(self, name="step", log_text="", exit_code=0, kind="test",
            runner="runner.py", namespace=""):
        step = orch.Step(name, "fixture step", [], kind=kind, runner=runner,
                         namespace=namespace)
        log = self.out / f"01-{name}.log"
        log.write_text(log_text, encoding="utf-8")
        return orch.StepContext(step, self.out, log, exit_code, time.time() - 5,
                                self.out / "failed-tests" / name)

    def write_json(self, relative, payload):
        path = self.out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path


class SharedInitializerAudit(Fixture):
    """The source guard follows actual shared ownership and rejects corruption."""
    def test_analytic_history_guard_and_early_return_are_required(self):
        root=_HERE.parent
        spec=importlib.util.spec_from_file_location('sep_source_contract_audit',
            root/'srcSEP/test/check_stage1_stage2_contracts.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        source=(root/'srcSEP/main_lib.cpp').read_text()
        signature='SEP::Background::AnalyticRadial::InitializeFixedNativeHistory()'
        body=module.function_body(module.code_without_comments(source),signature)
        module.check_analytic_history_owner(signature+body)
        # Corrupt the real guarded function independently. A wrong provider or
        # lost no-op must fail, even though other files still contain Analytic.
        wrong=body.replace('ConfiguredProvider()==BG::Provider::Analytic?1:0',
                           'ConfiguredProvider()==BG::Provider::Swmf?1:0')
        self.assertNotEqual(wrong,body)
        with self.assertRaises(AssertionError):module.check_analytic_history_owner(signature+wrong)
        wrong=body.replace('if(!maximumEnabled)return T::Status::Ok();','')
        self.assertNotEqual(wrong,body)
        with self.assertRaises(AssertionError):module.check_analytic_history_owner(signature+wrong)


class MaterialBackgroundReports(Fixture):
    """Software report controls; these fixtures make no physical claim."""
    def payload(self):
        return dict(schema='srcsep-material-background-v1',results=[
            dict(test_id=f'MBC{i:02d}',status='PASS') for i in range(1,21)])

    def test_complete_catalog_and_preserved_failures(self):
        payload=self.payload()
        payload['results'][3].update(status='FAIL',message='invalid supplied derivative')
        payload['results'][7].update(status='ERROR',message='compiler unavailable')
        self.write_json('summary.json',payload)
        rows=orch.reader_material_background(self.out)(self.ctx(exit_code=1))
        self.assertEqual(len(rows),20)
        self.assertEqual((rows[3].status,rows[7].status),('FAIL','ERROR'))

    def test_missing_foreign_duplicate_and_schema_are_errors(self):
        for change in ('missing','foreign','duplicate','schema','status'):
            payload=self.payload()
            if change=='missing':payload['results'].pop()
            elif change=='foreign':payload['results'][0]['test_id']='SBACK01'
            elif change=='duplicate':payload['results'][1]['test_id']='MBC01'
            elif change=='schema':payload['schema']='legacy-snapshot'
            else:payload['results'][0]['status']='SUCCESS'
            self.write_json('summary.json',payload)
            with self.assertRaises(orch.ResultsUnavailable):
                orch.reader_material_background(self.out)(self.ctx())

    def test_runner_build_error_keeps_every_case_and_actual_receipt(self):
        # A missing compiler fails the real maintained build immediately. The
        # structured report must retain every contract ID as ERROR and the
        # attempted make command/actual nonzero exit, never reuse old PASSes.
        root=_HERE.parent;out=self.out/'missing-compiler'
        result=subprocess.run([sys.executable,'-B',str(root/'srcSEP/test/run_material_background_tests.py'),
            '--output-dir',str(out)],cwd=root,env=dict(os.environ,CXX='/nonexistent/material-contract-compiler'),
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
        self.assertEqual(result.returncode,1)
        payload=json.loads((out/'summary.json').read_text())
        self.assertEqual([r['test_id'] for r in payload['results']],[f'MBC{i:02d}' for i in range(1,21)])
        self.assertTrue(all(r['status']=='ERROR' for r in payload['results']))
        self.assertEqual(payload['commands'][0]['argv'],['make','-j16','lib'])
        self.assertNotEqual(payload['commands'][0]['exit_code'],0)

    def test_existing_evidence_is_immutable(self):
        # Reusing an output directory is an invocation error. It must not
        # erase either successful or failed prior scientific evidence.
        root=_HERE.parent;out=self.out/'retained';out.mkdir()
        witness=out/'summary.json';witness.write_text('retained previous evidence')
        result=subprocess.run([sys.executable,'-B',str(root/'srcSEP/test/run_material_background_tests.py'),
            '--output-dir',str(out)],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
        self.assertEqual(result.returncode,1)
        self.assertEqual(witness.read_text(),'retained previous evidence')

    def test_both_registries_include_material_owner(self):
        root=_HERE.parent
        manifest=json.loads((root/'src/test/manifest.json').read_text())
        entry=next(t for t in manifest['tests'] if t['id']=='CORE-SEP-MATERIAL-BACKGROUND-001')
        self.assertEqual(entry['expected_cases'],[f'MBC{i:02d}' for i in range(1,21)])
        self.assertIn('material-background',[s.name for s in orch.srcsep_steps(self.out,argparse.Namespace())])


class BackendAndApplicationReports(Fixture):
    def application(self):
        return {'schema':'srcsep-stage6-application-v1','results':[
            dict(test_id=identity,mode=mode,status='PASS',physical_steps=400,
                 physical_dt_s=300.0,acceptance_relative_budget_bound=1e-10)
            for identity,mode in zip(('SAPP01','SAPP02','SAPP03'),('parker','fte-dmumu','fte-mfp'))]}

    def test_full_application_complete_catalog(self):
        self.write_json('verified-results.json',self.application())
        self.assertEqual(len(orch.reader_standalone_application(self.out)(self.ctx())),3)

    def test_full_application_rejects_partial_or_changed_window(self):
        for change in ('missing','mode','dt','steps','tolerance'):
            payload=self.application()
            if change=='missing':payload['results'].pop()
            elif change=='mode':payload['results'][0]['mode']='fte-mfp'
            elif change=='dt':payload['results'][0]['physical_dt_s']=299
            elif change=='steps':payload['results'][0]['physical_steps']=40
            else:payload['results'][0]['acceptance_relative_budget_bound']=1e-9
            self.write_json('verified-results.json',payload)
            with self.assertRaises(orch.ResultsUnavailable):orch.reader_standalone_application(self.out)(self.ctx())

    def test_full_application_preserves_actual_fail_and_error(self):
        payload=self.application();payload['results'][0].update(status='ERROR',message='actual MPI startup error')
        payload['results'][1].update(status='FAIL',message='signed ledger does not close')
        self.write_json('verified-results.json',payload)
        self.assertEqual([r.status for r in orch.reader_standalone_application(self.out)(self.ctx(exit_code=1))],['ERROR','FAIL','PASS'])

    def test_backend_rejection_catalog_and_scope_are_exact(self):
        payload={'schema':'srcsep-backend-admission-v1','native':True,'results':[
            dict(test_id=identity,status='PASS') for identity in ('SBACK06','SBACK07','SBACK08')]}
        self.write_json('summary.json',payload)
        self.assertEqual(len(orch.reader_backend_admission(self.out,native=True)(self.ctx())),3)
        with self.assertRaises(orch.ResultsUnavailable):orch.reader_backend_admission(self.out)(self.ctx())
        payload['results'].pop();self.write_json('summary.json',payload)
        with self.assertRaises(orch.ResultsUnavailable):orch.reader_backend_admission(self.out,native=True)(self.ctx())

    def test_stage8_9_shared_catalogs_are_registered(self):
        root=_HERE.parent;manifest=json.loads((root/'src/test/manifest.json').read_text())
        expected={'CORE-SEP-BACKEND-ADMISSION-001':['SBACK01','SBACK02','SBACK03','SBACK04','SBACK05'],
                  'CORE-SEP-BACKEND-NATIVE-001':['SBACK06','SBACK07','SBACK08'],
                  'CORE-SEP-FULL-APPLICATION-001':['SAPP01','SAPP02','SAPP03']}
        for name,ids in expected.items():
            self.assertEqual(next(t for t in manifest['tests'] if t['id']==name)['expected_cases'],ids)


class StandaloneIntegrationReports(Fixture):
    def test_shared_registry_matches_complete_numerical_catalog(self):
        # A newly added numerical case must remain visible through BOTH
        # validation entry points. The real shared run caught a stale wrapper
        # list even though all kernels passed; accepting a subset would hide
        # failures. Compare the actual frozen source catalogs, not counts.
        root=_HERE.parent
        source=ast.parse((root/'srcSEP/test/run_standalone_integration_tests.py').read_text())
        ids=[]
        for node in source.body:
            if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='IDS' for t in node.targets):
                ids.extend(ast.literal_eval(node.value))
            if isinstance(node,ast.AugAssign) and isinstance(node.target,ast.Name) and node.target.id=='IDS':
                ids.extend(ast.literal_eval(node.value))
        wrapper=ast.parse((root/'src/test/tests/models/cme/sep_standalone/run_test.py').read_text())
        assignment=next(node for node in ast.walk(wrapper) if isinstance(node,ast.Assign) and
            any(isinstance(t,ast.Name) and t.id=='expected' for t in node.targets))
        self.assertEqual(ast.literal_eval(assignment.value.orelse),ids)
        manifest=json.loads((root/'src/test/manifest.json').read_text())
        entry=next(test for test in manifest['tests'] if test['id']=='CORE-SEP-STANDALONE-INTEGRATION-001')
        self.assertEqual(entry['expected_cases'],ids)

    def test_incomplete_successful_report_rejects(self):
        self.write_json("native/summary.json", {"schema":"srcsep-native-standalone-integration-v1",
            "results":[{"test_id":"SINT08","status":"PASS"}]})
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_standalone_integration(self.out/"native",native=True)(self.ctx())

    def test_complete_report_preserves_infrastructure_error(self):
        ids=["SINT08","SINT09","SINT13N","SINT14N","SINT15N","SINT16N",
             "SINT25N","SINT26N","SINT27N","SINT28N","SINT29N","SINT30N","SINT31N","SINT32N","MCP08"]
        self.write_json("native/summary.json", {"schema":"srcsep-native-standalone-integration-v1",
            "results":[{"test_id":identity,"status":"ERROR","message":"MPI launch failed"} for identity in ids]})
        records=orch.reader_standalone_integration(self.out/"native",native=True)(self.ctx(exit_code=1))
        self.assertEqual([record.status for record in records],["ERROR"]*len(ids))

    def test_native_merger_remains_visible_in_both_registries(self):
        # Population reduction must preserve immutable birth energy in the
        # production adapter; a successful restart alone does not test a merge.
        root=_HERE.parent
        wrapper=ast.parse((root/'src/test/tests/models/cme/sep_standalone/run_test.py').read_text())
        assignment=next(node for node in ast.walk(wrapper) if isinstance(node,ast.Assign) and
            any(isinstance(t,ast.Name) and t.id=='expected' for t in node.targets))
        ids=ast.literal_eval(assignment.value.body)
        self.assertEqual(ids[-1],'MCP08')
        manifest=json.loads((root/'src/test/manifest.json').read_text())
        entry=next(t for t in manifest['tests'] if t['id']=='CORE-SEP-STANDALONE-NATIVE-001')
        self.assertEqual(entry['expected_cases'],ids)

    def test_material_report_requires_all_roles_and_preserves_errors(self):
        # Manufacture report structure only. These are isolated reader negative
        # controls, never claimed as native MPI or particle-history evidence.
        rows=[]
        for index in range(1,4):
            for ranks in (1,2,4,8):
                for role in ('continuous','write','resume'):
                    rows.append(dict(id=f'MCP{index:02d}',ranks=ranks,role=role,status='PASS',cwd='/fixture'))
                rows.append(dict(id=f'MCP{index+8:02d}',ranks=ranks,role='topology',status='PASS',cwd='/fixture'))
            for writer,reader in ((1,8),(8,1)):
                rows.append(dict(id=f'MCP{index+3:02d}',ranks=reader,writer_ranks=writer,role='repartition',status='PASS',cwd='/fixture'))
        for ranks in (1,8):
            rows.append(dict(id='MCP07',ranks=ranks,role='reject',status='PASS',cwd='/fixture'))
            rows.append(dict(id='MCP12',ranks=ranks,role='late-reject',status='PASS',cwd='/fixture'))
        payload=dict(schema='srcsep-material-checkpoint-v1',results=rows)
        self.write_json('material/summary.json',payload)
        self.assertEqual(len(orch.material_checkpoint_rows(self.out/'material')),58)
        for defect in ('missing','duplicate','wrong-donor','unknown-status'):
            corrupt=json.loads(json.dumps(payload))
            if defect=='missing':corrupt['results'].pop()
            elif defect=='duplicate':corrupt['results'][-1]=corrupt['results'][0]
            elif defect=='wrong-donor':
                next(r for r in corrupt['results'] if r['role']=='repartition')['writer_ranks']=2
            else:corrupt['results'][0]['status']='success-but-not-verified'
            self.write_json('material/summary.json',corrupt)
            with self.assertRaises(orch.ResultsUnavailable):
                orch.material_checkpoint_rows(self.out/'material')
        payload['results'][0]['status']='ERROR'
        payload['results'][0]['message']='actual MPI startup failed'
        self.write_json('material/summary.json',payload)
        records=orch.reader_material_checkpoint(self.out/'material')(self.ctx(exit_code=1))
        self.assertEqual(records[0].status,'ERROR')

    def test_stage7_requires_combined_thread_roles_and_individual_failures(self):
        # Structure-only controls. Omitting a resume or worker configuration
        # must not turn passing smoke output into restart qualification.
        rows=[]
        for identity in ('CP05C01','CP05C02','CP05C03','CP05C04','CP05C05',
                         'CP05O01','CP05O02','CP05O03','CP05O04','CP05O05'):
            rows.append(dict(id=identity,role='routine',threads=0,ranks=2,status='PASS',cwd='/fixture/'+identity))
        for identity in ('CP05F01','CP05F02','CP05F03','CP05OF01','CP05OF02','CP05OF03'):
            for role in ('write','resume'):
                rows.append(dict(id=identity,role=role,threads=0,ranks=2,status='PASS',cwd='/fixture/'+identity+role))
        rows.append(dict(id='MCP13',role='missing-fixture',threads=0,ranks=2,status='PASS',cwd='/fixture/guard'))
        for identity in ('TCR01','TCR02','TCR03'):
            for threads in (2,4):
                for ranks in ((1,8) if threads==4 and identity=='TCR02' else (1,)):
                    for role in ('continuous','write','resume'):
                        rows.append(dict(id=identity,role=role,threads=threads,ranks=ranks,status='PASS',cwd='/fixture/'+identity+role+str(threads)+str(ranks)))
        payload=dict(schema='srcsep-stage7-checkpoint-v1',fixtures_only=False,results=rows)
        self.write_json('stage7/summary.json',payload)
        self.assertEqual(len(orch.stage7_checkpoint_rows(self.out/'stage7')),44)
        for defect in ('missing-resume','wrong-workers','fixture-only'):
            altered=json.loads(json.dumps(payload))
            if defect=='missing-resume':altered['results'].pop()
            elif defect=='wrong-workers':altered['results'][-1]['threads']=1
            else:altered['fixtures_only']=True
            self.write_json('stage7/summary.json',altered)
            with self.assertRaises(orch.ResultsUnavailable):orch.stage7_checkpoint_rows(self.out/'stage7')
        payload['results'][-1]['status']='FAIL';payload['results'][-1]['message']='actual resumed history differs'
        self.write_json('stage7/summary.json',payload)
        records=orch.reader_stage7_checkpoint(self.out/'stage7')(self.ctx(exit_code=1))
        self.assertEqual(records[-1].status,'FAIL')
        steps={step.name:step for step in orch.srcsep_steps(self.out,argparse.Namespace())}
        self.assertIsNotNone(steps['stage7-checkpoint-native'].reader)


class StatusVocabulary(unittest.TestCase):
    def test_runner_spellings_map_to_four_statuses(self):
        cases = {"PASS": "PASS", "passed": "PASS", "ok": "PASS", "FAIL": "FAIL",
                 "failed": "FAIL", "SKIP": "SKIP", "skipped": "SKIP",
                 "ERROR": "ERROR", "errors": "ERROR", "MISSING": "ERROR",
                 "weird": "ERROR"}
        for raw, expected in cases.items():
            self.assertEqual(orch.normalize_status(raw), expected, raw)


class SrcSepRunTests(Fixture):
    def test_failed_test_section_is_copied_from_shared_log(self):
        report = self.out / "run-tests"
        self.write_json("run-tests/srcsep-tests.json", {"results": [
            {"id": "PARK01", "status": "PASS", "message": "ok"},
            {"id": "FTE01", "status": "ERROR", "message": "process exited 0"},
            {"id": "OV02", "status": "ERROR", "message": "TypeError gating"},
            {"id": "SCAT01", "status": "SKIP", "message": "no rule"}]})
        self.write_json("run-tests/individual/OV02/OV02.json", {})
        (report / "test-run.log").write_text(
            "header\n"
            "[2026-10-10T06:17:00+00:00] RUN /amps --test PARK01 --test-json x\n"
            "park output\n"
            "[2026-10-10T06:17:04+00:00] RUN /amps --test FTE01 --test-json y\n"
            "fte line 1\nAMPS: exit: focused_transport_dmumu.cpp:13\n"
            "[2026-10-10T06:18:00+00:00] RUN python3 run_case.py --case OV02 --output-dir z\n"
            "Traceback: TypeError gating\n", encoding="utf-8")
        records = orch.reader_srcsep_run_tests(report)(self.ctx("run-tests"))
        by_id = {r.test_id: r for r in records}
        self.assertEqual(orch.count(records),
                         {"PASS": 1, "FAIL": 0, "SKIP": 1, "ERROR": 2, "TOTAL": 4})
        fte = Path(by_id["FTE01"].log).read_text()
        self.assertIn("# source: ", fte)
        self.assertIn("focused_transport_dmumu.cpp:13", fte)
        self.assertNotIn("park output", fte)
        self.assertNotIn("TypeError", fte)
        self.assertIn("--case OV02", Path(by_id["OV02"].log).read_text())
        self.assertTrue(by_id["OV02"].evidence[0].endswith("individual/OV02/OV02.json"))
        self.assertEqual(by_id["PARK01"].log, str(report / "test-run.log"))

    def test_missing_report_is_unavailable(self):
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_srcsep_run_tests(self.out / "absent")(self.ctx())


class SrcSep3dRunTests(Fixture):
    def test_runner_log_is_used_and_fallback_locates_line(self):
        report = self.out / "rt"
        self.write_json("rt/srcsep3d-tests.json", {"results": [
            {"test_id": "R3D01", "status": "FAIL", "message": "deck", "log": "/logs/R3D01.log"},
            {"test_id": "LAY01", "status": "ERROR", "message": "crash"}]})
        ctx = self.ctx("run-tests-all", "[R3D01] FAIL\nnoise\n[LAY01] ERROR (0.1s) crash\n")
        by_id = {r.test_id: r for r in orch.reader_srcsep3d_run_tests(report)(ctx)}
        self.assertEqual(by_id["R3D01"].log, "/logs/R3D01.log")
        self.assertEqual((by_id["LAY01"].log, by_id["LAY01"].line), (str(ctx.step_log), 3))


class PdMover(Fixture):
    def test_missing_maps_to_error_and_suite_namespaces(self):
        report = self.write_json("pd.json", {"suites": [
            {"name": "library", "results": [
                {"test_id": "PD13-HOST-GATE", "status": "PASS", "suite": "library",
                 "detail": "ok", "log": "/l/library.log", "line": 9}]},
            {"name": "stage1", "results": [
                {"test_id": "CFG3D17", "status": "MISSING", "suite": "stage1",
                 "detail": "expected test did not report", "log": "/l/stage1.log", "line": 0}]}]})
        records = orch.reader_pd_mover(report, "pd.py", {"library": "lib"})(self.ctx("pd-mover"))
        self.assertEqual([r.status for r in records], ["PASS", "ERROR"])
        self.assertEqual(records[0].namespace, "lib")
        self.assertEqual(records[1].namespace, "pd-mover/stage1")
        self.assertEqual(records[0].runner, "pd.py [library]")
        self.assertEqual((records[0].log, records[0].line), ("/l/library.log", 9))


class NativeJson(Fixture):
    def test_results_point_to_step_log_line_and_artifacts(self):
        report = self.write_json("n/native.json", {"results": [
            {"id": "SCCM3D04", "status": "FAIL", "message": "epoch", "artifacts": ["/a/x"]}]})
        ctx = self.ctx("sccm3d", "init\n[SCCM3D04] FAIL epoch\n")
        record = orch.reader_native_json(report)(ctx)[0]
        self.assertEqual((record.status, record.line), ("FAIL", 2))
        self.assertEqual(record.evidence, [str(report), "/a/x"])

    def test_required_native_layout_cannot_pass_by_skip_or_missing_cases(self):
        report = self.write_json("spectral/results.json", {"results": [
            {"id": "PNW01", "status": "PASS", "message": "bins"},
            {"id": "PNW02", "status": "SKIP", "message": "spectral datum absent"},
            {"id": "PNW03", "status": "FAIL", "message": "restore changed helicity"}]})
        records = orch.reader_required_native_json(report, ["PNW01", "PNW02", "PNW03", "PNW04"])(
            self.ctx("polarized-wave-native"))
        self.assertEqual([(r.test_id, r.status) for r in records],
                         [("PNW01", "PASS"), ("PNW02", "ERROR"), ("PNW03", "FAIL"), ("PNW04", "ERROR")])
        self.assertIn("not exercised", records[1].message)


class NativeCliCampaign(Fixture):
    def test_individual_job_logs_and_failed_postchecks_are_authoritative(self):
        report = self.write_json("cli/summary.json", {
            "schema": "srcsep-native-cli-campaign-v1", "results": [
                {"id": "NCLI01", "status": "PASS", "message": "draws", "log": "/cli/seed.log",
                 "exit_code": 0, "expected_exit_code": 0, "postchecks": [{"status": "PASS"}]},
                {"id": "NCLI04", "status": "PASS", "message": "expected rank fault", "log": "/cli/fault.log",
                 "exit_code": 1, "expected_exit_code": 1, "postchecks": [{"status": "PASS"}]},
                {"id": "NCLI05", "status": "PASS", "message": "false green", "log": "/cli/invalid.log",
                 "postchecks": [{"status": "FAIL", "message": "missing failure provenance"}]}]})
        records = orch.reader_native_cli_campaign(report)(self.ctx("native-cli-campaign"))
        self.assertEqual([(r.test_id, r.status) for r in records],
                         [("NCLI01", "PASS"), ("NCLI04", "PASS"), ("NCLI05", "ERROR")])
        self.assertEqual(records[1].log, "/cli/fault.log")
        self.assertIn("missing failure provenance", records[2].message)
        self.assertEqual(records[0].evidence, [str(report)])

    def test_runner_records_launch_and_exit_failures_before_raising(self):
        script = orch.AMPS_ROOT / "srcSEP/test/run_native_cli_campaign_tests.py"
        for name, launcher, expected in (("missing", str(self.out / "no-launcher"), "ERROR"),
                                          ("bad-exit", "/bin/false", "FAIL")):
            output = self.out / name
            process = subprocess.run([sys.executable, str(script), "--amps", "/bin/true",
                                      "--mpiexec", launcher, "--output", str(output)],
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self.assertNotEqual(process.returncode, 0)
            payload = json.loads((output / "summary.json").read_text())
            self.assertEqual(len(payload["results"]), 9)
            self.assertEqual(payload["results"][0]["status"], expected)
            self.assertEqual(payload["results"][0]["id"], "NCLI01")
            self.assertEqual(payload["totals"]["SKIP"], 8)
            self.assertTrue(Path(payload["results"][0]["log"]).is_file())
            self.assertTrue(payload["results"][0]["postchecks"])

    def test_wave_and_cli_steps_have_explicit_readers(self):
        steps = {step.name: step for step in orch.srcsep_steps(self.out, argparse.Namespace())}
        for name, target, ids in (
                ("real-wave-normalization", "test-real-wave-normalization-unit", ["RWN01", "RWN02"]),
                ("polarized-wave-state", "test-polarized-wave-state-unit", ["PWS01", "PWS02", "PWS03", "PWS04", "PWS05", "PWS06"])):
            step = steps[name]
            self.assertEqual(step.commands, [["make", "-C", "srcSEP", target]])
            log = "".join("[%s] %s explicit gate\n" % (identifier, "FAIL" if index == 0 else "PASS")
                          for index, identifier in enumerate(ids))
            records = step.reader(self.ctx(name, log, exit_code=1))
            self.assertEqual([r.test_id for r in records], ids)
            self.assertEqual(records[0].status, "FAIL")
        self.assertTrue(steps["native-cli-campaign"].reader)
        self.assertEqual(steps["native-cli-campaign"].namespace, "srcsep-cli-infrastructure")


class FreshProcessCheckpoint(Fixture):
    """Literal retained files test the evidence boundary, without native PIC.

    A summary saying PASS is deliberately held constant in negative controls.
    Each control damages one actual prerequisite so that self-consistent
    summary labels cannot hide missing processes, reports, or final bytes.
    """

    def make_evidence(self):
        root = self.out / "fresh"
        executable = self.out / "fixture-native"
        executable.write_bytes(b"literal linked fixture executable")
        rows = []
        for identifier in ("CP05F01", "CP05F02", "CP05F03"):
            roles = []
            for name in ("write", "resume"):
                child = self.write_json("fresh/%s/%s/srcsep-tests.json" % (identifier, name),
                                        {"results": [{"id": identifier, "status": "PASS"}]})
                log = child.with_name("native.log")
                log.write_text("literal child process log\n")
                roles.append({"role": name, "status": "PASS", "exit_code": 0,
                              "report": str(child), "log": str(log)})
            (root / identifier / "physics.expected").write_bytes(b"complete encoded final state\0\xff")
            (root / identifier / "physics.resumed").write_bytes(b"complete encoded final state\0\xff")
            rows.append({"id": identifier, "status": "PASS", "message": "", "roles": roles,
                         "final_artifact_equality": True})
        payload = {"schema": "srcsep-fresh-process-checkpoint-v1", "executable": str(executable),
                   "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                   "results": rows}
        report = self.write_json("fresh/summary.json", payload)
        return report, payload, executable

    def test_success_requires_all_retained_native_roles_and_bytes(self):
        report, _, executable = self.make_evidence()
        records = orch.reader_fresh_process_checkpoint(report)(self.ctx("fresh-process-checkpoint"))
        self.assertEqual([(record.test_id, record.status) for record in records],
                         [("CP05F01", "PASS"), ("CP05F02", "PASS"), ("CP05F03", "PASS")])
        self.assertIn(str(executable), records[0].evidence)
        self.assertIn("/resume/", records[0].log)
        self.assertTrue(any(path.endswith("physics.expected") for path in records[0].evidence))

    def test_pass_summary_cannot_hide_missing_skip_exit_binary_or_byte_failures(self):
        for defect in ("missing-case", "duplicate-case", "skip-case", "missing-role", "skip-role",
                       "nonzero-exit", "boolean-exit", "missing-child", "skip-child", "wrong-child-id",
                       "malformed-child", "unequal-bytes", "missing-bytes", "false-equality",
                       "changed-binary", "missing-hash"):
            with self.subTest(defect=defect):
                report, payload, executable = self.make_evidence()
                row = payload["results"][0]
                role = row["roles"][1]
                if defect == "missing-case": payload["results"].pop(0)
                elif defect == "duplicate-case": payload["results"].append(row.copy())
                elif defect == "skip-case": row["status"] = "SKIP"
                elif defect == "missing-role": row["roles"].pop()
                elif defect == "skip-role": role["status"] = "SKIP"
                elif defect == "nonzero-exit": role["exit_code"] = 17
                elif defect == "boolean-exit": role["exit_code"] = False
                elif defect == "missing-child": Path(role["report"]).unlink()
                elif defect == "skip-child": Path(role["report"]).write_text(json.dumps({"results": [{"id": "CP05F01", "status": "SKIP"}]}))
                elif defect == "wrong-child-id": Path(role["report"]).write_text(json.dumps({"results": [{"id": "CP05F02", "status": "PASS"}]}))
                elif defect == "malformed-child": Path(role["report"]).write_text("[]")
                elif defect == "unequal-bytes": (report.parent / "CP05F01/physics.resumed").write_bytes(b"different")
                elif defect == "missing-bytes": (report.parent / "CP05F01/physics.resumed").unlink()
                elif defect == "false-equality": row["final_artifact_equality"] = False
                elif defect == "changed-binary": executable.write_bytes(b"changed executable")
                elif defect == "missing-hash": payload.pop("executable_sha256")
                report.write_text(json.dumps(payload))
                records = orch.reader_fresh_process_checkpoint(report)(self.ctx("fresh-process-checkpoint"))
                first = next(record for record in records if record.test_id == "CP05F01")
                self.assertEqual(first.status, "ERROR", first.message)
                self.assertTrue(first.message)

    def test_failed_child_stays_failed_and_schema_does_not_fall_back(self):
        report, payload, _ = self.make_evidence()
        payload["results"][0]["status"] = "FAIL"
        report.write_text(json.dumps(payload))
        records = orch.reader_fresh_process_checkpoint(report)(self.ctx("fresh-process-checkpoint", exit_code=1))
        self.assertEqual(records[0].status, "FAIL")
        payload["schema"] = "unknown-version"
        report.write_text(json.dumps(payload))
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_fresh_process_checkpoint(report)(self.ctx())


class NativeIsolatedCampaigns(Fixture):
    """A green aggregate cannot replace retained process/input/trace evidence."""

    def make_host(self):
        executable = self.out / "amps-pinned"
        executable.write_bytes(b"explicit linked executable fixture")
        fixture = self.write_json("host-fixture.json", {
            "schema": "srcsep-native-host-step-fixture-v1", "steps": 2, "source_macros": 4,
            "source_energy_mev": 1, "source_mu": 0.5, "constant_dmumu_per_s": 0.2})
        protocol = self.out / "fixture.args"
        protocol.write_text("SEP_NATIVE_HOST_STEP 1\nsteps 2\nsource_macros 4\nsource_energy_mev 1\nsource_mu 0.5\nconstant_dmumu_per_s 0.20000000000000001\n")
        guards = []
        for name, code, text in (("single-case", 1, "end-to-end test input requires exactly one explicit --test ID"),
                                 ("repeat-one", 2, "requires a fresh process")):
            log = self.out / (name + ".log");log.write_text(text)
            guards.append({"name": name, "status": "PASS", "exit_code": code, "log": str(log)})
        jobs = []
        for index, (identifier, mover) in enumerate((("HOST01", "parker"), ("HOST02", "fte-dmumu"), ("HOST03", "fte-mfp"))):
            child = self.write_json(identifier + "/srcsep-tests.json", {"results": [{"id": identifier, "status": "PASS", "metrics": [
                {"name": "assertion_failures", "value": 0}, {"name": "completed_host_steps", "value": 2},
                {"name": "global_first_birth_moves", "value": 2}]}]})
            for rank in range(2):
                prefix = identifier + "/artifacts/" + ("" if rank == 0 else "rank-1/")
                events = []
                for step in range(2):
                    event = {"step": step, "particle": 100 * (index + 1) + rank + 1000 * step,
                             "phase": "source", "clock_s": step, "parent_start_s": step, "parent_end_s": step + 1,
                             "first_duration_s": 1, "requested_s": 0, "source_radius_m": 1e10,
                             "pending_first": True, "cartesian_matches": True}
                    events.extend([event, dict(event, phase="mover", requested_s=1)])
                for sequence in range(2):
                    self.write_json(prefix + "host-step-trace-" + format(sequence, "06d") + ".json",
                        {"schema": "srcsep-native-host-step-trace-v2", "case_id": identifier,
                         "mover": mover, "rank": rank, "counter_at_end": True,
                         "snapshot_sequence": sequence, "complete": sequence == 1,
                         "requested_host_steps": 2, "completed_host_steps": 2 if sequence else 0,
                         "events": events if sequence else events[:1]})
            jobs.append({"id": identifier, "status": "PASS", "exit_code": 0, "report": str(child),
                         "log": str(child.parent / "native.log"), "postchecks": [{"status": "PASS"}]})
        payload = {"schema": "srcsep-native-host-step-campaign-v1", "executable": str(executable),
                   "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(), "fixture": str(fixture),
                   "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(), "mpi_ranks": 2,
                   "protocol_input": str(protocol), "protocol_input_sha256": hashlib.sha256(protocol.read_bytes()).hexdigest(),
                   "preflight_checks": guards, "results": jobs}
        return self.write_json("host-summary.json", payload), payload

    def test_host_reader_requires_real_guard_child_and_rank_trace_evidence(self):
        report, payload = self.make_host()
        rows = orch.reader_native_host_step(report)(self.ctx("host-step-native"))
        self.assertEqual([row.status for row in rows], ["PASS"] * 3)
        # Each defect starts from a newly materialized consistent campaign;
        # these are reader negative controls, never native physics evidence.
        for defect in ("binary", "fixture", "protocol", "protocol-mismatch", "guard-exit", "guard-diagnostic", "guard-mesh", "child-exit",
                       "child-skip", "extra-child", "missing-trace", "duplicate-first", "bad-x", "incomplete", "rewrite-prefix",
                       "failed-postcheck", "duplicate-job", "missing-job", "runner-exit"):
            with self.subTest(defect=defect):
                report, payload = self.make_host()
                job = payload["results"][0]
                if defect in ("binary", "fixture", "protocol", "protocol-mismatch"):
                    key = {"binary": "executable", "fixture": "fixture", "protocol": "protocol_input", "protocol-mismatch": "protocol_input"}[defect]
                    path = Path(payload[key]);path.write_bytes(b"changed")
                    if defect == "protocol-mismatch":payload[key + "_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
                elif defect == "guard-exit":payload["preflight_checks"][0]["exit_code"] = 2
                elif defect == "guard-diagnostic":Path(payload["preflight_checks"][0]["log"]).write_text("must be provided together")
                elif defect == "guard-mesh":
                    path = Path(payload["preflight_checks"][0]["log"]);path.write_text(path.read_text() + " The mesh is generated")
                elif defect == "child-exit":job["exit_code"] = 7
                elif defect in ("child-skip", "extra-child"):
                    child = Path(job["report"]);value = json.loads(child.read_text())
                    if defect == "child-skip":value["results"][0]["status"] = "SKIP"
                    else:value["results"].append(value["results"][0])
                    child.write_text(json.dumps(value))
                elif defect in ("missing-trace", "duplicate-first", "bad-x", "incomplete", "rewrite-prefix"):
                    trace = Path(job["report"]).parent / "artifacts/rank-1/host-step-trace-000001.json"
                    if defect == "missing-trace":trace.unlink()
                    else:
                        value = json.loads(trace.read_text())
                        if defect == "duplicate-first":value["events"].append(value["events"][-1])
                        elif defect == "bad-x":value["events"][0]["cartesian_matches"] = False
                        elif defect == "incomplete":value["complete"] = False
                        else:value["events"][0]["source_radius_m"] += 1e5
                        trace.write_text(json.dumps(value))
                elif defect == "failed-postcheck":job["postchecks"][0]["status"] = "FAIL"
                elif defect == "duplicate-job":payload["results"].append(job)
                elif defect == "missing-job":payload["results"].pop(0)
                report.write_text(json.dumps(payload))
                rows = orch.reader_native_host_step(report)(self.ctx("host-step-native", exit_code=1 if defect == "runner-exit" else 0))
                self.assertEqual(rows[0].status, "ERROR")
                self.assertTrue(rows[0].message)

    def make_radius(self):
        payload = {"schema": "srcsep-native-coefficient-radius-v1", "id": "RADIUSN01", "status": "PASS", "exit_code": 0}
        for key in ("executable", "template", "owned_input"):
            path = self.out / (key + ".input");path.write_text("retained " + key)
            payload[key] = str(path);payload[key + "_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        child = self.write_json("radius/srcsep-tests.json", {"results": [{"id": "RADIUSN01", "status": "PASS",
            "metrics": [{"name": "assertion_failures", "value": 0}]}]})
        payload["report"] = str(child)
        return self.write_json("radius-summary.json", payload), payload

    def test_radius_reader_rejects_prepared_changed_input_skip_and_hidden_child(self):
        report, _ = self.make_radius()
        self.assertEqual(orch.reader_native_coefficient_radius(report)(self.ctx())[0].status, "PASS")
        for defect in ("prepared", "binary", "template", "owned-input", "exit", "runner-exit", "skip", "extra-child", "assertion"):
            with self.subTest(defect=defect):
                report, payload = self.make_radius()
                if defect == "prepared":payload["status"] = "PREPARED"
                elif defect in ("binary", "template", "owned-input"):
                    key = {"binary": "executable", "template": "template", "owned-input": "owned_input"}[defect]
                    Path(payload[key]).write_text("changed")
                elif defect == "exit":payload["exit_code"] = 1
                elif defect in ("skip", "extra-child", "assertion"):
                    child = Path(payload["report"]);value = json.loads(child.read_text())
                    if defect == "skip":value["results"][0]["status"] = "SKIP"
                    elif defect == "extra-child":value["results"].append(value["results"][0])
                    else:value["results"][0]["metrics"][0]["value"] = 1
                    child.write_text(json.dumps(value))
                report.write_text(json.dumps(payload))
                row = orch.reader_native_coefficient_radius(report)(self.ctx(exit_code=1 if defect == "runner-exit" else 0))[0]
                self.assertEqual(row.status, "ERROR")

    def test_radius_reader_pins_selected_library_case_in_same_versioned_schema(self):
        report, payload = self.make_radius()
        payload["id"] = "RADIUSN02"
        child = Path(payload["report"]);value = json.loads(child.read_text())
        value["results"][0]["id"] = "RADIUSN02";child.write_text(json.dumps(value))
        report.write_text(json.dumps(payload))
        row = orch.reader_native_coefficient_radius(report, "RADIUSN02")(self.ctx())[0]
        self.assertEqual((row.test_id, row.status), ("RADIUSN02", "PASS"))
        self.assertEqual(orch.reader_native_coefficient_radius(report)(self.ctx())[0].status, "ERROR")


class Coupled(Fixture):
    def test_existing_diagnostic_used_else_output_copied(self):
        run_dir = self.out / "cr" / "runs" / "r1"
        diag = run_dir / "failures" / "native-amps" / "N01.log"
        diag.parent.mkdir(parents=True)
        diag.write_text("diagnostic")
        self.write_json("cr/summary.json", {"run_directory": str(run_dir), "results": [
            {"id": "N01", "scope": "native-amps", "status": "FAIL", "message": "m1"},
            {"id": "S01", "scope": "shared-model", "status": "ERROR", "message": "m2",
             "output": "captured traceback"},
            {"id": "S02", "scope": "shared-model", "status": "SKIP", "message": ""}]})
        records = {r.test_id: r for r in orch.reader_coupled(self.out / "cr")(self.ctx("coupled-runner"))}
        self.assertEqual(records["native-amps/N01"].log, str(diag))
        self.assertIn("captured traceback", Path(records["shared-model/S01"].log).read_text())
        self.assertEqual(records["shared-model/S02"].status, "SKIP")


class ReducedAndPositive(Fixture):
    def test_reduced_front_and_positive_shock_logs(self):
        self.write_json("rf/summary.json", {"results": [
            {"scope": "native", "test_id": "RSF07", "status": "FAIL", "message": "t",
             "log": "/rf/native.log", "report": "/rf/r.json"}]})
        self.write_json("ps/summary.json", {"checks": [
            {"run": "rank-4", "name": "arrival", "status": "FAIL", "detail": "late",
             "log": "/ps/execution.log"}]})
        rf = orch.reader_reduced_front(self.out / "rf")(self.ctx("reduced-front"))[0]
        ps = orch.reader_positive_shock(self.out / "ps")(self.ctx("positive-shock"))[0]
        self.assertEqual((rf.test_id, rf.log, rf.evidence), ("native/RSF07", "/rf/native.log", ["/rf/r.json"]))
        self.assertEqual((ps.test_id, ps.message, ps.log), ("rank-4/arrival", "late", "/ps/execution.log"))


class CoronalCme(Fixture):
    def test_fresh_report_copies_output_and_stale_report_is_rejected(self):
        report = self.write_json("cc/results.json", {"tests": [
            {"id": "SCCM01", "status": "FAIL", "message": "m", "output": "assert x"},
            {"id": "SCCM02", "status": "SKIP", "message": ""}]})
        records = orch.reader_coronal_cme(report)(self.ctx("coronal-cme"))
        self.assertIn("assert x", Path(records[0].log).read_text())
        old = time.time() - 3600
        os.utime(report, (old, old))
        with self.assertRaises(orch.ResultsUnavailable):
            orch.reader_coronal_cme(report)(self.ctx("coronal-cme"))


class UnittestText(Fixture):
    def test_fully_qualified_method_format_preserves_case_and_failure_log(self):
        # Python 3.14 prints the method in parentheses as well as before it.
        # This is actual stream syntax, not a change to the test's stable ID.
        text = ("test_a (test_mod.Viewer.test_a) ... ok\n"
                "test_b (test_mod.Viewer.test_b) ... FAIL\n"
                "FAIL: test_b (test_mod.Viewer.test_b)\n"
                "AssertionError: retained native failure\n")
        records = {r.test_id:r for r in orch.read_unittest(self.ctx("new-python", text))}
        self.assertEqual({k:r.status for k,r in records.items()},
                         {"Viewer.test_a":"PASS", "Viewer.test_b":"FAIL"})
        self.assertIn("retained native failure", Path(records['Viewer.test_b'].log).read_text())

    def test_launcher_step_requires_every_control_and_rejects_missing_or_skip(self):
        step = next(s for s in orch.srcsep_steps(self.out,argparse.Namespace())
                    if s.name=='standalone-application-runner')
        # Discover the real unittest names from syntax without importing its
        # fixtures, then require every name through the maintained step reader.
        source=ast.parse((_HERE.parent/'srcSEP/test/test_standalone_application_runner.py').read_text())
        methods=[node.name for owner in source.body if isinstance(owner,ast.ClassDef)
                 for node in owner.body if isinstance(node,ast.FunctionDef) and node.name.startswith('test_')]
        self.assertEqual(len(methods),11)
        def lines(selected,skip=False):
            return ''.join(f"{name} (__main__.ApplicationPlan.{name}) ... " +
                ("skipped 'missing prerequisite'" if skip and index==0 else 'ok')+'\n'
                for index,name in enumerate(selected))
        records=step.reader(self.ctx(log_text=lines(methods)))
        self.assertEqual({r.test_id for r in records},{'ApplicationPlan.'+m for m in methods})
        self.assertTrue(all(r.status=='PASS' for r in records))
        for text in (lines(methods[:-1]),lines(methods,skip=True)):
            self.assertTrue(any(r.status=='ERROR' for r in step.reader(self.ctx(log_text=text))))

    def test_verbose_lines_docstrings_and_failure_blocks(self):
        text = ("test_a (test_mod.Viewer) ... ok\n"
                "test_b (test_mod.Viewer)\nDocstring first line. ... FAIL\n"
                "test_c (test_mod.Viewer) ... skipped 'no display'\n"
                "\n" + "=" * 70 + "\n"
                "FAIL: test_b (test_mod.Viewer)\n" + "-" * 70 + "\n"
                "Traceback (most recent call last):\nAssertionError: 3 != 4\n\n"
                + "-" * 70 + "\nRan 3 tests in 0.1s\n\nFAILED (failures=1, skipped=1)\n")
        records = {r.test_id: r for r in orch.read_unittest(self.ctx("viewer-tests", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"Viewer.test_a": "PASS", "Viewer.test_b": "FAIL",
                          "Viewer.test_c": "SKIP"})
        block = Path(records["Viewer.test_b"].log).read_text()
        self.assertIn("FAIL: test_b", block)
        self.assertIn("# source:", block)

    def test_no_result_lines_is_unavailable(self):
        with self.assertRaises(orch.ResultsUnavailable):
            orch.read_unittest(self.ctx("viewer-tests", "ImportError: x\n"))


class BracketText(Fixture):
    def test_bracket_and_named_verdicts_worst_status_wins(self):
        text = ("[CMBGU01] PASS config\n[CMBGU05-CONTACT] flux_bound=1e-17\n"
                "BG3D-4 FAIL: leakage exceeds bound\n[ARCH01] PASS ok\n"
                "[CMBGU01] FAIL later\n")
        records = {r.test_id: r for r in orch.read_bracket_text(self.ctx("corona-swcme", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"CMBGU01": "FAIL", "BG3D-4": "FAIL", "ARCH01": "PASS"})
        self.assertEqual(records["BG3D-4"].line, 3)
        self.assertEqual(records["CMBGU01"].line, 5)


class MakeRecipe(Fixture):
    def test_every_binary_counted_and_hidden_failures_recovered(self):
        text = ("[ARCH01] PASS headers\n"
                "build/test_a\n[A1] PASS one\nBG3D-4 FAIL: leakage\n"
                "make: [makefile:122: test] Error 1 (ignored)\n"
                "build/test_b\nmetric=1\n"
                "build/test_c\nsegfault\n"
                "make: [makefile:122: test] Error 139 (ignored)\n"
                "build/test_d\n[D1] PASS ok\n"
                "make: [makefile:122: test] Error 2 (ignored)\n")
        reader = orch.reader_make_recipe(r"^build/(test_\w+)$")
        records = {r.test_id: r for r in reader(self.ctx("corona-swcme", text))}
        self.assertEqual({k: r.status for k, r in records.items()},
                         {"ARCH01": "PASS", "A1": "PASS", "BG3D-4": "FAIL",
                          "test_b": "PASS", "test_c": "ERROR", "D1": "PASS",
                          "test_d": "ERROR"})
        self.assertEqual(records["test_c"].message, "exited with status 139")
        self.assertEqual(records["test_c"].line, 8)


class RngLifecycleRegistration(Fixture):
    def test_finite_amplitude_step_requires_all_sixteen_actual_unittest_names(self):
        steps = {step.name: step for step in orch.srcsep_steps(self.out, argparse.Namespace())}
        step = steps["finite-amplitude-study"]
        self.assertEqual(step.commands, [["make", "-C", "srcSEP", "test-finite-amplitude-study-unit"]])
        # Inspect definitions only; never execute the mathematical study while
        # checking its discovery/report wiring. This catches misspelled class
        # or method selectors in the independently declared required-ID list.
        source = orch.AMPS_ROOT / "srcSEP/test/test_finite_amplitude_response.py"
        tree = ast.parse(source.read_text(encoding="utf-8"))
        names = [node.name + "." + method.name for node in tree.body if isinstance(node, ast.ClassDef)
                 for method in node.body if isinstance(method, ast.FunctionDef) and method.name.startswith("test_")]
        self.assertEqual(len(names), 16)
        def line(identifier, verdict):
            owner, method = identifier.split(".", 1)
            return method + " (__main__." + owner + ") ... " + verdict + "\n"
        records = step.reader(self.ctx("finite-amplitude-study", "".join(line(name, "ok") for name in names)))
        self.assertEqual({record.test_id for record in records}, set(names))
        self.assertTrue(all(record.status == "PASS" for record in records))
        broken = "".join(line(name, "skipped 'unavailable prerequisite'" if index == 0 else "ok")
                         for index, name in enumerate(names[:-1]))
        records = step.reader(self.ctx("finite-amplitude-study", broken))
        by_id = {record.test_id: record.status for record in records}
        self.assertEqual(by_id[names[0]], "ERROR")
        self.assertEqual(by_id[names[-1]], "ERROR")
        self.assertEqual(len(records), 16)

    def test_new_source_families_and_fresh_checkpoint_have_complete_readers(self):
        steps = {step.name: step for step in orch.srcsep_steps(self.out, argparse.Namespace())}
        for name, target, prefix, number in (
                ("analytic-radial", "test-analytic-radial-background-unit", "ARB", 6),
                ("analytic-shock-epoch", "test-analytic-shock-epoch-unit", "ASE", 5),
                ("material-tube-kinematics", "test-material-tube-kinematics-unit", "MTG", 6),
                ("standard-batsrus-waves", "test-standard-batsrus-wave-unit", "SBW", 5),
                ("flux-tube-state", "test-flux-tube-state-unit", "FLCP", 4),
                ("polarized-wkb", "test-polarized-wkb-unit", "WKB", 8),
                ("maxwell-juttner-flux", "test-maxwell-juttner-flux-unit", "MJ", 6)):
            with self.subTest(family=name):
                step = steps[name]
                self.assertEqual(step.commands, [["make", "-C", "srcSEP", target]])
                identifiers = ["%s%02d" % (prefix, index) for index in range(1, number + 1)]
                log = "".join("[%s] %s literal callback\n" % (identifier, "SKIP" if index == 0 else "PASS")
                              for index, identifier in enumerate(identifiers[:-1]))
                records = step.reader(self.ctx(name, log))
                by_id = {record.test_id: record.status for record in records}
                self.assertEqual(by_id[identifiers[0]], "ERROR")
                self.assertEqual(by_id[identifiers[-1]], "ERROR")
                self.assertEqual(len(records), number)
        fresh = steps["fresh-process-checkpoint"]
        self.assertEqual(fresh.commands[0][:3], ["python3", "-B", "srcSEP/test/run_fresh_process_checkpoint_tests.py"])
        command = fresh.commands[0]
        self.assertEqual(command[command.index("--mpi-ranks") + 1], "2")
        self.assertEqual(command[command.index("--output") + 1], str(self.out / "fresh-process-checkpoint"))
        radial = steps["analytic-radial-native"]
        command = radial.commands[0]
        self.assertEqual(command[:3], ["python3", "-B", "srcSEP/test/run_tests.py"])
        self.assertEqual(command[command.index("--group")+1], "analytic-radial-native")
        self.assertEqual(command[command.index("--mpi-np")+1], "2")
        report = self.write_json("analytic-radial-native/srcsep-tests.json", {"results": [
            {"id": "ARBN01", "status": "PASS"}, {"id": "ARBN02", "status": "SKIP"}]})
        records = radial.reader(self.ctx("analytic-radial-native"))
        self.assertEqual([(record.test_id, record.status) for record in records],
                         [("ARBN01", "PASS"), ("ARBN02", "ERROR"), ("ARBN03", "ERROR"), ("ARBN04", "ERROR")])
        for name, script in (("host-step-native", "run_native_host_step_tests.py"),
                             ("coefficient-radius-native", "run_native_coefficient_radius_tests.py"),
                             ("coefficient-radius-library-native", "run_native_coefficient_radius_tests.py")):
            command = steps[name].commands[0]
            self.assertEqual(command[:3], ["python3", "-B", "srcSEP/test/" + script])
            self.assertEqual(command[command.index("--mpi-ranks") + 1], "2")
            self.assertEqual(command[command.index("--output") + 1], str(self.out / name))
        host_command = steps["host-step-native"].commands[0]
        self.assertEqual(host_command[host_command.index("--fixture") + 1], "srcSEP/test/fixtures/host_step_smoke.json")
        library_command = steps["coefficient-radius-library-native"].commands[0]
        self.assertEqual(library_command[library_command.index("--case") + 1], "RADIUSN02")
        for name, group, required in (
                ("material-tube-native", "material-tube", ["MTG01", "MTG02", "MTG03", "MTG04", "MTG05", "MTG06"]),
                ("standard-batsrus-waves-native", "standard-batsrus-waves-native", ["SBWN01", "SBWN02"])):
            with self.subTest(native=name):
                step = steps[name]
                command = step.commands[0]
                self.assertEqual(command[:3], ["python3", "-B", "srcSEP/test/run_tests.py"])
                self.assertEqual(command[command.index("--group")+1], group)
                self.assertEqual(command[command.index("--mpi-np")+1], "2")
                self.write_json(name+"/srcsep-tests.json", {"results": [{"id": required[0], "status": "SKIP"}]})
                records = step.reader(self.ctx(name))
                self.assertEqual({record.test_id: record.status for record in records},
                                 {identifier: "ERROR" for identifier in required})

    def test_birth_warm_and_native_polarized_steps_retain_individual_gates(self):
        steps = {step.name: step for step in orch.srcsep_steps(self.out, argparse.Namespace())}
        birth = steps["birth-window"]
        self.assertEqual(birth.commands, [["make", "-C", "srcSEP", "test-birth-window-unit"]])
        records = birth.reader(self.ctx("birth-window", "[BIRTH01] PASS interval\n[BIRTH02] FAIL rounded-away positive remainder\n"))
        self.assertEqual([(r.test_id, r.status) for r in records], [("BIRTH01", "PASS"), ("BIRTH02", "FAIL")])
        warm = steps["warm-parallel-study"]
        self.assertEqual(warm.commands, [["make", "-C", "srcSEP", "test-warm-parallel-study-unit"]])
        warm_log = ("test_default_precision_serialized_frequency_roundtrip (__main__.WarmParallelTests) ... ok\n"
                    "test_winding_alias_negative_control (__main__.WarmParallelTests) ... FAIL\n")
        records = warm.reader(self.ctx("warm-parallel-study", warm_log, exit_code=1))
        self.assertEqual([(r.test_id, r.status) for r in records],
                         [("WarmParallelTests.test_default_precision_serialized_frequency_roundtrip", "PASS"),
                          ("WarmParallelTests.test_winding_alias_negative_control", "FAIL")])
        native = steps["polarized-wave-native"]
        command = native.commands[0]
        self.assertEqual(command[:4], ["python3", "-B", "srcSEP/test/run_tests.py", "--amps"])
        self.assertEqual(command[command.index("--mpi-np")+1], "2")
        self.assertEqual(command[command.index("--group")+1], "polarized-wave-native")
        self.assertEqual(command[command.index("--shock-model")+1], "swcme1d")
        self.assertEqual(command[command.index("--turbulence-model")+1], "wave-number-resolved")
        self.assertEqual(command[command.index("--turbulence-source")+1], "self-consistent-spectral")
        self.assertEqual(command[command.index("--coupling")+1], "off")
        birth_native=steps["birth-window-native"].commands[0]
        self.assertEqual(birth_native[:4],["python3","-B","srcSEP/test/run_tests.py","--amps"])
        self.assertEqual(birth_native[birth_native.index("--group")+1],"birth-window-native")
        self.assertEqual(birth_native[birth_native.index("--mpi-np")+1],"2")
        frozen=steps["frozen-collision-transaction"]
        self.assertEqual(frozen.commands,[["make","-C","srcSEP","test-frozen-collision-transaction-unit"]])
        records=frozen.reader(self.ctx("frozen-collision-transaction","[FCT01] PASS ledger\n[FCT03] FAIL budget\n",exit_code=1))
        self.assertEqual([(record.test_id,record.status) for record in records],[("FCT01","PASS"),("FCT03","FAIL")])

    def test_background_gradient_gate_retains_its_failed_physics_case(self):
        steps = {step.name: step for step in
                 orch.srcsep_steps(self.out, argparse.Namespace())}
        step = steps["background-derivatives"]
        self.assertEqual(step.commands,
                         [["make", "-C", "srcSEP", "test-background-derivatives-unit"]])
        log = ("[CP03G01] FAIL fixed-vertex expansion\n"
               "[CP03G02] PASS curved-field strain\n"
               "[CP03G03] PASS tensor validation\n")
        records = step.reader(self.ctx("background-derivatives", log, exit_code=1,
                                       namespace=step.namespace))
        self.assertEqual([(r.test_id, r.status) for r in records],
                         [("CP03G01", "FAIL"), ("CP03G02", "PASS"),
                          ("CP03G03", "PASS")])

    def test_source_lifecycle_gate_is_selectable_and_keeps_individual_failures(self):
        # The development-plan contract requires this source gate to remain
        # reachable from run_all_test.sh even before a new AMPS binary exists.
        # A literal child log verifies that one failed stochastic assertion is
        # retained individually instead of being hidden by five passing cases.
        steps = {step.name: step for step in
                 orch.srcsep_steps(self.out, argparse.Namespace())}
        step = steps["rng-lifecycle"]
        self.assertEqual(step.commands,
                         [["make", "-C", "srcSEP", "test-rng-lifecycle-unit"]])
        self.assertEqual(step.namespace, "srcsep-registry")
        log = ("[RNG01] PASS snapshot\n[RNG02] PASS cache\n"
               "[RNG03] FAIL repeated-call variance\n"
               "[RNG04] PASS event partitions\n[RNG05] PASS restored identity\n"
               "[RNG06] PASS invalid state rejected\n")
        records = step.reader(self.ctx("rng-lifecycle", log, exit_code=1,
                                       namespace=step.namespace))
        self.assertEqual(len(records), 6)
        self.assertEqual([(r.test_id, r.status) for r in records
                          if r.status != "PASS"], [("RNG03", "FAIL")])
        self.assertEqual(records[2].message, "repeated-call variance")


class StepGuards(Fixture):
    def run_step(self, command, reader):
        step = orch.Step("guard", "fixture", [command], reader=reader, runner="r")
        args = argparse.Namespace(stream=False)
        with contextlib.redirect_stdout(io.StringIO()):
            return orch.execute(step, 1, self.out, orch.Console(None), args)

    def test_nonzero_exit_without_failing_record_is_an_error(self):
        run = self.run_step(["sh", "-c", "exit 3"],
                            lambda ctx: [ctx.record("T1", "PASS")])
        self.assertEqual([r.test_id for r in run.records], ["T1", "guard:exit-status"])
        self.assertEqual(run.status, "FAIL")

    def test_missing_report_becomes_results_unavailable(self):
        def reader(ctx):
            raise orch.ResultsUnavailable("report not found: x")
        run = self.run_step(["sh", "-c", "echo hi"], reader)
        self.assertEqual([(r.test_id, r.status) for r in run.records],
                         [("guard:results-unavailable", "ERROR")])
        self.assertIn("hi", Path(run.log).read_text())


class Summary(Fixture):
    def test_counts_unique_dedupe_files_and_exit_status(self):
        def record(test_id, status, step, namespace, kind="test"):
            return orch.TestRecord(test_id, status, "msg", step, "r", log="/log",
                                   namespace=namespace, kind=kind)
        runs = [
            orch.StepRun(orch.Step("build", "", [], kind="check"), "PASS", 0, 1.0, "",
                         [record("build", "PASS", "build", "build", "check")]),
            orch.StepRun(orch.Step("a", "", []), "FAIL", 1, 1.0, "",
                         [record("X1", "PASS", "a", "shared"), record("X2", "FAIL", "a", "shared")]),
            orch.StepRun(orch.Step("b", "", []), "PASS", 0, 1.0, "",
                         [record("X1", "PASS", "b", "shared"), record("Y1", "SKIP", "b", "b")]),
        ]
        args = argparse.Namespace(dry_run=False, junit=self.out / "j.xml")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = orch.summarize("App", runs, orch.Console(None), self.out, args, "t0")
        text = buffer.getvalue()
        self.assertEqual(code, 1)
        self.assertIn("test executions: total=4 pass=2 fail=1 skip=1 error=0", text)
        self.assertIn("unique tests:    total=3 pass=1 fail=1 skip=1 error=0", text)
        self.assertIn("checks:          total=1 pass=1 fail=0", text)
        self.assertIn("FAIL   X2", text)
        report = json.loads((self.out / "run_all_test.json").read_text())
        self.assertEqual(report["counts"]["unique_tests"]["TOTAL"], 3)
        self.assertIn("X2", (self.out / "failed_tests.txt").read_text())
        self.assertIn('failures="1"', (self.out / "j.xml").read_text())

    def test_all_passing_returns_zero(self):
        runs = [orch.StepRun(orch.Step("a", "", []), "PASS", 0, 1.0, "",
                             [orch.TestRecord("X1", "PASS", "", "a", "r", namespace="a")])]
        args = argparse.Namespace(dry_run=False, junit=None)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(orch.summarize("App", runs, orch.Console(None), self.out, args, "t0"), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
