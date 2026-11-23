#!/usr/bin/env python3
"""Software boundary tests for maintained Routine/Birth entry points.

All fabricated reports are explicitly labelled software-only. They exercise
selection, export and process-status wiring; the real native/physical gates
remain separate and cannot be satisfied by this fixture.
"""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

SOURCE=Path(__file__).resolve().parents[1]
# The application is copied from AMPS/srcSEP to AMPS/build/main. Tests use the
# same two documented enclosing roots as the application makefile, so neither
# source nor generated-layout checks depend on the caller's working directory.
AMPS_ROOT=next(path for path in (SOURCE.parent,SOURCE.parent.parent)
               if (path/'ampsConfig.pl').is_file())


def load(path):
    spec=importlib.util.spec_from_file_location('srcsep_wiring_software_boundary',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


class WiringTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='sep-routine-wiring-software-')
        self.addCleanup(self.temporary.cleanup);self.directory=Path(self.temporary.name)
        self.runner=load(SOURCE/'test/run_tests.py');self.calls=[]
        self.ids=json.loads((SOURCE/'test/routine_fixture_manifest.json').read_text())['expected_ids']

    def report(self, ids, status='PASS'):
        rows=[{'id':i,'name':i,'group':'software-only','status':status,
               'message':'SOFTWARE BOUNDARY ONLY; NOT A MODEL CALCULATION',
               'configuration':['software_only=true'],'metrics':[],'artifacts':[],'elapsed_seconds':0,'seed':None} for i in ids]
        return {'schema':'srcsep-component-tests-v1','exit_code':2 if status=='ERROR' else 0,
                'results':rows,'totals':self.runner._report_totals(rows)}

    def dispatcher(self,status='PASS'):
        def run(args):
            self.calls.append(args)
            self.assertFalse(args.output.exists())
            args.output.mkdir(parents=True)
            report=self.report(getattr(args,"requested_ids",None) or self.ids,status)
            self.runner._write_json(args.output/'srcsep-tests.json',report)
            self.runner._write_json(args.output/'run_manifest.json',{'commands':[['SOFTWARE-ONLY',i] for i in self.ids]})
            return report['exit_code']
        return SimpleNamespace(run=run)

    def test_native_routine_command_cannot_collapse_to_all_tests(self):
        exe=self.directory/'labelled-executable';exe.write_text('software-only');exe.chmod(0o755)
        args=SimpleNamespace(amps=str(exe),all=False,routine=True)
        with self.assertRaisesRegex(self.runner.RunnerError,'descriptor-aware'):
            self.runner._native_command(args,self.directory/'x.json',self.directory/'x.xml')

    def test_routine_entry_preserves_full_exports_and_default_seed(self):
        out=self.directory/'out';export=self.directory/'explicit/report.json';xml=self.directory/'explicit/report.xml'
        with patch.object(self.runner,'_routine_protocol_module',return_value=self.dispatcher()),redirect_stdout(io.StringIO()):
            code=self.runner.main(['--routine','--amps','software-only-executable','--plot','none','--output-dir',str(out),
                                   '--report-json',str(export),'--report-junit',str(xml)])
        self.assertEqual(code,0);report=json.loads(export.read_text())
        self.assertEqual([r['id'] for r in report['results']],self.ids)
        self.assertIn('S13AD01', [r['id'] for r in report['results']])
        self.assertEqual(self.calls[0].native_args,[])
        self.assertIsNone(self.calls[0].mpi_np);self.assertTrue(xml.is_file())

    def test_repeat_seed_mpi_and_failure_status_are_not_changed(self):
        forwarded=['--test-seed=42','--test-repeat','2'];out=self.directory/'out'
        with patch.object(self.runner,'_routine_protocol_module',return_value=self.dispatcher('ERROR')),redirect_stdout(io.StringIO()):
            code=self.runner.main(['--routine','--amps','software-only','--plot','none','--mpi-np','2','--mpiexec','labelled-mpiexec',
                                   '--output-dir',str(out),'--',*forwarded])
        self.assertEqual(code,2);self.assertEqual(self.calls[0].native_args,forwarded)
        self.assertEqual((self.calls[0].mpi_np,self.calls[0].mpiexec),(2,'labelled-mpiexec'))
        self.assertEqual(json.loads((out/'srcsep-tests.json').read_text())['totals']['errors'],len(self.ids))

    def test_repeated_entry_retains_separate_raw_case_directories(self):
        out=self.directory/'out';out.mkdir();args=SimpleNamespace(amps='software-only',timeout=None,model_args=[],mpi_np=None,mpiexec='mpiexec')
        with patch.object(self.runner,'_routine_protocol_module',return_value=self.dispatcher()):
            self.runner._run_routine_tests(args,out);self.runner._run_routine_tests(args,out)
        self.assertNotEqual(self.calls[0].output,self.calls[1].output)
        self.assertTrue(all(c.output.is_dir() for c in self.calls))

    def test_birth_group_uses_all_six_complete_role_workflows(self):
        out=self.directory/'out';out.mkdir();ids=[f'BIRTHN0{i}' for i in range(1,7)]
        dispatcher=self.runner._routine_protocol_module()
        args=SimpleNamespace(amps='software-only',model_args=['--test-seed','123','--test-repeat=2'],mpi_np=2,mpiexec='mpiexec')
        commands=[]
        def execute(command,cwd,log,environment,timeout):
            commands.append(command);roleout=Path(command[command.index('--output')+1]);roleout.mkdir()
            self.runner._write_json(roleout/'srcsep-tests.json',self.report(ids));return 0
        with patch.object(self.runner,'_run_streaming',side_effect=execute):
            code,path,_=self.runner._run_birth_group(args,out,out/'log')
        self.assertEqual(code,0);self.assertEqual([r['id'] for r in json.loads(path.read_text())['results']],ids)
        cmd=commands[0];self.assertEqual(cmd[cmd.index('--mpi-ranks')+1],'2')
        self.assertEqual(cmd[cmd.index('--seed')+1],'123');self.assertEqual(cmd[cmd.index('--repeat-count')+1],'2')
        with self.assertRaisesRegex(ValueError,'does not support'):
            dispatcher.birth_campaign_arguments(['--particle-mover','parker'])

    def test_single_birth_id_uses_full_roles_and_default_two_rank_migration_contract(self):
        out=self.directory/'single';out.mkdir();commands=[]
        args=SimpleNamespace(amps='software-only',model_args=[],mpi_np=None,mpiexec='mpiexec',tests=['birthn05'])
        def execute(command,cwd,log,environment,timeout):
            commands.append(command);target=Path(command[command.index('--output')+1]);target.mkdir()
            self.runner._write_json(target/'srcsep-tests.json',self.report(['BIRTHN05']));return 0
        with patch.object(self.runner,'_run_streaming',side_effect=execute):
            code,path,_=self.runner._run_birth_group(args,out,out/'log')
        self.assertEqual(code,0);cmd=commands[0]
        self.assertEqual(cmd[cmd.index('--test-id')+1],'BIRTHN05')
        self.assertEqual(cmd[cmd.index('--mpi-ranks')+1],'2')
        self.assertEqual(len(json.loads(path.read_text())['results']),1)

    def test_selected_checkpoint_id_uses_same_reviewed_routine_protocol(self):
        out=self.directory/'cp'
        with patch.object(self.runner,'_routine_protocol_module',return_value=self.dispatcher()),redirect_stdout(io.StringIO()):
            code=self.runner.main(['--test','CP05C01','--amps','software-only','--plot','none','--output-dir',str(out)])
        self.assertEqual(code,0);self.assertEqual(self.calls[0].requested_ids,['CP05C01'])
        self.assertEqual([r['id'] for r in json.loads((out/'srcsep-tests.json').read_text())['results']],['CP05C01'])

    def test_modern_all_reuses_routine_fixture_and_keeps_extended_failure(self):
        out=self.directory/'all';out.mkdir();exe=self.directory/'software-only';exe.write_text('software-only');exe.chmod(0o755)
        args=SimpleNamespace(amps=str(exe),model_args=[],mpi_np=None,mpiexec='mpiexec',timeout=None)
        def discovery(executable,model_args,cwd):
            (cwd/'discovery.log').write_text('ID | group | class | initialization | build modes | description\n'
                'CP05C01 | checkpoint | routine | field-line-model | software-only | labelled fixture\n'
                'EXT99 | stress | extended | none | software-only | labelled extended fixture\n')
            return ['CP05C01','EXT99']
        routine_state={'status':'PASS','exit':0};extended_status={'value':'FAIL'}
        def routine(arguments,target,requested):
            self.assertEqual(requested,['CP05C01']);report=target/'srcsep-tests.json'
            self.runner._write_json(report,self.report(requested,routine_state['status']))
            return routine_state['exit'],report,[['SOFTWARE-ROUTINE','CP05C01']]
        calls=[]
        def extended(command,cwd,log,environment,timeout):
            calls.append(command);self.assertEqual(command[command.index('--test')+1],'EXT99')
            report=self.report(['EXT99']);report['results'][0]['status']=extended_status['value'];report['exit_code']=1 if extended_status['value']=='FAIL' else 0
            self.runner._write_json(Path(command[command.index('--test-json')+1]),report);return report['exit_code']
        with patch.object(self.runner,'_discover_all_ids',side_effect=discovery),patch.object(self.runner,'_run_routine_tests',side_effect=routine), \
             patch.object(self.runner,'_validation_case_ids',return_value=set()),patch.object(self.runner,'_run_streaming',side_effect=extended),redirect_stdout(io.StringIO()):
            code,report,_=self.runner._run_all_tests(args,out,out/'log')
            routine_state.update(status='SKIP',exit=2);extended_status['value']='PASS'
            second=self.directory/'required-skip';second.mkdir()
            skipped_code,skipped_report,_=self.runner._run_all_tests(args,second,second/'log')
        self.assertEqual(code,1);self.assertEqual(len(calls),2)
        self.assertEqual(skipped_code,2)
        self.assertEqual([(r['id'],r['status']) for r in json.loads(skipped_report.read_text())['results']],
                         [('CP05C01','SKIP'),('EXT99','PASS')])
        results=json.loads(report.read_text())['results']
        self.assertEqual([(r['id'],r['status']) for r in results],[('CP05C01','PASS'),('EXT99','FAIL')])

    def test_generated_layout_copies_required_module_manifest_and_resolves_root(self):
        # Match ampsConfig's unfiltered project-directory copy in an isolated
        # tree. This checks Python/JSON payloads, not merely C++ source files.
        project=self.directory/'srcSEP';target=self.directory/'build/main'
        files=['test/run_tests.py','test/routine_fixtures.py','test/routine_fixture_manifest.json',
               'test/run_birth_operator_stage_tests.py','validation/sep_native_execution.py']
        for rel in files:
            dest=project/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(SOURCE/rel,dest)
        shutil.copytree(project,target)
        relocated=load(target/'test/run_tests.py');dispatcher=relocated._routine_protocol_module()
        self.assertEqual(dispatcher.ROOT,target)
        self.assertEqual(json.loads((dispatcher.ROOT/'test/routine_fixture_manifest.json').read_text())['expected_ids'],self.ids)
        self.assertTrue(all((target/f).read_bytes()==(project/f).read_bytes() for f in files))
        config=(AMPS_ROOT/'ampsConfig.pl').read_text()
        self.assertIn('cp -r $ProjectSpecificSourceDirectory $ampsConfigLib::WorkingSourceDirectory/main',config)

    def test_make_and_orchestrator_preserve_required_gates(self):
        text=(SOURCE/'makefile').read_text()
        self.assertIn('@$(SEP_ROUTINE_COMMAND) --mpi-np "$(MPI_NP)" --mpiexec "$(MPIEXEC)"',text)
        self.assertEqual(text.count('@$(SEP_ROUTINE_COMMAND)'),3)
        self.assertIn('"$(SWCME_DIR)/test/run_tests.py"',text)
        self.assertIn('--routine --output-dir "$(abspath test_output/swcme-r1)"',text)
        self.assertIn('test-controlled-analytical:',text)
        orchestrator=(AMPS_ROOT/'tools/sep_test_orchestrator.py').read_text()
        self.assertIn('complete fixed-step role groups; all six native IDs required',orchestrator)
        for identity in ['BIRTHN01','BIRTHN02','BIRTHN03','BIRTHN04','BIRTHN05','BIRTHN06']:
            self.assertIn('"'+identity+'"',orchestrator)


if __name__=='__main__':unittest.main(verbosity=2)
