#!/usr/bin/env python3
"""Bounded real-process/status controls; no AMPS builds, MPI or observations.

Unittest discovery is supported for U16/global orchestration registration.
Direct --output-dir runs retain every controlled child/resource log. Tiny
temporary scripts are protocol fixtures only, never scientific evidence.
"""
from collections import Counter
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import m3_supporting_services as services
import run_all_test as global_runner

STATES=('PASS','FAIL','ERROR','SKIPPED')
EXITS={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}
EVIDENCE=None


class SupportingProtocolTests(unittest.TestCase):
    def setUp(self):
        if EVIDENCE is None:
            temporary=tempfile.TemporaryDirectory(prefix='amps-m3-support-protocol-')
            self.addCleanup(temporary.cleanup);self.root=Path(temporary.name)
        else:
            self.root=EVIDENCE/self._testMethodName;self.root.mkdir()

    def runner(self,name,body):
        path=self.root/(name+'.py')
        path.write_text('from pathlib import Path\nimport json,os,subprocess,sys,time\n'
            "out=Path(sys.argv[sys.argv.index('--output-dir')+1]);out.mkdir(parents=True)\n"+body+'\n')
        return path

    def receipt_runner(self,name,status,code=None,pause=0):
        return self.runner(name,"time.sleep("+repr(pause)+")\n"
            "(out/'result.json').write_text(json.dumps({'status':"+repr(status)+",'scope':'CONTROL ONLY'}))\n"
            'raise SystemExit('+repr(EXITS[status] if code is None else code)+')')

    def load(self,status,children,counts=None):
        case=self.root/('receipt-'+str(len(list(self.root.glob('receipt-*')))));case.mkdir()
        rows=[]
        for index,item in enumerate(children):
            if isinstance(item,str):
                log=case/('child-'+str(index)+'.log');log.write_text('CONTROL ONLY\n')
                item={'check':'CONTROL child '+str(index),'status':item,
                    'runner':str(self.root/'CONTROL-runner.py'),'log':str(log)}
            rows.append(item)
        if counts is None:counts={s:sum(isinstance(r,dict) and r.get('status')==s for r in rows) for s in STATES}
        payload={'id':'U15','name':'CONTROL ONLY','status':status,
            'supporting_checks':rows,'supporting_counts':counts}
        result=case/'result.json';result.write_text(json.dumps(payload));log=case/'runner.log';log.write_text('CONTROL ONLY\n')
        return global_runner.load_runner_result('U15',self.root/'CONTROL-runner.py',result,log,EXITS[status]),payload

    def test_real_child_status_and_exit_mappings(self):
        for status in STATES:
            with self.subTest(status=status):
                runner=self.receipt_runner('child-'+status,status)
                with patch.dict(services.RUNNERS,{'reaction-primitives':runner}):
                    result=services.run_services(self.root/status,('reaction-primitives',),10)
                self.assertEqual(result['status'],status)
                self.assertEqual(result['counts'],{s:int(s==status) for s in STATES})
                child=result['supporting_checks'][0]
                self.assertTrue(Path(child['log']).is_file())
                self.assertTrue(Path(child['resource_result']).is_file())
                self.assertTrue(Path(child['resource_history']).is_file())

    def test_child_status_exit_contradiction_is_error(self):
        runner=self.receipt_runner('bad-exit','PASS',1)
        with patch.dict(services.RUNNERS,{'reaction-primitives':runner}):
            result=services.run_services(self.root/'bad',('reaction-primitives',),10)
        self.assertEqual(result['status'],'ERROR')
        self.assertIn('contradiction',result['supporting_checks'][0]['reason'])

    def test_missing_or_nonobject_child_receipt_is_error(self):
        for name,body in [('missing','pass'),('nonobject',"(out/'result.json').write_text('[]')")]:
            with self.subTest(name=name):
                runner=self.runner(name,body)
                with patch.dict(services.RUNNERS,{'reaction-primitives':runner}):
                    result=services.run_services(self.root/(name+'-output'),('reaction-primitives',),10)
                self.assertEqual(result['status'],'ERROR')
                self.assertTrue(Path(result['supporting_checks'][0]['log']).is_file())

    def test_groups_share_one_clock_and_retain_each_resource_record(self):
        first=self.receipt_runner('first','PASS',pause=.15)
        second=self.receipt_runner('second','FAIL')
        with patch.dict(services.RUNNERS,{'reaction-primitives':first,'runtime-input':second}):
            value=services.run_services(self.root/'both',('reaction-primitives','runtime-input'),10)
        self.assertEqual(value['status'],'FAIL')
        a,b=value['supporting_checks']
        first_budget=json.loads(Path(a['resource_result']).read_text())
        second_history=Path(b['resource_history']).read_text().splitlines()
        second_start=float(second_history[1].split(',')[0])
        self.assertGreaterEqual(second_start,first_budget['wall_seconds'])
        self.assertNotEqual(a['resource_result'],b['resource_result'])
        self.assertTrue(Path(a['resource_history']).is_file())

    def test_timeout_reaps_grandchild_and_does_not_launch_remaining_group(self):
        if not Path('/proc/self/stat').is_file():self.skipTest('Linux owned-tree observer unavailable')
        first=self.runner('timeout',
            "child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'])\n"
            "fields=(Path('/proc')/str(child.pid)/'stat').read_text().rsplit(')',1)[1].split()\n"
            "(out/'owned-child.json').write_text(json.dumps({'pid':child.pid,'birth':int(fields[19])}))\n"
            'time.sleep(30)')
        second=self.receipt_runner('not-launched','PASS')
        child=None;begin=time.monotonic()
        try:
            with patch.dict(services.RUNNERS,{'reaction-primitives':first,'runtime-input':second}):
                value=services.run_services(self.root/'timeout-output',('reaction-primitives','runtime-input'),3)
            child=json.loads((self.root/'timeout-output/reaction-primitives/owned-child.json').read_text())
            record=Path('/proc')/str(child['pid'])/'stat'
            if record.is_file():
                fields=record.read_text().rsplit(')',1)[1].split()
                self.assertNotEqual(int(fields[19]),child['birth'],'owned grandchild was not reaped')
            self.assertLess(time.monotonic()-begin,5)
            self.assertEqual(value['status'],'ERROR')
            self.assertEqual(value['counts'],{'PASS':0,'FAIL':0,'ERROR':2,'SKIPPED':0})
            pending=value['supporting_checks'][1]
            self.assertEqual(pending['execution'],'not-run')
            self.assertTrue(Path(pending['log']).is_file())
            self.assertFalse((self.root/'timeout-output/runtime-input').exists())
        finally:
            # Even a regressed supervisor must not leave this protocol fixture
            # alive. Signal only the recorded PID with its exact birth identity.
            marker=self.root/'timeout-output/reaction-primitives/owned-child.json'
            if child is None and marker.is_file():child=json.loads(marker.read_text())
            if child:
                try:
                    fields=(Path('/proc')/str(child['pid'])/'stat').read_text().rsplit(')',1)[1].split()
                    if int(fields[19])==child['birth']:
                        os.kill(child['pid'],signal.SIGKILL)
                        try:os.waitpid(child['pid'],0)
                        except ChildProcessError:pass
                except FileNotFoundError:pass

    def test_formal_skip_does_not_promote_successful_or_unavailable_support(self):
        for children in (['PASS'],['PASS','SKIPPED'],[]):
            with self.subTest(children=children):
                value,payload=self.load('SKIPPED',children)
                self.assertEqual(value['status'],'SKIPPED')
                self.assertEqual(value['supporting_checks'],payload['supporting_checks'])

    def test_valid_nested_hierarchy_is_preserved(self):
        for status,children in [('PASS',['PASS']),('FAIL',['FAIL']),('FAIL',['PASS']),('ERROR',['ERROR']),('ERROR',['PASS'])]:
            with self.subTest(status=status,children=children):
                value,payload=self.load(status,children)
                self.assertEqual(value['status'],status)
                self.assertEqual(value['supporting_counts'],payload['supporting_counts'])

    def test_nested_failure_cannot_hide_and_prints_actual_runner_and_log(self):
        for status,child in [('SKIPPED','FAIL'),('SKIPPED','ERROR'),('PASS','FAIL'),('PASS','ERROR'),('PASS','SKIPPED'),('FAIL','ERROR')]:
            with self.subTest(status=status,child=child):
                value,payload=self.load(status,[child])
                self.assertEqual(value['status'],'ERROR')
                self.assertEqual(value['reported_status'],status)
                self.assertEqual(value['supporting_checks'],payload['supporting_checks'])
                self.assertEqual(value['supporting_counts'],{s:int(s==child) for s in STATES})
                if child in ('FAIL','ERROR'):
                    stream=io.StringIO()
                    with contextlib.redirect_stdout(stream):
                        global_runner.print_summary({'executed':1,'PASS':0,'FAIL':0,'SKIPPED':0,'ERROR':1},[value])
                    text=stream.getvalue();row=payload['supporting_checks'][0]
                    for field in ('check','runner','log'):self.assertIn(row[field],text)

    def test_incorrect_missing_boolean_or_extra_counters_are_error(self):
        for counts in [{'PASS':2,'FAIL':0,'ERROR':0,'SKIPPED':0},
                       {'PASS':True,'FAIL':0,'ERROR':0,'SKIPPED':0},
                       {'PASS':1}, {'PASS':1,'FAIL':0,'ERROR':0,'SKIPPED':0,'OTHER':0}]:
            with self.subTest(counts=counts):
                value,payload=self.load('SKIPPED',['PASS'],counts)
                self.assertEqual(value['status'],'ERROR')
                self.assertEqual(value['reported_supporting_counts'],counts)
                self.assertEqual(value['supporting_counts'],{'PASS':1,'FAIL':0,'ERROR':0,'SKIPPED':0})

    def test_exit_contradiction_still_retains_actual_child_failure(self):
        value,payload=self.load('SKIPPED',['FAIL'])
        value=global_runner.load_runner_result('U15',Path(value['runner']),
            Path(value['result_json']),Path(value['log']),1)
        self.assertEqual(value['status'],'ERROR')
        self.assertIn('exit',value['reason'])
        self.assertEqual(value['supporting_checks'],payload['supporting_checks'])
        self.assertEqual(value['supporting_counts'],{'PASS':0,'FAIL':1,'ERROR':0,'SKIPPED':0})

    def test_bad_row_fields_status_types_and_duplicates_are_error(self):
        good={'check':'CONTROL','status':'PASS','runner':'CONTROL.py','log':'CONTROL.log'}
        cases=[[dict(good,status=[])],[dict(good,status='UNKNOWN')],
               [{k:v for k,v in good.items() if k!='runner'}],
               [dict(good,log='')],[good,dict(good)]]
        for rows in cases:
            with self.subTest(rows=rows):
                value,payload=self.load('SKIPPED',rows)
                self.assertEqual(value['status'],'ERROR')
                self.assertEqual(value['reported_supporting_checks'],rows)


def main():
    global EVIDENCE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if args.output_dir:
        EVIDENCE=args.output_dir.resolve();EVIDENCE.mkdir(parents=True,exist_ok=False)
    else:
        EVIDENCE=Path(tempfile.mkdtemp(prefix='amps-m3-supporting-protocol-'))
    with (EVIDENCE/'checks.log').open('w') as log:
        result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SupportingProtocolTests))
    status='ERROR' if result.errors else 'FAIL' if result.failures else 'SKIPPED' if result.skipped else 'PASS'
    value={'status':status,'scope':'orchestration controls only; no native physics/observations',
        'executed':result.testsRun-len(result.skipped),'passed':result.testsRun-len(result.skipped)-len(result.failures)-len(result.errors),
        'failed':len(result.failures),'errors':len(result.errors),'skipped':len(result.skipped),'log':str(EVIDENCE/'checks.log')}
    (EVIDENCE/'result.json').write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]

if __name__=='__main__':raise SystemExit(main())
