#!/usr/bin/env python3
"""Protocol controls only: no physical events, source repairs or observations."""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import m3_protocol as protocol

class QualificationPreparationTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory(prefix='amps-m3-protocol-control-')
        self.addCleanup(temporary.cleanup);self.directory=Path(temporary.name)
        self.p=json.loads((protocol.TEST/'M3_QUALIFICATION_PROTOCOL.json').read_text())
        # Typed temporary metadata controls only. Refer to actual unchanged
        # specifications/contracts; no temporary file is a scientific oracle.
        self.p['physical_model_specification']['path']=str(protocol.TEST/self.p['physical_model_specification']['path'])
        for case in self.p['cases'].values():case['acceptance_path']=str(protocol.TEST/case['acceptance_path'])

    def evaluate(self):
        path=self.directory/'protocol.json';path.write_text(json.dumps(self.p))
        return protocol.validate(path)

    def test_actual_preparation_is_skip_with_all_four_physical_blocks(self):
        value=protocol.validate()
        self.assertEqual((value['status'],value['structural_status']),('SKIPPED','PASS'))
        self.assertEqual(set(value['cases']),set(protocol.IDS))
        self.assertTrue(all(r['status']=='SKIPPED' for r in value['cases'].values()))
        self.assertFalse(value['criteria_evaluated']);self.assertEqual(value['execution'],'not-run')

    def test_changed_threshold_is_error_not_relaxed_acceptance(self):
        for key,criteria in protocol.NUMERIC.items():
            for name in criteria:
                with self.subTest(key=key,name=name):
                    saved=self.p['cases'][key]['frozen_criteria'][name]
                    self.p['cases'][key]['frozen_criteria'][name]=saved+.01
                    self.assertEqual(self.evaluate()['status'],'ERROR')
                    self.p['cases'][key]['frozen_criteria'][name]=saved

    def test_changed_identity_is_error(self):
        self.p['physical_model_specification']['sha256']='0'*64
        self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_missing_evidence_remains_skip(self):
        self.p['physical_model_specification']['path']=str(self.directory/'missing.json')
        value=self.evaluate();self.assertEqual(value['status'],'SKIPPED')
        self.assertTrue(any('unavailable' in r for r in value['blocking_reasons']))

    def test_missing_or_repurposed_id_is_error(self):
        self.p['cases'].pop('I24');self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_fabricated_executor_boolean_is_error(self):
        self.p['formal_executor_installed']=True;self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_fabricated_physical_run_is_error(self):
        self.p['cases']['I21']['qualification_execution']='PASS';self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_metadata_cannot_approve_proposed_statistics(self):
        self.p['statistical_design']['status']='APPROVED';self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_frozen_contract_corruption_is_error(self):
        self.p['cases']['I22']['acceptance_sha256']='1'*64;self.assertEqual(self.evaluate()['status'],'ERROR')

    def test_clock_layout_resource_contracts_enforced(self):
        for obj,key,value in [(self.p['clock'],'mode','independent-history'),
            (self.p['execution'],'development_wall_seconds',7200),
            (self.p['execution'],'mpi_ranks',[8]),(self.p['execution'],'openmp_threads',2)]:
            with self.subTest(key=key):
                saved=copy.deepcopy(obj[key]);obj[key]=value
                self.assertEqual(self.evaluate()['status'],'ERROR');obj[key]=saved

    def test_hoeffding_design_satisfies_independent_bound(self):
        for d in (.02,.03):
            for k in (1,48,1000):
                row=protocol.cohort_bound(d,k,.01);n=row['global_original_parents']
                low,high=row['acceptance_probability_interval'];margin=min(.5-low,high-.5)
                self.assertEqual(n%8,0)
                self.assertLessEqual(2*k*math.exp(-2*n*margin*margin),.01)
                self.assertGreater(2*k*math.exp(-2*(n-8)*margin*margin),.01)

    def test_planning_rejects_invalid_family_sizes(self):
        for k in (0,-1,True,1.5):
            with self.subTest(k=k),self.assertRaises(ValueError):protocol.cohort_bound(.02,k,.01)

    def test_finer_precision_or_larger_family_needs_more_trials(self):
        n=lambda d,k:protocol.cohort_bound(d,k,.01)['global_original_parents']
        self.assertGreater(n(.02,48),n(.03,48));self.assertGreater(n(.02,48),n(.02,1))

    def test_malformed_json_returns_retained_error_artifacts(self):
        for index,raw in enumerate(('[]','null','{"schema":1,"schema":2}',
            '{"risk":NaN}','{"risk":Infinity}')):
            with self.subTest(raw=raw):
                path=self.directory/('bad-'+str(index)+'.json');path.write_text(raw)
                output=self.directory/('error-'+str(index))
                value=protocol.run(output,path)
                self.assertEqual(value['status'],'ERROR')
                for name in ('result.json','run_manifest.json','stdout.log','stderr.log'):
                    self.assertTrue((output/name).is_file())

    def test_boolean_layout_nonfinite_risk_and_duplicate_seeds_are_errors(self):
        for obj,key,value in ((self.p['execution'],'openmp_threads',True),
            (self.p['execution'],'mpi_ranks',[True,8]),
            (self.p['statistical_design'],'family_planning_risk',float('nan')),
            (self.p['statistical_design'],'family_planning_risk',True),
            (self.p['statistical_design'],'seeds',[731993,731993])):
            with self.subTest(key=key,value=value):
                saved=copy.deepcopy(obj[key]);obj[key]=value
                self.assertEqual(self.evaluate()['status'],'ERROR');obj[key]=saved

    def test_exposure_receipt_requires_every_declared_control_and_matching_exit(self):
        # Metadata transcript fixtures only, never executable/scientific proof.
        path=protocol.TEST/'reaction_execution/exposure_bridge/run.py'
        spec=importlib.util.spec_from_file_location('exposure_receipt_controls',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        good=''.join('PASS '+name+'\n' for name in sorted(module.EXPECTED))
        suffix='SUMMARY checks='+str(len(module.EXPECTED))+' failed=0\n'
        rows,count,failed=module.parse_receipt(good+suffix,0)
        self.assertEqual(count,len(module.EXPECTED));self.assertEqual(failed,0)
        for text,code in (('SUMMARY checks=0 failed=0\n',0),(good+suffix,1),
            (good.split('\n',1)[1]+suffix,0),(good+good.splitlines()[0]+'\n'+suffix,0)):
            with self.subTest(text=text,code=code),self.assertRaises(ValueError):module.parse_receipt(text,code)

    def test_common_time_refinements_share_integer_step_final_horizon(self):
        self.assertEqual(self.evaluate()['status'],'SKIPPED')
        self.p['numerical_design']['common_steps_to_final_time']=[8,16,31]
        self.assertEqual(self.evaluate()['status'],'ERROR')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=False)
    with (output/'checks.log').open('w') as stream:
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QualificationPreparationTests))
    status='ERROR' if result.errors else 'FAIL' if result.failures else 'SKIPPED' if result.skipped else 'PASS'
    value={'status':status,'scope':'protocol controls only; no physical qualification',
        'counts':{'PASS':result.testsRun-len(result.errors)-len(result.failures)-len(result.skipped),
            'FAIL':len(result.failures),'ERROR':len(result.errors),'SKIPPED':len(result.skipped)}}
    (output/'result.json').write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]

if __name__=='__main__':raise SystemExit(main())
