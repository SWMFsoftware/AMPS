#!/usr/bin/env python3
"""Production RNG and independent qualification-software controls only.

All cohort/profile/predecessor fixtures here are explicitly mathematical or
administrative. They are NEVER observations, real state profiles or formal
qualification evidence. The actual production RNG is compiled independently.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from decimal import Decimal,localcontext
import hashlib
import hmac
import io
import json
import math
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import unittest

import m3_qualification as q
import m3_qualification_executor as executor

ROOT=Path(__file__).resolve().parents[2]
CONTEXT=None
RNG_SAMPLES=50000
RNG_BIN_CLAIMS=10
RNG_FAMILY_ALPHA=.01
ARITHMETIC_BOUND=512*sys.float_info.epsilon


def python_reference(seed,parent,purpose,index,cohort,genealogy):
    u64=lambda value:struct.pack('>Q',value)
    text=cohort.encode('utf-8')
    message=b'AMPS-MOON-M3-HMAC-SHA256-COUNTER-v1'+u64(len(text))+text+u64(parent)+u64(len(genealogy))
    message+=b''.join(u64(value) for value in genealogy)+u64(purpose)+u64(index)
    digest=hmac.new(u64(seed),message,hashlib.sha256).hexdigest()
    return digest,(int(digest[:13],16)+.5)*2**-52


def fixture_plan(n=8000,identity='I22'):
    return {'schema':'srcMoon-M3-cohort-plan-v1','approved_addendum_sha256':q.APPROVED_SHA256,
        'test_id':identity,'cohort_id':'ANALYTICAL-MECHANICS-FIXTURE-NOT-OBSERVATIONS',
        'original_parents':n,'initial_represented_number':.125,'horizon_s':math.log(2),
        'reference_total_rate_s':1.,'source_channel_rates_s':[{'id':'analytic-a','rate_s':.25},{'id':'analytic-b','rate_s':.75}],
        'branch_claims':2,'seed':731993,'primary_estimator':'rate' if identity=='I22' else 'lifetime',
        'source_profile':{'path':'MISSING_REAL_PHYSICAL_PROFILE.json','sha256':'0'*64}}


def write_cohort(path,plan,survivors=None,event_time=None):
    n=plan['original_parents'];survivors=n//2 if survivors is None else survivors
    removed=n-survivors;boundary=removed//4
    with Path(path).open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(q.PARENT_COLUMNS)
        for i in range(n):
            survived=i<survivors
            writer.writerow([i,'SURVIVED' if survived else 'REACTED','' if survived else
                format(plan['horizon_s']/2 if event_time is None else event_time,'.17g'),
                '' if survived else 'analytic-a' if i-survivors<boundary else 'analytic-b',.125])


class RandomControls(unittest.TestCase):
    def test_independent_hmac_vectors_extreme_semantic_keys(self):
        cases=[]
        for seed in (0,731993,104729,130363,2**64-1):
            for parent in (0,1,2**64-1):
                for purpose in (1,2,3,4):
                    for genealogy in ([],[0],[1,0],[0,1],[2**64-1,7]):
                        cases.append((seed,parent,purpose,2**64-1,'explicit-α-cohort',genealogy))
        rows=CONTEXT.query(cases,'independent-vectors')
        for case,row in zip(cases,rows):
            expected,u=python_reference(*case)
            self.assertEqual(row[:2],[1,1]);self.assertEqual(row[2],expected);self.assertEqual(row[3],u)

    def test_retry_and_shard_order_preserve_same_original_parent(self):
        cases=[(731993,parent,1,index,'matched-dt-and-MPI-cohort',[7,2]) for parent in range(24) for index in range(3)]
        normal=CONTEXT.query(cases,'normal-order')
        order=list(range(len(cases)-1,-1,-1));reverse=CONTEXT.query([cases[i] for i in order],'reverse-order')
        self.assertEqual(reverse,[normal[i] for i in order])
        for rank in range(8):
            selected=[i for i,case in enumerate(cases) if case[1]%8==rank]
            self.assertEqual(CONTEXT.query([cases[i] for i in selected],'logical-shard-'+str(rank)),[normal[i] for i in selected])
        self.assertEqual(CONTEXT.query(cases,'unchanged-retry'),normal)

    def test_sibling_purpose_counter_and_cohort_keys_are_distinct(self):
        base=(731993,3,1,0,'cohort',[])
        cases=[base,(731993,3,1,0,'cohort',[0]),(731993,3,1,0,'cohort',[1]),
            (731993,3,2,0,'cohort',[]),(731993,3,1,1,'cohort',[]),(731993,3,1,0,'another-cohort',[])]
        rows=CONTEXT.query(cases,'semantic-domains')
        self.assertEqual(len({row[2] for row in rows}),len(rows))

    def test_invalid_key_preserves_digest_and_uniform_outputs(self):
        cases=[(731993,0,0,0,'cohort',[]),(731993,0,5,0,'cohort',[]),(731993,0,1,0,'',[])]
        for row in CONTEXT.query(cases,'invalid-key'):
            self.assertEqual(row[:2],[0,0]);self.assertEqual(row[2],'UNCHANGED');self.assertEqual(row[3],17.125);self.assertTrue(row[4])

    def test_exact_finite_grid_endpoints_are_strictly_open(self):
        minimum=(0+.5)*2**-52;maximum=((2**52-1)+.5)*2**-52
        self.assertEqual(minimum,2**-53);self.assertEqual(maximum,1-2**-53)
        self.assertGreater(minimum,0);self.assertLess(maximum,1)
        tail=q.finite_grid_tail_bound(RNG_SAMPLES)
        self.assertEqual(tail['omitted_tail_union_upper'],RNG_SAMPLES*2**-52)
        self.assertEqual(tail['maximum_optical_depth'],-math.log(minimum))

    def test_predeclared_distribution_bins_and_reference_plot(self):
        try:from scipy.stats import binom
        except ImportError:self.skipTest('SciPy fixed exact-binomial RNG controls unavailable')
        cases=[(731993,parent,1,0,'predeclared-rng-bin-control',[]) for parent in range(RNG_SAMPLES)]
        rows=CONTEXT.query(cases,'distribution')
        counts=Counter(min(RNG_BIN_CLAIMS-1,int(row[3]*RNG_BIN_CLAIMS)) for row in rows)
        tail=RNG_FAMILY_ALPHA/(2*RNG_BIN_CLAIMS)
        low=int(binom.ppf(tail,RNG_SAMPLES,1/RNG_BIN_CLAIMS));high=int(binom.isf(tail,RNG_SAMPLES,1/RNG_BIN_CLAIMS))
        CONTEXT.bins=[{'bin':i,'actual_count':counts[i],'independent_expected_count':RNG_SAMPLES/RNG_BIN_CLAIMS,
            'frozen_lower_count':low,'frozen_upper_count':high} for i in range(RNG_BIN_CLAIMS)]
        for row in CONTEXT.bins:self.assertLessEqual(low,row['actual_count']);self.assertLessEqual(row['actual_count'],high)
        self.assertTrue(all(0<row[3]<1 for row in rows))


class StatisticalControls(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='analytical-qualification-',dir=CONTEXT.output)
        self.addCleanup(self.temporary.cleanup);self.directory=Path(self.temporary.name)
        self.plan=fixture_plan();self.path=self.directory/'ANALYTICAL_NOT_OBSERVATIONS.csv';write_cohort(self.path,self.plan)

    def test_approved_addendum_and_original_bytes_are_bound(self):
        self.assertEqual(q.approved_protocol()['approval']['status'],'APPROVED')
        self.assertEqual(q.sha(q.TEST/'M3_QUALIFICATION_PROTOCOL.json'),q.BASE_SHA256)

    def test_changed_addendum_cannot_selfapprove_another_gate(self):
        path=self.directory/'changed-addendum.json';path.write_bytes(q.ADDENDUM.read_bytes()+b'\n')
        with self.assertRaises(ValueError):q.approved_protocol(path)

    def test_analytical_endpoint_and_complete_branch_scores(self):
        result=q.score_cohort(self.path,self.plan)
        self.assertEqual(result['status'],'PASS');self.assertEqual(result['metrics']['endpoint_rate_s'],1)
        self.assertEqual(result['metrics']['first_removals'],4000)
        self.assertEqual([r['count'] for r in result['branches']],[1000,3000])
        self.assertIn('NOT QUALIFIED',result['formal_status'])

    def test_event_time_diagnostic_cannot_replace_failing_primary(self):
        survivors=3200;removed=4800
        event=(removed-survivors*self.plan['horizon_s'])/removed
        write_cohort(self.path,self.plan,survivors,event)
        result=q.score_cohort(self.path,self.plan)
        self.assertEqual(result['primary_status'],'FAIL')
        self.assertLess(abs(result['metrics']['event_time_rate_diagnostic_s']-1),ARITHMETIC_BOUND)

    def test_rate_and_lifetime_planning_have_different_correct_intervals(self):
        for estimator in ('rate','lifetime'):
            row=q.survival_plan(1,math.log(2),.02,48,.01,estimator)
            n=row['independent_original_parents'];p=row['reference_survival_probability'];lo,hi=row['accepted_survival_interval']
            margin=min(p-lo,hi-p)
            self.assertEqual(n%8,0);self.assertLessEqual(2*48*math.exp(-2*n*margin*margin),.01)
            self.assertGreater(2*48*math.exp(-2*(n-8)*margin*margin),.01)
            expected=[math.exp(-math.log(2)*(1+.02)),math.exp(-math.log(2)*(1-.02))] if estimator=='rate' else [math.exp(-math.log(2)/(1-.02)),math.exp(-math.log(2)/(1+.02))]
            self.assertEqual([lo,hi],expected)

    def test_planning_uses_actual_represented_exposure(self):
        a=q.survival_plan(.4,1.1,.03,1,.01,'lifetime')
        self.assertEqual(a['represented_reference_cumulative_hazard'],.4*1.1)
        self.assertEqual(a['reference_survival_probability'],math.exp(-.4*1.1))

    def test_exact_clopper_pearson_against_decimal_binomial_polynomial(self):
        alpha=.05/17;n=11;k=3;low,high=q.clopper_pearson(k,n,alpha)
        with localcontext() as context:
            context.prec=70
            def beta_integer(x,a,b):
                x=Decimal.from_float(x);m=a+b-1
                return sum(Decimal(math.comb(m,j))*x**j*(1-x)**(m-j) for j in range(a,m+1))
            self.assertLess(abs(float(beta_integer(low,k,n-k+1))-alpha/2),ARITHMETIC_BOUND)
            self.assertLess(abs(float(beta_integer(high,k+1,n-k))-(1-alpha/2)),ARITHMETIC_BOUND)
        self.assertEqual(q.clopper_pearson(0,10,.05)[0],0)
        self.assertEqual(q.clopper_pearson(10,10,.05)[1],1)

    def test_rare_precision_is_not_replaced_by_confidence_containment(self):
        n=20;plan=fixture_plan(n);write_cohort(self.path,plan)
        result=q.score_cohort(self.path,plan)
        self.assertEqual(result['primary_status'],'PASS');self.assertEqual(result['status'],'FAIL')
        self.assertTrue(any(r['relative_halfwidth']>.1 for r in result['branches']))

    def test_zero_source_branch_is_not_assigned_fake_relative_error(self):
        self.plan['source_channel_rates_s'].append({'id':'source-zero','rate_s':0.})
        self.plan['branch_claims']=3
        result=q.score_cohort(self.path,self.plan)
        row=result['branches'][-1];self.assertEqual(row['status'],'PASS');self.assertIsNone(row['relative_halfwidth'])

    def test_exact_precision_plan_and_explicit_work_refusal(self):
        plan=q.branch_precision_plan(.25,4000,2)
        result=q.certify_branch_precision(plan,10000)
        self.assertEqual(result['status'],'PASS')
        with self.assertRaises(q.MissingPrerequisite):q.certify_branch_precision(plan,1)

    def test_complete_fixed_cohort_plan_preserves_rare_branches_and_bias_budget(self):
        plan=q.complete_cohort_plan('I22',1,math.log(2),[.001,.999],48,96,0,100000000,20)
        self.assertEqual(plan['status'],'PASS');self.assertTrue(plan['history'][-1]['adequate'])
        self.assertEqual(plan['unchanged_final_relative_gate'],.02)
        self.assertAlmostEqual((1+.002)*(1+plan['stochastic_relative_planning_margin'])-1,.02)
        self.assertEqual(len(plan['history'][-1]['relative_halfwidth_upper_bounds']),2)
        limited=q.complete_cohort_plan('I22',1,math.log(2),[.001,.999],48,96,0,8,20)
        self.assertEqual(limited['status'],'SKIPPED')

    def test_source_probability_roundoff_is_preserved_not_renormalized(self):
        probabilities=[.1,.2,.3,.4-2**-53]
        result=q.complete_cohort_plan('I22',1,math.log(2),probabilities,1,4,0,10000000,20)
        self.assertEqual(result['status'],'PASS');self.assertEqual(result['source_probabilities'],probabilities)
        self.assertEqual(result['source_probability_sum'],math.fsum(probabilities))
        self.assertFalse(result['source_renormalization'])
        with self.assertRaises(ValueError):q.complete_cohort_plan('I22',1,math.log(2),[.2,.799],1,2,0,10000000,20)

    def test_duplicate_parent_unequal_weight_and_censoring_refused(self):
        with self.path.open(newline='') as stream:rows=list(csv.reader(stream))
        for index,value in ((0,rows[1][0]),(1,'ESCAPED'),(4,'0.25')):
            changed=[r[:] for r in rows];changed[2][index]=value
            with self.path.open('w',newline='') as stream:csv.writer(stream).writerows(changed)
            with self.assertRaises(ValueError):q.score_cohort(self.path,self.plan)

    def test_no_channel_drop_renormalization_or_estimator_switch(self):
        for key,value in (('reference_total_rate_s',2),('primary_estimator','lifetime'),('branch_claims',1)):
            copied=dict(self.plan);copied[key]=value
            with self.assertRaises(ValueError):q.validate_plan(copied)

    def test_empty_or_all_surviving_endpoint_is_evaluated_fail(self):
        for survived in (0,self.plan['original_parents']):
            write_cohort(self.path,self.plan,survived)
            self.assertEqual(q.score_cohort(self.path,self.plan)['primary_status'],'FAIL')

    def test_zero_hazard_uses_exact_invariant_without_fake_finite_lifetime(self):
        self.plan.update(test_kind='zero-hazard-invariant',primary_estimator='zero-hazard',reference_total_rate_s=0.)
        for row in self.plan['source_channel_rates_s']:row['rate_s']=0.
        write_cohort(self.path,self.plan,self.plan['original_parents'])
        result=q.score_cohort(self.path,self.plan)
        self.assertEqual(result['status'],'PASS');self.assertNotIn('endpoint_lifetime_s',result['metrics'])
        write_cohort(self.path,self.plan,self.plan['original_parents']-4)
        self.assertEqual(q.score_cohort(self.path,self.plan)['status'],'FAIL')

    def test_numerical_budget_requires_more_than_decreasing_errors(self):
        rows=[{'steps':n,'actual':x,'reference':0,'scale':1} for n,x in zip((8,16,32),(.1,.05,.01))]
        self.assertEqual(q.score_refinement('I21',rows)['status'],'FAIL')
        rows[-1]['actual']=.001;self.assertEqual(q.score_refinement('I21',rows)['status'],'PASS')
        for row in rows:row['actual']=0
        self.assertEqual(q.score_refinement('I21',rows)['status'],'PASS')

    def test_independent_reference_precision_cannot_be_lost_before_subtraction(self):
        rows=[{'steps':n,'actual':str(10**20+error),'reference':str(10**20),'scale':1}
              for n,error in zip((8,16,32),(3,2,1))]
        result=q.score_refinement('I21',rows)
        self.assertEqual(result['status'],'FAIL')
        self.assertEqual(result['relative_errors'],[3.,2.,1.])

    def test_invalid_numeric_plans_and_statistical_counts_refused(self):
        for value in (True,0,-1,float('nan'),float('inf')):
            copied=dict(self.plan);copied['horizon_s']=value
            with self.assertRaises(ValueError):q.validate_plan(copied)
        with self.assertRaises(ValueError):q.clopper_pearson(11,10,.05)
        with self.assertRaises(ValueError):q.clopper_pearson(True,10,.05)
        with self.assertRaises(ValueError):q.finite('-1e-9999','actual event time')

    def test_unavailable_producer_remains_skip_even_with_green_fixture_gates(self):
        # Administrative controls ONLY; none is submitted as science evidence.
        profile=self.directory/'ADMINISTRATIVE_PROFILE_NOT_SCIENCE.json'
        profile.write_text(json.dumps({'schema':'srcMoon-M3-state-profile-v1','physical_activation':False}))
        self.plan['source_profile']={'path':str(profile),'sha256':q.sha(profile)}
        path=self.directory/'ANALYTICAL_PLAN.json';path.write_text(json.dumps(self.plan))
        for identity in tuple('I%02d'%i for i in range(13,21))+('U15','U24'):
            directory=self.directory/identity;directory.mkdir();(directory/'result.json').write_text('{"status":"PASS","scope":"ADMINISTRATIVE FIXTURE NOT QUALIFICATION"}')
        result=executor.run(path,self.directory,self.directory/'attempt')
        self.assertEqual(result['status'],'SKIPPED');self.assertFalse(result['criteria_evaluated'])
        self.assertTrue(any('not installed' in message for message in result['blocking_reasons']))
        self.assertFalse((self.directory/'attempt/receipt.txt').exists())

    def test_corruption_cannot_hide_behind_missing_producer(self):
        self.plan['source_profile']={'path':str(self.path),'sha256':'0'*64}
        path=self.directory/'plan.json';path.write_text(json.dumps(self.plan))
        result=executor.run(path,self.directory,self.directory/'bad-attempt')
        self.assertEqual(result['status'],'ERROR')
        for name in ('result.json','run_manifest.json','stdout.log','stderr.log','cohort_plan.json'):
            self.assertTrue((self.directory/'bad-attempt'/name).is_file())

    def administrative_native_evidence(self):
        # THESE ARE ADMINISTRATIVE RECEIPT FIXTURES, NOT A NATIVE RUN. The
        # positive guard below proves that their Boolean/string claims cannot
        # register an actual physical producer or generate qualification PASS.
        profile=self.directory/'ADMINISTRATIVE_PROFILE.json'
        profile.write_text(json.dumps({'schema':'srcMoon-M3-state-profile-v1','physical_activation':False}))
        self.plan['source_profile']={'path':str(profile),'sha256':q.sha(profile)}
        plan=self.directory/'ADMINISTRATIVE_PLAN.json';plan.write_text(json.dumps(self.plan))
        value={'schema':'srcMoon-M3-native-cohort-evidence-v1','test_id':'I22',
            'plan_sha256':q.sha(plan),'approved_addendum_sha256':q.APPROVED_SHA256,
            'source_profile_sha256':q.sha(profile),'producer_id':'ADMINISTRATIVE-NOT-A-PHYSICAL-PRODUCER',
            'execution_token':'ADMINISTRATIVE-NOT-A-LAUNCH','clock':'common-time','mpi_ranks':8,
            'openmp_threads':1,'native_returncode':0,'artifacts':{}}
        for role in executor.ARTIFACT_ROLES:
            path=self.directory/('ADMINISTRATIVE-'+role+'.txt')
            if role=='run_manifest':
                document={key:value[key] for key in ('execution_token','test_id','plan_sha256','source_profile_sha256')}
                document.update(execution='native-production-cohort',scope='ADMINISTRATIVE FIXTURE, NOT NATIVE EXECUTION')
                path.write_text(json.dumps(document))
            else:path.write_text('ADMINISTRATIVE PATH/HASH FIXTURE; NOT SCIENCE, NATIVE DATA OR EXECUTABLE\n')
            value['artifacts'][role]={'path':str(path),'sha256':q.sha(path)}
        path=self.directory/'ADMINISTRATIVE-EVIDENCE.json';path.write_text(json.dumps(value))
        return plan,path,value

    def test_hash_bound_receipt_metadata_cannot_install_a_physical_producer(self):
        plan,path,value=self.administrative_native_evidence()
        result=executor.run(plan,self.directory,self.directory/'artifact-attempt',native_evidence=path,expected_evidence_sha256=q.sha(path))
        self.assertEqual(result['status'],'SKIPPED');self.assertFalse(result['criteria_evaluated'])
        self.assertEqual(len(result['identities']),1+len(executor.ARTIFACT_ROLES))

    def test_changed_actual_artifact_is_error_even_if_another_is_missing(self):
        plan,path,value=self.administrative_native_evidence()
        Path(value['artifacts']['event_ledger']['path']).unlink()
        Path(value['artifacts']['input_snapshot']['path']).write_text('CHANGED ADMINISTRATIVE BYTES')
        result=executor.run(plan,self.directory,self.directory/'corrupt-artifact',native_evidence=path,expected_evidence_sha256=q.sha(path))
        self.assertEqual(result['status'],'ERROR')

    def test_run_manifest_token_contradiction_is_error(self):
        plan,path,value=self.administrative_native_evidence()
        manifest=Path(value['artifacts']['run_manifest']['path']);record=json.loads(manifest.read_text())
        record['execution_token']='WRONG ADMINISTRATIVE TOKEN';manifest.write_text(json.dumps(record))
        value['artifacts']['run_manifest']['sha256']=q.sha(manifest);path.write_text(json.dumps(value))
        result=executor.run(plan,self.directory,self.directory/'wrong-token',native_evidence=path,expected_evidence_sha256=q.sha(path))
        self.assertEqual(result['status'],'ERROR')


class Context:
    def __init__(self,output):
        self.output=output;self.commands=[];self.bins=[];self.executable=output/'qualification-random-probe'
        compiler=shutil.which('g++')
        if not compiler:raise q.MissingPrerequisite('C++17 compiler unavailable')
        self.sources=[ROOT/p for p in ('srcMoon/species/QualificationRandom.h','srcMoon/species/QualificationRandom.cpp',
            'src/models/chemical_database/DataIO.h','src/models/chemical_database/DataIO.cpp',
            'srcMoon/test/qualification/random_probe.cpp','srcMoon/test/m3_qualification.py',
            'srcMoon/test/m3_qualification_executor.py','srcMoon/test/test_m3_qualification.py',
            'srcMoon/test/M3_QUALIFICATION_ADDENDUM_20261118.json','srcMoon/test/qualification/CONTRACT.md')]
        self.hashes={str(path):q.sha(path) for path in self.sources}
        command=[compiler,'-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic',
            str(ROOT/'srcMoon/test/qualification/random_probe.cpp'),str(ROOT/'srcMoon/species/QualificationRandom.cpp'),
            str(ROOT/'src/models/chemical_database/DataIO.cpp'),'-o',str(self.executable)]
        result=subprocess.run(command,capture_output=True,text=True,timeout=180)
        (output/'compile.stdout.log').write_text(result.stdout);(output/'compile.stderr.log').write_text(result.stderr)
        self.commands.append({'command':command,'returncode':result.returncode})
        if result.returncode:raise RuntimeError('actual production RNG compilation failed')
        self.compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]

    def query(self,cases,label):
        def quote(value):return '"'+value.replace('\\','\\\\').replace('"','\\"')+'"'
        lines=[' '.join([str(seed),str(parent),str(purpose),str(index),quote(cohort),str(len(genealogy)),*(str(v) for v in genealogy)])
            for seed,parent,purpose,index,cohort,genealogy in cases]
        directory=self.output/'queries'/('%03d-'%len(self.commands)+label);directory.mkdir(parents=True)
        (directory/'stdin.txt').write_text('\n'.join(lines)+'\n')
        process=subprocess.run([str(self.executable)],input='\n'.join(lines)+'\n',capture_output=True,text=True,timeout=120)
        (directory/'stdout.log').write_text(process.stdout);(directory/'stderr.log').write_text(process.stderr)
        self.commands.append({'command':[str(self.executable)],'label':label,'rows':len(cases),'returncode':process.returncode,
                              'log':str(directory/'stdout.log')})
        if process.returncode:raise RuntimeError('production random probe abnormal exit')
        parsed=[shlex.split(row) for row in process.stdout.splitlines()]
        if len(parsed)!=len(cases) or any(len(row)!=5 for row in parsed):raise ValueError('native RNG row count/schema mismatch')
        return [[int(row[0]),int(row[1]),row[2],float(row[3]),row[4]] for row in parsed]


class Recorded(unittest.TextTestResult):
    def __init__(self,*a,**kw):super().__init__(*a,**kw);self.records=[]
    def addSuccess(self,test):super().addSuccess(test);self.records.append({'test':test.id(),'status':'PASS'})
    def addFailure(self,test,error):super().addFailure(test,error);self.records.append({'test':test.id(),'status':'FAIL'})
    def addError(self,test,error):super().addError(test,error);self.records.append({'test':test.id(),'status':'ERROR'})
    def addSkip(self,test,reason):super().addSkip(test,reason);self.records.append({'test':test.id(),'status':'SKIPPED','reason':reason})


def main():
    global CONTEXT
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();records=[];manifest={'scope':'analytical/software controls only; NOT observations/physical qualification'}
    try:
        CONTEXT=Context(output)
        suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c) for c in (RandomControls,StatisticalControls))
        with (output/'checks.log').open('w') as stream:
            result=unittest.TextTestRunner(stream=stream,verbosity=2,resultclass=Recorded).run(suite)
        records=result.records
        if any(q.sha(path)!=expected for path,expected in CONTEXT.hashes.items()):raise RuntimeError('qualification production/test source changed during run')
        manifest.update(source_sha256=CONTEXT.hashes,compiler=CONTEXT.compiler,commands=CONTEXT.commands,
            executable_sha256=q.sha(CONTEXT.executable),rng_samples=RNG_SAMPLES,rng_family_alpha=RNG_FAMILY_ALPHA)
        if CONTEXT.bins:
            with (output/'comparison.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(CONTEXT.bins[0]));writer.writeheader();writer.writerows(CONTEXT.bins)
            try:
                import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
                fig,ax=plt.subplots(figsize=(8,4));ax.bar([r['bin'] for r in CONTEXT.bins],[r['actual_count'] for r in CONTEXT.bins],label='actual production stream')
                ax.axhline(RNG_SAMPLES/RNG_BIN_CLAIMS,color='black',label='independent uniform reference')
                ax.set(xlabel='Fixed uniform bin',ylabel='Count',title='M3 semantic RNG software control — no observations');ax.legend()
                for extension in ('png','eps'):fig.savefig(output/('comparison.'+extension),bbox_inches='tight')
                plt.close(fig)
            except ImportError:manifest['visualization']={'status':'SKIPPED','reason':'optional Matplotlib unavailable; numerical CSV retained'}
    except q.MissingPrerequisite as error:records.append({'test':'infrastructure-prerequisite','status':'SKIPPED','reason':str(error)})
    except Exception as error:
        import traceback
        (output/'runner-error.log').write_text(traceback.format_exc());records.append({'test':'infrastructure','status':'ERROR','reason':str(error)})
    counts={status:sum(row['status']==status for row in records) for status in ('PASS','FAIL','ERROR','SKIPPED')}
    status='ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'SKIPPED' if counts['SKIPPED'] else 'PASS'
    value={'schema':'srcMoon-M3-qualification-software-controls-v1','status':status,'counts':counts,'checks':records,
        'scope':'production RNG and analytical qualification infrastructure only','formal_M3_status':'SKIPPED / NOT QUALIFIED',
        'observations_used':False,'native_physical_cohort_executed':False,'elapsed_seconds':time.monotonic()-started}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps(value,indent=2));return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]


if __name__=='__main__':raise SystemExit(main())
