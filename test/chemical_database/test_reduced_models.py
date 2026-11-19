#!/usr/bin/env python3
"""Compile/evaluate shared explicit-input models against independent arithmetic.

All cases are labelled mathematical fixtures, not observations, species defaults
or qualified physical source profiles. Decimal/Fraction references use independent
mass/energy/CDF/inventory calculations; the adapter calls only shipped C++ APIs.
The numerical comparison budget is fixed at 1e-12 relative to the named quantity
scale before execution; it is not an observational acceptance threshold.
"""
import argparse
import csv
from decimal import Decimal as D, localcontext
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

import test_protocol as common

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=None
EXPECTED_GUARDS=('ion_null_refused','ion_negative_energy_transactional','ion_no_positive_reference_mass',
 'ion_nan_and_infinite_refused','radiative_null_refused','radiative_zero_total_outside_domain',
 'radiative_negative_transactional','radiative_empty_nan_inf_refused','radiative_weak_positive_underflow_not_deleted',
 'opal_null_refused','opal_threshold_transactional','opal_width_and_quantile_refused','opal_nan_refused',
 'opal_zero_endpoint','opal_one_endpoint','two_body_null_refused','two_body_zero_direction_transactional',
 'two_body_negative_COM_no_clamp','two_body_invalid_masses','two_body_invalid_bounds','two_body_nonfinite_refused',
 'two_body_exact_zero_COM','two_body_lost_positive_mass_refused','two_body_underflowed_bulk_not_certified_zero',
 'inherit_null_refused','inherit_invalid_mass_transactional','inherit_nonfinite_charge_refused',
 'inherit_signed_negative_reservoir','inherit_small_Q_survives_bulk_cancellation','thermal_null_refused',
 'thermal_zero_temperature_transactional','thermal_negative_stock_refused','thermal_nan_inf_refused',
 'thermal_zero_attempt_explicit','thermal_zero_duration','thermal_rate_underflow_exposure_preserved',
 'thermal_prefactor_rescues_rate','thermal_large_exposure_explicit_overflow_diagnostic')


def _atan_series(x):
    total=term=x;power=x*x;n=1
    while abs(term)>D('1e-100'):
        term=-term*power
        total+=term/D(2*n+1);n+=1
        if n>1000:raise RuntimeError('independent Decimal atan did not converge')
    return total


def atan(x):
    # Machin's identity obtains pi independently of the production libm.
    pi=16*_atan_series(D(1)/5)-4*_atan_series(D(1)/239)
    if x>1:return pi/2-atan(1/x)
    multiplier=D(1)
    while x>D('.125'):
        x=x/(1+(1+x*x).sqrt());multiplier*=2
    return multiplier*_atan_series(x)


def tan(x):
    s,st,c,ct=x,x,D(1),D(1);xx=x*x;n=1
    while max(abs(st),abs(ct))>D('1e-100'):
        st=-st*xx/D((2*n)*(2*n+1));ct=-ct*xx/D((2*n-1)*(2*n))
        s+=st;c+=ct;n+=1
        if n>1000:raise RuntimeError('independent Decimal trigonometric reference did not converge')
    return s/c


def loss_fraction(x):
    # Direct series retains tiny losses that 1-Decimal.exp(-x) would cancel.
    if x>D('.1'):return 1-(-x).exp()
    total=term=x;n=1
    while abs(term)>abs(total)*D('1e-100'):
        n+=1;term=-term*x/n;total+=term
        if n>1000:raise RuntimeError('independent finite-inventory series did not converge')
    return total


def fixture_cases():
    cases=[]
    def add(op,name,*values):cases.append({'operation':op,'name':name,'values':[str(v) for v in values]})
    add('ION','ion_si_scale','4e-26','1e-30','1e-18','3e8')
    add('ION','ion_exact_binary_difference','1.00000095367431640625','1','0','1')
    add('ION','ion_separate_small_binding','1e200','1e190','1e200','1e100')
    add('RAD','radiative_closed_set',4,'1','2','0','3')
    add('RAD','radiative_weak_channel',3,'1','1e-20','2')
    for r in ('1e-8','.1','1','10','1e8'):
        for u in ('0','.001','.25','.5','.75','.999','1'):
            add('OPAL','opal_'+r+'_'+u,str(2*D(r)),'0','1',u)
    add('TWO','two_rest',3,5,0,7,0,0,0,1,0,0)
    add('TWO','two_moving',3,5,5,2,1,2,0,0,'.6','.8')
    add('TWO','two_endothermic',2,7,9,-3,1,0,2,0,0,1)
    add('TWO','two_SI_scale','4e-26','1e-30','1e-18','2e-19','1e-23','-2e-23','0','0','.6','.8')
    add('INHERIT','inherit_signed',5,4,0,1,-1,10,-2,3,1,-2,4,7,3)
    add('INHERIT','inherit_SI','4e-26','3.9e-26','0','1e-19','-1e-19',2000,-300,10,'1e-24','0','-2e-24','3e-19','-1e-20')
    for i,vals in enumerate((('2','1','1','1','10','3'),('1e-200','0','1','1','1e200','1e-100'),
        ('1e200','400','1','1','3','1e-100'),('1','0','1','1','10','0'),
        ('0','1','1','1','10','1'),('1','10','1','2','0','3'),
        ('10','0','1','1','10','100'),('1e-20','1','1','100','3e20','1e-5'))):
        add('THERM','thermal_'+str(i),*vals)
    return cases


class ReducedModelControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=OUTPUT or (common.CONTEXT.output/'reduced-models' if common.CONTEXT else Path(tempfile.mkdtemp(prefix='amps-reduced-models-')))
        cls.out.mkdir(parents=True,exist_ok=False) if not cls.out.exists() else None
        cxx=list(common.CONTEXT.cxx) if common.CONTEXT else [shutil.which('g++') or 'g++']
        if not shutil.which(cxx[0]):raise unittest.SkipTest('C++17 compiler unavailable')
        cls.sources=[ROOT/'src/models/chemical_database'/n for n in ('ReducedModels.h','Runtime.h','DataIO.h')]+[Path(__file__).resolve(),Path(__file__).resolve().with_name('reduced_models_probe.cpp')]
        # Preserve actual compiled source bytes, not just paths in a checkout
        # which can subsequently change. The parent shared campaign records
        # these sources too; this also makes the standalone artifact portable.
        for source in cls.sources:
            target=cls.out/'sources'/source.relative_to(ROOT)
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,target)
        cls.exe=cls.out/'reduced-models-probe'
        cmd=cxx+['-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(cls.sources[-1]),'-o',str(cls.exe)]
        (cls.out/'compile_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        compiled=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=180)
        (cls.out/'compile.stdout.log').write_text(compiled.stdout);(cls.out/'compile.stderr.log').write_text(compiled.stderr)
        if compiled.returncode:raise RuntimeError('actual shared reduced-model compilation failed; see '+str(cls.out))
        cls.cases=fixture_cases()
        text=''.join(' '.join([c['operation'],c['name']]+c['values'])+'\n' for c in cls.cases)
        (cls.out/'fixture_input.txt').write_text(text)
        (cls.out/'fixture_input.json').write_text(json.dumps({'scope':'DECLARED ANALYTICAL FIXTURES; NOT OBSERVATIONS','cases':cls.cases},indent=2)+'\n')
        process=subprocess.run([str(cls.exe),str(cls.out/'native_checks.csv'),str(cls.out/'native_metrics.csv')],input=text,text=True,capture_output=True,timeout=60)
        (cls.out/'stdout.log').write_text(process.stdout);(cls.out/'stderr.log').write_text(process.stderr)
        if process.returncode not in (0,1):raise RuntimeError('native adapter infrastructure failure: '+str(process.returncode))
        with(cls.out/'native_checks.csv').open(newline='')as f:cls.guards=list(csv.DictReader(f))
        if tuple(r['name'] for r in cls.guards)!=EXPECTED_GUARDS or any(r['passed'] not in ('0','1')for r in cls.guards):
            raise RuntimeError('native fixed guard membership/output contract mismatch')
        failed=sum(r['passed']=='0'for r in cls.guards)
        if process.stdout.strip()!=f'SUMMARY checks={len(EXPECTED_GUARDS)} failed={failed}' or process.returncode!=(1 if failed else 0):
            raise RuntimeError('native outcome/exit/summary contradiction')
        cls.metrics={};cls.comparisons=[]
        with(cls.out/'native_metrics.csv').open(newline='')as f:
            for r in csv.DictReader(f):
                key=(r['case'],r['key'])
                if key in cls.metrics:raise RuntimeError('duplicate native metric')
                value=D(r['value'])
                if not value.is_finite():raise RuntimeError('nonfinite native metric')
                cls.metrics[key]=value
        for c in cls.cases:
            if (c['name'],'status')not in cls.metrics:raise RuntimeError('missing case outcome')
        identity={'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'compiler':subprocess.check_output(cxx+['--version'],text=True).splitlines()[0],
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in cls.sources},
            'executable_sha256':hashlib.sha256(cls.exe.read_bytes()).hexdigest(),
            'scope':'PIC-independent analytical mathematics; no species defaults, source admission or observations',
            'relative_numerical_budget':'1e-12','native_guard_count':len(cls.guards),'fixture_count':len(cls.cases)}
        (cls.out/'run_manifest.json').write_text(json.dumps(identity,indent=2)+'\n')

    @classmethod
    def tearDownClass(cls):
        if not hasattr(cls,'comparisons'):return
        with(cls.out/'independent_comparisons.csv').open('w',newline='')as f:
            w=csv.DictWriter(f,fieldnames=('case','metric','native','reference','scale','relative_error','passed'));w.writeheader();w.writerows(cls.comparisons)

    def close(self,name,key,reference,scale=None):
        self.assertEqual(self.metrics[(name,'status')],0,name)
        value=self.metrics[(name,key)];reference=D(reference)
        scale=abs(reference)if scale is None else abs(D(scale))
        error=abs(value-reference)
        passed=error<=D('1e-12')*scale if scale else error==0
        self.comparisons.append({'case':name,'metric':key,'native':str(value),'reference':str(reference),'scale':str(scale),
            'relative_error':str(error/scale)if scale else str(error),'passed':int(passed)})
        self.assertTrue(passed,(name,key,value,reference,error,scale))

    def test_fixed_native_domains_transactionality_and_extremes(self):
        self.assertEqual([r['name']for r in self.guards if r['passed']!='1'],[])

    def test_ion_mass_independent_decimal(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='ION':continue
                n,e,q,s=map(D,c['values']);self.close(c['name'],'reference',n-e)
                self.close(c['name'],'binding',q/s/s);self.close(c['name'],'physical',n-e+q/s/s)

    def test_radiative_coherent_set_independent_decimal(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='RAD':continue
                values=list(map(D,c['values'][1:]));total=sum(values)
                self.close(c['name'],'total',total);self.close(c['name'],'lifetime',1/total)
                for i,a in enumerate(values):self.close(c['name'],'branch'+str(i),a/total)

    def test_opal_independent_high_precision_inverse_and_CDF(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='OPAL':continue
                e,i,w,u=map(D,c['values']);maximum=(e-i)/2
                expected=D(0)if u==0 else maximum if u==1 else w*tan(u*atan(maximum/w))
                self.close(c['name'],'energy',expected,maximum)
                self.close(c['name'],'support_native',maximum)
                actual=self.metrics[(c['name'],'energy')]
                # Domain membership is tested against the ACTUAL native
                # binary input support. Decimal '.1' is not exactly the value
                # read as long double; treating that decimal as an exact native
                # boundary would falsely fail its correctly returned endpoint.
                self.assertGreaterEqual(actual,0)
                self.assertLessEqual(actual,self.metrics[(c['name'],'support_native')])
                recovered=atan(actual/w)/atan(maximum/w)
                self.assertLessEqual(abs(recovered-u),D('1e-12'))

    def test_two_body_independent_decimal_momenta_Q_and_energy(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='TWO':continue
                m1,m2,k,q,*rest=map(D,c['values']);p=rest[:3];direction=rest[3:];total=m1+m2
                com=k+q-sum(x*x for x in p)/(2*total);mu=m1*m2/total
                norm=sum(x*x for x in direction).sqrt();relative=(2*mu*com).sqrt()
                product=[[mass/total*p[j]+sign*relative*direction[j]/norm for j in range(3)]for mass,sign in ((m1,1),(m2,-1))]
                momentum_scale=max(abs(x)for row in product for x in row)
                self.close(c['name'],'com',com);self.close(c['name'],'total',k+q)
                for index,mass in enumerate((m1,m2)):
                    self.close(c['name'],'k'+str(index),sum(x*x for x in product[index])/(2*mass),k+q)
                    for j in range(3):
                        self.close(c['name'],'p'+str(index)+str(j),product[index][j],momentum_scale)
                        self.close(c['name'],'v'+str(index)+str(j),product[index][j]/mass,momentum_scale/mass)

    def test_inherit_velocity_independent_signed_charge_momentum_energy(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='INHERIT':continue
                parent,heavy,pc,hc,bc,*r=map(D,c['values']);v,p=r[:3],r[3:6];bath,q=r[6:]
                speed2=sum(x*x for x in v);self.close(c['name'],'heavy_k',heavy*speed2/2)
                self.close(c['name'],'bath_e',(parent-heavy)*speed2/2+bath+q)
                self.close(c['name'],'bath_c',pc+bc-hc)
                for j in range(3):
                    self.close(c['name'],'v'+str(j),v[j]);self.close(c['name'],'p'+str(j),heavy*v[j])
                    self.close(c['name'],'bath_p'+str(j),(parent-heavy)*v[j]+p[j])

    def test_thermal_rate_and_finite_inventory_independent_decimal(self):
        with localcontext()as ctx:
            ctx.prec=110
            for c in self.cases:
                if c['operation']!='THERM':continue
                nu,e,k,t,n,dt=map(D,c['values']);rate=nu*(-e/(k*t)).exp();exposure=rate*dt
                self.close(c['name'],'rate',rate);self.close(c['name'],'hazard',exposure)
                self.close(c['name'],'remaining',n*(-exposure).exp(),n)
                self.close(c['name'],'released',n*loss_fraction(exposure))


def main():
    global OUTPUT
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();OUTPUT=args.output_dir.resolve();OUTPUT.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    with(OUTPUT/'unittest.log').open('w')as log:r=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReducedModelControls))
    errors=len(r.errors);failed=len(r.failures);skipped=len(r.skipped);passed=max(0,r.testsRun-errors-failed-skipped)
    status='ERROR'if errors else'FAIL'if failed else'PASS'if passed else'SKIPPED'
    value={'schema':'amps-shared-reduced-model-tests-v1','status':status,'counts':{'PASS':passed,'FAIL':failed,'ERROR':errors,'SKIPPED':skipped},
        'scope':'explicit-input analytical software verification only; NOT physical source/model qualification',
        'native_guard_count':len(getattr(ReducedModelControls,'guards',[])),
        'independent_metric_count':len(getattr(ReducedModelControls,'comparisons',[])),
        'elapsed_s':time.monotonic()-started,'runner':str(Path(__file__).resolve()),'log':str(OUTPUT/'unittest.log')}
    (OUTPUT/'result.json').write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]


if __name__=='__main__':raise SystemExit(main())
