#!/usr/bin/env python3
"""P09 actual production kernel/source controls; no lunar retention validation."""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal,localcontext
import hashlib
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time
import unittest

from helium_silica.prepare_source import prepare,PDF_SHA256
from chemical_database_build import shared_sources, shared_dependencies

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent/'helium_silica'
DEFAULT_PDF=ROOT/'test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z/tellez2012_silica_neutralization.pdf'
RELATIVE_BOUND=64*sys.float_info.epsilon # predeclared arithmetic bound, not a physical error budget
FILES=['srcMoon/species/HeliumSilicaChargeSurvival.h','srcMoon/species/HeliumSilicaChargeSurvival.cpp',
 'srcMoon/species/Photochemistry.h','srcMoon/species/Photochemistry.cpp',
 'srcMoon/test/test_m3_helium_silica.py','srcMoon/test/helium_silica/probe.cpp',
 'srcMoon/test/helium_silica/prepare_source.py','srcMoon/test/helium_silica/CONTRACT.md',
 'srcMoon/test/helium_silica/README.md']
FILES += [str(p.relative_to(ROOT)) for p in shared_dependencies(ROOT,'photochemistry')]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def adjacent_positive_binary64(value,toward_larger):
    """Return the adjacent representable positive energy, without an epsilon.

    Python3.8 does not provide math.nextafter. The source boundary controls
    need the exact neighbouring IEEE754 binary64 value, rather than an
    arbitrary physical tolerance. Positive finite binary64 encodings are
    ordered by their unsigned integer representation.
    """
    if not math.isfinite(value) or value<=0:
        raise ValueError('boundary control requires a positive finite value')
    bits=struct.unpack('>Q',struct.pack('>d',value))[0]
    return struct.unpack('>d',struct.pack('>Q',bits+(1 if toward_larger else -1)))[0]

class Context:
    def __init__(self,output,pdf,executable,table,deadline):
        self.output=output;self.pdf=pdf;self.executable=executable;self.table=table;self.deadline=deadline
        self.commands=[];self.comparisons=[];self.maximum_relative_error=0
        with table.open(newline='') as stream:self.reference={r['fit_id']:r for r in csv.DictReader(stream)}
    def query(self,mode,*arguments,pdf=None,table=None):
        command=[str(self.executable),str(pdf or self.pdf),str(table or self.table),mode]+[str(a) for a in arguments]
        remaining=self.deadline-time.monotonic()
        if remaining<=0:raise RuntimeError('whole component deadline exhausted; no further query launched')
        child=subprocess.run(command,capture_output=True,text=True,timeout=min(10,remaining))
        name=f'query-{len(self.commands):04d}';directory=self.output/'queries';directory.mkdir(exist_ok=True)
        (directory/(name+'.stdout.log')).write_text(child.stdout);(directory/(name+'.stderr.log')).write_text(child.stderr)
        self.commands.append({'command':command,'returncode':child.returncode,'stdout':str(directory/(name+'.stdout.log')),'stderr':str(directory/(name+'.stderr.log'))})
        value=json.loads(child.stdout)
        expected={'Evaluated':0,'Unavailable':77,'FAIL':1,'ERROR':2,
                  'InvalidInput':2,'NumericalError':2}.get(value.get('status'))
        if expected is None:
            raise RuntimeError('unrecognized native result status: '+name)
        if child.returncode!=expected:raise RuntimeError('native status/exit contradiction: '+name)
        return value

class KernelTests(unittest.TestCase):
    context=None
    def query(self,*args,**kwargs):return self.context.query(*args,**kwargs)
    def require_refusal(self,target,energy,x,status):
        value=self.query('evaluate',target,energy,x)
        self.assertEqual(value['status'],status);self.assertTrue(value['unchanged'])
    def test_complete_published_table_metadata(self):
        actual=self.query('catalog')['fits'];self.assertEqual(len(actual),7)
        self.assertEqual({r['fit_id'] for r in actual},set(self.context.reference))
        for row in actual:
            with self.subTest(row=row['fit_id']):
                ref=self.context.reference[row['fit_id']]
                for key,column in [('vc_m_s','characteristic_velocity_m_s'),('dvc_m_s','vc_reported_uncertainty_m_s'),('a','intercept'),('da','intercept_reported_uncertainty'),('energy_min_eV','energy_min_eV'),('energy_max_eV','energy_max_eV')]:self.assertEqual(row[key],float(ref[column]))
                self.assertEqual(row['minimum_included'],ref['energy_min_inclusive']=='true')
                self.assertEqual(row['maximum_included'],ref['energy_max_inclusive']=='true')
                for key,column in [('r_squared','r_squared'),('surface_ratio','O_Si_ratio'),('surface_ratio_uncertainty','O_Si_ratio_reported_uncertainty')]:self.assertEqual(row[key],float(ref[column]) if ref[column] else None)
    def test_published_law_and_reported_parameter_envelope(self):
        # These x queries are explicit analytical controls, not an inferred
        # Ei-to-speed conversion or synthetic experimental measurements.
        for identity,ref in self.context.reference.items():
            for x in (2e-6,5e-6,1e-5,1.5e-5,2e-5):
                with self.subTest(fit=identity,x=x):
                    energy=700 if ref['mechanism']=='AN' else 3000
                    actual=self.query('evaluate',ref['target'],energy,format(x,'.17g'))
                    self.assertEqual(actual['status'],'Evaluated');self.assertEqual(actual['fit']['fit_id'],identity)
                    vc=Decimal(ref['characteristic_velocity_m_s']);unc=Decimal(ref['vc_reported_uncertainty_m_s'])
                    with localcontext() as context:
                        context.prec=70;coordinate=Decimal.from_float(x)
                        expected=float((-vc*coordinate).exp());lower=float((-(vc+unc)*coordinate).exp());upper=float((-(vc-unc)*coordinate).exp())
                        logarithm=float(-vc*coordinate)
                    for key,value in [('probability',expected),('probability_at_vc_plus_uncertainty',lower),('probability_at_vc_minus_uncertainty',upper),('log_probability',logarithm)]:
                        relative=abs(actual[key]-value)/abs(value);self.context.maximum_relative_error=max(self.context.maximum_relative_error,relative)
                        self.assertLessEqual(relative,RELATIVE_BOUND)
                    self.assertLessEqual(actual['probability_at_vc_plus_uncertainty'],actual['probability'])
                    self.assertLessEqual(actual['probability'],actual['probability_at_vc_minus_uncertainty'])
                    self.context.comparisons.append({'fit_id':identity,'target':ref['target'],'mechanism':ref['mechanism'],'incident_energy_eV':energy,'source_coordinate_s_per_m':x,'production_probability':actual['probability'],'independent_probability':expected,'probability_at_reported_vc_plus':actual['probability_at_vc_plus_uncertainty'],'probability_at_reported_vc_minus':actual['probability_at_vc_minus_uncertainty'],'relative_error':abs(actual['probability']-expected)/expected})
    def test_oxygen_gaps_and_exact_endpoints_remain_unavailable(self):
        for target in ('O16','O18','O_mean'):
            for energy in (800,1000,1800,2000):
                with self.subTest(target=target,energy=energy):self.require_refusal(target,energy,1e-5,'Unavailable')
    def test_source_adjacent_energy_boundaries(self):
        for target in ('O16','O18','O_mean'):
            for energy,mechanism in [(600,'AN'),(adjacent_positive_binary64(800,False),'AN'),(adjacent_positive_binary64(2000,True),'CIN'),(7000,'CIN')]:
                with self.subTest(target=target,energy=energy):
                    value=self.query('evaluate',target,format(energy,'.17g'),'1e-5');self.assertEqual(value['status'],'Evaluated');self.assertEqual(value['fit']['mechanism'],mechanism)
            for energy in (adjacent_positive_binary64(600,False),adjacent_positive_binary64(7000,True)):self.require_refusal(target,format(energy,'.17g'),1e-5,'Unavailable')
    def test_silicon_complete_study_range_only(self):
        for energy in (600,800,1000,2000,7000):self.assertEqual(self.query('evaluate','Si',energy,1e-5)['status'],'Evaluated')
        for energy in (599,7001):self.require_refusal('Si',energy,1e-5,'Unavailable')
    def test_no_unpublished_silicon_AN_fit(self):
        value=self.query('missing-fit');self.assertEqual(value['status'],'Unavailable');self.assertTrue(value['unchanged'])
    def test_unloaded_registry_does_not_publish_probability(self):
        value=self.query('unloaded');self.assertEqual(value['status'],'Unavailable');self.assertTrue(value['unchanged'])
    def test_invalid_target_energy_or_coordinate_refused(self):
        self.require_refusal('He',3000,1e-5,'InvalidInput')
        for energy in ('0','-1','nan','inf'):self.require_refusal('Si',energy,1e-5,'InvalidInput')
        for x in ('0','-1','nan','inf'):self.require_refusal('Si',3000,x,'InvalidInput')
    def test_probability_underflow_is_not_fabricated_zero(self):self.require_refusal('Si',3000,'1e308','NumericalError')
    def test_velocity_projection_remains_unavailable(self):
        for target in ('Si','O16','O18','O_mean'):
            value=self.query('velocities',target,3000,1e-5)
            self.assertEqual(value['status'],'Unavailable');self.assertTrue(value['unchanged']);self.assertFalse(value['physical_velocity_mapping_available'])
    def test_changed_pdf_hash_is_error(self):
        value=self.query('catalog',pdf=self.context.output/'changed.pdf');self.assertEqual(value['status'],'ERROR')
    def test_changed_processed_source_table_is_error(self):
        value=self.query('catalog',table=self.context.output/'changed.csv');self.assertEqual(value['status'],'ERROR')
    def test_caller_hash_cannot_self_approve_another_version(self):self.assertEqual(self.query('caller-hash')['status'],'ERROR')
    def test_bad_reload_preserves_previous_registry(self):
        for pdf,table in [(self.context.output/'changed.pdf',self.context.table),(self.context.pdf,self.context.output/'changed.csv')]:
            value=self.query('reload',pdf,table);self.assertTrue(value['previous_registry_preserved'])
    def test_query_is_repeatable_and_uses_no_RNG(self):
        self.assertEqual(self.query('evaluate','O18',3000,1e-5),self.query('evaluate','O18',3000,1e-5))

def run(output,source_pdf=DEFAULT_PDF,compiler='g++',wall_seconds=180):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    result={'schema':'srcMoon-helium-silica-local-controls-v1','status':'ERROR',
      'scope':'P09 published conditional laboratory coordinate law and source/transcription controls only',
      'formal_lunar_ion_surface_status':'SKIPPED / NOT QUALIFIED','physical_velocity_mapping_status':'UNAVAILABLE',
      'observational_agreement':'NOT EVALUATED','counts':{}}
    manifest={'execution':'not-run','MPI':'none','random_seed':None,'stochastic_controls':'none',
      'canonical_data_modified':False,'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      'test_contract':str(HERE/'CONTRACT.md'),'test_contract_sha256':sha(HERE/'CONTRACT.md')}
    try:
        if not math.isfinite(wall_seconds) or wall_seconds<=0:raise ValueError('positive finite whole-component wall budget required')
        source_pdf=Path(source_pdf).resolve()
        missing=[tool for tool in (compiler,'pdftotext') if shutil.which(tool) is None]
        if not source_pdf.is_file() or missing:
            result.update(status='SKIPPED',reason='source publication or required native dependency unavailable',missing=missing,source_pdf=str(source_pdf),counts={'PASS':0,'FAIL':0,'ERROR':0,'SKIPPED':1});return result
        try:
            import matplotlib
            matplotlib.use('Agg');import matplotlib.pyplot as plt
        except ImportError:
            result.update(status='SKIPPED',reason='Matplotlib unavailable for required numerical comparison figures',counts={'PASS':0,'FAIL':0,'ERROR':0,'SKIPPED':1});return result
        manifest['source_hashes']={name:sha(ROOT/name) for name in FILES}
        evidence=prepare(source_pdf,output/'source_evidence');table=Path(evidence['output_table']['path'])
        manifest['source_evidence']=evidence
        executable=output/'helium-silica-probe';resolved=shutil.which(compiler)
        manifest['compiler']=subprocess.run([resolved,'--version'],capture_output=True,text=True,timeout=10).stdout
        command=[resolved,'-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT),str(ROOT/'srcMoon/species/HeliumSilicaChargeSurvival.cpp'),str(ROOT/'srcMoon/species/Photochemistry.cpp'),str(HERE/'probe.cpp'),*(str(p) for p in shared_sources(ROOT,'photochemistry')),'-o',str(executable)]
        manifest['compile_command']=command
        with (output/'compile.log').open('w') as log:
            child=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=min(60,wall_seconds))
        if child.returncode:raise RuntimeError('production kernel compilation failed; see compile.log')
        manifest['executable']={'path':str(executable),'sha256':sha(executable)};manifest['execution']='actual compiled production registry/conditional law'
        context=Context(output,source_pdf,executable,table,started+wall_seconds)
        (output/'changed.pdf').write_bytes(source_pdf.read_bytes()+b'\nSOURCE-CORRUPTION-CONTROL\n')
        (output/'changed.csv').write_bytes(table.read_bytes().replace(b'291000',b'291001',1))
        KernelTests.context=context
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(KernelTests)
        outcomes={case.id():'PASS' for case in suite}
        with (output/'checks.log').open('w') as log:
            tests=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
        # Several subtests can fail in one test method. Count each declared
        # method once, with ERROR taking priority over FAIL/SKIPPED, while
        # retaining every actual subtest failure below and in checks.log.
        for status,records in [('SKIPPED',tests.skipped),('FAIL',tests.failures),('ERROR',tests.errors)]:
            for case,_detail in records:
                outcomes[getattr(case,'test_case',case).id()]=status
        manifest['commands']=context.commands
        result.update(status='ERROR' if tests.errors else 'FAIL' if tests.failures else 'SKIPPED' if tests.skipped else 'PASS',
          counts={status:sum(value==status for value in outcomes.values()) for status in ('PASS','FAIL','ERROR','SKIPPED')},
          test_outcomes=outcomes,
          errors=[{'test':str(test),'detail':detail} for test,detail in tests.errors],failures=[{'test':str(test),'detail':detail} for test,detail in tests.failures],
          science_metrics={'compared_law_queries':len(context.comparisons),'maximum_relative_numerical_error':context.maximum_relative_error,'relative_arithmetic_bound':RELATIVE_BOUND,'published_fits':len(context.reference)})
        if not context.comparisons:raise RuntimeError('no actual law comparison artifact was produced')
        with (output/'comparison.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(context.comparisons[0]));writer.writeheader();writer.writerows(context.comparisons)
        fig,axes=plt.subplots(2,1,figsize=(8,8),sharex=True)
        for identity in sorted(context.reference):
            rows=[r for r in context.comparisons if r['fit_id']==identity]
            if not rows:continue
            axes[0].plot([r['source_coordinate_s_per_m']*1e6 for r in rows],[r['production_probability'] for r in rows],marker='o',label=identity)
            axes[0].plot([r['source_coordinate_s_per_m']*1e6 for r in rows],[r['independent_probability'] for r in rows],linestyle=':',color='black')
            axes[1].plot([r['source_coordinate_s_per_m']*1e6 for r in rows],[r['relative_error'] for r in rows],marker='o',label=identity)
        axes[0].set(ylabel='Conditional charged fraction',title='Published coordinate law: C++ versus Decimal\nNumerical controls only; no measured-data or lunar-retention comparison')
        axes[0].legend(fontsize=8,ncol=2,framealpha=1);axes[1].set(xlabel='Explicit source inverse-velocity coordinate [10^-6 s/m]',ylabel='Relative numerical error')
        for suffix in ('png','eps'):fig.savefig(output/('comparison.'+suffix),bbox_inches='tight')
        plt.close(fig)
        if {name:sha(ROOT/name) for name in FILES}!=manifest['source_hashes']:raise RuntimeError('production/test source changed during component run')
        if time.monotonic()-started>wall_seconds:raise RuntimeError('whole component wall budget exceeded')
        result['artifacts']={name:str(output/name) for name in ['compile.log','checks.log','comparison.csv','comparison.png','comparison.eps','source_evidence/qa.json']}
    except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as error:
        result.update(status='ERROR',reason=str(error))
        # A failed compile/export/source-stability check is an additional
        # infrastructure error. Preserve any already-evaluated kernel counts
        # rather than hiding the error behind their successful test results.
        if not result.get('counts'):
            result['counts']={'PASS':0,'FAIL':0,'ERROR':0,'SKIPPED':0}
        result['counts']['ERROR']+=1
    finally:
        result['elapsed_seconds']=time.monotonic()-started
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
        (output/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--source-pdf',type=Path,default=DEFAULT_PDF);parser.add_argument('--cxx',default='g++');parser.add_argument('--wall-seconds',type=float,default=180)
    args=parser.parse_args();value=run(args.output_dir,args.source_pdf,args.cxx,args.wall_seconds);print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]
if __name__=='__main__':raise SystemExit(main())
