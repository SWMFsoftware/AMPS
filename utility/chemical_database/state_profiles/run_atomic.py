#!/usr/bin/env python3
"""Bounded real-source Na ion preparation/native-loader controls, not activation.

Uses actual retained sources and the same shared mass kernel, profile resolver
and FormationQ as production. Administrative corruption controls never serve as
observations. Existing partial profile and its fourteen unavailable IDs remain.
"""
import argparse
import copy
import csv
from decimal import Decimal as D,localcontext
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import time

import prepare_atomic as production
import prepare as base
import run as original_runner

ROOT=Path(__file__).resolve().parents[3]


def run(args):
    start=time.monotonic();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    checks=[];manifest={'schema':'srcMoon-M3-atomic-profile-run-v1','commands':[],
      'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      'scope':'original-research profile/source/production-resolver supporting checks only',
      'MPI':'not used; dependency-free source preparation and PIC-independent production kernels',
      'physical_activation':False,'managed_database_modified':False}
    def check(name,passed):checks.append({'name':name,'status':'PASS'if passed else'FAIL'})
    try:
        inputs=[Path(args.raw_dir)/n for n in base.SOURCE_HASHES]+[args.energy_file,args.label_file,args.help_file,args.decision_register]
        before={str(p):(base.digest(p),p.stat().st_mtime_ns)for p in inputs}
        profile=production.prepare_atomic(args.raw_dir,args.decision_register,args.energy_file,args.label_file,args.help_file,args.canonical_root,out/'profile',out/'work')
        partial=base.strict_json(out/'work/partial/profile.json')
        check('default partial profile remains seven available/fourteen unavailable',partial['available_state_count']==7 and partial['unavailable_state_count']==14)
        check('opt-in profile adds exactly one available state',profile['available_state_count']==8 and profile['unavailable_state_count']==14)
        with Path(profile['states']['path']).open(newline='')as f:states=list(csv.DictReader(f))
        with Path(partial['states']['path']).open(newline='')as f:old=list(csv.DictReader(f))
        old_unavailable={r['state_id']:r for r in old if r['status']=='unavailable'}
        check('all fourteen original unavailable declarations byte-field equal',old_unavailable=={r['state_id']:r for r in states if r['status']=='unavailable'})
        check('no Na ion added to default partial output',not any(r['state_id']=='23Na+:ground'for r in old))
        compiler=shutil.which('g++')
        if not compiler:raise FileNotFoundError('g++ unavailable')
        sources=[ROOT/'srcMoon/species/PhysicalStateProfile.cpp',ROOT/'srcMoon/species/ReducedKinematics.cpp',
          ROOT/'src/models/chemical_database/DataIO.cpp',Path(__file__).with_name('atomic_profile_probe.cpp')]
        identities=sources+[Path(__file__).with_name(n)for n in ('probe.cpp','prepare.py','prepare_atomic.py','atomic_mass_probe.cpp','run_atomic.py')]+[
          ROOT/'srcMoon/species/PhysicalStateProfile.h',ROOT/'srcMoon/species/ReducedKinematics.h',
          ROOT/'srcMoon/species/ReactionExecution.h',ROOT/'src/models/chemical_database/DataIO.h',
          ROOT/'src/models/chemical_database/ReducedModels.h']
        manifest['source_sha256']={str(p.resolve()):base.digest(p)for p in identities}
        for p in identities:
            target=out/'source_snapshot'/p.resolve().relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        exe=out/'actual-atomic-profile-probe'
        def execute(command,label):
            manifest['commands'].append({'argv':command,'log':str(out/(label+'.log'))})
            r=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=90)
            (out/(label+'.log')).write_text(r.stdout+r.stderr);return r
        compiled=execute([compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT),*[str(p)for p in sources],'-o',str(exe)],'compile-native-profile')
        if compiled.returncode:raise RuntimeError('native actual profile/FormationQ compile failed')
        manifest['compiler']=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]
        manifest['executable_sha256']=base.digest(exe)
        actual=execute(original_runner.arguments(exe,profile),'native-atomic-profile')
        if actual.returncode not in (0,1):raise RuntimeError('actual native profile loader refused source-bound Na profile')
        lines=[line.split('\t')for line in actual.stdout.splitlines()]
        native_checks=[x for x in lines if x[0]=='CHECK']
        for _,name,status in native_checks:check('native:'+name,status=='PASS')
        summaries=[r for r in lines if r[0]=='ATOMIC-SUMMARY']
        check('native exact count/exit receipt',len(summaries)==1 and int(summaries[0][1])==len(native_checks)and int(summaries[0][2])==sum(x[2]!='PASS'for x in native_checks)and actual.returncode==(1 if any(x[2]!='PASS'for x in native_checks)else 0))
        atom={r[1]:r[2]for r in lines if r[0]=='ATOMIC'}
        receipt=base.strict_json(out/'profile/atomic_mass_derivation.json')
        cross=base.strict_json(out/'profile/na_source_energy_crosswalk.json')
        epsilon=D(next(r[1]for r in lines if r[0]=='EPSILON'))
        with localcontext()as context:
            context.prec=90
            references=receipt['independent_Decimal']
            comparison=[]
            for key in ('reference_mass_kg','physical_mass_kg'):
                expected=D(references[key]);relative=abs(D(atom[key])-expected)/expected
                check('independent-Decimal:'+key,relative<=2*epsilon)
                comparison.append({'quantity':key,'native':atom[key],'independent':str(expected),'relative_error':str(relative),'limit':str(2*epsilon)})
            ie=D(receipt['physical_ionization_energy_J'])
            check('independent-Decimal:formation',abs(D(atom['formation_J'])-ie)/ie<=2*epsilon)
            check('independent-Decimal:production-Q',abs(D(atom['Q_J'])+ie)/ie<=2*epsilon)
            check('crosswalk signed energy difference exact',D(cross['available_translational_energy_correction_eV'])==D(cross['source_threshold_eV'])-D(cross['physical_threshold_eV']))
        (out/'independent_comparisons.json').write_text(json.dumps(comparison,indent=2)+'\n')
        check('crosswalk exact native source identity',cross['physical_threshold_source_identity']==atom['formation_source'])
        with (out/'profile/na_source_energy_crosswalk.csv').open(newline='')as f:
            reader=csv.DictReader(f);cross_rows=list(reader);header=reader.fieldnames
        check('eight-column explicit crosswalk schema',header==production.CROSSWALK_COLUMNS and len(cross_rows)==1)
        check('CSV and JSON exact matched fields',cross_rows[0]=={k:str(cross[k])for k in production.CROSSWALK_COLUMNS})
        check('absolute rate bias remains unbounded',cross['absolute_rate_bias_bound'].startswith('UNAVAILABLE'))
        check('no numerical or physical source activation',not cross['physical_activation']and not cross['source_rates_modified'])
        # Corruption is exercised on a labelled temporary COPY, never an
        # original. The expected guard acts before any profile publication.
        bad=out/'CORRUPT_ADMINISTRATIVE_COPY.html';bad.write_bytes(Path(args.energy_file).read_bytes()+b'\nCHANGED FOR INTEGRITY GUARD ONLY\n')
        refused=False
        try:production.prepare_atomic(args.raw_dir,args.decision_register,bad,args.label_file,args.help_file,args.canonical_root,out/'forbidden-profile',out/'forbidden-work')
        except ValueError:refused=True
        check('changed authoritative source refused before outputs',refused and not(out/'forbidden-profile').exists()and not(out/'forbidden-work').exists())
        missing=False
        try:production.prepare_atomic(args.raw_dir,args.decision_register,out/'MISSING_SOURCE',args.label_file,args.help_file,args.canonical_root,out/'missing-profile',out/'missing-work')
        except FileNotFoundError:missing=True
        check('missing source never manufactured',missing and not(out/'missing-profile').exists())
        check('original sources bytes/timestamps unchanged',before=={str(p):(base.digest(p),p.stat().st_mtime_ns)for p in inputs})
        check('native/software source basis unchanged',manifest['source_sha256']=={str(p.resolve()):base.digest(p)for p in identities})
        check('research operation never declares managed import',profile['atomic_extension']['operation'].startswith('original-research-read-only'))
        failures=sum(x['status']=='FAIL'for x in checks)
        result={'status':'FAIL'if failures else'PASS','passed':len(checks)-failures,'failed':failures,'errors':0,'skipped':0}
    except FileNotFoundError as e:result={'status':'SKIPPED','reason':str(e),'passed':sum(x['status']=='PASS'for x in checks),'failed':sum(x['status']=='FAIL'for x in checks),'errors':0,'skipped':1}
    except (OSError,ValueError,KeyError,StopIteration,RuntimeError,ArithmeticError,subprocess.TimeoutExpired)as e:
        result={'status':'ERROR','reason':str(e),'passed':sum(x['status']=='PASS'for x in checks),'failed':sum(x['status']=='FAIL'for x in checks),'errors':1,'skipped':0}
    result.update(schema='srcMoon-M3-atomic-state-profile-support-v1',checks=checks,
      counts={k:result[v]for k,v in [('PASS','passed'),('FAIL','failed'),('ERROR','errors'),('SKIPPED','skipped')]},
      formal_M3_status='SKIPPED / NOT QUALIFIED',physical_activation=False,
      elapsed_s=time.monotonic()-start,peak_child_rss_kb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
      runner=str(Path(__file__).resolve()),output_dir=str(out))
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    source=ROOT/'test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/raw'
    p.add_argument('--raw-dir',type=Path,default=ROOT/'test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/raw')
    p.add_argument('--decision-register',type=Path,default=ROOT/'srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json')
    p.add_argument('--energy-file',type=Path,default=source/'nist_na_ionization_uncertainty.html')
    p.add_argument('--label-file',type=Path,default=source/'nist_na_handbook.html')
    p.add_argument('--help-file',type=Path,default=source/'nist_ie_help.html')
    p.add_argument('--canonical-root',type=Path,default=Path('/data/vtenishe/moon_validation_data/photochemistry'))
    a=p.parse_args()
    try:r=run(a)
    except(OSError,ValueError)as e:print(json.dumps({'status':'ERROR','reason':str(e)}));return 2
    print(json.dumps(r,indent=2));return{'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']]


if __name__=='__main__':raise SystemExit(main())
