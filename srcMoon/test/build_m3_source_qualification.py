#!/usr/bin/env python3
"""Bounded source-bound Na/Na+ or Na/Na+/H build in a disposable checkout.

No particle execution or source qualification is inferred. Public invocation
is from the active checkout; its private worker uses the normal maintained
Config/remove-build/make and actual producer-evidence functions. Read the
source_qualification_build contract before using this opt-in builder.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256

ROOT=Path(__file__).resolve().parents[2]
EXIT={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}
WORKFLOW=('srcMoon/test/build_m3_source_qualification.py','srcMoon/test/build_m2_qualification.py',
 'srcMoon/test/run_m3_development.py','srcMoon/test/run_integration_tests.py',
 'srcMoon/test/campaign_budget.py','utility/chemical_database/state_profiles/compile_recipe/prepare.py',
 'utility/chemical_database/state_profiles/compile_recipe/export.cpp')
EXPECTED_MAPPING={'Na':'23Na:ground','NA_PLUS':'23Na+:ground'}
GROUND_CAPTURE_MAPPING = {**EXPECTED_MAPPING, 'H': '1H:ground'}


def strict(path):
    def pairs(items):
        result={}
        for k,v in items:
            if k in result:raise ValueError('duplicate JSON key: '+k)
            result[k]=v
        return result
    def invalid(token):raise ValueError('nonfinite JSON: '+token)
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=invalid)


def validate_mapping(path):
    value=strict(path)
    # The third tracked particle is the fast H atom formed from the incident
    # plasma proton. Only the collision fixture requests this ABI. The order
    # fixes the species indices shared by the native PIC allocator and mover.
    if not isinstance(value, dict) or list(value.items()) not in (
            list(EXPECTED_MAPPING.items()), list(GROUND_CAPTURE_MAPPING.items())):
        raise ValueError('expected ordered ground-state mapping for Na, NA_PLUS and optionally H')
    return value


def checkout_identity(checkout):
    """Invoke the identical maintained read-only source inventory at that root."""
    path=Path(checkout)/'srcMoon/test/build_m2_qualification.py'
    spec=importlib.util.spec_from_file_location('source_build_identity_'+str(id(checkout)),path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.production_source_identity()


def validate_checkout(active,checkout,expected=None):
    active,checkout=Path(active).resolve(),Path(checkout).resolve()
    if checkout==active or active in checkout.parents or checkout in active.parents:
        raise ValueError('disposable checkout must be outside and disjoint from the active source tree')
    if not checkout.is_dir()or not(checkout/'Config.pl').is_file():
        raise ValueError('explicit disposable AMPS checkout unavailable')
    for relative in WORKFLOW:
        a,b=active/relative,checkout/relative
        if not a.is_file()or not b.is_file()or sha256(a)!=sha256(b):
            raise ValueError('disposable workflow source mismatch: '+relative)
        if active in b.resolve().parents:
            raise ValueError('private workflow aliases active checkout: '+relative)
    expected=production_source_identity()if expected is None else expected
    actual=checkout_identity(checkout)
    if actual!=expected:
        differing=sorted(k for k in set(actual)|set(expected)if actual.get(k)!=expected.get(k))
        raise ValueError('disposable production bytes differ from active working source: '+','.join(differing[:12]))
    for relative in actual:
        if active in (checkout/relative).resolve().parents:
            raise ValueError('private production source aliases active checkout: '+relative)
    # Absolute unre-based site paths can direct generated output back into the
    # active checkout even though cwd is private. Refuse before Config.pl.
    for name in ('Makefile','Makefile.def','Makefile.local','Makefile.conf','.amps.conf','.general.conf'):
        p=checkout/name
        if p.is_file()and str(active)in p.read_text(errors='replace'):
            raise ValueError('private site configuration still names active checkout: '+name)
    definition=checkout/'Makefile.def'
    if definition.is_file()and 'Makefile.test'in definition.read_text()and not(checkout/'Makefile.test').is_file():
        raise FileNotFoundError('private generated Makefile.test unavailable; run ./Config.pl -install in the DISPOSABLE checkout before this build')
    return actual


def worker(args):
    """Only called by the supervised public owner inside the disposable tree."""
    from run_integration_tests import clean_configure_build
    from run_m3_development import normal_restore_evidence
    out=Path(args.output_dir).resolve();active=Path(args.active_root).resolve()
    if ROOT!=Path(args.checkout).resolve()or ROOT==active or active in ROOT.parents:
        raise ValueError('private worker root/disposable boundary mismatch')
    validation=strict(out/'preflight.json')
    before=production_source_identity()
    if before!=validation['production_source_hashes']:raise ValueError('private source changed before worker staging')
    validate_mapping(args.mapping)
    recipe_path=ROOT/'utility/chemical_database/state_profiles/compile_recipe/prepare.py'
    spec=importlib.util.spec_from_file_location('source_bound_compile_recipe',recipe_path)
    recipe=importlib.util.module_from_spec(spec);spec.loader.exec_module(recipe)
    staged=out/'source-recipe'
    recipe.prepare(ROOT,args.profile,args.mapping,ROOT/'input/moon.input',staged,
      reference_problem=args.reference_problem,packet_layout='conservative-products')
    made=[]
    for suffix in ('.input','.input.pl','.input.spicekernels'):
        source=staged/('moon_source_profile'+suffix);target=ROOT/'input'/source.name
        if target.exists():raise ValueError('private derived application already exists: '+str(target))
        made.append((source,target))
    for source,target in made:shutil.copy2(source,target)
    # This function removes build first, configures the explicit derived
    # application and executes the existing bounded16-way -O3/-g make command.
    # The OUTER CampaignBudget supervises this worker and all compiler children.
    clean_configure_build('moon_source_profile',args.spice_path,args.spice_kernels,out,'moon-spice')
    proof=normal_restore_evidence(out,before)
    actual_recipe=strict(staged/'run_manifest.json')
    if actual_recipe.get('packet_layout')!='conservative-products':
        raise ValueError('source qualification requires explicitly staged conservative product packet ABI')
    # Check the effective generated macro after the ACTUAL successful clean
    # build, rather than inferring layout from a requested input string. This
    # reads the same generated definitions retained with the library receipt.
    probe='#include "picGlobal.dfn"\n#if !defined(_INDIVIDUAL_PARTICLE_WEIGHT_MODE_) || !defined(_INDIVIDUAL_PARTICLE_WEIGHT_ON_)\n#error missing individual particle weight definitions\n#endif\n#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_\n#error conservative product packet storage is disabled\n#endif\n'
    command=['g++','-std=c++17','-E','-P','-x','c++','-I',str(ROOT/'build/pic'),'-']
    checked=subprocess.run(command,input=probe,text=True,capture_output=True,timeout=30)
    (out/'packet_layout.stdout.log').write_text(checked.stdout)
    (out/'packet_layout.stderr.log').write_text(checked.stderr)
    if checked.returncode:raise ValueError('actual generated conservative product weight storage unavailable; see packet_layout.stderr.log')
    proof['source_qualification_recipe']={'path':str(staged/'run_manifest.json'),
      'sha256':sha256(staged/'run_manifest.json'),'reference_problem':args.reference_problem,
      'packet_layout':actual_recipe['packet_layout'],
      'packet_layout_proof':{'command':command,'returncode':checked.returncode,
        'generated_header':{'path':str(ROOT/'build/pic/picGlobal.dfn'),'sha256':sha256(ROOT/'build/pic/picGlobal.dfn')},
        'stdout':{'path':str(out/'packet_layout.stdout.log'),'sha256':sha256(out/'packet_layout.stdout.log')},
        'stderr':{'path':str(out/'packet_layout.stderr.log'),'sha256':sha256(out/'packet_layout.stderr.log')}},
      'mapping':actual_recipe['mapping'],'profile':actual_recipe['profile'],
      'compiled_charge_unit_C':actual_recipe['compiled_charge_unit_C'],
      'charge_source_identity':actual_recipe['charge_source_identity'],
      'derived_inputs':{str(p):sha256(p)for p in staged.iterdir()if p.is_file()}}
    proof['particle_execution']=False;proof['physical_activation']=False
    proof['formal_M3_status']='SKIPPED / NOT QUALIFIED'
    (out/'normal_build_evidence.json').write_text(json.dumps(proof,indent=2)+'\n')
    return proof


def run(args):
    started=time.monotonic();out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'schema':'srcMoon-M3-source-qualification-build-v1','status':'ERROR',
      'scope':'source-bound disposable build only; no particle execution or physical qualification',
      'physical_activation':False,'formal_M3_status':'SKIPPED / NOT QUALIFIED'}
    budget=None
    try:
        if not 0<args.wall_seconds<=3600:raise ValueError('one total wall allowance in(0,3600] seconds required')
        budget=CampaignBudget(args.wall_seconds,min(120,args.wall_seconds/10),memory_bytes=200_000_000_000)
        checkout=Path(args.checkout).resolve();validate_mapping(args.mapping)
        identity=validate_checkout(ROOT,checkout)
        for p in (args.profile,args.mapping):
            if not Path(p).is_file():raise FileNotFoundError('explicit source profile/mapping unavailable: '+str(p))
        preflight={'schema':'srcMoon-M3-source-build-preflight-v1','active_root':str(ROOT),
          'checkout':str(checkout),'production_source_hashes':identity,
          'workflow_sha256':{p:sha256(ROOT/p)for p in WORKFLOW},
          'profile':{'path':str(Path(args.profile).resolve()),'sha256':sha256(args.profile)},
          'mapping':{'path':str(Path(args.mapping).resolve()),'sha256':sha256(args.mapping)}}
        (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
        command=[sys.executable,'-B',str(checkout/'srcMoon/test/build_m3_source_qualification.py'),
          '--_worker','--active-root',str(ROOT),'--checkout',str(checkout),'--output-dir',str(out),
          '--profile',str(Path(args.profile).resolve()),'--mapping',str(Path(args.mapping).resolve()),
          '--reference-problem',args.reference_problem,'--spice-path',args.spice_path,'--spice-kernels',args.spice_kernels]
        execution=budget.run(command,checkout,out/'worker.log',output=out/'resources',allocated_cores=16,
          environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        result['execution']=execution
        if execution.get('reason'):raise RuntimeError('aggregate resource guard: '+execution['reason'])
        if execution['returncode']==77:
            result.update(status='SKIPPED',reason='source-bound recipe state/dependency unavailable; see worker.log')
        elif execution['returncode']!=0:raise RuntimeError('private source-bound build failed; see worker.log')
        else:
            evidence=out/'normal_build_evidence.json';proof=strict(evidence)
            if proof.get('schema')!='srcMoon-actual-normal-restore-producer-v1'or proof.get('status')!='PASS':
                raise ValueError('actual normal producer evidence unavailable/malformed')
            if identity!=production_source_identity():raise ValueError('active production source changed during private build')
            if sha256(args.profile)!=preflight['profile']['sha256']or sha256(args.mapping)!=preflight['mapping']['sha256']:
                raise ValueError('source profile/mapping changed during build')
            result.update(status='PASS',producer_manifest={'path':str(evidence),'sha256':sha256(evidence)})
    except FileNotFoundError as e:result.update(status='SKIPPED',reason=str(e))
    except Exception as e:
        (out/'exception.log').write_text(traceback.format_exc());result.update(status='ERROR',reason=str(e))
    result['counts']={s:int(result['status']==s)for s in EXIT}
    result['elapsed_wall_seconds']=time.monotonic()-started
    result['runner']=str(Path(__file__).resolve());(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('profile','mapping','checkout','output-dir'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--reference-problem',choices=('homogeneous-zero-acceleration','normal-template'),required=True)
    p.add_argument('--wall-seconds',type=float,default=3600)
    p.add_argument('--spice-path',default='/home/vtenishe/SPICE/cspice')
    p.add_argument('--spice-kernels',default='/home/vtenishe/SPICE/Kernels')
    p.add_argument('--_worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--active-root',type=Path,help=argparse.SUPPRESS)
    a=p.parse_args()
    if a._worker:
        try:worker(a)
        except (LookupError,FileNotFoundError)as e:print(str(e));return 77
        except Exception:traceback.print_exc();return 2
        return 0
    r=run(a);print(json.dumps(r,indent=2));return EXIT[r['status']]


if __name__=='__main__':raise SystemExit(main())
