#!/usr/bin/env python3
"""Bounded supporting normal-startup/semantic-PIC composition verification.

The analytical SI recipe keeps repository nominal masses and changes only the
declared charge-token conversion. It is not a physical sodium profile. Normal
startup consumes an authenticated existing ordinary build; mechanics owns an
explicit disposable checkout and a new clean normal application build.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import time
import traceback
import uuid
from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256 as _sha
from integration_common import REPO
from run_m3_particle_backend import derive_input
from run_m3_moon_force_tracers import native_sources
from test_m3_runtime_rates import change
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command

COHORT='DECLARED-SEMANTIC-PIC-COMPOSITION-NOT-LUNAR-PHYSICS'
STARTUP={
 'normal_amps_init_executes_semantic_registration_before_real_buffer',
 'normal_prebuffer_photon_field_uses_actual_furnished_spice_epoch_and_positions',
 'normal_prebuffer_callback_cannot_run_twice',
 'normal_prebuffer_callback_cannot_be_replaced_after_allocation',
 'normal_initialized_pic_buffer_accepts_explicit_source_semantic_state'}
MECHANICS={
 'variable_source_hazard_refused_without_activation','missing_explicit_source_state_resolver_refused',
 'complete_actual_core_hook_installation','actual_compiled_ions_match_explicit_charge_unit',
 'wrong_charge_unit_refused_before_reservation','staged_partial_child_failure_rolls_back_complete_packet',
 'same_plan_retry_commits_without_parent_rng_redraw',
 'staged_daughters_use_complete_child_genealogy_and_clock_counter_one',
 'staged_clocks_equal_original_prospective_child_draws',
 'staged_products_keep_represented_number_and_original_local_interval',
 'real_pic_ownership_clone_preserves_semantic_bytes',
 'recycled_source_slot_explicitly_starts_fresh_semantic_history',
 'transaction_controls_leave_no_native_slots_or_published_list',
 'real_staged_identity_chain_retains_full_genealogy_to_capacity',
 'actual_genealogy_capacity_refuses_before_clock_or_particle_publication',
 'actual_movers_continue_complete_semantic_daughter_packet',
 'actual_mpi_migration_preserves_full_daughter_genealogy',
 'semantic_mover_parent_is_deleted_once_without_slot_leak',
 'collector_rejects_mismatched_compiled_charge_unit_before_install',
 'collector_accepts_same_explicit_source_charge_as_controller_and_transaction',
 'actual_sphere_collector_preserves_owner_and_uses_source_SI_charge_ledger'}

def sha(path):return _sha(Path(path))
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def _start(output,wall_seconds,scope):
    output=Path(output).resolve()
    if output.exists() and any(output.iterdir()):raise ValueError('evidence output must be new or empty')
    output.mkdir(parents=True,exist_ok=True)
    result={'schema':'srcMoon-M3-normal-semantic-cohort-v1','status':'ERROR','scope':scope,
        'formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','observational_scoring':False,'commands':[]}
    identity={'execution':'not-run','source_root':str(REPO),'source_revision':subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),'random_seed':23,'cohort':COHORT,
        'genealogy_capacity':4,'mpi_layouts':[1,8],'orders':[1,2],'omp_threads':1,
        'physical_source_activated':False,'external_observations':[]}
    return output,result,identity,CampaignBudget(wall_seconds,min(180,wall_seconds*.1))
def _finish(output,result,identity,started):
    result['elapsed_wall_seconds']=time.monotonic()-started
    identity['artifact_sha256']={str(p.relative_to(output)):sha(p) for p in output.rglob('*')
        if p.is_file() and p.name not in ('result.json','run_manifest.json')}
    dump(output/'run_manifest.json',identity);dump(output/'result.json',result);return result
def _error(output,result,error):
    result.update(status='ERROR',reason=str(error),exception_type=type(error).__name__)
    (output/'infrastructure_exception.log').write_text(traceback.format_exc())
def _execute(budget,result,output,command,label,cwd=REPO,cores=1,native=False):
    value=budget.run(command,cwd,output/(label+'.log'),output=output,
        allocated_cores=cores,environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    result['commands'].append(value)
    if value.get('reason') or value['returncode'] not in ((0,1) if native else (0,)):
        raise RuntimeError(label+' infrastructure failure: '+str(value.get('reason') or value['returncode']))
    return value
def _checks(directory,required,pattern,execution):
    if execution.get('returncode') not in (0,1) or execution.get('reason'):
        raise ValueError('native signal/arbitrary exit is not an evaluated check outcome')
    with (directory/'semantic.checks.csv').open(newline='') as stream:
        reader=csv.DictReader(stream)
        if reader.fieldnames!=['check','status']:raise ValueError('native named-check header mismatch')
        rows=list(reader)
    summaries=re.findall(pattern,(directory/'native.log').read_text())
    if len(summaries)!=1 or len(rows)!=len(required) or {r['check'] for r in rows}!=required or any(
            r['status'] not in ('PASS','FAIL') for r in rows):
        raise ValueError('native completed named controls missing/duplicate/malformed')
    count,status=summaries[0]
    failed=any(r['status']=='FAIL' for r in rows)
    if int(count)!=len(required) or status not in ('PASS','FAIL') or (status=='FAIL')!=failed or bool(execution['returncode'])!=failed:
        raise ValueError('native completion summary/exit contradicts evaluated controls')
    return rows
def _compile(build,source,exe):
    command=_probe_compile_command(build,exe);command.insert(1,'-DMPI_ON')
    command[command.index(str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp'))]=str(source)
    return command
def _test_identity(identity,source,output):
    paths=[source,Path(__file__).resolve(),Path(__file__).parent/'semantic_cohort/CONTRACT.md']
    if source.name=='m3_native_semantic_cohort_probe.cpp':paths.append(Path(__file__).parent/'m3_native_reaction_controller_probe.cpp')
    identity['test_sources']={str(p.resolve()):sha(p) for p in paths}
    # Retain actual adapter bytes as well as hashes: a later justified test
    # correction must not make a failed attempt impossible to reconstruct.
    for path in paths:
        target=output/'test-sources'/path.resolve().relative_to(REPO)
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
def _verify_tests(identity):
    if any(sha(p)!=value for p,value in identity['test_sources'].items()):raise ValueError('test source changed during execution')
def python_uniform(parent,ordinal):
    """Independent RFC2104 HMAC and fixed52-bit midpoint conversion."""
    u64=lambda x:struct.pack('>Q',x)
    cohort=COHORT.encode();message=b'AMPS-MOON-M3-HMAC-SHA256-COUNTER-v1'+u64(len(cohort))+cohort
    message+=u64(parent)+u64(1)+u64(ordinal)+u64(3)+u64(0)
    digest=hmac.new(u64(23),message,hashlib.sha256).hexdigest()
    return (int(digest[:13],16)+.5)*2**-52
def verify_keys(directory,ranks):
    rows=[]
    for rank in range(ranks):
        with (directory/f'semantic.rank{rank}.keys.csv').open(newline='') as stream:
            reader=csv.DictReader(stream)
            if reader.fieldnames!=['context','original_parent','product_ordinal','prospective_uniform','optical_depth','product_clock_counter','owner']:
                raise ValueError('semantic key evidence header mismatch')
            current=list(reader)
        direct=[r for r in current if r['context']=='direct']
        if len(direct)!=2 or {int(r['product_ordinal']) for r in direct}!={0,1}:raise ValueError('per-rank complete staged key evidence missing')
        for r in current:
            if int(r['owner'])!=rank or int(r['original_parent'])!=73 or int(r['product_ordinal']) not in (0,1) or int(r['product_clock_counter'])!=1:
                raise ValueError('semantic key identity/counter/owner malformed')
            expected=python_uniform(73,int(r['product_ordinal']));actual=float(r['prospective_uniform']);depth=float(r['optical_depth'])
            if not math.isfinite(depth) or depth<=0:raise ValueError('semantic depth missing/nonfinite/nonpositive')
            if r['context']=='direct' and actual!=expected:raise ValueError('actual opaque child draw differs from independent HMAC oracle')
            if r['context'] not in ('direct','mover'):raise ValueError('unrecognized semantic key context')
            r.update(expected_uniform=expected)
        rows+=current
    movers=[r for r in rows if r['context']=='mover']
    if len(movers)!=2 or {int(r['product_ordinal']) for r in movers}!={0,1}:raise ValueError('global complete migrated daughter evidence missing')
    return rows
def _plots(output,rows):
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    with (output/'key_comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    selected=[r for r in rows if r['context']=='direct'];expected=[r['expected_uniform'] for r in selected]
    actual=[float(r['prospective_uniform']) for r in selected]
    figure,axes=plt.subplots(1,2,figsize=(10,4));axes[0].plot(expected,expected,label='independent HMAC identity')
    axes[0].scatter(expected,actual,s=16,label='actual staged PIC keys');axes[0].set_xlabel('Independent uniform');axes[0].set_ylabel('PIC prospective uniform')
    axes[1].plot(range(len(actual)),[a-b for a,b in zip(actual,expected)],'.');axes[1].set_xlabel('Retained native comparison');axes[1].set_ylabel('Actual minus independent uniform')
    for axis in axes:axis.grid(True)
    axes[0].legend(framealpha=1);figure.suptitle('Analytical semantic storage controls — no observations');figure.tight_layout()
    for suffix in ('png','eps'):figure.savefig(output/('key_comparison.'+suffix))
    plt.close(figure)

def run_startup(output,wall_seconds=600,build_evidence=None):
    """Authenticate and use a real normal SPICE build; never configure here."""
    started=time.monotonic();output,result,identity,budget=_start(output,wall_seconds,'actual normal pre-buffer registration and furnished SPICE; no physical provider')
    try:
        if build_evidence is None:
            result.update(status='SKIPPED',reason='authentic normal SPICE build proof unavailable');return _finish(output,result,identity,started)
        evidence=json.loads(Path(build_evidence).read_text()) if isinstance(build_evidence,(str,Path)) else build_evidence
        if not isinstance(evidence,dict) or evidence.get('status')!='PASS' or Path(evidence['source_root']).resolve()!=REPO:
            raise ValueError('normal build proof is not a successful matching checkout')
        if native_sources(evidence['production_source_hashes'])!=native_sources(production_source_identity()):raise ValueError('normal proof native source changed')
        library=Path(evidence['library']['path']);build=library.parent
        for role in ('library','executable','build_log'):
            if sha(evidence[role]['path'])!=evidence[role]['sha256']:raise ValueError('normal '+role+' changed')
        if re.findall(r'^EXIT: (-?\d+)$',Path(evidence['build_log']['path']).read_text(),re.M)!=['0','0']:
            raise ValueError('normal build proof lacks actual successful configure/make')
        for path,digest in evidence['generated_snapshots'].items():
            if sha(path)!=digest:raise ValueError('normal generated ABI changed')
        source=Path(__file__).with_name('m3_normal_semantic_registration_probe.cpp');_test_identity(identity,source,output)
        exe=output/'normal-semantic-registration';_execute(budget,result,output,_compile(build,source,exe),'compile-probe')
        counts=collections.Counter();layouts=[];external={}
        for order in ('first','second'):
            for ranks in (1,8):
                folder=output/f'{order}-np{ranks}';folder.mkdir();runtime=folder/'runtime';(runtime/'PT/plots').mkdir(parents=True)
                text=(REPO/f'srcMoon/examples/m3_tracing_{order}_order.in').read_text()
                for key,suffix in (('surface_cea_file','cea'),('surface_tecplot_file','dat')):text=change(text,key,str(folder/('sphere.'+suffix)))
                application=folder/'application.in';application.write_text(text)
                execution=_execute(budget,result,output,['mpiexec','-n',str(ranks),str(exe),str(application),str(folder/'semantic'),'1' if order=='first' else '2'],
                    f'{order}-np{ranks}/native',runtime,ranks,True)
                checks=_checks(folder,STARTUP,r'M3 normal semantic registration checks=(\d+) status=(PASS|FAIL|ERROR)',execution)
                counts.update(r['status'] for r in checks);layouts.append({'order':order,'mpi_ranks':ranks,'checks':len(checks)})
                paths=(folder/'semantic.spice_files.txt').read_text().splitlines()
                if not paths or len(paths)!=len(set(paths)):raise ValueError('furnished SPICE identities missing/duplicate')
                for p in paths:
                    p=Path(p);p=(p if p.is_absolute() else runtime/p).resolve();external[str(p)]=sha(p)
        _verify_tests(identity)
        if native_sources(evidence['production_source_hashes'])!=native_sources(production_source_identity()):raise ValueError('normal source changed during startup')
        identity.update(execution='actual amps_init prebuffer hook, both orders, MPI1/8',normal_build_evidence=evidence,
            executable={'path':str(exe),'sha256':sha(exe)},external_data_sha256=external)
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',counts={s:counts[s] for s in ('PASS','FAIL','ERROR','SKIPPED')},layouts=layouts)
    except Exception as error:_error(output,result,error)
    return _finish(output,result,identity,started)

def run(output,*,isolated_checkout=False,wall_seconds=1800,build_jobs=8,
        spice_path='/home/vtenishe/SPICE/cspice',spice_kernels='/home/vtenishe/SPICE/Kernels'):
    """Clean-build the explicitly analytical SI-charge eight-species recipe."""
    started=time.monotonic();output,result,identity,budget=_start(output,wall_seconds,'analytical SI-charge actual PIC transactions and semantic movers; NOTphysicalNa/CEX')
    created=[];application='moon_m3_semantic_'+uuid.uuid4().hex[:12]
    try:
        if not isolated_checkout:
            result.update(status='SKIPPED',reason='explicit disposable-checkout ownership required');return _finish(output,result,identity,started)
        identity['production_source_hashes']=production_source_identity()
        identity['compiler']={'command':['mpicxx','--version'],'output':subprocess.check_output(['mpicxx','--version'],text=True)}
        identity['mpi']={'command':['mpiexec','--version'],'output':subprocess.check_output(['mpiexec','--version'],text=True)}
        source=Path(__file__).with_name('m3_native_semantic_cohort_probe.cpp');_test_identity(identity,source,output)
        legacy=(REPO/'input/species.input').read_text();text=legacy
        constants=(REPO/'src/general/constants.h').read_text()
        match=re.search(r'(?m)^\s*#define\s+ElectronCharge\s+([0-9.eE+-]+)\s*$',constants)
        if not match:raise ValueError('legacy elementary-charge declaration not found')
        token=1.602176634e-19/float(match.group(1))
        # Each ion is explicit; no chemical identity is guessed from a suffix.
        for symbol in ('Na_PLUS','He_PLUS','Ne_PLUS','Ar_PLUS'):
            pattern=re.compile(r'(?im)(^\s*#component\s*=\s*'+symbol+r'\s*$)(.*?)(?=^\s*#component|^\s*#endspecies)',re.S)
            matches=list(pattern.finditer(text))
            if len(matches)!=1:raise ValueError('explicit species section missing/duplicate: '+symbol)
            m=matches[0];body,n=re.subn(r'(?im)^(\s*charge\s*=)[^\n]+',lambda q:q.group(1)+repr(token),m.group(2))
            if n!=1:raise ValueError('explicit species charge row missing/duplicate')
            text=text[:m.start(2)]+body+text[m.end(2):]
        recipe=output/'analytical_si_species.input';recipe.write_text(text)
        probe=output/'charge-expression.cpp';probe.write_text('#include "'+str(REPO/'src/general/constants.h')+'"\nint main(){return ('+repr(token)+'*ElectronCharge)==1.602176634e-19?0:1;}\n')
        _execute(budget,result,output,['g++','-std=c++17',str(probe),'-o',str(output/'charge-expression')],'compile-charge-expression')
        _execute(budget,result,output,[str(output/'charge-expression')],'verify-charge-expression')
        derived=derive_input((REPO/'input/moon.input').read_text(),application)
        derived,n=re.subn(r'(?im)^\s*#include\s+species\.input\s*$',lambda _: '#include '+str(recipe),derived)
        if n!=1:raise ValueError('canonical species include missing/duplicate')
        for suffix,content in {'.input':derived,'.input.pl':(REPO/'input/moon.input.pl').read_text(),
                '.input.spicekernels':(REPO/'input/moon.input.spicekernels').read_text()}.items():
            path=REPO/'input'/(application+suffix);path.write_text(content);created.append(path)
            snapshot=output/('compile'+suffix);snapshot.write_text(content)
        identity['compile_recipe']={'scope':'ANALYTICAL-NOMINAL-MASSES-EXACT-SI-CHARGE-NOT-PHYSICAL-PROFILE',
            'charge_unit_C':1.602176634e-19,'legacy_charge_unit_C':float(match.group(1)),
            'charge_token':token,'mass_policy':'unchanged repository nominal masses','species_input_sha256':sha(recipe),
            'source_species_input_sha256':sha(REPO/'input/species.input')}
        build=REPO/'build'
        if build.exists():shutil.rmtree(build)
        _execute(budget,result,output,[str(REPO/'Config.pl'),'-application='+application,'-amps-test=off',
            '-spice-path='+str(spice_path),'-spice-kernels='+str(spice_kernels)],'configure',cores=build_jobs)
        _execute(budget,result,output,['make','-j'+str(build_jobs),'OPT3=-O3 -g'],'build',cores=build_jobs)
        exe=output/'semantic-cohort';_execute(budget,result,output,_compile(build,source,exe),'compile-probe')
        generated={}
        for rel in ('pic/pic.h','pic/picGlobal.dfn','pic/picSpeciesMacro.dfn','pic/picParticleDataMacro.h','general/rnd.h'):
            target=output/'generated'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(build/rel,target);generated[str(target)]=sha(target)
        counts=collections.Counter();layouts=[];keyrows=[]
        for order in (1,2):
            for ranks in (1,8):
                folder=output/f'order{order}-np{ranks}';folder.mkdir();runtime=folder/'runtime';(runtime/'PT').mkdir(parents=True)
                execution=_execute(budget,result,output,['mpiexec','-n',str(ranks),str(exe),'mechanics',str(order),'si',str(folder/'semantic')],
                    f'order{order}-np{ranks}/native',runtime,ranks,True)
                checks=_checks(folder,MECHANICS,r'M3 semantic cohort mode=mechanics order='+str(order)+r' checks=(\d+) status=(PASS|FAIL|ERROR)',execution)
                counts.update(r['status'] for r in checks);layouts.append({'order':order,'mpi_ranks':ranks,'checks':len(checks)})
                keyrows+=verify_keys(folder,ranks)
        _plots(output,keyrows);_verify_tests(identity)
        if native_sources(identity['production_source_hashes'])!=native_sources(production_source_identity()):raise ValueError('native source changed during producer/run')
        identity.update(execution='fresh actual SI-charge normal application and PIC semantic composition',
            library={'path':str(build/'libAMPS.a'),'sha256':sha(build/'libAMPS.a')},
            executable={'path':str(exe),'sha256':sha(exe)},generated_snapshots=generated)
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',counts={s:counts[s] for s in ('PASS','FAIL','ERROR','SKIPPED')},layouts=layouts)
    except ImportError as error:result.update(status='SKIPPED',reason='plotting dependency unavailable: '+str(error))
    except Exception as error:_error(output,result,error)
    finally:
        for path in created:path.unlink(missing_ok=True)
        for suffix in ('.input','.input.pl','.input.spicekernels'):(REPO/(application+suffix)).unlink(missing_ok=True)
    return _finish(output,result,identity,started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--normal-startup-only',action='store_true');parser.add_argument('--build-evidence',type=Path)
    parser.add_argument('--isolated-checkout',action='store_true');parser.add_argument('--wall-seconds',type=float,default=1800)
    parser.add_argument('--build-jobs',type=int,default=8)
    parser.add_argument('--spice-path',default='/home/vtenishe/SPICE/cspice')
    parser.add_argument('--spice-kernels',default='/home/vtenishe/SPICE/Kernels');args=parser.parse_args()
    if not 0<args.wall_seconds<=3600 or args.build_jobs<1:parser.error('positive wall<=3600 and build jobs required')
    result=run_startup(args.output_dir,args.wall_seconds,args.build_evidence) if args.normal_startup_only else run(
        args.output_dir,isolated_checkout=args.isolated_checkout,wall_seconds=args.wall_seconds,build_jobs=args.build_jobs,
        spice_path=args.spice_path,spice_kernels=args.spice_kernels)
    print(json.dumps(result,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[result['status']])
