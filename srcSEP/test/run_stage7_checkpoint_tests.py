#!/usr/bin/env python3
"""Native Stage-7 fixture repair and combined threaded production restart gate.

Prepare explicit immutable jobs, execute real MPI, and collect actual receipts.
Separate processes supply the continuation reference; no physical reference is
generated or adjusted. Missing roles, stale builds, startup failures and absent
native reports are ERROR. Evaluated byte/budget mismatches are FAIL. The short
64-step thread campaign supplements the separate complete-output campaign.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import signal
import subprocess
import time
import xml.etree.ElementTree as ET

from run_native_standalone_integration import pic_dependency_closure
from verify_standalone_application import rows, balance, require, GateFailure
from run_fresh_process_checkpoint_tests import quoted

REPO=Path(__file__).resolve().parents[2]
MOVERS=('parker','fte-dmumu','fte-mfp')
NATIVE=('CP05C01','CP05C02','CP05C03','CP05C04','CP05C05')
FRESH=('CP05F01','CP05F02','CP05F03')
OPERATORS=('CP05O01','CP05O02','CP05O03','CP05O04','CP05O05')
OPFRESH=('CP05OF01','CP05OF02','CP05OF03')
THREADED=('TCR01','TCR02','TCR03')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def demand(condition,message):
    """Infrastructure requirements are never converted to numerical PASS."""
    if not condition:raise ValueError(message)


def request(directory,role,prefix,operator=False):
    # Test input is the established quoted C++ protocol, not production argv.
    # All physical assets are explicit absolute paths, checked and hashed.
    text='SEP_CHECKPOINT_FRESH_PROCESS '+('1' if operator else '2')+'\nrole '+quoted(role)+'\nprefix '+quoted(str(prefix))+'\n'
    if not operator:
        text+='startup_input '+quoted(str(REPO/'srcSEP/examples/reduced-shock/field_line.in'))+'\nevent '+quoted(str(REPO/'srcSEP/examples/reduced-shock/positive_1au.event'))+'\n'
    path=directory/'fixture.args';path.write_text(text);return path


def build_identity(build):
    """Audit the real PIC consumer closure and the changed application units.

    PIC template definitions may come from any archive member. Checking just
    a selected mover misses a stale weak pool method; reuse the complete MMD
    audit and additionally compare the exact newly modified generated sources.
    """
    files=pic_dependency_closure(build)
    for name in ('main.cpp','physics_checkpoint.cpp','physics_checkpoint.h',
                 'util/sep_checkpoint_native_validation.cpp','util/sep_checkpoint_native_validation.h'):
        primary=REPO/'srcSEP'/name;generated=build/'build/main'/name
        demand(primary.read_bytes()==generated.read_bytes(),'generated checkpoint source differs: '+str(generated))
        files += [primary,generated]
    files += [build/'amps',build/'build/libAMPS.a',build/'build/main/main.a',build/'Makefile.conf',
              build/'Makefile.def',build/'Makefile.local',build/'build/pic/picGlobal.dfn']
    return {str(path):digest(path) for path in files}


def prepare(args,output):
    hashes={str(Path(__file__).resolve()):digest(__file__)}
    # Freeze imported acceptance/quoting/dependency helpers too. Replacing an
    # imported budget tolerance after execution cannot certify an old campaign.
    for name in ('run_native_standalone_integration.py','run_fresh_process_checkpoint_tests.py',
                 'verify_standalone_application.py'):
        path=Path(__file__).resolve().parent/name;hashes[str(path)]=digest(path)
    assets=[REPO/'srcSEP/examples/reduced-shock'/name for name in ('field_line.in','positive_1au.event')]
    for path in assets:hashes[str(path)]=digest(path)
    build=args.build_root.resolve();hashes.update(build_identity(build))
    operator=output/'checkpoint-operators'
    env=dict(os.environ,AMPS_INTEGRATION_REPO=str(REPO),AMPS_INTEGRATION_OUTPUT=str(operator))
    command=['make','-f',str(REPO/'srcSEP/test/native_standalone_integration.mk'),'native-checkpoint-operators']
    with (output/'operator-build.log').open('w') as log:
        result=subprocess.run(command,cwd=build,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
    demand(result.returncode==0,'operator fixture compile/link failed; see operator-build.log')
    for path in [operator,REPO/'srcSEP/test/native_checkpoint_operators.cpp',REPO/'srcSEP/test/native_standalone_integration.mk']:
        hashes[str(path)]=digest(path)
    jobs=[]
    def add(identifier,role,executable,extra,ranks=2,thread_count=0,mode='',depends=()):
        directory=output/(identifier+'-'+role+('-threads'+str(thread_count)+'-ranks'+str(ranks) if thread_count else ''))
        directory.mkdir();(directory/'PT/plots').mkdir(parents=True)
        # This immutable geometry-measure cache is optional and affects only
        # startup cost. Its bytes are hashed if reused; it is never a solution
        # reference and cannot change a physical source or acceptance bound.
        cache=args.mesh_cache
        if cache:
            import shutil
            shutil.copyfile(cache,directory/cache.name);hashes[str(cache.resolve())]=digest(cache)
        job=dict(id=identifier,role=role,ranks=ranks,threads=thread_count,mover=mode,cwd=str(directory),
                 argv=['mpiexec','-n',str(ranks),'--bind-to','none',str(executable)]+extra,
                 timeout_s=300,depends=list(depends))
        jobs.append(job);return job,directory
    for identity in NATIVE:
        job,directory=add(identity,'routine',build/'amps',[])
        inp=request(directory,'routine',directory/'physics')
        job['argv'] += ['--test',identity,'--test-input',str(inp),'--test-output-dir',str(directory/'artifacts'),
                        '--keep-test-artifacts','--test-json',str(directory/'results.json'),'--test-junit',str(directory/'junit.xml')]
        hashes[str(inp)]=digest(inp)
    for operator_only,identities,executable in [(False,FRESH,build/'amps'),(True,OPFRESH,operator)]:
        for identity in identities:
            prefix=output/(identity+'-physics');writer=None
            for role in ('write','resume'):
                job,directory=add(identity,role,executable,[],depends=[] if writer is None else [writer['cwd']])
                inp=request(directory,role,prefix,operator_only)
                job['argv'] += ['--test',identity,'--test-input',str(inp),'--test-output-dir',str(directory/'artifacts'),
                                '--keep-test-artifacts','--test-json',str(directory/'results.json'),'--test-junit',str(directory/'junit.xml')]
                job['prefix']=str(prefix);hashes[str(inp)]=digest(inp)
                if role=='write':writer=job
    for identity in OPERATORS:
        job,directory=add(identity,'routine',operator,[])
        job['argv'] += ['--test',identity,'--test-json',str(directory/'results.json'),
                        '--test-junit',str(directory/'junit.xml')]
    job,directory=add('MCP13','missing-fixture',build/'amps',[])
    job['expected_exit_code']=2
    job['argv'] += ['--test','CP05C01','--test-json',str(directory/'results.json'),
                    '--test-junit',str(directory/'junit.xml')]
    if not args.fixtures_only:
        demand(args.thread_build2 and args.thread_build4,'supply both actual two- and four-thread builds')
        for count,path in [(2,args.thread_build2),(4,args.thread_build4)]:
            path=path.resolve();hashes.update(build_identity(path))
            flags=(path/'build/pic/picGlobal.dfn').read_text()
            demand(re.search(r'^#define\s+_PIC_FL_MOVER_MULTITHREAD_\s+_PIC_MODE_ON_',flags,re.M) and
                   re.search(r'^#define\s+_PIC_NUMBER_STD_THREADS_\s+'+str(count)+r'\b',flags,re.M),
                   'thread build flags disagree with requested worker count')
            for index,mode in enumerate(MOVERS,1):
                # All three movers on both worker counts; an additional
                # 8-rank Dmumu pair exercises MPI and threads together.
                for ranks in ([1,8] if count==4 and mode=='fte-dmumu' else [1]):
                    identity=THREADED[index-1];physical=['--input',str(assets[0]),'--reduced-shock-event',str(assets[1]),
                        '--reduced-shock-particles','on','--shock-release-distance-m','1e8',
                        '--shock-return-policy','delete','--shock-source-normalization','net-upstream-escape',
                        '--total-iterations','64','--no-test-manager','--turbulence-source','prescribed',
                        '--turbulence-coupling-policy','disabled','--coupling','off','--cascade','off','--reflection','off',
                        '--turbulence-advection','off','--shock-injection','off','--coefficient-source','prescribed',
                        '--pitch-angle-diffusion-provider','constant','--constant-dmumu','0.001',
                        '--spatial-diffusion-provider','from-dmumu','--mean-free-path-provider','from-spatial',
                        '--particle-mover',mode]
                    writer=None
                    for role in ('continuous','write','resume'):
                        job,directory=add(identity,role,path/'amps',list(physical),ranks,count,mode,
                            depends=[] if role!='resume' else [writer['cwd']])
                        job['argv'] += ['--checkpoint-save',str(directory/'final.checkpoint'),'--checkpoint-at-step','45' if role=='write' else '64']
                        if role=='resume':job['argv'] += ['--resume-checkpoint',str(Path(writer['cwd'])/'final.checkpoint')]
                        if role=='write':writer=job
    plan=dict(schema='srcsep-stage7-checkpoint-plan-v1',hashes=hashes,jobs=jobs,fixtures_only=args.fixtures_only,
              acceptance=dict(physical_dt_s=300,budget_relative_tolerance=1e-10,continuation='byte-exact full checkpoint',
                              source='approved prescribed reduced-core; absorbing return; fixed physical release distance'))
    (output/'prepared.json').write_text(json.dumps(plan,indent=2)+'\n');return plan


def audit_budget(text,epochs):
    """Independent fsum SI balances, including exact zero at zero scale."""
    ledgers=rows(text,'REDUCED_SOURCE_LEDGER')
    require([v['epoch_s'] for v in ledgers]==list(epochs),'missing/wrong physical endpoint sequence')
    maximum=0.0
    for v in ledgers:
        n,e=v['born'],v['born_kinetic_j']
        maximum=max(maximum,balance([n,-v['retained'],-v['returned'],-v['domain_exit'],-v['geometry_retired']],n,'number'),
            balance([e,-v['retained_birth_j'],-v['returned_birth_j'],-v['domain_birth_j'],-v['retired_birth_j']],e,'birth energy'),
            balance([e,v['background_work_j'],v['wave_work_j'],-v['retained_kinetic_j'],-v['returned_event_j'],
                -v['domain_event_j'],-v['retired_event_j']],abs(e)+abs(v['background_work_j'])+abs(v['wave_work_j']),'signed energy'))
    require(ledgers and ledgers[-1]['born']>0,'empty source cannot qualify continuation')
    if ledgers[-1]['epoch_s']==19200:
        require(ledgers[-1]['returned']>0,'complete short driver did not exercise absorbing return')
    return maximum


def expected_job_keys(fixtures_only=False):
    """Freeze the complete role/decomposition contract independently of jobs.

    A writer PASS alone cannot establish fresh continuation. Likewise, worker
    smoke without resume or one omitted mover cannot establish threaded restart.
    Compact IDs retain their established registry identity across this repair.
    """
    keys={(identity,'routine',0,2) for identity in NATIVE+OPERATORS}
    keys.update((identity,role,0,2) for identity in FRESH+OPFRESH for role in ('write','resume'))
    keys.add(('MCP13','missing-fixture',0,2))
    if not fixtures_only:
        for identity in THREADED:
            for count in (2,4):
                for ranks in ((1,8) if identity=='TCR02' and count==4 else (1,)):
                    keys.update((identity,role,count,ranks) for role in ('continuous','write','resume'))
    return keys


def collect(plan,output):
    for path,expected in plan['hashes'].items():demand(digest(path)==expected,'source/build identity changed: '+path)
    keys=[(j.get('id'),j.get('role'),j.get('threads'),j.get('ranks')) for j in plan['jobs']]
    expected=expected_job_keys(plan.get('fixtures_only',False))
    demand(len(keys)==len(expected) and set(keys)==expected,'incomplete or duplicate stage7 role/decomposition catalog')
    results=[]
    for job in plan['jobs']:
        directory=Path(job['cwd']);result=dict(job,status='ERROR',metrics={})
        try:
            receipt=json.loads((directory/'process-result.json').read_text())
            demand(receipt.get('argv')==job['argv'] and type(receipt.get('exit_code')) is int and
                   not receipt.get('launch_error'),'missing/failed actual MPI receipt')
            if not job['threads']:
                payload=json.loads((directory/'results.json').read_text());native=payload.get('results',[])
                if job['id']=='MCP13':
                    demand(receipt['exit_code']==2,'prerequisite rejection has the wrong process status')
                    demand(len(native)==1 and native[0].get('id')=='CP05C01','missing explicit prerequisite report')
                    require(native[0].get('status')=='ERROR' and 'requires --test-input protocol2' in native[0].get('message',''),
                            'unsupported fixture did not preserve an explicit infrastructure ERROR')
                    result['status']='PASS';result['metrics']['expected_native_error']=True;results.append(result);continue
                expected=job.get('ids',[job['id']])
                demand([r.get('id') for r in native]==expected,'missing/foreign native case report')
                demand(all(r.get('status') in ('PASS','FAIL') for r in native),'native checkpoint reported unavailable/error prerequisites')
                demand(receipt['exit_code']==int(any(r.get('status')=='FAIL' for r in native)),
                       'native process status disagrees with its assertions')
                require(all(r.get('status')=='PASS' for r in native),'native checkpoint assertion failed')
                for item in native:
                    for metric in item.get('metrics',[]):
                        if metric['name'] in ('retained_source_particles_at_save','native_scattering_events','changed_native_wave_cells'):
                            require(metric['value']>0,'required particle/scattering/wave exercise was empty')
                if job['role']=='resume':
                    a=Path(job['prefix']+'.expected');b=Path(job['prefix']+'.resumed')
                    demand(a.is_file() and b.is_file() and a.stat().st_size>0 and b.stat().st_size>0,'missing fresh-process final checkpoint artifact')
                    require(a.read_bytes()==b.read_bytes(),'fresh-process complete native continuation bytes differ')
                    result['metrics']['exact_final_bytes']=a.stat().st_size
            else:
                demand(receipt['exit_code']==0,'production MPI process failed or timed out')
                text=(directory/'application.rank0.log').read_text()
                steps=range(1,46) if job['role']=='write' else range(46,65) if job['role']=='resume' else range(1,65)
                result['metrics']['maximum_relative_budget_error']=audit_budget(text,[300*s for s in steps])
                if job['role']=='resume':
                    reference=next(j for j in plan['jobs'] if j['id']==job['id'] and j['role']=='continuous' and
                                   j['threads']==job['threads'] and j['ranks']==job['ranks'])
                    a=Path(reference['cwd'])/'final.checkpoint';b=directory/'final.checkpoint'
                    demand(a.is_file() and b.is_file() and a.stat().st_size>0 and b.stat().st_size>0,'missing threaded final checkpoint artifact')
                    require(a.read_bytes()==b.read_bytes(),'threaded full-driver checkpoint continuation differs')
                    result['metrics']['exact_final_bytes']=a.stat().st_size
            result['status']='PASS'
        except (OSError,ValueError,KeyError,TypeError) as error:result['message']=str(error)
        except (GateFailure,AssertionError) as error:result['status']='FAIL';result['message']=str(error)
        results.append(result)
    report=dict(schema='srcsep-stage7-checkpoint-v1',results=results,hashes=plan['hashes'],acceptance=plan['acceptance'],
                fixtures_only=plan.get('fixtures_only',False))
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    suite=ET.Element('testsuite',name='stage7-checkpoint',tests=str(len(results)))
    for row in results:
        case=ET.SubElement(suite,'testcase',classname=row['id'],name=Path(row['cwd']).name)
        if row['status']!='PASS':ET.SubElement(case,'failure' if row['status']=='FAIL' else 'error',message=row.get('message','missing evidence'))
    ET.ElementTree(suite).write(output/'junit.xml',encoding='utf-8',xml_declaration=True)
    if os.environ.get('SRCSEP_REPORT_DIR'):
        # The existing source-suite merger consumes one top-level compatibility
        # report. Keep every required role in the authoritative child summary;
        # the worst status propagates for each individually selectable ID.
        directory=Path(os.environ['SRCSEP_REPORT_DIR']);directory.mkdir(parents=True,exist_ok=True)
        imported=[]
        for identity in sorted({r['id'] for r in results}):
            related=[r for r in results if r['id']==identity]
            status='ERROR' if any(r['status']=='ERROR' for r in related) else 'FAIL' if any(r['status']=='FAIL' for r in related) else 'PASS'
            imported.append(dict(id=identity,status=status,group='stage7-checkpoint-native',
                message='; '.join(r.get('message','') for r in related if r.get('message')) or 'all required native roles verified',
                artifacts=[str(output/'summary.json'),str(output/'junit.xml')]))
        (directory/'stage7-checkpoint-results.json').write_text(json.dumps(dict(schema='srcsep-component-tests-v1',results=imported),indent=2)+'\n')
    print({status:sum(r['status']==status for r in results) for status in ('PASS','FAIL','SKIP','ERROR')})
    return int(any(r['status']!='PASS' for r in results))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-root',type=Path,default=os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD'))
    parser.add_argument('--thread-build2',type=Path,default=os.environ.get('AMPS_FIELD_LINE_THREAD2_BUILD'))
    parser.add_argument('--thread-build4',type=Path,default=os.environ.get('AMPS_FIELD_LINE_THREAD4_BUILD'))
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--mesh-cache',type=Path);parser.add_argument('--fixtures-only',action='store_true')
    parser.add_argument('--prepare-only',action='store_true');parser.add_argument('--collect',action='store_true')
    args=parser.parse_args();output=args.output.resolve()
    if args.collect:return collect(json.loads((output/'prepared.json').read_text()),output)
    demand(args.build_root,'supply the matching configured --build-root')
    demand(not output.exists() or not any(output.iterdir()),'choose a fresh output directory')
    output.mkdir(parents=True,exist_ok=True);plan=prepare(args,output)
    if args.prepare_only:print('Prepared actual native roles; no execution/PASS claimed:',output/'prepared.json');return 0
    for job in plan['jobs']:
        directory=Path(job['cwd']);error=''
        name='application.rank0.log' if job['threads'] else 'native.rank0.log'
        with (directory/name).open('w') as log:
            try:
                process=subprocess.Popen(job['argv'],cwd=directory,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:code=process.wait(timeout=job['timeout_s'])
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                    code=process.returncode;error='declared 300-s MPI deadline exceeded'
            except OSError as exception:code=None;error=str(exception)
        (directory/'process-result.json').write_text(json.dumps(dict(argv=job['argv'],exit_code=code,launch_error=error))+'\n')
    return collect(plan,output)


if __name__=='__main__':raise SystemExit(main())
