#!/usr/bin/env python3
"""Bounded M3 supporting verification under the approved development exception.

This owner runs the unchanged full M0 prerequisite, registered supporting
services and a real Moon rate callback. The explicit disposable-checkout flag
also selects real PIC packet, Moon force/tracer and constant-history controller
fixtures against one matching producer. Analytical controller events do not
activate a lunar physical kinetic policy or assign observational PASS. M1/M2
predecessors remain unchanged. All builds belong to the selected checkout.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from campaign_budget import CampaignBudget
from integration_common import REPO
from run_integration_tests import clean_configure_build
from run_m3_rate_callback import run as run_callback
from run_m3_particle_backend import run as run_particle_backend
from run_m3_moon_force_tracers import run as run_moon_force_tracers
from run_m3_reaction_controller import run as run_reaction_controller
from run_m3_trajectory_exposure import run as run_trajectory_exposure
from run_m3_transport_clock import run as run_transport_clock
from run_m3_transport_clock import run_startup as run_transport_startup
from run_m3_normal_cohort import run as run_semantic_cohort
from run_m3_normal_cohort import run_startup as run_semantic_startup
from run_m3_coefficient_stage import run as run_coefficient_stage
from build_m2_qualification import production_source_identity,sha256
from run_m3_moon_force_tracers import native_sources


def normal_restore_evidence(output,source_before):
    """Bind an ACTUALLY successful clean restore, never infer it from files.

    Caller invokes this only after clean_configure_build returned normally.
    Source bytes are recorded on both sides of that build; generated headers,
    archive, retained executable and complete successful command log are bound
    together. No compilation/configuration is performed by this function.
    """
    output=Path(output);after=production_source_identity()
    if native_sources(source_before)!=native_sources(after):
        raise ValueError('native source changed during ordinary Moon restore')
    log=output/'build_logs/moon-spice.log'
    lines=log.read_text().splitlines()
    statuses=[line for line in lines if line.startswith('EXIT: ')]
    commands=[line[len('COMMAND: '):] for line in lines if line.startswith('COMMAND: ')]
    if len(statuses)!=2 or statuses!=['EXIT: 0','EXIT: 0'] or len(commands)!=2:
        raise ValueError('ordinary restore lacks two successful configuration/build commands')
    current={};retained={}
    for relative in ('pic/pic.h','pic/picGlobal.dfn','pic/picSpeciesMacro.dfn',
            'pic/picParticleDataMacro.h','general/rnd.h'):
        source=REPO/'build'/relative;target=output/'generated'/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        current[str(source)]=sha256(source);retained[str(target)]=sha256(target)
    library=REPO/'build/libAMPS.a';executable=output/'executables/moon-spice/amps'
    value={'schema':'srcMoon-actual-normal-restore-producer-v1','status':'PASS',
        'source_root':str(REPO),'production_source_hashes':after,
        'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        'prebuild_production_source_hashes':source_before,
        'native_source_before_after_match':True,
        'library':{'path':str(library),'sha256':sha256(library)},
        'executable':{'path':str(executable),'sha256':sha256(executable)},
        'build_log':{'path':str(log),'sha256':sha256(log)},'build_commands':commands,
        'generated_snapshots':current,'retained_generated_snapshots':retained}
    (output/'normal_build_evidence.json').write_text(json.dumps(value,indent=2)+'\n')
    return value


def run_formal_owners(output,mpi_ranks):
    """Invoke the maintained per-ID gates, never fabricate four skip rows."""
    from integration_m3 import run
    rows=[]
    for number in range(21,25):
        identity='I'+str(number)
        value=run(identity,REPO/'amps',REPO/'srcMoon/examples/m3_reaction_rates.in',
            Path(output)/'formal',mpi_ranks,'regression')
        directory=Path(output)/'formal'/identity
        rows.append({**value,'runner':str(REPO/'srcMoon/test/integration-tests'/identity/'test.py'),
            'log':str(directory/'stdout.log'),'result_json':str(directory/'result.json'),
            'run_manifest':str(directory/'run_manifest.json')})
    return rows


class SupportingStopped(Exception):
    """An evaluated FAIL/SKIP stops dependents without inventing an ERROR."""


def worker(output,mpi_ranks,wall_seconds,spice_path,spice_kernels,cleanup_reserve_seconds=600,
        include_particle_backend=False):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    records=[];support=[];error=None;restore_required=False;restored_evidence=None
    def logged(command,path):
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w') as stream:
            stream.write('COMMAND: '+str(command)+'\n');stream.flush()
            process=subprocess.run(command,cwd=REPO,stdout=stream,stderr=subprocess.STDOUT,
                timeout=max(1,wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)))
        return process.returncode
    try:
        baseline=output/'m0'
        command=[sys.executable,'-B',str(REPO/'srcMoon/test/run_integration_tests.py'),
                 '--output-dir',str(baseline),'--spice-path',str(spice_path),
                 '--spice-kernels',str(spice_kernels)]
        code=logged(command,output/'m0.log')
        value=json.loads((baseline/'baseline_summary.json').read_text())
        passed=code==0 and value.get('m0_gate_status')=='PASS'
        # M0's aggregate gate is FAIL for BOTH evaluated failures and build/
        # runtime errors. Inspect its retained records/counts before assigning
        # our outcome; a failed build is never scientific disagreement.
        baseline_error=(value.get('counts',{}).get('ERROR',0)>0 or any(
            row.get('status')=='ERROR' for row in
            value.get('results',[])+value.get('local_prerequisite_results',[])))
        baseline_status=('ERROR' if baseline_error else 'PASS' if passed
            else 'FAIL' if value.get('m0_gate_status')=='FAIL' else 'ERROR')
        support.append({'check':'M0-linked-prerequisite','status':baseline_status,
            'runner':command[2],'log':str(output/'m0.log'),'summary':str(baseline/'baseline_summary.json')})
        if baseline_status!='PASS':
            if baseline_status=='FAIL':raise SupportingStopped()
            raise RuntimeError('M0 prerequisite did not PASS')
        # U16 was actually invoked by the complete registered U-series. Retain
        # its original record, never manufacture a local PASS from file presence.
        local=json.loads((baseline/'u-series/run_all_summary.json').read_text())
        entry=next(row for row in local['results'] if row['id']=='U16')
        records.append(entry)
        # The complete M0 local series now actually invokes documented photon
        # selection, event primitives, input controls and CEX rate services
        # through U15/U17. Preserve those original records and child artifacts;
        # their component PASS does not change their unqualified full scope.
        for identity in ('U15','U17'):
            component=next(row for row in local['results'] if row['id']==identity)
            records.append(component)
            support.extend(component.get('supporting_checks',[]))
        if entry['status']!='PASS':
            if entry['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            raise RuntimeError('registered U16 did not PASS')
        remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
        if remaining<1:raise RuntimeError('insufficient budget for native callback and restore')
        restore_required=True
        callback=run_callback(output/'native-rate-callback',
            REPO/'srcMoon/examples/m3_atomic_rates.in',mpi_ranks=mpi_ranks,
            wall_seconds=remaining,spice_path=spice_path,spice_kernels=spice_kernels)
        support.append({'check':'U16-native-production-callback','status':callback['status'],
            'runner':str(REPO/'srcMoon/test/run_m3_rate_callback.py'),
            'log':str(output/'native-rate-callback/stdout.log'),
            'result':str(output/'native-rate-callback/result.json')})
        if callback['status']=='ERROR':raise RuntimeError('native rate callback errored')
        if callback['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
        if include_particle_backend:
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for particle backend and normal restore')
            backend=run_particle_backend(output/'native-particle-backend',
                isolated_checkout=True,wall_seconds=remaining,build_jobs=8,
                spice_path=spice_path,spice_kernels=spice_kernels)
            support.append({'check':'M3-real-PIC-packet-backend','status':backend['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_particle_backend.py'),
                'log':str(output/'native-particle-backend/build.log'),
                'result':str(output/'native-particle-backend/result.json')})
            if backend['status']=='ERROR':raise RuntimeError('native particle backend errored')
            if backend['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for actual Moon force tracers and restore')
            force=run_moon_force_tracers(output/'native-moon-force-tracers',
                output/'native-particle-backend/run_manifest.json',wall_seconds=remaining)
            support.append({'check':'M3-actual-Moon-force-tracer-integration','status':force['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_moon_force_tracers.py'),
                'log':str(output/'native-moon-force-tracers/compile-probe.log'),
                'result':str(output/'native-moon-force-tracers/result.json')})
            if force['status']=='ERROR':raise RuntimeError('actual Moon force tracer checks errored')
            if force['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            # Reuse only this exact fresh producer, before ordinary restoration.
            # The controller runner independently checks library/source/ABI
            # identity and uses the actual normal Moon/PIC dispatch. A new
            # child never receives a fresh one-hour work allowance.
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for reaction controller and restore')
            controller=run_reaction_controller(output/'native-reaction-controller',
                output/'native-particle-backend/run_manifest.json',wall_seconds=remaining)
            support.append({'check':'M3-normal-path-constant-history-controller','status':controller['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_reaction_controller.py'),
                'log':str(output/'native-reaction-controller/compile-probe.log'),
                'result':str(output/'native-reaction-controller/result.json')})
            if controller['status']=='ERROR':raise RuntimeError('native reaction controller errored')
            if controller['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for actual trajectory exposure and restore')
            exposure=run_trajectory_exposure(output/'native-trajectory-exposure',
                output/'native-particle-backend/run_manifest.json',wall_seconds=remaining)
            support.append({'check':'M3-normal-trajectory-exposure-controller','status':exposure['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_trajectory_exposure.py'),
                'log':str(output/'native-trajectory-exposure/compile.log'),
                'result':str(output/'native-trajectory-exposure/result.json')})
            if exposure['status']=='ERROR':raise RuntimeError('native trajectory exposure checks errored')
            if exposure['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            # This final fixture uses that exact producer and normal Moon
            # source/mover/Begin/Commit adapters. It proves common chronology
            # separately from frozen environmental ET, without enabling any
            # physical chemistry or borrowing a new per-child wall allowance.
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for common transport clock and restore')
            clock=run_transport_clock(output/'native-common-transport-clock',
                output/'native-particle-backend/run_manifest.json',wall_seconds=remaining)
            support.append({'check':'M3-normal-common-transport-clock','status':clock['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_transport_clock.py'),
                'log':str(output/'native-common-transport-clock/compile-native.log'),
                'result':str(output/'native-common-transport-clock/result.json')})
            if clock['status']=='ERROR':raise RuntimeError('native common transport clock errored')
            if clock['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
            # Explicit source-SI charge admission owns a DIFFERENT analytical
            # normal application build. It must run only after every consumer
            # of the legacy producer above has finished, and before the one
            # ordinary restoration below. It never borrows another hour or
            # enables the unavailable physical sodium source/state profile.
            remaining=wall_seconds-cleanup_reserve_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('insufficient budget for semantic PIC composition and restore')
            semantic=run_semantic_cohort(output/'native-semantic-cohort',
                isolated_checkout=True,wall_seconds=remaining,build_jobs=8,
                spice_path=spice_path,spice_kernels=spice_kernels)
            support.append({'check':'M3-source-SI-semantic-PIC-composition','status':semantic['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_normal_cohort.py'),
                'log':str(output/'native-semantic-cohort/build.log'),
                'result':str(output/'native-semantic-cohort/result.json')})
            if semantic['status']=='ERROR':raise RuntimeError('native semantic PIC composition errored')
            if semantic['status'] in ('FAIL','SKIPPED'):raise SupportingStopped()
        else:
            support.append({'check':'M3-real-PIC-packet-backend','status':'SKIPPED',
                'runner':str(REPO/'srcMoon/test/run_m3_particle_backend.py'),
                'reason':'use --isolated-checkout in a disposable checkout to include the derived local-weight/time build'})
    except SupportingStopped:
        pass
    except (OSError,ValueError,RuntimeError,StopIteration,subprocess.TimeoutExpired) as exception:
        error=str(exception) or type(exception).__name__
    finally:
        if restore_required:
            try:
                before=production_source_identity()
                clean_configure_build('moon',str(spice_path),
                    str(spice_kernels),output/'restore','moon-spice')
                support.append({'check':'ordinary-Moon-restore','status':'PASS',
                    'runner':str(REPO/'srcMoon/test/run_integration_tests.py'),
                    'log':str(output/'restore/build_logs/moon-spice.log')})
                if include_particle_backend:
                    restored_evidence=normal_restore_evidence(output/'restore',before)
            except (OSError,ValueError,RuntimeError) as exception:
                error=(error+'; ' if error else '')+'restore failed: '+str(exception)
    # Full normal startup is a distinct integration layer after the manual
    # native hook fixture. Reuse the ordinary restored binary/archive; do not
    # reconfigure or start a second build. An upstream failure stops it, while
    # all actual earlier outcomes and formal skips remain preserved.
    if include_particle_backend and restored_evidence is not None and not error and all(row['status']=='PASS' for row in support):
        try:
            remaining=wall_seconds-(time.monotonic()-started)-min(60,cleanup_reserve_seconds/4)
            if remaining<=0:raise RuntimeError('insufficient original package budget for normal common-clock startup')
            startup=run_transport_startup(output/'normal-common-clock-startup',
                wall_seconds=remaining,build_evidence=restored_evidence)
            support.append({'check':'M3-full-normal-common-clock-startup','status':startup['status'],
                'runner':str(REPO/'srcMoon/test/run_m3_transport_clock.py'),
                'log':str(output/'normal-common-clock-startup/compile-startup.log'),
                'result':str(output/'normal-common-clock-startup/result.json')})
            if startup['status']=='ERROR':error='full normal common-clock startup errored'
            if startup['status']=='PASS':
                # The SAME authenticated ordinary archive/headers now prove
                # the actual pre-buffer callback and explicit-epoch SPICE
                # field in both orders and MPI1/8. This is a thin driver;
                # it does not configure/build or install a physical provider.
                remaining=wall_seconds-(time.monotonic()-started)-min(60,cleanup_reserve_seconds/4)
                if remaining<=0:raise RuntimeError('insufficient original budget for normal semantic registration')
                registration=run_semantic_startup(output/'normal-semantic-registration',
                    wall_seconds=remaining,build_evidence=restored_evidence)
                support.append({'check':'M3-full-normal-semantic-registration','status':registration['status'],
                    'runner':str(REPO/'srcMoon/test/run_m3_normal_cohort.py'),
                    'log':str(output/'normal-semantic-registration/compile-probe.log'),
                    'result':str(output/'normal-semantic-registration/result.json')})
                if registration['status']=='ERROR':error='full normal semantic registration errored'
                if registration['status']=='PASS':
                    # This later Moon-only extraction needs its OWN fresh
                    # matching normal archive. Never attach it to an older
                    # phase-two producer that predates the stage API. The
                    # wrapper checks the same source/header/library proof and
                    # runs actual normal SPICE stages without another build.
                    remaining=wall_seconds-(time.monotonic()-started)-min(60,cleanup_reserve_seconds/4)
                    if remaining<=0:raise RuntimeError('insufficient original budget for normal coefficient stage')
                    stage=run_coefficient_stage(output/'normal-coefficient-stage',
                        restored_evidence,wall_seconds=remaining)
                    support.append({'check':'M3-full-normal-coefficient-stage','status':stage['status'],
                        'runner':str(REPO/'srcMoon/test/run_m3_coefficient_stage.py'),
                        'log':str(output/'normal-coefficient-stage/compile.log'),
                        'result':str(output/'normal-coefficient-stage/result.json')})
                    if stage['status']=='ERROR':error='full normal coefficient stage errored'
        except (OSError,ValueError,RuntimeError,KeyError) as exception:
            error='full normal common-clock startup: '+str(exception)
    try:
        records.extend(run_formal_owners(output,mpi_ranks))
    except (OSError,ValueError,RuntimeError,KeyError) as exception:
        error=(error+'; ' if error else '')+'formal chemistry owner failed: '+str(exception)
    summary={'schema':'srcMoon-initial-M3-development-v1','results':records,
        'supporting_checks':support,'error':error,'formal_M3_status':'SKIPPED / NOT QUALIFIED',
        'wall_seconds':time.monotonic()-started,'observational_scoring':False,
        'counts':{s:sum(row['status']==s for row in records) for s in ('PASS','FAIL','SKIPPED','ERROR')},
        'supporting_counts':{s:sum(row['status']==s for row in support) for s in ('PASS','FAIL','SKIPPED','ERROR')}}
    summary['status']=('ERROR' if error or summary['counts']['ERROR'] or summary['supporting_counts']['ERROR']
        else 'FAIL' if summary['counts']['FAIL'] or summary['supporting_counts']['FAIL'] else 'SKIPPED')
    (output/'development_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2));return 2 if summary['status']=='ERROR' else 1 if summary['status']=='FAIL' else 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--mpi-ranks',type=int,default=8)
    parser.add_argument('--wall-seconds',type=float,default=3600)
    parser.add_argument('--cleanup-reserve-seconds',type=float,default=600)
    parser.add_argument('--spice-path',default='/home/vtenishe/SPICE/cspice')
    parser.add_argument('--spice-kernels',default='/home/vtenishe/SPICE/Kernels')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    parser.add_argument('--isolated-checkout',action='store_true',
        help='include actual PIC packet, Moon-force and reaction-controller fixtures; acknowledge this disposable checkout owns derived build/config changes')
    args=parser.parse_args()
    if args.mpi_ranks<1:parser.error('positive MPI rank count required')
    if not 0<args.cleanup_reserve_seconds<args.wall_seconds:
        parser.error('cleanup reserve must be positive and smaller than whole wall budget')
    if args.worker:return worker(args.output_dir,args.mpi_ranks,args.wall_seconds,
        args.spice_path,args.spice_kernels,args.cleanup_reserve_seconds,args.isolated_checkout)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error('output must be new or empty')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker',
        '--output-dir',str(args.output_dir.resolve()),'--mpi-ranks',str(args.mpi_ranks),
        '--wall-seconds',str(args.wall_seconds),'--spice-path',args.spice_path,
        '--spice-kernels',args.spice_kernels,'--cleanup-reserve-seconds',str(args.cleanup_reserve_seconds)]
    if args.isolated_checkout:command.append('--isolated-checkout')
    # The worker stops new work at wall-cleanup_reserve, then restores the
    # ordinary application inside that reserve. Giving the supervisor that
    # SAME reserve would kill the worker before restoration could use it.
    # Retain a smaller final termination/report allowance within the one
    # original deadline; the worker still owns the declared cleanup allowance.
    termination_reserve=min(60.0,args.cleanup_reserve_seconds/2)
    budget=CampaignBudget(args.wall_seconds,termination_reserve)
    execution=budget.run(command,REPO,args.output_dir/'campaign.log',output=args.output_dir)
    (args.output_dir/'execution_budget.json').write_text(json.dumps(execution,indent=2,sort_keys=True)+'\n')
    # Supervisor failure must be visible even if a killed child could not
    # finish its ordinary result. Never synthesize scientific success here.
    summary_path=args.output_dir/'development_summary.json'
    problem=None
    try:
        child=json.loads(summary_path.read_text())
        if not isinstance(child,dict) or child.get('status') not in ('PASS','FAIL','SKIPPED','ERROR'):
            raise ValueError('invalid M3 worker status/schema')
        for rows_key,counts_key in (('results','counts'),('supporting_checks','supporting_counts')):
            rows=child[rows_key]
            if not isinstance(rows,list) or any(not isinstance(row,dict) or
                    row.get('status') not in ('PASS','FAIL','SKIPPED','ERROR') for row in rows):
                raise ValueError('invalid M3 worker records')
            if child[counts_key]!={s:sum(row['status']==s for row in rows)
                    for s in ('PASS','FAIL','SKIPPED','ERROR')}:
                raise ValueError('M3 worker counts contradict retained records')
        if child.get('schema')!='srcMoon-initial-M3-development-v1':
            raise ValueError('unexpected M3 worker summary schema')
        derived=('ERROR' if child.get('error') or child['counts']['ERROR'] or child['supporting_counts']['ERROR']
            else 'FAIL' if child['counts']['FAIL'] or child['supporting_counts']['FAIL'] else 'SKIPPED')
        if child['status']!=derived:
            raise ValueError('M3 worker aggregate status hides a retained outcome')
        formal={row.get('id'):row['status'] for row in child['results']
                if row.get('id') in ('I21','I22','I23','I24')}
        if formal!={key:'SKIPPED' for key in ('I21','I22','I23','I24')}:
            raise ValueError('initial M3 supporting campaign must retain four formal chemistry SKIPs')
        expected={'PASS':0,'SKIPPED':0,'FAIL':1,'ERROR':2}[child['status']]
        if execution.get('reason') or execution['returncode']!=expected:
            raise ValueError('M3 worker exit/status or supervisor budget contradiction')
    except (OSError,ValueError,KeyError,TypeError) as exception:
        problem=str(exception)
    if problem is not None:
        # Preserve an existing malformed/contradictory summary for diagnosis
        # before publishing the aggregate infrastructure failure.
        if summary_path.is_file():
            summary_path.rename(args.output_dir/'invalid-worker-summary.json')
        summary_path.write_text(json.dumps({
            'schema':'srcMoon-initial-M3-development-v1','status':'ERROR',
            'error':execution.get('reason') or problem,
            'results':[],'supporting_checks':[{'check':'M3-supervisor','status':'ERROR',
                'reason':execution.get('reason') or problem,'log':str(args.output_dir/'campaign.log')}],
            'formal_M3_status':'SKIPPED / NOT QUALIFIED',
            'counts':{'PASS':0,'FAIL':0,'SKIPPED':0,'ERROR':0},
            'supporting_counts':{'PASS':0,'FAIL':0,'SKIPPED':0,'ERROR':1}},indent=2)+'\n')
        return 2
    return expected


if __name__=='__main__':raise SystemExit(main())
