#!/usr/bin/env python3
"""Bounded initial M3 rate verification under the approved development exception.

This owner runs the unchanged full M0 prerequisite, registered U16 production
kernel/runtime controls, and an eight-rank real Moon rate callback. It never
enables particle reactions, assigns chemistry observational PASS, or bypasses
M1/M2. All builds belong to the selected checkout; an isolated clone is the
recommended execution location to preserve a user's active generated tree.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from campaign_budget import CampaignBudget
from integration_common import REPO
from run_integration_tests import clean_configure_build
from run_m3_rate_callback import run as run_callback


class SupportingStopped(Exception):
    """An evaluated FAIL/SKIP stops dependents without inventing an ERROR."""


def worker(output,mpi_ranks,wall_seconds,spice_path,spice_kernels,cleanup_reserve_seconds=600):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    records=[];support=[];error=None;restore_required=False
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
        if not passed:
            if baseline_status=='FAIL':raise SupportingStopped()
            raise RuntimeError('M0 prerequisite did not PASS')
        # U16 was actually invoked by the complete registered U-series. Retain
        # its original record, never manufacture a local PASS from file presence.
        local=json.loads((baseline/'u-series/run_all_summary.json').read_text())
        entry=next(row for row in local['results'] if row['id']=='U16')
        records.append(entry)
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
    except SupportingStopped:
        pass
    except (OSError,ValueError,RuntimeError,StopIteration,subprocess.TimeoutExpired) as exception:
        error=str(exception) or type(exception).__name__
    finally:
        if restore_required:
            try:
                clean_configure_build('moon',str(spice_path),
                    str(spice_kernels),output/'restore','moon-spice')
                support.append({'check':'ordinary-Moon-restore','status':'PASS',
                    'runner':str(REPO/'srcMoon/test/run_integration_tests.py'),
                    'log':str(output/'restore/build_logs/moon-spice.log')})
            except (OSError,ValueError,RuntimeError) as exception:
                error=(error+'; ' if error else '')+'restore failed: '+str(exception)
    for number in range(21,25):
        records.append({'id':'I'+str(number),'status':'SKIPPED',
            'reason':'formal M1/M2 predecessors and chemistry product/ensemble qualification unavailable; local rates do not substitute',
            'scope':'formal-chemistry-not-qualified'})
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
    args=parser.parse_args()
    if args.mpi_ranks<1:parser.error('positive MPI rank count required')
    if not 0<args.cleanup_reserve_seconds<args.wall_seconds:
        parser.error('cleanup reserve must be positive and smaller than whole wall budget')
    if args.worker:return worker(args.output_dir,args.mpi_ranks,args.wall_seconds,
        args.spice_path,args.spice_kernels,args.cleanup_reserve_seconds)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error('output must be new or empty')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--worker',
        '--output-dir',str(args.output_dir.resolve()),'--mpi-ranks',str(args.mpi_ranks),
        '--wall-seconds',str(args.wall_seconds),'--spice-path',args.spice_path,
        '--spice-kernels',args.spice_kernels,'--cleanup-reserve-seconds',str(args.cleanup_reserve_seconds)]
    budget=CampaignBudget(args.wall_seconds,args.cleanup_reserve_seconds)
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
