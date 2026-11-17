#!/usr/bin/env python3
"""Explicit, scheduler-neutral long M2 planning/qualification entry.

Default invocation exports the fixed job registry only. --execute requires
real same-campaign numerical/environment predecessors, READY exact data,
frozen sites/normalizations/statistical support and a measured cost forecast.
No scheduler, cluster address, observation epoch or missing datum is guessed.
"""
from __future__ import annotations
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from campaign_budget import CampaignBudget
from integration_common import REPO,sha256,parse_receipt,moon_setting
from build_m2_qualification import verify_retained_executable
from m2_protocol import cases,execute_design,load_protocol
from m2_equivalence import score_campaign
from m2_initial_state_evidence import (Unavailable as InitialStateUnavailable,
    document as initial_state_document, load_qualification,
    verify_initialization_forecast,bind_frozen_sites,verify_normalization_physics)


def preflight(args):
    """Unavailable evidence SKIPPED; altered bytes or malformed proof ERROR."""
    gaps=[];evidence={};errors=[]
    if args.execute:
        for key in ('memory_limit_bytes','artifact_limit_bytes','host_free_floor_bytes','disk_free_floor_bytes'):
            value=getattr(args,key)
            if value is None:gaps.append('explicit cluster '+key+' required')
            elif value<0 or ('limit' in key and value==0):errors.append('invalid cluster '+key)

    required=('I13','I14','I15','I16','I17','I18')
    for test in required:
        if args.qualification_root is None:
            gaps.append(test+' same-campaign PASS unavailable');continue
        directory=args.qualification_root/test
        try:
            result=json.loads((directory/'result.json').read_text())
            manifest=json.loads((directory/'run_manifest.json').read_text())
            if result.get('status')!='PASS' or result.get('execution')=='not-run':
                gaps.append(test+' is not a linked PASS');continue
            # Supporting adapters intentionally reuse the frozen IDs. Their
            # PASS cannot be copied into a prerequisite directory and become
            # formal evidence merely because a receipt/executable also exists.
            scope=' '.join(str(value) for value in (
                result.get('scope',''),manifest.get('execution_scope',''))).lower()
            if (manifest.get('formal_test_status')=='SKIPPED' or
                any(word in scope for word in ('supporting','development','preflight','source-only'))):
                gaps.append(test+' evidence is explicitly supporting/preflight, not formal linked qualification');continue
            receipt=directory/'receipt.txt'
            if not receipt.is_file() or not manifest.get('executable',{}).get('sha256'):
                gaps.append(test+' linked receipt/build identity unavailable');continue
            native=parse_receipt(receipt)
            if native.get('test_id')!=test or native.get('mpi_ranks')!=manifest.get('mpi_ranks'):
                errors.append(test+' native receipt identity/layout disagrees with its manifest');continue
            evidence[test]={'result_sha256':sha256(directory/'result.json'),
                            'manifest_sha256':sha256(directory/'run_manifest.json'),
                            'receipt_sha256':sha256(receipt)}
        except FileNotFoundError:gaps.append(test+' prerequisite files unavailable')
        except (OSError,ValueError,KeyError) as error:errors.append(str(error))
    # Existing authority checks verify real raw and processed hashes, without
    # reading any model-observation residual or changing frozen membership.
    sys.path.insert(0,str(Path(__file__).parent/"data_processing"))
    from audit_m2_packages import inspect
    audit=inspect()
    for record in audit['records']:
        if record['status']=='ERROR':errors.append(record['reason'])
        elif record['status']!='PASS':gaps.append(record['check']+': '+record['reason'])
    documents={}
    for key in ('sites_file','forecast_file','normalizations_file','statistical_support_file'):
        path=getattr(args,key)
        if path is None:gaps.append(key+' required before long execution');continue
        try:
            documents[key]=initial_state_document(path);evidence[key]={'path':str(path),'sha256':sha256(path)}
        except FileNotFoundError:gaps.append(key+' unavailable')
        except (OSError,ValueError) as error:errors.append(str(error))
    if args.amps is None:gaps.append('explicit retained M2 executable required')
    elif args.execute:
        try:evidence['build']=verify_retained_executable(args.amps.resolve(),require_current_generated_tree=False)
        except (OSError,ValueError,KeyError) as error:errors.append(str(error))
    forecast=documents.get('forecast_file')
    if forecast:
        if forecast.get('status')!='READY' or not forecast.get('measured_evidence_sha256'):
            gaps.append('measured forecast/provenance is not READY')
        else:
            cost=float(forecast['allocated_core_hours'])
            if not math.isfinite(cost) or not 0<cost<=args.core_hours:
                gaps.append('measured forecast exceeds selected aggregate core-hour allowance')
            if forecast.get('protocol_sha256')!=sha256(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')):
                errors.append('forecast protocol identity differs from approved design')
            if args.amps and forecast.get('executable_sha256')!=sha256(args.amps):
                errors.append('forecast executable identity differs')
    sites=documents.get('sites_file')
    if sites:
        if sites.get('status')!='READY' or not sites.get('selection_frozen') or not sites.get('collocation_contract_sha256'):
            gaps.append('native-site selection/collocation proof not READY/frozen')
        elif not 1<=len(sites.get('sites',[]))<=48:errors.append('frozen site registry must have 1..48 sites')
        elif any('logical_surface_element' not in r or 'native_product_id' not in r for r in sites['sites']):
            errors.append('site registry lacks native identity/logical element mapping')
        else:
            evidence['shared_model_support']=(len({r['logical_surface_element'] for r in sites['sites']})!=len(sites['sites']))
    norms=documents.get('normalizations_file')
    if norms:
        if norms.get('status')!='READY' or not norms.get('freeze_before_residuals'):
            gaps.append('normalization reference must be READY and frozen before residuals')
        else:
            for case in ('source-fed','cohort'):
                for species in ('Na','Ar'):
                    value=float(norms[case][species])
                    if not math.isfinite(value) or value<=0:errors.append('invalid fixed physical normalization')
            if norms.get('source_configuration_sha256')!=sha256(args.input):errors.append('normalization input identity differs')
            if norms.get('protocol_sha256')!=sha256(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')):errors.append('normalization protocol identity differs')
            if args.amps and norms.get('executable_sha256')!=sha256(args.amps):errors.append('normalization executable identity differs')
            if args.sites_file and norms.get('sites_sha256')!=sha256(args.sites_file):errors.append('normalization frozen site identity differs')
            if args.amps is None:
                gaps.append('explicit executable required for native normalization physics binding')
            else:
                try:
                    evidence['normalization_physics']=verify_normalization_physics(norms,args.input,args.amps,
                        base=Path(args.normalizations_file).resolve().parent,mpi_ranks=args.mpi_ranks)
                except InitialStateUnavailable as error:gaps.append(str(error))
                except (OSError,ValueError,KeyError,TypeError) as error:errors.append(str(error))
    if forecast and forecast.get('mpi_ranks')!=args.mpi_ranks:
        gaps.append('cost forecast must cover the selected MPI layout')
    support=documents.get('statistical_support_file')
    if support and (support.get('status')!='READY' or not support.get('marginal_intervals_justified') or not support.get('rng_stream_policy_reviewed')):
        gaps.append('statistical interval/stream assumptions lack required qualification')
    # A generic epoch-level READY claim can accompany a completely different
    # applied t0 state. Require the real pilot's native FULL state, method,
    # canonical physical input, binary/source and explicit reviewed criteria.
    # This is model-initial-condition authentication, never observations or
    # independent thermal validation. No qualification is generated here.
    if args.amps is None or args.statistical_support_file is None:
        gaps.append('exact applied native initial-state qualification unavailable')
    else:
        try:
            qualified=load_qualification(support or {},args.input,args.amps,
                base=Path(args.statistical_support_file).resolve().parent)
            qualified=bind_frozen_sites(qualified,(sites or {}).get('sites'))
            evidence['initial_state_qualification']=qualified
            if forecast is None:
                raise InitialStateUnavailable('measured initialization cost forecast unavailable')
            evidence['initialization_cost']=verify_initialization_forecast(forecast,qualified,
                mpi_ranks=args.mpi_ranks,case_count=len(cases('cluster')))
        except InitialStateUnavailable as error:gaps.append(str(error))
        except (OSError,ValueError,KeyError,TypeError) as error:errors.append(str(error))
    return {'status':'ERROR' if errors else 'SKIPPED' if gaps else 'PASS',
            'scope':'preflight-only','gaps':gaps,'errors':errors,'evidence':evidence,'documents':documents}


def main():
    entered=time.monotonic()
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--input',type=Path,default=REPO/'srcMoon/examples/lola_surface.in')
    p.add_argument('--amps',type=Path);p.add_argument('--qualification-root',type=Path)
    for key in ('sites','forecast','normalizations','statistical-support'):p.add_argument('--'+key+'-file',type=Path)
    for key in ('memory-limit','artifact-limit','host-free-floor','disk-free-floor'):
        p.add_argument('--'+key+'-bytes',type=int)
    p.add_argument('--execute',action='store_true');p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--mpi-ranks',type=int,default=8);p.add_argument('--workers',type=int,choices=(1,2,4),default=2)
    p.add_argument('--allocated-cores',type=int);p.add_argument('--core-hours',type=float,default=300)
    a=p.parse_args();output=a.output_dir.resolve()
    if output.exists() and any(output.iterdir()):p.error('output directory must be new or empty')
    if not math.isfinite(a.core_hours) or a.core_hours<=0 or not 1<=a.mpi_ranks<=1009:
        p.error('positive finite allocation and 1..1009 ranks required by frozen seed stride')
    allocated=a.allocated_cores or a.mpi_ranks*a.workers
    if allocated<a.mpi_ranks*a.workers:p.error('allocated cores cannot be smaller than simultaneous MPI slots')
    output.mkdir(parents=True,exist_ok=True)
    plan={'schema':'srcMoon-M2-cluster-job-plan-v1','planning_complete':True,'execution':'not-run',
          'protocol_sha256':sha256(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')),
          'cases':cases('cluster'),'planned_cases':288,'allocated_cores':allocated,
          'aggregate_allocated_core_hour_budget':a.core_hours,'executor':'scheduler-neutral MPI; no remote submission'}
    (output/'job_plan.json').write_text(json.dumps(plan,indent=2,sort_keys=True)+'\n')
    try:gate=preflight(a)
    except (OSError,ValueError,KeyError,RuntimeError) as error:
        gate={'status':'ERROR','scope':'preflight-only','gaps':[],'errors':[str(error)],'evidence':{},'documents':{}}
    (output/'preflight.json').write_text(json.dumps(gate,indent=2,sort_keys=True)+'\n')
    summary={'schema':'srcMoon-M2-cluster-summary-v1','results':[], 'supporting_checks':[],
             'planned_cases':288,'formal_m2_status':'NOT QUALIFIED','execution':'not-run'}
    if not a.execute or gate['status']!='PASS':
        status='ERROR' if gate['status']=='ERROR' else 'SKIPPED'
        reason='; '.join(gate['errors']+gate['gaps']) if a.execute else 'cluster execution not selected; plan only'
        summary['results']=[{'id':test,'status':status,'execution':'not-run','scope':'formal-qualification',
            'reason':reason,'runner':str(Path(__file__)),'log':str(output/'preflight.json')} for test in ('I19','I20')]
    elif not a.worker:
        # The charged allocation is the WHOLE selected cluster bundle, including
        # coordinator/setup/plots/retries inside it, not a free budget per seed.
        wall=a.core_hours/allocated*3600
        budget=CampaignBudget(wall,max(.1,min(600,wall*.05)),core_hours=a.core_hours,
            memory_bytes=a.memory_limit_bytes,artifact_bytes=a.artifact_limit_bytes,
            host_free_floor_bytes=a.host_free_floor_bytes,disk_free_floor_bytes=a.disk_free_floor_bytes)
        budget.started=entered # Include preflight/coordinator time in the allocation.
        command=[sys.executable,'-B',str(Path(__file__).resolve()),'--execute','--worker',
            '--output-dir',str(output/'campaign'),'--input',str(a.input.resolve()),
            '--amps',str(a.amps.resolve()),'--qualification-root',str(a.qualification_root.resolve()),
            '--mpi-ranks',str(a.mpi_ranks),'--workers',str(a.workers),
            '--allocated-cores',str(allocated),'--core-hours',str(a.core_hours)]
        for key in ('sites_file','forecast_file','normalizations_file','statistical_support_file','memory_limit_bytes','artifact_limit_bytes','host_free_floor_bytes','disk_free_floor_bytes'):
            value=getattr(a,key)
            command.extend(['--'+key.replace('_','-'),str(value.resolve() if isinstance(value,Path) else value)])
        environment=dict(os.environ,AMPS_MOON_WORK_DEADLINE=str(budget.started+budget.wall_s-budget.reserve_s))
        execution=budget.run(command,REPO,output/'campaign.log',output=output,
                             allocated_cores=allocated,environment=environment)
        child=output/'campaign/cluster_summary.json'
        if child.is_file():summary=json.loads(child.read_text())
        if execution['returncode']!=0:
            summary['supporting_checks'].append({'check':'cluster-execution','status':'ERROR','scope':'infrastructure',
                'runner':str(Path(__file__)),'log':str(output/'campaign.log'),'reason':execution['reason'] or 'cluster worker failed'})
        summary['allocation_accounting']=execution
    else:
        documents=gate['documents'];elements=tuple(r['logical_surface_element'] for r in documents['sites_file']['sites'])
        partial=[]
        def checkpoint(result):
            partial.append(result);(output/'case_results.json').write_text(json.dumps(partial,indent=2)+'\n')
        records=execute_design(a.amps.resolve(),a.input.resolve(),output/'cases',tier='cluster',
            mpi_ranks=a.mpi_ranks,workers=a.workers,
            remaining_seconds=lambda:max(.01,float(os.environ['AMPS_MOON_WORK_DEADLINE'])-time.monotonic()),
            on_result=checkpoint,elements=elements,horizon_s=load_protocol()['physical_duration_s'],
            memory_limit_bytes=a.memory_limit_bytes)
        score=score_campaign(output/'cases',documents['normalizations_file'],documents['statistical_support_file'],
            application_input=a.input,executable=a.amps,support_base=Path(a.statistical_support_file).resolve().parent,
            sites=documents['sites_file']['sites'])
        (output/'numerical_score.json').write_text(json.dumps(score,indent=2,sort_keys=True)+'\n')
        summary.update(execution='executed',supporting_checks=records+[dict(check='full-numerical-equivalence',**score)])
        # Formal I20 additionally needs its independently collocated physical
        # comparison; a numerical score is never promoted to that validation.
        summary['results']=[{'id':'I19','status':score['status'],'scope':'linked-full-inventory-numerical-design',
            'runner':str(Path(__file__)),'log':str(output/'numerical_score.json')},
            {'id':'I20','status':'SKIPPED','scope':'formal-qualification',
             'reason':'numeric cohort score does not itself supply the independent PSR/thermal observation comparison',
             'runner':str(Path(__file__)),'log':str(output/'numerical_score.json')}]
    if a.execute and 'allocation_accounting' not in summary:
        summary['allocation_accounting']={'allocated_cores':allocated,
            'wall_seconds':time.monotonic()-entered,
            'allocated_core_hours':allocated*(time.monotonic()-entered)/3600,
            'scope':'selected bundle coordinator/preflight or worker allocation; retries must retain all receipts'}
    (output/'cluster_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    all_rows=summary['results']+summary['supporting_checks']
    print(json.dumps({'planned_cases':288,'execution':summary['execution'],'summary':str(output/'cluster_summary.json')},indent=2))
    return 2 if any(r['status']=='ERROR' for r in all_rows) else 1 if any(r['status']=='FAIL' for r in all_rows) else 0

if __name__=='__main__':raise SystemExit(main())
