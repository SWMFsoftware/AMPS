#!/usr/bin/env python3
"""Fail-closed M3 cohort preparation and conditional normal producer boundary.

The installed I21 adapter can execute only when source/state, normal compiled
producer, explicit inputs and predecessors actually exist. Current Na+ source
data remain incomplete. I22-I24 have no complete daughter provider. No guessed
CLI, Boolean activation, cached result or fixture can substitute for these.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import os

import m3_qualification as q

ROOT=q.TEST.parents[1]
LOCAL={'I21':('U14',),'I22':('U15','U24'),'I23':('U16',),'I24':('U17',)}
# A source-code integration, not a configurable JSON flag, must implement the
# normal production native-cohort adapter with complete physical providers.
# Current rate/zero-Q/controller fixtures are deliberately NOT in this list.
# This is a source-installed adapter identity, NOT a physical-readiness flag.
# It still requires source admission, actual execution and independent proofs.
PHYSICAL_PRODUCERS={'I21':'normal-source-bound-sodium-cohort-v1'}
ARTIFACT_ROLES=('run_manifest','parent_records','event_ledger','numerical_evidence',
    'budget_evidence','executable','build_manifest','input_snapshot')


def native_request(plan,base):
    """Authenticate a frozen six-or-more-launch homogeneous experiment.

    The request carries paths/hashes, not executable shell text or arbitrary
    activation flags. Both MPI layouts use the same predeclared >=3 increasing
    physical-step counts. Each native input independently binds its include
    graph and actual source/model quantities at runtime.
    """
    value=plan.get('native_adapter')
    if value is None:raise q.MissingPrerequisite('explicit frozen native adapter/inputs unavailable')
    if not isinstance(value,dict) or set(value)!={'schema','producer_id','producer_manifest','state_mapping','runs'}:
        raise ValueError('native adapter request schema/fields mismatch')
    if value['schema']!='srcMoon-M3-normal-sodium-execution-v1' or value['producer_id']!=PHYSICAL_PRODUCERS.get(plan['test_id']):
        raise ValueError('native adapter is not installed for this test ID')
    identities=[]
    def file(record,label):
        if not isinstance(record,dict) or set(record)!={'path','sha256'} or not isinstance(record['path'],str) or not record['path'] or not q.digest_syntax(record['sha256']):
            raise ValueError('explicit native path/SHA256 required: '+label)
        path=Path(record['path']);path=(path if path.is_absolute() else base/path).resolve()
        if not path.is_file():raise q.MissingPrerequisite('native '+label+' unavailable: '+str(path))
        if q.sha(path)!=record['sha256']:raise ValueError('native '+label+' identity changed')
        identities.append({'role':label,'path':str(path),'sha256':record['sha256']});return path
    producer=file(value['producer_manifest'],'compiled producer manifest')
    mapping=value['state_mapping']
    if not isinstance(mapping,dict) or set(mapping)!={'parent','ion','electron'} or any(not isinstance(v,str) or not v for v in mapping.values()):
        raise ValueError('explicit three-state source mapping required')
    runs=value['runs'];layouts={1:[],8:[]};resolved=[]
    if not isinstance(runs,list) or len(runs)<6:raise ValueError('at least three frozen refinements in BOTH1/8MPI layouts required')
    for run in runs:
        if not isinstance(run,dict) or set(run)!={'mpi_ranks','physical_steps','input','input_documents'}:raise ValueError('native launch requires explicit layout, steps and complete input identities')
        ranks=run['mpi_ranks'];steps=q.integer(run['physical_steps'],'frozen physical steps')
        if type(ranks) is not int or ranks not in layouts:raise ValueError('native layout must be1or8MPI,OMP1')
        main_input=file(run['input'],'native input');documents=run['input_documents']
        if not isinstance(documents,list) or not documents:raise ValueError('frozen complete input include graph required')
        paths=[file(record,'native input document') for record in documents]
        if len(paths)!=len(set(paths)) or main_input not in paths:raise ValueError('input include graph must be unique and include root input')
        layouts[ranks].append(steps);resolved.append({**run,'input_path':main_input,'document_paths':paths})
    if any(len(s)!=len(set(s)) or len(s)<3 for s in layouts.values()) or sorted(layouts[1])!=sorted(layouts[8]):
        raise ValueError('native matched layouts require the same distinct >=3 frozen refinements')
    return {'producer':producer,'mapping':mapping,'runs':resolved,'identities':identities}


def bind_native_evidence(path,expected_sha256,plan,plan_sha256):
    """Authenticate optional ACTUAL producer artifacts, without activation.

    This is an explicit forward integration seam; it does not recognize an
    arbitrary program as a physical producer or assign an absent physics
    criterion. A future production adapter must be installed in source and
    supply its model-specific ledger/numerical verifiers before formal use.
    """
    path=Path(path).resolve()
    if not q.digest_syntax(expected_sha256):raise ValueError('native evidence requires an explicit lowercase SHA256')
    if not path.is_file():raise q.MissingPrerequisite('actual native producer evidence manifest is unavailable')
    if q.sha(path)!=expected_sha256:raise ValueError('native producer evidence identity mismatch')
    value=q.strict_object(path)
    required={'schema','test_id','plan_sha256','approved_addendum_sha256','source_profile_sha256',
        'producer_id','execution_token','clock','mpi_ranks','openmp_threads','native_returncode','artifacts'}
    if set(value)!=required or value['schema']!='srcMoon-M3-native-cohort-evidence-v1':raise ValueError('native cohort evidence schema/field mismatch')
    if (value['test_id']!=plan['test_id'] or value['plan_sha256']!=plan_sha256 or
        value['approved_addendum_sha256']!=q.APPROVED_SHA256 or
        value['source_profile_sha256']!=plan['source_profile']['sha256']):
        raise ValueError('native cohort evidence is not bound to this exact plan/profile/protocol')
    if (not isinstance(value['producer_id'],str) or not value['producer_id'] or
        not isinstance(value['execution_token'],str) or not value['execution_token'] or
        value['clock']!='common-time' or type(value['mpi_ranks']) is not int or value['mpi_ranks'] not in (1,8) or
        type(value['openmp_threads']) is not int or value['openmp_threads']!=1 or
        type(value['native_returncode']) is not int or value['native_returncode']!=0):
        raise ValueError('native cohort completion/layout/token fields invalid')
    artifacts=value['artifacts']
    if not isinstance(artifacts,dict) or set(artifacts)!=set(ARTIFACT_ROLES):raise ValueError('native cohort evidence must bind every required actual artifact role')
    identities=[];missing=[];resolved={}
    for role in ARTIFACT_ROLES:
        record=artifacts[role]
        if not isinstance(record,dict) or set(record)!= {'path','sha256'} or not isinstance(record['path'],str) or not record['path'] or not q.digest_syntax(record['sha256']):
            raise ValueError('native artifact requires explicit path/SHA256: '+role)
        target=Path(record['path']);target=target if target.is_absolute() else path.parent/target
        if not target.is_file():missing.append(role+' artifact unavailable: '+str(target));continue
        actual=q.sha(target)
        if actual!=record['sha256']:raise ValueError('native '+role+' artifact hash mismatch')
        identities.append({'role':role,'path':str(target.resolve()),'sha256':actual});resolved[role]=target
    # An independently retained run manifest prevents a loose parent CSV from
    # being interpreted as this native launch merely because filenames agree.
    if 'run_manifest' in resolved:
        manifest=q.strict_object(resolved['run_manifest'])
        for key in ('execution_token','test_id','plan_sha256','source_profile_sha256'):
            if manifest.get(key)!=value[key]:raise ValueError('actual native run manifest contradicts '+key)
        if manifest.get('execution')!='native-production-cohort':raise ValueError('run manifest does not record native cohort execution')
    return {'manifest':value,'identities':identities,'missing':missing,
        'scope':'artifact identity only; named production provider/model verifiers still required'}


def readiness(plan_path,campaign_root,native_evidence=None,expected_evidence_sha256=None):
    """Read-only actual prerequisite assessment, preserving corrupt evidence."""
    plan_path=Path(plan_path).resolve();campaign_root=Path(campaign_root).resolve()
    protocol=q.approved_protocol()
    plan=q.strict_object(plan_path);q.validate_plan(plan)
    identity=plan['test_id'];blocks=[];identities=[]
    def bind(path,expected,label):
        path=Path(path)
        if not q.digest_syntax(expected):raise ValueError(label+' requires lowercase SHA256')
        if not path.is_file():blocks.append(label+' is unavailable: '+str(path));return None
        actual=q.sha(path)
        if actual!=expected:raise ValueError(label+' hash mismatch: '+str(path))
        identities.append({'role':label,'path':str(path),'sha256':actual});return path
    profile=plan.get('source_profile')
    if not isinstance(profile,dict) or set(profile)!= {'path','sha256'} or not isinstance(profile['path'],str) or not profile['path']:
        raise ValueError('cohort requires explicit source_profile path/SHA256 binding')
    path=Path(profile['path'])
    if not path.is_absolute():path=plan_path.parent/path
    path=bind(path,profile['sha256'],'physical source/state profile')
    if path:
        value=q.strict_object(path)
        if value.get('schema')!='srcMoon-M3-state-profile-v1':raise ValueError('unsupported actual state-profile schema')
        # These profiles authenticate available isotope/constants/states and
        # preserve unavailable ones; they provide NO complete physical channel.
        if value.get('physical_activation') is not False:
            raise ValueError('source-only state profile cannot self-activate physical chemistry')
        if plan.get('native_adapter') is None:
            blocks.append('source-state profile alone does not supply complete physical reaction providers')
    request=None
    if identity in PHYSICAL_PRODUCERS:
        try:
            request=native_request(plan,plan_path.parent);identities+=request['identities']
            # Independent source inspection precedes any costly compilation or
            # launch. It cannot infer an unavailable Na+ state or physical mass.
            if path:
                from sodium_cohort.verify_ledger import load_profile,Unavailable
                source_ids={}
                try:load_profile(path,{**request['mapping'],'profile_sha256':profile['sha256']},source_ids)
                except Unavailable as error:blocks.append(str(error))
                identities += [{'role':'source profile dependency','path':p,'sha256':h} for p,h in source_ids.items()]
        except q.MissingPrerequisite as error:blocks.append(str(error))
    if (native_evidence is None)!=(expected_evidence_sha256 is None):raise ValueError('native evidence requires BOTH explicit manifest path and expected SHA256')
    if native_evidence is not None:
        try:
            evidence=bind_native_evidence(native_evidence,expected_evidence_sha256,plan,q.sha(plan_path))
            identities+=evidence['identities'];blocks+=evidence['missing']
        except q.MissingPrerequisite as error:blocks.append(str(error))
    predecessors=[]
    for predecessor in tuple('I%02d'%i for i in range(13,21))+LOCAL[identity]:
        target=campaign_root/predecessor/'result.json'
        entry={'id':predecessor,'path':str(target),'status':'UNAVAILABLE'}
        if target.is_file():
            value=q.strict_object(target);status=value.get('status')
            if status not in ('PASS','FAIL','ERROR','SKIPPED'):raise ValueError('invalid actual predecessor status: '+str(target))
            entry.update(status=status,sha256=q.sha(target))
        if entry['status']!='PASS':blocks.append(predecessor+' has not qualified in this campaign')
        predecessors.append(entry)
    # All-green administrative predecessor fixtures CANNOT create a native
    # physical producer. Actual state/Q/kinetic/numerical/ledger evidence is
    # still required and no aggregate source/profile metadata closes it.
    if identity not in PHYSICAL_PRODUCERS:blocks.append('source-complete physical native cohort producer is not installed for '+identity)
    if request is None:blocks.append('actual independent transport/numerical-bias and physical ledger proofs remain required')
    return {'schema':'srcMoon-M3-qualification-readiness-v1','status':'SKIPPED' if blocks else 'PASS',
        'test_id':identity,'execution':'not-run','criteria_evaluated':False,
        'approved_addendum_sha256':q.APPROVED_SHA256,'plan_sha256':q.sha(plan_path),
        'identities':identities,'predecessors':predecessors,'blocking_reasons':blocks,
        'formal_M3_status':'SKIPPED / NOT QUALIFIED','observational_scoring':False,
        'resource_policy':protocol['execution']}


def execute_cohort(plan_path,plan,output,remaining_seconds):
    """Run actual installed native owner and independent auditors under ONE clock.

    Normal-source matching/build identities are rechecked by the maintained
    owner on every launch. Each nested runner receives only the remaining
    package allowance; this never resets a failed/exhausted campaign clock.
    A successful cohort does not close the whole formal I21 protocol matrix.
    """
    from campaign_budget import CampaignBudget
    budget=CampaignBudget(remaining_seconds,min(120,remaining_seconds*.1))
    request=native_request(plan,Path(plan_path).parent);frozen_plan_sha=q.sha(plan_path)
    records=[];commands=[];started=time.monotonic()
    journal=output/'execution_journal.json'
    def retain_progress():
        # This durable attempted-launch receipt survives exceptions/timeouts.
        # It records only actual commands, never fabricated science outcomes.
        journal.write_text(json.dumps({'execution':'native-attempted','native_runs':records,
            'commands':commands,'mpi_openmp':'declared1/8MPI,OMP1','plan_sha256':q.sha(plan_path)},indent=2)+'\n')
    def launch(argv,label,cores):
        directory=output/(label+'-resources');directory.mkdir()
        commands.append({'label':label,'command':argv,'state':'launched'})
        retain_progress()
        result=budget.run(argv,ROOT,output/(label+'.log'),output=directory,allocated_cores=cores,
            environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1'))
        commands[-1].update(state='completed',receipt=result);retain_progress()
        if result.get('reason'):raise RuntimeError('cohort resource termination: '+str(result['reason']))
        return result
    for index,run in enumerate(request['runs']):
        if q.sha(plan_path)!=frozen_plan_sha:raise ValueError('frozen cohort plan changed before native launch')
        for path,item in zip(run['document_paths'],run['input_documents']):
            if q.sha(path)!=item['sha256']:raise ValueError('frozen input document changed before native launch')
        elapsed=time.monotonic()-started;available=remaining_seconds-elapsed-30
        if available<=0:raise RuntimeError('original cohort deadline exhausted before next native launch; declare a cluster package')
        target=output/('native-%02d'%index)
        command=[sys.executable,'-B',str(q.TEST/'run_m3_sodium_cohort.py'),'--output-dir',str(target),
            '--producer-manifest',str(request['producer']),'--input',str(run['input_path']),
            '--mode','run','--mpi-ranks',str(run['mpi_ranks']),'--wall-seconds',str(available)]
        completion=launch(command,'launch-%02d'%index,run['mpi_ranks'])
        result=q.strict_object(target/'result.json');status=result.get('status')
        expected={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}
        if status not in expected or completion['returncode']!=expected[status]:raise ValueError('native owner status/exit contradiction')
        records.append({'output':str(target/'native-output'),'mpi_ranks':run['mpi_ranks'],'physical_steps':run['physical_steps'],
            'input':{'path':str(run['input_path']),'sha256':run['input']['sha256']},
            'owner_manifest':{'path':str(target/'run_manifest.json'),'sha256':q.sha(target/'run_manifest.json')},
            'owner_result':{'path':str(target/'result.json'),'sha256':q.sha(target/'result.json')},
            'input_documents':[{'path':str(path),'sha256':item['sha256']} for path,item in zip(run['document_paths'],run['input_documents'])]})
        retain_progress()
        if status!='PASS':return {'status':status,'reason':result.get('reason'),'execution':'native-attempted',
            'criteria_evaluated':False,'native_runs':records,'commands':commands}
    registry=output/'native_runs.json';registry.write_text(json.dumps(records,indent=2)+'\n')
    audited=launch([sys.executable,'-B',str(q.TEST/'sodium_cohort/audit.py'),'--plan',str(plan_path),
        '--runs',str(registry),'--output-dir',str(output/'independent-audit')],'independent-audit',1)
    result=q.strict_object(output/'independent-audit/result.json')
    if q.sha(plan_path)!=frozen_plan_sha:raise ValueError('frozen cohort plan changed during native execution/audit')
    if audited['returncode']!={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}.get(result.get('status')):
        raise ValueError('independent auditor status/exit contradiction')
    return {**result,'execution':'actual-native-cohort-and-independent-audit','criteria_evaluated':True,
        'native_runs':records,'commands':commands,'formal_M3_status':'NOT QUALIFIED BY ONE COHORT; full frozen I21-I24 matrix remains'}


def run(plan_path,campaign_root,output,amps=None,native_evidence=None,expected_evidence_sha256=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();plan_path=Path(plan_path).resolve()
    manifest={'schema':'srcMoon-M3-qualification-attempt-v1','created_utc':datetime.now(timezone.utc).isoformat(),
        'execution':'not-run','criteria_evaluated':False,'random_seed':None,
        'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'runner':str(Path(__file__).resolve()),'runner_sha256':q.sha(__file__),
        'scorer_sha256':q.sha(q.__file__),'approved_addendum_sha256':q.APPROVED_SHA256,
        'requested_plan':str(plan_path),'campaign_root':str(Path(campaign_root).resolve()),
        'MPI_OpenMP':'not launched; formal layouts remain1/8MPI,OMP1',
        'executable':None,'observational_scoring':False}
    if amps is not None:
        path=Path(amps).resolve();manifest['executable']={'path':str(path),'sha256':q.sha(path) if path.is_file() else None}
    try:
        if not plan_path.is_file():raise q.MissingPrerequisite('explicit frozen cohort plan is unavailable')
        shutil.copyfile(plan_path,output/'cohort_plan.json');manifest['plan_sha256']=q.sha(plan_path)
        value=readiness(plan_path,campaign_root,native_evidence,expected_evidence_sha256)
        if q.sha(plan_path)!=q.sha(output/'cohort_plan.json'):
            raise ValueError('frozen cohort plan changed during readiness assessment')
        manifest['source_data_identities']=value['identities']
        value['reason']='; '.join(value['blocking_reasons'])
        if value['status']=='PASS':
            elapsed=time.monotonic()-started
            if elapsed>=3600:raise RuntimeError('original cohort budget exhausted during prerequisite assessment')
            value=execute_cohort(plan_path,q.strict_object(plan_path),output,3600-elapsed)
    except q.MissingPrerequisite as error:
        value={'status':'SKIPPED','reason':str(error),'execution':'not-run','criteria_evaluated':False}
    except (OSError,ValueError,TypeError,KeyError,ArithmeticError,RuntimeError) as error:
        value={'status':'ERROR','reason':str(error),'execution':'not-run','criteria_evaluated':False}
    journal=output/'execution_journal.json'
    if journal.is_file():
        progress=q.strict_object(journal)
        if value['execution']=='not-run':value['execution']=progress['execution']
        manifest['native_runs']=progress['native_runs'];manifest['commands']=progress['commands']
        manifest['MPI_OpenMP']=progress['mpi_openmp']
    manifest['execution']=value['execution'];manifest['criteria_evaluated']=value['criteria_evaluated']
    if value['execution']!='not-run' and plan_path.is_file():
        actual_plan=q.strict_object(plan_path)
        manifest.update(random_seed=actual_plan['seed'],cohort_id=actual_plan['cohort_id'])
    value.update(scope='conditional per-cohort production verification; whole milestone closure remains separate',
        counts={s:int(value['status']==s) for s in ('PASS','FAIL','ERROR','SKIPPED')},
        wall_seconds=time.monotonic()-started)
    (output/'stdout.log').write_text(json.dumps(value,indent=2)+'\n')
    (output/'stderr.log').write_text(value['reason']+'\n' if value['status']=='ERROR' else '')
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (output/'result.json').write_text(json.dumps(value,indent=2)+'\n')
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-plan',type=Path,required=True)
    parser.add_argument('--campaign-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--amps',type=Path)
    parser.add_argument('--producer-evidence',type=Path,help='actual native receipt manifest; never an activation flag')
    parser.add_argument('--producer-evidence-sha256',help='explicit expected digest of that exact producer evidence')
    args=parser.parse_args();value=run(args.case_plan,args.campaign_root,args.output_dir,args.amps,
        args.producer_evidence,args.producer_evidence_sha256)
    print(json.dumps(value,indent=2));return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]


if __name__=='__main__':raise SystemExit(main())
