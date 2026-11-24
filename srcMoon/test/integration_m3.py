#!/usr/bin/env python3
"""Own the existing I21--I24 formal gates and retain their actual scope.

The I21 owner can execute a pinned fixed or moving cohort. I24 can exercise a
specified sodium event-plan cohort and a separate linked PIC mover matrix.
Each owner retains the original formal prerequisites and reports the status
of each supporting execution.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
from campaign_budget import CampaignBudget
from m3_qualification import MissingPrerequisite, digest_syntax, strict_object

ROOT=Path(__file__).resolve().parents[2]
TEST=Path(__file__).resolve().parent
CAPABILITIES=TEST/'M3_TRANSPORT_CAPABILITIES.json'
STATUSES=('PASS','FAIL','ERROR','SKIPPED')


def sha(path):
    """Byte identity only; matching bytes are not scientific readiness."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def i24_supporting(directory, result, manifest, lock, inventory, producer, ranks, deadline,refinement=False):
    """Supervise each requested I24 CLI and bind its independent artifacts.

    The event-plan cohort and the native PIC matrix measure different things.
    One missing dependency or failed check must not erase evidence from the
    other. The formal qualification clock started before predecessor checks;
    both children share the remaining seconds and the same cleanup reserve.
    """
    remaining = deadline - time.monotonic()
    budget = (CampaignBudget(remaining, min(4.0, remaining/10)) if remaining > 0 else None)
    choices = [('i24_sodium_provider', 'I24 sodium MPI provider', 'run_ground_capture.py',
                'sodium-ground-capture', '--test-runtime-seconds', ['--mpi-ranks', str(ranks)])]
    if producer is not None:
        pic_args = ['--producer-manifest', str(producer.resolve())]
        if refinement:
            pic_args.append('--refinement')
        choices.append(('i24_pic_commissioning', 'I24 sodium PIC commissioning', 'run_pic_ground_capture.py',
                        'sodium-pic-ground-capture', '--wall-seconds', pic_args))
    rows, failed_receipts = [], []
    for key, label, script, name, wall_flag, extra in choices:
        runner = TEST/'integration-tests/I24'/script
        child = directory/name
        log = directory/(name+'.log')
        log.touch()
        entry = {'runner':str(runner), 'log':str(log), 'directory':str(child)}
        manifest[key] = entry
        try:
            allowance = budget.remaining() if budget else 0
            if allowance <= 0:
                raise RuntimeError('I24 shared development deadline exhausted before this launch')
            command = [sys.executable, '-B', str(runner), '--lock', str(lock.resolve()),
                       '--output-dir', str(child), wall_flag, str(allowance), *extra]
            if inventory is not None:
                command += ['--mass-inventory', str(inventory.resolve())]
            receipt = budget.run(command, ROOT, log, output=directory/(name+'-supervision'))
            entry['supervision'] = receipt
            result_path, manifest_path = child/'result.json', child/'run_manifest.json'
            if receipt['reason']:
                raise RuntimeError(receipt['reason'])
            value = json.loads(result_path.read_text())
            source = json.loads(manifest_path.read_text())
            entry.update(result=str(result_path), result_sha256=sha(result_path),
                         manifest=str(manifest_path), manifest_sha256=sha(manifest_path))
            status = value['status']
            if {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}.get(status) != receipt['returncode'] or status not in STATUSES:
                raise ValueError('I24 child result status contradicts its process exit code')
            entry['scope'] = value.get('scope', 'native PIC first/second-order single-event commissioning')
            if status == 'FAIL':
                # Preserve the genuine native evidence for an evaluated failure.
                # The PIC manifest indexes each rank's emitted CSV, while the
                # independent ensemble has its own per-parent history index.
                native = child/'native_receipt.json' if key=='i24_sodium_provider' else manifest_path
                if key=='i24_pic_commissioning' and not any('-rank-' in f and f.endswith('.csv')
                        for f in source.get('result_files',{})):
                    raise ValueError('evaluated PIC failure has no retained native rank tables')
                failed_receipts.append({'check':label,'native_receipt':str(native),'sha256':sha(native)})
            reason = value.get('reason') or ('all native PIC checks met' if status=='PASS' else
                    'see per-layout checks in '+str(result_path))
        except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError) as error:
            status, reason = 'ERROR', str(error)
            with log.open('a') as stream:
                stream.write('I24 owner infrastructure error: '+reason+'\n')
        rows.append({'check':label, 'status':status, 'runner':str(runner), 'log':str(log), 'reason':reason})
    result.update(formal_status='SKIPPED / NOT QUALIFIED', supporting_checks=rows,
                  supporting_counts={s:sum(row['status']==s for row in rows) for s in STATUSES})
    if failed_receipts:
        (directory/'receipt.txt').write_text(json.dumps(failed_receipts,indent=2)+'\n')
    for status in ('ERROR','FAIL'):
        if result['supporting_counts'][status]:
            result.update(status=status,reason='; '.join(r['check']+' '+r['status']+': '+r['reason']
                          for r in rows if r['status'] in ('FAIL','ERROR')))
            break
    (directory/'stdout.log').write_text(''.join(r['check']+' '+r['status']+': '+r['reason']+'\n'
                                      for r in rows)+'Formal I24 remains SKIPPED.\n')
    (directory/'stderr.log').write_text(''.join(r['check']+': '+r['log']+'\n'
                                      for r in rows if r['status'] in ('FAIL','ERROR')))


def run(test_id,amps,application,output,mpi_ranks=8,mode='regression',*,
        i21_plan=None,i21_plan_sha256=None,wall_seconds=3600,
        i24_lock=None,i24_mass_inventory=None,i24_pic_producer=None,i24_pic_refinement=False):
    """Evaluate explicit prerequisites and produce an honest formal result.

    No native receipt is emitted when execution is not possible. The aggregate
    recognizes manifest execution=not-run and preserves absence rather than
    requiring a fabricated receipt.txt or copying another test's output.
    """
    started=time.monotonic()
    output=Path(output).resolve();directory=output/test_id
    directory.mkdir(parents=True,exist_ok=False)
    stdout=directory/'stdout.log';stderr=directory/'stderr.log'
    stdout.write_text('Formal chemistry execution has not started.\n')
    stderr.write_text('')
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    manifest={'schema':'srcMoon-M3-formal-owner-v1','test_id':test_id,'execution':'not-run',
        'created_utc':datetime.now(timezone.utc).isoformat(),'source_revision':revision,
        'mpi_ranks':mpi_ranks,'openmp_threads':1,'campaign_mode':mode,
        'observational_scoring':False,'random_seed':None,
        'requested_input_file':str(application.resolve()) if application else None,
        'input_snapshot_available':False,'predecessors':[],
        'executable':{'path':str(amps.resolve()),'sha256':sha(amps) if amps.is_file() else None}}
    result={'id':test_id,'name':test_id,'status':'ERROR','scope':'formal-M3-qualification',
        'metrics':{},'criteria_evaluated':False,'blocking_reasons':[]}
    try:
        if not math.isfinite(wall_seconds) or not 0 < wall_seconds <= 3600:
            raise ValueError('one total development allowance in (0, 3600] seconds required')
        if (i21_plan is None) != (i21_plan_sha256 is None) or (i21_plan is not None and test_id!='I21'):
            raise ValueError('explicit plan and SHA-256 must be supplied together for I21 only')
        if i24_lock is not None and (test_id!='I24' or mode!='regression'):
            raise ValueError('the controlled sodium provider case belongs to I24 regression mode')
        if i24_mass_inventory is not None and i24_lock is None:
            raise ValueError('I24 mass inventory requires an explicit I24 sodium release lock')
        if i24_pic_producer is not None and i24_lock is None:
            raise ValueError('I24 PIC producer requires an explicit I24 sodium release lock')
        if i24_pic_refinement and i24_pic_producer is None:
            raise ValueError('I24 PIC refinement requires an explicit I24 PIC producer')
        contract_path=TEST/'integration-tests'/test_id/'reference/acceptance.json'
        contract=json.loads(contract_path.read_text());result['name']=contract['name']
        manifest['acceptance']={'path':str(contract_path),'sha256':sha(contract_path)}
        # Bind the concrete physical/protocol specification even for a run
        # that cannot launch. A consistent preparation result is still SKIP;
        # corrupt identities are ERROR and cannot hide behind missing physics.
        from m3_protocol import validate
        preparation=validate()
        (directory/'protocol_preparation.json').write_text(json.dumps(preparation,indent=2)+'\n')
        manifest['protocol_preparation']=str(directory/'protocol_preparation.json')
        if preparation['status']=='ERROR':
            raise ValueError('M3 protocol preparation error: '+preparation.get('reason','unknown error'))
        # The original frozen preparation protocol is deliberately retained:
        # previous attempts must keep their original interpretation. The user
        # subsequently approved the numerical/statistical choices in a NEW
        # hash-bound addendum. Only those two superseded approval blockers are
        # removed from the EFFECTIVE current readiness assessment; physical,
        # native, source, budget and same-campaign gates remain unchanged.
        from m3_qualification import approved_protocol,ADDENDUM,APPROVED_SHA256
        effective=approved_protocol()
        shutil.copyfile(ADDENDUM,directory/'approved_qualification_addendum.json')
        manifest['approved_qualification_addendum']={'path':str(ADDENDUM),'sha256':APPROVED_SHA256,
            'snapshot':str(directory/'approved_qualification_addendum.json'),
            'scope':'approved protocol only; no physical/native qualification'}
        superseded={'numerical error budget requires review before execution',
            'full statistical family/stream/rare-branch/rate-estimator design requires review'}
        # Preserve the preparation snapshot verbatim. Its original statement
        # about absent executors predates the linked Na photon and capture
        # providers. Describe today's *whole-family* qualification gap while
        # the per-case support entries retain precisely what executed.
        old_execution='physical production providers/executor/scorers not installed; preparatory services are not qualification'
        for reason in preparation['blocking_reasons']:
            if reason not in superseded:
                result['blocking_reasons'].append(
                    'complete physical family, linked execution and numerical qualification remain required'
                    if reason==old_execution else reason)
        result['approved_numerical_bias_limit']=effective['numerical_relative_bias_budgets'][test_id]
        result['statistical_protocol_status']='APPROVED; actual frozen family/cohort/native evidence still required'
        result['blocking_reasons']+=preparation.get('cases',{}).get(test_id,{}).get('missing_physical_inputs',[])
        register=json.loads(CAPABILITIES.read_text())
        if register.get('schema')!='srcMoon-M3-transport-capability-register-v1':
            raise ValueError('unsupported M3 capability-register schema')
        capability=register['capabilities'][test_id]
        manifest['capability_register']={'path':str(CAPABILITIES),'sha256':sha(CAPABILITIES),
            'declaration':capability,'scope':'source capability declaration, not native qualification'}
        if application and application.is_file():
            shutil.copyfile(application,directory/'application.in')
            manifest['input_snapshot_available']=True
            manifest['input_sha256']=sha(directory/'application.in')
        else:
            result['blocking_reasons'].append('required CLI application input is unavailable')
        # Each per-ID contract explicitly lists same-campaign gates according
        # to unchanged formal ordering. Malformed/contradictory data are ERROR;
        # a missing or non-PASS prerequisite is a scientific SKIPPED condition.
        for predecessor in contract['linked_predecessor_gates']:
            path=output/predecessor/'result.json'
            record={'id':predecessor,'path':str(path),'status':'UNAVAILABLE'}
            if path.is_file():
                value=json.loads(path.read_text())
                if value.get('status') not in STATUSES:
                    raise ValueError('invalid predecessor status in '+str(path))
                record.update(status=value['status'],sha256=sha(path))
            manifest['predecessors'].append(record)
            if record['status']!='PASS':
                result['blocking_reasons'].append(predecessor+' has not qualified in this campaign')
        result['blocking_reasons']+=capability['missing_formal_capabilities']
        if not capability['missing_formal_capabilities']:
            # A capability-register edit is not an implementation of the
            # production executor/kinematics and cannot open a PASS path.
            raise ValueError('capability register claims completion without an installed formal execution protocol')
        result.update(status='SKIPPED',reason='; '.join(result['blocking_reasons']))
        if i24_lock is not None:
            i24_supporting(directory,result,manifest,i24_lock,i24_mass_inventory,i24_pic_producer,
                           mpi_ranks,started+wall_seconds,i24_pic_refinement)
        elif i21_plan is None or not manifest['input_snapshot_available']:
            stdout.write_text('SKIPPED '+test_id+': '+result['reason']+'\n'
                'No particle run, lifetime fit, branch score or observation comparison was performed.\n')
        else:
            from m3_qualification_executor import run as run_cohort
            plan=Path(i21_plan).resolve()
            if not digest_syntax(i21_plan_sha256):
                raise ValueError('I21 plan needs its explicit lowercase SHA-256')
            if not plan.is_file():
                raise MissingPrerequisite('I21 cohort plan unavailable: '+str(plan))
            if sha(plan)!=i21_plan_sha256:
                raise ValueError('I21 cohort plan SHA-256 mismatch')
            planned=strict_object(plan)
            if planned.get('test_id')!='I21':
                raise ValueError('the I21 owner requires an I21 cohort plan')
            case_name='trajectory' if planned.get('test_kind')=='trajectory-integrated-hazard' else 'homogeneous'
            mode_name='development' if mode=='regression' else 'qualification-candidate'
            child=directory/case_name
            manifest['i21_cohort']={'mode':mode_name,'plan':str(plan),
                'plan_sha256':i21_plan_sha256,'result':str(child/'result.json')}
            remaining=wall_seconds-(time.monotonic()-started)
            if remaining<=0:raise RuntimeError('I21 owner exhausted its original wall allowance during preflight')
            # The explicit development selection can exercise one numerical
            # matrix under the approved early-M3 exception. The candidate
            # qualification selection still checks U14 and I13-I20 before any
            # launch. Neither selection covers the complete I21 family.
            value=run_cohort(plan,output,child,development=mode=='regression',wall_seconds=remaining,
                application=application,expected_plan_sha256=i21_plan_sha256)
            if sha(plan)!=i21_plan_sha256:
                raise ValueError('I21 cohort plan changed during execution')
            child_status=value['status']
            manifest['i21_cohort'].update(execution=value['execution'],result_sha256=sha(child/'result.json'),
                manifest_sha256=sha(child/'run_manifest.json'))
            result.update(formal_status='SKIPPED / NOT QUALIFIED',
                supporting_checks=[{'check':'I21 '+case_name+' '+mode_name,'status':child_status,
                    'runner':str(TEST/'m3_qualification_executor.py'),'log':str(child/'stdout.log'),
                    'reason':value['reason']}],
                supporting_counts={status:int(child_status==status) for status in STATUSES})
            if child_status in ('FAIL','ERROR'):
                # Preserve the stronger supporting result in the aggregate.
                # The per-ID formal status still identifies an incomplete gate.
                failed=[entry['name'] for entry in value.get('checks',[])
                        if entry['status'] in ('FAIL','ERROR')]
                reason=(case_name+' I21 check failed: '+', '.join(failed) if failed else value['reason'])
                result.update(status=child_status,reason=reason)
                if child_status=='FAIL':
                    registry=child/'native_runs.json'
                    if not registry.is_file():
                        raise ValueError('evaluated I21 cohort failure has no retained native launch registry')
                    # The receipt indexes the actual completed native runs. It
                    # never exists for prerequisite-only or failed launches.
                    shutil.copyfile(registry,directory/'receipt.txt')
            else:
                reason=('one '+case_name+' cohort completed; full I21 coverage remains open'
                    if child_status=='PASS' else value['reason'])
                result['reason']=reason
            stdout.write_text((child/'stdout.log').read_text())
            stderr.write_text((child/'stderr.log').read_text())
    except MissingPrerequisite as error:
        result.update(status='SKIPPED',reason=str(error))
        stdout.write_text('SKIPPED '+test_id+': '+str(error)+'\nNo native execution was performed.\n')
    except (OSError,ValueError,KeyError,TypeError,RuntimeError) as error:
        result.update(status='ERROR',reason=str(error))
        stderr.write_text(str(error)+'\n')
    result['artifacts']={'stdout':str(stdout),'stderr':str(stderr),
        'run_manifest':str(directory/'run_manifest.json')}
    if manifest['input_snapshot_available']:
        result['artifacts']['input_snapshot']=str(directory/'application.in')
    if (directory/'receipt.txt').is_file():result['artifacts']['receipt']=str(directory/'receipt.txt')
    (directory/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (directory/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('test_id',nargs='?',default=Path(sys.argv[0]).resolve().parent.name,
        choices=['I21','I22','I23','I24'])
    parser.add_argument('--amps',type=Path,required=True)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--mpi-ranks',type=int,default=8)
    parser.add_argument('--campaign-mode',choices=('regression','qualification'),default='regression')
    parser.add_argument('--i21-cohort-plan',type=Path)
    parser.add_argument('--i21-cohort-plan-sha256')
    parser.add_argument('--i24-sodium-lock',type=Path,help='run the controlled r14 MPI provider in regression mode')
    parser.add_argument('--i24-pic-producer',type=Path,help='source-bound producer manifest for both native PIC mover orders at 1 and 8 ranks')
    parser.add_argument('--i24-pic-refinement',action='store_true',help='refine native timestep and reference weights at fixed physical horizon and population')
    parser.add_argument('--i24-mass-inventory',type=Path,help='relocated NIST reference mass inventory for that case')
    parser.add_argument('--test-runtime-seconds',type=float,default=3600)
    args=parser.parse_args()
    if args.mpi_ranks<1:parser.error('MPI ranks must be positive')
    value=run(args.test_id,args.amps,args.input,args.output_dir,args.mpi_ranks,args.campaign_mode,
        i21_plan=args.i21_cohort_plan,i21_plan_sha256=args.i21_cohort_plan_sha256,
        wall_seconds=args.test_runtime_seconds,i24_lock=args.i24_sodium_lock,i24_mass_inventory=args.i24_mass_inventory,
        i24_pic_producer=args.i24_pic_producer,i24_pic_refinement=args.i24_pic_refinement)
    print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]


if __name__=='__main__':raise SystemExit(main())
