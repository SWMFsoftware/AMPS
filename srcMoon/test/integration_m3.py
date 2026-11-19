#!/usr/bin/env python3
"""Own the existing I21--I24 formal gates and retain their actual scope.

The rate and transactional supporting services are executable; full physical
daughter policies/combined event dispatch remain unavailable. This owner
preserves the individual contract, nearest predecessor outcomes, actual input
and capability evidence. It never manufactures a survival run or interprets
a parser/rate/native-storage PASS as a qualified particle-chemistry result.
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
from m3_qualification import MissingPrerequisite

ROOT=Path(__file__).resolve().parents[2]
TEST=Path(__file__).resolve().parent
CAPABILITIES=TEST/'M3_TRANSPORT_CAPABILITIES.json'
STATUSES=('PASS','FAIL','ERROR','SKIPPED')


def sha(path):
    """Byte identity only; matching bytes are not scientific readiness."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(test_id,amps,application,output,mpi_ranks=8,mode='regression'):
    """Evaluate explicit prerequisites and produce an honest formal result.

    No native receipt is emitted when execution is not possible. The aggregate
    recognizes manifest execution=not-run and preserves absence rather than
    requiring a fabricated receipt.txt or copying another test's output.
    """
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
        result['blocking_reasons'] += [reason for reason in preparation['blocking_reasons'] if reason not in superseded]
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
        stdout.write_text('SKIPPED '+test_id+': '+result['reason']+'\n'
            'No particle run, lifetime fit, branch score or observation comparison was performed.\n')
    except MissingPrerequisite as error:
        result.update(status='SKIPPED',reason=str(error))
        stdout.write_text('SKIPPED '+test_id+': '+str(error)+'\nNo native execution was performed.\n')
    except (OSError,ValueError,KeyError,TypeError) as error:
        result.update(status='ERROR',reason=str(error))
        stderr.write_text(str(error)+'\n')
    result['artifacts']={'stdout':str(stdout),'stderr':str(stderr),
        'run_manifest':str(directory/'run_manifest.json')}
    if manifest['input_snapshot_available']:
        result['artifacts']['input_snapshot']=str(directory/'application.in')
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
    args=parser.parse_args()
    if args.mpi_ranks<1:parser.error('MPI ranks must be positive')
    value=run(args.test_id,args.amps,args.input,args.output_dir,args.mpi_ranks,args.campaign_mode)
    print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]


if __name__=='__main__':raise SystemExit(main())
