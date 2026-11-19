#!/usr/bin/env python3
"""Validate M3 protocol preparation; never execute or promote physical chemistry."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

TEST=Path(__file__).resolve().parent
ROOT=TEST.parents[1]
DATA=Path('/data/vtenishe/moon_validation_data')
IDS=('I21','I22','I23','I24')
# These are the unchanged roadmap numbers, independently checked against the
# protocol's machine-readable declarations. A user-edited looser criterion is
# infrastructure corruption, not permission to alter a physical exit gate.
NUMERIC={'I21':{'ensemble_lifetime_relative_error_max':.02,'inverse_square_rate_relative_error_max':.001},
    'I22':{'ensemble_rate_relative_error_max':.02,'linked_branch_confidence':.95,'local_branch_confidence':.99},
    'I23':{'direct_coefficient_relative_error_max':.02,'ensemble_lifetime_relative_error_max':.03},
    'I24':{'direct_rate_relative_error_max':.01,'ensemble_lifetime_relative_error_max':.03}}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_object(path):
    """Reject duplicate keys/nonfinite numbers instead of silently repairing JSON."""
    def pairs(items):
        value={}
        for key,item in items:
            if key in value:raise ValueError('duplicate JSON key: '+key)
            value[key]=item
        return value
    def nonfinite(value):raise ValueError('nonfinite JSON constant: '+value)
    value=json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=nonfinite)
    if not isinstance(value,dict):raise ValueError('JSON root must be an object: '+str(path))
    return value

def positive_integer(value):return type(value) is int and value>0

def cohort_bound(relative_error,comparisons,risk):
    """Sufficient Hoeffding N for the EXPLICIT half-survival planning design.

    This is not an observed count, exact-binomial minimum, acceptance interval,
    rare-branch design or an authorization to launch a physical case. Round up
    to a multiple of eight while retaining the original-parent trial unit.
    """
    if isinstance(comparisons,bool) or not isinstance(comparisons,int) or comparisons<1:
        raise ValueError('positive integer comparison-family size required')
    if not 0<relative_error<1 or not 0<risk<1:raise ValueError('invalid planning risk/error')
    a=math.log(2.)
    lower=math.exp(-a/(1-relative_error));upper=math.exp(-a/(1+relative_error))
    margin=min(.5-lower,upper-.5)
    count=8*math.ceil(math.log(2*comparisons/risk)/(2*margin*margin)/8)
    if not math.isfinite(count) or count<8:raise ValueError('planning bound is unrepresentable')
    return {'relative_lifetime_error':relative_error,'comparisons':comparisons,
        'family_planning_risk':risk,'global_original_parents':count,
        'survival_probability':.5,'acceptance_probability_interval':[lower,upper],
        'method':'two-sided Hoeffding sufficient bound; PROPOSAL ONLY'}

def validate(protocol_path=TEST/'M3_QUALIFICATION_PROTOCOL.json',comparisons=None):
    """Read-only identity/contract preparation with explicit missing evidence.

    No executable readiness claim is accepted from a boolean/string edit. This
    module has no physical executor and therefore deliberately no formal PASS
    return path. A future physical executor needs its own reviewed contract.
    """
    protocol_path=Path(protocol_path).resolve();base=protocol_path.parent
    result={'schema':'srcMoon-M3-protocol-preparation-v1','status':'ERROR',
        'structural_status':'ERROR','criteria_evaluated':False,'execution':'not-run',
        'formal_M3_status':'SKIPPED / NOT QUALIFIED','blocking_reasons':[],
        'identities':[],'cases':{},'observational_scoring':False}
    def identity(path,expected):
        path=Path(path)
        if not path.is_file():
            result['blocking_reasons'].append('required evidence unavailable: '+str(path));return False
        actual=sha(path);result['identities'].append({'path':str(path),'sha256':actual})
        if actual!=expected:raise ValueError('evidence hash mismatch: '+str(path))
        return True
    try:
        p=read_object(protocol_path);result['protocol_sha256']=sha(protocol_path)
        if p.get('schema')!='srcMoon-M3-qualification-protocol-v1' or p.get('activation')!='preparation-only':
            raise ValueError('unsupported or falsely activated qualification protocol')
        if p.get('formal_executor_installed') is not False or p.get('observational_scoring') is not False:
            raise ValueError('preparation cannot claim a physical executor or observation scoring')
        if set(p['cases'])!=set(IDS):raise ValueError('frozen I21-I24 family changed')
        if any(not positive_integer(p['execution'][name]) for name in (
                'development_wall_seconds','aggregate_memory_bytes','openmp_threads')):
            raise ValueError('resource/layout fields require positive integers, not booleans')
        if p['execution']['development_wall_seconds']!=3600 or p['execution']['aggregate_memory_bytes']!=200000000000:
            raise ValueError('approved development resource policy changed')
        if not isinstance(p['execution']['mpi_ranks'],list) or any(not positive_integer(v) for v in p['execution']['mpi_ranks']) or p['execution']['mpi_ranks']!=[1,8] or p['execution']['openmp_threads']!=1:
            raise ValueError('declared MPI/OpenMP verification layouts changed')
        if p['clock']['mode']!='common-time':raise ValueError('initial qualification clock changed')
        levels=p['numerical_design']['common_steps_to_final_time']
        fractions=p['numerical_design']['common_step_fractions_of_final_time']
        if not isinstance(levels,list) or len(levels)<3 or any(not positive_integer(n) for n in levels) or levels!=sorted(set(levels)) or fractions!=[1/n for n in levels] or p['numerical_design']['same_physical_final_time'] is not True:
            raise ValueError('common-time refinements require at least three declared integer-step levels to one final horizon')
        if p['numerical_design']['numerical_error_budget'] is None:
            result['blocking_reasons'].append('numerical error budget requires review before execution')
        if p['statistical_design']['status']!='PROPOSED_REQUIRES_PREEXECUTION_REVIEW':
            raise ValueError('statistical proposal cannot self-approve by metadata edit')
        risk=p['statistical_design']['family_planning_risk']
        if type(risk) not in (int,float) or not math.isfinite(risk) or not 0<risk<1:
            raise ValueError('planning risk must be a finite probability')
        seeds=p['statistical_design']['seeds']
        if not isinstance(seeds,list) or not seeds or any(not positive_integer(s) for s in seeds) or len(set(seeds))!=len(seeds):
            raise ValueError('proposed seed list must contain distinct positive integers')
        if p['statistical_design']['U17_minimum_branch_events']!=100000:
            raise ValueError('existing U17 event minimum changed')
        result['blocking_reasons'].append('full statistical family/stream/rare-branch/rate-estimator design requires review')
        spec_record=p['physical_model_specification'];spec_path=base/spec_record['path']
        if identity(spec_path,spec_record['sha256']):
            spec=read_object(spec_path)
            if spec.get('schema')!='srcMoon-M3-physical-model-spec-v1' or spec.get('physics_activation')!='disabled-pending-qualified-evidence':
                raise ValueError('unsupported or falsely activated physical-model specification')
            if set(spec['processes'])!=set(IDS):raise ValueError('physical-model process family changed')
            for key in IDS:
                s=spec['processes'][key]
                if s['status']!='BLOCKED_MISSING_DATA' or not s['missing_inputs'] or not s['required_policy_ids']:
                    raise ValueError('physical-model metadata cannot qualify an unimplemented joint product law: '+key)
                result['cases'][key]={'status':'SKIPPED','execution':'not-run',
                    'missing_physical_inputs':s['missing_inputs'],'required_policy_ids':s['required_policy_ids']}
            for record in spec['processed_source_files']:
                # Absolute canonical paths are deliberate; source-tree copies
                # and model outputs are not interchangeable validation inputs.
                identity(record['path'],record['sha256'])
        source_roadmap=ROOT/'srcMoon/development_and_validation.md'
        identity(source_roadmap,p['roadmaps']['source'])
        identity(DATA/'development_and_validation.md',p['roadmaps']['external'])
        for key in IDS:
            case=p['cases'][key]
            contract_path=base/case['acceptance_path']
            if identity(contract_path,case['acceptance_sha256']):
                contract=read_object(contract_path)
                if contract['id']!=key or contract['linked_predecessor_gates']!=case['linked_predecessors']:
                    raise ValueError('frozen ID/predecessor contract mismatch: '+key)
            for name,expected in NUMERIC[key].items():
                actual=case['frozen_criteria'][name]
                if isinstance(actual,bool) or actual!=expected:raise ValueError('frozen acceptance changed: '+key+'/'+name)
            if case['qualification_execution']!='not-run' or case['physical_policy_status']!='BLOCKED_MISSING_DATA':
                raise ValueError('preparation contains a fabricated physical execution: '+key)
        result['blocking_reasons']+=['same-campaign formal I13-I20 and applicable local gates remain required',
            'physical production providers/executor/scorers not installed; preparatory services are not qualification']
        if comparisons is not None:
            result['lifetime_planning']=[cohort_bound(d,comparisons,p['statistical_design']['family_planning_risk']) for d in (.02,.03)]
        result.update(status='SKIPPED',structural_status='PASS')
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as error:
        result.update(status='ERROR',reason=str(error))
    return result

def run(output,protocol_path=TEST/'M3_QUALIFICATION_PROTOCOL.json',comparisons=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    value=validate(protocol_path,comparisons)
    (output/'stdout.log').write_text(json.dumps(value,indent=2)+'\n');(output/'stderr.log').write_text(value.get('reason',''))
    manifest={'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'runner':str(Path(__file__).resolve()),'runner_sha256':sha(__file__),
        'protocol_path':str(Path(protocol_path).resolve()),'execution':'preparation-only; no MPI/particle/scoring run',
        'source_data_identities':value['identities'],'random_seed':None}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    value['counts']={s:int(value['status']==s) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    (output/'result.json').write_text(json.dumps(value,indent=2)+'\n')
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--plan-comparisons',type=int);args=parser.parse_args()
    value=run(args.output_dir,comparisons=args.plan_comparisons);print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]

if __name__=='__main__':raise SystemExit(main())
