#!/usr/bin/env python3
"""Freeze deterministic physical normalization BEFORE a cluster residual score.

This executes a thin linked production source-rate integral, not an exosphere
simulation or observation. No atom supply is inferred from random accepted
particles or fitted to model agreement. A source-free cohort scale is the
explicit initial physical stock, independent of refined packet weights.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys

from campaign_budget import CampaignBudget
from integration_common import (REPO,manifest,parse_receipt,sha256,moon_setting,
                                replace_moon_setting,MissingPrerequisite)
from build_m2_qualification import verify_retained_executable
from m2_protocol import write_case
from m2_initial_state import physics_identity_sha256,assignments


def _candidate_initial_state(application,identity,receipt):
    """Retain an UNQUALIFIED actual native state, never manufacture review.

    The ordered full boundary/layer CSV and actual furnished CSPICE/DEM
    identities use the same validators as formal evidence. This candidate
    merely describes what ran; method qualification and investigator review
    remain separate mandatory inputs, even for a uniform-temperature state.
    """
    from m2_initial_state_evidence import read_state,physical_external_data
    required=('initial_thermal_state_csv','initial_thermal_state_rows',
              'initial_thermal_state_element_count','initial_thermal_state_layer_count',
              'thermal_initialization_mode','thermal_initialization_et_s')
    if any(key not in receipt for key in required):
        raise ValueError('requested native I19 initial thermal-state export/receipt unavailable')
    state=Path(str(receipt['initial_thermal_state_csv']))
    if not state.is_absolute() or not state.is_file():
        raise ValueError('requested actual native initial thermal-state CSV is unavailable')
    shape=read_state(state,receipt['initial_thermal_state_element_count'],
                     receipt['initial_thermal_state_layer_count'])
    if (float(receipt['initial_thermal_state_rows'])!=shape['rows'] or
        float(moon_setting(application,'thermal_layer_count'))!=shape['layer_count']):
        raise ValueError('native initial thermal-state count/configuration contradicts its complete CSV')
    mode=str(receipt['thermal_initialization_mode'])
    selected=assignments(application).get('thermal_initialization_mode',{}).get('value','declared-uniform').lower()
    epoch=float(receipt['thermal_initialization_et_s'])
    if mode not in ('declared-uniform','periodic-spinup') or mode!=selected or not math.isfinite(epoch):
        raise ValueError('actual initial thermal method/epoch differs from the requested native input')
    if mode=='declared-uniform':
        temperature=float(moon_setting(application,'thermal_initial_temperature_k'))
        if shape['minimum_kinetic_temperature_K']!=temperature or shape['maximum_kinetic_temperature_K']!=temperature:
            raise ValueError('actual native uniform initial state differs from declared temperature')
    return {'schema':'srcMoon-M2-unqualified-initial-state-candidate-v1',
        'qualified':False,'scope':'actual model initial condition; not a reviewed qualification or observation',
        'native_state':{'path':str(state.resolve()),'sha256':sha256(state)},
        'state_shape':shape,'initialization_mode':mode,'ephemeris_time_tdb_s':epoch,
        'simulation_start_utc':moon_setting(application,'simulation_start_utc',preserve_case=True),
        'physics_identity_sha256':physics_identity_sha256(application),
        'physical_external_data':physical_external_data(identity,application,receipt)}


def prepare(executable,template,output,*,elements=None,sites_sha256=None,horizon_s=2550240.0,mpi_ranks=8,wall_s=3600,tier='development',record_initial_state=False):
    """Generate model-input normalization with explicit physical staging.

    Development preserves the existing bounded I18 mesh/Ar fixture. Formal
    cluster staging preserves every reviewed AMR/source/thermal setting and
    refuses unavailable physical prerequisites rather than enabling them.
    Both perform the same positive-H source quadrature. No transport occurs;
    any explicitly selected production periodic initialization occurs BEFORE
    the native integral and is not bypassed for a cheaper reference weight.
    """
    executable,template,output=Path(executable).resolve(),Path(template).resolve(),Path(output).resolve()
    if tier not in ('development','cluster'):
        raise ValueError('normalization tier must be development or cluster')
    output.mkdir(parents=True,exist_ok=True);(output/'runtime').mkdir()
    # This owner launches ONLY the retained binary. Its preserved generated
    # headers bind that binary; an unrelated canonical tree restored by the
    # global runner must not prevent this focused post-restore execution.
    verify_retained_executable(executable,require_current_generated_tree=False)
    case={'case':'normalization','seed':731993,'dt_divisor':1,'weight_divisor':1}
    application=output/'application.in'
    try:
        write_case(template,application,case,tier,elements=elements,horizon_s=horizon_s)
    except MissingPrerequisite as error:
        result={'status':'SKIPPED','scope':'model-input-reference','execution':'not-run',
                'production_input_tier':tier,'reason':str(error)}
        (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        return result
    if tier=='development':
        # The historical bounded source-only normalizer already integrates
        # the full positive H. Its native profile must therefore remain
        # 'cluster' (the parser reserves development for short transport),
        # while production_input_tier truthfully records bounded staging.
        text=application.read_text()
        text=replace_moon_setting(text,'m2_verification_profile','cluster',application)
        text=replace_moon_setting(text,'m2_verification_horizon_s',format(horizon_s,'.17g'),application)
        application.write_text(text)
    if record_initial_state:
        # Diagnostic request only: retain the actual applied native columns
        # without choosing an initial temperature, epoch, method or criterion.
        text=application.read_text()
        key='m2_verification_record_initial_thermal_state'
        if key in assignments(application):
            text=replace_moon_setting(text,key,'true',application)
        else:
            text=text.replace('#end application',key+' = true\n#end application')
        application.write_text(text)
    physics_hash=physics_identity_sha256(application)
    template_physics_hash=physics_identity_sha256(template)
    if tier=='cluster' and physics_hash!=template_physics_hash:
        raise ValueError('formal normalization changed reviewed physical configuration')
    # Count the actual derived receipt, including intentionally shared logical
    # supports, rather than assuming the default fixture or requested mapping
    # survived input transformation. One initial physical stock exists per tag.
    site_count=len(moon_setting(application,'m2_verification_surface_elements').split(','))
    command=['mpiexec','-n',str(mpi_ranks),str(executable),'--moon-integration-test','I19',
             '--m2-budget-check','--integration-output',str(output/'receipt.txt'),'-input',str(application)]
    identity=manifest('I19',executable,command,731993,application,mpi_ranks)
    identity.update(execution_scope='deterministic-model-input-source-integral-only',
                    production_input_tier=tier,physics_identity_sha256=physics_hash,
                    formal_test_status='SKIPPED',formal_reason='no linked transport/observational qualification')
    (output/'run_manifest.json').write_text(json.dumps(identity,indent=2,sort_keys=True)+'\n')
    budget=CampaignBudget(wall_s,min(60,wall_s*.1))
    execution=budget.run(command,output/'runtime',output/'stdout.log',output=output)
    (output/'stderr.log').write_text('Combined MPI stdout/stderr is retained in stdout.log.\n')
    if execution['returncode']!=0:
        result={'status':'ERROR','scope':'model-input-reference','reason':execution['reason'] or 'native integral failed'}
    else:
        r=parse_receipt(output/'receipt.txt');duration=float(r['duration_s'])
        expected_ar=float(moon_setting(application,'argon_source_rate_atoms_s'))*duration
        packets=float(moon_setting(application,'m2_verification_cohort_packets'))
        initial={'Na':float(r['initial_na_cohort_atoms_per_site']),
                 'Ar':float(r['initial_ar_cohort_atoms_per_site'])}
        weights={species:value/packets for species,value in initial.items()} if packets>0 else {}
        external={'Na':float(r['expected_external_na_atoms']),
                  'Ar':float(r['expected_external_ar_atoms'])}
        passed=(math.isfinite(duration) and duration>0 and
                math.isclose(duration,horizon_s,rel_tol=1e-14) and
                math.isclose(float(r['expected_external_ar_atoms']),expected_ar,rel_tol=1e-12) and
                math.isfinite(packets) and packets>0 and
                all(math.isfinite(value) and value>0 for value in initial.values()) and
                all(math.isfinite(value) and value>0 for value in weights.values()) and
                all(math.isfinite(value) and value>=0 for value in external.values()))
        result={'status':'READY' if passed else 'ERROR','freeze_before_residuals':True,
            'source_configuration_sha256':sha256(template),
            'production_input_tier':tier,'physics_identity_sha256':physics_hash,
            'template_physics_identity_sha256':template_physics_hash,
            'actual_input':{'path':str(application),'sha256':sha256(application)},
            'native_manifest':{'path':str(output/'run_manifest.json'),'sha256':sha256(output/'run_manifest.json')},
            'native_receipt':{'path':str(output/'receipt.txt'),'sha256':sha256(output/'receipt.txt')},
            'sites_sha256':sites_sha256,'protocol_sha256':sha256(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')),
            'executable_sha256':sha256(executable),'duration_s':duration,
            'source-fed':external,
            'cohort':{'Na':site_count*float(r['initial_na_cohort_atoms_per_site']),
                      'Ar':site_count*float(r['initial_ar_cohort_atoms_per_site'])},
            'reference_method':'converged Simpson integral of production MIV+Ar rates with SPICE heliocentric radius; independent of stochastic transfers',
            'native_weight_provenance':{
                'method':'native ordinary-initialization atoms per site divided by declared cohort packets, before numerical dt/weight overrides',
                'cohort_packets':packets,'initial_atoms_per_site':initial,
                'reference_weight_atoms':weights,
                'actual_input_sha256':sha256(application),
                'native_receipt_sha256':sha256(output/'receipt.txt')},
            'producer':{'path':str(Path(__file__)),'sha256':sha256(Path(__file__))},
            'native_receipt_sha256':sha256(output/'receipt.txt'),'execution':execution,
            'formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED'}
        # Cluster staging already requests this capability. Development can
        # explicitly opt in without turning its model-input READY result into
        # initial-state qualification. Missing/corrupt requested output is an
        # infrastructure ERROR, never a silently omitted candidate.
        requested=assignments(application).get('m2_verification_record_initial_thermal_state',{}).get('value','false').lower()
        if requested=='true':
            try:
                result['candidate_initial_state']=_candidate_initial_state(application,identity,r)
            except (OSError,ValueError,KeyError,RuntimeError) as error:
                result.update(status='ERROR',reason='requested actual initial-state evidence invalid: '+str(error))
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--amps',required=True,type=Path);p.add_argument('--input',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--sites-file',type=Path)
    p.add_argument('--horizon-s',type=float,default=2550240);p.add_argument('--mpi-ranks',type=int,default=8)
    p.add_argument('--tier',choices=('development','cluster'),default='development',
                   help='development keeps bounded fixture; cluster preserves reviewed physical template for formal references')
    p.add_argument('--record-initial-state',action='store_true',
                   help='retain actual native I19 initial columns as an explicitly UNQUALIFIED candidate')
    a=p.parse_args()
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error('output must be new or empty')
    if not math.isfinite(a.horizon_s) or a.horizon_s<=0:p.error('positive finite SI horizon required')
    elements=None;site_hash=None
    if a.sites_file:
        registry=json.loads(a.sites_file.read_text())
        if registry.get('status')!='READY' or not registry.get('selection_frozen'):p.error('site registry must be READY and frozen')
        elements=tuple(r['logical_surface_element'] for r in registry['sites']);site_hash=sha256(a.sites_file)
    result=prepare(a.amps.resolve(),a.input.resolve(),a.output_dir.resolve(),elements=elements,
                   sites_sha256=site_hash,horizon_s=a.horizon_s,mpi_ranks=a.mpi_ranks,tier=a.tier,
                   record_initial_state=a.record_initial_state)
    print(json.dumps({'status':result['status'],'normalization':str(a.output_dir/'result.json')},indent=2))
    return 0 if result['status']=='READY' else 77 if result['status']=='SKIPPED' else 2

if __name__=='__main__':raise SystemExit(main())
