#!/usr/bin/env python3
"""Run an explicit Stage-14 offline JSON job, separate from schema-5 AMPS input.

This CLI executes the same bounded producers used by the owning verification
registry. It cannot enable an unqualified native mover or convert a synthetic
campaign into observational evidence. Each output binds the entire input job.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from preprocessing.core import load_json, require, freeze_record, canonical, digest, PreprocessingError
from preprocessing.research import (research_configuration,validate_wave_asset,advance_resonant_waves,
    streaming_limit,family_bundle,transfer_campaign)
from preprocessing.research_products import renew,drift_exposure
from preprocessing.impulsive import impulsive_plan
from preprocessing.wind_envelope import gate_wind
from preprocessing.protocols import compare_backgrounds


def execute(job):
    require(job['schema']=='sep-stage14-offline-job-v1','wrong offline-job schema')
    p=job['parameters'];operation=job['operation']
    if operation=='configuration':result=research_configuration(p['capabilities'],p['authorities'])
    elif operation=='wave-step':result=advance_resonant_waves(p['state'],p['controls'],p['dt_s'])
    elif operation=='wave-asset':result=validate_wave_asset(p['asset'],p['active_bindings'])
    elif operation=='streaming-limit':result=streaming_limit(p['rows'],p['registration'])
    elif operation=='reference-family':result=family_bundle(p['members'],p['calibration'],p['target_peclet'])
    elif operation=='renewal':result=renew(p['checkpoint'])
    elif operation=='impulsive-source':result=impulsive_plan(p['asset'],p['physical_number'],p['front_generation'],p['birth_side'])
    elif operation=='wind-envelope':result=gate_wind(p['asset'],p['profiles'])
    elif operation=='event-transfer':result=transfer_campaign(p['protocol'],p['run'],p['observations'],p['predictions'],p['metric_registration'])
    elif operation=='mhd-comparison':result=compare_backgrounds(p['analytic'],p['mhd'],p.get('allow_combined_discrepancy',False))
    elif operation=='drift-report':result=drift_exposure(p['segments'],p['denominators'],p['frame'],p['generation'])
    else:raise PreprocessingError('unknown offline operation: '+str(operation))
    return freeze_record(dict(schema='sep-stage14-offline-result-v1',input_job_sha256=digest(job),operation=operation,
        result=result,production_adapter_qualified=False,observational_campaign_qualified=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        result=execute(load_json(args.job));require(not args.output.exists(),'output already exists')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('xb') as stream:stream.write(canonical(result))
        print('research_result='+str(args.output.resolve()))
        return 0 if result['result'].get('passed',True) else 1
    except (ValueError,KeyError,OSError) as error:
        print('research_error='+str(error));return 2
if __name__=='__main__':raise SystemExit(main())
