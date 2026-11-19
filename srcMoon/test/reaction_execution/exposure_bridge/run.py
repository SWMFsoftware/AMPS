#!/usr/bin/env python3
"""Compile and verify the production exposure bridge; analytical controls only."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
# Keep an explicit source dependency inventory. The whole generated application
# tree is deliberately untouched by this stand-alone production-kernel runner.
FILES=['srcMoon/species/'+name for name in (
    'ReactionExecution.h','ReactionExecution.cpp','PiecewiseHazard.h',
    'PiecewiseHazard.cpp','ReactionExposure.h','ReactionExposure.cpp')]+[
    'src/models/chemical_solver/ReactionNetwork.h',
    str(HERE.relative_to(ROOT)/'probe.cpp'),str(HERE.relative_to(ROOT)/'run.py'),
    str(HERE.relative_to(ROOT)/'README.md'),'srcMoon/test/M3_INDEPENDENT_SOFTWARE_CONTRACT.md']
EXPECTED={
    'origin_clock','pure_first_crossing','quadratic_analytic_event_time',
    'preview_preserves_origin','instantaneous_not_integrated_channel_CDF',
    'accepted_copy_consumed_only','immutable_retry_same_event',
    'no_event_ignores_channel_draw','analytic_integrated_debit',
    'partition_preserves_absolute_crossing','gap_refusal_transactional',
    'support_preflight_independent_of_early_draw','invalid_event_draw_transactional',
    'endpoint_uses_declared_left_limit','zero_rate_endpoint_preview_explicit',
    'zero_rate_endpoint_refuses_CDF_no_time_shift','zero_exposure_keeps_pending_depth',
    'duplicate_channel_refused','consumed_origin_refused','sub_ULP_exposure_compensation',
    'shared_selector_zero_total_transactional','shared_selector_constant_control',
    'precision_counterexample_clock_constructed','compensated_pending_depth_no_early_event',
    'public_integrated_clock_refuses_lossy_pending_depth'}

def parse_receipt(transcript,returncode):
    """A successful empty/incomplete transcript is ERROR, never local PASS."""
    rows=[{'status':status,'name':name} for status,name in re.findall(r'^(PASS|FAIL) (.+)$',transcript,re.M)]
    summary=re.findall(r'^SUMMARY checks=(\d+) failed=(\d+)$',transcript,re.M)
    if len(summary)!=1 or not rows or len({r['name'] for r in rows})!=len(rows) or {r['name'] for r in rows}!=EXPECTED:
        raise ValueError('missing/duplicate/incomplete probe receipt or required control')
    count,failed=map(int,summary[0])
    if count!=len(rows) or failed!=sum(r['status']=='FAIL' for r in rows) or returncode!=int(failed>0):
        raise ValueError('probe status/exit/count contradiction')
    return rows,count,failed

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(output,compiler='g++'):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    result={'schema':'srcMoon-reaction-exposure-bridge-v1','status':'ERROR',
        'scope':'local production numerical controls; no lunar kinetic qualification',
        'formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','counts':{}}
    manifest={'external_data':[],
        'random_seed':None,'stochastic_controls':'explicit uniforms, no RNG',
        'MPI':'none','time':'relative physical seconds, no absolute epoch'}
    try:
        manifest.update(source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            source_hashes={name:sha(ROOT/name) for name in FILES})
        resolved=shutil.which(compiler)
        if not resolved:
            result.update(status='SKIPPED',reason='C++17 compiler unavailable');return result
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            result.update(status='SKIPPED',reason='matplotlib required for declared diagnostics');return result
        version=subprocess.check_output([resolved,'--version'],text=True,timeout=15)
        (output/'compiler-version.log').write_text(version);manifest['compiler']=version
        executable=output/'exposure-bridge-probe'
        command=[resolved,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT)]+[
            str(ROOT/'srcMoon/species'/name) for name in
            ('ReactionExecution.cpp','PiecewiseHazard.cpp','ReactionExposure.cpp')]+[
            str(HERE/'probe.cpp'),'-o',str(executable)]
        manifest['compile_command']=command
        with (output/'compile.log').open('w') as stream:
            code=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=120).returncode
        if code:raise RuntimeError('production bridge compilation failed; see compile.log')
        manifest['executable']={'path':str(executable),'sha256':sha(executable)}
        proc=subprocess.run([str(executable)],capture_output=True,text=True,timeout=60)
        (output/'stdout.log').write_text(proc.stdout);(output/'stderr.log').write_text(proc.stderr)
        rows,count,failed=parse_receipt(proc.stdout,proc.returncode)
        metrics=re.findall(r'^METRIC actual_event_time=(\S+) actual_instantaneous_channel_7=(\S+) reference_event_time=(\S+)$',proc.stdout,re.M)
        if len(metrics)!=1:raise ValueError('missing native event diagnostic')
        actual_time,actual_probability,expected_time=map(float,metrics[0])
        if not all(math.isfinite(v) for v in (actual_time,actual_probability,expected_time)):
            raise ValueError('nonfinite native event diagnostic')
        if {name:sha(ROOT/name) for name in FILES}!=manifest['source_hashes']:
            raise RuntimeError('production or test source changed during run')
        # Independent closed-form diagnostics over a predeclared mathematical
        # depth grid; these are not model-observation residuals or solar data.
        table=[]
        for i in range(1,101):
            depth=i/100;time=2*depth/(1+math.sqrt(1+4*depth))
            table.append((depth,time,1/(1+2*time),time/depth))
        with (output/'analytic_reference.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['optical_depth','analytic_event_s',
                'instantaneous_channel_7_probability','integrated_fraction_NOT_event_probability']);writer.writerows(table)
        fig,ax=plt.subplots();ax.plot([r[0] for r in table],[r[2] for r in table],label='Instantaneous event probability')
        ax.plot([r[0] for r in table],[r[3] for r in table],label='Integrated fraction (incorrect event CDF)')
        ax.plot([1],[actual_probability],'o',label='Compiled production event')
        ax.set(xlabel='Optical depth',ylabel='Channel 7 fraction',title='Analytical control only: rates 1 and 2t');ax.legend()
        for suffix in ('png','eps'):fig.savefig(output/('analytic_reference.'+suffix),bbox_inches='tight')
        plt.close(fig)
        with (output/'production_comparison.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['quantity','production','independent_reference','units'])
            writer.writerow(['first_event_time',actual_time,expected_time,'s'])
            writer.writerow(['instantaneous_channel_7_probability',actual_probability,1/(1+2*expected_time),'1'])
        result.update(status='FAIL' if failed else 'PASS',results=rows,
            counts={'PASS':count-failed,'FAIL':failed,'ERROR':0,'SKIPPED':0},
            artifacts={name:str(output/name) for name in ('stdout.log','stderr.log','compile.log',
                'analytic_reference.csv','analytic_reference.png','analytic_reference.eps','production_comparison.csv')})
    except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as error:
        result.update(status='ERROR',reason=str(error))
    finally:
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--cxx',default='g++');args=parser.parse_args()
    value=run(args.output_dir,args.cxx);print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]

if __name__=='__main__':raise SystemExit(main())
