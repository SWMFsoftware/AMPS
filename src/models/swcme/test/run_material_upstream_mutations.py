#!/usr/bin/env python3
"""Prove the D10-1 numerical gates reject independent physical regressions.

Each variant is an evidence-owned COPY of the actual canonical header. The
production worktree, numerical tolerances, tests and reference values stay
unchanged. A control passes only when compilation succeeds and the intended
native test returns FAIL/exit1; a compiler error, crash or missing metric is
not evidence that a physical gate works.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parent
MODEL=ROOT.parent
# Exact, syntax-preserving defects chosen from the governing invariants, not
# selected after observing a candidate's residual. They exercise signed fields,
# independent gradients/force/work, real shocks, atomicity and fixed duration.
MUTATIONS=(
    ('mass-r-cubed','MAT1D02','const double density_ratio=ratio*ratio;',
     'const double density_ratio=ratio*ratio*ratio;'),
    ('isothermal-not-adiabatic','MAT1D06','const double thermal_ratio=std::pow(density_ratio,config_.gamma-1.0);',
     'const double thermal_ratio=1.0;'),
    ('reversed-parker-winding','MAT1D01','const double k=config_.rotation_rad_s/config_.wind_m_s;',
     'const double k=-config_.rotation_rad_s/config_.wind_m_s;'),
    ('transposed-magnetic-gradient','MAT1D05','value.magnetic_gradient_t_m_inv[3*i+j]=',
     'value.magnetic_gradient_t_m_inv[3*j+i]='),
    ('omitted-transverse-lorentz','MAT1D07','+gravity_density*er[i]-lorentz[i];',
     '+gravity_density*er[i];'),
    ('discarded-mechanical-work','MAT1D06','value.external_work_w_m3=detail::dot(value.velocity_m_s,value.external_force_n_m3);',
     'value.external_work_w_m3=0.0;'),
    ('rounded-endpoint-duration','MAT1D08','return profile_.move(position,fixed_dt_s,out);',
     'return profile_.move(position,end-begin,out);'),
    ('false-ambient-shock-success','MAT1D09','if(!candidate.jump.solver_converged)',
     'candidate.jump=shock::JumpResult{};candidate.jump.upstream=upstream;candidate.jump.downstream=upstream;\n    if(!candidate.jump.solver_converged)'),
    ('partial-sample-write','MAT1D10','ModelStatus sample(const Vector& position,Sample& out) const {',
     'ModelStatus sample(const Vector& position,Sample& out) const {\n    out.position_m=position;'),
)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--cxx',default=os.environ.get('CXX','g++'))
    args=parser.parse_args()
    if args.output.exists() or args.jobs<1:parser.error('fresh output and positive job count required')
    output=args.output.resolve();output.mkdir(parents=True)
    header=MODEL/'swcme_material_upstream.hpp';original=header.read_text();source_hash=sha(header)
    owned=[ROOT/'core/test_material_upstream.cpp',ROOT/'core/material_upstream_oracle.hpp']
    hashes={str(p):sha(p) for p in owned}
    def run(spec):
        name,identity,old,new=spec;directory=output/name;directory.mkdir()
        if original.count(old)!=1:
            return {'id':name,'status':'ERROR','reason':'source mutation marker is absent/ambiguous; no regression was tested'}
        target=directory/header.name;target.write_text(original.replace(old,new,1))
        binary=directory/'test';compile_command=[args.cxx,'-std=c++17','-O2','-Wall','-Wextra','-Wpedantic',
            '-I'+str(directory),'-isystem',str(MODEL),'-I'+str(ROOT/'core'),str(owned[0]),str(MODEL/'swcme.a'),'-o',str(binary)]
        started=time.monotonic()
        with (directory/'build.log').open('w') as stream:
            built=subprocess.run(compile_command,stdout=stream,stderr=subprocess.STDOUT)
        command=[str(binary),'--test',identity];code=None;stdout=''
        if built.returncode==0:
            with (directory/'native.log').open('w') as stream:
                tested=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
            code=tested.returncode;stdout=(directory/'native.log').read_text()
        rejected=built.returncode==0 and code==1 and f'[{identity}] FAIL' in stdout
        result={'id':name,'native_id':identity,'status':'PASS' if rejected else 'FAIL',
                'compile_argv':compile_command,'actual_compile_exit':built.returncode,
                'argv':command,'actual_native_exit':code,'duration_s':time.monotonic()-started,
                'mutated_header_sha256':sha(target),'scope':'validly compiled physical regression must fail its independent native gate'}
        (directory/'process.json').write_text(json.dumps(result,indent=2)+'\n')
        print(name,result['status'],'compile',built.returncode,'native',code,flush=True)
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        records=list(executor.map(run,MUTATIONS))
    unchanged=sha(header)==source_hash and all(sha(Path(p))==value for p,value in hashes.items())
    report={'schema':'swcme-material-negative-controls-v1','production_header_sha256':source_hash,
            'test_source_sha256':hashes,'worktree_sources_unchanged':unchanged,'results':records,
            'status':'PASS' if unchanged and all(r['status']=='PASS' for r in records) else 'FAIL'}
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0 if report['status']=='PASS' else 1

if __name__=='__main__':raise SystemExit(main())
