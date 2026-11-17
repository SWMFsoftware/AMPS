#!/usr/bin/env python3
"""Native M3 rate-configuration controls, never particle/observational evidence.

The executable compiles the maintained production MoonInput/ElectronImpact
sources, not a copied parser/formula. Inputs are explicit stationary-bath
fixtures derived from the complete bounded regression receipt. The published
native coefficients independently supply expected rates. Missing physical
dependencies are SKIPPED; compile/run failures ERROR; evaluated misses FAIL.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
NATIVE=Path('/data/vtenishe/moon_validation_data/photochemistry/raw/electron_impact/voronov/cfit.dat')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def change(text,key,value):
    """Edit an analytical input fixture; native code interprets the grammar."""
    pattern=r'(?m)^\s*'+re.escape(key)+r'\s*=[^\n]*$'
    if re.search(pattern,text):return re.sub(pattern,lambda _:key+' = '+str(value),text)
    return text.replace('#end application',key+' = '+str(value)+'\n#end application')


def run_checks(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    result={'schema':'srcMoon-M3-native-rate-controls-v1','scope':'native configuration/rate service only; no particles or observations',
        'status':'ERROR','checks':[],'formal_I23_status':'SKIPPED'}
    compiler=shutil.which('g++')
    dependencies=[NATIVE,ROOT/'srcMoon/examples/m1_regression.in']
    missing=[str(p) for p in dependencies if not p.is_file()]
    if not compiler or missing:
        result.update(status='SKIPPED',reason='required native compiler/input unavailable',missing=missing)
    else:
        executable=output/'runtime-rate-probe'
        command=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-pedantic',
                 '-I',str(ROOT/'srcMoon'),str(ROOT/'srcMoon/test/m3_runtime_rates_probe.cpp'),
                 str(ROOT/'srcMoon/MoonInput.cpp'),str(ROOT/'srcMoon/species/ElectronImpact.cpp'),
                 '-o',str(executable)]
        process=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
        (output/'compile.stdout.log').write_text(process.stdout)
        (output/'compile.stderr.log').write_text(process.stderr)
        result['compile']={'command':command,'returncode':process.returncode}
        if process.returncode:
            result['reason']='native production configuration adapter did not compile'
        else:
            template=(ROOT/'srcMoon/examples/m1_regression.in').read_text()
            # Fixture-only: these paths are not read in sphere mode, but their
            # complete receipt must remain explicit and independently parsed.
            base=template
            for key in ('surface_cea_file','surface_tecplot_file'):
                base=change(base,key,str(output/(key+'.unused')))
            enabled=change(base,'electron_impact_rate_model','maxwellian-voronov')
            enabled=change(enabled,'electron_impact_coefficients_file',str(NATIVE))
            enabled=change(enabled,'electron_density_m3','2500000')
            enabled=change(enabled,'electron_temperature_ev','10')
            corrupt=output/'malformed-cfit-not-authoritative.dat'
            corrupt.write_text('parser rejection fixture only; NOT scientific data\n')
            variants=[('legacy-unchanged',base,0,None),
                ('explicit-maxwellian',enabled,0,None),
                ('zero-density',change(enabled,'electron_density_m3',0),0,None),
                ('domain-minimum',change(enabled,'electron_temperature_ev',1),0,None),
                ('domain-maximum',change(enabled,'electron_temperature_ev',20000),0,None),
                ('below-domain',change(enabled,'electron_temperature_ev',0.999),4,None),
                ('above-domain',change(enabled,'electron_temperature_ev',20001),4,None),
                ('negative-density',change(enabled,'electron_density_m3',-1),3,None),
                ('nan-density',change(enabled,'electron_density_m3','nan'),3,None),
                ('infinite-temperature',change(enabled,'electron_temperature_ev','inf'),3,None),
                ('zero-temperature',change(enabled,'electron_temperature_ev',0),3,None),
                ('missing-density',re.sub(r'(?m)^electron_density_m3[^\n]*\n','',enabled),3,None),
                ('missing-native-table',change(enabled,'electron_impact_coefficients_file',str(output/'missing.dat')),3,None),
                ('native-table-corrupt',change(enabled,'electron_impact_coefficients_file',str(corrupt)),4,None),
                ('inactive-settings',change(enabled,'electron_impact_rate_model','legacy'),3,None),
                ('unknown-closure',change(enabled,'electron_impact_rate_model','kappa-assumed'),3,None),
                ('transactional-install',enabled,0,'rollback')]
            for name,text,expected,mode in variants:
                directory=output/name;directory.mkdir(exist_ok=False)
                application=directory/'application.in';application.write_text(text)
                invocation=[str(executable),str(application)]+([mode] if mode else [])
                p=subprocess.run(invocation,cwd=directory,text=True,capture_output=True)
                (directory/'stdout.log').write_text(p.stdout);(directory/'stderr.log').write_text(p.stderr)
                status='PASS' if p.returncode==expected else 'ERROR' if p.returncode<0 else 'FAIL'
                values=dict(line.split('=',1) for line in p.stdout.splitlines() if '=' in line)
                if status=='PASS':
                    if expected==3 and values.get('parse_preserved')!='1':status='FAIL'
                    if expected==4 and values.get('install_preserved')!='1':status='FAIL'
                    if expected==0:
                        active=name!='legacy-unchanged'
                        if values.get('has_table')!=str(int(active)):status='FAIL'
                        if active:
                            density=float(values['density_m3']);temperature=float(values['temperature_ev'])
                            # Independent oracle for wiring/unit checks. The
                            # dedicated U16 scorer uses high precision and
                            # native Fortran as the stronger formula reference.
                            coefficients={'He':(24.6,0,1.75e-14,.18,.35),
                                'Na':(5.1,1,1.01e-13,.275,.23),
                                'Ne':(21.6,1,1.5e-14,.0329,.43),
                                'Ar':(15.8,1,5.99e-14,.136,.26)}
                            for atom,(energy,power,a,x,k) in coefficients.items():
                                u=energy/temperature
                                reference=density*a*(1+power*math.sqrt(u))/(x+u)*u**k*math.exp(-u)
                                actual=float(values[atom+'_rate_s'])
                                if not math.isclose(actual,reference,rel_tol=1e-12,abs_tol=0):status='FAIL'
                        if mode and (values.get('failed_install_preserved')!='1' or
                                     values.get('repeat_install_rejected')!='1'):status='FAIL'
                result['checks'].append({'name':name,'status':status,'command':invocation,
                    'expected_returncode':expected,'returncode':p.returncode,
                    'input_sha256':sha(application),'log':str(directory/'stdout.log')})
            statuses=[row['status'] for row in result['checks']]
            result['status']='ERROR' if 'ERROR' in statuses else 'FAIL' if 'FAIL' in statuses else 'PASS'
            result['executable']={'path':str(executable),'sha256':sha(executable)}
    result['counts']={s:sum(row['status']==s for row in result['checks']) for s in ('PASS','FAIL','SKIPPED','ERROR')}
    result['wall_seconds']=time.monotonic()-started
    # Preserve authoritative production/adapter identity even for a compile
    # failure. These local controls are not linked particle simulations.
    result['source_files']={str(p.relative_to(ROOT)):sha(p) for p in (
        ROOT/'srcMoon/MoonInput.cpp',ROOT/'srcMoon/MoonInput.h',
        ROOT/'srcMoon/species/ElectronImpact.cpp',ROOT/'srcMoon/species/ElectronImpact.h',
        ROOT/'srcMoon/test/m3_runtime_rates_probe.cpp',Path(__file__).resolve())}
    result['data_files']={str(p):sha(p) for p in dependencies if p.is_file()}
    result['source_revision']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error('output must be new or empty')
    result=run_checks(args.output_dir)
    print(json.dumps({'status':result['status'],'counts':result['counts'],
                      'result':str(args.output_dir/'result.json')},indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[result['status']])
