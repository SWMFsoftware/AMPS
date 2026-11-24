#!/usr/bin/env python3
"""U16 standalone production-kernel checks; no observations or AMPS build.

Only the production ElectronImpact.cpp implements coefficients. This module
independently evaluates native coefficients with70-digit Decimal arithmetic,
cross-checks native DATA/CSV identities, and preserves every command/log/hash.
Temporary altered-byte tables are parser/arithmetic fixtures, never D05 data.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import tempfile
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'srcMoon/test'))
from chemical_database_build import shared_sources, shared_dependencies
DATA_ROOT=Path('/data/vtenishe/moon_validation_data')
NAMES={'He':2,'Na':11,'Ne':10,'Ar':18}
RELATIVE_TOLERANCE=1.0e-12

class MissingPrerequisite(RuntimeError):
    pass

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def document(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def preflight(data_root):
    package=data_root/'photochemistry'
    required=[package/x for x in ('README.md','provenance.json','qa_report.json',
                                  'processed/electron_impact/processing_summary.json')]
    if any(not p.is_file() for p in required):
        raise MissingPrerequisite('D05 package README/provenance/QA/processing summary missing')
    qa=document(package/'qa_report.json')
    if qa.get('status')!='PASS':
        raise MissingPrerequisite('D05 package QA is not PASS')
    if qa.get('record_counts',{}).get('voronov_coefficients')!=406:
        raise ValueError('D05 QA record count is inconsistent')
    provenance=document(package/'provenance.json')
    summary=document(required[-1])
    frozen=document(HERE/'reference/source_manifest.json')
    hashes={str(p):sha(p) for p in required}
    for relative,expected in frozen['files'].items():
        path=package/relative
        if not path.is_file():raise MissingPrerequisite('D05 required file missing: '+str(path))
        actual=sha(path)
        if actual!=expected:raise ValueError('D05 frozen source hash mismatch: '+str(path))
        hashes[str(path)]=actual
    for entry in summary['input_files']+summary['output_files']:
        if sha(package/entry['path'])!=entry['sha256']:
            raise ValueError('D05 processing-summary input/output hash mismatch')
    for entry in provenance['outputs']:
        if entry['path'].startswith('processed/electron_impact/'):
            if sha(package/entry['path'])!=entry['sha256']:
                raise ValueError('D05 provenance output hash mismatch')
    for name in ('cfit.dat','cfit.txt','cfit.f'):
        path=package/'raw/electron_impact/voronov'/name
        metadata=Path(str(path)+'.metadata.json')
        if not metadata.is_file():raise MissingPrerequisite('D05 native metadata missing')
        if document(metadata)['sha256']!=sha(path):raise ValueError('D05 raw native metadata hash mismatch')
        hashes[str(metadata)]=sha(metadata)
    return package,hashes

def independent_rows(package):
    """Parse native fixed columns independently of the C++ production parser."""
    raw=package/'raw/electron_impact/voronov/cfit.dat'
    rows={}
    for number,line in enumerate(raw.read_text().splitlines(),1):
        z,n=int(line[6:8]),int(line[9:11])
        rows[z,n]=dict(z=z,n=n,line=number,energy=Decimal(line[12:19].strip()),
            p=Decimal(line[20:21]),a=Decimal(line[22:31].strip()),
            x=Decimal(line[32:38].strip()),k=Decimal(line[39:43].strip()),
            minimum=Decimal(line[44:48].strip()),maximum=1000*Decimal(line[49:53].strip()))
    # Fortran DATA rows are a distinct source representation of all406 fits.
    text=(package/'raw/electron_impact/voronov/cfit.f').read_text()
    blocks=re.findall(r'DATA\(CF\(I,\s*(\d+),\s*(\d+)\),I=1,5\)/([^/]+)/',text)
    if len(rows)!=406 or len(blocks)!=406:raise ValueError('D05 incomplete native/Fortran fit set')
    for z,n,values in blocks:
        actual=[Decimal(v.strip()) for v in values.split(',')]
        row=rows[int(z),int(n)]
        if actual!=[row[k] for k in ('energy','p','a','x','k')]:
            raise ValueError('native/Fortran coefficient disagreement')
    selected={name:rows[z,z] for name,z in NAMES.items()}
    with (package/'processed/electron_impact/lunar_neutral_coefficients.csv').open(newline='') as stream:
        csv_rows=list(csv.DictReader(stream))
    if len(csv_rows)!=4:raise ValueError('selected CSV count disagreement')
    for record in csv_rows:
        row=selected[record['element']]
        if int(record['initial_charge'])!=0 or int(record['atomic_number'])!=row['z'] or int(record['electron_count'])!=row['n']:
            raise ValueError('D05 processed neutral stage mismatch')
        for field,key in (('ionization_energy_eV','energy'),('fit_P','p'),('fit_A_cm3_s','a'),('fit_X','x'),('fit_K','k'),('temperature_min_eV','minimum')):
            if Decimal(record[field])!=row[key]:raise ValueError('D05 processed coefficient mismatch')
        if Decimal(record['temperature_max_keV'])*1000!=row['maximum']:
            raise ValueError('D05 native Tmax conversion mismatch')
    return selected

def reference(row,temperature):
    """Independent direct Decimal expression; production uses a log expression."""
    with localcontext() as context:
        context.prec=70
        t=Decimal(str(temperature));u=row['energy']/t
        return row['a']*(1+row['p']*u.sqrt())/(row['x']+u)*(row['k']*u.ln()).exp()*(-u).exp()*Decimal('1e-6')

def reference_log_slope(row, temperature):
    """Independent log-difference of the direct expression, with Decimal headroom."""
    with localcontext() as context:
        context.prec=70
        t=Decimal(temperature)
        delta=Decimal('1e-18')
        # This evaluates the reference expression analytically across each
        # boundary. The actual production query still obeys the native range.
        return (reference(row,t*delta.exp()).ln()-reference(row,t*(-delta).exp()).ln())/(2*delta)

def grid(row):
    with localcontext() as context:
        context.prec=70
        lower=row['minimum'];upper=row['maximum'];span=(upper/lower).ln()
        points=[float(lower*(span*Decimal(i)/32).exp()) for i in range(33)]
    points[0]=float(lower);points[-1]=float(upper)
    return points

def run_command(command,directory,label,stdin=None):
    process=subprocess.run(command,input=stdin,text=True,capture_output=True,check=False)
    (directory/(label+'.stdout.log')).write_text(process.stdout)
    (directory/(label+'.stderr.log')).write_text(process.stderr)
    return process,{'command':command,'returncode':process.returncode,
                    'stdout':str(directory/(label+'.stdout.log')),'stderr':str(directory/(label+'.stderr.log'))}

def values(stdout):
    return dict(line.split('=',1) for line in stdout.splitlines() if '=' in line)

def adjacent_positive_float(value,increase):
    """One binary64 ULP for positive boundary controls (also on Python3.8)."""
    bits=struct.unpack('>Q',struct.pack('>d',value))[0]
    return struct.unpack('>d',struct.pack('>Q',bits+(1 if increase else -1)))[0]

def render_comparison(comparisons,work):
    """Plot unchanged mathematical comparison rows; never observations.

    Matplotlib is optional for coefficient execution. Missing plotting support
    is an explicit visualization SKIP, not a change to the numerical result.
    A symmetric-log error axis retains exact zeros without plotting a made-up
    positive residual. PNG/EPS use the identical declared grid and CSV values.
    """
    try:
        os.environ.setdefault('MPLCONFIGDIR',str(work/'matplotlib-config'))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as error:
        return {'status':'SKIPPED','reason':'optional matplotlib unavailable: '+str(error)}
    figure,axes=plt.subplots(2,1,figsize=(8,8),sharex=True)
    for name in NAMES:
        selected=[row for row in comparisons if row['species']==name]
        temperature=[float(row['temperature_ev']) for row in selected]
        line,=axes[0].plot(temperature,[row['independent_coefficient_m3_s'] for row in selected],label=name+' Decimal reference')
        axes[0].plot(temperature,[row['production_coefficient_m3_s'] for row in selected],
            'o',color=line.get_color(),markersize=3,fillstyle='none',label=name+' production')
        axes[1].plot(temperature,[row['coefficient_relative_error'] for row in selected],
            'o-',color=line.get_color(),markersize=3,label=name)
    axes[0].set_xscale('log');axes[0].set_yscale('log')
    axes[0].set_ylabel('Rate coefficient [m$^3$ s$^{-1}$]')
    axes[0].set_title('U16 mathematical reference comparison — no observations')
    axes[0].legend(fontsize=8,ncol=2,framealpha=1.0)
    axes[1].set_yscale('symlog',linthresh=1e-16)
    axes[1].axhline(RELATIVE_TOLERANCE,color='black',linestyle='--',label='declared 1e-12 numerical limit')
    axes[1].set_ylabel('Coefficient relative error');axes[1].set_xlabel('Electron temperature [eV]')
    axes[1].legend(fontsize=8,framealpha=1.0)
    for axis in axes:axis.grid(True,color='0.85')
    figure.tight_layout()
    artifacts=[]
    for extension in ('png','eps'):
        path=work/('coefficient_comparison.'+extension)
        figure.savefig(path,dpi=180)
        artifacts.append({'path':str(path),'sha256':sha(path)})
    plt.close(figure)
    return {'status':'PASS','scope':'visualization-of-mathematical-reference-not-observations','artifacts':artifacts}

def run_u16(build_dir,output_root):
    """Public common-runner hook; standalone compile does not use generated AMPS.

    A unique evidence directory preserves previous numerical results even when
    the common runner updates its per-ID summary. build_dir is provenance only;
    this test compiles exactly the production pure service and a thin adapter.
    """
    base=Path(output_root).resolve()/'U16';base.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='maxwellian-kernel-',dir=base))
    result={'status':'ERROR','scope':'local-production-kernel-not-linked-chemistry-or-validation',
            'artifact_directory':str(work),'checks':[],'probes':[],
            'generated_build_not_used':str(Path(build_dir).resolve())}
    try:
        revision=subprocess.run(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True,
                                capture_output=True,check=False)
        result['source_revision']=revision.stdout.strip() if revision.returncode==0 else 'unavailable'
        data_root=Path(os.environ.get('MOON_VALIDATION_DATA_ROOT',str(DATA_ROOT))).resolve()
        package,hashes=preflight(data_root)
        rows=independent_rows(package)
        result['reference_identity']={'native_atomic_stage_rows':406,
            'independent_Fortran_DATA_rows':406,'selected_neutral_rows':4,
            'all_five_fit_fields_match':True,
            'domain_ev':{name:[str(row['minimum']),str(row['maximum'])] for name,row in rows.items()}}
        raw=package/'raw/electron_impact/voronov/cfit.dat'
        compiler=shlex.split(os.environ.get('CXX','g++'))
        if not compiler or shutil.which(compiler[0]) is None:raise MissingPrerequisite('C++ compiler unavailable')
        executable=work/'electron-impact-probe'
        command=compiler+['-std=c++17','-O2','-Wall','-Wextra','-pedantic',
            '-I',str(ROOT/'srcMoon/species'),str(HERE/'probe.cpp'),
            str(ROOT/'srcMoon/species/ElectronImpact.cpp'),
            *(str(p) for p in shared_sources(ROOT,'electron-impact')),'-o',str(executable)]
        compiled,record=run_command(command,work,'compile')
        result['probes'].append(record)
        if compiled.returncode:raise RuntimeError('production standalone compile failed')
        version,_=run_command(compiler+['--version'],work,'compiler-version')
        result['compiler']=version.stdout.splitlines()[0]
        source_paths=[ROOT/'srcMoon/species/ElectronImpact.cpp',ROOT/'srcMoon/species/ElectronImpact.h',
                      HERE/'probe.cpp',Path(__file__),HERE/'CONTRACT.md',HERE/'reference/source_manifest.json']
        source_paths += shared_dependencies(ROOT,'electron-impact')
        result['source_sha256']={str(path):sha(path) for path in source_paths}
        result['roadmap_sha256']={str(path):sha(path) for path in (
            ROOT/'srcMoon/development_and_validation.md',data_root/'development_and_validation.md')}
        result['roadmap_control']='source: U16/API1e-12; external: immutable D05 acquisition; differences retained, not silently merged'
        result['data_hashes']=hashes
        (work/'data_hashes.json').write_text(json.dumps(hashes,indent=2,sort_keys=True)+'\n')
        def check(name,ok,details=None):
            result['checks'].append({'name':name,'status':'PASS' if ok else 'FAIL','details':details or {}})
        def query(mode,path,name,temperature,density,label):
            process,rec=run_command([str(executable),mode,str(path),name,str(temperature),str(density)],work,label)
            return process,values(process.stdout),rec
        comparisons=[]
        for name,row in rows.items():
            for index,temperature in enumerate(grid(row)):
                text=format(temperature,'.17g');expected=float(reference(row,text))
                process,actual,rec=query('evaluate',raw,name,text,1e6,f'{name}-{index:02d}')
                if process.returncode:raise RuntimeError('production query failed: '+str(rec))
                coefficient=float(actual['coefficient_m3_s']);rate=float(actual['rate_s'])
                error=abs(coefficient-expected)/expected
                rate_error=abs(rate-expected*1e6)/(expected*1e6)
                mapping=int(actual['atomic_number'])==row['z']==int(actual['electron_count'])
                check(f'{name}-{index:02d}',error<=RELATIVE_TOLERANCE and rate_error<=RELATIVE_TOLERANCE and mapping,
                      {'coefficient_relative_error':error,'rate_relative_error':rate_error,'log':rec['stdout']})
                slope=float(actual['log_temperature_derivative'])
                expected_slope=float(reference_log_slope(row,text))
                slope_error=abs(slope-expected_slope)/max(1,abs(expected_slope))
                check(f'log-temperature-slope-{name}-{index:02d}',slope_error<=RELATIVE_TOLERANCE,
                      {'actual':slope,'reference':expected_slope,'error':slope_error,'log':rec['stdout']})
                comparisons.append(dict(species=name,temperature_ev=text,density_m3=1e6,
                    production_coefficient_m3_s=coefficient,independent_coefficient_m3_s=expected,
                    production_rate_s=rate,coefficient_relative_error=error,rate_relative_error=rate_error))
            for density in (0,1e7):
                p,a,rec=query('evaluate',raw,name,10,density,f'{name}-density-{density}')
                expected=float(reference(row,'10'))*density
                check(f'{name}-density-{density}',p.returncode==0 and
                    (float(a['rate_s'])==0 if density==0 else abs(float(a['rate_s'])-expected)/expected<=RELATIVE_TOLERANCE),rec)
        with (work/'rate_vs_temperature.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
        result['visualization']=render_comparison(comparisons,work)
        result['metrics']={'point_count':len(comparisons),'maximum_coefficient_relative_error':max(x['coefficient_relative_error'] for x in comparisons),
                           'maximum_rate_relative_error':max(x['rate_relative_error'] for x in comparisons)}
        # All altered inputs below are local protocol/numerical robustness
        # fixtures. Authoritative preflight would reject their different hash.
        native_lines=raw.read_text().splitlines()
        fixtures=work/'control-fixtures-not-data';fixtures.mkdir()
        variants={
            'missing-row':native_lines[:-1],
            'duplicate-row':native_lines[:-1]+[native_lines[0]],
            'extra-row':native_lines+[native_lines[0]],
            'wrong-stage':native_lines.copy(),
            'malformed-width':native_lines.copy(),
            'invalid-domain':native_lines.copy(),
            'noninteger-bound':native_lines.copy(),
        }
        variants['wrong-stage'][1]='Na   '+native_lines[1][5:]
        variants['malformed-width'][0]=native_lines[0]+' '
        variants['invalid-domain'][1]=native_lines[1][:44]+'   0'+native_lines[1][48:]
        variants['noninteger-bound'][1]=native_lines[1][:44]+' 1.5'+native_lines[1][48:]
        bad_paths=[]
        for label,lines in variants.items():
            path=fixtures/(label+'.dat');path.write_text('\n'.join(lines)+'\n');bad_paths.append(path)
            p,rec=run_command([str(executable),'load',str(path)],work,'control-'+label)
            check('reject-'+label,p.returncode==4 and values(p.stdout).get('loaded')=='0',rec)
        p,rec=run_command([str(executable),'reload',str(raw),str(bad_paths[0])],work,'control-transactional-reload')
        a=values(p.stdout)
        check('transactional-reload',p.returncode==0 and a.get('reload_ok')=='0' and a.get('loaded')=='1' and
              abs(float(a['coefficient_m3_s'])-float(reference(rows['He'],'10')))/float(reference(rows['He'],'10'))<=RELATIVE_TOLERANCE,rec)
        for index,(name,temperature,density) in enumerate((
                ('He',0,1),('He',-1,1),('He','nan',1),('He','inf',1),
                ('He',adjacent_positive_float(1.0,False),1),('He',adjacent_positive_float(20000.0,True),1),
                ('He',10,-1),('He',10,'nan'),('He',10,'inf'),('He',0,0),
                ('he',10,1),('He+',10,1),('H2O',10,1))):
            p,a,rec=query('evaluate',raw,name,temperature,density,'control-invalid-'+str(index))
            check('invalid-input-'+str(index),p.returncode==4 and
                  (a.get('rate_s','-456')=='-456'),rec)
        for mode in ('uninitialized','enum-invalid'):
            p,a,rec=query(mode,raw,'He',10,1,'control-'+mode)
            check(mode,p.returncode==4 and a.get('coefficient_m3_s')=='-123' and a.get('rate_s')=='-456'
                  and a.get('log_temperature_derivative')=='-789',rec)
        for index,temperature in enumerate((11604.518121550082,0,-1,'nan','inf')):
            p,rec=run_command([str(executable),'temperature-k',str(temperature)],work,'kelvin-'+str(index))
            a=values(p.stdout)
            if index==0:
                expected=float(Decimal(str(temperature))*Decimal('1.380649e-23')/Decimal('1.602176634e-19'))
                check('exact-SI-Kelvin-conversion',p.returncode==0 and abs(float(a['temperature_ev'])-expected)/expected<=RELATIVE_TOLERANCE,rec)
            else:check('invalid-Kelvin-'+str(index),p.returncode==4 and a['temperature_ev']=='-123',rec)
        # U>80 is still mathematically positive and representable here. This
        # separate control catches accidental adoption of cfit.f's early-zero
        # shortcut; the deep-underflow control alone could not distinguish it.
        tail=native_lines.copy();tail[1]=tail[1][:12]+f'{100.0:7.1f}'+tail[1][19:]
        path=fixtures/'representable-tail-arithmetic-only.dat';path.write_text('\n'.join(tail)+'\n')
        p,a,rec=query('evaluate',path,'He',1,1e6,'control-representable-tail')
        tail_row=dict(rows['He'],energy=Decimal('100.0'))
        expected=float(reference(tail_row,'1'))
        check('no-artificial-U80-cutoff',p.returncode==0 and float(a['coefficient_m3_s'])>0 and
              abs(float(a['coefficient_m3_s'])-expected)/expected<=RELATIVE_TOLERANCE,rec)
        underflow=native_lines.copy();underflow[1]=underflow[1][:12]+f'{9999.9:7.1f}'+underflow[1][19:]
        path=fixtures/'underflow-arithmetic-only.dat';path.write_text('\n'.join(underflow)+'\n')
        p,a,rec=query('evaluate',path,'He',1,1e6,'control-natural-underflow')
        check('natural-representability-underflow',p.returncode==0 and float(a['coefficient_m3_s'])==0 and float(a['rate_s'])==0,rec)
        result['fortran']=run_fortran(package,rows,work)
        if result['fortran']['status']=='FAIL':check('native-Fortran-reference',False,result['fortran'])
        result['status']='FAIL' if any(x['status']=='FAIL' for x in result['checks']) else 'PASS'
    except MissingPrerequisite as error:
        result.update(status='SKIPPED',reason=str(error))
    except (OSError,ValueError,KeyError,RuntimeError,ArithmeticError) as error:
        result.update(status='ERROR',reason=str(error))
    result['counts']={state:sum(x['status']==state for x in result['checks']) for state in ('PASS','FAIL','ERROR','SKIPPED')}
    (work/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

def run_fortran(package,rows,work):
    compiler=shlex.split(os.environ.get('FC','gfortran'))
    if not compiler or shutil.which(compiler[0]) is None:
        return {'status':'SKIPPED','reason':'optional original-Fortran compiler unavailable'}
    driver=work/'native_fortran_driver.f90'
    driver.write_text('program verify\nimplicit none\ninteger z,n,ios\nreal t,c,te\ndo\nread(*,*,iostat=ios) z,n,t\nif(ios/=0)exit\ncall cfit(z,n,t,c)\nte=t*8.617385e-05\nwrite(*,\'(I3,1X,I3,1X,ES24.16,1X,ES24.16)\') z,n,te,c\nend do\nend program\n')
    executable=work/'original-cfit-fortran'
    p,compiled=run_command(compiler+['-O0','-ffixed-line-length-none',str(package/'raw/electron_impact/voronov/cfit.f'),str(driver),'-o',str(executable)],work,'fortran-compile')
    if p.returncode:raise RuntimeError('available original-Fortran compile failed')
    requests=[(name,row,t) for name,row in rows.items() for t in grid(row)]
    stdin=''.join(f"{row['z']} {row['n']} {t/8.617385e-5:.17g}\n" for name,row,t in requests)
    p,execution=run_command([str(executable)],work,'fortran-execute',stdin)
    lines=p.stdout.splitlines()
    if p.returncode or len(lines)!=len(requests):raise RuntimeError('original Fortran reference driver failed')
    records=[]
    for (name,row,t),line in zip(requests,lines):
        z,n,effective,actual=line.split();effective=float(effective);actual=float(actual)*1e-6
        expected=float(reference(row,str(effective)));relative=abs(actual-expected)/expected
        # Conservative a-priori single-precision operation/conditioning budget,
        # separate from the unchanged1e-12 production-vs-Decimal gate. Original
        # cfit constants/arithmetic are REAL, and Te is its actual rounded value.
        bound=64.0*2.0**-23*(1+float(row['energy'])/effective)
        records.append({'species':name,'requested_temperature_ev':t,'fortran_effective_temperature_ev':effective,
                        'fortran_coefficient_m3_s':actual,'independent_coefficient_m3_s':expected,
                        'relative_error':relative,'declared_single_precision_bound':bound,
                        'status':'PASS' if int(z)==row['z'] and int(n)==row['n'] and math.isfinite(actual) and relative<=bound else 'FAIL'})
    (work/'fortran_comparison.json').write_text(json.dumps(records,indent=2)+'\n')
    return {'status':'PASS' if all(x['status']=='PASS' for x in records) else 'FAIL','scope':'original-single-precision-supporting-oracle-not-production-tolerance',
            'maximum_relative_error':max(x['relative_error'] for x in records),'comparison_count':len(records),
            'temperature_constant_ev_per_k':8.617385e-5,'bound':'64*2^-23*(1+U); separate from1e-12 production Decimal gate',
            'compile':compiled,'execute':execution,'artifact':str(work/'fortran_comparison.json')}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--build-dir',type=Path,default=ROOT/'build')
    args=parser.parse_args();result=run_u16(args.build_dir,args.output_dir)
    print(json.dumps(result,indent=2,sort_keys=True))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[result['status']])
