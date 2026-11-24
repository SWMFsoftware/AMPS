"""Source-bound first-order surface residence comparison for the U24 owner."""
import csv
from decimal import Decimal as D, localcontext
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
MOON = HERE.parents[1]
ROOT = MOON.parent
sys.path.insert(0,str(ROOT/'utility/chemical_database/python'))
from amps_chemical_database import Database, DatabaseError

DEFAULT_LOCK = Path('/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14/database-lock.json')
LOCK_SHA = '9fd6272f652fa14e33517f1354d9dc3ee6bf47948a07b3438b3b9f4c7b14f0d8'
SOURCE_SHA = '9128965e06ed2fa21b39f2df3f97210ef719fbcbf7b627481ac50cc859c3c10c'
DATASET = 'lunar-oh-surface-model-precedents'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(lock, output,manifest):
    """Return the measured numerical status; the outer U24 gate stays open."""
    manifest.update({'scope':'conditional bound OH/H2O residence; published smooth-surface scenario',
              'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'input_lock':str(lock),'numerical_relative_tolerance':'1e-12','randomness':'none'})
    if not lock.is_file():
        return manifest, 'SKIPPED', 'explicit r14 chemistry lock unavailable'
    if digest(lock)!=LOCK_SHA:
        return manifest, 'ERROR', 'selected lock differs from reviewed r14 identity'
    db=Database.from_lock(lock)
    source,=(entry for entry in db.describe(DATASET)['files'] if entry['id']==DATASET+':laferriere-hops')
    data=source['quantity_metadata']
    if (source['sha256']!=SOURCE_SHA or set(data['species'])!={'OH','H2O'}
            or data['residence_time_s']!='1 / (f A exp(-E_d/(k_B T)))'):
        return manifest, 'ERROR', 'source identity, species or documented residence law changed'
    # f*A is the inverse-time prefactor in this specific published law. Do
    # not drop f or combine it with another model's trapping probability.
    f,attempt,barrier=(D(str(data[key])) for key in
                        ('coverage_fraction','attempt_frequency_s_inv','desorption_barrier_eV'))
    program=output/'residence-probe'
    sources=[HERE/'residence_probe.cpp',MOON/'MoonSurfacePhysics.cpp']
    compiler=shlex.split(os.environ.get('CXX','c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        return manifest,'SKIPPED','required C++17 compiler unavailable'
    argv=[*compiler,'-std=c++17','-O2','-I',str(MOON),*map(str,sources),'-o',str(program)]
    manifest.update(lock_sha256=LOCK_SHA,source_id=source['id'],source_sha256=source['sha256'],
                    numerical_source_sha256={str(p):digest(p) for p in [*sources,MOON/'MoonSurfacePhysics.h',Path(__file__)]},
                    compile_command=argv,execution_command=[str(program)])
    with (output/'compile.log').open('w') as log:
        subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,check=True)
    manifest['executable_sha256']=digest(program)
    with localcontext() as context:
        context.prec=60
        si_joule_per_ev=D('1.602176634e-19')
        kb=D('1.380649e-23')
        barrier_j=barrier*si_joule_per_ev
        rows=[]
        for temperature in map(D,('150','250','400')):
            tau=(barrier_j/(kb*temperature)).exp()/(f*attempt)
            # Figshare v1 lunar_processes.py uses these rounded constants.
            # Preserve that separate source convention as a cross-model
            # diagnostic. The lunar-kernel reference uses exact SI values.
            author_tau=(barrier*D('1.60218e-19')/(D('1.38064e-23')*temperature)).exp()/(f*attempt)
            for ratio in map(D,('.01','1','5')):
                rows.append({'T_K':temperature,'dt_s':ratio*tau,'tau_SI_s':tau,
                             'tau_published_constants_s':author_tau,'P_SI':1-(-ratio).exp()})
        inputs=''.join(f"{f} {attempt} {barrier_j} {row['T_K']} {row['dt_s']}\n" for row in rows)
        (output/'input.txt').write_text(inputs)
        result=subprocess.run([str(program)],input=inputs,capture_output=True,text=True,check=True)
        (output/'stdout.log').write_text(result.stdout)
        (output/'stderr.log').write_text(result.stderr)
        lines=result.stdout.splitlines()
        if len(lines)!=len(rows):
            raise ValueError('incomplete production residence output')
        for row,line in zip(rows,lines):
            tau,probability=map(D,line.split())
            row.update(tau_AMPS_s=tau,P_AMPS=probability)
            row['tau_relative_error']=abs(tau/row['tau_SI_s']-1)
            row['P_relative_error']=abs(probability/row['P_SI']-1)
            row['status']='PASS' if max(row['tau_relative_error'],row['P_relative_error'])<=D('1e-12') else 'FAIL'
    with (output/'comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    # Plot inverse temperature so the Arrhenius prediction is linear on the
    # logarithmic ordinate. Open markers reveal the numerical results even
    # where they coincide exactly with the independent reference.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure,(lifetime,release)=plt.subplots(1,2,figsize=(8,3.5),constrained_layout=True)
    samples=rows[::3]
    inverse=[1000/float(row['T_K']) for row in samples]
    lifetime.semilogy(inverse,[float(row['tau_SI_s']) for row in samples],label='Exact-SI model')
    lifetime.semilogy(inverse,[float(row['tau_AMPS_s']) for row in samples],'o',fillstyle='none',label='AMPS')
    lifetime.set(xlabel='1000 K / T',ylabel='Residence time (s)')
    ratios=[float(row['dt_s']/row['tau_SI_s']) for row in rows[:3]]
    release.semilogx(ratios,[float(row['P_SI']) for row in rows[:3]],'k+',label='Exact-SI model')
    release.semilogx(ratios,[float(row['P_AMPS']) for row in rows[:3]],'o',fillstyle='none',label='AMPS, 150 K')
    release.set(xlabel='Step / residence time',ylabel='Release probability',ylim=(-0.05,1.05))
    for axis in (lifetime,release):
        axis.legend(framealpha=1)
    for extension in ('png','eps'):
        figure.savefig(output/('comparison.'+extension))
    plt.close(figure)
    manifest['checks']={s:sum(row['status']==s for row in rows) for s in ('PASS','FAIL')}
    return manifest, 'FAIL' if manifest['checks']['FAIL'] else 'PASS', 'conditional residence time and finite-step release'


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    child=Path(tempfile.mkdtemp(prefix='residence-',dir=output))
    lock=Path(os.environ.get('AMPS_CHEMISTRY_LOCK',str(DEFAULT_LOCK))).resolve()
    manifest={'input_lock':str(lock)}
    try:
        manifest,status,reason=verify(lock,child,manifest)
    except DatabaseError as error:
        status='SKIPPED' if error.code=='Unavailable' else 'ERROR'
        reason=str(error)
    except ModuleNotFoundError as error:
        status,reason='SKIPPED',str(error)
    except (OSError,ValueError,KeyError,TypeError,ArithmeticError,subprocess.SubprocessError) as error:
        status,reason='ERROR',str(error)
    manifest.update(status=status,reason=reason)
    (child/'result.json').write_text(json.dumps(manifest,indent=2)+'\n')
    log=child/'runner.log'
    log.write_text(f'{status} {reason}\nLock: {lock}\nResult: {child/"result.json"}\n')
    return {'check':'published smooth-surface OH/H2O residence kernel','status':status,'reason':reason,
            'runner':str(Path(__file__).resolve()),'result':str(child/'result.json'),'log':str(log)}
