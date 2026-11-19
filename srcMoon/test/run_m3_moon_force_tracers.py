#!/usr/bin/env python3
"""Verify both actual Moon force tracer dispatches using a retained native build.

No application configure/build is performed. The successful real PIC backend
producer receipt, actual library/headers and current production source must
match. RK4 shares the force and qualifies integration only, not physical truth.
"""
from __future__ import annotations
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256
from integration_common import REPO
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command


def native_sources(records):
    """Compare actual production inputs; separately retain new thin test code.

    Adding a probe cannot alter already compiled library objects. Exclude the
    test directory only from this source-delta comparison, retaining the full
    original producer identity and independently hashing this test's sources.
    """
    return {key:value for key,value in records.items() if not key.startswith('srcMoon/test/')}


def run(output,producer_manifest,wall_seconds=600):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    producer_manifest=Path(producer_manifest).resolve();started=time.monotonic()
    result={'schema':'srcMoon-M3-actual-Moon-force-tracers-v1','status':'ERROR',
        'scope':'supporting frozen-snapshot production force/mover integration only',
        'formal_I21_I24_status':'NOT EVALUATED','observational_scoring':False,'commands':[]}
    manifest={'execution':'not-run','source_root':str(REPO),'external_observations':[],
        'reference_scope':'independently integrated RK4 using actual production force; not independent force truth',
        'mpi_layouts':[1,8],'orders':[1,2],'omp_threads':1,'random_seed':None}
    budget=CampaignBudget(wall_seconds,min(60,wall_seconds*.1))
    def execute(command,label,cwd=REPO,native=False,ranks=1):
        receipt=budget.run(command,cwd,output/(label+'.log'),output=output,allocated_cores=ranks,
            environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        result['commands'].append(receipt)
        if receipt['returncode']!=0 and not(native and receipt['returncode']==1 and not receipt.get('reason')):
            raise RuntimeError(label+' failed: '+str(receipt.get('reason') or receipt['returncode']))
        return receipt
    try:
        missing=[tool for tool in (os.environ.get('MPICXX','mpicxx'),'mpiexec') if shutil.which(tool) is None]
        if missing:
            result.update(status='SKIPPED',reason='required native tool unavailable',missing=missing);return result
        try:
            import matplotlib
        except ImportError:
            result.update(status='SKIPPED',reason='Matplotlib unavailable for required comparison artifacts');return result
        if not producer_manifest.is_file():
            result.update(status='SKIPPED',reason='successful native producer build receipt unavailable');return result
        producer=json.loads(producer_manifest.read_text());producer_result=producer_manifest.with_name('result.json')
        if json.loads(producer_result.read_text()).get('status')!='PASS':
            result.update(status='SKIPPED',reason='real PIC backend producer campaign is not PASS');return result
        if Path(producer['source_root']).resolve()!=REPO:
            raise ValueError('producer receipt belongs to another checkout; use its retained source/build')
        library=Path(producer['library']['path']);build=library.parent
        if sha256(library)!=producer['library']['sha256']:raise ValueError('producer library bytes changed')
        for path,digest in producer['generated_snapshots'].items():
            if sha256(Path(path))!=digest:raise ValueError('producer generated snapshot corrupted')
            relative=Path(path).relative_to(producer_manifest.parent/'generated')
            if sha256(build/relative)!=digest:raise ValueError('current generated ABI differs from producer')
        if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):
            raise ValueError('current production source differs from actual compiled producer')
        manifest.update(producer={'manifest':str(producer_manifest),'sha256':sha256(producer_manifest),
            'result':str(producer_result),'result_sha256':sha256(producer_result)},library=producer['library'],
            production_source_hashes=producer['production_source_hashes'],generated_snapshots=producer['generated_snapshots'])
        for source in (producer_manifest,producer_result):shutil.copy2(source,output/('producer_'+source.name))
        sources=[Path(__file__).resolve(),Path(__file__).resolve().with_name('m3_native_moon_force_tracer_probe.cpp'),
            Path(__file__).resolve().with_name('M3_MOON_FORCE_TRACER_CONTRACT.md')]
        manifest['test_sources']={str(path):sha256(path) for path in sources}
        for path in sources:shutil.copy2(path,output/path.name)
        executable=output/'native-moon-force-tracer-probe';command=_probe_compile_command(build,executable)
        command.insert(1,'-DMPI_ON') # same production preprocessing for inline mesh/channel APIs
        command[command.index(str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp'))]=str(sources[1])
        execute(command,'compile-probe');manifest['executable']={'path':str(executable),'sha256':sha256(executable)}
        outcomes=[];states={};reference_metrics=[];total=collections.Counter()
        expected=collections.Counter({'all_applicable_force_terms_nonzero':2,'actual_total_force_initial_wiring':2,
            'global_species_weight_clock_and_physical_interval':6,'normal_total_acceleration_macro_dispatch':6,
            'position_and_velocity_errors_decrease':2,'declared_first_or_second_order_convergence':2,
            'configured_boundary_birth_clock_and_step_endpoint':1,
            'configured_boundary_face_birth_clock_and_step_endpoint':1})
        for ranks in (1,8):
            for order in (1,2):
                name=f'np{ranks}-order{order}';directory=output/name;directory.mkdir();runtime=directory/'runtime';runtime.mkdir();(runtime/'PT').mkdir()
                manifest['execution']='actual normal Moon::MoveParticle/TotalParticleAcceleration'
                command=['mpiexec','-n',str(ranks),str(executable),str(order),str(directory/'moon-force')]
                execution=execute(command,name+'/stdout',runtime,True,ranks)
                (directory/'stderr.log').write_text('Combined MPI stdout/stderr is retained in stdout.log.\n')
                transcript=(directory/'stdout.log').read_text();summaries=re.findall(
                    r'M3 actual Moon force tracer order=(\d+) checks=(\d+) status=(PASS|FAIL|ERROR)',transcript)
                if len(summaries)!=1 or int(summaries[0][0])!=order:raise ValueError('missing actual native completion summary')
                if summaries[0][2]=='ERROR':raise ValueError('native completion reports ERROR')
                with (directory/'moon-force.checks.csv').open(newline='') as stream:checks=list(csv.DictReader(stream))
                if collections.Counter(row['check'] for row in checks)!=expected or len(checks)!=int(summaries[0][1]):
                    raise ValueError('incomplete native force/integrator controls')
                if any(row['status'] not in ('PASS','FAIL') for row in checks):raise ValueError('invalid native check status')
                failed=any(row['status']=='FAIL' for row in checks)
                if (execution['returncode']==1)!=failed or (summaries[0][2]=='FAIL')!=failed:
                    raise ValueError('native exit/summary/evaluated check disagreement')
                with (directory/'moon-force.states.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
                if len(rows)!=6 or {(int(row['species']),float(row['h_s'])) for row in rows}!={(s,h) for s in (0,1) for h in (1,.5,.25)}:
                    raise ValueError('incomplete two-species three-resolution state evidence')
                for row in rows:
                    if int(row['order'])!=order or not all(math.isfinite(float(v)) for v in row.values()):raise ValueError('nonfinite/invalid state evidence')
                with (directory/'moon-force.reference.csv').open(newline='') as stream:proof=list(csv.DictReader(stream))
                if len(proof)!=2 or {int(row['species']) for row in proof}!={0,1}:
                    raise ValueError('missing two-species reference-resolution evidence')
                for row in proof:
                    if (int(row['coarse_steps'])!=4096 or int(row['fine_steps'])!=8192 or
                        not all(math.isfinite(float(value)) for value in row.values()) or
                        not 0<=float(row['position_error_fraction'])<=.001 or
                        not 0<=float(row['velocity_error_fraction'])<=.001):
                        raise ValueError('reference-resolution evidence contradicts native completion')
                    for difference,finest,fraction in (
                        ('position_difference_m','finest_position_error_m','position_error_fraction'),
                        ('velocity_difference_m_s','finest_velocity_error_m_s','velocity_error_fraction')):
                        if not math.isclose(float(row[difference])/float(row[finest]),float(row[fraction]),rel_tol=1e-14,abs_tol=0):
                            raise ValueError('reference-resolution fraction/retained metric mismatch')
                    reference_metrics.append(dict(row,mpi_ranks=ranks,order=order))
                states[ranks,order]=rows;total.update(row['status'] for row in checks)
                outcomes.append({'mpi_ranks':ranks,'order':order,'status':'FAIL' if failed else 'PASS','checks':len(checks),
                    'path':str(directory)})
        layout_ok=True
        for order in (1,2):
            for a,b in zip(states[1,order],states[8,order]):
                for key in a:
                    av,bv=float(a[key]),float(b[key])
                    layout_ok=layout_ok and abs(av-bv)<=16*sys.float_info.epsilon*max(1,abs(av),abs(bv))
        total['PASS' if layout_ok else 'FAIL']+=1
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figure,axes=plt.subplots(1,2,figsize=(10,4))
        for order in (1,2):
            for species in (0,1):
                rows=[row for row in states[1,order] if int(row['species'])==species]
                for axis,key,label in ((axes[0],'position_error_m','Position error [m]'),(axes[1],'velocity_error_m_s','Velocity error [m/s]')):
                    axis.loglog([float(row['h_s']) for row in rows],[float(row[key]) for row in rows],marker='o',label=f'{"Na" if species==0 else "Na+"}, order{order}')
                    axis.set_xlabel('Maximum physical substep [s]');axis.set_ylabel(label);axis.grid(True)
        axes[1].legend();figure.suptitle('Supporting integration only: actual Moon force vs independent RK4 integration')
        figure.tight_layout()
        for extension in ('png','eps'):figure.savefig(output/('moon-force-convergence.'+extension))
        plt.close(figure)
        artifacts={str(path.relative_to(output)):sha256(path) for path in output.rglob('*') if path.is_file()}
        manifest['retained_artifacts']=artifacts
        result.update(status='FAIL' if total['FAIL'] else 'PASS',layouts=outcomes,
            mpi_layout_comparison={'status':'PASS' if layout_ok else 'FAIL','bound':'16 binary64 eps times max(1,absolute state value)'},
            reference_resolution=reference_metrics,
            counts={key:total[key] for key in ('PASS','FAIL','ERROR','SKIPPED')})
    except (OSError,ValueError,RuntimeError,KeyError,TypeError,ImportError) as error:
        result.update(status='ERROR',reason=str(error))
    finally:
        result['elapsed_wall_seconds']=time.monotonic()-started
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
        (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--producer-manifest',type=Path,required=True);p.add_argument('--wall-seconds',type=float,default=600)
    a=p.parse_args()
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error('output must be new or empty')
    if not 0<a.wall_seconds<=3600:p.error('wall seconds must be positive and at most3600')
    r=run(a.output_dir,a.producer_manifest,a.wall_seconds);print(json.dumps(r,indent=2,sort_keys=True))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
