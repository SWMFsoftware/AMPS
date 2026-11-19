#!/usr/bin/env python3
"""Actual Moon/PIC constant-history reaction controller, supporting controls.

Consumes the SAME successful fresh all-species backend producer; never changes
configuration/builds or activates a physical lunar channel. Mechanical fixtures
and expected capacity refusal are not Na/CEX kinetic validation or observations.
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
from run_m3_moon_force_tracers import native_sources
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command


def expected_capacity_refusal(returncode,transcript):
    """Only the probe's documented exit(1) is the intended resource control.

    A signal, launcher error or arbitrary nonzero exit remains infrastructure
    ERROR even if buffered stdout happens to contain an earlier refusal text.
    """
    return type(returncode) is int and returncode==1 and isinstance(transcript,str) and all(text in transcript for text in (
        'complete reaction packet cannot fit available PIC particle slots',
        'Increase particle-buffer capacity','without discarding products'))


def run(output,producer_manifest,wall_seconds=600):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    producer_manifest=Path(producer_manifest).resolve()
    result={'schema':'srcMoon-M3-constant-controller-controls-v1','status':'ERROR',
        'scope':'supporting actual PIC/controller/mover; analytical constant hazards and nominal zero-Q mechanics only',
        'formal_I21_I24_status':'NOT EVALUATED','observational_scoring':False,'commands':[]}
    manifest={'execution':'not-run','source_root':str(REPO),'external_observations':[],
        'constant_hazards':'explicit .25/.75 s^-1 analytical controls, not lunar capture/photo rates',
        'source_clocks':'explicit open uniforms from exp(-.4), exp(-.6), exp(-1)',
        'random_seed':None,'mpi_layouts':[1,8],'orders':[1,2],'omp_threads':1}
    budget=CampaignBudget(wall_seconds,min(60,wall_seconds*.1))
    def execute(command,label,cwd=REPO,ranks=1,expected_failure=False):
        receipt=budget.run(command,cwd,output/(label+'.log'),output=output,allocated_cores=ranks,
            environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        result['commands'].append(receipt)
        if receipt.get('reason') or (receipt['returncode'] not in (0,1) and not expected_failure):
            raise RuntimeError(label+' infrastructure failure: '+str(receipt.get('reason') or receipt['returncode']))
        return receipt
    try:
        missing=[tool for tool in (os.environ.get('MPICXX','mpicxx'),'mpiexec') if shutil.which(tool) is None]
        if missing:result.update(status='SKIPPED',reason='native dependency unavailable',missing=missing);return result
        try:
            import matplotlib
        except ImportError:
            result.update(status='SKIPPED',reason='Matplotlib unavailable for required analytical comparison figures');return result
        if not producer_manifest.is_file():result.update(status='SKIPPED',reason='successful all-species backend producer unavailable');return result
        producer=json.loads(producer_manifest.read_text());producer_result=producer_manifest.with_name('result.json')
        if json.loads(producer_result.read_text()).get('status')!='PASS':
            result.update(status='SKIPPED',reason='actual backend predecessor is not PASS');return result
        if Path(producer['source_root']).resolve()!=REPO:raise ValueError('producer belongs to another source checkout')
        library=Path(producer['library']['path']);build=library.parent
        if sha256(library)!=producer['library']['sha256']:raise ValueError('retained production library changed')
        for path,digest in producer['generated_snapshots'].items():
            relative=Path(path).relative_to(producer_manifest.parent/'generated')
            if sha256(Path(path))!=digest or sha256(build/relative)!=digest:raise ValueError('producer/current generated ABI mismatch')
        if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):
            raise ValueError('controller/production source is not the actual compiled producer')
        manifest.update(producer={'path':str(producer_manifest),'sha256':sha256(producer_manifest),
            'result':str(producer_result),'result_sha256':sha256(producer_result)},library=producer['library'],
            production_source_hashes=producer['production_source_hashes'],generated_snapshots=producer['generated_snapshots'])
        sources=[Path(__file__).resolve(),Path(__file__).resolve().with_name('m3_native_reaction_controller_probe.cpp'),
            Path(__file__).resolve().with_name('M3_REACTION_CONTROLLER_CONTRACT.md')]
        manifest['test_sources']={str(path):sha256(path) for path in sources}
        for path in sources:shutil.copy2(path,output/path.name)
        executable=output/'native-reaction-controller-probe';command=_probe_compile_command(build,executable);command.insert(1,'-DMPI_ON')
        command[command.index(str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp'))]=str(sources[1])
        if execute(command,'compile-probe')['returncode']!=0:raise RuntimeError('thin native controller probe compilation failed')
        manifest['executable']={'path':str(executable),'sha256':sha256(executable)};manifest['execution']='actual normal PIC walker + Moon mover + combined controller'
        required={'variable_source_hazard_refused_without_activation','missing_explicit_source_state_resolver_refused',
            'complete_actual_core_hook_installation','recycled_native_slot_starts_exactly_one_fresh_source_clock',
            'actual_face_boundary_macro_prepares_source_and_partial_children','actual_boundary_macro_prepares_source_and_partial_children',
            'strict_nextafter_future_deadline_does_not_fire_at_endpoint',
            'exact_internal_face_endpoint_event_commits_zero_time_children',
            'surviving_tiny_clock_continues_next_outer_interval_without_source_redraw',
            'repeated_refinement_previews_change_no_clock_or_rng','explicit_state_and_source_clock_once_per_birth',
            'combined_hazard_and_product_draw_counts_no_reseed_after_mpi','selected_source_actual_event_epoch_and_remaining_interval',
            'accepted_exact_material_and_represented_number_budget','all_daughters_complete_original_physical_interval_without_birth_dt_factor',
            'actual_pending_clock_crosses_mpi_before_event','no_unfinished_child_or_parent_continuation','no_parent_double_delete_or_native_slot_leak'}
        counts=collections.Counter();layouts=[];comparisons=[]
        for ranks in (1,8):
            for order in (1,2):
                label=f'np{ranks}-order{order}';directory=output/label;directory.mkdir();runtime=directory/'runtime';runtime.mkdir();(runtime/'PT').mkdir()
                command=['mpiexec','-n',str(ranks),str(executable),str(order),'normal',str(directory/'controller')]
                receipt=execute(command,label+'/stdout',runtime,ranks)
                transcript=(directory/'stdout.log').read_text();summaries=re.findall(r'M3 native reaction controller order=(\d+) checks=(\d+) status=(PASS|FAIL|ERROR)',transcript)
                if len(summaries)!=1 or int(summaries[0][0])!=order or summaries[0][2]=='ERROR':raise ValueError('incomplete/error native controller summary')
                with (directory/'controller.checks.csv').open(newline='') as stream:checks=list(csv.DictReader(stream))
                if len(checks)!=len(required) or {row['check'] for row in checks}!=required or int(summaries[0][1])!=len(checks):
                    raise ValueError('missing/duplicate native controller contracts')
                if any(row['status'] not in ('PASS','FAIL') for row in checks):raise ValueError('invalid native outcome status')
                failed=any(row['status']=='FAIL' for row in checks)
                if (receipt['returncode']==1)!=failed or (summaries[0][2]=='FAIL')!=failed:raise ValueError('exit/summary/check mismatch')
                candidates=[];particles=[]
                for rank in range(ranks):
                    for suffix,rows in [('candidates',candidates),('particles',particles)]:
                        with (directory/f'controller.rank{rank}.{suffix}.csv').open(newline='') as stream:rows.extend(csv.DictReader(stream))
                boundary=[row for row in candidates if row['case'] in ('boundary','boundary-face')];candidates=[row for row in candidates if row['case']=='main']
                if not failed and (len(boundary)!=2 or {row['case'] for row in boundary}!={'boundary','boundary-face'}):
                    raise ValueError('actual source-boundary macro diagnostics incomplete')
                if not failed and (len(candidates)!=3 or {int(row['source_tag']) for row in candidates}!={0,1,2} or len(particles)!=6):
                    raise ValueError('native candidate/product evidence incomplete')
                for row in candidates:
                    tag=int(row['source_tag']);delay=(.4,.6,1)[tag]
                    if not all(math.isfinite(float(value)) for key,value in row.items() if key!='case'):raise ValueError('nonfinite candidate diagnostic')
                    comparisons.append({'mpi_ranks':ranks,'order':order,'source_tag':tag,'expected_delay_s':delay,
                        'actual_delay_s':float(row['event_epoch_s'])-1000,'expected_channel':2*tag+(tag==1),
                        'actual_channel':int(row['channel']),'remaining_s':float(row['remaining_s'])})
                counts.update(row['status'] for row in checks);layouts.append({'mpi_ranks':ranks,'order':order,'status':'FAIL' if failed else 'PASS','checks':len(checks)})
        # Dedicated one-rank resource-negative run. An unrelated MPI/build
        # failure is ERROR; only the exact controller complete-packet refusal
        # with actionable instructions establishes this negative control.
        directory=output/'capacity-refusal';directory.mkdir();runtime=directory/'runtime';runtime.mkdir();(runtime/'PT').mkdir()
        command=['mpiexec','-n','1',str(executable),'2','capacity',str(directory/'controller')]
        receipt=execute(command,'capacity-refusal/stdout',runtime,1,True);transcript=(directory/'stdout.log').read_text()
        refusal=expected_capacity_refusal(receipt['returncode'],transcript)
        if not refusal:raise ValueError('expected resource ERROR was not the specific complete-packet controller refusal')
        counts['PASS']+=1
        with (output/'event_comparison.csv').open('w',newline='') as stream:
            fields=['mpi_ranks','order','source_tag','expected_delay_s','actual_delay_s','expected_channel','actual_channel','remaining_s']
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(comparisons)
        # Deterministic analytical event-time comparison, not observations or
        # fitted data. The second panel exposes floating clock representation
        # rather than hiding it beneath the equality line in the first panel.
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        figure,axes=plt.subplots(1,2,figsize=(10,4))
        axes[0].plot([0,1.1],[0,1.1],color='black',label='analytical equality')
        for ranks in (1,8):
            for order in (1,2):
                rows=[row for row in comparisons if row['mpi_ranks']==ranks and row['order']==order]
                x=[row['expected_delay_s'] for row in rows];actual=[row['actual_delay_s'] for row in rows]
                label=f'{ranks}MPI, order{order}'
                axes[0].scatter(x,actual,label=label)
                axes[1].scatter(x,[a-b for a,b in zip(actual,x)],label=label)
        axes[0].set_ylabel('Actual physical event delay [s]');axes[1].set_ylabel('Actual minus analytical delay [s]')
        for axis in axes:axis.set_xlabel('Analytical H / combined rate [s]');axis.grid(True)
        axes[1].legend(framealpha=1);figure.suptitle('Analytical controller controls only: no observations or physical CEX validation')
        figure.tight_layout()
        for suffix in ('png','eps'):figure.savefig(output/('event_comparison.'+suffix))
        plt.close(figure)
        manifest['artifacts']={str(path.relative_to(output)):sha256(path) for path in output.rglob('*') if path.is_file()}
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',layouts=layouts,capacity_refusal={'status':'PASS','expected_child_status':'ERROR','log':str(directory/'stdout.log')},
            counts={name:counts[name] for name in ('PASS','FAIL','ERROR','SKIPPED')},event_comparisons=comparisons)
    except (OSError,ValueError,RuntimeError,KeyError,TypeError,ImportError) as error:result.update(status='ERROR',reason=str(error))
    finally:
        result['elapsed_wall_seconds']=time.monotonic()-started
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
        (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--producer-manifest',type=Path,required=True);p.add_argument('--wall-seconds',type=float,default=600)
    a=p.parse_args()
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error('output must be new or empty')
    if not 0<a.wall_seconds<=3600:p.error('wall seconds must be in(0,3600]')
    r=run(a.output_dir,a.producer_manifest,a.wall_seconds);print(json.dumps(r,indent=2,sort_keys=True))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
