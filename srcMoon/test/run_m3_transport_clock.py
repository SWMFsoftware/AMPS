#!/usr/bin/env python3
"""Supporting common-physical-time clock controls, never observations.

--local-only compiles the actual PIC-independent clock and native production
parser/installer. Linked run() consumes a successful, hash-bound all-species
producer and compiles a thin normal Moon/PIC clock/source/mover probe, without
configuring a checkout or launching its own application build.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import traceback

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256 as _file_sha256
from integration_common import REPO
from run_m3_moon_force_tracers import native_sources
from test_m3_runtime_rates import change
from chemical_database_build import shared_sources, shared_dependencies
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command


def sha256(path):
    """Normalize JSON path strings before the shared Path-only file hasher.

    Retained manifests necessarily decode paths as strings; treating those
    strings as Path objects raised AttributeError before ERROR could be saved.
    Missing/invalid paths now follow the ordinary caught OSError/TypeError
    protocol and cannot turn source/ABI verification into a successful run.
    """
    return _file_sha256(Path(path))


ARITHMETIC_CONTRACTS={
    'uninitialized Begin preserves caller output',
    'nonfinite origin refused without initialization',
    'explicit origin initializes empty clock',
    'one-shot origin cannot silently reinitialize',
    'zero outer step refused',
    'negative outer step refused',
    'infinite outer step refused',
    'immutable first physical step has exact relative and projected endpoints',
    'nested Begin refused without token mutation',
    'zero offset is beginning birth with full interval',
    'source birth uses relative offset before absolute projection',
    'endpoint source has zero unused seconds',
    'negative source offset refused transactionally',
    'strict future source offset cannot be clamped to endpoint',
    'nonfinite departure offset refused',
    'remaining-interval source has same declared birth',
    'zero remaining source is exact endpoint',
    'full remaining source is exact beginning',
    'negative or oversized source remaining interval refused',
    'unrepresentable positive source residual does not become endpoint birth',
    'forged sequence cannot publish elapsed time',
    'forged duration cannot abort another active step',
    'PRE-WORK abort retains original accepted chronology',
    'duplicate abort refused',
    'fresh accepted step publishes exactly once',
    'duplicate commit cannot advance physical time',
    'source outside active physical step refused',
    'next step starts at accepted transport endpoint instead of frozen origin',
    'stale preceding token cannot commit next interval',
    'repeated source birth advances while external origin remains fixed',
    'external origin is immutable across transport commits',
    '100000 intervals match independent binary-dt product within long-double arithmetic bound',
    'relative compensated clock avoids repeated large-ET addition drift',
    'outer step below absolute epoch representation is explicit refusal',
    'negative explicit TDB epoch is valid and chronologically increasing',
    'derived two-addition epoch projection envelope accepted',
    'larger or nonfinite endpoint discrepancy cannot silently align',
}


def _start(output,wall_seconds,scope):
    output=Path(output).resolve()
    if output.exists() and any(output.iterdir()):raise ValueError('clock evidence output must be new or empty')
    output.mkdir(parents=True,exist_ok=True)
    result={'schema':'srcMoon-M3-common-transport-clock-v1','status':'ERROR','scope':scope,
        'formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','observational_scoring':False,'commands':[]}
    manifest={'execution':'not-run','source_root':str(REPO),'random_seed':None,
        'external_observations':[],'kinetic_providers_activated':False}
    # A direct invocation must retain revision as well as the complete native
    # hashes. Git metadata can be absent in a source archive; record that fact
    # explicitly instead of inventing a commit or borrowing another campaign.
    try:
        revision=subprocess.run(['git','rev-parse','HEAD'],cwd=REPO,
            text=True,capture_output=True,check=False)
        if revision.returncode or not re.fullmatch(r'[0-9a-fA-F]{40,64}',revision.stdout.strip()):
            manifest.update(source_revision='unknown',source_revision_error=revision.stderr.strip() or 'git did not return a revision')
        else:manifest['source_revision']=revision.stdout.strip()
    except OSError as error:
        manifest.update(source_revision='unknown',source_revision_error=str(error))
    return output,result,manifest,CampaignBudget(wall_seconds,min(60,wall_seconds*.1))


def _execute(budget,result,output,command,label,cwd=REPO,ranks=1):
    receipt=budget.run(command,cwd,output/(label+'.log'),output=output,allocated_cores=ranks,
        environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    result['commands'].append(receipt)
    if receipt.get('reason') or receipt.get('returncode') not in (0,1,3,4):
        raise RuntimeError(label+' infrastructure failure: '+str(receipt.get('reason') or receipt.get('returncode')))
    return receipt


def _finish(output,result,manifest,started):
    result['elapsed_wall_seconds']=time.monotonic()-started
    def path_only(value):
        if isinstance(value,Path):return str(value)
        raise TypeError('non-JSON metadata object: '+type(value).__name__)
    try:
        manifest['artifacts']={str(p.relative_to(output)):sha256(p) for p in output.rglob('*') if p.is_file() and p.name not in ('result.json','run_manifest.json')}
        manifest_text=json.dumps(manifest,indent=2,sort_keys=True,default=path_only)
    except Exception as error:
        # Artifact hashing/metadata serialization is infrastructure too. Keep
        # evaluated FAIL details separately if finalization itself fails; never
        # turn a partial or unserializable native receipt into PASS/no receipt.
        result['evaluated_status_before_finalization']=result.get('status')
        _infrastructure_error(output,result,error)
        (output/'metadata_before_finalization_error.txt').write_text(repr(manifest))
        manifest_text=json.dumps({'execution':str(manifest.get('execution','unknown')),
            'source_root':str(REPO),'source_revision':str(manifest.get('source_revision','unknown')),
            'metadata_finalization_error':str(error)},indent=2,sort_keys=True)
    (output/'run_manifest.json').write_text(manifest_text+'\n')
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def _infrastructure_error(output,result,error):
    """Retain programming/parser/infrastructure exceptions at runner boundary.

    Normal child FAIL is evaluated before this path and remains FAIL. This
    handler never catches process signals/SystemExit or invents a criterion;
    it preserves the traceback needed to diagnose e.g. malformed JSON shapes.
    """
    result.update(status='ERROR',reason=str(error),exception_type=type(error).__name__)
    (output/'infrastructure_exception.log').write_text(traceback.format_exc())


def run_local(output,wall_seconds=600):
    """Real clock/parser production code; no generated tree or data read."""
    started=time.monotonic();output,result,manifest,budget=_start(output,wall_seconds,
        'supporting production clock arithmetic and parser/installer; no linked particle/physical qualification')
    try:
        compiler=shutil.which('g++')
        if not compiler:result.update(status='SKIPPED',reason='C++17 compiler unavailable');return _finish(output,result,manifest,started)
        paths=['srcMoon/MoonTransportClock.cpp','srcMoon/MoonTransportClock.h',
            'srcMoon/test/m3_transport_clock_probe.cpp','srcMoon/test/M3_COMMON_TRANSPORT_CLOCK_CONTRACT.md',
            'srcMoon/MoonInput.cpp','srcMoon/MoonInput.h','srcMoon/test/m3_runtime_rates_probe.cpp',
            'srcMoon/test/run_m3_transport_clock.py','srcMoon/species/ElectronImpact.h',
            'srcMoon/species/Photochemistry.h','srcMoon/species/ChargeExchange.h','srcMoon/MoonThermalModel.h',
            'srcMoon/species/ElectronImpact.cpp','srcMoon/species/Photochemistry.cpp','srcMoon/species/ChargeExchange.cpp']
        paths += [str(p.relative_to(REPO)) for p in shared_dependencies(REPO,'runtime')]
        manifest['source_files']={name:sha256(REPO/name) for name in paths}
        exe=output/'clock-probe';command=[compiler,'-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic','-I',str(REPO/'srcMoon'),
            str(REPO/'srcMoon/MoonTransportClock.cpp'),str(REPO/'srcMoon/test/m3_transport_clock_probe.cpp'),'-o',str(exe)]
        if _execute(budget,result,output,command,'compile-clock')['returncode']!=0:raise RuntimeError('production clock controls did not compile')
        receipt=_execute(budget,result,output,[str(exe)],'clock-controls')
        text=(output/'clock-controls.log').read_text();checks=re.findall(r'^(PASS|FAIL) (.+)$',text,re.M)
        summary=re.findall(r'^SUMMARY checks=(\d+) failed=(\d+)$',text,re.M)
        if len(summary)!=1 or len(checks)!=len(ARITHMETIC_CONTRACTS) or {name for _,name in checks}!=ARITHMETIC_CONTRACTS or int(summary[0][0])!=len(checks) or int(summary[0][1])!=sum(s=='FAIL' for s,_ in checks):
            raise ValueError('incomplete native clock arithmetic contracts')
        if receipt['returncode'] not in (0,1) or bool(receipt['returncode'])!=any(s=='FAIL' for s,_ in checks):raise ValueError('clock status/exit mismatch')
        counts=collections.Counter(s for s,_ in checks)
        parser=output/'runtime-parser-probe'
        command=[compiler,'-std=c++17','-O2','-g','-I',str(REPO/'srcMoon'),str(REPO/'srcMoon/test/m3_runtime_rates_probe.cpp'),
            *(str(REPO/name) for name in ('srcMoon/MoonInput.cpp','srcMoon/species/ElectronImpact.cpp','srcMoon/species/Photochemistry.cpp','srcMoon/species/ChargeExchange.cpp')),
            *(str(p) for p in shared_sources(REPO,'runtime')),'-o',str(parser)]
        if _execute(budget,result,output,command,'compile-parser')['returncode']!=0:raise RuntimeError('native production parser/installer did not compile')
        variants=[]
        for order in ('first','second'):
            base=(REPO/f'srcMoon/examples/m3_tracing_{order}_order.in').read_text()
            variants += [(order+'-default',re.sub(r'(?m)^tracing_clock_mode[^\n]*\n','',base),0,'independent-history'),
                (order+'-explicit-independent',change(base,'tracing_clock_mode','independent-history'),0,'independent-history'),
                (order+'-explicit-common',change(base,'tracing_clock_mode','common-time'),0,'common-time'),
                (order+'-unknown-clock',change(base,'tracing_clock_mode','assumed-time'),3,None),
                (order+'-common-m2-refused',change(change(base,'tracing_clock_mode','common-time'),'m2_environment','enabled'),3,None),
                (order+'-duplicate-clock',change(base,'tracing_clock_mode','common-time').replace('#end application','tracing_clock_mode = independent-history\n#end application'),3,None)]
        legacy=(REPO/'srcMoon/examples/m1_regression.in').read_text()
        variants += [('legacy-absent-clock',legacy,0,'independent-history'),
            ('legacy-common-refused',change(legacy,'tracing_clock_mode','common-time'),3,None)]
        parsed=[]
        for name,text,expected,clock in variants:
            directory=output/name;directory.mkdir();application=directory/'application.in';application.write_text(text)
            receipt=_execute(budget,result,output,[str(parser),str(application)],name+'/native',directory)
            transcript=(directory/'native.log').read_text();values=dict(line.split('=',1) for line in transcript.splitlines() if '=' in line)
            good=receipt['returncode']==expected and (values.get('tracing_clock_mode')==clock if expected==0 else values.get('parse_preserved')=='1')
            status='PASS' if good else 'FAIL';counts[status]+=1;parsed.append({'name':name,'status':status,'expected_exit':expected,'actual_exit':receipt['returncode']})
        if any(sha256(REPO/name)!=expected for name,expected in manifest['source_files'].items()):
            raise ValueError('production source changed during local clock execution')
        manifest.update(execution='actual production clock and parser/installer executables',executable={'path':str(exe),'sha256':sha256(exe)},parser_executable={'path':str(parser),'sha256':sha256(parser)},mpi_layouts=[],omp_threads=1)
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',counts={s:counts[s] for s in ('PASS','FAIL','ERROR','SKIPPED')},
            arithmetic_checks=[{'name':name,'status':status} for status,name in checks],configuration_checks=parsed)
    except Exception as error:_infrastructure_error(output,result,error)
    return _finish(output,result,manifest,started)


def run(output,producer_manifest,wall_seconds=600):
    """Actual normal Moon/PIC hooks on the same fresh generated producer."""
    started=time.monotonic();output,result,manifest,budget=_start(output,wall_seconds,
        'supporting actual Moon common-time hooks/source-boundary/mover; analytical chronology, no observations or physical reaction activation')
    try:
        if any(shutil.which(tool) is None for tool in ('mpicxx','mpiexec')):
            result.update(status='SKIPPED',reason='native compiler/MPI unavailable');return _finish(output,result,manifest,started)
        import matplotlib
        producer_manifest=Path(producer_manifest).resolve()
        if not producer_manifest.is_file():result.update(status='SKIPPED',reason='successful all-species producer unavailable');return _finish(output,result,manifest,started)
        producer=json.loads(producer_manifest.read_text());predecessor=producer_manifest.with_name('result.json')
        if json.loads(predecessor.read_text()).get('status')!='PASS':
            result.update(status='SKIPPED',reason='all-species real PIC backend predecessor is not PASS');return _finish(output,result,manifest,started)
        if Path(producer['source_root']).resolve()!=REPO:raise ValueError('producer belongs to another actual checkout')
        library=Path(producer['library']['path']);build=library.parent
        if sha256(library)!=producer['library']['sha256']:raise ValueError('actual producer library bytes changed')
        for path,expected in producer['generated_snapshots'].items():
            relative=Path(path).relative_to(producer_manifest.parent/'generated')
            if sha256(path)!=expected or sha256(build/relative)!=expected:raise ValueError('retained/current generated ABI changed')
        if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):
            raise ValueError('common clock/native sources differ from compiled producer')
        manifest.update(producer= {'path':str(producer_manifest),'sha256':sha256(producer_manifest),'result':str(predecessor),'result_sha256':sha256(predecessor)},
            library=producer['library'],production_source_hashes=producer['production_source_hashes'],generated_snapshots=producer['generated_snapshots'],
            mpi_layouts=[1,8],orders=[1,2],omp_threads=1)
        source=Path(__file__).with_name('m3_native_transport_clock_probe.cpp').resolve();executable=output/'native-common-clock-probe'
        manifest['test_sources']={str(p):sha256(p) for p in (source,Path(__file__).resolve(),Path(__file__).with_name('M3_COMMON_TRANSPORT_CLOCK_CONTRACT.md'))}
        for path in manifest['test_sources']:shutil.copy2(path,output/Path(path).name)
        command=_probe_compile_command(build,executable);command.insert(1,'-DMPI_ON')
        command[command.index(str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp'))]=str(source)
        if _execute(budget,result,output,command,'compile-native')['returncode']!=0:raise RuntimeError('native common-clock thin probe did not compile')
        counts=collections.Counter();layouts=[];births=[];daughters=[]
        required={
            'explicit_common_mode_and_real_tracer_selected','common_storage_origin_and_frozen_coefficients_initialized',
            'unequal_actual_dt_refused_before_physical_step','changed_external_epoch_refused_before_begin',
            'first_common_step_begins_with_immutable_relative_endpoints','nested_begin_cannot_replace_active_step',
            'scoped_offset_source_uses_actual_advancing_transport_birth','first_normal_mpi_continuations_complete',
            'bad_endpoint_commit_refused_without_clock_publication','first_commit_accepts_only_completed_endpoint_population',
            'variable_common_duration_refused_without_advancing_clock','second_common_step_begins_after_accepted_transport_time',
            'recycled_later_source_uses_fresh_clock_target_and_face_birth','second_normal_population_completes_common_endpoint',
            'third_step_preserves_surviving_relative_chronology','all_packets_finish_common_physical_endpoint',
            'external_epoch_and_prescribed_fields_remain_frozen','no_unfinished_physical_continuation','no_native_particle_slot_leak'}
        required.add('real_conservative_daughters_inherit_common_target_and_finish_unused_interval')
        for ranks in (1,8):
            for order in (1,2):
                directory=output/f'np{ranks}-order{order}';directory.mkdir();runtime=directory/'runtime';runtime.mkdir();(runtime/'PT').mkdir()
                receipt=_execute(budget,result,output,['mpiexec','-n',str(ranks),str(executable),str(order),str(directory/'clock')],
                    str(directory.relative_to(output))+'/native',runtime,ranks)
                transcript=(directory/'native.log').read_text();summary=re.findall(r'M3 common transport clock order=(\d+) checks=(\d+) status=(PASS|FAIL|ERROR)',transcript)
                with (directory/'clock.checks.csv').open(newline='') as stream:checks=list(csv.DictReader(stream))
                if len(summary)!=1 or int(summary[0][0])!=order or int(summary[0][1])!=len(required) or summary[0][2]=='ERROR':
                    raise ValueError('native common-clock summary incomplete/error')
                if len(checks)!=len(required) or {r['check'] for r in checks}!=required or any(r['status'] not in ('PASS','FAIL') for r in checks):
                    raise ValueError('native clock contracts missing/duplicate/malformed')
                failed=any(r['status']=='FAIL' for r in checks)
                if receipt['returncode'] not in (0,1) or (receipt['returncode']==1)!=failed or (summary[0][2]=='FAIL')!=failed:
                    raise ValueError('native common-clock exit/check disagreement')
                count_before=len(births)
                for rank in range(ranks):
                    with (directory/f'clock.rank{rank}.births.csv').open(newline='') as stream:
                        for row in csv.DictReader(stream):births.append(dict(mpi_ranks=ranks,order=order,**row))
                daughter_before=len(daughters)
                for rank in range(ranks):
                    with (directory/f'clock.rank{rank}.daughters.csv').open(newline='') as stream:
                        reader=csv.DictReader(stream)
                        if reader.fieldnames!=['catalog_state','compiled_species','epoch_s','remaining_s','represented_number','event_epoch_s','event_remaining_s']:
                            raise ValueError('common daughter state diagnostics malformed')
                        for row in reader:daughters.append(dict(mpi_ranks=ranks,order=order,**row))
                if not failed and len(births)-count_before!=3:raise ValueError('successful native chronology omitted declared source births')
                if not failed and len(daughters)-daughter_before!=2:raise ValueError('successful common target composition omitted actual daughter evidence')
                if not failed:
                    actual_daughters=daughters[daughter_before:]
                    if {row['catalog_state'] for row in actual_daughters}!={'2','3'} or any(
                        not all(math.isfinite(float(value)) for value in row.values()) or
                        float(row['epoch_s'])!=1004 or float(row['remaining_s'])!=0 or
                        abs(float(row['represented_number'])-1)>=1e-13 for row in actual_daughters):
                        raise ValueError('successful common daughter CSV contradicts accepted target/state/number controls')
                counts.update(r['status'] for r in checks);layouts.append({'mpi_ranks':ranks,'order':order,'status':'FAIL' if failed else 'PASS','checks':len(checks)})
        with (output/'birth_comparison.csv').open('w',newline='') as stream:
            fields=['mpi_ranks','order','source_tag','step','relative_birth_s','expected_birth_et_s','actual_birth_et_s','remaining_s','source_owner']
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(births)
        with (output/'daughter_comparison.csv').open('w',newline='') as stream:
            fields=['mpi_ranks','order','catalog_state','compiled_species','epoch_s','remaining_s','represented_number','event_epoch_s','event_remaining_s']
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(daughters)
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        figure,axes=plt.subplots(1,2,figsize=(10,4))
        for ranks in (1,8):
            for order in (1,2):
                rows=[r for r in births if r['mpi_ranks']==ranks and r['order']==order]
                expected=[float(r['expected_birth_et_s'])-1000 for r in rows];actual=[float(r['actual_birth_et_s'])-1000 for r in rows]
                axes[0].scatter(expected,actual,label=f'{ranks}MPI/order{order}')
                axes[1].scatter(expected,[a-b for a,b in zip(actual,expected)])
        axes[0].plot([0,2],[0,2],color='black');axes[0].legend(framealpha=1)
        axes[0].set_ylabel('Actual source birth minus origin [s]');axes[1].set_ylabel('Actual minus analytical birth [s]')
        for axis in axes:axis.set_xlabel('Analytical source birth minus origin [s]');axis.grid(True)
        figure.suptitle('Common-time analytical wiring controls — no observations');figure.tight_layout()
        for suffix in ('png','eps'):figure.savefig(output/('birth_comparison.'+suffix))
        plt.close(figure)
        if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):
            raise ValueError('native production source changed during clock execution')
        if any(sha256(path)!=expected for path,expected in manifest['test_sources'].items()):
            raise ValueError('native clock test/contract source changed during execution')
        manifest.update(execution='actual normal Moon common Begin/source-boundary/mover/Commit hooks',executable={'path':str(executable),'sha256':sha256(executable)})
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',counts={s:counts[s] for s in ('PASS','FAIL','ERROR','SKIPPED')},layouts=layouts)
    except ImportError as error:result.update(status='SKIPPED',reason='required plotting dependency unavailable: '+str(error))
    except Exception as error:_infrastructure_error(output,result,error)
    return _finish(output,result,manifest,started)


def run_startup(output,wall_seconds=600,build_evidence=None):
    """Two actual normal Na-only startup/step runs, one per tracing order.

    The parent must supply the REAL successful ordinary-SPICE restore receipt.
    This function never configures/builds an application or manufactures that
    provenance from source existence. No evidence means SKIPPED, not PASS.
    """
    started=time.monotonic();output,result,manifest,budget=_start(output,wall_seconds,
        'supporting full normal amps_init plus two actual amps_time_step calls at8MPI, both orders; no observations or physical kinetic qualification')
    try:
        if build_evidence is None:
            result.update(status='SKIPPED',reason='authentic successful normal SPICE restore build evidence unavailable')
            return _finish(output,result,manifest,started)
        if isinstance(build_evidence,(str,Path)):
            evidence_path=Path(build_evidence).resolve();evidence=json.loads(evidence_path.read_text())
            manifest['build_evidence_file']={'path':str(evidence_path),'sha256':sha256(evidence_path)}
        else:evidence=build_evidence
        if not isinstance(evidence,dict) or evidence.get('status')!='PASS':
            raise ValueError('normal restore build evidence is not a successful actual build')
        if Path(evidence['source_root']).resolve()!=REPO:raise ValueError('normal restore belongs to another checkout')
        expected_sources=evidence['production_source_hashes']
        if native_sources(expected_sources)!=native_sources(production_source_identity()):raise ValueError('normal restored source does not match current native sources')
        library=Path(evidence['library']['path']).resolve();build=library.parent
        if library!=REPO/'build/libAMPS.a' or sha256(library)!=evidence['library']['sha256']:raise ValueError('normal restored library bytes/path changed')
        executable_proof=Path(evidence['executable']['path']).resolve()
        if sha256(executable_proof)!=evidence['executable']['sha256']:raise ValueError('normal restored executable bytes changed')
        build_log=Path(evidence['build_log']['path']).resolve()
        if sha256(build_log)!=evidence['build_log']['sha256']:raise ValueError('authentic normal build log changed')
        if not evidence.get('build_commands'):raise ValueError('normal successful build did not retain commands')
        log=build_log.read_text();statuses=re.findall(r'^EXIT: (-?\d+)$',log,re.M)
        if len(statuses)<2 or any(int(status)!=0 for status in statuses):raise ValueError('normal build log does not establish successful configuration and compilation')
        snapshots=evidence.get('generated_snapshots')
        if not isinstance(snapshots,dict) or not snapshots:raise ValueError('normal build lacks generated configuration identity')
        required_generated={'pic.h','picGlobal.dfn','picSpeciesMacro.dfn','picParticleDataMacro.h','rnd.h'}
        if not required_generated.issubset({Path(name).name for name in snapshots}):raise ValueError('normal build generated identity incomplete')
        for path,expected in snapshots.items():
            if sha256(path)!=expected:raise ValueError('normal generated configuration changed')
        import matplotlib
        manifest.update(normal_build_evidence=evidence,production_source_hashes=expected_sources,
            library=evidence['library'],mpi_layouts=[8],orders=[1,2],omp_threads=1)
        source=Path(__file__).with_name('m3_normal_transport_clock_startup_probe.cpp').resolve();binary=output/'normal-common-clock-startup'
        tests=[source,Path(__file__).resolve(),Path(__file__).with_name('M3_COMMON_TRANSPORT_CLOCK_CONTRACT.md').resolve()]
        manifest['test_sources']={str(path):sha256(path) for path in tests}
        for path in tests:shutil.copy2(path,output/path.name)
        command=_probe_compile_command(build,binary);command.insert(1,'-DMPI_ON')
        command[command.index(str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp'))]=str(source)
        if _execute(budget,result,output,command,'compile-startup')['returncode']!=0:raise RuntimeError('normal startup thin driver did not compile')
        required={
            'canonical_single_na_common_mode_after_normal_initializer','normal_spice_utc_origin_matches_explicit_transport_origin',
            'normal_startup_has_no_accepted_or_provisional_transport_step','normal_common_dt_is_positive_and_two_steps_do_not_cross_output_window',
            'first_actual_time_step_commits_full_common_duration','second_actual_time_step_commits_full_common_duration',
            'first_normal_sampling_precedes_actual_source_callbacks','second_normal_sampling_precedes_actual_source_callbacks',
            'first_actual_population_finishes_common_endpoint','second_actual_population_finishes_common_endpoint',
            'first_external_geometry_coefficient_epoch_remains_frozen','second_external_geometry_coefficient_epoch_remains_frozen',
            'first_normal_source_cohort_recorded','second_normal_source_cohort_recorded',
            'actual_normal_source_callbacks_were_exercised','actual_furnished_spice_dependencies_recorded'}
        counts=collections.Counter();layouts=[];external={};rows=[]
        for order in ('first','second'):
            directory=output/(order+'-order');directory.mkdir();runtime=directory/'runtime';runtime.mkdir();(runtime/'PT/plots').mkdir(parents=True)
            example=REPO/f'srcMoon/examples/m3_common_time_{order}_order.in'
            text=example.read_text()
            for key,suffix in (('surface_cea_file','cea'),('surface_tecplot_file','dat')):
                text=change(text,key,str(directory/('sphere_surface.'+suffix)))
            application=directory/'application.in';application.write_text(text)
            manifest.setdefault('applications',{})[order]={'path':str(application),'sha256':sha256(application),'template':str(example),'template_sha256':sha256(example)}
            receipt=_execute(budget,result,output,['mpiexec','-n','8',str(binary),str(application),str(directory/'startup')],
                order+'-order/native',runtime,8)
            transcript=(directory/'native.log').read_text();summary=re.findall(r'M3 normal common startup checks=(\d+) status=(PASS|FAIL|ERROR)',transcript)
            with (directory/'startup.checks.csv').open(newline='') as stream:checks=list(csv.DictReader(stream))
            if len(summary)!=1 or summary[0][1]=='ERROR' or int(summary[0][0])!=len(required) or len(checks)!=len(required) or {row['check'] for row in checks}!=required:
                raise ValueError('normal startup clock named controls incomplete/error')
            if any(row['status'] not in ('PASS','FAIL') for row in checks):raise ValueError('normal startup check status malformed')
            failed=any(row['status']=='FAIL' for row in checks)
            if receipt['returncode'] not in (0,1) or (receipt['returncode']==1)!=failed or (summary[0][1]=='FAIL')!=failed:raise ValueError('normal startup exit/criterion mismatch')
            with (directory/'startup.steps.csv').open(newline='') as stream:steps=list(csv.DictReader(stream))
            if not failed and (len(steps)!=2 or [row['step'] for row in steps]!=['1','2']):raise ValueError('normal startup step chronology diagnostics incomplete')
            for row in steps:rows.append(dict(order=order,**row))
            dependencies=(directory/'startup.spice_files.txt').read_text().splitlines()
            if not dependencies or len(dependencies)!=len(set(dependencies)):raise ValueError('actual furnished SPICE dependency paths incomplete/duplicate')
            for name in dependencies:
                path=Path(name);path=path if path.is_absolute() else runtime/path
                path=path.resolve();external[str(path)]=sha256(path)
            counts.update(row['status'] for row in checks);layouts.append({'order':order,'mpi_ranks':8,'checks':len(checks),'status':'FAIL' if failed else 'PASS'})
        if native_sources(expected_sources)!=native_sources(production_source_identity()) or sha256(library)!=evidence['library']['sha256']:
            raise ValueError('actual normal source/archive identity changed during startup execution')
        if any(sha256(path)!=expected for path,expected in manifest['test_sources'].items()):raise ValueError('startup test/contract changed during execution')
        for path,expected in snapshots.items():
            if sha256(path)!=expected:raise ValueError('normal generated headers changed during startup execution')
        with (output/'normal_startup_chronology.csv').open('w',newline='') as stream:
            fields=['order','step','relative_elapsed_s','expected_relative_elapsed_s','projected_endpoint_et_s','external_et_s','sampling_counter','emission_count','gas_count']
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        figure,axes=plt.subplots(1,2,figsize=(10,4))
        for order in ('first','second'):
            selected=[row for row in rows if row['order']==order]
            steps=[int(row['step']) for row in selected]
            expected=[float(row['expected_relative_elapsed_s']) for row in selected]
            actual=[float(row['relative_elapsed_s']) for row in selected]
            axes[0].plot(steps,expected,label=order+' analytical n*dt')
            axes[0].scatter(steps,actual,label=order+' actual normal clock')
            axes[1].scatter(steps,[a-b for a,b in zip(actual,expected)],label=order)
        axes[0].set_ylabel('Accepted relative transport time [s]');axes[1].set_ylabel('Actual minus analytical n*dt [s]')
        for axis in axes:axis.set_xlabel('Actual normal particle step');axis.grid(True);axis.legend(framealpha=1)
        figure.suptitle('Normal startup chronology — analytical wiring, no observations');figure.tight_layout()
        for suffix in ('png','eps'):figure.savefig(output/('normal_startup_chronology.'+suffix))
        plt.close(figure)
        manifest.update(execution='actual normal Na-only SPICE startup and two normal particle steps per tracer order',
            executable={'path':str(binary),'sha256':sha256(binary)},external_data_sha256=external,uses_spice_kernels=True,
            random_seed='each complete input explicitly declares731993; normal application installs it')
        result.update(status='FAIL' if counts['FAIL'] else 'PASS',counts={name:counts[name] for name in ('PASS','FAIL','ERROR','SKIPPED')},layouts=layouts)
    except ImportError as error:result.update(status='SKIPPED',reason='required native startup dependency unavailable: '+str(error))
    except Exception as error:_infrastructure_error(output,result,error)
    return _finish(output,result,manifest,started)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--local-only',action='store_true');p.add_argument('--normal-startup-only',action='store_true');p.add_argument('--build-evidence',type=Path);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--producer-manifest',type=Path);p.add_argument('--wall-seconds',type=float,default=600)
    a=p.parse_args()
    if not 0<a.wall_seconds<=3600:p.error('wall seconds must be in(0,3600]')
    if a.local_only and a.normal_startup_only:p.error('choose only one local/startup mode')
    if not a.local_only and not a.normal_startup_only and a.producer_manifest is None:p.error('linked mode requires --producer-manifest')
    value=run_local(a.output_dir,a.wall_seconds) if a.local_only else run_startup(a.output_dir,a.wall_seconds,a.build_evidence) if a.normal_startup_only else run(a.output_dir,a.producer_manifest,a.wall_seconds)
    print(json.dumps({k:value[k] for k in ('status','counts','reason') if k in value},indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
