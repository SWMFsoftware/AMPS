#!/usr/bin/env python3
"""Compile and exercise the maintained M3 parser/immutable registry installer.

These are input/rate-service controls, not linked particle events or measured
plasma. Every input/output/log is retained. Frozen raw/processed science files
are read only; deliberately damaged copies are rejection fixtures under the
requested output directory and are never observation/reference products.
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

from test_m3_runtime_rates import ROOT, change
from m2_initial_state import assignments
from chemical_database_build import shared_sources, shared_dependencies

DATA=Path('/data/vtenishe/moon_validation_data/photochemistry/processed/phidrates')
EXAMPLE=ROOT/'srcMoon/examples/m3_reaction_rates.in'
_research=assignments(EXAMPLE)['charge_exchange_research_directory']
_research_path=Path(_research['value'])
# The explicit receipt owns this external source path. Resolving it against a
# temporary checkout's test_output would silently lose the real publications
# in isolated builds; include-owner-relative paths follow the native parser.
RESEARCH=(_research_path if _research_path.is_absolute() else
    Path(_research['file']).parent/_research_path).resolve()


def run_checks(output):
    """Return actual evaluated controls, preserving ERROR versus SKIPPED."""
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    result={'schema':'srcMoon-M3-runtime-configuration-v1','status':'ERROR',
        'scope':'production parser/installer controls only; no U/I particle qualification',
        'checks':[],'observational_scoring':False}
    source_names=('MoonInput.cpp','species/ElectronImpact.cpp',
                  'species/Photochemistry.cpp','species/ChargeExchange.cpp')
    dependencies=[DATA/name for name in ('reaction_channels.csv','photo_rates_1au.csv',
                                         'thresholds.csv','cross_sections.csv')]
    dependencies += [RESEARCH/name for name in ('ornl_redbook_v1_osti.pdf','nist_morgan_1985_retry2.pdf')]
    package=DATA.parents[1]
    metadata=[package/name for name in ('README.md','provenance.json','qa_report.json')]
    dependencies += metadata
    dependencies.append(EXAMPLE)
    missing=[str(p) for p in dependencies if not p.is_file()]
    compiler=shutil.which('g++')
    metadata_issue=None;metadata_error=False
    if not missing:
        try:
            # This package's NATIVE readiness vocabulary is extraction-QA
            # PASS, not a fabricated READY key. It supplies rate tables, not
            # complete daughter kinetics or a calibration/hold-out selection.
            metadata[0].read_text()
            provenance=json.loads(metadata[1].read_text())
            qa=json.loads(metadata[2].read_text())
            if qa.get('status')!='PASS':
                metadata_issue='native photochemistry extraction QA is not PASS'
            else:
                for item in provenance['outputs']:
                    path=package/item['path']
                    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
                        raise ValueError('frozen processed identity mismatch: '+str(path))
            result['data_preflight']={'native_qa_status':qa.get('status'),
                'scope':'authoritative rate extraction only; full daughter kinetics unavailable',
                'solar_selection':'separate native quiet/active optical-depth-zero endpoints',
                'observational_training_holdout':'not applicable to analytical rate/input verification',
                'processed_identities_checked':len(provenance['outputs'])}
        except (OSError,ValueError,KeyError,TypeError) as error:
            metadata_issue=str(error);metadata_error=True
    if missing or not compiler:
        result.update(status='SKIPPED',reason='native compiler or required frozen source bytes unavailable',missing=missing)
    elif metadata_issue:
        result.update(status='ERROR' if metadata_error else 'SKIPPED',reason=metadata_issue)
    else:
        executable=output/'runtime-probe'
        command=[compiler,'-std=c++17','-O2','-g','-Wall','-Wextra','-pedantic',
            '-I',str(ROOT/'srcMoon'),str(ROOT/'srcMoon/test/m3_runtime_rates_probe.cpp'),
            *(str(ROOT/'srcMoon'/name) for name in source_names),
            *(str(p) for p in shared_sources(ROOT,'runtime')),'-o',str(executable)]
        process=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
        (output/'compile.stdout.log').write_text(process.stdout)
        (output/'compile.stderr.log').write_text(process.stderr)
        result['compile']={'command':command,'returncode':process.returncode}
        if process.returncode:
            result['reason']='production runtime adapter compilation failed'
        else:
            base=(ROOT/'srcMoon/examples/m1_regression.in').read_text()
            photon=change(change(change(change(base,'chemistry_mode','rate-services'),
                'photochemistry_rate_model','phidrates-v1'),
                'photochemistry_data_directory',DATA),'photochemistry_solar_case','quiet')
            both=change(change(photon,'charge_exchange_rate_model','total-capture-v1'),
                'charge_exchange_research_directory',RESEARCH)
            # These copies intentionally fail identity checks. They are not
            # replacements for native observations or edited science values.
            bad=output/'damaged-source-rejection-fixture';bad.mkdir()
            for source in dependencies[:4]:shutil.copyfile(source,bad/source.name)
            with (bad/'photo_rates_1au.csv').open('a') as stream:stream.write('\n')
            variants=[
                ('legacy',base,0,None),('quiet-photons',photon,0,None),
                ('active-photons',change(photon,'photochemistry_solar_case','active'),0,None),
                ('all-rate-registries',both,0,None),
                ('registry-install-rollback',both,0,'rollback-m3'),
                ('unknown-chemistry-mode',change(photon,'chemistry_mode','conservative-unqualified'),3,None),
                ('missing-explicit-mode',re.sub(r'(?m)^chemistry_mode[^\n]*\n','',photon),3,None),
                ('unknown-solar-spectrum',change(photon,'photochemistry_solar_case','intermediate'),3,None),
                ('missing-solar-case',re.sub(r'(?m)^photochemistry_solar_case[^\n]*\n','',photon),3,None),
                ('disabled-photon-settings',change(photon,'photochemistry_rate_model','legacy'),3,None),
                ('damaged-frozen-input',change(photon,'photochemistry_data_directory',bad),4,None),
                ('missing-data-directory',change(photon,'photochemistry_data_directory',output/'absent'),3,None),
                ('unknown-capture-model',change(both,'charge_exchange_rate_model','scaled-argon-assumed'),3,None),
                ('disabled-capture-settings',change(both,'charge_exchange_rate_model','disabled'),3,None),
            ]
            # These dimensioned, exactly representable test values distinguish
            # swapped fields from correct parsing. They are parser fixtures,
            # not recommended cell-crossing speeds or physical species weights.
            speed_key = 'hydrogen_characteristic_speed_m_s'
            ratio_key = 'hydrogen_to_sodium_reference_weight_ratio'
            hydrogen = change(change(base,speed_key,'1.5625e4'),ratio_key,'0.03125')
            variants.extend([
                ('hydrogen-complete',hydrogen,0,None),
                ('hydrogen-programmatic-incomplete',hydrogen,0,'reject-incomplete-hydrogen'),
                ('hydrogen-speed-alone',change(base,speed_key,'15625'),3,None),
                ('hydrogen-ratio-alone',change(base,ratio_key,'0.03125'),3,None),
            ])
            hydrogen_diagnostics = {'hydrogen-speed-alone':ratio_key, 'hydrogen-ratio-alone':speed_key}
            # Explicit zero cannot mean "omitted". A partial group, a signed
            # value, a non-finite value and a trailing unit are parser errors.
            for role, key in (('speed',speed_key),('ratio',ratio_key)):
                for label, invalid in (('zero','0'),('negative','-1'),('infinite','inf'),
                                       ('nan','nan'),('trailing-unit','4 m/s')):
                    name = 'hydrogen-'+role+'-'+label
                    variants.append((name,change(hydrogen,key,invalid),3,None))
                    hydrogen_diagnostics[name] = key
            for order in ('first','second'):
                tracing=(ROOT/f'srcMoon/examples/m3_tracing_{order}_order.in').read_text()
                variants += [(order+'-order-tracing',tracing,0,None),
                    (order+'-common-time-selection',change(tracing,'tracing_clock_mode','common-time'),0,None),
                    (order+'-absent-clock-default',re.sub(r'(?m)^tracing_clock_mode[^\n]*\n','',tracing),0,None),
                    (order+'-unknown-clock',change(tracing,'tracing_clock_mode','assumed-clock'),3,None),
                    (order+'-m2-parse-refusal',change(tracing,'m2_environment','enabled'),3,None),
                    (order+'-m2-install-refusal',tracing,0,'reject-dynamic-m2'),
                    (order+'-missing-tolerance',re.sub(r'(?m)^tracing_relative_tolerance[^\n]*\n','',tracing),3,None),
                    (order+'-negative-step',change(tracing,'tracing_minimum_substep_s',-1),3,None),
                    (order+'-bad-step-range',change(tracing,'tracing_minimum_substep_s',100),3,None),
                    (order+'-bad-field-shape',change(tracing,'tracing_magnetic_field_t','1,2'),3,None),
                    (order+'-unknown-frame',change(tracing,'tracing_field_frame','assumed-GSE'),3,None),
                    (order+'-evolving-provider-unavailable',change(tracing,'simulation_environment_mode','evolving'),3,None),
                    (order+'-inactive-controls',change(tracing,'particle_tracing_mover','legacy'),3,None)]
            references={'quiet':{'He':5.64e-8,'Na':7.26e-6,'Ne':1.84e-7,'Ar':3.42e-7},
                'active':{'He':1.68e-7,'Na':7.91e-6,'Ne':5.87e-7,'Ar':7.74e-7}}
            for name,text,expected,mode in variants:
                folder=output/name;folder.mkdir()
                application=folder/'application.in';application.write_text(text)
                invocation=[str(executable),str(application)]+([mode] if mode else [])
                native=subprocess.run(invocation,cwd=folder,text=True,capture_output=True)
                (folder/'stdout.log').write_text(native.stdout);(folder/'stderr.log').write_text(native.stderr)
                values=dict(line.split('=',1) for line in native.stdout.splitlines() if '=' in line)
                errors=[]
                if native.returncode!=expected:errors.append('unexpected process status')
                elif expected==3 and values.get('parse_preserved')!='1':errors.append('failed parse mutated output')
                elif expected==4 and values.get('install_preserved')!='1':errors.append('failed install published a registry')
                elif expected==0:
                    if name in ('legacy','hydrogen-complete','hydrogen-programmatic-incomplete'):
                        pair = tuple(float(values.get(key,'nan')) for key in (speed_key,ratio_key))
                        if pair != ((0,0) if name=='legacy' else (15625,0.03125)):
                            errors.append('installed hydrogen speed or reference-weight ratio differs from input')
                    if mode=='reject-incomplete-hydrogen' and values.get('incomplete_hydrogen_install_preserved')!='1':
                        errors.append('incomplete direct install was not refused before the successful retry')
                    if name in ('quiet-photons','active-photons','all-rate-registries','registry-install-rollback'):
                        solar='active' if name=='active-photons' else 'quiet'
                        if values.get('photon_catalogue_rows')!='17' or values.get('photon_reactive_rows')!='15':
                            errors.append('incomplete reaction family installed')
                        for atom,rate in references[solar].items():
                            for distance in (1,2):
                                actual=float(values.get(f'{atom}_photon_rate_{distance}au_s','nan'))
                                if not math.isclose(actual,rate/distance**2,rel_tol=1e-13):
                                    errors.append(atom+' endpoint/distance wiring')
                    if mode=='rollback-m3' and values.get('failed_m3_install_preserved')!='1':errors.append('rollback was not evaluated')
                    if mode=='reject-dynamic-m2' and values.get('incompatible_m2_install_preserved')!='1':errors.append('programmatic temporal-coupling refusal was not evaluated')
                    if name.endswith('-order-tracing'):
                        expected_mover='force-cut-surface-'+name[:-len('-tracing')]
                        if values.get('tracing_mover')!=expected_mover or values.get('tracing_field_frame')!='lso':
                            errors.append('tracing selection/frame mismatch')
                        if float(values.get('tracing_maximum_substep_s','nan'))!=10:errors.append('physical segment limit mismatch')
                    if 'common-time-selection' in name and values.get('tracing_clock_mode')!='common-time':
                        errors.append('common clock selection was not preserved')
                    if ('absent-clock-default' in name or name.endswith('-order-tracing')) and values.get('tracing_clock_mode')!='independent-history':
                        errors.append('existing independent clock default changed')
                if name.endswith('-m2-parse-refusal') and 'm2_environment=disabled' not in native.stderr:
                    errors.append('wrong refusal: unsupported physical-time coupling was not identified')
                if name in hydrogen_diagnostics and hydrogen_diagnostics[name] not in native.stderr:
                    errors.append('expected hydrogen setting was absent from the parser diagnostic')
                status='ERROR' if native.returncode<0 else 'FAIL' if errors else 'PASS'
                result['checks'].append({'name':name,'status':status,'command':invocation,
                    'expected_returncode':expected,'returncode':native.returncode,'errors':errors,
                    'log':str(folder/'stdout.log'),'input_sha256':hashlib.sha256(application.read_bytes()).hexdigest()})
            result['status']='FAIL' if any(r['status']=='FAIL' for r in result['checks']) else 'ERROR' if any(r['status']=='ERROR' for r in result['checks']) else 'PASS'
    result['counts']={s:sum(r['status']==s for r in result['checks']) for s in ('PASS','FAIL','SKIPPED','ERROR')}
    result['wall_seconds']=time.monotonic()-started
    result['source_files']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [*(ROOT/'srcMoon'/n for n in source_names),
                  *(ROOT/'srcMoon'/Path(n).with_suffix('.h') for n in source_names),
                  ROOT/'srcMoon/test/m3_runtime_rates_probe.cpp',Path(__file__).resolve(),
                  *shared_dependencies(ROOT,'runtime')]}
    result['data_files']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies if p.is_file()}
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();value=run_checks(args.output_dir)
    print(json.dumps({k:v for k,v in value.items() if k in ('status','counts','reason')},indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
