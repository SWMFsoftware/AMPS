#!/usr/bin/env python3
"""Build an explicit eight-species rate fixture and check the real callback.

This is a supporting U16 wiring runner. It enables neither electron-impact
particle reactions nor an ion-transport model. Normal production input and
regression references are never edited. Use an isolated checkout: like other
native AMPS builders this changes the selected checkout's generated tree.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import sys
import uuid

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,replace_once,sha256
from integration_common import REPO,manifest
from m2_initial_state import assignments
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command


def rate_manifest(executable,command,application,mpi_ranks):
    """Describe actual rate-only inputs without claiming SPICE consumption.

    This adapter initializes MPI, installs an explicit atomic coefficient
    table and queries scalar callbacks. It furnishes no ephemeris kernels,
    loads no AMR/terrain, advances no particles and scores no observations.
    Generic manifest inventory remains useful build/configuration context,
    but those installed files cannot be called physical runtime dependencies.
    """
    application=Path(application).resolve()
    records=assignments(application)
    if records.get('electron_impact_rate_model',{}).get('value')!='maxwellian-voronov':
        raise ValueError('rate-only manifest requires explicit maxwellian-voronov input')
    if 'electron_impact_coefficients_file' not in records:
        raise ValueError('rate-only manifest requires an explicit coefficient file')
    record=records['electron_impact_coefficients_file']
    coefficients=Path(record['value'])
    if not coefficients.is_absolute():
        coefficients=Path(record['file']).parent/coefficients
    coefficients=coefficients.resolve()
    if not coefficients.is_file():
        raise ValueError('configured rate coefficient file is unavailable: '+str(coefficients))
    # U16 has no generic observational-package dependency. Use the ordinary
    # metadata constructor, then narrow its installed-file inventory explicitly;
    # require_scientific_data=False is reserved there for I15/I20 operators.
    # The generic constructor's application-input branch assumes flat scalar
    # extraction and may inventory unused terrain. This callback consumes no
    # geometry; pass no application there and record the actual native include
    # owners below, preserving included-input support without staging/rewrite.
    identity=manifest('U16',Path(executable),command,None,None,mpi_ranks)
    identity['configuration_preflight_inventory']={
        'role':'installed build/configuration-preflight files; not consumed by the rate-only callback',
        'files':dict(identity.get('external_data',{})),
    }
    identity.update(uses_spice_kernels=False,
        configuration='native rate-only callback; SPICE may be linked but is not initialized or used',
        species_fixture='Na,NA_PLUS,He,HE_PLUS,Ne,NE_PLUS,Ar,AR_PLUS',
        execution_scope='production rate callback only; no particle events or observations',
        scientific_data_preflight_required=False,formal_test_status='SKIPPED',
        formal_I23_status='SKIPPED',observational_scoring=False,
        production_source_hashes=production_source_identity(),
        external_data={str(coefficients):sha256(coefficients)})
    # Included scalar ownership controls relative paths. Preserve every owner
    # and the specific coefficient assignment so the native input can be
    # reproduced without rebasing it against the launcher/root working folder.
    inputs=identity.setdefault('inputs',{})
    inputs[str(application)]=sha256(application)
    for path in {Path(item['file']).resolve() for item in records.values()}:
        inputs[str(path)]=sha256(path)
    identity['coefficient_input']={'configured_value':record['value'],
        'assignment_file':record['file'],'assignment_line':record['line'],
        'resolved_file':str(coefficients),'sha256':sha256(coefficients)}
    return identity


def run(output,input_file,*,mpi_ranks=8,wall_seconds=1800,
        spice_path='/home/vtenishe/SPICE/cspice',spice_kernels='/home/vtenishe/SPICE/Kernels'):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    # Keep clickable native stream paths even when configuration, linking or
    # provenance fails before MPI launch. CampaignBudget replaces stdout.log
    # with actual combined streams once execution starts; an unlaunched marker
    # is never mistaken for runtime/scientific evidence.
    (output/'stdout.log').write_text('Native callback not launched yet; see result.json and configure/build/probe logs.\n')
    (output/'stderr.log').write_text('')
    budget=CampaignBudget(wall_seconds,min(120,wall_seconds*.1))
    result={'status':'ERROR','scope':'supporting U16 production-rate callback; no reaction events',
        'formal_I23_status':'SKIPPED','commands':[],'mpi_ranks':mpi_ranks}
    application='moon_m3_rates_'+uuid.uuid4().hex[:12]
    created=[]
    def execute(command,label,cwd=REPO,*,allow_completed_failure=False):
        execution=budget.run(command,cwd,output/(label+'.log'),output=output)
        result['commands'].append(execution)
        if execution['returncode']!=0 and not (allow_completed_failure and
                execution['returncode']==1 and not execution.get('reason')):
            raise RuntimeError(label+' failed: '+str(execution.get('reason') or
                'child exit '+str(execution['returncode'])+'; see '+str(output/(label+'.log'))))
        return execution
    try:
        text=(REPO/'input/moon.input').read_text()
        text=replace_once(text,r'^(\s*SpeciesList\s*=\s*)Na(\s*!.*)$',
            r'\1Na,NA_PLUS,He,HE_PLUS,Ne,NE_PLUS,Ar,AR_PLUS\2','canonical Na-only species list')
        text=replace_once(text,r'^#block[ \t]+\.\./moon\.input\.pl[ \t]*[^\r\n]*$',
            '#block ../'+application+'.input.pl ! derived rate-only fixture','Moon parser block')
        text=replace_once(text,r'^\s*#include\s+moon\.input\.spicekernels\s*$',
            '#include '+application+'.input.spicekernels','Moon kernel include')
        staged={'.input':text,'.input.pl':(REPO/'input/moon.input.pl').read_text(),
                '.input.spicekernels':(REPO/'input/moon.input.spicekernels').read_text()}
        for suffix,content in staged.items():
            path=REPO/'input'/(application+suffix);path.write_text(content);created.append(path)
            # Preserve the exact derived compile receipt before unique input
            # cleanup; an investigator can inspect the actual built settings.
            snapshot=output/'compile-inputs'/(application+suffix)
            snapshot.parent.mkdir(exist_ok=True);snapshot.write_text(content)
        result['derived_compile_input']={suffix:sha256(REPO/'input'/(application+suffix)) for suffix in staged}
        build=REPO/'build'
        if build.exists():shutil.rmtree(build)
        execute([str(REPO/'Config.pl'),'-application='+application,'-amps-test=off',
            '-spice-path='+str(spice_path),'-spice-kernels='+str(spice_kernels)],'configure')
        execute(['make','-j16','OPT3=-O3 -g'],'build')
        executable=output/'native-rate-callback-probe'
        command=_probe_compile_command(build,executable)
        old=str(Path(__file__).resolve().parent/'common/production_kernel_probe.cpp')
        command[command.index(old)]=str(Path(__file__).resolve().parent/'m3_native_rate_callback_probe.cpp')
        execute(command,'compile-probe')
        # Preserve the exact generated ABI and species map actually linked.
        for relative in ('pic/pic.h','pic/picGlobal.dfn','pic/picSpeciesMacro.dfn','general/rnd.h'):
            destination=output/'generated'/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(build/relative,destination)
        runtime=output/'runtime';runtime.mkdir();(runtime/'PT').mkdir()
        command=['mpiexec','-n',str(mpi_ranks),str(executable),str(Path(input_file).resolve()),str(output/'callback')]
        identity=rate_manifest(executable,command,Path(input_file).resolve(),mpi_ranks)
        (output/'run_manifest.json').write_text(json.dumps(identity,indent=2,sort_keys=True)+'\n')
        shutil.copy2(input_file,output/'application.in')
        native_execution=execute(command,'stdout',runtime,allow_completed_failure=True)
        (output/'stderr.log').write_text('Combined MPI stdout/stderr is preserved in stdout.log.\n')
        rows=[];failures=[]
        # Independent physical coefficient oracle. D05's native values and
        # the high-precision U16 scorer qualify the formula separately; here
        # the expected 10-eV/2.5e6-m^-3 fixture additionally tests wiring.
        reference={'He':(24.6,0,1.75e-14,.18,.35),'Na':(5.1,1,1.01e-13,.275,.23),
                   'Ne':(21.6,1,1.5e-14,.0329,.43),'Ar':(15.8,1,5.99e-14,.136,.26)}
        from integration_common import moon_setting
        density=float(moon_setting(Path(input_file),'electron_density_m3'))
        temperature=float(moon_setting(Path(input_file),'electron_temperature_ev'))
        for rank in range(mpi_ranks):
            with (output/('callback.rank'+str(rank)+'.csv')).open(newline='') as stream:
                actual=list(csv.DictReader(stream))
            if len(actual)!=4 or {row['atom'] for row in actual}!=set(reference):
                raise ValueError('native rank output does not cover all four actual neutral species')
            for row in actual:
                energy,p,a,x,k=reference[row['atom']];u=energy/temperature
                expected=density*a*(1+p*math.sqrt(u))/(x+u)*u**k*math.exp(-u)
                masses=[float(row[key]) for key in ('neutral_mass_kg','ion_mass_kg')]
                neutral_charge=float(row['neutral_charge_c']);ion_charge=float(row['ion_charge_c'])
                if (int(row['rank'])!=rank or row['expected_ion_index']!=row['actual_ion_index'] or
                    not math.isclose(float(row['callback_rate_s']),expected,rel_tol=1e-12,abs_tol=0) or
                    row['callback_rate_s']!=row['scalar_rate_s'] or
                    not all(math.isfinite(mass) and mass>0 for mass in masses) or
                    neutral_charge!=0 or not math.isfinite(ion_charge) or ion_charge<=0):failures.append(row)
                rows.append(row)
        # A completed native evaluated miss returns one. Builds, signals,
        # launcher failures or incomplete/corrupt rows are infrastructure
        # ERROR instead. Require the native rank diagnostic as corroboration;
        # MPI's generic nonzero-exit message alone is not scientific evidence.
        transcript=(output/'stdout.log').read_text()
        native_failed=' status=FAIL' in transcript
        if native_execution['returncode']==1 and not native_failed:
            raise ValueError('nonzero MPI exit lacks a completed native FAIL diagnostic')
        if native_execution['returncode']==0 and native_failed:
            raise ValueError('native success exit contradicts retained rank FAIL diagnostic')
        if native_failed and not failures:
            failures.append({'check':'native ion exclusion/product mapping invariant',
                             'log':str(output/'stdout.log')})
        result.update(status='FAIL' if failures else 'PASS',rows=len(rows),failures=failures,
            executable={'path':str(executable),'sha256':sha256(executable)})
    except (OSError,ValueError,RuntimeError,KeyError,TypeError) as error:
        result.update(status='ERROR',reason=str(error))
        with (output/'stderr.log').open('a') as stream:
            stream.write('Runner infrastructure error: '+str(error)+'\n')
    finally:
        # Remove only the unique temporary inputs and their Config.pl copies.
        # No wildcard cleanup may touch another campaign or user configuration.
        for path in created:
            path.unlink(missing_ok=True)
        for path in REPO.glob(application+'*'):
            if path.is_file():path.unlink()
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--input',required=True,type=Path)
    parser.add_argument('--mpi-ranks',type=int,default=8)
    parser.add_argument('--wall-seconds',type=float,default=1800)
    parser.add_argument('--spice-path',default='/home/vtenishe/SPICE/cspice')
    parser.add_argument('--spice-kernels',default='/home/vtenishe/SPICE/Kernels')
    args=parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error('output must be new or empty')
    if args.mpi_ranks<1:parser.error('positive MPI rank count required')
    result=run(args.output_dir,args.input,mpi_ranks=args.mpi_ranks,wall_seconds=args.wall_seconds,
        spice_path=args.spice_path,spice_kernels=args.spice_kernels)
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['status']=='PASS' else 1 if result['status']=='FAIL' else 2)
