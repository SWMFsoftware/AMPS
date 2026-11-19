#!/usr/bin/env python3
"""Source-bound production composer against independent Decimal mechanics.

REAL frozen D04 rates/means are consumed without editing. Masses/states/frame
bindings are EXPLICIT ANALYTICAL FIXTURES, not a real Na threshold crosswalk,
observations or an admitted production profile. No PIC/AMPS build is invoked.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from decimal import Decimal,localcontext
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import traceback

from chemical_database_build import shared_sources,shared_dependencies
import m3_qualification as q

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent/'physical_reaction_events'
DATA=Path('/data/vtenishe/moon_validation_data/photochemistry')
ARITHMETIC_BOUND=64*sys.float_info.epsilon
GROUPS=[a+'-'+b for a in ('quiet','active') for b in ('rest','moving','affine')]
SUFFIXES=('compose','conservative-plan','one-ion-plus-bath-electron','number-no-timestep-factor',
    'original-unused-interval','fresh-daughter-clock','frozen-basis-distinct-physical-epoch','source-and-policy-identities')
EXPECTED={a+'-source-admission' for a in ('quiet','active')}|{g+'-'+s for g in GROUPS for s in SUFFIXES}|set('''
missing-crosswalk-preserves-provider missing-physical-ion-mass-refused missing-independent-error-budget-refused
wrong-state-threshold-refused actual-catalog-mass-mismatch-refused PIC-electron-substitution-refused
loose-numerical-bound-cannot-admit-physical-gap
endpoint-interpolation-refused missing-caller-admission-refused physical-model-error-budget-not-widened
omitted-solar-endpoint-refused omitted-state-transform-derivative-refused
inconsistent-physical-event-epoch-refused supported-absolute-projection-retains-actual-context
wrong-eV-conversion-refused unsupported-molecular-family-not-renormalized event-in-shadow-refused-no-output
closed-polar-draw-refused-no-output closed-azimuth-draw-refused-no-output missing-photon-ray-refused-no-output
improper-frame-refused-no-output wrong-parent-state-energy-refused-no-output outside-classical-parent-domain-refused
invalid-source-distance-refused missing-daughter-clock-refused finite-disc-or-cosine-not-substituted
point-shadow-hazard-exact-zero same-source-ion-collector-state
'''.split())


def preflight():
    for name in ('README.md','provenance.json','qa_report.json'):
        if not (DATA/name).is_file():raise q.MissingPrerequisite('required real D04 metadata unavailable: '+name)
    (DATA/'README.md').read_text();provenance=q.strict_object(DATA/'provenance.json');qa=q.strict_object(DATA/'qa_report.json')
    if qa.get('status')!='PASS':raise q.MissingPrerequisite('native D04 extraction QA is not PASS')
    hashes={str(DATA/name):q.sha(DATA/name) for name in ('README.md','provenance.json','qa_report.json')}
    for item in provenance['outputs']:
        path=DATA/item['path']
        if not path.is_file():raise q.MissingPrerequisite('required frozen processed source unavailable: '+str(path))
        if q.sha(path)!=item['sha256']:raise ValueError('frozen processed D04 identity changed: '+str(path))
        hashes[str(path)]=item['sha256']
    processed=DATA/'processed/phidrates'
    with (processed/'photo_rates_1au.csv').open(newline='') as stream:
        rates=[r for r in csv.DictReader(stream) if r['species']=='Na']
    with (processed/'thresholds.csv').open(newline='') as stream:
        thresholds=[r for r in csv.DictReader(stream) if r['species']=='Na']
    if len(rates)!=2 or len(thresholds)!=1 or {r['solar_activity'] for r in rates}!={'0.0','1.0'}:
        raise ValueError('independent source Na endpoint/threshold membership changed')
    if any(r['optical_depth']!='0.0' or r['distance_au']!='1.0' for r in rates):raise ValueError('independent source rate convention changed')
    # Match the documented binary64 source-registry storage, then perform
    # ALL independent photon/mechanics calculations in70-digit Decimal.
    # Conversion preserves the actual source parse; it is not fitted to output.
    ev=Decimal('1.602176634e-19')
    source={key:{'mean-J':Decimal.from_float(float(r['excess_energy_eV']))*ev,
                 'source-rate-1au':Decimal.from_float(float(r['rate_1au_s-1']))}
            for key in ('quiet','active') for r in rates if r['solar_activity']==('0.0' if key=='quiet' else '1.0')}
    threshold=Decimal.from_float(float(thresholds[0]['threshold_eV']))*ev
    return processed,hashes,source,threshold


def oracle(group,raw,source,threshold):
    """Closed-form photon Doppler, two-body COM and affine bath reference.

    This does not call any production kinematic function or infer real sodium
    states. The chosen +y relative direction is the exact mathematical image
    of the adapter's explicit(.5,.25) uniforms; finite native sin/cos roundoff
    is measured against the predeclared arithmetic bound, never retuned.
    """
    with localcontext() as context:
        context.prec=70
        D=Decimal;zero=D(0);c=D(299792458)
        m=D('4e-26');electron=D('1e-30');ion=m-electron;mass=ion+electron
        vector=lambda prefix:[raw[prefix+str(i)] for i in range(3)]
        rotation=[raw['input-R'+str(i)] for i in range(9)];offset=vector('input-offset-v');parent_event=vector('input-parent-event-v')
        transpose=lambda v:[sum(rotation[3*j+i]*v[j] for j in range(3)) for i in range(3)]
        rotate=lambda v:[sum(rotation[3*i+j]*v[j] for j in range(3)) for i in range(3)]
        dot=lambda a,b:sum(x*y for x,y in zip(a,b));norm=lambda v:dot(v,v).sqrt()
        parent=transpose([v-w for v,w in zip(parent_event,offset)])
        direction=vector('input-photon-n');direction=[x/norm(direction) for x in direction]
        gamma=1/(1-dot(parent,parent)/(c*c)).sqrt();doppler=gamma*(1-dot(parent,direction)/c)
        endpoint=source[group.split('-')[0]];rest=threshold+endpoint['mean-J'];energy=rest/doppler
        photon_p=[energy*x/c for x in direction];total=[m*v+p for v,p in zip(parent,photon_p)]
        boost=(gamma*gamma/(gamma+1)*dot(parent,photon_p)-gamma*energy)/(c*c)
        photon_rest_p=[p+boost*v for p,v in zip(photon_p,parent)]
        before=dot(parent,parent)*m/2;com=dot(total,total)/(2*mass)
        available=before+energy-threshold-com
        if available<0:raise ArithmeticError('analytical fixture unexpectedly has negative COM energy')
        reduced=ion*electron/mass;relative=(2*reduced*available).sqrt()
        pi=[ion/mass*p+(relative if i==1 else zero) for i,p in enumerate(total)]
        pe=[electron/mass*p-(relative if i==1 else zero) for i,p in enumerate(total)]
        vi=[p/ion for p in pi];ve=[p/electron for p in pe]
        vie=[v+w for v,w in zip(rotate(vi),offset)];vee=[v+w for v,w in zip(rotate(ve),offset)]
        bath_p=[p-incoming for p,incoming in zip(pe,photon_p)]
        bath_energy=dot(ve,ve)*electron/2-energy
        rotated=rotate(bath_p);bath_event_p=[p+electron*w for p,w in zip(rotated,offset)]
        bath_event_energy=bath_energy+dot(offset,rotated)+electron*dot(offset,offset)/2
        expected={'mass-parent':m,'mass-ion':ion,'mass-electron':electron,'threshold-J':threshold,
            'mean-J':endpoint['mean-J'],'source-Q-J':-threshold,'photon-rest-J':rest,
            'photon-inertial-J':energy,'doppler':doppler,'bath-inertial-J':bath_energy,
            'bath-event-J':bath_event_energy,'bath-mass-kg':electron,
            'source-rate-1au':endpoint['source-rate-1au'],'event-rate':endpoint['source-rate-1au']/4,
            'ion-reference-weight':D.from_float(.4),'ion-individual-correction':D('6.5625'),
            'ion-number':D('2.625'),'ion-internal-J':threshold,'remaining-s':D('.625'),
            'physical-epoch':D('100.375'),'basis-epoch':D(90)}
        scales={key:abs(value) for key,value in expected.items()}
        for prefix,value,scale in (
            ('parent-inertial-v',parent,max(norm(parent),D(1))),
            ('photon-inertial-p',photon_p,energy/c),
            ('photon-rest-p',photon_rest_p,rest/c),
            ('ion-inertial-v',vi,max(norm(vi),relative/ion)),
            ('electron-inertial-v',ve,max(norm(ve),relative/electron)),
            ('ion-event-v',vie,max(norm(vi),norm(offset),relative/ion)),
            ('electron-event-v',vee,max(norm(ve),norm(offset),relative/electron)),
            ('bath-inertial-p',bath_p,max(norm(pe),norm(photon_p))),
            ('bath-event-p',bath_event_p,max(norm(pe),norm(photon_p),electron*norm(offset)))):
            for i,v in enumerate(value):expected[prefix+str(i)]=v;scales[prefix+str(i)]=scale
        return expected,scales


def run(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();checks=[];comparisons=[];manifest={'scope':'actual source-backed numerical composition with ANALYTICAL STATES; NOT physical Na crosswalk/validation'}
    try:
        processed,hashes,source,threshold=preflight();manifest['external_data_sha256']=hashes
        compiler=shutil.which('g++')
        if not compiler:raise q.MissingPrerequisite('native C++17 compiler unavailable')
        fixture={'schema':'srcMoon-M3-analytical-composition-fixture-v1','scope':'ANALYTICAL ONLY, NOT PHYSICAL Na STATE DATA OR OBSERVATIONS',
            'neutral_reference_mass_kg':'4e-26','electron_reference_mass_kg':'1e-30','ion_formation':'exact native D04 threshold timesSIeV, fixture only',
            'rotation':'identity or exactz90deg','velocity_offset_m_s':[100,-50,25],
            'polar_uniform':'.5','azimuth_uniform':'.25','relative_arithmetic_bound':ARITHMETIC_BOUND,
            'source_is_real_frozen_D04':True,'real_Na_state_crosswalk_available':False}
        fixture_path=output/'ANALYTICAL_FIXTURE_NOT_PHYSICAL_PROFILE.json';fixture_path.write_text(json.dumps(fixture,indent=2)+'\n')
        paths=[ROOT/p for p in ('srcMoon/species/PhysicalReactionEvents.cpp','srcMoon/species/PhysicalReactionEvents.h',
            'srcMoon/species/PiecewiseHazard.cpp','srcMoon/species/PiecewiseHazard.h','srcMoon/species/ExactArithmetic.h',
            'srcMoon/species/ReducedKinematics.cpp','srcMoon/species/ReducedKinematics.h','srcMoon/species/ReactionExecution.cpp',
            'srcMoon/species/ReactionExecution.h','srcMoon/species/IonCollector.cpp','srcMoon/species/IonCollector.h',
            'srcMoon/species/Photochemistry.cpp','srcMoon/species/Photochemistry.h',
            'srcMoon/MoonTransportClock.cpp','srcMoon/MoonTransportClock.h',
            'srcMoon/test/physical_reaction_events/probe.cpp','srcMoon/test/physical_reaction_events/CONTRACT.md',
            'srcMoon/test/test_m3_physical_reaction_events.py','srcMoon/test/m3_qualification.py')]+shared_dependencies(ROOT,'photochemistry')
        manifest['source_sha256']={str(path):q.sha(path) for path in paths}
        cpp=[p for p in paths if p.suffix=='.cpp']
        executable=output/'physical-event-probe'
        command=[compiler,'-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT),*(str(p) for p in cpp),'-o',str(executable)]
        manifest['compile_command']=command
        compiled=subprocess.run(command,capture_output=True,text=True,timeout=180)
        (output/'compile.stdout.log').write_text(compiled.stdout);(output/'compile.stderr.log').write_text(compiled.stderr)
        if compiled.returncode:raise RuntimeError('actual production composer compile failed')
        command=[str(executable),str(processed),q.sha(fixture_path),str(output/'native')]
        native=subprocess.run(command,capture_output=True,text=True,timeout=120)
        manifest['native_command']=command
        (output/'stdout.log').write_text(native.stdout);(output/'stderr.log').write_text(native.stderr)
        with (output/'native.checks.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
        names=[r['name'] for r in rows]
        if len(names)!=len(set(names)) or set(names)!=EXPECTED or any(r['passed'] not in ('0','1') for r in rows):raise ValueError('actual native composer control registry/schema changed')
        failed=sum(r['passed']=='0' for r in rows)
        if native.returncode!=int(bool(failed)) or native.stdout.count('SUMMARY checks='+str(len(EXPECTED))+' failed='+str(failed)+'\n')!=1:
            raise ValueError('actual native summary/exit/complete registry mismatch')
        checks += [{'name':r['name'],'status':'PASS' if r['passed']=='1' else 'FAIL'} for r in rows]
        groups={}
        with (output/'native.metrics.csv').open(newline='') as stream:
            reader=csv.DictReader(stream)
            if reader.fieldnames!=['case','key','value']:raise ValueError('native metric schema changed')
            for row in reader:
                group=groups.setdefault(row['case'],{})
                if row['key'] in group:raise ValueError('duplicate native metric')
                value=Decimal(row['value'])
                if not value.is_finite():raise ValueError('nonfinite native composition metric')
                group[row['key']]=value
        if set(groups)!=set(GROUPS):raise ValueError('declared source endpoint/frame cohort changed')
        for group in GROUPS:
            expected,scales=oracle(group,groups[group],source,threshold)
            for key,reference in expected.items():
                actual=groups[group][key];scale=scales[key]
                normalized=abs(actual-reference)/scale if scale else abs(actual-reference)
                passed=normalized<=Decimal.from_float(ARITHMETIC_BOUND)
                checks.append({'name':group+':independent:'+key,'status':'PASS' if passed else 'FAIL'})
                comparisons.append({'case':group,'quantity':key,'production':str(actual),'independent_Decimal':str(reference),
                    'normalization_scale':str(scale),'normalized_error':str(normalized),'predeclared_arithmetic_bound':ARITHMETIC_BOUND})
        with (output/'comparison.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
        try:
            import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
            figure,axes=plt.subplots(1,2,figsize=(13,4))
            for key,axis,label in (('photon-inertial-J',axes[0],'Photon energy[J]'),('bath-event-J',axes[1],'Signed event-coordinate reservoir energy[J]')):
                selected=[r for r in comparisons if r['quantity']==key]
                axis.plot(range(len(selected)),[float(r['independent_Decimal']) for r in selected],label='independent Decimal',marker='x')
                axis.scatter(range(len(selected)),[float(r['production']) for r in selected],label='actual production')
                axis.set_xticks(range(len(selected)));axis.set_xticklabels([r['case'] for r in selected],rotation=35,ha='right');axis.set_ylabel(label);axis.legend(framealpha=1)
            figure.suptitle('Actual source/composer, ANALYTICAL masses/states — no physical Na validation');figure.tight_layout()
            for extension in ('png','eps'):figure.savefig(output/('comparison.'+extension),bbox_inches='tight')
            plt.close(figure)
        except ImportError:manifest['visualization']={'status':'SKIPPED','reason':'Matplotlib unavailable; independent numerical CSV retained'}
        if any(q.sha(path)!=expected for path,expected in manifest['source_sha256'].items()):raise RuntimeError('production/test source changed during composer verification')
        manifest.update(compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],
            executable_sha256=q.sha(executable),fixture_sha256=q.sha(fixture_path),
            source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    except q.MissingPrerequisite as error:checks.append({'name':'source/compiler-prerequisite','status':'SKIPPED','reason':str(error)})
    except Exception as error:
        (output/'runner-error.log').write_text(traceback.format_exc());checks.append({'name':'infrastructure','status':'ERROR','reason':str(error)})
    counts={status:sum(r['status']==status for r in checks) for status in ('PASS','FAIL','ERROR','SKIPPED')}
    status='ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'SKIPPED' if counts['SKIPPED'] else 'PASS'
    result={'schema':'srcMoon-M3-physical-event-composition-controls-v1','status':status,'counts':counts,'checks':checks,
        'scope':'source-bound numerical mechanics/transaction controls; ANALYTICAL states, NOT physical Na source crosswalk',
        'formal_M3_status':'SKIPPED / NOT QUALIFIED','physical_profile_activated':False,'observations_used':False,
        'comparison_count':len(comparisons),'elapsed_seconds':time.monotonic()-started}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();value=run(args.output_dir);print(json.dumps(value,indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
