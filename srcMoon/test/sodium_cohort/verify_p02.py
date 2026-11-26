#!/usr/bin/env python3
"""Independent source/actual-event P02 mechanical admission checks.

Only parsing/authentication is shared with the ledger audit. No production
kinematics/model-error evaluator is imported. Source decimals and native hex
values stay exact rationals; Decimal values are display diagnostics only.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import re
import sys
import traceback

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE));sys.path.insert(0, str(HERE.parent))
import verify_ledger as ledger
import atomic_targets
import audit
from prepare_physical_case import application_documents

BUDGET_SCHEMA = 'srcMoon-M3-P02-approved-mechanical-budget-v1'
BINDING_SCHEMA = 'srcMoon-M3-P02-owner-binding-v1'
CAPS = ('relative_reference_inertia_error_max','relative_acceleration_error_max',
        'relative_nonrelativistic_kinetic_error_max','particle_speed_max_m_s')
EXIT = {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}


def number(value):
    return ledger.decimal(value)


def shown(value):
    with localcontext() as context:
        context.prec = 60
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def load_budget(descriptor, base, identities):
    try:path = ledger.authenticate(descriptor, base, identities)
    except FileNotFoundError as error:raise ledger.Unavailable('approved budget file unavailable')from error
    value = ledger.strict_json(path)
    if value.get('schema') != BUDGET_SCHEMA or value.get('approval_status') != 'APPROVED':
        raise ValueError('explicit approved P02 budget schema/authority required')
    caps = {key:number(value[key]) for key in CAPS}
    if any(caps[key] < 0 or caps[key] >= 1 for key in CAPS[:3]) or caps[CAPS[3]] <= 0:
        raise ValueError('P02 relative caps require [0,1), positive explicit speed cap')
    return caps


def mass_errors(reference, physical):
    """Different denominators deliberately match the declared P02 quantities."""
    if reference <= 0 or physical <= 0:
        raise ValueError('positive source reference/physical masses required')
    difference = abs(reference-physical)
    return difference/physical, difference/reference


def kinetic_pass(speed2, c, limit):
    """Exact rational inequality for (Krel-Kclass)/Krel <= limit.

    With b=v²/c² and s=sqrt(1-b), the relative discrepancy is
    (1+b-s)/2. Compare s >= 1+b-2L by squaring ONLY when the right
    side is nonnegative. This avoids subtracting nearly equal energies,
    and cannot underflow positive subnormal-speed errors to a fake zero.
    The same reference mass cancels; reference/physical mass error remains
    a separate criterion with its own denominator.
    """
    if speed2 < 0 or c <= 0 or not 0 <= limit < 1:
        raise ValueError('nonnegative speed², positive c, limit in [0,1) required')
    b = speed2/(c*c)
    if b >= 1:
        return False
    rhs = 1+b-2*limit
    return rhs <= 0 or 1-b >= rhs*rhs


def kinetic_diagnostic(speed2, c):
    # Rationalized diagnostic expression retains tiny positive beta². It is
    # never used to decide the inequality above, including at a boundary.
    b = speed2/(c*c)
    if not 0 <= b < 1:
        return None
    with localcontext() as context:
        context.prec = 80
        beta2 = Decimal(b.numerator)/Decimal(b.denominator)
        s = (1-beta2).sqrt()
        return str(beta2*(2+s)/(2*(1+s)))


def invert_affine(matrix, offset, event_velocity):
    """Exact inverse of the ACTUAL recorded affine component transformation.

    No transpose=exact-inverse assumption is made for a rounded rotation.
    Proper-rotation admission is checked separately using its native bound.
    """
    if len(matrix) != 9 or len(offset) != 3 or len(event_velocity) != 3:
        raise ValueError('complete recorded affine transformation required')
    a,b,c,d,e,f,g,h,i = matrix
    det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    if det == 0:
        raise ValueError('singular actual inertial/event transformation')
    cof = [e*i-f*h,c*h-b*i,b*f-c*e,f*g-d*i,a*i-c*g,c*d-a*f,d*h-e*g,b*g-a*h,a*e-b*d]
    shifted = [x-y for x,y in zip(event_velocity,offset)]
    return [sum((cof[3*r+k]*shifted[k] for k in range(3)),F())/det for r in range(3)]


def native_controls(path):
    """Read ONLY existing owned scalar grammar; no physics values are chosen.

    Native ! comments and physical continuations are preserved. Each include
    must balance its blocks; includes resolve against the declaring file.
    Nested source blocks are inventoried but do not overwrite owner scalars.
    """
    values={};documents=[];chain=[];stack=[];owners=0
    def read(p):
        nonlocal owners
        p=Path(p).resolve()
        if p in chain or len(chain)>=64:raise ValueError('native include cycle/depth overflow')
        if p not in documents:documents.append(p)
        chain.append(p);entry=list(stack);pending=''
        try:
            for physical in p.read_text().splitlines():
                fragment=physical.split('!',1)[0].strip();more=fragment.endswith('\\')
                pending+=(' ' if pending and fragment else '')+(fragment[:-1].strip() if more else fragment)
                if more:continue
                text=pending;pending=''
                if not text:continue
                if text.lower().startswith('#include'):
                    match=re.fullmatch(r'#include\s+(.+)',text,re.I)
                    if not match:raise ValueError('malformed native include')
                    name=match[1].strip()
                    if name.startswith(('"','<')):
                        closing='"' if name[0]=='"' else '>'
                        if not name.endswith(closing):raise ValueError('malformed native include path')
                        name=name[1:-1]
                    child=Path(name);read(child if child.is_absolute() else p.parent/child);continue
                match=re.fullmatch(r'#begin\s+([A-Za-z_][A-Za-z_0-9]*)\s*:\s*(.+)',text,re.I)
                if match:
                    kind=match[1].lower()
                    if not stack:
                        if kind!='sodium_qualification' or owners:raise ValueError('one sodium owner required')
                        owners+=1
                    elif len(stack)!=1 or kind!='source':raise ValueError('unsupported nested native owner')
                    stack.append(kind);continue
                match=re.fullmatch(r'#end\s+([A-Za-z_][A-Za-z_0-9]*)',text,re.I)
                if match:
                    if not stack or len(stack)<=len(entry) or stack[-1]!=match[1].lower():raise ValueError('mismatched native owner end')
                    stack.pop();continue
                match=re.fullmatch(r'([A-Za-z_][A-Za-z_0-9.-]*)\s*=\s*(.+)',text)
                if not match or not stack:raise ValueError('malformed native owned assignment')
                if len(stack)==1:
                    key,value=match[1],match[2].strip()
                    if key in values:raise ValueError('duplicate native owner scalar '+key)
                    values[key]={'value':value,'file':p}
            if pending or stack!=entry:raise ValueError('unbalanced native include or continuation')
        finally:chain.pop()
    read(path)
    if owners!=1 or stack:raise ValueError('missing native sodium owner')
    return values,documents


def rotation_checks(matrix,absolute,relative):
    """Proper-rotation contract with explicit unit-vector norm scaleONE.

    Production declares columns of R as unit orthogonal vectors: R^T R=I.
    An off-diagonal target is zero, but its dimensionless normalization is
    still the product of UNIT vector norms, not the small measured dot product.
    Scaling by that dot product accidentally demands exact orthogonality.
    The determinant is likewise compared with the dimensionless unit scale.
    No caller bound, numerical epsilon or coefficient is changed here.
    """
    if len(matrix)!=9 or absolute<0 or relative<0:
        raise ValueError('complete matrix and nonnegative explicit rotation bounds required')
    limit=absolute+relative
    rows=[]
    for a in range(3):
        for b in range(3):
            dot=sum((matrix[3*k+a]*matrix[3*k+b]for k in range(3)),F())
            rows.append(('proper-rotation-'+str(a)+str(b),dot,dot-F(a==b),limit))
    r=matrix
    det=r[0]*(r[4]*r[8]-r[5]*r[7])-r[1]*(r[3]*r[8]-r[5]*r[6])+r[2]*(r[3]*r[7]-r[4]*r[6])
    rows.append(('proper-rotation-determinant',det,det-1,limit))
    return rows


def native_budget_checks(controls,caps):
    """Compare every actually parsed native cap with the pinned approval.

    Returning evaluated rows separates a genuinely changed cap (FAIL) from
    an unavailable native field (SKIPPED). The native value cannot supply its
    own independent approval, even if it would admit the observed velocities.
    """
    def scalar(key):
        if key not in controls:raise ledger.Unavailable('actual native P02 cap absent: '+key)
        return controls[key]['value']
    rows=[('all','native-speed-cap-equals-approved-budget',number(scalar('maximum_classical_speed_ms')),caps[CAPS[3]])]
    for state in ('parent','ion','electron'):
        for field,cap in [('inertial_mass',CAPS[0]),('acceleration',CAPS[1]),('relativistic_kinetic',CAPS[2])]:
            key=state+'_maximum_relative_'+field+'_error'
            rows.append((state,'native-'+field+'-cap-equals-approved-budget',number(scalar(key)),caps[cap]))
    return rows


def source_masses(profile_path, mapping, identities):
    states,constants=ledger.load_profile(profile_path,mapping,identities)
    profile=ledger.strict_json(profile_path)
    state_path=ledger.authenticate(profile['states'],Path(profile_path).parent,identities)
    rows={row['state_id']:row for row in ledger.csv_rows(state_path,ledger.STATE_HEADER)}
    # The P02 limits cover one admitted atomic target: the I21 Na/Na+ pair or,
    # under the same approved limits, the I22 He/He+, Ne/Ne+ or Ar/Ar+ pilot.
    parent_nuclei=states['parent']['nuclei']
    if len(parent_nuclei)!=1:raise ValueError('P02 budget scope requires a single-isotope atomic parent')
    nucleus,=parent_nuclei
    atomic_targets.element_of_nucleus(nucleus)
    isotope={nucleus:1}
    for key,charge,nuclei,kind in [('parent',0,isotope,'ground'),('ion',1,isotope,'ground'),('electron',-1,{},'free-electron')]:
        state=states[key];row=rows[state['id']]
        if state['charge']!=charge or state['nuclei']!=nuclei or row['state_kind']!=kind or number(row['excitation_J'])!=0:
            raise ValueError('P02 budget scope is explicitly an admitted atom/ion ground pair and free electron')
        state['physical_mass']=number(row['physical_rest_mass_kg'])
    c=constants['speed of light in vacuum'][0]
    return states,c


def verify(binding_path):
    binding_path=Path(binding_path).resolve();checks=[];ids={}
    result={'schema':'srcMoon-M3-independent-P02-v1','scope':'source-bound controlled Na event mechanical errors only',
            'formal_M3_status':'SKIPPED / NOT QUALIFIED BY THIS SUPPORTING AUDIT','observational_scoring':False}
    events_evaluated=0;owners_evaluated=0
    def check(owner,parent,state,name,passed,actual='',limit='',frame=''):
        checks.append(dict(owner=owner,original_parent_id=parent,state=state,criterion=name,
                           status='PASS' if passed else 'FAIL',actual=str(actual),limit=str(limit),frame=frame))
    try:
        if not binding_path.is_file():raise ledger.Unavailable('explicit P02 binding unavailable')
        ids[str(binding_path)]=ledger.digest(binding_path)
        value=ledger.strict_json(binding_path)
        if set(value)!={'schema','profile','budget','owners'} or value['schema']!=BINDING_SCHEMA:
            raise ValueError('explicit P02 owner binding schema required')
        caps=load_budget(value['budget'],binding_path.parent,ids)
        try:profile_path=ledger.authenticate(value['profile'],binding_path.parent,ids)
        except FileNotFoundError as error:raise ledger.Unavailable('source profile unavailable')from error
        owners=value['owners']
        if not isinstance(owners,list):raise ValueError('P02 owners must be an explicit list')
        if not owners:raise ledger.Unavailable('no completed physical native owner supplied')
        seen=set()
        for owned in owners:
            if not isinstance(owned,dict) or set(owned)!={'owner_manifest','owner_result','mpi_ranks'} or type(owned['mpi_ranks']) is not int or owned['mpi_ranks'] not in (1,8):
                raise ValueError('explicit completed1/8MPI owner descriptors required')
            mp=ledger.authenticate(owned['owner_manifest'],binding_path.parent,ids)
            rp=ledger.authenticate(owned['owner_result'],binding_path.parent,ids)
            if mp in seen:raise ValueError('duplicate native owner is not additional event evidence')
            seen.add(mp);manifest=ledger.strict_json(mp)
            launch={**owned,'owner_manifest':{'path':str(mp),'sha256':ids[str(mp)]},
                    'owner_result':{'path':str(rp),'sha256':ids[str(rp)]},
                    'input':manifest['input'],'output':str(mp.parent/'native-output')}
            native=audit.authenticate_launch(launch);receipt=ledger.strict_json(native/'input_receipt.json')
            if receipt.get('schema')!='srcMoon-native-sodium-input-receipt-v1' or receipt.get('phase')!='explicit-experimental-physical-input':
                raise ledger.Unavailable('actual physical native input receipt unavailable')
            mapping=receipt['mapping']
            if mapping.get('profile_sha256')!=value['profile']['sha256'] or receipt['profile']['sha256']!=value['profile']['sha256']:
                raise ValueError('actual event profile differs from frozen P02 profile')
            if mapping.get('event_frame')!='LSO' or mapping.get('photon_energy_policy')!='source-effective-photon':
                raise ValueError('actual source/frame convention is outside approved P02 controlled scope')
            try:states,c=source_masses(profile_path,mapping,ids)
            except FileNotFoundError as error:raise ledger.Unavailable('required profile dependency unavailable')from error
            if c!=F(299792458):raise ValueError('P02 SI profile must retain exactSI c=299792458m/s')
            if caps[CAPS[3]]>=c:raise ValueError('declared speed cap must be subluminal')
            actual_docs=audit.input_documents(receipt['input_documents'])
            for path,digest in actual_docs.items():
                if ledger.digest(path)!=digest:raise ValueError('native input graph changed after execution')
                ids[path]=digest
            controls,docs=native_controls(manifest['input']['path'])
            def scalar(key):
                if key not in controls:raise ledger.Unavailable('actual native scalar unavailable: '+key)
                return controls[key]['value']
            if scalar('inertial_origin')!='SUN':raise ValueError('P02 audit requires explicitlySUN-origin inertial coordinates')
            frame='SUN-origin:'+scalar('inertial_frame')
            record=controls.get('application_input')
            if record is None:raise ledger.Unavailable('actual application input identity absent')
            name=record['value']
            if name.startswith('"') and name.endswith('"'):name=name[1:-1]
            app=Path(name);app=(app if app.is_absolute() else record['file'].parent/app).resolve()
            app_values,app_docs=application_documents(app)
            if app_values.get('tracing_clock_mode',{}).get('value')!='common-time':raise ValueError('P02 controlled cohort requires actual common-time input')
            complete={str(p) for p in docs}|{r['path']for r in app_docs}
            if complete!=set(actual_docs):raise ValueError('actual input receipt omits/adds an include document')
            for state,label,actual,limit in native_budget_checks(controls,caps):
                check(str(mp),'*',state,label,actual==limit,shown(actual),shown(limit),frame)
            for state in ('parent','ion','electron'):
                inertia,acceleration=mass_errors(states[state]['mass'],states[state]['physical_mass'])
                for label,error,cap in [('source-reference-inertia',inertia,CAPS[0]),('source-reference-acceleration',acceleration,CAPS[1])]:
                    check(str(mp),'*',state,label,error<=caps[cap],shown(error),shown(caps[cap]),'source state scalar')
                check(str(mp),'*',state,'declared-maximum-speed-kinetic-error',kinetic_pass(caps[CAPS[3]]**2,c,caps[CAPS[2]]),kinetic_diagnostic(caps[CAPS[3]]**2,c),shown(caps[CAPS[2]]),frame)
            files=[native/('event_metrics.rank%d.hex.csv'%rank)for rank in range(owned['mpi_ranks'])]
            for path in [native/'input_receipt.json',native/'parent_records.csv']+files:ids[str(path)]=ledger.digest(path)
            events=ledger.read_events(files)
            parents=ledger.csv_rows(native/'parent_records.csv',list(audit.q.PARENT_COLUMNS))
            if type(receipt['original_parents'])is not int or receipt['original_parents']<=0 or {r['original_parent_id']for r in parents}!={str(i)for i in range(receipt['original_parents'])} or len(parents)!=receipt['original_parents']:
                raise ValueError('incomplete/duplicate actual original-parent cohort')
            reacted={r['original_parent_id']for r in parents if r['status']=='REACTED'}
            if any(r['status']not in ('REACTED','SURVIVED')for r in parents) or set(events)!=reacted:
                raise ValueError('P02 event set does not reconcile complete original-parent removals')
            def vector(row,key):
                try:return [row[key+str(k)]for k in range(3)]
                except KeyError as error:raise ledger.Unavailable('actual event velocity/transform missing: '+str(error))
            rotation=receipt['bounds']['rotation'];ra,rr=number(rotation['absolute']),number(rotation['relative'])
            if ra<0 or rr<0:raise ValueError('negative native rotation bound')
            for parent,row in events.items():
                try:r=[row['inertial_to_event_r'+str(k)]for k in range(9)]
                except KeyError as error:raise ledger.Unavailable('actual event matrix incomplete')from error
                rotation_ok=True
                for label,actual,residual,limit in rotation_checks(r,ra,rr):
                    passed=abs(residual)<=limit;rotation_ok=rotation_ok and passed
                    check(str(mp),parent,'all',label,passed,shown(actual),shown(limit),'recorded inertial-to-LSO; dimensionless unit scale1')
                if not rotation_ok:continue
                actual_ion=vector(row,'actual_ion_v');offset=vector(row,'velocity_offset_event')
                inertial={'parent':vector(row,'parent_inertial_v'),'ion':invert_affine(r,offset,actual_ion),'electron':vector(row,'electron_inertial_v')}
                event={'parent':vector(row,'parent_event_v'),'ion':actual_ion,'electron':vector(row,'electron_event_v')}
                for state in ('parent','ion','electron'):
                    for label,velocity in [('inertial',inertial[state]),('event-LSO',event[state])]:
                        v2=sum((v*v for v in velocity),F())
                        check(str(mp),parent,state,'actual-'+label+'-speed-cap',v2<=caps[CAPS[3]]**2,shown(v2),shown(caps[CAPS[3]]**2),frame if label=='inertial'else'LSO event coordinates; speed guard only')
                    speed2=sum((v*v for v in inertial[state]),F())
                    check(str(mp),parent,state,'actual-inertial-classical-kinetic-error',kinetic_pass(speed2,c,caps[CAPS[2]]),kinetic_diagnostic(speed2,c),shown(caps[CAPS[2]]),frame+('; actual ion affine reconstruction'if state=='ion'else''))
                    if state in ('parent','ion'):
                        key='actual_'+state+'_compiled_mass_kg'
                        if key not in row:raise ledger.Unavailable('actual compiled mass diagnostic missing')
                        ie,ae=mass_errors(row[key],states[state]['physical_mass'])
                        for label,error,cap in [('actual-PIC-inertia',ie,CAPS[0]),('actual-PIC-acceleration',ae,CAPS[1])]:
                            check(str(mp),parent,state,label,error<=caps[cap],shown(error),shown(caps[cap]),'actual compiled mass versus source physical mass')
                events_evaluated+=1
            owners_evaluated+=1
        if not events_evaluated:raise ledger.Unavailable('no actual event velocities evaluated; source/guard checks alone do not verify physical events')
    except ledger.Unavailable as error:
        checks.append(dict(owner='',original_parent_id='*',state='',criterion='required-evidence',status='SKIPPED',actual=str(error),limit='',frame=''))
    except Exception as error:
        result['exception']=traceback.format_exc()
        checks.append(dict(owner='',original_parent_id='*',state='',criterion='source-owner-infrastructure',status='ERROR',actual=str(error),limit='',frame=''))
    counts={s:sum(r['status']==s for r in checks)for s in EXIT}
    result.update(status=next(s for s in ('ERROR','FAIL','SKIPPED','PASS')if counts[s]),counts=counts,checks=checks,input_sha256=ids,events_evaluated=events_evaluated,owners_evaluated=owners_evaluated)
    return result


def write_artifacts(result,output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    with(output/'comparison.csv').open('w',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=['owner','original_parent_id','state','criterion','status','actual','limit','frame']);writer.writeheader();writer.writerows(result['checks'])
    (output/'checks.log').write_text('\n'.join(r['status']+' '+r['criterion']+' '+r['actual']for r in result['checks'])+'\n')
    (output/'run_manifest.json').write_text(json.dumps({'schema':'srcMoon-M3-P02-verifier-manifest-v1','command':sys.argv,'cwd':str(Path.cwd()),'verifier_sha256':{str(p):ledger.digest(p)for p in [Path(__file__),Path(ledger.__file__),Path(audit.__file__),Path(audit.q.__file__),Path(application_documents.__code__.co_filename),HERE.parent/'m2_initial_state.py']},'input_sha256':result['input_sha256'],'scope':result['scope']},indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--binding',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();result=verify(args.binding);write_artifacts(result,args.output_dir)
    print(json.dumps(result,indent=2));raise SystemExit(EXIT[result['status']])
