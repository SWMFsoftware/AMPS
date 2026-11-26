#!/usr/bin/env python3
"""Independent exact-rational audit of actual committed Na event records.

This module does not import production chemistry, mass, kinematic or ledger
code. It authenticates the supplied profile and independently reconstructs its
bookkeeping. It is a supporting verifier, never a source/readiness authority.
All native scalars are hexadecimal exact binary values; no long-double value
passes through Python's binary64 float. SI source decimals remain rationals.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
from decimal import Decimal,InvalidOperation
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import traceback
import sys
import subprocess
from datetime import datetime,timezone

sys.path.insert(0,str(Path(__file__).resolve().parent))
import atomic_targets


class Unavailable(Exception):
    """Required source/capability evidence is absent, rather than evaluated."""


def strict_json(path):
    def pairs(items):
        result={}
        for key,value in items:
            if key in result:raise ValueError('duplicate JSON key: '+key)
            result[key]=value
        return result
    def invalid(token):raise ValueError('nonfinite JSON token: '+token)
    value=json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=invalid)
    if not isinstance(value,dict):raise ValueError('JSON object required: '+str(path))
    return value


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def decimal(value):
    if not isinstance(value,str) or not value:raise ValueError('explicit decimal string required')
    number=Decimal(value)
    if not number.is_finite():raise ValueError('nonfinite decimal')
    return F(number)


def binary(value):
    """Decode native std::hexfloat exactly, including long-double exponents."""
    match=re.fullmatch(r'([+-]?)0x([0-9a-fA-F]+)(?:\.([0-9a-fA-F]*))?[pP]([+-]?\d+)',value)
    if not match:raise ValueError('exact native hexadecimal scalar required: '+value)
    sign,whole,fraction,exponent=match.groups();fraction=fraction or ''
    numerator=int(whole+fraction,16)*(-1 if sign=='-' else 1)
    power=int(exponent)-4*len(fraction)
    return F(numerator*(1<<power),1) if power>=0 else F(numerator,1<<(-power))


def rational_text(value):return str(value.numerator)+'/'+str(value.denominator)

def authenticate(descriptor,base,identities):
    if not isinstance(descriptor,dict) or not isinstance(descriptor.get('path'),str) or not re.fullmatch('[0-9a-f]{64}',descriptor.get('sha256','')):
        raise ValueError('source descriptor needs explicit path and lowercase SHA-256')
    path=Path(descriptor['path']);path=(path if path.is_absolute() else base/path).resolve()
    actual=digest(path)
    if actual!=descriptor['sha256']:raise ValueError('source SHA-256 changed: '+str(path))
    identities[str(path)]=actual
    return path


def csv_rows(path,header):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream)
        if reader.fieldnames!=header:raise ValueError('exact CSV schema mismatch: '+str(path))
        rows=list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):raise ValueError('incomplete CSV row: '+str(path))
    return rows

STATE_HEADER=['state_id','state_label','state_kind','nuclei','charge','formation_J','formation_source_id','formation_source_record','excitation_J','excitation_source_id','excitation_source_record','physical_rest_mass_kg','physical_mass_source_id','physical_mass_source_record','status','missing_reason']


def load_profile(profile_path,mapping,identities):
    profile_path=Path(profile_path).resolve();profile=strict_json(profile_path)
    expected=mapping.get('profile_sha256')
    if not re.fullmatch('[0-9a-f]{64}',expected or ''):raise ValueError('mapping must pin the profile JSON SHA-256')
    if digest(profile_path)!=expected:raise ValueError('source profile JSON identity changed')
    identities[str(profile_path)]=expected
    if profile.get('schema')!='srcMoon-M3-state-profile-v1':raise ValueError('unsupported physical state profile schema')
    sources={}
    for descriptor in profile['sources']:
        if descriptor['id'] in sources:raise ValueError('duplicate profile source identity')
        authenticate(descriptor,profile_path.parent,identities);sources[descriptor['id']]=descriptor
    isotope_path=authenticate(profile['isotopes'],profile_path.parent,identities)
    constant_path=authenticate(profile['constants'],profile_path.parent,identities)
    state_path=authenticate(profile['states'],profile_path.parent,identities)
    def source(row,prefix,role):
        sid=row[prefix+'source_id'];record=row[prefix+'source_record']
        if sid not in sources or not record or sources[sid]['role'] not in role:raise ValueError('missing/incompatible source record: '+sid)
    masses={}
    for row in csv_rows(isotope_path,['isotope','neutral_mass_kg','source_id','source_record']):
        if row['isotope'] in masses:raise ValueError('duplicate isotope mass')
        source(row,'',{'isotope-mass'});masses[row['isotope']]=decimal(row['neutral_mass_kg'])
        if masses[row['isotope']]<=0:raise ValueError('isotope mass must be positive')
    constants={}
    for row in csv_rows(constant_path,['name','value','unit','source_id','source_record']):
        if row['name'] in constants:raise ValueError('duplicate constant')
        source(row,'',{'constants'});constants[row['name']]=(decimal(row['value']),row['unit'])
    for name,unit in [('electron mass','kg'),('elementary charge','C'),('speed of light in vacuum','m s^-1')]:
        if name not in constants or constants[name][1]!=unit or constants[name][0]<=0:raise ValueError('missing positive SI constant '+name)
    states={}
    for row in csv_rows(state_path,STATE_HEADER):
        if row['state_id'] in states:raise ValueError('duplicate state identity')
        states[row['state_id']]=row
    selected={}
    for key in ('parent','ion','electron'):
        sid=mapping.get(key)
        if sid not in states:raise ValueError('explicit mapped state not in profile: '+str(sid))
        row=states[sid]
        if row['status']=='unavailable':raise Unavailable('required '+key+' state unavailable: '+row['missing_reason'])
        if row['status']!='available':raise ValueError('invalid source state status')
        if not re.fullmatch(r'-?(0|[1-9][0-9]*)',row['charge']) or row['charge']=='-0':raise ValueError('canonical integer charge required')
        charge=int(row['charge']);nuclei={}
        if row['nuclei']:
            for term in row['nuclei'].split(';'):
                atom,count=term.rsplit(':',1)
                if atom not in masses or atom in nuclei or not re.fullmatch('[1-9][0-9]*',count):raise ValueError('invalid isotope count')
                nuclei[atom]=int(count)
        if not nuclei and not (key=='electron' and row['state_kind']=='free-electron' and charge==-1):raise ValueError('only explicit free electron can omit nuclei')
        formation,excitation=decimal(row['formation_J']),decimal(row['excitation_J'])
        source(row,'formation_',{'energy','reference-definition'});source(row,'excitation_',{'energy','reference-definition'})
        mass=sum((masses[a]*n for a,n in nuclei.items()),F())-charge*constants['electron mass'][0]
        if mass<=0:raise ValueError('nonpositive independently reconstructed reference mass')
        if not row['physical_rest_mass_kg']:raise Unavailable('independent physical rest mass missing: '+sid)
        source(row,'physical_mass_',{'isotope-mass','constants','energy'})
        if decimal(row['physical_rest_mass_kg'])<=0:raise ValueError('nonpositive physical rest mass')
        selected[key]={'id':sid,'nuclei':nuclei,'charge':charge,'mass':mass,'internal':formation+excitation,
            'formation_source_identity':row['formation_source_id']+'@sha256:'+sources[row['formation_source_id']]['sha256']+'#'+row['formation_source_record']}
    return selected,constants



def load_photon_source(mapping,base,identities,constants,scope):
    """Independently select one explicitly pinned PHIDRATES endpoint row.

    This is not a threshold crosswalk or a claim that its daughter laws are
    complete. No native receipt value supplies its own external reference.
    """
    source=mapping.get('photon_source')
    if source is None:
        if scope=='source-bound-experimental':raise Unavailable('independent pinned PHIDRATES rate/threshold endpoint unavailable')
        return None
    if not isinstance(source,dict) or not source.get('reaction_id'):raise ValueError('explicit photon reaction ID required')
    activity=decimal(source.get('solar_activity'))
    if activity not in (0,1):raise ValueError('only explicit native quiet/active endpoint supported')
    rates=authenticate(source['rates'],base,identities);thresholds=authenticate(source['thresholds'],base,identities)
    rateHeader=['species','reaction_id','reaction_type','products','rate_1au_s-1','excess_energy_eV','solar_activity','optical_depth','distance_au','source_reference','source_file','source_record_or_index']
    thresholdHeader=['species','reaction_id','products','threshold_wavelength_angstrom','threshold_wavelength_nm','threshold_eV','source_reference','source_file','source_record_or_index']
    matches=[row for row in csv_rows(rates,rateHeader) if row['reaction_id']==source['reaction_id'] and decimal(row['solar_activity'])==activity and decimal(row['optical_depth'])==0 and decimal(row['distance_au'])==1]
    limits=[row for row in csv_rows(thresholds,thresholdHeader) if row['reaction_id']==source['reaction_id']]
    if len(matches)!=1 or len(limits)!=1:raise ValueError('nonunique/missing native source endpoint or threshold')
    row,limit=matches[0],limits[0]
    # One admitted atomic target (I21 Na; I22 atomic pilot He/Ne/Ar), from the
    # independent verifier table; the selected states must carry its isotope.
    element=row['species']
    if (row['reaction_type']!='photoionization' or element not in atomic_targets.TARGETS or
            row['products']!=atomic_targets.products(element) or limit['products']!=row['products']):
        raise ValueError('source is not the complete single reactive endpoint of an admitted atomic target')
    unit=constants['elementary charge'][0]
    return {'element':element,'rate':decimal(row['rate_1au_s-1']),'mean':decimal(row['excess_energy_eV'])*unit,
            'threshold':decimal(limit['threshold_eV'])*unit,'electron_volt_J':unit}


def verify_energy_crosswalk(mapping,states,source,base,identities):
    """Authenticate an explicit alternative source model without fitting a gap.

    Default strict receipts retain their original source=state criterion. For
    source-effective photons the independent verifier reads the pinned native
    crosswalk and separately binds state IDs, source threshold and the actual
    resolved ion formation-energy evidence. No native precomputed residual is
    an independent reference and no energy correction is user-supplied.
    """
    policy=mapping.get('photon_energy_policy','source-equals-state')
    if policy=='source-equals-state':return policy
    if policy!='source-effective-photon':raise ValueError('unsupported photon energy policy')
    if source is None:raise Unavailable('source-effective photon requires independent native photon-source evidence')
    descriptor=mapping.get('source_crosswalk')
    if descriptor is None:raise Unavailable('source-effective photon crosswalk unavailable')
    path=authenticate(descriptor,base,identities)
    header=['channel_id','parent_state_id','ion_state_id','electron_state_id','energy_reference_identity',
            'photon_energy_policy','source_threshold_eV','physical_threshold_source_identity']
    rows=csv_rows(path,header)
    if len(rows)!=1:raise ValueError('one complete Na source crosswalk required')
    row=rows[0]
    if [row[k]for k in header[:5]]!=[mapping['photon_source']['reaction_id'],states['parent']['id'],
            states['ion']['id'],states['electron']['id'],mapping.get('energy_reference_identity')]:
        raise ValueError('crosswalk channel/state/common-energy identity differs from the independently loaded profile')
    if row['photon_energy_policy']!=policy or row['physical_threshold_source_identity']!=states['ion']['formation_source_identity']:
        raise ValueError('crosswalk policy or physical formation-energy source identity differs')
    if decimal(row['source_threshold_eV'])*source['electron_volt_J']!=source['threshold']:
        raise ValueError('crosswalk changed the native source threshold')
    return policy


def _bound(bounds,key):
    value=bounds.get(key)
    if not isinstance(value,dict) or set(value)!={'absolute','relative'}:raise ValueError('explicit frozen absolute/relative bound required: '+key)
    a,r=decimal(value['absolute']),decimal(value['relative'])
    if a<0 or r<0:raise ValueError('negative numerical bound: '+key)
    return a,r


def read_events(paths):
    if not paths:raise Unavailable('no actual per-rank event files supplied')
    events={};owners={}
    for path in paths:
        path=Path(path).resolve()
        for row in csv_rows(path,['original_parent_id','key','value']):
            sid=row['original_parent_id']
            if not re.fullmatch('0|[1-9][0-9]*',sid):raise ValueError('canonical original parent ID required')
            event=events.setdefault(sid,{})
            if sid in owners and owners[sid]!=path:raise ValueError('one event original parent appears on multiple owners')
            owners[sid]=path
            if row['key'] in event:raise ValueError('duplicate native event scalar: '+sid+'/'+row['key'])
            event[row['key']]=binary(row['value'])
    return events


def verify(paths,profile_json,mapping,bounds):
    """Return independent criteria; inputs/data are never modified.

    `mapping` pins profile_sha256, parent/ion/electron IDs, expected compiled
    indices and event_frame. `bounds` has a schema, scope, explicit numerical
    bound objects and source/run receipt descriptors retained by the caller.
    Caller owns actual pre-execution freezing/producer proof; this component
    never converts a declared bound or source hash into a qualification claim.
    """
    result={'schema':'srcMoon-M3-independent-event-ledger-v1','status':'ERROR','scope':'independent committed-event bookkeeping only; NOT physical M3 qualification','formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','counts':{s:0 for s in ('PASS','FAIL','ERROR','SKIPPED')},'checks':[],'input_sha256':{}}
    try:
        if not isinstance(mapping,dict) or not isinstance(bounds,dict):raise ValueError('explicit mapping/bounds objects required')
        if bounds.get('schema')!='srcMoon-M3-event-ledger-bounds-v1' or not mapping.get('event_frame') or mapping['event_frame']!=bounds.get('event_frame'):raise ValueError('explicit matching event frame/bounds schema required')
        if bounds.get('scope') not in ('analytical-software-fixture','source-bound-experimental'):raise ValueError('explicit limited evidence scope required')
        result['evidence_scope']=bounds['scope']
        for native,si in [('mass','mass_kg'),('momentum','momentum_kg_m_s'),('energy','energy_j')]:
            if si not in bounds and native in bounds:bounds={**bounds,si:bounds[native]}
        states,constants=load_profile(profile_json,mapping,result['input_sha256'])
        for name in ('mass_kg','momentum_kg_m_s','energy_j','number','rotation'):_bound(bounds,name)
        for path in paths:result['input_sha256'][str(Path(path).resolve())]=digest(path)
        source=load_photon_source(mapping,Path(profile_json).resolve().parent,result['input_sha256'],constants,bounds['scope'])
        if source is not None:
            # The source parent's isotope must be the selected parent and ion.
            isotope={atomic_targets.nucleus(source['element']):1}
            if states['parent']['nuclei']!=isotope or states['ion']['nuclei']!=isotope:
                raise ValueError('selected parent/ion states do not carry the '+source['element']+' source isotope')
        energy_policy=verify_energy_crosswalk(mapping,states,source,Path(profile_json).resolve().parent,result['input_sha256'])
        if source is not None:_bound(bounds,'rate_s')
        events=read_events(paths);c=constants['speed of light in vacuum'][0]
        def check(parent,name,passed,residual=F(),limit=F()):
            row={'original_parent_id':parent,'criterion':name,'status':'PASS' if passed else 'FAIL','residual_exact':rational_text(residual),'limit_exact':rational_text(limit)}
            result['checks'].append(row)
        def close(parent,name,actual,expected,key,scale=None):
            a,r=_bound(bounds,key);scale=max(abs(actual),abs(expected)) if scale is None else abs(scale);limit=a+r*scale
            check(parent,name,abs(actual-expected)<=limit,actual-expected,limit)
        def vector(row,name):
            try:return [row[name+str(k)] for k in range(3)]
            except KeyError as error:raise Unavailable('required native vector/readback absent: '+str(error)) from error
        def dot(a,b):return sum((x*y for x,y in zip(a,b)),F())
        def kinetic(m,v):return m*dot(v,v)/2
        def mat(r,v):return [dot(r[3*k:3*k+3],v) for k in range(3)]
        required={'parent_number','ion_number','parent_reference_weight','parent_individual_correction','ion_reference_weight','ion_individual_correction','parent_internal_j','ion_internal_j','photon_inertial_j','bath_event_j','bath_mass_kg','electron_state_j','source_threshold_j','source_mean_excess_j','source_q_j','ledger_parent_mass_kg','ledger_ion_mass_kg','ledger_bath_mass_kg','ledger_tracked_de_j','ledger_bath_de_j','tracked_charge_change_e','bath_charge_change_e','exposure_error_bound','clock_storage_error_bound','actual_ion_reference_weight','actual_ion_individual_correction','actual_ion_internal_j','actual_ion_compiled_index','actual_ion_child_count','event_epoch_et_s','actual_ion_epoch_et_s','planned_interval_start_s','planned_interval_duration_s','planned_interval_elapsed_s','actual_ion_interval_start_s','actual_ion_interval_duration_s','actual_ion_interval_elapsed_s','actual_parent_compiled_mass_kg','actual_ion_compiled_mass_kg','actual_parent_charge_c','actual_ion_charge_c'}
        p,i,e=(states[k] for k in ('parent','ion','electron'))
        if energy_policy=='source-effective-photon':
            required.update(('physical_threshold_j','physical_minus_source_threshold_j'))
        composition=Counter(i['nuclei']);composition.update(e['nuclei'])
        source_q=p['internal']-i['internal']-e['internal']
        for parent,row in sorted(events.items(),key=lambda item:int(item[0])):
            missing=required-set(row)
            if missing:raise Unavailable('actual committed-child/source metrics unavailable: '+','.join(sorted(missing)))
            check(parent,'declared-Na-like-charge-topology',p['charge']==0 and i['charge']==1 and e['charge']==-1)
            check(parent,'integer-nuclei-stoichiometry',dict(composition)==p['nuclei'])
            check(parent,'integer-charge-stoichiometry',p['charge']==i['charge']+e['charge'])
            if source is not None:
                if 'source_rate_1au_s' not in row:raise Unavailable('native source rate endpoint absent')
                close(parent,'independent-PHIDRATES-rate-endpoint',row['source_rate_1au_s'],source['rate'],'rate_s')
                close(parent,'independent-PHIDRATES-threshold',row['source_threshold_j'],source['threshold'],'energy_j')
                close(parent,'independent-PHIDRATES-mean-excess',row['source_mean_excess_j'],source['mean'],'energy_j')
            close(parent,'actual-PIC-parent-source-reference-mass',row['actual_parent_compiled_mass_kg'],p['mass'],'mass_kg')
            close(parent,'actual-PIC-ion-source-reference-mass',row['actual_ion_compiled_mass_kg'],i['mass'],'mass_kg')
            for name,state in [('parent',p),('ion',i)]:
                expected_charge=F.from_float(float(state['charge']*constants['elementary charge'][0]))
                check(parent,'actual-PIC-source-SI-charge-'+name,row['actual_'+name+'_charge_c']==expected_charge,row['actual_'+name+'_charge_c']-expected_charge)
            close(parent,'source-formation-Q',row['source_q_j'],source_q,'energy_j')
            if energy_policy=='source-equals-state':
                close(parent,'source-threshold-Q',row['source_threshold_j'],-source_q,'energy_j')
            else:
                close(parent,'physical-threshold-Q',row['physical_threshold_j'],-source_q,'energy_j')
                close(parent,'separate-physical-minus-source-threshold',row['physical_minus_source_threshold_j'],
                      -source_q-row['source_threshold_j'],'energy_j')
            close(parent,'parent-source-internal-energy',row['parent_internal_j'],p['internal'],'energy_j')
            close(parent,'ion-source-internal-energy',row['ion_internal_j'],i['internal'],'energy_j')
            close(parent,'electron-source-internal-energy',row['electron_state_j'],e['internal'],'energy_j')
            n=row['parent_reference_weight']*row['parent_individual_correction'];ni=row['ion_reference_weight']*row['ion_individual_correction']
            check(parent,'positive-represented-weights',all(row[k]>0 for k in ('parent_reference_weight','parent_individual_correction','ion_reference_weight','ion_individual_correction')))
            close(parent,'parent-represented-number',row['parent_number'],n,'number')
            close(parent,'ion-represented-number',row['ion_number'],ni,'number')
            close(parent,'represented-stoichiometry',ni,n,'number')
            for key in ('reference_weight','individual_correction','internal_j'):
                check(parent,'actual-PIC-'+key,row['actual_ion_'+key]==row['ion_'+key],row['actual_ion_'+key]-row['ion_'+key])
            check(parent,'actual-PIC-child-count',row['actual_ion_child_count']==1)
            compiled=mapping.get('ion_compiled_species')
            if not isinstance(compiled,int) or isinstance(compiled,bool) or compiled<0:raise ValueError('explicit nonnegative ion compiled index required')
            check(parent,'actual-PIC-species-index',row['actual_ion_compiled_index']==compiled)
            catalog_index=mapping.get('ion_catalog_species')
            if not isinstance(catalog_index,int) or isinstance(catalog_index,bool) or catalog_index<0:raise ValueError('explicit nonnegative ion catalog index required')
            if 'actual_ion_catalog_species' not in row:raise Unavailable('actual committed ion catalog state absent')
            check(parent,'actual-PIC-catalog-state',row['actual_ion_catalog_species']==catalog_index)
            check(parent,'actual-PIC-event-epoch',row['actual_ion_epoch_et_s']==row['event_epoch_et_s'])
            for part in ('start','duration','elapsed'):
                check(parent,'actual-PIC-original-interval-'+part,row['actual_ion_interval_'+part+'_s']==row['planned_interval_'+part+'_s'])
            check(parent,'valid-original-interval',row['planned_interval_duration_s']>0 and 0<=row['planned_interval_elapsed_s']<=row['planned_interval_duration_s'])
            pv,iv,ev=vector(row,'parent_event_v'),vector(row,'ion_event_v'),vector(row,'electron_event_v')
            pa,ia=vector(row,'parent_event_x'),vector(row,'actual_ion_x')
            actualv=vector(row,'actual_ion_v')
            for axis in range(3):
                check(parent,'actual-PIC-velocity-'+str(axis),actualv[axis]==iv[axis],actualv[axis]-iv[axis])
                check(parent,'actual-PIC-position-'+str(axis),ia[axis]==pa[axis],ia[axis]-pa[axis])
            pin,ein=vector(row,'parent_inertial_v'),vector(row,'electron_inertial_v')
            photon=vector(row,'photon_inertial_p');photonE=row['photon_inertial_j']
            offset=vector(row,'velocity_offset_event')
            try:rotation=[row['inertial_to_event_r'+str(k)] for k in range(9)]
            except KeyError as error:raise Unavailable('actual frame matrix missing') from error
            # Orthogonality concerns UNIT columns, including zero expected
            # cross-dot products. Its dimensionless reference scale is one,
            # as in the declared R^T R identity. A generic relative comparison
            # with the tiny off-diagonal value would incorrectly require that
            # rounded SPICE coefficients have exactly zero cross-dot products.
            # Keep the caller's original rotation tolerances unchanged.
            for a in range(3):
                for b in range(3):close(parent,'proper-frame-orthogonality-'+str(a)+str(b),
                    dot(rotation[a::3],rotation[b::3]),F(a==b),'rotation',scale=F(1))
            determinant=rotation[0]*(rotation[4]*rotation[8]-rotation[5]*rotation[7])-rotation[1]*(rotation[3]*rotation[8]-rotation[5]*rotation[6])+rotation[2]*(rotation[3]*rotation[7]-rotation[4]*rotation[6])
            close(parent,'proper-frame-determinant',determinant,F(1),'rotation',scale=F(1))
            for label,actual,sourcev,mass in [('parent',pv,pin,p['mass']),('electron',ev,ein,e['mass'])]:
                rotated=mat(rotation,sourcev)
                for axis in range(3):close(parent,label+'-event-frame-velocity-'+str(axis),mass*actual[axis],mass*(rotated[axis]+offset[axis]),'momentum_kg_m_s')
            # Photon energy/momentum and rest-frame Doppler are checked without
            # square roots: compare exact squared positive intervals. The bound
            # remains the caller's energy tolerance, not a chosen squared gate.
            a,r=_bound(bounds,'energy_j');energyLimit=a+r*abs(photonE)
            photonNorm2=dot(photon,photon)*c*c
            check(parent,'positive-massless-photon',photonE>0 and max(F(),photonE-energyLimit)**2<=photonNorm2<=(photonE+energyLimit)**2)
            rest=row['source_threshold_j']+row['source_mean_excess_j'];doppler=photonE-dot(photon,pin);beta=1-dot(pin,pin)/(c*c)
            restLimit=a+r*abs(rest)
            check(parent,'parent-rest-source-photon-energy',rest>0 and doppler>0 and beta>0 and max(F(),rest-restLimit)**2*beta<=doppler*doppler<=(rest+restLimit)**2*beta)
            bathMass=e['mass'];bathPin=[bathMass*ein[k]-photon[k] for k in range(3)]
            bathEin=kinetic(bathMass,ein)+e['internal']-photonE
            rotatedBath=mat(rotation,bathPin)
            bathP=[rotatedBath[k]+bathMass*offset[k] for k in range(3)]
            bathE=bathEin+dot(offset,rotatedBath)+kinetic(bathMass,offset)
            close(parent,'independent-signed-bath-mass',row['bath_mass_kg'],bathMass,'mass_kg')
            close(parent,'independent-signed-bath-energy',row['bath_event_j'],bathE,'energy_j')
            trackedP=[i['mass']*iv[k]-p['mass']*pv[k] for k in range(3)]
            trackedE=kinetic(i['mass'],iv)+i['internal']-kinetic(p['mass'],pv)-p['internal']
            close(parent,'per-event-additive-mass-closure',i['mass']-p['mass']+bathMass,F(),'mass_kg',max(i['mass'],p['mass'],bathMass))
            close(parent,'per-event-energy-closure',trackedE+bathE,F(),'energy_j',max(abs(trackedE),abs(bathE),abs(kinetic(p['mass'],pv)+p['internal'])))
            actualN=row['actual_ion_reference_weight']*row['actual_ion_individual_correction']
            weightedE=actualN*(kinetic(i['mass'],actualv)+row['actual_ion_internal_j'])-n*(kinetic(p['mass'],pv)+row['parent_internal_j'])
            close(parent,'represented-energy-closure',weightedE+n*bathE,F(),'energy_j',max(abs(weightedE),abs(n*bathE)))
            close(parent,'ledger-parent-mass',row['ledger_parent_mass_kg'],n*p['mass'],'mass_kg')
            close(parent,'ledger-ion-mass',row['ledger_ion_mass_kg'],actualN*i['mass'],'mass_kg')
            close(parent,'ledger-bath-mass',row['ledger_bath_mass_kg'],n*bathMass,'mass_kg')
            close(parent,'ledger-tracked-energy',row['ledger_tracked_de_j'],weightedE,'energy_j')
            close(parent,'ledger-bath-energy',row['ledger_bath_de_j'],n*bathE,'energy_j')
            ledgerP,ledgerBath=vector(row,'ledger_tracked_dp'),vector(row,'ledger_bath_dp')
            declaredBath=vector(row,'bath_event_p')
            for axis in range(3):
                close(parent,'independent-signed-bath-momentum-'+str(axis),declaredBath[axis],bathP[axis],'momentum_kg_m_s')
                close(parent,'per-event-momentum-closure-'+str(axis),trackedP[axis]+bathP[axis],F(),'momentum_kg_m_s',max(abs(i['mass']*iv[axis]),abs(p['mass']*pv[axis]),abs(bathP[axis])))
                weightedP=actualN*i['mass']*actualv[axis]-n*p['mass']*pv[axis]
                close(parent,'represented-momentum-closure-'+str(axis),weightedP+n*bathP[axis],F(),'momentum_kg_m_s',max(abs(weightedP),abs(n*bathP[axis])))
                close(parent,'ledger-tracked-momentum-'+str(axis),ledgerP[axis],weightedP,'momentum_kg_m_s')
                close(parent,'ledger-bath-momentum-'+str(axis),ledgerBath[axis],n*bathP[axis],'momentum_kg_m_s')
            check(parent,'tracked-integer-charge-change',row['tracked_charge_change_e']==i['charge']-p['charge'])
            check(parent,'bath-integer-charge-change',row['bath_charge_change_e']==e['charge'])
            for atom in p['nuclei']:
                close(parent,'represented-nuclei-'+atom,actualN*i['nuclei'].get(atom,0)+n*e['nuclei'].get(atom,0),n*p['nuclei'][atom],'number')
            close(parent,'represented-elementary-charge',actualN*i['charge']+n*e['charge'],n*p['charge'],'number')
            check(parent,'nonnegative-exposure-storage-error',row['exposure_error_bound']>=0 and row['clock_storage_error_bound']>=0)
        if not events:
            check('*','valid-empty-event-set-identity',True)
            result['event_validation_evaluated']=False
            result['empty_event_scope']='empty event set only; owning audit must reconcile all original parent histories'
        else:result['event_validation_evaluated']=True
        result['events_evaluated']=len(events)
        for row in result['checks']:result['counts'][row['status']]+=1
        result['status']='FAIL' if result['counts']['FAIL'] else 'PASS'
    except Unavailable as error:
        result.update(reason=str(error));result['checks'].append({'original_parent_id':'*','criterion':'required-evidence-availability','status':'SKIPPED','residual_exact':'0/1','limit_exact':'0/1'})
    except (OSError,ValueError,KeyError,TypeError,OverflowError,InvalidOperation) as error:
        result.update(reason=str(error));result['checks'].append({'original_parent_id':'*','criterion':'input-infrastructure','status':'ERROR','residual_exact':'0/1','limit_exact':'0/1'})
    result['counts']={s:sum(row['status']==s for row in result['checks']) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    result['status']=next(s for s in ('ERROR','FAIL','SKIPPED','PASS') if result['counts'][s])
    return result


def write_artifacts(result,output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    source=Path(__file__).resolve();repo=source.parents[3]
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=repo,capture_output=True,text=True).stdout.strip()
    manifest={'schema':'srcMoon-M3-independent-event-ledger-manifest-v1','created_utc':datetime.now(timezone.utc).isoformat(),'source_revision':revision,'source_sha256':{str(source):digest(source)},'input_sha256':result.get('input_sha256',{}),'command':list(sys.argv),'cwd':str(Path.cwd()),'scope':result['scope'],'qualification':'supporting reconstruction only; actual native source/input/executable/predecessor proof remains caller-owned'}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if not result['checks']:return
    with (output/'comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(result['checks'][0]));writer.writeheader();writer.writerows(result['checks'])
    # Plot only evaluated ledger residuals; this is never a generated
    # observational curve or an alteration to the exact CSV comparisons.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axis=plt.subplots(figsize=(12,4));x=list(range(len(result['checks'])))
    # Binary pass/fail markers stay meaningful for exact-zero bounds and
    # arbitrary rational exponents; plots never narrow the science arithmetic.
    outcomes=[1 if row['status']=='PASS' else 0 for row in result['checks']]
    axis.scatter(x,outcomes,s=12,label='independently evaluated criterion')
    axis.set_yticks([0,1]);axis.set_yticklabels(['not met / unavailable','met'])
    axis.set_xlabel('Evaluated event criterion (exact residuals and SI limits in CSV)');axis.set_ylabel('Criterion outcome');axis.legend(framealpha=1)
    fig.suptitle('Independent event bookkeeping — no observational validation');fig.tight_layout()
    for suffix in ('png','eps'):fig.savefig(output/('comparison.'+suffix))
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--events',nargs='+',required=True,type=Path);p.add_argument('--profile',required=True,type=Path);p.add_argument('--mapping',required=True,type=Path);p.add_argument('--bounds',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path)
    a=p.parse_args()
    if a.output_dir.exists():p.error('output must not exist; existing evidence is never overwritten')
    try:r=verify(a.events,a.profile,strict_json(a.mapping),strict_json(a.bounds))
    except (OSError,ValueError,TypeError,InvalidOperation) as error:
        r={'schema':'srcMoon-M3-independent-event-ledger-v1','status':'ERROR','scope':'input infrastructure only','formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','counts':{'PASS':0,'FAIL':0,'ERROR':1,'SKIPPED':0},'checks':[],'input_sha256':{},'reason':str(error)}
    try:write_artifacts(r,a.output_dir)
    except (OSError,ImportError,ValueError,OverflowError) as error:
        r.update(status='ERROR',reason='artifact generation failed: '+str(error));r['counts']['ERROR']+=1
        if a.output_dir.is_dir():(a.output_dir/'result.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
