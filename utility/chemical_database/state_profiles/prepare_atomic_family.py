#!/usr/bin/env python3
"""Prepare original-research 4He, 20Ne and 40Ar ground-ion state profiles."""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path
import shutil
import subprocess

from prepare import audit_canonical, digest, strict_json
from run import arguments

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LEDGER_SHA256 = 'ac784bd1b8fdd9c6f44f99733c7ee386c88dc2d7098b3154ef094186c05ad834'
ELEMENTS = {'4He': 'He', '20Ne': 'Ne', '40Ar': 'Ar'}
EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}


def pin(path, expected, identities):
    path = path.resolve()
    actual = digest(path)
    if actual != expected:
        raise ValueError('source identity mismatch: '+str(path))
    identities[str(path)] = actual
    return path


def read_rows(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        rows=list(reader)
        if (not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames)
                or any(None in row or None in row.values() for row in rows)):
            raise ValueError('malformed CSV: '+str(path))
        return rows, reader.fieldnames


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def prepare(args, manifest):
    """Reuse existing state/constant rows and append only supported grounds."""
    identities = manifest['input_sha256']
    ledger = strict_json(pin(args.threshold_ledger, LEDGER_SHA256, identities))
    profile_path = pin(args.base_profile, args.base_profile_sha256, identities)
    profile = strict_json(profile_path)
    if profile['schema'] != 'srcMoon-M3-state-profile-v1' or profile['physical_activation'] is not False:
        raise ValueError('explicit source research profile required')
    for row in profile['sources']:
        row['path'] = str(pin(profile_path.parent/row['path'], row['sha256'], identities))
    files = {key: pin(profile_path.parent/value['path'], value['sha256'], identities)
             for key, value in profile.items() if key in ('isotopes','constants','states')}
    states, columns = read_rows(files['states'])
    by_state = {row['state_id']: row for row in states}
    if len(by_state) != len(states):
        raise ValueError('duplicate base-profile states')
    constants = {row['name']: row['value'] for row in read_rows(files['constants'])[0]}
    d04 = audit_canonical(args.d04)
    for row in d04['metadata']+d04['verified_processed_products']:
        pin(Path(row['path']),row['sha256'],identities)
    output = args.output_dir.resolve()
    write_json(output/'d04_audit.json', d04)
    table = args.d04/'processed/phidrates'
    thresholds = {row['reaction_id']: row for row in read_rows(table/'thresholds.csv')[0]}
    rates = read_rows(table/'photo_rates_1au.csv')[0]
    crosswalk = []
    photon_candidates = []
    with localcontext() as context:
        context.prec = 90
        electron = Decimal(constants['electron mass'])
        e = Decimal(constants['elementary charge'])
        c = Decimal(constants['speed of light in vacuum'])
        for source in ledger['isotope_mappings']:
            isotope = source['isotope']
            if isotope not in ELEMENTS or (source['parent_charge_e'],source['ion_charge_e']) != (0,1):
                raise ValueError('unexpected isotope or charge in threshold ledger')
            original = pin(Path(source['source_path']), source['source_sha256'], identities)
            neutral = by_state[isotope+':ground']
            ion = isotope+'+:ground'
            if (ion in by_state or neutral['status'] != 'available' or neutral['nuclei'] != isotope+':1'
                    or neutral['charge'] != '0' or Decimal(neutral['formation_J']) != 0
                    or Decimal(neutral['excitation_J']) != 0):
                raise ValueError('base state cannot define the named ground reference')
            conversion = source['conversion']
            if (Decimal(conversion['h_J_s']) != Decimal('6.62607015e-34')
                    or Decimal(conversion['c_m_s']) != c or Decimal(conversion['e_C']) != e):
                raise ValueError('retained exact SI conversion disagrees with the base constants')
            factor = 100*Decimal(conversion['h_J_s'])*c
            energy = Decimal(source['threshold_cm_inverse'])*factor
            reference_mass = Decimal(neutral['physical_rest_mass_kg'])-electron
            physical_mass = reference_mass+energy/(c*c)
            channel = source['native_D04_reaction_id']
            products = ELEMENTS[isotope]+'+ + e'
            endpoint = [row for row in rates if row['reaction_id']==channel]
            native_threshold = thresholds[channel]
            if (Decimal(native_threshold['threshold_eV']) != Decimal(source['native_D04_threshold_eV'])
                    or native_threshold['products'] != products or len(endpoint)!=2
                    or {row['solar_activity'] for row in endpoint}!={'0.0','1.0'}
                    or {row['reaction_id'] for row in rates if row['species']==ELEMENTS[isotope]}!={channel}
                    or any((row['products'],row['distance_au'],row['optical_depth'])!=(products,'1.0','0.0')
                           for row in endpoint)):
                raise ValueError('D04 endpoint or threshold differs from source record')
            source_id = 'GROUND-ION-'+isotope
            profile['sources'].append({'id': source_id,'path':str(original),
                'sha256':source['source_sha256'],'role':'energy'})
            crosswalk.append({'isotope':isotope,'parent_state':neutral['state_id'],'ion_state':ion,
                'primary':source,'physical_threshold_J':str(energy),
                'reported_threshold_uncertainty_J':str(Decimal(source['reported_uncertainty_cm_inverse'])*factor),
                'physical_threshold_eV':str(energy/e),'reference_mass_kg':str(reference_mass),
                'physical_mass_kg':str(physical_mass),'channel':channel,'D04_threshold':native_threshold,
                'physical_minus_D04_threshold_eV':str(energy/e-Decimal(native_threshold['threshold_eV'])),
                'D04_endpoints':endpoint,'ground_ion_yield':'UNAVAILABLE'})
            # This worksheet aligns the printed photon threshold/mean with
            # one possible physical product state. It never supplies a yield,
            # an approval, or a substitute photon distribution to the loader.
            threshold_ev = Decimal(native_threshold['threshold_eV'])
            for row in endpoint:
                rate = Decimal(row['rate_1au_s-1'])
                mean_ev = Decimal(row['excess_energy_eV'])
                if not rate.is_finite() or rate<=0 or not mean_ev.is_finite() or mean_ev<0:
                    raise ValueError('invalid frozen D04 rate or mean excess')
                photon_j = (threshold_ev+mean_ev)*e
                photon_candidates.append({
                    'isotope':isotope,'parent_state':neutral['state_id'],'ion_state':ion,
                    'source_channel':channel,'solar_activity':row['solar_activity'],
                    'rate_1au_s-1':str(rate),'source_mean_excess_eV':str(mean_ev),
                    'printed_threshold_eV':str(threshold_ev),'physical_threshold_J':str(energy),
                    'candidate_effective_photon_J':str(photon_j),
                    'energy_above_physical_threshold_J':str(photon_j-energy),
                    'energy_admissible_before_recoil':photon_j>=energy,
                    'native_rate_record':row['source_file']+'#'+row['source_record_or_index'],
                    'ground_ion_yield_status':'UNAVAILABLE','transport_status':'SKIPPED'})
            # A ground ion has zero excitation on the same reference as its
            # ground neutral. Its formation energy and physical mass carry
            # the primary source locators and are independently traceable.
            states.append({**neutral,'state_id':ion,'state_label':isotope+'+ '+source['ion_ground_state'],
                'charge':'1','formation_J':str(energy),'formation_source_id':source_id,
                'formation_source_record':source['source_locator'],'physical_rest_mass_kg':str(physical_mass),
                'physical_mass_source_id':'DERIVED-GROUND-ION-REFERENCE',
                'physical_mass_source_record':isotope+'/physical_mass_kg'})
    if len(crosswalk)!=len(ELEMENTS) or {row['isotope'] for row in crosswalk}!=set(ELEMENTS):
        raise ValueError('missing or duplicate ground-ion threshold records')
    path = output/'states.csv'
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(states)
    files['states'] = path
    for name, path in files.items():
        profile[name] = {'path':str(path),'sha256':digest(path)}
    profile.update(profile_scope='reference ground states only', recipe={'path':str(Path(__file__).resolve()),'sha256':digest(__file__)},
        available_state_count=sum(row['status']=='available' for row in states),
        unavailable_state_count=sum(row['status']=='unavailable' for row in states))
    crosswalk_path=output/'ground_ion_crosswalk.json'
    write_json(crosswalk_path, crosswalk)
    write_json(output/'source_effective_photon_candidates.json',{
        'schema':'amps-atomic-source-effective-photon-review-v1',
        'energy_convention':'D04-printed-threshold-plus-matched-endpoint-mean',
        'scope':'state-crosswalk preparation; ground-ion yield and observable bounds require source evidence',
        'physical_activation':False,'candidates':photon_candidates})
    profile['sources'].append({'id':'DERIVED-GROUND-ION-REFERENCE','path':str(crosswalk_path),
        'sha256':digest(crosswalk_path),'role':'energy'})
    return profile, crosswalk


def verify(profile, crosswalk, output, manifest):
    """Call the maintained production resolver with the actual prepared bytes."""
    compiler = shutil.which('g++')
    if compiler is None:
        raise FileNotFoundError('g++ is unavailable')
    sources = [HERE/'probe.cpp', ROOT/'srcMoon/species/PhysicalStateProfile.cpp',
               ROOT/'srcMoon/species/ReducedKinematics.cpp', ROOT/'src/models/chemical_database/DataIO.cpp']
    manifest['source_sha256'] = {str(path):digest(path) for path in [Path(__file__).resolve(), *sources]}
    exe = output/'native-probe'
    compile_command = [compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',
                       '-I',str(ROOT),*[str(path) for path in sources],'-o',str(exe)]
    def launch(command, label):
        path = output/(label+'.log')
        manifest['commands'].append({'argv':command,'log':str(path)})
        with path.open('x') as stream:
            return subprocess.run(command, cwd=ROOT, timeout=120, stdout=stream, stderr=subprocess.STDOUT)
    if launch(compile_command,'compile').returncode != 0:
        raise RuntimeError('ground-ion profile probe build failed')
    manifest['executable_sha256'] = digest(exe)
    status = launch(arguments(exe,profile),'native').returncode
    if status not in (0,1):
        raise RuntimeError('production profile resolver failed')
    rows = [line.split('\t') for line in (output/'native.log').read_text().splitlines()]
    checks = [{'name':row[1],'status':row[2]} for row in rows if row[0]=='CHECK']
    summary=next(row for row in rows if row[0]=='SUMMARY')
    if (int(summary[1])!=len(checks) or int(summary[2])!=sum(row['status']=='FAIL' for row in checks)
            or status!=int(int(summary[2])>0)):
        raise ValueError('native resolver check records contradict its summary or exit')
    resolved = {row[1]:row[2:] for row in rows if row[0]=='STATE'}
    masses = {row[1]:row[2] for row in rows if row[0]=='PHYSICAL'}
    epsilon = Decimal(next(row[1] for row in rows if row[0]=='EPSILON'))
    for row in crosswalk:
        ion = row['ion_state']
        if ion not in resolved or ion not in masses:
            checks.append({'name':ion+'/actual-resolved-state-and-mass','status':'FAIL'})
            continue
        mass, charge, formation, excitation = resolved[ion]
        for quantity, actual, expected in [('reference_mass_kg',mass,row['reference_mass_kg']),
                ('physical_mass_kg',masses[ion],row['physical_mass_kg']),
                ('formation_J',formation,row['physical_threshold_J'])]:
            with localcontext() as context:
                context.prec=90
                error=abs(Decimal(actual)/Decimal(expected)-1)
            checks.append({'name':ion+'/'+quantity,'status':'PASS' if error<=2*epsilon else 'FAIL',
                'actual':actual,'expected':expected,'relative_error':str(error),'limit':str(2*epsilon)})
        checks.append({'name':ion+'/charge-and-ground-state','status':'PASS' if charge=='1' and Decimal(excitation)==0 else 'FAIL'})
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('base-profile','threshold-ledger','d04','output-dir'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--base-profile-sha256',required=True)
    args=parser.parse_args()
    output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    manifest={'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'input_sha256':{},'commands':[],'physical_activation':False}
    try:
        profile, crosswalk=prepare(args,manifest)
        checks=verify(profile,crosswalk,output,manifest)
        status='FAIL' if any(row['status']=='FAIL' for row in checks) else 'PASS'
        if status=='PASS':
            write_json(output/'profile.json',profile)
    except FileNotFoundError as error:
        status='SKIPPED';checks=[{'name':'prerequisite','status':status,'reason':str(error)}]
    except (OSError,ValueError,ArithmeticError,KeyError,IndexError,StopIteration,RuntimeError,subprocess.TimeoutExpired) as error:
        status='ERROR';checks=[{'name':'infrastructure','status':status,'reason':str(error)}]
    result={'status':status,'scope':'ground-ion reference state preparation', 'formal_M3_status':'SKIPPED / NOT QUALIFIED',
            'counts':{status:sum(row['status']==status for row in checks) for status in EXIT},'checks':checks}
    write_json(output/'run_manifest.json',manifest)
    write_json(output/'result.json',result)
    print(json.dumps({'status':status,'counts':result['counts'],'result':str(output/'result.json')}))
    return EXIT[status]


if __name__=='__main__':
    raise SystemExit(main())
