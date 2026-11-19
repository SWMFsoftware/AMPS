#!/usr/bin/env python3
"""Opt-in source-bound Na I/Na II original-research profile preparation.

This never changes the old partial preparation, canonical data or an installed
database. Original ASD research use is distinct from the unresolved managed-
storage grant. No unsupported ion/molecule is inferred and no physical source
or observational campaign is activated. Read ATOMIC_EXTENSION_CONTRACT.md.
"""
from __future__ import annotations
import argparse
import copy
import csv
from decimal import Decimal as D, localcontext
import hashlib
from html import unescape
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import prepare as base

ROOT=Path(__file__).resolve().parents[3]
ENERGY_SHA='4f5750dbc6fb526406b31c5eb18c5994b706a2e9961b3eb3b5498dfe5f4278fc'
LABEL_SHA='dcdea4fbf9c06df35d0551a60d678113428f917dc3ebe720af380f0f9d92e35c'
HELP_SHA='1963ba2edf436e66bbd118cf15ba84643ce127b27ffa4b8134c4b27228604b5f'
CROSSWALK_COLUMNS=['channel_id','parent_state_id','ion_state_id','electron_state_id',
 'energy_reference_identity','photon_energy_policy','source_threshold_eV','physical_threshold_source_identity']


class ArithmeticMismatch(RuntimeError):
    pass


def _native(path,expected):
    path=Path(path).resolve()
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=expected:
        raise ValueError('authoritative original source SHA changed: '+str(path))
    return data.decode('utf-8')


def parse_sodium(energy_file,label_file,help_file):
    """Select the actual neutral first IE, not the next ion's ionization value.

    No brackets/theoretical flag, isotope population, unit or uncertainty kind
    is guessed. The SAME handbook explicitly documents23Na100% and both ground
    labels. Its old rounded masses/energies are deliberately not extracted.
    """
    energy=unescape(re.sub(r'<[^>]+>','',_native(energy_file,ENERGY_SHA)))
    labels=unescape(re.sub(r'<[^>]+>','',_native(label_file,LABEL_SHA)))
    help_text=unescape(re.sub(r'<[^>]+>','',_native(help_file,HELP_SHA)))
    if 'Ionization Energy (eV)'not in energy or 'Uncertainty (eV)'not in energy or '5.12'not in energy:
        raise ValueError('ASD version/energy/uncertainty units unavailable')
    selected=[]
    for line in energy.splitlines():
        fields=[x.strip()for x in line.split('|')]
        if len(fields)==8 and fields[0]=='11'and fields[2]=='Na':selected.append((line,fields))
    if len(selected)!=1:raise ValueError('neutral Na first ionization row missing/ambiguous')
    line,f=selected[0]
    if f[3]!='[Ne].3s'or f[4]!='2S<1/2>'or not re.fullmatch(r'\d+\.\d+',f[5])or not re.fullmatch(r'\d+\.\d+',f[6]):
        raise ValueError('native neutral Na state/value/uncertainty is unsupported')
    if not re.search(r'23Na[^\n]*100%',labels)or not re.search(r'Na I\s+Ground State\s+1s22s22p63s\s+2S1/2',labels)or not re.search(r'Na II\s*Ground State\s+1s22s22p6\s+1S0',labels):
        raise ValueError('explicit23Na abundance or neutral/ion ground-state mapping unavailable')
    if 'Uncertainty'not in help_text or 'conversion factors'not in help_text:
        raise ValueError('source uncertainty/units explanation unavailable')
    return {'value_eV':f[5],'reported_uncertainty_eV':f[6],
      'uncertainty_confidence':'not assigned a common sigma/confidence level by this preparation',
      'native_row':line,'native_references':f[7], 'native_neutral_configuration':f[3],
      'native_neutral_level':f[4],'ion_label':'Na II 1s2.2s2.2p6 1S0 electronic ground',
      'isotope_mapping':'23Na; SAME source handbook explicitly lists100% abundance and NaI/NaII ground labels',
      'energy_record':'ASD5.12:Z11;isoelectronic-sequence=Na;ground=[Ne].3s 2S1/2;units=eV;references='+f[7]}


def _reported_uncertainty(token):
    if token.strip()=='(exact)':return D(0)
    value=D(token.replace(' ',''))
    if not value.is_finite()or value<0:raise ValueError('native constant uncertainty unavailable')
    return value


def _write_csv(path,header,rows):
    with Path(path).open('w',newline='')as stream:
        writer=csv.writer(stream);writer.writerow(header);writer.writerows(rows)


def prepare_atomic(raw_dir,decision_register,energy_file,label_file,help_file,canonical_root,output_dir,work_dir):
    """Authenticate/derive fully before publishing a new research profile.

    The work directory is explicit/new and retains native compiler/evidence
    even on failure. The actual profile directory is created only after all
    authoritative source, semantic and independent arithmetic checks succeed.
    """
    output,work=Path(output_dir).resolve(),Path(work_dir).resolve()
    if output.exists()or work.exists():raise FileExistsError('profile/work output must be new')
    # Original-source use only. This code does not fetch/copy data into the
    # managed CHEMICAL_DATABASE and does not alter its source permission holds.
    sodium=parse_sodium(energy_file,label_file,help_file)
    isotopes,constants=base.extract_permitted(Path(raw_dir))
    na=[x for x in isotopes if x['identity']=='23Na']
    if len(na)!=1:raise ValueError('source23Na isotope missing/ambiguous')
    na=na[0]
    for name in ('elementary charge','speed of light in vacuum'):
        if constants[name]['native_uncertainty']!='(exact)':raise ValueError('native exactSI constant flag missing: '+name)
    audit=base.audit_canonical(canonical_root)
    canonical=Path(canonical_root).resolve()
    def csv_rows(name):
        p=canonical/'processed/phidrates'/name
        with p.open(newline='')as f:return p,list(csv.DictReader(f))
    threshold_path,thresholds=csv_rows('thresholds.csv')
    rate_path,rates=csv_rows('photo_rates_1au.csv')
    channel='PHIDRATES:Na:00:na_e'
    selected=[r for r in thresholds if r['reaction_id']==channel]
    endpoints=[r for r in rates if r['reaction_id']==channel]
    if len(selected)!=1 or len(endpoints)!=2 or {r['solar_activity']for r in endpoints}!={'0.0','1.0'}:
        raise ValueError('actual D04 Na threshold/quiet-active support unavailable')
    if any(r['products']!='Na+ + e'or r['optical_depth']!='0.0'or r['distance_au']!='1.0'for r in endpoints):
        raise ValueError('D04 Na product/normalization metadata changed')
    work.mkdir(parents=True,exist_ok=False)
    partial=base.prepare(raw_dir,decision_register,work/'partial',canonical)
    compiler=shutil.which('g++')
    if not compiler:raise FileNotFoundError('C++17 compiler unavailable')
    code=[Path(__file__).resolve(),Path(__file__).with_name('atomic_mass_probe.cpp').resolve(),
          ROOT/'src/models/chemical_database/ReducedModels.h',ROOT/'src/models/chemical_database/Runtime.h',
          ROOT/'src/models/chemical_database/DataIO.h',Path(base.__file__).resolve()]
    code_hashes={str(p):base.digest(p)for p in code}
    for p in code:
        target=work/'source_snapshot'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    executable=work/'atomic-mass-probe'
    command=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT),str(code[1]),'-o',str(executable)]
    (work/'compile_command.json').write_text(json.dumps(command,indent=2)+'\n')
    compiled=subprocess.run(command,text=True,capture_output=True,timeout=120)
    (work/'compile.stdout.log').write_text(compiled.stdout);(work/'compile.stderr.log').write_text(compiled.stderr)
    if compiled.returncode:raise RuntimeError('actual shared mass-kernel compile failed: '+str(work))
    with localcontext()as context:
        context.prec=90
        neutral=D(na['neutral_mass_kg']);electron=D(constants['electron mass']['value'])
        e=D(constants['elementary charge']['value']);c=D(constants['speed of light in vacuum']['value'])
        ie_ev=D(sodium['value_eV']);ie_j=ie_ev*e
        given=' '.join(map(str,(neutral,electron,ie_j,c)))+'\n'
        (work/'mass_input.txt').write_text(given)
        process=subprocess.run([str(executable)],input=given,text=True,capture_output=True,timeout=30)
        (work/'native.stdout.log').write_text(process.stdout);(work/'native.stderr.log').write_text(process.stderr)
        if process.returncode:raise RuntimeError('actual shared mass kernel refused authoritative inputs')
        native=dict(line.split('\t')for line in process.stdout.splitlines())
        if set(native)!={'reference_mass_kg','binding_mass_kg','physical_mass_kg','machine_epsilon'}:
            raise RuntimeError('shared mass adapter output contract changed')
        references={'reference_mass_kg':neutral-electron,'binding_mass_kg':ie_j/c/c,
                    'physical_mass_kg':neutral-electron+ie_j/c/c}
        comparisons=[]
        epsilon=D(native['machine_epsilon'])
        for key,value in references.items():
            relative=abs(D(native[key])-value)/abs(value)
            comparisons.append({'quantity':key,'native':native[key],'independent_Decimal':str(value),
                                'relative_error':str(relative),'limit':str(2*epsilon)})
            if relative>2*epsilon:raise ArithmeticMismatch('native/Decimal mass conversion criterion failed: '+key)
        match=re.fullmatch(r'(\d+\.\d+)\((\d+)\)',na['native_mass_with_uncertainty'])
        if not match:raise ValueError('native23Na mass uncertainty unavailable')
        ar=D(match[1]);uar=D(match[2]).scaleb(-len(match[1].split('.')[1]))
        mu=D(constants['atomic mass constant']['value'])
        contributions={'isotope_mass':mu*uar,
          'atomic_mass_constant':ar*_reported_uncertainty(constants['atomic mass constant']['native_uncertainty']),
          'electron_mass':_reported_uncertainty(constants['electron mass']['native_uncertainty']),
          'ionization_energy':D(sodium['reported_uncertainty_eV'])*e/c/c}
        source_threshold=selected[0]['threshold_eV'];delta=ie_ev-D(source_threshold)
        receipt={'schema':'srcMoon-M3-atomic-mass-derivation-v1','operation':'original-research-read-only',
          'sources':partial['sources']+[
            {'id':'NIST-ASD-NA-IE','path':str(Path(energy_file).resolve()),'sha256':ENERGY_SHA,'role':'energy'},
            {'id':'NIST-NA-STATE-LABELS','path':str(Path(label_file).resolve()),'sha256':LABEL_SHA,'role':'energy'},
            {'id':'NIST-ASD-IE-HELP','path':str(Path(help_file).resolve()),'sha256':HELP_SHA,'role':'energy'}],
          'sodium_source':sodium,'physical_ionization_energy_J':str(ie_j),
          'mass_equation':'mion=mneutral-me+IE/c^2; reference/binding contributions remain separate',
          'actual_shared_kernel':'AMPS::ChemicalDatabase::ReducedModels::IonPhysicalMass',
          'native':native,'independent_Decimal':{k:str(v)for k,v in references.items()},
          'source_uncertainty_contributions_kg':{k:str(v)for k,v in contributions.items()},
          'first_order_sum_of_source_reported_uncertainty_contributions_kg':str(sum(contributions.values())),
          'uncertainty_semantics':'linear first-order sum; no assumed independence, covariance, common confidence/sigma level or physical acceptance budget',
          'source_code_sha256':code_hashes,'compiler_command':command,
          'compiler':subprocess.check_output([compiler,'--version'],text=True).splitlines()[0],
          'executable_sha256':base.digest(executable),'independent_comparisons':comparisons,
          'physical_activation':False,'managed_database_modified':False,'original_sources_copied':False}
    if code_hashes!={str(p):base.digest(p)for p in code}:raise RuntimeError('source changed during preparation')
    output.mkdir(parents=True,exist_ok=False)
    receipt_path=output/'atomic_mass_derivation.json'
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n');receipt_sha=base.digest(receipt_path)
    profile=copy.deepcopy(partial)
    profile['sources']=receipt['sources']+[{'id':'DERIVED-NA-ATOMIC','path':str(receipt_path),'sha256':receipt_sha,'role':'energy'}]
    with Path(partial['states']['path']).open(newline='')as f:rows=list(csv.reader(f))
    for row in rows[1:]:
        if row[0]=='23Na:ground':row[1]='23Na I [Ne].3s 2S1/2 electronic ground'
    zero=['0','APPROVED-P02','P02:separated-neutral-ground-isotope-atoms-and-free-electron-reference-zero']
    rows.append(['23Na+:ground','23Na II 1s2.2s2.2p6 1S0 electronic ground','ground','23Na:1','1',
      receipt['physical_ionization_energy_J'],'DERIVED-NA-ATOMIC','physical_ionization_energy_J']+zero+
      [native['physical_mass_kg'],'DERIVED-NA-ATOMIC','native.physical_mass_kg','available',''])
    _write_csv(output/'states.csv',rows[0],rows[1:])
    for kind in ('isotopes','constants'):shutil.copyfile(partial[kind]['path'],output/(kind+'.csv'))
    for kind in ('isotopes','constants','states'):
        p=output/(kind+'.csv');profile[kind]={'path':str(p),'sha256':base.digest(p)}
    audit_path=output/'canonical_gap_audit.json';audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    profile['canonical_gap_audit']={'path':str(audit_path),'sha256':base.digest(audit_path)}
    energy_reference='APPROVED-P02@sha256:'+base.digest(decision_register)+'#'+zero[2]
    formation_identity='DERIVED-NA-ATOMIC@sha256:'+receipt_sha+'#physical_ionization_energy_J'
    crosswalk={'schema':'srcMoon-M3-Na-source-energy-crosswalk-v1','channel_id':channel,
      'parent_state_id':'23Na:ground','ion_state_id':'23Na+:ground','electron_state_id':'electron:free',
      'energy_reference_identity':energy_reference,'photon_energy_policy':'source-effective-photon',
      'source_threshold_eV':source_threshold,'physical_threshold_eV':sodium['value_eV'],
      'physical_threshold_J':receipt['physical_ionization_energy_J'],
      'physical_threshold_source_identity':formation_identity,
      'physical_minus_source_threshold_eV':str(delta),
      'available_translational_energy_correction_eV':str(-delta),
      'effective_photon_energy_equation':'Egamma=source_threshold+source_mean_excess; physical Q/recoil determines available translation',
      'solar_endpoints':endpoints,'source_files':[{'path':str(p),'sha256':base.digest(p)}for p in (threshold_path,rate_path)],
      'derivation':{'path':str(receipt_path),'sha256':receipt_sha},
      'absolute_rate_bias_bound':'UNAVAILABLE: requires original spectral/cross-section convolution; no rounding interval assumed',
      'source_rates_modified':False,'physical_activation':False}
    cross_path=output/'na_source_energy_crosswalk.json';cross_path.write_text(json.dumps(crosswalk,indent=2)+'\n')
    cross_csv=output/'na_source_energy_crosswalk.csv'
    _write_csv(cross_csv,CROSSWALK_COLUMNS,[[crosswalk[k]for k in CROSSWALK_COLUMNS]])
    profile.update(available_state_count=8,unavailable_state_count=14,
      profile_scope='opt-in source-bound Na I/Na II ground reference states plus unchanged partial inventory; no physical channel qualification',
      atomic_extension={'derivation':{'path':str(receipt_path),'sha256':receipt_sha},
        'crosswalk':{'path':str(cross_path),'sha256':base.digest(cross_path)},
        'crosswalk_csv':{'path':str(cross_csv),'sha256':base.digest(cross_csv)},
        'operation':'original-research-read-only; not a managed ASD import'},
      recipe={'path':str(Path(__file__).resolve()),'sha256':base.digest(__file__)},
      source_uncertainties_propagated=True,
      uncertainty_propagation_kind='first-order sum of source-reported contributions; no assumed common confidence/covariance',
      uncertainty_note='Native ASD uncertainty confidence is not assigned; the contribution sum is not a physical acceptance budget.')
    (output/'profile.json').write_text(json.dumps(profile,indent=2)+'\n')
    return profile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('raw-dir','decision-register','energy-file','label-file','help-file','canonical-root','output-dir','work-dir'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    try:value=prepare_atomic(a.raw_dir,a.decision_register,a.energy_file,a.label_file,a.help_file,a.canonical_root,a.output_dir,a.work_dir)
    except FileNotFoundError as e:print(json.dumps({'status':'SKIPPED','reason':str(e)}));return 77
    except ArithmeticMismatch as e:print(json.dumps({'status':'FAIL','reason':str(e)}));return 1
    except (OSError,ValueError,RuntimeError,KeyError,ArithmeticError,subprocess.TimeoutExpired)as e:
        print(json.dumps({'status':'ERROR','reason':str(e)}));return 2
    print(json.dumps({'status':'PASS','scope':'source-bound research preparation; no physical qualification','profile':value},indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
