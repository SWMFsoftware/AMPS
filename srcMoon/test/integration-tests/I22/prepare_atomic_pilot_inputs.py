#!/usr/bin/env python3
"""Prepare the source-bound inputs of the I22 atomic photoionization pilot.

Decision 2, part A (approved 2026-11-26; reference/atomic_pilot_approval_20261126.json):
He, Ne and Ar photoionization X -> X+ + e- to the ground ion under the named
source-effective-photon convention. For each parent this writes the same
eight-column energy crosswalk as the I21 sodium crosswalk, and a derived
state profile that lists the three crosswalks beside the unchanged Na
``atomic_extension``. Nothing physical is chosen here:

* the states, masses and formation-energy sources come from the pinned
  atomic-family profile (ground-ion references from CODATA 2018 for 4He,
  Eikema et al. 1994 for 20Ne and Velchev et al. 1999 for 40Ar);
* the channel, printed threshold and endpoint rows come from the authenticated
  D04 PHIDRATES package, whose QA must be PASS and whose bytes must match its
  provenance;
* the photon energy is the D04 printed threshold plus the D04 endpoint mean
  excess (applied by the production provider); the products use the physical
  ground-ion Q. Their difference is recorded, never corrected.

Inputs: --profile (atomic-family profile.json) and its expected --profile-sha256;
--output-dir must not exist. Outputs: <X>_source_energy_crosswalk.{json,csv},
profile.json, result.json, run_manifest.json. Status: PASS 0, ERROR 2 (changed
bytes, contradictory or malformed source), SKIPPED 77 (a required source or
state is absent). This verifies input construction only; formal I22 stays
SKIPPED / NOT QUALIFIED.
"""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
APPROVAL = HERE / 'reference/atomic_pilot_approval_20261126.json'
D04 = Path('/data/vtenishe/moon_validation_data/photochemistry')
TABLES = {'thresholds': 'processed/phidrates/thresholds.csv', 'rates': 'processed/phidrates/photo_rates_1au.csv'}
# D04 parent -> isotope of the approved pilot; the ion and product follow.
PILOT = {'He': '4He', 'Ne': '20Ne', 'Ar': '40Ar'}
ATOMIC_SCHEMA = 'srcMoon-M3-atomic-source-energy-crosswalk-v1'
COLUMNS = ('channel_id', 'parent_state_id', 'ion_state_id', 'electron_state_id', 'energy_reference_identity',
           'photon_energy_policy', 'source_threshold_eV', 'physical_threshold_source_identity')
EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}


class Missing(Exception):
    """A required source, state or capability is absent (SKIPPED)."""


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strict(path):
    def pairs(items):
        keys = [k for k, _ in items]
        if len(keys) != len(set(keys)):
            raise ValueError('duplicate JSON key in ' + str(path))
        return dict(items)
    def constant(token):
        raise ValueError('nonfinite JSON number in ' + str(path))
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs, parse_constant=constant)


def pin(record):
    """Authenticate a {path, sha256} descriptor; absent file is SKIPPED."""
    path = Path(record['path'])
    if not path.is_file():
        raise Missing('pinned source unavailable: ' + str(path))
    if sha(path) != record['sha256']:
        raise ValueError('pinned source bytes changed: ' + str(path))
    return {'path': str(path.resolve()), 'sha256': record['sha256']}


def rows(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def d04_tables():
    """Return authenticated D04 table pins after the AGENTS.md package checks."""
    for name in ('README.md', 'provenance.json', 'qa_report.json'):
        if not (D04 / name).is_file():
            raise Missing('D04 package metadata unavailable: ' + str(D04 / name))
    if strict(D04 / 'qa_report.json').get('status') != 'PASS':
        raise Missing('D04 QA is not PASS')
    outputs = {item['path']: item['sha256'] for item in strict(D04 / 'provenance.json')['outputs']}
    pins = {}
    for key, relative in TABLES.items():
        if relative not in outputs:
            raise ValueError('D04 provenance does not list ' + relative)
        pins[key] = pin({'path': str(D04 / relative), 'sha256': outputs[relative]})
    return pins


def prepare(profile_path, expected, output):
    profile_path = Path(profile_path).resolve()
    if not profile_path.is_file():
        raise Missing('atomic-family profile unavailable: ' + str(profile_path))
    if sha(profile_path) != expected:
        raise ValueError('atomic-family profile SHA-256 differs from the explicit expected value')
    profile = strict(profile_path)
    if profile.get('schema') != 'srcMoon-M3-state-profile-v1' or profile.get('physical_activation') is not False:
        raise ValueError('source-only srcMoon-M3-state-profile-v1 required')
    if 'atomic_photoionization_crosswalks' in profile:
        raise ValueError('input profile already lists atomic photoionization crosswalks')
    sources = {}
    for item in profile['sources']:
        if item['id'] in sources:
            raise ValueError('duplicate profile source')
        sources[item['id']] = {**pin(item), 'role': item['role']}
    states = {}
    for row in rows(pin(profile['states'])['path']):
        if row['state_id'] in states:
            raise ValueError('duplicate profile state')
        states[row['state_id']] = row
    constants = {row['name']: row for row in rows(pin(profile['constants'])['path'])}
    charge = Decimal(constants['elementary charge']['value'])
    if charge != Decimal('1.602176634E-19') or constants['elementary charge']['unit'] != 'C':
        raise ValueError('profile elementary charge is not the exact SI value')
    if 'DERIVED-GROUND-ION-REFERENCE' not in sources:
        raise Missing('ground-ion reference derivation is not a profile source')
    derivation = strict(sources['DERIVED-GROUND-ION-REFERENCE']['path'])
    tables = d04_tables()
    thresholds, rates = rows(tables['thresholds']['path']), rows(tables['rates']['path'])
    approval = {'path': str(APPROVAL), 'sha256': sha(APPROVAL)}

    def identity_of(state, component):
        # Same construction as PhysicalStateProfile: id@sha256:<source digest>#record.
        sid = state[component + '_source_id']
        if sid not in sources:
            raise ValueError('state source is not in the profile: ' + sid)
        return sid + '@sha256:' + sources[sid]['sha256'] + '#' + state[component + '_source_record']

    output.mkdir(parents=True, exist_ok=False)
    crosswalks, summary = {}, {}
    for element, isotope in PILOT.items():
        ids = {'parent': isotope + ':ground', 'ion': isotope + '+:ground', 'electron': 'electron:free'}
        chosen = {}
        for key, sid in ids.items():
            state = states.get(sid)
            if state is None or state['status'] != 'available':
                raise Missing('required state unavailable: ' + sid)
            chosen[key] = state
        parent, ion, electron = chosen['parent'], chosen['ion'], chosen['electron']
        if [(parent['nuclei'], parent['charge']), (ion['nuclei'], ion['charge']), (electron['nuclei'], electron['charge'])] \
                != [(isotope + ':1', '0'), (isotope + ':1', '1'), ('', '-1')]:
            raise ValueError(element + ': parent/ion/electron topology differs')
        if parent['state_kind'] != 'ground' or ion['state_kind'] != 'ground' or electron['state_kind'] != 'free-electron':
            raise ValueError(element + ': approved scope is the ground parent, ground ion and free electron')
        reference = identity_of(parent, 'formation')
        for state in (parent, electron):
            if Decimal(state['formation_J']) != 0 or identity_of(state, 'formation') != reference:
                raise ValueError(element + ': parent and electron must share the zero energy reference')
        for state in (parent, ion, electron):
            if Decimal(state['excitation_J']) != 0 or identity_of(state, 'excitation') != reference:
                raise ValueError(element + ': nonzero or differently referenced excitation')
        if ion['formation_source_id'] != 'GROUND-ION-' + isotope:
            raise ValueError(element + ': ion formation energy is not the retained ground-ion source')
        threshold_rows = [r for r in thresholds if r['species'] == element]
        family = [r for r in rates if r['species'] == element]
        products = element + '+ + e'
        if len(threshold_rows) != 1 or threshold_rows[0]['products'] != products:
            raise ValueError(element + ': D04 threshold family is not one X -> X+ + e channel')
        channel = threshold_rows[0]['reaction_id']
        if (len(family) != 2 or {r['reaction_id'] for r in family} != {channel}
                or {r['products'] for r in family} != {products} or {r['reaction_type'] for r in family} != {'photoionization'}
                or sorted(r['solar_activity'] for r in family) != ['0.0', '1.0']
                or {r['optical_depth'] for r in family} != {'0.0'} or {r['distance_au'] for r in family} != {'1.0'}):
            raise ValueError(element + ': D04 rate family is not the complete quiet/active 1 AU endpoint pair')
        recorded = [r for r in derivation if r.get('isotope') == isotope]
        if len(recorded) != 1 or Decimal(recorded[0]['D04_threshold']['threshold_eV']) != Decimal(threshold_rows[0]['threshold_eV']) \
                or recorded[0]['channel'] != channel or Decimal(recorded[0]['physical_threshold_J']) != Decimal(ion['formation_J']):
            raise ValueError(element + ': retained ground-ion derivation disagrees with the profile or D04')
        with localcontext() as context:
            context.prec = 60
            physical_ev = Decimal(ion['formation_J']) / charge
            difference = physical_ev - Decimal(threshold_rows[0]['threshold_eV'])
        meta = {'schema': ATOMIC_SCHEMA, 'channel_id': channel,
                'parent_state_id': ids['parent'], 'ion_state_id': ids['ion'], 'electron_state_id': ids['electron'],
                'energy_reference_identity': reference, 'photon_energy_policy': 'source-effective-photon',
                'source_threshold_eV': threshold_rows[0]['threshold_eV'],
                'physical_threshold_eV': str(physical_ev), 'physical_threshold_J': ion['formation_J'],
                'physical_threshold_source_identity': identity_of(ion, 'formation'),
                'physical_minus_source_threshold_eV': str(difference),
                'available_translational_energy_correction_eV': str(-difference),
                'effective_photon_energy_equation': 'Egamma=source_threshold+source_mean_excess; physical Q/recoil determines available translation',
                'solar_endpoints': family, 'source_files': [tables['thresholds'], tables['rates']],
                'derivation': {k: sources['DERIVED-GROUND-ION-REFERENCE'][k] for k in ('path', 'sha256')},
                'approval': approval,
                'product_state_scope': 'every event to the ground ion; admitted observables are heavy number, charge and mass budgets and event chronology only',
                'absolute_rate_bias_bound': 'UNAVAILABLE: requires original spectral/cross-section convolution; no rounding interval assumed',
                'source_rates_modified': False, 'physical_activation': False}
        meta_path = output / (element + '_source_energy_crosswalk.json')
        meta_path.write_text(json.dumps(meta, indent=2) + '\n')
        csv_path = output / (element + '_source_energy_crosswalk.csv')
        with csv_path.open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(COLUMNS)
            writer.writerow([meta[k] for k in COLUMNS])
        crosswalks[element] = {'crosswalk': {'path': str(meta_path), 'sha256': sha(meta_path)},
                               'crosswalk_csv': {'path': str(csv_path), 'sha256': sha(csv_path)}}
        summary[element] = {'channel_id': channel, 'source_threshold_eV': meta['source_threshold_eV'],
                            'physical_threshold_eV': meta['physical_threshold_eV'],
                            'physical_minus_source_threshold_eV': meta['physical_minus_source_threshold_eV']}
    derived = dict(profile)
    derived.update(
        atomic_photoionization_crosswalks=crosswalks,
        derived_from={'path': str(profile_path), 'sha256': expected},
        recipe={'path': str(Path(__file__).resolve()), 'sha256': sha(__file__)},
        approval=approval,
        profile_scope=profile.get('profile_scope', '') + '; I22 atomic-pilot He/Ne/Ar source-effective-photon '
                      'crosswalks (approved 2026-11-26); no physical channel qualification')
    (output / 'profile.json').write_text(json.dumps(derived, indent=2) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--profile', required=True, type=Path)
    parser.add_argument('--profile-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        print('ERROR output exists; preserve earlier artifacts: ' + str(output), file=sys.stderr)
        return 2
    result = {'schema': 'srcMoon-I22-atomic-pilot-inputs-v1', 'status': 'ERROR',
              'scope': 'I22 atomic-pilot input construction only; no particle run',
              'formal_I22_status': 'SKIPPED / NOT QUALIFIED'}
    try:
        result['crosswalks'] = prepare(args.profile, args.profile_sha256, output)
        result.update(status='PASS', reason='three source-bound crosswalks and the derived profile written')
    except Missing as error:
        result.update(status='SKIPPED', reason=str(error))
    except (OSError, ValueError, KeyError, ArithmeticError) as error:
        result.update(status='ERROR', reason=str(error))
    output.mkdir(parents=True, exist_ok=True)
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'source_revision': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True,
                                                  text=True).stdout.strip(),
                'command': [sys.executable] + sys.argv, 'script': {'path': str(Path(__file__).resolve()), 'sha256': sha(__file__)},
                'outputs': {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()}}
    (output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return EXIT[result['status']]


if __name__ == '__main__':
    sys.exit(main())
