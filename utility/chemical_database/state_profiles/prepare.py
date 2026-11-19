#!/usr/bin/env python3
"""Prepare a source-bound PARTIAL reference-state profile, without activation.

Only the two frozen, permission-cleared NIST isotope/CODATA raw responses are
used numerically. The approved P02 record supplies an energy-coordinate zero,
not measured atomic/molecular energies. Ion, excited and molecular entries
remain unavailable. No data are fetched or installed into a managed archive.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
from html import unescape
import json
from pathlib import Path
import re

SOURCE_HASHES = {
    'isotopes.html': '15de684d4be01faa70b0d72a19ddd5051cf330cf3c246bbcf36e423451a36983',
    'codata2022.txt': '77fb90e66c40db3e6eb16630bc9c88e4c7c8beddbe5e71be406f2f26e3f67e67',
}
ISOTOPES = (('H', 1), ('He', 4), ('O', 16), ('Ne', 20), ('Na', 23), ('Ar', 40))
CONSTANTS = {'atomic mass constant': 'kg', 'electron mass': 'kg',
             'elementary charge': 'C', 'speed of light in vacuum': 'm s^-1'}
HEADERS = {
    'isotopes': ['isotope', 'neutral_mass_kg', 'source_id', 'source_record'],
    'constants': ['name', 'value', 'unit', 'source_id', 'source_record'],
    'states': ['state_id', 'state_label', 'state_kind', 'nuclei', 'charge',
        'formation_J', 'formation_source_id', 'formation_source_record',
        'excitation_J', 'excitation_source_id', 'excitation_source_record',
        'physical_rest_mass_kg', 'physical_mass_source_id', 'physical_mass_source_record',
        'status', 'missing_reason'],
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strict_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate metadata key: ' + key)
            result[key] = value
        return result

    def invalid(token):
        raise ValueError('nonfinite JSON: ' + token)

    return json.loads(Path(path).read_text(), object_pairs_hook=unique,
                      parse_constant=invalid)


def extract_permitted(raw):
    """Use explicit native fields/units; uncertainty strings are never zeroed.

    The source-specific extraction follows the existing NIST preparation
    convention, but deliberately excludes its ASD threshold acquisition and
    mixed inventory. Decimal converts central neutral masses from u to kg;
    it is not an uncertainty propagation or an isotope population selection.
    """
    snapshots = {}
    for filename, expected in SOURCE_HASHES.items():
        payload = (raw / filename).read_bytes()
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError('frozen native source changed: ' + filename)
        snapshots[filename] = payload.decode('utf-8')
    native_constants = {}
    table = snapshots['codata2022.txt']
    if '2022 CODATA adjustment' not in table:
        raise ValueError('native constants version unavailable')
    for name, unit in CONSTANTS.items():
        rows = [(i + 1, line) for i, line in enumerate(table.splitlines())
                if line[:60].strip() == name]
        if len(rows) != 1:
            raise ValueError('nonunique native constant: ' + name)
        number, line = rows[0]
        native_value, uncertainty = line[60:85].strip(), line[85:110].strip()
        value = Decimal(native_value.replace(' ', ''))
        if not value.is_finite() or value <= 0 or line[110:].strip() != unit:
            raise ValueError('native constant value/unit changed: ' + name)
        native_constants[name] = dict(value=str(value), unit=unit,
            native_value=native_value, native_uncertainty=uncertainty,
            record='line:%d' % number)
    plain = unescape(re.sub(r'<[^>]+>', '', snapshots['isotopes.html']))
    isotope_records = []
    for symbol, number in ISOTOPES:
        matches = []
        for block_number, block in enumerate(plain.split('Atomic Number =')[1:], 1):
            fields = dict(re.findall(r'^\s*([^=\n]+?)\s*=\s*([^\n]*)', block, re.M))
            if fields.get('Atomic Symbol', '').strip() != symbol or fields.get('Mass Number', '').strip() != str(number):
                continue
            native = fields.get('Relative Atomic Mass', '').strip()
            found = re.fullmatch(r'(\d+(?:\.\d+)?)\((\d+)\)', native)
            if not found:
                raise ValueError('native isotope mass/uncertainty format unavailable')
            relative = Decimal(found.group(1))
            with localcontext() as context:
                context.prec = 80
                kilograms = relative * Decimal(native_constants['atomic mass constant']['value'])
            matches.append(dict(identity=str(number) + symbol, element=symbol,
                mass_number=number, relative_mass_u=str(relative), neutral_mass_kg=str(kilograms),
                native_mass_with_uncertainty=native, record='Atomic Number block:%d;Mass Number:%d' % (block_number, number)))
        if len(matches) != 1:
            raise ValueError('native isotope record missing/ambiguous: %d%s' % (number, symbol))
        isotope_records.extend(matches)
    return isotope_records, native_constants


def audit_canonical(canonical):
    """Read actual D04/D05 metadata and hashes, never repair physical gaps."""
    canonical = Path(canonical)
    required = [canonical / name for name in ('README.md', 'provenance.json', 'qa_report.json')]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError('canonical source audit prerequisite: ' + str(path))
    provenance, qa = strict_json(required[1]), strict_json(required[2])
    if qa.get('errors'):
        raise ValueError('canonical photochemistry QA contains errors')
    identities = []
    for descriptor in provenance['outputs']:
        path = canonical / descriptor['path']
        if digest(path) != descriptor['sha256']:
            raise ValueError('canonical processed identity changed: ' + str(path))
        identities.append(dict(path=str(path), sha256=descriptor['sha256']))
    def rows(name):
        with (canonical / 'processed/phidrates' / name).open(newline='') as stream:
            return list(csv.DictReader(stream))
    channels, thresholds, samples = rows('reaction_channels.csv'), rows('thresholds.csv'), rows('cross_sections.csv')
    threshold_by_id = {r['reaction_id']: Decimal(r['threshold_eV']) for r in thresholds}
    reactive = [r for r in channels if r['reaction_type'] != 'total_photoabsorption']
    long_wavelength = [dict(species=r['species'], reaction_id=r['reaction_id'],
        wavelength_nm=r['wavelength_nm'], photon_energy_eV=r['photon_energy_eV'],
        cross_section_m2=r['cross_section_m2'], source_file=r['source_file'],
        source_record_or_index=r['source_record_or_index']) for r in samples
        if r['reaction_id'] in threshold_by_id and Decimal(r['cross_section_m2']) > 0 and
        Decimal(r['photon_energy_eV']) < threshold_by_id[r['reaction_id']]]
    return dict(scope='source-audit-only-not-physical-qualification',
        metadata=[dict(path=str(p), sha256=digest(p)) for p in required],
        verified_processed_products=identities, reactive_channel_count=len(reactive),
        source_threshold_count=len(thresholds),
        reactive_channels_without_source_threshold=[r['reaction_id'] for r in reactive
                                                    if r['reaction_id'] not in threshold_by_id],
        positive_samples_below_listed_atomic_threshold=long_wavelength,
        source_reaction_labels=[{k:r[k] for k in ('reaction_id', 'species', 'products')}
                                for r in reactive],
        no_threshold_replacement=True, no_state_assignment=True,
        no_rate_or_sample_modification=True)


def prepare(raw_dir, decision_register, output_dir, canonical=None):
    raw, policy, output = Path(raw_dir).resolve(), Path(decision_register).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError('profile output must not exist')
    decisions = strict_json(policy)
    p02 = [d for d in decisions.get('decisions', []) if d.get('label') == 'P02']
    if len(p02) != 1 or p02[0].get('approval_status') != 'APPROVED':
        raise ValueError('approved P02 energy-reference definition not established')
    isotopes, constants = extract_permitted(raw)
    source_rows = [dict(id='NIST-ISOTOPES', path=str(raw / 'isotopes.html'),
        sha256=SOURCE_HASHES['isotopes.html'], role='isotope-mass'),
        dict(id='CODATA-2022', path=str(raw / 'codata2022.txt'),
        sha256=SOURCE_HASHES['codata2022.txt'], role='constants'),
        dict(id='APPROVED-P02', path=str(policy), sha256=digest(policy), role='reference-definition')]
    masses = [[r['identity'],r['neutral_mass_kg'],'NIST-ISOTOPES',r['record']] for r in isotopes]
    constant_rows = [[name,r['value'],r['unit'],'CODATA-2022',r['record']] for name,r in constants.items()]
    states = []
    zero = ['0', 'APPROVED-P02', 'P02:separated-neutral-ground-isotope-atoms-and-free-electron-reference-zero']
    for isotope in isotopes:
        states.append([isotope['identity'] + ':ground',
            isotope['identity'] + ' neutral ground reference cohort', 'ground',
            isotope['identity'] + ':1', '0'] + zero + zero +
            [isotope['neutral_mass_kg'], 'NIST-ISOTOPES', isotope['record'], 'available', ''])
    states.append(['electron:free', 'free electron reference state', 'free-electron', '', '-1'] +
        zero + zero + [constants['electron mass']['value'], 'CODATA-2022',
        constants['electron mass']['record'], 'available', ''])
    missing = []
    for isotope in isotopes:
        missing.append((isotope['identity'] + '+:unresolved', 'source-bound ion ground/residual state',
            isotope['identity'] + ':1', '1',
            'isotope/state ionization energy, uncertainty and source-threshold compatibility unresolved'))
    missing.extend([
        ('H2-16O:unresolved', 'water reference parent', '1H:2;16O:1', '0', 'isotope/rovibrational formation energy and source parent-state crosswalk missing'),
        ('16OH:unresolved', 'hydroxyl source/parent product', '16O:1;1H:1', '0', 'OH_exp compilation is not a source-state population; formation/rovibrational mapping missing'),
        ('1H2:unresolved', 'molecular hydrogen product', '1H:2', '0', 'isotope and ortho/para/rovibrational formation-state mapping missing'),
        ('16O:1D-unresolved', 'oxygen 1D term product', '16O:1', '0', 'fine-state/isotope population, complete radiative exits/quenching and threshold crosswalk unresolved'),
        ('16O:1S-unresolved', 'oxygen 1S term product', '16O:1', '0', 'isotope level uncertainty, complete radiative exits/quenching and threshold crosswalk unresolved'),
        ('16O:unlabelled-product', 'unlabelled O in water breakup', '16O:1', '0', 'native unlabelled O must not be silently assigned the ground term'),
        ('H2-16O+:unresolved', 'water ion product', '1H:2;16O:1', '1', 'molecular-ion state/formation energy and correlated electron partition missing'),
        ('16OH+:unresolved', 'hydroxyl ion product', '16O:1;1H:1', '1', 'molecular-ion state/formation energy and correlated electron partition missing'),
    ])
    for identifier, label, nuclei, charge, reason in missing:
        states.append([identifier,label,'unspecified',nuclei,charge] + [''] * 9 + ['unavailable',reason])
    audit = audit_canonical(canonical) if canonical is not None else None
    # All source extraction/semantic checks finish before creating outputs.
    output.mkdir(parents=True, exist_ok=False)
    manifest = dict(schema='srcMoon-M3-state-profile-v1', sources=source_rows,
        physical_activation=False, formal_M3_status='SKIPPED / NOT QUALIFIED',
        profile_scope='partial reference states only; no complete channel/event mapping',
        recipe=dict(path=str(Path(__file__).resolve()), sha256=digest(__file__)),
        native_isotope_metadata=isotopes, native_constant_metadata=constants,
        source_uncertainties_propagated=False,
        uncertainty_note='Native uncertainty strings retained, not converted into an assumed covariance or activation budget.',
        managed_database_modified=False, source_files_copied=False,
        available_state_count=7, unavailable_state_count=len(missing))
    for kind, rows in [('isotopes',masses),('constants',constant_rows),('states',states)]:
        path = output / (kind + '.csv')
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(HEADERS[kind]); writer.writerows(rows)
        manifest[kind] = dict(path=str(path), sha256=digest(path))
    if audit is not None:
        (output / 'canonical_gap_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
        manifest['canonical_gap_audit'] = dict(path=str(output / 'canonical_gap_audit.json'),
                                               sha256=digest(output / 'canonical_gap_audit.json'))
    (output / 'profile.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, required=True)
    parser.add_argument('--decision-register', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--canonical-root', type=Path)
    args = parser.parse_args(argv)
    try:
        value = prepare(args.raw_dir,args.decision_register,args.output_dir,args.canonical_root)
    except FileNotFoundError as error:
        print(json.dumps(dict(status='SKIPPED',reason=str(error)))); return 77
    except (OSError, ValueError, KeyError, ArithmeticError) as error:
        print(json.dumps(dict(status='ERROR',reason=str(error)))); return 2
    print(json.dumps(dict(status='PASS',scope='source-preparation-only',profile=value),indent=2)); return 0


if __name__ == '__main__':
    raise SystemExit(main())
