"""Lossless new molecular state/transition source views, never a Q provider.

The ExoMol mandatory four state fields and four transition fields have a
published schema. Additional state tokens are deliberately left positional:
the acquired neutral OH definition disagrees with its actual quantum-label
order. Preserving those tokens is useful source access; guessing their state
assignment would create unsupported science. Declared count discrepancies are
reported, with no dropped/created record and no complete-cascade claim.
"""
from __future__ import annotations
import argparse
import bz2
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strict_json(path):
    def pairs(rows):
        value = {}
        for key, item in rows:
            if key in value: raise ValueError('duplicate definition JSON key: ' + key)
            value[key] = item
        return value
    def bad(value): raise ValueError('nonfinite definition JSON: ' + value)
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs, parse_constant=bad)


def number(token):
    try: value = Decimal(token.replace('D', 'E').replace('d', 'e'))
    except InvalidOperation as error: raise ValueError('invalid native numerical token: ' + token) from error
    if not value.is_finite(): raise ValueError('nonfinite native numerical token: ' + token)
    return value


def integer(token):
    if not re.fullmatch('[0-9]+', token) or int(token) <= 0:
        raise ValueError('positive native integer required: ' + token)
    return int(token)


def tokens(path):
    raw = Path(path).read_bytes()
    if Path(path).suffix == '.bz2': raw = bz2.decompress(raw)
    lines = raw.decode('ascii').splitlines()
    if not lines or any(not line.strip() for line in lines):
        raise ValueError('empty/blank native row is not silently omitted')
    return [line.split() for line in lines]


def inspect(definition, definition_text, states, transitions):
    """Inspect actual rows and declared metadata without publishing anything.

    A completeness energy in the text definition is NOT a maximum allowed
    state energy. JSON calls the same quantity max_energy; both labels and the
    exact source line are retained so callers cannot infer a false defect.
    """
    value = strict_json(definition)
    metadata = value['dataset']; fields = metadata['states']['states_file_fields']
    names = [field['name'] for field in fields]
    if len(names) < 4 or names[0] not in ('ID', 'StateID') or names[1] != 'E' \
            or names[2] not in ('gtot', 'g_tot') or names[3] != 'J':
        raise ValueError('source definition does not declare the supported mandatory state fields')
    state_rows = tokens(states); transition_rows = tokens(transitions)
    widths = {len(row) for row in state_rows}
    if len(widths) != 1 or next(iter(widths)) < 4:
        raise ValueError('ragged native state fields')
    width = next(iter(widths)); by_id = {}; energies = []
    for row in state_rows:
        identifier = integer(row[0])
        if identifier in by_id: raise ValueError('duplicate numerical native state ID')
        integer(row[2]); j = number(row[3])
        if j < 0 or 2*j != (2*j).to_integral_value():
            raise ValueError('native J is not a nonnegative integer/half integer')
        energy = number(row[1]); energies.append(energy); by_id[identifier] = row
    for row in transition_rows:
        if len(row) != 4: raise ValueError('this source view requires native four-column transitions')
        upper, lower = integer(row[0]), integer(row[1])
        if upper not in by_id or lower not in by_id:
            raise ValueError('transition points to an unavailable native state')
        if number(row[2]) < 0 or number(row[3]) < 0:
            raise ValueError('negative native Einstein coefficient/transition wavenumber')
    completeness = []; versions = []
    for line_number, line in enumerate(Path(definition_text).read_text().splitlines(), 1):
        if '#' not in line: continue
        literal, label = (x.strip() for x in line.split('#', 1))
        if label.lower().startswith('version number'):
            versions.append(literal)
        if 'higher energy with complete set of transitions' in label.lower():
            # NA is an actual source gap (not zero or an infinite bound).
            # Lexical access remains possible; a physical completeness-based
            # provider would have to refuse this unavailable declaration.
            known = literal.upper() not in ('NA', 'NAN', 'NONE')
            if known: number(literal)
            completeness.append({'line': line_number, 'native_value': literal,
                                 'native_label': label, 'numeric_value_available': known})
    if len(completeness) != 1:
        raise ValueError('exact native transition-completeness energy declaration required')
    if versions != [str(metadata['version'])]:
        raise ValueError('native text and JSON definition versions differ or are unavailable')
    limitations = [
        'First four mandatory state fields alone are semantically labelled; additional native tokens remain positional.',
        'Native text completeness energy is not a maximum source-state energy.',
        'No absolute formation/ionization/dissociation energy or product-state crosswalk is supplied.',
        'No A-coefficient uncertainty, omitted-exit bound, full-lifetime or cascade completeness is invented.',
        'No interpolation, normalization, state creation, missing-record repair or physical activation occurs.']
    count_matches = len(state_rows) == metadata['states']['number_of_states']
    transitions_match = len(transition_rows) == metadata['transitions']['number_of_transitions']
    if not count_matches: limitations.append('Declared native state count differs from actual retained rows; both remain explicit.')
    if not transitions_match: limitations.append('Declared native transition count differs from actual retained rows; both remain explicit.')
    if len(names) != width: limitations.append('Declared state-column count differs from actual width; optional field mapping is unresolved.')
    columns = names[:4] + ['native_token_' + str(i+1) for i in range(4, width)]
    report = {'schema': 'amps-new-molecular-exomol-lexical-qa-v1',
        'status': 'PASS', 'scope': 'source byte/token/invariant QA only; not physical readiness',
        'scientific_status': 'PARTIAL_SOURCE_INPUT_NOT_QUALIFIED',
        'isotopologue': value['isotopologue'], 'dataset_name': metadata['name'],
        'dataset_version': metadata['version'], 'dataset_doi': metadata.get('doi'),
        'state_columns': columns, 'declared_state_columns': names,
        'declared_state_count': metadata['states']['number_of_states'], 'actual_state_count': len(state_rows),
        'declared_transition_count': metadata['transitions']['number_of_transitions'],
        'actual_transition_count': len(transition_rows), 'declared_state_count_matches': count_matches,
        'declared_transition_count_matches': transitions_match,
        'actual_maximum_state_energy_cm-1': str(max(energies)),
        'json_field_dataset_states_max_energy': metadata['states'].get('max_energy'),
        'native_transition_completeness_energy': completeness[0],
        'zero_energy_native_rows': [row for row, e in zip(state_rows, energies) if e == 0],
        'missing_transition_state_references': 0, 'limitations': limitations}
    return report, columns, state_rows, transition_rows


def project(definition, definition_text, states, transitions, output, *, expected_hashes):
    """Authenticate all four original inputs and publish new lexical CSVs.

    Every expected SHA is explicit, normally from the acquisition receipt.
    Output remains absent for invalid/truncated/changed inputs. Native tokens
    and row order pass to CSV unchanged; original bzip2 bytes are not rewritten.
    """
    paths = dict(zip(('definition', 'definition_text', 'states', 'transitions'),
                     map(lambda p: Path(p).resolve(), (definition, definition_text, states, transitions))))
    output = Path(output)
    if output.exists() or output.is_symlink(): raise ValueError('preserve existing source view')
    if set(expected_hashes) != set(paths): raise ValueError('all four explicit original SHA identities required')
    for key, path in paths.items():
        if not re.fullmatch('[0-9a-f]{64}', expected_hashes[key]) or sha(path) != expected_hashes[key]:
            raise ValueError('original molecular source identity mismatch: ' + key)
    report, columns, state_rows, transition_rows = inspect(**paths)
    for key, path in paths.items():
        if sha(path) != expected_hashes[key]: raise ValueError('source changed while preparing: ' + key)
    output.mkdir(parents=True, exist_ok=False)
    for name, header, rows in [('states.csv', columns, state_rows),
        ('transitions.csv', ['upper_state_id', 'lower_state_id', 'Einstein_A_s-1', 'wavenumber_cm-1'], transition_rows)]:
        with (output/name).open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(header); writer.writerows(rows)
    report['inputs'] = {key: {'path': str(path), 'sha256': expected_hashes[key]} for key, path in paths.items()}
    report['outputs'] = {name: {'path': str((output/name).resolve()), 'sha256': sha(output/name)}
                         for name in ('states.csv', 'transitions.csv')}
    (output/'qa.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('definition', 'definition-text', 'states', 'transitions', 'hashes-json', 'output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(project(args.definition, args.definition_text, args.states, args.transitions,
                            args.output_dir, expected_hashes=strict_json(args.hashes_json)), indent=2))
