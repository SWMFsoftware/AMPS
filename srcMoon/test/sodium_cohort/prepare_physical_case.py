#!/usr/bin/env python3
"""Prepare an explicitly controlled, source-bound EXPERIMENTAL sodium case.

This is an administrative input producer, not a physical qualification or a
particle runner. In particular it never supplies an error budget, source rate,
state-energy correction, timestep, seed or particle population. The maintained
native owner validates the actual linked configuration again before execution.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from profile_projection import project, pairs, sha
from m2_initial_state import assignments

SCHEMA = 'srcMoon-native-sodium-case-controls-v1'
CROSSWALK_COLUMNS = ('channel_id', 'parent_state_id', 'ion_state_id',
                     'electron_state_id', 'energy_reference_identity',
                     'photon_energy_policy', 'source_threshold_eV',
                     'physical_threshold_source_identity')
# These names are the native sodium owner fields, rather than a second model
# vocabulary. The six source/state fields below come from the authenticated CSV;
# the remaining scientific controls must ALL appear in the caller's JSON.
TEXT_KEYS = frozenset({'case_kind', 'inertial_frame', 'inertial_origin',
    'source_reference_distance_identity', 'earth_occultation', 'model_identity',
    'solar_endpoint', 'transport_policy_identity', 'reference_weight_mode'})
WEIGHT_KEYS = frozenset({'parent_reference_weight', 'ion_reference_weight'})
UINT_KEYS = frozenset({'campaign_seed', 'genealogy_capacity', 'original_parents',
    'physical_steps', 'parent_compiled_species', 'ion_compiled_species',
    'exposure_maximum_segments', 'exposure_maximum_geometry_intervals',
    'hazard_maximum_root_refinements'})
NUMERIC_KEYS = frozenset({
    'physical_step_s', 'physical_horizon_s', 'represented_number_per_parent',
    'source_reference_distance_m', 'exposure_root_time_tolerance_s',
    'exposure_absolute_depth_tolerance', 'exposure_relative_depth_tolerance',
    'maximum_classical_speed_ms', 'unit_vector_tolerance',
    *(f'source_position_lso_m_{i}' for i in range(3)),
    *(f'source_velocity_lso_ms_{i}' for i in range(3)),
    *(f'{s}_maximum_relative_{e}_error' for s in ('parent', 'ion', 'electron')
      for e in ('inertial_mass', 'acceleration', 'relativistic_kinetic')),
    *(f'{b}_{c}_bound' for b in ('mass', 'momentum', 'energy', 'threshold',
       'rotation', 'number', 'rate', 'planner_mass', 'planner_momentum',
       'planner_energy') for c in ('absolute', 'relative'))})
FILE_KEYS = frozenset({'transport_policy_evidence_file',
                       'transport_policy_evidence_sha256'})
REQUIRED = TEXT_KEYS | UINT_KEYS | NUMERIC_KEYS | FILE_KEYS
DECIMAL_TOKEN = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')


class MissingConfiguration(ValueError):
    """Required configuration/state is unavailable: SKIPPED, never fabricated."""


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'),
        object_pairs_hook=pairs,
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON ' + x)))


def safe(value, name):
    """Reject grammar injection; AMPS treats ! as comment even inside quotes."""
    if not isinstance(value, str) or not value or any(ord(c) < 32 for c in value) \
            or any(c in value for c in '!"\\'):
        raise ValueError(name + ': nonempty native scalar without quotes/comments/control bytes required')
    return value


def absolute(path, owner=None):
    value = safe(str(path), 'path')
    p = Path(value)
    return (p if p.is_absolute() else (owner or Path.cwd()) / p).resolve()


def pinned(record, owner):
    if not isinstance(record, dict) or not isinstance(record.get('path'), str) \
            or not isinstance(record.get('sha256'), str) \
            or not re.fullmatch('[0-9a-f]{64}', record['sha256']):
        raise ValueError('source descriptor requires path and lowercase SHA-256')
    p = absolute(record['path'], owner)
    if not p.is_file() or sha(p) != record['sha256']:
        raise ValueError('source byte identity mismatch: ' + str(p))
    return {'path': str(p), 'sha256': record['sha256']}


def csv_rows(path):
    reader = csv.DictReader(io.StringIO(Path(path).read_text(encoding='utf-8'), newline=''))
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError('missing/duplicate CSV columns')
    rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('ragged CSV record')
    return tuple(reader.fieldnames), rows


def application_documents(path):
    """Retain the complete include graph, including empty scalar fragments.

    The existing strict parser first validates block/assignment semantics. This
    second traversal only inventories includes with the same comment and physical
    continuation convention; it does not interpret or replace application physics.
    """
    parsed = assignments(path)
    result = []; seen = set(); active = set()
    def visit(p):
        p = absolute(p)
        if p in active:
            raise ValueError('application include cycle')
        if p in seen:
            return
        active.add(p); seen.add(p); result.append({'path': str(p), 'sha256': sha(p)})
        pending = ''
        for line in p.read_text(encoding='utf-8').splitlines():
            text = line.split('!', 1)[0].strip()
            more = text.endswith('\\')
            pending += (' ' if pending and text else '') + (text[:-1].strip() if more else text)
            if more:
                continue
            match = re.fullmatch(r'#include\s+(.+)', pending, re.IGNORECASE)
            if match:
                child = match[1].strip()
                if child.startswith(('"', '<')):
                    child = child[1:-1]
                visit(absolute(child, p.parent))
            pending = ''
        active.remove(p)
    visit(path)
    return parsed, result


def decimal(value, key):
    # Decimal STRINGS preserve the caller's precision. JSON floats have already
    # lost lexical information, so accepting them would disguise that loss.
    if not isinstance(value, str) or not DECIMAL_TOKEN.fullmatch(value):
        raise ValueError(key + ': finite decimal string required')
    number = Decimal(value)
    if not number.is_finite():
        raise ValueError(key + ': finite decimal string required')
    return number


def validate_controls(case_path):
    case_path = absolute(case_path)
    value = load_json(case_path)
    if not isinstance(value, dict) or set(value) != {'schema', 'cohort_id', 'controls'} \
            or value['schema'] != SCHEMA or not isinstance(value['controls'], dict):
        raise ValueError('case requires exact schema/cohort_id/controls object')
    cohort = safe(value['cohort_id'], 'cohort_id')
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', cohort):
        raise ValueError('cohort_id must be one unquoted native owner name')
    controls = value['controls']
    unknown = set(controls) - REQUIRED - WEIGHT_KEYS - {'optical_clock_precision'}
    if unknown:
        raise ValueError('unknown case controls: ' + ', '.join(sorted(unknown)))
    missing = [k for k in sorted(REQUIRED) if controls.get(k) is None]
    if missing:
        raise MissingConfiguration('explicit controls unavailable: ' + ', '.join(missing))
    result = {}
    for key in TEXT_KEYS:
        result[key] = safe(controls[key], key)
        if result[key] != result[key].strip():
            raise ValueError(key + ': leading/trailing scalar whitespace is not preserved by native parsing')
    # Omitted selectors preserve the complete historical native input. An
    # explicit policy must survive writing native.in; treating a misspelling
    # as the legacy clock would compare different numerical experiments.
    if 'optical_clock_precision' in controls:
        selected=controls['optical_clock_precision']
        if selected not in ('binary64-log','compensated-extended-log'):
            raise ValueError('optical_clock_precision must be binary64-log or compensated-extended-log')
        result['optical_clock_precision']=selected
    weight_mode = result['reference_weight_mode']
    if weight_mode == 'explicit':
        if any(controls.get(k) is None for k in WEIGHT_KEYS):
            raise MissingConfiguration('explicit reference weight mode requires both parent and ion weights')
        import math
        for key in WEIGHT_KEYS:
            number = decimal(controls[key], key)
            if number <= 0 or not math.isfinite(float(number)) or float(number) == 0:
                raise ValueError(key + ': positive representable binary64 required')
            result[key] = controls[key]
    elif weight_mode == 'production':
        if set(controls) & WEIGHT_KEYS:
            raise ValueError('production reference weights forbid explicit weight assignments')
    else:
        raise ValueError('reference_weight_mode must explicitly be production or explicit')
    for key in UINT_KEYS:
        v = controls[key]
        if isinstance(v, bool) or not isinstance(v, (str, int)) \
                or not re.fullmatch('0|[1-9][0-9]*', str(v)) or int(v) > 2**64-1:
            raise ValueError(key + ': canonical uint64 required')
        result[key] = str(v)
    for key in NUMERIC_KEYS:
        v = decimal(controls[key], key)
        if key.startswith(('source_position_lso_m_', 'source_velocity_lso_ms_')):
            pass  # Signed LSO coordinates and velocities, not positive rates.
        elif v < 0:
            raise ValueError(key + ': nonnegative bound/control required')
        result[key] = controls[key]
    for key in ('physical_step_s', 'physical_horizon_s', 'represented_number_per_parent',
                'source_reference_distance_m', 'maximum_classical_speed_ms'):
        if Decimal(result[key]) <= 0:
            raise ValueError(key + ': positive required')
    for key in UINT_KEYS - {'campaign_seed', 'parent_compiled_species', 'ion_compiled_species'}:
        if int(result[key]) == 0:
            raise ValueError(key + ': positive required')
    if max(int(result[k]) for k in ('parent_compiled_species', 'ion_compiled_species')) > 2**31-1 \
            or result['parent_compiled_species'] == result['ion_compiled_species']:
        raise ValueError('distinct native int-range compiled species required')
    if Decimal(result['unit_vector_tolerance']) >= 1:
        raise ValueError('unit_vector_tolerance must be below one')
    for key in NUMERIC_KEYS:
        if '_maximum_relative_' in key and Decimal(result[key]) >= 1:
            raise ValueError(key + ': explicit relative model error limit must be below one')
    if Decimal(result['exposure_root_time_tolerance_s']) <= 0 or (
            Decimal(result['exposure_absolute_depth_tolerance']) == 0 and
            Decimal(result['exposure_relative_depth_tolerance']) == 0):
        raise ValueError('positive root-time tolerance and at least one positive exposure bound required')
    if result['inertial_origin'] != 'SUN':
        raise ValueError('native physical cohort requires explicitly SUN-origin inertial coordinates')
    if result['case_kind'] not in ('homogeneous-fixed-position', 'trajectory-integrated-hazard') \
            or result['solar_endpoint'] not in ('quiet', 'active') \
            or result['earth_occultation'] not in ('on', 'off'):
        raise ValueError('unsupported explicit case/solar/occultation selection')
    if result['case_kind'] == 'homogeneous-fixed-position' and any(
            Decimal(result[f'source_velocity_lso_ms_{i}']) != 0 for i in range(3)):
        raise ValueError('homogeneous case requires explicitly zero initial velocity')
    # The native owner converts these values to binary64. Exact dyadic equality
    # is a conservative admission check: never truncate a final physical step.
    doubles = {'physical_step_s', 'physical_horizon_s', 'represented_number_per_parent',
               'source_reference_distance_m', 'exposure_root_time_tolerance_s',
               'number_absolute_bound', 'number_relative_bound',
               *(f'source_position_lso_m_{i}' for i in range(3)),
               *(f'source_velocity_lso_ms_{i}' for i in range(3)),
               *(f'planner_{b}_{c}_bound' for b in ('mass','momentum','energy')
                 for c in ('absolute','relative'))}
    import math
    for key in doubles:
        v = float(result[key])
        if not math.isfinite(v) or (Decimal(result[key]) != 0 and v == 0):
            raise ValueError(key + ': not representable as non-underflowing binary64')
    if Fraction.from_float(float(result['physical_step_s'])) * int(result['physical_steps']) \
            != Fraction.from_float(float(result['physical_horizon_s'])):
        raise ValueError('steps times binary64 step must equal exact binary64 horizon')
    evidence = pinned({'path': controls['transport_policy_evidence_file'],
                       'sha256': controls['transport_policy_evidence_sha256']}, case_path.parent)
    result['transport_policy_evidence_file'] = evidence['path']
    result['transport_policy_evidence_sha256'] = evidence['sha256']
    return cohort, result, evidence


def authenticate_crosswalk(profile, profile_path, crosswalk_path):
    extension = profile.get('atomic_extension')
    if not isinstance(extension, dict):
        raise MissingConfiguration('source-bound atomic crosswalk metadata unavailable')
    descriptor = pinned(extension['crosswalk_csv'], profile_path.parent)
    if absolute(crosswalk_path) != Path(descriptor['path']):
        raise ValueError('crosswalk must be the exact profile-bound CSV path')
    meta_pin = pinned(extension['crosswalk'], profile_path.parent)
    meta = load_json(meta_pin['path'])
    if not isinstance(meta, dict) or meta.get('schema') != 'srcMoon-M3-Na-source-energy-crosswalk-v1':
        raise ValueError('unsupported atomic crosswalk metadata schema')
    columns, rows = csv_rows(crosswalk_path)
    if columns != CROSSWALK_COLUMNS or len(rows) != 1:
        raise ValueError('one exact eight-column sodium crosswalk required')
    row = rows[0]
    for key in CROSSWALK_COLUMNS:
        safe(row[key], key)
        if row[key] != meta.get(key):
            raise ValueError('CSV and source-bound crosswalk metadata differ: ' + key)
    if row['photon_energy_policy'] != 'source-effective-photon':
        raise ValueError('eight-column owner requires explicit source-effective-photon')
    _, states = csv_rows(pinned(profile['states'], profile_path.parent)['path'])
    by_id = {}
    for state in states:
        if state['state_id'] in by_id:
            raise ValueError('duplicate profile state ID')
        by_id[state['state_id']] = state
    sources = {x['id']: x for x in profile['sources']}
    def identity(state, component):
        source = sources[state[component + '_source_id']]
        return source['id'] + '@sha256:' + source['sha256'] + '#' + state[component + '_source_record']
    selected = []
    for key in ('parent_state_id', 'ion_state_id', 'electron_state_id'):
        state = by_id.get(row[key])
        if state is None or state['status'] != 'available':
            raise MissingConfiguration('requested source state unavailable: ' + row[key])
        selected.append(state)
    parent, ion, electron = selected
    if [(s['nuclei'], s['charge']) for s in selected] != [('23Na:1','0'), ('23Na:1','1'), ('','-1')]:
        raise ValueError('explicit Na ground/ion/electron state topology differs')
    if identity(ion, 'formation') != row['physical_threshold_source_identity']:
        raise ValueError('physical threshold is not the exact ion formation source identity')
    for state in (parent, electron):
        if Decimal(state['formation_J']) != 0 or identity(state, 'formation') != row['energy_reference_identity']:
            raise ValueError('parent/electron must bind the declared common ground reference')
    for state in selected:
        if Decimal(state['excitation_J']) != 0 or identity(state, 'excitation') != row['energy_reference_identity']:
            raise ValueError('non-ground or differently referenced excitation is unsupported')
    files = [pinned(record, Path(meta_pin['path']).parent) for record in meta.get('source_files', [])]
    if not files:
        raise MissingConfiguration('exact native threshold/rate evidence unavailable')
    matches = []
    for item in files:
        cols, data = csv_rows(item['path'])
        if 'threshold_eV' in cols and 'reaction_id' in cols:
            matches += [r for r in data if r['reaction_id'] == row['channel_id']]
    if len(matches) != 1 or float(decimal(matches[0]['threshold_eV'], 'native threshold')) \
            != float(decimal(row['source_threshold_eV'], 'crosswalk threshold')):
        raise ValueError('crosswalk threshold differs from exact native binary64 source threshold')
    return row, [descriptor, meta_pin, *files]


def template():
    """An explicitly INCOMPLETE template; null values cannot produce native.in."""
    return {'schema': SCHEMA, 'cohort_id': 'REPLACE_WITH_EXPLICIT_COHORT_ID',
            'controls': {key: None for key in sorted(REQUIRED)}}


def prepare(profile, crosswalk, application_input, case_json, output_dir):
    if Path(output_dir).is_symlink():
        raise ValueError('output path is a symlink; preserve existing path')
    output = absolute(output_dir)
    if output.exists() or output.is_symlink():
        raise ValueError('output exists; preserve earlier artifacts: ' + str(output))
    if not output.parent.is_dir():
        raise ValueError('output parent must already exist')
    profile_path = absolute(profile); crosswalk_path = absolute(crosswalk)
    application_path = absolute(application_input); case_path = absolute(case_json)
    cohort, controls, evidence = validate_controls(case_path)
    parsed, documents = application_documents(application_path)
    expected = {'tracing_clock_mode': {'common-time'}, 'simulation_environment_mode': {'frozen'},
        'particle_tracing_mover': {'force-cut-surface-first-order', 'force-cut-surface-second-order'}}
    for key, allowed in expected.items():
        if key not in parsed:
            raise MissingConfiguration('application must explicitly declare ' + key)
        if parsed[key]['value'] not in allowed:
            raise ValueError('unsupported native cohort application setting: ' + key)
    # Never infer an omitted M2 selection; this is an explicit asynchronous model
    # compatibility guard, not a substitute for full native application parsing.
    if parsed.get('m2_environment', {}).get('value') not in ('off', 'disabled'):
        raise MissingConfiguration('application must explicitly declare m2_environment = off')
    with tempfile.TemporaryDirectory(prefix='moon-physical-case-') as temporary:
        stage = Path(temporary)
        projection = project(profile_path, stage / 'profile.native.in')
        value = load_json(profile_path)
        if sha(profile_path) != projection['profile_json_sha256']:
            raise ValueError('profile changed while projecting source identities')
        row, crosswalk_pins = authenticate_crosswalk(value, profile_path, crosswalk_path)
        all_pins = [{'path': str(p), 'sha256': sha(p)} for p in (profile_path, case_path)]
        all_pins += documents + [evidence] + crosswalk_pins
        source_lines = []
        native = {'schema': 'srcMoon-native-sodium-cohort-v1', 'application_input': str(application_path),
            'profile_metadata_file': str(profile_path), 'profile_metadata_sha256': sha(profile_path),
            'profile_native_manifest_file': str(output / 'profile.native.in'),
            'profile_native_manifest_sha256': projection['sha256'],
            'crosswalk_file': str(crosswalk_path), 'crosswalk_sha256': sha(crosswalk_path),
            **{key: row[key] for key in CROSSWALK_COLUMNS[:6]}, **controls}
        for key in ('isotopes', 'constants', 'states'):
            record = pinned(value[key], profile_path.parent); all_pins.append(record)
            native[key + '_file'] = record['path']; native[key + '_sha256'] = record['sha256']
        for item in value['sources']:
            record = pinned(item, profile_path.parent); all_pins.append(record)
            role = safe(item['role'], 'source role')
            if role != role.strip():
                raise ValueError('source role whitespace is not preserved by native parsing')
            source_lines += ['  #begin source: ' + item['id'],
                '    path = "' + record['path'] + '"', '    sha256 = ' + record['sha256'],
                '    role = ' + role, '  #end source']
        lines = ['! PREPARATION ONLY: EXPERIMENTAL; no qualification or budget approval implied.',
                 '#begin sodium_qualification: ' + cohort]
        for key in sorted(native):
            val = safe(native[key], key)
            lines.append('  ' + key + ' = ' + ('"' + val + '"' if key.endswith('_file') or key == 'application_input' else val))
        lines += source_lines + ['#end sodium_qualification', '']
        (stage / 'native.in').write_text('\n'.join(lines), encoding='utf-8')
        case = {'schema': SCHEMA, 'cohort_id': cohort, 'controls': controls,
                'source_crosswalk': row, 'scope': 'explicit input preparation only; NOT QUALIFIED'}
        (stage / 'case.json').write_text(json.dumps(case, indent=2) + '\n')
        unique = {}
        for record in all_pins:
            if record['path'] in unique and unique[record['path']] != record['sha256']:
                raise ValueError('contradictory input identity')
            unique[record['path']] = record['sha256']
        # Catch concurrent mutation before publication; no original is modified.
        for path, digest in unique.items():
            if sha(path) != digest:
                raise ValueError('source changed during preparation: ' + path)
        import subprocess
        version = subprocess.run(['git', '-C', str(HERE.parents[2]), 'rev-parse', 'HEAD'],
                                  text=True, capture_output=True)
        if version.returncode:
            raise ValueError('actual source revision unavailable: ' + version.stderr.strip())
        revision = version.stdout.strip()
        manifest = {'schema': 'srcMoon-sodium-case-preparation-manifest-v1',
            'source_revision': revision, 'physical_activation': False,
            'formal_M3_status': 'SKIPPED / NOT QUALIFIED', 'scope': 'input/source preflight only',
            'inputs': [{'path': p, 'sha256': h} for p, h in unique.items()],
            'application_documents': documents,
            'preparer_dependencies': [{'path': str(p), 'sha256': sha(p)} for p in
                (Path(__file__).resolve(), HERE / 'profile_projection.py',
                 HERE.parent / 'm2_initial_state.py', HERE / 'native.cpp',
                 HERE.parents[2] / 'src/general/input/amps_input.h')],
            'outputs': [{'path': str(output / name), 'sha256': sha(stage / name)}
                        for name in ('native.in', 'profile.native.in', 'case.json')]}
        (stage / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        result = {'schema': 'srcMoon-sodium-case-preparation-result-v1', 'status': 'PASS',
            'scope': 'explicit source-bound preparation only; no particles executed',
            'formal_M3_status': 'SKIPPED / NOT QUALIFIED', 'output_dir': str(output),
            'native_input': str(output / 'native.in'), 'passed': 1, 'failed': 0, 'errors': 0, 'skipped': 0}
        (stage / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        output.mkdir()  # Atomic exclusive directory creation; never replace an earlier case.
        try:
            for p in stage.iterdir():
                shutil.copyfile(p, output / p.name)
        except BaseException:
            shutil.rmtree(output)
            raise
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-template', metavar='NEW_JSON')
    for flag in ('profile', 'crosswalk', 'application-input', 'case-json', 'output-dir'):
        parser.add_argument('--' + flag)
    args = parser.parse_args(argv)
    try:
        if args.write_template:
            if any(getattr(args, k) for k in ('profile','crosswalk','application_input','case_json','output_dir')):
                raise ValueError('template mode cannot be combined with preparation')
            with absolute(args.write_template).open('x', encoding='utf-8') as stream:
                stream.write(json.dumps(template(), indent=2) + '\n')
            print(json.dumps({'status': 'SKIPPED', 'scope': 'INCOMPLETE template, NOT RUNNABLE'}))
            return 77
        if any(getattr(args, k) is None for k in ('profile','crosswalk','application_input','case_json','output_dir')):
            raise MissingConfiguration('all five explicit preparation arguments required')
        result = prepare(args.profile, args.crosswalk, args.application_input, args.case_json, args.output_dir)
        print(json.dumps(result, indent=2)); return 0
    except MissingConfiguration as error:
        print(json.dumps({'status': 'SKIPPED', 'reason': str(error), 'scope': 'no runnable case published'})); return 77
    except (ValueError, OSError, KeyError, TypeError, InvalidOperation, OverflowError) as error:
        print(json.dumps({'status': 'ERROR', 'reason': str(error), 'scope': 'no runnable case published'})); return 2


if __name__ == '__main__':
    raise SystemExit(main())
