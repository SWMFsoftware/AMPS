"""Canonical formal-case configuration identity, never thermal qualification.

The native .in parser remains authoritative. This deliberately narrow owner
implements its documented Moon scalar/block/include grammar for preflight:
comments start at EVERY '!', continuation joins physical lines, paths belong
to the assigning file, values are not shell-expanded or escape-decoded.
External DEM/kernel/state bytes are authenticated by the existing manifest
and initial-state evidence validators, not rehashed for each of 288 cases.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


# These are the ONLY settings excluded from physical configuration identity.
# In particular, an unknown m2_verification_* key remains significant. The
# approved design independently binds each listed case/seed/control in its
# native receipts; no unknown future physical switch receives that exemption.
NUMERICAL_KEYS = frozenset({
    'random_seed', 'm2_verification_case', 'm2_verification_profile',
    'm2_verification_coarse_steps', 'm2_verification_dt_divisor',
    'm2_verification_weight_divisor', 'm2_verification_horizon_s',
    'm2_verification_coarse_timestep_s', 'm2_verification_na_weight_atoms',
    'm2_verification_ar_weight_atoms', 'm2_verification_cohort_packets',
    'm2_verification_surface_elements', 'm2_verification_source_integral_tolerance',
    'm2_verification_record_initial_thermal_state',
})
OUTPUT_KEYS = frozenset({'surface_cea_file', 'surface_tecplot_file',
                         'mesh_output_tecplot'})
# These completion adapters select approval transport or diagnostic queries;
# production initial thermal/source coefficients do not read their contents.
# Their exact files/hashes are checked by the separate I15/I19/I20 contracts.
# Do not generalize this exemption to arbitrary '*_file' future settings.
TRANSPORT_KEYS = frozenset({'m2_formal_contract_file',
                            'm2_illumination_queries_file', 'm2_polar_queries_file'})
PATH_KEYS = frozenset({
    'spice_path', 'lola_image_file', 'lola_label_file',
    'surface_cea_file', 'surface_tecplot_file',
    'diviner_temperature_table_file', 'diviner_eclipse_history_file',
    'diviner_eclipse_selection_contract_file', 'diviner_eclipse_freeze_file',
    'diviner_eclipse_authorization_file', 'argon_source_map_file',
    'surface_inventory_checkpoint_file', 'm2_formal_contract_file',
    'm2_illumination_queries_file', 'm2_polar_queries_file',
})


def _unquote(value):
    """Match MoonInput.cpp Unquote: no shell/backslash/environment expansion."""
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
        return value[1:-1]
    return value


def assignments(input_path):
    """Read one flat application:moon owner, retaining assignment locations.

    Includes may supply scalar fragments, but must balance their own blocks,
    exactly as amps_input.h requires. Unknown directives/nested owners, cycles,
    duplicates, dangling continuations and malformed scalar lines are ERROR.
    No settings/defaults are inferred from the observed model output.
    """
    result = {}; chain = []; active = False; blocks = 0

    def read(path):
        nonlocal active, blocks
        path = Path(path).resolve()
        if path in chain or len(chain) >= 64:
            raise ValueError('formal input include cycle/depth overflow: ' + str(path))
        chain.append(path); entry = active; pending = ''; first = 0; continuing = False
        try:
            logical = []
            for number, physical in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
                # Native grammar treats ! as a comment even inside quotes.
                # Do not substitute shlex's different quoted-comment rules.
                fragment = physical.split('!', 1)[0].strip()
                if not continuing:
                    first = number
                more = fragment.endswith('\\')
                if more:
                    fragment = fragment[:-1].strip()
                pending += (' ' if pending and fragment else '') + fragment
                continuing = more
                if more:
                    continue
                if pending:
                    logical.append((first, pending))
                pending = ''
            if continuing:
                raise ValueError('dangling input continuation: ' + str(path))
            for number, text in logical:
                where = str(path) + ':' + str(number)
                directive = re.fullmatch(r'#([A-Za-z]+)\s*(.*)', text)
                if directive:
                    name, body = directive[1].lower(), directive[2].strip()
                    if name == 'include':
                        if body.startswith(('"', '<')):
                            closing = '"' if body[0] == '"' else '>'
                            end = body.find(closing, 1)
                            if end < 0 or body[end+1:].strip():
                                raise ValueError(where + ': malformed #include path')
                            body = body[1:end]
                        if not body:
                            raise ValueError(where + ': missing #include path')
                        target = Path(body)
                        read(target if target.is_absolute() else path.parent/target)
                    elif name == 'begin':
                        match = re.fullmatch(r'([A-Za-z][A-Za-z0-9_]*)\s*:\s*moon', body)
                        if active or not match or match[1].lower() != 'application' or blocks:
                            raise ValueError(where + ': expected exactly one flat application: moon owner')
                        active = True; blocks += 1
                    elif name == 'end':
                        if not active or entry or body.lower() != 'application':
                            raise ValueError(where + ': mismatched application end')
                        active = False
                    else:
                        raise ValueError(where + ': unsupported native directive')
                    continue
                match = re.fullmatch(r'([A-Za-z][A-Za-z0-9_.-]*)\s*=\s*(.+)', text)
                if not active or not match:
                    raise ValueError(where + ': expected native Moon scalar assignment')
                key, value = match[1], match[2].strip()
                if key in result:
                    raise ValueError(where + ': duplicate Moon setting ' + key)
                result[key] = {'value': _unquote(value), 'file': str(path), 'line': number}
            if active != entry:
                raise ValueError(str(path) + ': application block must balance within its file')
        finally:
            chain.pop()

    read(input_path)
    if blocks != 1 or active:
        raise ValueError('formal input requires exactly one application: moon block')
    return result


def canonical_value(key, record):
    """Resolve only documented path keys; unknown values remain literal."""
    value = record['value']
    if key in PATH_KEYS and value.lower() != 'none':
        path = Path(value)
        return str((path if path.is_absolute() else Path(record['file']).parent/path).resolve())
    return value


def physics_identity(input_path):
    """Return strict physical configuration, including unknown future keys.

    The explicit native initialization default is the sole default inserted:
    omitted and explicitly declared uniform initialization denote the SAME
    production algorithm. Periodic mode and all cycle/material/thermal-step
    controls remain identity-bearing. Checkpoint path/interval are outputs
    only for zero initial inventory; checkpoint initialization retains them.
    """
    parsed = assignments(input_path)
    excluded = set(NUMERICAL_KEYS | OUTPUT_KEYS | TRANSPORT_KEYS)
    if parsed.get('surface_inventory_initial_state', {}).get('value', '').lower() == 'zero':
        excluded.update(('surface_inventory_checkpoint_file',
                         'surface_inventory_checkpoint_interval_s'))
    values = {key: canonical_value(key, record) for key, record in parsed.items()
              if key not in excluded}
    values['thermal_initialization_mode'] = values.get(
        'thermal_initialization_mode', 'declared-uniform').lower()
    return {'schema': 'srcMoon-M2-physical-configuration-v1', 'settings': values}


def physics_identity_sha256(input_path):
    """SHA-256 of deterministic UTF-8 canonical configuration JSON."""
    raw = json.dumps(physics_identity(input_path), sort_keys=True,
                     separators=(',', ':'), ensure_ascii=False).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def _native_scalar(value):
    """Quote space-bearing values without introducing unsupported escapes."""
    if not value:
        return '""'
    if '\n' in value or '\r' in value or '!' in value:
        raise ValueError('native scalar cannot contain a newline or comment delimiter')
    if (any(character.isspace() for character in value) or
            (value[0] == value[-1] and value[0] in ('"', "'"))):
        for quote in ('"', "'"):
            if quote not in value:
                return quote + value + quote
        if value != value.strip():
            raise ValueError('cannot preserve whitespace-bearing scalar with both quote delimiters')
    return value


def stage_cluster_template(template, destination):
    """Flatten/rebase formal physics without selecting a physical parameter.

    The explicit existing physical prerequisites are validated, never enabled.
    Only mesh-file ownership and checkpoint/output suppression change here.
    Numerical controls/seed are added by m2_protocol.write_case afterwards.
    """
    from integration_common import MissingPrerequisite
    template, destination = Path(template).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError('refusing to overwrite an existing formal input receipt: ' + str(destination))
    parsed = assignments(template)
    required = {'argon_source_enabled': {'true'},
                'simulation_environment_mode': {'evolving'},
                'surface_inventory_initial_state': {'zero'}}
    valid = {'argon_source_enabled': {'true','false'},
             'simulation_environment_mode': {'evolving','frozen'},
             'surface_inventory_initial_state': {'zero','checkpoint','prescribed-uniform'}}
    for key, allowed in required.items():
        value=parsed.get(key, {}).get('value', '').lower()
        if value and value not in valid[key]:
            raise ValueError('invalid native value for formal input '+key)
        if value not in allowed:
            raise MissingPrerequisite('formal I19 requires explicitly declared ' + key +
                                      '; the runner does not change physical switches')
    values = {key: canonical_value(key, record) for key, record in parsed.items()}
    values['surface_inventory_checkpoint_file'] = 'none'
    for key, name in (('surface_cea_file', 'lola_surface.cea'),
                      ('surface_tecplot_file', 'lola_surface.dat')):
        if key in values:
            values[key] = str(destination.parent/name)
    if 'mesh_output_tecplot' in values:
        values['mesh_output_tecplot'] = 'false'
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Values are native trimmed remainders, not shell/JSON strings. In
    # particular, do not introduce backslash escape processing into paths.
    text = ('! Formal I19 receipt: physical settings preserved from ' + str(template) +
            '\n#begin application: moon\n' +
            ''.join(key + ' = ' + _native_scalar(value) + '\n' for key, value in values.items()) +
            '#end application\n')
    destination.write_text(text, encoding='utf-8')
    if physics_identity(destination) != physics_identity(template):
        raise ValueError('formal I19 staging changed physical configuration identity')

