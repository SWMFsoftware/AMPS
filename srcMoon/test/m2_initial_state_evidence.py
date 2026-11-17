"""Bind the APPLIED native t0 thermal state to an explicitly reviewed proof.

This module contains no thermal evolution, observation or residence law. The
full temperature grid is a model initial condition, never independent truth.
The production initializer must actually apply the reviewed method before
source/weight queries. A generic epoch-level READY claim cannot establish it.
"""
from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import sys

SCHEMA = 'srcMoon-M2-initial-thermal-state-qualification-v1'
STATE_HEADER = ['logical_surface_element', 'depth_index', 'depth_m', 'kinetic_temperature_K']
STATE_SUFFIX = '.I19_initial_thermal_state.csv'
HERE = Path(__file__).resolve().parent


class Unavailable(RuntimeError):
    """Exact authoritative qualification or native capability is unavailable."""


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def document(path):
    def pairs(values):
        result = {}
        for key, value in values:
            if key in result:
                raise ValueError('duplicate initial-state evidence key: '+key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('nonfinite initial-state JSON constant: '+value)
    value = json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=invalid)
    if not isinstance(value, dict):
        raise ValueError('initial-state evidence must be a JSON object')
    return value


def _number(value, label, *, minimum=None):
    if isinstance(value, bool):
        raise ValueError(label+' must not be Boolean')
    try:
        number = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(label+' must be finite numeric') from error
    if not math.isfinite(number) or (minimum is not None and number < minimum):
        raise ValueError(label+' is nonfinite or outside its declared domain')
    return number


def _integer(value, label, *, minimum=0):
    number = _number(value, label, minimum=minimum)
    if number != int(number):
        raise ValueError(label+' must be an integer')
    return int(number)


def _stamp(value):
    try:
        stamp = datetime.fromisoformat(str(value).replace('Z', '+00:00'))
    except ValueError as error:
        raise ValueError('invalid initial-state review/campaign timestamp') from error
    if stamp.tzinfo is None:
        raise ValueError('initial-state chronology must explicitly declare a timezone')
    return stamp


def _file(record, base, label):
    if record is None:
        raise Unavailable(label+' path/SHA-256 evidence unavailable')
    if not isinstance(record, dict) or set(record) != {'path', 'sha256'}:
        raise ValueError(label+' requires exactly path and sha256')
    if not isinstance(record['path'], str) or not record['path'] or not re.fullmatch('[0-9a-f]{64}', str(record['sha256'])):
        raise ValueError(label+' has malformed path/SHA-256')
    path = Path(record['path'])
    path = (path if path.is_absolute() else Path(base)/path).resolve()
    if not path.is_file():
        raise Unavailable(label+' file unavailable: '+str(path))
    if sha256(path) != record['sha256']:
        raise ValueError(label+' SHA-256 mismatch: '+str(path))
    return path


def read_state(path, elements, layers):
    """Validate the COMPLETE native boundary/layer grid in SI, no averaging.

    Native output is ordered by element and then boundary(-1)/layer centres.
    Requiring that exact complete sequence rejects duplicate/missing cells
    without keeping hundreds of thousands of dictionaries in RAM. One model
    depth grid is shared by the configured production columns.
    """
    elements = _integer(elements, 'logical surface element count', minimum=1)
    layers = _integer(layers, 'thermal layer count', minimum=2)
    count = 0;grid = None;current_grid = [];minimum = math.inf;maximum = 0.0
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != STATE_HEADER:
            raise ValueError('native initial thermal-state header differs from the declared full-grid schema')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('malformed native initial-state CSV row')
            expected_element, offset = divmod(count, layers+1)
            expected_layer = offset-1
            element = _integer(row['logical_surface_element'], 'native state element')
            layer = _integer(row['depth_index'], 'native state layer index', minimum=-1)
            depth = _number(row['depth_m'], 'native layer depth [m]', minimum=0)
            temperature = _number(row['kinetic_temperature_K'], 'native kinetic temperature [K]', minimum=0)
            if expected_element >= elements or element != expected_element or layer != expected_layer:
                raise ValueError('native initial-state element/layer grid is incomplete, duplicated or out of order')
            if temperature == 0 or (layer == -1 and depth != 0) or (layer >= 0 and depth <= 0):
                raise ValueError('native initial-state boundary/depth/positive-temperature identity failed')
            if offset == 0:
                current_grid = []
            if current_grid and depth <= current_grid[-1]:
                raise ValueError('native thermal layer centres must increase strictly in metres')
            current_grid.append(depth)
            if offset == layers:
                if grid is not None and current_grid != grid:
                    raise ValueError('native initial-state layer grid changes between logical elements')
                grid = current_grid[:]
            minimum = min(minimum, temperature);maximum = max(maximum, temperature);count += 1
    if count != elements*(layers+1):
        raise ValueError('native initial thermal-state file omits required elements/layers')
    return {'rows': count, 'logical_surface_element_count': elements, 'layer_count': layers,
            'boundary_and_layer_depth_m': grid, 'minimum_kinetic_temperature_K': minimum,
            'maximum_kinetic_temperature_K': maximum}


def _source_signature():
    from build_m2_qualification import production_source_identity
    return hashlib.sha256(json.dumps(production_source_identity(), sort_keys=True,
                         separators=(',', ':')).encode()).hexdigest()


def physical_external_data(manifest,application,receipt,*,expected=None):
    """Bind the ACTUAL CSPICE pool and configured DEM, ignoring observation extras.

    The ordered furnished-kernel list comes from native ktotal/kdata, not a
    guessed installation directory. Order/duplicates are retained because
    kernel load priority can matter. The manifest may use its established
    relative aliases; each must describe the actual named file's bytes. The
    returned map records the authoritative native absolute filenames/hashes.
    """
    from m2_initial_state import assignments,canonical_value
    declared=manifest.get('external_data')
    if not isinstance(declared,dict) or not declared:
        raise Unavailable('native physical SPICE/DEM external-data identities unavailable')
    if 'thermal_initialization_spice_kernel_count' not in receipt:
        raise Unavailable('native actual furnished-SPICE provenance capability unavailable')
    count=_integer(receipt['thermal_initialization_spice_kernel_count'],'actual furnished kernel count',minimum=1)
    keys={key for key in receipt if re.fullmatch(r'thermal_initialization_spice_kernel_\d+_file',key)}
    if len(keys)!=count or keys!={'thermal_initialization_spice_kernel_'+str(index)+'_file' for index in range(count)}:
        raise ValueError('actual furnished kernel receipt count/list is incomplete or contradictory')
    def checked(path,*,allow_relative_alias):
        if not path.is_absolute():raise ValueError('native physical file provenance must use absolute paths')
        path=path.resolve()
        if not path.is_file():raise Unavailable('actual physical initial-state file unavailable: '+str(path))
        actual=sha256(path)
        exact=[digest for name,digest in declared.items() if Path(name).is_absolute() and Path(name).resolve()==path]
        aliases=exact or ([digest for name,digest in declared.items() if
            not Path(name).is_absolute() and Path(name).name==path.name] if allow_relative_alias else [])
        if not aliases:raise Unavailable('native manifest omits actual physical file hash: '+str(path))
        if any(not isinstance(digest,str) or not re.fullmatch('[0-9a-f]{64}',digest) for digest in aliases):
            raise ValueError('malformed physical external SHA-256')
        if set(aliases)!={actual}:
            raise ValueError('actual furnished SPICE/DEM bytes contradict native manifest: '+str(path))
        return {'path':str(path),'sha256':actual}
    kernels=[checked(Path(str(receipt['thermal_initialization_spice_kernel_'+str(index)+'_file'])),
                     allow_relative_alias=True) for index in range(count)]
    parsed=assignments(application);terrain=[]
    if 'surface_geometry' not in parsed:raise Unavailable('explicit physical surface geometry unavailable')
    if canonical_value('surface_geometry',parsed['surface_geometry']).lower()=='lola':
        for key in ('lola_image_file','lola_label_file'):
            if key not in parsed:raise Unavailable('configured LOLA dependency missing from physical input')
            terrain.append(checked(Path(canonical_value(key,parsed[key])),allow_relative_alias=False))
    result={'ordered_actual_spice_kernels':kernels,'configured_lola_files':terrain}
    if expected is not None and result!=expected:
        raise ValueError('actual native case physical SPICE/DEM identities differ from qualified pilot')
    return result


def _receipt_state(receipt, application, state, shape, expected_mode, expected_epoch):
    """Check real native initialization diagnostics against explicit input."""
    from m2_illumination import setting
    required = ('thermal_initialization_mode', 'thermal_initialization_et_s',
                'initial_thermal_state_csv', 'initial_thermal_state_rows',
                'initial_thermal_state_element_count', 'initial_thermal_state_layer_count')
    if any(key not in receipt for key in required):
        raise Unavailable('native applied initial-state exporter/receipt capability unavailable')
    if receipt['thermal_initialization_mode'] != expected_mode or (
            _number(receipt['thermal_initialization_et_s'], 'applied initial-state TDB epoch') != expected_epoch):
        raise ValueError('native applied initialization mode/epoch differs from qualified state')
    if _number(receipt['initial_epoch_tdb_s'], 'native numerical initial TDB epoch') != expected_epoch:
        raise ValueError('numerical clock does not start at its applied initial thermal epoch')
    # The normal runner passes an absolute output prefix. Bind the path in
    # that actual receipt, rather than accepting an unrelated convenient CSV.
    recorded = Path(str(receipt['initial_thermal_state_csv']))
    recorded = (recorded if recorded.is_absolute() else Path(state).parent/recorded).resolve()
    if recorded != Path(state).resolve():
        raise ValueError('native receipt points to a different initial thermal-state CSV')
    for key, expected in (('initial_thermal_state_rows', shape['rows']),
                          ('initial_thermal_state_element_count', shape['logical_surface_element_count']),
                          ('initial_thermal_state_layer_count', shape['layer_count'])):
        if _integer(receipt[key], 'native '+key) != expected:
            raise ValueError('native initial-state count differs from complete qualified grid: '+key)
    mode = setting(application, 'thermal_initialization_mode', 'declared-uniform')
    if mode != expected_mode:
        raise ValueError('actual input selects a different thermal initialization method')
    if setting(application, 'surface_temperature_mode') != 'thermal-model':
        raise ValueError('applied depth-resolved initial thermal qualification requires thermal-model')
    if expected_mode == 'periodic-spinup':
        if _integer(receipt.get('thermal_initialization_spinup_completed_cycles'), 'native spin-up cycles', minimum=1) < 1:
            raise ValueError('periodic initialization did not complete a production spin-up cycle')
        end = _number(receipt.get('thermal_initialization_spinup_cycle_end_et_s'), 'native spin-up end TDB')
        begin = _number(receipt.get('thermal_initialization_spinup_cycle_start_et_s'), 'native spin-up begin TDB')
        if end != expected_epoch or begin >= end:
            raise ValueError('applied periodic forcing interval does not end at the initial simulation epoch')
        for key, control in (
                ('thermal_initialization_spinup_maximum_node_temperature_change_k', 'thermal_spinup_temperature_tolerance_k'),
                ('thermal_initialization_spinup_maximum_column_enthalpy_change_j_m2', 'thermal_spinup_energy_tolerance_j_m2')):
            if _number(receipt.get(key), 'native '+key, minimum=0) > _number(setting(application, control), control, minimum=0):
                raise ValueError('claimed qualified periodic state contradicts its unchanged production convergence criterion')
    elif expected_mode == 'declared-uniform':
        uniform = _number(setting(application, 'thermal_initial_temperature_k'), 'declared uniform initialization K', minimum=0)
        if uniform == 0 or shape['minimum_kinetic_temperature_K'] != uniform or shape['maximum_kinetic_temperature_K'] != uniform:
            raise ValueError('native state is not the explicitly selected uniform kinetic temperature')
    else:
        raise ValueError('unsupported explicit thermal initialization mode')


def load_qualification(support, application, executable, *, base,allow_surface_pilot=False):
    """Authenticate exact reviewed pilot state WITHOUT any numerical residual.

    A state may be reviewed after its real pilot ran, but that review must
    precede every numerical campaign case. No qualifying artifact is created
    here. Missing old/generic proof is unavailable; contradictory new proof
    is corruption. Review/source/method records remain distinct from science.
    """
    q = support.get('actual_epoch_thermal_qualification')
    if not isinstance(q, dict) or q.get('schema') != SCHEMA:
        raise Unavailable('exact applied native initial-state qualification unavailable; generic epoch proof is insufficient')
    review = q.get('review', {})
    if q.get('status') != 'READY' or not isinstance(review, dict) or review.get('approved') is not True or (
            not review.get('reviewer') or review.get('frozen_before_numerical_campaign') is not True):
        raise Unavailable('applied initial-state qualification is not explicitly reviewed/frozen')
    approved = _stamp(review.get('approved_utc'))
    from m2_initial_state import physics_identity_sha256
    from m2_illumination import setting
    from integration_common import parse_receipt
    application = Path(application).resolve();executable = Path(executable).resolve()
    if not application.is_file() or not executable.is_file():
        raise Unavailable('actual physical template or retained executable for initial-state proof unavailable')
    mode = q.get('initialization_mode')
    if mode not in ('declared-uniform', 'periodic-spinup'):
        raise ValueError('initial-state qualification names an unsupported initialization method')
    if q.get('simulation_start_utc') != setting(application, 'simulation_start_utc'):
        raise ValueError('qualified initial-state UTC differs from actual physical input')
    epoch = _number(q.get('ephemeris_time_tdb_s'), 'qualified initial TDB epoch')
    if q.get('physics_identity_sha256') != physics_identity_sha256(application):
        raise ValueError('qualified initial-state physical input identity differs from actual template')
    if q.get('executable_sha256') != sha256(executable) or q.get('production_source_signature_sha256') != _source_signature():
        raise ValueError('qualified initial-state executable/native-source identity differs')
    paths = {key: _file(q.get(key), base, 'qualified '+key)
             for key in ('pilot_input', 'pilot_manifest', 'pilot_receipt', 'native_state')}
    if physics_identity_sha256(paths['pilot_input']) != q['physics_identity_sha256']:
        raise ValueError('actual pilot applied a different physical input than the qualified campaign')
    native = parse_receipt(paths['pilot_receipt']);manifest = document(paths['pilot_manifest'])
    allowed={'I19','I20'} if allow_surface_pilot else {'I19'}
    pilot_test=manifest.get('test_id')
    if (pilot_test not in allowed or native.get('test_id') != pilot_test or
            manifest.get('executable', {}).get('sha256') != q['executable_sha256'] or
            manifest.get('inputs', {}).get(str(paths['pilot_input'])) != sha256(paths['pilot_input']) or
            manifest.get('mpi_ranks') != native.get('mpi_ranks') or manifest.get('openmp_threads') != 1):
        raise ValueError('actual initial-state pilot manifest/receipt/input/executable identity differs')
    physical_data=physical_external_data(manifest,paths['pilot_input'],native)
    shape_spec = q.get('state_shape', {})
    shape = read_state(paths['native_state'], shape_spec.get('logical_surface_element_count'), shape_spec.get('layer_count'))
    if pilot_test=='I19':
        _receipt_state(native, paths['pilot_input'], paths['native_state'], shape, mode, epoch)
    else:
        # I20 explicitly evolves a separately declared thermal history. Its
        # applied epoch may be its own spin-up end, not a cohort's t0. Accept
        # that only through this explicitly selected REAL pilot and its own
        # reviewed proof, never by substituting I19's initial condition.
        verify_surface_initial_state(native,manifest,paths['pilot_input'],paths['native_state'],
            {'state_sha256':q['native_state']['sha256'],'shape':shape,'initialization_mode':mode,
             'ephemeris_time_tdb_s':epoch,'physics_identity_sha256':q['physics_identity_sha256'],
             'executable_sha256':q['executable_sha256'],'qualification_sha256':'pilot-proof-being-authenticated',
             'physical_external_data':physical_data})
    method = q.get('method_qualification', {})
    if not isinstance(method, dict) or method.get('status') != 'READY' or not isinstance(method.get('criteria'),dict) or not method.get('criteria') or not method.get('evidence_files'):
        raise Unavailable('explicit actual-state method criteria/qualification evidence unavailable')
    if not q.get('source_evidence_files'):
        raise Unavailable('actual-state native-source proof unavailable')
    for label, items in (('method', method['evidence_files']), ('source', q['source_evidence_files'])):
        if not isinstance(items, list):
            raise ValueError('initial-state '+label+' evidence must be an explicit descriptor array')
        for index, descriptor in enumerate(items):
            _file(descriptor, base, label+' evidence '+str(index))
    digest = hashlib.sha256(json.dumps(q, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    return {'schema': SCHEMA, 'qualification_sha256': digest, 'qualification': q,
            'approved_utc': approved.isoformat(), 'paths': {key: str(path) for key, path in paths.items()},
            'state_sha256': q['native_state']['sha256'], 'shape': shape,
            'initialization_mode': mode, 'ephemeris_time_tdb_s': epoch,
            'pilot_test_id':pilot_test,
            'physical_external_data':physical_data,
            'physics_identity_sha256': q['physics_identity_sha256'], 'executable_sha256': q['executable_sha256'],
            # These two controls are approved recipe defaults already used by
            # write_case, not physical-law assumptions. They are excluded from
            # the physics projection but must NOT vary silently in a matrix.
            'numerical_controls':{key:_number(setting(application,key,default),key,minimum=0)
                for key,default in (('m2_verification_cohort_packets','32'),
                                    ('m2_verification_source_integral_tolerance','1e-8'))},
            'frozen_surface_elements':None}


def bind_frozen_sites(qualification, sites):
    """Bind the separately reviewed ordered native-site recipe, duplicates kept."""
    if sites is None:
        raise Unavailable('explicit frozen site registry needed to authenticate numerical case controls')
    if not isinstance(sites,list) or not 1<=len(sites)<=48:
        raise ValueError('frozen initial-state/case site recipe must contain 1..48 native sites')
    elements=[_integer(row['logical_surface_element'],'frozen native site element') for row in sites]
    if any(element>=qualification['shape']['logical_surface_element_count'] for element in elements):
        raise ValueError('frozen site element lies outside complete qualified initial thermal grid')
    return {**qualification,'frozen_surface_elements':elements}


def verify_case_initial_state(directory, qualification):
    """Bind one actual native case to the complete reviewed applied t0 state.

    SHA equality to the already fully parsed qualified grid proves the same
    complete shape/depths/temperatures for each file. Every case is re-hashed;
    there is no path/mtime hash cache or physical state-loading shortcut.
    """
    from m2_initial_state import physics_identity_sha256
    from m2_illumination import setting
    from integration_common import parse_receipt
    directory = Path(directory).resolve();application = directory/'application.in'
    receipt_path = directory/'receipt.txt';manifest_path = directory/'run_manifest.json'
    if not all(path.is_file() for path in (application,receipt_path,manifest_path)):
        raise Unavailable('selected native initial-state input/receipt/manifest unavailable: '+str(directory))
    receipt = parse_receipt(receipt_path);manifest = document(manifest_path)
    if 'thermal_initialization_mode' not in receipt or 'initial_thermal_state_csv' not in receipt:
        raise Unavailable('selected native case lacks applied initial-state receipt/export capability: '+str(directory))
    if not qualification.get('physical_external_data'):
        raise Unavailable('qualified pilot physical external-data identities unavailable')
    physical_external_data(manifest,application,receipt,expected=qualification['physical_external_data'])
    if setting(application, 'm2_verification_record_initial_thermal_state', 'false') != 'true':
        raise ValueError('formal native case did not request the actual initial thermal-state record')
    controls=qualification.get('numerical_controls')
    elements=qualification.get('frozen_surface_elements')
    if not isinstance(controls,dict) or elements is None:
        raise Unavailable('reviewed numerical recipe and ordered frozen sites unavailable for applied-state binding')
    for key,value in controls.items():
        if _number(setting(application,key),'actual native '+key,minimum=0)!=value:
            raise ValueError('native case changes reviewed numerical recipe: '+key)
    actual_elements=setting(application,'m2_verification_surface_elements')
    if actual_elements is None or [_integer(value,'actual native site element') for value in actual_elements.split(',')]!=elements:
        raise ValueError('native case changes ordered frozen surface-element recipe')
    state = directory/('receipt.txt'+STATE_SUFFIX)
    if not state.is_file():
        raise Unavailable('actual selected native initial-state file unavailable: '+str(state))
    actual_hash = sha256(state)
    if actual_hash != qualification['state_sha256']:
        raise ValueError('selected native case applied a different complete initial thermal state: '+str(directory))
    physical_hash = physics_identity_sha256(application)
    if physical_hash != qualification['physics_identity_sha256']:
        raise ValueError('selected native staged physics changed from its qualified initial-state problem')
    if manifest.get('physics_identity_sha256') != physical_hash:
        raise ValueError('native case manifest does not bind the actual staged physical input identity')
    if (manifest.get('test_id') != 'I19' or manifest.get('executable', {}).get('sha256') != qualification['executable_sha256'] or
            manifest.get('inputs', {}).get(str(application)) != sha256(application)):
        raise ValueError('selected native initial-state manifest/input/executable identity differs')
    if 'started_utc' not in manifest:
        raise Unavailable('native case chronology needed for pre-campaign state freeze unavailable')
    if _stamp(manifest['started_utc']) <= _stamp(qualification['approved_utc']):
        raise ValueError('native numerical case started before the applied initial-state qualification was frozen')
    _receipt_state(receipt, application, state, qualification['shape'],
                   qualification['initialization_mode'], qualification['ephemeris_time_tdb_s'])
    return {'path': str(state), 'sha256': actual_hash, 'shape': qualification['shape'],
            'shape_verification': 'exact bytes match independently complete-parsed qualified model initial condition',
            'physics_identity_sha256': physical_hash, 'initialization_mode': qualification['initialization_mode'],
            'ephemeris_time_tdb_s': qualification['ephemeris_time_tdb_s'],
            'qualification_sha256': qualification['qualification_sha256']}


def verify_surface_initial_state(receipt, manifest, application, state, qualification):
    """Bind I20's ACTUAL pre-query state to its EXPLICITLY selected proof.

    I20 has a thermal-history exporter rather than a numerical budget receipt,
    so its native field names differ. Its complete grid still must be the
    declared applied initial condition, not an extra unreported spin-up.
    That proof can have a different physical epoch from the cohort experiment;
    it is never chosen automatically from I19. No annual conclusion follows.
    """
    from m2_initial_state import physics_identity_sha256
    required=('initial_thermal_state_file','initial_thermal_surface_elements',
              'initial_thermal_layer_count','initial_thermal_state_tdb_s',
              'thermal_initialization_mode','thermal_initialization_et_s')
    if any(key not in receipt for key in required):
        raise Unavailable('native I20 applied initial-state export capability unavailable')
    state=Path(state).resolve()
    if not state.is_file():
        raise Unavailable('actual I20 initial-state CSV unavailable')
    if not qualification.get('physical_external_data'):
        raise Unavailable('explicit I20 proof lacks physical SPICE/DEM identities')
    physical_external_data(manifest,application,receipt,expected=qualification['physical_external_data'])
    if (Path(str(receipt['initial_thermal_state_file'])).resolve()!=state or
            sha256(state)!=qualification['state_sha256'] or
            receipt['thermal_initialization_mode']!=qualification['initialization_mode'] or
            _number(receipt['thermal_initialization_et_s'],'I20 applied initial epoch')!=qualification['ephemeris_time_tdb_s'] or
            _number(receipt['initial_thermal_state_tdb_s'],'I20 exported initial epoch')!=qualification['ephemeris_time_tdb_s']):
        raise ValueError('fresh I20 state/method/epoch differs from its explicitly qualified applied initial condition')
    shape=qualification['shape']
    if (_integer(receipt['initial_thermal_surface_elements'],'I20 initial surface element count')!=shape['logical_surface_element_count'] or
            _integer(receipt['initial_thermal_layer_count'],'I20 initial thermal layer count')!=shape['layer_count'] or
            physics_identity_sha256(application)!=qualification['physics_identity_sha256'] or
            manifest.get('executable',{}).get('sha256')!=qualification['executable_sha256']):
        raise ValueError('fresh I20 initial-state grid/physical problem/binary differs from qualification')
    # I20's explicit existing spin-up is allowed to end at its separately
    # declared epoch even for legacy uniform startup. Require its ACTUAL
    # convergence metrics, without assigning its epoch to the cohort run.
    if qualification['initialization_mode']=='periodic-spinup':
        from m2_illumination import setting
        if (_integer(receipt.get('thermal_initialization_spinup_completed_cycles'),'I20 native spin-up cycles',minimum=1)<1 or
                _number(receipt.get('thermal_initialization_spinup_cycle_end_et_s'),'I20 spin-up end')!=qualification['ephemeris_time_tdb_s'] or
                _number(receipt.get('thermal_initialization_spinup_cycle_start_et_s'),'I20 spin-up begin')>=qualification['ephemeris_time_tdb_s']):
            raise ValueError('I20 periodic applied state lacks matching actual spin-up interval')
        for key,control in (
                ('thermal_initialization_spinup_maximum_node_temperature_change_k','thermal_spinup_temperature_tolerance_k'),
                ('thermal_initialization_spinup_maximum_column_enthalpy_change_j_m2','thermal_spinup_energy_tolerance_j_m2')):
            if _number(receipt.get(key),'I20 '+key,minimum=0)>_number(setting(application,control),control,minimum=0):
                raise ValueError('I20 initial-state claim contradicts its unchanged production convergence criterion')
    return {'path':str(state),'sha256':qualification['state_sha256'],
            'qualification_sha256':qualification['qualification_sha256'],
            'scope':'actual applied model initial condition, not independent observations'}


def verify_initialization_forecast(forecast, qualification, *, mpi_ranks, case_count):
    """Require measured actual initialization in the total allocation budget.

    The production implementation re-executes initialization for every case;
    no free cached state exists. Measured timing evidence is reviewed input,
    not inferred from numerical residuals or assigned a convenient duration.
    """
    cost = forecast.get('initialization_cost')
    if not isinstance(cost, dict) or cost.get('included_in_allocated_core_hours') is not True or (
            cost.get('reexecuted_for_every_case') is not True) or not cost.get('evidence_files'):
        raise Unavailable('cost forecast does not explicitly include actual per-case thermal initialization')
    if (cost.get('qualification_sha256') != qualification['qualification_sha256'] or
            cost.get('native_state_sha256') != qualification['state_sha256'] or
            cost.get('initialization_mode') != qualification['initialization_mode'] or
            _integer(cost.get('case_count'), 'initialization forecast case count', minimum=1) != case_count or
            _integer(cost.get('measurement_mpi_ranks'), 'initialization measurement MPI ranks', minimum=1) != mpi_ranks):
        raise ValueError('initialization cost forecast refers to a different method/state/layout/case count')
    measured = _number(cost.get('measured_per_case_wall_seconds'), 'measured initialization wall seconds', minimum=0)
    allocated = _number(cost.get('allocated_initialization_core_hours'), 'allocated initialization core-hours', minimum=0)
    total = _number(forecast.get('allocated_core_hours'), 'total allocated core-hours', minimum=0)
    lower = measured*mpi_ranks*case_count/3600
    if measured == 0 or allocated < lower*(1-16*sys.float_info.epsilon) or total < allocated:
        raise ValueError('forecast omits measured initialization cost from the complete allocation')
    base = Path(qualification['paths']['pilot_manifest']).parent
    # Require absolute timing descriptors: a forecast file can be staged far
    # from its pilot. Relative guessing would silently attach another proof.
    for descriptor in cost['evidence_files']:
        if not Path(descriptor.get('path', '')).is_absolute():
            raise ValueError('initialization timing evidence must use explicit absolute paths')
        _file(descriptor, base, 'measured initialization timing evidence')
    return {'measured_initialization_core_hours_lower_bound': lower,
            'allocated_initialization_core_hours': allocated, 'case_count': case_count,
            'qualification_sha256': qualification['qualification_sha256']}


def verify_normalization_physics(normalizations, application, executable, *, base, mpi_ranks):
    """Authenticate the native normalization pilot's ACTUAL production inputs.

    A legacy bounded development pilot can be internally correct but solves a
    different physical problem. Only a cluster-tier preserved-physics pilot
    can supply this formal matrix's normalization. No normalization value is
    fitted or changed here; these are identity and native provenance checks.
    """
    from m2_initial_state import physics_identity_sha256
    from integration_common import parse_receipt
    if normalizations.get('production_input_tier') != 'cluster':
        raise Unavailable('normalization needs a preserved-physics cluster-tier native pilot')
    expected = physics_identity_sha256(application)
    if (normalizations.get('physics_identity_sha256') != expected or
            normalizations.get('template_physics_identity_sha256') != expected):
        raise ValueError('normalization actual/template physics differs from reviewed formal problem')
    files = {key:_file(normalizations.get(key), base, 'normalization '+key)
             for key in ('actual_input','native_manifest','native_receipt')}
    if physics_identity_sha256(files['actual_input']) != expected:
        raise ValueError('native normalization pilot changed actual staged AMR/source/thermal physics')
    manifest=document(files['native_manifest']);receipt=parse_receipt(files['native_receipt'])
    if (manifest.get('test_id') != 'I19' or receipt.get('test_id') != 'I19' or
            manifest.get('mpi_ranks') != mpi_ranks or receipt.get('mpi_ranks') != mpi_ranks or
            manifest.get('openmp_threads') != 1 or
            manifest.get('executable',{}).get('sha256') != sha256(executable) or
            manifest.get('inputs',{}).get(str(files['actual_input'])) != sha256(files['actual_input'])):
        raise ValueError('normalization native manifest/receipt/input/binary/layout identity differs')
    physical_data=physical_external_data(manifest,files['actual_input'],receipt)
    weight=normalizations.get('native_weight_provenance')
    if not isinstance(weight,dict) or not all(key in weight for key in (
            'method','cohort_packets','initial_atoms_per_site','reference_weight_atoms',
            'actual_input_sha256','native_receipt_sha256')):
        raise Unavailable('normalization native initialized-weight provenance unavailable')
    from m2_illumination import setting
    packets=_number(weight['cohort_packets'],'native normalization cohort packet scale',minimum=0)
    if (not isinstance(weight['method'],str) or not weight['method'] or packets==0 or
            packets!=_number(setting(files['actual_input'],'m2_verification_cohort_packets','32'),'declared normalization packets',minimum=0) or
            weight['actual_input_sha256']!=sha256(files['actual_input']) or
            weight['native_receipt_sha256']!=sha256(files['native_receipt'])):
        raise ValueError('normalization initialized-weight provenance is not bound to the actual input/receipt')
    for species,prefix in (('Na','na'),('Ar','ar')):
        atoms=_number(receipt['initial_'+prefix+'_cohort_atoms_per_site'],'native initial normalization atoms',minimum=0)
        recorded=_number(weight['initial_atoms_per_site'][species],'recorded normalization atoms',minimum=0)
        reference=_number(weight['reference_weight_atoms'][species],'recorded ordinary initialized weight',minimum=0)
        if atoms==0 or recorded!=atoms or not math.isclose(reference,atoms/packets,rel_tol=16*sys.float_info.epsilon,abs_tol=0):
            raise ValueError('normalization reference weight is not the actual ordinary initialized stock/packet scale')
    return {'physics_identity_sha256':expected,'physical_external_data':physical_data,
            'native_manifest_sha256':sha256(files['native_manifest']),
            'native_receipt_sha256':sha256(files['native_receipt']),
            'actual_input_sha256':sha256(files['actual_input'])}
