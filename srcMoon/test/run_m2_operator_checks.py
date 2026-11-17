#!/usr/bin/env python3
"""Bounded REAL production I15/I20 exporter checks, not lunar validation.

I15 samples predeclared Cartesian/latitude/epoch fixture queries against the
actual LOLA/SPICE terrain. I20 replays a deliberately 1 ms forcing interval
until the EXISTING thermal boundary/column tolerances hold, then samples two
actual production columns. That is a controlled API/wiring fixture, not an
approximation to a synodic cycle, qualified lunar initial state, annual polar
temperature or independent observation. Formal I15/I20 remain SKIPPED.

The ordinary executable and callbacks perform all geometry and thermal work.
There is no fake predecessor receipt, observation table, or test-only law.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time

from build_m2_qualification import verify_retained_executable
from integration_common import (REPO, manifest, parse_receipt, prepare_na_ar_input,
                                 replace_moon_setting, sha256)
from m1_execution import execute_levels
from m2_illumination import setting
from m2_protocol import without_verification_controls

# Numerical formula-evaluation allowance, fixed before executing the native
# kinetics exporter. This is not a laboratory or observational error budget.
RESIDENCE_REFERENCE_RTOL = 1e-11
INITIAL_STATE_ROUNDOFF_RTOL = 64*sys.float_info.epsilon


class MissingDependency(Exception):
    """A required compiler/runtime/data capability was not available."""


class FixtureFailure(AssertionError):
    """A completed native output missed a declared mathematical invariant."""


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def _put(text, key, value, source):
    """Add an optional fixture control inside the unique application block."""
    if re.search(r'(?m)^\s*'+re.escape(key)+r'\s*=', text):
        return replace_moon_setting(text, key, value, source)
    endings = list(re.finditer(r'(?m)^\s*#end application[^\n]*', text))
    if len(endings) != 1:
        raise ValueError('fixture needs one application block end')
    at = endings[0].start()
    return text[:at]+key+' = '+value+'\n'+text[at:]


def _utc(text):
    # The .in key explicitly declares UTC. This limited calendar helper only
    # constructs the millisecond wiring fixture; CSPICE remains authoritative
    # for native UTC/TDB conversion. Unsupported leap-second spelling is an
    # explicit setup error, never a substituted epoch or observation correction.
    instant = datetime.fromisoformat(text.replace('Z', '+00:00'))
    if instant.tzinfo is None:
        instant = instant.replace(tzinfo=timezone.utc)
    return instant.astimezone(timezone.utc)


def _utc_text(instant):
    return instant.isoformat(timespec='microseconds').replace('+00:00', 'Z')


def _base_input(template, destination):
    """Derive an explicitly bounded fixture; retain original material laws."""
    prepare_na_ar_input(template, destination)
    text = without_verification_controls(destination.read_text())
    # Paths moved into the campaign must retain source-relative interpretation.
    for key in ('spice_path', 'lola_image_file', 'lola_label_file'):
        value = Path(setting(template, key))
        if not value.is_absolute():
            value = template.parent/value
        text = replace_moon_setting(text, key, str(value.resolve()), template)
    for key, value in (('m2_formal_contract_file', 'none'),
                       ('diviner_eclipse_selection', 'training'),
                       ('diviner_eclipse_freeze_file', 'none'),
                       ('m2_illumination_queries_file', 'none'),
                       ('m2_polar_queries_file', 'none'),
                       ('illumination_terrain_solar_model', 'uniform-disc'),
                       ('illumination_disc_radial_points', '4'),
                       ('illumination_disc_azimuth_points', '16'),
                       ('illumination_aberration_correction', 'NONE'),
                       ('surface_temperature_mode', 'thermal-model'),
                       ('thermal_initialization_mode', 'declared-uniform'),
                       ('mesh_output_tecplot', 'false'),
                       ('surface_inventory_checkpoint_file', 'none')):
        text = _put(text, key, value, template)
    return '! SUPPORTING OPERATOR WIRING FIXTURE; NO OBSERVATIONS OR FORMAL QUALIFICATION.\n'+text


def prepare_jobs(amps, template, output, mpi_ranks):
    """Freeze the complete two-fixture workload before any native run."""
    jobs = {};contracts = {}
    for test in ('I15', 'I20'):
        directory = output/test;directory.mkdir();(directory/'runtime').mkdir()
        application = directory/'application.in'
        text = _base_input(template, application)
        query = directory/'geometry_queries.csv'
        if test == 'I15':
            fields = ('query_id', 'site_id', 'latitude_deg', 'longitude_east_deg',
                      'epoch', 'time_system', 'observer_height_m', 'height_direction', 'include_earth')
            rows = [{'query_id': f'geometry-{i}-{j}', 'site_id': f'fixture-site-{i}-{j}',
                     'latitude_deg': latitude, 'longitude_east_deg': longitude,
                     'epoch': 0, 'time_system': 'TDB', 'observer_height_m': 0,
                     'height_direction': 'radial', 'include_earth': 'true'}
                    for i, latitude in enumerate((-80, -30, 30, 80))
                    for j, longitude in enumerate((0, 180))]
            # Keep the eight original TDB=0 geometry fixtures unchanged. A
            # ninth query supplies an INDEPENDENT native UTC->TDB conversion
            # of the actual configured startup epoch. This catches an
            # installed provider whose epoch stayed at its default zero;
            # Python never guesses leap-second/TDB offsets for the check.
            source_start_utc = setting(template, 'simulation_start_utc')
            source_query_id = 'source-provider-startup-epoch'
            rows.append({'query_id': source_query_id, 'site_id': 'fixture-source-start',
                         'latitude_deg': -80, 'longitude_east_deg': 0,
                         'epoch': source_start_utc, 'time_system': 'UTC',
                         'observer_height_m': 0, 'height_direction': 'radial', 'include_earth': 'true'})
            text = _put(text, 'm2_illumination_queries_file', str(query), application)
            contracts[test] = {'scope': 'actual production geometry API; no known-observation truth',
                               'queries': rows, 'ray_count': 64, 'fixed_tdb_geometry_queries': 8,
                               'source_start_utc': source_start_utc,
                               'source_epoch_query_id': source_query_id}
        else:
            # Keep the original 0.01 K/10000 J m^-2 criteria. A millisecond
            # replay is deliberately NOT a lunar spin-up: it exercises the
            # same state reset/convergence/export path at bounded cost only.
            temperature_limit = float(setting(template, 'thermal_spinup_temperature_tolerance_k'))
            energy_limit = float(setting(template, 'thermal_spinup_energy_tolerance_j_m2'))
            if temperature_limit != .01 or energy_limit != 10000:
                raise ValueError('operator fixture requires unchanged approved 0.01 K / 10000 J m^-2 spin-up tolerances')
            start = _utc(setting(template, 'simulation_start_utc'))
            end = start+timedelta(milliseconds=1)
            times = [_utc_text(end), _utc_text(end+timedelta(milliseconds=1))]
            for key, value in (('thermal_spinup_cycle_start_utc', _utc_text(start)),
                               ('thermal_spinup_cycle_end_utc', _utc_text(end)),
                               ('simulation_start_utc', _utc_text(end)),
                               ('thermal_initialization_mode', 'periodic-spinup'),
                               ('thermal_time_step_s', '0.001')):
                text = _put(text, key, value, application)
            fields = ('query_id', 'site_id', 'logical_surface_element', 'epoch',
                      'time_system', 'depth_m', 'sample_weight', 'time_support_id')
            rows = [{'query_id': f'column-{element}-time-{i}-depth-{j}',
                     'site_id': f'fixture-element-{element}', 'logical_surface_element': element,
                     'epoch': epoch, 'time_system': 'UTC', 'depth_m': depth,
                     'sample_weight': 1, 'time_support_id': 'millisecond-wiring-fixture-not-a-lunar-period'}
                    for element in (0, 3000) for i, epoch in enumerate(times)
                    for j, depth in enumerate((0, .02))]
            text = _put(text, 'm2_polar_queries_file', str(query), application)
            contracts[test] = {'scope': 'controlled millisecond state/wiring fixture, not lunar thermal qualification',
                               'queries': rows, 'spinup_cycle_start_utc': _utc_text(start),
                               'spinup_cycle_end_utc': _utc_text(end), 'thermal_time_step_s': .001,
                               'thermal_initialization_mode': 'periodic-spinup',
                               'simulation_start_utc': _utc_text(end),
                               'unchanged_temperature_limit_K': temperature_limit,
                               'unchanged_column_energy_limit_J_m-2': energy_limit,
                               'no_cycle_limit_or_threshold_relaxation': True}
        with query.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
        application.write_text(text)
        command = ['mpiexec', '-n', str(mpi_ranks), str(amps), '--moon-integration-test', test,
                   '--integration-output', str(directory/'receipt.txt'), '-input', str(application)]
        identity = manifest(test, amps, command, int(setting(application, 'random_seed')),
                            application, mpi_ranks, require_scientific_data=False)
        identity.update(execution_scope='supporting-native-operator-wiring', formal_test_status='SKIPPED',
                        formal_reason='explicit analytical/control fixture; I13/observational qualification not supplied',
                        geometry_queries={'path': str(query), 'sha256': sha256(query)},
                        fixture=contracts[test], actual_source_files={str(path): sha256(path) for path in
                            (REPO/'srcMoon/MoonSurfaceEnvironment.cpp', REPO/'srcMoon/MoonPolarDiagnostics.cpp',
                             REPO/'srcMoon/MoonSurfaceTemperature.cpp')})
        identity['inputs'][str(query)] = sha256(query)
        _write(directory/'run_manifest.json', identity)
        _write(directory/'fixture_contract.json', contracts[test])
        jobs[test] = {'directory': directory, 'command': command}
    return jobs, contracts


def _rows(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('native CSV has an absent or duplicate header')
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError('native CSV has missing or excess columns')
        return rows


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError('native operator export requires a numerical value') from error
    if not math.isfinite(value):
        raise ValueError('native operator export contains a nonfinite value')
    return value


def evaluate_illumination(directory, contract):
    """Check native identities, masks, actual PSD hook and startup clock.

    Native source-start TDB is compared to a separate geometry query converted
    from the frozen template UTC by CSPICE itself. The eight TDB0 geometry
    probes remain independent of that initial source epoch.
    """
    receipt = parse_receipt(directory/'receipt.txt')
    rows = [row for path in sorted(directory.glob('receipt.txt.I15.rank*.csv')) for row in _rows(path)]
    expected = {row['query_id']: row for row in contract['queries']}
    actual = {row['query_id']: row for row in rows}
    if len(rows) != len(expected) or len(actual) != len(rows) or set(actual) != set(expected):
        raise ValueError('native I15 did not emit every predeclared geometry fixture once')
    for key, row in actual.items():
        source = expected[key]
        epoch = _finite(row['ephemeris_time_tdb_s'])
        if row['site_id'] != source['site_id'] or row['input_time_system'] != source['time_system'] or (
                _finite(row['latitude_deg']) != source['latitude_deg'] or
                _finite(row['longitude_east_deg']) != source['longitude_east_deg'] or
                int(row['ray_count']) != 64):
            raise ValueError('native I15 geometry/ray-count identity differs from fixture')
        if source['time_system'] == 'TDB':
            if _finite(row['input_epoch']) != float(source['epoch']) or epoch != float(source['epoch']):
                raise ValueError('native I15 explicit TDB geometry epoch differs from fixture')
        elif source['time_system'] == 'UTC':
            if row['input_epoch'] != source['epoch']:
                raise ValueError('native I15 UTC epoch identity differs from frozen template start')
        else:
            raise ValueError('operator fixture requires an explicit UTC/TDB epoch')
        terrain, earth, joint = [_finite(row[name]) for name in
            ('terrain_visible_area_fraction', 'earth_visible_area_fraction', 'jointly_visible_area_fraction')]
        rounding = contract['ray_count']*sys.float_info.epsilon
        if not all(0 <= value <= 1 for value in (terrain, earth, joint)) or joint > min(terrain, earth)+rounding:
            raise FixtureFailure('native I15 joint masks violate probability/set bounds')
        for name, area in (('sampled_any_terrain_visible', terrain), ('sampled_any_jointly_visible', joint)):
            flag = int(row[name])
            if flag not in (0, 1):
                raise ValueError('native I15 visibility flag is not binary')
            if flag != int(area > 0):
                raise FixtureFailure('native I15 occurrence and projected area diagnostics disagree')
        if not 0 <= _finite(row['projected_incidence']) <= 1+rounding:
            raise FixtureFailure('native I15 projected irradiance is outside its geometric bounds')
    hook = receipt['runtime_psd_geometry_hook']
    if hook not in (0, 1):
        raise ValueError('native runtime PSD geometry hook flag is not binary')
    if hook != 1:
        raise FixtureFailure('actual Moon production PSD geometry provider was not installed')
    source_query = actual[contract['source_epoch_query_id']]
    if source_query['input_time_system'] != 'UTC' or source_query['input_epoch'] != contract['source_start_utc']:
        raise ValueError('source startup epoch is not bound to the frozen template UTC')
    expected_epoch = _finite(source_query['ephemeris_time_tdb_s'])
    source_epoch = _finite(receipt['source_begin_epoch_tdb_s'])
    rounding_seconds = 16*sys.float_info.epsilon*max(1, abs(expected_epoch))
    if abs(source_epoch-expected_epoch) > rounding_seconds:
        raise FixtureFailure('actual PSD provider epoch differs from native SPICE conversion of template startup UTC')
    return {'native_geometry_rows': len(rows), 'ray_count_per_query': 64,
            'runtime_psd_geometry_hook': hook, 'source_begin_epoch_tdb_s': source_epoch,
            'expected_source_begin_epoch_tdb_s': expected_epoch,
            'source_epoch_roundoff_bound_s': rounding_seconds,
            'scope_limit': 'bounds/identities only; independent horizon/eclipse oracles remain U09'}


def evaluate_thermal(directory, contract):
    """Check actual production thermal exports without inventing a solution."""
    receipt = parse_receipt(directory/'receipt.txt')
    rows = _rows(directory/'receipt.txt.I20_thermal_model.csv')
    expected = {row['query_id']: row for row in contract['queries']}
    actual = {row['query_id']: row for row in rows}
    if len(rows) != 8 or len(actual) != len(rows) or set(actual) != set(expected) or receipt.get('thermal_model_rows') != 8:
        raise ValueError('native I20 did not emit the declared production state queries once')
    emissivity = float(setting(directory/'application.in', 'thermal_infrared_emissivity'))
    application = directory/'application.in'
    reporting_timescale = float(setting(application, 'cold_trap_retention_timescale_s'))
    activation_ev = float(setting(application, 'sodium_activation_energy_ev'))
    attempt_hz = float(setting(application, 'sodium_attempt_frequency_hz'))
    argon_energy_j_mol = float(setting(application, 'argon_desorption_energy_j_mol'))
    if receipt.get('surface_kinetics_export_schema') != 'native-species-residence-retention-v1':
        raise ValueError('native I20 does not identify the required production surface-kinetics exporter')
    kinetics_fields = ('surface_retention_timescale_s', 'na_surface_residence_time_s',
                       'na_surface_thermal_retention', 'ar_surface_residence_time_s',
                       'ar_surface_thermal_retention')
    for key, row in actual.items():
        query = expected[key]
        if row['site_id'] != query['site_id'] or int(row['logical_surface_element']) != query['logical_surface_element'] or (
                _finite(row['depth_m']) != query['depth_m'] or _finite(row['sample_weight']) != 1 or
                row['time_support_id'] != query['time_support_id'] or row['frame'] != 'MOON_ME_DE421'):
            raise ValueError('native I20 state identity/depth/weight differs from fixture')
        temperature = _finite(row['model_physical_temperature_K'])
        if temperature <= 0 or _finite(row['absorbed_solar_flux_W_m-2']) < 0:
            raise FixtureFailure('native I20 temperature/solar flux is outside its physical domain')
        if query['depth_m'] == 0:
            brightness = _finite(row['model_bolometric_brightness_temperature_K'])
            if not math.isclose(brightness, temperature*emissivity**.25, rel_tol=32*sys.float_info.epsilon):
                raise FixtureFailure('native surface grey observation operator is not epsilon^(1/4) T')
            if _finite(row['surface_retention_timescale_s']) != reporting_timescale:
                raise ValueError('native retention-reporting timescale differs from the actual input')
            # Independent closed-form references use the actual input laws
            # and SI constants (2019 exact h, k, e and N_A). No temperature,
            # residence or source-strength observation is invented here.
            log_references = {
                'na': activation_ev*1.602176634e-19/(1.380649e-23*temperature)-math.log(attempt_hz),
                'ar': math.log(6.62607015e-34/(1.380649e-23*temperature))+
                      argon_energy_j_mol/(8.31446261815324*temperature)}
            for species, logarithm in log_references.items():
                reference = (math.inf if logarithm > math.log(sys.float_info.max)
                             else 0 if logarithm < math.log(sys.float_info.min) else math.exp(logarithm))
                residence = float(row[species+'_surface_residence_time_s'])
                retained = row[species+'_surface_thermal_retention']
                if math.isnan(residence) or residence < 0 or retained not in ('0', '1'):
                    raise ValueError('native surface residence/classification is malformed')
                if not math.isclose(residence, reference, rel_tol=RESIDENCE_REFERENCE_RTOL, abs_tol=0):
                    raise FixtureFailure('native '+species+' residence does not match its independent input-law reference')
                if int(retained) != int(residence >= reporting_timescale):
                    raise FixtureFailure('native surface retention does not match its explicit kinetic timescale')
        elif row['model_bolometric_brightness_temperature_K'] != '':
            raise FixtureFailure('native I20 falsely labels a subsurface kinetic state as brightness')
        elif any(row[name] != '' for name in kinetics_fields):
            raise FixtureFailure('native I20 falsely assigns surface retention to a subsurface layer')
        if not -90 <= _finite(row['model_latitude_deg']) <= 90 or not 0 <= _finite(row['model_longitude_east_deg']) < 360:
            raise FixtureFailure('native I20 model body coordinates are invalid')
    if receipt.get('thermal_spinup_completed_cycles', 0) < 1:
        raise ValueError('native I20 did not execute its production spin-up convergence path')
    initial = _check_applied_initial_columns(directory, receipt, rows)
    return {'native_thermal_rows': len(rows), 'spinup_completed_cycles': receipt['thermal_spinup_completed_cycles'],
            **initial,
            'scope_limit': 'millisecond replay/state export only; not synodic spin-up or annual thermal qualification'}


def _check_applied_initial_columns(directory, receipt, samples):
    """Prove this bounded fixture actually uses its recorded complete t0 state.

    This is a model-state/wiring identity, not an independent physical truth.
    Compare only query samples at the exact native initial epoch; later
    temperatures must evolve normally rather than being replaced by a cache.
    """
    if receipt.get('thermal_initialization_mode') != 'periodic-spinup':
        raise FixtureFailure('native operator did not apply the requested periodic initial state')
    epoch = _finite(receipt['thermal_initialization_et_s'])
    if epoch != _finite(receipt['thermal_spinup_end_tdb_s']):
        raise ValueError('native initial state and declared replay end epochs differ')
    path = directory/'receipt.txt.I20_initial_thermal_state.csv'
    if Path(str(receipt['initial_thermal_state_file'])).resolve() != path.resolve():
        raise ValueError('native initial-state path differs from the selected exporter artifact')
    count_value = _finite(receipt['initial_thermal_surface_elements'])
    layers_value = _finite(receipt['initial_thermal_layer_count'])
    count, layers = int(count_value), int(layers_value)
    if count != count_value or layers != layers_value:
        raise ValueError('native initial-state dimensions must be exact integers')
    if count <= 0 or layers != int(setting(directory/'application.in', 'thermal_layer_count')) or layers <= 0:
        raise ValueError('native initial-state dimensions differ from the actual input')
    rows = _rows(path)
    if len(rows) != count*(layers+1):
        raise ValueError('native initial-state export is not the complete logical grid')
    profiles = {}
    for position, row in enumerate(rows):
        element, offset = divmod(position, layers+1)
        if int(row['logical_surface_element']) != element or int(row['depth_index']) != offset-1:
            raise ValueError('native initial-state rows are duplicated, incomplete or unordered')
        depth, temperature = _finite(row['depth_m']), _finite(row['kinetic_temperature_K'])
        profile = profiles.setdefault(element, [])
        if temperature <= 0 or (depth != 0 if offset == 0 else depth <= profile[-1][0]):
            raise ValueError('native initial profile has invalid kinetic temperature or depth ordering')
        profile.append((depth, temperature))
    compared = 0
    for sample in samples:
        if _finite(sample['ephemeris_time_tdb_s']) != epoch:
            continue
        profile = profiles[int(sample['logical_surface_element'])]
        depth = _finite(sample['depth_m'])
        if depth == 0:
            expected = profile[0][1]
        else:
            bracket = next(((a, b) for a, b in zip(profile, profile[1:]) if a[0] <= depth <= b[0]), None)
            if bracket is None:
                raise ValueError('initial query lies outside the recorded native model depth support')
            a, b = bracket
            expected = a[1]+(b[1]-a[1])*(depth-a[0])/(b[0]-a[0])
        if not math.isclose(_finite(sample['model_physical_temperature_K']), expected,
                            rel_tol=INITIAL_STATE_ROUNDOFF_RTOL, abs_tol=0):
            raise FixtureFailure('native t0 query differs from the actual exported/applied model initial state')
        compared += 1
    if compared == 0:
        raise ValueError('operator fixture did not sample its exact actual initial epoch')
    return {'initial_thermal_rows': len(rows), 'initial_queries_compared': compared,
            'initial_thermal_state_sha256': sha256(path), 'initialization_mode': 'periodic-spinup'}


def run_checks(amps, template, output, *, mpi_ranks=8, wall_seconds=None):
    """Run both actual exporters under the caller's remaining resource budget."""
    amps, template, output = Path(amps).resolve(), Path(template).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('operator-check output must be new or empty')
    output.mkdir(parents=True, exist_ok=True)
    results = []
    try:
        if mpi_ranks < 1:
            raise ValueError('positive MPI layout required')
        if not amps.is_file() or not template.is_file() or shutil.which('mpiexec') is None:
            raise MissingDependency('actual executable, explicit input or MPI launcher unavailable')
        for key in ('lola_image_file', 'lola_label_file'):
            configured = Path(setting(template, key))
            if not configured.is_absolute():
                configured = template.parent/configured
            if not configured.is_file():
                raise MissingDependency('configured native LOLA dependency unavailable: '+str(configured))
        if setting(template, 'surface_geometry') != 'lola':
            raise ValueError('operator fixture must use the actual LOLA production terrain')
        verify_retained_executable(amps)
        deadline = os.environ.get('AMPS_MOON_WORK_DEADLINE')
        remaining = max(.01, float(deadline)-time.monotonic()) if deadline else 3600.0
        if wall_seconds is not None:
            if not math.isfinite(wall_seconds) or wall_seconds <= 0:
                raise ValueError('positive finite remaining runtime required')
            remaining = min(remaining, wall_seconds)
        jobs, contracts = prepare_jobs(amps, template, output, mpi_ranks)
        monitor = output/'monitor';monitor.mkdir()
        codes = execute_levels(jobs, monitor, 1, remaining, memory_gib=200000000000/1024**3,
                               render_previews=False)
        for test, evaluate in (('I15', evaluate_illumination), ('I20', evaluate_thermal)):
            directory = jobs[test]['directory']
            result = {'id': test, 'check': test+'-production-operator-wiring', 'status': 'ERROR',
                      'scope': 'supporting-native-operator-wiring', 'formal_test_status': 'SKIPPED',
                      'formal_reason': 'analytical/control fixture cannot qualify M1/M2 or observations',
                      'runner': str(Path(__file__)), 'log': str(directory/'stdout.log')}
            try:
                if codes.get(test) != 0:
                    raise RuntimeError('native '+test+' exited '+str(codes.get(test))+'; inspect retained streams/resource monitor')
                native = parse_receipt(directory/'receipt.txt')
                if native.get('test_id') != test or native.get('mpi_ranks') != mpi_ranks:
                    raise ValueError('native test/MPI receipt identity disagrees with actual launch')
                result.update(status='PASS', metrics=evaluate(directory, contracts[test]),
                              reason='actual production exporter completed and its declared wiring checks passed')
            except FixtureFailure as error:
                result.update(status='FAIL', reason=str(error))
            except (OSError, ValueError, KeyError, RuntimeError) as error:
                result['reason'] = str(error)
            _write(directory/'result.json', result);results.append(result)
    except MissingDependency as error:
        results.append({'check': 'operator-wiring-dependencies', 'status': 'SKIPPED',
                        'scope': 'supporting-native-operator-wiring', 'reason': str(error),
                        'runner': str(Path(__file__)), 'log': str(output)})
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        results.append({'check': 'operator-wiring-infrastructure', 'status': 'ERROR',
                        'scope': 'supporting-native-operator-wiring', 'reason': str(error),
                        'runner': str(Path(__file__)), 'log': str(output)})
    counts = {status: sum(row['status'] == status for row in results)
              for status in ('PASS', 'FAIL', 'SKIPPED', 'ERROR')}
    summary = {'schema': 'srcMoon-M2-native-operator-checks-v1',
               'status': 'ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'SKIPPED' if counts['SKIPPED'] else 'PASS',
               'scope': 'supporting-native-operator-wiring', 'counts': counts, 'results': results,
               'formal_i15_status': 'SKIPPED', 'formal_i20_status': 'SKIPPED',
               'observational_scoring': False, 'assumed_observations': False}
    _write(output/'operator_checks_summary.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--amps', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--mpi-ranks', type=int, default=8)
    parser.add_argument('--test-runtime-seconds', type=float)
    args = parser.parse_args()
    result = run_checks(args.amps, args.input, args.output_dir, mpi_ranks=args.mpi_ranks,
                        wall_seconds=args.test_runtime_seconds)
    print(json.dumps({'status': result['status'], 'counts': result['counts'],
                      'summary': str(args.output_dir/'operator_checks_summary.json')}))
    return {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}[result['status']]


if __name__ == '__main__':
    raise SystemExit(main())
