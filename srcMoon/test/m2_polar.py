#!/usr/bin/env python3
"""Independent I20 polar collocation, aggregation and comparison.

This module owns evidence handling and arithmetic, not lunar physics. Native
AMPS exports the actual thermal state and passive cohorts. D02/D03 provide
independent, labelled evidence. Their time support, depth, temperature kind,
selection and acceptance must be explicitly reviewed before residuals exist.
Missing metadata stays SKIPPED; no default annual interval, polar-temperature
tolerance, uncertainty, model renormalization or approval is supplied here.

The pure selection/aggregation functions are also exercised by infrastructure
tests with explicitly named arithmetic fixtures. Those fixtures never enter
the formal data path and never provide observational validation evidence.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
CANONICAL_DATA = Path('/data/vtenishe/moon_validation_data')
SCHEMA = 'srcMoon-I20-polar-contract-v1'
FRAME = 'MOON_ME_DE421'
PREDECESSORS = ('I13', 'I14', 'I15', 'I16', 'I17', 'I18', 'I19')
# These are definitions of the native products, not interchangeable aliases.
QUANTITIES = {
    'any_solar_disc_visibility_occurrence_fraction': ('1', 0.0),
    'permanent_shadow_classification': ('1', 0.0),
    'annual_maximum_surface_temperature_K': ('K', 0.0),
    'annual_average_temperature_at_2cm_K': ('K', 0.02),
}
COMPARISON_FIELDS = ('comparison_id', 'x_value', 'x_unit', 'reference_value',
                     'model_value', 'y_unit', 'uncertainty_value', 'selection')
KINETICS_EXPORT_SCHEMA = 'native-species-residence-retention-v1'
SURFACE_KINETICS_FIELDS = ('surface_retention_timescale_s',
                          'na_surface_residence_time_s', 'na_surface_thermal_retention',
                          'ar_surface_residence_time_s', 'ar_surface_thermal_retention')


class MissingEvidence(Exception):
    """A scientific prerequisite is absent; this is not an evaluated failure."""


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _number(value, name, *, minimum=None, maximum=None):
    if isinstance(value, bool):
        raise ValueError(name + ' must be a finite numeric value, not Boolean')
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(name + ' must be numeric') from error
    if not math.isfinite(result) or (minimum is not None and result < minimum) or (
            maximum is not None and result > maximum):
        raise ValueError(name + ' is nonfinite or outside its explicit domain')
    return result


def _integer(value, name, *, minimum=0):
    number = _number(value, name, minimum=minimum)
    if number != int(number):
        raise ValueError(name + ' must be an integer')
    return int(number)


def _timestamp(text):
    # Naive times are never guessed to be UTC. Review/export chronology is
    # provenance, distinct from SPICE science-epoch interpretation.
    result = datetime.fromisoformat(str(text).replace('Z', '+00:00'))
    if result.tzinfo is None:
        raise ValueError('provenance timestamp must declare a timezone')
    return result


def _json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('duplicate JSON evidence key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('nonfinite JSON evidence constant: ' + value)
    with Path(path).open(encoding='utf-8') as stream:
        result = json.load(stream, object_pairs_hook=unique, parse_constant=invalid)
    if not isinstance(result, dict):
        raise ValueError('JSON evidence must be an object: ' + str(path))
    return result


def _descriptor(record, base, label):
    """Resolve and hash one immutable descriptor without interpreting values."""
    if record is None:
        raise MissingEvidence(label + ' path/SHA-256 descriptor unavailable')
    if not isinstance(record, dict) or not isinstance(record.get('path'), str) or not record['path']:
        raise ValueError(label + ' descriptor has an invalid explicit path')
    expected = record.get('sha256')
    if not isinstance(expected, str) or not re.fullmatch('[0-9a-f]{64}', expected):
        raise ValueError(label + ' needs a literal lowercase SHA-256')
    path = Path(record['path'])
    if not path.is_absolute():
        path = Path(base) / path
    path = path.resolve()
    if not path.is_file():
        raise MissingEvidence(label + ' file unavailable: ' + str(path))
    if sha256(path) != expected:
        raise ValueError(label + ' SHA-256 mismatch: ' + str(path))
    return path


def _review(document, label):
    review = document.get('review', {})
    if not isinstance(review, dict):
        raise ValueError(label + ' review must be an object')
    if document.get('status') not in ('APPROVED', 'READY') or not (
            document.get('selection_frozen') is True and
            review.get('approved') is True and review.get('reviewer') and
            review.get('approved_utc') and
            review.get('before_model_evaluation') is True):
        raise MissingEvidence(label + ' is not explicitly reviewed/frozen before model evaluation')
    return _timestamp(review['approved_utc'])


def resolve_contract(path):
    """Use the same reviewed direct/bundle dispatch as other formal M2 IDs."""
    from m2_contracts import resolve_contract as resolve_dispatch
    return resolve_dispatch(Path(path).resolve(), 'I20', unavailable=MissingEvidence)


def _csv(path, *, allow_empty=False):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('CSV has missing or duplicate field names: ' + str(path))
        rows = list(reader)
    if (not rows and not allow_empty) or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('CSV is empty or has malformed rows: ' + str(path))
    return rows


def select_sites(candidates, *, maximum_distance_m):
    """Apply the APPROVED metadata-only native-pixel stratified selection.

    The caller's collocation radius must be independently reviewed. No model
    value/residual is accepted as an input field for the selection. Real gaps
    and shared logical elements survive: three sites per available stratum do
    not imply three independent model cells or a representative global sample.
    """
    distance_limit = _number(maximum_distance_m, 'collocation distance', minimum=0)
    protocol = _json(HERE / 'M2_I19_I20_PROTOCOL.json')['proposed_site_selection']
    bands = protocol['absolute_latitude_band_edges_deg']
    per_stratum = protocol['sites_per_band_class']
    groups = {(h, i, p): [] for h in ('N', 'S') for i in range(len(bands)-1) for p in (0, 1)}
    seen = set()
    for row in candidates:
        if any('residual' in key.lower() or key.startswith('model_value') for key in row):
            raise ValueError('native-site selection cannot consume model residual/value fields')
        latitude = _number(row['latitude_deg'], 'native latitude', minimum=-90, maximum=90)
        _number(row['longitude_east_deg'], 'native east longitude', minimum=0, maximum=360)
        hemisphere = 'N' if latitude >= 0 else 'S'
        if row['hemisphere'] != hemisphere:
            raise ValueError('native latitude/hemisphere conflict')
        native = (str(row['native_product_id']),
                  _integer(row['native_line'], 'native line', minimum=1),
                  _integer(row['native_sample'], 'native sample', minimum=1))
        if native in seen:
            raise ValueError('duplicate native pixel in selection candidates')
        seen.add(native)
        psr = _integer(row['permanent_shadow'], 'native PSR')
        edge = _integer(row['edge_ambiguous'], 'native edge ambiguity')
        if psr not in (0, 1) or edge not in (0, 1):
            raise ValueError('PSR/edge values must be native binary flags')
        _integer(row['logical_surface_element'], 'logical model element')
        distance = _number(row['collocation_distance_m'], 'native collocation distance', minimum=0)
        if edge or row['d03_quality_flag'] != 'GOOD' or distance > distance_limit:
            continue
        band = next((i for i in range(len(bands)-1)
                     if bands[i] <= abs(latitude) < bands[i+1] or
                     (i == len(bands)-2 and abs(latitude) == bands[-1])), None)
        if band is None:
            continue  # Preserve unavailable latitude coverage, do not extend it.
        digest = hashlib.sha256(json.dumps(list(native), separators=(',', ':'),
                                         ensure_ascii=False).encode()).hexdigest()
        groups[(hemisphere, band, psr)].append((digest, native, dict(row)))
    selected, coverage = [], []
    for stratum, group in groups.items():
        group.sort(key=lambda item: (item[0], item[1]))
        chosen = group[:per_stratum]
        coverage.append({'hemisphere': stratum[0], 'latitude_band_index': stratum[1],
                         'native_psr': stratum[2], 'available': len(group),
                         'selected': len(chosen), 'missing': max(0, per_stratum-len(chosen))})
        for digest, native, row in chosen:
            row['selection_sha256'] = digest
            selected.append(row)
    if len(selected) > protocol['maximum_sites']:
        raise ValueError('native site registry exceeds the approved maximum')
    return {'sites': selected, 'strata': coverage,
            'shared_model_support': len({str(r['logical_surface_element']) for r in selected}) < len(selected),
            'interpretation': 'coverage probes; shared support is retained, not independent global sampling'}


def prepare_site_registry(candidate_file, collocation_contract, destination):
    """Create a REVIEW CANDIDATE only; never grant scientific approval.

    Candidates contain native metadata and an explicit mapping only, with no
    AMPS residuals. A later investigator-reviewed registry must preserve the
    exact selected row identities, input hashes and genuine empty strata.
    """
    candidate_file = Path(candidate_file).resolve()
    collocation_contract = Path(collocation_contract).resolve()
    contract = _json(collocation_contract)
    _review(contract, 'collocation contract')
    selection = select_sites(_csv(candidate_file), maximum_distance_m=contract['maximum_distance_m'])
    result = dict(selection, schema='srcMoon-I20-native-site-registry-v1',
                  status='REVIEW_REQUIRED_NOT_FROZEN', selection_frozen=False,
                  observational_scoring=False,
                  candidates={'path': str(candidate_file), 'sha256': sha256(candidate_file)},
                  collocation_contract={'path': str(collocation_contract), 'sha256': sha256(collocation_contract)})
    _write_new_json(Path(destination), result)
    return result


def aggregate_thermal_samples(samples, supports, sites):
    """Aggregate MODEL samples under explicitly qualified native operators.

    Input samples come from MoonPolarDiagnostics. A support supplies exact
    query IDs/weights per site, interval and depth, quantity/temperature kind,
    and a qualified sampling proof. There is no default annual window or
    interpolation of observations. Missing a required model sample is ERROR.
    """
    index = {}
    for row in samples:
        query = row['query_id']
        if query in index:
            raise ValueError('duplicate native thermal query_id')
        index[query] = row
    result = []
    for support in supports:
        quantity = support['quantity']
        if quantity not in ('annual_maximum_surface_temperature_K', 'annual_average_temperature_at_2cm_K'):
            raise ValueError('unsupported thermal comparison quantity')
        if support.get('coverage_qualified') is not True or not support.get('sampling_evidence'):
            raise MissingEvidence('thermal temporal sampling/coverage is not independently qualified')
        expected_aggregation = 'maximum' if quantity.startswith('annual_maximum') else 'weighted-mean'
        expected_depth = QUANTITIES[quantity][1]
        if support.get('aggregation') != expected_aggregation or _number(
                support['depth_m'], 'operator depth', minimum=0) != expected_depth:
            raise ValueError('PRP quantity/aggregation/depth definition mismatch')
        kind = support['temperature_kind']
        if kind not in ('physical_temperature_K', 'bolometric_brightness_temperature_K') or (
                expected_depth > 0 and kind != 'physical_temperature_K'):
            raise ValueError('subsurface kinetic and surface brightness temperature cannot be interchanged')
        low = _number(support['start_tdb_s'], 'support start')
        high = _number(support['stop_tdb_s'], 'support stop')
        if high <= low or support['endpoint_convention'] not in ('closed', 'left-closed-right-open'):
            raise ValueError('invalid explicit thermal temporal support')
        for site_id, queries in support['queries_by_site'].items():
            if site_id not in sites or not queries or len(queries) != len(set(queries)):
                raise ValueError('unknown site or empty/duplicate query membership')
            values, weights = [], []
            for query_id in queries:
                if query_id not in index:
                    raise ValueError('required production thermal query missing: ' + query_id)
                row = index[query_id]
                epoch = _number(row['ephemeris_time_tdb_s'], 'model sample epoch')
                if not (low <= epoch <= high) or (
                        support['endpoint_convention'] == 'left-closed-right-open' and epoch == high):
                    raise ValueError('model sample outside the frozen temporal support')
                if row['site_id'] != site_id or row['time_support_id'] != support['time_support_id'] or (
                        row['frame'] != FRAME or _number(row['depth_m'], 'sample depth') != expected_depth or
                        _integer(row['logical_surface_element'], 'sample element') !=
                        _integer(sites[site_id]['logical_surface_element'], 'registry element')):
                    raise ValueError('model thermal sample site/frame/depth/support identity mismatch')
                values.append(_number(row['model_' + kind], 'native model temperature', minimum=0))
                actual_weight = _number(row['sample_weight'], 'native sample weight', minimum=0)
                frozen_weight = _number(support['sample_weights'][query_id], 'frozen sample weight', minimum=0)
                if actual_weight != frozen_weight:
                    raise ValueError('model sample weight differs from the frozen operator')
                weights.append(actual_weight)
            if not sum(weights) > 0:
                raise ValueError('thermal temporal operator has no positive measure')
            value = (max(v for v, w in zip(values, weights) if w > 0)
                     if expected_aggregation == 'maximum' else
                     math.fsum(v*w for v, w in zip(values, weights))/math.fsum(weights))
            result.append({'site_id': site_id, 'quantity': quantity, 'units': 'K',
                           'frame': FRAME, 'time_support_id': support['time_support_id'],
                           'depth_m': expected_depth, 'model_value': value,
                           'temperature_kind': kind,
                           'logical_surface_element': sites[site_id]['logical_surface_element']})
    return result


def compare_rows(reference, model, sites, criteria):
    """Pure matched-row arithmetic; it grants no data/runtime qualification.

    Absolute comparison limits MUST be supplied by an externally reviewed
    pre-residual contract. Native missing uncertainties remain empty and are
    not replaced by zero or a quantization-derived observational error bar.
    """
    def identity(row):
        quantity = row['quantity']
        if quantity not in QUANTITIES:
            raise ValueError('unknown polar observable: ' + quantity)
        units, depth = QUANTITIES[quantity]
        if row['units'] != units or row['frame'] != FRAME or _number(row['depth_m'], 'row depth') != depth:
            raise ValueError('quantity/unit/frame/depth mismatch')
        if row['site_id'] not in sites or not row['time_support_id']:
            raise ValueError('unknown native site or missing explicit temporal support')
        return row['site_id'], quantity
    references, predictions = {}, {}
    for row in reference:
        key = identity(row)
        if key in references:
            raise ValueError('duplicate independent reference site/quantity')
        if row['quality_flag'] != 'GOOD' or not row.get('source_product_id') or not row.get('source_record'):
            raise ValueError('reference row lacks a native good-quality provenance identity')
        references[key] = row
    for row in model:
        key = identity(row)
        if key in predictions:
            raise ValueError('duplicate production model site/quantity')
        if _integer(row['logical_surface_element'], 'model element') != _integer(
                sites[key[0]]['logical_surface_element'], 'registered element'):
            raise ValueError('production logical support differs from frozen native-site mapping')
        predictions[key] = row
    required = {(site_id, quantity) for site_id in sites for quantity in QUANTITIES}
    if set(references) != required or set(predictions) != required:
        raise ValueError('missing/extra native site/quantity rows; no subset or gap filling is allowed')
    comparisons, diagnostics, failures = [], [], []
    for index, key in enumerate(sorted(required)):
        observed, predicted = references[key], predictions[key]
        if observed['time_support_id'] != predicted['time_support_id']:
            raise ValueError('model and reference time support differ')
        if QUANTITIES[key[1]][0] == 'K' and observed.get('temperature_kind') != predicted.get('temperature_kind'):
            raise ValueError('model and native temperature kind differ')
        actual = _number(predicted['model_value'], 'model value', minimum=0)
        expected = _number(observed['reference_value'], 'reference value', minimum=0)
        if key[1] in ('permanent_shadow_classification', 'any_solar_disc_visibility_occurrence_fraction'):
            if actual > 1 or expected > 1:
                raise ValueError('illumination/classification must be bounded in [0,1]')
        if key[1] == 'permanent_shadow_classification' and (actual not in (0, 1) or expected not in (0, 1)):
            raise ValueError('PSR classifications must be binary')
        criterion = criteria.get(key[1])
        if not criterion or 'maximum_absolute_error' not in criterion:
            raise MissingEvidence('reviewed acceptance criterion missing for ' + key[1])
        limit = _number(criterion['maximum_absolute_error'], 'reviewed comparison limit', minimum=0)
        uncertainty = observed.get('uncertainty_value', '')
        if uncertainty != '':
            uncertainty = _number(uncertainty, 'documented native uncertainty', minimum=0)
        residual = actual-expected
        passed = abs(residual) <= limit
        diagnostics.append({'site_id': key[0], 'quantity': key[1], 'time_support_id': observed['time_support_id'],
                            'reference_value': expected, 'model_value': actual, 'residual': residual,
                            'maximum_absolute_error': limit, 'passed': passed,
                            'source_product_id': observed['source_product_id'], 'source_record': observed['source_record'],
                            'logical_surface_element': sites[key[0]]['logical_surface_element']})
        if not passed:
            failures.append({'site_id': key[0], 'quantity': key[1], 'reason': 'frozen comparison criterion missed'})
        comparisons.append({'comparison_id': 'I20_' + key[1], 'x_value': index,
                            'x_unit': 'frozen-site-quantity-row', 'reference_value': expected,
                            'model_value': actual, 'y_unit': observed['units'],
                            'uncertainty_value': uncertainty, 'selection': 'verification'})
    return {'status': 'FAIL' if failures else 'PASS', 'scope': 'matched-row arithmetic only; provenance gates external',
            'comparisons': comparisons, 'diagnostics': diagnostics, 'failures': failures,
            'evaluated_rows': len(comparisons)}


def _load_native_illumination(illumination_rows, sites):
    """Validate distinct benchmark visibility and TERRAIN-only PSR evidence.

    An Earth eclipse can suppress the benchmark's joint visibility without
    turning a terrain-visible site into a geometric permanent-shadow region.
    The new terrain occurrence is mandatory; inferring it from a joint mask
    would erase precisely the distinction this diagnostic must preserve.
    """
    observed = {}
    for row in illumination_rows:
        site_id = row['site_id']
        if site_id not in sites or site_id in observed:
            raise ValueError('unknown/duplicate I15 model site in polar assembly')
        if row['frame'] != FRAME or not row['time_support_id'] or _integer(
                row['logical_surface_element'], 'I15 logical support') != _integer(
                    sites[site_id]['logical_surface_element'], 'frozen logical support'):
            raise ValueError('I15 model geometry differs from the polar native-site registry')
        count = _integer(row['epoch_count'], 'I15 sampled epoch count', minimum=1)
        occurrence = _number(row['any_disc_occurrence_fraction'], 'I15 occurrence', minimum=0, maximum=1)
        if 'terrain_any_disc_occurrence_fraction' not in row:
            raise MissingEvidence('native I15 terrain-only occurrence export unavailable')
        terrain = _number(row['terrain_any_disc_occurrence_fraction'],
                          'I15 terrain-only occurrence', minimum=0, maximum=1)
        psr = _integer(row['geometric_psr'], 'I15 geometric PSR')
        if psr not in (0, 1) or bool(psr) != (terrain == 0) or occurrence > terrain:
            raise ValueError('I15 geometric PSR/benchmark visibility disagrees with terrain-only model rays')
        # A nonzero occurrence must arise from at least one actual sample;
        # this is an arithmetic count identity, not an invented error bound.
        for fraction in (occurrence, terrain):
            if not math.isclose(fraction*count, round(fraction*count), rel_tol=0,
                                abs_tol=32*sys.float_info.epsilon*count):
                raise ValueError('I15 occurrence is not an equal-epoch frequency')
        observed[site_id] = row
    if set(observed) != set(sites):
        raise ValueError('I15 model export is missing required frozen native sites')
    return observed


def assemble_model_rows(thermal_samples, illumination_rows, supports, sites):
    """Build separate model observables without equating shadow and retention.

    Benchmark occurrence retains its declared Earth policy. Geometric PSR
    derives solely from native terrain rays, not an external PSR label or an
    eclipse. Annual thermal values use independently qualified time supports.
    """
    rows = aggregate_thermal_samples(thermal_samples, supports, sites)
    illumination = _load_native_illumination(illumination_rows, sites)
    for site_id, row in illumination.items():
        occurrence = _number(row['any_disc_occurrence_fraction'], 'I15 occurrence')
        psr = _integer(row['geometric_psr'], 'I15 geometric PSR')
        for quantity, value in (('any_solar_disc_visibility_occurrence_fraction', occurrence),
                                ('permanent_shadow_classification', psr)):
            rows.append({'site_id': site_id, 'quantity': quantity, 'units': '1',
                         'frame': FRAME, 'time_support_id': row['time_support_id'],
                         'depth_m': 0.0, 'model_value': value,
                         'temperature_kind': '',
                         'logical_surface_element': sites[site_id]['logical_surface_element']})
    return rows


def surface_retention_samples(samples, supports, sites, illumination_rows, *,
                              retention_timescale_s, native_export_schema=None):
    """Report native surface kinetics; do not implement a Python release law.

    Residence times originate in the SAME production kernel as desorption.
    The only evaluated kinetic invariant here is its diagnostic comparison to
    the explicit timescale. A finite instantaneous classification is never a
    statement about irreversible capture, annual retention, or an observation.
    """
    timescale = _number(retention_timescale_s, 'surface retention timescale', minimum=0)
    if timescale == 0:
        raise ValueError('surface retention timescale must be strictly positive')
    if native_export_schema is not None and native_export_schema != KINETICS_EXPORT_SCHEMA:
        raise ValueError('unsupported declared native surface-kinetics export schema')
    illumination = _load_native_illumination(illumination_rows, sites)
    required = {}
    for support in supports:
        if support['quantity'] != 'annual_maximum_surface_temperature_K':
            continue
        for site_id, queries in support['queries_by_site'].items():
            if site_id not in sites:
                raise ValueError('surface kinetics support contains an unknown native site')
            for query in queries:
                identity = (site_id, support['time_support_id'],
                            _number(support['start_tdb_s'], 'surface support start'),
                            _number(support['stop_tdb_s'], 'surface support stop'),
                            support['endpoint_convention'])
                if identity[3] <= identity[2] or identity[4] not in ('closed', 'left-closed-right-open'):
                    raise ValueError('invalid explicit surface-kinetics temporal support')
                if query in required and required[query] != identity:
                    raise ValueError('surface query has conflicting native temporal support')
                required[query] = identity
    if not required:
        raise MissingEvidence('no explicit native surface queries for thermal-retention reporting')
    rows, failures, seen = [], [], set()
    for sample in samples:
        absent = set(SURFACE_KINETICS_FIELDS)-set(sample)
        if absent:
            if native_export_schema == KINETICS_EXPORT_SCHEMA:
                raise ValueError('declared native surface-kinetics CSV is missing fields: '+','.join(sorted(absent)))
            raise MissingEvidence('native surface-residence/retention export unavailable')
        query = sample['query_id']
        depth = _number(sample['depth_m'], 'native kinetics sample depth', minimum=0)
        if depth != 0:
            if any(sample[key] != '' for key in SURFACE_KINETICS_FIELDS):
                raise ValueError('subsurface samples must not claim surface residence/retention')
            continue
        if query not in required or query in seen:
            raise ValueError('unexpected or duplicate native surface-kinetics query')
        seen.add(query)
        site_id, support_id, low, high, endpoints = required[query]
        epoch = _number(sample['ephemeris_time_tdb_s'], 'native kinetics TDB epoch')
        if (sample['site_id'] != site_id or sample['time_support_id'] != support_id or
                sample['frame'] != FRAME or _integer(sample['logical_surface_element'], 'kinetics element') !=
                _integer(sites[site_id]['logical_surface_element'], 'frozen kinetics element') or
                not low <= epoch <= high or (endpoints == 'left-closed-right-open' and epoch == high)):
            raise ValueError('surface-kinetics site/frame/element/time/depth identity differs')
        if _number(sample['surface_retention_timescale_s'], 'native retention timescale', minimum=0) != timescale:
            raise ValueError('native surface-retention timescale differs from the reviewed input/protocol')
        temperature = _number(sample['model_physical_temperature_K'], 'native kinetic surface temperature', minimum=0)
        if temperature == 0:
            raise ValueError('native kinetic surface temperature must be positive')
        for species, prefix in (('Na', 'na'), ('Ar', 'ar')):
            text = sample[prefix+'_surface_residence_time_s']
            # CSPICE/CSV and JSON have different representations: retain the
            # literal string `inf` in artifacts so no nonfinite JSON number
            # is emitted. +infinity is the physical zero-release limit.
            residence = math.inf if text == 'inf' else _number(text, 'native species residence time', minimum=0)
            retained = _integer(sample[prefix+'_surface_thermal_retention'], 'native thermal-retention classification')
            if retained not in (0, 1):
                raise ValueError('native thermal-retention classification is not binary')
            if bool(retained) != (residence >= timescale):
                failures.append({'query_id': query, 'species': species,
                                 'reason': 'native thermal-retention label disagrees with residence/time-scale invariant'})
            sky = illumination[site_id]
            rows.append({'query_id': query, 'site_id': site_id, 'species': species,
                         'logical_surface_element': sample['logical_surface_element'],
                         'ephemeris_time_tdb_s': epoch, 'depth_m': 0, 'frame': FRAME,
                         'thermal_time_support_id': support_id,
                         'sample_weight': _number(sample['sample_weight'], 'native surface sample weight', minimum=0),
                         'kinetic_surface_temperature_K': temperature,
                         'surface_residence_time_s': 'inf' if math.isinf(residence) else residence,
                         'surface_retention_timescale_s': timescale, 'thermal_retention': retained,
                         'illumination_time_support_id': sky['time_support_id'],
                         'terrain_only_geometric_psr': _integer(sky['geometric_psr'], 'geometric PSR'),
                         'terrain_any_disc_occurrence_fraction': _number(sky['terrain_any_disc_occurrence_fraction'], 'terrain occurrence'),
                         'benchmark_any_disc_occurrence_fraction': _number(sky['any_disc_occurrence_fraction'], 'benchmark occurrence'),
                         'permanent_retention_established': False})
    if seen != set(required):
        raise ValueError('missing required native surface-kinetics query')
    return {'status': 'FAIL' if failures else 'PASS', 'rows': rows, 'failures': failures,
            'scope': 'native instantaneous surface kinetics and separate finite-support terrain PSR; not permanent trapping'}


def bind_cohort_sites(rows, hazards, sites, *, declared_elements, horizon_s,
                      initial_epoch_s, seed_index, expected_initial_atoms=None):
    """Bind actual cohort labels to the ordered frozen native-site registry.

    Each original site has a cohort, even if several sites share one logical
    model element. Native elapsed times and Na/Ar identities are validated
    before the existing independent hazard integrator is called. No hazard,
    temperature, missing interval, or observational value is manufactured.
    """
    ordered = list(sites)
    elements = [_integer(sites[site]['logical_surface_element'], 'frozen cohort element') for site in ordered]
    if [_integer(value, 'declared cohort element') for value in declared_elements] != elements:
        raise ValueError('ordered native I19 cohort elements differ from frozen I20 native sites')
    horizon = _number(horizon_s, 'native cohort horizon', minimum=0)
    start = _number(initial_epoch_s, 'native cohort initial TDB epoch')
    seed = _integer(seed_index, 'native cohort seed index')
    if not rows or horizon == 0:
        raise ValueError('native cohort history/horizon is empty')
    seen, groups, diagnostics, initial_by_group = set(), {}, [], {}
    for row in rows:
        species = row['species']
        cohort = _integer(row['cohort'], 'native cohort label')
        if species not in ('Na', 'Ar') or cohort >= len(ordered):
            raise ValueError('native cohort species/label is outside frozen Na/Ar site family')
        if _integer(row['element'], 'native cohort element') != elements[cohort]:
            raise ValueError('native cohort label maps to a different frozen site element')
        elapsed = _number(row['elapsed_s'], 'native cohort elapsed seconds', minimum=0, maximum=horizon)
        key = (species, cohort, elapsed)
        if key in seen:
            raise ValueError('duplicate native site/species/cohort checkpoint')
        seen.add(key);groups.setdefault((species, cohort), []).append(elapsed)
        initial = _number(row['initial_atoms'], 'native cohort initial atoms', minimum=0)
        if initial == 0:
            raise ValueError('native cohort initial atoms must be positive')
        group = (species, cohort)
        if group in initial_by_group and initial_by_group[group] != initial:
            raise ValueError('native cohort initial physical stock changes between checkpoints')
        initial_by_group[group] = initial
        if expected_initial_atoms is not None and not math.isclose(initial,
                _number(expected_initial_atoms[species], 'configured initial cohort atoms', minimum=0),
                rel_tol=16*sys.float_info.epsilon, abs_tol=0):
            raise ValueError('native cohort stock differs from configured packet scale and actual initialized particle weight')
        diagnostic = {'seed_index': seed, 'site_id': ordered[cohort], 'cohort': cohort,
                      'species': species, 'logical_surface_element': elements[cohort],
                      'initial_epoch_tdb_s': start, 'elapsed_s': elapsed,
                      'ephemeris_time_tdb_s': start+elapsed, 'origin_depth_m': 0,
                      'frame': FRAME, 'initial_atoms': initial}
        for field in ('first_survival', 'restricted_first_residence_s', 'surface_fraction',
                      'gas_fraction', 'escape_fraction', 'numerical_loss_fraction'):
            diagnostic[field] = _number(row[field], 'native cohort '+field)
        diagnostics.append(diagnostic)
    expected = {(species, cohort) for species in ('Na', 'Ar') for cohort in range(len(ordered))}
    if set(groups) != expected:
        raise ValueError('native cohort history omits a required Na/Ar frozen site')
    checkpoint_times = None
    for group, times in groups.items():
        times.sort()
        if times[0] != 0 or times[-1] != horizon:
            raise ValueError('native cohort does not cover the exact initial/final elapsed interval')
        if checkpoint_times is not None and times != checkpoint_times:
            raise ValueError('native cohort checkpoint clocks differ between species/sites')
        checkpoint_times = times
    seen_hazards = set()
    for row in hazards:
        species = row['species'];cohort = _integer(row['cohort'], 'native hazard cohort')
        if (species, cohort) not in expected or _integer(row['element'], 'native hazard element') != elements[cohort]:
            raise ValueError('native hazard species/cohort/element differs from frozen site')
        begin = _number(row['begin_elapsed_s'], 'native hazard begin', minimum=0)
        duration = _number(row['dt_s'], 'native hazard duration', minimum=0)
        if duration == 0 or begin+duration > horizon+16*sys.float_info.epsilon*max(1, horizon):
            raise ValueError('native hazard interval lies outside its declared cohort horizon')
        key = (species, cohort, begin)
        if key in seen_hazards:
            raise ValueError('duplicate native site/species/hazard interval')
        seen_hazards.add(key)
        if _number(row['temperature_K'], 'native hazard kinetic temperature', minimum=0) == 0:
            raise ValueError('native hazard kinetic temperature must be positive')
        for field in ('total_release_hazard_s-1', 'thermal_release_hazard_s-1', 'photon_release_hazard_s-1'):
            _number(row[field], 'native '+field, minimum=0)
    return diagnostics


def validate_kinetic_run(record, base, sites, protocol, runtime, application, i19_manifest):
    """Require actual, qualified I19 bytes for one site-bound finest seed.

    A standalone PASS label or a different, scientifically plausible cohort
    file cannot replace the cohort files included in the formal I19 evidence
    list. Exact binary/build identity connects both diagnostic exporters to
    the same production implementation; per-input native time/seed/site
    controls are checked independently below.
    """
    paths = {name: _descriptor(record.get(name), base, 'kinetic '+name)
             for name in ('cohorts', 'hazards', 'application', 'receipt', 'manifest', 'budget', 'initial_state')}
    receipts = i19_manifest.get('native_evidence')
    source = i19_manifest.get('native_source_evidence')
    if not isinstance(receipts, list) or not receipts or not isinstance(source, list) or not source:
        raise MissingEvidence('formal I19 native evidence list required for site-cohort binding')
    qualified_files = set()
    for descriptor in receipts+source:
        path = Path(descriptor['path'])
        if not path.is_absolute():
            raise ValueError('I19 aggregate native evidence must use absolute paths')
        qualified_files.add((str(path.resolve()), descriptor['sha256']))
    if any((str(path), sha256(path)) not in qualified_files for path in paths.values()):
        raise ValueError('I20 cohort artifact was not part of the qualified same-campaign I19 evidence')
    from integration_common import parse_receipt
    from m2_illumination import setting
    from m2_protocol import cases
    seed_index = _integer(record['seed_index'], 'kinetic seed index')
    expected_case = next((case for case in cases('cluster')
                          if case['case'] == 'cohort' and case['dt_divisor'] == 4 and
                          case['weight_divisor'] == 4 and case['seed_index'] == seed_index), None)
    if expected_case is None:
        raise ValueError('kinetic seed is outside the fixed sixteen-seed finest design')
    applied = i19_manifest.get('applied_initial_state')
    qualification = i19_manifest.get('predecessor_and_build_evidence', {}).get('initial_state_qualification')
    if not isinstance(applied, dict) or not isinstance(qualification, dict):
        raise MissingEvidence('formal I19 applied native initial-state qualification unavailable to I20')
    from m2_initial_state_evidence import verify_case_initial_state,Unavailable as InitialStateUnavailable
    try:
        state = verify_case_initial_state(paths['application'].parent, qualification)
    except InitialStateUnavailable as error:
        raise MissingEvidence(str(error)) from error
    retained_state = applied.get('cases', {}).get(expected_case['key'])
    if (applied.get('native_case_count') != 288 or not retained_state or
            retained_state.get('path') != str(paths['initial_state']) or
            retained_state.get('sha256') != state['sha256'] or
            applied.get('qualification_sha256') != state['qualification_sha256']):
        raise ValueError('I20 cohort applied initial state differs from full qualified I19 native source evidence')
    if any(_integer(site['logical_surface_element'],'frozen cohort element') >=
           state['shape']['logical_surface_element_count'] for site in sites.values()):
        raise ValueError('frozen I20 cohort site is outside the qualified full initial-state grid')
    identity = _json(paths['manifest']);receipt = parse_receipt(paths['receipt'])
    expected_seed = protocol['seed_base']+protocol['seed_stride']*seed_index
    if (identity.get('test_id') != 'I19' or identity.get('execution_scope') != 'cluster-numerical-design' or
            identity.get('case') != expected_case or receipt.get('test_id') != 'I19' or
            receipt.get('verification_case') != 'cohort' or receipt.get('seed') != expected_seed or
            identity.get('random_seed') != expected_seed or
            identity.get('protocol', {}).get('sha256') != sha256(HERE/'M2_I19_I20_PROTOCOL.json')):
        raise ValueError('native I19 case/seed/protocol identity differs from fixed design')
    if (identity.get('mpi_ranks') != runtime.get('mpi_ranks') or
            receipt.get('mpi_ranks') != identity.get('mpi_ranks') or
            identity.get('openmp_threads') != 1 or runtime.get('openmp_threads') != 1 or
            identity.get('executable', {}).get('sha256') != runtime.get('executable', {}).get('sha256') or
            identity.get('inputs', {}).get(str(paths['application'])) != sha256(paths['application'])):
        raise ValueError('native I19 layout/executable/input differs from bound production run')
    build = identity.get('build_configuration_hashes')
    if not build or build != runtime.get('build_configuration_hashes') or not identity.get('source_revision'):
        raise ValueError('native I19 source/build identity is absent or differs from the thermal exporter')
    for key, expected in (('m2_verification_profile', 'cluster'), ('m2_verification_case', 'cohort'),
                          ('m2_verification_dt_divisor', '4'), ('m2_verification_weight_divisor', '4')):
        if setting(paths['application'], key) != expected:
            raise ValueError('native I19 input control differs: '+key)
    for key in ('simulation_start_utc', 'cold_trap_retention_timescale_s'):
        own, other = setting(paths['application'], key), setting(application, key)
        if own is None or other is None or (
                (own != other) if key == 'simulation_start_utc' else
                (_number(own, 'I19 retention timescale') != _number(other, 'I20 retention timescale'))):
            raise ValueError('native I19 and I20 epoch/retention configuration differ: '+key)
    horizon = _number(receipt['duration_s'], 'native I19 duration', minimum=0)
    if not math.isclose(horizon, protocol['physical_duration_s'], rel_tol=1e-14, abs_tol=0):
        raise ValueError('native I19 duration differs from fixed physical horizon')
    start = _number(receipt['initial_epoch_tdb_s'], 'native I19 initial TDB epoch')
    budgets = _csv(paths['budget'])
    clock_rows = {'Na': [], 'Ar': []}
    for row in budgets:
        if row['species'] not in clock_rows:
            raise ValueError('native I19 budget has a species outside Na/Ar')
        elapsed = _number(row['elapsed_s'], 'native budget elapsed time', minimum=0, maximum=horizon)
        if _number(row['species_elapsed_s'], 'native species elapsed time') != elapsed:
            raise ValueError('native cohort budget species clocks disagree')
        epoch = start+elapsed
        for key in ('epoch_tdb_s', 'actual_epoch_tdb_s'):
            if abs(_number(row[key], 'native budget '+key)-epoch) > 16*sys.float_info.epsilon*max(1, abs(epoch)):
                raise ValueError('native cohort absolute clock differs from declared start plus elapsed time')
        clock_rows[row['species']].append(elapsed)
    if (any(not times or min(times) != 0 or max(times) != horizon or len(times) != len(set(times))
            for times in clock_rows.values()) or sorted(clock_rows['Na']) != sorted(clock_rows['Ar'])):
        raise ValueError('native cohort absolute clock evidence omits or duplicates required endpoints/species')
    elements = setting(paths['application'], 'm2_verification_surface_elements')
    if elements is None:
        raise ValueError('native I19 input has no explicit cohort site elements')
    packet_scale = _number(setting(paths['application'], 'm2_verification_cohort_packets'),
                           'configured physical cohort packet scale', minimum=0)
    if packet_scale == 0:
        raise ValueError('configured physical cohort packet scale must be positive')
    # Production initializes fixed PHYSICAL stocks using its saved, pre-sweep
    # weights, then varies emitted packet weights for numerical convergence.
    # Do not multiply by the refined injection weight: that would change the
    # physical initial-value problem with the numerical weight divisor.
    initial_atoms = {species: packet_scale*_number(receipt['initialized_'+prefix+'_weight_atoms'],
                        'native initialized particle weight', minimum=0)
                     for species, prefix in (('Na', 'na'), ('Ar', 'ar'))}
    diagnostic = bind_cohort_sites(_csv(paths['cohorts']), _csv(paths['hazards'], allow_empty=True), sites,
        declared_elements=elements.split(','), horizon_s=horizon, initial_epoch_s=start, seed_index=seed_index,
        expected_initial_atoms=initial_atoms)
    return paths, diagnostic, {'seed_index': seed_index, 'source_revision': identity['source_revision'],
        'build_configuration_hashes': build, 'initial_epoch_tdb_s': start,
        'simulation_start_utc': setting(paths['application'], 'simulation_start_utc'),
        'horizon_s': horizon, 'mpi_ranks': identity['mpi_ranks'], 'openmp_threads': 1,
        'manifest': record['manifest']}


def publish_model_export(contract_path, *, thermal_samples, illumination_model,
                         numerical_qualification, kinetic_runs, native_evidence,
                         application, executable, source_revision, started_utc,
                         output_dir, application_template=None,
                         runtime_manifest=None, runtime_receipt=None,initial_state=None):
    """Publish a reproducible model-only bundle after real native execution.

    This helper does not approve the contract, invent a qualified temporal
    history, launch AMPS, or read reference temperatures. The owning runner
    must provide actual native artifacts; preflight/score verify them again.
    Descriptor paths are absolute so staging cannot change their meaning.
    """
    path = Path(contract_path).resolve()
    contract = _json(path)
    _review(contract, 'I20 contract')
    if application_template is None or sha256(application_template) != contract.get('application_template_sha256'):
        raise ValueError('production export requires the exact predeclared application template')
    if runtime_manifest is None or runtime_receipt is None or initial_state is None:
        raise MissingEvidence('native I20 runtime manifest/receipt required for site-kinetics binding')
    registry_path = _descriptor(contract.get('site_registry'), path.parent, 'native site registry')
    registry = _json(registry_path)
    _review(registry, 'native-site registry')
    sites = {str(row['site_id']): row for row in registry['sites']}
    thermal_samples = Path(thermal_samples).resolve()
    illumination_model = Path(illumination_model).resolve()
    rows = assemble_model_rows(_csv(thermal_samples), _csv(illumination_model),
                               contract['thermal_supports'], sites)
    output = Path(output_dir).resolve()
    if output == CANONICAL_DATA or CANONICAL_DATA in output.parents:
        raise ValueError('native model bundles cannot be written into validation data')
    output.mkdir(parents=True, exist_ok=True)
    target = output/'polar_model_rows.csv'
    if target.exists():
        raise ValueError('refusing to overwrite an earlier polar model export')
    fields = ('site_id', 'quantity', 'units', 'frame', 'time_support_id',
              'depth_m', 'model_value', 'temperature_kind', 'logical_surface_element')
    with target.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    def descriptor(value):
        resolved = Path(value).resolve()
        return {'path': str(resolved), 'sha256': sha256(resolved)}
    manifest = {'schema': 'srcMoon-I20-production-export-v1', 'status': 'PASS',
                'execution': 'executed', 'source_revision': source_revision,
                'started_utc': started_utc, 'comparison_contract_sha256': sha256(path),
                'application_template_sha256': sha256(application_template),
                'executable': descriptor(executable), 'input': descriptor(application),
                'runtime_manifest': descriptor(runtime_manifest), 'runtime_receipt': descriptor(runtime_receipt),
                'initial_state':descriptor(initial_state),
                'model_rows': descriptor(target), 'thermal_samples': descriptor(thermal_samples),
                'illumination_model': descriptor(illumination_model),
                'numerical_qualification': descriptor(numerical_qualification),
                'kinetic_runs': kinetic_runs, 'native_evidence': native_evidence,
                'scope': 'model-only artifact publication; formal qualification assigned by scorer'}
    manifest_path = output/'polar_model_manifest.json'
    _write_new_json(manifest_path, manifest)
    return manifest_path


def _write_new_json(path, value):
    text = json.dumps(value, indent=2, sort_keys=True) + '\n'
    if path.exists():
        if path.read_text(encoding='utf-8') != text:
            raise ValueError('refusing to overwrite earlier scientific/protocol evidence: ' + str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')


def _write_new_csv(path, rows):
    """Preserve the exact evaluated MODEL diagnostics without overwriting."""
    if not rows:
        raise ValueError('refusing to write an empty native diagnostic CSV')
    path = Path(path)
    if path.exists():
        raise ValueError('refusing to overwrite an earlier native diagnostic CSV')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def _predecessors(root):
    evidence = {}
    for test in PREDECESSORS:
        directory = Path(root) / test
        try:
            result = _json(directory/'result.json')
            manifest = _json(directory/'run_manifest.json')
        except FileNotFoundError as error:
            raise MissingEvidence(test + ' same-campaign linked evidence unavailable') from error
        scope = ' '.join(str(value) for value in (result.get('scope', ''), manifest.get('execution_scope', ''))).lower()
        if result.get('status') != 'PASS' or result.get('execution') == 'not-run' or (
                manifest.get('formal_test_status') == 'SKIPPED' or
                any(word in scope for word in ('supporting', 'preflight', 'development', 'source-only'))):
            raise MissingEvidence(test + ' is not formal linked PASS evidence')
        receipt = directory/'receipt.txt'
        if receipt.is_file():
            from integration_common import parse_receipt
            native = parse_receipt(receipt)
            if native.get('test_id') != test or native.get('mpi_ranks') != manifest.get('mpi_ranks'):
                raise ValueError(test + ' native receipt identity/layout differs from manifest')
        elif test == 'I19' and manifest.get('native_evidence'):
            for index, descriptor in enumerate(manifest['native_evidence']):
                _descriptor(descriptor, directory, 'I19 native evidence ' + str(index))
        else:
            raise MissingEvidence(test + ' actual native runtime receipt unavailable')
        if not manifest.get('executable', {}).get('sha256'):
            raise ValueError(test + ' native executable identity absent')
        evidence[test] = {'result_sha256': sha256(directory/'result.json'),
                          'manifest_sha256': sha256(directory/'run_manifest.json')}
    return evidence


def _qualified_initial_state(files,application=None,executable=None):
    """Authenticate the EXPLICIT I20 proof; never borrow a cohort epoch.

    Before a fresh run, existing real pilot artifacts authenticate the proof.
    After execution the actual input/binary are checked against it again. The
    contract may explicitly reference an I19 proof only when those physical
    settings and the actual I20 mode/epoch/state match; there is no fallback.
    """
    from m2_initial_state_evidence import load_qualification,Unavailable,SCHEMA
    path=Path(files['initial_state_qualification']);q=_json(path)
    if q.get('schema')!=SCHEMA:
        raise MissingEvidence('explicit I20 applied initial-state qualification unavailable')
    pilot_input=_descriptor(q.get('pilot_input'),path.parent,'I20 qualified pilot input')
    manifest_path=_descriptor(q.get('pilot_manifest'),path.parent,'I20 qualified pilot manifest')
    pilot_manifest=_json(manifest_path)
    pilot_executable=_descriptor(pilot_manifest.get('executable'),manifest_path.parent,'I20 qualified pilot executable')
    try:
        return load_qualification({'actual_epoch_thermal_qualification':q},
            pilot_input if application is None else application,
            pilot_executable if executable is None else executable,
            base=path.parent,allow_surface_pilot=True)
    except Unavailable as error:
        raise MissingEvidence(str(error)) from error


def preflight(contract_path, qualification_root, *, model_manifest_path=None, data_root=None,
              require_model=True):
    """Check metadata and reviewed model t0 state; no residual is scored.

    The frozen contract precedes model evaluation. The separate post-runtime
    manifest binds to its hash and contains the actual model/native file
    descriptors, avoiding a circular contract/model-manifest digest.
    Only the explicitly reviewed pilot's model INITIAL state is read; no
    observed row, annual-model comparison or December residual is consumed.
    """
    result = {'id': 'I20', 'status': 'SKIPPED', 'scope': 'preflight-only',
              'observational_scoring': False, 'gaps': [], 'errors': [], 'evidence': {}}
    try:
        if contract_path is None or not Path(contract_path).is_file():
            raise MissingEvidence('explicit reviewed/frozen I20 comparison contract unavailable')
        path, dispatch = resolve_contract(contract_path)
        result['evidence']['dispatch'] = dispatch
        contract = _json(path)
        if contract.get('schema') != SCHEMA:
            raise ValueError('unsupported I20 polar contract schema')
        approved = _review(contract, 'I20 contract')
        if contract.get('frame') != FRAME or contract.get('coordinate_convention') != 'planetocentric-east-positive':
            raise ValueError('I20 contract must preserve the native ME/DE421 coordinates')
        from integration_common import dataset_package_preflight
        for dataset in ('D02', 'D03'):
            # The established package checker verifies native raw/processed
            # hashes and READY status, not merely a directory or CSV header.
            if data_root is not None:
                # The existing checker uses an environment-independent canonical
                # authority. A different data root is unsupported rather than
                # bypassing that authority with caller-supplied READY flags.
                if Path(data_root).resolve() != CANONICAL_DATA:
                    raise MissingEvidence('noncanonical I20 data authority requires a separately reviewed package adapter')
            ready, reason, proof = dataset_package_preflight(dataset)
            result['evidence'][dataset] = proof
            if not ready:
                raise MissingEvidence(reason)
        result['evidence']['predecessors'] = _predecessors(qualification_root)
        files = {}
        for key in ('site_registry', 'candidates', 'collocation_contract', 'reference_rows',
                    'operator_qualification', 'thermal_queries','initial_state_qualification'):
            files[key] = _descriptor(contract.get(key), path.parent, key)
        registry = _json(files['site_registry'])
        _review(registry, 'native-site registry')
        if registry.get('schema') != 'srcMoon-I20-native-site-registry-v1':
            raise ValueError('unsupported native-site registry schema')
        collocation = _json(files['collocation_contract'])
        _review(collocation, 'collocation contract')
        if registry.get('candidates', {}).get('sha256') != sha256(files['candidates']) or (
                registry.get('collocation_contract', {}).get('sha256') != sha256(files['collocation_contract'])):
            raise ValueError('native-site freeze refers to different candidates/collocation')
        selected = select_sites(_csv(files['candidates']), maximum_distance_m=collocation['maximum_distance_m'])
        if registry.get('sites') != selected['sites'] or registry.get('strata') != selected['strata'] or (
                registry.get('shared_model_support') != selected['shared_model_support']):
            raise ValueError('frozen native sites differ from the approved preselection rule')
        sites = {str(row['site_id']): row for row in selected['sites']}
        if not sites or len(sites) != len(selected['sites']):
            raise ValueError('native-site registry is empty or has duplicate site_id')
        operator = _json(files['operator_qualification'])
        if operator.get('status') != 'PASS' or operator.get('execution') != 'executed' or (
                operator.get('site_registry_sha256') != sha256(files['site_registry']) or
                operator.get('collocation_contract_sha256') != sha256(files['collocation_contract']) or
                operator.get('scorer_sha256') != sha256(Path(__file__)) or
                not operator.get('evidence_files')):
            raise MissingEvidence('independent polar observation/temporal operator qualification unavailable')
        for index, descriptor in enumerate(operator['evidence_files']):
            _descriptor(descriptor, files['operator_qualification'].parent, 'operator proof ' + str(index))
        # No numerical polar tolerance is made up. Review must supply these
        # criteria before a model row is consumed, distinct from native sigma.
        if not isinstance(contract.get('criteria', {}), dict):
            raise ValueError('polar criteria must be an object keyed by quantity')
        for quantity in QUANTITIES:
            criterion = contract.get('criteria', {}).get(quantity)
            if not criterion or 'maximum_absolute_error' not in criterion:
                raise MissingEvidence('reviewed polar criterion unavailable: ' + quantity)
            _number(criterion['maximum_absolute_error'], 'reviewed polar limit', minimum=0)
        if not re.fullmatch('[0-9a-f]{64}', str(contract.get('application_template_sha256', ''))):
            raise ValueError('I20 contract must bind the existing reviewed application-template SHA-256')
        if not isinstance(contract.get('thermal_supports', []), list):
            raise ValueError('polar thermal_supports must be an explicit array')
        for support in contract.get('thermal_supports', []):
            if not isinstance(support, dict):
                raise ValueError('thermal support must be an explicit object')
            if support.get('coverage_qualified') is not True:
                raise MissingEvidence('native polar thermal temporal support is unqualified')
            _descriptor(support.get('sampling_evidence'), path.parent, 'thermal sampling proof')
        if not contract.get('thermal_supports'):
            raise MissingEvidence('explicit polar thermal operators unavailable')
        if not require_model:
            result['evidence']['initial_state_qualification']=_qualified_initial_state(files)
            result.update(status='PASS', reason='reference/geometry prerequisites available; native model not launched',
                          contract=contract, sites=sites,
                          files={key: str(value) for key, value in files.items()},
                          contract_path=str(path))
            result['evidence']['comparison_contract_sha256'] = sha256(path)
            return result
        if model_manifest_path is None or not Path(model_manifest_path).is_file():
            raise MissingEvidence('actual production I20 export manifest unavailable')
        model_manifest_path = Path(model_manifest_path).resolve()
        native = _json(model_manifest_path)
        if native.get('schema') != 'srcMoon-I20-production-export-v1' or native.get('status') != 'PASS' or (
                native.get('execution') != 'executed' or native.get('comparison_contract_sha256') != sha256(path)):
            raise ValueError('production polar export does not match its frozen comparison contract')
        if _timestamp(native['started_utc']) <= approved:
            raise ValueError('polar contract was not frozen before production evaluation')
        if not native.get('source_revision') or not native.get('executable', {}).get('sha256') or (
                native.get('application_template_sha256') != contract['application_template_sha256']):
            raise ValueError('polar model source/build/input provenance differs from frozen configuration')
        files['input'] = _descriptor(native.get('input'), model_manifest_path.parent, 'actual native input')
        _descriptor(native.get('executable'), model_manifest_path.parent, 'actual native executable')
        for key in ('model_rows', 'thermal_samples', 'illumination_model', 'numerical_qualification',
                    'runtime_manifest', 'runtime_receipt','initial_state'):
            files[key] = _descriptor(native.get(key), model_manifest_path.parent, key)
        numerical = _json(files['numerical_qualification'])
        if numerical.get('status') != 'PASS' or numerical.get('scope') != 'fixed full numerical design; not observational agreement':
            raise MissingEvidence('full fixed-design numerical qualification unavailable')
        if not native.get('native_evidence') or not native.get('kinetic_runs'):
            raise MissingEvidence('actual production native/kinetic evidence unavailable')
        if not isinstance(native['native_evidence'], list) or not isinstance(native['kinetic_runs'], list):
            raise ValueError('native/kinetic evidence must be explicit arrays')
        for index, descriptor in enumerate(native['native_evidence']):
            _descriptor(descriptor, model_manifest_path.parent, 'production runtime proof ' + str(index))
        result['evidence']['initial_state_qualification']=_qualified_initial_state(files,
            files['input'],_descriptor(native['executable'],model_manifest_path.parent,'actual I20 executable'))
        result.update(status='PASS', reason='all declared metadata/prerequisites available; no residual scored',
                      contract=contract, native_manifest=native, sites=sites,
                      files={key: str(value) for key, value in files.items()},
                      contract_path=str(path), model_manifest_path=str(model_manifest_path))
        result['evidence']['comparison_contract_sha256'] = sha256(path)
    except MissingEvidence as error:
        result['gaps'].append(str(error));result['reason'] = str(error)
    except (OSError, ValueError, TypeError, KeyError) as error:
        result.update(status='ERROR', reason=str(error));result['errors'].append(str(error))
    return result


def prepare_input(template, destination, campaign_root):
    """Stage geometry-only input after metadata/applied-state authentication.

    The fresh model outputs do not exist yet and are never required by this
    hook. Missing reference/predecessor evidence raises the shared typed
    MissingPrerequisite exception so the existing global runner writes its
    ordinary honest SKIPPED receipt. Exact runtime input and geometry hashes
    are bound in a separate pre-launch record, not written into the frozen
    reference contract and never used as observational truth.
    """
    from integration_common import MissingPrerequisite, replace_moon_setting
    import shutil
    template, destination = Path(template).resolve(), Path(destination).resolve()
    from m2_illumination import setting
    # Only absence is SKIPPED. Duplicate or malformed scalar syntax remains
    # ValueError/ERROR and cannot be hidden by an optional-input fallback.
    name = setting(template, 'm2_formal_contract_file')
    if not name or name.lower() == 'none':
        raise MissingPrerequisite('I20 reviewed polar comparison contract unavailable')
    contract_path = Path(name)
    if not contract_path.is_absolute():
        contract_path = template.parent/contract_path
    gate = preflight(contract_path, campaign_root, require_model=False)
    if gate['status'] == 'SKIPPED':
        raise MissingPrerequisite(gate['reason'])
    if gate['status'] != 'PASS':
        raise ValueError(gate['reason'])
    if sha256(template) != gate['contract']['application_template_sha256']:
        raise ValueError('I20 application template differs from the frozen contract')
    from m2_initial_state import physics_identity_sha256
    initial=gate.get('evidence',{}).get('initial_state_qualification')
    if not isinstance(initial,dict):
        raise MissingPrerequisite('explicit I20 applied initial-state proof unavailable')
    if physics_identity_sha256(template)!=initial.get('physics_identity_sha256'):
        raise ValueError('I20 template physical initial-value problem differs from its explicit initial-state proof')
    queries = Path(gate['files']['thermal_queries'])
    destination.parent.mkdir(parents=True, exist_ok=True)
    staged_queries = destination.parent/'polar_thermal_queries.csv'
    if staged_queries.exists() and sha256(staged_queries) != sha256(queries):
        raise ValueError('refusing to overwrite changed I20 geometry queries')
    if not staged_queries.exists():
        shutil.copy2(queries, staged_queries)
    # Formal I20 is not the bounded I18 distribution fixture. Preserve every
    # physical AMR/source/thermal setting in the frozen template rather than
    # inheriting prepare_na_ar_input's intentional volume coarsening/Ar switch.
    # Only path ownership and geometry-only diagnostics are staged here.
    if setting(template, 'surface_temperature_mode') != 'thermal-model':
        raise MissingPrerequisite('I20 template must explicitly select thermal-model; the formal runner does not change physics')
    text = template.read_text(encoding='utf-8')
    for key in ('spice_path', 'lola_image_file', 'lola_label_file',
                'diviner_temperature_table_file', 'diviner_eclipse_history_file',
                'diviner_eclipse_selection_contract_file', 'diviner_eclipse_freeze_file'):
        original = setting(template, key)
        if original is None or original.lower() == 'none':
            continue
        absolute = Path(original)
        if not absolute.is_absolute():
            absolute = template.parent/absolute
        text = replace_moon_setting(text, key, str(absolute.resolve()), template)
    for key, name in (('surface_cea_file', 'lola_surface.cea'),
                      ('surface_tecplot_file', 'lola_surface.dat')):
        text = replace_moon_setting(text, key, str((destination.parent/name).resolve()), template)
    for key, value in (('m2_formal_contract_file', str(Path(gate['contract_path']).resolve())),
                       ('m2_polar_queries_file', str(staged_queries))):
        # Optional completion keys may be absent in a legacy template. Add
        # them only inside the unique application block; mandatory existing
        # physics settings retain the normal shared replacement routine.
        if re.search(r'(?m)^\s*'+re.escape(key)+r'\s*=', text):
            text = replace_moon_setting(text, key, value, template)
        else:
            endings = list(re.finditer(r'(?m)^\s*#end application[^\n]*', text))
            if len(endings) != 1:
                raise ValueError('I20 input needs one application block end')
            at = endings[0].start();text = text[:at]+key+' = '+value+'\n'+text[at:]
    destination.write_text(text, encoding='utf-8')
    binding = {'schema': 'srcMoon-I20-prelaunch-input-v1', 'observational_scoring': False,
               'comparison_contract_sha256': sha256(Path(gate['contract_path'])),
               'application_template': {'path': str(template), 'sha256': sha256(template)},
               'actual_input': {'path': str(destination), 'sha256': sha256(destination)},
               'thermal_queries': {'path': str(staged_queries), 'sha256': sha256(staged_queries)},
               'dispatch': gate['evidence'].get('dispatch')}
    _write_new_json(destination.parent/'polar_input_binding.json', binding)
    _write_new_json(destination.parent/'polar_reference_manifest.json', gate)
    return gate


def finish_run(work, receipt):
    """Join a fresh native I20 run to existing same-campaign I15/I19 outputs.

    The shared integration runner calls this only after the ordinary native
    process has returned. All descriptors are generated from actual files;
    the original frozen contract contains no future output checksum. This is
    the normal fresh-run route, not an import-only or hand-assembled shortcut.
    """
    work = Path(work).resolve()
    metrics = {'native_thermal_export_executed': True}
    try:
        proof = _json(work/'polar_reference_manifest.json')
        binding = _json(work/'polar_input_binding.json')
        run = _json(work/'run_manifest.json')
        if proof.get('status') != 'PASS' or proof.get('scope') != 'preflight-only':
            raise ValueError('native polar run lacks its successful metadata preflight')
        contract_path = Path(proof['contract_path'])
        if sha256(contract_path) != binding['comparison_contract_sha256']:
            raise ValueError('polar contract changed after pre-launch binding')
        for key in ('input_document', 'selected_document'):
            if key in (binding.get('dispatch') or {}):
                _descriptor(binding['dispatch'][key], work, 'reviewed M2 dispatch '+key)
        application = _descriptor(binding.get('actual_input'), work, 'bound native polar input')
        template = _descriptor(binding.get('application_template'), work, 'bound application template')
        _descriptor(binding.get('thermal_queries'), work, 'bound geometry-only queries')
        if run.get('test_id') != 'I20' or receipt.get('test_id') != 'I20' or (
                receipt.get('mpi_ranks') != run.get('mpi_ranks') or
                run.get('inputs', {}).get(str(application)) != sha256(application)):
            raise ValueError('fresh native I20 receipt/input/layout disagrees with run manifest')
        native_csv = work/'receipt.txt.I20_thermal_model.csv'
        if not native_csv.is_file() or _integer(receipt.get('thermal_model_rows', 0), 'native row count') < 1:
            raise ValueError('fresh I20 thermal model export is absent or empty')
        executable = _descriptor(run.get('executable'), work, 'fresh I20 executable')
        campaign_root = work.parent
        i19 = _json(campaign_root/'I19/run_manifest.json')
        directory = Path(i19.get('native_campaign_directory', ''))
        if not directory.is_absolute():
            raise ValueError('I19 formal manifest lacks an absolute native campaign directory')
        native_runs = []
        def descriptor(file):
            file = Path(file).resolve()
            if not file.is_file():
                raise MissingEvidence('actual native polar predecessor artifact unavailable: ' + str(file))
            return {'path': str(file), 'sha256': sha256(file)}
        for seed in range(16):
            case = directory/'cases'/f'cohort-d4-w4-s{seed:02d}'
            native_runs.append(dict(seed_index=seed, **{
                name: descriptor(case/filename) for name, filename in
                (('cohorts', 'm2_cohorts.csv'), ('hazards', 'm2_hazards.csv'),
                 ('application', 'application.in'), ('receipt', 'receipt.txt'),
                 ('manifest', 'run_manifest.json'), ('budget', 'm2_budget.csv'),
                 ('initial_state', 'receipt.txt.I19_initial_thermal_state.csv'))}))
        export = publish_model_export(contract_path,
            thermal_samples=native_csv,
            illumination_model=campaign_root/'I15/illumination_model.csv',
            numerical_qualification=campaign_root/'I19/independent_numerical_score.json',
            kinetic_runs=native_runs,
            native_evidence=[descriptor(work/'receipt.txt'), descriptor(work/'run_manifest.json')],
            application=application, application_template=template, executable=executable,
            runtime_manifest=work/'run_manifest.json', runtime_receipt=work/'receipt.txt',
            initial_state=work/'receipt.txt.I20_initial_thermal_state.csv',
            source_revision=run['source_revision'], started_utc=run['started_utc'],
            output_dir=work/'polar-comparison')
        result = score(contract_path, campaign_root, work/'polar-comparison', model_manifest_path=export)
        metrics.update(result)
        metrics['artifacts'] = result.get('artifacts', {})
        return result['status'], metrics, result['reason']
    except MissingEvidence as error:
        return 'SKIPPED', metrics, str(error)
    except (OSError, ValueError, TypeError, KeyError) as error:
        return 'ERROR', metrics, str(error)


def score(contract_path, qualification_root, output_dir, *, model_manifest_path=None,
          render=True, data_root=None):
    """Evaluate formal I20 only after the unchanged metadata/runtime gates.

    The enclosing global/native runner owns executable launch, resource budget
    and full run_manifest.json. This function preserves a separate preflight
    and row-level score, so a missing/corrupt contract cannot become an empty
    green plot. Actual output rows are never subset-selected by their errors.
    """
    output = Path(output_dir).resolve()
    if output == CANONICAL_DATA or CANONICAL_DATA in output.parents:
        raise ValueError('I20 artifacts must not be written into canonical validation data')
    output.mkdir(parents=True, exist_ok=True)
    gate = preflight(contract_path, qualification_root, model_manifest_path=model_manifest_path, data_root=data_root)
    _write_new_json(output/'polar_preflight.json', gate)
    result = {'id': 'I20', 'status': gate['status'], 'scope': 'formal-polar-qualification',
              'reason': gate.get('reason', ''), 'observational_scoring': False,
              'runner': str(Path(__file__)), 'log': str(output/'polar_preflight.json'),
              'artifacts': {'preflight': str(output/'polar_preflight.json')}}
    if gate['status'] != 'PASS':
        _write_new_json(output/'polar_result.json', result)
        return result
    try:
        contract, native = gate['contract'], gate['native_manifest']
        files = {key: Path(value) for key, value in gate['files'].items()}
        rows = _csv(files['model_rows'])
        samples = _csv(files['thermal_samples'])
        illumination_rows = _csv(files['illumination_model'])
        from integration_common import parse_receipt
        from m2_illumination import setting
        runtime = _json(files['runtime_manifest'])
        thermal_receipt = parse_receipt(files['runtime_receipt'])
        if (runtime.get('test_id') != 'I20' or thermal_receipt.get('test_id') != 'I20' or
                runtime.get('mpi_ranks') != thermal_receipt.get('mpi_ranks') or
                runtime.get('executable', {}).get('sha256') != native['executable']['sha256'] or
                runtime.get('inputs', {}).get(str(files['input'])) != sha256(files['input']) or
                runtime.get('source_revision') != native['source_revision'] or
                runtime.get('started_utc') != native['started_utc']):
            raise ValueError('native I20 runtime receipt/manifest differs from bound model export')
        export_schema = thermal_receipt.get('surface_kinetics_export_schema')
        if export_schema is None:
            raise MissingEvidence('native I20 surface-kinetics exporter capability unavailable')
        if _integer(thermal_receipt.get('thermal_model_rows'), 'native thermal row count', minimum=1) != len(samples):
            raise ValueError('native I20 receipt count differs from actual thermal/kinetics export')
        initial=gate.get('evidence',{}).get('initial_state_qualification')
        if not isinstance(initial,dict):
            raise MissingEvidence('explicit I20 applied initial-state proof unavailable; no cohort-phase fallback')
        if _timestamp(native['started_utc'])<=_timestamp(initial['approved_utc']):
            raise ValueError('I20 native history began before its applied initial-state proof was frozen')
        from m2_initial_state_evidence import verify_surface_initial_state,Unavailable as InitialStateUnavailable
        try:
            applied=verify_surface_initial_state(thermal_receipt,runtime,files['input'],files['initial_state'],initial)
        except InitialStateUnavailable as error:
            raise MissingEvidence(str(error)) from error
        _write_new_json(output/'initial_state_binding.json',applied)
        result['artifacts']['initial_state_binding']=str(output/'initial_state_binding.json')
        protocol = _json(HERE/'M2_I19_I20_PROTOCOL.json')
        timescale = _number(setting(files['input'], 'cold_trap_retention_timescale_s'),
                            'explicit native retention timescale', minimum=0)
        if timescale != protocol['source_and_inventory_controls']['retention_reporting_timescale_s']:
            raise ValueError('native retention-reporting timescale differs from reviewed M2 numerical protocol')
        # Thermal summaries must be reproducible from the actual native
        # per-epoch export; a caller cannot inject a convenient annual value.
        actual_model = assemble_model_rows(samples, illumination_rows,
                                           contract['thermal_supports'], gate['sites'])
        actual_thermal = [row for row in actual_model if QUANTITIES[row['quantity']][0] == 'K']
        supplied_thermal = [row for row in rows if row['quantity'] in (
            'annual_maximum_surface_temperature_K', 'annual_average_temperature_at_2cm_K')]
        def key(row): return row['site_id'], row['quantity']
        by_key = {key(row): row for row in supplied_thermal}
        if set(by_key) != {key(row) for row in actual_thermal}:
            raise ValueError('polar model thermal rows differ from actual native support')
        for row in actual_thermal:
            other = by_key[key(row)]
            if _number(other['model_value'], 'exported thermal value') != row['model_value'] or (
                    other['time_support_id'] != row['time_support_id'] or
                    other.get('temperature_kind') != row['temperature_kind']):
                raise ValueError('polar annual quantity is not the declared native-sample aggregation')
        all_keys = {(row['site_id'], row['quantity']): row for row in rows}
        if set(all_keys) != {(row['site_id'], row['quantity']) for row in actual_model}:
            raise ValueError('polar environment export changed native site/quantity membership')
        for row in actual_model:
            other = all_keys[(row['site_id'], row['quantity'])]
            if _number(other['model_value'], 'exported model value') != row['model_value'] or (
                    other['time_support_id'] != row['time_support_id']):
                raise ValueError('polar model row differs from the actual native environment export')
        # Existing numerical qualification covers the complete fixed design.
        # Recheck conditional law/warming and tag accounting for every supplied
        # finest-seed cohort using existing independent production diagnostics.
        retention = surface_retention_samples(samples, contract['thermal_supports'], gate['sites'],
            illumination_rows, retention_timescale_s=timescale, native_export_schema=export_schema)
        _write_new_csv(output/'surface_retention.csv', retention['rows'])
        _write_new_json(output/'surface_retention.json', retention)
        result['artifacts'].update(surface_retention_csv=str(output/'surface_retention.csv'),
                                   surface_retention_json=str(output/'surface_retention.json'))
        if retention['status'] != 'PASS':
            result.update(status='FAIL', reason='native surface-retention invariant failed; no observational residuals scored',
                          failures=retention['failures'])
            _write_new_json(output/'polar_result.json', result)
            return result
        from m2_protocol import evaluate_cohorts
        runs = native['kinetic_runs']
        if sorted(_integer(row['seed_index'], 'kinetic seed') for row in runs) != list(range(16)):
            raise ValueError('I20 kinetic evidence must preserve all sixteen predeclared finest-grid seeds')
        i19_path = Path(qualification_root)/'I19/run_manifest.json'
        if sha256(i19_path) != gate['evidence']['predecessors']['I19']['manifest_sha256']:
            raise ValueError('formal I19 manifest changed after polar preflight')
        i19_manifest = _json(i19_path)
        kinetics, site_kinetics, kinetic_identities = [], [], []
        for row in runs:
            paths, diagnostics, identity = validate_kinetic_run(row, Path(gate['model_manifest_path']).parent,
                gate['sites'], protocol, runtime, files['input'], i19_manifest)
            if kinetic_identities and identity['initial_epoch_tdb_s'] != kinetic_identities[0]['initial_epoch_tdb_s']:
                raise ValueError('native cohort initial absolute epochs differ between finest seeds')
            kinetic_identities.append(identity);site_kinetics.extend(diagnostics)
            checked = evaluate_cohorts(paths['cohorts'], paths['hazards'], paths['application'])
            checked.update(seed_index=int(row['seed_index']))
            kinetics.append(checked)
        _write_new_csv(output/'site_kinetics.csv', site_kinetics)
        _write_new_json(output/'conditional_kinetics.json', {
            'scope': 'site-bound conditional kinetics; independent annual thermal validation separate',
            'results': kinetics, 'native_run_identities': kinetic_identities,
            'permanent_retention_established': False,
            'shared_model_support': len({str(site['logical_surface_element']) for site in gate['sites'].values()}) < len(gate['sites'])})
        result['artifacts'].update(site_kinetics_csv=str(output/'site_kinetics.csv'),
                                   conditional_kinetics=str(output/'conditional_kinetics.json'))
        if any(row['status'] != 'PASS' for row in kinetics):
            # A resolved production-kinetic failure is the first failing
            # invariant. Preserve it without exposing observational residuals
            # that could distract from or substitute for the earlier layer.
            result.update(status='FAIL', reason='conditional production kinetic/retention identity failed; no observational residuals scored',
                          failures=[row for row in kinetics if row['status'] != 'PASS'])
            result['artifacts']['conditional_kinetics'] = str(output/'conditional_kinetics.json')
            _write_new_json(output/'polar_result.json', result)
            return result
        comparison = compare_rows(_csv(files['reference_rows']), rows, gate['sites'], contract['criteria'])
        result.update(status=('FAIL' if any(row['status'] != 'PASS' for row in kinetics) else comparison['status']),
                      reason='all declared independent polar and conditional kinetic criteria evaluated',
                      observational_scoring=True, evaluated_rows=comparison['evaluated_rows'],
                      failures=comparison['failures'], shared_model_support=(
                          len({str(row['logical_surface_element']) for row in gate['sites'].values()}) < len(gate['sites'])))
        _write_new_json(output/'polar_rows.json', {'rows': comparison['diagnostics']})
        csv_path = output/'comparison.csv'
        if csv_path.exists():
            raise ValueError('refusing to overwrite an earlier evaluated polar comparison')
        with csv_path.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=COMPARISON_FIELDS, lineterminator='\n')
            writer.writeheader();writer.writerows(comparison['comparisons'])
        result['artifacts'].update(comparison_csv=str(csv_path), diagnostics=str(output/'polar_rows.json'),
                                   conditional_kinetics=str(output/'conditional_kinetics.json'))
        if render:
            from common.write_comparison_artifacts import render as render_rows
            render_rows(csv_path, output/'comparison.png', output/'comparison.eps', 'I20 independent polar comparison')
            result['artifacts'].update(comparison_png=str(output/'comparison.png'), comparison_eps=str(output/'comparison.eps'))
    except MissingEvidence as error:
        result.update(status='SKIPPED', reason=str(error))
    except (OSError, ValueError, TypeError, KeyError, ImportError) as error:
        result.update(status='ERROR', reason=str(error))
    _write_new_json(output/'polar_result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path)
    parser.add_argument('--qualification-root', type=Path)
    parser.add_argument('--model-manifest', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--prepare-sites', type=Path, metavar='CANDIDATES_CSV')
    parser.add_argument('--collocation-contract', type=Path)
    args = parser.parse_args()
    if args.prepare_sites:
        if not args.collocation_contract:
            parser.error('--prepare-sites requires an explicit reviewed --collocation-contract')
        result = prepare_site_registry(args.prepare_sites, args.collocation_contract,
                                       args.output_dir/'site_registry_candidate.json')
        print(json.dumps({'status': result['status'], 'selected': len(result['sites']),
                          'observational_scoring': False}))
        return 0
    if args.qualification_root is None:
        parser.error('--qualification-root is required for formal comparison')
    result = score(args.contract, args.qualification_root, args.output_dir,
                   model_manifest_path=args.model_manifest)
    print(json.dumps({'status': result['status'], 'reason': result['reason'],
                      'result': str(args.output_dir/'polar_result.json')}))
    return {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}[result['status']]


if __name__ == '__main__':
    raise SystemExit(main())
