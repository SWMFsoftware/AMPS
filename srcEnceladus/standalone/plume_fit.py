#!/usr/bin/env python3
"""Standalone CLI. Run `python3 plume_fit.py --help`; no AMPS build needed."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from fitting import SOURCE_CASES, fit_density, load_config


def json_value(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    raise TypeError(type(x).__name__)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def load_density(path, metadata_path):
    meta = json.loads(Path(metadata_path).read_text())
    required = ('observable', 'frame', 'position_units', 'density_units', 'uncertainty_units',
                'uncertainty_model', 'calibration_reference', 'sticking_treatment',
                'selection_reference', 'data_kind', 'time_sampling')
    if any(not meta.get(k) for k in required):
        raise ValueError('Missing explicit observation metadata: '+', '.join(required))
    expected = dict(observable='number_density', frame='IAU_ENCELADUS', position_units='m',
                    density_units='m^-3', uncertainty_units='m^-3',
                    uncertainty_model='independent_gaussian', time_sampling='instantaneous')
    if any(meta[k] != v for k, v in expected.items()):
        raise ValueError('Unsupported units/observable/error model: raw counts, correlated errors and exposure averaging require dedicated operators')
    if meta['data_kind'] not in ('synthetic', 'calibrated_observation'):
        raise ValueError('Explicit synthetic or calibrated_observation data_kind required')
    with Path(path).open(newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError('No observations')
    x = np.array([[float(r[k]) for k in ('x_m', 'y_m', 'z_m')] for r in rows])
    y = np.array([float(r['density_m3']) for r in rows])
    sd = np.array([float(r['sigma_m3']) for r in rows])
    return x, y, sd, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, help='Explicit spherical-model JSON')
    parser.add_argument(
        '--source-case', choices=SOURCE_CASES,
        help='Override source_case in the model JSON')
    parser.add_argument(
        '--jet-catalog',
        help=('External Porco/Portyankina catalog JSON; alternatively set '
              'jet_catalog_file or ENCELADUS_JET_CATALOG'))
    parser.add_argument('--observations', required=True, help='Qualified density CSV; never L1A counts')
    parser.add_argument('--metadata', required=True, help='Observation provenance/units/response JSON')
    parser.add_argument('--output', required=True, help='New result JSON; existing files are refused')
    args = parser.parse_args()
    try:
        config = load_config(args.config, args.source_case, args.jet_catalog)
        x, y, sigma, metadata = load_density(args.observations, args.metadata)
        result = fit_density(config, x, y, sigma)
        result['observation_metadata'] = metadata
        provenance_files = [
            args.config, args.observations, args.metadata, __file__,
            Path(__file__).with_name('plume_model.py'),
            Path(__file__).with_name('fitting.py')]
        if config.get('_stripe_geometry_path'):
            provenance_files.append(config['_stripe_geometry_path'])
        if config.get('_jet_catalog_path'):
            provenance_files.append(config['_jet_catalog_path'])
        result['provenance'] = {
            str(path): digest(path) for path in provenance_files}
        result['qualification'] = (
            'stage-0 density-profile fit; not validation of gravity, collisions, '
            'chemistry, time variability, or the digitized observation operator')
        with Path(args.output).open('x') as f:
            json.dump(result, f, default=json_value, indent=2, allow_nan=False)
            f.write('\n')
    except (ValueError, KeyError, OSError, RuntimeError, ArithmeticError) as exc:
        parser.exit(2, 'ERROR: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
