#!/usr/bin/python3
"""Prepare the CME regression fit and verify the native shock trace."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

from pack_observations import interpolate_hci, pointing

HERE = Path(__file__).resolve().parent
AMPS = HERE.parents[2]
sys.path.insert(0, str(AMPS / 'src/models/sep_corona_swcme/forecast'))

import replay
from manufactured_event import AMBIENT, REFERENCE, TRUTH, synthetic, timestamp

STEP_S = 1200
END_S = 18000
HANDOFF_S = 10800
WIND_M_S = 400_000
DRAG_M_INV = 5e-12
HANDOFF_MIN_M = 20 * replay.SOLAR_RADIUS_M
EPOCHS = np.arange(0, END_S + STEP_S, STEP_S)
HEADER = ('elapsed_s', 'ray_id', 'x_m', 'y_m', 'z_m', 'n_x', 'n_y', 'n_z',
          'normal_speed_m_s', 'status', 'fast_mach', 'density_compression',
          'theta_bn_rad', 'max_rh_residual')
POSITIONS = ('x_m', 'y_m', 'z_m')
NORMALS = ('n_x', 'n_y', 'n_z')
SOLVED = 'solved-fast-shock'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fixture(path):
    # Generate an editable example outside the regression directory. The
    # legacy Table reads the frozen observations.json checked into the repo.
    event = synthetic()
    event.pop('rays_hci')
    event.pop('observers_hci_m')
    event['propagation']['drag_m_inv'] = DRAG_M_INV
    event['propagation']['horizon_s'] = END_S - HANDOFF_S
    path.parent.mkdir(parents=True, exist_ok=True)
    replay.write_json(path, event)


def ambient_path(event, data, issue, selected):
    # Use one coherent upstream realization. Automatic selection is safe
    # only when precisely one member was available at the issue time.
    members = event.get('ambient_members', [])
    eligible = [i for i, member in enumerate(members) if replay.utc(member['available_utc']) <= issue]
    if selected is None:
        require(eligible, 'a native trace requires a qualified upstream member available at issue time')
        require(len(eligible) == 1, 'supply --ambient-member when more than one upstream member is available')
        selected = eligible[0]
    require(selected in eligible, 'selected ambient member is missing or unavailable at issue time')
    return (data.parent / members[selected]['event_file']).resolve(strict=True)


def prepare(args):
    raw = args.data.read_bytes()
    event = json.loads(raw)
    require(event['schema'] == 'shock-front-replay-v1' and event['frame'] == 'HCI',
            'expected shock-front-replay-v1 data in HCI')
    for key in ('light_time_policy', 'provenance'):
        require(isinstance(event[key], str) and event[key].strip(), f'document {key}')
    reference = replay.utc(event['reference_utc'])
    require(event['reference_utc'] == reference.strftime('%Y-%m-%dT%H:%M:%SZ'),
            'the native event reference needs whole-second UTC in YYYY-MM-DDTHH:MM:SSZ form')
    # Check the geometry before plasma admission. An out-of-domain observed
    # front should report the geometric failure even with no upstream member.
    fit, _, epochs, diagnostics = replay.fit_tracks(event, args.issue, set(args.withhold_view))
    ambient = ambient_path(event, args.data, args.issue, args.ambient_member)
    first, last = epochs[0], epochs[-1]
    start = first if args.start_elapsed_s is None else replay.finite(args.start_elapsed_s, 'trace start')
    require(start >= first, 'the trace may start at or after the first usable front epoch')
    end = start + END_S
    outer = None
    if end > last:
        require(event.get('propagation'), 'a trace beyond the fit requires explicit propagation')
        outer = replay.resolve_propagation(event['propagation'])
        initial_radius, _ = replay.kinematics(fit, last / replay.HOUR_S)
        require(end - last <= outer[2], 'five-hour trace exceeds the declared propagation horizon')
        require(initial_radius * replay.SOLAR_RADIUS_M >= outer[3],
                'final fitted apex is below the minimum drag handoff radius')

    axis = replay.direction(*fit[3:5])
    rays = ([replay.vector(ray, 'sample ray', unit=True) for ray in event['rays_hci']]
            if 'rays_hci' in event else replay.default_rays(fit))
    require(len(rays) > 0, 'the trace needs at least one HCI ray')
    # Fitter units are solar radii and hours. Every value passed to the C++
    # evaluator is in SI with HCI Cartesian directions. Coefficients are
    # sufficient to reconstruct the trajectory continuously.
    fit_si = [fit[0] * replay.SOLAR_RADIUS_M, fit[1] * replay.SOLAR_RADIUS_M / 3600,
              fit[2] * replay.SOLAR_RADIUS_M / 3600**2, *axis, fit[5]]
    if args.legacy:
        # Preserve the original v1 regression input and its stronger nose
        # gate. A user event gets its own versioned provenance and results.
        require(first == start == 0 and last == HANDOFF_S and len(rays) == 49
                and event['reference_utc'] == REFERENCE and ambient == AMBIENT.resolve()
                and args.issue == replay.utc(timestamp(11400))
                and diagnostics['excluded'] == {'not_yet_available': 1}
                and diagnostics['picks'] == 332
                and outer == (WIND_M_S, DRAG_M_INV, END_S - HANDOFF_S, HANDOFF_MIN_M),
                'legacy fixture does not match its frozen input contract')
        header = ['cme-shock-evolution-v1', REFERENCE, json.dumps(str(ambient)),
                  ' '.join(format(x, '.17g') for x in fit_si),
                  f'{HANDOFF_S} {HANDOFF_MIN_M:.17g} {WIND_M_S} {DRAG_M_INV:.17g} {END_S - HANDOFF_S}']
    else:
        # The trace consumes only this receipt and its ambient file. Keep the
        # issue time and original observation checksum for later auditing.
        header = ['cme-shock-evolution-v2', event['reference_utc'],
                  args.issue.isoformat().replace('+00:00', 'Z'), hashlib.sha256(raw).hexdigest(),
                  json.dumps(str(args.data.resolve()), ensure_ascii=False),
                  json.dumps(str(ambient), ensure_ascii=False),
                  ' '.join(format(x, '.17g') for x in (first, last, start)),
                  ' '.join(format(x, '.17g') for x in fit_si)]
        header += ['fitted-only'] if outer is None else [
            'quadratic-drag ' + ' '.join(format(x, '.17g') for x in outer)]
    payload = '\n'.join([
        *header, str(len(rays)),
        *(' '.join(format(x, '.17g') for x in ray) for ray in rays),
    ]) + '\n'
    args.path.parent.mkdir(parents=True, exist_ok=True)
    temp = None
    try:
        with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=args.path.parent, delete=False) as stream:
            temp = Path(stream.name)
            stream.write(payload)
        os.replace(temp, args.path)
    finally:
        if temp is not None and temp.exists():
            temp.unlink()
    print(f'prepared {len(rays)} rays from {diagnostics["picks"]} calibrated picks for {start:g}–{end:g} s')


def check_archived_observations(args):
    # Feed the actual pack_observations.py CLI output to the fitter only
    # after comparison with the frozen HCI reference. This exercises the
    # filename and serialization path as well as the coordinate transform.
    raw = args.data.read_bytes()
    require(raw == args.reference.read_bytes(), 'packed HCI observations differ from the reference')
    event = json.loads(raw)
    sources = json.loads(args.sources.read_text())
    require(args.issue == replay.utc(sources['issue_utc']), 'archived observation issue time changed')
    fit_frames = {row['frame_file'] for row in event['tracks']}
    held_frames = {row['frame_file'] for row in event['holdout_tracks']}
    require(len(event['tracks']) == 37 and len(fit_frames) == 4
            and len(event['holdout_tracks']) == 22 and len(held_frames) == 3
            and fit_frames.isdisjoint(held_frames), 'archived exposure selection changed')
    require(all(replay.utc(row['available_utc']) <= args.issue for row in event['tracks']),
            'archived input contains tracks unavailable at the frozen cutoff')

    # The spacecraft header and independent Horizons solution should
    # place STEREO-A at the same HCI longitude, latitude and distance.
    # A 1,000 km discrepancy corresponds to about 1.5 arcsec here and
    # readily detects a sign, axis, or kilometre-to-metre mistake.
    discrepancies = []
    for frame in sources['frames']:
        if frame['instrument'] == 'STEREO-A/COR2':
            epoch, position, _ = pointing(frame, sources['ephemerides'])
            reference = interpolate_hci(sources['ephemerides'][frame['ephemeris']]['knots'], epoch)
            discrepancies.append(float(np.linalg.norm(position - reference)))
    require(len(discrepancies) == 3 and max(discrepancies) < 1e6,
            'STEREO-A archive and Horizons positions disagree in HCI')

    try:
        replay.fit_tracks(event, args.issue, set())
    except replay.FitBoundaryError as error:
        require(set(error.limits) == {'half_width'}
                and error.limits['half_width']['side'] == 'upper',
                f'archived event failed for a different reason: {error}')
        # Save an admission diagnostic only. No covariance, geometry map,
        # upstream realization or C++ event is produced from this minimum.
        report = {
            'status': 'REJECTED',
            'reason': 'SSE half-width reached its upper bound',
            'fit_boundary': error.limits,
            'trial_rms_angular_arcsec': error.rms_angular_rad * 180 * 3600 / math.pi,
            'observation_sha256': hashlib.sha256(raw).hexdigest(),
            'issue_utc': sources['issue_utc'],
            'fit_picks': len(event['tracks']),
            'fit_frames': sorted(fit_frames),
            'held_out_frames': sorted(held_frames),
            'max_stereo_position_difference_m': max(discrepancies),
        }
        replay.write_json(args.path, report)
        print(f'archived 2020-11-29 CME: {error}; '
              f'trial RMS {report["trial_rms_angular_arcsec"]:.1f} arcsec')
    else:
        raise ValueError('archived event unexpectedly fits within the reduced SSE domain; review admission')


def read_frame(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == list(HEADER), f'{path.name}: CSV schema mismatch')
        rows = list(reader)
    require(len(rows) == 49 and all(None not in row and None not in row.values() for row in rows),
            f'{path.name}: expected 49 complete ray records')
    return rows


def verify_manifest(path):
    # Recompute the native transactional artifact's FNV-1a checksum from its
    # bytes. This detects an incomplete write or a stale sidecar before CSV
    # metrics are read. It is a consistency check, not a physics reference.
    payload = path.read_bytes()
    lines = path.with_suffix('.csv.manifest').read_text().splitlines()
    checksum = 14695981039346656037
    for byte in payload:
        checksum = ((checksum ^ byte) * 1099511628211) & ((1 << 64) - 1)
    require(lines[:5] == ['SEP_ARTIFACT 1', 'schema cme-shock-evolution-v1',
                          'records 49', f'bytes {len(payload)}', f'checksum {checksum:016x}']
            and len(lines) == 7 and lines[6] == 'complete true',
            f'{path.name}: invalid completion manifest or payload checksum')
    return lines[5]


def independent_motion():
    # The oracle uses the manufactured truth rather than the prepared fitted
    # parameters. After 3 h it integrates the drag ODE independently of the
    # production closed form. Time, distance, and speed are SI throughout.
    r0, v0, a0 = TRUTH[:3]
    r0 *= replay.SOLAR_RADIUS_M
    v0 *= replay.SOLAR_RADIUS_M / 3600
    a0 *= replay.SOLAR_RADIUS_M / 3600**2
    clipped = np.minimum(EPOCHS, HANDOFF_S)
    radii = r0 + clipped * (v0 + 0.5 * clipped * a0)
    speeds = v0 + clipped * a0
    outer = EPOCHS > HANDOFF_S

    def rhs(_, y):
        relative = y[1] - WIND_M_S
        return y[1], -DRAG_M_INV * relative * abs(relative)

    ode = solve_ivp(rhs, (HANDOFF_S, END_S), [radii[~outer][-1], speeds[~outer][-1]],
                    t_eval=EPOCHS[outer], max_step=300, rtol=1e-11, atol=1e-7)
    require(ode.success, 'independent drag integration failed')
    radii[outer], speeds[outer] = ode.y
    return radii, speeds


def verify_geometry(rows, radius_m, speed_m_s):
    # Intersect the truth's generating sphere analytically, and compute the
    # outward normal velocity from the center and radius derivatives. No
    # projection fit, C++ geometry, or reported jump is used as this oracle.
    axis = replay.direction(*TRUTH[3:5])
    rays = np.asarray(replay.default_rays(TRUTH))
    expansion = math.sin(TRUTH[5])
    center = radius_m / (1 + expansion) * axis
    sphere_radius = radius_m * expansion / (1 + expansion)
    center_projection = rays @ center
    ray_distance = center_projection + np.sqrt(
        sphere_radius**2 - (center @ center - center_projection**2))
    positions = rays * ray_distance[:, None]
    normals = (positions - center) / sphere_radius
    normal_speeds = speed_m_s / (1 + expansion) * (normals @ axis + expansion)
    observed_positions = [[float(row[c]) for c in POSITIONS] for row in rows]
    observed_normals = [[float(row[c]) for c in NORMALS] for row in rows]
    observed_speeds = [float(row['normal_speed_m_s']) for row in rows]
    # Bounds cover numerical fitting and ODE integration roundoff. Even at
    # five hours a 2 m geometry error is below 1e-10 of the apex distance.
    np.testing.assert_allclose(observed_positions, positions, rtol=0, atol=2.0)
    np.testing.assert_allclose(observed_normals, normals, rtol=0, atol=2e-9)
    np.testing.assert_allclose(observed_speeds, normal_speeds, rtol=0, atol=0.005)


def verify_physics(rows):
    # A physical no-shock outcome carries no downstream state. Every solved
    # sample must be superfast, compressive, and conservative to the
    # prescribed ambient profile's 1e-8 RH bound at every snapshot.
    require(rows[0]['status'] == SOLVED, 'the nose must remain a fast shock')
    for row in rows:
        mach = float(row['fast_mach'])
        theta = float(row['theta_bn_rad'])
        require(math.isfinite(mach) and 0 <= theta <= math.pi / 2 + 1e-14,
                'nonfinite Mach number or invalid obliquity')
        if row['status'] == SOLVED:
            require(mach > 1 and 1 < float(row['density_compression']) < 4
                    and 0 <= float(row['max_rh_residual']) <= 1e-8,
                    'invalid fast-branch jump or RH residual')
            continue
        physical_absence = ((row['status'] == 'subfast-front' and 0 < mach <= 1)
                            or (row['status'] == 'non-forward-inflow' and mach <= 0))
        require(physical_absence and row['density_compression'] == ''
                and row['max_rh_residual'] == '', 'invalid no-shock classification or downstream state')


def compare_final(actual, reference):
    # Compare the five-hour numerical snapshot with the frozen reference.
    # RH residual magnitudes can vary at machine precision and already face
    # the 1e-8 analysis gate at every epoch.
    for observed, frozen in zip(actual, reference):
        for column in ('elapsed_s', 'ray_id', 'status'):
            require(observed[column] == frozen[column], f'final {column} differs from reference')
        for column in POSITIONS + NORMALS + ('normal_speed_m_s', 'fast_mach',
                                             'density_compression', 'theta_bn_rad'):
            if frozen[column] == '':
                require(observed[column] == '', f'final {column} must be absent')
            else:
                require(math.isclose(float(observed[column]), float(frozen[column]),
                                     rel_tol=2e-8, abs_tol=2e-9),
                        f'final ray {observed["ray_id"]} {column} differs from reference')


def check(directory, reference):
    names = [f'shock-{t:05d}.csv' for t in EPOCHS]
    require(sorted(file.name for file in directory.glob('shock-*.csv')) == names,
            'expected one initial snapshot and 15 successive 20-minute snapshots')
    radii, speeds = independent_motion()
    configuration = None
    for epoch, name, radius, speed in zip(EPOCHS, names, radii, speeds):
        output = directory / name
        fingerprint = verify_manifest(output)
        if configuration is None:
            configuration = fingerprint
        require(fingerprint == configuration, f'{name}: prepared event changed during the trace')
        rows = read_frame(output)
        require(all(row['elapsed_s'] == str(epoch) and row['ray_id'] == str(i)
                    for i, row in enumerate(rows)), f'{name}: epoch or ray order mismatch')
        verify_geometry(rows, radius, speed)
        verify_physics(rows)
    compare_final(rows, read_frame(reference))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    demo = sub.add_parser('fixture', help='write the manufactured observation JSON')
    demo.add_argument('path', type=Path)
    make_input = sub.add_parser('prepare')
    make_input.add_argument('path', type=Path)
    make_input.add_argument('--data', type=Path, required=True, help='calibrated observation-event JSON')
    make_input.add_argument('--issue', type=replay.utc, required=True, help='observation cutoff in UTC')
    make_input.add_argument('--ambient-member', type=int, help='zero-based index in ambient_members')
    make_input.add_argument('--withhold-view', action='append', default=[], help='omit a viewpoint')
    make_input.add_argument('--start-elapsed-s', type=float, help='seconds from reference_utc for first frame')
    make_input.add_argument('--legacy', action='store_true', help='apply the frozen legacy fixture contract')
    archived = sub.add_parser('check-archived-observations',
                              help='verify the frozen 2020-11-29 observations and the SSE admission result')
    archived.add_argument('path', type=Path, help='write the observational admission report')
    archived.add_argument('--data', type=Path, required=True)
    archived.add_argument('--sources', type=Path, required=True)
    archived.add_argument('--reference', type=Path,
                          default=HERE / 'observed_2020_11_29' / 'observations.json')
    archived.add_argument('--issue', type=replay.utc, required=True)
    validate = sub.add_parser('check')
    validate.add_argument('directory', type=Path)
    validate.add_argument('--reference', type=Path, default=HERE / 'final-reference.csv')
    args = parser.parse_args()
    try:
        if args.action == 'fixture':
            fixture(args.path)
        elif args.action == 'prepare':
            prepare(args)
        elif args.action == 'check-archived-observations':
            check_archived_observations(args)
        else:
            check(args.directory, args.reference)
    except (OSError, ValueError, AssertionError, KeyError, TypeError) as failure:
        parser.exit(1, f'CME shock evolution: {failure}\n')


if __name__ == '__main__':
    main()
