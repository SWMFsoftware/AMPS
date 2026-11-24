#!/usr/bin/python3
"""Replay time-gated coronagraph front picks with the AMPS finite SSE core."""
import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import NamedTuple, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

SOLAR_RADIUS_M = 695_700_000.0
AU_M = 149_597_870_700.0
# Fit in solar radii and hours, with angles in radians. Scaled parameters keep
# the optimizer and its rank test from comparing SI metres against radians.
HOUR_S = 3600.0
SHOCK = 'solved-fast-shock'
NONSHOCK = {'subfast-front', 'non-forward-inflow', 'outside-front-support',
            'below-physical-inner-boundary'}


class FrontState(NamedTuple):
    parameters: Optional[np.ndarray]
    # The core wire format uses SI; the fitter retains solar radii and hours.
    apex: Optional[Tuple[float, float]]
    status: str


def utc(text):
    """Use explicit UTC and a caller-supplied elapsed SI time for dynamics."""
    if not text.endswith('Z'):
        raise ValueError('UTC timestamps must end in Z')
    return datetime.fromisoformat(text[:-1] + '+00:00')


def finite(value, name, positive=False):
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'invalid {name}')
    return value


def vector(value, name, unit=False):
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f'invalid {name} vector')
    if unit and abs(np.linalg.norm(result) - 1) > 1e-10:
        raise ValueError(f'{name} must be a unit vector')
    return result


def direction(longitude, latitude):
    return np.array([np.cos(latitude) * np.cos(longitude),
                     np.cos(latitude) * np.sin(longitude), np.sin(latitude)])


def kinematics(parameters, hour):
    r, v, a = parameters[:3]
    return r + hour * (v + 0.5 * hour * a), v + hour * a


def projection(parameters, hours, observers, sightlines):
    """Exact angular tangency to the finite-distance SSE generating sphere.

    The second returned term is the signed radial-normal cosine at the
    predicted tangent. A negative value selects the sphere's rear surface,
    which is outside the physical leading SSE cap.
    """
    longitude, latitude, width = parameters[3:]
    d = direction(longitude, latitude)
    r, _ = kinematics(parameters, hours)
    sine_width = np.sin(width)
    centers = r[:, None] * d / (1 + sine_width)
    radii = r * sine_width / (1 + sine_width)
    separation = centers - observers
    distance = np.linalg.norm(separation, axis=1)
    axis = separation / np.maximum(distance, 1e-16)[:, None]
    scalar = np.einsum('ij,ij->i', sightlines, axis)
    lateral = sightlines - scalar[:, None] * axis
    lateral_norm = np.linalg.norm(lateral, axis=1)
    # Limit only the trial iterate passed to asin. The associated domain
    # penalty below rejects inside-camera and nonphysical radius solutions.
    sine = np.clip(radii / np.maximum(distance, 1e-16), 0, 1)
    ring_normal = (-sine[:, None] * axis
                   + np.sqrt(1 - sine**2)[:, None] * lateral
                   / np.maximum(lateral_norm, 1e-16)[:, None])
    tangent = centers + radii[:, None] * ring_normal
    outward = (np.einsum('ij,ij->i', tangent, ring_normal)
               / np.maximum(np.linalg.norm(tangent, axis=1), 1e-16))
    domain = np.minimum.reduce([r, distance - radii, lateral_norm])
    return np.arctan2(lateral_norm, scalar) - np.arcsin(sine), outward, domain


def fit_tracks(event, issue, withheld):
    tracks = []
    rejected = Counter()
    epoch_utc = utc(event['reference_utc'])
    if issue < epoch_utc:
        raise ValueError('issue time precedes the event reference')
    for row in event['tracks']:
        observed = utc(row['observed_utc'])
        available = utc(row['available_utc'])
        if available < observed:
            raise ValueError('track availability precedes acquisition')
        if available > issue:
            rejected['not_yet_available'] += 1
        elif row['view'] in withheld:
            rejected['withheld_view'] += 1
        else:
            tracks.append(row)
    elapsed = np.array([finite(row['elapsed_s'], 'track elapsed seconds') for row in tracks])
    for row, seconds in zip(tracks, elapsed):
        if seconds < 0 or seconds > (utc(row['observed_utc']) - epoch_utc).total_seconds():
            raise ValueError('front epoch must lie between the reference and acquisition time')
    fit_epochs = np.unique(elapsed)
    hours = elapsed / HOUR_S
    # Observation UTC is receive/exposure time. Elapsed SI time is the caller's
    # epoch for the front model after applying the declared light-time policy.
    # Inputs deliberately carry both so issue-time gating cannot leak a shot.
    if len(fit_epochs) < 3 or len(hours) < 8:
        raise ValueError('fit needs at least eight picks spanning three epochs')
    observers = np.array([vector(row['observer_hci_m'], 'camera') / SOLAR_RADIUS_M
                          for row in tracks])
    sightlines = np.array([vector(row['los_hci'], 'sightline', unit=True) for row in tracks])
    sigma = np.array([finite(row['sigma_rad'], 'angular uncertainty', positive=True)
                      for row in tracks])
    seed = np.array(event['initial_fit'], dtype=float)
    if seed.shape != (6,) or not np.isfinite(seed).all():
        raise ValueError('initial_fit must contain six finite SSE parameters')
    fit_acceleration = event['fit_acceleration']
    if not isinstance(fit_acceleration, bool):
        raise ValueError('fit_acceleration must be Boolean')
    active = np.array([0, 1, 2, 3, 4, 5] if fit_acceleration else [0, 1, 3, 4, 5])
    if not fit_acceleration:
        seed[2] = 0.0
    lower = np.array([1.0, 0, -np.inf, -np.inf, -np.pi / 2, 1e-4])
    upper = np.array([np.inf, np.inf, np.inf, np.inf, np.pi / 2, np.pi / 2 - 1e-4])

    def expand(active_parameters):
        result = seed.copy()
        result[active] = active_parameters
        return result

    def residual(active_parameters):
        p = expand(active_parameters)
        angular, outward, domain = projection(p, hours, observers, sightlines)
        _, speed = kinematics(p, hours[[np.argmin(hours), np.argmax(hours)]])
        # Inequality penalties guide trial iterates back to the supported cap.
        # The final fit is checked explicitly before any physics is evaluated.
        return np.concatenate([angular / sigma, 1000 * np.minimum(outward, 0),
                               1000 * np.minimum(domain, 0), 1000 * np.minimum(speed, 0)])

    solution = least_squares(residual, seed[active], bounds=(lower[active], upper[active]),
                             x_scale='jac', max_nfev=250)
    center = expand(solution.x)
    angular, outward, domain = projection(center, hours, observers, sightlines)
    limits = np.array([hours.min(), hours.max()])
    radii, speeds = kinematics(center, limits)
    if (not solution.success or np.any(solution.active_mask) or np.min(outward) < -1e-8
            or np.min(domain) <= 0 or np.min(speeds) <= 0 or np.min(radii) <= 1):
        raise ValueError('no regular advancing fit on the leading SSE cap')

    # Only genuine measurements contribute to the information matrix. An
    # active domain penalty cannot make an otherwise unidentifiable fit look
    # constrained. No fitted covariance is published for ill-conditioned data.
    _, singular, axes = np.linalg.svd(solution.jac[:len(tracks)], full_matrices=False)
    if singular[-1] <= singular[0] * 1e-6:
        raise ValueError('geometry is underdetermined by the supplied tracks')
    chi_squared = float(np.sum((angular / sigma)**2))
    reduced_chi_squared = chi_squared / (len(tracks) - len(active))
    covariance_scale = math.sqrt(max(1.0, reduced_chi_squared))
    rng = np.random.default_rng(event['seed'])
    samples = event['ensemble_samples']
    if type(samples) is not int or not 1 <= samples <= 1024:
        raise ValueError('ensemble_samples must be an integer between 1 and 1024')
    draws = rng.standard_normal((samples, len(active))) @ (
        covariance_scale * axes / singular[:, None])
    ensemble = []
    for draw in draws:
        p = expand(solution.x + draw)
        radii, speeds = kinematics(p, limits)
        _, valid_cap, valid_distance = projection(p, hours, observers, sightlines)
        physical = (np.all(p[active] > lower[active]) and np.all(p[active] < upper[active])
                    and min(radii) > 1.0 and min(speeds) > 0
                    and min(valid_cap) >= 0 and min(valid_distance) > 0)
        ensemble.append(p if physical else None)

    return center, ensemble, fit_epochs, {
        'picks': len(tracks),
        'views': sorted(set(row['view'] for row in tracks)),
        'excluded': dict(rejected),
        'condition_number': float(singular[0] / singular[-1]),
        'rms_angular_rad': float(np.sqrt(np.mean(angular**2))),
        'reduced_chi_squared': reduced_chi_squared,
        'parameters': center.tolist(),
    }


class Core:
    """Keep the production ambient and shock kernels behind a small batch CLI."""

    def __init__(self, executable, paths):
        self.command = [str(executable.resolve()), *(str(path) for path in paths)]

    def batch(self, rows):
        if not rows:
            return []
        body = ''.join(' '.join(map(str, row)) + '\n' for row in rows)
        proc = subprocess.run(self.command, input=body, capture_output=True,
                              text=True, check=False)
        if proc.returncode:
            raise ValueError(proc.stderr.strip() or 'forecast core failed')
        if proc.stderr:
            # A failed upstream member remains unresolved, with its original
            # provider diagnostic available in the replay log.
            sys.stderr.write(proc.stderr)
        result = [line.split() for line in proc.stdout.splitlines()]
        if len(result) != len(rows):
            raise ValueError('forecast core emitted an incomplete batch')
        return result


def ray_request(parameters, elapsed_s, ray, member, offset, outer=None):
    r, v = (kinematics(parameters, elapsed_s / HOUR_S) if outer is None else outer)
    return ['R', r * SOLAR_RADIUS_M, v * SOLAR_RADIUS_M / HOUR_S,
            *direction(*parameters[3:5]), parameters[5], *ray,
            member, elapsed_s + offset]


def quantiles(values):
    if not values:
        return None
    return np.quantile(values, [0.05, 0.50, 0.95]).tolist()


def summarize(replies):
    counts = Counter(reply[0] for reply in replies)
    solved = counts[SHOCK]
    known_absent = sum(counts[name] for name in NONSHOCK)
    # Keep inaccessible/failed draws in the denominator. These lower/upper
    # fractions express what the ensemble can resolve, rather than assigning
    # every unresolved member to either a hit or a miss.
    unresolved = len(replies) - solved - known_absent
    physical = [reply for reply in replies if len(reply) == 11]
    return {
        'members': len(replies), 'statuses': dict(counts),
        'fast_shock_fraction_bounds': [solved / len(replies),
                                       (solved + unresolved) / len(replies)],
        'unresolved_fraction': unresolved / len(replies),
        'normal_speed_m_s_p05_p50_p95': quantiles([float(r[7]) for r in physical]),
        'fast_mach_given_solved_p05_p50_p95': quantiles(
            [float(r[8]) for r in physical if r[0] == SHOCK]),
        'compression_given_solved_p05_p50_p95': quantiles(
            [float(r[9]) for r in physical if r[0] == SHOCK]),
        'obliquity_rad_given_solved_p05_p50_p95': quantiles(
            [float(r[10]) for r in physical if r[0] == SHOCK and math.isfinite(float(r[10]))]),
    }


def report_map(core, states, elapsed_s, rays, members):
    backgrounds = list(enumerate(members)) or [(-1, 0)]
    best, *ensemble = states
    jobs = []
    for ray in rays:
        if best.apex is not None:
            jobs.append(ray_request(best.parameters, elapsed_s, ray, -1, 0, outer=best.apex))
        for state in ensemble:
            if state.apex is not None:
                jobs.extend(ray_request(state.parameters, elapsed_s, ray, i, offset,
                                        outer=state.apex)
                            for i, offset in backgrounds)
    replies = iter(core.batch(jobs))
    result = []
    for ray in rays:
        nominal = next(replies) if best.apex is not None else [best.status]
        # Preserve the same set of draws and weights at every epoch, including
        # those with no valid handoff or with no valid fit to the observations.
        values = [next(replies) if state.apex is not None else [state.status]
                  for state in ensemble for _ in backgrounds]
        result.append({'ray_hci': ray.tolist(), 'nominal_geometry': {
            'status': 'FRONT_GEOMETRY' if len(nominal) == 11 else nominal[0],
            'position_hci_m': list(map(float, nominal[1:4])) if len(nominal) == 11 else None,
            'normal_hci': list(map(float, nominal[4:7])) if len(nominal) == 11 else None,
            'normal_speed_m_s': float(nominal[7]) if len(nominal) == 11 else None,
        }, **summarize(values)})
    return result


def resolve_propagation(propagation):
    wind = finite(propagation['wind_m_s'], 'effective wind', positive=True)
    drag = finite(propagation['drag_m_inv'], 'drag coefficient')
    horizon = finite(propagation['horizon_s'], 'forecast horizon', positive=True)
    min_handoff = finite(propagation['minimum_handoff_apex_m'], 'minimum handoff radius',
                         positive=True)
    if drag < 0 or horizon > 14 * 86400:
        raise ValueError('drag must be nonnegative and horizon at most 14 days')
    return wind, drag, horizon, min_handoff


def drag_states(core, parameters, elapsed_s, last_s, propagation):
    wind, drag, _, min_handoff = propagation
    # Match every realization at the final fitted epoch. The launch-fit
    # polynomial ends there. The same r and v initialize the outer history.
    handoff = [kinematics(p, last_s / HOUR_S) if p is not None else None for p in parameters]
    jobs = [['T', elapsed_s - last_s, r * SOLAR_RADIUS_M,
             v * SOLAR_RADIUS_M / HOUR_S, wind, drag]
            for r, v in (apex for apex in handoff if apex is not None)
            if r * SOLAR_RADIUS_M >= min_handoff]
    replies = iter(core.batch(jobs))
    result = []
    for p, apex in zip(parameters, handoff):
        if apex is None:
            result.append(FrontState(p, None, 'UNRESOLVED_FIT_DRAW'))
        elif apex[0] * SOLAR_RADIUS_M < min_handoff:
            result.append(FrontState(p, None, 'BELOW_DRAG_HANDOFF'))
        else:
            label, r, v = next(replies)
            if label != 'DRAG_STATE':
                raise ValueError('forecast core returned an invalid drag state')
            result.append(FrontState(p, (float(r) / SOLAR_RADIUS_M,
                                         float(v) * HOUR_S / SOLAR_RADIUS_M),
                                     'quadratic-drag'))
    return result


def report_evolution(core, parameters, fit_epochs, requested, rays, members, propagation):
    first_s, last_s = fit_epochs[0], fit_epochs[-1]
    if requested is None:
        epochs = fit_epochs.tolist()
    else:
        if not isinstance(requested, list) or not requested or any(
                type(value) not in (float, int) for value in requested):
            raise ValueError('evolution_elapsed_s must be a nonempty list of elapsed seconds')
        epochs = sorted({finite(t, 'evolution elapsed seconds') for t in requested})
    if epochs[0] < first_s:
        raise ValueError('evolution request precedes the first usable front epoch')
    if epochs[-1] > last_s:
        if not propagation:
            raise ValueError('evolution beyond the fitted interval requires configured propagation')
        propagation = resolve_propagation(propagation)
        if epochs[-1] - last_s > propagation[2]:
            raise ValueError('evolution request exceeds the declared drag horizon')

    frames = {}
    # The same HCI rays and fit draws are used throughout this issue. Evaluate
    # the actual surface intersection and upstream state at every requested
    # epoch. Interpolating already-classified jump values would hide transitions
    # between absent, subfast, fast and unresolved states.
    for elapsed_s in sorted(set(epochs) | {last_s}):
        if elapsed_s <= last_s:
            source = 'fitted-history'
            states = [FrontState(p, kinematics(p, elapsed_s / HOUR_S), source)
                      if p is not None else FrontState(p, None, 'UNRESOLVED_FIT_DRAW')
                      for p in parameters]
        else:
            source = 'quadratic-drag'
            states = drag_states(core, parameters, elapsed_s, last_s, propagation)
        best = states[0]
        frames[elapsed_s] = {
            'elapsed_s': float(elapsed_s), 'source': source, 'trajectory_status': best.status,
            'apex_radius_m': best.apex[0] * SOLAR_RADIUS_M if best.apex is not None else None,
            'apex_speed_m_s': (best.apex[1] * SOLAR_RADIUS_M / HOUR_S
                               if best.apex is not None else None),
            'map': report_map(core, states, elapsed_s, rays, members),
        }
    return [frames[t] for t in epochs], frames[last_s]['map']


def report_arrival(core, ensemble, last_s, issue_elapsed, observers, members, propagation):
    if not observers:
        return {}
    if not propagation:
        return {name: {'status': 'PROPAGATION_UNAVAILABLE'} for name in observers}
    wind, drag, horizon, min_handoff = resolve_propagation(propagation)
    backgrounds = list(enumerate(members)) or [(-1, 0)]
    result = {}
    for name, xyz in observers.items():
        position = vector(xyz, 'fixed forecast observer')
        radius = float(np.linalg.norm(position))
        # A norm of an exactly 1-AU Cartesian position can round just above
        # AU_M. Accept that representation error without moving the observer.
        if not SOLAR_RADIUS_M < radius <= AU_M + 8 * np.spacing(AU_M):
            raise ValueError('forecast observer radius must lie above the Sun and at most 1 AU')
        ray = position / radius
        geometries = iter(core.batch([ray_request(p, last_s, ray, -1, 0)
                                     for p in ensemble if p is not None]))
        outcomes = []
        queries = []
        pending = []
        for p in ensemble:
            if p is None:
                outcomes.append(['UNRESOLVED_FIT_DRAW'])
                continue
            geometry = next(geometries)
            if len(geometry) == 1:
                outcomes.append(geometry)
                continue
            r, v = kinematics(p, last_s / HOUR_S)
            if r * SOLAR_RADIUS_M < min_handoff:
                outcomes.append(['BELOW_DRAG_HANDOFF'])
                continue
            ray_radius = np.linalg.norm(np.asarray(geometry[1:4], dtype=float))
            # Fixed SSE scales linearly with apex radius on any supported
            # heliocentric ray. Reach the observer on that ray, not at the nose.
            target = r * radius / ray_radius
            if target <= r:
                outcomes.append(['ALREADY_PASSED_OR_IN_RANGE'])
                continue
            queries.append(['A', r * SOLAR_RADIUS_M, v * SOLAR_RADIUS_M / HOUR_S,
                            wind, drag, target * SOLAR_RADIUS_M, horizon])
            pending.append((p, target))
        physical_jobs = []
        crossing_epochs = []
        for reply, (p, target) in zip(core.batch(queries), pending):
            if reply[0] != 'ARRIVAL':
                outcomes.append(reply)
                continue
            elapsed = last_s + float(reply[1])
            if elapsed <= issue_elapsed:
                outcomes.append(['CROSSING_BEFORE_ISSUE'])
                continue
            v = float(reply[2]) * HOUR_S / SOLAR_RADIUS_M
            for i, offset in backgrounds:
                physical_jobs.append(ray_request(p, elapsed, ray, i, offset,
                                                 outer=(target, v)))
                crossing_epochs.append(elapsed)
        physical = core.batch(physical_jobs)
        all_members = [outcome for outcome in outcomes for _ in backgrounds] + physical
        hits = [epoch for epoch, reply in zip(crossing_epochs, physical)
                if reply[0] == SHOCK]
        result[name] = {'arrival_elapsed_s_given_solved_p05_p50_p95': quantiles(hits),
                        'front_crossing_elapsed_s_p05_p50_p95': quantiles(crossing_epochs),
                        **summarize(all_members)}
    return result


def default_rays(parameters):
    # Sample heliocentric directions uniformly in solid angle on concentric
    # rings. This is a point map. It carries no implicit curved-area weights.
    d = direction(*parameters[3:5])
    axis = np.eye(3)[np.argmin(np.abs(d))]
    u = np.cross(d, axis)
    u /= np.linalg.norm(u)
    v = np.cross(d, u)
    rays = [d]
    for fraction in (0.25, 0.50, 0.75, 0.99):
        cosine = 1 - fraction * (1 - np.cos(parameters[5]))
        sine = np.sqrt(1 - cosine**2)
        rays.extend(cosine * d + sine * (np.cos(phi) * u + np.sin(phi) * v)
                    for phi in np.arange(12) * (2 * np.pi / 12))
    return rays


def replay(event, issue, withheld, executable, directory):
    if event['schema'] != 'shock-front-replay-v1' or event['frame'] != 'HCI':
        raise ValueError('expected the shock-front-replay-v1 schema and HCI')
    if not event['light_time_policy'] or not event['provenance']:
        raise ValueError('document the elapsed-time convention and observation provenance')
    center, draws, fit_epochs, fit_metrics = fit_tracks(event, issue, withheld)
    selected = []
    excluded_ambient = []
    for member in event.get('ambient_members', []):
        availability = utc(member['available_utc'])
        path = (directory / member['event_file']).resolve()
        if availability <= issue:
            selected.append((path, availability))
        else:
            excluded_ambient.append({'path': str(path), 'available_utc': availability.isoformat()})
    paths = [path for path, _ in selected]
    if len(paths) != len(set(paths)):
        raise ValueError('duplicate ambient event paths would change ensemble weights')
    core = Core(executable, paths)
    event_epoch = utc(event['reference_utc'])
    receipts = []
    members = []
    definitions = core.batch([['M', i] for i in range(len(paths))])
    for (path, availability), (epoch, fingerprint) in zip(selected, definitions):
        # Event time is independent of the epoch of every upstream realization.
        # Use its resolved UTC origin, never the reference of a different
        # member. Replay intervals are limited to dates without leap seconds.
        offset = (event_epoch - utc(epoch)).total_seconds()
        members.append(offset)
        receipts.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'physics_fingerprint': fingerprint, 'reference_utc': epoch,
                         'available_utc': availability.isoformat(),
                         'event_to_ambient_offset_s': offset})
    rays = ([vector(ray, 'sample ray', unit=True) for ray in event['rays_hci']]
            if 'rays_hci' in event else default_rays(center))
    evolution, latest_map = report_evolution(core, [center, *draws], fit_epochs,
                                             event.get('evolution_elapsed_s'), rays,
                                             members, event.get('propagation'))
    return {
        'schema': 'shock-front-replay-result-v1', 'status': 'EXPERIMENTAL',
        'reference_utc': event['reference_utc'],
        'issue_utc': issue.isoformat(), 'map_elapsed_s': float(fit_epochs[-1]),
        'fit_interval_elapsed_s': [float(fit_epochs[0]), float(fit_epochs[-1])],
        'light_time_policy': event['light_time_policy'], 'provenance': event['provenance'],
        'fit': fit_metrics, 'propagation': event.get('propagation'),
        'ambient_members': receipts, 'ambient_not_available_at_issue': excluded_ambient,
        'evolution': evolution, 'map': latest_map,
        'observers': report_arrival(core, draws, fit_epochs[-1],
                                    (issue - utc(event['reference_utc'])).total_seconds(),
                                    event.get('observers_hci_m', {}), members,
                                    event.get('propagation')),
    }


def write_json(path, data):
    """Publish only a complete, finite JSON result, including missing values."""
    tmp = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as stream:
            tmp = stream.name
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write('\n')
        os.replace(tmp, path)
    finally:
        if tmp is not None and os.path.exists(tmp):
            os.unlink(tmp)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('event', type=Path, help='calibrated HCI picks and replay configuration')
    parser.add_argument('--issue', type=utc, required=True, help='observation cutoff in UTC')
    parser.add_argument('--withhold-view', action='append', default=[], help='omit this camera')
    parser.add_argument('--core', type=Path, required=True, help='built shock_forecast_core')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        raw = args.event.read_bytes()
        result = replay(json.loads(raw), args.issue, set(args.withhold_view),
                        args.core, args.event.parent)
        result['event_sha256'] = hashlib.sha256(raw).hexdigest()
        result['core_sha256'] = hashlib.sha256(args.core.read_bytes()).hexdigest()
        write_json(args.output, result)
    except (ValueError, OSError, TypeError, KeyError, OverflowError) as error:
        parser.exit(2, f'shock replay: {error}\n')


if __name__ == '__main__':
    main()
