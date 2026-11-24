"""Manufactured coronagraph tangencies for C++ and Python shock regressions."""
from datetime import timedelta
from pathlib import Path

import numpy as np

import replay

AMPS = Path(__file__).resolve().parents[4]
AMBIENT = AMPS / 'srcSEP/examples/reduced-shock/positive_1au.event'
REFERENCE = '2026-10-06T00:00:00Z'
TRUTH = np.array([8.0, 7.50, 0.25, 0.30, 1.06, 0.82])


def timestamp(t):
    return (replay.utc(REFERENCE) + timedelta(seconds=t)).isoformat().replace('+00:00', 'Z')


def synthetic(noise=False):
    """Generate measured rays with the independent sphere tangent construction.

    Each ray is built from a surface point whose normal has a dot product of
    zero with the camera ray, using a right triangle in physical space. This
    fixture deliberately omits later frames and includes a delayed outlier.
    """
    r0, v0, a, lon, lat, width = TRUTH
    nose = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    rng = np.random.default_rng(7218)
    tracks = []
    for name, theta in (('A', 0.0), ('B', 1.7), ('C', 3.6)):
        observer = np.array([np.cos(theta), np.sin(theta), 0]) * replay.AU_M
        for hour in (0.0, 0.33, 1.2, 2.0, 3.0):
            if name == 'B' and hour == 1.2:
                continue  # A real imaging gap leaves the other views intact.
            r = (r0 + v0 * hour + a * hour**2 / 2) * replay.SOLAR_RADIUS_M
            radius = r * np.sin(width) / (1 + np.sin(width))
            center = (r - radius) * nose
            view = (observer - center) / np.linalg.norm(observer - center)
            across = np.cross(view, [0, 0, 1])
            across /= np.linalg.norm(across)
            up = np.cross(across, view)
            term = radius / np.linalg.norm(observer - center)
            sigma = 1.0e-6
            for phi in np.arange(32) * (2 * np.pi / 32):
                transverse = np.cos(phi) * across + np.sin(phi) * up
                normal = term * view + np.sqrt(1 - term**2) * transverse
                point = center + radius * normal
                if np.dot(point, normal) < 0.12 * np.linalg.norm(point):
                    continue
                los = (point - observer) / np.linalg.norm(point - observer)
                if noise:
                    error = rng.normal(size=3)
                    los += sigma * (error - np.dot(error, los) * los)
                    los /= np.linalg.norm(los)
                tracks.append({'view': name, 'observed_utc': timestamp(hour * 3600 + 500),
                               'available_utc': timestamp(hour * 3600 + 545),
                               'elapsed_s': hour * 3600, 'observer_hci_m': observer.tolist(),
                               'los_hci': los.tolist(), 'sigma_rad': sigma})
    tracks.append({**tracks[-1], 'elapsed_s': 4 * 3600, 'observed_utc': timestamp(14900),
                   'available_utc': timestamp(14945), 'los_hci': [1, 0, 0]})
    return {
        'schema': 'shock-front-replay-v1', 'frame': 'HCI',
        'reference_utc': REFERENCE, 'light_time_policy': 'manufactured 500 s delay',
        'provenance': 'independent manufactured tangent rays', 'seed': 1729,
        'ensemble_samples': 16, 'fit_acceleration': True,
        'initial_fit': [8.3, 7.3, 0.38, 0.23, 0.99, 0.75], 'tracks': tracks,
        'rays_hci': [nose.tolist(), (-nose).tolist()],
        'ambient_members': [{'event_file': str(AMBIENT), 'available_utc': REFERENCE}],
        'observers_hci_m': {'nose': (nose * replay.AU_M).tolist(),
                            'far_side': (-nose * replay.AU_M).tolist()},
        'propagation': {'wind_m_s': 400_000, 'drag_m_inv': 0, 'horizon_s': 140_000,
                        'minimum_handoff_apex_m': 20 * replay.SOLAR_RADIUS_M},
    }
