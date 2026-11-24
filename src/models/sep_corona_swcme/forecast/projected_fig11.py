#!/usr/bin/python3
"""Reproduce the projected 17 May 2012 shock-track replay from Figure 11."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import label

from replay import HOUR_S, SOLAR_RADIUS_M, write_json

FIGURE_SHA256 = '726c5f35075fdcacafcd83a1d9701d661d318a262ded688243aa427e7bd5fb60'
DOI = 'https://doi.org/10.1007/s11207-025-02449-5'


def projected_replay(path):
    """Read only the authors' COR2-A white-light shock overlay.

    This plot contains a red dashed guide along the shock feature, plus black
    CME and radio markers. Segmenting the red dashes avoids accidental use of
    ejecta positions. The dashes are samples of a drawn line and are correlated.
    """
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != FIGURE_SHA256:
        raise ValueError('unexpected Figure 11 edition or raster size')
    with Image.open(path) as image:
        red, green, blue = np.asarray(image.convert('RGB')).transpose(2, 0, 1)
    red_track = (red > 150) & (green < 75) & (blue < 75)
    components, count = label(red_track, np.ones((3, 3), int))
    points = []
    for index in range(1, count + 1):
        rows, cols = np.nonzero(components == index)
        if len(rows) > 20:
            # Figure 11 axes are time from UTC midnight (0 to 8 hours) and
            # projected height from Sun centre (0 to 30 nominal solar radii).
            # Half-pixel origin values locate the centres of the axis strokes.
            hour = 8.0 * (cols.mean() - 97.5) / (785.5 - 97.5)
            radius = 30.0 * (514.5 - rows.mean()) / (514.5 - 2.5)
            points.append((hour, radius))
    if len(points) != 5:
        raise ValueError(f'expected five red shock-track dashes, found {len(points)}')
    points.sort()
    hours, radii = np.asarray(points).T
    speed, intercept = np.polyfit(hours, radii, 1)
    rms = np.sqrt(np.mean((speed * hours + intercept - radii)**2))
    return {
        'status': 'OBSERVED_PROJECTED_TRACK_ONLY',
        'reference_utc': '2012-05-17T00:00:00Z',
        'source': DOI, 'figure': 11, 'figure_sha256': digest,
        'view': 'STEREO-A SECCHI/COR2',
        'track': [{'elapsed_s': hour * HOUR_S, 'projected_radius_rs': height}
                  for hour, height in points],
        'projected_speed_m_s': float(speed * SOLAR_RADIUS_M / HOUR_S),
        'line_rms_rs': float(rms),
        'three_dimensional_geometry': 'SKIP_NO_ASTROMETRIC_FRONT_PICKS',
        'shock_strength': 'SKIP_NO_QUALIFIED_UPSTREAM',
        'observer_forecast': 'SKIP_NO_3D_GEOMETRY',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figure', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        write_json(args.output, projected_replay(args.figure))
    except (ValueError, OSError) as error:
        parser.exit(2, f'projected replay: {error}\n')


if __name__ == '__main__':
    main()
