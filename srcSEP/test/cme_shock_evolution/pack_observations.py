#!/usr/bin/python3
"""Turn rectified SECCHI/COR2 and LASCO pixel marks into AMPS replay rays."""

import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ARCSEC = math.pi / (180 * 3600)
C_M_S = 299792458.0
# Horizons geometric ICRF vectors use the IAU76 mean ecliptic. Rotate the
# IAU solar pole (ICRF RA and Dec) into that same frame before forming HCI.
J2000_OBLIQUITY = 84381.448 * ARCSEC
SOLAR_POLE_RA, SOLAR_POLE_DEC = np.deg2rad([286.13, 63.87])
KEYS = ('NAXIS1 NAXIS2 CTYPE1 CTYPE2 CUNIT1 CUNIT2 CRPIX1 CRPIX2 '
        'CRVAL1 CRVAL2 CDELT1 CDELT2 DATE-OBS EXPTIME RECTIFY').split()
SECCHI_KEYS = ('DATE-AVG DATE-END PC1_1 PC1_2 PC2_1 PC2_2 HCIX_OBS HCIY_OBS '
               'HCIZ_OBS').split()
LASCO_KEYS = ['TIME-OBS', 'CROTA2']


def utc(value):
    if not value.endswith('Z'):
        raise ValueError(f'explicit UTC required: {value}')
    return datetime.fromisoformat(value[:-1] + '+00:00')


def unit(vector):
    return vector / np.linalg.norm(vector)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primary_header(stream):
    """Read FITS cards and leave the stream at the start of the data block."""
    values = {}
    for block in iter(lambda: stream.read(2880), b''):
        for offset in range(0, len(block), 80):
            card = block[offset:offset + 80].decode('ascii')
            name = card[:8].strip()
            if name == 'END':
                return values
            if card[8:10] == '= ':
                field = card[10:].strip()
                values[name] = (field.split("'")[1] if field.startswith("'")
                                else field.split('/')[0]).strip()
    raise ValueError(f'missing FITS END card: {stream.name}')


def archive_header(path, instrument):
    """Preserve archive precision for the cards used to pack the sight lines."""
    with path.open('rb') as stream:
        values = primary_header(stream)
    extra = SECCHI_KEYS if instrument == 'STEREO-A/COR2' else LASCO_KEYS
    return {key: values[key] for key in KEYS + extra}


def horizons_rows(path):
    """Preserve the geometric centre-of-Sun ecliptic positions and UTC knots."""
    content = path.read_text()
    if ('Ecliptic of J2000.0' not in content or 'GEOMETRIC' not in content.upper()
            or 'KM-S' not in content or 'Sun (10)' not in content):
        raise ValueError('expected Horizons geometric Sun-centred J2000 ecliptic KM-S')
    table = content.split('$$SOE')[1].split('$$EOE')[0]
    rows = []
    for fields in csv.reader(table.strip().splitlines()):
        epoch = datetime.strptime(fields[1].strip(), 'A.D. %Y-%b-%d %H:%M:%S.%f')
        rows.append([epoch.isoformat() + 'Z', *map(float, fields[2:5])])
    return rows


def ecliptic_to_hci(xyz_km):
    # HCI X lies along the ascending node of the solar equator on the
    # J2000 ecliptic. Completing the right-handed triad fixes the longitude
    # zero and avoids treating an Earth-relative longitude as inertial.
    c, s = math.cos(J2000_OBLIQUITY), math.sin(J2000_OBLIQUITY)
    equatorial_to_ecliptic = np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
    pole = np.array([math.cos(SOLAR_POLE_DEC) * math.cos(SOLAR_POLE_RA),
                     math.cos(SOLAR_POLE_DEC) * math.sin(SOLAR_POLE_RA),
                     math.sin(SOLAR_POLE_DEC)])
    z = equatorial_to_ecliptic @ pole
    x = unit(np.cross([0, 0, 1], z))
    y = np.cross(z, x)
    return 1000 * (np.array([x, y, z]) @ xyz_km)


def interpolate_hci(knots, epoch):
    seconds = np.array([utc(row[0]).timestamp() for row in knots])
    if np.any(np.diff(seconds) <= 0) or not seconds[0] <= epoch.timestamp() <= seconds[-1]:
        raise ValueError(f'ordered Horizons UTC knots must bracket {epoch.isoformat()}')
    positions = np.array([row[1:] for row in knots])
    return ecliptic_to_hci([np.interp(epoch.timestamp(), seconds, values)
                            for values in positions.T])


def pointing(frame, ephemerides):
    """Resolve the pointing and exposure epoch without reorienting the raster."""
    h = frame['header']
    axes = h['CTYPE1'], h['CTYPE2']
    if h['RECTIFY'] not in ('T', 'TRUE') or any(
            h[key].lower() != 'arcsec' for key in ('CUNIT1', 'CUNIT2')):
        raise ValueError(f'expected rectified arcsecond WCS for {frame["file"]}')
    if frame['instrument'] == 'STEREO-A/COR2' and axes == ('HPLN-TAN', 'HPLT-TAN'):
        epoch = utc(h['DATE-AVG'] + 'Z')
        o = np.array([float(h[f'HCI{axis}_OBS']) for axis in 'XYZ'])
        pc = [[float(h[f'PC{i}_{j}']) for j in (1, 2)] for i in (1, 2)]
    elif (frame['instrument'] in ('SOHO/LASCO-C2', 'SOHO/LASCO-C3')
          and axes == ('SOLAR-X', 'SOLAR-Y')):
        start = utc(h['DATE-OBS'].replace('/', '-') + 'T' + h['TIME-OBS'] + 'Z')
        epoch = start + timedelta(seconds=0.5 * float(h['EXPTIME']))
        o = interpolate_hci(ephemerides[frame['ephemeris']]['knots'], epoch)
        # These legacy LASCO headers use equal scales and CROTA2 for the
        # roll. Restrict this adapter to those headers so unequal scales
        # cannot be rotated in the wrong order.
        if float(h['CDELT1']) != float(h['CDELT2']):
            raise ValueError('unequal LASCO scales require an explicit WCS matrix')
        roll = math.radians(float(h['CROTA2']))
        c, s = math.cos(roll), math.sin(roll)
        pc = [[c, -s], [s, c]]
    else:
        raise ValueError(f'unsupported WCS for {frame["file"]}: {axes}')
    scale = ARCSEC * np.array([float(h['CDELT1']), float(h['CDELT2'])])
    matrix = np.diag(scale) @ pc
    return epoch, o, matrix


def sightline(h, pixel, observer, matrix):
    # Review images and marks use zero-based pixel centres with origin at
    # lower left. FITS CRPIX is one-based. CRVAL is a spherical angle pair,
    # so invert TAN in its reference-point basis before rotating into HCI.
    crpix = np.array([float(h['CRPIX1']), float(h['CRPIX2'])])
    xi, eta = matrix @ (pixel + 1 - crpix)
    lon, lat = ARCSEC * np.array([float(h['CRVAL1']), float(h['CRVAL2'])])
    centre = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon),
                       math.sin(lat)])
    lon_tangent = np.array([-math.sin(lon), math.cos(lon), 0])
    lat_tangent = np.cross(centre, lon_tangent)
    local = centre + xi * lon_tangent + eta * lat_tangent
    # Helioprojective positive X points solar west. The other two local
    # axes are the projected north pole and the line from observer to Sun.
    radial = unit(observer)
    west = unit(np.cross([0, 0, 1], radial))
    north = np.cross(radial, west)
    return unit(np.column_stack((-radial, west, north)) @ local)


def pack(source, picks):
    event = dict(source['replay'])
    epoch = utc(event['reference_utc'])
    frames = {item['file']: item for item in source['frames']}
    ephemerides = source['ephemerides']
    if len(frames) != len(source['frames']):
        raise ValueError('duplicate source frame')
    for key in ('max_light_time_error_s', 'speed_budget_m_s'):
        if not math.isfinite(source[key]) or source[key] < 0:
            raise ValueError(f'invalid uncertainty term: {key}')
    pointings = {name: pointing(frame, ephemerides) for name, frame in frames.items()}
    tracks, holdout = [], []
    seen = set()
    with picks.open(newline='') as table:
        for pick in csv.DictReader(table):
            frame = frames[pick['frame']]
            h = frame['header']
            t, o, matrix = pointings[pick['frame']]
            pixel = np.array([float(pick['x_px']), float(pick['y_px'])])
            sigma_px = float(pick['sigma_px'])
            dimensions = [int(h['NAXIS1']), int(h['NAXIS2'])]
            key = (frame['file'], *pixel)
            if (key in seen or not np.all((0 <= pixel) & (pixel < dimensions))
                    or not math.isfinite(sigma_px) or sigma_px <= 0):
                raise ValueError(f'invalid or duplicate input mark: {pick}')
            seen.add(key)
            available = utc(frame['available_utc'])
            # Use the centre-of-Sun light time because a surface brightness
            # edge supplies a ray but no independent line-of-sight depth.
            # Its unresolved depth contributes an angular error budget.
            distance = np.linalg.norm(o)
            elapsed = (t - epoch).total_seconds() - distance / C_M_S
            if elapsed < 0 or available < t:
                raise ValueError(f'inconsistent observation or availability time: {frame["file"]}')
            marking = sigma_px * np.linalg.norm(matrix, ord=2)
            depth = source['max_light_time_error_s'] * source['speed_budget_m_s'] / distance
            row = {
                'view': frame['instrument'], 'frame_file': frame['file'],
                'observed_utc': t.isoformat(timespec='microseconds').replace('+00:00', 'Z'),
                'available_utc': frame['available_utc'], 'elapsed_s': round(elapsed, 9),
                'observer_hci_m': [round(float(x), 3) for x in o],
                'los_hci': sightline(h, pixel, o, matrix).round(14).tolist(),
                'sigma_rad': round(math.hypot(marking, depth), 14),
            }
            # Whole exposures are held out together. Neighbouring marks
            # on a single exposure share pointing and feature uncertainty.
            if frame['role'] == 'holdout':
                holdout.append(row)
            elif frame['role'] == 'fit':
                tracks.append(row)
            else:
                raise ValueError(f'frame role must be fit or holdout: {frame["file"]}')
    return {**event, 'tracks': tracks, 'holdout_tracks': holdout}


def sync_archive(source, directory, record):
    # FITS bytes and their header are immutable inputs. A new Horizons
    # response can change its generation date or descriptive text. Numeric
    # UTC and geometric positions still have to match the archived knots.
    for frame in [*source['frames'], *source['backgrounds']]:
        path = directory / frame['file']
        derived = {'sha256': sha256(path), 'header': archive_header(path, frame['instrument'])}
        if record:
            frame.update(derived)
        elif any(frame[key] != value for key, value in derived.items()):
            raise ValueError(f'archive hash or FITS header mismatch: {path}')
    for info in source['ephemerides'].values():
        path = directory / info['file']
        knots = horizons_rows(path)
        if record:
            info.update(sha256=sha256(path), knots=knots)
        elif info['knots'] != knots:
            raise ValueError(f'Horizons position or epoch mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--picks', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    archives = parser.add_mutually_exclusive_group()
    archives.add_argument('--verify-archive-dir', type=Path,
                          help='verify downloaded FITS files and Horizons numeric knots')
    archives.add_argument('--record-archive-dir', type=Path,
                          help='fill WCS, source hashes and ephemeris knots for a new event manifest')
    args = parser.parse_args()
    source = json.loads(args.sources.read_text())
    if args.verify_archive_dir is not None:
        sync_archive(source, args.verify_archive_dir, record=False)
    if args.record_archive_dir is not None:
        sync_archive(source, args.record_archive_dir, record=True)
    packed = pack(source, args.picks)
    if args.record_archive_dir is not None:
        args.sources.write_text(json.dumps(source, indent=2, allow_nan=False) + '\n')
    args.output.write_text(json.dumps(packed, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
