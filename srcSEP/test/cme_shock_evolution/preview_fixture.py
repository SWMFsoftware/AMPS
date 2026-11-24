#!/usr/bin/python3
"""Materialize and check small archival FITS cutouts for the preview regression."""

import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path
from statistics import median

import numpy as np
from PIL import Image

from pack_observations import primary_header, sha256


OBS = Path(__file__).resolve().parent / 'observed_2020_11_29'
FIXTURE = OBS / 'preview-cutouts.json'
REFERENCE = OBS / 'preview-reference.json'
# The asymmetric windows surround existing picks and allow both axes and
# the image origin to be tested. Use one detector for each supported archive
# adapter and the display stretch and median window of the real case.
CUTS = [
    ('22793508.fts', 97, 452, 13, 9, 103, 455),
    ('32643517.fts', 360, 540, 13, 9, 365, 544),
    ('A135400.fts', 601, 957, 13, 9, 605, 962),
]
CALIBRATION = ('BSCALE', 'BZERO', 'EXPTIME')


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def extract(path, digest, bounds):
    # Slice the original signed 16-bit samples before applying either the
    # offset or the count-rate conversion. Preserve the original header
    # scale, zero point and exposure values as decimal strings.
    if sha256(path) != digest:
        raise ValueError(f'archive checksum mismatch: {path}')
    x, y, width, height = bounds
    with path.open('rb') as stream:
        h = primary_header(stream)
        shape = int(h['NAXIS2']), int(h['NAXIS1'])
        if h['BITPIX'] != '16' or x + width > shape[1] or y + height > shape[0]:
            raise ValueError(f'unsupported cutout: {path}')
        raw = np.fromfile(stream, dtype='>i2', count=shape[0] * shape[1]).reshape(shape)
    return {'fits_sha256': digest, 'calibration': {k: h[k] for k in CALIBRATION if k in h},
            'pixels': raw[y:y + height, x:x + width].tolist()}


def expected_panel(panel):
    # A scalar, exact-rational oracle for the plotted values. The division
    # and exposure normalization are independent of NumPy and the median
    # window below mirrors at half-pixel boundaries without using SciPy.
    def rates(cut):
        h = cut['calibration']
        scale = Fraction(h.get('BSCALE', '1'))
        zero = Fraction(h.get('BZERO', '0'))
        exptime = Fraction(h['EXPTIME'])
        return [[(scale * v + zero) / exptime for v in row] for row in cut['pixels']]

    pre = rates(panel['background'])
    cur = rates(panel['exposure'])
    difference = [[(after - before) / max(before, 1) for before, after in zip(p, c)]
                  for p, c in zip(pre, cur)]
    height, width = len(pre), len(pre[0])
    half = panel['frame']['preview']['median_px'] // 2

    def reflect(i, n):
        i %= 2 * n
        return min(i, 2 * n - 1 - i)

    result = [[float(median(difference[reflect(y + j, height)][reflect(x + i, width)]
                            for j in range(-half, half + 1) for i in range(-half, half + 1)))
               for x in range(width)] for y in range(height)]
    style = panel['frame']['preview']
    return {'difference': result, 'picks': [panel['pick']], 'origin': 'lower',
            'extent': [-0.5, width - 0.5, -0.5, height - 0.5],
            'color_limits': [style['vmin'], style['vmax']]}


def record(archive, fixture, reference):
    manifest = json.loads((OBS / 'sources.json').read_text())
    images = {item['file']: item for item in manifest['frames'] + manifest['backgrounds']}
    with (OBS / 'picks.csv').open(newline='') as stream:
        picks = {(row['frame'], float(row['x_px']), float(row['y_px']))
                 for row in csv.DictReader(stream)}
    cutouts = []
    for file, x, y, width, height, mark_x, mark_y in CUTS:
        image = images[file]
        previous = image['difference_frame']
        if (file, mark_x, mark_y) not in picks:
            raise ValueError(f'archive pick moved: {file}, {mark_x}, {mark_y}')
        bounds = x, y, width, height
        frame = {k: image[k] for k in ('file', 'difference_frame', 'instrument', 'preview')}
        cutouts.append({'frame': frame, 'archive_origin_px': [x, y],
                        'pick': [mark_x - x, mark_y - y],
                        'exposure': extract(archive / file, image['sha256'], bounds),
                        'background': extract(archive / previous, images[previous]['sha256'], bounds)})
    dump(fixture, cutouts)
    dump(reference, {p['frame']['file']: expected_panel(p) for p in cutouts})


def write_fits(path, cutout):
    # Write the minimal FITS primary image consumed by the production
    # preview CLI. The cutout has its own zero-based grid and no astrometric
    # WCS. Restore the original stored bytes and radiometric calibration.
    pixels = np.asarray(cutout['pixels'], dtype='>i2')
    height, width = pixels.shape
    h = {'SIMPLE': 'T', 'BITPIX': 16, 'NAXIS': 2, 'NAXIS1': width, 'NAXIS2': height,
         **cutout['calibration']}
    cards = ''.join(f'{key:<8}= {value!s:>20}'.ljust(80) for key, value in h.items()) + 'END'.ljust(80)
    data = pixels.tobytes()
    path.write_bytes(cards.encode('ascii') + b' ' * (-len(cards) % 2880)
                     + data + bytes(-len(data) % 2880))


def prepare(directory):
    # Materialize in the case run directory. The fixture contains 3 pairs
    # of subimages of the original files, so routine runs need no network
    # or full-sized images. Use a fresh directory to avoid stale outputs.
    cutouts = json.loads(FIXTURE.read_text())
    directory.mkdir(parents=True)
    dump(directory / 'sources.json', {'frames': [p['frame'] for p in cutouts]})
    with (directory / 'picks.csv').open('w', newline='') as stream:
        table = csv.writer(stream)
        table.writerow(('frame', 'x_px', 'y_px'))
        for panel in cutouts:
            frame = panel['frame']
            write_fits(directory / frame['file'], panel['exposure'])
            write_fits(directory / frame['difference_frame'], panel['background'])
            table.writerow((frame['file'], *panel['pick']))


def check(directory, reference):
    expected = json.loads(reference.read_text())
    for suffix in ('.png', '.npz'):
        names = {file.name for file in directory.glob(f'*{suffix}')}
        if names != {name + suffix for name in expected}:
            raise ValueError(f'unexpected preview {suffix} artifact set: {sorted(names)}')
    for name, reference_panel in expected.items():
        with np.load(directory / (name + '.npz'), allow_pickle=False) as actual:
            # Only the floating-point arithmetic can differ from the
            # exact-rational reference. A 2e-15 absolute tolerance is
            # round-off at these values. Geometry and marks are exact.
            if set(actual.files) != set(reference_panel):
                raise ValueError(f'{name}: plotted data schema differs from reference')
            for field in ('origin', 'extent', 'picks', 'color_limits'):
                np.testing.assert_array_equal(actual[field], reference_panel[field],
                                              err_msg=f'{name}: {field}')
            np.testing.assert_allclose(actual['difference'], reference_panel['difference'],
                                       rtol=0, atol=2e-15, equal_nan=False,
                                       err_msg=f'{name}: filtered radiometric difference')
        # Decode the visual artifact as well. Orange pixels verify that
        # the scatter markers reached the rendered PNG. Exact raster
        # colors, font glyphs and compression can vary across installations.
        with Image.open(directory / (name + '.png')) as preview:
            if preview.format != 'PNG':
                raise ValueError(f'{name}: expected a PNG preview')
            rgb = np.asarray(preview.convert('RGB'))
        red, green, blue = rgb.transpose(2, 0, 1)
        if not np.any((red > 230) & (green > 125) & (green < 210) & (blue < 90)):
            raise ValueError(f'{name}: PNG missing the orange pixel mark')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('prepare').add_argument('directory', type=Path)
    checker = sub.add_parser('check')
    checker.add_argument('directory', type=Path)
    checker.add_argument('--reference', type=Path, default=REFERENCE)
    recorder = sub.add_parser('record')
    recorder.add_argument('archive', type=Path)
    recorder.add_argument('--fixture', type=Path, required=True)
    recorder.add_argument('--reference', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == 'prepare':
            prepare(args.directory)
        elif args.action == 'check':
            check(args.directory, args.reference)
        else:
            record(args.archive, args.fixture, args.reference)
    except (OSError, ValueError, AssertionError, KeyError) as error:
        parser.exit(1, f'CME observation preview: {error}\n')


if __name__ == '__main__':
    main()
