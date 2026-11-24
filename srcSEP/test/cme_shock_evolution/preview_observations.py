#!/usr/bin/python3
"""Display archive running differences in the original FITS pixel coordinates."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter

from pack_observations import primary_header


def count_rate(path):
    # These archive products have a 2-D signed 16-bit primary array. FITS
    # stores it big-endian, with BSCALE/BZERO applied before the exposure
    # normalization. A different product needs its own reader and review.
    with path.open('rb') as stream:
        h = primary_header(stream)
        if (h['NAXIS'], h['BITPIX']) != ('2', '16'):
            raise ValueError(f'expected a 16-bit primary image: {path}')
        shape = int(h['NAXIS2']), int(h['NAXIS1'])
        image = np.frombuffer(stream.read(2 * np.prod(shape)), dtype='>i2').reshape(shape)
    return (float(h.get('BSCALE', 1)) * image + float(h.get('BZERO', 0))) / float(h['EXPTIME'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--archive-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--picks', type=Path, help='overlay existing pixel marks')
    parser.add_argument('--save-panel-data', action='store_true',
                        help='write the plotted array, axes and marks beside each PNG as an NPZ')
    args = parser.parse_args()
    source = json.loads(args.sources.read_text())
    picks = []
    if args.picks is not None:
        with args.picks.open(newline='') as stream:
            picks = list(csv.DictReader(stream))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for frame in source['frames']:
        previous = count_rate(args.archive_dir / frame['difference_frame'])
        exposure = count_rate(args.archive_dir / frame['file'])
        style = frame['preview']
        # Use the same detector pixels for both images. Registration,
        # vignetting and the median filter are outside the astrometric
        # calculation. These dimensionless panels serve feature selection.
        difference = median_filter((exposure - previous) / np.maximum(previous, 1),
                                   size=style['median_px'])
        fig, ax = plt.subplots(figsize=(9, 8))
        raster = ax.imshow(difference, origin='lower', cmap='gray', interpolation='nearest',
                           vmin=style['vmin'], vmax=style['vmax'])
        points = [(float(row['x_px']), float(row['y_px'])) for row in picks
                  if row['frame'] == frame['file']]
        offsets = np.empty((0, 2))
        if points:
            marks = ax.scatter(*np.asarray(points).T, facecolors='none', edgecolors='orange', s=40)
            offsets = np.asarray(marks.get_offsets())
        ax.set(xlabel='FITS column, zero based', ylabel='FITS row, zero based',
               title=f'{frame["instrument"]} {frame["file"]} minus {frame["difference_frame"]}')
        fig.savefig(args.output_dir / (frame['file'] + '.png'), dpi=140, bbox_inches='tight')
        if args.save_panel_data:
            # Read the actual image and scatter artists. This records what
            # was plotted and permits a numerical reference independent of
            # fonts, antialiasing or PNG compression.
            np.savez_compressed(args.output_dir / (frame['file'] + '.npz'),
                                difference=np.asarray(raster.get_array()), picks=offsets,
                                extent=raster.get_extent(), origin=raster.origin,
                                color_limits=raster.get_clim())
        plt.close(fig)


if __name__ == '__main__':
    main()
