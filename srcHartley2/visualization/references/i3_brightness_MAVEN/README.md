# MAVEN/IUVS 3I brightness references

`maven_iuvs_brightness.csv` contains the H Ly-alpha and OH(0,0) radial
brightness points extracted from the 12 PNG panels supplied in `Archive.zip`.
Distance is in units of `10^5 km`; brightness is in kR.  The table is the
default observation source for `../../plot_i3_brightness_comparison.py`.

The original numerical measurements were not supplied.  These values are a
raster digitization and must not be represented as higher-precision source
data.  `digitization_pixels.csv` preserves every reviewed center-row anchor
and every vertical uncertainty bar that is visibly separable from the marker
and connecting curve.  Blank vertical uncertainty columns mean that no
separate vertical error bar could be resolved in the raster; they do not mean
zero physical uncertainty.

The point centers and horizontal bounds are not guessed from black-pixel
widths.  Each center is `point_index * resolution_1e5_km`, using the spatial
resolution independently transcribed from `srcHartley2/i3/reference/i3_cases.csv`.
The lower and upper bounds are one half resolution on either side, matching
the contiguous horizontal bins visible in every panel.  Brightness centers
and visible vertical uncertainty endpoints are mapped from reviewed pixels
through the detected logarithmic plot spines and the common source-panel
limits (`0.2--20` in `10^5 km`, `0.09--20 kR`).

`metadata.json` records the archive, case-table, anchor-table and output-table
SHA-256 digests; every member PNG digest, image size and detected plot box;
axis limits; units; extraction policy; and row count.  Rebuild and validate
the CSV from the supplied archive with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  srcHartley2/visualization/extract_i3_brightness_maven.py
```

The extractor requires Python 3, NumPy, and Pillow.  Its regression test is:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  srcHartley2/tests/test_i3_brightness_reference.py
```

The extraction ignores the red legacy-model curves.  `Archive.zip` is needed
only to audit or regenerate these committed numerical references; ordinary
comparison plotting reads this directory and does not require the archive.
