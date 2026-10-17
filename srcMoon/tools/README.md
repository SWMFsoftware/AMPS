# srcMoon visualization tools

## Surface visualization

`visualize_surface.py` reads either surface representation written by the
production `srcMoon` geometry path:

- CEA long format: node/face counts, Cartesian nodes in metres, then one-based
  triangle connectivity;
- Tecplot point-packed `FETRIANGLE`: Cartesian nodes and stored radius in
  metres, followed by one-based triangle connectivity.

The tool validates declared record counts, finite coordinates, connectivity,
the Tecplot radius column, non-degenerate faces, and outward orientation before
rendering. It does not modify, interpolate, smooth, repair, or resample the
mesh. Every figure is marked `VISUALIZATION ONLY - NOT V&V EVIDENCE`.

From the repository root, visualize the Tecplot output produced by the LOLA
demonstration runner:

```sh
python3 -B srcMoon/tools/visualize_surface.py \
  test_output/srcMoon/lola-surface-demo/meeting/lola_surface.dat \
  --color-by elevation-m \
  --reference-radius-m 1737400 \
  --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \
  --output test_output/srcMoon/lola-surface-demo/meeting/lola_surface_view.png \
  --output test_output/srcMoon/lola-surface-demo/meeting/lola_surface_view.eps
```

The reference radius and frame in this example come from the native
`LDEM_4.LBL`; the mesh formats themselves do not encode either value. Do not
reuse those options for a different dataset without reading its metadata.
Omitting `--color-by elevation-m` uses centre-relative radius in kilometres and
therefore requires no assumed reference sphere. The default output is
`<input-stem>_visualization.png` next to the input.
Existing figures are not replaced unless `--force` is supplied.

The figure contains the requested 3-D view, its antipodal view, an
east-longitude/planetocentric-latitude triangle-centre map, and the complete
vertex radial distribution. `--view-latitude-deg` and
`--view-longitude-deg` change the first camera; `--show-edges` exposes the
triangulation.

### Rotating the surface and saving the selected view

Use the dedicated interactive mode when an X display is available:

```sh
python3 -B srcMoon/tools/visualize_surface.py \
  test_output/srcMoon/lola-surface-demo/meeting/lola_surface.cea \
  --interactive-3d \
  --color-by elevation-m \
  --reference-radius-m 1737400 \
  --vertical-exaggeration 20 \
  --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \
  --selected-view-output test_output/srcMoon/lola-surface-demo/meeting/selected_view.png
```

The command first writes the normal four-panel visualization and then opens a
standalone rotatable globe. Drag the globe with the mouse. Press `S` or click
**Save selected view** to write a clean single-panel figure at the current
elevation and azimuth; press `Q` to close the window. The selected output may
be repeated to save, for example, both PNG and EPS. The JSON receipt is printed
after the window closes and records the final camera angles and the hashes of
any views saved through the button or `S` key.

On a remote host, interactive mode requires an X display (for example an SSH
session created with X forwarding). Without a display, choose the view using
`--view-latitude-deg` and `--view-longitude-deg` and generate only static
figures. Existing normal or selected-view files require `--force`; repeated
saves within one interactive session update that session's selected view.

### Displaying small topographic relief

The physical LOLA relief is only a small fraction of the lunar radius and can
look nearly spherical at global scale. `--vertical-exaggeration F` applies the
following transform to 3-D display coordinates only:

```text
displayed_radius = reference_radius + F * (true_radius - reference_radius)
```

The default `F=1` is the true geometry. Any other positive factor requires an
explicit `--reference-radius-m`; the tool never invents the datum. Vertex
colors, the longitude/latitude map, histogram, geometry checks, source mesh,
and reported physical radius range remain unexaggerated. Static, interactive,
and saved selected-camera figures visibly state `VERTICAL EXAGGERATION xF -
DISPLAY GEOMETRY ONLY`, and the JSON receipt records the factor. This option is
for presentation and cannot be used as model or V&V evidence.

By default every face is rendered. Very fine surfaces can exceed practical
Matplotlib memory. An explicit `--max-render-faces N` selects a deterministic
display subset while topology checks and the histogram still use the complete
mesh. Both the title and JSON receipt disclose that subset, so a reduced
render cannot be mistaken for a changed production surface. `--show` opens an
interactive Matplotlib window in addition to writing requested files.

Dependencies are Python 3, NumPy, and Matplotlib. The script prints a JSON
receipt containing the input SHA-256, format, complete geometry counts and
radius range, rendering choices, and hashes of saved figures.

Run the utility tests with:

```sh
python3 -B srcMoon/tools/test_visualize_surface.py
```

These are tool-infrastructure tests, not a new U/I test ID and not lunar-model
verification or validation.
