# Running the srcMoon demonstration runners

This document describes two independent demonstration runners:

1. `environment_demo/run.py` generates illumination, shadowing, and thermal
   state products from one production-kernel run.
2. `lola_surface_demo/run.py` generates and checks the triangulated LOLA
   surface without advancing the illumination or thermal models.

The runners have different purposes and output directories. A successful
execution of one runner does not imply that the other runner, or a formal
I-series validation test, has passed.

## Runner 1 — Illumination, shadowing, and thermal state

### Description

`srcMoon/test/environment_demo/run.py` produces a self-contained lunar-surface
demonstration from the normal `srcMoon` production kernels. It reconstructs a
triangulated lunar surface from the D01 `LDEM_4` LOLA digital elevation model,
evaluates SPICE Sun/Earth geometry, separates geometric night, terrain shadow,
and Earth occultation, and advances the production depth-resolved thermal
columns. The runner then writes complete CSV tables, CEA/Tecplot surface files,
PNG/EPS figures, logs, hashes, and a run manifest.

The Python runner does not contain an alternative illumination or thermal
model. A thin compiled adapter calls the same production geometry, absorbed
solar-flux, and thermal-column routines used by `srcMoon`. Before exporting a
case, it requires local tests U08, U09, and U11 to pass.

This is an `EXPERIMENTAL` meeting-demonstration runner, not a formal validation
test. A successful run does not establish I15, I17, or I20, and its model
output must not be represented as D02 illumination or D03 Diviner observation
data.

### Commands: one-hour demonstration

The single `environment_demo/run.py` invocation below generates all three
products together: illumination, geometric/terrain shadowing, and the lunar
thermal state. Separate physics runs are not required. Static PNG/EPS figures
are created automatically; the later visualization command is needed only for
an additional rotatable view.

Run these commands from a shell. The example uses a 100 km target surface
resolution and advances the thermal model by 3,600 seconds from the explicit
initial condition in `srcMoon/examples/lola_surface.in`.

```sh
cd /data/vtenishe/Moon-test/AMPS

test -r /home/vtenishe/SPICE/cspice/include/SpiceUsr.h
test -d /home/vtenishe/SPICE/Kernels
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL

rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j

RUN_OUTPUT=test_output/srcMoon/environment-demo/one-hour-100km

python3 -B srcMoon/test/environment_demo/run.py \
  --input srcMoon/examples/lola_surface.in \
  --build-dir build \
  --surface-resolution-m 100000 \
  --advance-seconds 3600 \
  --initialization input-initial \
  --output-dir "$RUN_OUTPUT"

python3 -m json.tool "$RUN_OUTPUT/result.json"
```

The output directory must either not exist or be empty. The runner deliberately
refuses to overwrite an earlier result. Choose a new `RUN_OUTPUT` name for each
case instead of deleting a scientifically meaningful run.

The one-hour case is an initial-condition transient. It is useful for checking
the complete data and visualization workflow, but it is not a converged lunar
day/night thermal solution. For an observation-independent periodic thermal
initialization, use the input-controlled production spin-up:

```sh
cd /data/vtenishe/Moon-test/AMPS

RUN_OUTPUT=test_output/srcMoon/environment-demo/periodic-spinup

python3 -B srcMoon/test/environment_demo/run.py \
  --input srcMoon/examples/lola_surface.in \
  --build-dir build \
  --surface-resolution-m 100000 \
  --initialization periodic-spinup \
  --snapshot-utc 2011-06-15T16:50:00.000Z \
  --output-dir "$RUN_OUTPUT"
```

Periodic spin-up can be substantially longer than the one-hour demonstration.
It uses the convergence criteria and cycle definition recorded in the input
file and has no arbitrary cycle cap. The requested snapshot must be at or
after the converged cycle end.

### Expected results and artifacts

For a successful run, the terminal prints a JSON object containing:

- `execution_status: PASS` for artifact generation and integrity checks;
- `local_prerequisites` with U08, U09, and U11 all equal to `PASS`;
- `surface_element_count: 6000`;
- `thermal_profile_row_count: 480000`, corresponding to 80 thermal layers for
  every logical surface cell;
- separate counts for terrain-visible and terrain-occluded cells;
- minimum and maximum modeled surface temperatures;
- I15, I17, and I20 recorded separately as `SKIPPED / NOT_VALIDATED`.

The exact shadow counts and temperatures depend on the selected epoch,
resolution, input parameters, and thermal initialization. They must not be
replaced with expected constants or adjusted to match a previous run.

With the `RUN_OUTPUT` value above, results are written under:

```text
test_output/srcMoon/environment-demo/one-hour-100km/
```

The principal files are:

| File | Contents |
| --- | --- |
| `surface_environment.csv` | One row per logical surface cell: location, terrain normal, local solar time, incidence, night/terrain/eclipse flags, direct-beam fraction, absorbed flux, and surface temperature. |
| `thermal_profiles.csv` | Every finite-volume thermal layer for every logical surface cell, including depth and temperature. |
| `lola_surface.cea` | LOLA-adjusted triangulated surface in an AMPS-loadable format. |
| `lola_surface.dat` | The same triangulated surface for Tecplot. |
| `moon_environment_illumination_shadowing.png` and `.eps` | Static illumination, terrain-shadow, eclipse, and three-dimensional direct-beam summary. |
| `moon_environment_thermal_state.png` and `.eps` | Static surface-temperature, local-time, distribution, and subsurface-profile summary. |
| `application.in` | Exact input snapshot used by the adapter. |
| `stdout.log`, `stderr.log` | Compiled-adapter output. A successful checked run should have an empty `stderr.log`. |
| `visualizer_stdout.log`, `visualizer_stderr.log` | Visualization receipt output and diagnostics. |
| `run_manifest.json` | Source revision, dirty paths, commands, build identity, compiler, roadmap hashes, and external-data hashes. |
| `result.json` | Machine-readable execution result, metrics, artifact hashes, and formal V&V status. |

### Visualization

Static figures are generated automatically. To open a rotatable 3-D view of
the completed case without replacing its original figures, write the custom
rendering into a separate directory:

```sh
cd /data/vtenishe/Moon-test/AMPS

RUN_OUTPUT=test_output/srcMoon/environment-demo/one-hour-100km
VIEW_OUTPUT="$RUN_OUTPUT/custom-temperature-view"

python3 -B srcMoon/tools/visualize_environment.py \
  "$RUN_OUTPUT/surface_environment.csv" \
  "$RUN_OUTPUT/thermal_profiles.csv" \
  --output-dir "$VIEW_OUTPUT" \
  --interactive-3d \
  --color-by surface_temperature_k \
  --selected-view-output "$VIEW_OUTPUT/selected_temperature_view.png"
```

Drag with the mouse to rotate the globe, press `S` to save the selected view,
and press `Q` to close it. Other supported interactive quantities are listed
by:

```sh
python3 -B srcMoon/tools/visualize_environment.py --help
```

### V&V status

All generated plots are explicitly labelled `EXPERIMENTAL MODEL OUTPUT - NOT
V&V EVIDENCE`. Formal comparison with illumination or Diviner observations
must follow the predecessor, convergence, data-readiness, and frozen-selection
requirements in `srcMoon/development_and_validation.md`.

## Runner 2 — Standalone LOLA surface map

### Description

The separate `srcMoon/test/lola_surface_demo/run.py` runner creates a
presentation-oriented surface package without running the illumination or
thermal model. It uses the production Moon input parser, native LOLA decoder,
icosphere generator, radial vertex adjustment, and CEA/Tecplot writers. It
also independently decodes the native D01 elevation selected for every mesh
vertex and records the comparison in `lola_dem_comparison.csv`.

### Commands: 100 km surface demonstration

Starting from the repository root, configure and build the SPICE-enabled Moon
application and run the 100 km surface case as follows:

```sh
cd /data/vtenishe/Moon-test/AMPS

test -r /home/vtenishe/SPICE/cspice/include/SpiceUsr.h
test -d /home/vtenishe/SPICE/Kernels
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
test -r /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL

rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j

LOLA_OUTPUT=test_output/srcMoon/lola-surface-demo/one-hundred-km

python3 -B srcMoon/test/lola_surface_demo/run.py \
  --build-dir build \
  --surface-resolution-m 100000 \
  --output-dir "$LOLA_OUTPUT"

python3 -m json.tool "$LOLA_OUTPUT/result.json"
```

The output directory must not contain an earlier run. The 100 km target is the
maximum requested pre-topography great-circle edge length, not an elevation
sampling interval. In the checked example it produced 10,242 vertices and
20,480 triangles; these counts are consequences of the icosphere refinement
and should be read from each run's `result.json`, not treated as acceptance
constants.

### Expected results and artifacts

The surface-runner output includes:

| File | Contents |
| --- | --- |
| `lola_surface.cea` | AMPS-loadable triangulated surface with radially adjusted vertex coordinates. |
| `lola_surface.dat` | Tecplot `FETRIANGLE` representation of the same surface. |
| `surface_vertices.csv` | Complete generated vertex table. |
| `lola_dem_comparison.csv` | Production radius/elevation and independently decoded D01 value for every vertex. |
| `lola_surface_demonstration.png` and `.eps` | Static surface maps, production-versus-D01 comparison, and residual distribution. |
| `application.in` | Exact derived input used for the run. |
| `verification/U08/result.json` | Required local LOLA geometry verification result. |
| `stdout.log`, `stderr.log`, `run_manifest.json`, `result.json` | Execution logs, provenance, hashes, metrics, and status. |

### Visualization

To inspect the generated geometry interactively, run:

```sh
cd /data/vtenishe/Moon-test/AMPS

LOLA_OUTPUT=test_output/srcMoon/lola-surface-demo/one-hundred-km

python3 -B srcMoon/tools/visualize_surface.py \
  "$LOLA_OUTPUT/lola_surface.cea" \
  --interactive-3d \
  --color-by elevation-m \
  --reference-radius-m 1737400 \
  --vertical-exaggeration 20 \
  --frame-label "MEAN EARTH/POLAR AXIS OF DE421" \
  --selected-view-output "$LOLA_OUTPUT/selected_view.png"
```

The reference radius and frame label above come from the native `LDEM_4.LBL`
metadata and must not be reused for another product without checking its
label. Drag with the mouse to rotate the surface, press `S` to save the chosen
camera view, and press `Q` to close it. The vertical exaggeration affects only
visualization; it does not alter the CEA/Tecplot surface coordinates.

### V&V status

This runner records artifact generation as `execution_status=PASS` only after
U08 and its independent per-vertex comparison pass. It remains an
`EXPERIMENTAL / MEETING_DEMONSTRATION` product: formal I14 is recorded as
`SKIPPED / NOT_VALIDATED` until its linked-production and predecessor gates
are satisfied.
