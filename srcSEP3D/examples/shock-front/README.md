# Reduced shock-front production examples

This directory contains two scientifically distinct zero-particle examples:

- `positive_1au.in` / `positive_1au.event` is the runnable production-style
  synthetic case documented below. Its axis reaches the declared 1-AU
  observer as an accepted fast shock.
- `corona_to_1au.in` / `corona_to_1au.event` is the preserved negative
  reference. Its prescribed geometry reaches 1 AU but its observer state is
  `non-forward-inflow`; do not change its parameters or reinterpret it as a
  shock arrival.
- `handoff_smoke.in` is the short native integration case. It is not a 1-AU
  propagation result.

These are reduced-model examples. The 3-D AMPS volume contains the maintained
undisturbed ambient plasma/IMF reference. The separately published surface
contains the prescribed front and its immediate Rankine--Hugoniot limits. No
sheath, ejecta, wake, or downstream CME volume is implemented or claimed, and
none of these results qualifies deferred BG3D-4.

## Positive synthetic case

All dimensional values are SI and all vectors use the heliocentric inertial
(`HCI`) Cartesian frame. The Sun is at the origin. The launch and observer are
on `+Z`, an open polar line of the bundled synthetic axial dipole.

The finite front is a fixed-width self-similar-expansion (SSE) surface with
45-degree half width. Its apex starts at 804,000 km (`1.1557 R_sun`) with a
speed of 100 km/s. A quintic smoothstep changes the speed to 1400 km/s over
1800 s. Thus velocity and acceleration join continuously at both ends; the
maximum acceleration is about 1.354 km/s2, representative of a fast synthetic
low-coronal launch rather than a fitted event. The apex reaches the exact
matched handoff at `20 R_sun = 13,914,000 km` at 10,200 s.

Beyond handoff, the existing quadratic-drag trajectory uses an effective wind
of 400 km/s and `gamma = 5.037350032041981e-12 m-1`. Radius and speed are
continuous at handoff. The coefficient is in the weak-drag fast-CME proxy
regime and gives an axis speed of about 1035.22 km/s at 1 AU. Its final digits
place the exact analytical 1-AU passage at 124,200 s, on a 600-s native time
boundary; this is output synchronization, not observational calibration. The
canonical ambient and the drag-law effective wind are reported separately.

An independent provider evaluation at the observer gives an accepted
`solved-fast-shock` state (`M_f ~= 2.211`, density compression `~= 2.479`) for
the committed inputs. This establishes synthetic model plausibility only. The
magnetic asset is not a magnetogram, the kinematics are not reconstructed from
coronagraph data, and the result is not observational validation.

## Ambient plasma and magnetic field

`assets/synthetic_dipole.harmonics.csv` is repository-local and its SHA-256 is
frozen in each event. The maintained coronal/SWCME composite evaluates PFSS
inside the source surface (`2.5 R_sun`), then uses its Parker continuation.
The reference state at 1 AU is:

- electron density `7.0e6 m-3`;
- electron and proton temperatures `1.4e6 K`;
- no alpha particles and no electron contribution to mass density;
- spherical radial mass-flux density continuation;
- solar rotation `2.86533e-6 rad/s`.

Pressure is the total declared species pressure. The mesh output adds mass
density using exactly `rho = n_e m_p` for this hydrogen case; the general
output code uses the event's alpha abundance and electron-mass switch. It does
not infer density from magnetic-field magnitude and applies no density, speed,
or Mach-number floor.

The simple dipole has an exact exterior equatorial null. The `+Z` case avoids
that null deliberately. It is useful for a reproducible open-field numerical
example but is not a global map of a real Carrington rotation.

## Domain, resolution, and cadence

The deck allocates the complete Sun-centred domain through `1.55e11 m`, not a
Parker-tube pruning mask. The full finite front and the 1-AU observer therefore
remain inside the represented domain. Outer cells are 20 Gm; through the first
five solar radii the requested cell scale is 1.25 Gm (about 1.8 solar radii),
which the cubic hierarchy realizes as approximately 1.21 Gm. No
heliospheric transport tube is refined because particles are disabled and the
front is represented by its own surface quadrature. The fine spherical
transition deliberately ends below handoff: the volume is an ambient
reference, while front geometry and RH limits come from the independent
surface quadrature. Refining ambient-only cells throughout a 20-solar-radius
sphere or at this scale all the way to 1 AU would add cost without
constructing a sheath. Blocks contain `4^3` cells and at most six AMR
refinements. The declared memory budget is 16 GiB.

The host step and background cadence are 600 s. Regular volume/surface output
is every 60 steps (ten hours). The application additionally writes committed
snapshots near:

- launch, pulse midpoint, and pulse end;
- one step before, at, and one step after the exact 10,200-s handoff;
- one step before, at, and one step after the exact 124,200-s 1-AU passage.

Change `[output].cadence_steps` to change regular output. Landmark snapshots
are retained because they diagnose phase continuity. Change launch/history,
width, ambient, handoff, or endpoint physics in `positive_1au.event`; every
event key and magnetic-asset checksum contributes to the event fingerprint.
Change domain, AMR, cadence, and output paths in `positive_1au.in`; these enter
the application configuration fingerprint. A changed event must be rechecked
for fast-mode admissibility rather than made to pass with a floor.

## Build and run

From the AMPS root, the required clean native build is:

```bash
pwd
ps -C make -C gmake -C g++ -C gcc -C cc1plus -C mpiexec -C mpirun -C amps -o pid=,stat=,cmd=
rm -rf -- build
./Config.pl -application=sep3d
./ampsConfig.pl -input sep3d.input -no-compile
make -C srcSEP3D prepare-production
make -j16 amps
```

The one-line four-rank production run is:

```bash
mpiexec -n 4 ./amps --input srcSEP3D/examples/shock-front/positive_1au.in --output-dir test_output/reduced-front/positive-1au-run
```

Use a fresh `--output-dir`: propagation telemetry is deliberately
non-overwriting. The initialization files keep the paths declared in the deck;
change those three paths as well when retaining multiple initialization
products. No test-registry interpretation or test-only option is required.

For a parser/provider preflight without AMPS allocation:

```bash
./amps --input srcSEP3D/examples/shock-front/positive_1au.in --dry-run
```

The dependency-light example gate is:

```bash
make -C srcSEP3D -j16 test/stage1
(cd srcSEP3D && ./test/stage1 --test RSHAPP04)
```

The complete native qualification (one rank, four ranks, and a second
four-rank run at twice the output frequency) is a single independent command:

```bash
python3 srcSEP3D/test/validate_positive_shock_example.py \
  --output-root test_output/reduced-front/positive-1au-qualification
```

The runner creates relocated decks so the three runs cannot overwrite one
another, invokes the same public `amps --input ... --output-dir ...` CLI shown
above, and writes `summary.txt` and `summary.json`.  It reads an actual valid
native volume row and every front node, in addition to auditing receipts,
epochs, owner/ghost hashes, particle counts and cadence.  A failed check is
named in both summaries and points to its exact `execution.log`; the runner
does not reuse old receipts or turn a missing native run into PASS.

## Products and physical meaning

For each scheduled tick `NNNNNNNN`, the production output directory contains:

- `positive-1au-tick-NNNNNNNN-ambient.dat`: maintained distributed AMPS
  FEBRICK/tetrahedron volume output. Application columns include the actual
  installed `B`, `U`, electron number density, mass density, pressure,
  temperature, time, ambient generation, and front generation. Validity flags
  distinguish physical cells from interpolation padding.
- `positive-1au-tick-NNNNNNNN-front.dat`: FEQUADRILATERAL visualization mesh
  over the provider's supported equal-area surface samples. It includes
  position, outward normal, normal speed, typed status code, accepted flag,
  fast Mach number, `theta_Bn`, density/magnetic compression, and immediate
  upstream/downstream `rho,p,U,B`. The Tecplot encoding is deliberately
  finite. When `shock_accepted=0`, it writes `fast_mach=0`, both compression
  ratios as one, and copies the upstream ambient primitive into the downstream
  display columns. `downstream_valid` remains zero: these copied values are
  not an RH state or downstream CME plasma. `theta_Bn_valid` and
  `magnetic_compression_valid` distinguish finite diagnostic placeholders from
  physical values in magnetic-null cases.
- `...-front.json`: the provider's area/count/identity ledger.
- `...-receipt.json`: committed time/generations, event identity, owner and
  received-ghost readback results, decomposition-independent owner hashes,
  exact observer-arrival status, and actual particle/injection counts.

`shock-history.csv` records every native boundary. `shock_active=1` means the
apex is an accepted fast shock; radius is still recorded when the prescribed
geometric front exists but is not a shock. `native-runtime.json` records the
completed runtime. Surface status codes are the declaration order in
`shock_front/provider.h`; use `shock_accepted`, not a numeric-code guess, for
plot selection.  A surface point has `shock_accepted=1` only when the shared
provider classified it `solved-fast-shock` and supplied a complete valid RH
downstream state.  Every non-forward, sub-fast, unavailable, or unresolved
point has `shock_accepted=0`, irrespective of its finite plotting placeholders.

The volume files are ambient reference values even geometrically behind the
surface. The surface's downstream columns are infinitesimal RH limits only.
They are such limits only where both `shock_accepted=1` and
`downstream_valid=1`; elsewhere they are the documented ambient visualization
fallback. Do not interpolate those values through the volume or label volume
cells as downstream CME plasma.

## VisIt and Tecplot recipe

In VisIt, open one `*-ambient.dat` as Tecplot, make a pseudocolor plot of
`number_density_m-3`, `pressure_Pa`, or `temperature_K`, and add vector plots
from `(B_x_T,B_y_T,B_z_T)` or `(U_x_m_per_s,U_y_m_per_s,U_z_m_per_s)`. Open the
matching `*-front.dat` as a second database, add a Mesh plot plus pseudocolor
of `fast_mach` or `shock_accepted`, and use the same HCI axes. Create a time
series by grouping filenames in tick order.

In Tecplot 360, load the ambient and matching front file into one frame. Use
the ambient FEBRICK zones for slices/isosurfaces and show the front
FEQUADRILATERAL zone as a translucent mesh. Apply a value blanking condition
`shock_accepted < 0.5` when only accepted patches are desired. Do not blank on
`fast_mach` alone because non-forward and numerical states have separate
meaning.

### Standalone Python front viewer

If VisIt cannot open a standalone `*-front.dat` surface, use the included
Python viewer.  It requires Python 3, NumPy and Matplotlib (the same packages
used by the maintained plotting validations):

```bash
python3 srcSEP3D/examples/shock-front/view_front.py --help
```

The built-in help contains copyable examples for discovering and selecting
variables, choosing sequential or diverging color maps, fixing color limits,
filtering accepted shocks, adding the Sun/ruler, saving headless images,
selecting a camera, and encoding GIF/MP4 movies.

```bash
python3 srcSEP3D/examples/shock-front/view_front.py \
  test_output/reduced-front/positive-1au-finite-tecplot-20261006-final/rank-4-cadence-60/products/positive-1au-tick-00000207-front.dat \
  --variable fast_mach
```

The window supports normal Matplotlib 3-D rotation and zoom.  Toggle the
`accepted shock only` checkbox to hide non-shock surface areas.  The equivalent
command-line initial state is:

```bash
python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \
  --variable density_compression --accepted-only
```

Show the Sun and a heliocentric distance ruler through the middle of the front
with:

```bash
python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \
  --variable fast_mach --accepted-only \
  --show-sun --show-distance-axis
```

The Sun is a one-solar-radius sphere at the HCI origin.  The viewer infers the
front centerline from its area-weighted surface centroid; for the equal-area
SSE output this converges to the analytical symmetry axis.  The ruler follows
that Sun-to-front direction and its labels are heliocentric solar radii even
when the Cartesian axes use metres or AU.  This inferred axis is display-only:
it does not replace or feed back into provider geometry, normals, or shock
classification.

The filter uses `shock_accepted`, not Mach or the integer status code.  A quad
is accepted only when all four vertices have `shock_accepted=1`.  In the full
view, quads with four non-shock vertices use their finite output values and
mixed classification-boundary quads are neutral gray.  In accepted-only mode,
both kinds are hidden.  This conservative rule avoids inventing a physical
value by averaging across the accepted/non-shock boundary.

List every plottable column with:

```bash
python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat --list-variables
```

For a compute node without a display, save an image without opening a window:

```bash
python3 srcSEP3D/examples/shock-front/view_front.py FRONT.dat \
  --variable theta_Bn_rad --accepted-only \
  --save front-theta-bn.png --no-show
```

Axes default to solar radii and remain in the file's HCI frame.  Select metres
or AU with `--length-unit m` or `--length-unit au`; use `--vmin`, `--vmax`,
`--cmap`, `--elev`, and `--azim` to control the presentation.  Node values are
averaged over homogeneous quads for visualization only.  They remain
front-local samples and must not be interpreted as a volumetric downstream
sheath or ejecta reconstruction.

### Movies and a fixed user-selected viewpoint

Pass a directory, multiple files, or a quoted shell pattern and select a GIF
destination.  For example, this encodes every surface in one production run:

```bash
python3 srcSEP3D/examples/shock-front/view_front.py \
  test_output/reduced-front/positive-1au-finite-tecplot-20261006-final/rank-4-cadence-60/products \
  --variable fast_mach --accepted-only \
  --show-sun --show-distance-axis \
  --elev 24 --azim -58 --fps 3 \
  --movie positive-1au-front.gif
```

Movie files are ordered by the embedded `time_s`, never merely by filename.
The camera elevation/azimuth, equal-aspect HCI bounding cube, radial-ruler
extent, and color normalization are computed once and then frozen for every
frame.  Consequently, neither viewpoint, zoom, ruler ticks, nor color scale
can jump as the CME front propagates.  The global low-corona-to-1-AU scale
necessarily makes the earliest front appear small; this is the physical size
ratio rather than per-frame auto-zoom.

To choose a camera interactively, first open any representative snapshot,
rotate it to the desired view, and press `v`.  The terminal prints reusable
arguments such as:

```text
selected_view=--elev 24.5 --azim -61
```

Copy those `--elev` and `--azim` values into the movie command.  This keeps the
selected HCI viewpoint identical in every encoded frame.  `--fps` controls
playback speed and `--dpi` controls raster resolution.  GIF output uses the
installed Pillow writer.  An `.mp4` destination is also supported when the
host has the external `ffmpeg` executable; otherwise the viewer fails with an
explicit message and recommends GIF rather than silently changing format.

The parser, classification boundary, Sun/ruler render, fixed-view time sorting,
headless PNG path, and GIF encoder are tested with:

```bash
python3 srcSEP3D/test/test_shock_front_viewer.py
```

## Validation and resource expectations

Qualification requires both 1- and 4-rank runs to agree on event identity,
committed epochs, front geometry/classification, owner fingerprints, and exact
observer passage, with zero particles/injections in every receipt. The exact
provider root is bracketed by the 600-s native time boundaries. Output-cadence
convergence compares otherwise identical four-rank runs writing every 60 and
30 steps; output scheduling must not change the analytical root, accepted
state, native owner fingerprint, or surface product. Actual measured wall time,
block/cell counts, output volume, command lines, and receipts are recorded in
`CODEX_REDUCED_SHOCK_PLAN.md` after execution.

The qualified full-domain mesh allocated 736 leaf blocks and 106,208 physical
cells.  On the qualification host, AMPS reported a peak of about 278 MB for
the one-rank run and 669 MB total (about 167 MB/rank) for the four-rank run.
One-/four-rank default-cadence campaigns took 444/145 s and occupied 2.3/2.4
GiB; an arrival ambient file was 167,920,092/176,709,641 bytes respectively.
The doubled-output-cadence run took 151 s and occupied 2.9 GiB.  These are
measurements, not portable limits: plan for four MPI ranks, the deck's
conservative 16-GiB memory budget, at least 3 GiB per default run, and about 8
GiB for the complete three-run qualification.  A one-rank run is supported
for agreement evidence but is slower and cannot exercise received MPI ghosts.

## Limitations

- The front trajectory and width are prescribed. Drag is a kinematic proxy,
  not a solved CME momentum equation or eruption-instability prediction.
- Local RH states do not supply a spatial downstream sheath/ejecta solution.
- The ambient is a time-independent reference in this selected event; epochs
  and native generations still evolve coherently.
- Independent surface patches do not construct a solenoidal downstream
  volume, and the axial-dipole field is synthetic.
- A positive synthetic shock is not evidence that a particular observed CME
  would retain a shock at 1 AU.
- Particle injection, particle transport, srcSEP coupling, and BG3D-4 remain
  outside this example and are unchanged/deferred.
