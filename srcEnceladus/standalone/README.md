# Enceladus Stage 0: standalone neutral-plume inference

This directory implements Stage 0 of the Enceladus development plan: fit a
steady, collisionless multi-source H2O plume to the Cassini INMS E3 and E5
density profiles.  It is a separate Python application and neither builds nor
links the AMPS numerical machinery.

See [`model.md`](model.md) for the complete physics, observation operators,
inference method, assumptions, configuration reference, and run procedures.

The delivered reference configuration contains three selectable source cases.
Its default legacy/hybrid case contains all three original neutral components:

- eight discrete jets from Tables 2 and 3 of Tenishev et al. (2014);
- four normal and three high-speed/extra tiger-stripe line sources, represented
  by the 102 locations in the existing 2014 application table; and
- an empirical background `c0 + a (r/R)^alpha`, with `alpha = -2` by default.

The fit is useful for producing a Stage-1 source prescription and for exposing
what E3/E5 can and cannot constrain.  It is **not** an independent validation:
E3/E5 are development/fitting flybys.  Per the project split, E7/E13 are held
for development validation and E14/E17/E18/E21 for final validation; this
application does not fit those holdouts.

## Quick start

The application requires Python 3, NumPy, SciPy, and Pillow.  The tested local
environment is Python 3.8.10, NumPy 1.24.4, SciPy 1.10.1, and Pillow 10.4.0.

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
ENCELADUS_PLUME_FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
ENCELADUS_JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  python3 -m unittest discover -s tests -v

python3 plume_fit.py \
  --config config/tenishev2014_stage0.json \
  --observations /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --metadata /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.metadata.json \
  --output /tmp/e3_stage0_fit.json
```

`--output` must name a new file.  Refusing to overwrite evidence is deliberate.
The equivalent E5 command substitutes `e5` in the observation, metadata, and
output names.  Reference outputs are retained in `results/e3_fit.json` and
`results/e5_fit.json`.

## External observation data

Observation data are deliberately not stored in the AMPS source tree. The
compact, verified release is stored at
`/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/`. It
contains:

- the 79-point E3 and 79-point E5 Teolis sticking-corrected H2O profiles used
  for absolute local-density fitting;
- 518 E3 and 519 E5 PDS L1A closed-source mass-44 samples, retained in count
  space with counter, quality, pointing, state, and lineage fields;
- 97 samples from the 2007 zeta Orionis stellar occultation and 1,133 samples
  from the 2010 E10 solar occultation, expressed as count attenuation,
  optical-depth proxy, tangent geometry, and normalized column shape; and
- two separate published H2O column constraints for the UVIS amplitudes.

The raw products, calibrated UVIS spectra, SPICE kernels, archive checksum
manifests, and reproducibility scripts are also under
`/data/vtenishe/cassini_validation_data`. The release's
`release_manifest.json` records checksums and the exact fit role of each table.

The code never infers an observation directory. Supply the selected CSV and
its metadata explicitly with `plume_fit.py --observations PATH --metadata
PATH`. This makes alternate releases and machines selectable without embedding
a filesystem location in the model. Use `python3 plume_fit.py --help` for the
complete CLI.

The present `plume_fit.py` interface consumes the H2O-density tables. It does
not yet implement the UVIS line-of-sight likelihood or the INMS mass-44 forward
instrument response. Until those pieces are implemented, the UVIS normalized
shapes must not be interpreted as pointwise absolute columns, and mass-44
counts must not be relabeled as CO2 or H2O density.

### E3/E5 and UVIS fitting study

Run the complete independent-fit, E3+E5-combined, UVIS-combined,
all-data-combined, source-family ablation, and cross-prediction study with:

```bash
python3 joint_fitting_study.py \
  --config config/tenishev2014_stage0.json \
  --data-root /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
  --output-dir /data/vtenishe/cassini_validation_data/studies/my_fit \
  --quadrature-order 192
```

The output directory must not already exist. The checked study is stored at
`/data/vtenishe/cassini_validation_data/studies/e3_e5_uvis_fit/`; its
`fitting_study.md` gives the scientific interpretation and its JSON/CSV files
retain every fitted coefficient, per-dataset combined-objective contribution,
diagnostic, input hash, and cross-prediction metric. `fit_profiles.png` and
`.pdf` compare each independent fit with its instrument-combined fit and the
all-data-combined diagnostic.

UVIS tangent tracks come from the SwRI release. For the distant 2007 stellar
event, the line direction is reconstructed as the unit vector normal to both
the tangent radius and its track derivative, matching the legacy SwRI
construction. For the 2010 solar event, `spice_observer_ranges.cpp` supplies
the finite Cassini-to-Sun ray direction and observer coordinate. This matters
after the infinite-line tangent moves behind Cassini: those 308 rows are
integrated outward from the spacecraft rather than through a fictitious chord.
The code fits resolved UVIS shape in optical-depth space with a signed
reference-count offset; it applies the separately published peak/representative
column once to establish the physical scale. Thus the undocumented legacy
scalar calibration is not used.

Rebuild the 2010 finite-observer geometry with:

```bash
make spice-observer-ranges CSPICE_ROOT=/home/vtenishe/SPICE/cspice
cd /data/vtenishe/cassini_validation_data/data/raw/spice/metakernels
/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/build/spice_observer_ranges \
  cas_uvis_2010_e10_solar.tm \
  /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/uvis_2010_e10_solar_optical_depth.csv \
  /tmp/uvis_2010_e10_solar_observer_range.csv
```

The checked geometry and its provenance metadata are already stored with the
external fit release. Generate into a new path and compare checksums before
replacing that evidence product.

## Architecture

| Path | Responsibility |
| --- | --- |
| `plume_model.py` | SI-coordinate transforms, velocity moments, flux normalization, point/line density kernels, local-density matrices, and UVIS line-of-sight column matrices |
| `fitting.py` | configuration validation, source assembly, nonnegative weighted least squares, optional bounded shape optimization, and identifiability diagnostics |
| `plume_fit.py` | strict data/metadata checks, command-line interface, hashes, and JSON result writing |
| `joint_fitting_study.py` | E3/E5 and UVIS independent fits, instrument-separated and all-data combined shape fits, fixed-shape/per-observation rate fits, ablations, cross-predictions, plots, and machine-readable report |
| `config/tenishev2014_stage0.json` | reference H2O physics, three source-case definitions, empirical background, and fit controls |
| `data/published_stripes_2014.json` | seven IAU_ENCELADUS polylines, 102 nodes in metres |
| `/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json` | external 98-jet geometry, uncertainty metadata, and published 41-jet selection |
| `/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/` | external verified E3/E5/UVIS joint-fit tables and release manifest |
| `tools/extract_stripe_geometry.py` | checked extraction of the legacy 2014 C++ geometry table |
| `tools/build_teolis_profiles.py` | reproducible figure digitization and attachment of PDS V2 trajectory geometry |
| `tools/spice_ground_tracks.cpp` | reusable CSPICE computation of Cassini states and ellipsoidal ground tracks |
| `tools/spice_observer_ranges.cpp` | CSPICE finite-observer geometry for the 2010 Cassini-to-Sun UVIS rays |
| `config/e3_e5_spice_ground_tracks.json` | exact E3/E5 event definitions, kernel/library hashes, and output checksums |
| `figure_data/e[35]_spice_ground_track.csv` | 1 s SPICE states and sub-spacecraft coordinates over +/-1600 s |
| `tools/extract_figure3_vector.py` | recover the published Figure 3 curve centerlines from PDF vector paths |
| `tools/plot_reconstructed_figures.py` | polar-map reprojection and deterministic Figure 2/3 rendering |
| `figures/figure2_ground_tracks_sources.*` | reconstructed source/stripe/ground-track map in PNG and PDF |
| `figures/figure3_published_density_comparison.*` | reconstructed published density comparison in PNG and PDF |
| `figures/figure3_current_stage0_fit.*` | paper-style comparison generated from the current Stage-0 fits |
| `tests/` | independent numerical/physics and input-contract tests |

All model calculations use SI units.  Source and observation coordinates are
in `IAU_ENCELADUS`.  Latitudes are signed planetocentric latitudes and west
longitude is positive westward.  The paper tabulates positive magnitudes for
south latitude; those values are therefore stored as negative latitudes in the
configuration.  Jet azimuth is clockwise from local north and tilt is measured
from the outward radial normal.

## Recomputing Cassini ground tracks with CSPICE

`tools/spice_ground_tracks.cpp` is independent of AMPS and is deliberately
flyby-agnostic.  It links directly to the supplied CSPICE N0067 static library;
adding a flyby requires a UTC reference, interval, event PCK, and reconstructed
SPK rather than a new trajectory implementation.

Build it with:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone
make spice-ground-tracks \
  CSPICE_ROOT=/home/vtenishe/SPICE/cspice
```

The command interface is:

```text
build/spice_ground_tracks LSK PCK SPK TAG UTC_REFERENCE T0 T1 DT OUTPUT.csv
```

The checked E3 invocation is:

```bash
build/spice_ground_tracks \
  /home/vtenishe/SPICE/Kernels/NAIF/naif0012.tls \
  /home/vtenishe/SPICE/Kernels/CASSINI/cpck10Mar2008.tpc \
  /home/vtenishe/SPICE/Kernels/CASSINI/080228AP_SCPSE_08050_08076.bsp \
  E3 2008-03-12T19:06:12 -1600 1600 1 /tmp/e3_ground_track.csv
```

For E5, use `cpck16Sep2008.tpc`,
`081010AP_SCPSE_08282_08315.bsp`, and reference time
`2008-10-09T19:06:43`.  The complete machine-readable invocations and SHA-256
kernel identities are in `config/e3_e5_spice_ground_tracks.json`.  The program
uses `spkezr_c` for the geometric (`NONE`) Cassini state and `subpnt_c` with
`Intercept: ellipsoid` for the surface point.  Output includes:

- spacecraft and subpoint Cartesian coordinates in metres;
- observer-to-ellipsoid-intercept distance in metres;
- positive-west longitude in `[0, 360)`; and
- both planetocentric latitude (appropriate for the USGS mosaic) and
  planetographic latitude (used by the legacy Figure 2 code).

The reference UTC is retained from the existing 2014 application and is not
silently relabeled as exact closest approach.  A 0.1 s SPICE scan resolves the
minimum ellipsoid-intercept distance at `-0.5 s` / `52.431 km` for E3 and
`-3.3 s` / `24.960 km` for E5 relative to those references.  Each committed
CSV has 3,201 samples plus its header.  Output files are opened exclusively, so
an existing evidence product cannot be overwritten accidentally.

The event PCKs reproduce the longitude frame used by the legacy 2014 model.
The current USGS controlled mosaic applies the later IAU convention, shifted
3.5 degrees west.  The stored SPICE results remain in their event-PCK frame;
the 3.5-degree correction belongs only in a plot overlay and must be recorded
there.  This separation also leaves the same state/geometry pipeline usable
for later UVIS line-of-sight work, where CK, SCLK, FK, and IK kernels will be
added rather than approximating instrument pointing from a ground track.

## Reconstructed Figures 2 and 3

The requested figures are retained in both raster and vector form:

- `figures/figure2_ground_tracks_sources.png` and `.pdf` show the E3/E5 SPICE
  ground tracks, eight discrete sources, four normal tiger stripes, and three
  extra model fissures over the clean USGS Cassini ISS controlled mosaic.
- `figures/figure3_published_density_comparison.png` and `.pdf` show the E3/E5
  observed and simulated density traces recovered from Figure 3 of Tenishev et
  al. (2014).

The clean map source is the USGS controlled 100 m/pixel high-pass-filtered
global mosaic,
`Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif` (SHA-256
`1ff7fd5526ebf873e62a1daa8d63eef9b99408f7d37cc2d862992393d01a2a2e`), from
<https://astrogeology.usgs.gov/search/map/enceladus_cassini_iss_global_mosaic_hpf_110m>.
The plotter reprojects the south-polar part of that equirectangular mosaic to an
azimuthal-equidistant polar view.  Source, stripe, and event-PCK trajectory
longitudes receive the documented `+3.5 degrees west` shift for alignment with
the current-IAU mosaic; archived numerical coordinates are not altered.

Figure 3 is explicitly a reconstruction of the published curves, not a rerun
of the 2014 simulation and not the new Stage-0 fit.  Its machine-readable data
and extraction metadata are in `figure_data/figure3_published_curves.csv` and
`.metadata.json`.  The extractor converts PDF page 5 to SVG, identifies the
filled vector stroke polygons, and recovers each centerline.  It records the
paper hash, axis calibration, series interpretation, and precision limitation.

Rebuild into an empty output directory with:

```bash
python3 tools/extract_figure3_vector.py \
  /home/vtenishe/Downloads/tiger-stripes.pdf \
  /tmp/figure3_published_curves.csv \
  /tmp/figure3_published_curves.metadata.json

python3 tools/plot_reconstructed_figures.py \
  --figure3-curves /tmp/figure3_published_curves.csv \
  --output-dir /tmp/enceladus-reconstructed-figures
```

Both tools refuse to overwrite existing evidence products.

The separate `figure3_current_stage0_fit.png` and `.pdf` products are generated
from this application's current inference rather than copied from the paper.
Black is the actual fitting observable (the digitized Teolis sticking-corrected
profile), red is the complete fitted model, and blue is an independently
optimized no-stripe fit.  The plotted numerical inputs, hashes, and chi-square
values are recorded in `figure3_current_stage0_fit.metadata.json`.  Rebuild it
into a new `--output-stem` while explicitly selecting the external release:

```bash
python3 tools/plot_current_fits.py \
  --data-root /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
  --output-stem /tmp/figure3_current_stage0_fit
```

The downloaded `CO-S-INMS-3-L1A-U-V2.0` archive cannot replace the black curve
directly: it is an uncalibrated Level-1A count-rate archive and has no ambient
H2O-density field.  Conversion requires species sensitivity, fragmentation,
counter behavior, background removal, ram-angle transmission, and—in the
Enceladus water plume—time-dependent inlet adsorption/desorption.  The old PDS
routine named `inms_compute_density.pro` explicitly calls itself deprecated,
targets Titan, and does not implement H2O.  The V2 errata also says a new
calibration file is still awaiting delivery.  Therefore raw mass-18 rates must
not be labeled observed density.  A reusable later-flyby validation product
will require either a qualified Enceladus inlet-response inversion or a
forward instrument model that compares simulated ambient density in count
space.

## Physical model

### Ballistic velocity law

Each point source emits a drifting Maxwellian

```text
f(v) proportional to exp[-|v - u a|^2 / v_th^2],
v_th = sqrt(2 k_B T / m),     s = u / v_th,
```

where `a` is the unit jet axis, `u` is bulk speed, `T` is temperature, and `m`
is the H2O molecular mass.  For a source-to-observer unit ray `e`, define
`c = e dot a` and

```text
K_p(c,s) = integral_0^infinity t^p exp(-t^2 + 2 s c t - s^2) dt.
```

`K_2` determines number density and `K_3` determines radial particle flux.
They are evaluated analytically with a stable negative-tail expression; the
tests compare them against independent adaptive quadrature.

For the default Tenishev (2010) point-source law, the normalization is

```text
N = integral_outward K_3(e dot a, s) dOmega,
n(x | Q) = Q K_2(e dot a, s) / [N v_th d^2],
F_r(x | Q) = Q K_3(e dot a, s) / [N d^2],
```

where `d` is distance from the source and `Q` is the **total physical source
rate in molecules/s**.  Thus integrating `F_r d^2` over the outward hemisphere
returns `Q`; fitted amplitudes are not arbitrary historical scale factors.
Rays below the local horizon (`e dot n < 0`) are excluded.  Gravity,
collisions, chemistry, and time dependence are not included at this stage.

### Selectable source interpretation

Two emission laws are implemented:

- `tenishev2010_point` (default): the point-source/ray law above;
- `surface_flux_maxwellian`: a locally planar surface-crossing Maxwellian,
  which multiplies both density and normalization integrands by
  `mu = max(e dot n, 0)`.

Select a law by editing the top-level field in a copied input file:

```json
{
  "source_law": "surface_flux_maxwellian"
}
```

The top-level value applies to every jet and stripe.  A source dictionary may
override it with its own `source_law`.  The default remains
`tenishev2010_point`, as requested.  These laws express different physical
definitions of a source and must not be compared using an unnormalized
amplitude; both are normalized here to the same total `Q`.

### Jets, tiger stripes, and background

Each discrete jet is one point source with its published position, axis,
temperature, and bulk speed.  A tiger stripe is a polyline of radially directed
point sources.  Trapezoidal arc-length weights make its configured `Q` a total
rate uniformly distributed per unit fissure length, independent of node
spacing.  The four normal stripes use 180 K and 350 m/s; the three extra
stripes use 180 K and 950 m/s, following the 2014 model.

Three mutually exclusive source cases are selected with the top-level
`source_case` field or the `--source-case` CLI override. The source case fixes
the spatial basis used by the inversion. It does not prescribe the fitted
production rates: every listed jet or stripe component receives its own
nonnegative rate `Q`, which is estimated again for the selected observation.

| `source_case` | Fitted spatial components | Direction and gas kinematics | Total fitted columns with the reference background |
|---|---:|---|---:|
| `tenishev2014_hybrid_8` | 8 jets + 7 stripe components | Original per-source values from the 2014 model | 17 |
| `porco2014_98` | 98 catalog jets; no separate stripe components | Published tilts/azimuths; controlled-comparison `T = 180 K`, `u = 1.5 km/s` | 100 |
| `portyankina2022_solar_41` | 41 selected catalog jets; no separate stripe components | Local-normal axes, `T = 180 K`, `u = 1.5 km/s` | 43 |

The last column includes the constant and inverse-square background amplitudes
enabled by the reference input. Changing `background.model` changes that
count.

#### Case 1: `tenishev2014_hybrid_8`

This is the historical Tenishev et al. (2014) source architecture and the
default case. It contains:

- eight discrete point sources with the published location, direction,
  temperature, and speed of each source;
- four normal line sources on Alexandria, Baghdad, Cairo, and Damascus,
  represented by radially directed nodes at 180 K and 350 m/s; and
- three extra high-speed line components on Baghdad, Cairo, and Damascus at
  180 K and 950 m/s.

The seven line components contain 102 geometry nodes in total. A component
has one fitted total production rate regardless of its number of nodes;
arc-length quadrature distributes that rate uniformly along the component.
The model therefore fits 15 physical source rates plus the configured
backgrounds, rather than fitting 102 independent node rates. Use this case to
reproduce and extend the original tiger-stripe-plus-eight-source model.

This case is self-contained in the source tree: it uses the discrete-source
records in `config/tenishev2014_stage0.json` and the checked stripe geometry in
`data/published_stripes_2014.json`. It does not require the external jet
catalog.

#### Case 2: `porco2014_98`

This case materializes Porco et al. (2014) Table 2 jets 1 through 98. Their
catalog locations, tilts, and azimuths are used, with 26 sources associated
with Cairo, 37 with Baghdad, 29 with Damascus, and 6 with Alexandria. These
are individual point sources; the seven line sources from Case 1 are removed.
The fact that catalog jets lie on the tiger-stripe system must not be confused
with adding a separate continuous stripe component.

Porco et al. constrained particulate-jet geometry, not a unique H2O gas
temperature and bulk speed for every jet. The supplied 180 K and 1.5 km/s
values are therefore explicit controlled-comparison settings borrowed from
the 2010 UVIS source setup. They are assumptions of this fit, not Porco
measurements. Alternative values must be recorded by changing
`source_case_options.porco2014_98` in a copied configuration.

With the reference background this fit has 100 nonnegative columns but each
E3 or E5 profile has only 79 samples. It is necessarily underdetermined even
before accounting for nearby jets with almost indistinguishable signatures.
The fit can constrain combinations of sources and a model-dependent total;
it cannot uniquely measure 98 individual production rates.

#### Case 3: `portyankina2022_solar_41`

This case uses the 41 Porco Table-2 sources corresponding to strictly positive
coefficients in the Portyankina et al. (2022) 2010 solar-occultation solution.
There are 16 Cairo, 9 Baghdad, 11 Damascus, and 5 Alexandria sources. Following
that calculation, every axis is the local outward normal (the catalog tilt and
azimuth are retained only as provenance), `T = 180 K`, and `u = 1.5 km/s`. No
separate continuous stripe component is included.

The 41-jet list is a source *selection*, not a frozen Portyankina source map.
This application independently refits all 41 nonnegative production rates to
E3, E5, or each UVIS observation. The catalog retains the published 2010
coefficient for traceability, but the coefficient is not used as a constraint
or initial fitted amplitude.

The 41 Porco Table-2 IDs are:

```text
1, 4, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 21, 24, 25,
27, 28, 32, 35, 36, 37, 40, 54, 60, 61, 62, 63, 66, 67, 68,
73, 84, 85, 86, 89, 90, 91, 92, 93, 94
```

Portyankina Appendix A labels its table from zero (`jet 00`); the external
catalog performs and records the translation to the one-based Porco Table-2
IDs. The nonzero selection is one possible nonnegative solution after removal
of linearly dependent columns, not evidence that the other 57 vents were
physically inactive.

#### Selecting a case and its data

The Porco/Portyankina catalog deliberately remains outside the source tree at:

```text
/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json
```

The recommended selection mechanism is the command line, which leaves the
reference configuration unchanged:

```bash
python3 plume_fit.py \
  --config config/tenishev2014_stage0.json \
  --source-case porco2014_98 \
  --jet-catalog /data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  --observations /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --metadata /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.metadata.json \
  --output /tmp/e3_porco98_fit.json
```

`jet_catalog_file` in a copied input JSON and the environment variable
`ENCELADUS_JET_CATALOG` are equivalent alternatives to `--jet-catalog`.  The
result records the resolved case, materialized source list, catalog path, and
catalog SHA-256.  The empirical background setting is unchanged across the
three cases; only the surface source architecture changes.

Selection precedence is:

1. `--source-case` overrides the top-level `source_case` in the JSON input;
2. for Cases 2 and 3, `--jet-catalog` overrides `jet_catalog_file`, which in
   turn overrides `ENCELADUS_JET_CATALOG`; and
3. Case 1 ignores the catalog because all of its geometry is local.

An input-only setup is also supported. Copy the reference JSON, set one of the
three exact `source_case` names, and, for Cases 2 or 3, set an absolute catalog
path:

```json
{
  "source_case": "portyankina2022_solar_41",
  "jet_catalog_file": "/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json"
}
```

Do not combine the 98- or 41-jet selections with the legacy line stripes in
the same named case. A combined catalog-plus-continuous-stripe hypothesis
would be a fourth physical model and should receive a new case name and tests.

### Running the fitting procedure for all three cases

There are two fitting entry points:

- `plume_fit.py` fits one qualified local-density profile, currently E3 or E5;
- `joint_fitting_study.py` runs independent E3, E5, 2007 UVIS, and 2010 UVIS
  fits; combined E3+E5 and combined 2007+2010 UVIS fits; an all-data balanced
  fit; cross-predictions, component ablations, and comparison figures.

In both programs the source geometry, direction, temperature, and speed are
fixed by the selected case. For that fixed basis, the code constructs one
density or column-density response per component and solves for nonnegative
component rates and backgrounds. `shape_parameters` is empty in the reference
configuration, so these commands do not optimize jet locations or kinematics.

The joint study produces three distinct fitting modes:

1. **Individual:** four different source vectors, one fitted independently to
   E3, E5, the 2007 UVIS occultation, and the 2010 UVIS occultation.
2. **Combined by instrument:** one source vector shared by E3+E5 and a separate
   source vector shared by the two UVIS occultations.
3. **All-data combined:** one source vector shared by all four observations.

Combined fits minimize the sum of the mean weighted-square residual of each
included dataset. Equivalently, every point uncertainty in a dataset is
multiplied by `sqrt(N_dataset)`. This prevents the 1,133-point 2010 UVIS event
from dominating the 97-point 2007 event only because it contains more rows.
E3 and E5 both have 79 samples, so this balancing leaves their relative weight
unchanged. These are balanced diagnostics rather than formal joint likelihoods
because the correlated calibration systematics and complete instrument
response covariances are unavailable.

#### Fit E3 and E5 independently for every case

Run the following from the standalone directory. Choose a new `RUN_TAG` for
each execution because result files are opened exclusively and are never
overwritten.

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

MODEL_CONFIG=config/tenishev2014_stage0.json
FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1
JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json
RUN_TAG=source_cases_001
FIT_OUTPUT=/data/vtenishe/cassini_validation_data/studies/${RUN_TAG}/inms
mkdir -p "$FIT_OUTPUT"

for SOURCE_CASE in \
  tenishev2014_hybrid_8 \
  porco2014_98 \
  portyankina2022_solar_41
do
  for FLYBY in e3 e5
  do
    python3 plume_fit.py \
      --config "$MODEL_CONFIG" \
      --source-case "$SOURCE_CASE" \
      --jet-catalog "$JET_CATALOG" \
      --observations "$FIT_DATA/${FLYBY}_h2o_sticking_corrected.csv" \
      --metadata "$FIT_DATA/${FLYBY}_h2o_sticking_corrected.metadata.json" \
      --output "$FIT_OUTPUT/${FLYBY}_${SOURCE_CASE}.json"
  done
done
```

This produces six result files. Each contains the resolved source case and
materialized configuration, component labels, nonnegative amplitudes, source
rates in molecules/s, fitted profile, normalized residuals, chi-square,
rank/singular-value diagnostics, warnings, and SHA-256 provenance. Compare
total production and predicted profiles across cases; do not compare an
individual Case-2 or Case-3 jet coefficient as though it were uniquely
identified.

To run just one dataset, use one iteration of the inner command. The
`--jet-catalog` option may be omitted for `tenishev2014_hybrid_8`.

#### Run the E3/E5/UVIS study for every case

Choose a new study root. `joint_fitting_study.py` requires each case output
directory not to exist, preventing accidental replacement of an earlier
study.

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

MODEL_CONFIG=config/tenishev2014_stage0.json
FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1
JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json
RUN_TAG=source_cases_001
STUDY_OUTPUT=/data/vtenishe/cassini_validation_data/studies/${RUN_TAG}/joint
mkdir -p "$STUDY_OUTPUT"

for SOURCE_CASE in \
  tenishev2014_hybrid_8 \
  porco2014_98 \
  portyankina2022_solar_41
do
  python3 joint_fitting_study.py \
    --config "$MODEL_CONFIG" \
    --source-case "$SOURCE_CASE" \
    --jet-catalog "$JET_CATALOG" \
    --data-root "$FIT_DATA" \
    --output-dir "$STUDY_OUTPUT/$SOURCE_CASE" \
    --quadrature-order 192
done
```

Each case directory contains:

- `study_results.json`, with all independent, instrument-combined, and
  all-data-combined fitted vectors, source totals, diagnostics, catalog hash,
  cross-predictions, and quadrature checks;
- `component_amplitudes.csv`, with the independent, E3+E5-combined,
  UVIS-combined, and all-data-combined coefficients;
- `model_comparison.csv`, with full/no-stripe/no-jet/background-only
  ablations;
- `cross_prediction_metrics.csv`;
- `source_rate_shape_strategies.csv`, with one independently fitted total rate
  per observation under each of four spatial-shape inference strategies;
- `fit_profiles.png` and `fit_profiles.pdf`; and
- `fitting_study.md`, the generated interpretation report.

The four source-rate strategies are designed for the orbital/tidal comparison:

1. `independent_each_observation`: every observation independently determines
   both its relative jet/stripe shape and its total source rate.
2. `shared_shape_all_observations`: E3, E5, 2007 UVIS, and 2010 UVIS determine
   one normalized relative source shape; one total rate is then fitted to each
   observation separately.
3. `shared_shape_inms_only`: E3 and E5 determine the normalized shape, which is
   held fixed while every INMS and UVIS observation receives its own total-rate
   fit.
4. `shared_shape_uvis_only`: the two UVIS occultations determine the normalized
   shape, which is held fixed while every observation receives its own
   total-rate fit.

In strategies 2--4, "shape" means only the unit-sum relative amplitudes of the
physical jet and stripe components. The empirical constant and radial
backgrounds are excluded from that normalization and are refitted as nuisance
amplitudes for each observation. Thus the reported rate is always the sum of
jet and stripe production, not plume plus background. UVIS rate fits use the
documented physical-column proxy derived from the independently fitted signed
optical-depth baseline and the published absolute column constraint; they do
not treat raw optical depth as an absolutely calibrated column density.

After all three directories are complete, build the cross-case CSV, JSON,
Markdown report, and PNG/PDF comparison figure with:

```bash
python3 compare_source_case_studies.py \
  --study-root /data/vtenishe/cassini_validation_data/studies/${RUN_TAG}/joint
```

The comparison tool verifies the three case identities, hashes each input
`study_results.json`, refuses to overwrite an existing synthesis, treats a
zero-residual-degree-of-freedom fit as saturated rather than assigning it a
reduced chi-square, and retains both chi-square per sample and rank-adjusted
diagnostics.

For Cases 2 and 3, `no_stripes` is intentionally identical to `full` because
those cases contain no line-stripe components. For Case 1, `no_stripes` tests
removal of all seven line components. In every case, `no_jets` removes all
discrete point sources. These definitions must be considered when comparing
the ablation tables across source cases.

The UVIS portion fits the resolved optical-depth *shape* and applies one
separately published H2O column constraint to set each event's physical scale.
It is not a calibrated point-by-point UVIS column-density inversion. For the
UVIS-combined and all-data-combined fits, each independent signed reference
offset and published peak column first define a nonnegative physical-column
proxy; one source vector is then fitted to both proxy profiles. The JSON stores
the reference offsets, tau-to-column scale factors, clipping counts, and
dataset weights. Each combined-fit object also stores `per_dataset` raw
chi-square, chi-square per sample, and its contribution to the balanced
objective. Independent observation-specific fits remain the appropriate
outputs for measuring how the preferred source map changes between epochs;
combined fits test whether one map can describe the selected epochs.

The empirical background is selectable as `none`, `constant`, `power_law`, or
`both`.  In the reference input,

```text
n_background(r) = c0 + a (r / R)^-2,
R = 252.1 km.
```

The fitted `a` is density at `R`.  The result also reports the equivalent
coefficient `c1 = a R^2` for the paper's `c0 + c1/r^2` notation.

## Observation construction and calibration boundary

The public INMS L1A V2 files contain counts and Cassini target-frame geometry,
not a ready-to-fit sticking-corrected H2O density.  Teolis et al. (2010) used a
1293-element inlet/surface response model, including adsorption and
re-emission.  The published description gives the sticking law and conversion
efficiencies, but the element geometry and complete response operator needed
to reproduce that inversion are not in the material available to this work.
Applying one scalar count-to-density factor would therefore be physically
incorrect.  `plume_fit.py` explicitly rejects raw counts.

The approved alternative is a traceable reconstruction of the published,
sticking-corrected density observable:

1. Digitize the red “Sticking Model” E3/E5 traces reprinted in Figure 5 of
   Tenishev et al. (2010), over outbound `+10` to `+400 s` at 5 s spacing.
2. Normalize each trace to the exact Teolis et al. (2010) reported maximum:
   `7.6e13 m^-3` for E3 and `3.4e14 m^-3` for E5.
3. Interpolate `sc_pos_t_x/y/z` versus `time_ca` from the official PDS
   `CO-S-INMS-3-L1A-U-V2.0` tables.  PDS positions are target-centered
   IAU_ENCELADUS kilometres and are converted to metres.

The independent fitting sigma is the density change corresponding to two
rendered pixels: 5.157% for E3 and 5.174% for E5.  The 23% INMS neutral
calibration uncertainty reported by Teolis et al. (2015) is a correlated common
scale uncertainty.  It is recorded in metadata and intentionally not added to
each point as if independent.  Consequently, result `chi2` values assess fit to
the digitized shape under those pixel weights; they are not complete
instrument-likelihood statistics.

Every profile metadata file records selection, units, frame, response status,
uncertainty semantics, trajectory source, PDF SHA-256, render calibration, and
normalization.  The CLI refuses missing or incompatible metadata.

## Inference method and current fit

For fixed source shapes, density is linear in nonnegative source rates and
background amplitudes.  Columns are weighted by the supplied independent
sigmas and normalized before SciPy bounded least squares.  Source columns
multiply `Q / 1e28`; output `source_rate_s` is restored to molecules/s.  If
explicit `shape_parameters` are supplied, bounded nonlinear least squares wraps
the amplitude solution (variable projection).  No shape parameters are enabled
in the reference fit because E3/E5 alone cannot support them without additional
priors.

| Metric | E3 | E5 |
| --- | ---: | ---: |
| Samples | 79 | 79 |
| Total fitted source rate (molecules/s) | `2.061e28` | `4.520e28` |
| Constant background (m^-3) | `4.091e11` | `5.150e11` |
| Power-law density at R (m^-3) | effectively zero | `7.291e13` |
| chi-square / nominal degrees of freedom | `404.06 / 62` | `83.05 / 62` |
| Nominal reduced chi-square | `6.517` | `1.339` |
| Scaled condition number | `3.53e7` | `3.02e8` |
| Effective rank at relative SVD cutoffs `1e-2 / 1e-3 / 1e-4` | `8 / 10 / 11` | `7 / 9 / 10` |

There are 17 fitted columns (8 jets, 7 stripes, and 2 backgrounds).  Although
both matrices have formal floating-point rank 17, their very large condition
numbers and much smaller effective ranks show that the flight tracks do not
uniquely distinguish individual jets and stripes.  The nonnegative solver
therefore selects one of many sensitive partitions.  **Individual component
rates must not be interpreted as uniquely measured source strengths.**  Even
the total rate is conditional on this geometry, source law, and background.
The result JSON carries this warning and the full singular spectrum so a
downstream stage cannot silently lose the limitation.

E3's reduced residual is also substantially larger than one under the narrow
digitization-only weights.  This is evidence of model/data-shape mismatch or
unrepresented error, not a reason to inflate errors or declare validation.

## Verification and acceptance tests

Run:

```bash
ENCELADUS_PLUME_FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
ENCELADUS_JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  python3 -m unittest discover -s tests -v
python3 -m py_compile plume_model.py fitting.py plume_fit.py \
  tools/extract_stripe_geometry.py tools/build_teolis_profiles.py
```

The tests require all of the following:

- analytic `K_2` and `K_3` agree with independent numerical quadrature;
- integrated outward radial flux equals configured `Q` for both laws;
- density has exact inverse-square scaling along a fixed ray;
- the infinite-chord integral of an `r^-2` background matches its analytic
  value and is invariant to ray-direction reversal;
- below-horizon density is zero, and the surface law has its projected-area
  zero at a tangent ray while the point law does not;
- signed latitude and west-longitude conversion agrees with an independent
  known vector;
- line-source weights integrate uniform production along physical length;
- the reference configuration resolves 8 jets plus 7 stripes/102 nodes;
- all three source cases resolve to exactly 8+7, 98, and 41 spatial components,
  and the 41-source case matches the published 2010 ID selection;
- both source-law selections propagate through all configured sources;
- nonnegative inference recovers an independent synthetic linear fixture and
  reports unobservable columns; and
- raw counts are rejected and digitized E3/E5 maxima match the published
  values exactly; and
- the E3/E5 SPICE subpoint geometry, output convention, and closest-approach
  samples agree with regression fixtures derived from the pinned kernels.

These tests verify the mathematical implementation and internal conservation.
They do not validate the collisionless assumptions, published source geometry,
digitization as an instrument response, or applicability to held-out flybys.

## Rebuilding source data

Stripe geometry is reproducibly extracted from the existing application table:

```bash
python3 tools/extract_stripe_geometry.py \
  ../EnceladusMultiPlume_TigerStripes3D.cpp /tmp/published_stripes_2014.json
```

The extractor fails unless it finds exactly seven named fissures and 102 nodes.
The qualified generated profiles live in
`/data/vtenishe/cassini_validation_data/data/derived/inms_sticking_corrected/`.
To rebuild them into a new staging directory for comparison:

```bash
mkdir -p /tmp/enceladus-stage0-profiles
python3 tools/build_teolis_profiles.py \
  --pdf /home/vtenishe/Downloads/Tenishev-2010-JGR.pdf \
  --inms-root /data/vtenishe/cassini_validation_data/data/raw/inms \
  --output-dir /tmp/enceladus-stage0-profiles
```

The digitizer deliberately checks the rendered page dimensions and refuses a
changed layout, because its axis calibration would then be invalid.

## Known limits and handoff rules

- The point-source solution is singular at a source location and is intended
  for spacecraft points outside the adopted 252.1 km reference sphere.
- Enceladus gravity, source-surface topography, neutral-neutral collisions,
  condensation, photochemistry, electron impact, ionization, pickup ions,
  plasma feedback, dust, and temporal source variability are absent.
- A spherical radial normal is used for every fissure node.  Local terrain
  normals are not available in the adopted source table.
- The steady solution has no exposure-time convolution.  The digitized traces
  are treated as instantaneous published densities.
- The 23% calibration term can rescale source amplitudes coherently and must be
  propagated as a systematic in later comparisons.
- E3/E5 constrain a profile along two trajectories, not a globally unique
  source map.  Stage 1 should carry alternative regularized/source-family
  ensembles rather than treating the NNLS partition as ground truth.
- The held-out flybys must remain untouched until their declared validation
  stage.  Any tuning to E7/E13 or E14/E17/E18/E21 changes their role and must be
  recorded before results are interpreted.

## Primary provenance

- Tenishev et al. (2010), *An approach to numerical simulation of the gas
  distribution in the atmosphere of Enceladus*, JGR, DOI
  `10.1029/2009JA015223` (local source:
  `/home/vtenishe/Downloads/Tenishev-2010-JGR.pdf`).
- Teolis et al. (2010), *Detection and measurement of ice grains and gas
  distribution in the Enceladus plume by Cassini's Ion Neutral Mass
  Spectrometer*, DOI `10.1029/2009JA015192`.
- Tenishev et al. (2014), *Effect of the Tiger Stripes on the water vapor
  distribution in Enceladus' exosphere*, DOI `10.1002/2014JE004700` (local source:
  `/home/vtenishe/Downloads/tiger-stripes.pdf`).
- Porco, DiNino, and Nimmo (2014), *How the Geysers, Tidal Stresses, and
  Thermal Emission across the South Polar Terrain of Enceladus Are Related*,
  DOI `10.1088/0004-6256/148/3/45`, Table 2 source geometry.
- Portyankina et al. (2022), *Modeling the complete set of Cassini's UVIS
  occultation observations of Enceladus' plume*, DOI
  `10.1016/j.icarus.2022.114918`, Appendix A 2010 selection and source setup.
- Teolis et al. (2015), INMS sensitivity/calibration update, DOI
  `10.1007/s11214-014-0133-8`.
- Cassini INMS calibrated archive `CO-S-INMS-3-L1A-U-V2.0`, used here only for
  time-tagged IAU_ENCELADUS trajectory geometry.
