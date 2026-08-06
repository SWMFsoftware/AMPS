# C19A — GOES EPEAD East–West directional-access validation

C19A is the public-data directional observational validation for the AMPS
geomagnetic cutoff calculator. It compares modeled broad-aperture proton access
with simultaneous eastward- and westward-looking GOES-13 and GOES-15 EPEAD
measurements during the 17 May 2012 SEP/GLE71 event.

The test is intentionally distinct from C17:

- C17 verifies an exact internal charge-sign/velocity-reversal symmetry.
- C19A tests agreement with directional spacecraft observations.

The primary observable is

```text
log10[(background-subtracted physical EAST flux) /
      (background-subtracted physical WEST flux)]
```

for the P4 and P5 proton channels. The model folds an AMPS directional cutoff
map through a uniform elliptical top-hat approximation to each detector field
of view and a common incident spectrum `J(E) proportional to E^-gamma`.

C19A is a broad-aperture validation, not a full detector-response simulation.
It does not claim to reproduce absolute count rates or the complete EPEAD
energy-angle response matrix.

## Directory contents

```text
srcEarth/test/C19/
├── README.md
├── run_C19.py
├── build_goes_reference.py
├── AMPS_PARAM_C19_mode3d.in
├── AMPS_PARAM_C19_gridless.in
├── event_C19_may2012.json
└── tests/
    └── run_self_tests.py
```

Only source, configuration, documentation, and self-test files are committed.
The `data/` directory is created automatically when the observational reference
is generated; downloaded NOAA files, generated reference tables, provenance
records, and user-supplied magnetic-field drivers are runtime products and are
not included in this commit package.

Generated observational products are written under `data/` by
`build_goes_reference.py`. AMPS results are written under
`test_output/C19_goes_epead_ew/` by default.

## Scientific configuration

### Event

The default manifest selects:

```text
analysis interval:   2012-05-17 06:00 UTC through 2012-05-18 06:00 UTC
background interval: 2012-05-16 00:00 UTC through 2012-05-16 12:00 UTC
spacecraft:          GOES-13 and GOES-15
channels:            P4 and P5
```

The prompt event onset is excluded from the default analysis interval to reduce
contamination by interplanetary directional anisotropy.

### EPEAD channel model

The manifest uses the published nominal channel bounds and fields of view:

| Channel | Nominal energy range | Nominal elliptical half angles |
|---|---:|---:|
| P4 | 15–40 MeV | 45° north–south × 25° equatorial |
| P5 | 38–82 MeV | 60° north–south × 30° equatorial |

For the May 2012 event, the fixed telemetry-head mapping is:

| Spacecraft | Telemetry W head | Telemetry E head |
|---|---|---|
| GOES-13 | physical EAST | physical WEST |
| GOES-15 | physical WEST | physical EAST |

The W/E telemetry labels are invariant instrument labels, not guaranteed
physical look directions. The event mapping follows the GOES-13/15 orientation
reported for May 2012 by Rodriguez et al. (2014).

### Model configurations

The runner supports:

```text
solver:      GRIDDED, GRIDLESS, or BOTH
field model: T96, T05, or both
```

`GRIDDED` initializes the empirical field on the standalone Mode3D AMR mesh and
traces through the production interpolation stencil. `GRIDLESS` evaluates the
same empirical model directly along each trajectory.

The user supplies the May 2012 AMPS-format magnetic-field driver. No event
driver is included in this package.

## Requirements

- Python 3.8 or newer.
- An AMPS executable containing:
  - standalone `gridless` and/or `3d` cutoff modes;
  - `DIRECTIONAL_MAP` support;
  - T96 and T05 field evaluation;
  - GEO trajectory-file conversion through SPICE;
  - the Mode3D field-mesh interpolation path for `GRIDDED`.
- `mpirun` available in `PATH`, or another launcher supplied with `--mpirun`.
- `matplotlib` for comparison plots.
- Network access only when `build_goes_reference.py --download` is used.

Install the plotting dependency when it is not already available:

```bash
python3 -m pip install matplotlib
```

The reference builder itself uses only the Python standard library. A separate
per-test `requirements.txt` is intentionally not included because matplotlib is
the only non-standard Python dependency and is commonly shared by the other
validation runners.

## 1. Generate the observational reference

### Automatic NOAA download

From the AMPS codebase root:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download
```

The script downloads the May 2012 five-minute GOES-13 and GOES-15
`epead_p17ew` CSV files from the NOAA/NCEI GOES Space Environment Monitor
archive. It then:

1. reads P4E/P4W and P5E/P5W directional proton fluxes;
2. applies the event-specific telemetry-head-to-physical-direction mapping;
3. calculates a separate pre-event median background for each spacecraft,
   channel, and telemetry head;
4. subtracts those backgrounds;
5. rejects invalid, flagged, nonpositive, or low-signal samples;
6. forms the physical EAST/WEST ratio;
7. writes a compact compressed reference table and provenance record.

Generated files:

```text
srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz
srcEarth/test/C19/data/reference_C19_goes_epead_ew_provenance.json
```

Downloaded source files are cached under:

```text
srcEarth/test/C19/data/cache/
```

### Build from already downloaded NOAA files

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/g13_epead_p17ew_5m_20120501_20120531.csv \
  --goes15-particle /path/to/g15_epead_p17ew_5m_20120501_20120531.csv
```

### Optional ephemeris

Without ephemeris files, the reference uses the nominal GOES-East and GOES-West
geostationary slots recorded in the event manifest. For a publication run, pass
public one-minute ephemeris CSV files when available:

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/goes13_particle.csv \
  --goes15-particle /path/to/goes15_particle.csv \
  --goes13-ephemeris /path/to/goes13_ephemeris.csv \
  --goes15-ephemeris /path/to/goes15_ephemeris.csv
```

The ephemeris parser searches common timestamp, longitude, latitude, and
altitude column names and uses the nearest record within 180 seconds.

### Reference-builder quality control

The default signal gate requires both physical look directions to satisfy:

```text
(raw flux - background) / background >= 3
```

Override it only for a documented sensitivity study:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download \
  --min-signal-to-background 2
```

Every generated reference has a JSON provenance file containing source paths,
source URLs, SHA-256 checksums, direction mapping, background values, processing
settings, and output checksum.

## 2. Supply the magnetic-field driver

C19A does not generate or bundle the May 2012 T96/T05/TS05 driver. Pass it to
the runner with `--driver`.

The expected AMPS format is one UTC timestamp plus 19 numerical values:

```text
# YYYY-MM-DDTHH:MM:SS Bx By Bz Vx Vy Vz Np Temp SYM-H IMFflag SWflag Tilt Pdyn W1 W2 W3 W4 W5 W6
```

No example driver is committed. The required row format is documented above,
and the real driver is supplied explicitly with `--driver`.

The real driver must:

- cover every selected C19 reference epoch;
- have strictly increasing timestamps;
- have a five-minute median cadence;
- contain no gap larger than ten minutes;
- contain all 19 numerical columns on every data row.

AMPS consumes the complete driver through `DRIVER_FILE`. The Python
postprocessor additionally reads `Tilt` to transform the GEO/GSM observation
position into the SM frame used to label the directional cutoff map.

## 3. Validate the package before running AMPS

```bash
python3 srcEarth/test/C19/tests/run_self_tests.py
```

This executes:

- Python syntax compilation;
- reference-builder synthetic parsing/orientation/background test;
- runner directional-map/aperture/metrics/plot test;
- rendering checks for both GRIDLESS and GRIDDED input templates;
- command-line help checks.

The individual self-tests are also available:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --self-test
python3 srcEarth/test/C19/run_C19.py --self-test
```

## 4. Preview commands and generated inputs

A dry run validates the reference and supplied driver, creates every run
directory, renders `AMPS_PARAM_C19.in` and `C19_trajectory.txt`, and prints the
exact AMPS commands without launching AMPS:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile SMOKE \
  --solver GRIDDED \
  --models T96,T05 \
  --driver /path/to/may2012_driver.txt \
  --amps ./amps \
  -np 4 -nt 16 \
  --dry-run
```

Inspect the generated files under:

```text
test_output/C19_goes_epead_ew/
```

## 5. Run C19A

### Recommended first run

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile SMOKE \
  --solver GRIDDED \
  --models T05 \
  --driver /path/to/may2012_driver.txt \
  --amps ./amps \
  -np 4 -nt 16
```

### Routine T96/T05 Mode3D comparison

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --reference srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz \
  --driver /path/to/may2012_driver.txt \
  --amps ./amps \
  -np 4 -nt 16
```

### Compare direct and mesh-interpolated field evaluation

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver BOTH \
  --models T96,T05 \
  --driver /path/to/may2012_driver.txt \
  --amps ./amps \
  -np 4 -nt 16
```

### Enable POSIX-thread field initialization

This option affects only the `GRIDDED` branch:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --driver /path/to/may2012_driver.txt \
  --mode3d-parallel-field-init \
  --amps ./amps \
  -np 4 -nt 16
```

### Alternative MPI launcher

`mpirun` is the default, consistent with the other test runners. Override it as
needed:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --driver /path/to/may2012_driver.txt \
  --mpirun mpiexec_mpt
```

## Profiles

| Profile | Reference cadence | Intended use |
|---|---:|---|
| `SMOKE` | a small fixed subset | software, parser, and geometry check |
| `ROUTINE` | 60 minutes per spacecraft | routine model/data validation |
| `FULL` | all retained five-minute samples | publication and sensitivity studies |

The cadence can be overridden with `--time-step-minutes`; `0` keeps every
reference epoch.

## Directional calculation

For each selected `(UTC, spacecraft, field model, solver)` combination, the
runner:

1. writes a one-sample GEO trajectory file;
2. launches AMPS with `DIRECTIONAL_MAP T`;
3. reads the SM-labeled global directional cutoff map;
4. constructs local physical EAST and WEST detector boresights at the
   spacecraft position;
5. selects map cells inside the nominal elliptical aperture;
6. converts each directional cutoff to a P4/P5 channel transmission assuming
   `J(E) proportional to E^-gamma`;
7. averages transmission over solid angle;
8. forms the modeled EAST/WEST ratio;
9. compares modeled and observed `log10(E/W)`.

The default angular grid is 10° × 10°. Refine it with:

```text
--dir-lon-res-deg
--dir-lat-res-deg
```

The default spectrum is `gamma=3`; change it with `--spectral-index` and archive
the sensitivity.

## Output products

Top-level files include:

```text
C19_commands.json
C19_reference_used.csv
C19_model.csv
C19_comparison.csv
C19_metrics.csv
C19_aperture_samples.csv
C19_result.json
C19_summary.txt
```

Per-run directories contain:

```text
AMPS_PARAM_C19.in
C19_trajectory.txt
C19_amps.log
cutoff_gridless_dir_map_point_0000.dat
or
cutoff_3d_dir_map_loc_000000.dat
```

Plots are generated for each solver/field-model pair:

```text
C19_comparison_<solver>_<model>.png
C19_scatter_<solver>_<model>.png
C19_transmission_<solver>_<model>.png
C19_aperture_diagnostic.png
```

## Acceptance metrics

The runner reports, separately for each solver, field model, and channel:

- valid modeled fraction;
- east/west sign agreement fraction;
- mean log-ratio bias;
- mean absolute log-ratio error;
- log-ratio RMSE;
- Pearson correlation.

Current thresholds are provisional:

```text
valid fraction       >= 0.85
sign agreement       >= 0.90
correlation          >= 0.60
mean absolute error  <= 0.20 in log10(E/W)
RMSE                 <= 0.30 in log10(E/W)
```

By default, these metrics are reported but do not change the process exit code.
Use `--enforce-acceptance` after the thresholds have been reviewed against the
actual reference and model sensitivity studies.

Exit codes:

```text
0  numerical processing completed; and, when requested, acceptance passed
1  observational acceptance failed with --enforce-acceptance
2  missing input, AMPS failure, parse failure, or incomplete numerical result
```

## Recommended validation sequence

1. Run both Python self-tests.
2. Generate the NOAA reference and inspect its provenance and time series.
3. Validate the supplied magnetic driver and run `--dry-run`.
4. Run `SMOKE`, `GRIDDED`, `T05`.
5. Repeat with and without `--mode3d-parallel-field-init`; results should agree.
6. Run `ROUTINE` for T96 and T05.
7. Run `BOTH` to quantify GRIDDED–GRIDLESS differences.
8. Refine the directional grid and Mode3D mesh.
9. Test spectral-index sensitivity.
10. Use exact ephemeris and all five-minute samples for publication results.

## Known limitations

- The P4/P5 response is a uniform elliptical top-hat, not a measured
  energy-angle response matrix.
- No relative calibration correction is currently applied between the two
  EPEAD heads.
- The source is isotropic and has a single power-law spectral index.
- Prompt-onset interplanetary anisotropy is not modeled.
- Nominal GEO positions are used unless ephemeris files are supplied.
- The physical W/E mapping is event-specific and must not be reused blindly for
  another event or yaw state.
- T96/T05 driver uncertainty and field-model uncertainty are not included in the
  observational error bars.

These limitations should appear in any publication using C19A.

## Public data and references

### NOAA data

- NOAA/NCEI GOES 1–15 Space Weather Instruments:
  https://www.ncei.noaa.gov/products/goes-1-15/space-weather-instruments
- GOES Space Environment Monitor direct archive:
  https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/
- GOES 1–15 SEM L1B/L2 dataset metadata:
  https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncei.swx%3Agoes_sem-l1b-l2_ncei

### Instrument and East–West references

1. Rodriguez, J. V., Onsager, T. G., and Mazur, J. E. (2010),
   “The east–west effect in solar proton flux measurements in geostationary
   orbit: A new GOES capability,” *Geophysical Research Letters*, 37, L07109.
   https://doi.org/10.1029/2010GL042531
2. Rodriguez, J. V., Krosschell, J. C., and Green, J. C. (2014),
   “Intercalibration of GOES 8–15 solar proton detectors,” *Space Weather*, 12,
   92–109. https://doi.org/10.1002/2013SW000996
3. Bruno, A. (2017), “Calibration of the GOES-13/15 high-energy proton
   detectors based on PAMELA solar energetic particle observations,”
   *Space Weather*, 15. https://doi.org/10.1002/2017SW001672

## Reproducibility record

For every scientific run, archive:

- AMPS code revision;
- compiler and build configuration;
- command line;
- generated `AMPS_PARAM_C19.in` files;
- supplied driver plus SHA-256;
- reference table and provenance JSON;
- event manifest plus SHA-256;
- MPI launcher and rank count;
- thread count and field-initialization mode;
- solver, field model, mover, mesh, and directional resolution;
- spectral index;
- all raw directional maps and final C19 CSV/JSON/PNG products.
