# C19A — GOES EPEAD East–West directional-access validation

## 1. Purpose

C19A is the first observational validation of the **directional** geomagnetic-access capability of the AMPS cutoff calculator. Internal symmetry tests such as C17 can verify charge-sign and velocity-reversal consistency, but they do not establish that the calculated directional access agrees with measurements. C19A compares AMPS directional cutoff maps with simultaneous eastward- and westward-looking proton measurements from the GOES-13 and GOES-15 EPEAD instruments.

The implemented public-data event is the 17 May 2012 SEP/GLE71 event. The default analysis interval is the event decay from 2012-05-17 06:00 UTC through 2012-05-18 06:00 UTC. The prompt onset is excluded because interplanetary beam anisotropy can imitate or obscure a geomagnetic East–West effect.

The primary observable is

```text
log10[(physical EAST background-subtracted flux) /
      (physical WEST background-subtracted flux)]
```

for EPEAD P4 and P5. A negative value means that the eastward-looking detector measured less flux than the westward-looking detector.

C19A is a **broad-aperture validation**, not a complete detector simulation. It uses nominal energy intervals, a uniform elliptical top-hat angular response, and a common isotropic power-law incident spectrum. It is intended to test the sign, energy dependence, and temporal variability of directional geomagnetic shielding.

## 2. Implemented comparison

| Item | C19A implementation |
|---|---|
| Spacecraft | GOES-13 and GOES-15 |
| Event | 17 May 2012 SEP/GLE71 decay |
| Observation cadence | Public NOAA/NCEI 5-minute EPEAD product |
| Primary channels | P4: 15–40 MeV; P5: 38–82 MeV |
| P4 nominal FOV | ±45° north–south, ±25° equatorial |
| P5 nominal FOV | ±60° north–south, ±30° equatorial |
| Observation | Background-subtracted physical East/West flux ratio |
| Field models | IGRF + T96 and IGRF + T05/TS05 |
| Solvers | GRIDDED, GRIDLESS, or BOTH |
| AMPS product | Global SM directional cutoff map at each GOES position and epoch |
| Source spectrum | `J(E) ∝ E^-gamma`, default `gamma=3` |
| Instrument model | Uniform elliptical top-hat response inside nominal FOV |

### Event-specific detector orientation

The NOAA files use invariant telemetry labels `E` and `W`; these labels are not themselves guaranteed to identify physical east and west because GOES-13–15 can change yaw orientation. The May 2012 mapping follows Rodriguez et al. (2014), Table 2:

| Spacecraft | Telemetry W | Telemetry E |
|---|---|---|
| GOES-13 | physical EAST | physical WEST |
| GOES-15 | physical WEST | physical EAST |

The mapping is stored in `event_C19_may2012.json`. It must be changed when C19A is extended to a different event.

## 3. Directory contents

```text
srcEarth/test/C19/
├── README.md
├── AMPS_PARAM_C19_gridless.in
├── AMPS_PARAM_C19_mode3d.in
├── build_goes_reference.py
├── event_C19_may2012.json
├── references.bib
├── requirements.txt
├── run_C19.py
├── data/
│   ├── README.md
│   └── reference_C19_goes_epead_ew_schema.csv
└── tools/
    └── prepare_official_ts05_driver.py
```

Generated data and result files are intentionally not included in the source archive. They are reconstructed from the public source products with recorded SHA-256 provenance.

## 4. Requirements

### AMPS executable

C19A requires an AMPS executable with:

- standalone Earth `gridless` and/or `3d` cutoff calculations;
- `DIRECTIONAL_MAP T` output;
- T96 and T05 field models;
- SPICE enabled and the required kernels loaded;
- trajectory input in `TRAJ_FRAME GEO`;
- MPI; and
- POSIX-thread field initialization only when `--mode3d-parallel-field-init` is requested.

SPICE is required because each GOES latitude/longitude/altitude sample is transformed from ITRF93/GEO to GSM and because the directional-map labels are defined in SM and rotated to GSM at the selected epoch.

### Python

Python 3.8 or newer is supported. Install the plotting dependency from the repository root:

```bash
python3 -m pip install -r srcEarth/test/C19/requirements.txt
```

The data-download and TS05-driver scripts otherwise use only the Python standard library.

### MPI launcher

The runner uses `mpirun` by default, as do the other AMPS test runners. Use `--mpirun <launcher>` only when a different command is required by the local system.

## 5. Obtain and build the observational reference

### 5.1 Public NOAA files

The builder uses these fixed public monthly files from the historical NOAA/SWPC operational-average tree. These exact files preserve the directional `p17ew` E/W-head schema used by C19A. NCEI also publishes a newer Version 1 collection for general GOES 1–15 use; replacing the fixed C19A sources requires revalidating variable definitions, detector-head mapping, quality flags, and the resulting reference.

```text
GOES-13:
https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/2012/05/goes13/csv/g13_epead_p17ew_5m_20120501_20120531.csv

GOES-15:
https://www.ncei.noaa.gov/data/goes-space-environment-monitor/access/avg/2012/05/goes15/csv/g15_epead_p17ew_5m_20120501_20120531.csv
```

From the AMPS repository root, download both files and generate the compact reference table:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download
```

This writes:

```text
srcEarth/test/C19/data/cache/g13_epead_p17ew_5m_20120501_20120531.csv
srcEarth/test/C19/data/cache/g15_epead_p17ew_5m_20120501_20120531.csv
srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz
srcEarth/test/C19/data/reference_C19_goes_epead_ew_provenance.json
```

The provenance file records the source URLs and SHA-256 values, manifest checksum, background method, directional mapping, and generated-reference checksum.

To use files downloaded separately:

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/g13_epead_p17ew_5m_20120501_20120531.csv \
  --goes15-particle /path/to/g15_epead_p17ew_5m_20120501_20120531.csv
```

### 5.2 Reference construction

For every spacecraft, channel, and telemetry head, the script calculates the median valid flux over the manifest background interval:

```text
2012-05-16 00:00–12:00 UTC
```

It retains a sample only when:

- both physical directions have finite positive raw flux;
- available quality flags are zero;
- both background-subtracted fluxes are positive; and
- `(raw-background)/background` is at least 3 for both directions by default.

The threshold can be changed explicitly:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download \
  --min-signal-to-background 2.0
```

Changing this value changes the observational reference and is recorded in the provenance file.

### 5.3 Spacecraft position

The routine reference uses the nominal GOES operational slots from the event manifest:

```text
GOES-13: 75° W, geostationary altitude
GOES-15: 135° W, geostationary altitude
```

NOAA also lists public one-minute ephemeris products for GOES 6–15. For publication runs, download the appropriate May 2012 ephemeris CSV files through the NOAA GOES 1–15 Space Weather Instruments page and pass them to the builder:

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/g13_particle.csv \
  --goes15-particle /path/to/g15_particle.csv \
  --goes13-ephemeris /path/to/goes13_ephemeris.csv \
  --goes15-ephemeris /path/to/goes15_ephemeris.csv
```

The ephemeris parser recognizes common UTC, longitude, latitude, and altitude column names. It uses the nearest position within 180 seconds and otherwise falls back to the nominal slot. Every output row records `position_source`.

### 5.4 Reference-builder self-test

```bash
python3 srcEarth/test/C19/build_goes_reference.py --self-test
```

This uses synthetic NOAA-format files to test header discovery, P4/P5 parsing, quality/background filtering, the event-specific E/W mapping, gzip output, and provenance generation. It does not contact NOAA.

## 6. Obtain and prepare the T05/TS05 driver

The official Tsyganenko archive provides yearly five-minute input files containing IMF, solar-wind, SYM-H, dipole tilt, pressure, and W1–W6. Build the C19A driver with:

```bash
python3 srcEarth/test/C19/tools/prepare_official_ts05_driver.py
```

The default source is:

```text
https://geo.phys.spbu.ru/~tsyganenko/TS05_data_and_stuff/2012_OMNI_5m_with_TS05_variables.zip
```

The default selected interval is:

```text
2012-05-16 00:00 UTC through 2012-05-18 06:00 UTC, inclusive
```

The output is:

```text
srcEarth/test/C19/data/ts05_driver_may2012.txt
```

The script copies all physical values and quality flags without interpolation or smoothing and records source URL, archive member, source SHA-256, and selected coverage in comments.

Use an archive downloaded locally:

```bash
python3 srcEarth/test/C19/tools/prepare_official_ts05_driver.py \
  --source /path/to/2012_OMNI_5m_with_TS05_variables.zip
```

Driver-script self-test:

```bash
python3 srcEarth/test/C19/tools/prepare_official_ts05_driver.py --self-test
```

Before AMPS is launched, `run_C19.py` independently verifies the timestamp-plus-19-value schema, strict time ordering, five-minute median cadence, maximum ten-minute internal gap, finite values, and coverage of all selected observation epochs.

## 7. Verify the package before a production run

Run all non-network self-tests:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --self-test
python3 srcEarth/test/C19/tools/prepare_official_ts05_driver.py --self-test
python3 srcEarth/test/C19/run_C19.py --self-test
```

The runner self-test checks directional-map parsing, detector-aperture selection, energy folding, E/W sign, template rendering for both solvers, driver validation, CSV output, and plot generation.

After generating the public reference and driver, preview the exact AMPS commands and rendered inputs:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile SMOKE \
  --solver GRIDDED \
  --models T96,T05 \
  --dry-run \
  --amps ./amps \
  -np 4 -nt 16
```

## 8. Run C19A

### 8.1 Quick smoke run

The SMOKE profile selects the first, middle, and last retained observation epoch for each spacecraft:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile SMOKE \
  --solver GRIDDED \
  --models T96,T05 \
  --amps ./amps \
  -np 4 -nt 16
```

### 8.2 Routine regression

The ROUTINE profile samples the generated five-minute reference at 60-minute spacing:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --amps ./amps \
  -np 4 -nt 16
```

The default reference and driver paths are used automatically. Explicit equivalents are:

```bash
--reference srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz
--driver srcEarth/test/C19/data/ts05_driver_may2012.txt
```

### 8.3 Parallel Mode3D field initialization

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --mode3d-parallel-field-init \
  --amps ./amps \
  -np 4 -nt 16
```

The same `-nt` value controls the Mode3D cutoff thread backend and the number of temporary POSIX workers requested for background-field initialization. The caller also participates in the field initialization, following the implementation documented by C18.

### 8.4 GRIDLESS/GRIDDED cross-solver run

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver BOTH \
  --models T96,T05 \
  --amps ./amps \
  -np 4 -nt 16
```

### 8.5 Full five-minute comparison

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile FULL \
  --solver GRIDDED \
  --models T96,T05 \
  --amps ./amps \
  -np 4 -nt 16
```

This can require many AMPS launches: one launch per selected `(epoch, spacecraft, solver, field model)` group. P4 and P5 share the same directional map for a given group and are folded in postprocessing.

### 8.6 Custom cadence or time interval

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile FULL \
  --time-step-minutes 15 \
  --start 2012-05-17T08:00:00Z \
  --end   2012-05-17T20:00:00Z \
  --solver GRIDDED \
  --models T05 \
  --amps ./amps
```

### 8.7 Direction-mapping diagnostic

The AMPS cutoff implementation and the GOES detector geometry use two different
vector meanings that must not be conflated:

- the AMPS directional-map vector is the **incoming particle arrival/velocity direction** at the observation point; the cutoff code backtraces by launching the reversed vector, and its vertical-arrival direction points toward Earth;
- the GOES EPEAD EAST/WEST direction is a telescope **look direction**. A telescope looking toward `+d` receives particles whose forward-time velocity is approximately `-d`.

Therefore the production C19 mapping is fixed as:

```text
detector look direction = - AMPS directional-map arrival direction
```

This is implemented in one audited helper in `run_C19.py`; it is not a fit parameter.
The runner also reproduces the **old direct-vector comparison** in postprocessing
as `LEGACY_DIRECT_DIAGNOSTIC` and writes both results to
`C19_direction_sense_diagnostic.csv`. The legacy diagnostic never contributes to
acceptance. Its purpose is to make an east/west reversal obvious and to preserve a
direct comparison with older C19 output without rerunning AMPS.

## 9. Numerical calculation

For each selected spacecraft and epoch, the runner writes a one-line GEO trajectory file:

```text
UTC latitude_deg longitude_deg_east altitude_km
```

AMPS produces a directional cutoff map on a regular global SM longitude/latitude grid. The default map resolution is 10° × 10°. For each P4/P5 detector direction, the runner:

1. rotates the observation position from GSM to SM using the driver dipole tilt;
2. constructs local physical eastward and westward detector boresights;
3. converts the AMPS incoming arrival/velocity direction to the opposite EPEAD telescope look direction;
4. selects map cells inside the channel's nominal elliptical aperture;
5. converts cutoff rigidity to proton kinetic energy;
6. integrates `E^-gamma` over the accessible part of the nominal energy interval; and
7. calculates modeled broad-aperture East and West transmissions and their ratio.

The directional map is independent of the scalar `CUTOFF_SAMPLING` result. The input templates retain a single vertical scalar sample while requesting the additional full directional map with `DIRECTIONAL_MAP T`.

### 9.1 AMPS arrival direction versus detector look direction

The production mapping is now explicit and fixed. In the AMPS cutoff source, the
directional-map vector is called the physical **ARRIVAL direction**. The same code
uses a vertical arrival direction that points toward Earth and initializes the
backward trajectory with the opposite vector. Consequently, the map vector is the
forward-time incoming particle-velocity direction.

The EPEAD EAST/WEST aperture, however, is described by the direction the telescope
faces. These two vectors are antiparallel. C19 therefore negates every directional-
map unit vector before testing membership in the physical EAST or WEST aperture.

The runner self-test builds a synthetic asymmetric map using this AMPS convention:
at a point on +SM-X, an EAST-looking telescope faces +SM-Y but particles entering
that telescope have AMPS arrival direction -SM-Y. The self-test verifies that the
production mapping gives the expected E/W sign and that the legacy direct mapping
flips it.

`C19_direction_sense_diagnostic.csv` contains both mappings for every observation:

```text
AMPS_ARRIVAL_TO_DETECTOR_LOOK   production; used for all metrics and PASS/FAIL
LEGACY_DIRECT_DIAGNOSTIC        old direct comparison; diagnostic only
```

This diagnostic is especially useful when comparing new results with an older C19
run that showed a systematic sign reversal.

### 9.2 Zero transmission is not aperture failure

Earlier C19 revisions classified a row as `INSUFFICIENT_APERTURE_COVERAGE` whenever the modeled west transmission was zero. That was misleading: a broad aperture could contain dozens of valid cells while the model predicted that the entire channel was geomagnetically blocked.

The runner now uses explicit statuses:

```text
VALID
ZERO_EAST_TRANSMISSION
ZERO_WEST_TRANSMISSION
ZERO_BOTH_TRANSMISSION
NO_EAST_APERTURE_CELLS
NO_WEST_APERTURE_CELLS
UNRESOLVED_EAST_APERTURE
UNRESOLVED_WEST_APERTURE
NEGATIVE_TRANSMISSION
NONFINITE_MODELED_RATIO
```

`ZERO_EAST_TRANSMISSION` represents `log10(E/W) -> -infinity`; `ZERO_WEST_TRANSMISSION` represents `log10(E/W) -> +infinity`. C19 does **not** replace either case by an arbitrary epsilon. They therefore do not receive a finite MAE/RMSE value, but they do contribute to the east/west **sign-agreement** metric because their sign is physically unambiguous.

### 9.3 Metrics are reported by spacecraft and in aggregate

Metrics are now produced for each:

```text
(solver, field model, spacecraft, channel)
```

and, when both spacecraft are present, an additional `spacecraft=ALL` aggregate is written for each channel. This prevents a GOES-13 orientation problem from being hidden by mixing it with GOES-15, or vice versa.

For each metric row:

- `valid_fraction` is the fraction with a **finite** modeled log ratio suitable for MAE/RMSE/correlation;
- `saturated_fraction` is the fraction with one-sided zero transmission (`ZERO_EAST_TRANSMISSION` or `ZERO_WEST_TRANSMISSION`);
- `sign_evaluable_fraction` includes both finite and one-sided saturated results;
- `sign_agreement_fraction` is evaluated over all sign-evaluable results; and
- MAE, RMSE, bias, and correlation use only finite `VALID` rows.

This design keeps quantitative metrics mathematically well defined while preventing strong zero-transmission predictions from silently disappearing from the directional validation.

### 9.4 Sensitivity and convergence controls

Important sensitivity controls include:

```text
--spectral-index
--dir-lon-res-deg
--dir-lat-res-deg
--cutoff-scan-n
--cutoff-emin-mev
--cutoff-emax-mev
--dt-trace
--max-trace-time
--mode3d-mesh-res-earth-re
```

Publication results should include convergence tests for angular resolution, rigidity scan, mesh resolution, and assumed spectral index. A useful initial convergence matrix is 10°, 5°, and 2.5° directional sampling together with `CUTOFF_UPPER_SCAN_N` values such as 120, 240, and 480. In narrow energy channels, a modest cutoff change can move an entire direction from transmitting to completely blocked, so convergence of the zero-transmission fraction is as important as convergence of finite ratios.

## 10. Outputs

The top-level output directory defaults to:

```text
test_output/C19_goes_epead_ew/
```

Per-run directories are arranged as:

```text
<solver>/<field-model>/<spacecraft>/<UTC-token>/
```

Each contains the generated input, trajectory file, AMPS log, and directional-map product.

Top-level machine-readable products:

| Product | Contents |
|---|---|
| `C19_commands.json` | Every launch command, working directory, epoch, spacecraft, solver, and field model |
| `C19_reference_used.csv` | Exact selected observational rows |
| `C19_model.csv` | Modeled E/W ratio, transmissions, aperture counts, production direction mapping, status, and map provenance |
| `C19_comparison.csv` | Observation/model rows and finite log-ratio residuals |
| `C19_metrics.csv` | Per-spacecraft and aggregate finite/saturated fractions, sign agreement, bias, MAE, RMSE, correlation, and gate result |
| `C19_direction_sense_diagnostic.csv` | Diagnostic-only comparison of selected and opposite directional-vector conventions |
| `C19_aperture_samples.csv` | Cell-level aperture membership, production direction mapping, cutoff, and transmission diagnostics for one representative case |
| `C19_result.json` | Complete result, thresholds, file hashes, direction-sense summary, failures, limitations, and overall scientific status |
| `C19_summary.txt` | Compact human-readable result with numerical, observational, and overall status separated |

Generated plots for every selected solver/model pair:

```text
C19_comparison_<solver>_<model>.png
C19_scatter_<solver>_<model>.png
C19_parity_<solver>_<model>.png
C19_residual_<solver>_<model>.png
C19_transmission_<solver>_<model>.png
C19_aperture_diagnostic.png
```

Plot interpretation:

- `C19_comparison_*` shows observed and finite modeled `log10(E/W)` versus time. Exact one-sided zero transmission is displayed explicitly at the top/bottom of the panel instead of disappearing as a NaN.
- `C19_scatter_*` is the **zoomed** comparison. Observed and modeled axes are scaled independently from their own finite data ranges, so the points occupy the plotting area even when the model is far from the observations. If the 1:1 line is outside that zoomed window, the plot says so explicitly.
- `C19_parity_*` uses a common x/y range and the 1:1 line. This is the correct figure for visually assessing absolute agreement.
- `C19_residual_*` shows finite `modeled - observed` log-ratio residuals versus time.
- `C19_transmission_*` shows broad-aperture East and West modeled transmission and is particularly useful for diagnosing one-sided zero transmission.
- `C19_aperture_diagnostic.png` visualizes the directional cells entering one representative detector response.

The zoomed scatter and parity plot intentionally answer different questions: the zoomed scatter reveals structure in the finite data, while the parity plot shows how far the model lies from equality.

## 11. Acceptance behavior

The initial observational thresholds remain provisional:

```text
finite modeled fraction      >= 0.85
correct E/W sign fraction    >= 0.90
correlation                  >= 0.60
mean absolute log10 error    <= 0.20
RMS log10 error              <= 0.30
```

The important change is that **scientific PASS/FAIL is no longer controlled by `--enforce-acceptance`**.

C19 now defines:

```text
numerical_complete     = all requested AMPS runs/postprocessing completed
observational_passed   = every reported acceptance metric passed
overall passed         = numerical_complete AND observational_passed
```

Therefore, if an individual observational metric says `FAIL`, the scientific `overall:` line also says `FAIL`. The old behavior could print individual failures while still reporting `overall: PASS` merely because acceptance enforcement was disabled; that ambiguity has been removed.

`--enforce-acceptance` now controls **only the program exit status**, which is useful for deciding whether exploratory validation failures should stop a larger automated workflow:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE --solver GRIDDED --models T96,T05 \
  --enforce-acceptance --amps ./amps -np 4 -nt 16
```

The summary explicitly separates:

```text
numerical calculation: PASS/FAIL
observational validation: PASS/FAIL
overall: PASS/FAIL
acceptance enforcement: ON/OFF
```

Exit codes remain:

```text
0  numerical calculation completed; observational failure is allowed only when
   --enforce-acceptance is OFF
1  observational gates failed and --enforce-acceptance was requested
2  input, launch, output, or postprocessing failure
```

Thus an exploratory run may still return shell exit code 0 while its scientifically meaningful `overall` result is `FAIL`.

## 12. Interpretation and limitations

C19A supports the following claim when successful:

> The AMPS cutoff calculator reproduces the sign and broad temporal/energy behavior of the GOES EPEAD East–West proton-access asymmetry for the selected event.

A systematic sign disagreement deserves special attention before model tuning. If an older C19 result has the opposite sign, inspect `C19_direction_sense_diagnostic.csv`: the production `AMPS_ARRIVAL_TO_DETECTOR_LOOK` mapping should be compared with the legacy direct-vector diagnostic. The production minus sign follows from the distinction between incoming particle direction and telescope look direction and should not be treated as a tunable convention.

C19A alone does not support a claim of exact detector count-rate prediction because:

- the full energy–angle response matrix is not used;
- out-of-aperture and secondary responses are not modeled;
- the two detector heads may retain relative calibration differences;
- the external SEP spectrum is represented by a common isotropic power law;
- the prompt event onset is unsuitable when interplanetary anisotropy is large;
- nominal GEO locations are used unless ephemeris files are supplied; and
- provisional acceptance thresholds require refinement from multiple events.

A publication-quality extension should include exact ephemeris, detector-response folding, detector-head intercalibration, spectral-slope uncertainty, additional events/yaw states, direction-mapping regression checks, angular/rigidity convergence, and T96/T05/SWMF comparisons.

## 13. Troubleshooting

### Reference file is missing

```text
C19A reference is missing ...
```

Run:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download
```

### Driver file is missing

Run:

```bash
python3 srcEarth/test/C19/tools/prepare_official_ts05_driver.py
```

### `mpirun` is not available

Load the system MPI environment or override the launcher:

```bash
--mpirun mpiexec_mpt
```

### GEO trajectory transformation fails

Verify that AMPS was built with SPICE, the required kernels are available, and the epoch is covered. `TRAJ_FRAME GEO` cannot be converted correctly without SPICE.

### Directional-map file is missing

Verify that the executable supports `DIRECTIONAL_MAP T`, that `DIRMAP_LON_RES` and `DIRMAP_LAT_RES` are positive, and that the selected UPPER_SCAN calculation completed.

### Modeled E/W has the opposite sign at nearly every epoch

1. Inspect `C19_direction_sense_diagnostic.csv`.
2. Compare `AMPS_ARRIVAL_TO_DETECTOR_LOOK` with `LEGACY_DIRECT_DIAGNOSTIC` in `C19_summary.txt`.
3. Confirm that the production result uses `AMPS_ARRIVAL_TO_DETECTOR_LOOK`.
4. Check the event-specific telemetry-head-to-physical-look-direction mapping independently.
5. If an older C19 result matches only `LEGACY_DIRECT_DIAGNOSTIC`, that is evidence that the older postprocessor compared AMPS incoming particle direction directly with the telescope look direction.

Do not swap the GOES EAST/WEST labels to compensate. The AMPS-arrival-to-look minus
sign is part of the detector geometry and is now exercised by the runner self-test.

### Many rows report `ZERO_WEST_TRANSMISSION` or `ZERO_EAST_TRANSMISSION`

This is no longer treated as missing aperture coverage. Check `C19_transmission_*`, `C19_metrics.csv`, and the cutoff-map resolution. A large saturated fraction can be physical, but it can also indicate a directional convention error, an overly coarse angular/rigidity scan, or a cutoff-energy range that brackets the detector channel poorly. Perform the convergence tests described in Section 9.4 before interpreting it as a robust model prediction.

### Too many launches

Use `--profile SMOKE`, a larger `--time-step-minutes`, one spacecraft, one model, or one solver during debugging.

## 14. References and public sources

1. NOAA National Centers for Environmental Information, *NOAA GOES 1–15 Space Environment Monitor (SEM) L1B & L2 Data, Version 0 (superseded)*, operational 5-minute GOES-13/15 EPEAD `p17ew` subset used by C19A, dataset identifier `gov.noaa.ncei.swx:goes_sem-l1b-l2_swpc`, DSI `2086_01`. NCEI recommends Version 1.0 for general use; C19A retains the fixed historical directional files for reproducibility.
2. Rodriguez, J. V., T. G. Onsager, and J. E. Mazur (2010), “The east–west effect in solar proton flux measurements in geostationary orbit: A new GOES capability,” *Geophysical Research Letters*, 37, L07109, doi:10.1029/2010GL042531.
3. Rodriguez, J. V., J. C. Krosschell, and J. C. Green (2014), “Intercalibration of GOES 8–15 solar proton detectors,” *Space Weather*, 12, 92–109, doi:10.1002/2013SW000996.
4. Kress, B. T., J. V. Rodriguez, J. E. Mazur, and M. Engel (2013), “Modeling solar proton access to geostationary spacecraft with geomagnetic cutoffs,” *Advances in Space Research*, 52(11), 1939–1948, doi:10.1016/j.asr.2013.08.019.
5. Tsyganenko, N. A., and M. I. Sitnov (2005), “Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms,” *Journal of Geophysical Research*, 110, A03208, doi:10.1029/2004JA010798.

Machine-readable BibTeX entries are provided in `references.bib`.
