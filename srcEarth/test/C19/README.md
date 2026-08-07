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

## 9. Numerical calculation

For each selected spacecraft and epoch, the runner writes a one-line GEO trajectory file:

```text
UTC latitude_deg longitude_deg_east altitude_km
```

AMPS produces a directional cutoff map on a regular global SM longitude/latitude grid. The default map resolution is 10° × 10°. For each P4/P5 detector direction, the runner:

1. rotates the observation position from GSM to SM using the driver dipole tilt;
2. constructs local physical eastward and westward boresights;
3. selects SM directional-map cells inside the channel's nominal elliptical aperture;
4. converts the cutoff rigidity to proton kinetic energy;
5. integrates `E^-gamma` over the accessible part of the nominal energy interval; and
6. calculates modeled broad-aperture East and West transmissions and their ratio.

The directional map is independent of the scalar `CUTOFF_SAMPLING` result. The input templates retain a single vertical scalar sample while requesting the additional full directional map with `DIRECTIONAL_MAP T`.

### Mode3D cutoff task parallelism

C19 launches one spacecraft position per AMPS run.  In the former Mode3D cutoff
scheduler, the smallest MPI/thread work item was one observation location.  Therefore a
C19 GRIDDED launch with `N_locations=1` reduced a requested 16-thread pool to one active
worker and executed all directional-map trajectories serially inside that worker.

The Mode3D cutoff scheduler now flattens independent trajectory products into tasks.  At
the default 10° × 10° directional resolution:

```text
longitude cells       = 360/10 = 36
latitude cells        = 180/10 + 1 = 19
directional cells     = 36*19 = 684
primary scalar cutoff = 1
total cutoff tasks    = 685 per C19 location
```

The MPI scheduler distributes these 685 tasks across ranks, and the selected shared-memory
backend distributes each rank's fetched task range across its workers.  Consequently
`-np 4 -nt 16` can expose up to 64 concurrent trajectory workers even though the input
contains only one observation location.

This is safe for the GRIDDED/MESH path because the Mode3D magnetic-field snapshot is a
complete read-only spatial field `B(x,y,z)`, not one mutable background value.  Each
worker owns a private `cMode3DMeshFieldEval` (including its tree-search hint and
interpolation state) and only reads the shared compact field arrays.  Different workers
can therefore interpolate different field values at different trajectory positions
simultaneously.

For Mode3D cutoff, `MODE3D_MPI_DYNAMIC_CHUNK` / `-mode3d-mpi-dynamic-chunk` now counts
**flattened cutoff tasks**, not locations.  Recent C19 runners use a 32-task dynamic chunk
by default; older runners may derive the GRIDDED chunk from `-nt`.  Either is sufficient
to feed 16 direct workers when `-nt 16` is selected.  If a manually selected dynamic
chunk is smaller than the thread count, the Mode3D startup banner prints a warning.

A correctly rebuilt C19 GRIDDED executable should report a task-level banner similar to:

```text
Cutoff backend : THREADS
Cutoff workers : 16 per rank
Work unit      : flattened cutoff trajectory task
Tasks/location : 685
Global tasks   : 685
MPI scheduler  : DYNAMIC
MPI dyn chunk  : 32 global cutoff task(s) per atomic fetch
```

If a one-location C19 run still reports `global location(s) per atomic fetch` or clamps
the dynamic chunk to one, that executable still contains the old location-level cutoff
scheduler.

`Task` is now the authoritative progress counter.  For task-level products the progress
line labels the derived location-equivalent counter as `LocEq`; independent tasks from the
same physical location may finish out of order on different MPI ranks.

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

Publication results should include convergence tests for angular resolution, rigidity scan, mesh resolution, and assumed spectral index.

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
| `C19_model.csv` | Modeled E/W ratio, transmissions, aperture counts, status, and map provenance |
| `C19_comparison.csv` | Observation/model rows and log-ratio residuals |
| `C19_metrics.csv` | Valid fraction, E/W sign agreement, bias, MAE, RMSE, correlation, and gate result |
| `C19_aperture_samples.csv` | Cell-level aperture membership and cutoff diagnostics for one representative case |
| `C19_result.json` | Complete result, thresholds, file hashes, failures, limitations, and overall status |
| `C19_summary.txt` | Compact human-readable result |

Generated plots for every selected solver/model pair:

```text
C19_comparison_<solver>_<model>.png
C19_scatter_<solver>_<model>.png
C19_transmission_<solver>_<model>.png
C19_aperture_diagnostic.png
```

The comparison plot shows observed and modeled `log10(E/W)` versus time. The scatter plot shows modeled versus observed log ratio with the 1:1 line. The transmission plot shows modeled broad-aperture East and West transmissions. The aperture diagnostic visualizes which directional-map cells enter the nominal detector response.

## 11. Acceptance behavior

The initial observational thresholds are provisional:

```text
valid modeled fraction       >= 0.85
correct E/W sign fraction    >= 0.90
correlation                  >= 0.60
mean absolute log10 error    <= 0.20
RMS log10 error              <= 0.30
```

By default, the runner reports these gates but returns success when the calculation completed numerically. Add `--enforce-acceptance` to return exit code 1 when an observational gate fails:

```bash
python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE --solver GRIDDED --models T96,T05 \
  --enforce-acceptance --amps ./amps -np 4 -nt 16
```

Exit codes:

```text
0  numerical calculation completed; and observational gates passed when enforced
1  observational gates failed and --enforce-acceptance was requested
2  input, launch, output, or postprocessing failure
```

## 12. Interpretation and limitations

C19A supports the following claim when successful:

> The AMPS cutoff calculator reproduces the sign and broad temporal/energy behavior of the GOES EPEAD East–West proton-access asymmetry for the selected event.

C19A alone does not support a claim of exact detector count-rate prediction because:

- the full energy–angle response matrix is not used;
- out-of-aperture and secondary responses are not modeled;
- the two detector heads may retain relative calibration differences;
- the external SEP spectrum is represented by a common isotropic power law;
- the prompt event onset is unsuitable when interplanetary anisotropy is large;
- nominal GEO locations are used unless ephemeris files are supplied; and
- provisional acceptance thresholds require refinement from multiple events.

A publication-quality extension should include exact ephemeris, detector-response folding, detector-head intercalibration, spectral-slope uncertainty, additional events/yaw states, and T96/T05/SWMF comparisons.

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

### Too many launches

Use `--profile SMOKE`, a larger `--time-step-minutes`, one spacecraft, one model, or one solver during debugging.

## 14. References and public sources

1. NOAA National Centers for Environmental Information, *NOAA GOES 1–15 Space Environment Monitor (SEM) L1B & L2 Data, Version 0 (superseded)*, operational 5-minute GOES-13/15 EPEAD `p17ew` subset used by C19A, dataset identifier `gov.noaa.ncei.swx:goes_sem-l1b-l2_swpc`, DSI `2086_01`. NCEI recommends Version 1.0 for general use; C19A retains the fixed historical directional files for reproducibility.
2. Rodriguez, J. V., T. G. Onsager, and J. E. Mazur (2010), “The east–west effect in solar proton flux measurements in geostationary orbit: A new GOES capability,” *Geophysical Research Letters*, 37, L07109, doi:10.1029/2010GL042531.
3. Rodriguez, J. V., J. C. Krosschell, and J. C. Green (2014), “Intercalibration of GOES 8–15 solar proton detectors,” *Space Weather*, 12, 92–109, doi:10.1002/2013SW000996.
4. Kress, B. T., J. V. Rodriguez, J. E. Mazur, and M. Engel (2013), “Modeling solar proton access to geostationary spacecraft with geomagnetic cutoffs,” *Advances in Space Research*, 52(11), 1939–1948, doi:10.1016/j.asr.2013.08.019.
5. Tsyganenko, N. A., and M. I. Sitnov (2005), “Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms,” *Journal of Geophysical Research*, 110, A03208, doi:10.1029/2004JA010798.

Machine-readable BibTeX entries are provided in `references.bib`.
