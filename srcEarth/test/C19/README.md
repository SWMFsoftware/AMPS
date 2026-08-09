# C19A — GOES EPEAD East–West directional-access validation

## 1. Purpose

C19A is the first observational validation of the **directional** geomagnetic-access capability of the AMPS cutoff calculator. Internal symmetry tests such as C17 can verify charge-sign and velocity-reversal consistency, but they do not establish that the calculated directional access agrees with measurements. C19A compares the AMPS directional access function with simultaneous eastward- and westward-looking proton measurements from the GOES-13 and GOES-15 EPEAD instruments. The production path uses the explicit three-state `A(E,Ω)` access cube; a full directional cutoff/penumbra map is generated only when `--cutoff-search PENUMBRA_SCAN` is requested.

The implemented public-data event is the 17 May 2012 SEP/GLE71 event. The default analysis interval is the event decay from 2012-05-17 06:00 UTC through 2012-05-18 06:00 UTC. The prompt onset is excluded because interplanetary beam anisotropy can imitate or obscure a geomagnetic East–West effect.

The primary observable is

```text
log10[(physical EAST background-subtracted flux) /
      (physical WEST background-subtracted flux)]
```

for EPEAD P4 and P5. A negative value means that the eastward-looking detector measured less flux than the westward-looking detector.

C19A is a **broad-aperture observational validation**. The current science chain is explicit: the reference records the exact NOAA flux variables and ephemeris provenance; a time-dependent event spectrum is estimated from the physical-WEST P4/P5 measurements (or supplied independently); both GRIDDED Mode3D and GRIDLESS evaluate the same direct three-state access function `A(E,Ω)` on the directional sky grid; and the postprocessor folds that access through a file-defined detector energy response and the nominal elliptical angular aperture. The current runner also includes fine angular/mesh defaults, external detector-attitude support, a bounded upstream-anisotropy option, explicit finite-rigidity-grid uncertainty, and a high-energy response extension appropriate to the default uncorrected GOES product. The response is still factorized in energy and angle rather than a complete calibrated energy-angle matrix, so the test remains a controlled approximation rather than a full instrument simulator.

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
| AMPS product | Default: `DIRECT_ACCESS` three-state `A(E,Ω)` cube only. Optional diagnostic: `PENUMBRA_SCAN` map plus the same direct `A(E,Ω)` companion cube |
| Source spectrum | Default: epoch-dependent power law fit to physical-WEST background-subtracted P4/P5 flux; explicit spectrum CSV and fixed-gamma sensitivity modes are supported |
| Instrument model | File-defined piecewise energy response (`data/epead_response_C19_uncorrected_extended.csv` by default) × nominal elliptical angular aperture; identical direct `A(E,Ω)` fold for GRIDDED and GRIDLESS |

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
├── C19_trajectory.txt
├── build_goes_reference.py
├── event_C19_may2012.json
├── run_C19.py
├── tests/
│   └── run_self_tests.py
└── data/
    ├── reference_C19_goes_epead_ew.csv.gz
    ├── reference_C19_goes_epead_ew_provenance.json
    ├── epead_response_C19_uncorrected_extended.csv
    ├── epead_response_C19_nominal.csv
    └── ts05_driver_may2012.txt
```

The compact observational reference and the five-minute T05 driver used by this
snapshot are committed with the test so a routine C19 run is reproducible without a
network download.  `build_goes_reference.py` remains the supported way to regenerate
the GOES reference from the public NOAA files; its provenance JSON records the inputs
and hashes.

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

Python 3.8 or newer is supported.  The C19 runner and CSV reference-building path use
the Python standard library for their core workflow.  Reading the native NOAA/NCEI
`goes-l2-orb1m` NetCDF ephemeris requires `xarray`; converting that ephemeris to CSV
removes the NetCDF dependency.  Plot generation is optional and requires Matplotlib:

```bash
python3 -m pip install xarray matplotlib
```

If Matplotlib is unavailable, C19 still writes the CSV/JSON diagnostics and reports that
plot generation was skipped.

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

From the AMPS repository root, download the particle files and supply the official NOAA/NCEI one-minute GOES ephemeris for both spacecraft. The reference builder intentionally does **not** silently replace a missing ephemeris with a nominal slot:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download \
  --goes13-ephemeris /path/to/goes13_goes-l2-orb1m.nc \
  --goes15-ephemeris /path/to/goes15_goes-l2-orb1m.nc
```

Native NOAA NetCDF (`geo_llr`) and CSV ephemeris files are supported. For an explicit legacy reproduction only, `--allow-nominal-ephemeris` restores the nominal-slot fallback and labels every such row `NOMINAL_GEO_SLOT_LEGACY_OVERRIDE`.

This writes:

```text
srcEarth/test/C19/data/cache/g13_epead_p17ew_5m_20120501_20120531.csv
srcEarth/test/C19/data/cache/g15_epead_p17ew_5m_20120501_20120531.csv
srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz
srcEarth/test/C19/data/reference_C19_goes_epead_ew_provenance.json
```

The provenance file records source URLs and SHA-256 values, manifest checksum, background method, directional mapping, exact flux-product policy, exact source flux variables actually selected, ephemeris policy, and generated-reference checksum. The same flux-variable/correction-state fields are also embedded in every reference row, so a sub-selected `C19_reference_used.csv` remains auditable without the sidecar JSON.

To use files downloaded separately, provide both the directional particle products and
the real one-minute ephemeris products:

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/g13_epead_p17ew_5m_20120501_20120531.csv \
  --goes15-particle /path/to/g15_epead_p17ew_5m_20120501_20120531.csv \
  --goes13-ephemeris /path/to/goes13_orb1m.nc \
  --goes15-ephemeris /path/to/goes15_orb1m.nc
```

### 5.2 Reference construction

The current reference builder removes the previous silent `UNCOR_FLUX → FLUX → COR_FLUX` fallback. The event manifest now declares `goes_flux_product_policy`, and the builder accepts an explicit override:

```bash
--flux-product UNCORRECTED|CORRECTED|AUTO
```

`UNCORRECTED` accepts only `*_UNCOR_FLUX`; `CORRECTED` never falls back to an uncorrected variable; `AUTO` exists only for deliberate legacy reproduction. Each output row records `east_flux_variable`, `west_flux_variable`, and their correction-state labels.

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

Science reference generation requires real spacecraft ephemeris. The builder accepts the official NOAA/NCEI GOES 6–15 one-minute orbit product in native NetCDF form (`goes-l2-orb1m`, variable `geo_llr = [latitude, longitude, geocentric radius]`) or a CSV containing UTC/longitude/latitude and altitude or radius. The nearest sample within 180 s is used.

```bash
python3 srcEarth/test/C19/build_goes_reference.py \
  --goes13-particle /path/to/g13_particle.csv \
  --goes15-particle /path/to/g15_particle.csv \
  --goes13-ephemeris /path/to/goes13_orb1m.nc \
  --goes15-ephemeris /path/to/goes15_orb1m.nc
```

If a selected particle epoch has no real ephemeris within 180 s, reference generation fails. `--allow-nominal-ephemeris` is an explicit legacy override only. The runner can additionally enforce science provenance with `--require-real-ephemeris`.

### 5.4 Reference-builder self-test

```bash
python3 srcEarth/test/C19/build_goes_reference.py --self-test
```

This uses synthetic NOAA-format files to test header discovery, P4/P5 parsing, quality/background filtering, the event-specific E/W mapping, gzip output, and provenance generation. It does not contact NOAA.

## 6. T05/TS05 driver used by C19

This source snapshot contains the validated five-minute event driver directly at:

```text
srcEarth/test/C19/data/ts05_driver_may2012.txt
```

`run_C19.py` uses that file by default.  It independently verifies the
timestamp-plus-19-value schema, strict time ordering, five-minute median cadence,
maximum ten-minute internal gap, finite values, and coverage of all selected
observation epochs before launching AMPS.  The source/archive provenance is retained
in the comments at the top of the committed driver and the runner records its SHA-256
in `C19_result.json`.

If the event driver is regenerated externally, pass the replacement explicitly with
`--driver /path/to/driver.txt`; the same validation is applied before any run.

## 7. Verify the package before a production run

Run all non-network package self-tests:

```bash
python3 srcEarth/test/C19/tests/run_self_tests.py
```

The package test compiles both Python programs, runs their unit/self-tests, performs a
normal BOTH-solver dry run, and performs a P0 dry run that verifies the A/B/C
algorithm-policy matrix, the 400/800/1600-Re trace-budget scaling, and the centered-
dipole Störmer-anchor input.  No AMPS executable or network access is required for
these checks.

The component tests can also be run separately:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --self-test
python3 srcEarth/test/C19/run_C19.py --self-test
```

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


### 8.8 Single current workflow

C19 no longer exposes P0/P1/P2 as alternate execution modes. Those names describe the
development history only. The validated changes are integrated into every ordinary run.
In particular, the runner always renders:

```text
CUTOFF_SEARCH_ALGORITHM   PENUMBRA_SCAN
CUTOFF_TRACE_LIMIT_POLICY UNRESOLVED
```

and both GRIDDED Mode3D and GRIDLESS always request the direct three-state `A(E,Ω)` companion cube used
for the synthetic detector fold. The current production defaults also adopt the finest
settings from the former P2 convergence implementation:

```text
directional grid:       2.5° × 2.5°
near-Earth mesh:        0.025 Re
boundary mesh:          1.0 Re
mover:                  RK4
upstream model:         ISOTROPIC (unless explicitly changed)
```

There is therefore no `--p0-diagnostic` or `--p2-diagnostic` runner flag. A normal
SMOKE/ROUTINE/FULL invocation follows the same code path and, after successful
post-processing, always produces the standard comparison plots. Numerical studies can
still be performed by explicitly changing physical/numerical inputs such as angular
resolution, mesh resolution, trace budget, detector orientation, or anisotropy, but
those changes do not select a different runner implementation.

## 9. Numerical calculation and validity architecture

For each selected spacecraft and epoch, the runner writes a one-line GEO trajectory file:

```text
UTC latitude_deg longitude_deg_east altitude_km
```

AMPS represents detector directions on a regular global SM longitude/latitude grid.
The current C19 default angular resolution is 2.5° × 2.5°, with the grid pruned to the
union of requested instrument apertures unless `FULL_SPHERE` is selected. In the default
`DIRECT_ACCESS` mode AMPS writes only the explicit three-state access cube on those cells.
`PENUMBRA_SCAN` additionally writes a long-form directional cutoff map. The first four
historical variables in that optional map are preserved:

```text
lon_deg lat_deg Rc_GV Emin_MeV
```

and the following diagnostics are appended:

```text
Rc_lower_GV
Rc_effective_GV
Rc_upper_GV
n_transitions
n_allowed_intervals
n_unresolved_samples
lower_bracket_unresolved
upper_bracket_unresolved
lower_below_range / lower_above_range
upper_below_range / upper_above_range
n_trajectory_evaluations
n_outer_boundary_allowed
n_inner_boundary_forbidden
n_magnetically_trapped_forbidden
n_time_limit
n_step_limit
n_distance_limit
max_trace_time_s
max_trace_distance_Re
max_trace_steps
```

Mode3D DIPOLE maps can additionally contain `Rc_stormer_GV` for independent standalone Størmer regression checks.

`Rc_GV` remains the scalar value consumed by existing readers. For a directional
`PENUMBRA_SCAN` map it is `Rc_effective_GV`; the lower and upper cutoffs and topology are
retained explicitly instead of being discarded after the scan.

For each P4/P5 physical detector direction, `run_C19.py`:

1. rotates the observation position from GSM to SM using the driver dipole tilt;
2. constructs the physical detector boresight (or reads the supplied attitude vector);
3. converts the AMPS incoming arrival/velocity direction to the opposite EPEAD telescope look direction;
4. selects regular-grid cells inside the detector's elliptical aperture;
5. reads the three-state access value at each explicitly requested detector-response rigidity;
6. folds `A(E,Ω)` with the event spectrum and detector response, propagating unresolved and finite-rigidity-grid uncertainty; and
7. calculates East and West transmission and their E/W ratio only when the response bounds satisfy the configured validity criteria.

`DIRECT_ACCESS` is independent of the scalar `CUTOFF_SAMPLING` cutoff result even though
the common input syntax retains `CUTOFF_SAMPLING VERTICAL`. In production mode no scalar
cutoff trajectory and no directional cutoff-map/PENUMBRA_SCAN task is scheduled. The
runner constructs the geometry-only directional map needed for aperture folding directly
from the `A(E,Ω)` cube. By default only regular-grid cells needed by the requested EPEAD
apertures are scheduled; `FULL_SPHERE` remains available as an explicit reference coverage.
The retained cells have exactly the same SM lon/lat coordinates as the corresponding cells
in a full-sphere run. `PENUMBRA_SCAN` remains an explicit diagnostic mode and writes the
full lower/effective/upper cutoff topology in addition to the same direct-access cube.

### 9.1 Three-state cutoff classification

There are two cutoff-classification paths in the current solver. The historical
`UPPER_SCAN` path calls the Boolean `TraceAllowed3D`/gridless equivalent. That Boolean
interface intentionally maps every recognized cutoff-forbidden termination, including
configured time/step/distance caps, to `false`. It does not use the three-state trace
policy.

Both `DIRECT_ACCESS` and `PENUMBRA_SCAN` use the structured three-state classifier for
the samples that feed C19. A trajectory that escapes through the outer boundary is
`ALLOWED`; a physical inner-boundary loss or validated magnetic-trap termination is
`PHYSICAL_FORBIDDEN`; and a configured `TIME_LIMIT`, `STEP_LIMIT`, or `DISTANCE_LIMIT` is
`UNRESOLVED` when the input requests `CUTOFF_TRACE_LIMIT_POLICY UNRESOLVED`.

C19 therefore does **not** attempt to make `UPPER_SCAN` pass by increasing a timeout or by
relabeling its result after the calculation. The production runner defaults to
`DIRECT_ACCESS`; `--cutoff-search PENUMBRA_SCAN` requests the more expensive full cutoff-band
diagnostic while retaining the same direct `A(E,Ω)` companion samples. The historical
Boolean `UPPER_SCAN` path is not selectable from the C19 runner CLI.

### 9.2 Unresolved aperture cells are never silently removed

Earlier C19 folding skipped an unresolved directional cell and divided by the solid angle
of the remaining cells. A mostly unresolved aperture could therefore look like a precise
finite detector response.

The new fold keeps all in-aperture solid-angle weights in the denominator and writes
bounds:

```text
T_min : every unresolved cell is assigned transmission 0
T_max : every unresolved cell is assigned transmission 1
```

The model row contains:

```text
east_transmission_min / east_transmission_max
west_transmission_min / west_transmission_max
unresolved_east_fraction
unresolved_west_fraction
```

A central transmission is reported only when the unresolved solid-angle fraction is at
or below `--max-unresolved-aperture-fraction` (default 0.05). Above that threshold the
row is explicitly invalidated as:

```text
EXCESSIVE_UNRESOLVED_EAST_APERTURE
EXCESSIVE_UNRESOLVED_WEST_APERTURE
```

This tolerance is explicit and auditable; unresolved sky cells can no longer disappear
from the detector normalization.

### 9.3 Termination-reason forensics

A cutoff value alone cannot distinguish a physical shielding transition from a finite
trajectory-budget artifact. `PENUMBRA_SCAN` therefore accumulates raw termination counts
for every unique rigidity evaluation used by the coarse scan, transition refinement, and
effective-cutoff integration.

The long-form map records, separately:

```text
OUTER_BOUNDARY_ALLOWED
INNER_BOUNDARY_FORBIDDEN
MAGNETICALLY_TRAPPED_FORBIDDEN
TIME_LIMIT
STEP_LIMIT
DISTANCE_LIMIT
```

The `n_unresolved_samples` topology field refers to unresolved samples on the regular
cutoff-band grid; the termination counters cover all unique trajectory evaluations,
including refinement/integration evaluations. Use the latter when diagnosing why a
particular direction exhausted the trace budget.

### 9.4 Trace-budget and frozen-field guardrail

Trace-budget studies, when needed, are performed by rerunning the same current workflow with explicitly changed `--max-trace-distance-re`, `--max-trace-time`, and `--max-steps`. The purpose is **not** to keep increasing the path length until every trajectory escapes.
Each C19 epoch uses a frozen background field. A trajectory requiring a very long
integration can enter a regime where the static-field approximation itself deserves
scrutiny. The runner therefore carries a separate guardrail:

```text
--frozen-field-warning-seconds 300
```

If any contributing directional scan reports a longer individual trajectory, the row is
flagged `STATIC_FIELD_TRACE_GUARDRAIL` and is excluded from quantitative E/W validation.
The trajectory is not reclassified as physically forbidden.

### 9.5 Independent centered-dipole Størmer regression

GRIDDED-versus-GRIDLESS agreement is useful but is only an internal consistency test; a shared sign or trajectory-classification error could affect both. For independent solver regression, `FIELD_MODEL DIPOLE` directional output can include the analytic Størmer cutoff `Rc_stormer_GV` at the same position and incoming direction as the numerical sky cell. Such a standalone regression checks trajectory/cutoff machinery independently of T05 and GOES and must not be tuned using the May-2012 observations.

### 9.6 AMPS arrival direction versus detector look direction

The production mapping remains explicit and fixed. In the AMPS cutoff source, the
directional-map vector is the physical **ARRIVAL direction**. The backward trajectory is
initialized with the opposite vector. The EPEAD EAST/WEST aperture is described by the
direction the telescope faces, so C19 uses:

```text
detector look direction = - AMPS directional-map arrival direction
```

The runner self-test uses a synthetic asymmetric map to verify the sign. The legacy
direct-vector comparison is still written as `LEGACY_DIRECT_DIAGNOSTIC`, but never enters
acceptance.

### 9.7 Zero transmission is a resolved saturation state, not missing coverage

Exact one-sided zero transmission remains a valid model state when the aperture itself is
resolved:

```text
ZERO_EAST_TRANSMISSION  -> log10(E/W) -> -infinity
ZERO_WEST_TRANSMISSION  -> log10(E/W) -> +infinity
ZERO_BOTH_TRANSMISSION
```

C19 does not manufacture an epsilon. One-sided zeros are excluded from finite MAE/RMSE
but retained in saturation and sign diagnostics.

Additional non-saturation statuses now include:

```text
NO_EAST_APERTURE_CELLS
NO_WEST_APERTURE_CELLS
EXCESSIVE_UNRESOLVED_EAST_APERTURE
EXCESSIVE_UNRESOLVED_WEST_APERTURE
UNRESOLVED_EAST_APERTURE
UNRESOLVED_WEST_APERTURE
STATIC_FIELD_TRACE_GUARDRAIL
NEGATIVE_TRANSMISSION
NONFINITE_MODELED_RATIO
```

### 9.8 Staged validity gates

C19 no longer calls a run numerically complete merely because all processes returned
successfully. Four stages are reported separately:

```text
execution_complete
trajectory_resolution_passed
instrument_fold_valid
observational_passed
```

`execution_complete` means the requested AMPS runs and postprocessing completed.
`trajectory_resolution_passed` additionally requires each EAST/WEST aperture to satisfy
the unresolved-solid-angle threshold and the frozen-field duration guardrail.
`instrument_fold_valid` requires usable aperture geometry and transmission bounds.
`observational_passed` is the provisional GOES agreement gate.

The scientific overall state is:

```text
passed = execution_complete
         AND trajectory_resolution_passed
         AND instrument_fold_valid
         AND observational_passed
```

For compatibility, `C19_result.json` still contains `numerical_complete`, but it is marked
deprecated and now means `execution_complete AND trajectory_resolution_passed`.

### 9.9 Sensitivity and convergence controls

The production trajectory classification defaults to `DIRECT_ACCESS + UNRESOLVED`.
`--cutoff-search PENUMBRA_SCAN` explicitly enables the more expensive full cutoff-band
diagnostic while retaining the same direct `A(E,Ω)` companion cube. Historical Boolean
cutoff modes are not exposed by the C19 runner. Numerical
inputs that remain intentionally configurable for controlled studies are:

```text
--max-unresolved-aperture-fraction
--frozen-field-warning-seconds
--spectral-index
--dir-lon-res-deg
--dir-lat-res-deg
--cutoff-scan-n
--cutoff-emin-mev
--cutoff-emax-mev
--dt-trace
--max-steps
--max-trace-time
--max-trace-distance-re
--mode3d-mesh-res-earth-re
--mode3d-mesh-res-boundary-re
```

The current workflow integrates three-state trajectory classification, source-variable provenance, real-ephemeris support, event-derived spectra, direct directional rigidity access, detector-response folding, finite-rigidity-grid uncertainty bounds, documented high-energy P4/P5 response support, fine angular/mesh settings, exact-file detector geometry support, and explicit upstream-anisotropy controls. Replacement of the factorized response by a fully calibrated energy-angle response matrix remains a separate instrument-calibration task.

## 9.10 Observational-commensurability implementation

### Exact GOES source-variable provenance

`build_goes_reference.py` uses an explicit flux-product policy and records the exact NOAA column used for each physical EAST/WEST measurement. No correction state can change silently. Legacy committed reference rows that predate exact provenance are read as `LEGACY_UNRECORDED`; regenerate the reference to obtain full provenance.

### Real GOES ephemeris

The reference builder supports the native NOAA/NCEI one-minute NetCDF orbit product and CSV equivalents. Real ephemeris is required by default when rebuilding the reference. The old manifest slot is available only through `--allow-nominal-ephemeris`. For science runs use `run_C19.py --require-real-ephemeris` so a legacy nominal-slot reference cannot accidentally enter publication output.

### Event-derived spectrum

The default `--spectrum-source OBSERVED_WEST` fits

```text
J(E,t) = J0(t) [E/E0]^-gamma(t)
```

from the physical-WEST background-subtracted P4/P5 intensities at each epoch, using channel effective energies stored in the event manifest. If quality filtering removes one channel at an isolated epoch, gamma is interpolated from neighboring two-channel fits and the available WEST channel sets the normalization; this is marked explicitly in `C19_spectrum_used.csv`. An independent upstream spectrum is preferable when available:

```bash
--spectrum-source FILE --spectrum-file spectrum.csv
```

with `utc,gamma[,j0,e0_mev]`. `--spectrum-source FIXED --spectral-index ...` is retained only for sensitivity/legacy reproduction.

### Direct GRIDDED/GRIDLESS `A(E,Ω)` product

C19 supports two current cutoff products selected with `--cutoff-search`:

```text
DIRECT_ACCESS   (default)
    Trace only the requested (direction,rigidity) samples needed by the detector fold.
    No scalar cutoff and no directional PENUMBRA_SCAN map are computed.

PENUMBRA_SCAN
    Compute the full lower/effective/upper cutoff topology for every selected direction
    and, in addition, trace the same requested direct-access rigidity list used by the
    detector fold.
```

With `DIRECT_ACCESS`, `DIRECTIONAL_MAP T`, and a non-empty `CUTOFF_RIGIDITY_LIST_GV`,
both solvers write the direct cube without first calculating a directional cutoff map.
With `PENUMBRA_SCAN`, the same cube is written as a companion product. The solver-specific
file names are: The solver-specific file names are:

```text
GRIDDED : cutoff_3d_dir_access_loc_000000.dat
GRIDLESS: cutoff_gridless_dir_access_point_0000.dat
```

Each row is one `(lon,lat,rigidity)` trajectory and contains `rigidity_GV`, corresponding proton `energy_MeV`, `access_state`, `allowed`, and `unresolved`. States use the same three-state classifier as the production penumbra path:

```text
0 = PHYSICAL_FORBIDDEN
1 = ALLOWED
2 = UNRESOLVED
```

The GRIDDED and GRIDLESS direct-access files are required to carry the same sky-cell set and the same ordered energy/rigidity grid for a given C19 case. In `PENUMBRA_SCAN` mode the runner additionally checks the cube against the companion directional map. In `DIRECT_ACCESS` mode no cutoff map exists by design; the runner constructs a geometry-only map directly from the cube's frame, observation position, and `(lon,lat)` cells. Missing or truncated direct-access output fails post-processing rather than falling back to a scalar-cutoff proxy.

The task scheduler flattens `(location, sky-cell, rigidity)` work so those trajectories are distributed across MPI ranks rather than serialized inside one map cell. Both Mode3D and GRIDLESS may additionally use their configured intra-rank thread backend. `DIRECT_ACCESS` is the default science product because the detector fold consumes the rigidity-list cube itself; `PENUMBRA_SCAN` is retained when full cutoff topology is explicitly required. The two files intentionally have the same columns and three-state semantics so `run_C19.py` uses one parser and one folding implementation.

A terminology detail is important: `RIGIDITY_LIST` remains the historical shell-oriented Mode3D product. C19 uses the distinct `DIRECT_ACCESS` token for point/trajectory directional access. GRIDDED and GRIDLESS implement the same `DIRECT_ACCESS` semantics and schedule one independent three-state trajectory for every requested `(direction,rigidity)` pair.

### Runtime implication and execution commands

For a representative 2.5-degree C19 aperture selection with 1,760 retained directions
and 55 requested detector-response rigidities, `DIRECT_ACCESS` schedules

```text
1760 x 55 = 96,800 trajectory tasks
```

per observation location. The old/full `PENUMBRA_SCAN` path schedules those same 96,800
direct-access trajectories **plus** 1,760 directional penumbra tasks. Each penumbra task
internally evaluates roughly `CUTOFF_UPPER_SCAN_N` trajectories (120 by default), so the
physical trajectory count is approximately 308,000 before any extra refinement. This is
why `DIRECT_ACCESS` is now the production default.

One-line GRIDDED commands (4 MPI ranks, 16 threads/rank, T05, SMOKE profile) are:

```bash
python3 srcEarth/test/C19/run_C19.py --profile SMOKE --solver GRIDDED --models T05 --cutoff-search DIRECT_ACCESS --amps ./amps -np 4 -nt 16 --keep
python3 srcEarth/test/C19/run_C19.py --profile SMOKE --solver GRIDDED --models T05 --cutoff-search PENUMBRA_SCAN --amps ./amps -np 4 -nt 16 --keep
```

The same switch works with `--solver GRIDLESS` or `--solver BOTH`.

### Response fold of direct access

For either solver, `run_C19.py` builds a response-energy grid over the **complete positive support of the configured response file**, writes the same strictly increasing rigidity list into the input, reads that solver's direct access cube, and evaluates the head transmission schematically as

```text
T_head(t) = ∫dΩ w(Ω) ∫dE J(E,t) G_channel(E) A(E,Ω,t)
            ------------------------------------------------
             ∫dΩ w(Ω) ∫dE J(E,t) G_channel(E)
```

where `w(Ω)` is the nominal elliptical angular-aperture weight and `G_channel(E)` is read from `--detector-response`.

The energy quadrature no longer linearly interpolates the binary access state. For every adjacent pair of checked rigidities:

```text
ALLOWED   -> ALLOWED      whole interval transmitted
FORBIDDEN -> FORBIDDEN    whole interval blocked
ALLOWED   -> FORBIDDEN    transition lies somewhere in the interval
FORBIDDEN -> ALLOWED      transition lies somewhere in the interval
anything  -> UNRESOLVED   physical state of the interval is unresolved
```

For a resolved state change the whole response-weighted interval is carried as a finite-grid uncertainty `[0, full interval]`; it is **not** replaced by the artificial 0.5 contribution produced by trapezoidal interpolation of the 0/1 endpoints. `C19_model.csv` records `discrete_transition_east_fraction` and `discrete_transition_west_fraction`. The default `--max-discrete-transition-fraction 0.05` invalidates a quantitative E/W value when too much detector response lies inside such unresolved transition brackets. Increasing `--access-energy-points` narrows those brackets and is therefore a direct rigidity-grid convergence test.

`UNRESOLVED` trajectory states remain a separate uncertainty source and continue to use `--max-unresolved-aperture-fraction`. The lower/upper signal and transmission bounds include both effects, while their diagnostic fractions remain separate so a coarse rigidity grid cannot be confused with a trace-limit problem.

`J(E)G(E)` is integrated analytically over every energy interval for the current power-law spectrum and piecewise-constant response. Consequently exact detector-response edges are inserted as grid nodes but the old `E*(1±10^-8)` edge-bracketing trajectories have been removed.

The default `epead_response_C19_uncorrected_extended.csv` matches the default **uncorrected** GOES reference more closely than the old primary-only response. It contains the nominal primary bands plus the documented GOES 8–15 secondary proton energy ranges from Rodriguez et al. (2017), Table A2 / the GOES-N EPS calibration handbook:

| Channel | Component | Energy range (MeV) | Published geometrical factor (cm² sr) | C19 relative weight |
|---|---|---:|---:|---:|
| P4 | primary | 15–40 | 0.21 | 1.0 |
| P4 | secondary | 80–115 | 0.038 | 0.180952 |
| P4 | secondary | 115–150 | 0.25 | 1.190476 |
| P5 | primary | 38–82 | 0.36 | 1.0 |
| P5 | secondary | 80–110 | 0.091 | 0.252778 |
| P5 | secondary | 110–150 | 0.57 | 1.583333 |
| P5 | secondary | 150–190 | 0.21 | 0.583333 |

The relative weights are each secondary geometrical factor divided by that channel's primary factor; absolute normalization cancels in the same-channel E/W ratio. With this file the direct-access grid extends automatically to **190 MeV** instead of ending at 82 MeV. The primary-only `epead_response_C19_nominal.csv` remains available for corrected-flux or controlled sensitivity studies.

Important limitation: the published secondary factors are integrated side/rear penetrating-particle responses, not a resolved angular response matrix. C19 currently treats them as an **energy-response extension inside the same factorized angular model**. That is better than assigning zero response above 82 MeV for an uncorrected product, but it is not a substitute for a measured energy-angle response. If a corrected NOAA flux product is selected, the runner warns when a response containing `SECONDARY` components is used because that would tend to double-count the correction.

## 9.11 Current production numerical and physics settings

The former P2 work is now part of the normal C19 implementation rather than a separate
one-epoch diagnostic suite. The runner uses a single production path.

### Angular resolution

The default directional map is `2.5° × 2.5°`, the finest level from the previous
10°/5°/2.5° convergence ladder. `--dir-lon-res-deg` and `--dir-lat-res-deg` remain
available only when intentionally performing a new resolution study.

### Directional coverage: arbitrary instrument apertures versus full sphere

The regular 2.5° sky grid contains 144 longitude values and 73 latitude values, or
10,512 nominal cells. Earlier C19 versions traced every one of those directions even
though the synthetic observable subsequently uses only directions viewed by the
instrument heads at that epoch. The optimized implementation is now **instrument-
neutral**: AMPS is given one or more physical detector-look vectors and schedules only
the regular-grid cells inside the union of those apertures.

The current C19 default is:

```text
DIRMAP_COVERAGE      VECTOR_APERTURES
DIRMAP_APERTURE_FILE C19_directional_apertures.dat
```

`run_C19.py` generates `C19_directional_apertures.dat` separately for every spacecraft
and epoch. For an authoritative attitude file, each head can point in an arbitrary
direction and the heads do **not** have to be antipodal. The May-2012 observational
reference still reports the conventional physical EAST/WEST flux ratio, but it also
records `telemetry_head_east` and `telemetry_head_west`. C19 uses those fields to map
each numerator/denominator stream back to the **actual instrument head** (for example
raw head `E` or `W`) and then looks up that head's time-dependent boresight. Thus the
words EAST/WEST do not define the geometry used by AMPS.

Each aperture record has the form:

```text
# name frame bx by bz upx upy upz horizontal_half_angle_deg vertical_half_angle_deg
```

The fields define **one finite instrument field of view**, not one particle trajectory.
Their meanings are:

| Field | Meaning |
|---|---|
| `name` | Arbitrary aperture/head identifier such as `E`, `W`, `HEAD_A`, or `TELESCOPE_1`. The name carries **no geometric meaning**; the vectors determine the pointing. |
| `frame` | Coordinate frame in which both vector triplets are supplied: `SM`, `GSM`, or `LOCAL_SM`. |
| `bx by bz` | Components of the detector **LOOK boresight** vector. This is the direction in the sky toward which the instrument points. |
| `upx upy upz` | Components of a roll/up reference vector. Together with the boresight, this fixes the rotation of a non-circular/elliptical FOV about the boresight. |
| `horizontal_half_angle_deg` | Angular semi-width of the elliptical FOV along its derived horizontal axis, in degrees. |
| `vertical_half_angle_deg` | Angular semi-width of the elliptical FOV along its derived vertical/up axis, in degrees. |

For example,

```text
HEAD_A SM  0.0 1.0 0.0   0.0 0.0 1.0   30.0 60.0
```

defines an aperture whose detector looks toward SM `+Y`, whose nominal up direction is
SM `+Z`, and whose FOV has horizontal and vertical half-angles of 30 degrees and 60
degrees, respectively. The vector magnitudes are not used as physical quantities: the
solver normalizes the boresight and orthonormalizes the supplied up reference. The
boresight must be nonzero, and the up reference must not be parallel (or numerically
nearly parallel) to the boresight.

The supplied up vector need not be exactly perpendicular to the boresight. After both
vectors have been transformed into the directional-map frame, the solver constructs the
aperture basis as

```text
b = normalize(boresight)
v = normalize(up - (up dot b) b)
h = normalize(v cross b)
```

where `b` is the detector look boresight, `v` is the vertical/up aperture axis, and `h`
is the horizontal aperture axis. Thus a small parallel component in an attitude-file up
vector is removed automatically rather than rotating the requested FOV.

A directional-grid cell is converted from the AMPS particle-arrival direction to the
detector-look direction,

```text
look = -arrival
```

and is first required to lie in the forward hemisphere (`look dot b > 0`). The solver
then forms the horizontal and vertical angular offsets

```text
ah = atan2(look dot h, look dot b)
av = atan2(look dot v, look dot b)
```

and retains the grid cell when

```text
(ah / horizontal_half_angle_deg)^2
  + (av / vertical_half_angle_deg)^2 <= 1
```

with `ah` and `av` expressed in degrees. Therefore the two half-angle fields are the
semi-axes of an **elliptical angular aperture**; they are not independent rectangular
`|ah|`/`|av|` cuts.

The supported frames are:

- `SM`: `boresight` and `up` are Cartesian components in the global Solar Magnetic
  frame used to label the directional map.
- `GSM`: the vectors are Cartesian components in GSM. AMPS transforms them into the SM
  directional-map frame at the run epoch before constructing the aperture basis.
- `LOCAL_SM`: the three components are coefficients in the local orthonormal
  `(radial, local-east, local-north)` basis at the observation point. This mode is mainly
  for proxy/engineering configurations and synthetic tests; authoritative instrument
  attitude should normally be supplied as global SM or GSM vectors.

The boresight is explicitly a detector **LOOK** vector. It is **not** the velocity used
to start the backward particle trace. For an incoming particle at the detector, the
AMPS directional-map arrival vector is antiparallel to the telescope look vector. Thus
an aperture row with `b = (0,1,0)` means that the detector looks toward `+Y` in the
selected frame, while the corresponding incoming-particle arrival/velocity direction is
`-Y`. C19 performs this sign conversion internally; aperture files must contain the
physical detector look vector and must not pre-reverse it.

For GOES, identifiers such as `E` and `W` should be interpreted as telemetry/instrument
head names only. They do **not** tell AMPS that a head points geographically east or
west. At every epoch, the actual pointing is determined only by the corresponding
boresight/up vectors. The two heads may be non-antipodal and may change orientation with
time.

The solver does **not** create a new angular lattice for optimized coverage. It first
defines exactly the same regular SM longitude/latitude grid as `FULL_SPHERE`, then keeps
only cells whose detector-look vector is inside at least one configured ellipse. AMPS
directional-map vectors are incoming particle **arrival** directions, so the selector
uses `look = -arrival`, the same convention used by the C19 postprocessor. Retained
cells therefore have exactly the same direction and trajectory as the corresponding
cells in a full-sphere calculation.

For current P4/P5 C19 runs, the pruning envelope uses the widest relevant P5 response
(`30°` horizontal and `60°` vertical half-angle) for each observed detector head. The
actual P4/P5 detector fold still uses the channel-specific response/field-of-view data;
the 30°×60° values only determine which AMPS directions can safely be omitted. With two
roughly opposite GEO heads this typically retains about 1,700--1,900 of the 10,512
2.5° sky cells, but the exact number now follows the actual attitude vectors and need
not be symmetric.

The runner exposes the choice directly:

```bash
# Fast/default science coverage: use the actual per-epoch instrument look vectors
python3 run_C19.py ... --direction-coverage INSTRUMENT_APERTURES

# Complete sky for diagnostic/full-map comparisons
python3 run_C19.py ... --direction-coverage FULL_SPHERE
```

The same generic selector is available directly in an AMPS input file. Apertures can be
provided through a file:

```text
#CUTOFF_RIGIDITY
DIRECTIONAL_MAP T
DIRMAP_LON_RES  2.5
DIRMAP_LAT_RES  2.5
DIRMAP_COVERAGE VECTOR_APERTURES
DIRMAP_APERTURE_FILE instrument_apertures.dat
```

or inline with repeatable records:

```text
DIRMAP_COVERAGE VECTOR_APERTURES
DIRMAP_APERTURE HEAD_A SM  0.10 0.97 0.22   0.00 0.22 0.98  30 60
DIRMAP_APERTURE HEAD_B GSM 0.20 -0.94 0.27  0.01 0.27 0.96  30 60
```

The AMPS command line provides equivalent controls:

```bash
-cutoff-dirmap-coverage VECTOR_APERTURES \
-cutoff-dirmap-aperture-file instrument_apertures.dat
```

and a repeatable inline form:

```bash
-cutoff-dirmap-aperture "HEAD_A SM 0.10 0.97 0.22 0.00 0.22 0.98 30 60"
```

CLI aperture-file selection is resolved before the file is opened, so it can safely
override a path present in a generic input template. `FULL_SPHERE` in either interface
recovers the complete historical sky-map calculation.

For multiple observation locations in one AMPS launch, the solver retains the union of
all cells required by all configured apertures at all locations. This preserves the
constant-size flattened MPI task layout while remaining conservative: no direction that
can contribute at any requested location is pruned.

The legacy `SM_PROXY` orientation remains available in C19 only as a fallback data
source. The runner converts that approximation into ordinary `LOCAL_SM` vector-aperture
records before launching AMPS; there is no EAST/WEST special case in the C++ selector.
With `--detector-orientation-source FILE`, the runner instead writes the actual per-head,
per-epoch boresight vectors to the generated aperture file, so optimized coverage is
fully compatible with time-varying spacecraft attitude.

### GRIDDED and GRIDLESS

GRIDDED and GRIDLESS are now observationally equivalent C19 science paths. Both emit the same three-state `A(E,Ω)` values on the same requested rigidity list and selected directional cells, and both are folded by the identical spectrum/response/aperture postprocessor. The difference is intentionally confined to **magnetic-field evaluation along the trajectory**: GRIDDED interpolates the precomputed Mode3D mesh while GRIDLESS evaluates the configured background model directly. `--solver BOTH` is therefore a true apples-to-apples solver comparison rather than a direct-access-versus-effective-cutoff comparison. The scalar `PENUMBRA_SCAN` maps remain useful diagnostics but no longer define the GRIDLESS observational observable.

### GRIDLESS mesh-free execution and intra-rank threading

Standalone C19 `GRIDLESS` is intentionally **mesh-free**.  The early `-mode gridless`
branch in `srcEarth/main.cpp` initializes MPI and SPICE, parses the input, and calls the
gridless solver directly.  It returns before the historical `amps_init_mesh()` path.
Therefore a standalone T05/T96/T01/TA15/TA16/IGRF/DIPOLE GRIDLESS run does **not**:

- construct/read the AMR tree for the background field;
- allocate Mode3D cell-centered magnetic/electric-field arrays;
- populate a field mesh;
- gather the compact Mode3D field snapshot; or
- interpolate a precomputed field during particle tracing.

Instead, each trajectory step evaluates the selected background field directly at the
particle position.  `PIC::InitMPI()` is still required because it initializes the MPI
runtime; it is not a mesh initialization call.  The live SWMF-coupled build is a separate
case: its "gridless" field evaluator samples the SWMF field stored on the coupled AMPS
mesh and is therefore not the standalone mesh-free path used by C19.

C19 now uses the same two-level parallel design for GRIDLESS cutoff/direct-access work
that Mode3D uses for GRIDDED work:

```text
MPI ranks
  +-- rank/main thread: fetches MPI work chunks; performs MPI/progress/output updates
  |     +-- worker 0: independent trajectory task, private cFieldEvaluator
  |     +-- worker 1: independent trajectory task, private cFieldEvaluator
  |     +-- ...
  |     +-- worker N-1: independent trajectory task, private cFieldEvaluator
  +-- next MPI rank ...
```

The input-file controls are:

```text
GRIDLESS_PARALLEL             THREADS
GRIDLESS_THREADS              16
GRIDLESS_MPI_SCHEDULER        DYNAMIC
GRIDLESS_MPI_DYNAMIC_CHUNK    0
```

and the equivalent CLI is:

```bash
-gridless-parallel THREADS \
-gridless-threads 16 \
-gridless-mpi-scheduler DYNAMIC
```

`GRIDLESS_MPI_DYNAMIC_CHUNK 0` means **automatic**.  The common scheduler receives the
actual number of local workers and chooses a chunk proportional to that count (currently
about four trajectory tasks per worker).  An explicit CLI value may still be supplied,
for example `-gridless-mpi-dynamic-chunk 64`.  A chunk smaller than the worker count can
under-fill the local thread team and is normally a poor choice.

Only the rank/main thread calls MPI.  Worker threads compute complete flattened cutoff
or `(direction,rigidity)` direct-access tasks and write their result into private batch
slots.  After all workers in the batch join, the rank/main thread updates `RcMin`, the
directional-map arrays, the direct `A(E,Omega)` cube, and the MPI progress counter.
Consequently this implementation does **not** require `MPI_THREAD_MULTIPLE` and does not
put locks/atomics around the large result arrays.

The direct-field model interfaces have process-global Geopack/Tsyganenko snapshot state.
For the normal C19 architecture (one spacecraft/epoch per AMPS process), the epoch and
driver snapshot are frozen.  Worker evaluators are created serially, refreshed to that
same snapshot, and the shared model parameters are installed once before the worker team
starts.  If a general GRIDLESS `TRAJECTORY` input contains multiple distinct epochs in
one process, the code conservatively falls back to **SERIAL intra-rank field evaluation**
while keeping the selected MPI scheduler active.  This prevents multiple worker threads
from trying to represent different process-global Geopack epochs simultaneously.

The normal C19 runner requests `THREADS` explicitly for GRIDLESS.  `-nt N` becomes
`-gridless-threads N`.  With the runner's default `--dynamic-chunk 0`, GRIDLESS leaves
chunk sizing to the C++ auto resolver; an explicit positive `--dynamic-chunk` overrides
that behavior.

### Mode3D mesh

The current GRIDDED defaults are:

```text
--mode3d-mesh-res-earth-re 0.025
--mode3d-mesh-res-boundary-re 1.0
```

These are the finest settings from the previous mesh-convergence ladder. They may be
overridden explicitly for a new convergence study, but there is no separate convergence
runner.

### Detector attitude

The historical `SM_PROXY` basis remains the fallback when no authoritative attitude
file is available. For publication geometry use:

```bash
--detector-orientation-source FILE \
--detector-orientation-file /path/to/epead_orientation.csv
```

Preferred columns are one row per detector head and epoch:

```text
utc,spacecraft,detector,frame,
boresight_x,boresight_y,boresight_z,
aperture_north_x,aperture_north_y,aperture_north_z[,source]
```

`frame` is `SM` or `GSM`; `detector` is the actual instrument/telemetry head ID. For
regenerated GOES-13/15 references, the runner reads the required IDs from
`telemetry_head_east`/`telemetry_head_west`, so an attitude file normally contains rows
named `E` and `W`, not rows whose names assert a geometric direction. The two boresights
are independent and need not be opposite. This matters because the instrument geometry
is whatever each head was actually looking at at that epoch, not a hard-coded
geographic-east/geographic-west pair. Legacy references without telemetry-head fields
fall back to `EAST`/`WEST` IDs, and the old one-row `east_boresight_*` attitude schema is
still accepted only for backward compatibility.

Small explicit `--orientation-yaw-deg` and `--orientation-pitch-deg` offsets can be
applied to quantify pointing sensitivity while keeping the same production workflow.
The same perturbed vectors are used both for AMPS aperture pruning and for the detector
fold, preventing the optimization from selecting a different geometry than the science
calculation.

### Upstream anisotropy

The default source is isotropic. A bounded first-order directional model may be selected
explicitly:

```text
J(E,Ω) = J0(E) [1 + A (u · Ω)],    |A| < 1
```

using `--anisotropy-model DIPOLE`, `--anisotropy-amplitude`, and the directional-map
axis longitude/latitude options. The modeled E/W ratio is formed from the synthetic
detector signals, not merely normalized transmissions, so a selected anisotropic source
is represented consistently.

## 10. Outputs

The normal top-level output directory defaults to:

```text
test_output/C19_goes_epead_ew/
```

Per-run directories remain:

```text
<solver>/<field-model>/<spacecraft>/<UTC-token>/
```

Normal machine-readable products include:

| Product | Contents |
|---|---|
| `C19_commands.json` | Exact launch commands and working directories |
| `C19_reference_used.csv` | Selected observational rows plus actual detector IDs, exact flux-variable/correction-state, and ephemeris provenance |
| `C19_spectrum_used.csv` | Epoch-dependent gamma/J0/E0 and measured/interpolated/file source |
| `C19_detector_response_used.csv` | Exact response intervals/components used by direct-response |
| `C19_access_energy_grid.csv` | Energy and proton rigidity values requested identically from GRIDDED and GRIDLESS for direct `A(E,Ω)` |
| per-run `cutoff_3d_dir_access_loc_000000.dat` / `cutoff_gridless_dir_access_point_0000.dat` | Solver-native direct three-state `A(E,Ω)` cubes consumed by the common detector fold |
| `C19_model.csv` / `C19_comparison.csv` | E/W results, transmission bounds, spectrum source, response model, `DIRECT_A_E_OMEGA` access-product label, unresolved fractions, maximum trace time, search algorithm/policy, and map provenance |
| `C19_metrics.csv` | Per-spacecraft and aggregate finite/saturated fractions, sign agreement, bias, MAE, RMSE, correlation, and provisional gate |
| `C19_direction_sense_diagnostic.csv` | Production arrival→look convention and legacy opposite-convention diagnostic |
| `C19_aperture_samples.csv` | Representative aperture cells including lower/effective/upper cutoff, topology, raw termination counts, trace maxima, and transmission |
| `C19_result.json` | Staged validity gates, thresholds, hashes, failures, limitations, and overall result |
| `C19_summary.txt` | Human-readable execution/trajectory/fold/observation/overall status |

Every completed C19 run also writes the standard visual products without any special
flag:

| Plot | Contents |
|---|---|
| `C19_comparison_<solver>_<field>.png` | Observed and modeled log10(E/W) versus time |
| `C19_scatter_<solver>_<field>.png` | Data-ranged observed-versus-modeled scatter |
| `C19_parity_<solver>_<field>.png` | Common-range parity view with the 1:1 line |
| `C19_residual_<solver>_<field>.png` | Model-minus-observation residual versus time |
| `C19_transmission_<solver>_<field>.png` | EAST/WEST aperture transmission diagnostics |
| `C19_aperture_diagnostic.png` | Representative aperture-cell cutoff/access diagnostic |

Exact zero transmission remains a dedicated saturation state rather than an arbitrary
finite substitute. Unresolved cells are carried through lower/upper bounds and the
configured unresolved-aperture validity gate.

## 11. Acceptance behavior

The observational thresholds remain provisional:

```text
finite modeled fraction      >= 0.85
correct E/W sign fraction    >= 0.90
correlation                  >= 0.60
mean absolute log10 error    <= 0.20
RMS log10 error              <= 0.30
```

They are evaluated **after** the trajectory/fold validity stages. Sign agreement remains a direction-
convention/coarse-physics diagnostic and cannot compensate for a poorly resolved or
strongly saturated magnitude comparison.

The summary now reports:

```text
execution complete:       PASS/FAIL
trajectory resolution:    PASS/FAIL
instrument fold:          PASS/FAIL
observational validation: PASS/FAIL
overall:                  PASS/FAIL
acceptance enforcement:   ON/OFF
```

`--enforce-acceptance` controls only shell exit policy; it never changes the scientific
labels. Exit status 2 is reserved for input/launch/output/postprocessing failures. A run
that executes successfully but fails trajectory resolution, detector folding, or GOES
agreement is a scientific FAIL and becomes shell exit status 1 only when
`--enforce-acceptance` is requested.



### GRIDLESS live progress reporting

The threaded GRIDLESS solver uses the same **completed-task** progress semantics as
Mode3D.  Progress is not based on tasks merely fetched from the dynamic MPI queue.
Each successful trajectory worker increments a rank-local atomic completion counter;
the rank/main thread polls that counter every 200 ms and transfers newly completed work
to the global MPI RMA completion counter.  Worker threads never call MPI.

The 200-ms polling cadence is intentionally **not** the stdout cadence.  Rank 0 suppresses
identical progress lines and normally prints only after about 0.1% of the global task set
has newly completed.  For unusually slow trajectories, a line may be printed after 10 s
provided at least one additional task has completed.  Consequently a stalled count such
as `Task 0/98897` is printed once, not once per second.  The final 100% line is also
printed exactly once.

ETA is withheld until at least 0.1% of the work (and at least 32 tasks) has completed and
at least 10 s have elapsed.  This avoids meaningless startup estimates based on one or a
few exceptionally expensive trajectories.

A GRIDLESS run therefore prints an initial line immediately and then lines of the form

```text
[Gridless cutoff TRAJECTORY] [rank 0/global over 4 MPI ranks] [####--------------------------------] 12.3%  (LocEq 0/1, Task 1234/9999)  ETA 00:08:42
```

The important detail is that updates occur **while a threaded MPI chunk is still
running**.  With 16 workers and an automatic chunk of roughly 64 tasks, the display no
longer waits for all ~64 trajectories to finish before advancing.  If rank 0 happens to
be tracing an unusually long trajectory, its main thread still polls the global RMA
counter and displays completions reported by the other MPI ranks.  After the assignment
queue is exhausted, rank 0 continues polling until the global completed-task count
reaches 100%, matching the slow-tail behavior of the GRIDDED Mode3D progress display.

This progress mechanism changes only reporting granularity; it does not alter task
scheduling, trajectory physics, result accumulation, or MPI/thread ownership rules.

## 12. Interpretation and limitations

C19A should be interpreted in layers even though it has only one execution path. The
trajectory/cutoff product must first be resolved and numerically trustworthy; the
synthetic detector fold must then be valid; only after those gates should the GOES
observational comparison be interpreted. A successful GRIDDED C19A supports the more
specific claim:

> With the documented event spectrum, detector response, real spacecraft position, and direct three-state geomagnetic access `A(E,Ω)`, AMPS reproduces the measured GOES EPEAD East–West proton-access ratio for the selected event within the stated validation tolerances.

A systematic sign disagreement deserves special attention before model tuning. If an older C19 result has the opposite sign, inspect `C19_direction_sense_diagnostic.csv`: the production `AMPS_ARRIVAL_TO_DETECTOR_LOOK` mapping should be compared with the legacy direct-vector diagnostic. The production minus sign follows from the distinction between incoming particle direction and telescope look direction and should not be treated as a tunable convention.

C19A alone does not support a claim of exact detector count-rate prediction because:

- the committed default response now includes the documented high-energy secondary proton *energy* response for the uncorrected P4/P5 products, but it is still a factorized proxy rather than the full calibrated energy-angle response matrix;
- the published secondary geometrical factors represent side/rear penetrating-particle response integrated over incidence geometry; C19 currently applies those energy weights with the same nominal angular-aperture model, so component-specific out-of-aperture angular response is not yet reconstructed;
- the two detector heads may retain relative calibration differences;
- `OBSERVED_WEST` is an event-derived spectral approximation; an independent upstream spectrum is preferred when available;
- the production default source is isotropic; anisotropy quantifies a bounded dipole-anisotropy sensitivity, but a measured upstream angular distribution is still preferable when available;
- the prompt event onset remains unsuitable when interplanetary anisotropy is large; and
- provisional observational acceptance thresholds require refinement from multiple events.


A publication-quality extension should replace the current factorized primary+secondary response CSV with the best available calibrated **energy-angle** EPEAD response, use an independent upstream spectrum when possible, quantify detector-head intercalibration and spectral uncertainty, add more events/yaw states, and recheck numerical convergence when changing the selected production settings and additional field models (T96/T05/SWMF as appropriate). Exact ephemeris ingestion, direct response folding, angular/mesh convergence machinery, solver cross-checks, attitude-file support, and anisotropy sensitivity are now implemented rather than deferred; the committed legacy reference remains a nominal-slot artifact and must be regenerated with NOAA ephemeris for a science run.

## 13. Troubleshooting

### Reference file is missing

```text
C19A reference is missing ...
```

Run:

```bash
python3 srcEarth/test/C19/build_goes_reference.py --download \
  --goes13-ephemeris /path/to/goes13_orb1m.nc \
  --goes15-ephemeris /path/to/goes15_orb1m.nc
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

With the default `DIRECT_ACCESS` mode, a directional cutoff-map file is intentionally **not** written and is not required by the runner; aperture geometry is reconstructed from the direct-access cube. If `--cutoff-search PENUMBRA_SCAN` was explicitly requested, verify that the executable supports `DIRECTIONAL_MAP T`, that `DIRMAP_LON_RES` and `DIRMAP_LAT_RES` are positive, and that the cutoff calculation completed.

### Modeled E/W has the opposite sign at nearly every epoch

1. Inspect `C19_direction_sense_diagnostic.csv`.
2. Compare `AMPS_ARRIVAL_TO_DETECTOR_LOOK` with `LEGACY_DIRECT_DIAGNOSTIC` in `C19_summary.txt`.
3. Confirm that the production result uses `AMPS_ARRIVAL_TO_DETECTOR_LOOK`.
4. Check the event-specific telemetry-head-to-physical-look-direction mapping independently.
5. If an older C19 result matches only `LEGACY_DIRECT_DIAGNOSTIC`, that is evidence that the older postprocessor compared AMPS incoming particle direction directly with the telescope look direction.

Do not swap the GOES EAST/WEST labels to compensate. The AMPS-arrival-to-look minus
sign is part of the detector geometry and is now exercised by the runner self-test.

### Many rows report `ZERO_WEST_TRANSMISSION` or `ZERO_EAST_TRANSMISSION`

This is no longer treated as missing aperture coverage. Inspect the per-cell termination diagnostics and
inspect `C19_aperture_samples.csv` and the extended directional-map termination counters. A high apparent cutoff should not be interpreted until `DISTANCE_LIMIT`/`TIME_LIMIT`/`STEP_LIMIT` terminations have been separated from physical loss/trapping and trace-budget sensitivity has been reviewed. If the production
PENUMBRA/UNRESOLVED map is well resolved, then continue with angular/rigidity/spectrum
and instrument-response sensitivity rather than relabeling trace limits as forbidden.

### Many rows report `EXCESSIVE_UNRESOLVED_*`

Inspect the EAST/WEST unresolved solid-angle fractions and the per-cell termination
counts. A large `n_distance_limit`, `n_time_limit`, or `n_step_limit` means the current
finite trace budget did not establish physical access for those rigidity samples. Rerun the same workflow with a controlled larger trace budget to determine whether the unresolved fraction decreases robustly.
Do not convert the unresolved cells to forbidden merely to obtain a finite E/W ratio.

### Rows report `STATIC_FIELD_TRACE_GUARDRAIL`

At least one directional penumbra scan required a trajectory longer than
`--frozen-field-warning-seconds`. This is separate from geomagnetic access. Review the
trace-budget experiment and the physical validity of a frozen T05 epoch before increasing
the guardrail.

### Too many launches

Use `--profile SMOKE`, a larger `--time-step-minutes`, one spacecraft, one model, or one solver during debugging.

## 14. References and public sources

1. NOAA National Centers for Environmental Information, *NOAA GOES 1–15 Space Environment Monitor (SEM) L1B & L2 Data, Version 0 (superseded)*, operational 5-minute GOES-13/15 EPEAD `p17ew` subset used by C19A, dataset identifier `gov.noaa.ncei.swx:goes_sem-l1b-l2_swpc`, DSI `2086_01`. NCEI recommends Version 1.0 for general use; C19A retains the fixed historical directional files for reproducibility.
2. Rodriguez, J. V., T. G. Onsager, and J. E. Mazur (2010), “The east–west effect in solar proton flux measurements in geostationary orbit: A new GOES capability,” *Geophysical Research Letters*, 37, L07109, doi:10.1029/2010GL042531.
3. Rodriguez, J. V., J. C. Krosschell, and J. C. Green (2014), “Intercalibration of GOES 8–15 solar proton detectors,” *Space Weather*, 12, 92–109, doi:10.1002/2013SW000996.
4. Kress, B. T., J. V. Rodriguez, J. E. Mazur, and M. Engel (2013), “Modeling solar proton access to geostationary spacecraft with geomagnetic cutoffs,” *Advances in Space Research*, 52(11), 1939–1948, doi:10.1016/j.asr.2013.08.019.
5. Rodriguez, J. V. et al. (2017), “On the Use of NOAA GOES Energetic Particle Detector Integral Fluxes,” *Space Weather*, 15, 290–309, doi:10.1002/2016SW001533. Table A2 documents the P4/P5 high-energy secondary proton geometrical factors used by the committed uncorrected-response proxy.
6. Tsyganenko, N. A., and M. I. Sitnov (2005), “Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms,” *Journal of Geophysical Research*, 110, A03208, doi:10.1029/2004JA010798.

Machine-readable BibTeX entries are provided in `references.bib`.
