# C9 - PAMELA public-data global-shell cutoff-latitude validation

## 1. Purpose

C9 is the routinely executable, public-data form of the PAMELA December 2006
storm validation. It tests whether the AMPS time-dependent IGRF+T05/TS05
backtracing calculation reproduces the observed equatorward motion and recovery
of the low-Earth-orbit proton cutoff boundary. The same observable is evaluated
through two explicit solver branches: GRIDLESS direct field evaluation and the
GRIDDED standalone Mode3D mesh/interpolation path.

The observational reference is **Supporting Information Table S1** from:

> Adriani, O., et al. (2016), *PAMELA's measurements of geomagnetic cutoff
> variations during the 14 December 2006 storm*, Space Weather, 14, 210-220,
> doi:10.1002/2016SW001364.

Table S1 states that its values are geomagnetic cutoff **AACGM latitudes**, with
statistical uncertainties, and that each listed date is the midpoint of a
94-minute interval. The checked-in reference CSV contains all 37 intervals and
all seven published rigidity bins; the single missing table cell remains
missing and is never interpolated.

C9 intentionally avoids a claim that public information is sufficient to
recreate PAMELA's exact event geometry. The public test does **not** need:

- definitive Resurs-DK1 ephemeris;
- spacecraft attitude quaternions;
- the PAMELA-to-spacecraft mounting matrix;
- individual PAMELA track directions.

Instead, it compares the published orbit-averaged product with an AMPS global
shell evaluated at the representative PAMELA altitude.

## 2. Fixed physical defaults

The agreed routine defaults are:

| Quantity | C9 default |
|---|---:|
| Shell altitude | 475 km |
| Incidence | vertical protons |
| Rigidity representation | geometric center of each Table-S1 bin |
| Field model | IGRF + T05/TS05 |
| Driver cadence | 5 minutes |
| Cutoff method | `PENUMBRA_SCAN` |
| Backtrace convention | reversed velocity and charge |
| Trace-limit policy | forbidden |
| Geographic shell spacing | 30 deg longitude x 2 deg latitude |
| Routine temporal sampling | one field snapshot at each Table-S1 midpoint |
| Solver selection | `BOTH` (GRIDLESS and GRIDDED branches) |
| Mode3D near-Earth mesh target | 0.02 Re |
| Mode3D boundary mesh target | 2.0 Re |
| Mode3D coarsening | `LINEAR`, exponent 1.0 |

The seven geometric-center rigidities are approximately:

```text
0.423556372 GV
0.498397432 GV
0.587877538 GV
0.692820323 GV
0.817006732 GV
0.962081078 GV
1.131017241 GV
```

The runner calculates the exact centers from the reference CSV rather than
copying rounded values from this README.


## 3. Two AMPS solver branches

C9 applies identical event times, TS05 driver values, shell points, vertical
incidence, rigidity scan, trajectory limits, AACGM conversion, contour
extraction, and PAMELA acceptance metrics to both branches.

### `GRIDLESS`

The executable is launched with `-mode gridless`. IGRF+T05 is evaluated
directly at trajectory positions. The numerical template is:

```text
AMPS_PARAM_C9_gridless.in
```

The raw penumbra product is:

```text
cutoff_gridless_shells_penumbra.dat
```

### `GRIDDED`

The executable is launched with `-mode 3d`, `-mode3d-field-eval MESH`, and the
Mode3D thread/MPI scheduler controls. IGRF+T05 is first sampled onto the
standalone AMR mesh; trajectory field values then come from the Mode3D
interpolation stencil. The numerical template is:

```text
AMPS_PARAM_C9_mode3d.in
```

The raw penumbra product is:

```text
cutoff_3d_shells_penumbra.dat
```

The production starting mesh is 0.02 Re near Earth and 2.0 Re at the outer
boundary with linear coarsening. A 0.01-Re near-Earth rerun is recommended for
mesh-convergence work.

`--solver BOTH` is the default. Each branch is independently judged against
PAMELA. When both are run together, C9 also writes a GRIDDED-minus-GRIDLESS
comparison as a diagnostic of mesh/interpolation effects; that diagnostic does
not replace either observational pass/fail decision.

## 4. What AMPS calculates

Each solver branch evaluates a complete 475-km geographic shell. AMPS returns:

```text
Rc_lower_GV
Rc_effective_GV
Rc_upper_GV
penumbra topology and unresolved-bracket diagnostics
```

in `cutoff_gridless_shells_penumbra.dat` for GRIDLESS or
`cutoff_3d_shells_penumbra.dat` for GRIDDED.

The historical scalar file `cutoff_gridless_shells.dat` contains only one cutoff
number and is not used by C9. The external comparison is based on
`Rc_effective_GV`, while lower and upper cutoffs remain available for diagnosis.

The two penumbra files do not have identical column counts. GRIDLESS includes
three additional Størmer diagnostic columns and normally has 21 columns, while
GRIDDED omits those GRIDLESS-only values and normally has 18. The C9
postprocessor reads the quoted Tecplot `VARIABLES` record and locates the common
fields by name; it must not assume a fixed positional schema shared by both
branches.

### Geographic-to-AACGM conversion

For every AMPS shell point, `run_C9.py` calls `aacgmv2` at the exact snapshot
UTC and 475-km altitude. This is required because Table S1 is in AACGM latitude,
not geographic or GSM latitude.

Install the dependency with:

```bash
python3 -m pip install -r srcEarth/test/C9/requirements.txt
```

The runner stops with a clear error when scientific postprocessing is attempted
without `aacgmv2`. Reference validation and `--dry-run` do not require it.

### Boundary extraction at one rigidity

For each sampled geographic longitude and each AACGM hemisphere, the model rows
are sorted from low to high absolute AACGM latitude. The code finds the first
poleward crossing satisfying:

```text
Rc_effective > R  ->  Rc_effective <= R
```

and linearly interpolates the crossing latitude. The first equator-to-pole
transition is used as the operational cutoff boundary. Additional transitions
are counted as nonmonotonic-profile diagnostics.

For a single snapshot, the primary modeled cutoff latitude is the **median** of
all valid longitude/hemisphere crossings. The output also retains:

- arithmetic mean and standard deviation;
- minimum and maximum crossing latitude;
- separate north and south medians;
- north-south difference;
- number of requested and valid crossings;
- number of longitude/hemisphere profiles with multiple crossings.

For `--interval-samples N > 1`, the model result for a PAMELA interval is the
arithmetic mean of the N snapshot medians. This approximates the 94-minute
observational averaging while giving equal weight to each field snapshot.

## 5. Reference data

`reference_C9_pamela_table_s1.csv` contains 259 rows:

```text
37 interval midpoints x 7 rigidity bins = 259
```

There are 258 numerical observations and one missing value. Columns include:

```text
interval_midpoint_utc
interval_start_utc
interval_end_utc
rigidity_min_gv
rigidity_max_gv
rigidity_geometric_center_gv
pamela_cutoff_aacgm_deg
sigma_plus_deg
sigma_minus_deg
missing
source
source_pdf_sha256
notes
```

The 94-minute start/end fields are midpoint +/-47 minutes. They are metadata for
temporal averaging; the original table prints only the midpoint.

Validate the complete transcription without AMPS:

```bash
python3 srcEarth/test/C9/run_C9.py --validate-references
```

For audit/reproduction, `tools/extract_table_s1.py` recreates the CSV from
`swe20314-sup-0001-supinfo.pdf`. Runtime execution does not depend on Camelot or
on the PDF.

## 6. TS05 driver bundled with C9

C9 uses the checked-in five-minute event driver:

```text
srcEarth/test/C9/data/ts05_driving.txt
```

Its columns are:

```text
UTC Bx By Bz Vx Vy Vz Np Temp SYM-H IMFflag SWflag Tilt Pdyn W1 W2 W3 W4 W5 W6
```

The file covers `2006-12-14T00:00:00` through
`2006-12-17T00:00:00` at five-minute cadence and includes the abrupt solar-wind
and IMF change near 14 December 14:10 UT.  The exact file supplied for C9 is
protected by this SHA-256 digest:

```text
cb3f3f1959763660beb1e26e5a49489b132708944fb91c4e1ee37cfc3a6c4317
```

The runner uses this file automatically; `--driver` is not required.  It checks:

1. timestamp plus 19 numerical fields on every data row;
2. strictly increasing UTC times;
3. nominal five-minute median cadence;
4. no gap larger than ten minutes;
5. complete coverage of every requested field snapshot; and
6. agreement with the checked-in SHA-256 digest.

Validate the bundled driver without running AMPS:

```bash
python3 srcEarth/test/C9/run_C9.py --validate-driver
```

An alternate driver may be selected with either:

```bash
--driver /path/to/driver.txt
```

or:

```bash
export C9_TS05_DRIVER=/path/to/driver.txt
```

An alternate driver is scientifically eligible when it either matches the
bundled checksum or contains provenance produced by
`tools/prepare_official_ts05_driver.py`.  A structurally valid but unverified
alternate driver is accepted only with `--allow-unverified-driver`, and the
result JSON is marked `scientific_validation_eligible=false`.

## 7. Execution profiles

### `SMOKE`

Two Table-S1 intervals:

```text
2006-12-14T14:31:00Z  shock interval
2006-12-15T03:03:00Z  cutoff minimum
```

This verifies input generation, driver interpolation, shell output, AACGM
conversion, contour extraction, and comparison logic. It is not the default
scientific regression.

### `ROUTINE` - default

Seven representative intervals:

```text
2006-12-14T09:49:00Z  pre-shock reference
2006-12-14T14:31:00Z  shock interval
2006-12-14T16:05:00Z  early storm response
2006-12-14T23:55:00Z  main-phase development
2006-12-15T03:03:00Z  observed minimum cutoff latitude
2006-12-15T09:19:00Z  recovery
2006-12-16T04:08:00Z  late recovery
```

One TS05 field snapshot is used at each published midpoint unless
`--interval-samples` is increased.

### `FULL`

All 37 Table-S1 intervals. Use this for publication-quality reproduction or
sensitivity work rather than frequent source regression.

Example with five snapshots per 94-minute interval:

```bash
python3 srcEarth/test/C9/run_C9.py \
  --profile FULL \
  --interval-samples 5 \
  -np 32
```

Custom Table-S1 midpoints can be selected with a comma-separated list:

```bash
--timestamps 2006-12-14T14:31:00Z,2006-12-15T03:03:00Z
```

## 8. Routine command

From the AMPS repository root, run both branches:

```bash
python3 srcEarth/test/C9/run_C9.py \
  --solver BOTH \
  --profile ROUTINE \
  --cutoff-scan-n 160 \
  --shell-lon-res-deg 30 \
  --shell-lat-res-deg 2 \
  -np 4 -nt 16
```

Run branches separately when different MPI layouts are preferred:

```bash
# Direct empirical-field evaluation
python3 srcEarth/test/C9/run_C9.py --solver GRIDLESS \
  --profile ROUTINE --cutoff-scan-n 160 \
  --shell-lon-res-deg 30 --shell-lat-res-deg 2 -np 16

# Mode3D mesh/interpolation evaluation
python3 srcEarth/test/C9/run_C9.py --solver GRIDDED \
  --profile ROUTINE --cutoff-scan-n 160 \
  --shell-lon-res-deg 30 --shell-lat-res-deg 2 \
  --mode3d-mesh-res-earth-re 0.02 \
  --mode3d-mesh-res-boundary-re 2.0 \
  --mode3d-mesh-coarsening LINEAR -np 4 -nt 16
```

Preview all generated commands and numerical inputs without running AMPS:

```bash
python3 srcEarth/test/C9/run_C9.py \
  --solver BOTH --profile ROUTINE \
  --dry-run -np 4 -nt 16
```

By default the runner expects `./amps`. Select another executable with:

```bash
--amps /path/to/amps
```

## 9. Numerical controls

Useful options are:

```text
--solver GRIDLESS|GRIDDED|BOTH
-np 4
-nt 16
--scheduler DYNAMIC
--dynamic-chunk 0
--mode3d-mesh-res-earth-re 0.02
--mode3d-mesh-res-boundary-re 2.0
--mode3d-mesh-coarsening LINEAR
--mode3d-mesh-exponent 1.0
--cutoff-scan-n 160
--rigidity-min-gv 0.30
--rigidity-max-gv 1.35
--shell-lon-res-deg 30
--shell-lat-res-deg 2
--max-trace-time 20
--mover BORIS
```

The rigidity bracket must cover every Table-S1 geometric center with margin.
The default 0.30-1.35 GV bracket gives a linear scan spacing of about 0.0066 GV
for 160 samples.

Recommended convergence checks are:

```bash
# Latitude-grid sensitivity
--shell-lat-res-deg 2
--shell-lat-res-deg 1

# Longitude sampling
--shell-lon-res-deg 30
--shell-lon-res-deg 15

# Rigidity scan
--cutoff-scan-n 160
--cutoff-scan-n 320

# Temporal averaging
--interval-samples 1
--interval-samples 5
--interval-samples 19
```

## 10. Acceptance metrics

The routine defaults are intentionally looser than PAMELA's statistical
uncertainties because C9 is a global-shell approximation, not an exact
spacecraft-acceptance reproduction.

| Metric | Default requirement |
|---|---:|
| Valid modeled/reference fraction | >= 0.85 |
| Overall latitude RMSE | <= 3.0 deg |
| Absolute mean bias | <= 2.0 deg |
| PAMELA/AMPS correlation | >= 0.80 |
| Lowest-bin storm suppression | 4-9 deg |
| Time of lowest-bin minimum | within 100 min |

All thresholds are CLI-configurable. The runner also reports:

- mean absolute and maximum absolute error;
- uncertainty-normalized residuals;
- fraction inside PAMELA 1-sigma and 2-sigma intervals;
- observed and modeled low-rigidity suppression;
- observed and modeled time of minimum;
- longitude/hemisphere and temporal scatter.

The statistical uncertainties are diagnostics, not the sole pass criterion:
they do not include global-shell geometry, rigidity-bin, magnetic-field-model,
or temporal-sampling systematics.

## 11. Output products

Top-level output defaults to `test_output/C9`:

```text
C9_reference_used.csv
C9_driver_info.json
C9_commands.json
C9_result.json
driver/ts05_driver.txt
gridless/
gridded/
C9_cross_solver_comparison.csv   # written for --solver BOTH
C9_cross_solver_result.json      # written for --solver BOTH
```

Each selected branch contains its own `C9_model.csv`, `C9_comparison.csv`,
`C9_comparison.png`, `C9_result.json`, command inventory, interval directories,
and snapshot directories. A snapshot contains:

```text
AMPS_PARAM_C9.in
C9_amps.log
cutoff_gridless_shells_penumbra.dat  # GRIDLESS
cutoff_3d_shells_penumbra.dat        # GRIDDED
C9_snapshot_boundaries.csv
```

The exact driver, generated input, command line, and comparison settings are
therefore archived with every run.

## 12. Interpretation and limitations

C9 answers this question:

> Do the GRIDLESS direct-field and GRIDDED Mode3D-mesh implementations of the
> time-dependent AMPS IGRF+TS05 global vertical-access boundary at 475 km reproduce the rigidity dependence, timing, magnitude, and recovery of
> the published PAMELA AACGM cutoff-latitude product?

It does not answer the stricter question:

> Can AMPS reproduce every individual PAMELA event using exact Resurs-DK1
> position, attitude, detector acceptance, and reconstructed track direction?

That stricter future test is C9C. A map-comparison extension using published
PAMELA global cutoff maps can be treated as C9B.

Other known approximations are:

- one geometric-center rigidity represents each finite PAMELA rigidity bin;
- the primary boundary is the `Rc_effective = R` contour;
- default temporal sampling uses the interval midpoint only;
- longitude and both hemispheres are aggregated instead of following the orbit;
- geodetic vertical is used rather than the exact instrument boresight;
- the GRIDDED branch adds finite mesh resolution and interpolation error, which
  must be assessed with the supplied mesh controls.

These limitations are recorded in `C9_result.json` and should be retained in
reports that use the test.

## 13. Source files

```text
AMPS_PARAM_C9_gridless.in
AMPS_PARAM_C9_mode3d.in
run_C9.py
reference_C9_pamela_table_s1.csv
requirements.txt
tools/extract_table_s1.py
tools/prepare_official_ts05_driver.py
data/README.md
data/ts05_driving.txt
```
