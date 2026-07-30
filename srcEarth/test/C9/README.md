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
| Cutoff evaluation | `FULL_SCAN` by default; optional GRIDDED `DIRECT_ACCESS` |
| Backtrace convention | reversed velocity and charge |
| Trace-limit policy | forbidden |
| Geographic shell spacing | 30 deg longitude x 2 deg latitude |
| Routine temporal sampling | one field snapshot at each Table-S1 midpoint |
| Solver selection | `BOTH` for `FULL_SCAN`; `GRIDDED` for `DIRECT_ACCESS` |
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



## 3. Numerical products and solver branches

C9 separates the **scientific boundary definition/product** from the magnetic-
field evaluation branch.

### `FULL_SCAN` - backward-compatible default

`FULL_SCAN` preserves the original C9 implementation. At every node of the
complete configured shell, AMPS evaluates a `PENUMBRA_SCAN` with
`--cutoff-scan-n` rigidity samples. The postprocessor locates the
`Rc_effective = R` contour for every PAMELA rigidity.

`FULL_SCAN` supports `GRIDLESS`, `GRIDDED`, and `BOTH`. It is the reference
product for penumbra topology, `Rc_lower`, `Rc_effective`, and `Rc_upper`
sensitivity, and cross-solver mesh comparisons.

Raw products:

```text
GRIDLESS: cutoff_gridless_shells_penumbra.dat
GRIDDED : cutoff_3d_shells_penumbra.dat
```

### `DIRECT_ACCESS` - efficient GRIDDED validation product

`DIRECT_ACCESS` launches `RIGIDITY_LIST` in standalone Mode3D. It traces exactly
the seven geometric-center rigidities derived from the checked reference CSV,
rather than constructing a 160-node rigidity scan at every location. It also
retains only shell nodes satisfying

```text
access_abs_lat_min <= |geodetic latitude| <= access_abs_lat_max
```

while keeping both hemispheres and all configured longitudes. The default band
is 35-75 degrees. For the default 30-degree x 2-degree shell this retains 480
of 1092 locations. The nominal trajectory count per snapshot changes from

```text
FULL_SCAN:     1092 locations x 160 scan nodes = 174,720 trajectories
DIRECT_ACCESS:  480 locations x   7 rigidities =   3,360 trajectories
```

which is approximately a 50-fold reduction before differences in trajectory
termination time are considered.

For each longitude, rigidity, and hemisphere, rows are ordered by absolute
AACGM latitude. The modeled boundary is the first resolved transition moving
poleward from `PHYSICAL_FORBIDDEN` to `ALLOWED`. Because the two sampled
transmission values are 0 and 1, linear interpolation places the `T=0.5`
boundary halfway between the adjacent latitude nodes. Unresolved nodes are not
reclassified and cannot form a boundary bracket.

`DIRECT_ACCESS` is intentionally GRIDDED-only. Selecting it with `GRIDLESS` or
`BOTH` is rejected explicitly. Its raw product is:

```text
cutoff_3d_shells_access.dat
```

The file uses one row per `(shell node, rigidity)` and records numeric access
state (`0=PHYSICAL_FORBIDDEN`, `1=ALLOWED`, `2=UNRESOLVED`) together with
redundant `allowed` and `unresolved` flags that the C9 parser cross-checks.

### Field-evaluation branches

For a selected product, C9 applies the same event times, driver, trajectory limits,
AACGM conversion, aggregation, and PAMELA acceptance metrics to the supported
field-evaluation branches.

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

`GRIDDED` writes `cutoff_3d_shells_penumbra.dat` for `FULL_SCAN` or
`cutoff_3d_shells_access.dat` for `DIRECT_ACCESS`.

The production starting mesh is 0.02 Re near Earth and 2.0 Re at the outer
boundary with linear coarsening. A 0.01-Re near-Earth rerun is recommended for
mesh-convergence work.

`--solver BOTH` remains the default for `FULL_SCAN`. Each branch is independently
judged against PAMELA. When both are run together, C9 also writes a
GRIDDED-minus-GRIDLESS diagnostic; it does not replace either observational
pass/fail decision. `DIRECT_ACCESS` requires `--solver GRIDDED`.

## 4. Postprocessing and boundary extraction

### FULL_SCAN boundary

For `FULL_SCAN`, each branch evaluates the complete 475-km geographic shell and
returns `Rc_lower_GV`, `Rc_effective_GV`, `Rc_upper_GV`, and penumbra topology.
The C9 comparison uses `Rc_effective_GV`; lower and upper values remain
available for sensitivity diagnostics. The postprocessor reads columns by name
from the Tecplot `VARIABLES` record because GRIDLESS and GRIDDED have different
numbers of solver-specific diagnostic columns.

For each longitude and AACGM hemisphere, rows are sorted from low to high
absolute AACGM latitude. C9 finds the first poleward crossing

```text
Rc_effective > R  ->  Rc_effective <= R
```

and linearly interpolates its latitude.

### DIRECT_ACCESS boundary

For `DIRECT_ACCESS`, the raw rows already contain one binary/three-state access
classification at each requested rigidity. For each rigidity, longitude, and
hemisphere, C9 finds the first resolved poleward transition

```text
PHYSICAL_FORBIDDEN (T=0)  ->  ALLOWED (T=1)
```

and places the `T=0.5` boundary halfway between those adjacent AACGM latitude
nodes. A state of `UNRESOLVED` breaks the bracket. More than one resolved state
change in a profile is counted as a nonmonotonic diagnostic.

### Geographic-to-AACGM conversion

For every retained AMPS shell point, `run_C9.py` calls `aacgmv2` at the exact
snapshot UTC and 475-km altitude. This is required because Table S1 is in AACGM
latitude, not geographic or GSM latitude. AACGM conversion failures near the
magnetic equator are ignored; the validation boundary is in the mid/high-
latitude band.

Install the dependency with:

```bash
python3 -m pip install -r srcEarth/test/C9/requirements.txt
```

For either product, the primary modeled cutoff latitude for one snapshot is the
**median** of all valid longitude/hemisphere crossings. Output also retains the
mean, standard deviation, range, north and south medians, north-south
difference, crossing coverage, and nonmonotonic-profile count.

For `--interval-samples N > 1`, the modeled PAMELA interval is the arithmetic
mean of the N snapshot medians. Each generated command includes an explicit
`--epoch`, and each sample also has a self-contained input file with the same
`EPOCH` directive.

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

## 8. Commands

Run the efficient GRIDDED product with five samples per PAMELA interval:

```bash
python3 srcEarth/test/C9/run_C9.py --solver GRIDDED --cutoff-evaluation DIRECT_ACCESS --profile ROUTINE --interval-samples 5 --access-abs-lat-min-deg 35 --access-abs-lat-max-deg 75 --output-root test_output/C9_direct --amps /home/vtenishe/T11/AMPS/amps --shell-lon-res-deg 30 --shell-lat-res-deg 2 --dynamic-chunk 64 -np 4 -nt 16
```

Run the backward-compatible complete GRIDDED penumbra scan:

```bash
python3 srcEarth/test/C9/run_C9.py --solver GRIDDED --cutoff-evaluation FULL_SCAN --profile ROUTINE --interval-samples 1 --output-root test_output/C9_full --amps /home/vtenishe/T11/AMPS/amps --cutoff-scan-n 160 --shell-lon-res-deg 30 --shell-lat-res-deg 2 --dynamic-chunk 64 -np 4 -nt 16
```

Run both field-evaluation branches with the full scan:

```bash
python3 srcEarth/test/C9/run_C9.py --solver BOTH --cutoff-evaluation FULL_SCAN --profile ROUTINE --cutoff-scan-n 160 --shell-lon-res-deg 30 --shell-lat-res-deg 2 -np 4 -nt 16
```

Preview the exact per-sample commands and generated inputs without launching
AMPS:

```bash
python3 srcEarth/test/C9/run_C9.py --solver GRIDDED --cutoff-evaluation DIRECT_ACCESS --profile SMOKE --interval-samples 5 --dry-run --output-root test_output/C9_direct_dry -np 4 -nt 16
```

Use `--skip-run --keep` with the same product and output root to reprocess
existing raw files without recalculating trajectories.

## 9. Numerical controls

Useful options are:

```text
--solver GRIDLESS|GRIDDED|BOTH
--cutoff-evaluation FULL_SCAN|DIRECT_ACCESS
-np 4
-nt 16
--scheduler DYNAMIC
--dynamic-chunk 0
--mode3d-mesh-res-earth-re 0.02
--mode3d-mesh-res-boundary-re 2.0
--mode3d-mesh-coarsening LINEAR
--mode3d-mesh-exponent 1.0
--cutoff-scan-n 160                 # FULL_SCAN only
--access-abs-lat-min-deg 35          # DIRECT_ACCESS only
--access-abs-lat-max-deg 75          # DIRECT_ACCESS only
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

# FULL_SCAN rigidity scan
--cutoff-scan-n 160
--cutoff-scan-n 320

# DIRECT_ACCESS latitude band and resolution
--access-abs-lat-min-deg 35 --access-abs-lat-max-deg 75
--shell-lat-res-deg 2
--shell-lat-res-deg 1

# Product-definition comparison
--cutoff-evaluation DIRECT_ACCESS
--cutoff-evaluation FULL_SCAN

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
and snapshot directories.  In `C9_comparison.png`, PAMELA and AMPS use the
same color for the same rigidity; dashed circles identify PAMELA and solid
x-marked curves identify AMPS.  PAMELA Table S1 statistical uncertainties are
shown as error bars.  A snapshot contains:

```text
AMPS_PARAM_C9.in
C9_amps.log
cutoff_gridless_shells_penumbra.dat  # GRIDLESS FULL_SCAN
cutoff_3d_shells_penumbra.dat        # GRIDDED FULL_SCAN
cutoff_3d_shells_access.dat          # GRIDDED DIRECT_ACCESS
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
- `FULL_SCAN` uses `Rc_effective = R`; `DIRECT_ACCESS` uses an adjacent-node `T=0.5` crossing;
- default temporal sampling uses the interval midpoint only;
- longitude and both hemispheres are aggregated instead of following the orbit;
- `DIRECT_ACCESS` assumes the selected geodetic latitude band brackets every cutoff;
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
