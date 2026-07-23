# C6 — global IGRF effective vertical-cutoff validation

C6 compares the AMPS gridless cutoff solver with published global tables of
**effective vertical cutoff rigidity**.  It contains three independently
selectable subtests:

| Subtest | Reference | Bundled coverage |
|---|---|---|
| `INITIAL` | Smart & Shea epoch-2000 printed table | 444 points: 5° latitude × 30° longitude |
| `COMPLETE` | FAA CARI-7 cutoff tables | 386,640 points: six epochs, 1° output grid |
| `MODERN` | Gerontidou et al. 2010/2015/2020 grids | author-supplied 2015 grid bundled; loader supports authoritative 2010/2020 files |

The test uses `FIELD_MODEL IGRF`; no T96, T01, T05, TA15, or TA16 external
field is added.  Every case is run at 20 km geodetic altitude with vertical
incidence.

## Why C6 needs `Rc_effective_GV`

A geomagnetic penumbra may alternate between allowed and forbidden rigidity
bands.  Three cutoff values are therefore distinct:

- `Rc_lower_GV`: first resolved forbidden-to-allowed transition;
- `Rc_upper_GV`: last resolved forbidden-to-allowed transition;
- `Rc_effective_GV`: upper cutoff minus the total width of all allowed bands
  between the lower and upper cutoffs.

C6 compares only `Rc_effective_GV`, because that is the quantity tabulated by
all three references.  Comparing a published effective cutoff with
`UPPER_SCAN` would produce a systematic definition error even if every
trajectory were integrated correctly.

The implementation evaluates the complete `PENUMBRA_SCAN` access sequence,
refines every resolved access transition, and integrates the allowed rigidity
width.  Any `TIME_LIMIT`, `STEP_LIMIT`, or `DISTANCE_LIMIT` sample leaves the
effective cutoff unresolved rather than silently counting it as forbidden.

## Required model support added for C6

C6 relies on these input options:

```text
#BACKGROUND_FIELD
FIELD_MODEL       IGRF
EPOCH             2000-01-01T00:00:00

#CUTOFF_RIGIDITY
CUTOFF_SEARCH_ALGORITHM PENUMBRA_SCAN
CUTOFF_SCAN_SPACING     LINEAR

#OUTPUT_DOMAIN
OUTPUT_MODE       SHELLS
SHELL_ALTS_KM     20
SHELL_LON_RES_DEG 30
SHELL_LAT_RES_DEG 5
SHELL_GEOMETRY    GEODETIC
```

`IGRF` evaluates only `Geopack::IGRF::GetMagneticField`.  `GEODETIC` constructs
WGS-84 positions and local ellipsoid-normal verticals, then rotates both into
GSM at the selected epoch.  `SPHERICAL` remains the default for older tests.

`CUTOFF_SCAN_SPACING` defaults to `LOG`, preserving all earlier behavior.
C6 explicitly selects `LINEAR` because effective-cutoff tables are formed from
rigidity-band widths.  With the defaults `R=[0.01,20] GV` and 2000 vertices,
the spacing is exactly 0.01 GV.

## References and bundled data

### INITIAL: Smart–Shea epoch 2000

`reference_C6_smart_shea_2000.csv` is transcribed from Table 1 of:

D. F. Smart and M. A. Shea, “World Grid of Calculated Cosmic Ray Vertical
Cutoff Rigidities for Epoch 2000.0,” *Proceedings of the 30th International
Cosmic Ray Conference*, Vol. 1, pp. 737–740, 2008.

The printed table has 37 latitudes from +90° to −90° in 5° increments and 12
east longitudes from 0° to 330° in 30° increments.

### COMPLETE: CARI-7

`reference_C6_cari7_1965_2010.csv.gz` contains the numerical appendices from:

K. Copeland, *CARI-7 Documentation: Geomagnetic Cutoff Rigidity Calculations
and Tables for 1965–2010*, FAA report DOT/FAA/AM-19/4, 2019.

It includes epochs 1965, 1980, 1990, 1995, 2000, and 2010, latitudes −89°…+89°,
and all 360 east longitudes.  A few printed/interpolated 1965 entries are
`−0.01 GV`; because a cutoff rigidity cannot be negative, the bundled CSV
clamps those entries to zero and marks their `source` field accordingly.

The archive is complete, but a full calculation is intentionally not the
default.  Six epochs × 64,440 locations × 2000 rigidity nodes would require
hundreds of millions of trajectories.  The default `COMPLETE` run samples a
10° subset.  Use `--full-grid` for the complete one-degree comparison.

### MODERN: Gerontidou et al.

`reference_C6_gerontidou_2015.csv` is converted from the author-supplied 2015
5° × 15° spreadsheet associated with:

M. Gerontidou, N. Katzourakis, H. Mavromichalaki, V. Yanke, and E. Eroshenko,
“World grid of cosmic ray vertical cut-off rigidity for the last decade,”
*Advances in Space Research*, 67, 2231–2240, 2021,
DOI: 10.1016/j.asr.2021.01.011.

The paper describes 2010, 2015, and 2020 grids.  Only an authoritative 2015
machine-readable table was available for bundling during implementation.  C6
does **not** fill the missing years from CARI or interpolate them.  After
obtaining authoritative 2010/2020 tables, convert them to the common schema and
pass them with repeated `--modern-reference` options.

Required columns:

```text
epoch_year,latitude_deg,longitude_deg_east,altitude_km,rc_effective_gv,source
```

## Running the subtests

Initial printed table:

```bash
python srcEarth/test/C6/run_C6.py --subtest INITIAL -np 8 -nt 16
```

CARI-7, one epoch on a 5° subset:

```bash
python srcEarth/test/C6/run_C6.py \
  --subtest COMPLETE \
  --complete-epochs 2010 \
  --complete-grid-step-deg 5 \
  -np 16 -nt 16
```

Complete one-degree CARI epoch:

```bash
python srcEarth/test/C6/run_C6.py \
  --subtest COMPLETE \
  --complete-epochs 2000 \
  --full-grid \
  -np 64 -nt 16
```

Bundled modern table:

```bash
python srcEarth/test/C6/run_C6.py --subtest MODERN --modern-epochs 2015
```

All three modern epochs after obtaining the missing tables:

```bash
python srcEarth/test/C6/run_C6.py \
  --subtest MODERN \
  --modern-epochs 2010,2015,2020 \
  --modern-reference reference_C6_gerontidou_2010.csv \
  --modern-reference reference_C6_gerontidou_2020.csv
```

A mover may be selected explicitly:

```bash
python srcEarth/test/C6/run_C6.py --subtest INITIAL --mover RK4
```

When `--mover` is omitted, the runner passes no mover override and preserves the
current AMPS default.

## Reference-only and dry-run checks

Validate file schemas, uniqueness, and coverage without AMPS:

```bash
python srcEarth/test/C6/run_C6.py --validate-references
```

Render inputs and commands without executing them:

```bash
python srcEarth/test/C6/run_C6.py --subtest INITIAL --dry-run
```

## Output

Each expanded epoch case writes:

- `AMPS_PARAM_C6.in` — rendered, self-describing input;
- `reference_C6_selected.csv` — exact rows used in that case;
- `C6_amps.log` — combined AMPS output;
- `C6_comparison.csv` — point-by-point reference/model comparison;
- `C6_result.json` — case metrics and pass/fail;
- `C6_comparison.png` — optional diagnostic plot when Matplotlib is available.

The root directory additionally contains `C6_summary.csv`, `C6_result.json`,
`C6_commands.json`, and `C6_reference_inventory.json`.

## Acceptance criteria

For each valid point, the allowed error is

```text
max(abs_tol_gv, rel_tol * max(reference_GV, relative_floor_GV))
```

The default thresholds are deliberately validation-scale rather than
bit-for-bit regression thresholds because the references differ in IGRF
generation, trajectory code, Earth geometry, rigidity resolution, and in some
cases spatial interpolation.  C6 requires by default:

- at least 98% matched and resolved points;
- at least 85% of resolved points inside the point tolerance;
- global RMSE no greater than 0.60 GV.

Use the CLI tolerance options to establish tighter project-specific baselines
after the first controlled production runs.  Do not relax tolerances merely to
hide unresolved trajectory terminations; those are reported separately.

## Rebuilding the reference CSVs

`tools/build_reference_tables.py` documents the extraction workflow and can
recreate the CSVs from locally downloaded source PDF/text/XLSX files.  Source
files are not redistributed by this test.
