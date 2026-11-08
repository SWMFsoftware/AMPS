# DIM1 comet-result visualization

## 3I/ATLAS brightness comparison

`plot_i3_brightness_comparison.py` overlays the in-run, particle-resolved AMPS
H Ly-alpha or OH(0,0) brightness on the corresponding numerical MAVEN/IUVS
reference in `references/i3_brightness_MAVEN/maven_iuvs_brightness.csv`.
That table contains all six cases and is the default; ordinary plotting no
longer requires `Archive.zip`.

The observations were supplied as raster plots rather than original numerical
tables.  The reference directory therefore also retains reviewed source-pixel
anchors, source/member hashes, detected plot boxes, axis calibration, and an
explicit uncertainty policy.  Horizontal bounds use the independently
tabulated case resolution.  Visible vertical uncertainty bars are retained;
blank vertical bounds mean that the raster did not resolve a separate bar,
not that physical uncertainty is zero.  The archive's red legacy-model curve
is not part of the extracted reference.

If the AMPS curve falls below the source panels' 0.09-kR lower edge, the new
plot automatically extends the Y axis by up to two useful decades and marks
that source-raster boundary.  This makes a physically low model result visible
instead of silently clipping it.

Case 1 examples using the first completed sampling window are:

```sh
python3 srcHartley2/visualization/plot_i3_brightness_comparison.py \
  validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/pic.H.brightness.s=0.out=0.dat \
  --case 1 --species H \
  --output validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/i3_case1_H_brightness_comparison.png

python3 srcHartley2/visualization/plot_i3_brightness_comparison.py \
  validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/pic.OH.brightness.s=3.out=0.dat \
  --case 1 --species OH \
  --output validation_runs/i3_case1_rays200_sampling10000_max30000_indexfix_20261107_12mpi/production/PT/plots/i3_case1_OH_brightness_comparison.png
```

The case can normally be inferred from a standard run-directory name such as
`i3_case1_*`, and species can be inferred from the Tecplot variables.  Supply
`--case` and `--species` explicitly for renamed or relocated outputs.  Run
`--help` for axis overrides.  After saving the PNG, the script opens the
comparison in an interactive Matplotlib window by default; close that window
to return to the shell.  Use `--no-show` on a batch node.  If neither
`DISPLAY` nor `WAYLAND_DISPLAY` is available, the script saves the PNG and
reports that no window can be opened.

An alternate reference location must contain the same
`maven_iuvs_brightness.csv` schema and is selected with:

```sh
python3 srcHartley2/visualization/plot_i3_brightness_comparison.py \
  PT/plots/pic.H.brightness.s=0.out=0.dat \
  --case 1 --species H \
  --reference-dir /full/path/to/i3_brightness_MAVEN
```

Reproduce the committed table from the original supplied archive with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  srcHartley2/visualization/extract_i3_brightness_maven.py \
  --archive /full/path/Archive.zip
```

Extraction details and limitations are in
`references/i3_brightness_MAVEN/README.md`.  `Archive.zip` is needed for that
audit/regeneration command only, not by the comparison plotter.

The comparison plotter requires Python 3, NumPy, and matplotlib.  The
reference extractor additionally requires Pillow.  Validate the committed
table, its recorded hashes, alternate-reference loading, and byte-exact
archive regeneration with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  srcHartley2/tests/test_i3_brightness_reference.py
```

Without `--output`, the PNG is written under `PT/plots` relative to the
current working directory, using
`i3_caseN_SPECIES_brightness_comparison.png`.  Thus a command launched from
the AMPS repository root writes to that root's `PT/plots`, even when its input
Tecplot table is under `validation_runs`.  Use `--output` when the comparison
must remain beside a frozen validation run instead.

`plot_1d.py` reads the quoted variable names in AMPS Tecplot output instead of
assuming fixed column positions. It plots one selected physical quantity from
one or more files on common axes. The X axis is always the logarithmic
cometocentric radius. The Y scale must be selected explicitly as `linear`,
`log`, or `symlog`.

The script requires Python 3, NumPy, and matplotlib. Saving with `--output`
selects matplotlib's non-interactive backend automatically, so no display is
needed on a batch node.

List the available physical quantities:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  PT/plots/pic.H2O.s=13.out=0.dat --list-quantities
```

Plot water number density with logarithmic axes and radius in kilometres:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  PT/plots/pic.H2O.s=13.out=0.dat \
  -q dens --yscale log --radius-unit km --rmax 1000 \
  --output PT/plots/h2o_number_density
```

Compare the same physical quantity from several output files:

```sh
python3 srcHartley2/visualization/plot_1d.py \
  run-a/pic.H2O.s=13.out=0.dat run-b/pic.H2O.s=13.out=4.dat \
  -q temp --yscale linear --output temperature-comparison
```

`-q` is the short form of `--quantity`. The standard aliases are:

| Alias | AMPS quantity |
| --- | --- |
| `r` | Radius `X` |
| `level`, `id`, `rank` | AMR level, temporary mesh ID, and MPI rank/thread |
| `dens`, `part` | Number density and model-particle number |
| `vel` or `vx`, `vy`, `vz` | Radial/X, Y, and Z velocity components |
| `speed` | Velocity magnitude `\|V\|` |
| `temp`, `tpar`, `ttan` | Total, parallel, and tangential translational temperature |
| `dt`, `weight` | Local timestep and local particle weight |
| `coll-SPECIES` | Collision frequency with `SPECIES`, such as `coll-h2o` or `coll-co2` |

Run `--help` for the complete species-specific alias list and copyable examples
for listing quantities, plotting one file, comparing runs, and batch output.
Exact quantity names remain case-insensitive, and an unambiguous prefix is
accepted. Ambiguous selections are rejected and report the matching columns.

Limit the plotted cometocentric distance with one of these mutually exclusive
controls:

- `--max-distance MAX` or `--rmax MAX` plots from the available inner radius
  through `MAX`;
- `--distance-range MIN MAX` or `--rrange MIN MAX` plots an inclusive radial
  interval.

All distance limits use the unit selected by `--radius-unit`. For example,
`--radius-unit km --rrange 6 100` selects cometocentric radii from 6 km through
100 km. Both limits must be positive and a range must have `MIN < MAX`.

Each curve's legend label is the input file's basename. With logarithmic Y
scale, the script skips non-positive points, reports their count on standard
error, and breaks the curve at each skipped point. Use
`--yscale symlog --linthresh VALUE` when zero or signed values must remain
visible. Without a distance limit, the full radial extent remains visible when
logarithmic Y scaling skips zero-valued outer cells.

If a file has no positive values left for a logarithmic Y plot, its entire
curve is skipped with a diagnostic and the remaining input files are still
plotted. The command fails only when none of the requested files has plottable
data. This permits shell patterns such as `PT/plots/*.dat` when some species
are absent in an early sampling window.

## Rubin et al. (2011) Halley comparison

`compare_rubin2011.py` overlays the calculated neutral-species profiles with
the published Rubin et al. DSMC curves from Figs. 1--17 and overlays the two
water radial-velocity profiles corresponding to Fig. 18. It selects the latest
output index that is complete for every configured neutral species unless
`--output-index` is supplied. Number densities are converted from the AMPS
output unit of m^-3 to the paper's cm^-3, and the paper's species groupings and
neutral-curve line styles are retained.

```sh
python3 srcHartley2/visualization/compare_rubin2011.py \
  PT/plots --run-log halley1d.run.log
```

The default `--paper-data` table is
`srcHartley2/halley2011/reference/rubin2011_neutral_profiles.csv`. It was
extracted from vector paths in the local journal PDF, including the paper's
piecewise logarithmic distance axes. The accompanying metadata JSON records
the exact PDF SHA-256, plot bounds, source pages, curve colors, and extraction
method. Regenerate the table from the same PDF with:

```sh
python3 srcHartley2/halley2011/digitize_rubin2011.py \
  ~/Downloads/Rubin-2011-icarus.pdf
```

The default destination is `PT/plots/rubin2011_comparison`. It contains a PNG
and EPS for each figure, a combined multipage PDF,
`comparison_profiles.csv` at 10, 100, 1,000, 10,000, and 100,000 km, and a
JSON provenance manifest. The CSV contains side-by-side AMPS and Rubin values
at the tabulated radii. The manifest records the exact AMPS files, output
index, nominal iteration, density conversion, convergence record, reference
table path, and reference-table digest.

When `--run-log` is supplied, the script parses the selected window's
`consecutive=N/M` record. Every PNG, EPS, and combined-PDF page is visibly
marked `CONVERGED` or `NONCONVERGED`; missing or unrecognized records are
also labeled rather than silently treated as converged. The manifest stores
the same result in `scientific_status` and adds an explicit diagnostic-only
limitation for every state other than `converged`.

Figs. 2--17 in the paper also contain measured and simulated Giotto NMS mass
channel signals. The current comet application does not model the NMS ion
source or detector, so those different observables are not overlaid with
neutral density. The comparison contains both AMPS and Rubin neutral DSMC
curves and states the NMS limitation in the manifest. Methanol likewise
remains an empty AMPS curve when it is disabled in `halley1d.input`.

## Stage 0 minor-species statistics

`analyze_minor_species.py` is the read-only, machine-readable counterpart to
the plotting tools. It summarizes sampled particle occupancy, zero-cell
fraction, first and last populated radii, block-local `W`, `dt`, and `W/dt`,
window-to-window changes at requested radii, the recorded convergence state,
and Rubin comparison values. It writes both JSON and long-form CSV:

```sh
python3 srcHartley2/visualization/analyze_minor_species.py \
  PT/plots --species C2 C2H C2H2 \
  --output-index 3 --sampling-length 10000 \
  --run-log run.log \
  --comparison-csv PT/plots/rubin2011_comparison/comparison_profiles.csv \
  --output-prefix PT/plots/minor_species_stage0
```

The label `log10(Local Particle Weight)` is retained verbatim because that is
the historical AMPS output-column name; the radial application writes the
actual base weight in that column, and the report records the source label.

The integrated density change uses exactly the same definition as
`RadialSimulation.cpp`: the L1 norm of `n*cell_measure`, normalized by the
larger adjacent-window norm. Legacy plot files do not include the core cell
measure. For those files the JSON/CSV explicitly reports the metric as
unavailable; it does not reconstruct a volume from plotting coordinates.
Future diagnostic output may add `Cell Measure [m^3]`, or an exact aligned CSV
can be supplied with `--cell-measures` using the columns `radius_m` and
`measure_m3`.

The preserved high-statistics C2 acceptance values and artifact checksums are
in `halley2011/reference/halley1d_stage0_baseline.json`. Run all Stage 0 tests
with:

```sh
srcHartley2/tests/run_stage0_tests.sh
```

The runner checks deterministic counter accumulation, identical reductions at
1, 4, and 12 MPI ranks, exact-measure L1 analysis, legacy missing-measure
handling, and the external archived C2 final window. Set
`HARTLEY2_STAGE0_BASELINE_ROOT` when the validation archive is not at its
documented repository-relative location; only the archive test is skipped
when neither location exists.

For an instrumented run, set `HARTLEY2_STAGE0_PRODUCTION_ROOT` to its
production directory. The runner then invokes
`test_stage0_production_report.py --require-observed-c2`. That checker reads
`PT/diagnostics/cell_measures.csv` and every `comet_stage0.out=N.{csv,json}`,
rejects shifted CSV columns or inconsistent reaction/product/flux identities,
requires both configured C2 production paths and an actual C2 carrier, and
recomputes the `max(5 sigma, 1%)` acceptance instead of trusting the stored
boolean. A zero C2 destruction or escape rate is retained as an observed
zero, not replaced by a synthetic event.

`--output NAME` always saves `NAME.png` and `NAME.eps`. If the argument has a
different extension, such as `.pdf`, that format is saved in addition to PNG
and EPS. After saving, the script opens the PNG in the system's desktop image
viewer. Use `--no-show` for a batch job. When no graphical display is
available, the script reports that condition, saves the files, and skips the
viewer automatically. If `--output` is omitted, matplotlib's interactive
viewer is used directly.

## Stage 3 exposure-controlled radial bins

`adaptive_radial_bins.py` reads one species table plus the exact
`PT/diagnostics/cell_measures.csv` written by the radial model. It merges
contiguous cells from the nucleus outward until the requested exposure is
reached, where exposure is mean sampled model-particle occupancy times the
completed sampling length. It never changes the original AMPS table.
The underlying radial mesh already expands logarithmically; aggregate edges
remain exact existing cell faces rather than interpolated synthetic log edges.

Use `stage3_cells.SPECIES.out=N.csv`, not `pic.*.dat`. The latter is a
node-interpolated Tecplot representation with repeated AMR interfaces and
cannot be joined exactly to center-cell shell exposure. The dedicated Stage-3
table is emitted for complete diagnostic sampling windows and contains the
authoritative center-cell density and mean model-particle occupancy. It also
embeds cell index and shell geometry, sampling length, output index, random
seed plus an explicit-seed flag, residual mode, and MPI-rank count. Current
tables additionally embed the block index, local cell index, cells per block,
resolved population target, lock flag, derived population class, and winning
policy source. The tool checks the index arithmetic and requires one
internally consistent policy per block. Legacy tables without any of the
policy columns remain readable for exploratory plots, but the permanent
Stage-3 science evaluator requires the complete authoritative metadata and
the unlocked 50--200 `CHEMICAL_ONLY` policy for every important daughter.
The tool cross-checks embedded geometry against the supplied cell-measure
table and rejects inconsistent provenance instead of interpolating it.

Every CSV/JSON row records the inner/outer support, geometric representative
radius, cell count, volume, exposure, volume-weighted density, integrated
number, mean model particles per cell, nonpositive-cell count, and whether the
last bin is underexposed or completely empty. Zeros are retained.

Repeat `--replicate TABLE` for independent seeds. All replicas must have the
same exact mesh, species, output index, sampling length, residual mode, and MPI
metadata. Multi-replica analysis also requires distinct explicitly configured
seeds; duplicate paths and duplicate seeds are fatal. The CLI sampling length
must match the embedded value. Common bin boundaries are selected from mean
exposure; density standard deviation, standard error, CV, and
nonpositive-replica fraction are then evaluated across independent runs. Four
or more total tables are marked `decision_sample_status=ready`; this only
means enough replicas are present, not that an unconfigured scientific
threshold has passed.

Run the focused tests with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_adaptive_bins.py
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_diagnostic_semantics.py
```

The first command currently runs eleven cases, including duplicate-seed,
sampling-length, optional block-policy preservation, and corrupt block-index
rejection. The second runs four cases proving that the
Bernoulli five-sigma criterion is not reused for systematic residuals while
pre-Stage-3 Bernoulli reports remain readable.

The permanent all-eight-run gate uses
`srcHartley2/validation/evaluate_stage3_science_gate.py`, not eight separate
invocations of this plotting utility. That evaluator derives one common
support from the minimum exposure over all four Bernoulli and all four
Systematic runs, preventing a low-exposure seed from being hidden by an
ensemble mean. See `srcHartley2/validation/README.md` for the fixed thresholds
and archive contract.
