# Change log

## 2026-09-11 — parallel dynamics analysis

- Added `--analysis-workers AUTO|N` to the study orchestrator and standalone
  dynamics reducer, with affinity-aware AUTO selection and an eight-worker
  memory guard.
- Parallelized the expensive lag/bootstrap calculation by complete
  altitude/rigidity/hemisphere series and hysteresis matching by complete
  altitude/rigidity/hemisphere/MLT series.
- Added process-local driver interpolation caches so repeated physical series
  do not repeatedly scan the five-minute TS05 driver table.
- Made bootstrap streams key-derived and output sorting parent-controlled, so
  serial and parallel scientific tables are exactly reproducible.
- Added live phase/series progress, ETA reporting, `analysis_timings.csv`,
  manifest fields, package checks, and serial-versus-parallel regression tests.

## 2026-09-10 — parallel epoch postprocessing

- Parallelized the formerly silent serial interval after each batched AMPS
  launch by assigning complete epoch products, including both shells, to a
  bounded local process pool.
- Added `--postprocess-workers AUTO|N`; AUTO respects CPU affinity, epoch count,
  and an eight-worker memory guard, while `1` preserves a serial reference path.
- Added per-epoch progress messages and `postprocessing_timings.csv` so long
  FULL reductions are observable and worker throughput can be diagnosed.
- Cached AACGM transforms by unique geographic location within each epoch and
  altitude, indexed access rows once for all MLT/rigidity boundary fits, and
  reused the strict multi-shell parse for R50 calculation.
- Preserved deterministic product ordering, process-isolated output paths, the
  existing scientific predicates, and GEO-only global-map behavior.
- Added process-pool, cache-call-count, indexed-versus-legacy parity, worker
  policy, and source-contract regression coverage.

## 2026-09-10 — relative cutoff-erosion maps

- Added cell-resolved maximum relative cutoff decrease, expressed as a
  fraction and percent of the positive BRACKETED quiet-cell R50.
- Added independent per-shell relative extrema to the JSON summary and a
  filled, coastlined publication map in PNG, vector EPS, and PDF formats.
- Kept below-range event estimates explicitly marked as conservative lower
  bounds and masked invalid denominators instead of emitting unstable ratios.
- Extended figure-stage contracts, documentation, and regression tests for the
  new product without requiring another AMPS trajectory calculation.

## 2026-09-10 — GEO-only global-map postprocessing

- Added an explicit `--geo-only` mode to the shared morphology engine and made
  the dedicated global runner select it for SMOKE, ROUTINE, and FULL.
- Removed unnecessary AACGM conversion and ACCESS_T50 boundary fitting from
  complete-shell global-map production, eliminating repeated expected errors
  where AACGM is undefined near the magnetic equator.
- Preserved the observation-facing path unchanged and rejected combinations of
  GEO-only mode with C9/C10 observation staging.
- Added result-contract fields and regression checks proving every global
  profile bypasses AACGM while retaining native batched mesh reuse.

## 2026-09-10 — vector global cutoff maps

- Added a matching vector EPS file for every per-epoch multi-shell global
  cutoff-rigidity PNG.
- Extended the figure manifest and result contract with explicit PNG/EPS paths
  and made a missing epoch EPS a hard figure-stage failure.
- Reused the quiet PostScript backend wrapper so long FULL sequences do not
  flood batch logs with transparency warnings.
- Added regression coverage for EPS creation, nonempty files, and manifest
  traceability while retaining the established PNG-returning Python API.

## 2026-09-10 — live DYNAMIC multi-shell progress

- Enabled once-per-second completed-task progress for Mode3D `DYNAMIC + SHELLS`
  calculations instead of suppressing all intermediate shell updates.
- Kept remote-counter polling outside the hot scheduler-fetch path and deferred
  per-shell details until their final global reduction.
- Prevented nonterminal rounding from displaying a full bar or 100.0%, made the
  terminal ETA exactly zero, and removed duplicate terminal progress lines.
- Added a source-contract regression test covering the dynamic shell path and
  invariants shared with POINTS, TRAJECTORY, STATIC, and BLOCK_CYCLIC runs.

## 2026-09-09 — continuous geographic cutoff maps

- Replaced point-marker global cutoff panels with cyclic filled contours on a
  conventional 180 W--180 E longitude axis.
- Added continental outlines from the existing AMPS
  `earth-continental-map.dat` asset, avoiding new GIS dependencies and runtime
  downloads.
- Preserved masked invalid cells, explicit rigidity-range censoring, canonical
  single-cell poles, and date-line continuity in the visualization contract.
- Applied the same filled-map treatment to the maximum cutoff-decrease figure
  and added regression coverage for the cyclic seam and pole expansion.

## 2026-09-09 — bounded global-map SMOKE workload

- Added an isolated global-map SMOKE override: two quiet/main-phase epochs,
  both shells, a complete 30-degree-by-10-degree geographic grid, and a
  17-point 0.025--20-GV bracket.
- Reduced SMOKE from 1,389,024 to 15,504 fixed-rigidity trajectories while
  retaining native cross-epoch/cross-shell mesh reuse and all downstream
  postprocessing and figure contracts.
- Added pre-launch grid and trajectory-count reporting plus regression guards
  proving that ROUTINE/FULL retain the publication-resolution grid.

## 2026-09-09 — dedicated global cutoff maps

- Added a spacecraft-independent global cutoff-map runner with mandatory native
  Mode3D `SNAPSHOT_LIST` batching and explicit cross-shell/cross-epoch mesh-reuse
  verification.
- Added a small configuration overlay covering the complete GEO sphere at two
  altitudes with a 0.025--20-GV rigidity bracket suitable for polar through
  equatorial LEO cutoffs.
- Added restartable postprocessing that validates grid closure, canonicalizes
  duplicate pole longitudes, preserves censored and unresolved states, and
  derives quiet-relative spatial and temporal erosion diagnostics.
- Added per-epoch two-shell cutoff maps plus maximum-decrease and erosion-
  evolution figures in publication formats.
- Added synthetic full-sphere regression tests, package validation, a data
  dictionary extension, and a dedicated workflow README.

## 2026-09-09

- Added quality-controlled spatial R50 inversion for every shell/epoch
  DIRECT_ACCESS product, including explicit below-range, above-range,
  unbracketed, incomplete, and unresolved-sample diagnostics.
- Added one multi-shell cutoff-rigidity PNG per modeled epoch and a panel/source
  manifest suitable for assembling an event animation.
- Added quiet-referenced spatial cutoff-change and time-evolution products that
  locate the largest erosion, report its epoch/AACGM/MLT context, and preserve
  below-range event values as clearly flagged conservative lower bounds.
- Added publication PNG/EPS/PDF figures for the maximum cutoff-decrease map and
  the time evolution/spatial extent of erosion.
- Expanded dynamics postprocessing with cell-resolved quiet-reference erosion,
  second MLT harmonics, spherical-shell-equivalent accessible area, exact-key
  altitude response, storm extrema, recovery times, and compact best-lag tables.
- Added a machine-readable three-state analysis-availability contract so SMOKE
  spatial diagnostics cannot be mistaken for FULL temporal inference.
- Suppressed bootstrap confidence intervals below 24 paired epochs while
  retaining sparse correlations and matched pairs for deterministic QA.
- Added required publication figures for MLT evolution, altitude dependence,
  and accessible-area expansion in PNG, EPS, and PDF.
- Extended top-level dynamics and figure output checks, unit tests, the data
  dictionary, package manifest, and README for the new scientific products.

## 2026-09-08

- Replaced the unnecessary `TwoSlopeNorm` dependency with an exactly
  equivalent symmetric `Normalize` scale so publication figures run on older
  system Matplotlib installations.
- Converted pandas Series to NumPy/native-datetime plot inputs, fixing the
  old-Matplotlib/new-pandas `value[:, None]` incompatibility, and suppressed the
  repeated harmless PostScript transparency messages during EPS export.
- Fixed BATCHED morphology postprocessing for the current Mode3D
  `fixed_rigidity_access` format, whose generic Tecplot zone identifies each
  shell through the row-level zero-based `shell_index` column.
- Added strict validation and cross-checking of column-, altitude-, and
  zone-based shell identities; malformed or contradictory indices remain hard
  failures.
- Added top-level `--keep` recovery so complete deterministic AMPS snapshot
  batches can be reused after a postprocessing-only failure while all derived
  tables, validation decisions, and figures are regenerated.
- Made an empty optional SMOKE hysteresis table a documented figure skip rather
  than a pipeline exception; required cutoff-degradation PNG/EPS products remain
  enforced by the figure manifest.
- Enabled native Mode3D `SNAPSHOT_LIST` for static `SHELLS` domains and restored
  true cross-epoch mesh-topology reuse. `BATCHED` now runs up to 16 epochs and
  both shells per AMPS process, with deterministic epoch-suffixed products.
- Kept IGRF+TS05/Geopack field reinitialization inside the snapshot loop; only
  time-invariant mesh allocation is reused across epochs.
- Added source-contract and preparation tests that reject mixed installations,
  per-epoch pseudo-batching, missing temporal decks, or ambiguous output names.
- Made RK4 the primary mover in morphology and both independent or staged C9/C10
  command paths.
- Added immediate morphology error reporting and top-level output-contract
  checks so a zero exit without boundary/staged products cannot appear valid.
- Added publication cutoff-degradation heat maps and peak-degradation curves in
  300-dpi PNG, vector EPS, and PDF, plus a machine-readable figure manifest.
- Extended multi-shell Tecplot parsing to accept explicit zero-based shell-zone
  labels while continuing to reject row-order inference.
- Added `BATCHED`, `PER_EPOCH`, and `STANDALONE` morphology layouts. The
  production default combines 475/850-km shells and up to 16 epochs per mesh.
- Added the commented multi-shell input template, strict altitude-zone
  splitting, epoch-suffixed output addressing, batch-complete resume checks, and
  separate logical-product versus MPI-launch counts.
- Added shared-product staging so C9 and C10 reuse the matching 475-km and
  850-km raw access maps without losing their independent observation operators
  or acceptance decisions. Exact observation midpoints are merged into the
  production workset.
- Added `verify_mesh_reuse.py` for exact access-state, unresolved-fraction, key
  closure, and derived-boundary equivalence against `STANDALONE` calculations.
- Added `--independent-observation-runs` for a fully separate C9/C10 numerical
  baseline and documented the new execution and output contracts.
- Added live, line-buffered terminal output for every study stage and for the
  underlying C9, C10, and morphology AMPS/MPI launches while preserving the
  complete per-stage and per-case logs.
- Added completed/remaining/starting counters before every C9, C10, and
  morphology AMPS launch.
- Changed fail-fast handling so PAMELA/C9 and POES/C10 both run as independent
  observational validations; any failure still blocks production morphology
  and inference unless diagnostic continuation is explicitly enabled.
- Added stage start/pass/fail banners, commands, working directories, elapsed
  times, exit codes, a final status table, and timing fields in the run manifest.
- Added early checks for a missing/non-executable AMPS binary and unavailable
  MPI launcher.
- Documented the observational provenance and interpretation of the PAMELA
  Table S1, NOAA POES/MetOp boundary product, and five-minute TS05 driver,
  including the independent archive-rebuild and OMNI comparison audits.
- Set the shared default output root to
  `test_output/dec2006_ts05_cutoff_erosion/`.
- Added repository-root discovery for installation under
  `srcEarth/studies/dec2006_ts05_cutoff_erosion/`.
- Made explicit relative `--output-root` values resolve from the launch
  directory and preserved absolute output paths in generated subcommands.
- Applied the shared default to morphology, comparison, dynamics, figures, and
  TS05 sensitivity runners.
- Added a regression test that enforces the stable output directory name.
