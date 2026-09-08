# Change log

## 2026-09-08

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
