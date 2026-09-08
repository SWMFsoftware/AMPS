# Change log

## 2026-09-08

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
