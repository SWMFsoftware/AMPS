# Generic comet test suite

Run all portable generic-comet tests from the repository root with one public
entry point:

```sh
srcHartley2/tests/run_all_stage_tests.sh
```

Use `--plan` to preflight tools and data without compiling or launching MPI.
Use `--log-root NEW_DIRECTORY` to retain every child log, the human summary,
and the checksum-bound JSON result. The directory must not already exist.

The orchestrator runs the application-independent runtime-input grammar, I3
six-case generation, MAVEN/IUVS reference regeneration, resonance-fluorescence
and brightness tests, and every Stage 0--10 chemistry/transport test. The
`run_stage*.sh` scripts and Python/C++ checkers remain small implementation
components, but they are not separate operator steps. This prevents a partial
selection of individually successful tests from being reported as the model
test result.

## Portable data

All default test evidence is under `data/`. The compact `data/archived` tree
preserves only records consumed by strict checkers: protocol logs, diagnostic
sidecars, focused numerical tables, generated-policy provenance, and compact
contract fixtures. It deliberately excludes executables, particle dumps, full
species output trees, plots, and unrelated production files. The original
multi-gigabyte validation archives remain external scientific evidence.

`data/i3/Archive.zip` is the reviewed source archive for the twelve MAVEN/IUVS
H/OH panels. Keeping it with the tests makes raster-to-table regeneration
mandatory rather than an optional test that silently skips on another host.
`data/SHA256SUMS` authenticates every tracked fixture; the orchestrator checks
it before using any archived record.

Portable protocol fixtures prove that readers, policies, invariants, and
cross-rank result contracts remain compatible. They do not prove a newly built
executable produced those historical records and do not replace a fresh,
converged comet science run. A complete original archive can still be audited:

```sh
srcHartley2/tests/run_all_stage_tests.sh \
  --evidence-root /absolute/path/to/validation_runs
```

The complete archive must have the dated directory layout printed by `--help`.

## Requirements

The default suite requires Bash, GNU make, Perl, Git, SHA-256 utilities,
Python 3 with NumPy and Pillow, C++11/C++17 compilers, `mpicxx`, and an MPI
launcher. The temporary Stage-0 diagnostic runs at 1, 4, and 12 ranks. No
configured AMPS `build/`, assembled root input, AMPS executable, network
download, or production `validation_runs/` directory is required.

The suite never configures AMPS, deletes `build/`, changes an input, or starts
a comet production run. Fresh Stage 4--9 executable modes remain available
through the advanced options printed by `--help`; they require caller-supplied
already-built executables and new output directories.

