# Steps 13–14 implementation and validation report

Date: 2026-09-13

## Scope completed

Step 13 adds a fail-closed acceptance path to the existing srcSEP test CLI.
Selected tests can emit JSON and JUnit, report-write failures alter the process
status, hidden assertion counters cannot be masked by callback PASS, and MPI
root retains complete rank-local pre-reduction evidence. Four small production
descriptors cover standalone background ownership and cross-mover limits.

Step 14 removes the legacy mover implementation/API surface. The production
archive now contains one adapter per canonical mover plus `mover_state.cpp`;
the mutable public dispatch pointer, retired declarations, old source
monoliths, mover-name aliases, and inactive globals are absent. A migration
manifest and automated documentation/source/archive gates make this boundary
auditable.

## Source-only evidence

The following commands completed successfully from a clean source tree:

```sh
for step in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
  ./test/run_step${step}_tests.sh
done
./test/run_stochastic_repeat_tests.sh
sh -n test/run_step*.sh
make -n test-acceptance-unit
make -n test-documentation-unit
```

Observed acceptance evidence included:

- all prior parser, background, geometry/source, three-mover, field-line scope,
  shared transport, Parker, both focused-transport, coefficient, turbulence,
  and deterministic-reduction checks passing;
- `BG01`, `BG02`, `CROSS01`, and `CROSS02` passing under ASan/UBSan;
- JSON and JUnit files containing every Step 13 acceptance ID;
- the intentional `HIDDEN01` fixture being converted from callback PASS to
  result FAIL while its parent regression test passed;
- `DOC01`, `DOC02`, `DOC03`, `WARN01`, `STATIC01`, and source-only `SAN01`
  passing;
- Steps 7, 8, 9, and 12 producing byte-identical fixed-seed evidence on
  consecutive executions;
- no object, static/shared library, coverage/profiling file, generated result,
  or nested archive present inside the source tree.

## Native integration boundary

The enclosing AMPS configuration is not present in this handoff
(`../../Makefile.conf` is absent), so this environment cannot compile PIC/MPI
adapter translation units or link `amps`. Consequently, the following remain
required in a configured native checkout and are not represented as passed:

```sh
make WARNINGS='-Wall -Wextra -Wpedantic -Werror' lib amps
make test SEP_EXECUTABLE=../amps JSON=results.json JUNIT=results.xml
make test-mpi MPI_NP=4 SEP_EXECUTABLE=../amps \
  JSON=mpi-results.json JUNIT=mpi-results.xml
make test-stress SEP_EXECUTABLE=../amps \
  JSON=stress-results.json JUNIT=stress-results.xml
```

Those commands validate generated PIC types, the real field-line particle
attachment path, MPI collectives/report evidence, SWMF coupling, and extended
legacy diagnostics. The source-only suite deliberately does not emulate or
claim those host-dependent results.

## Packaging contract

The downloadable archive is created one directory above `srcSEP`, preserving
the top-level `srcSEP/` directory. `.gitignore` and `DOC03` prevent build
products, test reports, output directories, and nested delivery archives from
entering future source handoffs.
