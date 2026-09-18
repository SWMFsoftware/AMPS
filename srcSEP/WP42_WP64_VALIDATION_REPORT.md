# WP42--WP64 B05 validation and evidence status

## Decision represented by this report

B05 replaces the historical overlay report with current, reproducible status.
The authority for every package is
[`WP42_WP64_DISPOSITION.json`](WP42_WP64_DISPOSITION.json). The four extension
implementations are experimental component APIs and are deliberately absent
from the production `MAINLIBOBJ`. WP59, WP60, WP63, and WP64 are external
qualification gates. None of these work packages is claimed as an active
end-to-end production feature by this source-only report.

Two supporting changes are accepted in current interfaces because they resolve
concrete compatibility/correctness defects without selecting new physics:

1. `SEP::Transport::ParkerMeasure` is defined once in canonical `sep_common`,
   giving both applications the same per-arc-length/per-volume vocabulary.
2. Event-driven focused transport records the selected wave branch and the
   actual momenta immediately before and after each wave-frame scatter. The
   record is deposited only after the full event succeeds, and the production
   adapter no longer substitutes enclosing-shell endpoints for a discrete
   event.

## Dependency-light command and interpretation

Run from `srcSEP`:

```sh
make test-wp42-wp64-experimental
```

The target compiles the experimental sources and their current dependencies
with C++11, `-Wall -Wextra -Wpedantic -Werror`, AddressSanitizer, and
UndefinedBehaviorSanitizer. It then runs one controlled contract for each
WP42--WP64 identifier and executes the machine-readable disposition audit plus
its negative documentation-claim control.

A PASS establishes only that the typed algorithms satisfy their controlled
fixtures and that source/build/runner classification is internally consistent.
It does not establish that the experimental objects are linked, selected, or
executed by AMPS. In particular, the test's WP59/WP60/WP63/WP64 component
checks validate fail-closed evidence schemas; they are not substitutes for the
external evidence named by those schemas.

The deterministic numerical reference metrics retained by the component test
are:

| Work package | Metric | Expected current record |
|---|---|---:|
| WP53 | Brownian parent-minus-children residual | `0` |
| WP53 | accepted / rejected / maximum depth | `3 / 2 / 2` |
| WP54 | coarse L1 error | approximately `0.00385425` |
| WP54 | fine L1 error | approximately `0.00122031` |
| WP54 | observed refinement order | approximately `1.6592` |

These values are component regression evidence, not a scientific validation
of the production model.

## What the component suite covers

- WP42--WP46: prototype turbulence ownership/checkpoint, atomic transaction,
  derivative, shock-state, and source-key contracts.
- WP47--WP52: explicit Parker measure, dynamic resonance, named 90-degree
  closure, wave invariant, conservative cascade, and versioned-profile APIs.
- WP53--WP56: same-path adaptive SDE, periodic MUSCL/SSPRK2 fixture,
  conservative overlap remap, and typed transaction/failure primitives.
- WP57--WP58: lineage-aware estimator and multidimensional detector response.
- WP59--WP64: validators for native traces, decomposition signatures,
  refinement fits, preregistered statistics, external manifests, and scaling
  records, including negative fixtures that reject insufficient evidence.

The disposition audit additionally requires all 23 ordered IDs, reasons and
promotion gates; verifies the experimental objects are outside the production
object list; verifies both suites are selectable from the unified runner; and
rejects unqualified historical production/pass claims.

## External evidence remains separate

Select the native/external gate explicitly:

```sh
python3 test/run_tests.py --suite wp59-wp64-native \
  --plot none --output-dir test_output/wp59-wp64-native
```

Without `SRCSEP_NATIVE_GATE`, the wrapper emits the exact intentional-skip
marker and the unified runner records `SUITE-WP59-WP64-NATIVE` as `SKIP`. That
is acceptable for a development source-only profile but leaves native,
parallel, coupled/observational, and scaling qualification incomplete.

At a configured site, run:

```sh
SRCSEP_NATIVE_GATE='/reviewed/site/runner --matrix --restart --mpi' \
python3 test/run_tests.py --suite wp59-wp64-native \
  --plot none --output-dir /evidence/wp59-wp64-native

make test-swmf-validation SWMF_MANIFEST=/evidence/swmf-replay.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS='/evidence/event-1.json /evidence/event-2.json'
```

Promotion requires executable/compiler checksums, resolved configuration,
callback traces, decomposition and restart signatures, real input lineage,
preprocessing/exclusion records, uncertainty, and frozen hardware/workload
metadata as applicable. Those artifacts must be reviewed before changing any
entry in the disposition manifest to `integrated`.
