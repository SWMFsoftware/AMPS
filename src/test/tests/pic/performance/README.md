# Order-preserving performance-kernel tests

`CORE-PIC-PERFORMANCE-001` is a dependency-light, nonexclusive owner: it compiles
private executables, never configures AMPS or changes a root build. Tool paths,
MPI launcher/rank option and work/artifact directories come from the mandatory
machine configuration. The MPI fixture reduces receiver failures across all
ranks. Every launch must exit successfully and publish exactly one result per
case; missing output is not a PASS.

| Case | What/how/why | Where |
|---|---|---|
| `PERF-PERIODIC-LOCAL-001` | Exact asymmetric byte law through local ordered copies | `test_ordered_exchange.cpp`, production `parallel/ordered_periodic_exchange.h` |
| `PERF-PERIODIC-REMOTE-001` | Every destination checks remote records and byte counts | Same |
| `PERF-PERIODIC-LAYOUT-001` | Reuse workspace after source epochs, pair order and record sizes change | Same |
| `PERF-PERIODIC-CHUNK-001` | Force 17-byte chunks across record boundaries | Same |
| `PERF-PERIODIC-ALIAS-001` | Serial write-to-later-source chain must propagate 17, not stale 99 | Same |
| `PERF-PERIODIC-EMPTY-001` | Empty/zero payload and invalid owner handling | Same |
| `PERF-PERIODIC-COMMUNICATOR-001` | Split/reversed rank groups detect hidden world-communicator assumptions | Same |
| `PERF-STENCIL-EXACT-001` | Nonmonotone IDs, zero/unequal weights and exact insertion order | `test_complete_stencil.cpp`, production `interpolation/complete_stencil.h` |
| `PERF-STENCIL-FALLBACK-001` | Missing/outside/negative support returns without partial writes | Same |
| `PERF-STENCIL-REFRESH-001` | Subsequent query replaces pointers/IDs/weights completely | Same |
| `PERF-SPMV-ORDER-001` | Independent indexed reference and cancellation-sensitive row order | `test_sparse_matvec.cpp`, production `general/ordered_sparse_matvec.h` |
| `PERF-SPMV-REFRESH-001` | Live RHS values and rebuilt coefficients/storage/topology | Same |
| `PERF-SPMV-EMPTY-INVALID-001` | Empty rows/snapshots, negative lengths, null RHS and invalid row handling | Same |

The seven transport identities each run at **1, 2 and 4 ranks**. The six
stencil/sparse identities use a separate one-process harness. They are not
silently counted as three-rank-matrix tests. Expected values are explicit byte
laws or independent indexed arithmetic, not production-generated references.

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --match CORE-PIC-PERFORMANCE-001 \
  --output-dir test_output/performance-kernel-tests
```

Choose a fresh output directory and overlay site paths as needed. The owner is
included in nightly/extended/distributed selection and `--all`. Reports retain
individual statuses plus `performance-kernels.log` with rank counts and child
results. The global owner's log retains exact compiler/launcher commands and
compilation diagnostics. Permanent new cases must update the manifest,
relevant component inventory, README and regenerated `src/test/TEST_CATALOG.md`.
Also register each new compiled case source in
`src/test/source_traceability.json` under `case_registry_sources`, and add its
case/suite identities to the correct ownership area. Run
`python3 -B src/test/validate_traceability.py` and
`python3 -B src/test/generate_test_catalog.py --check`; a case visible only in
the manifest is not completely registered. Agent contributors must follow
these same rules and retain the independent reference equations.
When adding a suite or changing campaign membership, update the deliberately
pinned stream counts and campaign sets in
`src/test/tests/harness/self/run_test.py`. Keep explicit expected sets; deriving
the expected selection from the implementation would erase its rejection
value. Run `CORE-HARNESS-001` after the registry edit.

These tests qualify transport/guarded arithmetic primitives, not all PIC
physics. The integration gates remain the configured periodic/mesh/MPI native
owners, ECSIM operator/boundary owners, new native particle-signature case, and
the unchanged full four-rank fast-wave comparison. See
[`test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md`](../../../../../test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md).
