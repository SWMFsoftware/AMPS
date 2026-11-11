# Fast-wave optimization implementation and evidence

Updated 2026-11-10 on passion. P0–P4 production changes are implemented and
the final native fast-wave regression passes. Complete core revalidation and
documented owner retries are finished; all previous passing identities remain
passing and all eighteen added cases pass. Initial registration/preflight
failures and their actual retries are retained below.
No Git state has been changed by this work.

## Scientific and numerical boundary

The driver still evolves 60 iterations with the existing four-rank deck,
species, particle counts, seed, sampling, deterministic ordering, validation
checks, output and solver controls. The existing proton/electron
`DiffNum.pl -a=2.0E-14` comparisons remain. Every successful candidate also
reproduces **all six full sampled species files byte-for-byte**, including
the column excluded by the historical reference comparison. References were
copied privately for measurement, their original hashes are unchanged, and
the original output symlink was restored. Timing is not a scientific gate.

## Production changes and contracts

| Change | Implementation and numerical reasoning | Documentation |
| --- | --- | --- |
| Periodic exchange | Aggregate exact full associated-data bytes by peer, post receives/sends, wait and unpack in original pair order. Eligibility requires default pure-copy callbacks and distinct physical/ghost leaves. Custom or uncertain sources use streaming reference semantics. Only buffer capacities persist; owners, offsets and payload sizes are rebuilt every call. | [`src/pic/parallel/README.md`](../../src/pic/parallel/README.md) |
| Interpolation | Complete eight-node uniform pointer stencils use the existing coefficient expressions, node order and normalization. Missing nodes, adaptive layouts and row stencils retain the general path. Rejected fast paths leave the destination unchanged. | [`src/pic/interpolation/README.md`](../../src/pic/interpolation/README.md) |
| Deposition | Per-call block metadata hoists ownership/periodic-face checks and cell geometry. Uniform-only magnetic-node pointers dereference live values; no field or particle-position weight survives as a cached state. Existing flat cell/particle order and synchronization remain. Other backends retain the original geometry fallback. | [`src/pic/README.md`](../../src/pic/README.md) |
| Signatures | Immutable host tables and four-byte slicing reproduce the historical unsigned-long CRC recurrence, including LP64 upper bits. Whole-record copy then zeroed list fields equals the prior canonical CloneParticle record. Sorting retains one signature per particle, comparator/tie behavior and list validation; descriptors retain capacity only. Device CRC retains the bytewise algorithm. | [`src/general/README.md`](../../src/general/README.md) |
| Sparse multiply | Compact scalar host row offsets, coefficients and live RHS pointers exist only within the linear-solver callback scope. Matrix/vector-generation or coefficient changes refresh the view. Scalar summation order, input copy for alias protection, MPI, tolerances and iteration ceilings remain. CUDA, AVX and std::thread branches keep their existing kernels. | [`srcInterface/README.md`](../../srcInterface/README.md) |

All new implementation paths have comments describing dependencies, arithmetic
order, lifetime, failure conditions and fallback boundaries. This bounded
optimization does not establish new CME, GPU, AVX or adaptive-model physical
qualification.

## Native measurements

Every candidate used the exact command `make test_fast-wave -j` from the
idle AMPS root. The maintained recipe verifies root/idle state, removes only
generated `build/`, configures the application through the normal workflow,
regenerates hooks and compiles through the top-level workflow with `-j16`.
Site configuration and libraries are preserved. The compiler policy remains
the previously measured GNU-only `-O1` fallback when no level is configured.

| Archive under `test_output/fast-wave-optimization-20261110/` | Complete command (s) | Simulation (s) | Native result |
| --- | ---: | ---: | --- |
| Prior optimized baseline, archived separately | 87.082 | 35.043 | PASS |
| `P1-batched` | 83.654 | 32.099 | PASS |
| `P2-geometry-stencils` (stencils/descriptor precursor only) | 82.827 | 30.524 | PASS |
| `P2b-block-geometry` (full geometry implementation) | 82.224 | 30.021 | PASS |
| `P3-signatures` | 80.589 | 28.040 | PASS |
| `P4-compact-spmv` | 80.263 | 26.573 | PASS |
| `final-fast-wave` (final guarded source) | 79.387 | 26.326 | PASS |

Each PASS includes completed iterations/rank logs, empty proton/electron/
combined diffs and six exact full-file hashes. Per-stage `result.json`,
`run.log`, `lines.jsonl`, `time.txt`, executable and output files retain the
measurements. The stencil precursor is not evidence for full block geometry.
Elapsed phase markers include scheduling and other system load. Small stage
differences should not be projected to other machines or workloads.

Relative to the already optimized baseline, final simulation time falls
24.9% and complete command time falls 8.8%. Against the original unoptimized
119.036 s command, the total reduction is 33.3%, combining the earlier compiler
policy with the current algorithm changes.

## Profile confirmation

The final frozen executable has SHA-256
`1bf560dd94a1b25de81c2ac98cf9aedfb81b08455410e85e592f35a4f39b4015`.
`final-profile/run-summary.json` records a second uninstrumented control
(26.410 s) and instrumented execution (27.626 s, 4.60% overhead). Both PASS
unchanged numerical comparisons and six exact file hashes. The profiler is
compiled separately; this confirmation does not rebuild AMPS.

The sampler collects 10,791 main-thread PCs at 100 Hz, 9,204 outside wrapped
MPI. Self-CPU shares in that category are deposition 17.29%, memory copy
14.69%, corner stencils 8.85%, mover 7.79%, center stencils 4.40%, signatures
3.76% and sparse multiply 3.73%. Changed percentages have changed denominators;
they do not measure individual elapsed-time savings or identify memory-copy
callers. MPI elapsed includes library CPU, polling and descheduling. See the
[`profiling/README.md`](profiling/README.md) for signal safety, interval,
symbols, call coverage and measurement limitations.

The original periodic send/receive call sites take 5.325 s mean elapsed per
rank in the earlier profile; the final periodic batch wait site takes 1.117 s.
The implementations and instrumented runs differ, so this is supporting
call-site evidence rather than a critical-path speedup. Final profile PNG/PDF
outputs are under `final-profile/`; `performance-comparison.png` and `.pdf`
and `performance-summary.json` summarize unchanged-workload stage timings
and independently rechecked full-file hashes in the parent evidence directory.

## Permanent independent tests

Eighteen new case identities are registered, bringing the complete registry
to 56 suites / 1,192 cases (core: 51 suites / 1,085 cases):

- Seven periodic transport cases, each at 1/2/4 MPI ranks: local, remote,
  changed layout/owners, fragmented messages, alias fallback, empty traffic
  and non-world communicator. Independent record bytes and ordered reference
  results check every participating rank.
- Three complete-stencil cases: exact pointer/ID/weight order, atomic fallback
  rejection and refreshed nodes/weights.
- Three sparse-view cases: cancellation-sensitive ordered arithmetic,
  coefficient/storage/topology refresh and empty/invalid inputs.
- Four CRC cases: independent bitwise stream oracle, split updates, copied
  state and multithreaded construction/use.
- Native `PIC-PBUFFER-SIGNATURE-001`: compiled record layout, independent CRC,
  historical canonical CloneParticle reference, link inclusion/exclusion,
  incremental state and source immutability.

The focused performance suite passes 13/13; the general numerics suite passes
14/14, including its ten previous cases. The complete native PIC owner passes
95/95 in the current core campaign, including the new signature case.
Commands, independent oracles and individual IDs are in
[`src/test/tests/pic/performance/README.md`](../../src/test/tests/pic/performance/README.md),
[`src/test/tests/general/numerics/README.md`](../../src/test/tests/general/numerics/README.md)
and [`test/srcCoreUnitTests/README.md`](../srcCoreUnitTests/README.md).

## Complete core campaign and registration repair

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --config test_output/fast-wave-optimization-20261110/machine-overlay.json \
  --match 'CORE-*' --jobs 8 --build-jobs 16 \
  --output-dir test_output/fast-wave-optimization-20261110/core-validation
```

Final reconciled outcome: **1,085 individual results = 843 PASS / 1 FAIL /
240 SKIP / 1 ERROR**; **51 suites = 44 PASS / 1 FAIL / 5 SKIP / 1 ERROR**.
All **825 previous PASS identities** remain PASS, no prior identity is missing,
and all **18 new cases** PASS. Exact remaining nonpasses are unchanged:

| Identity | State | Retained limitation |
| --- | --- | --- |
| `CSWC-FULL-CSWC0628-MASS-TIME` | ERROR | Routine refinement executable exceeds its unchanged 180 s ceiling; partial levels do not qualify BG3D-4/G1 |
| `T12-COVERAGE-STRICT-001` | FAIL | Components remain partial/planned; strict completeness remains open |

The first campaign takes 1,666.301 s and reports 797 PASS / 4 FAIL / 240 SKIP /
8 ERROR across 1,049 rows. Its harness catches three registration checks:
selection, traceability and coverage governance. New compiled case sources
were missing from the traceability ledger, interpolation lacked its ledger
entry, a suite reference belonged to the wrong component, and the selection
test still pinned the previous totals/campaign sets. These are corrected
explicitly; the expected sets still check every previous owner plus the new
one. Standalone traceability passes 29 areas, 218 translation units, 213 public
contracts, 25 configuration contracts, 56 suites and 1,192 cases.

`implementation-source-freeze.json` records the pre-campaign kernel and
registry hashes. Only registration/traceability metadata changed during this
run; do not interpret the old inventory hash as the corrected registry hash.
`source-freeze-audit.json` identifies the changed registration metadata and
records the separate selection-test correction. Physical kernels and the
native signature fixture retain their pre-campaign hashes. Scientific gates,
tolerances and solver controls are unchanged.

During the first campaign, a documentation-only `make PDF` was incorrectly
scheduled beside native testing. The conservative native process guard counts
every make invocation, so seven profile adapters (VSS, rotational Z=2/Z=10,
packed individual split/merge, aligned fixed/individual layout and local
weight/time) refused preflight before any compilation or physics execution.
These are execution errors, not passes or algorithm failures. All seven
subsequently pass real clean native reruns from an idle root. Do not run even
documentation make targets beside a native validation campaign.

The eight-owner retry takes 328.049 s: 55 PASS / 1 FAIL across 56 cases.
Its only failure is the still-stale selection expectation. After pinning the
new explicit totals/campaign sets, the final complete harness retry takes
11.513 s and passes **13/13**. The sum of execution intervals is 2,005.863 s;
the gap between campaigns is excluded. Every earlier attempt is preserved.

Exact retry command (same two `--config` options as above):

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --config test_output/fast-wave-optimization-20261110/machine-overlay.json \
  --test CORE-HARNESS-001 --test CORE-PIC-COLLISION-VSS-001 \
  --test CORE-PIC-IDF-ROTZ-Z2-001 --test CORE-PIC-IDF-ROTZ-Z10-001 \
  --test CORE-PIC-SPLIT-MERGE-001 \
  --test CORE-PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-001 \
  --test CORE-PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-001 \
  --test CORE-PIC-WEIGHT-TIME-LOCAL-001 --jobs 8 --build-jobs 16 \
  --output-dir test_output/fast-wave-optimization-20261110/core-validation-retry

python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --config test_output/fast-wave-optimization-20261110/machine-overlay.json \
  --test CORE-HARNESS-001 --jobs 1 --build-jobs 16 \
  --output-dir test_output/fast-wave-optimization-20261110/core-harness-final
```

Evidence root: `/home/vtenishe/Mars2/AMPS/test_output/fast-wave-optimization-20261110/`.
`core-validation-reconciled.json` and `.xml` publish the final dispositions
with attempt history, exact original/retry report paths, real owner logs,
source identities and a combined stream receipt. They do not overwrite the
raw reports in `core-validation/`, `core-validation-retry/` or
`core-harness-final/`. `core-validation-comparison.json` lists the eighteen
new cases and unchanged nonpasses. `reconcile_validation.py` checks all prior
identities, all prior PASS states, all new PASS states, the exact existing
FAIL/ERROR identities and physical source hashes before writing the report.
