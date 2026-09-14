# WP31–WP41 validation report

## Source-only result

`make test-wp31-wp41-unit` passes under C++11 strict warnings, AddressSanitizer,
and UndefinedBehaviorSanitizer. It reports eleven named PASS results plus the
WP35 generated-table and production source-integration checks.

| Work package | Executed evidence | Result |
| --- | --- | --- |
| WP31 | Shared planner, stiff limit, automatic temporal-order fit | PASS, observed order 0.969891 |
| WP32 | Typed positivity correction, signed rejection, nonfinite failure | PASS |
| WP33 | Split/merge moments and deterministic lineage | PASS |
| WP34 | Native-observation contract and false-promotion rejection | PASS at source-contract level |
| WP35 | 90-row Cartesian product, diagnostics, preflight, generated table | PASS |
| WP36 | Global closure, identity/restart round trip, checksum fault | PASS |
| WP37 | Versioned domain-separated seed panel and mean gate | PASS |
| WP38 | IEEE boundary cases, counterexample identity, fault hit | PASS |
| WP39 | Analytical two-bin instrument response | PASS at analytical-core level |
| WP40 | Exact work and environment-specific wall-time policy | PASS at contract level |
| WP41 | Evidence-level validation and single flush owner | PASS |

## External gates

| Gate | Status in this handoff | Required completion input |
| --- | --- | --- |
| Linked native AMPS mover/adapter execution | BLOCKED | Enclosing AMPS checkout, generated PIC types, executable |
| OpenMP scheduler and MPI decomposition equality | BLOCKED | Native executable and launch environment |
| Strong and weak scaling | BLOCKED | Stable hardware/compiler baseline and native workloads |
| Real SWMF replay | BLOCKED | Checksum-complete SWMF manifest and archived inputs |
| Held-out observations | BLOCKED | Instrument response, spacecraft products, event manifests |

These blockers do not downgrade the source-only PASS results, but those results
must not be cited as native, coupled, scaling, or observational evidence.
