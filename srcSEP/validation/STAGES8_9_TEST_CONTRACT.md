# Stages 8–9 test contract

Baseline commit: `b498050`. The user's Makefile and unrelated untracked assets
are preserved. Evidence root:
`/data/vtenishe/TEST/sep-stages8-9-20261123T003839Z`.
`baseline.json`, `baseline.patch`, `baseline-status.txt` and `baseline-Makefile`
retain exact source identities and starting state. The accepted physical scope
is [the backend support matrix](STANDALONE_BACKEND_SUPPORT.md).

## Invariants and oracles

The existing reduced core remains the physics authority. All 1-D displacement
is relative to the material line; U/grad(U) remain in geometry, source flux,
strain and work. Physical dt stays **300 s**, with **400 samples** per full
application output. Internal birth lifetimes/numerical subdivisions create no
extra samples. No reference, confidence bound or existing **1e-10** independent
relative number/energy budget limit changes.

| IDs / procedure | Production entry and independent acceptance |
| --- | --- |
| SBACK01 | Typed backend admission; reduced prerequisite available, invalid enum rejects, actual profile checks remain mandatory |
| SBACK02 | Actual canonical SWCME coefficients; independent algebraic `div(nU)` and radial-flux identity expose missing conservative evolution, which must reject |
| SBACK03 | Adapter vs direct canonical SI n/U/B/div(U), thermal availability and epoch/state identity; failed preparation preserves prior state |
| SBACK04 | Real neutral codec/importer; labelled software fixture preserves exchange, material request invokes no callback and retains identity, missing line/unknown use rejects |
| SBACK05 | Full regional backend returns unsupported with outstanding BG3D-4 qualification |
| SBACK06–08 | Actual standalone main, MPI=2, all movers: exit 1 with intended unsupported diagnostic before PT allocation/iteration; absent diagnostic/startup is ERROR |
| SAPP01–03 | Actual complete host/source/movers/sampler: original 400-step checker requires physical endpoints/cadence, independent `fsum` signed ledgers, identical source history, output schema and segment-age/exposure reconstruction |
| Existing numerical/native/material/stage-7 gates | Original coefficient/transport, geometry/source/merge, ownership and fresh restart assertions; exact same-decomposition/output comparison, 1↔8 keyed repartition, late failure rollback |
| Five existing core/geometry suites | Original reduced, SWCME, coronal-CME, full regional and PIC cases; unsupported/dormant/unqualified failures remain visible |
| Four existing coupling profiles | Europa, OH, PIC coupling/restart/divE/dynamic allocation and dynamic flux rope: original physical decks and DiffNum references/limits, actual MPI receipts and restart dependency order |

Strict source builds use warnings as errors. ASan/UBSan explicitly disable
LeakSanitizer where unavailable; no leak qualification is claimed. Native
tests use actual MPI and matching compiler-dependency/object/archive/executable
closure. Reports retain SI/HCI conventions, seed 0 or existing fixed seeds,
configuration/assets/hashes, metrics, commands and individual verdicts.

## Maintained execution and negative controls

The shared registry adds `CORE-SEP-BACKEND-ADMISSION-001`,
`CORE-SEP-BACKEND-NATIVE-001` and `CORE-SEP-FULL-APPLICATION-001` through the
existing standalone wrapper. The application orchestrator adds matching
`backend-admission`, `backend-admission-native` and
`standalone-full-application` steps with strict individual result readers.
Its `standalone-application-runner` step additionally requires all eleven
software controls, rather than accepting a successful subset of the unittest
stream. Both supported Python verbose method formats preserve the same IDs.
The application step is long: `--quick` is not release qualification. Make
targets and source-suite selection expose the same procedures. Catalog,
traceability and independently pinned harness totals change together.

The application launcher audits the matching configured native closure and
ordinary main, freezes a private executable, and writes exact commands/input
hashes. PREPARED is not PASS. Collection binds each actual receipt to planned
argv/cwd before invoking the existing independent verifier. Missing build,
output or MPI evidence yields every required case as ERROR and a nonzero exit;
JSON/JUnit remain available.

Negative controls cover missing movers, shortened windows, stale source or
executable identity, foreign receipts, changed dt/budget limits and failed
native exits. Existing physical gates retain wrong composition, stale epochs,
extra bulk displacement, unsupported roots and segment-lifetime denominators.
Verifier controls damage isolated software-output copies only; physical inputs
and reference observations are never edited. No successful subset hides a
FAIL, ERROR or missing required role.

## Completion rule

Stage 8 closes for the user-approved admission/support matrix, without claiming
transport on rejected models. Stage 9 closes for the supported reduced profile
after its required numerical/native/restart/application/coupling results are
inspected. Broad core-suite unqualified cases remain explicit and block their
dependent backend. Report actual counts, metrics, source/build identities,
commands and scope; preparation/source preflight cannot replace native gates.
