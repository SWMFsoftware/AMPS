# AMPS core-test status

Last updated: 2026-11-10

## Fast-wave P0–P4 implementation (2026-11-10)

Implemented ordered per-peer periodic byte transport with a streaming
fallback; complete uniform pointer stencils and per-call deposition geometry;
exact full-width legacy CRC/canonical particle-copy improvements; and a
solve-scoped compact scalar sparse operator. Source comments and adjacent
READMEs document equations, arithmetic order, ownership/lifetimes, failure
conditions, backend guards and fallback behavior. Physics settings, solver
tolerances, deterministic checks and regression references are unchanged.

The final exact `make test_fast-wave -j` passes: 79.387 s complete / 26.326 s
simulation versus the prior optimized 87.082 / 35.043 s. All measured stages
finish 60 iterations on four ranks, leave three empty diff files and reproduce
six full sampled species files byte-for-byte. A second final-binary control
passes in 26.410 s. See `test/srcFastWave/STATUS.md` for full evidence and
portable interpretation limits.

Eighteen new permanent cases bring the registry to 56 suites / 1,192 cases;
the core selection is 51 / 1,085. The focused performance owner passes 13/13,
general numerics 14/14 and native PIC 95/95. Complete core validation with
documented retries finishes at **1,085 cases: 843 PASS / 1 FAIL / 240 SKIP /
1 ERROR**; **51 suites: 44 PASS / 1 FAIL / 5 SKIP / 1 ERROR**. All **825
previous PASS identities** remain PASS, no previous identity is missing and
all eighteen new cases PASS. The final harness passes all 13 checks.

The initial campaign takes 1,666.301 s. Registration/traceability omissions
and stale explicit selection counts/sets cause three harness checks to fail.
A competing documentation-only make also causes seven native preflight
refusals, before physics execution. Metadata is repaired while physical
kernel hashes remain frozen. The eight-owner retry takes 328.049 s; all seven
native owners pass, with only the selection check still failing. Its corrected
explicit expectations then pass the complete 13-case harness retry in
11.513 s. No acceptance tolerance is relaxed. The source-freeze audit records
metadata hashes separately; sum of execution intervals is 2,005.863 s.

Evidence root: `test_output/fast-wave-optimization-20261110/`.
`core-validation-reconciled.json` and `.xml` retain final dispositions, actual
owner logs, the untouched first/retry reports and every attempt in provenance.
`core-validation-comparison.json` records exact case-identity preservation,
all eighteen new PASS identities and unchanged nonpasses. The raw report
directories are `core-validation/`, `core-validation-retry/` and
`core-harness-final/`; the reconciled report is explicitly a multi-attempt
evidence assembly, not a claim of one all-green run.

The only actual remaining nonpasses are unchanged:
`CSWC-FULL-CSWC0628-MASS-TIME` ERROR at its existing 180 s routine ceiling,
and `T12-COVERAGE-STRICT-001` FAIL for incomplete coverage. BG3D-4 and G1
remain unqualified. The earlier core campaign below is historical.

## Fast-wave individual-function profiling (2026-11-10)

Completed isolated control and CPU/PMPI-profile runs of the previously
validated frozen -O1 binary, with no AMPS rebuild or production algorithm
change. Both PASS: 60 iterations, four rank logs, both unchanged numeric
comparisons and all six sampled species files byte-identical to baseline.
Elapsed times: control 35.191 s, profile 36.419 s (3.49% overhead).
The sampler collected 14,313 main-thread CPU PCs at 100 Hz with no drops.

Top self CPU outside wrapped MPI: deposition 14.58%, memory copy 10.47%,
matrix-vector multiply 10.27%, particle signatures 9.10%, corner interpolation
7.71%, cell-center interpolation 6.81% and the Lapenta mover body 6.75%.
The periodic ghost-block MPI send/receive sites cost 5.325 s mean elapsed per
rank; timer-report gathers cost only 0.00318 s per rank. Timer percentages,
MPI elapsed and possible speedups are not interchangeable.

Added optional profiling tools and detailed methodology under
test/srcFastWave/profiling/. Proposed implementation stages, independent
validation requirements, source pointers and limitations are in
test/srcFastWave/PERFORMANCE_OPTIMIZATION_PLAN.md. The proposal preserves
nightly sorting/checks, physical work and tolerances; it does not implement
or claim any further speedup. Evidence and plots are under
test_output/fast-wave-function-profile-20261110/. The full-core campaign
result below remains historical.

## Legacy fast-wave performance validation (2026-11-10)

Timed and optimized the separate `make test_fast-wave -j` procedure. Its
GCC-only `-O1` fallback applies when the effective C/C++ flags lack an explicit
optimization level; it preserves configured levels, non-GNU/NVCC policies,
site files and all physical/numerical gates. Disable with
`FAST_WAVE_CXXFLAGS=`. The current host-generated recipe and maintained Table
both forward the fallback only for the native fast-wave build.

Four clean four-rank native runs PASS: baseline 119.036 s, explicit -O1
86.989 s, explicit -O2 91.794 s and the exact command with the new default
87.082 s. The final command is 26.8% faster; its simulation phase is 2.21x
faster. Every run completes 60 iterations and all three diff files are empty;
all six sampled species files are byte-identical across the four runs.
Twelve compiler-policy make probes PASS. Native builds retain the prescribed
preflight and -j16 top-level workflow. Reference hashes are unchanged and the
original output symlink was restored. An initial sandbox MPI socket failure
is retained separately and excluded from the successful native timings.

Detailed commands, methodology and limits: `test/srcFastWave/README.md`.
Evidence: `test_output/fast-wave-performance-20261110/summary.json`, policy
probe JSON and per-run logs/timing/results. The complete core result below
remains historical; this change does not claim a new all-core campaign.

## Complete core revalidation finished (2026-11-10)

The source-frozen full core campaign finished in **1958.293 s**:
**1,067 cases = 825 PASS / 1 FAIL / 240 SKIP / 1 ERROR**;
**50 suites = 43 PASS / 1 FAIL / 5 SKIP / 1 ERROR**. The return code is 1,
retaining the two actual nonpasses. All 792 previously passing identities are
still PASS, all 32 added cases PASS, no baseline case is missing, and the only
changed baseline result is CMBGU05 FAIL -> PASS. The original 785-pass
preintegration baseline remains protected by the retained earlier identity
comparison and this complete 792-pass preservation.

Exact executed command:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --config test_output/remaining-failures-repair-20261110/local-machine.json \
  --match 'CORE-*' --jobs 8 --build-jobs 16 \
  --output-dir test_output/remaining-failures-repair-20261110/core-current
```

The explicit selector excludes application suites. The AMPS-local overlay
keeps output/work below the permitted checkout. Native adapters recorded the
required root/idle, generated-build removal, configuration workflow and -j16
build preflight. All core native, 1/2/4-rank ownership, 1D/2D/3D mesh, physics,
layout and weight/time cases pass. Existing typed provider/tool/campaign SKIPs
remain. This is complete core regression evidence, not a native full-CME run.

Remaining exact nonpasses:

| Case | State | Reason / preserved gate |
| --- | --- | --- |
| `CSWC-FULL-CSWC0628-MASS-TIME` | ERROR | Unchanged full production refinement exceeds its original 180 s ceiling; no partial-level PASS is inferred. |
| `T12-COVERAGE-STRICT-001` | FAIL | Valid real inventory still has 2 complete / 26 incomplete components; structural contradictions remain separate ERROR conditions. |

The six T11 import/audit suites now report 668 rows: 434 PASS / 233 SKIP /
1 ERROR. T12 reports 13 PASS / 2 SKIP / 1 FAIL: Clang and TSan are skipped,
strict coverage fails, and the other bounded tool/ABI/sanitizer/backend/perf
cases pass. The full-volume owner has 46 rows: 45 PASS / 1 ERROR.

Full evidence directory:
`/home/vtenishe/Mars2/AMPS/test_output/remaining-failures-repair-20261110/core-current/`.
The combined result is `amps-core-tests.json`/`.xml`; stream progress/report is
`streams/core.log`/`.json`; suite logs and local numeric transcripts are linked
by each JSON suite row. `core-current-baseline-comparison.json` in the parent
evidence directory retains identity comparisons and every remaining nonpass.
`full-core-preservation.json` verifies no source/test change during the run,
unchanged HEAD, and unchanged Makefile.conf/.amps.conf. Makefile.local was
regenerated through the normal configuration workflow; no direct .amps.conf
edit or Git mutation occurred. Current Git-read-only instructions remain in
force; no Git-writing application test was selected.

## Documentation and subsequent diagnostic boundary

All new APIs/numerical repairs, units/frames, invariants, input/API boundaries,
exact commands, independent oracles, portable settings, failure/report policy
and unqualified physics are documented in the relevant READMEs:

- `src/models/sep_corona_swcme/README.md` and `README_CROSS_RAY.md`;
- `src/models/sep_corona_swcme/test/README.md`;
- `src/test/README.md`, `tests/models/cme/README.md` and
  `tests/models/cme/sep_corona_swcme_full/README.md`;
- `src/test/tests/harness/self/README.md`,
  `src/test/tests/portability/matrix/README.md`, `src/test/lib/README.md`, and
  `src/test/tests/models/comet_libraries/README.md`.

Documentation validation: all 60 local links across these ten READMEs
resolve; the generated catalog check passes at 55 suites / 1,174 cases;
traceability passes at 28 areas, 218 production translation units, 210 public
contracts and 25 configuration contracts. The material-transport C++ example
compiles with warnings treated as errors, links against the actual tested
shared archives and prints the expected `8 0.25`. The harness governance
recheck passes all 13 cases. `git diff --check` passes. Detailed receipts are
`readme-link-audit.json`, `readme-example-validation.json`,
`readme-governance-check/amps-core-tests.json` and
`readme-documentation-validation.json` in the evidence root. These checks
validate documentation, registration and governance; the full physics result
above remains the source-frozen core receipt.

The source-level physical implementation has not changed since this full run.
The new comparison PNG uses recorded contact/advection/reconstruction errors:
[bounded comparison plot](../../test_output/remaining-failures-repair-20261110/bounded-comparisons.png),
with data/provenance JSON beside it.

Private ambient-flow/variational diagnostics identify a numerical trace issue
before adopting a production analytical lateral map. Closed-corotation/Parker
probes conserve density/field to roundoff/about 1e-12. The open-PFSS 60 s probe
has a timestep-independent density mismatch about 4.88e-8; an independent
analytic dipolar invariant finds about 1.8e-6 pointwise relative density error
from the existing endpoint trace. A PRIVATE exact-dipolar-trace prototype
reduces the material-density mismatch to about 1.6e-10 and endpoint-density
error to about 1e-13. Its complete unchanged-gate CMBGU05 run passes in 80.137 s.
No such prototype or ambient-map production API is imported into source. The
adiabatic probes also expose the expected need for explicit heating when
matching the isothermal ambient; they do not qualify thermal/source balance.
See README_CROSS_RAY.md and the analytical design document for provenance and
next steps. These experiments do not replace the current production receipt.

The following implementation checkpoints retain their historical snapshots.

## 2026-11-10 analytical cross-ray implementation and validation

The user explicitly resumed G1, retained the full 3-D ambient, and selected
conservative lateral MATERIAL transport without an MHD runtime. The original
runtime prohibition remains in force. HLL/static-fluid operators are reference
conservation tools; they are not a production MHD time integrator. Further
user decisions select a field-threaded ideal material contact, exact startup
corotation, and a 1200-s quintic transition to zero angular rate. The final
orientation retains the agreed 0.0985-degree offset without a reset.

Implemented bounded components under `src/models/sep_corona_swcme`:

| Component | Independent new cases | Current evidence |
| --- | ---: | --- |
| SI MHD/ALE flux and conservative transfer references | 5 | Optimized, debug and ASan/UBSan/LeakSanitizer PASS |
| Exact static spherical sectors | 3 | Analytic octant, rotations, thin cells and 336-cell shell PASS |
| Shared face/edge topology | 3 | 792 internal interfaces, integer incidence, curved edges and authority PASS |
| Static face-flux CT reference | 3 | Independent Stokes/Faraday and deliberate monopoles/stale data PASS |
| Cell-field reconstruction | 3 | Constant rotated/scaled fields; affine-field errors 0.0152182/0.00391239/0.000984985 PASS |
| Atomic fluid/CT energy reference | 4 | Uniform states, opposite transfers, explicit work and whole-epoch rejection PASS |
| Corotating threaded-contact geometry | 4 | Full startup, independently differentiated motion, interface stress/E/work and exclusions PASS |
| Analytical material transport | 4 | Full Cauchy mass/flux/entropy, independent rates, nonlinear transport and map rejection PASS |

These are **29 new cross-ray cases**, plus independent contact-derivative and
stable-shell controls. No production lateral velocity/map has been selected
or guessed. No new runtime selector/application coupling is active. The
remaining regional closure must supply one compatible map, initial/boundary
states, EOS/entropy and force/heating/work, then satisfy the original computed-
shock, mass/time/ray, contact, magnetic and validity gates. BG3D-4/G1 remain
UNQUALIFIED. See the normative cross-ray, startup and material-transport docs.

Key retained outcomes:

- Complete `CMBGU05` PASS, retaining the original three derivative steps,
  physical leakage bounds and convergence requirement. The user approved the
  stronger distributed acceptance requiring EVERY level to pass unchanged
  limits. Physical force/work and cross-stage authority mismatch remain
  diagnostic; this software PASS does not qualify the old volume map.
- Original complete coarse refinement: 471.164 s and TWO detected shock zones.
  Algebraically stable exact shell moments: comparable coarse diagnostic
  127.902 s, same printed shock values and two zones. The unchanged full matrix
  retains ERROR at its 180 s routine limit; no smaller grid/ray/epoch set,
  extra timeout, changed branch filter or relaxed physical gate is accepted.
- Stable production shell geometry uses independent Gaussian polynomial
  moments as its test oracle. A cancellation-sensitive old centroid oracle
  initially exposed CSWC0625 failure and was corrected independently. Radial
  flux error is now 4.01794e-16 against the unchanged 2e-13 bound; CSWC0624/26
  also PASS. No PFSS optimization prototype was imported.
- The stationary legacy HCI contact has initial ambient normal speeds from
  -905 to +410 m/s on the four supported rays. The selected rotating contact
  reduces maximum startup mass flux to 1.10794e-26 kg/m^2/s, below the unchanged
  1e-20 bound, retaining abs(Bn)/abs(B) up to 0.323732. Actual two-sided ejecta
  matching remains a later G1 gate.
- Fresh complete 42-case model-owner plus 13-case harness snapshot:
  **55 total, 54 PASS, 0 FAIL, 0 SKIP, 1 ERROR; 541.326 s**. The sole ERROR is
  the original refinement ceiling. It predates the four subsequent material
  cases and later private-factory/excluded-frame controls; those have separate
  current build/debug/sanitizer receipts. Do not silently relabel this snapshot.
- A fresh current archive/neutral-header build passes all 25 earlier cross-ray
  cases, derivative/shell controls, CSWC0625 and architecture after validated-
  factory-only construction and HCI exclusion bookkeeping. The four new
  material cases independently pass optimized/debug and all sanitizers.
- Latest harness: **13/13 PASS**. In addition to the real validator and strict
  CLI, an injected ledger contradiction must produce ERROR=2 with its unique
  diagnostic; the legitimate incomplete inventory remains FAIL=1. The first
  new bootstrap fixture formatting error is retained, then corrected. No
  production coverage component was promoted by these governance controls.

Evidence root: `test_output/remaining-failures-repair-20261110/`. Important
receipts: `model-owner-with-corotation/amps-core-tests.json`,
`final-kernel-validation.json`, `cross-ray-sanitizer-factory-results.json`,
`contact-latest-profiles.json`, `contact-latest-sanitizers.log`,
`material-controls-current.json`, `material-profiles.json`,
`material-sanitizers.log`, and
`harness-final-controls-corrected/amps-core-tests.json`.

Current registry: **55 suites / 1,174 cases**, including **50 core suites /
1,067 core cases**. At that checkpoint a complete core campaign was running in
`core-current/`, selected with `--match 'CORE-*' --jobs 8 --build-jobs 16` and
explicit Passion plus AMPS-local path overlays. Registered applications have
zero selected suites. Native adapters perform and record root/idle checks,
removal of only `build`, Config.pl regeneration and top-level `make -j16`.
Source/test definitions are frozen during the campaign; documentation updates
continue. Its final counts and identity comparison are recorded in the completed section above.
The previous complete 1035-case receipt (792 PASS, 2 FAIL, 240 SKIP, 1 ERROR)
remains historical evidence, not a claimed result for the enlarged registry.

Strict coverage remains **2 complete / 26 incomplete** and exits 1. The
separate Hartley2 3-D geometry/runtime limitation remains open. No Git state
is modified; no test-output tree is imported as application source. No native
one-/four-rank full-CME or completed-coverage qualification is claimed.

The following checkpoints retain their historical counts and evidence.

## Remaining-failure repair implementation in progress

The full `CMBGU05` executable now passes with the unmodified coronal library
(223.624 s in the retained run, including concurrent diagnostic work). The
fourth-order position derivative uses equivalent paired differences, retaining
all frozen steps and physical/convergence bounds. The user approved changing
the old expectation that some distributed contact level fail into the stronger
requirement that every level satisfy the unchanged limits. All subsequent
inventory, material/differential, transactional, spatial-inverse and residual
checks execute and pass. This remains bounded software validation: the old
map still reports its physical force/work and cross-stage authority mismatch.

The new independent derivative/slip controls pass optimized/debug execution and
ASan/UBSan/LeakSanitizer. Restricted-sandbox LeakSanitizer failure is retained
separately; the same binary passes with every sanitizer enabled outside that
ptrace environment. A prototype-linked contact run is retained as superseded
diagnostic evidence; the cited final contact receipt verifies the original
coronal archive was rebuilt. The harness passes 13/13 after correcting its
private inventory fixture's required schema version. Strict coverage now
cross-checks the real source/configuration ledgers, rejects contradictory
coverage labels, and preserves setup ERROR versus valid incomplete FAIL.

Refinement remains open. The unchanged executable reproduces the 180-second
timeout. A private degree-one field prototype also times out and has not been
added to production. A separate algebraically stable shell-geometry diagnostic
completes the 40-cell/3600-s level in 127.902 s but detects two shock zones;
this single-level measurement is not an acceptance PASS. The unchanged
original executable's complete coarse measurement is being checked next.
At that earlier checkpoint there was no production piston-solver/geometry or PFSS change.
The later stable shell changes and validation are described above.

Evidence root: `test_output/remaining-failures-repair-20261110/`.
Registry: 55 suites / 1,144 cases, including 50 core suites / 1,037 cases.
No new complete core campaign or strict-coverage PASS is claimed.
The user explicitly resumed full-CME implementation and qualification; work
returns to G1/BG3D-4 in its authoritative ledger. G2 and later particle stages
remain outside the resumed scope. All original numerical/physical gates remain.

## Remaining-failure repair plan prepared

The user requested a plan that preserves the validation gates. The detailed
[remaining-failure repair plan](REMAINING_FAILURES_FIX_PLAN.md) covers
`CSWC-FULL-CMBGU05`, `CSWC-FULL-CSWC0628-MASS-TIME`, and
`T12-COVERAGE-STRICT-001` separately. Inspection identifies the immediate
contact failure as its finite-difference monotonicity clause: all three
printed pointwise fluxes are below the frozen physical bound. Roundoff remains
a hypothesis requiring independent measurement; the plan does not remove the
convergence condition. The refinement must complete every unchanged level
within the existing 180-second routine limit. Strict coverage must earn actual
closure of 26 incomplete components, without exemptions or label-only changes.
This is planning only; the current 792 PASS / 2 FAIL / 240 SKIP / 1 ERROR
receipt remains authoritative. The separate Hartley2 3-D geometry/runtime issue
and deferred full-CME qualification remain open.

## Native configuration-isolation repair validated

The planned repair is implemented in the exosphere, dust charging/injection,
and ECSIM stencil adapters. They explicitly select no-SPICE mode and verify
saved/effective make settings and real compiler commands, retaining configuration
receipts and rejecting contradictions before MPI. The new
`HARNESS-NATIVE-CONFIG-ISOLATION-001` passes 24 synthetic adapter scenarios:
eight valid enabled/disabled starts and sixteen expected rejection scenarios.
The complete harness passes 12/12. Real exosphere and dust charging runs both
PASS from independently seeded SPICE-enabled and disabled starts. All four
native matrix rows and the reverse-order pair PASS, as does reselecting a
SPICE-enabled application through Config.pl. Full core revalidation from that
enabled initial state completed in 1,744.066 s: 1,035 cases = 792 PASS / 2 FAIL /
240 SKIP / 1 ERROR; 50 suites = 43 PASS / 1 FAIL / 5 SKIP / 1 ERROR. All original
785 baseline PASS identities are preserved, every other original status is
unchanged, and all seven new cases PASS. Compared with the user's failing run,
exactly the exosphere and dust charging cases change FAIL -> PASS; all 789
previously passing cases remain PASS. The nonzero runner exit is caused by
the original two failures and one error, not a new configuration failure.
Evidence root:
`test_output/native-config-isolation-20261110/`.

Exact final full-core command (exit 1 retains the original nonpasses):

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --match 'CORE-*' --streams 1 --jobs 8 --output-dir test_output/native-config-isolation-20261110/core-final
```

| Evidence | Result |
| --- | --- |
| `harness-final/amps-core-tests.json` | PASS 12/12; new control case includes 24 synthetic scenarios |
| `native-matrix.json` | PASS: four independently seeded native rows plus the reverse-order pair (six actual physics-case executions) |
| `spice-boundary.json` | PASS: normal application SPICE selection remains available |
| `core-final/amps-core-tests.json` and `.xml` | Complete 50-suite / 1,035-case report, including exact retained nonpasses |
| `baseline-comparison.json` | PASS 785/785 baseline passes; no original status changes |
| `reported-run-comparison.json` | PASS; exactly the two requested FAIL -> PASS repairs |
| `invariant-audit.json` | PASS: site file, Config.pl, ampsConfig.pl and all four native input decks unchanged; no staged paths; HEAD unchanged |

Current registry: 55 suites / 1,142 cases globally; 50 suites / 1,035 cases in
the core stream. The nightly selection has 43 suites / 382 cases. Historical
1,034-case receipts below are retained without relabeling their counts.

## Historical configuration failures and separate legacy 3-D runtime timeout

The pre-repair user-reported run `run-20261110T092615Z-924820` under
`/data/vtenishe/AMPS-test-results/`. It reports 789 PASS / 4 FAIL / 240 SKIP /
1 ERROR, 50 suites, 1,723.884 seconds. Compared with the frozen baseline and
the earlier passing integration campaign below, precisely two cases changed
PASS -> FAIL:

- `CORE-EXOSPHERE-SPHERE-INTEGRATION-001`: SPICE leapsecond-kernel error at runtime.
- `CORE-DUST-CHARGING-DYNAMICS-001`: undefined `Exosphere::IAU_FRAME` at link time.

Both adapters inherited global SPICE build selection instead of explicitly
selecting their intended no-SPICE mode. The later Hartley2 build exposed this
configuration-isolation gap. Before repair, 783/785 original passes remained
PASS, with all six then-newly registered cases PASS. The accepted repair and
current complete comparison are recorded above.
The original two failed cases and timeout error are unchanged. Exact evidence:
`test_output/comet-integration-20261109/provenance/core-reported-20261110T092615Z-comparison.json`.

The [configuration-isolation repair plan](CONFIGURATION_ISOLATION_FIX_PLAN.md)
is implemented and validated. Its original steps remain as the acceptance
checklist; the current counts include the new harness case.

The approved `1d_comet_dev` integration adds `INP-05` for the generic case
selector and `CORE-COMET-LIBRARIES-001` for reaction-network, collision-data,
radiative-cooling, fluorescence, and reference-data checks. The expanded
comet-era registry was 55 suites / 1,141 cases; core alone was 50 suites / 1,034 cases.
That nightly selection had 43 suites / 381 cases. Catalog generation and
structural traceability pass. The earlier Passion integration campaign completed in
1,719.797 seconds: 1,034 cases = 791 PASS / 2 FAIL / 240 SKIP / 1 ERROR;
50 suites = 43 PASS / 1 FAIL / 5 SKIP / 1 ERROR. The exact baseline comparison
preserves all 785 original PASS identities and every other baseline status;
all six new cases PASS. The two failures and one error remain
`CSWC-FULL-CMBGU05`, `T12-COVERAGE-STRICT-001`, and
`CSWC-FULL-CSWC0628-MASS-TIME` (180-second routine ceiling).

Evidence is under `test_output/comet-integration-20261109/core-qualified/`;
the machine-readable identity comparison is
`test_output/comet-integration-20261109/provenance/final-core-comparison.json`.
Both final-tree legacy make regressions pass with zero-byte numerical diffs.
Hartley2's all-stage runner passes 46 jobs across 14 groups, including fresh
Stage-8/9 matrices at 1/4/12 MPI ranks; earlier fresh Stage-4–7 receipts are
retained separately. The user approved the remaining two-file 3-D build
exception on 2026-11-10. Clean DIM3/global-weight builds pass with the site
default flags (63.932 s) and a make-only `OPT3=-O2 -g` override (97.007 s).
The optimized 100-iteration/four-rank nightly run reached its local 1,800 s
ceiling during initial cut-cell volume construction, before particle stepping.
All owned MPI processes exited. The approved build fixes are verified, but
the 3-D runtime is unvalidated. The user-selected performance investigation is
complete; see [its diagnosis and proposed spatial-index validation plan](GEOMETRY_PERFORMANCE_INVESTIGATION.md).
The geometry algorithm was unchanged during the investigation. See
[the integration ledger](../../COMET_CORE_INTEGRATION_STATUS.md) for commands,
all evidence paths and approved changes. These software
checks do not close the existing coverage or physical-qualification gates.

## Core/application execution streams

The global runner now derives three independent streams from the complete
registered manifest:

| Stream | Suites | Cases | Ownership |
| --- | ---: | ---: | --- |
| `core` | 50 | 1,035 | Registered `src/...` core/model owners; coverage remains partial |
| `application-srcEnceladus` | 3 | 61 | Forward, native startup, and E3 Enceladus owners |
| `application-srcMoon` | 2 | 46 | Moon regression and M0 boundary owners |

`--streams N` sets the maximum number of concurrently active streams and
defaults to one. `--jobs N` remains a distinct per-stream limit for
consecutive nonexclusive suites. The make wrapper exposes the same settings as
`STREAMS=N JOBS=N`. A writer-priority checkout gate spans all streams:
nonexclusive readers may overlap, while every `exclusive: true` suite waits
for all readers and runs without any other stream inspecting or changing the
generated checkout. Native configuration and build operations therefore
remain globally serialized.

A bare `--list` now inventories the complete registry, summarizes all
registered applications, and prints separate core, Enceladus, and Moon stream
sections. Explicit campaign/component/test/case selectors still narrow the
selected stream contents. Execution writes `streams/<id>.log` and
`streams/<id>.json` for every selected stream; `selection.json` records the
plan and `amps-core-tests.json` references all stream receipts while retaining
the single combined result inventory.

Validation evidence:

| Command/evidence | Outcome |
| --- | --- |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --list --streams 3` | PASS: two registered applications plus three separate stream sections; 49/1,028 core, 3/61 Enceladus, and 2/46 Moon |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --all --dry-run --streams 3 --jobs 2 --output-dir /tmp/amps-stream-runner-dry-all` | PASS control-plane run: 54 explicit dry-run SKIPs, six per-stream log/JSON files, and one reconciled combined JSON/JUnit result |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-GENERAL-INPUT-001 --test APPLICATION-ENCELADUS-FORWARD-001 --test APPLICATION-MOON-M0-001 --streams 3 --jobs 1 --output-dir /tmp/amps-stream-runner-real-smoke` | PASS control plane and executed evidence: core 4/4 PASS, Enceladus 42/42 PASS, Moon ten explicit SKIPs; all three streams visibly overlap without a native rebuild |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-HARNESS-001 --streams 2 --output-dir /tmp/amps-stream-runner-harness-final` | PASS 11/11, including exact stream partition/list formatting, bounded stream scheduling, global exclusive-gate invariants, and per-stream receipts |

No native AMPS rebuild was performed while qualifying this orchestration
layer; the checkout gate is independently exercised with deterministic
concurrency fixtures, while the retained real cross-stream run uses only
dependency-light owners and the typed Moon M0 non-execution boundary.

## Required qualification-machine profile

The global CLI now refuses every invocation that omits `--config`, including
`--list` and `--coverage`, with exit 2. Its diagnostic names the resolved
`test/qualification-machine/` directory, the short-hostname profile it
expects, and `template.json`. This makes the choice of compiler, MPI launcher,
filesystem roots, rank count, timeout, and optional qualification
capabilities an explicit part of every test receipt rather than an ambient
default. Machine profiles are deliberately stored in the repository's
top-level `test/` tree rather than `src/test/` because they are operational
test configuration, not AMPS source code. The make wrapper supplies
`../../test/qualification-machine/$(hostname -s).json` explicitly and supports a
`CONFIG=...` override.

The new maintained files are:

- `test/qualification-machine/passion.json`: Passion's installed GNU, OpenMPI,
  MPICH, Fortran, Valgrind, and CUDA commands; `/data/vtenishe` durable roots;
  `/tmp/vtenishe` work root; local CSPICE/kernel paths; 16 build jobs; eight
  MPI ranks; and the Moon-SPICE/Enceladus-E3 capabilities;
- `test/qualification-machine/template.json`: schema-valid placeholders and all
  optional portability-tool keys for a new host; and
- `test/qualification-machine/README.md`: naming, copy/edit, validation, overlay,
  and execution instructions.

T12's recursive strict-coverage identity now writes a schema-limited profile
from its resolved parent context and passes that file explicitly. It therefore
obeys the public requirement without selecting a hostname file behind the
operator's back.

Validation evidence:

| Command/evidence | Outcome |
| --- | --- |
| `python3 src/test/run_tests.py --list` | Expected exit 2 before selection, with the Passion-profile and template paths in the error |
| `python3 src/test/run_tests.py --help` | PASS without a profile; help labels `--config` required and names the profile directory |
| JSON parsing and `py_compile` for the two profiles and three modified runners | PASS |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --full-validation --list` | PASS: all 54 suites and 1,135 cases selected |
| `make -C src/test list` | PASS: the wrapper visibly invokes the global runner with `../../test/qualification-machine/passion.json` |
| `python3 src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-HARNESS-001 --output-dir /tmp/amps-required-config-harness-relocated` | PASS 11/11 after relocation, including real schema loading of Passion/template and omission rejection even when `AMPS_TEST_CONFIG` is set |
| exact `T12-COVERAGE-STRICT-001` at `/tmp/amps-required-config-t12-nested` | Expected FAIL because coverage remains 2 covered/24 incomplete; its retained nested transcript contains the inventory, not a missing-config error |

Historical commands elsewhere in this ledger are retained exactly as they
were executed. Commands recorded before this change can omit `--config`; all
new runs and reproductions must supply a reviewed profile.

## Moon and Enceladus application integration

The global registry now contains 54 suites and 1,135 individual cases across
26 inventory/traceability areas. Five `application-smoke`/`extended` owners
add 107 application-level identities without changing the 42-suite, 375-case
deterministic core nightly selection. The routine `application-smoke` subset
contains four owners and 94 cases; the costly 13-case E3 campaign is retained
in `extended` behind an explicit machine capability.

Strict coverage remains intentionally open: 2 of 26 inventory components are
covered and 24 are partial/planned. Adding globally visible application
campaigns does not convert either application root to covered.

| Global owner | Cases | Current execution boundary |
| --- | ---: | --- |
| `APPLICATION-ENCELADUS-FORWARD-001` | 42 | Executes all 11 dependency-light local software checks and 31 independent numeric reference comparisons using the configured compiler argv |
| `APPLICATION-ENCELADUS-NATIVE-STARTUP-001` | 6 | Exclusive clean Enceladus build followed by exact import of CLI/default-deck, invalid-configuration, and generated-hook startup gates |
| `APPLICATION-ENCELADUS-E3-001` | 13 | Owns every 10,000-step E3 process/ledger/flyby assertion; defaults to individual SKIP unless `enceladus-e3` or global `--full-validation` explicitly activates it |
| `APPLICATION-MOON-REGRESSION-001` | 36 | Exclusive clean Moon/SPICE build, then exact import of U01-U26 and I11-I20 outcomes in routine regression or full-validation qualification mode |
| `APPLICATION-MOON-M0-001` | 10 | Audits every I01-I10 runner/README and emits explicit SKIP until the local multi-build orchestrator satisfies the current native preflight |

The Moon and native Enceladus adapters follow the mandatory native sequence:
they confirm the repository and absence of another native build/run, lock
`Config.pl`, remove only root `build` with `rm -rf -- build`, configure the
selected application, and execute top-level `make amps` with
`MAKEFLAGS=-j16`. The Moon owner then passes `--skip-build` to its local
aggregate, so there is one build authority. SPICE locations and the
`moon-spice` capability are portable machine-overlay settings or the
`AMPS_TEST_SPICE_PATH`/`AMPS_TEST_SPICE_KERNELS` environment overrides; no
application adapter embeds a site path. The E3 local launcher now accepts the
configured launcher argv and rank flag, including wrapper options encoded as
a JSON string array.

Validation and retained evidence:

| Command/evidence | Outcome |
| --- | --- |
| `python3 src/test/run_tests.py --campaign application-smoke --jobs 8 --output-dir /tmp/amps-application-smoke-integrated-final` | Complete integrated routine: 48 PASS / 46 explicit SKIP / 0 FAIL / 0 ERROR across all 94 cases; forward and native Enceladus pass, every unavailable Moon identity and reason is printed and retained |
| `python3 src/test/run_tests.py --test APPLICATION-ENCELADUS-FORWARD-001 --output-dir /tmp/amps-global-enceladus-integration` | PASS: 42/42 cases, zero FAIL/SKIP/ERROR; JSON/JUnit, local receipt, transcript and SVG retained |
| `python3 src/test/run_tests.py --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001 --output-dir /tmp/amps-global-enceladus-native-integration` | PASS: 6/6 after the exact clean native preflight; configured executable/hash, build transcript, local receipt, JSON and JUnit retained |
| `python3 src/test/run_tests.py --test APPLICATION-ENCELADUS-E3-001 --output-dir /tmp/amps-global-enceladus-e3-integration` | SKIP: all 13 E3 identities listed individually because the default profile deliberately leaves `enceladus-e3` unavailable; the 10,000-step simulation was not executed |
| `python3 src/test/run_tests.py --test APPLICATION-MOON-M0-001 --output-dir /tmp/amps-global-moon-m0-integration` | SKIP: all ten I01-I10 identities listed individually with the preflight-conversion blocker |
| `python3 src/test/run_tests.py --test APPLICATION-MOON-REGRESSION-001 --output-dir /tmp/amps-global-moon-registry-integration` | SKIP: all 36 U/I identities listed individually because the default profile deliberately leaves `moon-spice` unavailable |
| `python3 src/test/run_tests.py --case MOON-U01 --config /tmp/amps-moon-machine.json --output-dir /tmp/amps-global-moon-u01-integration` | Clean Moon build PASS; imported U01 ERROR because the restricted sandbox denied all OpenMPI network interfaces; this is retained environment failure evidence, not a test pass |
| unsandboxed `python3 srcMoon/test/U01_los_geometry/test.py --build-dir /home/vtenishe/Mars2/AMPS/build --output-dir /tmp/amps-moon-u01-unsandboxed` | PASS on the exact generated Moon build: all six U01 geometry/output checks true |
| `python3 src/test/run_tests.py --test CORE-HARNESS-001 --output-dir /tmp/amps-app-integration-harness-final4` | PASS: 11/11 after pinning the four-suite `application-smoke` and 37-suite `extended` selections |
| `python3 src/test/validate_traceability.py --output /tmp/amps-app-integration-traceability-final.json` | PASS: 26 areas, 209 translation units, 194 public contracts, 24 configuration contracts, 54 suites and 1,135 cases |
| `python3 src/test/generate_test_catalog.py --check` | PASS: generated catalog contains 54 suites and 1,135 cases |

The complete Moon regression campaign, Moon qualification mode, Moon M0
multi-build campaign, and 10,000-step Enceladus E3 simulation were not executed
by this integration. The earlier selected Moon case establishes its
adapter/build path and separately demonstrates that the sandbox error is
environmental; it does not qualify the other 45 Moon cases. The repository's
current generated tree is the clean Enceladus configuration left by the final
integrated application-smoke run.

## Full-validation selector

The global runner now accepts `--full-validation`. It selects all 54 suites
and 1,135 cases and records the complete campaign union in both
`selection.json` and `amps-core-tests.json`. Unlike ordinary `--all`, it also
sets `execution.full_validation=true` in every adapter context. Enceladus E3
therefore attempts the bounded 10,000-step run, while the Moon U/I aggregate
selects `--campaign-mode qualification`.

The selector is intentionally incompatible with campaign, suite, component,
glob, case, and `--all` selectors. Missing tools, SPICE paths, external data,
or unfinished application gates retain their typed outcomes; full validation
does not weaken prerequisites. Moon I01-I10 remain a visible SKIP until their
multi-build orchestrator satisfies the mandatory native preflight.

Validation evidence:

| Command/evidence | Outcome |
| --- | --- |
| `python3 src/test/run_tests.py --full-validation --list` | PASS: complete manifest-order selection of 54 suites and 1,135 cases |
| `python3 src/test/run_tests.py --full-validation --campaign nightly --list` | Expected exit 2 with a clear incompatible-selector diagnostic |
| `python3 src/test/run_tests.py --full-validation --dry-run --jobs 8 --output-dir /tmp/amps-full-validation-dry-run` | PASS control-plane exercise: one combined JSON/JUnit receipt for all 54 selected owners; every suite is explicitly marked dry-run SKIP and no native build runs |
| `python3 src/test/run_tests.py --test CORE-HARNESS-001 --output-dir /tmp/amps-full-validation-harness` | PASS: 11/11, including complete-selection and narrowing-rejection assertions |

The expensive scientific campaigns were not executed while adding this CLI
control. An actual run requires an explicitly provisioned machine overlay and
is started with `python3 src/test/run_tests.py --full-validation --jobs N
--config /path/to/qualification-machine.json`.

## Contributor-guide checkpoint

The `src/test/README.md` “Add tests to the global procedure” section is now a
standalone contributor workflow rather than a five-step summary. It documents
case-versus-suite ownership, ID and manifest schemas, exact registry/selection
rules, the schema-1 result protocol, PASS/FAIL/SKIP/ERROR semantics,
independent physical/numerical oracles, dependency-light C++ and configured
native tests, MPI/thread/accelerator evidence, portable paths and machine
profiles, coverage/traceability ledgers, exact verification commands, and a
final checklist.

The core-test development roadmap is now co-located with the runner at
`src/test/TEST.md` (moved from `test/TEST.md`); repository documentation and
handoff references use the new path.

`src/test/templates/run_test.py.template` supplies a commented, copyable,
fail-closed adapter. Its unmodified placeholders intentionally publish two
FAIL rows so unfinished tests cannot appear green. Syntax validation and a
synthetic schema-1 execution both pass; the latter returns the required exit 1
and exact two-case receipt. Relative README links, `git diff --check`, and
structural traceability pass. At that checkpoint the registry remained 49
suites and 1,028 cases.

`src/test/TEST_CATALOG.md` was introduced as a generated projection of the
expanded manifest and component inventory. Before the application integration
above, it contained a 49-suite index and 1,028 individual case
rows; every case has nonempty What, How, Why, and Where fields, with maintained
source descriptions used for 728 cases and conservative stable-ID prose for
the remaining 300. `generate_test_catalog.py --check` is enforced by both
structural traceability and `HARNESS-TRACEABILITY-001`. The focused harness
receipt at `/tmp/amps-test-catalog-harness-final` passes 11/11, traceability at
`/tmp/amps-test-catalog-traceability-final.json` passes, and all
catalog-relative links resolve.

## T11/T12 checkpoint (before application integration)

At this checkpoint the registry contained 49 suites and 1,028 cases. Structural
traceability passes with 24 areas, 209 production translation units, 194
public contracts, 24 configuration contracts, and the complete live suite/case
registry. T11's integration mechanism is implemented; T12's bounded matrix is
implemented, but T12 coverage closure is deliberately still open.

T11 contributes six extended suites and 637 identities:

| Imported owner | PASS | FAIL | SKIP | ERROR | Qualification boundary |
| --- | ---: | ---: | ---: | ---: | --- |
| SEP common | 91 | 0 | 52 | 0 | Long pitch-angle Monte-Carlo acceptance is not run by the routine import |
| SWCME | 4 | 0 | 143 | 0 | Four R1 checks execute; missing extended sources and observational runner remain typed |
| staged coronal CME | 219 | 0 | 3 | 0 | Software verification only; production and observational qualification remain false |
| full Corona--SWCME | 13 | 1 | 0 | 1 | `CMBGU05` rejects material-contact leakage; `CSWC0628-MASS-TIME` exceeds the 180-s routine ceiling |
| reduced front plus ambient | 62 | 0 | 35 | 0 | Dependency-light rows execute; native rows use explicit `--skip-native` and no historical PASS is reused |
| legacy disposition audit | 13 | 0 | 0 | 0 | Directory/owner governance only; complete legacy behavioral migration is not claimed |
| **T11 total** | **402** | **1** | **233** | **1** | No imported result qualifies BG3D-4, G1, a full CME, or an observational campaign |

All large registries are explicit in `case_registries.json`. Adapters require
exact ordered central/local equality before applying a one-case selection, and
the harness has a negative narrowed-registry control. The full-volume adapter
preserves typed nonzero markers and treats untyped crashes/setup exits as
ERROR; SEP-common reconciles the pitch runner's exit with its JSON statuses.

T12 adds 16 extended cases without a native AMPS rebuild. The retained global
run reports 10 PASS, five SKIP, one FAIL and zero ERROR. GNU debug/optimized,
Fortran ABI, ASan, UBSan, Valgrind, OpenMP 1/2/4 invariance, NVCC header/kernel
compilation and two broad array/I/O guardrails pass. Clang is unavailable,
TSan lacks a usable host runtime, both installed MPI launchers are denied
sockets by the sandbox, and CUDA has no usable device; these are five separate
SKIPs. `T12-COVERAGE-STRICT-001` truthfully FAILs because only 2 of 24
components are covered and 22 remain incomplete.

Retained evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t11-sep-common-final` | 91 PASS / 52 SKIP after exit/report and exact-registry hardening |
| `/tmp/amps-t11-swcme-global` | 4 PASS / 143 SKIP with complete dormant/observational identities |
| `/tmp/amps-t11-sccm-fixed` | 219 PASS / 3 SKIP, bounded software-verification metadata retained |
| `/tmp/amps-t11-cswc-full-final` | 13 PASS / 1 FAIL / 1 ERROR; architecture still executed and passed after earlier failure |
| `/tmp/amps-t11-reduced-final2` | 62 PASS / 35 native SKIP over the complete 97-row current registry |
| `/tmp/amps-t11-legacy-final` | 13 PASS governance/disposition rows; not behavioral migration qualification |
| `/tmp/amps-t12-portability-final2` | 10 PASS / 5 SKIP / 1 strict-coverage FAIL after compiler-runtime and CUDA evidence-isolation hardening |
| `/tmp/amps-t12-harness-final3` | 11/11 harness checks PASS, including malformed-selection fallback and case-preserving infrastructure ERROR receipts |
| `/tmp/amps-t12-traceability-final3.json` | PASS: 49 suites / 1,028 cases and all structural inventories reconcile |

The preceding bring-up receipt `/tmp/amps-t11-reduced-final` is retained as
failure evidence: exact-order validation rejected a central registry that put
restart rows before the maintained local runner's 1-AU phase. Reordering the
central identity registry to match that local authority produced the final
97-row result above; no reduced-front physics or expected value was changed.

T12 remains open. Required remaining evidence includes whole-AMPS
compiler/sanitizer profiles, a usable Clang and TSan profile, successfully
executed OpenMPI and MPICH (or two approved scheduler profiles), hybrid
MPI/OpenMP and concurrent particle-list behavior, actual supported accelerator
backends, calibrated production mesh/push/collision/halo/restart performance,
and closure or explicit reviewed exemption of all 22 incomplete components.

## T10 checkpoint

T10 is implemented and qualified for its bounded CPU, three-dimensional
Cartesian, fixed-weight neutral scope. The live registry contains 42 suites
and 375 cases. Structural traceability passes with 24 component areas, 209
production translation units, 194 public contracts, and 24 configuration
contracts. The new exclusive owner contributes 16 collective cases, each run
independently at 1, 2, and 4 MPI ranks; all 48 matrix elements require a local
PASS record from every participating rank rather than accepting rank-zero
output alone.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Ownership/ghosts | Replicated `Thread` owner map, exact current-rank load-list membership, allocated received-ghost lists, deterministic ownership checksums | PASS at 1/2/4 ranks; remote intrusive owner lists are intentionally not traversed because their links are reused by received-ghost lists |
| Halos/masks | Full-block center-node Replace/Add, actual owner/destination send relations, face/edge/corner and same/coarse-fine masks | PASS at 1/2/4 ranks on the bounded nonuniform 2:1 fixture; ECSIM remote-stencil consumption remains open |
| Serialization/migration | Associated center/corner data and particle pack/unpack, exact-once ghost-to-owner exchange, weighted block/particle repartition | PASS with exact data and particle identity checks; opaque-field semantic reduction and distributed source injection remain open |
| Boundaries/geometry | Periodic ghost-to-real transfer, outer loss, distributed cut-cell ray hit identities | PASS for Cartesian CPU profiles; non-Cartesian and larger/nonmanifold surfaces remain open |
| Restart/invariance | Same and changed decomposition particle restart; invariant count, physical weight/mass, momentum, kinetic energy, sampled affine rate and commutative checksum | PASS at 1/2/4 ranks; mesh-topology file restart, callback schemas and hybrid per-thread RNG state remain open |
| Portability | CPU-list parsing and rank-local online-CPU discovery | PASS without mutating process affinity; scheduler/cgroup-authorized affinity mutation is deliberately not claimed |

Final retained T10 evidence:

| Command/evidence | Result |
| --- | --- |
| `python3 src/test/run_tests.py --test CORE-MESH-MPI-OWNERSHIP-001 --output-dir /tmp/amps-t10-qualified --jobs 1` | PASS 16/16 cases, 48/48 rank-matrix invocations, 112 participating rank processes; both profiles rebuilt through the mandatory clean workflow |
| `python3 src/test/run_tests.py --test CORE-MESH-MPI-OWNERSHIP-001 --order reverse --output-dir /tmp/amps-t10-qualified-reverse --jobs 1` | PASS 16/16 in reverse order; matrix receipts are byte-identical to normal order and both retained executable SHA-256 values reproduce exactly |
| `/tmp/amps-t10-qualified/artifacts/CORE-MESH-MPI-OWNERSHIP-001/t10-matrix-receipts.json` | Exact 1/2/4-rank receipts, including identical global count/hash/weight/mass/momentum/energy/rate values |
| `/tmp/amps-t10-qualified/artifacts/CORE-MESH-MPI-OWNERSHIP-001/native-t10-mpi-output.txt` | Complete configure/build and rank-local run transcript |
| `/tmp/amps-t10-qualified/artifacts/CORE-MESH-MPI-OWNERSHIP-001/native-build/{baseline,periodic}` | Retained executables, generated headers, selected inputs, SHA-256 file inventory and build identity |
| `python3 src/test/run_tests.py --test CORE-HARNESS-001 --output-dir /tmp/amps-t10-harness-final3 --jobs 1` | PASS 11/11 governance checks after pinning T10's exact campaign membership, three alias-to-case mappings, and generated-configuration contract |
| `python3 src/test/validate_traceability.py --output /tmp/amps-t10-traceability-final.json` | PASS: 24 areas, 209 translation units, 194 public contracts, 24 configuration contracts, 42 suites, 375 cases |

The initial full bring-up at `/tmp/amps-t10-full1` passed 12/15 cases and
exposed three fixture/oracle defects: the send-mask callback had been invoked
outside its real owner/destination caller contract, block-pack validation had
not initialized unrelated owner blocks, and a rank-invariant floating receipt
required bitwise rather than tolerance-bounded reduction agreement. The
contract-corrected focused receipt `/tmp/amps-t10-fix3` passes all three cases
at 1/2/4 ranks. The later qualified run adds the CPU-affinity diagnostic and
is the superseding 16-case authority; no failed bring-up is counted as
qualification evidence.

Both native profiles use the required preflight: confirm repository/process
state, remove only `build`, select the application and SPICE mode through
`Config.pl`, and compile through top-level `make amps` with
`MAKEFLAGS=-j16`. The two test input decks are covered by the repository's
broad `*.input` ignore rule; Git control files were not changed, and a human
preparing a commit must explicitly force-add those decks. No Git state
operation was performed.

## T8/T9 checkpoint

T8 and T9 are implemented and qualified for their bounded one-rank host
scope. The live registry contains 41 suites and 359 cases. Structural
traceability passes with 24 component areas, 209 production translation
units, 194 public contracts, and 23 configuration contracts. T8 contributes
34 cases in seven canonical suites. T9 contributes 13 dependency-light cases
in three canonical suites plus one actual AMPS sampling/particle restart case
selected through `CORE-PIC-RESTART-NATIVE-001` from the baseline native owner.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| ECSIM operators | Native periodic curl, grad-div, Laplacian and curl-curl orders 2/4/6/8; exact constant/affine/quadratic/mixed moments; compatible div-curl and curl-grad identities | PASS for periodic interiors; remote halos are T10 |
| ECSIM boundaries | Periodic/no-row, fixed-field Dirichlet and homogeneous-Neumann increment rows | PASS for supported rows; perfect-conductor component physics is typed unsupported |
| Deposition and fields | CIC charge/qv totals and exact 1-D segment continuity; Poisson CG residual/order; Yee dispersion/staggered invariant; particle-field energy ledger | PASS for the dependency-light one-rank kernels; legacy multi-block GMRES assembly and feedback remain open |
| Field line and auxiliary | Straight/curved interpolation, exchange uniqueness, intrinsic split, conservative repartition, ray clipping, column quadrature, stopping limits and affine background dispatch/interpolation | PASS for local geometry/data contracts; remote exchange belongs to T10 |
| Callbacks and trajectory | Registration order, exactly-once particle processing, deterministic trajectory schema/records and finite-value rejection | PASS for process-local callbacks and base trajectory columns; application optional columns retain their own layout responsibility |
| Persistent envelope | AMPSPF1 magic/version/kind, byte order, pointer/scalar widths, layout ID, length, payload/header CRC, atomic publication and explicit native/canonical/legacy policy | PASS for envelope version 1 and particle payload revision 2 |
| Native restart | Real sampling bytes/counters and an AMR-owned H2O particle round-trip through production save/read, including owner cell, x/v/species, global weight/time step, internal energy, opaque extension bytes and scalar MPI RNG seed | PASS for one-rank MPI; changed-rank ownership, hybrid per-thread RNG and application callback schema identity remain open |
| Data/output | Affine background fixture, concatenated chunks consumed under a changed partition, kind/layout rejection, Tecplot variable order/units/dimensions and nonfinite rejection | PASS for canonical fixture and trajectory base schema; raw legacy background/trajectory files remain same-build-only |

Production corrections include periodic wavenumber and eighth-order curl-curl
coefficients in the native stencil application, checked ECSIM boundary-row
construction, shared CIC/interpolation/geometry/integration/field/field-line
kernels used at production call sites, and checked callback/trajectory
formatting. Sampling and particle restart publication now uses binary rank
parts and an atomically published AMPSPF1 envelope. Particle payload revision
2 restores the global time-step table/initialization flag and scalar MPI RNG
seed in addition to the existing weight and complete particle records.
Unversioned files remain readable through
`PIC::Restart::AllowLegacyRestartFiles`, which defaults true for baseline
compatibility and can be disabled before recovery.

Final retained T8/T9 evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t8-operators-final2` | PASS 5/5 native ECSIM operator/moment/identity cases after a clean configured rebuild |
| `/tmp/amps-t8-final` | PASS 29/29 remaining T8 cases and 6/6 suites in normal order |
| `/tmp/amps-t8-reverse` | PASS 29/29 remaining T8 cases and 6/6 suites in reverse order |
| `/tmp/amps-t9-qualified` | PASS 13/13 persistence/data/output cases and 3/3 suites in normal order |
| `/tmp/amps-t9-qualified-shuffle-809` | PASS 13/13 persistence/data/output cases and 3/3 suites in deterministic shuffled order |
| `/tmp/amps-t9-native-restart-state` | PASS actual one-rank native sampling/particle restart case after a clean configured rebuild |
| `/tmp/amps-t89-harness-qualified` | PASS 11/11 harness schema, campaign, selection, scheduler, protocol and traceability checks |
| `python3 src/test/validate_traceability.py` | PASS: 24 areas, 209 translation units, 194 public contracts, 23 configuration contracts, 41 suites, 359 cases |

The first post-integration ECSIM rebuild failed because the established user-
particle hook is a function-like configuration macro and was incorrectly
passed as a callable. The production bridge now invokes the macro inside an
exactly-once lambda; `/tmp/amps-t8-operators-final2` is the superseding clean
receipt. The earlier failure is discovery evidence, not a qualified result.

Every configured T8/T9 run performed the mandatory root/process preflight,
removed only `build`, selected its application through `Config.pl`, and built
through top-level `make amps` with `MAKEFLAGS=-j16`. No Git state operation was
performed.

The qualification boundary is explicit. The old diagnostic-only linear-
solver, fast-wave, light-wave and external-data reader programs were not
relabeled as passing applications; their reproducible numerical contracts are
now asserted in production-linked small fixtures. Component-aware conductor
physics, legacy multi-block assembly, particle-field feedback, legacy-format
migration, full native application continuation and distributed restart/
field-line/halo ownership remain later work.

## T7 checkpoint

T7 is complete for the bounded one-rank spherical scope selected for this
stage. The live registry contains 31 suites and 311 cases; structural
traceability passes with 24 component areas, 209 production translation
units, 186 public contracts, and 23 configuration contracts. The final T7
campaign reports 29 PASS, no FAIL/SKIP/ERROR, and six passing suites. Its 26
dependency-light cases also pass in normal, reverse, and deterministic
seed-707 shuffled order.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Spherical surface | Flat/index round trips, positive exact curved areas, total `4*pi*R^2`, unit normals, default inventory/sample-buffer semantics, deterministic accommodation and adsorption-minus-desorption exchange | PASS for spherical one-rank contracts; nonspherical surfaces, target chemistry/temperature laws and MPI reduction remain open |
| Exosphere sources | Area-weighted shapes, uniform and default Knudsen-cosine laws, stationary Poisson event times, cumulative process selection, physical/sampled ledgers, and production particle insertion | PASS for the five target-independent compiled spherical modes through native AMR ownership and the mover; target-specific, external-domain, background-plasma, SPICE and MPI modes remain open |
| Sputtering | Rubin/Teolis boundaries and interpolation; product catalogue; H2O/O2 speed/weight laws; selectable angles with Knudsen-cosine default; incident-rate, daughter-weight and number/mass/momentum/energy ledgers | PASS for the shared source layer; native sputtered-particle insertion and exact page-level provenance of inherited numeric tables remain open |
| Dust kernels | Log bins/power laws, SI potential/Lorentz limits, spherical direction/mass, Fink Mie interpolation, weighted moments, and dust-owned individual-weight split/merge within explicit exact-intrinsic size groups | PASS for the deterministic physical and conservation APIs; the two-output merge cannot generally preserve the pressure tensor, higher moments or raw speed moment |
| Dust charging and injection | Independent Horanyi-current fixtures at three potentials and one equilibrium root; transient convergence for four radii; translated spherical production injection, DUST records, native AMR owner/mover, sampled mass rate and signed reservoir closure | PASS for bounded one-rank profiles; multi-rank arbitration/migration/reduction, restart and long stochastic convergence remain open |

Production corrections include bounds-safe Rubin/Teolis interpolation, the
Teolis interval-offset correction, a single shared sputtering product
catalogue with idempotent initialization, dimensionally correct Mie bounds in
log size parameter, dimension-safe dust injection normalization, translated
sphere geometry, initialized particle records and checked AMR ownership, the
removable `p=1` dust-size limit, and shared checked dust bin/force/moment
kernels used by production call sites. The legacy exosphere injector now uses
the shared exact angular inverse CDFs, stationary Poisson scheduler and
cumulative compiled-process selector. New headers and tests document
equations, SI units, coordinate conventions, numerical bounds, ownership and
unsupported domains.

Final retained T7 evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t7-final` | PASS 29/29 cases and 6/6 suites; no FAIL/SKIP/ERROR |
| `/tmp/amps-t7-kernels-current` | PASS 26/26 dependency-light cases in normal order |
| `/tmp/amps-t7-kernels-reverse` | PASS 26/26 dependency-light cases in reverse order |
| `/tmp/amps-t7-kernels-shuffle-707` | PASS 26/26 dependency-light cases in deterministic shuffled order |
| `/tmp/amps-t7-exosphere-native-final` | PASS native five-mode spherical exosphere source, ownership and ledger case after final relative-rate closure hardening |
| `/tmp/amps-t7-dust-charging-native-current` | PASS native signed-current/equilibrium/transient charging case |
| `/tmp/amps-t7-dust-injection-native-current3` | PASS native translated-sphere dust injection and mass-ledger case |
| `/tmp/amps-t7-harness-final` | PASS 11/11 harness schema, campaign, selection, scheduler, protocol and traceability checks |
| `python3 src/test/validate_traceability.py` | PASS: 24 areas, 209 translation units, 186 public contracts, 23 configuration contracts, 31 suites, 311 cases |

The final native exosphere transcript emitted 618 numerical particles in four
source intervals: 122 impact-vaporization, 141 PSD, 112 thermal-desorption,
122 sputtering and 121 vertical-injection records. The accumulated calculated,
surface-sampled and global production ledgers each report
`6.1799999999999988e11 s^-1`; source occurrence counts are observability checks,
not rate estimators. The dust injection receipt emitted one represented grain
in the first half-grain interval and none in the second, closing the signed
reservoir from `-4.1887902047863901e-17 kg` to exactly zero while reporting
`8.3775804095727796e-16 kg s^-1` sampled mass rate.

Earlier dust receipts that reported two capability SKIPs are retained only as
diagnostics. They motivated the dedicated DUST-species native injection
profile and the dust-owned size-group split/merge API; the all-green final
campaign supersedes them. Likewise, the first native injection attempt exposed
radius-record initialization order and the second exposed a floating endpoint
at `exp(log(a_max))`; both defects were corrected before the retained receipt.

Every configured T7 run performed the mandatory root/process preflight,
removed only `build`, selected its application through `Config.pl`, and
compiled through top-level `make amps` with `MAKEFLAGS=-j16`. Focused
warnings-as-errors and ASan/UBSan runs passed the dependency-light dust and
sputtering/exosphere kernels; LeakSanitizer is unavailable under this
environment's ptrace restrictions. No Git state operation was performed.

## T6 checkpoint

T6 is qualified for the implemented, bounded reaction domains. The live
registry at this checkpoint contains 25 suites and 282 cases, and structural
traceability passes with 24 component areas, 209 production translation units,
181 public contracts, and 20 configuration contracts. The gate covers the
dependency-light chemical integrator, configured H2O photolysis execution,
one individual-weight daughter-production profile, Burger electron-impact
rates/products, and the two compiled charge-exchange rate laws.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Chemical integrator | `A -> B`, reversible `A <-> B`, stiff decay, coupled chain, positivity/conservation, invalid networks | PASS for the dependency-light finite network API; runtime deck assembly and cell coupling remain open |
| Photolysis execution | Populated-cell parent survival/decay, intrusive-list integrity, disabled/invalid lifetime paths, idempotent initialization | PASS for configured H2O quiet-Sun execution; application multi-parent networks remain open |
| Photolysis daughters | Individual-weight H2O channel, daughter identity/count, weight, charge, mass, momentum and kinetic-energy accounting | PASS for the selected one-channel native profile; MPI/threaded insertion and general runtime networks remain open |
| Charge exchange | Lindsay--Stebbings and Maher--Tinsley formulas, finite-domain values and invalid-domain sentinel | PASS for the exposed H/H+ lifetime functions; event selection, products and particle weights are not exposed by this module |
| Electron impact | Burger total rates, thresholds, eligible-channel selection and H2O product statistics | PASS for compiled parent-rate dispatch and H2O products; runtime electron coupling and other-parent product insertion remain open |

The populated-cell fixture now calls
`PIC::DomainBlockDecomposition::UpdateBlockTable()` after allocating owner
blocks. This is a substantive test correction: without rebuilding that table,
whole-mesh chemistry processors observe `nLocalBlocks == 0` and a nominally
passing parent-decay case is a no-op. The final receipts therefore supersede
the earlier diagnostic runs.

Final retained T6 evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t6-chemical-integrated` | PASS 7/7 chemical-solver cases |
| `/tmp/amps-t6-photochem-baseline-blocktable` | PASS 6/6 populated-cell photochemistry cases |
| `/tmp/amps-t6-photochem-daughter-blocktable` | PASS 1/1 individual-weight daughter case |
| `/tmp/amps-t6-reactions-build2` | PASS 7/7 charge-exchange/electron-impact cases |
| `/tmp/amps-t6-gate-reverse` | PASS 21/21 cases and 4/4 suites in reverse order |
| `/tmp/amps-t6-gate-shuffle-706` | PASS 21/21 cases and 4/4 suites in deterministic shuffled order |
| `/tmp/amps-t6-harness-prebuild-2` | PASS 11/11 harness self-tests |
| `python3 src/test/validate_traceability.py` | PASS: 24 areas, 209 translation units, 181 public contracts, 20 configuration contracts, 25 suites, 282 cases |

Each native T6 owner used the mandatory root/process preflight, removed only
`build`, selected its deck through `Config.pl`, and compiled through top-level
`make amps` with `MAKEFLAGS=-j16`. No Git state operation was performed.

## T4/T5 checkpoint

The live registry contains 22 suites and 261 cases. Structural traceability
passes with 24 component areas, 209 production translation units, 181 public
contracts, and 19 configuration contracts. T4 is implemented and qualified
for the bounded one-rank profiles below but remains partial against the full
roadmap acceptance gate; T5 satisfies its supported-capability gate while
retaining typed unsupported domains.

| Area | Implemented and directly tested | State / remaining boundary |
| --- | --- | --- |
| Charged movers | Uniform-E velocity/position convergence; uniform-B phase/energy convergence; crossed-field drift; relativistic gamma and gyrofrequency | PASS for constant-field Boris profiles; coarse/fine crossing and additional production movers remain open |
| Guiding center | Parallel electric/mirror force, `E x B`, grad-B and curvature drift; neutral, zero-B, negative-mu and non-finite statuses | PASS for the field-local RHS; full orbit integration, adiabatic validity and convergence remain open |
| Boundaries | Static sphere reflection/absorption/tangent/start-inside, outer loss, exact moving plane/sphere, periodic real/ghost wrap and data sync | PASS for registered modes; circle/rotation-body/Nastran mover dispatch, diffuse law, multiple boundaries and MPI migration remain open |
| Injection | Maxwell half-space flux; Poisson waits; uniform cell, triangle, and ball positions | PASS for distribution kernels; ownership, species-rate/weight normalization and distributed injection remain open |
| Sampling | Exact density, flux, full pressure tensor, energy, speed, temperature, rescaling invariance, collect/swap/flush | PASS for raw cell moments; histogram under/overflow and distributed reduction remain open |
| Collision selectors | Same independent rate/conservation fixture for NTC, simplified NTC, MF, simplified MF, and Yinsi; zero/edge/limiting paths; unequal global W/dt on advanced selectors | PASS in declared domains; simplified selectors remain fixed-weight-only and Yinsi record weights are typed unimplemented |
| Molecular laws | Generated VHS/VSS power tables; Maxwell quadrature; VSS alpha moments and elastic invariants | PASS for the two compile profiles; broader transport calibration and alternate mixture laws remain open |
| LB internal energy | RT hot/cold equilibrium, update flags, mode families, vibrational limits, generated Zrot 2/5/10 | PASS for RT; unequal rotating-pair Zrot is typed ambiguous, VT is unverified/terminating, VV is absent |

Production corrections made under these tests include the previously omitted
VHS/VSS power-index emission, VSS angular dispatch for every selector,
unbiased per-particle collision limiting without Yinsi global-threshold
mutation, generated rather than literal `Zrot`, a typed guiding-center RHS,
and shared checked Poisson/uniform-ball source samplers. New code and fixtures
carry equations, SI units, frames, stochastic tolerances, validity limits, and
typed-failure comments.

Final retained evidence:

| Command/evidence | Result |
| --- | --- |
| `/tmp/amps-t45-baseline-final-normal` | PASS 51/51, one clean baseline build |
| `/tmp/amps-t45-baseline-final-reverse` | PASS 51/51, reverse case order |
| `/tmp/amps-t45-baseline-final-shuffle-1729` | PASS 51/51, deterministic shuffled order |
| `/tmp/amps-t45-profile-variants-final` | PASS 7/7 cases and 6/6 suites: periodic, VHS, VSS, Zrot 2, Zrot 10, individual-weight NTC/MF |
| `/tmp/amps-t2-t5-profile-qualified` | PASS 194/194 cases and 12/12 suites on the corrected final tree, covering every T2-T5 configured profile |
| `python3 src/test/validate_traceability.py --output /tmp/amps-t45-traceability-pretest.json` | PASS: 24 areas, 209 translation units, 181 public contracts, 19 configuration contracts, 22 suites, 261 cases |

The native adapter performed the mandatory root/process preflight, removed
only `build`, selected each application with `Config.pl`, and compiled with
top-level `make amps` and `MAKEFLAGS=-j16` for every retained build. No Git
state operation was performed.

The first full-profile diagnostic at `/tmp/amps-t2-t5-profile-final`
reported 193 PASS and one FAIL: the pre-T4 species oracle still asserted a
seven-neutral deck after H+ and electrons were appended for the charged-mover
profile. The oracle now documents and checks all nine symbols, indices,
masses, charges, and gas types. `/tmp/amps-t45-species-nine-profile-fix`
passes the focused 5/5 correction, and the complete
`/tmp/amps-t2-t5-profile-qualified` rerun supersedes the diagnostic with
194/194 PASS. The failed receipt is retained as discovery evidence and is not
reported as qualification.

## Implemented

- Portable global runner with validated registry, isolated work/artifact
  directories, timeouts, typed PASS/FAIL/SKIP/ERROR states, JSON, and JUnit.
- The global runner now ends every execution with a stable
  `FINAL INDIVIDUAL TEST RESULTS` table.  It prints one PASS/FAIL/SKIP/ERROR
  row for each executed case, then separate individual-case and suite
  summaries.  The JSON receipt preserves the same normalized rows in
  `individual_results` and their totals in `individual_counts`; suites which
  cannot publish case evidence are explicitly labeled as suite-level rows.
- Layered machine configuration with external JSON profiles, environment
  overrides, relative-path resolution, and no hard-coded installation paths.
- Explicit component inventory with structural drift audit and a separate
  strict-completeness gate. `source_traceability.json` now classifies every
  discovered production translation unit and public compile contract, while
  `configuration_contracts.json` records the machine, Config.pl, generated-
  file, compile-profile, and runtime-mode boundaries.
- The dependency-free `validate_traceability.py` expands those rules, rejects
  unclassified/overlapping files and stale suite/case/local-registry owners,
  and cross-checks the 24-area component inventory. This is structural
  governance evidence, not behavioral physics coverage.
- Eleven registered harness self-tests cover the original configuration,
  registry, selection, protocol, scheduler, and report contracts plus schema,
  hostile-subprocess failures, real timeout and SIGINT process-tree cleanup,
  bounded concurrency/fail-fast, and the full source/configuration
  traceability contract.
- Generic runtime-input parser adapter (four existing independent fixtures).
- General utility tests for arrays, bit arrays, stack allocation/reuse, linked
  containers, string helpers, and sparse 4-D indexing.
- T1 adds 10 general-numerics cases for asymmetric arrays 1-D through 5-D,
  empty/copy/move/alias/resize/bounds behavior, managed-stack reuse,
  special-function identities, and independent quadrature exactness,
  convergence, and invalid-domain checks.
- T1 adds seven deterministic RNG cases for recurrence/reseeding, endpoint
  exclusion, independent interleaved streams, uniform/normal/gamma moments,
  and isotropic directions, plus a separate two-rank RNG-stream case.
- T1 adds four two-rank `CMPI_channel` payload cases, four diagnostic-logger
  cases, five SI/CGS/normalized-units cases, and four He/Na radiation-table
  cases. The associated production hardening covers safe multidimensional
  array ownership/extent arithmetic, complete hybrid RNG reseeding, checked
  MPI framing/alignment, finite and representable unit scales, Gamma/bulk-
  normal domains, GCD `INT_MIN`, finite radiation-table clamping, and invalid
  or unrepresentable heliocentric inputs.
- T1 adds `CORE-INTERFACE-ABI-001`: six executable C/Fortran/Geopack cases and
  five explicit optional-provider gates. ARMS, BATSRUS, CCMC, Kameleon, and
  SWMF remain `SKIP` until independent libraries and versioned fixtures are
  configured; these skips are not behavioral passes.
- New dependency-free C++ harness shared by standalone and configured tests.
- One combined `CORE-PIC-NATIVE-001` registered owner for the configured
  `test/srcCoreUnitTests` application. It replaces two executable native suite
  entries while preserving their old IDs as selection aliases, so compatible
  particle-buffer and physics cases require one configured build. T1 extends
  it from 24 to 31 cases and adds `CORE-SPECIES-001` as a third alias for
  duplicate-species rejection, compiled constants, generated mapping/tables,
  and unknown-name diagnostics. Two additional production-linked cases cover
  transactional/idempotent global AMU/charge initialization and conversion-
  table output.
- T2 extends the baseline exclusive native owner from 31 to 43 cases and adds
  five profile owners, for six suites and 85 native profile cases. The matrix
  covers packed/aligned-16 and fixed/individual-weight records,
  species-global/species-local/single-global weight and time-step modes, and
  packed/aligned split/merge. Production record requests now enforce typed
  power-of-two alignment and checked offset/stride arithmetic; record reuse
  clears stale payload; invalid host pointers are contained; deletion marks
  the addressed record; and host weight setters require finite positive
  values before accessing storage. Species-local distribution copies rescale
  physical density when weight or time step changes.
- T2 split/merge cases assert model count, physical weight, species, weighted
  position, momentum, translational energy, rotational energy, and each
  vibrational mode. The central-variance implementation avoids catastrophic
  cancellation. Opaque application bytes follow a clone/survivor policy
  because generic PIC code has no field semantics. Buffer exhaustion is
  isolated per MPI subprocess and must exit nonzero with the expected
  diagnostic.
- T3 adds `CORE-MESH-NATIVE-001`, a 58-case owner that performs clean 1-D,
  2-D, and 3-D configured builds. Fourteen topology/adaptation cases execute in
  each dimension; eight production interpolation and eight cut-cell/NASTRAN
  cases execute in 3-D. The former roadmap suite IDs remain selection aliases.
  Production hardening covers lower-dimensional AMR/PIC compilation, safe
  split/coarsen data transfer, full neighbor-slot interpolation selection,
  affine reproduction across all 256 legacy refinement masks, open-interior
  cut-surface contact, interface-error-driven integration depth, translated
  sphere recursion, ray/edge metrics, and sparse NASTRAN identifiers.
- Expanded 19-case native physics suite:
  - four photochemistry kernel cases for Huebner H2O rates/yield/radial
    scaling, exponential decay probability/event time, product momentum/
    energy and two-body angular isotropy, and the seven-channel branching
    distribution plus neutral-only compiled-product availability;
  - two hard-sphere collision-kernel cases for unequal-mass momentum/energy/
    relative-speed conservation and spherical isotropy;
  - one drifting-Maxwellian distribution case for component means, analytic
    mean thermal speed, and recovered kinetic temperature;
  - six Larsen-Borgnakke cases for optional particle storage, rotational
    equipartition, vibrational harmonic-oscillator means, beta energy
    partition, two-stage RT redistribution/conservation, and the analytic
    `1/Z_rot=1/5` relaxation occurrence rate;
  - three mover cases for neutral ballistic Simple/Boris limits, destination-
    cell insertion, and an analytic specular reflection at a registered sphere;
  - three populated 3-D-cell cases: NTC/LB H2O rotational equilibration with
    momentum/total-energy conservation, the independent Maxwellian H2O/O
    collision rate `n0*n1*sigma*<g>`, and 2,400-particle photolytic parent
    survival with no daughter production.
- Native physics statistics use fixed case-local seeds, documented sample
  counts, analytical expectations, and explicit confidence/tolerance bounds.
- The photoproduct angular test exposed and fixed uniform-polar-angle sampling
  in `GetProductVelocity_2`; the implementation now samples `cos(theta)`
  uniformly for the correct solid-angle measure.
- The neutral-only branching test exposed and fixed a second production defect:
  product IDs were compacted around unavailable charged daughters while their
  velocity vectors retained original channel slots. Both ordinary and
  given-product generators now compact IDs and velocities together, with a
  NaN-sentinel regression for a surviving slot-1 daughter.
- Native-run safety includes a repository-specific advisory lock, an active
  build/run preflight, exact validation of the generated `build` target,
  application selection through `Config.pl`, fixed `MAKEFLAGS=-j16`, and a
  retained configure/build/MPI transcript.
- The native lock now uses the shared checkout's read-only `Config.pl` inode,
  so cooperating jobs on different nodes or with different `TMPDIR` values
  cannot reconfigure the same checkout concurrently. Explicit output reuse
  invalidates the old per-suite receipt before launch; an adapter must publish
  fresh evidence. PASS receipts must contain at least one executed case.
- Exact order-independent case registries detect a missing, renamed,
  duplicated, or unlinked C++ test even when the remaining executable returns
  success, without relying on unspecified cross-translation-unit static
  initialization order.
- Exact `--case` selection and deterministic `--order normal|reverse|shuffle`
  with a recorded shuffle seed are implemented. The harness checks selection,
  aliasing, and deterministic order construction. The retained T0
  normal/reverse/shuffle campaigns each pass their historical 45-case
  registry. All 21 T0
  dependency-light cases also pass through one-case global selections; the 22
  executable native cases are launched as exact one-case MPI processes from
  one clean build, while the two configuration cases execute as individually
  identified adapter cases.
- The global runner accepts `--jobs N` (default `1`) and executes only
  consecutive nonexclusive suites concurrently. Every `exclusive: true`
  native suite has a leading and trailing barrier, so repository-wide
  configuration/build state cannot overlap another suite. Results are restored
  to manifest order; concurrent case-start lines carry their suite ID.
- The runner now streams `[CASE] RUN <ID> <description>` immediately before
  every Python, dependency-light C++, parser, and configured native case. The
  terminal shows case-level progress while full compiler/MPI output remains in
  the suite log; streaming retains timeout enforcement and kills the complete
  POSIX child process group on timeout or user interrupt. The harness exercises
  both paths with a real descendant process.
- The legacy `test/srcUnitTests` entry point now returns GoogleTest's result
  instead of reporting process success after failed assertions.
- `SRC_UNIT_TESTS_MIGRATION.md` maps every legacy source/case to its registered
  replacement or an explicit open item. It records the legacy out-of-bounds RT
  assertion, allocator-mismatched fixture teardown, non-reset statistical
  accumulators, and still-unported individual-weight split/merge profile.
- `array_1d` now has explicit deep-copy/resizing-copy and ownership-transferring
  move semantics. The container tests exposed that the previous implicit copy
  constructor aliased owned storage and could double-delete when copy elision
  was unavailable.
- CRC diagnostic formatting no longer uses overlapping source/destination
  `sprintf`; bounded append formatting preserves termination on long paths.
- The manifest distinguishes smoke, nightly, profile, extended, distributed,
  and application-smoke campaigns. The retained T1 qualification registry
  assigned eight suites and 51 cases to smoke and 11 suites/98 cases to
  nightly. At the T3 checkpoint the registry had 17 executable suites and 210
  cases: smoke contained 51, nightly/`--all` contained 210, `profile`
  contained seven exclusive native suites and 143 cases, and `extended`
  contained 125. The current T7 registry instead has 31/311 total and 17/211
  in `profile`; `distributed` still retains its five-case MPI selection. The
  harness pins the current distinct sets. Empty campaign selections do not
  imply implemented coverage; every native profile remains a scheduler
  barrier wherever selected.
- Framework, suite-adapter, native-fixture, and core-change comments document
  protocol/failure semantics, path containment, portability, concurrency,
  memory ownership, MPI reduction, test oracles, and validity limits. The
  framework/design and configured-application READMEs carry the matching
  maintainer and machine-configuration contracts.

## Not yet complete

Only `src/general/input` and `src/pic/units` are currently marked `covered`.
`src/species` has complete computational path coverage for this T1 scope but
remains `partial` until a separately versioned paper/digitization fixture
independently revalidates the retained He/Na tables. The test harness is
`partial`, and all other `partial` and
`planned` entries in
`component_inventory.json` remain open. Existing model/application runners
outside the adapted T7 spherical-exosphere and dust applications still remain
to be integrated;
distributed mesh ownership, restart, and rank-invariance coverage are open.
T2 closes the declared host particle-layout, full-buffer, weight/time-mode,
individual-weight, and split/merge contracts, but not every use of those
facilities by later physics.

The PIC and photolytic-reaction entries remain `partial`. Production NTC cell
selection, a cross-species analytic rate, and H2O RT equilibration now have
one-rank populated-cell evidence. Species-global, species-local, and
single-global weight/time modes, individual particle weights, local rescaling,
and conservative split/merge now have T2 host evidence. T5 covers all five
collision selectors in their declared weight domains, collision limiting, and
the configured VHS/VSS laws; Yinsi individual-record weights and simplified-
selector unequal-weight extensions remain open. T4 covers constant-field
charged movers plus bounded static/moving/periodic surfaces, but coarse/fine
crossings, additional production surface types and rank migration remain open.
T6 covers three-body photolysis angular statistics and one individual-weight
H2O daughter-insertion profile; quiet-Sun/runtime networks and MPI/thread-safe
daughter processing remain open. VT/VV relaxation is not claimed because the
current VT implementation explicitly exits as unverified.

The RT occurrence case deliberately qualifies the current LB closure's
hard-coded `Zrot=5`. The input generator emits and initialization validates
`RotationZnumber`, but `LB::RedistributeEnergy` does not currently use that
array. Selecting a cross-species combination rule is a production-physics
decision, so this test campaign records the discrepancy instead of inventing
one. A later correction must set a non-5 input and verify the configured rate.

`HARNESS-TRACEABILITY-001` imports and runs the same exhaustive validation used
by the standalone `validate_traceability.py` command. It expands and audits
every current production translation unit/public contract and validates the
configuration ledger, so structural drift is part of smoke/nightly execution.
The explicit standalone command remains useful for retaining the expanded
per-file JSON receipt.

T1 is qualified for its declared bounded contracts. The stale campaign oracle
from the first smoke run was corrected; `/tmp/amps-t1-smoke-qualified` passes
51/51. The final 98-case normal, reverse, and seed-1729 shuffled campaigns each
report 93 PASS, five typed optional-provider SKIPs, and no FAIL/ERROR. This
qualification does not promote the explicit `partial`, `planned`, or
`unavailable` contracts described above.

T2 is qualified for the declared one-rank host profile matrix. The baseline
owner passes 43/43, and the six-suite profile campaign passes 85/85 in reverse
and deterministic seed-1729 shuffled orders. The matrix covers four record
layouts, three supported weight/time organizations, individual finite-positive
weight validation, local rescaling, and packed/aligned conservative
split/merge. Earlier species-global sentinel restoration and the initial T2
layout/profile fixture failures remain retained diagnostics; each was
corrected without weakening the production invariant. CUDA/device diagnostics,
semantic reduction of opaque fields, unequal-weight collision physics,
injection/sampling normalization, restart, and distributed ownership remain
outside T2. The separately implemented T3 scope is described below.

T3 is qualified for its declared one-rank Cartesian scope. The first normal,
reverse, and seed-1729 shuffled 58/58 receipts were followed by an independent
audit that found and corrected one additional lower-dimensional array-bound
defect, so those receipts are retained only as superseded bring-up evidence.
The corrected tree then passed all 58 cases in normal, reverse, and seed-1729
shuffled order with identical per-profile executable hashes. The scope does not
include MPI owner/ghost exchange, restart, non-Cartesian coordinates, arbitrary
application-associated data, deep imbalance, or general nonmanifold/open/
self-intersecting and multi-rank surfaces. The legacy 3-D-hard-coded
`OutputLayerNodeMap` is not exercised or qualified in the lower-dimensional
profiles.

The lower-dimensional bring-up was substantive rather than a source-only
compile: dimension-bounded mesh-signature and associated-data loops, compact
neighbor access, the 2-D neighbor call signature, dimension-generic mesh and
node-ID types, mover local scope, and always-three-component field-line vectors
were corrected before 1-D/2-D cases could execute. The final audit also bounded
global-domain initialization and `buildMesh` reads/writes to active dimensions,
with explicitly zeroed three-component diagnostics. The per-dimension uniform
case now executes production `buildMesh` on a private topology-only mesh and
asserts its root split, generation/maximum level, and unique `2^D` bottom list.
The interpolation audit then corrected the coarse-surrounded-by-fine `<=`
selector, searched every neighbor
slot before selecting the coarsest support, and repaired the edge-1/edge-3
predicates. Cut-cell work distinguishes open-interior intersection from face
contact and selects recursive depth from the interface error `O(A*h)`; sphere
subdivision preserves unsplit bounds, and NASTRAN lookup supports sparse,
out-of-order IDs.

`make -C src/test strict-coverage` is expected to fail until those items are
implemented and qualified. A passing smoke run must not be described as
complete AMPS coverage.

## Current validation evidence

The live source tree registers 49 executable suites and 1,028 individual
cases. Its `profile` campaign selects 19 configured owners and 233 cases,
including the expanded baseline owner, historical T2 variants, the 58-case T3
dimensional mesh owner, T4/T5 compile variants, T6 reaction profiles, three
T7 applications, the T8 periodic ECSIM profile, and the T10 collective MPI
owner. The earlier seven-suite/143-case combined T2+T3
profile remains retained historical evidence rather than the current campaign
membership. The historical bounded T1
qualification remains 11 suites/98 cases with 31 native cases; its normal,
reverse, and shuffled native rebuilds retain the identical executable SHA-256
`17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`.
T2 has a 43/43 normal baseline receipt plus 85/85 reverse and deterministic
shuffled profile receipts. T3 has final corrected-tree 58/58 normal, reverse,
and seed-1729 shuffled receipts; its first three 58/58 receipts remain
superseded bring-up evidence. The final-tree combined profile regression passes
all 143 T2+T3 configured cases. T0's earlier 45-case receipts remain
historical evidence, and the retained T1 receipts remain the authority for T1
itself.

| Date | Command/evidence | Outcome | Scope |
|---|---|---|---|
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order normal --output-dir /tmp/amps-t2-t3-profile-qualified` | PASS 143/143, 7/7 suites, 0 FAIL/SKIP/ERROR (426.294 s) | Historical T2+T3 checkpoint of all six T2 configured owners (85 cases) plus the T3 mesh owner (58 cases); the live profile campaign has since expanded. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal --output-dir /tmp/amps-t3-qualified-normal` | PASS 58/58 | Final corrected-tree normal order. Executable SHA-256: 1-D `54e79ddf1c34e2910a059a23e38232bea6322ab7a1ee6d31e35204acd8dbf4e0`; 2-D `8d5311b99bfba3812599336b2cdbb4fbcd53ab16a3e625122b6db44a30b4aad9`; 3-D `96f13a4be39ed221b7a252bb8313db3371ed49b2fd464ea1b5d09a545bf42e34`. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse --output-dir /tmp/amps-t3-qualified-reverse` | PASS 58/58 | Final corrected-tree reverse order; all three executable hashes match the normal receipt exactly. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t3-qualified-shuffle-1729` | PASS 58/58 | Final corrected-tree deterministic shuffled order; all three executable hashes match normal and reverse exactly. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order normal --output-dir /tmp/amps-t3-native-normal` | Superseded bring-up PASS 58/58 | All three dimensional profiles passed, but an independent audit subsequently corrected a lower-dimensional array-bound defect; superseded by the qualified rerun above. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order reverse --output-dir /tmp/amps-t3-native-reverse` | Superseded bring-up PASS 58/58 | Reverse case order across the same three clean profile builds; not final authority after the post-campaign fix. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-MESH-NATIVE-001 --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t3-native-shuffle-1729` | Superseded bring-up PASS 58/58 | Deterministic shuffled order before the audit correction; final authority is the qualified receipt above. |
| 2026-11-07 | `/tmp/amps-t3-final-fixes-fourth` | PASS 4/4 focused | Final 3-D interpolation-mask and planar cut-cell regressions before the independent lower-D audit. |
| 2026-11-07 | `/tmp/amps-t3-1d-escalated-initial`, `/tmp/amps-t3-2d-second` | PASS 6/6 + 6/6 focused | Search, neighbors, coarsen, data transfer, and cycles after lower-dimensional compilation/neighbor hardening. These are focused bring-up receipts, not the final matrix. |
| 2026-11-07 | `/tmp/amps-t3-3d-initial`, `/tmp/amps-t3-3d-second` | Historical FAIL 16/21, then 3/5 | Exposed missing interpolation configuration, sphere-recursion bounds, stale/incorrect AMR interpolation selection, and planar cut-cell contact/depth errors; later focused and 58-case bring-up runs passed. |
| 2026-11-07 | `python3 src/test/validate_traceability.py --output /tmp/amps-t3-doc-audit.json` | PASS structural audit | 24 areas, 209 production translation units, 181 public contracts, 16 configuration contracts, 17 suites, and 210 cases. This structural receipt is distinct from the final behavioral campaigns. |
| 2026-11-07 | `python3 src/test/run_tests.py --test CORE-PIC-NATIVE-001 --output-dir /tmp/amps-t2-native-full-normal` | PASS 43/43 | Complete baseline packed/fixed, species-global native owner in normal order. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order reverse --output-dir /tmp/amps-t2-profile-reverse` | PASS 85/85, 6/6 suites | All T2 generated profiles and cases in reverse suite/case order; each owner performed its required clean rebuild. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign profile --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t2-profile-shuffle-1729` | PASS 85/85, 6/6 suites | All T2 generated profiles and cases in deterministic shuffled order. |
| 2026-11-07 | `/tmp/amps-t2-layout-packed-fixed-recheck` | PASS 7/7 | Packed fixed record offsets, clone/reuse, list, capacity, invalid access, and exhaustion after correcting test offset types to the checked `long int&` API. |
| 2026-11-07 | `/tmp/amps-t2-split-merge-packed-recheck` | PASS 13/13 | Packed individual layout and split/merge after making the two-particle bin reachable and keeping the high-bulk fixture below the debugger speed guard. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-fixed-initial` | PASS 7/7 | Aligned-16 fixed record profile with mixed-size anonymous extensions. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-individual-recheck` | PASS 13/13 | Aligned-16 individual layout and split/merge after the opaque-state oracle accepted an exact clone of either valid source tuple. |
| 2026-11-07 | `/tmp/amps-t2-weight-local-recheck` | PASS 5/5 | Species-local block-associated weight/time and rescaling after adding the target-free exosphere generator block required by mesh headers. |
| 2026-11-07 | `/tmp/amps-t2-weight-single-global-initial` | PASS 4/4 | Single-global round trips, cross-species slot aliasing, and invalid host inputs. |
| 2026-11-07 | `/tmp/amps-t2-harness-final-prebuild` | PASS 11/11 | Harness registry, campaign-count, schema, protocol, subprocess, and traceability regression before the final native profile campaign. |
| 2026-11-07 | `/tmp/amps-t2-layout-packed-fixed-initial` | Historical build failure | The fixture used `int` momentum offsets with a production API requiring checked `long int&`; the type mismatch was corrected and the seven cases passed in the recheck. |
| 2026-11-07 | `/tmp/amps-t2-split-merge-packed-initial` | Historical FAIL: 11 PASS, 2 FAIL | Independent fixture exposed an unreachable two-particle split bin and a merge bulk speed above the debugger guard; corrected with identical-velocity input and a `1e9` bulk offset. |
| 2026-11-07 | `/tmp/amps-t2-layout-aligned-individual-initial` | Historical FAIL: 12 PASS, 1 FAIL | Oracle over-specified which opaque tuple member was cloned; corrected to accept either exact existing record under the documented clone/survivor policy. |
| 2026-11-07 | `/tmp/amps-t2-weight-local-initial` | Historical build failure | Minimal local deck omitted the exosphere-generated internal-boundary base classes; corrected by adding the target-free exosphere block. |
| 2026-11-07 | `/tmp/amps-t2-weight-time-initial` | PASS 5/5 | Focused species-global weight/time API cases from one clean native build; executable SHA-256 `46678c6054daad2ff9b7e3a1aabce9792d12139515e02f93bfaa0296210f5ffe`. |
| 2026-11-07 | `/tmp/amps-t2-all-normal` | Historical FAIL: 96 PASS, 2 FAIL, 5 SKIP | First 103-case normal campaign. `PHYS-3D-COLLISION-RELAXATION-001` and `PHYS-3D-COLLISION-RATE-001` exposed a test-fixture defect: restoring the `-1` uninitialized sentinel conflicted with the new finite-positive production setter contract. |
| 2026-11-07 | `/tmp/amps-t2-weight-time-collision-recheck` | PASS 7/7 | Corrective focused run: the two populated-cell collision cases and all five T2 weight/time cases pass together after raw saved lifecycle-state restoration; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `/tmp/amps-t2-all-qualified-normal` | 98 PASS, 5 SKIP, 0 FAIL/ERROR | Qualified normal-order receipt for the live 103-case registry; all 36 native cases pass and the five optional providers remain typed SKIPs. |
| 2026-11-07 | `/tmp/amps-t2-native-qualified-reverse` | PASS 36/36 | Complete native owner in reverse case order after a clean rebuild; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `/tmp/amps-t2-native-qualified-shuffle-1729` | PASS 36/36 | Complete native owner in deterministic shuffled case order after a clean rebuild; executable SHA-256 `6329674e470d9cadda436aca2d410f1ef5ecad27f4918108d56bd3078a1f4925`. |
| 2026-11-07 | `python3 src/test/validate_traceability.py` | Historical T2 PASS | Then-current T2 structural classification: 24 areas, 209 production translation units, 181 public contracts, 15 configuration contracts, 16 suites, and 152 registered cases. No behavioral case execution. |
| 2026-11-07 | `python3 src/test/run_tests.py --campaign smoke --jobs 4 --output-dir /tmp/amps-t1-smoke-qualified` | PASS 51/51, 8/8 suites | Corrected campaign-selection oracle and complete dependency-light T1 smoke gate. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order normal --jobs 4 --output-dir /tmp/amps-t1-all-qualified-normal` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete 98-case normal order; native 31/31 PASS; provider SKIPs remain typed and explicit. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order reverse --jobs 4 --output-dir /tmp/amps-t1-all-qualified-reverse` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete reverse suite/case order with a mandatory clean native rebuild. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729 --jobs 4 --output-dir /tmp/amps-t1-all-qualified-shuffle-1729` | 93 PASS, 5 SKIP, 0 FAIL/ERROR | Complete deterministic shuffled order with a mandatory clean native rebuild. All three final native executables have SHA-256 `17b98e5ecd22915fc181cce51b8d71bd6abef517802a456ff2efba659c581919`. |
| 2026-11-07 | `/tmp/amps-t1-rng-hybrid-final`, `/tmp/amps-t1-mpi-integrated-final`, `/tmp/amps-t1-numerics-integrated-final`, `/tmp/amps-t1-units-integrated-final`, `/tmp/amps-t1-species-radiation-final-2` | PASS: 7 + 5 + 10 + 5 + 4 focused cases | Focused post-audit evidence for RNG, MPI channel/rank streams, arrays/numerics, units, and radiation tables. The complete campaigns above are the qualification authority. |
| 2026-11-07 | `/tmp/amps-t1-smoke-final` | Historical FAIL: 50 PASS, 1 FAIL | Preserved first live smoke receipt that exposed the stale `HARNESS-SELECTION-001` oracle; resolved by `/tmp/amps-t1-smoke-qualified`. |
| 2026-11-07 | `/tmp/amps-t1-interface-global` | 6 PASS, 5 SKIP | Executable C/Fortran/Geopack contracts pass; five external provider capabilities remain explicitly unavailable. |
| 2026-11-07 | `python3 src/test/validate_traceability.py --output /tmp/amps-core-traceability.json` | PASS | Structural classification: 24 areas, 209 production translation units, 181 public contracts, 12 configuration contracts, 4 suites, 45 registered cases. No behavioral case execution. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order normal --output-dir /tmp/amps-t0-all-final-normal` | PASS 45/45, 4/4 suites | Final-tree normal order; JSON/JUnit reconcile; 22 native executable cases each ran in a fresh exact-case MPI process. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order reverse --output-dir /tmp/amps-t0-all-final-reverse` | PASS 45/45, 4/4 suites | Final-tree reverse suite/case order; exact reverse receipt and identical native executable hash. |
| 2026-11-07 | `python3 src/test/run_tests.py --all --order shuffle --shuffle-seed 1729 --output-dir /tmp/amps-t0-all-final-complete-shuffle-1729` | PASS 45/45, 4/4 suites | Final-tree deterministic shuffled suite/case order; includes real timeout/SIGINT cleanup, terminal/JSON/JUnit reconciliation, distinct campaign assertions, effective generated-mode checks, and an explicit retained native preflight transcript. |
| 2026-11-07 | `/tmp/amps-t0-harness-final-complete` | PASS 11/11, 1/1 suite | Final harness-only receipt after campaign-set, interrupt, traceability, and JSON/JUnit reconciliation hardening. |
| 2026-11-07 | `/tmp/amps-t0-individual-<CASE-ID>` for all 21 dependency-light IDs | PASS 21/21 | Each dependency-light case selected and executed alone through the global `--case` interface. Native executable cases use exact one-case MPI subprocesses in every all-order campaign. |
| 2026-11-07 | `/tmp/amps-t0-individual-CONFIG-GENERATOR-001`, `/tmp/amps-t0-individual-CONFIG-INVALID-DECK-001` | PASS 1/1 each | Both adapter-owned configuration cases pass as exact one-case global selections after independent required clean rebuilds. Together with the 21 dependency-light probes and 22 one-case native MPI processes, every registered case has individual-selection evidence. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-normal` | PASS 21/21, 3/3 suites | Final dependency-light registry in normal order, including exact campaign-set assertions. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-reverse` | PASS 21/21, 3/3 suites | Same dependency-light registry in reverse suite/case order. |
| 2026-11-07 | `/tmp/amps-t0-smoke-final-shuffle-1729` | PASS 21/21, 3/3 suites | Same dependency-light registry with deterministic shuffle seed 1729. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b` | PASS 24/24, 1/1 suite | Combined native normal order from one clean configured build; includes both Config.pl/generated-deck cases, three particle-buffer cases, and 19 physics cases. |
| 2026-11-07 | `/tmp/amps-core-native-t0-b-reverse` | PASS 24/24, 1/1 suite | Combined native reverse order from one clean configured build; JSON case selection/results reconcile exactly. |
| 2026-11-07 | `test_output/core-tests/run-20261107T072523Z-3412060/` | PASS 38/38, 5/5 suites | Last complete pre-consolidation execution; historical registry only. |
| 2026-11-07 | `/tmp/amps-core-report-final` | PASS 16/16, 3/3 suites | Last retained dependency-light receipt before the 11-case harness expansion. |
| 2026-11-07 | `/tmp/amps-core-physics-audit-final` | PASS 19/19 | Last retained native-physics family receipt before suite consolidation. |
| 2026-11-07 | `/tmp/amps-core-pbuffer-expanded` | PASS 3/3 | Last retained particle-buffer family receipt before suite consolidation. |

The following records are retained as historical bring-up/progression. Their
different case counts are expected for their then-current registries; none
supersedes the consolidated final-tree 45-case evidence above.

## Historical validation evidence

- `/tmp/amps-t0-final-report-reconciliation` correctly failed
  `HARNESS-REPORT-001` while adding direct JUnit coverage because the synthetic
  self-test receipt omitted the real writer's required `seconds` field. The
  fixture was corrected without changing the writer; the immediate rerun at
  `/tmp/amps-t0-final-report-reconciliation-2` and the final 11-case harness
  receipt both pass.

- `python3 src/test/run_tests.py --all` — PASS on 2026-11-07 after both
  required clean configured native builds: 38/38 individual tests and 5/5
  suites passed with 0 FAIL/SKIP/ERROR.  The command printed the final 38-row
  individual result table and the two summaries.  JSON/JUnit receipts and all
  suite transcripts are retained in
  `test_output/core-tests/run-20261107T072523Z-3412060/`; the JSON report
  independently records `individual_counts={TOTAL:38,PASS:38,FAIL:0,SKIP:0,ERROR:0}`
  and suite `counts={TOTAL:5,PASS:5,FAIL:0,SKIP:0,ERROR:0}`.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-report-final` — PASS after strengthening the formatter
  regression: `HARNESS-REPORT-001` directly captures and checks the heading,
  per-case states, and typed suite-level fallback; all 16 dependency-light
  cases were printed and recorded, with 3/3 suites and 0 FAIL/SKIP/ERROR.

- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-smoke` — PASS: the three nonexclusive smoke suites
  executed concurrently, all 15 individual cases passed, suite-qualified case
  starts were streamed, and aggregate results remained in manifest order.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-4` — PASS after the required clean configured
  build and one-rank MPI run: the exact 19-case registry passed with 0
  FAIL/SKIP/ERROR in 54.629 seconds. JSON/JUnit reports and the full native
  transcript are retained in that directory. This was the expanded
  native-physics evidence for that registry state.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-expanded` — PASS after the required
  clean rebuild with RT and collision-frequency storage enabled: exact 3-case
  registry, 0 FAIL/SKIP/ERROR in 37.360 seconds.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final` — PASS after interruption/scratch/fail-fast
  hardening: 3 suites and 15 cases, 0 FAIL/SKIP/ERROR in 3.007 seconds.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-expanded-final` — PASS after the independent audit
  fixes and required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 54.818 seconds. This supersedes
  `/tmp/amps-core-physics-expanded-4` as the then-current native evidence.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-concurrent-final-2` — PASS after the final scheduler timeout,
  interrupt-cleanup, fail-fast, and scratch-identity changes: 3 suites, 15
  cases, 0 FAIL/SKIP/ERROR in 3.015 seconds. This supersedes the earlier
  concurrent smoke receipts.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-audit-final` — PASS after all final independent-audit
  fixes and the required clean rebuild: exact 19-case registry, 0
  FAIL/SKIP/ERROR in 55.407 seconds. This was the final native-family evidence
  before consolidation.
- `python3 src/test/run_tests.py --tier smoke --jobs 3 --output-dir
  /tmp/amps-core-smoke-audit-final` — PASS: final scheduler/protocol tree, 3
  suites and 15 cases, 0 FAIL/SKIP/ERROR in 2.987 seconds. This was a prior
  dependency-light receipt.
- `mpiexec -n 1 ./amps --filter PIC-PBUFFER-` — PASS on the executable from
  the final clean native build: 3/3 particle-buffer cases after hybrid-thread
  fixture restoration was added.

- `python3 src/test/run_tests.py --tier smoke --output-dir /tmp/amps-core-test-final-2`
  — PASS: 3 suites, 14 individual cases, 0 FAIL/SKIP/ERROR. JSON and JUnit
  receipts are in that output directory.
- After the detailed code/documentation audit,
  `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-test-final-3` — PASS: the same 3 suites and 14 individual
  cases, 0 FAIL/SKIP/ERROR. This was the final dependency-light receipt for
  that earlier registry state.
- `AMPS_TEST_OUTPUT_ROOT=/tmp/amps-core-top-level make core_test_nightly`
  — PASS: validated the top-level Make/nightly entry point with the same 14
  cases.
- `SCRATCH=/tmp python3 src/test/run_tests.py --config
  src/test/config/machine.example.json --test CORE-HARNESS-001 --output-dir
  /tmp/amps-core-machine-profile-2` — PASS and verified that transient work is
  routed through the configured scratch root.
- `python3 src/test/run_tests.py --coverage --require-complete-coverage` —
  expected exit 1: 2 complete and 22 incomplete inventory components.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-native-final-3` — PASS after a clean configured
  native build and one-rank MPI run: 1 suite, 3 cases, 0 FAIL/SKIP/ERROR. The
  JSON/JUnit reports, full build log, and particle-buffer transcript are in
  `/tmp/amps-core-native-final-3`.
- `python3 src/test/run_tests.py --test CORE-PIC-PHYSICS-001 --output-dir
  /tmp/amps-core-physics-5` — PASS after a clean configured build and one-rank
  MPI run: 1 suite, 11 cases, 0 FAIL/SKIP/ERROR. JSON/JUnit reports and the
  complete native transcript are retained in that directory.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-pbuffer-physics-final` — PASS after a second
  clean configured build: the pre-existing exact 3-case registry remains green
  with the expanded species and optional particle-data layout.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-smoke-physics-final` — PASS: 3 dependency-light suites, 14
  cases, 0 FAIL/SKIP/ERROR after registering the physics suite.
- `python3 src/test/run_tests.py --all --list` and `python3
  src/test/run_tests.py --coverage` — PASS registry/inventory audit: 5 suites,
  24 components, 2 covered, and 22 honestly partial/planned.
- `python3 src/test/run_tests.py --tier smoke --output-dir
  /tmp/amps-core-case-progress-final` — PASS: all 14 cases were announced by
  stable ID and description before execution, with unchanged JSON/JUnit
  results.
- `python3 src/test/run_tests.py --test CORE-PIC-PARTICLE-BUFFER-001
  --output-dir /tmp/amps-core-case-progress-native` — PASS after the required
  clean native build: all 3 MPI-launched cases were announced live and the
  exact registry remained green.
- `git diff --check` — PASS after the implementation and documentation audit.
- `PYTHONPYCACHEPREFIX=/tmp/amps-core-pycache python3 -m compileall -q
  src/test` — PASS for every Python runner/framework module without writing
  bytecode into the source tree.

## Native bring-up record

The first clean native attempt (`/tmp/amps-core-native-final`) failed while
compiling `meshAMR`: the new input deck had not invoked `exosphere.pl`, which
is the established generator for AMPS' default internal-boundary user-data
classes. The fix adds a target-free, SPICE-free parser block. It supplies type
definitions only; exosphere physics and source modules remain excluded.

The second clean attempt (`/tmp/amps-core-native-final-2`) compiled and linked
successfully but reported no cases. AMPS stores application objects in the
static archive `main.a`, and the linker did not load `particle_buffer.o` when
it contained only static registration constructors. The fixture now exports
`LinkParticleBufferTests()`, called by `main.cpp`, to make archive membership
deterministic. The third clean attempt passed the exact three-case registry.
Both failed attempts remain useful diagnostic evidence and were not described
as test passes.

## Native physics bring-up record

- `/tmp/amps-core-physics-1` failed during input generation because the
  collision table used a single backslash instead of AMPS' double-backslash
  continuation token.
- `/tmp/amps-core-physics-2` reached the IDF parser and rejected unsupported
  `RotZnum all` shorthand. The unused rotation-collision and pair-index values
  were removed because every relaxation switch is off in this kernel fixture.
- `/tmp/amps-core-physics-3` built and passed all nine non-mover cases. The
  mover fixture incorrectly selected an unallocated AMR internal/root node;
  the second case inherited it and triggered the expected signal handler.
- `/tmp/amps-core-physics-4` resolved the owner leaf and produced correct
  ballistic states, but inspected the starting leaf's temporary list for an
  endpoint exactly on a leaf face. The fixture now avoids an ambiguous face
  endpoint and resolves the destination owner explicitly.
- `/tmp/amps-core-physics-5` passed the exact 11-case registry. None of the
  failed bring-up runs is counted as validation evidence.
- `/tmp/amps-core-physics-expanded-1` executed all 19 expanded cases and
  reported 16 PASS / 3 FAIL. The failures were test-fixture/oracle defects:
  the neutral-only build intentionally omits unavailable ionic daughters, the
  RT oracle initially used one molecule's rotational exponent rather than the
  pair's six rotational degrees, and the Boris free-flight trajectory crossed
  the newly registered sphere. All three definitions were corrected without
  relaxing a production tolerance.
- `/tmp/amps-core-physics-expanded-2` was manually interrupted because the
  240,000-draw branching check also generated full two/three-body kinematics
  on every draw and was unsuitable for routine execution. Its orphaned adapter
  and MPI processes were explicitly cleaned up. The final test separates the
  high-statistics production channel selector from a smaller product-contract
  sample; product kinematics retain their dedicated conservation case.
- `/tmp/amps-core-physics-expanded-3` produced a typed `ERROR` because the
  interrupted adapter still held the repository-specific native lock. No test
  result from that invocation is counted as scientific validation.
- `/tmp/amps-core-physics-expanded-4` passed the final exact 19-case registry
  in 54.629 seconds. It supersedes the three diagnostic expanded bring-up runs;
  those failures/interruption remain recorded and are not described as passes.
