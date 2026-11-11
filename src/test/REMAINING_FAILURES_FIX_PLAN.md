# Plan to repair the remaining validation failures

> **2026-11-10 current checkpoint:** CMBGU05 passes unchanged physical bounds
> and the explicitly approved stronger all-level gate. Twenty-nine new bounded
> cross-ray/material/contact cases and derivative/shell controls pass. The
> user retained the 3-D ambient, selected a threaded contact and 1200-s
> corotation taper, and retained an analytical/semi-empirical closure WITHOUT
> an MHD runtime. Production lateral transport/qualification remains open.
> Full 1,067-case core revalidation is complete: 825 PASS / 1 FAIL /
> 240 SKIP / 1 ERROR; every old 792 PASS and all 32 new cases pass.
> CSWC0628's original timeout and T12's 26 incomplete components remain visible.
> See STATUS.md and the authoritative full-CME ledger for receipts and limits.


Prepared: 2026-11-10. Implementation is in progress. R1's full contact
executable passes with the original coronal library; its independent reference
and rejection controls pass in ordinary, debug and sanitizer profiles. The
user explicitly approved replacing the obsolete expectation of a failed
distributed contact level with the stronger requirement that every level
satisfy the original physical bounds. The harness passes 13/13, including new
strict-coverage rejection controls. R2 and complete coverage remain open.

The user subsequently explicitly authorized resuming full-CME implementation
and qualification. This supersedes the deferral mentioned in the original plan
below, with work resuming at G1/BG3D-4 and the reopened contact dependency.
G2 and later particle stages remain outside this authorization. Qualification
must still follow every original gate. Current receipts are under
`test_output/remaining-failures-repair-20261110/`; no complete repaired core
campaign is claimed yet.

The objective is reliable validation, with retained evidence that an incorrect
implementation fails. A green report obtained by reducing the problem,
changing its expected answer, or hiding an unfinished qualification is not
completion of this plan.

## 1. Scope and frozen starting evidence

This plan covers the two remaining FAIL cases and the one ERROR in the latest
complete core run. The exosphere and dust configuration failures are already
repaired and must stay repaired.

| Exact case | Current result | Observed failure | Repair track |
| --- | --- | --- | --- |
| `CSWC-FULL-CMBGU05` | FAIL | The combined pointwise contact check rejects its finite-difference sequence | A: independent contact derivative and complete sheath regression |
| `CSWC-FULL-CSWC0628-MASS-TIME` | ERROR | The refinement executable exceeds its 180-second routine execution limit | B: measure and reduce cost, then grade the complete unchanged refinement matrix |
| `T12-COVERAGE-STRICT-001` | FAIL | The real inventory contains 2 complete and 26 incomplete components | C: implement and qualify the missing coverage; preserve the strict audit |

Authoritative starting report:

- [Complete JSON](../../test_output/native-config-isolation-20261110/core-final/amps-core-tests.json)
  and [JUnit](../../test_output/native-config-isolation-20261110/core-final/amps-core-tests.xml):
  **1,035 cases: 792 PASS, 2 FAIL, 240 SKIP, 1 ERROR**;
  50 suites: 43 PASS, 1 FAIL, 5 SKIP, 1 ERROR; 1,744.066 seconds.
- [Full-volume model transcript](../../test_output/native-config-isolation-20261110/core-final/artifacts/CORE-MODEL-SEP-CORONA-SWCME-FULL-001/full-volume-local.txt).
- [T12 transcript, including the actual strict inventory](../../test_output/native-config-isolation-20261110/core-final/artifacts/CORE-T12-PORTABILITY-001/t12-portability-matrix.txt).
- [Original baseline comparison](../../test_output/native-config-isolation-20261110/baseline-comparison.json):
  all 785 original passes preserved. The new repair baseline also protects all
  792 currently passing individual identities.

The current registry has 55 suites and 1,142 cases across core and application
streams. Counts alone are insufficient: compare stable identities, statuses,
selected cases, skipped cases, and evidence provenance.

The Hartley2 four-rank initialization timeout is a **separate application
runtime issue**, outside these three core results. Its existing
[geometry investigation and optimization plan](GEOMETRY_PERFORMANCE_INVESTIGATION.md)
remain applicable. A core repair does not close that issue.

## 2. Requirements that cannot be relaxed

1. Preserve each existing case identity, registration, physical fixture,
   acceptance threshold, and truthful global outcome. Add new cases where
   needed; do not replace an existing required case with a smaller smoke test.
2. Do not increase contact-flux bounds, loosen convergence ratios, coarsen the
   refinement grids, shorten the physical integration, change the shock branch,
   drop a ray, or tune the physical controls against the failing answer.
3. Do not change FAIL/ERROR to SKIP, expected failure, warning, or PASS to obtain
   a green aggregate. A known defective state may be used in an additional
   negative control; it cannot stand in for the required positive validation.
4. Keep the refinement's **180-second routine ceiling** and the owner's current
   3,600-second outer timeout. A longer explicitly labeled diagnostic run may
   establish what happens after 180 seconds; it cannot satisfy the routine gate.
5. Keep strict coverage over the complete current inventory, including registered
   applications. A core-only execution selector does not narrow that inventory.
   Do not add exemptions, remove components, or mark partial work covered to
   repair T12. New discovered components remain visible.
6. Preserve numerical error and physical model error as separate quantities.
   A numerical improvement cannot qualify a physically incompatible closure.
7. Preserve all existing native build preflight, stream/build serialization,
   machine-profile selection, configuration receipts, and no-SPICE fixture
   isolation. Every native rebuild must confirm the AMPS root and idle state,
   remove only root `build/`, regenerate through Config.pl, and use top-level
   `make amps` with `MAKEFLAGS=-j16`.
8. Git remains read-only. Leave source and documentation changes unstaged and
   uncommitted. Preserve site settings, user changes, reference repositories,
   existing inputs, and all original failure evidence.

## 3. Execution order and decision points

| Stage | Work and deliverable | Exit condition |
| --- | --- | --- |
| R0 | Freeze case-level baseline, inputs, binaries, tools, controls, limits, and failure predicates | Reproductions distinguish a physical rejection, an unstable numerical measurement, a timeout, and incomplete implementation |
| R1 | Track A: contact diagnostics and independent derivative validation | Complete `CMBGU05` passes every existing assertion, with added independent and negative evidence; its qualification scope remains explicit |
| R2 | Track B: profile the unchanged full refinement workload and optimize measured costs | Every prescribed level and comparison completes within 180 seconds and passes unchanged gates |
| R3 | Track C: convert the 26 incomplete rows and their subordinate gaps into individually owned work packages | Every gap has an oracle, executable cases, configuration/rank coverage, evidence requirements, and a closure decision |
| R4 | Implement and qualify those work packages in dependency order | Components advance only on retained behavioral evidence; strict coverage remains FAIL until the last required component is complete |
| R5 | Run focused, full-core, application, order-isolation, and qualification regressions; compare identities | No current PASS regresses, no required identity disappears, and each claimed repaired case passes with full evidence |

R0 is partly informed by the inspected receipts below, but fresh execution and
profiling are implementation work, not evidence supplied by this plan.
Finish R1 and R2 before claiming the two local model problems repaired.
R3/R4 are a larger validation program and cannot honestly be closed by a small
runner patch. Their progress must remain visible even after R1/R2 pass.

The active repository program defers full-volume BG3D-4 production work. This
plan identifies the work required for that closure without silently resuming
the entire campaign. Numerical test repair and bounded shared-library
validation must not be reported as full-CME qualification. Decisions about
unimplemented physics or unavailable external references require clarification
before dependent implementation; no physical closure is invented here.

## 4. Track A: contact-leakage validation (`CSWC-FULL-CMBGU05`)

### 4.1 What is established, and what is not

Owned implementation and tests:

- [test_bg3d4.cpp](../models/sep_corona_swcme/test/test_bg3d4.cpp):
  `Fourth`, `TestMaterialContactAndStartup`, reference fixtures, material
  inventory, differential identities, distributed residuals, and spatial inverse.
- [sheath_model.cpp](../models/sep_corona_swcme/src/sheath_model.cpp):
  `ShockFedSheathModel::Evaluate` and `EvaluateContact`.
- [Frozen regional controls](../models/sep_corona_swcme/examples/bg3d1/ejecta.asset).
- [Authoritative correction contract](../models/sep_corona_swcme/docs/BG3D4_CODEX_REVIEW_AND_CORRECTIONS.md)
  and [selected production piston closure](../models/sep_corona_swcme/docs/BG3D4_PISTON_CLOSURE.md).

The retained pointwise diagnostic is:

```text
flux_bound = 1.001e-17 kg m^-2 s^-1
h_s        = 0.1,        0.05,        0.025
fd_fluxes  = 2.58964e-20, 3.45297e-20, 8.14684e-20 kg m^-2 s^-1
```

The frozen bound is

```text
F_limit = F_absolute + epsilon_relative * F_reference
        = 1e-20 + 1e-10 * 1e-7
        = 1.001e-17 kg m^-2 s^-1.
```

The printed finite-difference fluxes are below that bound. The immediate
nonpassing clause is `leakage[1] < leakage[0]`, within a combined assertion
whose error text describes excessive leakage. The API-reported flux is zero
by construction because its boundary velocity is copied from the primitive
velocity. That zero is not an independent impermeability measurement.

The current fourth-order derivative evaluates differences of heliocentric
positions of order `1e9 m`. Cancellation is a plausible explanation for the
nondecreasing small errors; it has **not** been established by a conditioning
study. Neither the small local flux nor a future stable derivative proves
global contact compatibility.

There are other deliberate limitations in this executable: it detects a
fixed-fraction/material-contact mismatch and requires the old diagnostic map
to remain explicitly unqualified in distributed contact diagnostics. A PASS
of this bounded software regression must retain that meaning. Do not change
its message to imply that the replacement production closure is qualified.

### 4.2 Diagnosis before choosing a fix

1. Preserve the original failing run. Reproduce the exact supported chart,
   material labels, epoch, positions, density, normal, primitive velocity,
   reported boundary velocity, and all stencil samples at full precision.
2. Report each predicate separately: reported flux, independent flux at every
   level, physical bound, and convergence inequality. Keep the final case
   failing if any required predicate fails. This corrects diagnosis without
   changing acceptance.
3. Perform a diagnostic step-size study extending both sides of the frozen
   sequence, only within declared temporal support. Record coordinate scale,
   cancellation estimate, integration/discretization error, and floating-point
   precision. Do not select the smallest measured error after seeing results.
4. Compare independent calculations: the original four-point formula,
   mathematically equivalent paired differences, and higher-precision stencil
   accumulation. Casting already rounded double positions to higher precision
   does not recover position information; quantify that limit explicitly.
5. Verify the same derivative method on independent moving planar and curved
   material surfaces with known position, normal, normal speed, density, and
   traction. Include a large fixed HCI translation and ordinary SI scales.
   Test fourth-order convergence where truncation error is resolvable.
6. If needed, independently evaluate the trajectory equations at higher
   precision or differentiate a separate reference implementation. Do not
   obtain expected derivatives from the returned production velocity, and do
   not add a new package or platform-dependent precision assumption silently.

Deliverable: a table distinguishing truncation, position-evaluation error,
roundoff, and physical slip, plus a justified numerical representation.

### 4.3 Permitted implementation and stop conditions

Prefer the smallest demonstrated correction: stable evaluation of the same
fourth-order stencil, or a proven correction to production trajectory
evaluation if that is the defect. Retain the frozen `0.1/0.05/0.025 s` checks,
their physical flux limit, and the convergence predicate. Do not simply remove
the monotonic clause, add an arbitrary floor, or accept increasing errors
because all three values are small.

If no trustworthy measurement at those controls can satisfy the existing
predicate, report that limitation before changing the definition of the test.
A mathematically justified redesign would need independent error certification
and the same physical discrimination; this plan does not preauthorize a
replacement gate. Until then the case remains FAIL. A diagnosed artifact is
not permission to mark the test passed.

If independent evidence finds real contact slip, fix the relevant trajectory,
normal, material map, or chosen physical closure. Do not force a zero flux by
assigning velocities equal, tuning the contact fraction, suppressing crossings,
changing retained inventory, or shrinking the supported domain. If the remedy
requires the deferred production closure, document the dependency and obtain
the necessary physical decisions before implementation.

After the first assertion is repaired, execute **all** later checks that its
early exit currently hides. Keep absolute area-time leakage integrals, typed
unsupported charts, curved-volume measures, mass/admission ledgers, no material
loss at subfast epochs, Cauchy induction, thermal evolution, smooth residuals,
map/inverse regularity, transactional failures, and traction checks. Another
failure becomes the next defect to investigate; do not stop after the first
line passes.

### 4.4 Prove the test still rejects a defect

Add registered tests that feed independent manufactured samples through the
actual derivative/flux grader. Prescribe normal slip safely above the frozen
limit and require rejection; also test tangential reparameterization, translated
coordinates, a bad normal, invalid samples, and unsupported stencil endpoints.
Ensure the defect reaches the measured positions or primitive state, rather
than only editing a printed flux. A separate deliberately incorrect material
history must fail the integrated absolute-leakage check without cancellation
between positive and negative flux.

These controls strengthen the existing test. They do not replace its positive
trajectory, full map, or physics assertions. Retain ordinary and optimized
compiler runs and ASan/UBSan evidence for modified shared code. Preserve all
existing planar, spherical, gas, MHD, contact, ambient, and architecture tests.

Exit condition: full original executable success, exact global case PASS,
independent derivative evidence, rejection controls PASS, and no unsupported
promotion of BG3D-4/G1. Physical production qualification remains a distinct
open gate wherever the actual closure still fails its contract.

## 5. Track B: complete the full refinement test within its budget

Sources:

- [test_bg3d4_piston_refinement.cpp](../models/sep_corona_swcme/test/test_bg3d4_piston_refinement.cpp).
- [piston_solver.cpp](../models/sep_corona_swcme/src/piston_solver.cpp),
  [piston_sheath.cpp](../models/sep_corona_swcme/src/piston_sheath.cpp), and
  the associated contact/ambient implementations.
- [Global full-volume adapter](tests/models/cme/sep_corona_swcme_full/run_test.py).

### 5.1 Freeze the complete numerical experiment

| Quantity | Required unchanged experiment |
| --- | --- |
| Physical epoch | `3600 s` for every graded level |
| Mass refinement | 40, 80, 160 initial cells, each at CFL `0.2` |
| Time refinement | 40 cells at CFL `0.4`, `0.2`, `0.1`; the same existing CFL-0.2 run is reused |
| Buffer controls | Existing resolution-scaled minimum/append cell counts; same physical maintenance lattice and append policy |
| Compared quantities | Shock radius, compression, speed, fast Mach, contact total pressure |
| Branch requirement | Exactly one common detected shock zone on all three mass levels |
| Mass convergence | Fine/coarse change comparison uses the existing `1.1` factor |
| Fine mass limits | Radius `<3%`; compression, speed, Mach `<8%`; contact pressure `<10%` |
| Time convergence | Existing `1.25` comparison factor; all fine changes `<1%` |
| Identity | Every distinct refinement control set changes the immutable event fingerprint |
| Routine execution | Entire executable, all levels and comparisons, finishes within `180 s` |

This is five distinct production advances, not one diagnostic level. The
executable's single-level mode emits no suite PASS and cannot qualify it.
The earlier `1800 s` fixture was rejected for differing shock branches; do
not restore it to obtain a faster answer. Do not assume that removing the
timeout will yield a numerical PASS: the current report contains no complete
matrix grading.

### 5.2 Measure cost before optimizing

1. Reproduce with the configured compiler/flags and no competing test/build.
   Record wall/CPU time, memory, ray and cell counts, accepted steps, minimum
   timestep, append events, source/ambient/contact evaluations, and time spent
   in each level. Retain partial output on timeout as diagnostic evidence.
2. Complete one explicitly budgeted diagnostic run if necessary to learn the
   full runtime and eventual numerical result. Label it nonqualifying for the
   180-second gate. If profiling cannot complete within available resources,
   record that limitation rather than inventing a budget or result.
3. Profile actual hot paths. Code inspection suggests repeated RK-stage
   allocation/copy, ambient/source interpolation, contact history evaluation,
   and per-maintenance state copies as candidates; these are hypotheses, not
   measured root causes.
4. Check timestep histories for unexpected collapse, excess appending, or
   wrong resolution scaling. Fixing a solver/domain defect takes precedence
   over making the same defective calculation faster.

### 5.3 Optimize implementation, preserve equations and ownership

Optimize only measured costs. Candidates include reusable per-instance RK
scratch, immutable projection lookup data, and reuse of exactly identical
boundary evaluations at the same immutable event/time. Every cache key must
include all relevant event/control identity and have bounded lifetime.

Preserve RK4 stage equations, CFL and thermal-step restrictions, artificial
viscosity, source work, magnetic invariants, append schedules, all rays, and
all-or-nothing commit behavior. Scratch reuse must not alias a committed state
or survive a rejected private advance. Avoid shared mutable caches and do not
enable fast-math or change precision as a shortcut.

Parallel independent-level execution is not the first remedy: it can hide
serial cost and oversubscribe the global `--jobs` budget. Consider it only
with explicit resource accounting and identical complete grading. It cannot
drop levels or introduce an unbounded internal worker pool.

After each change compare unoptimized and optimized states/ledgers on frozen
fixtures; bitwise equality is preferred when arithmetic order is unchanged.
For intentional arithmetic changes, freeze an independent error budget before
grading, within every existing acceptance limit. Preserve rollback tests,
well balance, planar/curved gas/MHD references, signed magnetic flux, energy
work, and shock detection.

Add deterministic work/allocation counters for a representative hot-path
fixture so a faster host cannot be the sole explanation of a performance
improvement. Counter limits must come from the justified algorithmic work,
not an arbitrary timing observed after tuning.

### 5.4 Reporting and exit condition

The global adapter must continue requiring the unique final PASS marker and
zero exit from the complete executable. Timeout, crash, missing final levels,
and malformed receipts remain ERROR; failed numerical criteria remain FAIL.
Retain the per-level metrics, exact binary/input hashes, compiler settings,
runtime measurements, and full comparison output.

Exit condition: repeated complete runs on the declared Passion profile within
the existing 180-second limit, all frozen checks passing, and a fresh global
`CSWC-FULL-CSWC0628-MASS-TIME` PASS. Repeatedly exceeding the budget leaves this
track open; an increased timeout is not an accepted repair under this plan.

## 6. Track C: earn strict coverage with actual tests

### 6.1 Why this is larger than a failing test implementation

[The T12 adapter](tests/portability/matrix/run_test.py) executes the real
`--coverage --require-complete-coverage` audit and imports its nonzero exit.
[The inventory](component_inventory.json) presently has **28 components,
2 covered and 26 partial**. It correctly fails. The current complete entries
are `src/general/input` and `src/pic/units`.

The inventory audit checks declarations and discovered ownership; it does not
execute physics merely because an entry says `covered`. Therefore changing
labels would bypass the intended validation. Preserve
[source traceability](source_traceability.json),
[configuration contracts](configuration_contracts.json),
[public algorithm inventory](t1_public_algorithm_inventory.json), local
qualification records, and fresh executed receipts as the closure evidence.
Registering or compiling a test alone is insufficient.

The following is the current closure backlog, extracted from the live inventory.
During R3 expand each row into a finite contract checklist with stable owning
cases, valid/invalid domains, independent oracles, profiles, and explicit
capability requirements. These summaries do not replace any subordinate gap.

| Incomplete component | Required work before coverage can advance | Independent evidence to require |
| --- | --- | --- |
| `src/general` | Logger lifecycle, CRC/debug output, synchronization/CUDA helpers, remaining grid/distribution containers, wider sanitizer and performance profiles | Independent checksums/data structures; lifecycle and fault injection; thread/race and device evidence where supported |
| `src/interface` | Remaining empirical fields and ARMS/BATSRUS/CCMC/Kameleon/SWMF production wrappers | Versioned external reference states, ABI/error/domain tests, configured real providers; synthetic echo alone is insufficient |
| `src/meshAMR` | Nonmanifold and large surfaces, non-Cartesian paths, deeper imbalance, topology-file restart | Analytic geometry, topology invariants, exact restart, owner/ghost and conservation checks |
| `src/models/charge_exchange` | Define and qualify the missing event/product/branch/velocity/weight interfaces where required | Physical network decisions first; independent rates, particle/charge/momentum/energy balances and product statistics |
| `src/models/chemical_solver` | Runtime network assembly, cell coupling, elemental/charge metadata, multi-network application integration | Analytic network solutions and separate stiff references, conservation, invalid networks, actual cell integration |
| `src/models/cometary_physics` | Native transport/collision coupling beyond the standalone contracts | Fresh Hartley application fixtures and independent collision/energy/transport references |
| `src/models/resonance_fluorescence` | Observation and transfer validation beyond published local formulas/tables | Versioned measurement/reference provenance, optically thin limits, declared transfer regime and uncertainty |
| `src/models/dust` | MPI injection/migration/reduction, restart, full stochastic convergence | Size-group physical ledgers, exact restart, rank invariance, preregistered sampling error; preserve dust-owned weight policy |
| `src/models/electron_impact` | Remaining parent/product distributions, runtime electrons, insertion/weights, thread-safe return storage | Independent branching/kinematics and elemental/charge/energy closure, distribution tests and race evidence |
| `src/models/exosphere` | MPI sources/reduction, SPICE/orbit inputs, target maps/user sources, external-domain and nonspherical modes | Independent source normalization, geometry/ephemeris references, physical inventory and rank invariance |
| `src/models/photolytic_reactions` | Quiet-Sun tables, runtime networks/cycles, MPI/thread-safe daughters and multi-parent callbacks | Published rate provenance, analytic survival/branching, independent daughter conservation and exact-once processing |
| `src/models/sep_common` | Qualify all 52 expensive pitch-angle Monte Carlo cases currently skipped in routine execution | Complete fresh campaigns, independent transport limits and registered statistical budgets; partial histories cannot qualify |
| `src/models/sep_corona_swcme` | Tracks A/B, reduced native cases, and the declared full-volume qualification gaps | Independent contact/conservation/residual/refinement evidence; actual native receipts; no reduced/full scope substitution |
| `src/models/sep_coronal_cme` | Remaining campaign evidence and physical/observational qualification | Independently sourced events, complete staged evidence and native campaigns required by the owning contract |
| `src/models/sputtering` | Native insertion, real incident-flux providers, exact provenance of inherited numeric tables | Independent yield/table references, angular/speed distributions, source and momentum/energy ledgers |
| `src/models/surface` | Nonspherical production surfaces, target temperature/chemistry laws, MPI ownership/reduction | Analytic/manufactured surfaces, accommodation/adsorption/desorption balances and rank-invariant inventory |
| `src/models/swcme` | Restore six absent extended test translation units and the missing observation runner; run 127 extended and 16 observational identities | Recover authoritative files from a read-only source or clarify missing material; independent event fixtures and complete receipts |
| `src/pic` | Remaining CUDA, opaque-extension semantics, coarse/fine and surface movers, distributed sources, histogram/callback/restart/legacy format paths, Yinsi weights, unequal Zrot and VT/VV | Independent trajectories and boundary outcomes; conservation/rates/relaxation and restart continuation across supported profiles/ranks |
| `src/pic/ecsim` | Conductor boundaries, multi-block assembly/GMRES and remote halo stencils | Manufactured fields/solutions, solver residual/order, interface and remote-consumption checks |
| `src/pic/gyro` | Full orbit integration, adiabatic validity, magnetic-moment conservation, convergence and terminal handling | Analytic uniform fields, independently integrated variable-field orbits, controlled scale separation |
| `src/pic/mesh` | Collector reductions and collector-specific physical fields | Analytic integrals/fluxes and multi-rank reduction invariants |
| `src/pic/parallel` | Affinity-policy behavior and ECSIM remote-stencil use | Typed scheduler constraints and actual remote readback; a successful send alone is insufficient |
| `src/species` | Independent scientific provenance/digitization fixtures for He/Na tables | Versioned bibliographic/page-level references separate from production arrays |
| `srcEnceladus` | Broader INMS/UVIS comparisons, MPI behavior, long-run source distribution and costly E3 validation | Independent measured data and uncertainty, native source/rank ledgers, complete selected observation campaigns |
| `srcMoon` | Bring I01-I10 multi-build orchestration into required preflight; complete M1/M2 and external validation | Fresh compliant builds, independent physical references and local/global identity reconciliation |
| `src/test` | Manifest-owned standalone traceability, legacy-oracle reviews, cross-platform fault injection, sharded aggregation, cancellation/resource stress and complete isolation evidence | Adversarial subprocess/receipt fixtures, exact identity reconciliation, real cross-platform and order tests |

The table describes completion obligations already present in the inventory;
it does not authorize inventing missing physics, measurements, or provider
packages. Missing references need an explicit acquisition/provenance decision.
Missing physical APIs need a specified contract before implementation. Deferred
CME stages retain their authorization and ordering boundary. T12 stays FAIL
while any of these required gaps remains unresolved.

### 6.2 Work package and promotion rules

For every component:

1. Enumerate public behaviors and supported configuration variants against
   actual implementation and all existing gap ledgers. Keep unsupported or
   unfinished capabilities typed and visible.
2. Define independent valid and invalid fixtures, analytical/reference values,
   units, conservation laws, convergence requirements, and stochastic confidence
   budgets before executing candidate implementations. Record provenance.
3. Exercise production entry points, including native allocation, ownership,
   dispatch, state update, and output where those are part of the contract.
   A test-only reimplementation cannot prove the production path.
4. Add negative controls that demonstrate failure for representative defects:
   wrong rates/signs, lost particles/weights, bad ownership, broken conservation,
   incorrect units, wrong boundary conditions, and corrupt/truncated receipts.
5. Register cases and exact local/global mapping; document what/how/why/where;
   update the README and regenerate `TEST_CATALOG.md`. Extend source and
   configuration traceability together.
6. Retain fresh required profile/rank/toolchain results and all typed nonpasses.
   Optional-unavailable evidence does not satisfy a required qualification.
   Keep deterministic seeds, tolerances and statistical trial budgets fixed;
   do not reroll stochastic tests until one passes.
7. Review all remaining gaps, then advance the component only when its complete
   declared contract is evidenced. No automatic promotion from case counts,
   coverage percentages, build success, or the absence of a reported failure.

Implement in dependency waves: harness/provenance and configuration first;
general/species/ABI next; mesh/PIC/movers/fields/parallel next; chemistry and
source/dust/surface modules next; model and application scientific campaigns
last. Honor the existing CME stage order within those dependencies.

### 6.3 Keep the strict audit itself a meaningful regression

Add isolated harness fixtures demonstrating that strict coverage rejects an
incomplete component, a newly discovered component, missing case ownership,
stale/unknown IDs, missing required evidence, and contradictory covered/gap
declarations. Preserve the existing structural checks; extend them only to
make evidence harder to fake. A fully evidenced synthetic inventory may pass
the harness audit, but that receipt must never mark the real AMPS inventory
complete.

When all real rows are legitimately complete, execute both the strict CLI and
`T12-COVERAGE-STRICT-001` with the real profile and retain their outputs.
The existing strict gate is metadata completeness; fresh required executions
remain an additional necessary validation step, not something that changing
metadata executes automatically.

Exit condition: no required inventory/traceability/configuration/qualification
gap remains, the strict audit exits zero, the T12 case passes, and every
required behavioral campaign has a valid fresh PASS receipt. If full closure
is not achieved, publish the remaining rows and retain the real FAIL.

## 7. Regression commands and retained evidence

The physics, native, coverage, and qualification campaigns below are planned
future execution from `/home/vtenishe/Mars2/AMPS`; they are not new test results
supplied by this document. Preparing this plan includes only documentation
consistency checks. Use a unique output directory per attempt and the appropriate
reviewed machine profile on another host.

Focused cases (the current full-volume adapter still executes its complete
local registry before filtering selected result rows; these commands are not
a promise of a cheap single-binary run):

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --case CSWC-FULL-CMBGU05 --output-dir test_output/remaining-failures-repair/contact
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --case CSWC-FULL-CSWC0628-MASS-TIME --output-dir test_output/remaining-failures-repair/refinement
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --case T12-COVERAGE-STRICT-001 --output-dir test_output/remaining-failures-repair/strict-case
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --coverage --require-complete-coverage
```

After focused acceptance, validate the complete owner and core selection:

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-MODEL-SEP-CORONA-SWCME-FULL-001 --output-dir test_output/remaining-failures-repair/full-volume-owner
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --test CORE-HARNESS-001 --output-dir test_output/remaining-failures-repair/harness
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --match 'CORE-*' --streams 1 --jobs 8 --output-dir test_output/remaining-failures-repair/core-normal
```

Use focused normal/reverse/seed-1729 order runs to investigate shared-state
isolation. Once those pass, retain a complete alternate-order regression if
stateful implementation changes justify it; no need to repeat an expensive
whole campaign after documentation-only edits. Recheck the no-SPICE fixtures
from both enabled and disabled starting states when build orchestration changes.

After production core changes, rerun the required Hartley2 stage suite and
the Moon/fast-wave make regressions through compliant clean native preparation:

```sh
bash srcHartley2/tests/run_all_stage_tests.sh
make -j16 test_Moon
make -j16 test_fast-wave
```

Do not run the three commands concurrently if they can configure/build the
same checkout. The command list does not waive native preflight or authorize
real-repository Git writes; Stage-10 synthetic commits retain only their
existing explicitly authorized private-repository scope. Production geometry
changes additionally require the canonical four-rank Hartley2 run to complete
all 100 steps under the separate geometry plan.

Full qualification, when its required dependencies and authorization are ready:

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json --full-validation --streams 3 --jobs 8 --output-dir test_output/remaining-failures-repair/full-qualification
python3 -B src/test/generate_test_catalog.py --check
python3 -B src/test/validate_traceability.py
git --no-optional-locks diff --check
```

The existing `--full-validation` selection is intentionally non-narrowable;
it includes costly, distributed, portability, and application qualification
owners. Missing tools/data or required physical rejection remain visible.
It must not be treated as synonymous with a routine core-only PASS.

Retain under each attempt:

- Exact argv, working directory, machine profile and resolved tool versions.
- Source/input/configuration/binary hashes, unchanged gate values, build logs,
  and the root/idle/clean/configure/top-level-build receipt for native builds.
- Contact samples, error/flux tables, independent references and defect controls.
- Every refinement level, fingerprint, shock branch, comparison, runtime and
  profiler/counter evidence; partial timed-out output is retained but ungraded.
- Component-by-component gap closure and bibliographic/measurement provenance.
- Individual test logs, suite/stream receipts, single combined JSON/JUnit, all
  failed and skipped identities, and exact case-level baseline comparison.

Preserve all 792 current PASS identities and separately all 785 original
pre-integration passes. New cases must pass their own defined gates. A
disappearing case or unexplained change to a baseline SKIP is not a successful
comparison. Explain intended repaired-status changes individually.

Update `src/test/STATUS.md`, relevant module READMEs, this plan's per-stage
evidence, and the integration ledger after each accepted stage. If CME
implementation resumes, also maintain its authoritative execution ledger and
preserve the reduced/full qualification boundary.

## 8. Completion claims

The first achievable milestone is **the two local full-volume software cases
repaired with no gate relaxation**. It does not close strict coverage.

The second is **strict coverage earned by complete required evidence**. This
is the end of the larger AMPS validation backlog, not a relabeling operation.

Hartley2 3-D runtime completion and full-CME physical qualification remain
separately named obligations. They may be claimed only from their own actual
required results. No implementation or new PASS is claimed by this plan.
