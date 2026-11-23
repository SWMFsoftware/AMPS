# Standalone core integration: implementation and qualification status

## Stage 7 COMPLETE, 2026-11-22

Stage 7 is complete for the prescribed reduced SSE/Parker material profile.
All 150 required real-MPI job checks pass against the final rebuilt executables.
The legacy checkpoint fixtures and the numerical/threading defects exposed by
their qualification are repaired. Stages 8–9 have not begun in this batch.
Current evidence root:
`/data/vtenishe/TEST/sep-stage7-completion-20261122T182436Z`.
The original dirty checkout and Makefile are retained in `baseline.json`,
`baseline.patch`, `baseline-status.txt` and `baseline-Makefile`. The earlier
stage-7 evidence under
`/data/vtenishe/TEST/sep-integration-stages7-9-20261122T091334Z` remains
historical evidence for its own frozen source snapshot.

The physical timestep remains **300 s** in the material campaigns. The source
normalization, absorbing return boundary, prescribed turbulence and **1e-10**
SI conservation bound are unchanged. No physical reference was regenerated.

### Repairs and preserved contracts

CP05C01–05 now run the actual host on the registered moving reduced SSE/Parker
line, with explicit startup/event input and a nonempty accelerated population.
CP05F01–03 retain their stable IDs and separate-process byte-exact comparison.
MFP must observe physical scattering. The old wave and codec assertions are
retained as CP05O01–05 / CP05OF01–03 in a separate native operator executable,
which never moves particles or advances a host on frozen nonzero-U geometry.
This separates software-state checks from physical transport qualification;
no material-owner or self-consistent wave-closure guard is bypassed. MCP13
requires missing material inputs to produce native ERROR/exit 2 before
allocation, with structured JSON rather than an uncaught abort.

A material MFP continuation test exposed inconsistent particle Cartesian X
following inner segment insertion. Adding an integer prefix to a binary64 line
coordinate can round its fractional part. Geometry preparation now derives X
from the exact stored translated coordinate. Native SINT08 reproduces **two
inconsistent records among 32** before the fix and zero afterward. Original
records and process results remain in `prefix-before64-ranks1/` and the final
native campaign. No momentum, weight, stochastic history, solar-wind-relative
particle equation or physical timestep changes in this repair.

The eight-rank/four-thread Dmumu driver then exposed a second actual failure:
at step 34 the worker tried to read dt through an unallocated Eulerian AMR
block (`pic_field_line.cpp:2438` in the retained old source). Field-line and
Cartesian MPI decompositions differ. All segment worker paths now select the
same configured global timestep as the serial segment sweep and pass a null
Eulerian node. Node-attached block-local stepping remains unchanged; local
segment timesteps remain explicitly unsupported. The failure's actual exit 1
and logs are in `checkpoint-final/TCR02-continuous-threads4-ranks8/`.
Its exact tested source/object/archive/executable closure is retained as 1422
content-addressed files in `before-thread-dispatch-fix/manifest.json`.

The repaired eight-rank run then completed all 64 steps but its fresh resume
failed the unchanged byte-exact comparison. All geometry, particle histories
and source ledgers were identical; four binary64 sampled vertex sums differed
by one ULP. Material sampling now reads a canonical physical-input index,
without changing particle/list state or sampling formulas. SINT08 reverses the
same 32 manufactured records: the old traversal reports `exact=0` with actual
exit 1; the corrected native run reports `exact=1` and exit 0. The diagnostic
restores every touched buffer, exposure, cycle counter and original list before
continuing the original geometry tests. This correction preserves a complete
physical sample at the fixed dt; it does not add opportunities or relax exact
continuation. The failing final images and all 1422 tested source/build files
are retained in `checkpoint-dispatch-qualified/` and
`before-sampling-order-fix/manifest.json`.

The full threaded collector subsequently found one MFP mismatch: two complete
432-byte particle records on segment 48 had opposite list order. Independent
comparison confirmed an exact record multiset, with all other state identical.
Material serialization now orders complete records by immutable particle ID;
it does not discard or round a field. Keyed material sampling, inventory and
population control already define their own deterministic order. Fixed-line
shared-stream operators retain their semantic list order. CP05C01–03 additionally
reverse actual material lists, require the complete image to remain identical,
and restore both link directions before further host steps; MFP requires a
nonempty multi-record segment. The failed 43-PASS/1-FAIL report remains in
`checkpoint-qualified/`, and its full build snapshot in
`before-canonical-checkpoint-order/manifest.json`.

A real Make negative control also exposed `MAINOBJ` being expanded before its
definition. Main/debug objects could retain stale weak PIC template methods
after worker-count/header changes. The object list is now defined before the
dependency stanza; the old ordering fails the independent Make test and the
new ordering passes. All configured builds audit the full PIC-header consumer
closure, rather than selected mover objects only.

MCP12 rejects a valid populated restart on the last MPI rank after detached
geometry, ownership, flux/source plans and particles have been prepared. It
requires complete live-state and actual pointer/counter equality, then an
ordinary successful retry. The rejection-only count probe is prohibited in
production run purpose, has no CLI/environment hook, and cannot modify staged
state. Intrinsic primitive validation now precedes exact core-source equality;
both guards remain mandatory with the same admitted checkpoint set.

### Final qualification evidence

The authoritative final campaigns are `checkpoint-stage7-final/` (44 jobs),
`material-stage7-final/` (58 jobs) and `native-stage7-final/` (48 jobs).
Preparation alone is not evidence of execution. Their collectors require
actual MPI receipts, unchanged tested source/build identities and the complete requested catalogs; evaluated
mismatches remain FAIL, missing/crashed/unavailable cases remain ERROR.

| Final executed gate | Result | Evidence below this root |
| --- | --- | --- |
| Material CP05C/F, isolated CP05O/OF, missing-fixture rejection and threaded full-driver restart | **44 PASS, 0 FAIL/SKIP/ERROR**; 20 logical IDs | `checkpoint-stage7-final/summary.json`, `junit.xml` |
| Moving material restart, sampling phases/topology and 1↔8 repartition, all three movers | **58 PASS, 0 FAIL/SKIP/ERROR** | `material-stage7-final/summary.json`, `junit.xml` |
| Native geometry, connectivity, transport, source and triggered merger | **48 MPI jobs PASS**; 15 IDs pass at each rank count 1/2/4/8 | `native-stage7-final/summary.json`, `junit.xml` |
| Applied production CLI before/after initialization changes | **PCLI01–12 PASS**, twice each on two ranks | `pcli-before/`, `pcli-stage7-final/artifacts/results.json`, `process-result.json` |
| Strict numerical ASan/UBSan suite, LSan disabled | **30 PASS** | `numerical-final/summary.json`, `numerical-final.log` |
| Shared harness plus field-line kernels | **24 PASS, 0 FAIL/SKIP/ERROR** | `shared-final/amps-core-tests.json` |
| Material-reader, dependency and stage-7-runner negative controls | **16 PASS** | `runner-units-qualified.log` |
| Shared orchestrator readers | **38 PASS** | `orchestrator-final.log` |
| CLI/registry focused controls | **10 PASS**, including CLI01–06 and expected hidden-failure rejection | `cli-python-regression.log` |
| Python orchestration and plot/report controls | **27 PASS** with system Python and its installed matplotlib | `python-runner-system-python-final.log` |
| Catalog and source traceability | **63 suites / 1277 cases**, checks PASS | `catalog-stage7-final.log`, `traceability-final.json` |

All 150 original planned argv/actual process receipts were checked again in
`stage7-final-verification.json`. MCP13 intentionally requires native ERROR and
actual MPI exit 2; every other final MPI job exits 0. It tests rejection, not
transport acceptance. Complete final same-decomposition checkpoint comparisons
remain byte-exact, including the **21,425,968-byte** eight-rank/four-thread image.
The maximum independent relative budget residual is **1.466478e-15** for
threaded production and **1.070965e-15** for the material campaign,
against the unchanged **1e-10** bound. MCP12 stages 19 real particle records
and rejects on rank 0/7 for MPI 1/8, then successfully reloads step 4.

The three actual production builds (`native/`, `native-thread2/`,
`native-thread4/`) compiled successfully without compiler warnings in the final
changed-unit build logs. Their exact source/object/archive/executable/helper
closure is retained as **1422 content-addressed files** in
`stage7-qualified-build/manifest.json`. The ordinary checkout Makefile remains
byte-identical to `baseline-Makefile`; no reference solution or Git index was
changed. Tests used these owned configured builds, not the unrelated checkout
`./amps` executable.

Two preliminary external readers were accidentally launched before their
writers completed. Their actual missing-checkpoint exit-2 failures are retained
in `premature-MPI-launches/`; the final campaign executes all dependency edges
before starting readers. No physics code or gate was changed for that error.
The first Python runner-unit attempt used a Python installation without
matplotlib. System Python supplies the required dependency; the final 27-case
run also verifies the repaired backward-compatible source-suite argument
default. Original failed attempts remain retained.

Exact prepared MPI argv and working directories are in each `prepared.json`.
The 44-job plan records its dependency edges directly; material-reader edges
were derived from their named writer paths and are retained in
`external-execution-plan.json`. Actual launched commands and return codes are
in per-job `process-result.json`. The final collectors were executed from AMPS as:

```sh
evidence_root=/data/vtenishe/TEST/sep-stage7-completion-20261122T182436Z
python3 srcSEP/test/run_stage7_checkpoint_tests.py \
  --output "$evidence_root/checkpoint-stage7-final" --collect
python3 srcSEP/test/run_native_standalone_integration.py \
  --build-root "$evidence_root/native" --ranks 1 2 4 8 \
  --output-dir "$evidence_root/native-stage7-final" --collect
python3 srcSEP/test/run_material_checkpoint_tests.py \
  --output "$evidence_root/material-stage7-final" --collect
```

The first command requires the complete 44-job catalog. Independent shared
orchestrator readers also imported **44/44** and **58/58** final role results
as PASS. The maintained [test instructions](../test/README.md) document fresh
preparation/execution through Make, including explicit two-/four-worker builds.

The implementation changes are in `src/pic/pic_field_line.cpp`,
`srcSEP/physics_checkpoint.{h,cpp}`, `srcSEP/main.cpp`,
`srcSEP/util/sep_checkpoint_native_validation.{h,cpp}` and the application
makefile. Maintained native fixtures, the 44-/58-job runners, negative controls,
shared registry/traceability entries and orchestrator readers accompany them.
Detailed code comments and the PIC/application/test READMEs explain the frame,
timestep, canonical ordering, rollback and numerical contracts.

Suggested commit: `Complete SEP material restart qualification across MPI and threads`.


The new cases are integral to maintained validation: Make targets
`test-stage7-checkpoint-native` and `test-material-checkpoint-native`, shared
suites `CORE-SEP-STAGE7-CHECKPOINT-001` and
`CORE-SEP-MATERIAL-CHECKPOINT-001`, and matching orchestrator steps.
Routine source tests require the negative-control/dependency tests. The old
bare fixed-line checkpoint-group orchestrator launch is replaced by the full
stage-7 gate; all stable CP05C/F assertions remain included. A required skipped
role cannot be hidden by another passing role.

The admitted physics scope is one moving prescribed reduced SSE/Parker line,
keyed upstream escape at explicit reference distance and absorbing returns.
Isolated two-branch operators do not qualify active particle-wave feedback or
four-polarization closure. Short 64-step threaded restart runs do not replace
the stage-9 400-step output and current AMPS–BATS-R-US coupling campaigns.
Other backends, downstream transport and broader release qualification remain
stages 8–9.

## Stages 5–6 update, 2026-11-22

**Stage 6 COMPLETE, 2026-11-22.** All three final real-MPI one-rank application
runs completed 400 physical steps and returned actual exit code zero in
`application-predictable400/`. The maintained complete-driver gate passed
SAPP01–03 and all eleven audit negative controls (`stage6-final-gate.log`).
The admitted scope is the declared synthetic, prescribed upstream reduced
SSE/Parker consumer; broader backends and moving restart remain later stages.

| Final gate | Result | Evidence below the root given below |
| --- | --- | --- |
| Strict numerical integration/source/boundary cases | 30 PASS | `predictable-boundary-after/summary.json` |
| ASan/UBSan, including changed core paths; LSan disabled | 30 PASS | `sanitizers-predictable30/summary.json` |
| Actual native geometry/source/three movers on MPI 1 and 4 | 14 IDs PASS at each rank count; 22 jobs exit 0 | `native-fixtures-predictable/summary.json` |
| Numerical suite plus orchestrator self-tests | 65 PASS | `orchestrated-predictable30/run_all_test.json` |
| Shared validation registry | 30 individual PASS | `shared-predictable30/amps-core-tests.json` |
| Complete-driver 400-step campaigns, all three movers | SAPP01–03 PASS; actual process exits 0 | `application-predictable400/verified-results.json` |
| Deliberately corrupted evidence controls | 11 PASS | `application-predictable400/negative-control-results.json` |

Every final application used the unchanged physical dt **300 s**, sampling
window **400**, deadline **1800 wall seconds** and conservation bound **1e-10**.
All reached physical epoch **120000 s**. The output correctly represents the
departure sampling epoch **119700 s**, before geometry commit 400. Its 357
segments retain independently verified lifetimes of **5–400 samples**, each
with fixed dt 300 s. All movers receive exactly the same source history.
The largest independently recomputed relative number/birth/current-energy
residual is **2.203e-15**, below the original bound. No reference was regenerated.

The complete gate was executed as:

```sh
make -C srcSEP test-standalone-integration-application \
  EVIDENCE_DIR=/data/vtenishe/TEST/sep-integration-stages5-6-20261122T023833Z/application-predictable400
```

`application-predictable400/launch-plan.json` contains each exact MPI command;
its `build-sha256.json` and per-mode `process-result.json` files retain tested
identities and actual process receipts. JSON/JUnit reports accompany the
numerical, native, shared-registry and complete-driver results. The checker
requires every recorded source/header/archive hash to match, so old successful
output cannot certify an edited checkout. Its source-identity negative control
failed before that protection and passes after it.

The exact iteration-84 trajectory inputs now form SINT30: the retained old
controller fails, while endpoint-preserving bisection passes with 1.27e-6 m
endpoint error against an independent converged midpoint reference and the
unchanged 0.008503654-m tolerance. The earlier iteration-52 MFP closed-endpoint
defect is covered by SINT28 (66/4096 old-kernel misses, zero after correction).
Both original failures and actual nonzero process receipts remain retained.

Profiling identified repeated vertex trajectories and identical PFSS source
queries. Exact source-point memoization and validated search hints preserve
canonical field values and root coordinates (SINT31 and native full-range/
full-polyline comparisons). Predictable boundary-refinement history avoids
repeatedly rejected coarse proposals, while retaining all actual geometry,
topology, spatial-error, minimum-step and 64-attempt checks before noise draws.
SINT32 fails before that hint and passes afterward: 149 checks replace 1024
cold checks for the same 128 accepted manufactured intervals at a 1-mm bound.
The physical dt, sampling schedule, 1800-s campaign deadline and conservation
bounds remain unchanged. The original unoptimized MFP campaign timed out with
actual exit 124; it is not reported as a numerical/physics PASS.

The maintained `test-standalone-integration-application` Make target checks
the actual complete-driver receipts, 400 clock endpoints, source equality,
segment sampling lifetimes reconstructed from topology commits, active output
metadata and independently recomputed number/birth/current-energy budgets.
Eleven isolated corruptions of a completed real synthetic run must be rejected,
including a changed source identity. All recorded source/build hashes are
checked; unchanged old output cannot validate a subsequently edited checkout.
The shared wrapper's stale case list was also corrected; its new catalog
self-test fails before the correction and passes afterward. No subset report
is accepted to hide a missing numerical case.

Earlier campaigns, including the actual MFP failures and runtime timeouts,
remain in the evidence history. The final revision passes the complete gate
above. This report makes no observational, moving-restart, self-consistent-wave
or downstream-volume qualification claim.

Current evidence root:

```text
/data/vtenishe/TEST/sep-integration-stages5-6-20261122T023833Z/
```

`baseline.json`, `baseline.patch` and `baseline-Makefile` retain the starting
dirty tree. This batch changes the source/boundary consumer, three native
mover integrations, source configuration/history, focused-escape CDF, closed
RK endpoint arithmetic, and exact search/caching optimizations. All physical
and statistical acceptance bounds remain unchanged. The existing root
Makefile and prior regression references are preserved.
`final-stage6-audit.json`, `final-stage6-tracked.patch` and
`final-stage6-status.txt` retain the final source identities and checkout review.

### Implemented physical contract

- The existing core local-compression DSA **number** law is
  `dN/dp proportional to p^(2-3r/(r-1))`. The parent angular distribution is
  isotropic in the upstream plasma frame. Focused escaping directions use
  positive normal flux; Parker uses the angular marginal and an isotropic
  diagnostic direction. No second shock-acceleration operator is applied.
- The selected strength is net upstream first passage at an explicit positive
  normal offset `L_ref`. Release points intersect the actual material line at
  the actual stratified birth time. There is no numerical production default
  for `L_ref`, no off-line normal placement, no electron/proton substitution,
  and no source renormalization after return loss.
- Default return is absorption/deletion. `pass` and gross-emission selections
  take the requested AMPS startup error path; no downstream volume or missing
  accelerator-to-reference escape law is fabricated.
- The native injection callback executes **after departure-time sampling**
  inside `PIC::TimeStep`. Each birth consumes its actual remaining duration;
  every physical AMPS step and ordinary sampling opportunity remain unchanged.
  Relative field-line displacement has no added solar-wind velocity. Actual
  core U remains in material motion, gradients/cooling, source flux and shock
  motion relative to the line.
- Parker retains its sampled endpoint and includes Brownian-bridge hidden
  shock returns. Dmumu stops its accepted streaming interpolant at a hit. MFP
  checks every between-event deterministic interval before drawing a later
  scatter/clock. Internal front-curvature refinement precedes stochastic draws.
- Signed accounting separates injection, retained particles, shock returns,
  domain exits, geometry retirement, deterministic background work and
  prescribed scattering-frame work. Immutable birth kinetic energy is carried
  in particle transport schema 5; incompatible old records reject.
- The dedicated `test/sep_reduced_front__field_line` configuration compiles a
  physical proton mass. The shared legacy H_PLUS=1-amu catalog is preserved;
  this source refuses to treat that mass as a proton silently.

### Earlier bounded qualification history

| Gate | Observed result | Evidence |
| --- | --- | --- |
| Strict source/adapter/geometry/transport suite | 25 PASS | `parallel-bracket-after/summary.json` |
| ASan/UBSan, including changed ambient/source/RH core paths | 25 PASS; LeakSanitizer disabled | `sanitizers-core-final/summary.json` |
| Actual native geometry/connectivity/source/three movers | 13 IDs PASS on both 1 and 4 ranks, 20 launched MPI jobs | `native-fixtures-cache/summary.json` |
| Actual full PIC driver, same source and all three movers | 48/48 physical steps each, actual exit 0, endpoint 14400 s | `application-cache/verified-results.json` |
| Applied production CLI | PCLI01–12 PASS, twice on 2 ranks | `pcli-cache/results.json` |
| Real CLI unsupported/missing-input traps | 3 expected rejections PASS | `negative-cli-results.json` |
| Core coronal-CME stage-14 registry | 219 PASS, 3 data-dependent SKIP, 0 FAIL | `core-parallel-stage14-complete.log` |
| Reduced-front regression after axis shortcut | 57 PASS, 0 FAIL/SKIP/ERROR | `polar-provider-after/core-polar-provider-reduced.log` |
| Background configuration/history | CMBGU01/02 PASS | `core-polar-provider.log` (earlier reduced-harness setup error retained separately) |

The full-driver postcheck independently recomputes number and signed energy
balances from each of the 48 endpoint rows, verifies the exact 300-s clock
sequence, and requires the original 400-step sampling counter on every step.
The largest relative ledger residual is below 1.3e-15 against the unchanged
1e-10 bound. All three movers receive identical candidate/committed source
sequences. These 48 steps do **not** complete a 400-step output interval.

The focused-escape CDF independent test failed before the threshold/inverse
correction (`source-before/`). The strict closed-interval RK case likewise
failed before explicit endpoint propagation (`closed-interval-before/`).
Those failures and their passing after-results are retained. Native MPI source
fixtures initially left outbound lists unjoined; the fixture now joins all
destination lists and executes the existing production migration. Its unchanged
production window-close check remains authoritative.

Earlier slow application attempts and CPU-placement timeouts remain in
`application/`, `application-polar-provider/` and `application-unbound/`.
They are not passes. Direct debugger sampling found repeated core RH scans in
the boundary path (`profile-source/profile.log`). Exactly parallel/normal
incident states now use the gas-dynamic identity only to select an existing
scan bracket, with the same bisection, weak-root cutoff, entropy, characteristic
and RH checks. Exact-axis winding removal and exact-epoch per-worker root
memoization introduce no substitute background, approximate timestamp or
relaxed accuracy. Native repeated queries remain bitwise identical and pass
the independent full-intersection 1-mm gate across material-owner changes.

The 400-step campaign is a distinct longer check with the same physical dt,
source, 400-step sampling cadence and gates; its budget is 1800 s versus 600 s
for the 48-step check. No failed 48-step deadline was increased. Its actual
receipts and complete output must be inspected before assigning a result.

### Remaining scope boundaries

The implementation admits one transverse supported leading connection on one
material line. Tangency, multiple absorbing fronts and occupied connectivity/
criticality creation/loss require additional qualified operators and reject.
Gross conversion, pass-through/downstream transport, local wave evolution,
moving restart and stage-7 population-merger accounting remain unsupported.
An actual conservative-merger trigger rejects before losing the new immutable
birth-energy account; no population-control threshold is increased. Exact
splitting retains that account. Prescribed read-only scattering is the admitted
first consumer; this is not self-consistent polarized-wave qualification.

The synthetic event, labelled monoenergy and `L_ref=1e8 m` verification input
are not observational calibration. Fresh AMPS–BATS-R-US coupling and additional
existing-backend qualification remain the separately authorized later stages.

## Historical stages 0–4 report

This report covers the user-authorized stages 0–4 of
[`INTEGRATION_ROADMAP_STANDALON.md`](../INTEGRATION_ROADMAP_STANDALON.md).
The selected physical authority is the **existing shared reduced SSE front
and Parker exterior**, with prescribed, read-only background. No new shock
model, source/acceleration law or downstream-volume closure is introduced.

The implementation has bounded native and numerical qualification. This is
not a claim that the entire repository regression campaign is green.

## Source and evidence

Baseline HEAD: `e65d7d16d89038019038f497f92040ae9d54deb1`, with the user's
pre-existing dirty tree preserved. Complete evidence is outside the source
package:

```text
/data/vtenishe/TEST/sep-integration-stages0-4-20261121T203723Z/
```

`baseline-source.json`, `baseline.patch` and `baseline-Makefile` retain the
starting state. `final-source-audit.json`, `final-native-campaign/summary.json`,
source/build hashes and exact launch receipts identify the tested final state.
The root Makefile remains byte-identical to its retained starting version.
No reference solution, physical timestep, sampling cadence, numerical gate or
statistical acceptance bound was relaxed. Nothing was staged or committed.

Native builds use the owned `native/` configured checkout and real MPI. One
initial command accidentally relinked the primary ignored executable from its
existing objects; that discarded relink is not qualification evidence. All
native results below use the separately rebuilt owned tree.

## Stage dispositions

| Stage | Implemented behavior | Observed bounded gate |
| --- | --- | --- |
| 0 | Frozen source/configuration; actual prerequisite registries and capability inventory | All five named suites executed; individual limitations below remain visible |
| 1 | Explicit species densities, core thermal/derivative state, immutable full-interval views and prepared publication | SINT01–04 PASS; legacy adapter checks PASS |
| 2 | Core-B initialization, actual-U material trajectories, both endpoints, B-traced empty insertion, tail loss/exposure, atomic native publication | SINT05–09 and SINT07R PASS; actual native 1-/4-rank container transactions PASS |
| 3 | Finite supported SSE roots, multiple/endpoint/tangent classification, material-time queries, entry/loss of connection | SINT10–13 and native SINT13N PASS |
| 4 | All three native solvers on the same material background, tube metric, gyrotropic coefficients and independent limits | SINT14–16 and native SINT14N–16N PASS; existing PARK/FTED/FTEM suites PASS |

Native stage-4 particles are explicit manufactured initial conditions. The
production reduced path remains zero-particle until stage-5 source and
unmodelled-downstream boundary choices are supplied. The stage-4 MFP native
characteristic is conditioned to survive without an event; existing FTEM
cases independently grade event clocks, residual hazards and wave-frame
conservation. This is not an adopted physical source or feedback model.

## Stage-0 capability inventory

The five-suite baseline produced **341 PASS, 181 SKIP, 1 ERROR, 0 FAIL** from
523 individual registered rows:

| Existing backend/suite | PASS | SKIP | ERROR | Admission consequence |
| --- | ---: | ---: | ---: | --- |
| Reduced front + ambient | 62 | 35 | 0 | Reduced physical capability supported; skipped srcSEP3D native/restart rows remain unqualified by this baseline |
| SWCME1D | 4 | 143 | 0 | Maintained adapter/backend retained; dormant extended rows are not admitted production physics |
| Coronal-CME neutral implementation | 219 | 3 | 0 | Maintained exchange capability retained; moving srcSEP admission is stage 8 |
| Full corona–SWCME volume | 45 | 0 | 1 | BG3D-4/G1 remains open; no full-volume transport claim |
| Generic PIC field line | 11 | 0 | 0 | Existing container/kernel gates preserved |

The full-volume error is `CSWC-FULL-CSWC0628-MASS-TIME`:
`test_bg3d4_piston_refinement` exceeded its existing **180-second** routine
limit. It occurred before this implementation. The limit was not changed,
and partial refinement output was not converted to PASS.

The stock example's self-consistent polarized production selection still
rejects its unfinished closure. A deliberately prescribed zero-particle
baseline reached the old unregistered-material-owner failure; stages 1–2
replace that failure with the registered complete native owner.

## Numerical corrections and independent evidence

Parker physical-number representatives carry `4 pi p^2 A f ds dp` with
`A |B|=Phi`. Their material diffusion drift is
`partial_s kappa + kappa partial_s ln A`; coordinate-to-arc conversion does not
supply this term. Native transport now passes `-partial_s ln|B|` explicitly,
while constant-area manufactured helper callers retain the documented zero
metric default. Cooling remains `p exp(-div(U) dt/3)`.

The first native focused-diffusion comparison failed: streaming used the
final cooled speed for the whole interval. At unchanged 300 s its displacement
was 179.536559843 m. Coupled midpoint integration of log(p)/mu, half-speed
coefficient prediction and trapezoidal parallel velocity corrected it to
179.882721554 m; the native MFP no-event path gives 179.883164057 m. The original
0.1-m displacement and momentum/pitch bounds were retained. SINT15 adds an
independent exact relativistic cooling/streaming reference and internal
second-order refinement; the physical step is never shortened.

Other independently graded results include:

- Actual-core radial/azimuthal trajectories against separately implemented
  scalar quadrature/inversion, inverse-map error, mass continuity and ideal
  magnetic induction along material worldlines.
- Curved material B-tangent errors 3.04159e-7 and 7.60985e-8 for 9/17 vertices
  at the same 600-s endpoint: approximately fourfold error reduction.
- Polar SSE position/speed, off-axis sphere roots, rear/support rejection,
  endpoint ownership, two genuine intersections and tangency.
- Core composition with alpha/proton abundance 0.1: `ne=1.2 np`, preventing
  a pure-H assumption from hiding incorrect native density mapping.
- Parker exact conditional drift and Gaussian moments with N=32768,
  seed 4289 and fixed 6-sigma bounds.
- A common diffusion comparison with N=2048, seed 5297 and equivalent
  `D_mumu=nu(1-mu^2)/2`, `lambda=v/nu`, `kappa=v^2/(3nu)`; fixed finite-time
  moment references and declared confidence/discretization bounds.
- Actual native conserved flux/volume, material particle coordinate remap,
  empty insertion, represented tail loss and distinct segment exposure.
- A late last-vertex fault on one rank rejects both state installation and
  motion without publishing native geometry, metadata or core cache changes.

Isolated negative controls remove the tube drift, conflate electron/proton
density, or restore final-speed streaming. Each triggers its required failing
physics case (SINT14, SINT02, SINT15). Their expected failures are retained in
`negative-controls/`; none changes production code or reference data.

## Executed validation

Commands run from the AMPS root unless a configured-build directory is named.
Exact argument arrays, logs, hashes and JSON/JUnit reports are retained with
each campaign.

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --build-jobs 4 --test CORE-SEP-STANDALONE-INTEGRATION-001 \
  --test CORE-MODEL-SEP-CORONA-SWCME-REDUCED-001 \
  --test CORE-PIC-FIELD-LINE-001 --output-dir <fresh-shared-output>
python3 srcSEP/test/run_standalone_integration_tests.py --sanitizers \
  --output-dir <fresh-sanitizer-output>
python3 srcSEP/test/run_tests.py --suite standalone-integration \
  --output-dir <fresh-unified-output> --plot none
python3 tools/sep_test_orchestrator.py --app srcsep --skip-build \
  --only standalone-core-integration --only orchestrator-selftest \
  --output-dir <fresh-orchestrator-output>
make -C srcSEP test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit \
  test-coefficients-unit test-transport-common-unit test-state-unit test-cli-unit
```

Observed results:

- Shared bounded campaign: **88 PASS, 35 SKIP, 0 FAIL, 0 ERROR**. The 35 skips
  are the existing explicit srcSEP3D native/restart omissions, not these native
  srcSEP fixtures.
- New integration catalog: **15/15 PASS**, including strict-warning and
  ASan/UBSan runs. Unified source runner: 16 PASS including its suite wrapper.
- Application orchestrator: **49/49 PASS** (15 physical/software cases plus
  34 orchestrator checks).
- Existing controlled solver catalogs: PARK 7/7, FTED 8/8, FTEM 8/8 PASS,
  including their strict-warning ASan/UBSan checks.
- Native diagnostic cases: all six IDs PASS separately on **1 and 4 real MPI
  ranks**. The full native application also exits zero after four complete
  300-s steps on both rank counts, at epochs 300, 600, 900, 1200 s, with global
  particle count zero at every endpoint.
- Native production CLI PCLI01–12: **12/12 PASS**, each repeated twice on two
  MPI ranks; these are applied-software-state checks, not feedback physics.
- Legacy reduced adapter and initialization checks PASS. Background derivative
  and flux-state sanitizer checks PASS with LeakSanitizer explicitly unavailable.
- Root Makefile byte comparison, source packaging required-member check and
  `git diff --check` PASS. The packaging required-member check is narrower than
  the broken Make target described below.

## Explicit infrastructure and outstanding boundaries

The Python-launched MPI campaign cannot create Open MPI interface sockets
inside this session's sandbox. It records six ERROR cases and a failing suite
wrapper; those receipts remain in `final-native-unified/`. The same rebuilt
native binary was separately launched with direct real `mpiexec -n 1` and
`-n 4` commands. Its successful rank logs, process exits and exact commands
are in `final-native-campaign/`. Direct receipts do not rewrite the automated
ERROR report. The registered runner remains executable on a host with the
required MPI capability.

LeakSanitizer reports that it cannot run under ptrace. The original failing
receipt is preserved; ASan and UBSan were run with `ASAN_OPTIONS=detect_leaks=0`
and **no leak-check qualification** is claimed.

An additional existing `test-package-hygiene-unit` Make target is broken:
it passes `--root --self-test` to the repository's `sep_package_hygiene.py`,
which accepts a positional root and has no self-test option. It remains an
explicit unrelated validation-interface error; its self-test requirement was
not removed to obtain a pass. The separate existing positional-root,
`--allow-generated` required-member check passed and is not a replacement for
that missing gate.

One existing native build warning remains in
`turbulence_wave_number_resolved.cpp:1580` (`-Wformat-truncation`); that source
was not changed by these stages. New pure adapter/numerical fixtures build
with warnings as errors.

Stages 5+ have not begun. No physical source/release-side/acceleration or
unmodelled-downstream policy is inferred. PFSS material evolution, moving
local-wave remap, full CME-volume admission, moving-geometry restart,
thread/decomposition production campaigns, complete 1-AU SEP transport,
observational validation and fresh AMPS–BATS-R-US coupling qualification remain
their separately specified later gates.
