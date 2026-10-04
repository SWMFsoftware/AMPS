# CME execution plan — full background before particles

Revision 2026-10-03 (America/Chicago), sequence 2.0. This replaces the earlier
CODEX_CME_PLAN.md interleaved L0–L7 execution order. Read the updated shared
roadmap revision 1.3, especially Section 0, and the maintained coronal and
Corona–SWCME composite physics specifications. This file records execution
and evidence; the specifications own the equations and closure contracts.

## Active user decision

1. Complete the spatial CME plasma/IMF background in **srcSEP3D**, from a
   prescribed low-coronal launch through continuous SWCME handoff to 1 AU.
2. Once that works, implement the same background coupling in **srcSEP**.
3. Only after both background milestones pass, implement particles first in
   **srcSEP3D**.
4. Once the srcSEP3D particle milestone works, implement particles in **srcSEP**.

**Current authorization ends at G1.** After BG3D-0 through BG3D-9 all pass,
finish the G1 plan/README/evidence/limitations report and stop. Do not begin
G2, G3 or G4 without a new explicit user instruction. A required SKIP or any
unresolved failure means G1 remains incomplete.

**A moving front with quiet ambient is insufficient.** The first deliverable
needs spatial plasma and magnetic fields in the ambient, sheath, ejecta and
other sampled CME-related regions, plus shock positions and local parameters.
Do not fill unavailable downstream volume with quiet ambient and call it a
CME solution. A local RH jump is only an interface state. Full spatial state
must come from the selected analytical regional closure and its validation.

The model prescribes eruption motion/expansion and selected analytical plasma
closures; it does not predict instability onset or solve global MHD. Inputs
for mass, composition, thermodynamics and magnetic structure are required by
the selected volume closure. A spheromak is not a mandatory choice. Required
force/heating/work and physical residual budgets are explicit limitations.

## Working directories and preserved behavior

- Work in `~/Mars2/AMPS`, containing the verified pre-L0–L7 baseline.
- `~/Mars1/AMPS` is a read-only reference. Preserve it and its uncommitted work.
- Check the actual baseline archive/commit and source fingerprints. This plan
  has not inspected the user's latest checkout or MPI receipts.
- First reproduce a baseline build and bounded relevant tests in Mars2.
- Transfer useful background kernels selectively, then revalidate them in
  Mars2. Restore particle-only files/sections to the baseline. Do not remove
  baseline particle tests, schemas, CLI behavior or legacy SWCME support.
- Keep shared physics under `src/models`, application adapters thin, and the
  generic PIC core functional for other applications. A composite facade may
  retain its existing public location; one component owns each physical state.
- Do not commit or push. Do not delegate unless explicitly authorized.

## Gate and stage ledger

The initial status below is an honest scheduling state, not a finding that all
existing implementations are absent. Audit actual reusable code and record
its new evidence before changing an item to implemented or qualified.

| Gate | Required stages | Implementation status | Validation status | May start next gate? |
|---|---|---|---|---|
| G1: srcSEP3D full background | BG3D-0 through BG3D-9 | BG3D-0/1/2 qualified; BG3D-3 shock subset qualified, contact subset reopened; BG3D-4 partial/blocked | OPEN | No |
| G2: srcSEP background | BGSEP-0 through BGSEP-3 | Deferred until G1 | OPEN | No |
| G3: srcSEP3D particles | P3D-0 through P3D-4 | Deferred until G2 | OPEN | No |
| G4: srcSEP particles | PSEP-0 through PSEP-3 | Deferred until G3 | OPEN | N/A |

## G1 ordered implementation

| Stage | Work | Required evidence before dependent work |
|---|---|---|
| BG3D-0 | Baseline/source audit; changed-file classification; reconcile coronal/composite physics and coverage; freeze supported closures | Baseline tests/build, public-library independence, particle-file hashes, compatibility and physical contract table |
| BG3D-1 | Background-only configuration, owning checksummed histories/assets, launch geometry and full-time/radial support | Actual parser/resolver and independent history/attachment/clipping/capability negatives |
| BG3D-2 | Canonical ambient plasma/IMF/EOS, thermodynamics, topology, valid smooth/one-sided derivatives | Independent EOS/wind/flux/gradient/coverage cases; declared plasma support through >1 AU |
| BG3D-3 | Full shock/contact/ejecta geometry, local normal motion and admissible shock states | Independent geometry and weak/oblique/subfast shock cases; compression, obliquity and conservation residuals |
| BG3D-4 | Spatial analytical sheath with shock admission, valid initial inventory, downstream material mapping and all exits | Planar/expanding spatial fixtures, shock-limit matching, independent inventory/interface/convergence checks |
| BG3D-5 | Spatial ejecta and remaining sampled regions, magnetic/thermal/material closure and interface matching | Positive maps, independent mass/flux/induction/EOS checks, thermal/work/force budgets and full requested coverage |
| BG3D-6 | Dense-time matched crossing; complete shape/material transition and outer SWCME continuation | Direct/finite-transition cases, derivative/reset/incompatible-shape negatives and continuous mass/flux/geometry state |
| BG3D-7 | Native runtime provider, DATAFILE slots, actual owner readback/received halos, private preparation and readiness | Actual evolving 1-/4-rank zero-particle runs, all-region probes, both slots/derivatives/categories and failure controls |
| BG3D-8 | Separate short handoff-smoke deck; multi-epoch receipts; background-only checkpoint/repartition | Actual pre/inside/post transition, rank agreement, no relaunch/inventory loss and tested late-write policy |
| BG3D-9 | Separate full low-corona-to-1-AU deck and campaign; diagnostics/comments/READMEs | Actual initial/crossing evidence, full spatial plasma, time/mesh/material/integrator refinement, 1-/4-rank and restart evidence |

Section 0.4 of the roadmap specifies the detailed algorithms, identities and
negative cases. BG3D-3 is an intermediate diagnostic checkpoint. G1 cannot
close without the regional plasma work in BG3D-4/5. Do not start srcSEP
adapters or particle algorithms merely because a front or initialization passes.

## Critical implementation contracts

1. **One physical state:** preserve event identity and complete geometry and
   material state. SWCME continuation starts from the matched state, without
   radius/speed/shape reset or a new ambient normalization. Retain one composite
   native authority and component provenance through the entire run.
2. **Full geometry:** an equal apex does not establish equal surfaces or
   velocities. Preserve center/axes/orientation/support and feature identity.
   Differentiate the transition; include its weight derivative. Test normal
   velocity, regularity and any claimed acceleration continuity.
3. **Full plasma:** explicitly implement every sampled supported region. Use
   a spatial sheath/ejecta closure, with mass/flux/thermodynamic/interface
   constraints. No arbitrary field-vector blending, compression floor,
   replacement Parker field or post hoc conservation renormalization.
4. **Honest shock:** determine local fast-shock admissibility from upstream
   relative normal speed and canonical EOS/field. A geometrical front can be
   subfast; preserve typed status and stable physical geometry. A genuine
   solver failure is not an excuse to replace the state with ambient.
5. **Independent background epoch:** geometry/plasma/regions/shock ownership
   and readiness do not require sources, particle counts or feedback state.
   Preserve any future richer epoch separately. A selected fixed prescribed
   wave pressure/energy is background closure, not implicit SEP feedback.
6. **Actual native evidence:** owner/ghost agreement uses actual received
   values and both relevant time slots. Unit fixtures and final-only receipts
   cannot qualify multi-epoch MPI behavior.
7. **Failure integrity:** private candidate rejection preserves committed
   identities and values. Late writes/halo errors require tested staging/
   restoration or declared fail-stop semantics. Never advance the clock with
   stale fields or advertise rollback that the memory layout cannot provide.

## Deferred G2, G3 and G4

G2 implements independent srcSEP dispatch/build, real traced line/node
coordinates, regional plasma/derivative reduction, signed B versus arc tangent,
current/previous native line data, metrics/ownership/readiness, zero-particle
MPI/restart and matched background probes. It consumes the shared provider,
not srcSEP3D outputs. The qualified frozen/moving line representation is stated
explicitly. Begin only after G1 passes.

G3 implements physical release/reference/source measures, event integration,
native nonzero allocation, srcSEP3D movers/scattering/interfaces/boundaries,
number/energy/loss accounting, statistical and numerical convergence and
particle MPI/restart. Requested wave feedback has its own closure and energy
tests. Begin only after G2 passes.

G4 implements finite-tube source projection, native srcSEP births and transport,
streaming/contact/coordinate handling, particle MPI/restart and matched
Parker/focused limits. Begin only after G3 passes. General 3-D perpendicular
transport is not certified by a field-aligned comparison.

Historical L0–L13 assertions retain their identity and full requirements.
Source-free checks do not close particle assertions. Background restart is
needed now where outer integration/material inventory requires it; particle
restart remains deferred. New CMBG IDs in the roadmap are proposed: reconcile
the current registry before implementing them, and report registrations and
execution separately.

## Build, test and campaign rules

- Before **every** native AMPS rebuild, execute and record the following exact
  preflight in order: (1) confirm the working directory is
  `/home/vtenishe/Mars2/AMPS` and verify that no build or test is running;
  (2) remove only the generated `build/` directory with `rm -rf -- build`;
  (3) regenerate the selected application's configuration and production
  hooks through the established workflow; and (4) compile through the
  established top-level workflow with `-j16`. Preserve source files,
  `Makefile.conf`, site settings/libraries and all of Mars1.
- Read existing site/application build instructions and inspect actual decks.
  Do not invent configuration switches or a missing srcSEP compile deck.
- Before **every native rebuild**, confirm `~/Mars2/AMPS`, remove only its
  generated `build/`, regenerate selected application configuration and
  production hooks, then compile through the established workflow, e.g.
  parallelism `-j16`. Never modify generated copies as the only permanent fix.
- Preserve site Makefiles, configuration scripts, libraries and baseline
  source files. Check for hard-coded Mars1 build/output paths. Save binary
  fingerprints/logs before switching the compiled application.
- Use an available allocation for native MPI. Required background runs use
  1 and 4 ranks, zero allocated/injected particles throughout, and real
  evolving epochs. Add relevant thread settings to evidence.
- Keep handoff smoke and full low-corona-to-1-AU decks separate. Choose
  physical duration from the trajectory. Ten short steps or an outer-only run
  cannot qualify a low-coronal launch-to-1-AU claim.
- Receipts retain pre/inside/post handoff and region coverage, actual shock
  parameters, owner/native/received ghost values, generation bindings and
  input/closure/binary identities. Export plasma/IMF and regional histories,
  not only a front outline.
- Compare against independent equations/data, not production-generated
  expected answers. Freeze tolerances and physical discrepancy budgets before
  grading. Use three refinement levels for convergence claims and distinguish
  numerical error, physical parameter sensitivity and prescribed-model error.
- Required SKIP is not PASS. Missing tools/data, incomplete implementation and
  physically rejected profiles are separate categories. Preserve failures.
- Every new or modified CME/background code path, including earlier code in
  this task, requires reasoning-level comments covering equations, assumptions,
  approximations, validity, SI units/frame/sign conventions, algorithms,
  tolerances and failure behavior. Relevant READMEs must give inputs and
  commented examples, exact build/run/test commands, validation criteria,
  outputs, limitations, implemented-versus-validated status, regional/shock
  semantics, handoff/outer evolution and (when implemented) epoch/native/MPI
  ownership. Do not document planned native storage or synchronization as
  implemented. This rule applies to every subsequent change and must preserve
  unrelated and baseline particle documentation.

## Required per-stage progress record

| Field | Record |
|---|---|
| Stage/current action | Exact BG3D/BGSEP/P3D/PSEP stage and next dependency |
| Implementation | Not started / reused pending verification / partial / implemented |
| Validation | Not run / blocked / failed / partial evidence / qualified |
| Changes | Files and functions, physical reason, source of reused algorithms |
| Tests | IDs, exact commands, PASS/FAIL/SKIP/ERROR, negative controls, refinement metrics |
| Provenance | Source/deck/asset/closure/binary hashes, MPI ranks/threads, output directory |
| Physics limits | Unsupported regions/operations, prescribed-force/heat/work residuals |
| Deferred items | Particle changes left in Mars1; unresolved future full assertions |

Maintain this ledger after each meaningful stage. A documentation revision is
not implementation evidence. The prior package's portable/initialization
observations remain historical; populate actual Mars2 results here.

## Mars2 execution ledger

### BG3D-0 — baseline identity and initial regression, 2026-10-03

| Field | Record |
|---|---|
| Stage/current action | `BG3D-0` QUALIFIED; `BG3D-1` is now active. Later BG3D stages remain sequential. |
| Implementation | Baseline/source audit, preservation sentinels, supported-profile matrix, selected analytical closures, and pre-implementation tolerance policy are frozen in `src/models/sep_corona_swcme/BG3D0_CONTRACT.md`. No Mars1 implementation has yet been transferred. Root instructions encode G1 -> G2 -> G3 -> G4 and the mandatory native rebuild preflight. |
| Validation | BG3D-0 qualified. Shared archives/tests pass; a compliant clean native rebuild and evolving zero-particle 1-/4-rank baseline controls pass. The initial srcSEP3D routine's only failure was repaired as a Matplotlib API compatibility defect and its focused 21-test case now passes. `V2D01` remains a declared skip because no srcSEP source was supplied; it is not required for the deliberately srcSEP3D-only BG3D-0 gate. The native controls begin at 20 solar radii and are not G1 low-corona/handoff/1-AU qualification. |
| Baseline identity | Mars2 Git `fde03997f408792a544291146c5f81f057aa32a8`; origin `git@github.com:SWMFsoftware/AMPS`; Mars1 reference has the same commit but a substantially different dirty tree. Initial Mars2 executable SHA-256 `5e872cc7e8f0c0db8ecbe34720686e95d419c2858c8b14a3c0d21f8ba40cece4`. `Makefile.conf` SHA-256 `219d843c1256bd9be2294147a712c32167c653527c344cb1daac77abff0340da`. Selected application: `sep3d`. |
| Commands/results | `make -C src/models/sep_common -j16 verify`: PASS (`SEP_COMMON01`). `make -C src/models/swcme -j16 test`: PASS 4, FAIL 0, SKIP 0, ERROR 0. `make -C src/models/sep_coronal_cme -j16 test`: PASS 219, FAIL 0, SKIP 3 (`SLM3D01`, `EVT3D01`, `XMD3D01` require independent campaigns), ERROR 0. `python3 srcSEP3D/test/run_tests.py --routine --amps-source . --output-dir test_output/cme-codex/bg3d-0-baseline-routine --rebuild`: PASS 166, FAIL 1 (`CME3D01`), SKIP 1 (`V2D01`), ERROR 0. Compliant rebuild: `pwd`; process audit; `rm -rf -- build`; `./Config.pl -application=sep3d`; `./ampsConfig.pl -input sep3d.input -no-compile`; `make -C srcSEP3D prepare-production`; `env MAKEFLAGS=-j16 make amps`; `make -C srcSEP3D audit-production-symbols`: PASS. Native source-off commands used `mpiexec -n {1,4} ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in --test-steps 2 --expect-mpi-ranks {1,4}` with the recorded JSON/artifact/output paths. Rank 1: PASS 9, FAIL 0, SKIP 1 (`SWBGAMPS03`, no remote blocks), ERROR 0. Rank 4: PASS 10, FAIL 0, SKIP 0, ERROR 0. |
| Resolved baseline failure | The original `CME3D01` renderer failure was reproduced and fixed by replacing the unsupported Matplotlib `layout=` keyword with the equivalent supported `constrained_layout=True`; model data and metrics are unchanged. `python3 srcSEP3D/test/run_tests.py --test CME3D01 --amps-source . --output-dir test_output/cme-codex/bg3d-0-cme3d01-compat --no-build --verbose`: PASS 1 registered case / 21 Python tests, FAIL/SKIP/ERROR 0. `V2D01` skipped because the baseline runner was not given an explicit srcSEP source. |
| Provenance | Coronal modular physics hashes: `physics.md` `91379645...`, `configuration_validation.md` `4cc1d51a...`, `architecture_exchange.md` `e6edf02b...`, `testing_validation.md` `02d4908f...`, `requirements.yaml` `f3dee57e...`. Composite design hash: `291c8887...`. Rebuilt `amps` SHA-256 `4c766b6e0aeeed1b85a1123b14a29e5991b4c244db02ea5b55600dff452f1b96`; `build/main/mainlib.a` `e2441a16...`; `build/main/main.a` `0ec5c448...`; MPI ranks 1/4, one OpenMP thread per rank. Evidence: `test_output/cme-codex/bg3d-0-baseline-routine`, `test_output/cme-codex/bg3d-0-native-r1`, and `test_output/cme-codex/bg3d-0-native-r4`. |
| Preservation audit | The nine particle-only files listed in `BG3D0_CONTRACT.md` are byte-identical to Git baseline; their ordered hash-manifest SHA-256 is `c52c6bdd658e2f285dd65e0e09b9c35785e0879bfc344e344db2a97bdf71d6c1`. Twenty-six mixed files were reviewed as hunk-only transfer boundaries; their ordered hash-manifest SHA-256 is `4eac9c069ce9761157ac3d1006600e3a7b1bccfde712484f6259dc2c19ac9389`. The independent Mars1 clean reference is the same commit and its recorded source archive SHA-256 is `8fa1e375...`; the archive bytes are absent, so this is recorded provenance rather than a recomputed digest. |
| Architecture/contract evidence | `make -C src/models/sep_coronal_cme -j16 check-architecture check-adapters`: PASS (`ARCHSCCM01`, five source-audit tests, external C++17 public-header consumer, both adapter compiles). The selected G1 profile is `pfss-parker-shock-fed-map-v1`: one PFSS/Parker ambient; complete ellipsoid and matched DBM history; RH-fed positive-J sheath map; distinct vector-potential/positive-J ejecta map; explicit contact/post-event disposition; frozen mass/flux/induction/interface/thermal/force/heat/work checks. Legacy `SHOCK_ONLY`, `FULL_ICME_DIAGNOSTIC`, ambient-plus-front, and spheromak/import modes are controls or unsupported, not G1. Contract SHA-256 `a6ed590cd563456cf84611740398222978baa84f61d3cbe02c568d9620cb756b`. |
| Physics status | BG3D-0 selected the required full spatial analytical closures but did not implement them. Existing code/tests cover ambient, finite front/local shocks and legacy SWCME regional fixtures. BG3D-1 must now implement owning background input/history/assets; BG3D-4/5 mass/flux/induction/thermodynamic/interface/force/heat/work evidence remains absent. G1 remains OPEN. |
| Build-process note | The baseline routine's `BLDL3D01` invoked a configured build before the user's reiterated exact clean-build sequence was recorded. A subsequent fresh rebuild followed the four-step root-verify/no-process, `rm -rf -- build`, regenerate-config/hooks, `-j16` sequence and passed. Every later native rebuild must repeat it. Configuration emitted nonfatal host/Schedule/dataARMS warnings; compiler wrappers repeatedly emitted nonfatal OpenMPI `opal_ifinit: socket() failed with errno=1` warnings. |

### BG3D-1 — owning input, assets, and continuous launch, 2026-10-03

| Field | Record |
|---|---|
| Stage/current action | `BG3D-1` QUALIFIED as a shared, background-only software stage; `BG3D-2` ambient/plasma implementation is active. This does not qualify spatial CME regions or G1. |
| Implementation | Added dependency-light `src/models/sep_corona_swcme` library. `cme_event.{h,cpp}` implements a closed configuration grammar, application-injected asset acquisition, SHA-256 byte verification, typed PFSS/Parker and regional-map inputs, immutable identity, exact cubic-Hermite between-knot positivity/outward checks, exact triaxial radial extrema, attached-then-detached validation, and one full-shape quintic-differentiated transition into sign-aware constant-wind DBM. Added public byte-checksum wrapper over the maintained coronal SHA-256 implementation. No application, MPI, source, wave, or particle type is referenced. |
| Assets/assumptions | Current frozen example under `src/models/sep_corona_swcme/examples/bg3d1`: history SHA `06f55a7...`; typed ambient asset SHA `b959c0c4...`; typed regional/ejecta asset SHA `679d7b56...`; event input SHA `e3894760...`. BG3D-4 re-identified the regional asset after its independent compatibility test selected `contact_apex_fraction=0.99`; the earlier 0.82 identity is retained as a rejected fold negative, not silently treated as prior evidence. The regional record owns shock-admission start, contact fraction, ejecta reference density/pressure, axial/poloidal fluxes, minimum J, zero added heating, and BG3D-0 force/work caps. Current event physics fingerprint `fa3dbef584ccc34511ab5f24b4d16beb1abd2eddd84b166fe09f1c83353a1876`. End state at `200000 s` has apex `1.59878e11 m`, beyond 1 AU and within 2-AU coverage. This is configuration/trajectory evidence, not a native propagation receipt. |
| Tests | `make -C src/models/sep_corona_swcme -j16 test`: `CMBGU01` PASS strict grammar/checksums/assets/capabilities/path-independent identity and corrupt/missing/unknown/unsupported negatives; `CMBGU02` PASS continuous history, exact/sampled extrema, solar attachment then detachment, transition phase/velocity derivative, no extrapolation and actual 1-AU trajectory coverage. `ELL3D01`–`ELL3D10`, each invoked through the canonical shared runner: PASS 10, FAIL/SKIP/ERROR 0 for independent ellipsoid geometry/normals/speeds/area/clipping/identity/history/basis/nesting. Full `make -C src/models/sep_coronal_cme -j16 test`: PASS 219, FAIL 0, SKIP 3 campaign prerequisites, ERROR 0. Optimized composite archive SHA `ae2d09bf...`; test binary SHA `c322f3e3...`. |
| Sanitizers/failures | Initial composite compile failed on one `-Werror` unused helper; removed. Initial link exposed a stale dependent archive rule; changed the composite makefile to invoke the coronal library build before linking. First sanitizer launch used root-relative paths from the module directory and failed before compilation; rerun from the root compiled successfully. ASan/UBSan execution with `detect_leaks=0` PASS. LeakSanitizer reports that it cannot run under the environment's ptrace setup, so leak detection is SKIP. No result is hidden or reclassified. |
| Validation boundary | BG3D-1 resolves and certifies inputs and motion only. The typed ambient asset is not yet an evaluated plasma/IMF provider; the regional asset is not yet a material map. BG3D-2 through BG3D-9, native 1-/4-rank publication, complete region coverage, residual budgets and the separate long campaign remain required. G2/G3/G4 remain deferred. |

### Current requested status checkpoint, 2026-10-03

| Question | Accurate status |
|---|---|
| Current milestone/stage | `G1`, `BG3D-4` partial and blocked at its material-contact/finite-inventory/evolution gates. BG3D-0/1/2 and the BG3D-3 shock/front subset remain qualified; the BG3D-3 shared contact/interface subset is reopened and unqualified. G1 is OPEN. BG3D-5 cannot start, and current authorization ends after G1; G2/G3/G4 must not start. |
| Implemented versus validated | BG3D-4 implements the experimental relaxing interior map, regular curved cohort reference metrics, Cauchy/adiabatic state, inventory/query APIs, independent manufactured fixtures, absolute contact-leakage measurement, transactional committed-material readback and distributed momentum/energy/work diagnostics. Those verification components pass. **The production stage is not validated:** it still lacks a justified finite initial volume and single cross-stage contact authority, distributed leakage exceeds its frozen bound, and the converged force/energy residual is persistent model error. |
| Exact accumulated relevant counts | Baseline `sep_common`: PASS 1, FAIL/SKIP/ERROR 0. Baseline SWCME top-level suite: PASS 4, FAIL/SKIP/ERROR 0. Current full maintained coronal suite: PASS 219, FAIL 0, SKIP 3 (`SLM3D01`, `EVT3D01`, `XMD3D01` campaign prerequisites), ERROR 0. Focused ellipsoid tests: PASS 10, FAIL/SKIP/ERROR 0. Current composite registered checks: PASS 6 (`CMBGU01`–`05`, `ARCHCSWC01`), FAIL 0, SKIP 0, ERROR 0. Focused Matplotlib compatibility regression: PASS 1 registered case comprising 21 Python tests, FAIL/SKIP/ERROR 0. Baseline native controls: rank 1 PASS 9, FAIL 0, SKIP 1 remote-ghost case, ERROR 0; rank 4 PASS 10, FAIL/SKIP/ERROR 0. ASan/UBSan `CMBGU05`: PASS 1, FAIL/SKIP/ERROR 0 with leak detection disabled; LeakSanitizer capability remains SKIP. Counts are reported per suite because focused cases overlap the maintained suite. |
| Remaining failures/blockers | The blocking physical prerequisite is a closed replacement construction with one compatible shock/contact history and finite initial sheath state. The current fixed-fraction geometry and oldest-cohort surface disagree by about `1.0e-2`. Independently differentiated leakage estimates exceed the frozen budget but are nonconvergent with comparable derivative uncertainty, so this is a numerical verification failure—not proof of a resolved physical leak and not a qualification from the API's identical-data zero. By contrast, the production force ratio about `0.646` and absolute energy residual about `1.26e16 W` are stable over the tested time refinement and are persistent model/approximation residuals. Sub-fast continuation is deliberately unsupported and rejects transactionally; no compression-wave state is claimed. LeakSanitizer remains SKIP under ptrace, and maintained campaign tests `SLM3D01`, `EVT3D01`, `XMD3D01` remain SKIP. |
| Spatial sheath/ejecta | **Sheath: partially implemented and internally tested, not qualified. Ejecta: not implemented or tested.** The spatial inverse and material identities pass, but the nonzero rear-boundary mass transfer fails the selected material-contact requirement. The vector-potential ejecta inputs remain frozen only. |
| Shock diagnostics | **Implemented and standalone-tested for BG3D-3.** At the reference epoch, all `1152` patches retain geometry/upstream state; `576` are fast and `576` sub-fast. Fast patches retain both-sided states, Mach number, obliquity and compression; maximum normalized RH residual is `9.93759e-11`, minimum compression is `1.46728`. These interface states do not supply the missing downstream sheath volume. |
| Continuous SWCME handoff | Full-shape kinematic transition/DBM continuation is implemented and standalone-tested in BG3D-1, including a numerical derivative oracle. Material/regional state continuity through handoff and native publication are **not** implemented or tested, so BG3D-6 and G1 remain open. |
| Actual propagation to 1 AU | The frozen analytical trajectory evaluates past 1 AU (`1.59878e11 m` at `200000 s`) and ambient samples/flux are tested at 1.0/1.1 AU. There is **no actual evolving native low-corona-to-1-AU spatial-CME run** yet. The existing 1-/4-rank controls start at 20 solar radii and run only two 60-s steps; they do not qualify G1. |
| Preservation | All nine particle-only file hashes still match their frozen baseline values; ordered manifest sentinel remains `c52c6bdd...`. The mixed shock-provider addition defaults to legacy source-enabled behavior, and the full maintained regression passes, so baseline particle functionality is preserved. Mars1 remains read-only at commit `fde03997...`; no command modified it, and its tracked-only (`-uno`) status digest remains `512baf7f...`. |
| Now/next | Remain in BG3D-4. Resolve the missing physical contact-history/initial-inventory prerequisite, then add the no-through-flow/normal-velocity interface test and requalify the complete stage. BG3D-5, native work, srcSEP and particles remain forbidden. |

### BG3D-2 — ambient and low-coronal plasma, qualified 2026-10-03

| Field | Record |
|---|---|
| Implementation | Added `ambient_state.{h,cpp}` in the shared composite library, selectively adapting the independently audited Mars1 ambient prototype to the particle-free event. It composes maintained PFSS, topology tracing, closed hydrostatic plasma, species EOS, transonic isothermal Parker wind, angular-momentum flow, and retarded Parker IMF kernels. Samples carry rho/U/p/T/B, region, signed sector, event identity, epoch/generation, gradients, stencil category and readiness. Magnetic nulls, branch-mixing derivatives, invalid epochs and radial over/under-coverage fail with typed status. A one-sided PFSS source-surface coordinate guard repairs floating-point shell overshoot without flooring B. |
| Validation | QUALIFIED for BG3D-2. `make -C src/models/sep_corona_swcme -j16 test`: PASS 4 (`CMBGU01`–`03`, `ARCHCSWC01`), FAIL/SKIP/ERROR 0. `CMBGU03` independently checks PFSS open/closed and exterior branches, north/south sector, corotation, electron EOS, Parker winding/angular momentum and isothermal momentum equations, exact radial mass-flux conservation, current-sheet-null/coverage/generation negatives, a separate five-point derivative oracle, normalized exterior divergence, and signed flux at 12x24, 24x48 and 48x96. Flux residuals: `6.11174e-17`, `5.62905e-16`, `2.66055e-15`; mass flux `9.18831e7 kg s^-1 sr^-1`. `ARCHCSWC01` passes source scan, external public-header compile and archive symbol audit. ASan/UBSan `CMBGU03` passes with leak detection disabled; the already recorded ptrace restriction leaves LeakSanitizer SKIP. Full maintained coronal regression remains PASS 219 / SKIP 3 / FAIL 0 / ERROR 0. Composite archive SHA `98e00500...`; BG3D-2 test binary SHA `3ed77fa0...`; ambient header/source SHA `7112ae3f...` / `83d9ee0d...`. |
| Preserved failures during development | Initial Parker-winding assertion used a unit-azimuth formula but compared the code's spherical component without the required `sin(theta)` factor; the independent oracle was corrected. The five-point oracle then exposed a real source-surface roundoff defect (`PFSS query is outside its spherical shell`); production evaluation now uses the explicit one-sided interior coordinate limit. A monotonic roundoff assertion on already `1e-15` signed-flux residuals was replaced with the predeclared `1e-12` absolute normalized target at all three resolutions; raw values remain recorded. |
| Exit/next boundary | BG3D-2 supplies ambient only. It does not provide a front, shock, sheath, ejecta or native state. BG3D-3 must bind complete event geometry and local both-sided shocks to this exact ambient identity. No native rebuild was needed for BG3D-2. |

### BG3D-3 — shock/front qualified; shared contact reopened 2026-10-04

| Field | Record |
|---|---|
| Stage/current action | The BG3D-3 front/local-shock subset remains QUALIFIED. Its fixed-fraction contact/interface subset is REOPENED and UNQUALIFIED because BG3D-4 uses a different material surface. `BG3D-4` remains active. G1 is OPEN, and no application or particle stage may start. |
| Implementation | Added `surface_shock.{h,cpp}` to the shared composite library. One transactional epoch owns the qualified event and ambient, complete tessellated front, rear-aligned nested contact/ejecta boundary, patch positions/normals/areas/normal speeds/upstream samples, canonical shock snapshot, generation and event identity. The contact span and all axes are derived from the same front rather than a second radius/apex. Maintained ellipsoid physical IDs now start at one. The maintained shock provider accepts the event's gamma and a source-enable flag whose default preserves legacy callers; this background explicitly disables source eligibility and all particle measures. Failed preparation preserves the committed epoch. |
| Validation | QUALIFIED for BG3D-3. `make -C src/models/sep_corona_swcme -j16 test`: PASS 5 (`CMBGU01`–`04`, `ARCHCSWC01`), FAIL 0, SKIP 0, ERROR 0. `CMBGU04` checks complete and distinct nested geometry, stable nonzero IDs, finite-difference normal velocity, three-level surface-area convergence, direct upstream binding, fast/sub-fast separation, both-sided RH state, Mach/compression/obliquity, zero source eligibility/rates/measures, and failed-candidate transaction integrity. At `10000 s`: `1152` patches, `576` fast, `576` sub-fast, maximum RH residual `9.93759e-11`, minimum compression `1.46728`; absolute area errors versus the `96x192` reference are `2.94955e17`, `7.00027e16`, `1.39894e16 m2` at `12x24`, `24x48`, `48x96`. A separate slow-event fixture remains wholly sub-fast. ASan/UBSan `CMBGU04` passes with `ASAN_OPTIONS=detect_leaks=0`; LeakSanitizer remains the recorded ptrace-environment SKIP. Full maintained coronal regression after the mixed-file compatibility additions: PASS 219, FAIL 0, SKIP 3, ERROR 0. |
| Commands/provenance | Shared tests: `make -C src/models/sep_corona_swcme -j16 test`. Sanitizer compile: `g++ -I src/models/sep_corona_swcme/include -I src/models/sep_coronal_cme/include -I src/models/sep_common -std=c++17 -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined src/models/sep_corona_swcme/src/cme_event.cpp src/models/sep_corona_swcme/src/ambient_state.cpp src/models/sep_corona_swcme/src/surface_shock.cpp src/models/sep_corona_swcme/test/test_bg3d3.cpp src/models/sep_coronal_cme/build/libsep_coronal_cme.a -o /tmp/sep-corona-swcme-bg3d3-sanitized`; execution `env ASAN_OPTIONS=detect_leaks=0 /tmp/sep-corona-swcme-bg3d3-sanitized`. Compatibility: `make -C src/models/sep_coronal_cme -j16 test`. Composite archive SHA `4b414daa...`; BG3D-3 binary `dcbc9859...`; surface header/source/test `a4fc2693...` / `e7ee7628...` / `0c4a9dd9...`; mixed provider header/source `300e21a2...` / `68ebba35...`. No native rebuild was performed. |
| Physics boundary | The local discontinuity is now real analytical state and diagnostics, but it is not a volume closure. No spatial sheath/ejecta plasma or IMF has been claimed. BG3D-4 must construct admitted downstream sheath material, initialize its inventory, evolve its spatial map, define every exit/disposition, and independently check planar/expanding limits, mass and interface residuals. BG3D-5 must provide the distinct ejecta magnetic/thermal/material closure and force/heat/work budgets. |

### BG3D-4 — shock-fed spatial sheath, partial/blocked 2026-10-03

| Field | Record |
|---|---|
| Stage/current action | `BG3D-4` remains active but is BLOCKED at the required material-contact interface. `BG3D-5` has not started. G1 remains OPEN; G2/G3/G4 remain unauthorized. |
| Implementation | Added `sheath_model.{h,cpp}`. Material labels are complete front coordinates plus unique shock-crossing time. The prescribed ballistic map exactly recovers canonical `U2` at birth; its crossing/current derivative bases form `F` and positive `J`, with `rho=rho2/J`, `p=p2 J^-gamma`, `B=F B2/J`, and analytic material velocity. Admission uses both equal RH mass fluxes and actual birth area. First contact/photosphere/outer exits are event-located and ledgered before any later fold. The committed inventory has unique patch/time cells, explicit fast/sub-fast area-time, mass dispositions and transaction integrity. A damped Newton inverse uses `dx/d(theta,phi,tau)` to return actual Cartesian spatial state; unsupported points fail instead of receiving ambient. Fourth-order same-branch deficit derivatives use explicit one-sided fallback only at physical support boundaries. |
| Validation | PARTIAL, not qualified. `make -C src/models/sep_corona_swcme -j16 test`: component PASS 6 (`CMBGU01`–`05`, `ARCHCSWC01`), FAIL 0, SKIP 0, ERROR 0. `CMBGU05` checks the exact RH shock limit, Cauchy determinant/field and adiabatic identities, unique cells, geometric-exit closure, off-grid Cartesian inversion and unsupported-point rejection. Admission-time quadrature errors against 16 bins are `6.63936e7`, `1.58266e7`, `3.16624e6 kg`; admitted mass is `6.46441e10 kg`. Minimum J is `0.859165`; maximum admission/inventory residuals are `1.0242e-15` / `3.9022e-16`. Independent Lagrangian mass/thermal/induction residuals are `2.99574e-12`, `3.67904e-12`, `1.60474e-11`; three material-time stencils give maxima `2.04692e-12`, `1.97002e-12`, `5.97736e-12`. However, the 200-s fixture transfers `2.08416e11 kg` through the nominal contact. This is a required physical interface failure even though the permeable-boundary ledger closes. The old 0.82 contact also folds before exit. Maintained exact moving-planar sheath and RH references remain PASS 219 / SKIP 3 / FAIL 0 / ERROR 0. ASan/UBSan passes with leak detection disabled; LeakSanitizer remains the ptrace SKIP. |
| Commands/provenance | Optimized/architecture command: `make -C src/models/sep_corona_swcme -j16 test`. Sanitizer compiles all four composite sources plus `test_bg3d4.cpp` with `-fsanitize=address,undefined` to `/tmp/sep-corona-swcme-bg3d4-sanitized`, then runs `env ASAN_OPTIONS=detect_leaks=0 ...`: PASS. Compatibility: `make -C src/models/sep_coronal_cme -j16 test`: PASS 219, SKIP 3. Composite archive SHA `80d4981c...`; BG3D-4 binary `259c9fef...`; sheath header/source/test `164b4716...` / `be11dda7...` / `822a65f5...`; regional asset `679d7b56...`; event input `e3894760...`; event fingerprint `fa3dbef5...`. No native AMPS rebuild was performed. |
| Preserved failures/limits | Contact fraction 0.82 fails the minimum-J gate before the rear boundary. Fraction 0.99 avoids the fold by allowing early rear-boundary transfer, but that transfer violates the selected material-contact/no-through-flow contract; it therefore does not rescue qualification. The missing contact history and initial inventory are physical inputs/state, not a numerical tolerance. BG3D-4 also does not claim ejecta, sub-fast compression-layer state, post-event recovery, force/work balance, native storage or 1-AU native propagation. |

### BG3D-4 authoritative correction intake, 2026-10-04

| Field | Record |
|---|---|
| Authority | Read the user-supplied review at its actual Mars2 location, `src/models/sep_coronal_cme/docs/BG3D4_CODEX_REVIEW_AND_CORRECTIONS.md` (the requested `src/models/sep_corona_swcme/docs/...` path does not exist). It now governs BG3D-4 over older shorthand in this plan, `BG3D0_CONTRACT.md`, `model.md`, code comments and tests. Root instructions and the contract/model references were updated before implementation. |
| Reproduction | From `/home/vtenishe/Mars2/AMPS`: `make -C src/models/sep_corona_swcme -j16 build/test_bg3d4 && cd src/models/sep_corona_swcme && ./build/test_bg3d4`. Executed current source: executable PASS 1, FAIL/SKIP/ERROR 0, but reproduced the required physical failure `permeable_contact_mass_kg=2.08416e11`, `material_contact_qualified=false`. This is diagnostic/component success, not BG3D-4 qualification. No native AMPS rebuild occurred. |
| Corrections accepted | Curved volume must use the full physical label-map metric; `rho d\ell` is not an exact curved measure. Cauchy/adiabatic identities require a regular three-dimensional reference map with `F_rel=A A0^-1`, not an implicit identity in angle/time labels. Contact leakage uses pointwise dimensional-plus-relative tolerances and accumulated absolute, not only signed, transfer. Zero mass/admission limits use absolute terms and no floors. Zero-volume startup is a declared limiting fixture only; inventory may not reset at each local first-fast classification. Each boundary uses its own point/normal, and cross-layer comparisons require a declared correspondence. |
| Active implementation boundary | Remove the fixed-fraction rear surface as a purported material contact and replace it with a fully defined compatible contact/startup construction. Add exact planar and independent expanding-curved oracles, curved-volume negative controls, reference-map/flux compatibility, pointwise and absolute contact leakage, zero-admission/startup behavior, global injectivity/coverage, smooth mass/induction/divergence/thermal/momentum diagnostics, and separately refined time/surface/material/timestep evidence. A velocity interpolation, closed mass ledger or positive local Jacobian alone cannot qualify the stage. |

### BG3D-4 bounded verification/diagnostics increment, criteria frozen 2026-10-04

| Field | Record |
|---|---|
| Cross-stage status | Reopen only the BG3D-3 contact/interface subset. The independently checked front geometry, normal motion, shock admission and both-sided RH diagnostics remain qualified. The fixed-fraction surface is retained only as an explicitly unqualified ejecta-reference geometry until a single contact authority replaces it. |
| Manufactured fixtures | Add (1) a finite initial planar slab between a material piston/contact and a constant-speed admissible MHD shock, with exact initial inventory plus subsequent admission, and (2) a finite curved spherical-shell material map under force-free homologous expansion. Each fixture owns explicit initial/boundary states and nonzero divergence-free magnetic topology. They are verification solutions, not physical CME qualification. |
| Numerical criteria frozen before implementation | Exact algebraic RH/map/state/energy identities: relative error `<=1e-12` plus explicit expected-zero absolute scales. Independent fourth-order derivatives: finest normalized error `<=1e-7` with three-level convergence until roundoff. Curved midpoint volume/mass quadrature: three decreasing errors with observed order at least `1.8` and finest relative error `<=1e-3`. Contact position agreement for a future shared authority will require `1e-10` relative geometry and the already frozen nonsingular flux/velocity budgets; the present cross-authority mismatch is expected to fail and must be reported, not tuned away. |
| Production diagnostics | Sample the actual `ShockFedSheathModel::Evaluate` material map at distributed supported labels. Report the dimensional momentum residual field, volume-integrated absolute residual and force scale, volume-weighted local P99, signed and absolute residual work, differential total-energy residual and resolved inertia/pressure/Lorentz/gravity work. Three derivative step sizes separate numerical convergence from the extrapolated physical/model residual. The existing event caps (`2`, `5`, `2`) are reporting limits only in this increment and cannot qualify the experimental evolution. |
| Admission/inventory failure policy | Shock admission and material queries are distinct. A sub-fast or otherwise unsupported candidate epoch is rejected transactionally; it is not relabeled as a solved compression wave. The last committed epoch, inventory handle and stored material cells remain queryable. No later compression-region capability is claimed. |
| Production-law boundary | This increment does not replace or qualify `rh-relaxing-material-map-v1`. A replacement requires a separate closed construction specifying its unknowns, equations, compatible finite initial state, shock/contact boundary data and solution procedure. Running the current empirical relaxation for an earlier interval is not accepted as prehistory. |
| Pre-change reproduction | `make -C src/models/sep_corona_swcme -j16 build/test_bg3d4`; `(cd src/models/sep_corona_swcme && ./build/test_bg3d4)`: PASS 1, FAIL 0, SKIP 0, ERROR 0 in approximately 154 s. Recorded `momentum_pa_m=2.75827e-13`, single-point `momentum_ratio=0.646356`; no distributed force or energy/work qualification existed. No native AMPS rebuild was performed. |
| Implementation increment | Added independent finite-inventory planar shock/piston and homologous curved-shell fixtures with explicit plasma/IMF states, maps, boundary conditions, topology and exact moving-boundary energy budgets. Added reusable public-map diagnostics for distributed momentum terms, conservative total-energy residual, resolved work and volume-weighted integration. Surface epochs now type the fixed-fraction contact as unqualified. Added committed-cell readback that does not re-run a rejected later shock admission. The empirical production evolution equation is unchanged and remains explicitly unqualified. |
| Focused post-change evidence | `make -C src/models/sep_corona_swcme -j16 build/test_bg3d4 && (cd src/models/sep_corona_swcme && ./build/test_bg3d4)`: PASS 1, FAIL 0, SKIP 0, ERROR 0. Manufactured curved volume errors `1.22956e-3,3.07164e-4,7.67770e-5` have orders `2.00106,2.00026`; energy/work errors `2.73638e-8,1.70983e-9,1.06838e-10`. Cross-stage contact mismatch is `1.00149e-2`, so the shared contract remains unqualified. Independently differentiated absolute leakage is `3031.26,1004.33,2645.37 kg` versus allowed `196.047,99.1243,170.685 kg`; every level is explicitly unqualified and the identically zero API self-report is shown separately. A sub-fast `60 s` candidate rejects while the committed `20 s` inventory and all `128` stored cells remain available. Eight production samples give stable integrated force ratio `0.645573`, local P99 `0.646356`, absolute energy residual `1.25746e16 W`, and force-work ratio `0.279886` at all three time stencils. This stability is evidence of persistent approximation/model residual, not disappearing temporal error. |
| Regression/sanitizer evidence | `make -C src/models/sep_corona_swcme -j16 test`: PASS 6 (`CMBGU01`–`05`, `ARCHCSWC01`), FAIL 0, SKIP 0, ERROR 0. `make -C src/models/sep_coronal_cme -j16 test`: PASS 219, FAIL 0, SKIP 3 (`SLM3D01`, `EVT3D01`, `XMD3D01` require campaign evidence), ERROR 0. ASan/UBSan compilation of all five composite sources plus `test_bg3d4.cpp`, followed by `ASAN_OPTIONS=detect_leaks=0`, is PASS 1, FAIL/SKIP/ERROR 0; LeakSanitizer itself remains the recorded ptrace-environment SKIP. No native AMPS rebuild was needed or performed. |
| Provenance | Mars2 HEAD at this increment is user-created commit `8b771a27c43101b3da0e4f8a9e0a50973a2e4842`; no commit/push was made by Codex. Composite archive SHA-256 `4c6fa79c41d1b02aa034979db709f06929c6e70be3dda9a384fe72190367922d`; focused binary `c832f2af91a209933045b5a9093d2349307d1252fab7c021600aae8898a37af0`. New diagnostics header/source hashes `db44464e...` / `20798a68...`; independent fixture `078d625b...`; focused test `d605ed99...`. All nine particle-only preservation hashes still exactly match the frozen table. Mars1 remains at `fde03997...` with unchanged tracked-status digest `512baf7f...`. |
| Remaining boundary | The fixtures and diagnostic machinery are validated components; the production map is not. No justified finite production initial inventory, single contact history, compression-region continuation, or closed replacement evolution has been implemented. The broad configured force/work caps are not used to qualify this increment. BG3D-4 and G1 remain OPEN. |

## Completion delivery

Return focused direct source changes, detailed comments/READMEs, documented
inputs, independent tests, complete logs/native receipts and scientific
residual/convergence figures. Include an inventory explaining every retained
change and a baseline particle/core compatibility report. Do not include
generated objects/dependencies/binaries/output in source archives. State which
gate is qualified and which later gates remain deferred. No commit/push.
